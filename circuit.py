#!/usr/bin/env python3
"""
circuit — un petit simulateur de circuits électriques façon SPICE, écrit de zéro.

Il lit une « netlist » (la liste des composants et de leurs nœuds), écrit les
équations de Kirchhoff sous forme matricielle (analyse nodale modifiée, MNA)
et les résout par élimination de Gauss. Pas de numpy : tout est fait main.

Composants : R, C, L, sources V et I (continu, SIN, PULSE, PWL, AC),
             diode D (non linéaire : Newton-Raphson), source commandée E
             (pratique pour un amplificateur opérationnel idéal).
Analyses   : .op (point de repos), .dc (balayage), .tran (temporel, règle
             des trapèzes), .ac (réponse en fréquence, nombres complexes).
"""
import cmath
import math
import os
import re
import sys

VT = 0.025852  # tension thermique kT/q à 27 °C


# =================================================================== valeurs

SUFFIXES = {"f": 1e-15, "p": 1e-12, "n": 1e-9, "u": 1e-6, "µ": 1e-6, "m": 1e-3,
            "k": 1e3, "meg": 1e6, "g": 1e9, "t": 1e12}


def valeur(texte):
    """'4.7k' -> 4700, '100n' -> 1e-7, '10meg' -> 1e7, '5V' -> 5 (comme SPICE, m = milli)."""
    m = re.fullmatch(r"([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)(meg|[fpnuµmkgt])?[a-zω]*", texte.strip().lower())
    if not m:
        raise ValueError(f"valeur illisible : {texte!r}")
    return float(m.group(1)) * SUFFIXES.get(m.group(2), 1.0)


def joli(x, unite=""):
    """4700 -> '4.7 k', 1e-7 -> '100 n'."""
    if x == 0 or not math.isfinite(x):
        return f"{x:g} {unite}".strip()
    for pref, f in (("T", 1e12), ("G", 1e9), ("M", 1e6), ("k", 1e3), ("", 1), ("m", 1e-3), ("µ", 1e-6), ("n", 1e-9), ("p", 1e-12)):
        if abs(x) >= f * 0.9995:
            return f"{x / f:.4g} {pref}{unite}".strip()
    return f"{x:.3g} {unite}".strip()


# ============================================================ algèbre linéaire

def resous(A, b):
    """Élimination de Gauss avec pivot partiel. Marche aussi en complexe."""
    n = len(b)
    A = [ligne[:] for ligne in A]
    b = b[:]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(A[i][k]))
        if abs(A[p][k]) < 1e-300:
            raise ArithmeticError("système singulier : un nœud est peut-être en l'air (relié à rien d'autre)")
        if p != k:
            A[k], A[p] = A[p], A[k]
            b[k], b[p] = b[p], b[k]
        piv = A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] / piv
            if f:
                Ai, Ak = A[i], A[k]
                for j in range(k, n):
                    Ai[j] -= f * Ak[j]
                b[i] -= f * b[k]
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / A[i][i]
    return x


# =================================================================== sources

class Source:
    """Valeur d'une source V ou I en fonction du temps, plus sa part « AC »."""

    def __init__(self, texte):
        t = texte.lower()
        self.dc, self.ac, self.phase, self.forme, self.args = 0.0, 0.0, 0.0, None, []
        m = re.search(r"\b(sin|pulse|pwl)\s*\(([^)]*)\)", t)
        if m:
            self.forme, self.args = m.group(1), [valeur(v) for v in m.group(2).replace(",", " ").split()]
            t = t[:m.start()] + t[m.end():]
        m = re.search(r"\bac\s+(\S+)(?:\s+(\S+))?", t)
        if m:
            self.ac = valeur(m.group(1))
            if m.group(2):
                try:
                    self.phase = valeur(m.group(2))
                except ValueError:
                    pass
            t = t[:m.start()] + t[m.end():]
        m = re.search(r"\bdc\s+(\S+)", t)
        if m:
            self.dc = valeur(m.group(1))
        elif t.split():
            self.dc = valeur(t.split()[0])
        if self.forme == "sin" and len(self.args) < 3:
            raise ValueError("SIN(décalage amplitude fréquence [retard])")
        if self.forme == "pulse" and len(self.args) < 7:
            raise ValueError("PULSE(v1 v2 retard montée descente largeur période)")

    def __call__(self, t):
        a = self.args
        if self.forme == "sin":
            vo, va, f = a[0], a[1], a[2]
            td = a[3] if len(a) > 3 else 0.0
            return vo if t < td else vo + va * math.sin(2 * math.pi * f * (t - td))
        if self.forme == "pulse":
            v1, v2, td, tr, tf, pw, per = a[:7]
            if t < td:
                return v1
            t = (t - td) % per if per > 0 else t - td
            if t < tr:
                return v1 + (v2 - v1) * t / tr if tr > 0 else v2
            if t < tr + pw:
                return v2
            if t < tr + pw + tf:
                return v2 + (v1 - v2) * (t - tr - pw) / tf if tf > 0 else v1
            return v1
        if self.forme == "pwl":
            pts = list(zip(a[::2], a[1::2]))
            if t <= pts[0][0]:
                return pts[0][1]
            for (t0, v0), (t1, v1) in zip(pts, pts[1:]):
                if t <= t1:
                    return v0 + (v1 - v0) * (t - t0) / (t1 - t0)
            return pts[-1][1]
        return self.dc

    def valeur_repos(self):
        return self(0.0) if self.forme else self.dc


# ================================================================== netlist

class Composant:
    def __init__(self, nom, noeuds, **kw):
        self.nom, self.noeuds = nom, noeuds
        self.__dict__.update(kw)

    @property
    def type(self):
        return self.nom[0].upper()


class Circuit:
    def __init__(self, texte, chemin="?"):
        self.titre, self.composants, self.analyses, self.sondes = None, [], [], []
        self.noeuds = {"0": -1}
        lignes = []
        for brute in texte.splitlines():
            l = brute.split(";")[0].rstrip()
            if l.startswith("+") and lignes:
                lignes[-1] = (lignes[-1][0], lignes[-1][1] + " " + l[1:])
            else:
                lignes.append((len(lignes) + 1, l))
        for num, l in lignes:
            l = l.strip()
            if not l:
                continue
            if l.startswith("*"):
                if self.titre is None:
                    self.titre = l.lstrip("* ").strip()
                continue
            try:
                self._ligne(l)
            except (ValueError, IndexError) as e:
                raise SystemExit(f"{chemin}, ligne {num} : {e}\n  {l}")
        self.n_noeuds = len(self.noeuds) - 1
        # une inconnue de courant par source de tension, source commandée et bobine
        k = self.n_noeuds
        for c in self.composants:
            if c.type in "VEL":
                c.k, k = k, k + 1
        self.taille = k

    def noeud(self, nom):
        nom = nom.lower()
        if nom in ("gnd", "masse"):
            nom = "0"
        if nom not in self.noeuds:
            self.noeuds[nom] = len(self.noeuds) - 1
        return self.noeuds[nom]

    def _ligne(self, l):
        mots = l.split()
        t = mots[0][0].upper()
        if mots[0].startswith("."):
            cmd = mots[0].lower()
            if cmd in (".op", ".tran", ".ac", ".dc"):
                self.analyses.append((cmd, mots[1:]))
            elif cmd in (".sonde", ".print", ".probe", ".plot"):
                self.sondes += [m for m in mots[1:] if re.match(r"[vi]\(", m, re.I)]
            elif cmd in (".end", ".fin"):
                pass
            else:
                raise ValueError(f"commande inconnue {mots[0]}")
            return
        nom = mots[0]
        if t in "RCL":
            a, b, v = self.noeud(mots[1]), self.noeud(mots[2]), valeur(mots[3])
            if v <= 0:
                raise ValueError("une résistance, un condensateur ou une bobine doit être > 0")
            self.composants.append(Composant(nom, (a, b), val=v))
        elif t in "VI":
            self.composants.append(Composant(nom, (self.noeud(mots[1]), self.noeud(mots[2])), src=Source(" ".join(mots[3:]))))
        elif t == "D":
            params = {k.lower(): valeur(v) for k, v in (p.split("=") for p in mots[3:] if "=" in p)}
            self.composants.append(Composant(nom, (self.noeud(mots[1]), self.noeud(mots[2])),
                                             Is=params.get("is", 1e-14), n=params.get("n", 1.0)))
        elif t == "E":
            self.composants.append(Composant(nom, tuple(self.noeud(m) for m in mots[1:5]), gain=valeur(mots[5])))
        else:
            raise ValueError(f"composant inconnu {nom} (R, C, L, V, I, D, E)")

    def nom_noeud(self, i):
        return next(n for n, j in self.noeuds.items() if j == i)


# ============================================================== mise en équation

class Systeme:
    """Matrice MNA en construction : on y « tamponne » chaque composant."""

    def __init__(self, n, complexe=False):
        z = 0j if complexe else 0.0
        self.A = [[z] * n for _ in range(n)]
        self.b = [z] * n

    def g(self, a, b, g):  # conductance entre a et b
        A = self.A
        if a >= 0:
            A[a][a] += g
        if b >= 0:
            A[b][b] += g
        if a >= 0 and b >= 0:
            A[a][b] -= g
            A[b][a] -= g

    def i(self, a, b, i):  # courant i qui va de a vers b à travers le composant
        if a >= 0:
            self.b[a] -= i
        if b >= 0:
            self.b[b] += i

    def v(self, a, b, k, v):  # source de tension : V(a) - V(b) = v, courant inconnu n° k
        A = self.A
        if a >= 0:
            A[a][k] += 1
            A[k][a] += 1
        if b >= 0:
            A[b][k] -= 1
            A[k][b] -= 1
        self.b[k] += v


def tension(x, a):
    return x[a] if a >= 0 else 0.0


def limite_diode(v, vieux, d):
    """Limitation de pas de SPICE (pnjlim) : évite que exp() explose quand
    Newton propose une tension de diode beaucoup trop grande."""
    nvt = d.n * VT
    vcrit = nvt * math.log(nvt / (math.sqrt(2) * d.Is))
    if v > vcrit and abs(v - vieux) > 2 * nvt:
        if vieux > 0:
            arg = 1 + (v - vieux) / nvt
            v = vieux + nvt * math.log(arg) if arg > 0 else vcrit
        else:
            v = nvt * math.log(v / nvt)
    return v


def diode(d, v):
    """Courant et conductance de la diode (loi de Shockley) à la tension v."""
    nvt = d.n * VT
    e = math.exp(min(v / nvt, 200))
    return d.Is * (e - 1), d.Is * e / nvt


def resout_non_lineaire(ckt, tamponne, x0=None, iterations=200):
    """Newton-Raphson : on remplace chaque diode par sa tangente au point
    courant (une conductance et une source de courant), on résout le système
    linéaire, on recommence jusqu'à ce que ça ne bouge plus."""
    diodes = [c for c in ckt.composants if c.type == "D"]
    x = x0[:] if x0 else [0.0] * ckt.taille
    vd = {d.nom: tension(x, d.noeuds[0]) - tension(x, d.noeuds[1]) for d in diodes}
    for it in range(iterations):
        S = Systeme(ckt.taille)
        tamponne(S)
        for d in diodes:
            idd, gd = diode(d, vd[d.nom])
            S.g(d.noeuds[0], d.noeuds[1], gd + 1e-12)  # + gmin : une vraie petite fuite en parallèle,
            S.i(d.noeuds[0], d.noeuds[1], idd - gd * vd[d.nom])  # comme dans SPICE
        nouveau = resous(S.A, S.b)
        # tolérances de SPICE : relative (reltol = 1e-3) + absolue (1 µV)
        converge = all(abs(a - b) <= 1e-3 * max(abs(a), abs(b)) + 1e-6 for a, b in zip(nouveau, x))
        x = nouveau
        if not diodes:
            return x
        for d in diodes:
            v = tension(x, d.noeuds[0]) - tension(x, d.noeuds[1])
            v2 = limite_diode(v, vd[d.nom], d)
            converge &= abs(v2 - vd[d.nom]) <= 1e-3 * abs(v2) + 1e-6
            vd[d.nom] = v2
        if converge:
            return x
    raise ArithmeticError("Newton ne converge pas")


def tamponne_lineaire(ckt, S, t=None, valeurs=None):
    """Composants linéaires en continu (t=None : point de repos)."""
    for c in ckt.composants:
        a, b = c.noeuds[0], c.noeuds[1]
        if c.type == "R":
            S.g(a, b, 1 / c.val)
        elif c.type == "L":
            S.v(a, b, c.k, 0.0)  # en continu une bobine est un fil
        elif c.type in "VI":
            v = (valeurs or {}).get(c.nom, c.src.valeur_repos() if t is None else c.src(t))
            if c.type == "V":
                S.v(a, b, c.k, v)
            else:
                S.i(a, b, v)
        elif c.type == "E":
            S.v(a, b, c.k, 0.0)
            cp, cn = c.noeuds[2], c.noeuds[3]
            if cp >= 0:
                S.A[c.k][cp] -= c.gain
            if cn >= 0:
                S.A[c.k][cn] += c.gain
        # C : circuit ouvert en continu


def point_de_repos(ckt, valeurs=None, x0=None):
    return resout_non_lineaire(ckt, lambda S: tamponne_lineaire(ckt, S, None, valeurs), x0)


# =================================================================== analyses

def mesure(ckt, x, sonde):
    """Valeur d'une sonde : V(n), V(a,b) ou I(V1)."""
    m = re.fullmatch(r"([vi])\(([^,)]+)(?:,([^)]+))?\)", sonde.strip().lower())
    if not m:
        raise SystemExit(f"sonde illisible : {sonde}")
    if m.group(1) == "v":
        a = ckt.noeuds.get(m.group(2))
        b = ckt.noeuds.get(m.group(3) or "0")
        if a is None or b is None:
            raise SystemExit(f"nœud inconnu dans {sonde}")
        return tension(x, a) - tension(x, b)
    for c in ckt.composants:
        if c.nom.lower() == m.group(2) and hasattr(c, "k"):
            return x[c.k]
    raise SystemExit(f"{sonde} : on ne peut mesurer que le courant d'une source V, E ou d'une bobine L")


def sondes(ckt):
    return ckt.sondes or [f"V({n})" for n in ckt.noeuds if n != "0"]


def analyse_op(ckt):
    x = point_de_repos(ckt)
    print("point de repos :")
    for nom, i in ckt.noeuds.items():
        if i >= 0:
            print(f"  V({nom}) = {joli(x[i], 'V')}")
    for c in ckt.composants:
        if c.type in "VEL":
            print(f"  I({c.nom}) = {joli(x[c.k], 'A')}")
        elif c.type == "D":
            v = tension(x, c.noeuds[0]) - tension(x, c.noeuds[1])
            print(f"  diode {c.nom} : {joli(v, 'V')}, {joli(diode(c, v)[0], 'A')}")
    return None


def analyse_dc(ckt, args):
    nom, debut, fin, pas = args[0], valeur(args[1]), valeur(args[2]), valeur(args[3])
    if not any(c.nom.lower() == nom.lower() and c.type in "VI" for c in ckt.composants):
        raise SystemExit(f".dc : pas de source {nom}")
    nom = next(c.nom for c in ckt.composants if c.nom.lower() == nom.lower())
    xs, series, x = [], {s: [] for s in sondes(ckt)}, None
    n = int(round((fin - debut) / pas))
    for k in range(n + 1):
        v = debut + k * pas
        x = point_de_repos(ckt, {nom: v}, x)  # on repart de la solution précédente
        xs.append(v)
        for s in series:
            series[s].append(mesure(ckt, x, s))
    return {"titre": f"balayage de {nom}", "x": xs, "xlabel": f"{nom} (V)", "series": series, "log": False}


def analyse_tran(ckt, args):
    """Intégration temporelle. Chaque condensateur et chaque bobine est remplacé,
    à chaque pas h, par son « modèle compagnon » : une conductance et une
    source de courant qui résument le passé. Règle des trapèzes (précise,
    d'ordre 2), sauf au premier pas où Euler implicite amortit le démarrage."""
    h, tfin = valeur(args[0]), valeur(args[1])
    x = point_de_repos(ckt, {c.nom: c.src(0.0) for c in ckt.composants if c.type in "VI"})
    etat = {}
    for c in ckt.composants:
        if c.type == "C":
            etat[c.nom] = [tension(x, c.noeuds[0]) - tension(x, c.noeuds[1]), 0.0]  # v, i
        elif c.type == "L":
            etat[c.nom] = [0.0, x[c.k]]
    # en temporel les bobines n'ont plus besoin d'inconnue de courant : on la garde à 0
    ts, series = [0.0], {s: [mesure(ckt, x, s)] for s in sondes(ckt)}
    n = int(round(tfin / h))
    for pas in range(1, n + 1):
        t = pas * h
        trapeze = pas > 1

        def tamponne(S):
            for c in ckt.composants:
                a, b = c.noeuds[0], c.noeuds[1]
                if c.type == "R":
                    S.g(a, b, 1 / c.val)
                elif c.type == "C":
                    v, i = etat[c.nom]
                    geq = (2 if trapeze else 1) * c.val / h
                    J = geq * v + (i if trapeze else 0.0)
                    S.g(a, b, geq)
                    S.i(a, b, -J)
                elif c.type == "L":
                    v, i = etat[c.nom]
                    geq = h / ((2 if trapeze else 1) * c.val)
                    K = i + (geq * v if trapeze else 0.0)
                    S.g(a, b, geq)
                    S.i(a, b, K)
                    S.A[c.k][c.k] += 1  # inconnue de courant inutilisée : équation x_k = i
                    S.b[c.k] += 0.0
                elif c.type == "V":
                    S.v(a, b, c.k, c.src(t))
                elif c.type == "I":
                    S.i(a, b, c.src(t))
                elif c.type == "E":
                    S.v(a, b, c.k, 0.0)
                    cp, cn = c.noeuds[2], c.noeuds[3]
                    if cp >= 0:
                        S.A[c.k][cp] -= c.gain
                    if cn >= 0:
                        S.A[c.k][cn] += c.gain

        x = resout_non_lineaire(ckt, tamponne, x)
        for c in ckt.composants:
            if c.type in "CL":
                v_nouv = tension(x, c.noeuds[0]) - tension(x, c.noeuds[1])
                v, i = etat[c.nom]
                if c.type == "C":
                    geq = (2 if trapeze else 1) * c.val / h
                    i_nouv = geq * v_nouv - (geq * v + (i if trapeze else 0.0))
                else:
                    geq = h / ((2 if trapeze else 1) * c.val)
                    i_nouv = geq * v_nouv + i + (geq * v if trapeze else 0.0)
                    x[c.k] = i_nouv  # pour pouvoir sonder I(L)
                etat[c.nom] = [v_nouv, i_nouv]
        ts.append(t)
        for s in series:
            series[s].append(mesure(ckt, x, s))
    return {"titre": "analyse temporelle", "x": ts, "xlabel": "temps (s)", "series": series, "log": False}


def analyse_ac(ckt, args):
    """Petits signaux : le circuit est linéarisé autour du point de repos
    (une diode devient une simple conductance), puis résolu en complexe
    pour chaque fréquence : un condensateur a l'admittance jωC, une bobine 1/(jωL)."""
    mode, n, f1, f2 = args[0].lower(), int(valeur(args[1])), valeur(args[2]), valeur(args[3])
    repos = point_de_repos(ckt)
    if mode == "dec":
        points = int(round(n * math.log10(f2 / f1)))
        freqs = [f1 * 10 ** (k / n) for k in range(points + 1)]
    elif mode == "lin":
        freqs = [f1 + (f2 - f1) * k / max(1, n - 1) for k in range(n)]
    else:
        raise SystemExit(".ac dec|lin points fdébut ffin")
    noms = sondes(ckt)
    series = {f"|{s}| (dB)": [] for s in noms}
    phases = {f"phase {s} (°)": [] for s in noms}
    for f in freqs:
        w = 2 * math.pi * f
        S = Systeme(ckt.taille, complexe=True)
        for c in ckt.composants:
            a, b = c.noeuds[0], c.noeuds[1]
            if c.type == "R":
                S.g(a, b, 1 / c.val)
            elif c.type == "C":
                S.g(a, b, 1j * w * c.val)
            elif c.type == "L":
                S.v(a, b, c.k, 0)
                S.A[c.k][c.k] -= 1j * w * c.val  # V(a) - V(b) - jωL·i = 0
            elif c.type in "VI":
                amp = cmath.rect(c.src.ac, math.radians(c.src.phase))
                S.v(a, b, c.k, amp) if c.type == "V" else S.i(a, b, amp)
            elif c.type == "D":
                v = tension(repos, a) - tension(repos, b)
                S.g(a, b, diode(c, v)[1] + 1e-12)
            elif c.type == "E":
                S.v(a, b, c.k, 0)
                cp, cn = c.noeuds[2], c.noeuds[3]
                if cp >= 0:
                    S.A[c.k][cp] -= c.gain
                if cn >= 0:
                    S.A[c.k][cn] += c.gain
        x = resous(S.A, S.b)
        for s in noms:
            z = mesure(ckt, x, s)
            series[f"|{s}| (dB)"].append(20 * math.log10(max(abs(z), 1e-30)))
            phases[f"phase {s} (°)"].append(math.degrees(cmath.phase(z)))
    return [{"titre": "réponse en fréquence (gain)", "x": freqs, "xlabel": "fréquence (Hz)", "series": series, "log": True},
            {"titre": "réponse en fréquence (phase)", "x": freqs, "xlabel": "fréquence (Hz)", "series": phases, "log": True}]


# ==================================================================== tracés

def trace_texte(r, largeur=72, hauteur=18):
    xs = r["x"]
    tx = [math.log10(v) for v in xs] if r["log"] else xs
    toutes = [v for s in r["series"].values() for v in s]
    ymin, ymax = min(toutes), max(toutes)
    if ymax - ymin < 1e-12:
        ymin, ymax = ymin - 1, ymax + 1
    grille = [[" "] * largeur for _ in range(hauteur)]
    marques = "●○■□▲△◆◇"
    for k, (nom, ys) in enumerate(r["series"].items()):
        for x, y in zip(tx, ys):
            c = int((x - tx[0]) / ((tx[-1] - tx[0]) or 1) * (largeur - 1))
            l = hauteur - 1 - int((y - ymin) / (ymax - ymin) * (hauteur - 1))
            grille[l][c] = marques[k % len(marques)]
    print(f"\n  {r['titre']}")
    for i, ligne in enumerate(grille):
        etiquette = joli(ymax - (ymax - ymin) * i / (hauteur - 1)) if i in (0, hauteur // 2, hauteur - 1) else ""
        print(f"{etiquette:>10} │{''.join(ligne)}")
    print(" " * 11 + "└" + "─" * largeur)
    print(" " * 11 + f"{joli(xs[0]):<{largeur // 2}}{joli(xs[-1]):>{largeur // 2}}  {r['xlabel']}")
    print("   " + "   ".join(f"{marques[k % len(marques)]} {nom}" for k, nom in enumerate(r["series"])))


def trace_svg(r, chemin, largeur=720, hauteur=360):
    """Un graphique SVG sans bibliothèque : axes, graduations, courbes, légende."""
    g, d, h_, b = 70, 20, 30, 50
    W, H = largeur - g - d, hauteur - h_ - b
    xs = r["x"]
    tx = [math.log10(v) for v in xs] if r["log"] else xs
    toutes = [v for s in r["series"].values() for v in s]
    ymin, ymax = min(toutes), max(toutes)
    marge = (ymax - ymin) * 0.08 or 1
    ymin, ymax = ymin - marge, ymax + marge
    X = lambda v: g + (v - tx[0]) / ((tx[-1] - tx[0]) or 1) * W
    Y = lambda v: h_ + (ymax - v) / (ymax - ymin) * H
    couleurs = ["#b45cff", "#ff8c42", "#3cc7d6", "#4cc76a", "#ff5c8a", "#f2d02b"]
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{largeur}" height="{hauteur}" viewBox="0 0 {largeur} {hauteur}" '
         f'font-family="sans-serif" font-size="12">',
         f'<rect width="{largeur}" height="{hauteur}" fill="#0b0710"/>',
         f'<text x="{g}" y="20" fill="#e8e0f0" font-size="14">{r["titre"]}</text>']

    def graduations(a, z, n=6):
        pas = 10 ** math.floor(math.log10((z - a) / n or 1))
        for m in (1, 2, 5, 10):
            if (z - a) / (pas * m) <= n:
                pas *= m
                break
        v = math.ceil(a / pas) * pas
        while v <= z + 1e-12:
            yield v
            v += pas

    for v in graduations(ymin, ymax):
        o.append(f'<line x1="{g}" x2="{g + W}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="#2a2133"/>')
        o.append(f'<text x="{g - 6}" y="{Y(v) + 4:.1f}" fill="#a89cb8" text-anchor="end">{joli(v)}</text>')
    if r["log"]:
        for e in range(math.floor(tx[0]), math.ceil(tx[-1]) + 1):
            if tx[0] - 1e-9 <= e <= tx[-1] + 1e-9:
                o.append(f'<line y1="{h_}" y2="{h_ + H}" x1="{X(e):.1f}" x2="{X(e):.1f}" stroke="#2a2133"/>')
                o.append(f'<text y="{h_ + H + 16}" x="{X(e):.1f}" fill="#a89cb8" text-anchor="middle">{joli(10 ** e)}</text>')
    else:
        for v in graduations(tx[0], tx[-1]):
            o.append(f'<line y1="{h_}" y2="{h_ + H}" x1="{X(v):.1f}" x2="{X(v):.1f}" stroke="#2a2133"/>')
            o.append(f'<text y="{h_ + H + 16}" x="{X(v):.1f}" fill="#a89cb8" text-anchor="middle">{joli(v)}</text>')
    o.append(f'<text x="{g + W / 2}" y="{hauteur - 10}" fill="#a89cb8" text-anchor="middle">{r["xlabel"]}</text>')
    o.append(f'<rect x="{g}" y="{h_}" width="{W}" height="{H}" fill="none" stroke="#4a3d5c"/>')
    for k, (nom, ys) in enumerate(r["series"].items()):
        pts = " ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in zip(tx, ys))
        c = couleurs[k % len(couleurs)]
        o.append(f'<polyline points="{pts}" fill="none" stroke="{c}" stroke-width="2" stroke-linejoin="round"/>')
        lx = g + W - 170 * (len(r["series"]) - k)  # légende en haut à droite, au-dessus du cadre
        o.append(f'<rect x="{lx}" y="14" width="12" height="3" fill="{c}"/>')
        o.append(f'<text x="{lx + 18}" y="19" fill="#e8e0f0">{nom}</text>')
    o.append("</svg>")
    with open(chemin, "w", encoding="utf-8") as f:
        f.write("\n".join(o))


# ======================================================================= main

def simule(chemin, dossier_svg=None, texte=True):
    ckt = Circuit(open(chemin, encoding="utf-8").read(), chemin)
    if ckt.titre:
        print(f"── {ckt.titre} ──")
    print(f"{len(ckt.composants)} composants, {ckt.n_noeuds} nœuds, système {ckt.taille}×{ckt.taille}")
    resultats = []
    for cmd, args in ckt.analyses or [(".op", [])]:
        r = {".op": lambda a: analyse_op(ckt), ".dc": lambda a: analyse_dc(ckt, a),
             ".tran": lambda a: analyse_tran(ckt, a), ".ac": lambda a: analyse_ac(ckt, a)}[cmd](args)
        if r:
            resultats += r if isinstance(r, list) else [r]
    base = os.path.splitext(os.path.basename(chemin))[0]
    for k, r in enumerate(resultats):
        if texte:
            trace_texte(r)
        if dossier_svg:
            os.makedirs(dossier_svg, exist_ok=True)
            nom = os.path.join(dossier_svg, f"{base}{'-' + str(k + 1) if len(resultats) > 1 else ''}.svg")
            trace_svg(r, nom)
            print(f"  -> {nom}")
    return ckt, resultats


def main():
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--aide"):
        print(__doc__)
        print("usage : python3 circuit.py netlist.cir [--svg dossier] [--sans-texte]")
        return
    svg = args[args.index("--svg") + 1] if "--svg" in args else None
    fichiers = [a for i, a in enumerate(args) if not a.startswith("--") and (i == 0 or args[i - 1] != "--svg")]
    for f in fichiers:
        simule(f, svg, "--sans-texte" not in args)


if __name__ == "__main__":
    main()
