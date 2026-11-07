#!/usr/bin/env python3
"""Compare le simulateur à des résultats calculés à la main."""
import math
import circuit as c

resultats = []


def verifie(nom, obtenu, attendu, tol):
    ok = abs(obtenu - attendu) <= tol
    resultats.append(ok)
    print(f"  {'ok   ' if ok else 'ÉCHEC'} {nom} : {obtenu:.5g} (attendu {attendu:.5g} ± {tol:g})")


def ckt(texte):
    return c.Circuit(texte)


def a_frequence(r, f):
    """Valeur interpolée d'une série AC à la fréquence f (axe log)."""
    xs = r["x"]
    k = min(range(len(xs)), key=lambda i: abs(math.log(xs[i] / f)))
    return {nom: s[k] for nom, s in r["series"].items()}, xs[k]


# 1. pont diviseur
x = c.point_de_repos(ckt("V1 a 0 12\nR1 a b 10k\nR2 b 0 4.7k"))
verifie("diviseur : V(b) = 12·4.7/14.7", x[1], 12 * 4.7 / 14.7, 1e-9)

# 2. RC : un tau après l'échelon, 63 % de la valeur finale
k = ckt("V1 e 0 PULSE(0 5 0 1n 1n 1 2)\nR1 e s 1k\nC1 s 0 1u\n.sonde V(s)")
r = c.analyse_tran(k, ["1u", "3m"])
i = r["x"].index(min(r["x"], key=lambda t: abs(t - 1e-3)))
verifie("RC : V(s) à t = τ vaut 5·(1 − 1/e)", r["series"]["V(s)"][i], 5 * (1 - math.exp(-1)), 5e-3)

# 3. RC en fréquence : -3 dB et -45° à fc = 1/(2πRC)
k = ckt("V1 e 0 AC 1\nR1 e s 1k\nC1 s 0 1u\n.sonde V(s)")
fc = 1 / (2 * math.pi * 1e3 * 1e-6)
gain, phase = c.analyse_ac(k, ["lin", "1", str(fc), str(fc)])
verifie("RC : gain à fc", gain["series"]["|V(s)| (dB)"][0], -10 * math.log10(2), 1e-6)
verifie("RC : phase à fc", phase["series"]["phase V(s) (°)"][0], -45, 1e-6)

# 4. RLC série : pic de résonance de la tension du condensateur
R, L, C = 20, 10e-3, 100e-9
Q, f0 = math.sqrt(L / C) / R, 1 / (2 * math.pi * math.sqrt(L * C))
k = ckt(f"V1 e 0 AC 1\nR1 e a {R}\nL1 a b {L}\nC1 b 0 {C}\n.sonde V(b)")
gain, _ = c.analyse_ac(k, ["dec", "2000", "1k", "20k"])
s = gain["series"]["|V(b)| (dB)"]
kmax = max(range(len(s)), key=s.__getitem__)
verifie("RLC : fréquence du pic f0·√(1 − 1/2Q²)", gain["x"][kmax], f0 * math.sqrt(1 - 1 / (2 * Q * Q)), 10)
verifie("RLC : hauteur du pic Q/√(1 − 1/4Q²) en dB", s[kmax], 20 * math.log10(Q / math.sqrt(1 - 1 / (4 * Q * Q))), 1e-3)

# 5. RLC en temporel : pseudo-période des oscillations
k = ckt(f"V1 e 0 PULSE(0 1 0 1n 1n 1 2)\nR1 e a {R}\nL1 a b {L}\nC1 b 0 {C}\n.sonde V(b)")
r = c.analyse_tran(k, ["0.5u", "1m"])
v, t = r["series"]["V(b)"], r["x"]
maxima = [t[i] for i in range(1, len(v) - 1) if v[i] > v[i - 1] and v[i] >= v[i + 1]]
w0, alpha = 1 / math.sqrt(L * C), R / (2 * L)
verifie("RLC : pseudo-période 2π/√(ω0² − α²)", maxima[1] - maxima[0], 2 * math.pi / math.sqrt(w0 ** 2 - alpha ** 2), 2e-6)

# 6. amplificateur inverseur
k = ckt("V1 e 0 AC 1\nR1 e m 1k\nR2 m s 10k\nE1 s 0 0 m 100k\n.sonde V(s)")
gain, phase = c.analyse_ac(k, ["lin", "1", "1k", "1k"])
verifie("inverseur : gain −R2/R1 = −10 → 20 dB", gain["series"]["|V(s)| (dB)"][0], 20, 1e-3)
verifie("inverseur : déphasage 180°", abs(phase["series"]["phase V(s) (°)"][0]), 180, 1e-3)

# 7. Sallen-Key Butterworth : -3 dB à fc, -40 dB/décade au-delà
k = c.Circuit(open("exemples/sallen_key.cir", encoding="utf-8").read())
gain, _ = c.analyse_ac(k, ["lin", "1", "1000.4", "1000.4"])
verifie("Sallen-Key : -3 dB à fc", gain["series"]["|V(sortie)| (dB)"][0], -10 * math.log10(2), 0.01)
g1, _ = c.analyse_ac(k, ["lin", "1", "100k", "100k"])
g2, _ = c.analyse_ac(k, ["lin", "1", "1meg", "1meg"])
verifie("Sallen-Key : pente en dB par décade", g2["series"]["|V(sortie)| (dB)"][0] - g1["series"]["|V(sortie)| (dB)"][0], -40, 0.5)

# 8. diode : loi de Shockley inversée à 1 mA
Is, n = 1e-14, 1.0
k = ckt("I1 0 a 1m\nD1 a 0")
x = c.point_de_repos(k)
verifie("diode : V à 1 mA = n·Vt·ln(1 + I/Is)", x[0], n * c.VT * math.log(1 + 1e-3 / Is), 1e-4)

# 9. redresseur : tension lissée ≈ crête − 2 diodes
k = c.Circuit(open("exemples/redresseur.cir", encoding="utf-8").read())
r = c.analyse_tran(k, ["50u", "60m"])
fin = r["series"]["V(plus)"][-400:]
verifie("pont de diodes : tension moyenne en sortie", sum(fin) / len(fin), 10.4, 0.6)

print(f"{sum(resultats)} / {len(resultats)}")
raise SystemExit(0 if all(resultats) else 1)
