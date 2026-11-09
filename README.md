<p align="center"><img src="docs/banner-circuit-spice.png" alt="circuit-spice" width="100%"></p>

<p align="center"><img src="https://img.shields.io/badge/Python-no%20dependencies-3776ab?style=flat-square" alt="Python"> <img src="https://img.shields.io/badge/license-MIT-610089?style=flat-square" alt="license"></p>

An electrical circuit simulator in the style of **SPICE**, written from scratch in pure Python: no numpy, no scipy. You describe the circuit in a text file, the program writes Kirchhoff's laws as a matrix and solves them.

<p align="center"><a href="README.fr.md">🇫🇷 Version française</a></p>

![Damped oscillations of a series RLC](images/rlc-1.svg)

![Bode plot of a 2nd-order Sallen-Key filter](images/sallen_key-1.svg)

![Diode bridge with a smoothing capacitor](images/redresseur.svg)

## Usage

```sh
python3 circuit.py exemples/rc.cir                 # curves in the terminal
python3 circuit.py exemples/rlc.cir --svg images   # + SVG charts
python3 tests.py                                   # 13 comparisons with theory
```

A netlist looks like this (`.sonde` = probe):

```
* Filtre RC passe-bas
V1 entree 0 PULSE(0 5 0.5m 1u 1u 10m 20m) AC 1
R1 entree sortie 1k
C1 sortie 0 1u
.tran 10u 6m
.ac dec 20 1 100k
.sonde V(entree) V(sortie)
```

| Line | Meaning |
|---|---|
| `R1 a b 4.7k` | resistor (also `C`, `L`); suffixes `f p n u m k meg g` |
| `V1 a b 5` | voltage source; also `SIN(offset amplitude frequency)`, `PULSE(v1 v2 delay rise fall width period)`, `PWL(t1 v1 t2 v2 ...)`, `AC amplitude` |
| `I1 a b 1m` | current source (from a to b through the source) |
| `D1 anode cathode IS=2.5n N=1.75` | diode (Shockley's law) |
| `E1 s 0 p m 100k` | voltage-controlled voltage source: V(s) = 100,000·(V(p) − V(m)), an ideal op-amp |
| `.op` `.dc V1 0 5 0.1` `.tran step end` `.ac dec n f1 f2` | the four analyses |
| `.sonde V(n) V(a,b) I(V1)` | what to plot |

Bundled examples: voltage divider, RC, RLC, diode bridge, diode characteristic, inverting amplifier, Sallen-Key filter.

## How it works

### Modified nodal analysis (MNA)

The unknowns are the node voltages (ground is 0) and the current of each voltage source. Each component adds its contribution, its "stamp", into a matrix `A` and a right-hand side `b`, then `A·x = b` is solved by Gaussian elimination with partial pivoting.

- a resistor between `a` and `b` adds `1/R` at `(a,a)` and `(b,b)`, and subtracts `1/R` at `(a,b)` and `(b,a)`: that's Kirchhoff's current law;
- a voltage source adds a row `V(a) − V(b) = E` and a column for its current.

### Transient: companion models

A capacitor obeys `i = C·dv/dt`. With a step `h`, the **trapezoidal rule** gives

```
i(n+1) = (2C/h)·(v(n+1) − v(n)) − i(n)
```

That's a conductance `2C/h` in parallel with a current source that only depends on the past: at every step, the circuit becomes a circuit of resistors and sources again. Same thing for an inductor. The first step uses backward Euler, which damps the start-up shock.

### Diodes: Newton-Raphson

`i = Is·(e^(v/nVt) − 1)` isn't linear. It is replaced by its tangent at the current point (a conductance plus a current source), solved, and repeated until nothing moves any more. Two details taken from SPICE:

- **step limiting** (`pnjlim`): if Newton proposes 5 V across a diode, `e^(5/0.026)` overflows; the step is shortened logarithmically;
- **gmin**: a 10⁻¹² S leak in parallel with each diode. It has to be a real conductance (in the current *and* in the matrix), otherwise Newton converges at a snail's pace: that's exactly what was stalling the diode bridge at first.

### AC: complex numbers

The circuit is linearised around its operating point, then solved for each frequency with complex admittances: `jωC` for a capacitor, `1/(jωL)` for an inductor. The same Gauss code works as is, since Python can divide complex numbers.

## Checks (`tests.py`)

The output is in French: `diviseur` = divider, `pont de diodes` = diode bridge.

```
ok  diviseur : V(b) = 12·4.7/14.7
ok  RC : V(s) à t = τ vaut 5·(1 − 1/e)
ok  RC : gain à fc                          -3.0103 dB
ok  RC : phase à fc                         -45°
ok  RLC : fréquence du pic f0·√(1 − 1/2Q²)
ok  RLC : hauteur du pic Q/√(1 − 1/4Q²) en dB
ok  RLC : pseudo-période 2π/√(ω0² − α²)
ok  inverseur : gain −R2/R1 = −10 → 20 dB
ok  inverseur : déphasage 180°
ok  Sallen-Key : -3 dB à fc
ok  Sallen-Key : pente en dB par décade     -40
ok  diode : V à 1 mA = n·Vt·ln(1 + I/Is)
ok  pont de diodes : tension moyenne en sortie
```

## Limits

No transistors, no adaptive time step, no sparse matrices: it's made for understanding, and for circuits of a few dozen components.
