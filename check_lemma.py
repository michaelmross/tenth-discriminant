"""Numerical checks of the zero-forcing lemma.

For h(-d) = 1, d > 4:
  zeta_K(s) = zeta(2s) + g(s) (sqrt(d)/2)^(1-2s) zeta(2s-1) + E(s),
  g(s) = sqrt(pi) Gamma(s-1/2)/Gamma(s),
  E(s) = (4 pi^s/Gamma(s)) sum_{y,m>=1} (m/c_y)^(s-1/2) K_{s-1/2}(2 pi m c_y) cos(pi m y),
  c_y = y sqrt(d)/2.
Check against zeta(s) L(s, chi_{-d}) computed independently, then locate the real
zero beta and compare 1 - beta with 6/(pi sqrt d) and the lemma's bound 2/sqrt d.
"""
from mpmath import mp, mpf, zeta, gamma, besselk, cos, pi, sqrt, dirichlet, findroot
from sympy import jacobi_symbol

mp.dps = 30

def chi_list(d):
    # Kronecker (-d/n) = Jacobi (n/d) for d prime, d = 3 mod 4
    return [int(jacobi_symbol(n, d)) if n % d else 0 for n in range(d)]

def zetaK_direct(s, d, chi):
    return zeta(s) * dirichlet(s, chi)

def zetaK_formula(s, d, Y=40, M=40):
    g = sqrt(pi) * gamma(s - mpf(1) / 2) / gamma(s)
    main = zeta(2 * s) + g * (sqrt(d) / 2) ** (1 - 2 * s) * zeta(2 * s - 1)
    E = mpf(0)
    for y in range(1, Y + 1):
        c = y * sqrt(d) / 2
        for m in range(1, M + 1):
            E += (m / c) ** (s - mpf(1) / 2) * besselk(s - mpf(1) / 2, 2 * pi * m * c) * cos(pi * m * y)
    E *= 4 * pi ** s / gamma(s)
    return main + E, E

for d in [11, 19, 43, 67, 163]:
    chi = chi_list(d)
    print(f"\nd = {d}")
    for s in [mpf("0.6"), mpf("0.9"), mpf("0.99"), mpf("1.5"), mpf(2)]:
        a = zetaK_direct(s, d, chi)
        b, E = zetaK_formula(s, d)
        print(f"  s={float(s):4.2f}  direct={float(a): .12e}  formula={float(b): .12e}  |diff|={float(abs(a-b)):.1e}  |E|={float(abs(E)):.1e}")
    # real zero of L(s, chi) in (1/2, 1)
    L = lambda s: dirichlet(s, chi)
    lo, hi = mpf("0.5001"), mpf("0.99999")
    if L(lo) * L(hi) > 0:
        print("  no sign change of L on (1/2, 1)")
        continue
    beta = findroot(L, (lo, hi), solver="anderson")
    print(f"  beta = {float(beta):.10f}   1-beta = {float(1-beta):.6f}   "
          f"6/(pi sqrt d) = {float(6/(pi*sqrt(d))):.6f}   2/sqrt d = {float(2/sqrt(d)):.6f}")
