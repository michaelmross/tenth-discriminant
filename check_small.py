"""Independent sanity checks for the sieve.
1. Brute force: d = 3 mod 8, d < 2e6, all primes p < (d+1)/4 have (-d/p) = -1.
   Expect exactly {3, 11, 19, 43, 67, 163}.
2. Class numbers by counting reduced forms for all fundamental -d, d < 2e4;
   h = 1 list should be the nine Heegner numbers.
"""
from sympy import primerange, legendre_symbol
from math import isqrt

PR = list(primerange(3, 10**6))

def all_inert(d):
    if d % 8 != 3:
        return False
    A = (d + 1) // 4
    for p in PR:
        if p >= A:
            return True
        if legendre_symbol((-d) % p, p) != -1:
            return False
    raise RuntimeError("prime list too short")

hits = [d for d in range(3, 2_000_000, 8) if all_inert(d)]
print("brute-force all-inert d < 2e6:", hits)

def h(D):  # D < 0, count reduced forms (a,b,c), b^2-4ac = D
    n = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a) == 0:
                c = (b * b - D) // (4 * a)
                if c >= a and not (b < 0 and a == c):
                    n += 1
        a += 1
    return n

def fundamental(D):
    if D % 4 == 1:
        m = -D
        return all(m % (p * p) for p in range(2, isqrt(m) + 1))
    if D % 4 == 0:
        m = -D // 4
        return m % 4 in (1, 2) and all(m % (p * p) for p in range(2, isqrt(m) + 1))
    return False

ones = [-D for D in range(-3, -20000, -1) if fundamental(D) and h(D) == 1]
print("h = 1 fundamental discriminants |D| < 2e4:", ones)
