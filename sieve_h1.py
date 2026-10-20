"""
Rabinowitsch sieve: no imaginary quadratic field Q(sqrt(-d)) with 163 < d < LIMIT
has class number one.

If h(-d) = 1 and d > 8, then D = -d is odd (a principal ideal of norm 2 would need
2 = x^2 + m y^2 with m >= 3 when D = -4m), and every prime p < (d+1)/4 is inert:
an ideal of norm p would be principal, p = x^2 + xy + ((d+1)/4) y^2, forcing
y != 0 and hence p >= (d+1)/4.  For p = 2 inertness means d = 3 (mod 8); for odd p
it means the Legendre symbol (-d / p) = -1.

We show every d = 3 (mod 8) with 163 < d < LIMIT has some prime p < (d+1)/4 with
(-d/p) != -1.  Candidates are enumerated by CRT over the first few odd primes and
then filtered by further primes, vectorised in numpy.
"""
import sys
import numpy as np
from math import prod
from sympy import primerange

LIMIT = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**10
D_MIN = 163

def nonres_classes_for_d(p):
    """Residues r mod p such that -r is a quadratic non-residue mod p."""
    squares = {(x * x) % p for x in range(1, p)}
    return [r for r in range(p) if r % p != 0 and ((-r) % p) not in squares]

# Stage 1: CRT enumeration over modulus 8 * prod(base primes)
base = [3, 5, 7, 11, 13, 17, 19, 23]
mod = 8
classes = np.array([3], dtype=np.int64)          # d = 3 mod 8
for p in base:
    allowed = nonres_classes_for_d(p)
    new = []
    # combine classes mod `mod` with classes mod p
    inv = pow(mod, -1, p)
    for r in allowed:
        # x = c + mod*t, x = r mod p  ->  t = (r - c)*inv mod p
        t = ((r - classes) * inv) % p
        new.append(classes + mod * t)
    classes = np.concatenate(new)
    mod *= p
classes.sort()
print(f"CRT modulus {mod}, admissible classes {len(classes)}")

# Stage 2: expand to all d < LIMIT in these classes
reps = LIMIT // mod + 1
cand = (classes[None, :] + mod * np.arange(reps, dtype=np.int64)[:, None]).ravel()
cand = cand[(cand > D_MIN) & (cand < LIMIT)]
print(f"candidates after base sieve: {len(cand)}")

# Small d where base primes may exceed (d+1)/4 are handled exactly in stage 3.
# Stage 3: filter by further primes; a candidate is eliminated by prime p only if
# p < (d+1)/4 and (-d/p) != -1.
small_mask = cand < 4 * 23 + 4   # none expected, since d > 163
assert not small_mask.any()

for p in primerange(29, 2000):
    if len(cand) == 0:
        break
    allowed = np.zeros(p, dtype=bool)
    allowed[nonres_classes_for_d(p)] = True
    applies = p < (cand + 1) // 4          # (d+1)/4 is an integer since d = 3 mod 8
    survive = allowed[cand % p] | ~applies
    cand = cand[survive]

print(f"survivors after primes < 2000: {len(cand)}")
if len(cand):
    print("survivors:", cand[:50])
