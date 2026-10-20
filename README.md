# tenth-discriminant

Class number one for imaginary quadratic fields, derived from W. D. Banks's explicit bound excluding Landau–Siegel zeros ([arXiv:2610.15483](https://arxiv.org/abs/2610.15483)), without Heegner points.

Written in response to Terence Tao's challenge, [*An experimental AI-math challenge: using tenth discriminant bounds to explore arguments ruling out Siegel zeroes*](https://terrytao.wordpress.com/2026/10/20/an-experimental-ai-math-challenge-using-tenth-discriminant-bounds-to-explore-arguments-ruling-out-siegel-zeroes/) (20 October 2026).

## Result

- **Zero forcing.** If $h(-d) = 1$ and $d \ge 10^{10}$, then $L(s,\chi_{-d})$ has a real zero $\beta$ with $1 - \beta < 2/\sqrt d$. The sharp asymptotic constant is $6/\pi$.
- **With Banks's theorem.** Banks's Theorem 1.1 gives $(1-\beta)\log q \ge 3.0 \times 10^{-3}$ for $q \ge 10^{10}$. The lemma gives $(1-\beta)\log d < 4.6 \times 10^{-4}$, so any tenth discriminant has $d < 10^{10}$.
- **Sieve.** No $d$ with $163 < d < 10^{10}$ satisfies Rabinowitsch's condition.

The result is conditional on Banks's Theorem 1.1 (first assertion). The full argument is in [NOTE.md](NOTE.md).

## Files

| File | Purpose |
|---|---|
| `NOTE.md` | The note: reductions, zero-forcing lemma with proof, combination with Banks, the sieve, verification |
| `sieve_h1.py` | Shows no $d \equiv 3 \pmod 8$ in $(163, 10^{10})$ has $(-d/p) = -1$ for all primes $p < (d+1)/4$ |
| `check_small.py` | Independent sanity checks: brute-force version of the sieve condition for $d < 2 \times 10^6$, and class numbers by counting reduced forms for $\lvert D \rvert < 2 \times 10^4$ |
| `check_lemma.py` | Checks the Chowla–Selberg expansion of $\zeta_K$ against $\zeta(s)L(s,\chi)$ for small Heegner $d$, and looks for real zeros of $L$ on $(\tfrac12, 1)$ |

## Running

Requires Python 3 with `numpy`, `sympy` and `mpmath`.

```bash
pip install numpy sympy mpmath
python3 sieve_h1.py 1e10
python3 check_small.py
python3 check_lemma.py
```

`sieve_h1.py` takes the upper limit as an optional argument (default $10^{10}$).

## Expected output

`sieve_h1.py 1e10` (about 0.4 s):

```
CRT modulus 892371480, admissible classes 142560
candidates after base sieve: 1597562
survivors after primes < 2000: 0
```

`check_small.py` (a few seconds):

```
brute-force all-inert d < 2e6: [3, 11, 19, 43, 67, 163]
h = 1 fundamental discriminants |D| < 2e4: [3, 4, 7, 8, 11, 19, 43, 67, 163]
```

`check_lemma.py`: for each $d \in \{11, 19, 43, 67, 163\}$, a table comparing $\zeta_K(s)$ computed directly and from the expansion, agreeing to about $10^{-29}$, followed by `no sign change of L on (1/2, 1)`. The absence of real zeros for these small $d$ is expected: the zero-forcing lemma is a statement about large $d$ (see NOTE.md, Remark 2.5).

## Dependencies of the argument

Used: Banks's Theorem 1.1 (first assertion), the Chowla–Selberg expansion, elementary bounds for $\zeta$, $\Gamma$, $\psi$ and $K_\nu$, Rabinowitsch's criterion, and the computation above.

Not used: Heegner points, Gross–Zagier, Baker's method, or the Lu–Zaman–Zhao computation that Banks's second assertion relies on.

## AI disclosure

The lemma, its proof, the note and the scripts were developed with Claude (Anthropic). The author has checked the argument and run the scripts.

## Author

Michael M. Ross · michaelmross@cantab.net · [michaelmross.github.io](https://michaelmross.github.io)
