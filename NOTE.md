# Class number one from Banks' explicit Siegel-zero bound

Michael M. Ross · michaelmross@cantab.net · 20 October 2026

## AI acknowledgment

The lemma, its proof, this note and the scripts were developed with Claude Opus (Anthropic). I have checked the argument and run the scripts.

## Summary

We show that the explicit exclusion of Landau–Siegel zeros given by Banks [1], following the OpenAI manuscript [4] (see also the expositions of Green and Kerr [2]), implies the class number one theorem for imaginary quadratic fields. The argument does not use Heegner points or the small-conductor computation of Lu, Zaman and Zhao.

It has three steps.

1. **Zero forcing (Lemma 2).** If $h(-d) = 1$ and $d \ge 10^{10}$, then $L(s,\chi_{-d})$ has a real zero $\beta$ with $1 - \beta < 2/\sqrt d$. The proof uses only the Chowla–Selberg expansion of $\zeta_K$ and two textbook inequalities for $\zeta$.
2. **Banks' theorem.** For conductors $q \ge 10^{10}$, every real zero satisfies $(1-\beta)\log q \ge 3.0 \times 10^{-3}$. Lemma 2 gives $(1-\beta)\log d < 4.6 \times 10^{-4}$, so any tenth discriminant has $d < 10^{10}$.
3. **Sieve.** Rabinowitsch's criterion forces $(-d/p) = -1$ for every prime $p < (d+1)/4$. A numpy sieve shows that no $d$ in $(163, 10^{10})$ satisfies this. The run takes under a second.

This is the route of the comment by twinheterodox on Tao's challenge post [6], which proves $1-\beta \le 10/\sqrt q$ by a theta-integral argument. The contribution here is the constant $2$ in place of $10$, with $6/\pi$ identified as the sharp asymptotic constant, and a short proof whose analytic inputs are all classical.

**Everything is conditional on Banks' Theorem 1.1 (first assertion).** That is the step that most needs independent checking.

## 1. Reductions

Let $D = -d$ be a negative fundamental discriminant with $d > 8$, $K = \mathbb{Q}(\sqrt{D})$, $\chi = \chi_D$ the associated Kronecker character, and $h = h(D)$.

**Lemma 1 (Rabinowitsch).** If $h = 1$ and $d > 8$, then $d \equiv 3 \pmod 8$, and $(-d/p) = -1$ for every odd prime $p < (d+1)/4$.

*Proof.* If $D \equiv 0 \pmod 4$, write $D = -4m$ with $m \equiv 1, 2 \pmod 4$ and $m \ge 3$. The prime $2$ ramifies, so there is an ideal of norm $2$. It is principal, so $2 = x^2 + my^2$, which is impossible for $m \ge 3$. Hence $D \equiv 1 \pmod 4$ and $d \equiv 3 \pmod 4$.

Then $\mathcal{O}_K = \mathbb{Z}[\omega]$ with $\omega = (1+\sqrt{-d})/2$, and the norm form is

$$Q(x,y) = x^2 + xy + \tfrac{d+1}{4}\,y^2 = \left(x + \tfrac{y}{2}\right)^2 + \tfrac{d}{4}\,y^2 .$$

Suppose a prime $p < (d+1)/4$ is not inert. Then there is an ideal of norm $p$, and it is principal, so $p = Q(x,y)$. If $y = 0$ then $p = x^2$, which is impossible. If $|y| = 1$ then $Q(x,y) = x^2 \pm x + (d+1)/4 \ge (d+1)/4$. If $|y| \ge 2$ then $Q(x,y) \ge d$. Each case contradicts $p < (d+1)/4$.

So every prime below $(d+1)/4$ is inert. For $p = 2$ (which lies below $(d+1)/4$ once $d > 7$) this means $D \equiv 5 \pmod 8$, that is $d \equiv 3 \pmod 8$. For odd $p$ it means $(-d/p) = -1$. ∎

## 2. Zero forcing

**Lemma 2.** Let $D = -d$ be a fundamental discriminant with $h(D) = 1$ and $d \ge 10^{10}$. Then $L(s,\chi_D)$ has a real zero $\beta$ with

$$1 - \beta < \frac{2}{\sqrt d}.$$

### 2.1 The Chowla–Selberg expansion

Since $h = 1$ and the unit group has order $w = 2$, for $`\mathrm{Re}(s) > 1`$

$$\zeta_K(s) = \sum_{\mathfrak a} N\mathfrak a^{-s} = \frac12 \sum_{(x,y) \ne (0,0)} Q(x,y)^{-s}.$$

The terms with $y = 0$ contribute $\zeta(2s)$. By the symmetry $(x,y) \mapsto (-x,-y)$, the remaining terms contribute

$$\sum_{y \ge 1} \sum_{x \in \mathbb{Z}} \left( (x + \alpha_y)^2 + c_y^2 \right)^{-s}, \qquad \alpha_y = \frac{y}{2}, \quad c_y = \frac{y\sqrt d}{2}.$$

For $c > 0$ and $`\mathrm{Re}(s) > \tfrac12`$, the Fourier transform of $t \mapsto (t^2 + c^2)^{-s}$ is

$$\int_{\mathbb R} (t^2+c^2)^{-s} e^{-2\pi i \xi t}\,dt = \begin{cases} \sqrt{\pi}\,\dfrac{\Gamma(s-\frac12)}{\Gamma(s)}\, c^{1-2s}, & \xi = 0, \\[2ex] \dfrac{2\pi^s}{\Gamma(s)} \left( \dfrac{|\xi|}{c} \right)^{s-\frac12} K_{s-\frac12}(2\pi|\xi| c), & \xi \ne 0. \end{cases}$$

Poisson summation in $x$, followed by summation over $y$, gives

$$\boxed{\;\zeta_K(s) = \zeta(2s) + g(s) \left( \frac{\sqrt d}{2} \right)^{1-2s} \zeta(2s-1) + E(s)\;}$$

where

$$g(s) = \sqrt{\pi}\,\frac{\Gamma(s - \frac12)}{\Gamma(s)}, \qquad E(s) = \frac{4\pi^s}{\Gamma(s)} \sum_{y \ge 1} \sum_{m \ge 1} \left( \frac{m}{c_y} \right)^{s-\frac12} K_{s-\frac12}(\pi m y \sqrt d)\, \cos(\pi m y).$$

The double sum defining $E$ converges absolutely and locally uniformly in $s$, because $K_\nu(z)$ decays like $e^{-z}$. So $E$ is entire, and the identity extends from $`\mathrm{Re}(s) > 1`$ to all $s \ne 1$ by analytic continuation.

As a check on normalization, the residue at $s = 1$ is $g(1) \cdot (2/\sqrt d) \cdot \tfrac12 = \pi/\sqrt d$. This equals $L(1,\chi)$, as the class number formula requires for $h = 1$, $w = 2$. The identity is also checked numerically in `check_lemma.py` (Section 5).

### 2.2 Bounding $E$ near $s = 1$

Let $\sigma \in [0.99, 1]$ and $\nu = \sigma - \tfrac12 \in [0.49, 0.5]$.

- For fixed $z > 0$, $K_\nu(z)$ is increasing in $\nu \ge 0$, so $K_\nu(z) \le K_{1/2}(z) = \sqrt{\pi/(2z)}\, e^{-z}$. Here $z = \pi m y \sqrt d \ge \pi$, so $\sqrt{\pi/(2z)} \le 1$.
- $c_y \ge \sqrt d / 2 \ge 1$, so $(m/c_y)^{\nu} \le m^{1/2}$.
- $\Gamma$ is decreasing on $(0, 1.46)$, so $\Gamma(\sigma) \ge \Gamma(1) = 1$ and $4\pi^\sigma/\Gamma(\sigma) \le 4\pi$.

With $a = \pi\sqrt d$, this gives

$$|E(\sigma)| \le 4\pi \sum_{m \ge 1} m^{1/2} \sum_{y \ge 1} e^{-a m y} \le 4\pi \cdot (1.001)^2\, e^{-a} \le 8\pi\, e^{-\pi \sqrt d},$$

since $a \ge \pi \cdot 10^5$ makes both geometric tails negligible. For $d \ge 10^{10}$ this is far below $10^{-100000}$.

### 2.3 The sign of $\zeta_K$ at $\sigma = 1 - 2/\sqrt d$

Put $\eta = 2/\sqrt d \le 2 \times 10^{-5}$ and $\sigma = 1 - \eta$. We bound each term of the expansion.

**(i) First term.** $\zeta$ is decreasing on $(1, \infty)$ and $2\sigma < 2$, so $\zeta(2\sigma) > \zeta(2) = \pi^2/6 > 1.6449$.

**(ii) $\zeta(2\sigma - 1)$.** For $0 < u < 1$,

$$\zeta(u) = \frac{u}{u-1} - u \int_1^\infty \{t\}\, t^{-u-1}\, dt .$$

The integral lies in $(0, 1/u)$, so $\dfrac{1}{u-1} \le \zeta(u) < \dfrac{u}{u-1} < 0$. With $u = 1 - 2\eta$ this gives $-\dfrac{1}{2\eta} \le \zeta(2\sigma - 1) < 0$.

**(iii) $g(\sigma)$.** We have $\dfrac{d}{d\eta} \log g(1-\eta) = \psi(1-\eta) - \psi(\tfrac12 - \eta)$, which is positive and, for $0 \le \eta \le 2 \times 10^{-5}$, at most $1.39$ (it equals $2\log 2 \approx 1.386$ at $\eta = 0$). Since $g(1) = \pi$,

$$\pi < g(\sigma) \le \pi e^{1.39\eta} \le \pi (1 + 3 \times 10^{-5}).$$

**(iv) The power of $d$.** $(\sqrt d/2)^{1-2\sigma} = \dfrac{2}{\sqrt d}\left(\dfrac d4\right)^{\eta}$, and

$$\left(\frac d4\right)^{\eta} \le \exp\!\left( \frac{2\log d}{\sqrt d} \right) \le \exp(4.61 \times 10^{-4}) \le 1.0005,$$

using that $\log d/\sqrt d$ is decreasing for $d \ge e^2$.

**Combining.** The second term of the expansion is negative, with absolute value at most

$$\pi(1 + 3\times 10^{-5}) \cdot \frac{2}{\sqrt d} \cdot 1.0005 \cdot \frac{1}{2\eta} = \frac{\pi \cdot 1.00053}{\eta \sqrt d} = \frac{\pi \cdot 1.00053}{2} < 1.5717 .$$

Therefore

$$\zeta_K(\sigma) \ge 1.6449 - 1.5717 - 8\pi e^{-\pi\sqrt d} > 0.07 .$$

### 2.4 Conclusion of the proof

On $(\tfrac12, 1)$, $\zeta_K(s) = \zeta(s) L(s,\chi)$ is real and continuous. As $s \to 1^-$, $\zeta(s) \to -\infty$ while $L(s,\chi) \to L(1,\chi) = \pi/\sqrt d > 0$, so $\zeta_K(s) \to -\infty$. Since $\zeta_K(\sigma) > 0$, the intermediate value theorem gives $\beta \in (\sigma, 1)$ with $\zeta_K(\beta) = 0$. By (ii), $\zeta(\beta) < 0$, so $L(\beta,\chi) = 0$, and $1 - \beta < \eta = 2/\sqrt d$. ∎

### 2.5 Remarks

**Sharpness.** Balancing the two main terms, $\zeta(2) \approx \pi/(\eta\sqrt d)$, gives $\eta \approx 6/(\pi\sqrt d) \approx 1.91/\sqrt d$. Evaluating the main terms at $s = 1 - K/\sqrt d$ (the term $E$ is negligible) confirms this:

| $d$ | $K = 1.8$ | $K = 1.9$ | $K = 1.92$ | $K = 2$ |
|---|---|---|---|---|
| $10^{7}$ | $-0.1143$ | $-0.0223$ | $-0.0051$ | $+0.0604$ |
| $10^{8}$ | $-0.1055$ | $-0.0136$ | $+0.0036$ | $+0.0691$ |
| $10^{10}$ | $-0.1010$ | $-0.0092$ | $+0.0080$ | $+0.0735$ |
| $10^{12}$ | $-0.1005$ | $-0.0086$ | $+0.0086$ | $+0.0741$ |
| $10^{20}$ | $-0.1004$ | $-0.0085$ | $+0.0087$ | $+0.0741$ |

The margin $+0.0735$ at $d = 10^{10}$, $K = 2$ matches the analytic lower bound $0.073$ from Section 2.3.

**The lemma is asymptotic.** For $d = 163$, $L(s,\chi_{-163})$ has no sign change on $(\tfrac12, 1)$ and is positive at every sampled point of $[0.3, 0.95]$. The quantitative inequalities of Section 2.3 are what make the lemma work, so a size hypothesis on $d$ is genuinely needed. The table suggests the sign condition already holds from about $d = 10^7$, but only $d \ge 10^{10}$ is proved here, and that is all Section 3 requires.

## 3. Combining with Banks' theorem

**Theorem (Banks [1], Theorem 1.1, first assertion).** For every primitive real Dirichlet character $\chi$ of conductor $q \ge 10^{10}$, every real zero $\beta \in (0,1)$ of $L(s,\chi)$ satisfies $(1 - \beta)\log q \ge \delta_0$, where $\delta_0 = 81/(5888\, e^{3/2}) = 3.069\ldots \times 10^{-3}$.

Banks proves this assertion himself. His second assertion, covering all $q \ge 3$, relies on a computation by Lu, Zaman and Zhao that he did not verify. We use only the first assertion.

**Corollary.** If $h(-d) = 1$, then $d < 10^{10}$.

*Proof.* The conductor of $\chi_{-d}$ is $d$. If $d \ge 10^{10}$, Lemma 2 gives a real zero with

$$(1-\beta)\log d < \frac{2\log d}{\sqrt d} \le \frac{2 \log(10^{10})}{10^5} < 4.61 \times 10^{-4},$$

which contradicts $(1-\beta)\log d \ge 3.0 \times 10^{-3}$. ∎

The margin is a factor of about $6.5$, so the argument would survive a correspondingly weaker $\delta_0$.

## 4. The range $163 < d < 10^{10}$

By Lemma 1, a class-number-one $d > 163$ satisfies $d \equiv 3 \pmod 8$ and $(-d/p) = -1$ for every odd prime $p < (d+1)/4$. `sieve_h1.py` checks that no $d$ in $(163, 10^{10})$ has this property.

1. **CRT stage.** For each prime $p \in \{3, 5, \ldots, 23\}$, keep the $(p-1)/2$ residues $r$ with $-r$ a non-residue mod $p$. Combined with $d \equiv 3 \pmod 8$, this leaves $142{,}560$ classes modulo $8 \cdot 3 \cdot 5 \cdots 23 = 892{,}371{,}480$, and $1{,}597{,}562$ candidates in $(163, 10^{10})$. Every candidate has $(d+1)/4 > 41$, so all of these primes are below $(d+1)/4$ and the conditions genuinely apply.
2. **Filtering stage.** For each prime $29 \le p < 2000$, remove the candidates with $p < (d+1)/4$ and $(-d/p) \ne -1$.

No candidate survives. The run takes about $0.4$ s.

## 5. Verification

| Script | What it checks | Result |
|---|---|---|
| `sieve_h1.py` | Section 4, for $163 < d < 10^{10}$ | 0 survivors |
| `check_small.py` | Brute-force Lemma 1 condition for $d < 2 \times 10^6$; class numbers by counting reduced forms for $\lvert D\rvert < 2 \times 10^4$ | $\{3, 11, 19, 43, 67, 163\}$; $h = 1$ exactly for $\lvert D\rvert \in \{3, 4, 7, 8, 11, 19, 43, 67, 163\}$ |
| `check_lemma.py` | The expansion of Section 2.1 against $\zeta(s)L(s,\chi)$ computed independently, for $d \in \{11, 19, 43, 67, 163\}$ and several $s$; real zeros of $L$ on $(\tfrac12, 1)$ | Agreement to $10^{-29}$; no real zero for these small $d$ (see Remark 2.5) |

The brute-force range of `check_small.py` overlaps the sieve range on $(163, 2 \times 10^6)$, and the two agree.

## 6. Inputs used

- Banks [1], Theorem 1.1, first assertion.
- The Chowla–Selberg expansion of $\zeta_K$ (Poisson summation and the Fourier transform of $(t^2+c^2)^{-s}$) [3].
- Elementary facts: $\zeta(u) = u/(u-1) - u\int_1^\infty \{t\}t^{-u-1}dt$, monotonicity of $\Gamma$, $\psi$ and $K_\nu$, and the class number formula.
- Rabinowitsch's criterion [5] and a finite computation.

Not used: Heegner points, modular forms, Baker's method, Gross–Zagier, or the Lu–Zaman–Zhao computation.

## References

1. W. D. Banks, *On the OpenAI proof that Landau–Siegel zeros do not exist*, arXiv:2610.15483 (2026).
2. B. Green, *On Landau–Siegel zeros*, hexagon:2610.00311v2 (2026). B. Kerr, *A simpler proof that there are no Landau–Siegel zeros*, arXiv:2610.21099 (2026).
3. S. Chowla and A. Selberg, *On Epstein's zeta-function*, J. Reine Angew. Math. 227 (1967), 86–110.
4. OpenAI, *Uniform exclusion of Landau–Siegel zeros*, manuscript dated 1 October 2026.
5. G. Rabinowitsch, *Eindeutigkeit der Zerlegung in Primzahlfaktoren in quadratischen Zahlkörpern*, J. Reine Angew. Math. 142 (1913), 153–164.
6. T. Tao, *An experimental AI-math challenge: using tenth discriminant bounds to explore arguments ruling out Siegel zeroes*, What's new, 20 October 2026, and the comment thread there.
