# Reverse-Cuckoo leakage against regular and stationary LPN

This note studies the ideal independent-hash Reverse-Cuckoo channel.  It is an
attack analysis, not a proof of the leakage-robust Ring-LPN assumption.  The
implementation additionally reuses a Goldreich-feature root across a batch;
that joint channel must be studied separately.

## Assessment

The leakage does not make the syndrome solution space dense.  It removes or
downweights support candidates while leaving the planted solution sparse.  The
security question is therefore whether it makes the sparse solution easier to
*find*, not whether it creates many solutions.

There is no evident direct linear attack.  Conditional on a candidate linear
test being a codeword, its bias is exactly the posterior probability that its
support avoids the hidden noise support.  Reverse-Cuckoo does not reveal noise
coefficients and therefore gives no direct cancellation equation.  However,
the test may be chosen as a function of the leakage.  The linear-test theorem
for Stationary Syndrome Decoding (SSD) does not cover that adaptive choice; it
reduces our linear question to a leakage-aware codeword-avoidance problem.

Finite-$K$ occupancy selection substantially suppresses the soft placement-
multiplicity signal.  It does not suppress the hard statement that the true
support admits a matching.  Reduced exact experiments show that the remaining
signal is overwhelmingly correlation information rather than single-coordinate
information.  This makes coordinate-weighted Prange attacks less concerning,
but leaves joint list decoding, belief propagation, and algebraic solving as
the important open attack classes.

The analogy with regular noise and SSD is encouraging but incomplete:

- In the sparse regime, the best known regular-ISD attacks have the same leading
  asymptotic exponent as ordinary ISD, even though regular supports form a tiny
  fraction of all supports.
- The SSD analysis proves resistance to its linear-test model and finds that an
  adapted XL attack is cheapest on one instance rather than many stationary
  instances.
- Reverse-Cuckoo adds candidate-dependent constraints, Ring-LPN uses structured
  circulant operators, and the real hash descriptors share algebraic features.
  None of those differences is covered by the cited results.

The current evidence supports treating a finite-$K$, $2N$ Reverse-Cuckoo path
as a heuristic LPN-optimized construction, not as an ordinary-DMPF construction
with a negligible statistical error.  Raising the regular weight while keeping
$d=t$ is not presently attractive: the independent-list leakage diagnostic has
grown faster than the nominal regular-Prange exponent in the explored rows.

## 1. Leakage channel

Let $A$ be a set of $n$ active addresses.  A two-choice descriptor $h$ gives
each address one bin in each of two partitions of size $d$.  Let $Z_A(h)$ be
the number of injective placements of $A$ under $h$, and let $U$ be the uniform
descriptor distribution.  Ordinary Reverse Cuckoo has density

$$
P_{\rm raw}(h\mid A)
=U(h)\frac{Z_A(h)}{\mu_n},
\qquad
\mu_n=\mathbb E_U[Z_A(H)].
$$

For the occupancy selector, let

$$
D_A(h)=2n-\left|\{h_0(a):a\in A\}\right|
          -\left|\{h_1(a):a\in A\}\right|
$$

be the endpoint-occupancy deficit and set

$$
a_n(D)=\frac{1}{\mathbb E_U[Z_A(H)\mid D_A(H)=D]}.
$$

Draw $K$ independent raw proposals and select one with probability proportional
to $a_n(D)$.  The exact selected-descriptor likelihood is

$$
P_K(h\mid A)
=P_{\rm raw}(h\mid A)g_{K,n}(D_A(h)),
$$

where

$$
g_{K,n}(D)
=K\,\mathbb E\left[
\frac{a_n(D)}{a_n(D)+\sum_{i=2}^K a_n(D_i)}
\right]
$$

and the $D_i$ are independent deficits under $P_{\rm raw}$.  Equivalently,
for $a=a_n(D)$,

$$
g_{K,n}(D)
=K\int_0^1
\left(\sum_{D'}P_{\rm raw}(D')
u^{a_n(D')/a}\right)^{K-1}du.
$$

The experiment evaluates this one-dimensional integral deterministically.

### Hard support is invariant under finite $K$

All occupancy rates are positive.  Therefore, for every finite $K$,

$$
P_K(h\mid A)>0 \quad\Longleftrightarrow\quad Z_A(h)>0.
$$

Finite $K$ changes likelihood ratios among feasible descriptors, but it does
not change the feasibility language.  In particular, any SAT, Gröbner-basis,
or combinatorial attack that uses only the constraints $Z_A(h)>0$ is completely
unaffected by $K$.

The best channel attainable while retaining precisely this support is the
uniform descriptor conditioned on feasibility,

$$
P_{\rm feas}(h\mid A)
=\frac{U(h)\mathbf 1[Z_A(h)>0]}
       {\Pr_U[Z_A(H)>0]}.
$$

Exact inverse-$Z$ races converge to this channel as $K$ grows.  Occupancy races
need not: even at infinite $K$ they can retain multiplicity variation inside a
fixed deficit class.  Thus the feasibility channel is a lower benchmark for
any positive-rate selector over planted proposals, not the automatic limit of
the occupancy selector.

At $t=d=16$, a uniform two-choice graph is feasible with probability about
$0.938$.  The feasibility-only KL divergence is only about $0.093$ bits per
list, but across the idealized 256 independent lists it still gives a 23.7-bit
exhaustive-ranking diagnostic.  That figure is not a 23.7-bit loss from the
best Ring-LPN attack.

## 2. Exact characterization of linear tests

Work over a field $\mathbb F$ and divide $n=tL$ coordinates into $t$ blocks of
length $L$.  Let the stationary support be
$S=(S_0,\ldots,S_{t-1})\in[L]^t$.  In instance $i\in[q]$, let
$e_i[j,S_j]=\eta_{i,j}$ and all other entries be zero, where the coefficients
$\eta_{i,j}$ are fresh and uniform in $\mathbb F$.  Consider LPN samples

$$
b_i=A_i s_i+e_i.
$$

After seeing the public matrices and leakage $\Lambda$, a linear-test
algorithm chooses vectors $v_i$.  If some $v_i^T A_i\ne0$, the secret $s_i$
makes $\sum_i v_i^Tb_i$ uniform.  The only relevant case is therefore
$v_i^T A_i=0$ for every $i$.

For each regular block define the forbidden coordinates

$$
T_j(v)=\{r\in[L]:\text{there exists }i\text{ with }v_i[j,r]\ne0\}.
$$

For any nontrivial additive character $\chi$ of $\mathbb F$, independence and
uniformity of the fresh coefficients give the exact conditional character
bias

$$
\left|\mathbb E\left[
\chi\!\left(\sum_i v_i^Tb_i\right)
\middle|A_1,\ldots,A_q,\Lambda
\right]\right|
=\Pr[S_j\notin T_j(v)\text{ for all }j\mid\Lambda].
\tag{1}
$$

Equation (1) is the cleanest description of what linear attacks can obtain
from Reverse-Cuckoo leakage:

1. the leakage cannot reveal or linearly cancel fresh coefficient values;
2. it can help only by guiding the choice of codeword supports away from the
   stationary locations;
3. all expansions enter through the union $T_j(v)$, so a fixed test does not
   benefit merely from having more independently weighted instances.

Without leakage, the $S_j$ are independent and uniform.  For one nonzero
codeword with block weights $d_j$, Equation (1) becomes

$$
\prod_{j=0}^{t-1}\left(1-\frac{d_j}{L}\right)
\leq\left(1-\frac{d_{\rm dual}}{n}\right)^t,
$$

which is the SSD linear-test bound when the code family has dual distance at
least $d_{\rm dual}$.  With $\Lambda$, the posterior of $S$ is no longer a
product distribution and the inequality does not follow.  A theorem for our
setting would need a *conditional dual-avoidance* statement: no efficient
algorithm given $(A_1,\ldots,A_q,\Lambda)$ can find nonzero codewords for which
the right side of Equation (1) is noticeably larger than its no-leakage value.

If coefficients are restricted to $\mathbb F^*$, a hit contributes a small
signed character term rather than exactly zero; over a large field it is of
size $1/(|\mathbb F|-1)$.  The binary fixed-coefficient case is materially
different and must use the binary version of the SSD analysis.  The exact
zero-on-hit formula above matches the stationary PCG convention in which
fresh coefficients are sampled from all of $\mathbb F$.

### What this says about concrete linear attacks

| attack | effect of $\Lambda$ | current conclusion |
|---|---|---|
| Fixed Fourier/parity test | Only changes the posterior avoidance probability in Equation (1). | No direct coefficient leakage. |
| BKW or linear sample combinations | A successful combination still corresponds to a codeword/test vector; $\Lambda$ may guide which combinations are sought. | Covered by conditional dual avoidance only when the choice is output-independent. |
| Prange/ISD | The information set can be selected from $\Lambda$ before the syndrome test. | Not covered by the published SSD bound; this is the primary linear-algebraic attack target. |
| Coordinate-weighted ISD | Uses one-coordinate posterior marginals as reliability scores. | Reduced experiments show small and decreasing marginal signal. |
| Correlation-aware ISD | Selects whole tuples or lists using matching/collision factors. | Plausible and not ruled out; most important attack to implement. |
| Output-span tests across expansions | Coefficients are eliminated after observing outputs. | Outside the fixed linear-test model.  Fresh, sufficiently independent code/mask instances are essential. |

The SSD paper gives a concrete warning about the last row.  If every stationary
instance reuses the same generator matrix, then $B=G'_S C$ has column span
contained in the fixed hidden submatrix $G'_S$; after roughly $t$ outputs an
output-span test distinguishes the samples.  SSD avoids this by requiring a
large, weakly correlated code family.  In the Ring-LPN application, fresh ring
masks are therefore security-critical.  Their matrices are circulant rather
than uniformly random, so the absence of an analogous cross-expansion span or
module attack must be audited rather than assumed.

## 3. Relation to regular SD and SSD

Esser and Santini's regular-ISD analysis gives the strongest positive analogy.
For $t=o(n)$, ordinary ISD and regular-ISD have the same leading asymptotic
running time, differing in lower-order terms.  Their sparse concrete rows
include $(n,k,t)=(245760,245460,15)$ and $(40960,40660,20)$, where ordinary
ISD and their best regular-ISD estimates are respectively about 127 and 126
bits.  Regularity can help significantly in denser regimes, so the conclusion
is parameter-dependent.  The key lesson is that the negative logarithm of a
structured-support event is not itself a bit-security loss.

SSD is closer to the stationary setting.  It reuses the $t$ support locations
across $q$ independent matrices while refreshing all nonzero coefficients.
Its linear-test theorem gives exactly the no-leakage specialization of
Equation (1).  It also analyzes a nonlinear XL attack by adding, across all
instances, quadratic equations that force different coordinates in the same
regular block never to be simultaneously nonzero.  For the paper's
$n=2^{16}$, $k=2^{15}$, $t=69$ experiment, the estimated XL-hybrid work grows
from 128 bits at $q=1$ to 132, 144, 156, and 168 bits at
$q=2,16,128,1024$.  The best attack in that model targets one instance.

This is meaningful evidence for robustness to stationary reuse, but it has
three limitations here:

1. SSD uses random or very large code families; Ring-LPN gives structured
   circulant multiplication matrices.
2. SSD has regular-support equations only; Reverse-Cuckoo adds a fixed,
   support-dependent matching CSP.
3. Its XL estimates rely on Hilbert-series/maximal-rank assumptions.  The 2025
   soundness analysis shows that the assumptions fail in full generality but
   supports important degree-two/asymptotic regimes and proves separate
   polynomial-time regimes.  Higher witness-degree behavior remains less
   settled.

The SSD experiment should therefore be cited as a successful attack attempt
that found no stationary advantage, not as evidence that all nonlinear attacks
are harmless.

## 4. Nonlinear attack surface

### 4.1 Likelihood-aware support search

Fix the corrupt party's offsets $x_{a,u}$ and one hidden regular polynomial
$y=(y_0,\ldots,y_{t-1})$.  For public descriptors
$h_{a,k}$, the ideal posterior factor is

$$
W(y)=\prod_{a=0}^{P-1}\prod_{k=0}^{t-1}
P_{\rm channel}\!\left(
h_{a,k}\middle|\{x_{a,u}+y_{k-u}:u\in[t]\}
\right).
$$

For independent ideal descriptors, the four hidden polynomials factor into
four independent $t$-variable problems.  At the production point, each such
problem has 16 variables in a domain of size $L=N/16$ and 64 arity-16 factors.
The full 256 descriptors constrain the same 64 hidden offsets, not 256
independent secrets.

For two choices, $Z$ is determined by connected components of a bipartite
multigraph: a tree contributes its vertex count, a unicyclic component
contributes two, and an overfull component contributes zero.  Most variance in
$\log Z$ is explained by endpoint collisions.  Expanding the likelihood into
pair collisions turns the problem into a dense-looking but weak, repulsive
pairwise constraint system on the offsets.  This suggests the following
attack families:

- belief propagation or mean-field inference using collision potentials;
- MCMC/coordinate descent with exact $Z$ as the acceptance score;
- branch-and-bound using feasibility and partial Hall obstructions;
- collision-weighted regular-ISD list generation and nearest-neighbor merging;
- a hybrid that guesses a few offsets, conditions all incident factors, and
  invokes an ordinary syndrome decoder on the remainder.

The syndrome must be included when evaluating any of these.  Leakage-only
candidate rank and Shannon information are diagnostics, not attack costs.

### 4.2 Quadratic system

Use one-hot variables $z_{b,j,r}$ for the hidden offsets and matching variables
$\rho_{\lambda,u,v}$ for the chosen bin of every item in every descriptor.
The constraints are quadratic:

- Boolean and one-hot equations for $z$ and $\rho$;
- $\rho_{\lambda,u,v}(1-C_{\lambda,u,v}(z))=0$ for candidate incidence;
- $\rho_{\lambda,u,v}\rho_{\lambda,u',v}=0$ for injectivity;
- expanded error equations $E_{b,j,r}=\alpha_{b,j}z_{b,j,r}$;
- the linear Ring-LPN syndrome equations in the expanded error variables.

At $P=4,t=d=16,B=256$, the direct encoding has about $4N$ support one-hot
variables, 131,072 matching variables, and 983,040 pairwise injectivity
equations before adding Boolean, one-hot, incidence, and syndrome equations.
For $N=2^{20}$, a naive degree-two Macaulay matrix over all support variables
would have on the order of $10^{13}$ quadratic monomials.  This rules out a
naive dense XL implementation; it is not evidence against sparse, structured,
elimination, SAT, or message-passing attacks.

Eliminating the matching variables leaves Hall constraints over subsets of
items.  That raises the apparent degree/arity and discards the placement-count
multiplicity.  Keeping the matching witnesses preserves degree two but makes
the system large and gives each support solution a number of lifts equal to
its Reverse-Cuckoo likelihood.  Finite $K$ changes those multiplicities and
soft likelihoods, but leaves the zero-versus-nonzero solution set unchanged.

### 4.3 Stationary expansions

Each additional expansion reuses $z$ and the matching constraints while adding
fresh coefficient variables and another public syndrome.  A compact support-
variable formulation therefore becomes more overdetermined with $q$, even
though the SSD XL formulation became more expensive as $q$ grew.  Small-scale
Gröbner/SAT experiments should test this exact combined system.  In particular,
we should not infer from the SSD table that extra expansions are harmless once
the fixed Reverse-Cuckoo CSP is present.

### 4.4 Actual shared-feature channel

The ideal model treats each descriptor as an independent random function.  The
implementation uses Goldreich-style quadratic features and reuses a public
feature root across the batch.  `rev_cuckoo_goldreich.py` now samples that
shared root, samples the revealed coefficient vectors uniformly from their
programmed affine fibers, and exposes the induced full-domain tables to the
attack.  At reduced sizes it also scores the exact rank-weighted coefficient-
vector likelihood rather than treating feature-rank dependence as an abort.

Conditioned on full row rank for every placed feature submatrix, the candidate
likelihood reduces exactly to the ideal placement-multiplicity likelihood.
The shared root nevertheless correlates the descriptor tables away from the
planted rows.  At the implementation-valid power-of-two setting
$P=2,t=3,d=4,L=8,K=4$, two independent comparisons found shared-minus-ideal
joint-rectangle differences of $+0.008$ and $-0.004$ bit.  Marginal and pair-
collision differences were also about one hundredth of a bit or less.  The
experiment therefore detects no stable additional shared-feature leakage at
this toy size.  Raising the matching size to $t=d=4$ gives the same result:
two independent comparisons found joint differences of $+0.035$ and $+0.005$
bit, with a pooled difference of about $+0.02$ bit below the sampling
resolution.  This is encouraging evidence as matching complexity grows, but
it is not a joint pseudorandomness proof.

At that $t=d=4$ point, the two-seed pooled exact joint gains for
$K=1,2,4,8$ are respectively $2.893$, $1.834$, $1.207$, and $0.810$ bits.
The reductions from each doubling are $1.059$, $0.627$, and $0.397$ bit,
while the corresponding production cost model is $1.00$, $1.35$, $1.73$,
and $2.50$ times the raw baseline.  This makes $K=4$ a reduced-attack Pareto
knee, but not a security threshold.

At $P=4,t=d=16,L=4096,K=4$, the matching-aware production screen still reached
only $1/256$ late rectangle survival and no positive gain.  The sampled batch
had no rank-deficient systems among 512 planted affine solves.  Production
scoring currently uses the good-rank likelihood after ten items, so a global
rank-aware or algebraic attack remains open.  The shared-feature channel is no
longer wholly untested, but a proof still needs a joint pseudorandomness/rank
argument or a stronger global attack study.

The identical production screen was also run for raw $K=1$ and occupancy
$K=2$.  Their local pair slopes are larger in magnitude, as expected, but both
rectangles again collapse to $1/256$ late survival.  Thus this production
attack does not justify $K=4$: it fails for every tested $K$.  Parameter choice
must currently use the reduced exact curve and MPC cost, explicitly under the
leakage assumption, until a stronger production optimizer supplies a usable
bound.

We also tested a full-batch conditional selector that removes the factorwise
approximation entirely.  Posterior SMC particles supply contexts, and every
candidate value of one offset is scored by the product of all exact list
likelihoods before the next rectangle restriction.  Independent SMC
populations score the resulting rectangle.  It matches the previous selector
on eight exact toys (mean paired difference about $-0.001$ bit), but both four-
and sixteen-context production runs at raw $K=1$ have zero late survival.
Thus even complete local Hall/matching information does not produce a useful
single rectangle on this instance.  A stronger attack must cover multiple
posterior modes and charge each attempted rectangle, or use a non-rectangular
information-set strategy.

A charged four-mode rectangle portfolio does not improve the exact toys.  We
also implemented a non-rectangular one-offset attack.  It selects an
exact-width value set from a posterior-particle histogram and uses an
independent SMC population to estimate its mass with exact full-conditionals;
the baseline for a width-$s$ set is $s/L$.  Independent scoring is necessary:
using the same few contexts for selection and evaluation spuriously reported
7--10 bit singleton gains after maximizing over $tL$ noisy candidates.

The holdout estimator matches exact enumeration within 0.03 bit on a reduced
width-two test.  At production `P=4,t=d=16,L=4096`, a width-256 set captures
mass `0.0569` for occupancy $K=4$, below the uniform baseline `0.0625`, hence
the demonstrated gain is zero.  At $K=1$ it also finds no positive gain, but
the minimum SMC ESS is below 0.09 in both independent populations, so that
case is a sampler failure rather than affirmative security evidence.  The
current evidence therefore keeps $K=4$ as a cost/heuristic knee, not a proven
security threshold, and supplies no defensible production claim that $K=1$ is
broken.

We therefore measured posterior entropy directly instead of requiring two
samplers to agree on one mode.  If $W_\Lambda(x)$ is the full transcript
likelihood, SMC estimates the Shannon loss
$E_{x\mid\Lambda}[\log_2 W_\Lambda(x)]-\log_2 E_U[W_\Lambda]$ and the Renyi-2
loss $\log_2 E_U[W_\Lambda^2]-2\log_2 E_U[W_\Lambda]$.  The latter subtracts
from the 192-bit prior offset space to give effective collision support.  The
estimator agrees with exact enumeration within 0.01 bit on a reduced test.

For production shared-Goldreich transcripts, residual collision entropy is
about 134--148 bits for $K=1$, 155--160 bits for $K=2$, and 161--165 bits for
$K=4$.  An ideal-hash screen is slightly more conservative: 130--134 bits for
$K=1$ and 149--154 bits for $K=2$.  Thus $K=1$ does leave an enormous and very
sparse posterior, matching the central heuristic, but its worst observed ideal
margin over 128 bits is only about two bits.  $K=2$ captures most of the support
improvement of $K=4$ at modeled cost 1.37 rather than 1.73 times raw.

Collision entropy gives a formal leakage-only list bound: if $H_2$ is the
posterior collision entropy, any set carrying posterior mass $\epsilon$ has
size at least $\epsilon^2 2^{H_2}$.  It does not rule out an algebraic decoder
that combines the descriptor with the LPN equations, so these measurements are
support evidence rather than a standard-assumption reduction.

The first collision estimates above use $E_{p_\Lambda}[W]$ from a fixed
posterior population and may under-sample heavy modes.  We therefore continue
SMC from target $W$ to $W^2$, adding each second-copy descriptor factor in
quarter-power steps.  This gives `0.553` bit on an exact toy whose true Renyi-2
loss is `0.547`.  At production size, two 1024-particle shared-$K=1$ runs leave
only `138.1` and `131.5` collision bits, and one lower-budget run reaches
`127.7`.  Conversely, three shared-$K=2$ runs leave `151.6--156.5` bits with
powered minimum ESS above `0.90`.  A tempered ideal $t=18,K=1$ stress test
leaves only `128.6` bits after correcting to fixed $N$.  Tail-sensitive support
estimation therefore favors $t=16,K=2$; $K=1$ remains structurally sparse but
does not retain a credible security margin.

## 5. Reduced finite-$K$ experiments

`rev_cuckoo_regular_noise.py` now samples and scores the raw, finite-$K$
occupancy, and feasibility-only channels exactly at reduced sizes.  The tables
report information about one hidden regular polynomial.  `P=2`; hence there
are $2t$ descriptor factors.  The Prange column is a one-shot rectangular
avoidance gain based only on one-coordinate marginals.

For $t=d=3$, $L=8$ (512 candidates, six lists, 512 trials):

| channel | full information | marginal information | correlation information | marginal Prange gain, 6/8 excluded |
|---|---:|---:|---:|---:|
| raw | 0.496 | 0.132 | 0.365 | 1.041 bits |
| occupancy $K=2$ | 0.215 | 0.048 | 0.167 | 0.547 bits |
| occupancy $K=4$ | 0.123 | 0.021 | 0.102 | 0.292 bits |
| occupancy $K=8$ | 0.094 | 0.011 | 0.083 | 0.175 bits |
| feasibility only | 0.080 | 0.004 | 0.076 | 0.071 bits |

For $t=d=4$, $L=6$ (1,296 candidates, eight lists, 256 trials):

| channel | full information | marginal information | correlation information | marginal Prange gain, 3/6 excluded |
|---|---:|---:|---:|---:|
| raw | 0.830 | 0.315 | 0.515 | 1.272 bits |
| occupancy $K=4$ | 0.227 | 0.054 | 0.173 | 0.463 bits |
| occupancy $K=8$ | 0.165 | 0.027 | 0.138 | 0.278 bits |
| feasibility only | 0.120 | 0.010 | 0.109 | 0.102 bits |

The numerical conclusion is stable across both toys:

- $K=4$ removes most raw information at moderate setup cost;
- $K=8$ gives a smaller second reduction;
- coordinate marginals shrink faster than total information;
- the residual approaches a correlated feasibility/CSP signal rather than
  independent noisy hints about each location.

The experiments are much too small to extrapolate a production bit-security
loss.  They do identify the attack that an extrapolation must model: a joint
correlation-aware decoder, not a coordinate-only reliability decoder.

### Exact basic Ring-LPN experiment

The support-only experiments above do not measure Ring-LPN security loss.  We
therefore extended the exact coefficient-aware enumerator to the finite-$K$
channel.  It samples two regular polynomials $e,f$ with independent nonzero
coefficients, samples a uniform public $r$, computes the complete negacyclic
syndrome $f+re$ in $\mathbb F_q[X]/(X^N+1)$, and reveals all $4t$ product-list
descriptors induced by the corrupt party's two known regular polynomials.  For
each product list it uses the implementation's grouping by output block,
applies the negacyclic sign to wrapped terms, deduplicates equal local
addresses, adds their product coefficients, and removes a row when that sum
is zero.  The enumerator then considers every regular support and every
nonzero coefficient assignment for $(e,f)$ having that same syndrome.

Thus the measured conditional information is exactly
$I((e,f);\Lambda\mid r,f+re)$ in the reduced model.  We report a separate
computational metric for exhaustive support decoding: order all regular
support guesses by their exact descriptor likelihood and test each guess
against the Ring-LPN syndrome, solving for its coefficients.  Its expected
rank gain is the base-two reduction relative to a random ordering.  This is a
direct leakage-aware guess-and-check attack, not a support-entropy proxy.

Across four exact toys---$(q,t,L)=(3,2,2),(3,2,3),(5,2,2),(3,3,2)$---the
results are:

| selector | conditional information (range) | support-enumeration rank gain (range) |
|---:|---:|---:|
| $K=1$ | $0.083$--$0.182$ bits | $0.296$--$0.560$ bits |
| $K=2$ | $0.023$--$0.115$ bits | $0.174$--$0.391$ bits |
| $K=4$ | $0.001$--$0.035$ bits | $0.033$--$0.151$ bits |
| $K=8$ | $0.002$--$0.031$ bits | $0.032$--$0.132$ bits |

Small negative point estimates were rounded to zero because conditional mutual
information is nonnegative and those estimates are within one standard error
of zero.  Raw Reverse Cuckoo and occupancy selection with $K=1$ are the same
channel distribution.  The detailed estimates, standard errors, candidate
counts, and seeds are recorded in
`rev_cuckoo_ring_lpn_product_results.csv`; the support-only control is in
`rev_cuckoo_ring_lpn_channel_results.csv`, and the reproducible driver is
`rev_cuckoo_ring_lpn_channel_sweep.py`.  Moving from the support-only control
to coefficient-aware product deduplication does not materially increase the
observed raw-channel attack: the largest $K=1$ support-rank gain moves from
$0.534$ to $0.560$ bits.  This comparison is deliberately conservative about
cancellation effects because the toy fields $q=3,5$ make zero sums far more
common than the implementation's 31- or 64-bit fields.  Indeed, exact
enumeration shows that a zero cancellation removes at least one deduplicated
row in $9.4\%$--$36.7\%$ of the toy product lists, depending on the
configuration.  The detailed row counts are in
`rev_cuckoo_product_cancellation_results.csv`.

This changes the interpretation of the earlier support-only results.  At
these exact sizes, both the additional information conditioned on the
syndrome and the gain of the explicit support-enumeration attack are sub-bit.
In particular, the experiment does not show that $K=1$ is broken, and it gives
no empirical reason to call $K=4$ a security threshold.  It also does not
establish production security: the enumerated rings have dimensions four or
six, use independent ideal hash functions, and cannot reveal a
large-parameter transition, a speedup of a non-enumerative Ring-LPN decoder,
or the effect of the shared Goldreich feature root.  The next basic-LPN task
is therefore to scale a syndrome-conditioned decoder, rather than further
refining the support-only entropy model.  SSD remains a separate later step.

#### Repeatable basic Ring-LPN Prange attack

We next instantiate the basic $P=2$ regular syndrome matrix
$[I\mid M_r]$.  A Prange iteration retains exactly $L/2$ coordinates in each
of the $2t$ regular blocks, so its uniform support success is $2^{-2t}$.  The
experiment keeps the complete coefficient-deduplicated product channel above
and audits the rank of every selected square matrix over $\mathbb F_{65537}$;
all selected matrices in the recorded sweep were full rank.

The marginal and joint selectors use the exact reduced posterior.  The pair
selector is scalable: it fits the descriptor likelihood by public collision
features, runs mean-field inference in polynomial time, and uses exact
posterior enumeration only to evaluate its success in the experiment.  At
$(t,L)=(2,4),(3,4),(3,8),(4,4)$, the best deterministic pair rectangle gives
one-shot gains of respectively $0.749,1.216,1.345,1.324$ bits for raw $K=1$.
For $K=4$ these become $0.216,0.431,0.463,0.584$ bits.  The exact joint
rectangle oracles are larger, reaching $0.924,1.617,1.625,2.232$ bits for
$K=1$.  These are real one-shot probability gains, but repeating the same
deterministic rectangle is not a decoding algorithm.

To obtain a repeatable attack, let the inferred block weights be $p_i$ and
sample every $L/2$-subset $C$ with probability proportional to
$\prod_{i\in C}p_i^\beta$.  The inclusion probability of coordinate $i$ is
computed exactly as
\[
 \Pr[i\in C]
 =\frac{p_i^\beta e_{L/2-1}(p_{-i}^\beta)}
        {e_{L/2}(p^\beta)},
\]
so independent iterations and the complete success curve can be evaluated
without observing rare support hits.  Optimizing $\beta$, the scalable pair
attack's gain in work to reach $50\%$ success is only $0.195,0.180,0.138,0.080$
bits for raw $K=1$ at the four points above.  For $K=4$ the corresponding
values are $0,0.038,0.045,0.060$ bits.  Its gain under the stricter average
expected-iterations metric is at most $0.076$ bits for $K=1$ and $0.020$ bits
for $K=4$ in this sweep.

This is the first attack here with an explicit repeatable information-set
distribution.  It demonstrates that the leakage is algorithmically usable,
but currently finds no material basic Ring-LPN work-factor reduction.  The
large gap between deterministic one-shot mass and repeated-attack work comes
from hard instances: aggressively favoring one rectangle makes supports
outside it much more expensive.  The detailed results are in
`rev_cuckoo_basic_ring_lpn_isd_results.csv`, generated by
`rev_cuckoo_basic_ring_lpn_isd_sweep.py`.  These reduced $t\leq4$ experiments
do not rule out growth at the production noise weight or a better joint
randomized selector.

#### Production-shaped basic Ring-LPN attack

The reduced conclusion does not survive increasing the regular noise weight.
We implemented a production-size repeatable selector that avoids the
$O(L^2)$ elementary-symmetric-polynomial calculation.  In each block it sorts
the pair mean-field beliefs and divides the ranking into $L/P$ groups.  Group
$i$ contains ranks $i,i+L/P,\ldots,i+(P-1)L/P$, and the selector chooses
exactly one coordinate from every group.  Within a group $G$, it selects
coordinate $j$ with probability
\[
 \frac{p_j^\beta}{\sum_{k\in G}p_k^\beta}.
\]
Thus every iteration contains exactly $L/P$ coordinates from each block,
every coordinate retains positive probability, and $\beta=0$ gives every
coordinate inclusion probability exactly $1/P$.  The planted-support
probability is consequently the exact regular-Prange baseline $P^{-Pt}$,
although the complete subset distribution is stratified rather than uniform.
The belief computation and sampling are polynomial in $P,t,L$ and do not
enumerate support candidates.

The experiment samples two known and two hidden regular polynomials, forms all
$4t$ coefficient-deduplicated negacyclic product descriptors, and evaluates
the planted support under independent ideal two-choice hash functions.  At
fixed ring dimension $N=tL=65536$, the measured support-stage work gains are:

| $t$ | $L$ | $K$ | trials | gain, fixed $\beta=.25$ | gain, fixed $\beta=.5$ | $50\%$ gain, best fixed $\beta$ |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 4096 | 1  | 24 | 0.087 | 0.139 | 0.279 |
| 16 | 4096 | 4  | 16 | 0.030 | 0.045 | 0.089 |
| 32 | 2048 | 1  | 8  | 0.884 | 1.041 | 2.423 |
| 64 | 1024 | 1  | 6  | 11.837 | 14.987 | 21.618 |
| 64 | 1024 | 4  | 6  | 8.584 | 10.024 | 12.959 |
| 64 | 1024 | 8  | 4  | 8.541 | 10.879 | 13.881 |
| 64 | 1024 | 16 | 4  | 6.508 | 7.261 | 11.117 |

The small $K=8$ versus $K=4$ reversal is not evidence that more proposals
leak more: those rows have only four and six trials.  The important result is
the scale.  For $t=64$, every recorded instance has positive one-iteration
log gain at $\beta=.25$, including all four $K=16$ instances.  At
$t=16,K=1$, the fixed-$\beta=.5$ expected-work estimate is 0.139 bits; a
nonparametric bootstrap over the 24 saved instances gives an approximate 95%
interval $[0.025,0.253]$.  At $t=16,K=4$, the analogous 0.045-bit estimate has
an interval spanning zero.  The $t=64$ sample counts are sufficient to show a
multi-bit effect, but not to give narrow confidence intervals.

This scaling has a direct security interpretation.  Basic $P=2$ regular
Prange has support success $2^{-2t}$, so $t=16$ supplies only a 32-bit support
exponent whereas $t=64$ supplies 128 bits.  Increasing $t$ adds baseline ISD
work, but it also increases both the number and arity of the leaked product
constraints.  Consequently the low-noise observation that $K=1$ has a
sub-bit concrete effect cannot justify $K=1$ at a 128-bit basic Ring-LPN
point.  The occupancy race also debiases too slowly: at $t=64$ the fitted
pair-collision coefficient moves from $-0.583$ for $K=1$ to only $-0.536$,
$-0.497$, $-0.458$, $-0.407$, and $-0.355$ for respectively
$K=4,8,16,32,64$.

We also tested a correlated-center construction.  It forces a
collision-model center into each block and fills the remaining $L/2-1$
positions uniformly, which reduces exactly to uniform regular Prange for a
uniform center.  Exact posterior sampling gains 0.089 expected-work bits in a
$(t,L,K)=(3,8,1)$ diagnostic, showing that this interface can use joint
correlation.  A polynomial coordinate-ascent center is already neutral by
$t=16,L=256$ and is much weaker than the balanced-half sampler at $t=64$.

These are support-stage gains.  Reduced large-field experiments found every
selected structured matrix full rank, but the $N=65536$ scorer does not yet
audit the rank event.  The experiment also intentionally remains in the
independent ideal-hash channel and makes no claim about shared implementation
features.  Subject to those caveats, it is a concrete leakage-aware decoding
attack, not an exhaustive-rank or Shannon-information heuristic.  The driver
is `rev_cuckoo_basic_ring_lpn_isd_production.py`; aggregate and per-instance
CSVs use the prefix `rev_cuckoo_basic_ring_lpn_isd_production_`.

The generalized selector permits a direct test of the implementation-shaped
$P=4,t=16,L=4096$ configuration.  It also has baseline support work
\[
 Pt\log_2P=4\cdot16\cdot2=128\text{ bits},
\]
but keeps $L/4$ coordinates in each of 64 blocks and leaks all
$P^2t=256$ product descriptors.  Six-instance ideal-channel screens give:

| $P$ | $t$ | $L$ | $K$ | gain, fixed $\beta=.25$ | gain, fixed $\beta=.5$ | $50\%$ gain, best fixed $\beta$ |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 16 | 4096 | 1 | 0.399 | 0.609 | 1.224 |
| 4 | 16 | 4096 | 4 | 0.168 | 0.257 | 0.463 |

These results sharply distinguish two ways of obtaining a 128-bit support
exponent.  The basic $P=2,t=64$ point loses roughly 10--15 support-stage bits
for $K=1$--4 under the same scalable attack, whereas the $P=4,t=16$ point
loses substantially less than one expected-work bit.  This supports using
more low-weight polynomials rather than increasing the product weight inside
a basic instance.  It does not prove the $P=4$ construction secure: the
sample is small, the full-rank caveat remains, and a stronger joint or
non-Prange decoder could perform better.  The saved $P=4$ rows use the prefix
`rev_cuckoo_ring_lpn_isd_production_p4_`.

#### Sequential collision-conditioned selector

The preceding selector samples every block independently after computing its
pair mean-field beliefs.  This discards a usable part of the joint signal:
once a subset has been selected in one block, its realized collision labels
change which coordinates are attractive in the remaining blocks.  We
therefore implemented a correlated, repeatable selector.  It orders the
blocks by marginal entropy and samples them sequentially.  Before sampling a
new block, it adjusts that block's log weights by the difference between the
collision-label mass of the subsets selected so far and the corresponding
mean-field mass.  It then applies the same fixed-cardinality stratified
sampler.  At $\beta=0$ the corrected ordering cannot change any coordinate's
conditional inclusion probability, so the construction recovers the exact
$P^{-Pt}$ support baseline.

Directly observing complete successes would require about $2^{128}$ trials.
The evaluator instead uses a bootstrap particle filter.  At each block it
samples the attack conditioned on including the planted coordinate,
multiplies by that exact conditional inclusion probability, and resamples
the partial selections according to their inclusion weights.  The planted
support is used only by this rare-event evaluator; the attack's transition
rule sees only the leaked descriptors, public multipliers, and its own prior
selections.  At the reduced point $P=4,t=3,L=16,K=1$, increasing the particle
count from 16 to 64 changes the sequential expected-work gain only from
$0.786$ to $0.799$ bits (with $0.807$ at 32 particles), providing a basic
stability check.

At $P=4,t=16,L=4096$ and fixed $\beta=.5$, the six-instance estimates are:

| $K$ | independent expected-work gain | sequential expected-work gain | independent $50\%$ gain | sequential $50\%$ gain |
|---:|---:|---:|---:|---:|
| 1 | 0.609 | 0.977 | 0.811 | 1.087 |
| 4 | 0.257 | approximately 0.37--0.40 | 0.323 | 0.434 |

The $K=4$ range reflects particle noise in one of the six instances: the
16-particle aggregate is $0.403$ bits, while a 64-particle re-evaluation of
that instance moves the aggregate to about $0.369$ bits.  The conclusion is
therefore not that the independent approximation was exact, but that the
first explicitly correlated repeatable attack still extracts less than one
expected-work bit from $K=4$.  Raw $K=1$ is measurably worse, at about one
bit under the same attack.

This remains a pair-collision surrogate.  It does not propagate exact
matching counts or Hall obstructions, and it does not optimize a global
non-Prange decoder.  Those are the natural directions for a stronger attack.
The implementation and per-instance rows are in
`rev_cuckoo_ring_lpn_sequential_isd.py` and
`rev_cuckoo_ring_lpn_sequential_isd_p4_*.csv`.

We next tested the matching/Hall direction directly while retaining a
repeatable, randomized information-set distribution.  For each hidden
polynomial, the new selector first runs SMC against the complete product of
exact Reverse-Cuckoo matching likelihoods.  Its posterior particles therefore
encode feasibility, component structure, and placement multiplicity across
all $Pt$ descriptor factors.  During an information-set iteration, the
selector blends the pair message with the posterior-particle histogram for
the next block.  After selecting a fixed-cardinality subset, it conditions and
resamples the posterior population, so later blocks follow the posterior mode
chosen by earlier subsets.  A nonzero pair component preserves full support.

At $P=4,t=3,L=16$, we can replace SMC by exact enumeration of the complete
matching posterior.  With 1,024 exact-posterior contexts and a $10\%$ matching
blend, the $K=1$ expected-work gains are $0.946$ bits for pair-only and $0.990$
bits for matching guidance.  This $0.044$-bit difference is below the roughly
$0.05$-bit null fluctuation observed when both evaluators use the pair-only
distribution.  For $K=4$, the corresponding gains are $0.0755$ and $0.0752$
bits.  Thus even an exact reduced posterior supplies no detectable improvement
through these sequential conditional histograms.

At $P=4,t=16,L=4096,K=4$, two production-shaped instances give pair-only
one-instance log gains of $0.118$ and $0.785$ bits, versus $-0.067$ and $0.625$
bits with a $10\%$ global matching blend.  The two-instance expected-work
aggregates are respectively $0.413$ and $0.237$ bits.  The posterior SMC has
minimum ESS fractions $0.263,0.303$, minimum unique fractions $0.188,0.215$,
and mutation acceptance $0.225,0.233$.  These diagnostics are imperfect but
do not indicate total particle collapse; moreover, the exact reduced control
has the same negative conclusion.  Since an attacker can always discard the
matching guide, this experiment does not lower the pair-only estimate.  It
instead suggests that the remaining exact matching information is multimodal
or otherwise poorly compressed into conditional coordinate histograms.

The sharper raw $K=1$ posterior is harder to sample.  Uniform-prior SMC on the
same two production instances has minimum ESS fractions $0.0063,0.0157$ and
minimum unique fraction $0.0195$ in both cases, so its matching-guided scores
are not credible posterior diagnostics.  We therefore tested a valid but more
aggressive attack distribution $q_{\rm pair}(x)W_{\rm match}(x)$.  Particles
start from the pair marginals, exact matching factors are introduced
sequentially, and coordinate Metropolis proposals use the same pair
distribution; the proposal ratio cancels the $q_{\rm pair}$ factor.  This
raises the minimum ESS to $0.043,0.064$ and uniqueness to $0.055,0.063$.
Pair-only one-instance gains are $0.573,1.447$ bits, while matching-guided
gains are $0.512,1.492$ bits.  Their two-instance expected-work aggregates are
$0.945$ and $0.920$ bits.  Thus matching guidance wins on one seed and loses
on the other, with no aggregate improvement.  The remaining diagnostics are
still weak, but the exact reduced control reaches the same neutral conclusion.

These experiments close the conditional-histogram refinement of regular
Prange for now: independent, collision-conditioned, exact-posterior, and
pair-proposal matching guides have all been tested, and the strongest measured
loss remains about one bit for $K=1$ and $0.4$ bits for $K=4$.  They do not
close attacks that retain whole posterior particles or couple matching
witnesses directly to the Ring-LPN equations.

This is still not an exhaustive global attack.  A non-Prange decoder could use
whole posterior particles, matching witnesses, or the Ring-LPN equations
directly rather than reduce them to block subsets.  The matching-particle
selector also has substantially more preprocessing than pair-only, a cost not
charged in the figures above.  The driver and recorded rows are
`rev_cuckoo_ring_lpn_matching_isd.py` and
`rev_cuckoo_ring_lpn_matching_isd_p4_*.csv`.

#### Equation-coupled four-list decoder

We next abandoned rectangular information sets and coupled the complete
matching likelihood directly to a Ring-LPN equation.  The reduced decoder uses
$P=4$ regular unit-coefficient error polynomials over a prime-field
negacyclic ring.  A field of size 17 makes the syndrome essentially unique at
tractable dimensions, avoiding the thousands of aliases in the earlier small
binary structural toy.  For each hidden polynomial it enumerates all $L^t$
regular supports, computes the polynomial's syndrome contribution, and assigns
either a uniform score, the fitted pair-collision score, or the product of all
$Pt$ exact Reverse-Cuckoo likelihood factors incident to that polynomial.

The decoder combines polynomials zero and one on the left and polynomials two
and three on the right.  Each Cartesian pair list is generated best-first by
the sum of its two leakage scores, without materializing the complete sorted
list.  The two sides are alternately inserted into hash tables and joined on
the complete Ring-LPN syndrome.  All recorded runs recover the planted
solution.  The baseline uses independently randomized pair-list orders.

We report two work models.  The optimistic model counts one operation per
per-polynomial candidate constructed and one per pair-list entry emitted.  The
factor-charged model additionally counts every complete descriptor-likelihood
evaluation as one operation.  It still omits sorting, heap, memory, and the
larger constant of exact matching, so it should not be read as a production
bit cost.  The average-expected-work gains of the exact selector are:

| $t$ | $L$ | $K$ | trials | optimistic gain | factor-charged gain | paired bootstrap 95% interval |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 8 | 1 | 16 | 2.139 | 1.773 | $[1.279,2.508]$ |
| 3 | 8 | 4 | 16 | 0.809 | 0.652 | $[0.287,1.082]$ |
| 4 | 4 | 1 | 32 | 2.515 | 1.474 | $[1.152,1.787]$ |
| 4 | 4 | 2 | 32 | 1.048 | 0.581 | $[0.301,0.871]$ |
| 4 | 4 | 4 | 32 | 0.925 | 0.491 | $[0.200,0.770]$ |
| 4 | 4 | 8 | 32 | 0.600 | 0.244 | $[-0.037,0.502]$ |

The exact matching score materially outperforms the pair surrogate.  For
example, at $(t,L)=(4,4)$ the pair selector's optimistic expected-work gains
are $1.204$ bits for $K=1$ and $0.451$ bits for $K=4$, versus $2.515$ and
$0.925$ bits for exact matching.  This is the first decoder here to turn the
residual higher-order matching structure into a clear equation-coupled work
reduction.  Finite-$K$ selection suppresses but does not eliminate the signal:
$K=4$ retains about half a charged bit in both tested geometries, whereas the
$K=8$ interval includes zero.

This result does not extrapolate directly to production.  The decoder first
enumerates $L^t$ candidates per polynomial, which is $4096^{16}=2^{192}$ at
the implementation-shaped point, before its two pair lists even begin.  The
prime-field unit-coefficient syndrome is a reduced structural model, the
descriptors use independent ideal hashes, and the cost model is optimistic.
The result instead identifies the next scaling problem precisely: replace
complete per-polynomial enumeration by an inner beam, Stern, or Wagner decoder
that uses partial collision bounds and evaluates exact matching factors when
they close.  The driver, deterministic tests, and raw rows are
`rev_cuckoo_ring_lpn_list_attack.py`,
`test_rev_cuckoo_ring_lpn_list_attack.py`, and
`rev_cuckoo_ring_lpn_list_attack_p4_*.csv`.

#### Blockwise inner-list search

We next replaced the complete per-polynomial enumeration by a bounded
blockwise beam.  For a partial assignment of the unknown regular offsets, the
pair-collision surrogate can already count every collision between active
items whose variables have been assigned.  Adding another offset can only add
collisions.  Since the fitted collision slope is negative, the partial score
is therefore an upper bound on the score of every completion.  At each level
the beam keeps the best partial assignments under this bound.  On completing
all $t$ offsets, it keeps twice the requested final list width, evaluates the
complete matching likelihood on those candidates, and reranks them exactly.
The incremental score was checked candidate-by-candidate against the original
complete pair surrogate on exhaustive small instances, including cases in
which duplicate product indices collapse in the active set.

The metric here is deliberately different from the preceding decoder's
stopping time.  At final width $B$, a uniform list contains one planted
polynomial with probability $B/L^t$.  We measure the planted-support retention
of each of the four leakage-guided lists.  To compare at equal success
probability, we take the geometric mean of the four measured coordinate
retention rates and choose the uniform width with that same retention.  The
work comparison charges the beam's generated children, one operation per
complete descriptor-likelihood factor, and two exhaustive $B^2$ outer pair
lists.  We do **not** interpret the fixed-$B$ success-probability ratio as an
independent-restart exponent: the leakage instance is fixed, so rerunning a
deterministic top-$B$ search produces the same lists.

At the 30-times-larger reduced support space $(t,L)=(5,6)$, a final width of
$B=256$, pair oversampling by two, and beam width 512 give:

| $K$ | trials | uniform retention | effective exact retention | equivalent uniform width | equal-success charged gain | bootstrap 95% interval |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 64 | 0.0329 | 0.2963 | 2304 | 5.80 bits | $[5.12,6.29]$ |
| 4 | 128 | 0.0329 | 0.1162 | 904 | 3.10 bits | $[2.17,3.73]$ |
| 8 | 128 | 0.0329 | 0.0814 | 633 | 2.08 bits | $[0.79,2.88]$ |

For $K=4$, the pair-only beam gives 2.12 charged bits
($[0.84,2.82]$), so exact reranking contributes about one additional bit.
A cutoff sweep confirms that the result is not monotone in list width:
$B=128,256,512$ give exact point estimates of approximately
$1.94,3.10,2.69$ bits, respectively, under the same accounting.  The middle
width best balances posterior concentration against beam and quadratic merge
cost in this reduced geometry.

This is meaningful evidence that avoiding $L^t$ enumeration does not erase
the leakage signal, and it reverses the earlier reduced impression that
$K=8$ makes the equation-coupled advantage statistically invisible.  It is
still not a production attack.  The four coordinate rates are combined under
an independence approximation; direct four-way success is too rare for the
current trial counts to estimate.  The model does not execute the outer
Ring-LPN join, it charges each exact graph-likelihood factor as a unit
operation, and a width-$W$ beam requires roughly $tLW$ child expansions.
With $L=4096$ even this polynomial-in-$L$ inner search can be prohibitive at
the widths needed for useful retention.  The next step is a global
Stern/Wagner-style search that uses syndrome projections while partial
supports are merged, followed by experiments at larger $t$, field size, and
$K$.

The driver and raw rows are `rev_cuckoo_ring_lpn_inner_beam.py` and
`rev_cuckoo_ring_lpn_inner_beam_t*_k*.csv`.

#### Global syndrome projections

We then tested whether a Stern/Wagner-style syndrome projection can prune the
outer search rather than merely improve four independent lists.  The first
construction is deliberately free of planted advice.  It projects every
left-pair and right-pair syndrome onto $r$ public field coordinates, computes
the posterior mass of all compatible projection buckets from the four exact
leakage scores, and searches buckets in decreasing product mass.  Polynomial
projection masses are combined by cyclic convolution over
$\mathbb F_q^r$; the implementation uses a dense multidimensional FFT and
charges $2q^r\log_2(q^r)$ convolution operations.  Within each searched
bucket, complete pair entries are ordered by exact leakage score and joined on
the full syndrome.  No true projection value is supplied to the attack.

At $(t,L,q,r)=(4,4,17,2)$, 16 paired trials give:

| $K$ | flat exact gain | projected exact gain | projected minus flat | bootstrap 95% interval for difference |
|---:|---:|---:|---:|---:|
| 1 | 1.335 bits | 0.792 bits | -0.543 bits | $[-1.074,-0.026]$ |
| 4 | 0.385 bits | 0.195 bits | -0.190 bits | $[-0.608,0.232]$ |
| 8 | 0.057 bits | -0.007 bits | -0.065 bits | $[-0.488,0.354]$ |

Thus this global projection gives no improvement over the existing flat
best-first pair merge.  It often locates the correct coarse projection bucket
earlier than a uniform bucket order, but a high-scoring pair in a lower-mass
bucket is delayed behind every entry of the earlier buckets.  The flat merge's
candidate-level order is strictly finer and avoids the FFT preprocessing.  A
one-coordinate $K=4$ control reaches the same conclusion: projected exact and
flat exact gains are approximately 0.10 and 0.48 bits over 32 trials.

There is also a generic obstruction to pushing this projection inside each
polynomial.  Assigning $r$ projected syndrome coordinates separately to $P$
polynomials requires guessing $P-1$ shares, hence $q^{r(P-1)}$ share choices.
Each share filters one $L^t$ polynomial list by only $q^r$, so enumerating the
share space costs a factor $q^{r(P-2)}$ more candidate work than the original
lists; for $P=4$ this is $q^{2r}$.  Allocating only one share between the two
outer polynomial pairs costs $q^r$ and filters pair work by $q^r$, which is
work-neutral.  The posterior bucket experiment above tests the remaining hope
of searching those pair shares nonuniformly, and does not find an improvement.

This does not rule out all Stern/Wagner attacks.  A useful construction would
need additional representations or a way to reuse partial-list work across
many projection shares; simply guessing shares or imposing intermediate zero
conditions silently discards the planted sparse solution and repays the
filtering gain in repetitions.  The driver, tests, and recorded rows are
`rev_cuckoo_ring_lpn_projected_bucket.py`,
`test_rev_cuckoo_ring_lpn_projected_bucket.py`, and
`rev_cuckoo_ring_lpn_projected_bucket_t4_l4_r*_k*.csv`.

#### Single-syndrome q-ary regular Stern

We next instantiated the enumeration-based regular-ISD construction of Esser
and Santini for one complete Ring-LPN syndrome.  Here the ambient vector has
$n=PtL$ coordinates, the parity-check matrix has $N=tL$ rows, and the error
has one unknown nonzero field element in each of $Pt$ blocks of length $L$.
The information set $J$ contains $k+\ell=(P-1)N+\ell$ coordinates.  Its
integer block widths $v_b$ are balanced, differ by at most one, and sum to
$k+\ell$.

For even $p$, the Stern iteration asks for exactly $p/2$ planted errors in
each half of $J$.  If $H$ is one half of the regular blocks, its leaf-list
size is
\[
  \mathcal L_H
  = \sum_{\substack{S\subseteq H\\|S|=p/2}}
      \prod_{b\in S} v_b(q-1).
\]
The factor $q-1$ is essential: unlike binary regular noise, each selected
Ring-LPN error has an unknown nonzero coefficient.  Matching on $\ell$ field
equations leaves an expected $\mathcal L_0\mathcal L_1/q^\ell$ collisions.
We charge the maximum of the two leaf lists and the collision list, while
omitting polynomial matrix-elimination and memory costs.  This is favorable
to the attack.

Reverse-Cuckoo leakage changes only the iteration-success probability.  For
each of the 64 regular blocks, the production mean-field selector gives the
probability that the planted coordinate lies in that block's chosen part of
$J$.  Two exact Poisson-binomial calculations then give the probability of
$p/2$ hits in each fixed half.  We average the inverse success probability
over leakage instances and charge the list work on every iteration.  This is
a one-syndrome attack: there is no stationary-support reuse or amortization.

At $(P,t,L,q)=(4,16,4096,65537)$, searching even $p\leq8$ and
$0\leq\ell\leq16$ gives:

| channel | trials | best $\beta$ | best total work | best forced $p>0$ | forced work | forced minus best |
|---:|---:|---:|---:|---:|---:|---:|
| uniform | analytic | -- | $2^{128.000}$ | $(p,\ell)=(2,3)$ | $2^{147.419}$ | 19.419 bits |
| occupancy $K=1$ | 6 | 1.0 | $2^{127.123}$ | $(2,3)$ | $2^{146.526}$ | 19.403 bits |
| occupancy $K=4$ | 4 | 1.0 | $2^{127.360}$ | $(2,3)$ | $2^{146.788}$ | 19.428 bits |
| occupancy $K=8$ | 4 | 0.5 | $2^{127.841}$ | $(2,3)$ | $2^{147.262}$ | 19.421 bits |

Every global optimum is $p=0,\ell=0$, namely the leakage-weighted regular
Prange attack already studied above.  The leakage advantage for the best
nontrivial Stern row is nearly the same sub-bit advantage as for $p=0$ and
does not bridge the roughly 19.4-bit list penalty.  Thus ordinary regular
Stern remains a valid nonlinear attack, but it is not competitive at the
current large-field parameters under this measured posterior model.  This
conclusion does not follow from the SSD linear-test theorem and does not rule
out representation-based MMT/BJMM: those algorithms require a separate
$q$-ary representation-count analysis.

The estimator, deterministic checks, and full search grids are
`rev_cuckoo_ring_lpn_stern.py`,
`test_rev_cuckoo_ring_lpn_stern.py`, and
`rev_cuckoo_ring_lpn_stern_p4_k*_t16.csv`.

#### q-ary MMT/BJMM representations

We then generalized the representation count from the binary regular-ISD
analysis to field-valued regular errors.  Write $W=Pt$ for the number of
regular blocks and let a target vector have regular weight $a$.  In a
representation as the sum of two regular vectors of equal weight $b$, an
active target coefficient may either be assigned wholly to one addend or
split into two nonzero coefficients.  If $s$ active blocks are split, each
has $q-2$ ordered coefficient choices.  An inactive block may contain a
cancelling pair at any of its $v$ coordinates with any of $q-1$ nonzero
values.  If there are $u$ such blocks, equal addend weight requires
$2b=a+s+2u$.  The complete ordered representation count is therefore
\[
 R_q(a,b)=
 \sum_{\substack{0\leq s\leq a,\ 0\leq u\leq W-a\\
                  2b=a+s+2u}}
 \binom as(q-2)^s
 \binom{a-s}{(a-s)/2}
 \binom{W-a}{u}\bigl(v(q-1)\bigr)^u.
\]
For $q=2$, all $s>0$ terms vanish and this reduces exactly to the binary
formula of Esser and Santini.  We also checked the formula by exhaustive
enumeration of every equal-weight regular decomposition in a small
$q=3$ instance.

The one-level MMT model sets $\ell_x=\log_q R_q(p,p_x)$, enumerates the two
regular addend lists, and joins them on the remaining $\ell-\ell_x$ field
equations.  The two-level BJMM model applies the same count recursively with
$\ell_y=\log_qR_q(p_x,p_y)$.  We use fractional filter lengths, sum all
decomposition types with the same addend weight, and omit polynomial,
elimination, and memory-access costs.  Thus the model is deliberately
optimistic for the attack.

At $(P,t,L,q)=(4,16,4096,65537)$, an exhaustive search through
$0\leq\ell\leq128$ gives:

| method | best nontrivial parameters | iteration | inverse success | total work | gap above Prange |
|---:|---:|---:|---:|---:|---:|
| Stern | $(p,\ell)=(2,3)$ | 32.585 bits | 114.834 bits | 147.419 bits | 19.419 bits |
| one-level MMT | $(p,p_x,\ell)=(2,2,5)$ | 32.585 bits | 113.860 bits | 146.445 bits | 18.445 bits |
| two-level BJMM | $(p,p_x,p_y,\ell)=(2,3,2,6)$ | 47.562 bits | 113.861 bits | 161.423 bits | 33.423 bits |

MMT recovers about one bit relative to Stern, principally by using all
placements of the two information-set errors rather than requiring one in
each fixed half.  Large-field coefficient splitting supplies many
representations, but cannot remove the $2^{32.585}$ leaf list containing one
unknown position and one unknown field coefficient.  The extra BJMM level is
counterproductive because its first merge grows to $2^{47.562}$.

As a lower bound on any benefit from the May--Ozerov final nearest-neighbor
step, we repeated the entire parameter search while declaring the final merge
free.  The optimum changes only from 146.4445 to 146.4432 bits: the leaf list,
not the final merge, is already dominant.  Nearest-neighbor acceleration of
that merge therefore cannot make this regular-representation family
competitive at the current parameters.

The measured Reverse-Cuckoo advantage in the corresponding $p=2$ Stern event
is at most 0.893 bit among $K=1,4,8$.  Even granting that complete gain to MMT
would leave about $2^{145.55}$ work for $K=1$, versus the measured
$2^{127.12}$ leakage-weighted $p=0$ attack.  This transfer is only a
diagnostic because the MMT success event omits Stern's fixed-half condition;
the roughly 18.4-bit separation is nevertheless too large for the distinction
to matter at the present precision.

The estimator, exhaustive-count tests, and compact search summary are
`rev_cuckoo_ring_lpn_bjmm.py`,
`test_rev_cuckoo_ring_lpn_bjmm.py`, and
`rev_cuckoo_ring_lpn_bjmm_p4_t16.csv`.  These calculations close the standard
single-syndrome regular Prange/Stern/MMT/BJMM progression under the current
random-matrix cost model.  They do not exclude an attack exploiting the
circulant Ring-LPN matrix or the full higher-order matching leakage.

#### Sparse-factor reduction and projected leakage

The random-matrix progression is not the controlling attack for the reducible
ring used by the protocol.  For a factor $X^n+a$ of $X^N+1$, reduction modulo
that factor folds every original exponent into its residue modulo $n$.  With
$P$ regular error polynomials, each containing $t$ points, the expected total
reduced weight is
\[
  \bar w_n=Pn\left(1-(1-1/n)^t\right).
\]
At the implementation-shaped point $(P,t,n)=(4,16,128)$ this is $60.383$.
This is exactly the $(c,w,i,w_i)=(4,64,7,60)$ row in the Ring-LPN parameter
analysis.  Its conservative statistical-decoding estimate is
\[
  n\left(\frac{Pn}{n-1}\right)^{\bar w_n},
\]
or $2^{128.450}$ field operations before rounding the reduced weight.  This
calibration resolves the apparent discrepancy with the full-length $128$-bit
Prange calculation: the relevant structured attack works on a length-$512$,
degree-$128$ reduced instance, not on the original $65536$-row matrix.

The reduced code has length $Pn=512$, dimension $(P-1)n=384$, and a dual
space of dimension $n=128$.  Consequently, every chosen set $Z$ of $n-1=127$
zero coordinates is realizable by a nonzero dual codeword: imposing
$z_j=0$ for $j\in Z$ gives only 127 homogeneous constraints on the
128-dimensional dual space.  If the reduced error support $T$ is contained
in $Z$, the resulting parity check is noiseless.  This turns global zero-set
selection into a concrete statistical-decoding attack under the paper's
favorable assumption that arbitrary dual checks may be precomputed.  Solving
the small homogeneous system on demand is also polynomial, although charging
that solve would only make the attack slower.

We project the scalable pair-collision posterior into this factor.  For block
$j$, its inferred distribution on the $L=4096$ original offsets is summed by
$x\mapsto(jL+x)\bmod128$.  The score $s_{p,r}$ is the inferred probability
that reduced coordinate $r$ of polynomial $p$ is occupied.  Flattening the
four score vectors gives 512 weights, and the attack samples a global
127-subset according to
\[
  \Pr_\beta[Z]\ \propto\ \prod_{j\in Z}s_j^\beta.
\]
For planted reduced support $T$, the exact containment probability is
\[
 \frac{\left(\prod_{j\in T}s_j^\beta\right)
 e_{127-|T|}((s_j^\beta)_{j\notin T})}
 {e_{127}((s_j^\beta)_{j\in[512]})}.
\]
This Rao--Blackwellizes over the randomized zero set, so no rare support hit
must be observed.  At $\beta=0$ it is the exact hypergeometric probability.
The earlier selector that forced 32 zeros per polynomial is retained only as
a stricter control.

Four production trials at each channel give:

| channel | uniform per-instance work | best $\beta$ | work gain at $\beta=1$ | signal gain at $\beta=1$ | work gain at $\beta=8$ | signal gain at $\beta=8$ |
|---:|---:|---:|---:|---:|---:|---:|
| occupancy $K=1$ | $147.927$ bits | 0 | $-0.022$ | $-0.052$ | $-0.236$ | $-0.456$ |
| occupancy $K=4$ | $147.927$ bits | 0 | $-0.015$ | $-0.016$ | $-0.151$ | $-0.151$ |

Here ``work'' is $\log_2\mathbb E[1/p]$ plus the factor-$128$ check cost,
while ``signal'' is $-\log_2\mathbb E[p]$ plus that cost.  Both averages choose
$\beta=0$, and all eight paired production trials have a negative ratio for
every positive $\beta$.  Thus the tested Reverse-Cuckoo signal does not
survive the modulo-128 folding in a form useful even to the directly
realizable global dual-check selector.  This closes the obvious
leakage-adaptive version of the paper's statistical decoder under the current
pair-mean-field posterior.  A better attack would need to preserve
higher-order matching information through factor reduction, rather than just
reweight reduced coordinates.

We tested that residual direction directly at the highest-leakage $K=1$
point.  The pair-marginal distribution was used as an SMC proposal, and all
$Pt$ complete Reverse-Cuckoo matching likelihoods were then applied to each
candidate.  The resulting aggressively leakage-weighted heuristic target is
$q_{\rm pair}W_{\rm match}$; it is not claimed to be an exact posterior.  We
evaluated two ways of converting its final
particles into a degree-128 dual check.  First, particle residue-occupancy
marginals were blended into the score $s_j$ above.  Second, a stronger joint
selector planted the complete reduced support of one particle for each of the
four error polynomials, then filled the remaining positions of $Z$ uniformly.
For a planted prediction $U$, its exact conditional success is
\[
 \frac{\binom{512-|U|-|T\setminus U|}
                   {127-|U|-|T\setminus U|}}
      {\binom{512-|U|}{127-|U|}},
\]
and we average this value exactly over the empirical particle populations.
Thus the joint experiment preserves every within-polynomial support
correlation and does not estimate a rare support hit by sampling.

With two production trials, 128 particles, and two mutation moves per factor,
the uniform selector costs $145.975$ bits on these sampled weights.  Planting
supports drawn independently from the pair marginals costs $151.262$ bits, a
$5.286$-bit loss.  Planting the complete matching-particle supports costs
$180.923$ bits, a $34.948$-bit loss.  For the marginal attack, the optimum over
mixes $0,.25,.5,.75$ and $\beta=0,.25,.5,1,2,4,8$ is again no guidance:
mix $0$ and $\beta=0$.  Even the weakest positive matching blend
$(\mathrm{mix},\beta)=(.25,.25)$ loses $0.252$ bit in expected work.

This is negative evidence against this direct higher-order factor attack, not
a proof that matching leakage is useless.  The complete-factor SMC is
difficult: the minimum ESS and unique-particle fractions were $0.0255$ and
$0.0234$, respectively.  Particle impoverishment can only make this heuristic
selector worse, so a fundamentally better Hall- or matching-aware decoder is
not excluded.  Nevertheless, neither exact factor messages projected to
residue marginals nor complete particle supports provide any measured gain at
the worst tested leakage point.  This closes the concrete attack families we
set out to test; further security claims should be stated as an assumption
against efficient leakage-aware decoding rather than as an empirical proof.

The absolute baseline also clarifies how conservative the published estimate
is.  For fixed reduced weight $W$, exact uniform zero-set containment costs
\[
 n\,\frac{\binom{Pn}{n-1}}
          {\binom{Pn-W}{n-1-W}}.
\]
At $W=60$ this is $2^{146.774}$, whereas replacing the without-replacement
probability by $((n-1)/(Pn))^W$ gives the paper's approximately $2^{127.679}$
integer-weight estimate.  The paper explicitly presents its statistical
decoder as a conservative lower bound, and this gap shows that the bound is
very favorable to the attacker at $W/n\approx0.47$.  The four-trial absolute
values are not tail estimates and should not replace a calculation over the
complete noise law; their paired leakage ratios are the relevant output.

The factor reduction exposes a separate, non-cuckoo implementation issue.
The parameter analysis explicitly rejects noise whose weight after reduction
by the smallest relevant sparse factor is below its expectation.  The current
`RingLpnTriple::genDpf` implementation instead samples each of the 64 regular
positions once and performs no such test.  The exact occupancy recurrence
gives
\[
  \Pr[W<\lceil\bar w_{128}\rceil]
  =\Pr[W<61]=0.498086.
\]
Averaging the paper's optimistic independent-coordinate bound over the
unconditioned weight distribution gives $122.482$ bits; conditioning on
$W\ge61$ raises it to $130.572$ bits.  Averaging the exact hypergeometric dual
check instead gives $136.431$ bits unconditioned and $150.526$ bits after
conditioning.  Both models show that the lower-weight tail matters even
though their absolute exponents differ substantially.  The rejection accepts
with probability $0.501914$, or $1.992$ samples on average, which is the
roughly one-bit entropy loss anticipated in the Ring-LPN paper.  It should be
added independently of any decision about Reverse Cuckoo: the position holder
can fold its locally sampled supports modulo 128 and resample before
constructing the DPF.  Nonzero factor scalars do not change support collisions,
and coefficient cancellation is negligible over the large implementation
fields.

> **Implementation TODO (security).** Add weak-factor rejection to
> `RingLpnTriple::genDpf` before product positions or DPF state are derived.
> Parameterize the smallest relevant 1-sparse factor degree and minimum folded
> weight. For the current $(P,t)=(4,16)$ row, repeatedly sample the complete
> regular support, fold absolute positions modulo 128 separately within each
> polynomial, and accept only when the total distinct folded weight is at
> least 61. Add a deterministic support-folding test and a sampling test for
> the predicted $0.501914$ acceptance probability. This TODO is independent
> of the selected DMPF construction and must therefore apply to both Reverse
> Cuckoo and Waterfall configurations.

The driver, deterministic tests, summaries, and paired trials are
`rev_cuckoo_ring_lpn_factor_attack.py`,
`test_rev_cuckoo_ring_lpn_factor_attack.py`, and
`rev_cuckoo_ring_lpn_factor_k{1,4}*.csv`.  The matching-aware follow-up is in
`rev_cuckoo_ring_lpn_factor_matching.py`,
`test_rev_cuckoo_ring_lpn_factor_matching.py`, and
`rev_cuckoo_ring_lpn_factor_matching_k1.csv`.

#### Structured-rank audit

The support scorer above does not construct the $65536$-square selected
matrix.  We therefore sampled the selector's actual fixed-cardinality sets and
audited
\[
 [I\mid M_{r_1}\mid\cdots\mid M_{r_{P-1}}]_C
\]
by finite-field Gaussian elimination at scaled dimensions.  Each comparison
uses the same independently sampled public multipliers for three selections:
ordinary uniform regular Prange, the stratified selector at $\beta=0$, and
the leakage-guided selector at $\beta=.5$.  The multipliers are independent
of the support leakage, as in the Ring-LPN game.

For $P=4,t=16,L=16$ ($N=256$), the full-rank counts are:

| rank field | samples/selector | uniform | stratified $\beta=0$ | guided $\beta=.5$ |
|---:|---:|---:|---:|---:|
| $\mathbb F_{17}$ | 400 | 371 | 387 | 375 |
| $\mathbb F_{257}$ | 1000 | 997 | 996 | 998 |
| $\mathbb F_{65537}$ | 400 | 400 | 400 | 400 |

Every observed failure had deficiency one.  In the paired
$\mathbb F_{257}$ audit, uniform alone failed three times and guided alone
failed twice; there were no shared failures.  Doubling to $L=32$
($N=512$) over $\mathbb F_{257}$ gave respectively $99/100,100/100,100/100$
full-rank matrices.  The implied guided-versus-baseline rank correction lies
between approximately $-0.046$ and $+0.016$ bits in these amplified-failure
experiments and has no detectable adverse direction.  Thus there is currently
no empirical reason to discount either the independent or sequential
$P=4,K=4$ support-stage gain for rank.

We repeated the audit with subsets sampled by the sequential correlated
selector.  Over $\mathbb F_{257}$ the uniform, stratified-$\beta=0$, guided,
and sequential selectors respectively produced $499,498,499,498$ full-rank
matrices out of 500.  Over $\mathbb F_{17}$ the corresponding counts were
$381,372,375,381$ out of 400.  Every failure again had deficiency one.  The
sequential selector therefore shows no detectable rank penalty relative to
uniform selection in these deliberately small-field tests.

This is not a production-rank proof.  Dense elimination cannot reach
$N=65536$, and a formal statement would require either a nonvanishing-minor
argument for the selected block-circulant matrix or a fast black-box rank test
using negacyclic convolution.  The audit and raw paired rows are in
`rev_cuckoo_ring_lpn_rank_audit.py` and
`rev_cuckoo_ring_lpn_rank_audit_p4_*.csv`.

## 6. Parameter implications

The current independent-list diagnostics at $t=d=16$ are:

| channel | modeled exhaustive-rank gain over 256 ideal lists |
|---|---:|
| raw | 166.8 bits |
| occupancy $K=2$ | 115.4 bits |
| occupancy $K=4$ | 85.3 bits |
| occupancy $K=8$ | 69.5 bits |
| feasibility only | 23.7 bits |

These are gains over uniform exhaustive enumeration of $L^{64}$ supports, not
losses from a 128-bit decoder.  The uniform expected-rank exponent is roughly
767 bits already at $N=2^{16}$, so even the raw diagnostic remains far above
128 bits on that irrelevant baseline.

Increasing $t$ while retaining a $2N$ table and $d=t$ has not improved this
diagnostic tradeoff.  Raising $t$ from 16 to 18 adds 16 nominal regular-Prange
bits for four hidden polynomials, but the measured occupancy-$K=4$ independent-
list diagnostic rises by roughly 26 bits.  This comparison is not an attack,
but it argues against “add slightly more noise” as the first mitigation.  At
fixed $2N$ expansion, debiasing the channel and analyzing the decoder are more
promising than increasing $t$.

If the goal is a conventional proof from ordinary Ring-LPN/SSD, finite $K$ is
insufficient because hard feasibility remains.  Waterfall Cuckoo or an exact
canonical-witness correction is the appropriate path.  The earlier heuristic
identification of $K=4$ as a general cost knee is superseded by the basic
$P=2,t=64$ attack above: no tested $K\leq16$ makes that 128-support-bit point
quiet.  In contrast, the direct $P=4,t=16$ experiment finds only a 0.257-bit
expected-work gain for $K=4$ under the current scalable selector.  Thus $K$
cannot be chosen independently of $(P,t,L)$; the measured evidence presently
favors the implementation-shaped $P=4$ row, while still falling short of a
concrete 128-bit security claim.

### Regular-Prange experiment

The first decoder experiment is now implemented in
`rev_cuckoo_regular_isd.py`.  It chooses exactly $N$ columns of the structured
matrix $[I\mid M_{r_1}\mid\cdots\mid M_{r_{P-1}}]$, verifies full rank, and
uses the exact conditional posterior mass of the selected regular rectangle as
its one-iteration success probability.  Over $\mathbb F_{65537}$ every sampled
matrix in the recorded rows was full rank.

At the production-shaped polynomial ratio $P=4$ but reduced regular weights,
the measured gains are:

| $(t,L)$ | channel | marginal selector | joint rectangle selector |
|---:|---|---:|---:|
| $(3,8)$ | raw | 7.43 bits | 9.78 bits |
| $(3,8)$ | occupancy $K=4$ | 2.15 bits | 3.61 bits |
| $(3,8)$ | occupancy $K=8$ | 1.37 bits | 2.27 bits |
| $(3,8)$ | feasibility only | 0.58 bits | 0.83 bits |
| $(3,16)$ | occupancy $K=8$ | 1.07 bits | 1.95 bits |
| $(3,16)$ | feasibility only | 0.48 bits | 0.80 bits |
| $(4,8)$ | occupancy $K=4$ | 4.13 bits | 6.54 bits |
| $(4,8)$ | occupancy $K=8$ | 2.53 bits | 4.63 bits |
| $(4,8)$ | feasibility only | 1.02 bits | 1.48 bits |

This establishes that leakage can improve a genuine regular-Prange iteration;
the gain is not confined to exhaustive support ranking.  It also confirms that
$K=8$ materially reduces the gain and that the hard feasibility floor remains
measurable.  The joint selector still enumerates $L^t$ support tuples, so these
numbers are attack diagnostics rather than net work-factor reductions.  The
full method and all reduced rows are recorded in
`rev_cuckoo_regular_isd.md`.

We then replaced the leakage-aware choice by a polynomial pair-collision
selector.  It fits the exact per-list log likelihood by the public number of
within-partition collision pairs and runs damped mean-field inference on the
induced offset constraints.  The exact posterior is used only to score its
chosen rectangle in the reduced experiments, not to construct it.  At
`P=4,t=3,L=8`, its gains are 7.52 bits for raw Reverse Cuckoo, 2.22 bits for
occupancy `K=4`, and 1.33 bits for occupancy `K=8`.  At `t=3,L=16`, the last
two fall to 1.64 and 1.05 bits.  At `t=4,L=8`, they rise to 4.17 and 2.35
bits.  Thus a demonstrated scalable heuristic does use the leakage, but its
current concrete advantage is only a few bits in these toys.  It recovers
little of the additional advantage seen by the exponential joint-rectangle
oracle.

Feasibility-only leakage has zero fitted pair slope because all feasible
graphs have the same likelihood in that channel.  The remaining feasibility
advantage is therefore a genuinely higher-order target involving Hall
obstructions, components, or matching counts.  This is both the main residual
attack surface and a reason not to interpret the exact joint numbers as an
efficient attack.

The scalable selector can now be run directly at the representative
$(t,L)=(16,4096)$ point using complete-factor aggregation.  A direct,
unbiased coordinate diagnostic shows the selected quarter containing 43.2\%
of planted coordinates for $K=4$ and 34.6\% for $K=8$ at $L=64$.  The signal
then decays toward the 25\% baseline: at $L=1024$ the corresponding small-run
estimates are 28.5\% and 26.6\%, while at $L=4096$ they are 27.7\% and 29.7\%.
The production-size runs have only four trials and do not establish a
nonzero coordinate advantage.  In particular, their reversed $K$ ordering
shows that more sampling is needed.

We also attempted to score the full selected rectangle without $L^t$
enumeration using the exact likelihood-ratio identity and defensive importance
sampling.  It reproduces the exact $t=3$ result, but its effective sample
fraction collapses below $10^{-3}$ at $t=16$.  The resulting provisional
30--50-bit gains are discarded as statistically unsupported.  This is an
important distinction: the large-$L$ experiment supports marginal smoothing,
but it neither proves nor refutes a constant joint-correlation loss.  This
failure motivates the sequential Monte Carlo scorer below.

A sequential Monte Carlo scorer introduces the exact $Pt$ list factors one at
a time, resamples, and applies exact Metropolis rejuvenation.  At
$(P,t,L,K)=(4,3,8,4)$ its separate-normalizer version matches the enumerated
per-polynomial rectangle gains within 0.026 bits.  At production parameters,
however, that version has a heavy factor-order tail even when per-stage ESS is
moderate.  In the shared-known-support batch its dominant component changed
from $+5.66$ to $-3.30$ bits when six replicates were added.  We therefore
retract the preliminary positive per-polynomial screen and its speculative
8.6-bit four-polynomial extrapolation.

A more appropriate rare-event version imposes the sixteen rectangle
restrictions sequentially and directly estimates each conditional probability.
It agrees with exact toy gains within about 0.11 bits.  At
$(P,t,L)=(4,16,4096)$, both the independently thresholded mean-field rectangle
and a new pair-collision coordinate-ascent rectangle reach late-stage survival
of at most $1/256$ and sometimes zero, even with diversity refresh and sixteen
exact-target mutation moves per restriction.  Their precise negative gains are
not resolved, but neither gives a positive attack over uniform regular Prange.

Pair-coordinate refinement is not vacuous: on one exact toy it raises true
gain from 0.810 to 0.956 bits.  Its production failure instead indicates that
the exact higher-order matching constraints are poorly represented by the
pair objective.  The batched experiment therefore provides no evidence for a
single-digit $K=4$ loss from the current scalable selectors.  A Hall-,
component-, or matching-aware rectangle optimizer remains the relevant attack
to try.

We therefore implemented sampled exact matching-factor messages.  For fixed
values of the other $t-1$ offsets, varying one offset produces at most $d^2$
endpoint pairs, so memoization evaluates feasibility, component structure, and
placement multiplicity for all $L$ candidates efficiently.  This is a real
higher-order improvement on the toy: two rounds raise exact gain from 0.956 to
1.083 bits.  On the production $K=4$ descriptor, however, the matching-aware
rectangle again reaches only $1/256$-scale late conditional survival and gives
no positive regular-Prange advantage.

The completed attack portfolio now contains marginal reliability, pairwise
mean field, pair-coordinate rectangles, and local exact matching-factor
messages.  None demonstrates a production-point gain.  A global Hall-aware
optimizer or a different non-rectangle decoder remains possible, so this is
evidence for the leakage-robust assumption rather than a reduction to ordinary
Ring-LPN.

## 7. Required cryptanalysis

The next study should proceed in this order:

1. **Ring-structure-aware decoding.** The random-matrix regular
   Prange/Stern/MMT/BJMM progression is now screened through production
   parameters.  The next decoding attack should exploit the block-circulant
   Ring-LPN operator or leakage-adaptive dual codewords rather than add another
   generic representation level.
2. **Global Hall-aware attack.** If further cryptanalysis is desired, move
   beyond the completed local matching-factor messages to a global Hall-aware
   optimizer or a non-rectangle decoder.  Use conditional-splitting SMC for
   evaluation and charge reusable preprocessing against Prange iterations.
3. **Conditional dual avoidance.** On reduced random-code instances, find the
   best leakage-adaptive low-weight dual codeword and measure the multiplicative
   increase in Equation (1), including the selected-matrix rank event.
4. **Pair-potential inference.** Fit the collision surrogate at $t=16$ and test
   BP/MCMC recovery of one hidden 16-offset polynomial before adding the
   syndrome.  This is a cheap screen for a dangerous planted-CSP transition.
5. **Combined algebraic scaling.** Generate the exact stationary
   Ring-LPN--Reverse-Cuckoo quadratic system for tiny parameters and compare
   solving time as the number of expansions grows.
6. **Actual feature family.** Repeat support-ranking and inference experiments
   with the shared Goldreich feature matrices used by the implementation.
7. **Tail metrics.** Measure smooth max-information and high likelihood-ratio
   quantiles; average KL alone is inadequate for a concrete failure claim.

## References

- V. Kolesnikov, S. Peceny, S. Raghuraman, and P. Rindal,
  [Stationary Syndrome Decoding for Improved PCGs](https://eprint.iacr.org/2025/295),
  ePrint 2025/295.  The formal definitions, linear-test theorem, output-span
  warning, algebraic model, and XL tables were cross-checked against Chapter 2
  of Peceny's 2026 dissertation.
- A. Esser and P. Santini,
  [Not Just Regular Decoding: Asymptotics and Improvements of Regular Syndrome Decoding Attacks](https://eprint.iacr.org/2023/1568),
  ePrint 2023/1568 / CRYPTO 2024.
- P. Briaud and M. Øygarden,
  [A New Algebraic Approach to the Regular Syndrome Decoding Problem and Implications for PCG Constructions](https://eprint.iacr.org/2023/176),
  ePrint 2023/176 / EUROCRYPT 2023.
- M. Cueto Noval, S.-P. Merz, P. Stählin, and A. Ünal,
  [On the Soundness of Algebraic Attacks against Code-based Assumptions](https://eprint.iacr.org/2025/415),
  ePrint 2025/415.
