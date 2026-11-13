# Research plan: reliability and cost of Lagrange-multiplier reduced contact models

Status: DRAFT for discussion. Built on the measured results of Chapter 14 of the book (`packages_book/book_v3`) and a literature check done on 13 Nov 2026 (abstracts, plus parts of two full texts; Google Scholar could not be searched).

## 1. Working title
A computable admissibility estimator and a cost model for reduced contact models with Lagrange multipliers.

## 2. Problem and gap
A reduced contact model with multipliers can satisfy its own reduced KKT conditions to 1e-16 and still return a lifted solution that is wrong by orders of magnitude (measured in the book: displacement error 3e2 to 4e2 with a signed SVD dual basis, against 3e-2 to 9e-2 with an NMF basis at the same sizes; lifted multipliers down to about -1.4e2).

What the literature covers (verified from abstracts or text):
- Positive dual bases: NMF (Balajewicz et al. 2016), cone-projected greedy (Benaceur, Ern, Ehrlacher 2020), greedy index selection (Sambataro, Ehrlacher 2026), restriction to a reduced integration domain (Fauque et al. 2018; Le Berre et al. 2022).
- Certified bounds for reduced variational inequalities (Zhang, Bader, Veroy 2014): they assume a bijective constraint operator B and give feasible dual approximations by construction.
- Full-order residual estimators for mixed contact (Hild, Nicaise 2007).
- Conditioning of the projected contact matrix linked to dual precision (Le Berre et al. 2022).

- Full-order estimators built by replacing the multiplier with a post-processed discrete one (Banz, Schönauer, Schröder, arXiv 2609.26001, Sept 2026; variational inequalities of the second kind, friction-type). The device is the same as the lambda+ projection used here, so the novelty claimed is the reduced-order lifted setting, not the technique.
- Greedy positive-cone compression with hyper-reduction for discrete contact (Sambataro, Ehrlacher 2026), which postpones an error indicator to future work (in the context of basis enrichment for inf-sup stability).

What I did not find (to be re-checked before submission): an estimator for the lifted solution of a reduced model with a rectangular constraint matrix; a measured case where the reduced residual passes and the lifted solution fails; one benchmark that compares the dual-basis constructions; a cost model for when rank repair plus reduced-domain enrichment removes the saving.

## 3. Research questions
RQ1. Can a computable a posteriori estimate computed from the lifted pair detect inadmissible reduced solutions that the reduced residual accepts?
RQ2. Which dual basis (NMF, cone-projected greedy, greedy index selection, signed SVD with repair, reduced-integration-domain restriction) is most reliable, at what size?
RQ3. When does the total cost (reduced unknowns after rank repair, reduced integration domain, contact rows) fall below the full-order cost?

## 4. Proposed method
Lifted pair: u~ = V u_r, lambda~ = U lambda_r. Define lambda+ = max(lambda~, 0), lambda- = min(lambda~, 0).

Estimator components:
- eta_sign = ||lambda-||, cost O(contact rows).
- eta_gap = ||(c - B u~)+||, cost O(contact rows).
- eta_comp = |lambda+^T (B u~ - c)|, cost O(contact rows).
- eta_eq = ||K u~ + B^T lambda+ - f||_{K^-1}, needs a full-order operation, so it is a validation or greedy-selection quantity, not an online one.

Online admissibility indicator (candidate, to be tested in E2): eta_adm = eta_sign + eta_gap + sqrt(eta_comp), with each term scaled to a common unit (for instance through the K and K^-1 norms) because force, length and work differ in units. The square root on eta_comp is a definition to test, not a proved property: the claim that all three terms scale linearly with the error is measured in E2. Use as an online circuit breaker: if eta_adm exceeds a tolerance, flag the prediction and trigger basis enrichment or a fall-back to the full-order solver. eta_eq stays an offline certificate (validation and greedy selection).

Target (to prove or to establish empirically): the energy-norm error is bounded by a constant times a combination of these terms, following the constant-sign multiplier idea of Hild and Nicaise. A proved bound is not assumed; if it fails, the paper reports the empirical effectivity instead.

Cost model, two tiers.
Tier 1, size: unknowns after enrichment = reduced displacement columns + unit vectors added to restore rank (14 on the book's chain) + reduced multipliers; compared with full-order unknowns, as a function of contact-zone size and the number of contact rows.
Tier 2, work: measured wall-clock time and fitted scaling per solve, including dense reduced-matrix factorisation (O(r^3)) against sparse full-order factorisation, and the number of active-set or LCP iterations k_iter per step (a poorly conditioned dual basis can need more iterations than the full-order solver, which removes the saving even when r is much smaller than N). The break-even is reported per solver and hardware from measured scaling; no closed-form crossover is assumed.

## 5. Experiments
E1. Reproduce the book's failure on the 324-unknown membrane and the 60-spring chain; compute all four estimator components for every basis and size.
E2. Effectivity: estimate against the true lifted error, for each basis, in three regimes: (i) in-range parameters; (ii) load-magnitude extrapolation (the contact patch grows beyond the training snapshots, where rank deficiency and negative multipliers are expected); (iii) contact-patch shift or sliding (rows of B inactive in training become active, testing whether eta_gap and eta_comp catch unseen active constraints). Also measure the scaling of each term and of eta_adm against the error. Report the failure cases.
E3. Comparison of the five dual-basis constructions on identical snapshots (accuracy, admissibility violation, online cost, offline cost).
E4. Cost model validation: measured against predicted break-even on the rope and on the membrane, then on a 3D case.
E5. Scale-up: 3D Hertzian contact (elastic hemisphere or cylinder on a block; analytical solution for verification), 1e4 to 1e5 unknowns, with the exact (K, B, c, f) exported from a scriptable FEM backend so that the ROM code is identical to the 1D and 2D code. Backend not yet chosen: FEniCSx contact libraries are marked experimental and I could not confirm two-body support; GetFEM was not checked. A one-week feasibility spike must export (K, B, c, f) for a small 3D case before a tool is fixed.
E6. Dynamics (optional): the reduced LCP with BDF2, with the estimator evaluated per step.

## 6. Success criteria (set before running)
- No false negative: every case with true error above a stated tolerance is flagged by eta_sign, eta_gap or eta_comp.
- Effectivity index between 1 and 10 on at least 90 % of cases for the best basis, or a clear statement of where it fails.
- Spearman rank correlation of the estimate with the true error above 0.9 across bases and sizes.
- Predicted break-even within a factor of 2 of the measured one.
- Online indicator eta_adm: (i) scaling: the slope of log eta_adm against log error, fitted over all cases, lies between 0.8 and 1.2, otherwise the form (units or square root) is changed and the change reported; (ii) circuit breaker: with the tolerance set on the training range only, no case in the extrapolation and patch-shift regimes of E2 with true error above the stated tolerance goes unflagged, and the false-alarm rate (flagged cases with true error below tolerance) is reported, with 10 % as the target; (iii) online cost of eta_adm stays below 10 % of one reduced solve step.
- Thresholds above are initial targets set before running and are revised only with the reason recorded.

## 7. Risks
- The bound may not be provable; mitigation: report it empirically.
- eta_eq is not online-cheap; mitigation: use it for validation and greedy selection and state the limitation.
- Overlap with Le Berre et al. on conditioning; mitigation: cite and show what the measured rank repair adds.
- A recent preprint (Sambataro, Ehrlacher 2026) may be followed by an indicator; mitigation: monitor arXiv before submission.
- The 3D benchmark needs a verified contact solver; this is the main cost. Mitigation: the feasibility spike in E5 before committing.
- The unit scaling and the square-root form of eta_adm may not give linear scaling with the error; mitigation: measure it, and report the form that works.

## 8. Timeline (about 8 months)
Months 1-2: literature recheck with Scholar, estimator derivation, E1-E2 on the existing benchmarks.
Months 3-4: E3, the comparison with the four other bases.
Months 5-6: E4-E5, the 3D solver and cost model.
Month 7: E6 if time allows; writing.
Month 8: internal review, submission.

## 9. Target journals
First choice Computer Methods in Applied Mechanics and Engineering or International Journal for Numerical Methods in Engineering if the estimator has a proof; Computers & Structures or Computational Mechanics if the contribution is the comparison and the cost model.

## 10. Reuse from the package
`ch13_assets.py` (`ObstacleMembrane`, `NonlinearChainContact`, `solve_contact_dynamics`, `ContactROM`, `nonneg_dual_basis`, `svd_dual_basis`, `greedy_contact_sampling`, `nonlinearity_indicator`); `ch14_intro.py`, `ch14_basis.py`, `ch14_more.py` for the measured evidence.

## 11. Still to verify before writing
- Full texts of Balajewicz et al. (2016) and Benaceur et al. (2020): is there an error estimate for the lifted multiplier? STILL OPEN (only abstracts read).
- Google Scholar searches: `"reduced basis" contact "Lagrange multiplier" "a posteriori"`; `"non-negative" "dual basis" "reduced" contact "error estimate"`; `"hyper-reduction" contact "inf-sup" "unilateral"`.
- Whether any paper shows a passing reduced residual with a failing lifted solution.

## 12. Sources
- Banz, Schönauer, Schröder, A posteriori error estimates for variational inequalities of the second kind in an abstract framework, arXiv 2609.26001 (2026).
- Zhang, Bader, Veroy, A Duality Approach to Error Estimation for Variational Inequalities, arXiv 1410.2095.
- Sambataro, Ehrlacher, arXiv 2605.04037 (2026).
- Hild, Nicaise, ESAIM M2AN 41(5):897-923 (2007).
- Balajewicz, Amsallem, Farhat, IJNME 106(8):644-663 (2016), arXiv 1503.01000.
- Benaceur, Ern, Ehrlacher, IJNME 121(6):1170-1197 (2020).
- Fauque, Ramiere, Ryckelynck, IJNME 115(1):117-139 (2018).
- Le Berre et al., Mathematics 10(9):1495 (2022).
- Niakh, Drouet, Ehrlacher, Ern, SMAI-JCM 10:29-54 (2024), arXiv 2307.11541.
- Manvelyan-Stroot et al., arXiv 2505.18050 (2025).
- Kollepara et al., arXiv 2406.11461 (2024).
- Pawar, Kulkarni, Computers & Structures 315:107838 (2025).

## 13. Data and reproducibility
The benchmark matrices (K, B, c, f), unified Python implementations of all five dual-basis builders, and the estimator verification scripts will be released in a public repository at submission, with the random seeds and solver versions recorded.
