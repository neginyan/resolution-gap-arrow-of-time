# Independent checks

These scripts were written separately from the study code in `../study/`. Methods (iii) and (iv) of
Appendix C of the paper.

## Entropy by finite differences (method iii)

No identity of Sec. II of the paper is used (no Stein lemma, no Tweedie formula, no posterior moments).
The fine-grained density is transported exactly by the flow over ±Δt (fourth-order Runge–Kutta for the
characteristics), blurred with the Gaussian window, and Σ = dS_C/dt is taken as the central difference of the
entropy; then R = −Σ/(σ² tr F).

| script | what it checks | usage |
|---|---|---|
| `indep.py` | core routine (FFT blur on a grid); run directly, it evaluates R for two-component pendulum mixtures of different widths (Appendix C) | `python indep.py` |
| `top.py` | the largest-R state stored in a results file | `python top.py ../study/results/<file>.json` |
| `topE.py` | the same for nested result files (used for the state with R = 0.98545), at grid size N | `python topE.py ../study/results/<file>.json N` |
| `topH.py` | a mixture that beats the best Gaussian state (the violations of R ≤ R\* at σ = 0.1, 0.3), chosen by K and σ | `python topH.py ../study/results/h2_adversarial.json K sigma N pad` |
| `twist_direct.py` | the best K = 8 twist-flow state; blur by summing kernels over transported quadrature nodes (for very thin components) | `python twist_direct.py ../study sigma i4_twist_chain501.json n_per_sigma dt` |

`top*.py` import `indep.py` from the current directory, so run them from this folder.
`twist_direct.py` reads the parametrisation of the stored state from `../study` (`i3_flows.unpack`); the evaluation itself is independent.

## Best Gaussian state (method iv)

| script | what it checks | usage |
|---|---|---|
| `rstar.py` | R\*(σ) by direct optimisation of R = λφ over all Gaussian states (four parameters) | `python rstar.py` |
| `rstar_check.py` | R\*(σ) over zero-width segments against the series of Theorem 5 and the closed form of Theorem 6; the leading coefficient 2√(2βκ) for the pendulum, the double well and the quartic oscillator | `python rstar_check.py` |
