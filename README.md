# resolution-gap-arrow-of-time

[![verify](https://github.com/neginyan/resolution-gap-arrow-of-time/actions/workflows/verify.yml/badge.svg)](https://github.com/neginyan/resolution-gap-arrow-of-time/actions/workflows/verify.yml)

Code, results and verification for the paper

> T. Namba, *A resolution gap for the coarse-grained arrow of time in nonlinear Hamiltonian flows*, preprint (2026).

The paper continues
T. Namba, *Sharp criterion for the coarse-grained arrow of time in Hamiltonian dynamics*, preprint v1.1.1, Zenodo (2026),
[doi:10.5281/zenodo.23163727](https://doi.org/10.5281/zenodo.23163727).

## Question

A Hamiltonian flow is observed through a Gaussian window of width σ that coarsens at rate r = σ̇/σ.
The coarse-grained entropy S_C satisfies dS_C/dt = σ² tr F (r − R), with F the Fisher matrix of the
coarse-grained density and R a *critical rate* that depends on the state. The earlier paper proved
R ≤ κ (κ = largest stretching rate of the flow) for linear flows and conjectured it for nonlinear ones.
This study asks how close a state can bring R to κ when the flow is nonlinear.

## Results

- **Exact decomposition.** For every state, R = L + C + D: a law term L ≤ κ, a covariance term C between
  the local focus and the local strain, and a Stein defect D. A violation of R ≤ κ needs C or D.
- **Gaussian states.** For every Gaussian state, in any dimension, R = tr(F Sym J̄)/tr F ≤ κ
  (statistical linearization), so the conjecture holds for all Gaussian states.
- **Best Gaussian state.** For one degree of freedom, R = λφ (strain felt × anisotropy of focus).
  For the pendulum the supremum reduces to one variable,
  R\*(σ) = sup_x ½(1 + e^{−x²/4}) x/√(x² + 4σ²), with 1 − R\* = σ − 7/8 σ² + 65/96 σ³ − …,
  attained by a slightly tilted segment of length ≈ 2√σ. Quadratic strain has a closed form, and
  near a non-degenerate maximum κ − R\* = 2√(2βκ) σ + O(σ²): the *resolution gap*.
- **Mixtures beat the best Gaussian state.** Fans of thin segments crossing at the hyperbolic point
  exceed R\* by up to 5 % of the gap, but the gap stays proportional to σ (smallest observed 0.946σ).
  Gaussian states are therefore not extremal.
- **State-wise lower bound (pendulum).** 1 − R ≥ Θ′σ/(1 + Θ′σ/2) − ε̃ − D for every state, from a
  weighted Cramér–Rao inequality along the contracting direction; at least 0.898σ on all stored states
  with σ ≤ 10⁻². A uniform bound would follow from control of the focus slope (third derivatives of
  ln ρ_C), which remains open.

## Repository layout

| path | contents |
|---|---|
| `paper/` | LaTeX source (`main.tex`, `sections/`, `figures/`) and the compiled PDF |
| `study/` | numerical study, stages A–K: scripts, raw results (`results/`), logs, figures, stage summaries |
| `study/verification/` | verification suite: one script per stage and `run_all_stages.py` |
| `independent_checks/` | checks written separately from the study code (entropy by finite differences, best Gaussian state) |

`study/README.md` describes every stage and every script; it also lists the commands that regenerate the raw results.

### Paper figures

| figure | file in `study/figures/` | made by |
|---|---|---|
| Fig. 1 (mechanism) | `stageC_spatial.png` | `study/make_figures_C.py` |
| Fig. 2 (best Gaussian state) | `stageH_rstar.png` | `study/h_analyze.py` |
| Fig. 3 (mixtures beat it) | `stageH_violation.png` | `study/h_analyze.py` |
| Fig. 4 (four flows) | `stageI_flows.png` | `study/i_analyze.py` |
| Fig. 5 (state-wise bound) | `stageJ_bound.png` | `study/j_analyze.py` |

## Verify

```bash
pip install -r requirements.txt
cd study
python verification/run_all_stages.py     # 281 checks
```

Each stage script re-derives every number it reports from the stored raw results and checks the identities and
inequalities used in the paper (for example the decomposition to 10⁻¹¹, the second-order Tweedie formula, and each
step of the state-wise bound on every state). Each check reports the fraction of its tolerance it uses; none uses
more than 20 %. The suite runs on every push (GitHub Actions badge above).

The independent checks are described in [`independent_checks/README.md`](independent_checks/README.md).

### Environment

The study and the verification used Python 3.13.16 with numpy 2.5.3, scipy 1.18.1 and matplotlib 3.11.2
(pinned in `requirements.txt`). No GPU and no machine-learning framework are used; all integrals are done by quadrature.

## Use of AI assistants

This work was carried out with the help of AI assistants: Google Gemini (Gemini 3.8 Flash) and Anthropic Claude.
Their roles are described in the acknowledgments of the paper. Every analytical statement, number and figure is
checked by the scripts in this repository.

## License

MIT (see `LICENSE`).
