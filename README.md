# Four-way decomposition of a longitudinal outcome slope

How much of an exposure contrast in **annual outcome change** can be attributed to mediation, interaction, or both? This worked example applies Tyler J. VanderWeele's four-way decomposition to a longitudinal outcome using a Bayesian two-level model in Mplus. It estimates the mediator regression and the outcome trajectory in one analysis, then decomposes the difference in mean slopes at each posterior draw. The example uses a continuous outcome, one continuous mediator measured after exposure, and one pre-exposure covariate. All 320 participants and 1,458 visits are simulated.

## Run the example

Download the repository and keep the input and data files together. Open `four_way.inp` in Mplus and choose **Run**, or run the following command from this directory with the Mplus executable on your path. The example was tested with Mplus 8.7 for Mac with multilevel modeling support.

```sh
mplus four_way.inp
```

In `four_way.out`, find **MODEL RESULTS → New/Additional Parameters**. The five rows from `CDE` through `TE` give posterior medians and equal-tailed 95% credible intervals. The run also produces `four_way.gh5` for Mplus plots and `posterior.dat` for draw-level calculations.

The supplied `synthetic.dat` is ready to use. To regenerate it or print the effect and magnitude-share tables, use base R (tested with R 4.5.2). The second command reads the outputs of a completed Mplus run.

```sh
Rscript --vanilla simulate.R
Rscript --vanilla summarize.R
```

## Expected results

The example compares exposure `a1 = 1` with `a0 = −0.5`, fixes the mediator at `mstar = −0.5` for the controlled direct effect, and evaluates the contrast at covariate value `cref = 0`. The simulation coefficients in `simulate.R` give a known population value for each effect. The last two columns come from the seeded Mplus example.

| Component | Interpretation | Simulation truth | Posterior median | 95% credible interval |
| --- | --- | ---: | ---: | ---: |
| CDE | Exposure contrast in slope with the mediator fixed at `mstar` | 0.3000 | 0.247 | 0.180 to 0.311 |
| INTref | Interaction attributable to the mediator level already present under `a0`, relative to `mstar` | 0.1650 | 0.178 | 0.147 to 0.211 |
| INTmed | Interaction attributable to the exposure-induced change in the mediator | 0.2250 | 0.229 | 0.188 to 0.273 |
| PIE | Slope difference from the exposure-induced mediator change, evaluated at `a0` | 0.1125 | 0.133 | 0.091 to 0.179 |
| TE | Total exposure contrast in mean slope | 0.8025 | 0.787 | 0.715 to 0.858 |

Here, a total effect of 0.787 means that the mean outcome rises 0.787 units more **per year** under the target exposure than under the reference exposure, allowing the mediator to take its natural value under each exposure. The four signed components add to TE within every posterior draw. The sum of their posterior medians can differ from the posterior median of TE.

The optional R summary also reports each component's absolute magnitude divided by the sum of the four absolute magnitudes, calculated within each draw. These percentages describe relative size. Opposing components can cancel in the signed total, so retain the signed effects when reporting direction.

## Model and estimand

For person $i$ at visit $j$, the continuous outcome follows $Y_{ij}=I_i+S_i t_{ij}+e_{ij}$. The random intercept $I_i$ and the residual of the random slope $S_i$ may covary. The example uses a common mean outcome at time zero across exposure and mediator levels. The labels below match `four_way.inp`.

```text
E[M | A, C]    = b0 + b1 A + bc C
E[S | A, M, C] = t0 + t1 A + t2 M + t3 A M + tc C
```

Let `da = a1 − a0`, `m0 = b0 + b1*a0 + bc*cref`, and `dm = b1*da`. Here `m0` is the expected mediator under the reference exposure and `dm` is its expected change under the contrast. Linearity in the mediator allows these expectations to enter the slope contrast directly.

```text
CDE    = (t1 + t3*mstar) * da
INTref = t3 * (m0 - mstar) * da
INTmed = t3 * da * dm
PIE    = (t2 + t3*a0) * dm
TE     = CDE + INTref + INTmed + PIE
```

As a check, TE also equals the predicted slope at `(a1, m0 + dm, cref)` minus the predicted slope at `(a0, m0, cref)`. Changing only `mstar` redistributes the effect between CDE and INTref while preserving TE. The decomposition follows VanderWeele's additive-scale definitions, applied here to the longitudinal slope. [VanderWeele, Epidemiology (2014)](https://doi.org/10.1097/EDE.0000000000000121) · [Mplus User's Guide, example 9.2](https://www.statmodel.com/usersguide/chapter9.shtml)

## Adapt the example to your data

1. Prepare a headerless, numeric, space-delimited file with one row per visit. Sort rows by participant and time. The columns in order are `id` (numeric participant identifier), `time` (elapsed years), `y` (observed outcome), `a` (exposure), `m` (mediator), and `c` (pre-exposure covariate). Repeat the same `a`, `m`, and `c` values on every row for a person. Use actual elapsed times for irregular visits. Align time zero and the follow-up window with the mediator measurement.
2. Point `DATA: FILE` to your file. The example assumes complete input values. For incomplete data, explicitly declare your missing-value code and check the retained participant and visit counts in Mplus. A missing exposure or mediator also makes the observed interaction product missing.
3. Set `a0`, `a1`, `mstar`, and `cref` in `MODEL CONSTRAINT` on the scales used in the file. For a binary exposure coded 0/1, set `a0 = 0` and `a1 = 1`. For a continuous exposure, either direction is supported. Choose a meaningful reference mediator value with empirical support.
4. For additional baseline confounders, add columns to `NAMES`, `USEVARIABLES`, and `BETWEEN`. Include each in both the slope and mediator regressions, label each coefficient in the mediator regression, and include its coefficient times its target value in `m0`. Keep the generated variable `am` last in `USEVARIABLES`. Centering covariates on participant-level means makes zero the mean profile. With the additive covariate terms shown, evaluating at that profile also averages the contrast over those participants.
5. If the expected outcome at time zero depends on the predictors, add `y ON a m am c;` in the between-person model and assess its fit. For nonlinear trajectories, time-varying exposures or mediators, or binary outcomes, the model and target contrast need different specifications.

## Inference and interpretation

The example uses Mplus 8.7 default priors, two chains, and 4,000 fixed iterations per chain. The R summary uses the last 2,000 draws from each chain. Inspect `TECH8` in the output and inspect trace plots for the slope coefficients and decomposition terms. Rerun with a longer chain length and compare estimates, intervals, and the potential scale reduction statistic. For other data, check the implied linear trajectories, residual variation, covariate overlap, and sensitivity to priors and the observation process.

A causal interpretation requires a well-defined exposure contrast, consistency, no interference, and support for the required exposure and mediator levels. The pre-exposure covariates must suffice to control exposure–outcome, exposure–mediator, and mediator–outcome confounding. The natural effects also require the cross-world identification assumption, including the absence of exposure-induced mediator–outcome confounding under the usual causal model. The timing, confounder set, and handling of missing visits require study-specific justification. [Identification assumptions and derivation](https://doi.org/10.1097/EDE.0000000000000121)
