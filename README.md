# Sparse control of resonance permutations around multiblock exceptional points

Code and reference outputs accompanying the paper by Fatemeh Mohammadi.
The programs construct sparse perturbations of nilpotent matrices, check
characteristic polynomials and discriminants, follow eigenvalues around
parameter loops, and calculate passive scattering responses.

## Run the calculations

The suite was tested with Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8.
No TeX installation or computer algebra system is needed.
From the repository directory, run:

```bash
python -m pip install -r requirements-lock.txt
python reproduce.py
```

Precomputed reports, examples and figures are collected in
`reference_outputs.zip`. Extract that archive into a separate directory
when comparing with the reported results. The driver does not read it: it
generates fresh results, figures and logs in `reproduced/`.
No external dataset or precomputed report is required.
To choose a different output directory:

```bash
python reproduce.py --output-root /path/to/fresh-output
```

Use a fresh directory for each run. To run only the five regression tests:

```bash
python -m unittest discover -s experiments -p 'test_*.py' -v
```

`requirements-lock.txt` pins the two tested primary packages, rather than
all transitive dependencies. `requirements.txt` gives broader version ranges.

## Files

| Path | Contents |
| --- | --- |
| `reproduce.py` | Driver for the complete calculation suite. |
| `verify_all_jordan.py`, `verify_higher_rank.py` | Exact full-matrix determinant checks. |
| `independent_all_jordan_check.py`, `global_group_exact_checks.py` | Separate polynomial and permutation-group checks. |
| `certificates/` | Matrix constructors, permutation compiler, rank/support checks, discriminants and worked examples. |
| `experiments/` | Regression tests, root continuation, response calculations, reference substitutions and figure scripts. |
| `reproduced/results/` | Generated reports and numerical records, including resolvents and coherent port vectors. |
| `reproduced/figures/` | The five regenerated figures in PDF, SVG and PNG formats. |
| `CITATION.cff` | Citation metadata for the code. |
| `reference_outputs.zip` | All 43 reference output files from the complete computational archive. |

Commands for individual calculations are documented in
`certificates/README.txt`, `certificates/GENERAL_FUSION_README.txt`,
`certificates/SUPPORT_OBSTRUCTION_README.txt`, `experiments/README.md`,
`experiments/CONTINUATION_README.txt` and
`experiments/TROPICAL_REFINEMENT_README.txt`.
The paper's Appendix E lists the verification ranges and counts.

## Interpreting the results

The exact checks use integer or rational arithmetic and, where stated,
complete interpolation grids. They check finite ranges; the general
statements and support bounds are proved in the paper. The support-pattern
audit tests three specified weight assignments, rather than every
complex coefficient choice.

Root continuation uses floating arithmetic. The stored permutations are
compared at two step sizes, and their generated groups are then computed
exactly. Agreement after refinement is a consistency check, not certified
tracking. The routines implement the explicit families in the paper;
they are not a general solver for arbitrary matrix pencils.

Response differences are evaluated through cofactor and rank-one resolvent
identities to avoid subtracting nearly equal matrices. The response plots
use zero drive frequency and a fixed control phase. They do not estimate
frequency suprema or simulate measurement noise.

Support and norm refer to the stated Jordan-chain basis. Rank does not
count independently adjustable controls. The joint support minimum at
rank `g` is determined for `N >= 3` when `gcd(N,g)` is one or two;
larger gcd values have lower and upper bounds. A zero trace certificate
is inconclusive.

The cone and higher Newton calculations are in
`experiments/cone_geometry.py` and
`experiments/check_tropical_refinement.py`. Their comparisons keep the
specified coupling resources fixed within the given families. They do
not optimize over all sparse profiles.

Matrix indices in the certificates are one-based; permutation images are
zero-based. Products act rightmost first, so exported chronological loop
commands reverse the mathematical word. Sheet labels follow the chosen
branch-cycle convention and may need to be matched to an eigensolver or
experiment by a separate continuation.

The calculations concern pole permutations, not full braid images,
dynamical state transfer or quantum-generator compatibility.

## Citation

Please cite the accompanying paper when using these constructions or
calculations. Add its publication DOI to `CITATION.cff` when available.
