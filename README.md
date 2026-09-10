# Tracial transport diagnostics

Reproducible companion to **Tracial Commutator Diagnostics for Incompressible
Transport: Exact Formulas, Fourier Saturation, and Limits of Regularity
Detection**, by Gautier-Edouard Filardo.

Software version: **0.3.0**. The supplied archive is a local delivery; it does
not itself establish that a GitHub repository or release has been published.
The requested destination is `gautierfilardo-efrei/tracial-transport-diagnostics`.

## Scientific scope

For cubic Fourier observables on the unit torus, the code evaluates

\[
S_K(\Phi)=\int\log\left(1+4\sum_{0<\|k\|_\infty\le K}
\sin^2(\pi k\cdot(\Phi(x)-x))\right)\,dx.
\]

The manuscript proves the bound `S_K <= log(1 + 4*m_K)`, where
`m_K = (2*K+1)**d - 1`. Hence `S_K/K**2` tends to zero for every map. The
fixed-map saturation theorem and an explicit family of shears characterize
information the diagnostic does not retain. A derivative probe recovers
mean-square deformation.

The examples are exact toral maps and particle trajectories of an explicit,
globally smooth, decaying Beltrami solution. This software does not solve a
generic Navier-Stokes PDE, demonstrate blow-up, or establish a regularity
criterion. Earlier purported Hou-Luo values are retained only in a clearly
identified audit CSV, with the five bound violations marked.

## Reproduce from a clean extraction

Use Python 3.12 and run the following commands in this directory:

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 -m unittest discover -s tests -v
python3 reproduce.py --output-dir build/reproduction
python3 verify_reproduction.py --generated build/reproduction
```

On Windows, activate the environment with `.venv\Scripts\activate` and use
`python` if `python3` is unavailable. Virtual-environment creation is optional
when an appropriate environment already exists.

The reference run uses Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8.
The output directory keeps new runs separate from the packaged reference.
To intentionally regenerate the reference files in place, run `python3
reproduce.py` without `--output-dir`.

The verifier compares all CSV tables and raw NPY arrays, using relative
tolerance `1e-9` and absolute tolerance `1e-10`, and checks the rounded TeX table.
Figure files are checked for existence; PDF bytes are not compared because
rendering metadata can vary. Analytical tests provide independent checks in
addition to this numerical regression comparison.

## Files

| Path | Purpose |
| --- | --- |
| `diagnostics.py` | Direct Fourier sum, Dirichlet formula, analytical bound and ABC particle flow |
| `reproduce.py` | Regenerate every table, sample dataset and figure |
| `verify_reproduction.py` | Compare a separate run against the reference data |
| `tests/` | Nine analytical/numerical checks and three groups of input checks |
| `results/` | Six CSV tables, five raw NPY datasets, environment and verification records |
| `figures/` | Three figures used in the article |
| `benchmark_table.tex` | Table generated from the computed benchmark values |
| `CITATION.cff` | Software citation metadata, without an invented DOI |
| `docs/` | Mathematical conventions, provenance, references and GitHub setup |

## GitHub checks

The included workflow runs the same tests, regenerates the results and compares
them with the reference on pushes and pull requests. It requires no repository
write permission. Local test results do not imply that this workflow has
already run on GitHub.

The workflow structure follows the official
[GitHub Python workflow documentation](https://docs.github.com/en/actions/tutorials/build-and-test-code/python).
See `docs/CREER_GITHUB.md` for the prepared account-specific creation procedure.

## Citation and status

Use `CITATION.cff` for the software metadata. The revised article is supplied
separately in the Overleaf archive. Its older preprint has different conclusions;
do not present it as an identical version. Add a repository URL, commit hash and
archive DOI to the manuscript only after the corresponding records exist.

No distribution license has been selected in this preparation. The supplied
creation procedure starts with a private repository. License selection and
public release are subsequent author decisions. The calculations and assistance
used in this revision are described in `docs/PROVENANCE.md`.
