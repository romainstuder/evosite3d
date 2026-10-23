# Engineering

Reference notes on the tools and conventions behind evosite3d: how the environment is managed, how code quality is
enforced, and which MCP connectors are useful for life-science work.

## Contents

| Document                                   | What it covers                                                                   |
| ------------------------------------------ | -------------------------------------------------------------------------------- |
| [dev-dependencies.md](dev-dependencies.md) | The `dev` dependency group: linting, type checking, testing, pre-commit hooks    |
| [mcp-listing.md](mcp-listing.md)           | MCP connectors for biology, pharma and cheminformatics, with endpoints and tools |

## Stack at a glance

- **Python 3.12** (`requires-python = ">=3.12,<3.13"`, pinned locally in `.python-version`).
  The project only targets Python releases in **security** maintenance status: not end-of-life, not
  still receiving bug-fix releases, and not pre-release. Check the
  [Python status of versions](https://devguide.python.org/versions/) before changing the pin.
- **[uv](https://docs.astral.sh/uv/)** for packages and environments. Never use `pip`, `poetry` or `conda`, and never
  create `requirements.txt` or `setup.py`. Dependencies live in [`pyproject.toml`](../pyproject.toml) and are locked in
  [`uv.lock`](../uv.lock).
- **Runtime libraries** (see `[project.dependencies]`): Biopython, RDKit, OpenMM and MDTraj for structure and
  simulation; NumPy, pandas, SciPy, statsmodels and scikit-learn for analysis; PyTorch (pinned), Transformers and PEFT
  for models; Matplotlib and seaborn for plotting; JupyterLab for notebooks.
- **Dev tools**: Ruff, ty, mypy, pylint, pytest, pytest-cov and pre-commit. See
  [dev-dependencies.md](dev-dependencies.md).
- **[Astral](https://astral.sh/) tools**: uv, Ruff and ty all come from Astral, whose tools for Python are written in
  Rust and are extremely fast. That speed is why ty can run on every commit.

## Code quality

Run everything through `uv run`, so the project's virtual environment is always used.

| Task                    | Command                                         |
| ----------------------- | ----------------------------------------------- |
| Install the environment | `uv sync`                                       |
| Lint and auto-fix       | `uv run ruff check --fix .`                     |
| Format                  | `uv run ruff format .`                          |
| Type check              | `uv run ty check` (or `uv run mypy .`, slower)  |
| Test                    | `uv run pytest`                                 |
| Coverage                | `uv run pytest --cov --cov-report=term-missing` |
| Install commit hooks    | `uvx pre-commit install`                        |
| Run hooks on every file | `uvx pre-commit run --all-files`                |

The pre-commit hooks run file hygiene checks, Prettier on Markdown (120 columns), Ruff and ty, so most problems are
caught before a commit. Tests are collected from `scripts/tests`, `software/tests` and `tutorials/tests`.

## MCP

[mcp-listing.md](mcp-listing.md) lists remote MCP servers that Claude and other AI tools can call for literature search (PubMed, bioRxiv),
compounds and targets (ChEMBL, Open Targets), structures (RCSB PDB, beta), structure prediction (Boltz) and clinical data (ClinicalTrials.gov). It
also describes gaps, such as UniProt, which has no connector, and example multi-step workflows.

## Adding to this folder

- Keep one topic per file and link it from the table above.
- Check version numbers and commands against `pyproject.toml`, `uv.lock` and `.pre-commit-config.yaml` before
  documenting them.
- Format Markdown with Prettier at 120 columns (the pre-commit hook does this).
