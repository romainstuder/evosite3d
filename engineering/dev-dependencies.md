# Development Dependencies

The tools in the `dev` dependency group of [`pyproject.toml`](../pyproject.toml). They are used to write, check and test
the code, and are never needed to run the scripts, tutorials or software.

```toml
[dependency-groups]
dev = [
    "mypy>=2.1.0",
    "pre-commit>=3.0.0",
    "pylint>=3.0.0",
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "ruff>=0.15.16",
    "scipy-stubs>=1.17.1",
    "ty>=0.0.46",
]
```

| Module                                              | Locked version | Purpose                                     |
| --------------------------------------------------- | -------------- | ------------------------------------------- |
| [ruff](https://docs.astral.sh/ruff/)                | 0.15.16        | Linting, import sorting, formatting         |
| [ty](https://docs.astral.sh/ty/)                    | 0.0.46         | Fast type checking, run on every commit     |
| [mypy](https://mypy-lang.org/)                      | 2.1.0          | Reference type checking                     |
| [scipy-stubs](https://github.com/scipy/scipy-stubs) | 1.17.1.5       | Type information for SciPy                  |
| [pylint](https://pylint.readthedocs.io/)            | 4.0.5          | Deeper static analysis and code smells      |
| [pytest](https://docs.pytest.org/)                  | 9.0.3          | Test runner                                 |
| [pytest-cov](https://pytest-cov.readthedocs.io/)    | 7.1.0          | Test coverage reports                       |
| [pre-commit](https://pre-commit.com/)               | 4.6.0          | Runs the checks automatically before commit |

Locked versions come from [`uv.lock`](../uv.lock); the lower bounds above are the minimum accepted.

## Installing

`uv sync` installs the `dev` group by default alongside the runtime dependencies. To install only what is needed to run
the code, leave it out:

```bash
uv sync                  # runtime + dev
uv sync --no-dev         # runtime only
uv add --dev <package>   # add a new dev tool
```

---

## Ruff

Linter and formatter in one, replacing Flake8, isort and Black. Configured in `pyproject.toml`:

```toml
[tool.ruff]
target-version = "py312"
line-length = 100

[tool.ruff.lint]
select = ["E", "F", "W", "I"]   # pycodestyle, Pyflakes, warnings, isort
ignore = ["E501"]               # line length is left to the formatter
```

```bash
uv run ruff check --fix .
uv run ruff format .
```

## ty

Astral's type checker, fast enough to run on every commit through pre-commit. Standalone scripts add `lib/` to
`sys.path` at runtime, so ty is told about it statically:

```toml
[tool.ty.environment]
extra-paths = ["lib"]
```

```bash
uv run ty check
```

## mypy

The established type checker, kept as a second opinion where ty (still in beta) does not yet cover a feature.

```bash
uv run mypy .
```

## scipy-stubs

Type stubs for SciPy, which ships without full type annotations. They let ty and mypy check code that calls
`scipy.stats`, `scipy.optimize` and similar modules. Nothing to run: the type checkers pick them up automatically.

## pylint

Slower and stricter than Ruff, it catches design issues Ruff does not cover (too many branches, unused arguments,
duplicated code). Used for occasional reviews rather than on every commit.

```bash
uv run pylint scripts software lib
```

## pytest

Test runner. Test paths are set in `pyproject.toml`:

```toml
[tool.pytest.ini_options]
testpaths = ["scripts/tests", "software/tests", "tutorials/tests"]
```

```bash
uv run pytest
```

## pytest-cov

Coverage plugin for pytest, showing which lines the tests exercise.

```bash
uv run pytest --cov --cov-report=term-missing
```

## pre-commit

Runs the checks from [`.pre-commit-config.yaml`](../.pre-commit-config.yaml) before each commit: file hygiene, Prettier
for Markdown, Ruff (lint and format) and ty. The Ruff hook pins its own version (`rev` in the config), separate from
the one locked in `uv.lock`; keep the two in step when you upgrade Ruff.

```bash
uvx pre-commit install
uvx pre-commit run --all-files
```

---

See the [engineering overview](README.md) for the runtime stack and the day-to-day commands.
