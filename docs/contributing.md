# Contributing

Contributions are welcome. This document describes how to set up a development environment and the conventions the
project follows.

---

## Development setup

```bash
git clone https://github.com/andreassag/clystere.git
cd clystere

# Create a Python virtual environment (Python 3.11+)
python -m venv .venv
source .venv/bin/activate

# Install development tools
pip install pre-commit pytest ruff
npm install

# Activate the pre-commit git hooks
git config core.hooksPath .githooks
```

---

## Code style

### Python (`bin/`)

Python scripts are linted and formatted with [Ruff](https://docs.astral.sh/ruff/). The target interpreter is Python
3.11.

```bash
ruff check .
ruff format --check .
```

### Nextflow (`.nf`, `.config`) & Documentation

Nextflow and config files follow EditorConfig rules and Prettier standards:

```bash
npx editorconfig-checker
npx prettier --check .
```

---

## Pre-commit hooks

Pre-commit git hooks run automatically before each `git commit` to guarantee clean code and formatting:

| Check                   | Tool                   | Purpose                          |
| :---------------------- | :--------------------- | :------------------------------- |
| Python Linting          | `ruff check`           | Code style and import order      |
| Python Formatting       | `ruff format`          | Code layout consistency          |
| Line Length & Standards | `editorconfig-checker` | Indentation and line endings     |
| Document Formatting     | `prettier`             | Markdown, JSON, and YAML styling |

Run all checks manually:

```bash
npm run check:all
```

---

## Testing

### Python unit tests

```bash
pytest tests/python/ -v
```

Tests live in `tests/python/` and exercise the pure Python functions in `bin/`.

### nf-test (Nextflow pipeline & module tests)

```bash
# Run all nf-test suites (integration + local modules)
nf-test test tests/ --verbose

# Run specific local module tests
nf-test test tests/modules/local/summary --verbose
nf-test test tests/modules/local/combgc --verbose
nf-test test tests/modules/local/bigscape --verbose
```

### Full pipeline test

```bash
nextflow run . \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile test,docker
```

---

## Branching and pull requests

- `main` is the stable release branch.
- Feature branches: `feature/<description>`
- Bug/hot fixes: `fix/<description>`
- Non-code tasks: `chore/<description>`
- Preparing a release: `release/<description>`
- Open a pull request against `main`; CI must be green before merge.

---

## Releases

Releases follow [Semantic Versioning](https://semver.org/). To cut a release:

1. Update the version in `nextflow.config` (`manifest.version`).
2. Commit with message `release: bump version to X.Y.Z`.
3. Push a `X.Y.Z` tag — the `release.yml` workflow creates the GitHub Release automatically.

```bash
git tag 1.1.0
git push origin 1.1.0
```
