# ProLIF Development Commands

## Essential Development Commands

### Code Quality and Formatting
- `poe format` - Format code using ruff and black
- `poe lint` - Lint code using ruff 
- `poe style-fix` - Fix formatting and linting issues (runs format + lint)
- `poe style-check` - Check formatting and linting without fixing

### Type Checking
- `poe type-check` - Run mypy type checking across the codebase
- `poe type-check FILE_PATHS` - Run mypy on specific files

### Testing
- `poe test` - Run test suite with pytest
- `pytest tests/test_fingerprint.py::TestFingerprint::test_method_name` - Run single test

### Documentation
- `poe docs` - Build documentation with Sphinx

### Comprehensive Checks
- `poe check` - Run all checks (style, type, test, docs) - **Use this before committing**

## Individual Check Commands
- `poe format-check` - Check if code needs formatting
- `poe lint-check` - Check if code needs linting

## Package Management
- `uv sync` - Install dependencies using uv package manager
- `uv run COMMAND` - Run commands in the uv environment

## System Commands (Darwin/macOS)
- `git` - Git version control
- `ls` - List directory contents
- `find` - Find files and directories
- `grep` - Search text patterns (prefer `rg` ripgrep if available)

## Key Files to Know
- `pyproject.toml` - Project configuration and dependencies
- `prolif/_version.py` - Version information
- `CLAUDE.md` - Project-specific development guidelines
- `README.rst` - Project documentation