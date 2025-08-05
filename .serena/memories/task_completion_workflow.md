# Task Completion Workflow

## Before Committing Changes

### Mandatory Pre-commit Checks
Run the comprehensive check command:
```bash
poe check
```

This runs all essential checks in sequence:
1. **Style checks** (`poe style-check`)
   - Code formatting check (`ruff format --diff`)
   - Notebook formatting check (`black --diff docs/notebooks/`)
   - Linting check (`ruff check --preview --diff`)

2. **Type checking** (`poe type-check`)
   - mypy type validation across codebase

3. **Testing** (`poe test`)
   - Full pytest suite execution

4. **Documentation** (`poe docs`)
   - Sphinx documentation build

### If Any Checks Fail
1. **Fix formatting/linting issues**:
   ```bash
   poe style-fix
   ```

2. **Fix type errors**: Address mypy complaints manually

3. **Fix test failures**: Debug and resolve failing tests

4. **Fix documentation issues**: Resolve Sphinx build errors

### Alternative: Individual Checks
If you prefer to run checks individually:
```bash
poe style-check    # Check formatting and linting
poe type-check     # Check types
poe test          # Run tests
poe docs          # Build documentation
```

## Development Workflow
1. Make code changes
2. Run `poe check` to validate all aspects
3. Fix any issues that arise
4. Re-run `poe check` until all checks pass
5. Commit changes

## Important Notes
- **Never skip the checks**: The `poe check` command is designed to catch issues before they reach the repository
- **Fix issues promptly**: Don't commit with failing checks
- **Use style-fix for automation**: Let tools fix formatting/linting automatically when possible
- **Test thoroughly**: Ensure all tests pass, especially after making changes to core functionality

## Git Workflow
- All changes should be on feature branches
- Use descriptive commit messages
- Run checks before pushing to remote
- Create pull requests for code review