# ProLIF Code Style and Conventions

## Code Formatting and Linting Tools
- **ruff**: Primary linter and formatter (version 0.11.2)
- **black**: Jupyter notebook formatting
- **mypy**: Type checking

## Python Version Support
- Requires Python >= 3.10
- Supports Python 3.10, 3.11, 3.12, 3.13

## Type Hints
- **Required**: Type hints are enforced with mypy
- **Settings**: `disallow_untyped_defs = true`
- **Style**: Use proper type hints for all functions and methods
- **No `any` types**: Avoid using `Any` or `unknown` types when possible

## Docstring Style
Based on the code examples, the project uses **comprehensive docstrings** with:
- **Google/NumPy style**: Parameter and return type documentation
- **Version changes**: Documented with `.. versionchanged::` and `.. versionadded::`
- **Examples**: Include practical usage examples
- **Cross-references**: Link to related classes and methods

Example docstring format:
```python
def method(self, param: str, optional_param: bool = False) -> ReturnType:
    """Brief description of the method.
    
    Parameters
    ----------
    param : str
        Description of the parameter
    optional_param : bool
        Description with default value
        
    Returns
    -------
    ReturnType
        Description of return value
        
    Example
    -------
    ::
    
        >>> obj.method("example")
        
    .. versionadded:: 2.0.0
    """
```

## Code Structure Patterns
- **Class organization**: Methods grouped logically (init, properties, main methods, private methods)
- **Error handling**: Comprehensive error messages with migration guidance
- **Overloads**: Use `@overload` for different return types based on parameters
- **Type checking**: Use `TYPE_CHECKING` for import-only type hints

## Import Organization
- Standard library imports first
- Third-party imports
- Local imports
- Type checking imports in `if TYPE_CHECKING:` block

## Naming Conventions
- **Classes**: PascalCase (e.g., `Fingerprint`, `Molecule`)
- **Methods/functions**: snake_case (e.g., `run_from_iterable`, `get_residues_near_ligand`)
- **Private methods**: Leading underscore (e.g., `_set_interactions`, `_run_serial`)
- **Constants**: UPPER_CASE (e.g., `DEFAULT_INTERACTIONS`)

## Special Patterns in Codebase
- **Conditional imports**: Handle optional dependencies gracefully
- **Dummy classes**: For removed dependencies (MDAnalysis replacement)
- **Metadata handling**: Rich interaction metadata instead of simple atom indices
- **Progress bars**: Use tqdm for long-running operations