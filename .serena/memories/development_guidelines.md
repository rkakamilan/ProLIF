# ProLIF Development Guidelines

## Key Development Principles

### 1. MDAnalysis Migration Context
This is a specialized version of ProLIF where **MDAnalysis dependency has been removed**:
- **Focus**: Docking results analysis only (PDB/SDF files)
- **Licensing**: Moved from LGPL3 to Apache 2.0 for broader usage
- **Breaking changes**: Some functionality deprecated with clear migration paths

### 2. Error Handling Philosophy
- **Informative errors**: All removed functions provide clear error messages
- **Migration guidance**: Errors include suggestions for alternatives
- **Graceful degradation**: Handle missing dependencies without crashing

### 3. Backward Compatibility
- **API preservation**: Maintain existing method signatures where possible
- **Deprecation warnings**: Use proper versioning for breaking changes
- **Documentation**: Clear version change notes in docstrings

## Coding Standards

### Type Safety
- **Mandatory type hints**: All functions must have type annotations
- **mypy compliance**: Code must pass mypy type checking
- **Generic types**: Use proper generic types for containers and callables
- **Literal types**: Use `Literal` for string constants and options

### Documentation Requirements
- **Comprehensive docstrings**: All public methods need full documentation
- **Parameter documentation**: Use NumPy/Google style parameter descriptions
- **Examples**: Include practical usage examples in docstrings
- **Version tracking**: Use `.. versionadded::` and `.. versionchanged::` annotations

### Testing Standards
- **Comprehensive coverage**: All new functionality must have tests
- **Conditional testing**: Handle optional dependencies in tests
- **Fixture usage**: Use pytest fixtures for common test data
- **Integration tests**: Test complete workflows, not just individual functions

## Architecture Guidelines

### Interaction Classes
- **Base class inheritance**: Inherit from Distance, SingleAngle, or DoubleAngle
- **SMARTS patterns**: Use clear, well-documented SMARTS for molecular matching
- **Metadata support**: Return rich metadata, not just boolean results
- **Performance**: Consider performance for large molecular datasets

### Data Structures
- **Type-safe collections**: Use proper typing for dictionaries and lists
- **Sparse representation**: Use sparse data structures for interaction fingerprints
- **Serialization**: Support pickle/dill for saving and loading results

### Error Handling
- **Custom exceptions**: Use specific exception types from `prolif.exceptions`
- **Clear messages**: Error messages should guide users toward solutions
- **Input validation**: Validate inputs early and provide meaningful feedback

## Performance Considerations

### Parallel Processing
- **Multiprocessing support**: Design methods to work with multiprocessing
- **Memory efficiency**: Consider memory usage for large datasets
- **Progress tracking**: Use tqdm for long-running operations

### RDKit Integration
- **Molecule handling**: Efficient RDKit molecule manipulation
- **SMARTS matching**: Optimize SMARTS pattern matching
- **Coordinate handling**: Efficient 3D coordinate processing

## Special Patterns in This Codebase

### Conditional Features
```python
try:
    from MDAnalysis import AtomGroup
    _HAS_MDANALYSIS = True
except ImportError:
    _HAS_MDANALYSIS = False
    class AtomGroup:
        pass
```

### Rich Metadata
- Store complete interaction metadata, not just atom indices
- Support both boolean and count-based fingerprints
- Enable detailed analysis of interaction patterns

### Migration Support
- Provide clear error messages for deprecated functionality
- Include suggestions for alternative approaches
- Maintain API compatibility where feasible

## Git Workflow
- **Feature branches**: Use descriptive branch names (e.g., `feature/hydrophobic-improvements`)
- **Atomic commits**: Each commit should represent a single logical change
- **Clear messages**: Use descriptive commit messages
- **Pull requests**: All changes should go through pull request review