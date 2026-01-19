# DRIFT Release Process

This document describes the process for creating and publishing new releases of the DRIFT framework.

## Versioning

DRIFT follows semantic versioning (SemVer): `MAJOR.MINOR.PATCH`

- **MAJOR**: Breaking changes that affect API compatibility
- **MINOR**: New features that maintain backward compatibility
- **PATCH**: Bug fixes that maintain backward compatibility

## Pre-Release Checklist

Before creating a release, ensure:

- [ ] All tests pass (`pytest tests/`)
- [ ] Code coverage is maintained or improved
- [ ] Documentation is up-to-date
- [ ] Examples work correctly
- [ ] CHANGELOG.md is updated with all notable changes
- [ ] Dependencies are current and secure
- [ ] README.md reflects the new version's features

## Creating a Release

### 1. Update Version Number

Update the version in `pyproject.toml`:

```toml
[project]
name = "drift-sim"
version = "X.Y.Z"  # Update this
```

### 2. Update Documentation

Make sure all documentation reflects the changes in the new version:

- Update any version-specific information in README.md
- Update examples if needed
- Update API documentation if interfaces changed

### 3. Commit and Tag

Commit the version changes:

```bash
git add .
git commit -m "Bump version to vX.Y.Z"
git tag -a vX.Y.Z -m "Release version X.Y.Z"
git push origin main
git push origin vX.Y.Z
```

### 4. Automated Release Process

Once tagged, the GitHub Actions release workflow will:

1. Build the package using `python -m build`
2. Create a GitHub release with changelog
3. Upload the package to PyPI (if credentials are configured)

## Manual Release Process

If needed, you can manually build and upload the package:

```bash
# Install build tools
pip install build twine

# Build the package
python -m build

# Check the built package
twine check dist/*

# Upload to TestPyPI first (optional but recommended)
twine upload --repository testpypi dist/*

# Upload to PyPI
twine upload dist/*
```

## Post-Release Tasks

After a successful release:

- [ ] Verify the PyPI package works correctly
- [ ] Update the project website/documentation if needed
- [ ] Announce the release on relevant channels
- [ ] Update any dependent projects to use the new version

## PyPI Configuration

To publish to PyPI, you need to configure credentials:

1. Create an API token on PyPI
2. Add it as a secret in GitHub repository settings named `PYPI_API_TOKEN`

## Release Branch Strategy

- Main development happens on the `main` branch
- Releases are tagged from `main`
- Hotfix releases can be branched from release tags if needed