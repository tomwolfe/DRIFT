# Contributing to DRIFT

Thank you for your interest in contributing to DRIFT! We appreciate your time and effort to help improve this research framework.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Pull Request Process](#pull-request-process)

## Code of Conduct

This project and everyone participating in it is governed by the DRIFT Code of Conduct. By participating, you are expected to uphold this code.

## How Can I Contribute?

### Reporting Bugs
- Use the issue tracker to report bugs
- Describe the problem in detail
- Include steps to reproduce the issue
- Specify your environment (OS, Python version, etc.)

### Suggesting Enhancements
- Open an issue with your enhancement idea
- Explain the problem you're trying to solve
- Describe your proposed solution
- Discuss potential alternatives

### Pull Requests
- Fork the repository
- Create a branch for your feature or bug fix
- Follow the coding standards
- Add tests for new functionality
- Update documentation as needed
- Ensure all tests pass

## Development Setup

1. Fork and clone the repository:
```bash
git clone https://github.com/YOUR_USERNAME/DRIFT.git
cd DRIFT
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the package in development mode:
```bash
pip install -e .
```

4. Install development dependencies:
```bash
pip install pytest pytest-cov flake8 black
```

## Coding Standards

- Follow PEP 8 style guidelines
- Use type hints for all function parameters and return values
- Write docstrings for all public functions, classes, and modules
- Keep functions focused and reasonably sized
- Use meaningful variable and function names

## Testing

- Write unit tests for all new functionality
- Ensure all tests pass before submitting a PR
- Aim for high test coverage (>80%)

Run the test suite:
```bash
pytest tests/
```

Run tests with coverage:
```bash
pytest tests/ --cov=drift --cov-report=html
```

## Documentation

- Update docstrings when modifying code
- Add examples to function docstrings when helpful
- Update the README if adding new features
- Add/update tutorial documentation for new functionality

## Pull Request Process

1. Ensure your PR addresses a single issue or adds a single feature
2. Update the README.md with details of changes if applicable
3. Add tests for new functionality
4. Update documentation as needed
5. Ensure all tests pass
6. Submit your PR with a clear title and description

Thank you for contributing to DRIFT!