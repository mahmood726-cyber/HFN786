# Contributing to CNMA

Thank you for considering contributing to CNMA! This document provides guidelines for contributing to the project.

## Code of Conduct

By participating in this project, you agree to maintain a respectful and collaborative environment.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues. When creating a bug report, include:

* **Clear title** - Descriptive summary of the issue
* **Steps to reproduce** - Minimal reproducible example
* **Expected behavior** - What you expected to happen
* **Actual behavior** - What actually happened
* **Environment** - R version, OS, package versions
* **Code example** - Minimal code to reproduce the issue

Example:

```r
# Bug: validate_cnma_input fails on valid data

library(cnma)

data <- data.frame(
  studlab = c("S1", "S2"),
  treat1 = c("A", "A"),
  treat2 = c("B", "C"),
  TE = c(0.5, 0.3),
  seTE = c(0.1, 0.2)
)

validate_cnma_input(data)  # Error: unexpected behavior

# R version 4.3.0
# cnma version 1.0.0
# OS: Ubuntu 22.04
```

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. Include:

* **Clear use case** - Why is this enhancement useful?
* **Proposed solution** - How should it work?
* **Alternatives** - Other approaches considered
* **Examples** - Code examples showing desired usage

### Pull Requests

1. **Fork the repository** and create a branch from `main`
2. **Make your changes** following the coding guidelines below
3. **Add tests** for new functionality
4. **Update documentation** (roxygen comments, README, NEWS.md)
5. **Run checks** (`R CMD check`, tests pass)
6. **Submit pull request** with clear description

## Development Setup

### Prerequisites

* R >= 4.0.0
* RStudio (recommended)
* Git
* Required packages:
  ```r
  install.packages(c(
    "devtools", "testthat", "roxygen2",
    "netmeta", "ggplot2", "dplyr", "tidyr", "purrr", "tibble"
  ))
  ```

### Getting Started

```bash
# Clone your fork
git clone https://github.com/your-username/HFN786.git
cd HFN786

# Create a branch for your changes
git checkout -b feature/my-new-feature
```

In R:

```r
# Load development version
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()

# Build documentation
devtools::document()
```

## Coding Guidelines

### Style Guide

Follow the [tidyverse style guide](https://style.tidyverse.org/):

* Use `snake_case` for function and variable names
* Use `<-` for assignment (not `=`)
* Limit lines to 80 characters
* Use 2 spaces for indentation (no tabs)
* Add spaces around operators (`x + y`, not `x+y`)
* Use explicit returns when returning values

Example:

```r
# Good
calculate_effect_size <- function(mean_diff, sd_pooled) {
  effect_size <- mean_diff / sd_pooled
  return(effect_size)
}

# Bad
calculateEffectSize=function(meanDiff,sdPooled){
return(meanDiff/sdPooled)
}
```

### Function Documentation

All exported functions must have roxygen2 documentation:

```r
#' Brief description (one line)
#'
#' Detailed description explaining the function's purpose,
#' behavior, and any important details.
#'
#' @param param1 Description of first parameter
#' @param param2 Description of second parameter
#' @return Description of return value
#' @export
#' @examples
#' # Example usage
#' result <- my_function(param1 = "value", param2 = 42)
my_function <- function(param1, param2) {
  # Implementation
}
```

### Testing

* Write tests for all new functions
* Aim for >80% code coverage
* Use `testthat` framework
* Include edge cases and error conditions

Example test:

```r
test_that("calculate_effect_size works correctly", {
  # Normal case
  expect_equal(calculate_effect_size(1, 2), 0.5)

  # Edge cases
  expect_equal(calculate_effect_size(0, 1), 0)
  expect_error(calculate_effect_size(1, 0), "division by zero")

  # Multiple values
  expect_equal(
    calculate_effect_size(c(1, 2), c(2, 4)),
    c(0.5, 0.5)
  )
})
```

### Error Handling

* Use informative error messages
* Use `.stop_hint()` helper to provide hints
* Validate inputs early
* Handle edge cases gracefully

```r
# Good
validate_input <- function(x) {
  if (!is.numeric(x)) {
    .stop_hint(
      "x must be numeric",
      "Convert your input using as.numeric()"
    )
  }
  if (length(x) == 0) {
    .stop_hint(
      "x must have at least one element",
      "Check that your data is not empty"
    )
  }
}
```

## Project Structure

```
cnma/
├── R/                          # R source files
│   ├── aaa-package.R          # Package documentation
│   ├── utilities.R            # Utility functions
│   ├── config.R               # Configuration
│   ├── validation.R           # Input validation
│   ├── data-simulation.R      # Data generation
│   ├── analysis.R             # Main analysis
│   ├── methods.R              # S3 methods
│   └── parallel.R             # Parallel processing
├── tests/                     # Test files
│   ├── testthat.R
│   └── testthat/
│       ├── test-utilities.R
│       ├── test-config.R
│       ├── test-validation.R
│       ├── test-data-simulation.R
│       ├── test-analysis.R
│       ├── test-methods.R
│       └── test-parallel.R
├── man/                       # Documentation (auto-generated)
├── DESCRIPTION               # Package metadata
├── NAMESPACE                 # Export declarations
├── README.md                 # Main documentation
├── NEWS.md                   # Changelog
└── LICENSE                   # License file
```

## Adding New Features

### 1. Design Phase

* Open an issue to discuss the feature
* Get feedback from maintainers
* Write a design document for complex features

### 2. Implementation

```r
# Create new file if needed (e.g., R/new-feature.R)

#' New Feature Function
#'
#' Detailed description of what this does.
#'
#' @param data Input data
#' @param ... Additional parameters
#' @return Results object
#' @export
#' @examples
#' \donttest{
#' result <- new_feature(data)
#' }
new_feature <- function(data, ...) {
  # Validate inputs
  stopifnot("data must be data.frame" = is.data.frame(data))

  # Implementation
  result <- process_data(data, ...)

  # Return
  return(result)
}
```

### 3. Testing

```r
# tests/testthat/test-new-feature.R

test_that("new_feature works with valid input", {
  data <- create_test_data()
  result <- new_feature(data)

  expect_s3_class(result, "expected_class")
  expect_true("required_field" %in% names(result))
})

test_that("new_feature validates input", {
  expect_error(new_feature("not a data frame"), "data must be")
})
```

### 4. Documentation

* Update roxygen comments
* Add examples to README if relevant
* Update NEWS.md
* Consider adding a vignette for complex features

### 5. Submit PR

* Run `devtools::check()` - should have 0 errors, 0 warnings
* Run `devtools::test()` - all tests should pass
* Update NEWS.md with your changes
* Submit pull request with clear description

## Release Process

(For maintainers)

1. Update version in DESCRIPTION
2. Update NEWS.md with release notes
3. Run `devtools::check()` locally
4. Run `devtools::test()` - all tests pass
5. Build package: `devtools::build()`
6. Submit to CRAN (if applicable)
7. Create GitHub release with tag
8. Update main branch

## Questions?

* Open an issue for questions
* Check existing documentation
* Review closed issues/PRs for similar questions

## Recognition

Contributors will be acknowledged in:
* README.md contributors section
* Package DESCRIPTION file (for significant contributions)
* Release notes

Thank you for contributing to CNMA! 🎉
