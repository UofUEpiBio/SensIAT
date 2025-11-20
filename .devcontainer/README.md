# SensIAT Development Container

This devcontainer provides a complete R development environment for the SensIAT package.

## Features

- **Base Image**: `rocker/tidyverse:latest`
  - R (>= 4.4.0)
  - tidyverse packages pre-installed
  - RStudio Server available on port 8787
  
- **Tools Included**:
  - Quarto CLI (latest)
  - Pandoc
  - Git
  - C++ compiler for Rcpp packages
  
- **VS Code Extensions**:
  - R language support
  - R debugger
  - Quarto extension
  - LaTeX Workshop

## Quick Start

1. **Open in Container**:
   - Open this folder in VS Code
   - Press `F1` and select "Dev Containers: Reopen in Container"
   - Wait for the container to build (first time will take several minutes)

2. **Development Workflow**:
   ```r
   # Load the package for development
   devtools::load_all()
   
   # Run tests
   devtools::test()
   
   # Check the package
   devtools::check()
   
   # Build documentation
   devtools::document()
   
   # Build and preview pkgdown site
   pkgdown::build_site()
   ```

3. **Access RStudio Server** (optional):
   - Navigate to `http://localhost:8787`
   - Login with username: `rstudio` and password: `rstudio`

## Installed Dependencies

The container automatically installs all package dependencies from `DESCRIPTION`:

### Imports
- assertthat, dplyr, generics, ggplot2, glue, KernSmooth, MASS, MAVE, methods
- orthogonalsplinebasis, pracma, purrr, Rcpp, rlang, splines, stats, survival
- tibble, tidyr, utils

### Suggests
- dfoptim, furrr, future, inline, knitr, ManifoldOptim, metR, progress
- rmarkdown, spelling, testthat, tidyverse

### Additional Development Tools
- devtools, roxygen2, pkgdown, remotes, BB

## Rebuilding the Container

If you modify the devcontainer configuration:

1. Press `F1`
2. Select "Dev Containers: Rebuild Container"

## Manual Dependency Installation

If needed, you can reinstall dependencies:

```bash
Rscript .devcontainer/install-deps.R
```

## Notes

- The container runs as the `rstudio` user
- Your workspace is mounted at `/workspaces/pcoriRPackage` or `/workspace`
- C++17 support is enabled for Rcpp compilation
- First build may take 10-15 minutes to install all dependencies
