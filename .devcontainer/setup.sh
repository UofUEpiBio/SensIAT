#!/bin/bash
set -e

echo "Setting up SensIAT development environment..."

# System dependencies (pandoc, qpdf) are installed via devcontainer features

echo "Installing R package dependencies..."
Rscript .devcontainer/install-deps.R

echo "Installing development tools..."
R -e "install.packages(c('devtools', 'roxygen2', 'pkgdown', 'remotes'), repos='https://cloud.r-project.org')"

echo "Building package to check for compilation issues..."
cd /workspaces/pcoriRPackage || cd /workspace || cd ~
R CMD build .
R CMD check --no-manual --no-build-vignettes *.tar.gz || echo "Note: Package check had warnings/notes. Review as needed."
rm -f *.tar.gz

echo "Setup complete! You can now work on the SensIAT package."
