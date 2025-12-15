# .Rprofile for SensIAT package development
# This file loads common development packages automatically

# Load development packages if available
if (interactive()) {
    # Load devtools for package development workflow
    if (requireNamespace("devtools", quietly = TRUE)) {
        library(devtools)
        cat("✓ devtools loaded\n")
    } else {
        cat("✗ devtools not available - install with: install.packages('devtools')\n")
    }

    # Load usethis for package setup and maintenance
    if (requireNamespace("usethis", quietly = TRUE)) {
        library(usethis)
        cat("✓ usethis loaded\n")
    } else {
        cat("✗ usethis not available - install with: install.packages('usethis')\n")
    }

    # Optional: Set up convenient aliases
    if (exists("load_all")) {
        la <- load_all
    }
    if (exists("document")) {
        doc <- document
    }
    if (exists("test")) {
        t <- test
    }

    cat("Package development environment ready!\n")
    cat("Shortcuts: la() = load_all(), doc() = document(), t() = test()\n")
}
