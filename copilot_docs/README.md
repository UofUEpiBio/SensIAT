# Copilot Documentation

This directory contains documentation generated during AI-assisted development sessions with GitHub Copilot.

## Contents

### Bug Analysis & Debugging

- **[VECTORIZED_INTEGRATION_BUG.md](VECTORIZED_INTEGRATION_BUG.md)**: Comprehensive analysis of the critical bug in the C++ adaptive Simpson integration implementation. Documents root cause, evidence, and required fixes.

- **[VECTORIZED_INTEGRATION_DEBUG_PLAN.md](VECTORIZED_INTEGRATION_DEBUG_PLAN.md)**: Systematic debugging plan for the vectorized integration feature. Outlines hypothesis testing, diagnostic tests, and verification strategies.

### Testing Strategies

- **[JACKKNIFE_TESTING_STRATEGY.md](JACKKNIFE_TESTING_STRATEGY.md)**: Multi-level testing strategy for jackknife resampling functionality. Includes fast tests for CI, medium-scale tests for validation, and comprehensive tests for major changes.

## Purpose

These documents serve as:

1. **Development History**: Track the debugging process and decision-making
2. **Knowledge Base**: Document complex bugs and their solutions
3. **Testing Guidance**: Provide strategies for testing computationally expensive features
4. **Onboarding**: Help new contributors understand implementation challenges

## Status

**Current State (2025-12-12)**: Vectorized integration is work-in-progress with a known bug in the C++ adaptive Simpson implementation. The integrand calculation is correct, but the integration algorithm needs to be rewritten to use per-segment convergence instead of global convergence flags.

## Next Steps

See [VECTORIZED_INTEGRATION_BUG.md](VECTORIZED_INTEGRATION_BUG.md) for the recommended fix implementation.
