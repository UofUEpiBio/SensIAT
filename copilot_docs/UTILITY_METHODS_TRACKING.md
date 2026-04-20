# S3 Utility Methods Tracking

Issue: https://github.com/UofUEpiBio/SensIAT/issues/15

## Status Legend
- ❌ Not implemented
- 🚧 In progress
- ✅ Implemented

## Primary Model Classes

| Class | print | summary | coef | vcov | confint |
|-------|:-----:|:-------:|:----:|:----:|:-------:|
| `SensIAT_within_group_model` | ✅ | ✅ | ✅ | ✅ | ❌ |
| `SensIAT_fulldata_model` | ✅ | ✅ | ✅ | ✅ | ❌ |
| `SensIAT_marginal_mean_model` | ✅ | ✅ | ✅ | ✅ | ❌ |
| `SensIAT_marginal_mean_model_generalized` | ✅ | ✅ | ✅* | ✅* | ❌ |

*Inherited from `SensIAT_marginal_mean_model` base class

## Outcome Model Classes

| Class | print | summary | coef | vcov | confint |
|-------|:-----:|:-------:|:----:|:----:|:-------:|
| `SensIAT::Single-index-outcome-model` | ✅ | ✅ | ✅ | ✅* | ❌ |

*Note: `vcov()` returns NULL with a warning directing users to use `jackknife()` for variance estimation.

## Jackknife Results Classes

| Class | print | summary | coef | vcov | confint |
|-------|:-----:|:-------:|:----:|:----:|:-------:|
| `SensIAT_withingroup_jackknife_results` | ❌ | ❌ | ❌ | ❌ | ❌ |
| `SensIAT_fulldata_jackknife_results` | ❌ | ❌ | ❌ | ❌ | ❌ |

## Existing Methods (for reference)

These methods already exist for various classes:

- **predict**: `SensIAT_within_group_model`, `SensIAT_fulldata_model`, `SensIAT::Single-index-outcome-model`
- **autoplot**: `SensIAT_within_group_model`, `SensIAT_fulldata_model`, `SensIAT_withingroup_jackknife_results`, `SensIAT_fulldata_jackknife_results`
- **jackknife**: `SensIAT_within_group_model`, `SensIAT_fulldata_model`
- **prune**: `SensIAT::Single-index-outcome-model`, `SensIAT_within_group_model`
- **formula**: `SensIAT::Single-index-outcome-model`
- **model.frame**: `SensIAT::Single-index-outcome-model`
- **model.matrix**: `SensIAT::Single-index-outcome-model`

## Implementation Order

Starting with methods that will be called by higher-level methods:

1. `SensIAT::Single-index-outcome-model` - outcome model (foundation layer) ✅
2. `SensIAT_marginal_mean_model` - marginal model (builds on outcome model) ✅
3. `SensIAT_marginal_mean_model_generalized` - generalized marginal model ✅
4. `SensIAT_within_group_model` - within-group model (user-facing) ✅
5. `SensIAT_fulldata_model` - full data model (user-facing) ✅
6. `SensIAT_withingroup_jackknife_results` - jackknife results
7. `SensIAT_fulldata_jackknife_results` - jackknife results

## Notes

- Class names with `::` require backtick quoting in method definitions
- Some classes are internal (not exported) but still benefit from utility methods for debugging
- `SensIAT_marginal_outcome_model` was renamed to `SensIAT_marginal_mean_model` for consistency
- `SensIAT_marginal_mean_model_generalized` inherits from `SensIAT_marginal_mean_model`
