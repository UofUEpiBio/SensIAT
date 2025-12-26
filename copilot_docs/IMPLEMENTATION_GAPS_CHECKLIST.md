# Implementation Gaps Checklist

Last Updated: December 26, 2025

## High Priority

### 1. Export Generalized Fitting Function
- [ ] Add `@export` tag to `fit_SensIAT_marginal_mean_model_generalized()`
- [ ] Re-document to update NAMESPACE
- [ ] Verify function appears in package exports
- [ ] Add to package documentation index

**Status:** Not Started  
**Impact:** Critical - Users cannot access the generalized fitting functionality  
**Estimated Effort:** 5 minutes  

---

### 2. Fix Simulation Outcome Model to Use Actual Splines
- [ ] Replace polynomial approximation with actual `splines::ns()` in `generate_outcome()`
- [ ] Ensure spline parameters match those used in fitting functions
- [ ] Update tests to verify spline consistency
- [ ] Update documentation to reflect accurate implementation

**Status:** Not Started  
**Impact:** High - Simulation should match the actual model structure  
**Estimated Effort:** 1-2 hours  

---

### 3. Add Dedicated Simulation Unit Tests
- [ ] Create `tests/testthat/test-simulate_SensIAT_data.R`
- [ ] Test each link function produces correct outcome types
  - [ ] Identity → continuous
  - [ ] Log → count (non-negative integers)
  - [ ] Logit → binary (0/1)
- [ ] Test reproducibility with seed
- [ ] Test edge cases (extreme parameters)
- [ ] Test that all subjects have follow-up
- [ ] Test two-group simulation
- [ ] Test treatment effect is properly applied

**Status:** Not Started  
**Impact:** High - Ensure simulation reliability  
**Estimated Effort:** 2-3 hours  

---

## Medium Priority

### 4. Complete Loss and Links Vignette
- [ ] Add practical examples using simulation
- [ ] Demonstrate `fit_SensIAT_marginal_mean_model_generalized()` usage
- [ ] Compare results across different link functions
- [ ] Show when to use each link function
- [ ] Add interpretation guidance
- [ ] Include performance comparisons

**Status:** Not Started  
**Impact:** Medium - Users need guidance on new features  
**Estimated Effort:** 3-4 hours  

---

### 5. Implement Proper Stratified Baseline Hazard in Simulation
- [ ] Create stratified baseline hazard function
- [ ] Support vector of baseline hazards by visit stratum
- [ ] Support baseline hazard as function of (time, visit_num)
- [ ] Update documentation
- [ ] Add tests for stratified simulation

**Status:** Not Started  
**Impact:** Medium - Better alignment with actual models  
**Estimated Effort:** 2-3 hours  

---

### 6. Resolve Quasi-Likelihood Identity Link API Inconsistency
- [ ] Investigate why quasi-likelihood + identity test is skipped
- [ ] Determine if generalized and original return structures should match
- [ ] Implement fix or document intentional difference
- [ ] Enable the skipped test
- [ ] Add comparison test similar to lp_mse identity test

**Status:** Not Started  
**Impact:** Medium - API consistency important for users  
**Estimated Effort:** 1-2 hours  

---

## Low Priority

### 7. Enhance Simulation Examples and Documentation
- [ ] Add more `@examples` to simulation functions
- [ ] Show various parameter configurations
- [ ] Demonstrate treatment effect specification
- [ ] Add to package README
- [ ] Create simulation-specific vignette if needed

**Status:** Not Started  
**Impact:** Low - Helpful but not critical  
**Estimated Effort:** 2-3 hours  

---

### 8. Document Best Practices for Link Function Selection
- [ ] Create decision guide for link function choice
- [ ] Document when to use each loss function
- [ ] Add statistical background
- [ ] Provide practical examples from literature
- [ ] Add to vignette or separate documentation

**Status:** Not Started  
**Impact:** Low - Educational, not functionality  
**Estimated Effort:** 2-4 hours  

---

## Additional Observations

### Testing Coverage
- Current: 382 tests passing
- Simulation functions used in: 1 test file
- Generalized function tested in: 1 test file
- **Recommendation:** Expand test coverage for new functionality

### Documentation Status
- Simulation functions: ✅ Documented and exported
- Generalized fitting: ❌ Not exported (documented but inaccessible)
- Vignettes: ⚠️ Partial (loss_and_links.Rmd exists but incomplete)

---

## Progress Summary

**Completed:** 0/8 items  
**In Progress:** 0/8 items  
**Not Started:** 8/8 items  

**Total Estimated Effort:** 15-25 hours

---

## Notes

- This checklist should be updated as work progresses
- Mark items complete with `[x]` when finished
- Update status and add notes as needed
- Prioritize based on user needs and project timeline
