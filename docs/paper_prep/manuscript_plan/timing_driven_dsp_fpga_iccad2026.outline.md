# ICCAD 2026 First-Pass Manuscript Outline

## Title Directions
1. Timing-Driven DSP Placement for Systolic Arrays on FPGAs
2. Timing-Driven Placement of DSP-Based Systolic Arrays on FPGAs
3. A Timing-Driven DSP Placement Flow for Systolic Arrays on FPGAs

## Preferred Title Direction
Use a conservative scope-first title. Right now option 1 is the cleanest because it matches the strongest writing pattern in the corpus.

## Abstract Beat List
1. Systolic-array-based FPGA accelerators are structurally regular but timing-sensitive.
2. Existing generic FPGA placement, and even earlier specialized placement, do not directly optimize the timing-critical DSP-placement problem.
3. We propose a timing-driven DSP placement flow for systolic arrays on FPGA.
4. The key mechanism is an iterative weighted placement loop that uses target-based slack or critical-edge information to update placement pressure.
5. The method preserves legality/structural constraints while driving the layout toward better timing behavior.
6. Experiments show improved timing-oriented quality over baselines on systolic-array benchmarks.

## Section Blueprint

### 1. Introduction
Goal:
- state the bottleneck quickly,
- explain why systolic-array DSP placement is a narrower and cleaner problem than generic FPGA placement,
- explain why timing must enter the placement loop,
- list contributions explicitly.

Contribution order:
1. problem framing: timing-driven DSP placement for systolic arrays on FPGA,
2. method: iterative timing-driven reweighting around a weighted placement core,
3. realization: legal structured placement flow with single-column/multicolumn handling,
4. evidence: timing-oriented improvement on completed experiments.

### 2. Background And Problem Formulation
Goal:
- give only the minimum needed background on systolic arrays, DSP columns, and the placement target,
- formalize the optimization target and the role of threshold `T` or the timing surrogate.

### 3. Proposed Method
Suggested subsection order:
3.1 Overall flow
3.2 Weighted placement core
3.3 Critical-edge or slack-driven reweighting
3.4 Termination/feasibility criterion
3.5 Multicolumn realization if it is a real part of the contribution

### 4. Experimental Setup
Goal:
- benchmark suite,
- baselines,
- implementation details,
- metrics,
- hardware/runtime notes if needed.

### 5. Results And Analysis
Suggested subsection order:
5.1 Main timing-oriented comparison
5.2 Wirelength/secondary metric tradeoff
5.3 Ablation on reweighting policy or target `T`
5.4 Case study or visualization of one representative placement

### 6. Related Work
Group by:
- systolic-array-specific FPGA placement,
- timing-driven FPGA placement,
- constrained/architecture-aware FPGA placement.

### 7. Conclusion
Goal:
- restate the problem-specific contribution,
- summarize why timing-driven DSP placement matters,
- avoid introducing new claims.

## Figure/Table Allocation
### Figures
1. Problem figure: systolic array mapped to FPGA DSP resources.
2. Method flow figure: solve -> evaluate critical edges -> update weights -> repeat.
3. Case-study figure: placement example before/after timing-driven optimization.

### Tables
1. Benchmark summary table.
2. Main comparison table with timing-oriented metrics.
3. Ablation/sensitivity table for `T`, `top_ratio`, or weighting policy.
4. Optional runtime/overhead table.

## Missing Evidence Slots
1. Final table metrics are still missing from the current outline.
2. The exact naming of the timing metric in the paper needs to match what your experiments can support honestly.
3. The role of multicolumn flow needs to be promoted or demoted depending on experiment evidence.

## Drafting Order Recommendation
1. finalize problem statement and contribution bullets,
2. finalize method overview figure and section 3 skeleton,
3. write abstract draft 0,
4. insert experiment tables when you provide them,
5. revise title and abstract after results are fixed.
