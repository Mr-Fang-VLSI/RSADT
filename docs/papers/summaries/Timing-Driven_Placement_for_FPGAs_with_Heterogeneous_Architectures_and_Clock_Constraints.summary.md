## Metadata
- title: Timing-Driven Placement for FPGAs with Heterogeneous Architectures and Clock Constraints
- authors: unknown
- venue/year: DATE 2021
- local_pdf: docs/papers/pdf/Timing-Driven_Placement_for_FPGAs_with_Heterogeneous_Architectures_and_Clock_Constraints.pdf

## One-Paragraph Abstract (rewritten)
The paper argues that modern FPGAs require timing-driven placement that simultaneously respects heterogeneous resources and clocking constraints. It proposes an analytical placement flow with a heterogeneity-aware delay model, clock-region refinement, and timing-driven legalization, showing that timing and wirelength can both improve relative to a commercial tool.

## Problem and Context
Modern FPGA timing closure cannot be separated from architecture heterogeneity and clock-resource legality. Placers that optimize only wirelength or ignore clock constraints can produce placements that remain infeasible or timing-poor.

## Method
- key idea: Couple timing modeling, timing-driven global placement, and timing-aware legalization under explicit clocking constraints.
- assumptions: Assumes heterogeneous FPGA fabric and structured clock-region limits; delay estimation must reflect both routing distance and block-type differences.
- model/formula: Uses a heterogeneity-aware delay model plus a memory-friendly delay lookup/interpolation scheme to drive analytical placement.

## Evidence
- result-1 (abstract): The abstract claims all timing violations are resolved while wirelength improves by 7.2% with 4.8% shorter runtime versus Vivado 2019.1 (abstract).
- result-2 (abstract/introduction): The contribution list explicitly combines timing modeling, clock refinement, and timing-based co-optimization (introduction).
- result-3 (abstract/introduction): The formulation includes legality, clock-region capacity, and non-negative slack constraints together (problem formulation).

## Limits and Threats to Validity
The method is tailored to explicit clock-resource rules and heterogeneous FPGA architectures, so portability to other fabrics or simpler accelerator-specific layouts may not be direct.

## Relevance to Current Project
- supports/contradicts/neutral: supports
- why: This is the closest technical writing sample for your current timing-driven direction. It shows how to write a placement paper whose novelty is not just timing weighting but timing-aware modeling plus constrained legalization.

## Actionable Next Step
Borrow this paper’s way of coupling timing objective, constraint story, and legalization story, then adapt it from clock constraints to DSP-column and systolic-array structure.
