# First-Pass Method Story

## Scope
This note is a first-pass story artifact derived from the current `src/` codebase, not a final paper claim set. It is intended to stabilize the manuscript narrative before the final result tables are inserted.

## Main Code Entry Points
- `src/main_td.cpp`
- `src/main_singlecol.cpp`
- `src/momentum_weighter.h`
- `src/weighter_policy.h`
- `src/lightOCShortest.cpp`
- `src/multicol_splitter.cpp`

## One-Paragraph Problem Statement
The project targets timing-driven placement for systolic-array-style DSP-heavy designs on FPGA. Generic placement objectives such as equal-weight HPWL are too coarse because the design class has highly structured local communication, strong DSP-column regularity, and a timing target that is more naturally expressed as a per-edge or path-length budget. The code suggests that the practical goal is to generate legal structured placements while iteratively emphasizing timing-critical edges whose current geometric length violates a target threshold `T`.

## One-Paragraph Method Statement
The current implementation appears to solve placement through an iterative weighted optimization loop. In each round, a base placement solver (`lightOCShortest`) produces a layout under current edge weights. A timing-oriented reweighting stage then computes edge lengths from the layout, derives slack-like quantities `s_e = T - L_e`, identifies the most critical edges, and increases their weights through a momentum-style update in log-weight space. The next round solves placement again under the updated weights. In effect, the method turns timing pressure into adaptive edge weights and repeatedly pushes the placement toward satisfying a target maximum edge-length budget while preserving legality.

## Pipeline In Ordered Steps
1. Initialize horizontal and vertical edge weights, optionally from external matrices.
2. Export current weights to integer form for the placement solver.
3. Solve a placement instance with `lightOCShortest` under the current weighted objective.
4. Compute equal-weight and weighted HPWL for monitoring.
5. Derive slack-like quantities from edge lengths with respect to target `T`.
6. Select the most critical edges, currently by top-ratio criticality.
7. Update edge weights through a momentum-style rule in log space.
8. Repeat until the layout becomes feasible under the target or the round budget ends.
9. Emit a single-column placement, with a multicolumn replication/splitting flow available separately.

## Mechanism Explanation

### A. Timing Signal Injection
The code injects timing pressure through the target threshold `T` and edge slack `s_e = T - L_e`. This means timing is not represented only as a post hoc metric; it directly changes weights for the next placement round.

### B. Adaptive Critical-Edge Emphasis
`momentum_weighter.h` maintains `logw` and `dlog` terms for horizontal and vertical edges. The update rule amplifies selected critical edges through a momentum-like recurrence, which makes the solver remember previous pressure instead of fully resetting every round.

### C. Weighted Placement Core
`lightOCShortest` receives weighted horizontal and vertical cost matrices. This suggests the actual placement engine is solving a weighted geometric optimization problem, while the timing-driven part lives in the outer reweighting loop.

### D. Feasibility-Style Termination
The loop stops early when `WNS >= 0`, which operationally means all tracked edges satisfy the target `T`. This is a clean story hook for the paper: the algorithm is not merely reducing average wirelength; it is trying to drive the placement into a target-feasible region.

## Core Novelty Split
- `core_publishable_mechanism`: iterative timing-driven reweighting under a target edge-length budget, with momentum-style updates and repeated weighted placement solves.
- `supporting_engineering`: file parsing, integer weight export, output formatting, and multicolumn placement replication/splitting.
- `must_not_oversell`: do not claim full signoff timing optimization from the current code reading alone. The present evidence supports a timing-inspired or slack-driven placement surrogate story, but not yet a full STA-coupled timing-closure claim.

## Figure Candidates
1. Problem figure: systolic-array communication pattern mapped onto FPGA DSP columns.
2. Flow figure: initialize weights -> weighted placement solve -> slack extraction -> critical-edge update -> repeat.
3. Mechanism figure: how `T - L_e` drives reweighting of critical edges.
4. Example figure: single-column base placement and multicolumn expansion/splitting.

## Open Ambiguities
1. The current code evidence is strongest at edge-level length/slack control; whether the final paper should say path-based timing or edge-budget timing still needs confirmation from experiments and any external notes.
2. The relation between single-column optimization and multicolumn realization needs a cleaner narrative link.
3. The exact novelty over the predecessor papers still needs to be aligned against your intended final claim set.
