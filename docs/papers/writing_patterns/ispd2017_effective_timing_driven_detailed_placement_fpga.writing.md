- paper_id: ispd2017_effective_timing_driven_detailed_placement_fpga
- title: An Effective Timing-Driven Detailed Placement Algorithm for FPGAs
- venue_year: ISPD 2017
- writing_role: method_narration_anchor

## 1. Title Pattern
Conservative and effective title. It emphasizes task and platform, not branding. This is still a safe style for ICCAD/ISPD placement papers.

## 2. Abstract Move Order
1. timing bottleneck
2. limits of previous optimization style
3. new formulation
4. why it explores more solution space
5. timing result

## 3. Introduction Framing
Strong motivation section. It explicitly contrasts net-based and path-based approaches and then bullet-lists contributions in reviewer-friendly form.

## 4. Method Narration
The method is narrated through the optimization object itself. That is useful when your paper’s novelty is a specific algorithmic formulation.

## 5. Experiment Narration
Results are likely centered on Fmax improvement with minimal side effects on wirelength and congestion. Good pattern for primary/secondary metric ordering.

## 6. Wording Habits
Uses formulate, optimize, control, stack up, and negligible effect. Good vocabulary for calibrated but confident claims.

## 7. Figure And Table Logic
A key intuition figure for the candidate-location or layered-network formulation should appear early; tables then carry timing gains.

## 8. Reusable Lessons For Our Paper
- imitate: I
- avoid: m
