- paper_id: date2021_timing_driven_fpga_clock_constraints
- title: Timing-Driven Placement for FPGAs with Heterogeneous Architectures and Clock Constraints
- venue_year: DATE 2021
- writing_role: method_narration_anchor

## 1. Title Pattern
Title states task, platform, and the two hard constraints directly. This is a strong pattern if your paper must foreground timing and DSP/resource constraints together.

## 2. Abstract Move Order
1. modern FPGA context
2. two missing considerations
3. proposed analytical flow
4. mechanism components
5. quantitative result

## 3. Introduction Framing
The introduction gets to the technical bottleneck early and uses a contribution list with separate bullets for model, constraint handling, co-optimization, and results.

## 4. Method Narration
Method is staged: timing modeling, global placement, legalization. This staged narration is strong for papers that combine several technical pieces without losing clarity.

## 5. Experiment Narration
Setup likely names architecture and commercial baseline early, then reports timing legality and wirelength together. Main metrics appear before secondary discussion.

## 6. Wording Habits
Good verbs for this genre: resolve, optimize, satisfy, preserve, improve. Claims are strong but still metric-grounded.

## 7. Figure And Table Logic
A flow figure and architecture figure arrive early. Tables should carry quantitative comparison; equations support delay-model credibility.

## 8. Reusable Lessons For Our Paper
- imitate: I
- avoid: m
