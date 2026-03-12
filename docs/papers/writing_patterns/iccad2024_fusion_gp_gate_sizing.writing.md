- paper_id: iccad2024_fusion_gp_gate_sizing
- title: Fusion of Global Placement and Gate Sizing with Differentiable Optimization
- venue_year: ICCAD 2024
- writing_role: venue_style_anchor

## 1. Title Pattern
Mechanism and method class are both in the title. It is modern, precise, and slightly more assertive than older FPGA placement titles.

## 2. Abstract Move Order
1. why the current flow split is limiting
2. new fused flow
3. differentiable mechanism
4. main metrics
5. runtime/process implication

## 3. Introduction Framing
Strong motivation arc from design-flow separation to optimization-space loss. Good sample for writing when the novelty is a fused flow.

## 4. Method Narration
Uses high-level flow change first, then the differentiable objective details. That separation is useful for readability.

## 5. Experiment Narration
Main results are metric-heavy and appear early. Good model for concise, quantitative abstracts.

## 6. Wording Habits
Useful phrases include shift-left, larger exploration space, jointly optimize, and differentiable objective.

## 7. Figure And Table Logic
A flow-comparison figure is important and likely appears very early.

## 8. Reusable Lessons For Our Paper
- imitate: Imitate how it sells a flow-level change with concrete metrics.
- avoid: Avoid borrowing ASIC-specific motivation that does not translate to FPGA placement.
