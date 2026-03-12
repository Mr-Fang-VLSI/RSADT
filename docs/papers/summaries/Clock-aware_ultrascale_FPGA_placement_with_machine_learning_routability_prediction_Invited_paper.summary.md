## Metadata
- title: Clock-Aware UltraScale FPGA Placement with Machine Learning Routability Prediction (Invited Paper)
- authors: unknown
- venue/year: ICCAD 2017
- local_pdf: docs/papers/pdf/Clock-aware_ultrascale_FPGA_placement_with_machine_learning_routability_prediction_Invited_paper.pdf

## One-Paragraph Abstract (rewritten)
This paper targets FPGA placement under complex clocking architecture constraints and augments the flow with machine-learning-based routability prediction. It is a good example of a constraint-heavy FPGA placement paper where architecture legality and predictive guidance are both part of the main story.

## Problem and Context
Modern FPGA clocking architecture and routability demands cannot be handled well by placement flows that optimize only generic wirelength. Constraint violations can dominate outcome quality.

## Method
- key idea: Integrate clock-aware placement with predictive routability guidance so the placer accounts for both architectural and routing considerations.
- assumptions: Assumes clock architecture and routability pressure materially shape good placement quality on UltraScale devices.
- model/formula: The extracted text emphasizes predictive guidance and architecture awareness rather than one compact mathematical model.

## Evidence
- result-1 (abstract): The title itself shows a two-part contribution: architecture awareness plus ML-based routability prediction (title).
- result-2 (abstract/introduction): The abstract positions the work around meeting timing while satisfying clock constraints and routing concerns (abstract).
- result-3 (abstract/introduction): It is a useful rhetorical sample for papers that combine a core optimization flow with an auxiliary prediction model (abstract).

## Limits and Threats to Validity
The paper is not centered on systolic arrays or DSP placement, and the predictive component may be more auxiliary than core for your own work.

## Relevance to Current Project
- supports/contradicts/neutral: neutral
- why: Use this primarily as a structure sample for combining hard architecture constraints with a supplementary modeling component.

## Actionable Next Step
Borrow only the multi-part storytelling pattern if your current method also has a main optimization core plus a supporting model or heuristic.
