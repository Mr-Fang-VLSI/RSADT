## Metadata
- title: Systolic Array Placement on FPGAs
- authors: unknown
- venue/year: ICCAD 2023
- local_pdf: docs/papers/pdf/Systolic_Array_Placement_on_FPGAs.pdf

## One-Paragraph Abstract (rewritten)
The paper frames systolic-array placement as a domain-specific FPGA CAD problem rather than a generic FPGA placement instance. It argues that the regular communication structure and DSP-column constraints of systolic arrays deserve a dedicated placement strategy, then presents a tailored placer and shows that the specialized flow improves placement quality on systolic-array-style designs.

## Problem and Context
General-purpose FPGA placers do not exploit the regular communication and resource structure of systolic arrays. That mismatch can leave wirelength and timing quality on the table, especially when the design is dominated by DSP-centered regular compute tiles.

## Method
- key idea: Build a systolic-array-aware placement method that respects the regular structure of the accelerator and the FPGA resource organization.
- assumptions: Assumes the target design has strong systolic regularity and that a placement flow can benefit from exploiting repeated local communication patterns.
- model/formula: No single formula is central in the extracted text; the contribution is a structure-aware placement flow rather than a compact analytical model.

## Evidence
- result-1 (abstract): The abstract claims specialized placement is needed because existing FPGA layout techniques are not matched to systolic arrays (abstract).
- result-2 (abstract/introduction): The paper positions itself as a domain-specific alternative to generic FPGA layout flows for CNN-style accelerators (abstract/introduction).
- result-3 (abstract/introduction): The evidence focus is on demonstrating that exploiting systolic regularity yields better placement outcomes than generic treatment (abstract).

## Limits and Threats to Validity
The method is closely tied to systolic-array structure, so transfer to more irregular accelerators may be limited. The extracted text does not yet expose a detailed limitation section or broader cross-architecture validation scope.

## Relevance to Current Project
- supports/contradicts/neutral: supports
- why: This is the closest predecessor and the best reference for how to motivate a specialized FPGA placement paper around systolic arrays. It also sets the baseline narrative that the current work should extend from placement quality toward timing-driven DSP placement.

## Actionable Next Step
Use this paper as the main template for problem framing, then explicitly state what the current timing-driven DSP placement method adds beyond it.
