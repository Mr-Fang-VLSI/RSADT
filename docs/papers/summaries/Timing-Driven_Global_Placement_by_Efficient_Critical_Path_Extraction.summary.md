## Metadata
- title: Timing-Driven Global Placement by Efficient Critical Path Extraction
- authors: unknown
- venue/year: DATE 2025
- local_pdf: docs/papers/pdf/Timing-Driven_Global_Placement_by_Efficient_Critical_Path_Extraction.pdf

## One-Paragraph Abstract (rewritten)
This work revisits timing optimization during global placement and argues that pin-level net reweighting misses path-level structure. It proposes an efficient critical-path extraction view for timing-driven global placement, combining path sensitivity with a placement-friendly optimization flow and reporting strong TNS and WNS gains.

## Problem and Context
Timing-driven global placement remains difficult because fast net-level proxies often miss path-level critical structure, while richer path-based methods can be too expensive. The paper targets that modeling gap.

## Method
- key idea: Introduce a more path-aware timing-driven global placement scheme centered on efficient critical path extraction.
- assumptions: Assumes path information can be extracted efficiently enough to guide placement without collapsing runtime.
- model/formula: The exact formulation is not expanded in the extracted text, but the method is clearly framed as path-aware timing-driven optimization.

## Evidence
- result-1 (abstract): The abstract explicitly contrasts simple pin-level weighting with path-based timing information (abstract).
- result-2 (abstract/introduction): The reported evidence emphasizes TNS and WNS gains, which is useful for learning result narration in timing papers (abstract).
- result-3 (abstract/introduction): The framing is current and stylistically relevant even though the domain is not FPGA-specific (abstract).

## Limits and Threats to Validity
This is not an FPGA paper, so architecture-specific resource constraints and DSP-placement structure are outside its scope. It is more useful for writing style than direct method borrowing.

## Relevance to Current Project
- supports/contradicts/neutral: neutral
- why: Use this paper mainly as a recent timing-driven writing sample: abstract pacing, metric ordering, and how to justify path-aware timing modeling. It is not the main related-work anchor for FPGA-specific claims.

## Actionable Next Step
Mine this paper for wording and result-presentation patterns, not for FPGA-specific problem framing.
