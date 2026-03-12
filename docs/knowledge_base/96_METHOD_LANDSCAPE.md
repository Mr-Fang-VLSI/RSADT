# Method Landscape For Timing-Driven DSP Placement On Systolic Arrays

## Main Method Families In The Corpus

### 1. Structure-aware specialized placement
Representative papers:
- Systolic Array Placement on FPGAs
- SysMix

Pattern:
Use domain structure to reshape placement rather than only reweight generic objectives.

### 2. Timing-aware optimization via explicit path or delay modeling
Representative papers:
- Timing-Driven Placement for FPGAs with Heterogeneous Architectures and Clock Constraints
- An Effective Timing-Driven Detailed Placement Algorithm for FPGAs
- Timing-Driven Global Placement by Efficient Critical Path Extraction

Pattern:
Move beyond plain net weights by injecting timing through delay models, critical-path handling, or path-aware extraction.

### 3. Constraint-aware placement with legality/resource reasoning inside the flow
Representative papers:
- RippleFPGA
- Clock-Aware UltraScale FPGA Placement with Machine Learning Routability Prediction
- Timing-Driven Placement for FPGAs with Heterogeneous Architectures and Clock Constraints

Pattern:
Treat congestion, clocking, or heterogeneous-site legality as co-equal concerns with quality optimization.

### 4. Flow reformulation or flow coupling
Representative papers:
- A New Paradigm for FPGA Placement Without Explicit Packing
- Fusion of Global Placement and Gate Sizing with Differentiable Optimization

Pattern:
Change where optimization authority lives in the flow when the standard stage split is too restrictive.

### 5. Scalability/runtime-centered flow engineering
Representative papers:
- UTPlaceF 3.0

Pattern:
Justify engineering design by empirical bottleneck analysis, then present targeted performance-preserving acceleration.

## Method Positioning Guidance For The Current Project
The strongest positioning for the current paper is likely specialized placement for systolic-array/DSP structure, plus explicit timing-driven guidance, while preserving FPGA legality and structural resource constraints. A weaker positioning would be to present the work as only another wirelength-focused specialized placer. The literature already supports the claim that timing needs its own modeling story.
