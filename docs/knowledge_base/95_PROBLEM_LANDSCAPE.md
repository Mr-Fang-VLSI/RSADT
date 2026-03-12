# Problem Landscape For Timing-Driven DSP Placement On Systolic Arrays

## Main Problem Clusters In The Corpus

### 1. Generic FPGA placement is not specific enough
Representative papers:
- Systolic Array Placement on FPGAs
- SysMix
- A New Paradigm for FPGA Placement Without Explicit Packing

Observation:
These papers argue that standard FPGA placement abstractions lose important structure. In your case, the missing structure is the combination of systolic communication regularity, DSP-column placement constraints, and timing sensitivity.

### 2. Timing closure is not equivalent to wirelength minimization
Representative papers:
- Timing-Driven Placement for FPGAs with Heterogeneous Architectures and Clock Constraints
- An Effective Timing-Driven Detailed Placement Algorithm for FPGAs
- Timing-Driven Global Placement by Efficient Critical Path Extraction

Observation:
The literature repeatedly shows that pure wirelength or simple net weighting is not enough. A convincing ICCAD paper in this area must explain exactly how the timing signal enters the placement loop.

### 3. FPGA architectural constraints reshape the optimization problem
Representative papers:
- RippleFPGA
- Clock-Aware UltraScale FPGA Placement with Machine Learning Routability Prediction
- Timing-Driven Placement for FPGAs with Heterogeneous Architectures and Clock Constraints

Observation:
Routing congestion, clock regions, heterogeneous resources, and legality are not afterthoughts. They define what a credible placement method must preserve.

### 4. Accelerator-specific hierarchy creates a second layer of structure
Representative papers:
- SysMix
- Systolic Array Placement on FPGAs

Observation:
The design class itself creates structure. For your project, the paper should frame systolic-array regularity and DSP usage pattern as problem-defining information, not just benchmark flavor.

## Position Of The Current Project
The current project sits at the intersection of accelerator-specific placement, FPGA heterogeneous-resource constraints, timing-driven optimization, and DSP-centric structural regularity. That position is narrower and cleaner than generic FPGA placement. The paper should exploit that narrowness instead of apologizing for it.
