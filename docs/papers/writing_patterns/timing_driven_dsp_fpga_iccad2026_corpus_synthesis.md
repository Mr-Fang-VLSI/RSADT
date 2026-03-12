# FPGA Timing/Placement Writing Pattern Synthesis

## Corpus
This corpus contains 11 papers spanning ICCAD 2016-2024, DATE 2021/2025, ISPD 2017, and TCAD 2019. The strongest direct anchors for the current paper are:
- `Systolic Array Placement on FPGAs` (ICCAD 2023)
- `SysMix: Mixed-Size Placement for Systolic-Array-Based Hierarchical Designs` (ICCAD 2024)
- `Timing-Driven Placement for FPGAs with Heterogeneous Architectures and Clock Constraints` (DATE 2021)
- `An Effective Timing-Driven Detailed Placement Algorithm for FPGAs` (ISPD 2017)

## Repeated Writing Patterns
1. Titles are usually scope-first and technical, not slogan-first.
2. Abstracts usually follow: application/problem context -> concrete CAD gap -> method idea -> key mechanism pieces -> quantitative outcome.
3. Introductions reach the bottleneck quickly and do not spend many paragraphs on generic FPGA background.
4. Contributions are commonly listed explicitly, often separating model, optimization engine, legalization/constraint handling, and empirical evidence.
5. Method sections are easier to follow when staged into 3-5 blocks instead of one monolithic optimization story.
6. Results sections usually show benchmark/setup first, then primary QoR metrics, then ablations or case studies.
7. Good papers maintain calibrated verbs: `propose`, `develop`, `optimize`, `satisfy`, `preserve`, `improve`, `resolve`, `demonstrate`.

## Narrative Guidance For The Current Paper
1. Start from the mismatch between generic placement and timing-driven DSP placement for systolic arrays on FPGA.
2. Present the method as the next controlled extension after `Systolic Array Placement on FPGAs` and `SysMix`, not as an unrelated new direction.
3. Separate the story into why timing matters here, what DSP/systolic structure changes, what mechanism is introduced, how legality/resource structure is preserved, and what experiments validate the mechanism.
4. Keep the title conservative. A good title should probably name `timing-driven`, `DSP placement`, `systolic arrays`, and `FPGA` directly.
5. In the abstract, avoid leading with implementation details. Lead with the bottleneck, then the method idea, then the main measured outcome.
6. In the introduction, list contributions explicitly. Do not hide the novelty split.

## Anti-Patterns To Avoid
1. Spending too much space on generic ML-accelerator background before stating the placement bottleneck.
2. Mixing the core method contribution with engineering cleanup in the same claim sentence.
3. Overselling with words like `new paradigm` unless the formulation truly changes at that level.
4. Reporting timing gains without clarifying what is modeled or optimized differently.
5. Treating result tables as a dump. Each table must have one rhetorical job.
