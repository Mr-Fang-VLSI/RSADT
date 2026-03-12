# ICCAD Writing Guidelines From Original Paper Text

## Scope

This note is derived from direct raw-text extraction of the following core papers, rather than from downstream summaries:

- [Systolic_Array_Placement_on_FPGAs.raw.txt](/mnt/research/Hu_Jiang/Students/Fang_Donghao/RSADT/RSADT/docs/papers/summaries/raw/Systolic_Array_Placement_on_FPGAs.raw.txt)
- [3676536.3676690.raw.txt](/mnt/research/Hu_Jiang/Students/Fang_Donghao/RSADT/RSADT/docs/papers/summaries/raw/3676536.3676690.raw.txt)
- [Timing-Driven_Placement_for_FPGAs_with_Heterogeneous_Architectures_and_Clock_Constraints.raw.txt](/mnt/research/Hu_Jiang/Students/Fang_Donghao/RSADT/RSADT/docs/papers/summaries/raw/Timing-Driven_Placement_for_FPGAs_with_Heterogeneous_Architectures_and_Clock_Constraints.raw.txt)
- [3036669.3036682.raw.txt](/mnt/research/Hu_Jiang/Students/Fang_Donghao/RSADT/RSADT/docs/papers/summaries/raw/3036669.3036682.raw.txt)

## What The Original Papers Actually Do

1. They reach the technical bottleneck quickly.
The opening paragraph gives only enough application context to justify the placement problem, then immediately states the missing capability in existing methods.

2. They frame novelty as a controlled delta, not as a project diary.
`RSAD` and `SysMix` both position the new method as a specific extension of a known design setting. They do not narrate the research process or implementation history.

3. They separate problem structure from algorithm.
The papers first explain why the problem has a special structure, and only then introduce the optimization engine that exploits it.

4. They use contribution lists as reviewer-facing contracts.
The bullet lists are not generic. Each bullet names one concrete technical piece or one concrete empirical outcome.

5. They avoid source-code-level narration.
The method sections describe formulations, models, constraints, and optimization flow. They do not mention helper functions, file names, or engineering switches.

6. They open method sections with the optimization object.
The strongest method openings first state what is being optimized and why previous formulations are insufficient. Algorithm mechanics come after that.

7. They keep proof sketches short but structurally complete.
A proof sketch names the claim, the case split, the recurring inequalities, and the conclusion. It does not retell exploratory false starts.

8. They give each table one rhetorical job.
Main tables establish the primary empirical claim. Secondary tradeoffs are discussed afterward, rather than mixed into the same sentence.

## Concrete Guidance For The Current Paper

1. The paper should not mention source files, function names, executables, or implementation switches.

2. `OC` should be introduced as a structural hypothesis from prior work, then justified by a theorem-style local-optimality argument, and only after that used for formulation.

3. The method opening should start from the mathematical object:
the weighted single-column placement problem under a justified OC structure.

4. The shortest-path contribution should be phrased as a formulation result:
weighted multi-net HPWL optimization is reduced to one shortest path on a state DAG.

5. The multi-column stage should be framed as decomposition, not implementation.
State that the column partition follows the `RSAD` principle using lower and upper wirelength bounds; do not narrate wrappers or binaries.

6. The `CTR` idea should be presented as a reduction theorem or reduction principle.
The important message is: for tall strips, the problem can be reduced from `m x h` to `h x h`, then expanded by inserting a row-sweep center region.

7. The result section should keep one primary claim:
the proposed method achieves higher `Fmax` and better `WNS/TNS` while keeping `HPWL` essentially at the `RSAD` level.

8. The paper should describe `RSAD` respectfully and precisely.
`RSAD` contributes the structural insight and the OC viewpoint; the present paper contributes the proof-oriented weighted formulation and the timing-driven extension.

## Anti-Patterns To Remove

1. Do not write sentences that sound like internal project notes.

2. Do not explain the method through code artifacts.

3. Do not use phrases like `in the implementation` unless the paper is explicitly discussing engineering overhead or complexity.

4. Do not narrate abandoned proof routes if a direct proof sketch already exists.

5. Do not mix theorem-level claims with empirical timing claims in the same sentence.

## Immediate Revision Checklist

1. Rewrite the abstract so it states bottleneck, formulation, and quantitative outcome with no implementation flavor.

2. Rewrite the first two paragraphs of the introduction so they mirror the cadence of `RSAD`, `SysMix`, and the DATE/ISPD timing-driven papers.

3. Rewrite the `OC` section as a proposition-style argument with a direct three-case split for arbitrary-distance two-swaps.

4. Rewrite the method opening so the paper first introduces the justified OC-structured problem and only then explains the shortest-path conversion.
