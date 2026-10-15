# CLAUDE.md — Veinik-based thermodynamic modeling library (Julia)

## What we are building

A Julia library for multi-domain thermodynamic modeling built on A. I. Veinik's generalized thermodynamics:

- A system is assembled from **bodies**. Each body has a set of **forms of motion** (domains), each a conjugate pair of generalized charge E and generalized force P.
- Within a body, forms of motion are **uncoupled or coupled** through the body's state function U(E).
- Bodies are connected by **conductors** (dissipative transport, possibly with coupled fluxes), **streams** (advection) and **converters** (reversible cross-domain transformation).
- Adding a new form of motion should mean declaring a new (E, P) pair, not writing a new solver.

The end goal, and the yardstick for every slice, is the intended functionality in DESIGN §0.1: power cycles, transient heat transfer, control-volume and control-mass systems, process systems, circuits, mass transfer and other practical thermodynamic simulations, all built from the same few element kinds.

It is used for engineering consulting (thermal and balance-of-plant analysis for fusion systems), so correctness, traceability and loss attribution (where entropy is produced) matter more than raw speed at first.

## Read before writing code

Everything needed from Veinik's books is in `docs/veinik/` (local only, git-ignored, D-37). **Do not ask for or try to OCR the source books.** Start with `docs/veinik/00_README.md`, then read in order:

1. `01_core_framework.md` — charge, force, first law, equations of state, coupling, capacity
2. `02_thermokinetics.md` — fluxes, dissipation, entropy production, Onsager, couples
3. `03_unsteady_and_cycles.md` — regimes, diffusion equation, lumped relaxation, cycles
4. `05_cautions_and_errata.md` — **mandatory**: where Veinik is wrong and what to do instead
5. `07_architecture_notes.md` — recommended architecture and the reasoning behind it
6. `08_validation_cases.md` — benchmarks and invariants

`04_domain_catalog.md`, `06_modern_mapping.md` and `09_equation_index.md` are for reference and lookup.

## Non-negotiable physics

1. Every non-thermal charge is conserved by every **connection between bodies**, except **converters**, which transform one charge into another at a fixed ratio while conserving power with zero entropy production. Changes *inside* a body (chemical reactions, radioactive decay and transmutation such as tritium decay) are explicit internal transformations. Reactions are affinity-driven; nuclear processes are prescribed source terms that conserve nucleon number and charge (see `07` §5). **Entropy *production* is ≥ 0**, and it is the only internal source of any charge. External sources (nuclear heating, prescribed cooling) are separate, explicit boundary or source terms, and their entropy can have either sign. Never conserve entropy.
2. Energy is conserved across every connection. Conductors use the exact discrete formulation in `07` §3.1, which guarantees Ṡ_gen ≥ 0 at any finite ΔT.
3. Forces are defined strictly as P = ∂U/∂E, so the volumetric force is −p. No per-domain sign exceptions.
4. Static reciprocity comes from U being one potential (its Hessian is symmetric). Kinetic coefficient matrices L have a positive semidefinite symmetric part. L is symmetric unless a magnetic field is present (Onsager–Casimir).
5. Use true forces across material boundaries: μ, not concentration (the contact rule).
6. SI units only. Veinik's kgf, kcal, A = 1/427 and γ (specific weight) never appear in code.

## Start of every session

1. Read `docs/STATUS.md`: where the project stands, what is next, known shortcuts.
2. Read `docs/DESIGN.md` (the target architecture) and `docs/DECISIONS.md` (why things are as they are) as needed for the task. Do not re-open a decision recorded there without new information; supersede it with a new entry instead.
3. Run the tests (`julia --project=. -e 'using Pkg; Pkg.test()'`) before changing code, so you know the baseline.

DESIGN.md describes the target. STATUS.md and the code say what exists. When they differ, believe the code and fix the docs.

## How to work

1. **Design is approved** (DESIGN.md revision 2, reviewed 2026-10-14; of the revision 3 additions, mechanics §18 was accepted 2026-10-15, while rockets §19 are still proposed). Build in the vertical slices of DESIGN §15. Any change to the design itself (new element kinds, new state representation, a new kind of solver) goes into DESIGN.md and DECISIONS.md and is put to the user for review before it is built.
2. **Tests first, alongside code.** For each slice, write its tests from `08_validation_cases.md` (§A invariants, then analytic cases B → C → D → E → G) before the code, then write code until they pass. Every integration test also checks `closure` (energy and entropy balance). When a slice delivers a use from DESIGN §0.1, add its **application test** to `test/test_applications.jl` and mark the row.
3. **Keep `src/` minimal**: only what the current slice needs.
4. When a book equation and a non-negotiable disagree, the non-negotiable wins. Note the disagreement in `docs/DECISIONS.md`.
5. **DECISIONS.md is append-only**: one entry per decision, with alternatives considered. Log any deviation from DESIGN.md there. Supersede old entries rather than rewriting them.
6. **No ModelingToolkit, anywhere** (D-2).
6a. **Generality over special cases** (D-34). Never bend the framework's rules (element kinds, the converter and conductor contracts, the bookkeeping) to fit one use case. Prefer an explicit, ledgered, special-purpose element (a "virtual" one-off component such as `GravityField`, a momentum source) that enables the case within the existing rules. Record each such element in DECISIONS.md, flag it as a candidate to revisit, and say when it would need to change.
7. `docs/veinik/` belongs to the user. Change it only with the user's approval, and say exactly which lines changed.

## End of every slice (do this without being asked)

1. Full test suite passes; note the count and run time.
2. Update `docs/STATUS.md`: rewrite it to the current state (not a log).
3. Mark the slice done in DESIGN §15. Update any DESIGN text the implementation proved wrong, logging the change in DECISIONS.md.
4. Update this file if the working rules changed.
5. Commit on a feature branch (never directly on `master`/`main`), one commit per slice, with a message that lists what the slice added and its tests. Ask before merging or pushing.

## Background: the author's earlier library

The author previously built a Julia thermal library (ThermalModeling / TSM), a conventional lumped-element code built on ModelingToolkit. That dependency proved a major and nearly irreversible mistake, which is why this project uses no ModelingToolkit (D-2). The old code is **not** a reference for this project and is not consulted. If a specific piece might help (a heat-transfer or friction correlation, a property wrapper, a validated case for the regression comparison in `08` §F), ask the author rather than looking for it.

## Network plotting reference

Network plotting (DESIGN §17, track P; not implemented yet) **must** use the compact orthogonal node-placement algorithm of Freivalds & Glagoļevs, *Graph Compact Orthogonal Layout Algorithm*: `_reference/graphing/Freivalds, Graph Compact Orthogonal Layout Algorith.pdf` (12 pages, with Algorithm 1 on p. 7). Read it, and DESIGN §17, before implementing any plotting.

The paper specifies node placement only. All plotting PDFs are local only (git-ignored, D-30). Its references for the other phases are also in `_reference/graphing/`:
- Freivalds–Kikusts (`layopt.pdf`): normalization;
- IPSep-CoLa (`An Incremental Procedure for.pdf`): the QP solver for compaction and normalization;
- Dwyer–Marriott–Stuckey (`Fast node overlap removal.pdf`): non-overlap constraints;
- Biedl–Madden–Tollis (`biedl2000.pdf`): port assignment;
- Lengauer: edge routing. Local only, not in git; PDF page = printed page + 18.

DESIGN §17.1 maps each phase to its paper, section and page.

## Source material (only when needed)

Raw Veinik sources live in `_reference/veinik_source/` (only `README.md` is in git; the section notes, the OCR texts and `page_images/` are local only, D-36; the plotting references are in `_reference/graphing/`, see above): OCR text of the whole 1964 book, page images, section notes, and the Russian 1966 and 1973 texts. Use them only to check what Veinik actually wrote or to verify an equation number; `docs/veinik/` remains the working reference. Read `_reference/veinik_source/README.md` first. OCR garbles equations, so open the page image when it matters.
