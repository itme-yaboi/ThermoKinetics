# CLAUDE.md — Veinik-based thermodynamic modeling library (Julia)

## What we are building

A Julia library for multi-domain thermodynamic modeling built on A. I. Veinik's generalized thermodynamics:

- A system is assembled from **bodies**. Each body has a set of **forms of motion** (domains), each a conjugate pair of generalized charge E and generalized force P.
- Within a body, forms of motion are **uncoupled or coupled** through the body's state function U(E).
- Bodies are connected by **conductors** (dissipative transport, possibly with coupled fluxes), **streams** (advection) and **converters** (reversible cross-domain transformation).
- Adding a new form of motion should mean declaring a new (E, P) pair, not writing a new solver.

It is used for engineering consulting (thermal and balance-of-plant analysis for fusion systems), so correctness, traceability and loss attribution (where entropy is produced) matter more than raw speed at first.

## Read before writing code

Everything needed from Veinik's books is in `docs/veinik/`. **Do not ask for or try to OCR the source books.** Start with `docs/veinik/00_README.md`, then read in order:

1. `01_core_framework.md` — charge, force, first law, equations of state, coupling, capacity
2. `02_thermokinetics.md` — fluxes, dissipation, entropy production, Onsager, couples
3. `03_unsteady_and_cycles.md` — regimes, diffusion equation, lumped relaxation, cycles
4. `05_cautions_and_errata.md` — **mandatory**: where Veinik is wrong and what to do instead
5. `07_architecture_notes.md` — recommended architecture and the reasoning behind it
6. `08_validation_cases.md` — benchmarks and invariants

`04_domain_catalog.md`, `06_modern_mapping.md` and `09_equation_index.md` are for reference and lookup.

## Non-negotiable physics

1. Every non-thermal charge is conserved by every **connection between bodies**. Changes *inside* a body (chemical reactions, radioactive decay and transmutation such as tritium decay) are explicit internal transformations that conserve elements or total mass as appropriate. **Entropy *production* is ≥ 0**, and it is the only internal source of any charge. External sources (nuclear heating, prescribed cooling) are separate, explicit boundary or source terms, and their entropy can have either sign. Never conserve entropy.
2. Energy is conserved across every connection. Dissipated energy re-enters as entropy at the temperature of the body where it is deposited.
3. Forces are defined strictly as P = ∂U/∂E, so the volumetric force is −p. No per-domain sign exceptions.
4. Static reciprocity comes from U being one potential (its Hessian is symmetric). Kinetic coefficient matrices L have a positive semidefinite symmetric part. L is symmetric unless a magnetic field is present (Onsager–Casimir).
5. Use true forces across material boundaries: μ, not concentration (the contact rule).
6. SI units only. Veinik's kgf, kcal, A = 1/427 and γ (specific weight) never appear in code.

## How to work

1. **Design first.** Before implementation, write `docs/DESIGN.md`: types, state representation (**charges-as-states is decided**, see `07` §1; design the Newton inversion for property backends), assembly, solver strategy, and how each non-negotiable above is guaranteed structurally. Stop and ask for review of DESIGN.md before building past a minimal prototype.
2. **Tests alongside code.** Implement the invariant tests in `08_validation_cases.md` §A first, then the analytic cases, in roughly the order B → C → D → E.
3. **Small vertical slices.** First slice: thermal bodies, fixed-temperature reservoirs, an energy-form thermal conductor (`07` §3), and a chain of cells for a discretized rod (B2, C1, C2). Then add the volumetric domain and an ideal gas (B3, B4). Then a stream (B5, C7). Then coupled conductors (C4, C5).
4. When a book equation and a non-negotiable disagree, the non-negotiable wins. Note the disagreement in `docs/DECISIONS.md`.
5. Keep a running `docs/DECISIONS.md`: one entry per design decision, with alternatives considered.

## Existing library (reference only)

The author has an earlier Julia thermal library (ThermalModeling / TSM). It is a conventional lumped-element code, not generalized the way this one should be. See `docs/EXISTING_LIBRARY.md` for when and how to use it. Short version: **do not read it until DESIGN.md is drafted, never mirror its architecture, and use it only for validated cases, correlations, property wrappers and lessons learned.**
