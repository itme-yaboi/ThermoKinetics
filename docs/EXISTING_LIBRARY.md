# Using the existing ThermalModeling / TSM library

The earlier library works and has been used on real consulting jobs. But it is a standard lumped-element (heat-rate-based) design, and the new library has a different foundation. It is valuable as **evidence and parts**, not as a **template**.

## Phase rules

| Phase | Access | What to produce |
|---|---|---|
| 1. Design | **No access.** Work only from `docs/veinik/`. | `docs/DESIGN.md` drafted from first principles |
| 2. Review | Read-only access to a **copy** of the old library | `docs/EXISTING_LIBRARY_REVIEW.md` (below) |
| 3. Build | Read-only, consult as needed | Port specific pieces only through the review's "reuse" list |
| 4. Validate | Run both libraries | Regression comparison (`08_validation_cases.md` §F) |

Why the delay: an agent that reads a working codebase first tends to reproduce its structure. Designing first and reading second keeps the old code from anchoring the new architecture while still mining it for value.

## What `EXISTING_LIBRARY_REVIEW.md` should contain

1. **Inventory**: modules, element types, solvers, property and correlation code, I/O, plotting, tests and examples.
2. **Reuse list** (port with attribution, adapted to the new types):
   - heat-transfer and friction correlations (Nusselt, friction factor, two-phase multipliers, …)
   - property-library wrappers and unit handling
   - validated example networks and any measured or benchmark data
   - solver settings that proved robust on stiff problems
   - I/O formats clients already receive
3. **Do-not-carry list**: architectural choices that conflict with the non-negotiables. Examples to look for: heat rate as the thermal "flow" with no entropy ledger; per-element sign conventions; convection terms that do not expose mixing generation; concentration-driven mass transfer across material boundaries.
4. **Lessons**: pain points found in the code (workarounds, TODOs, special cases). Each one is a requirement the new design should meet natively.
5. **Mapping table**: old element type → new construct (body, conductor, stream, converter, reservoir, source).

## Safety

- Give access to a separate clone or snapshot, never the working copy. For example, clone it into `reference/ThermalModeling/` and add that folder to `.gitignore`.
- Strip any client-confidential models, data or reports from the copy first.
