# DECISIONS

One entry per design decision: what was chosen, the alternatives considered, and why. Status is **accepted** (fixed by the docs, CLAUDE.md or the reviewer), **proposed** (not yet reviewed), or **superseded**. Disagreements between Veinik's book and the non-negotiables are listed at the end.

Review history: revision 1 of `DESIGN.md` was reviewed on 2026-10-14. The answers are recorded in DESIGN §16 and reflected below.

---

### D-1 State representation: charges as states — accepted
- **Chosen:** the state vector is every body's charges (S, V, m_k, q, p⃗, …). Forces and U are functions of state.
- **Alternatives:**
  - (a) forces as states (a DAE with mass matrix K(P));
  - (b) U instead of S as the thermal state ("conservative variables").
- **Why:** this is Veinik's parameters of state (`01` §2.4), decided in CLAUDE.md and `07` §1. Charge totals are linear invariants, conserved exactly by the right-hand side and preserved to round-off by Rosenbrock and BDF integrators.
  - (a) conserves charge only to solver tolerance.
  - (b) would make energy a linear invariant, but it demotes entropy to a derived quantity.
- **Consequence (accepted):** integrated energy is conserved to integrator tolerance, not round-off. It is reported as a closure residual on every transient run (D-10). 08 A1 and 07 §7 were reworded to match.

### D-2 Build the core directly; no ModelingToolkit, anywhere — accepted (firm)
- **Chosen:** a hand-written assembler with explicit bookkeeping on OrdinaryDiffEq and NonlinearSolve, with topology held in a graph (D-17).
- **Alternatives:**
  - ModelingToolkit acausal components (as a core, a backend or an export);
  - a bond-graph or port-Hamiltonian junction core;
  - a generic DAE builder.
- **Why:** the reviewer's experience is that the previous library's ModelingToolkit dependency was a major mistake and almost irreversible. A symbolic layer also rewrites the entropy bookkeeping that must stay auditable. A power-conserving junction core makes entropy production an afterthought (`06`, `07` §9). No MTK export is planned.

### D-3 Element contract: laws return fluxes; one assembler writes derivatives — accepted
- **Chosen:** elements are pure maps from forces to a `Transfer`. There are five `Transfer` kinds, each with exactly one `apply!`, and `apply!` is the only code that writes `du`.
- **Alternatives:**
  - each element writes its own `du` contributions (the typical lumped-element design);
  - a global conservation projection after each step.
- **Why:** conservation of charge and energy and the sign of entropy production become properties of five small, exhaustively tested functions, not of every law. A projection step would hide bugs instead.

### D-4 Conductors: one bookkeeping rule, two proven heat-distribution kernels — accepted
- **Chosen:** a single `apply!` with E_k moved a → b, heat deliveries Q_i to thermal targets, S_i += Q_i/T_i, and Σ Q_i = D enforced by computing the last Q_i as a remainder. Two kernels compute the Q_i:
  - **Kernel A** (direct, `07` §3.1 generalized to weights over {a, b, h}):
    - X_Q = 1/T_b − 1/T_a and Y_k = ΔP_k Σ w_i/T_i, giving Ṡ_gen = xᵀM_sym x ≥ 0 at any ΔT;
    - used for pure thermal conductors (energy form Q̇ = GΔT) and for coupled conductors with a thermal row (weights on a and b only).
  - **Kernel B** (zero-capacity internal thermal node, eliminated analytically):
    - T_m = (Σ G_i T_i + D)/Σ G_i and Q_i = G_i(T_m − T_i);
    - Ṡ_gen = D/T_m + Σ Q_i(1/T_i − 1/T_m) ≥ 0;
    - **the default for non-thermal dissipation** (Joule, viscous, diffusion), with required leg conductances.
  - `Into(h)` and `Split(w)` (kernel A) remain explicit options for non-thermal dissipation.
- **Alternatives:**
  - an entropy-form law J_s = L·ΔT (no unique reference T);
  - Onsager at a mean T̄ (**superseded**: not ≥ 0 at large ΔT);
  - `07` §3.2's series internal node carrying the non-thermal charge (needs a zero-capacity non-thermal DOF, i.e. a DAE);
  - kernel A only, with w = ½ default (revision 1).
- **Why:** kernel B was suggested by the reviewer. It puts dissipation where it physically occurs and attributes losses in two parts (generation at the node, conduction from the node). It also stays explicit, because exactly-mapped engineering laws (Ohm, Hagen–Poiseuille on μ) make D independent of T_m. Kernel A covers coupled thermal rows, where D can be negative, and gives the Kelvin relations without a T̄.

### D-5 Conductance parameterization: M = BBᵀ + W by type — accepted
- **Chosen:**
  - `Conductance(B, W)`, with W antisymmetric and settable only through a `MagneticField` context;
  - `ScalarConductance(g ≥ 0)`, which errors rather than clamping;
  - coupled data rescaled only by diagonal congruence.
- **Alternatives:**
  - accept any matrix and check PSD at runtime;
  - check at construction only (misses state-dependent laws);
  - clamp negatives (hides bugs).
- **Why:** non-negotiable 4 and the sign of Ṡ_gen become type invariants. Onsager–Casimir is supported without hard-coding L = Lᵀ.

### D-6 Converters: explicit antisymmetric (gyrator) form now — accepted
- **Chosen:** F = W P with constant antisymmetric W, so Pᵀ W P ≡ 0. A kinetic DOF on one side keeps the system an ODE. Real converters are an ideal converter plus loss conductors.
- **Alternatives:**
  - transformer form (a constraint, so a DAE) now;
  - a large-stiffness penalty;
  - modeling everything as conductors.
- **Why:** explicit, power-exact, and it matches the reversible half of GENERIC. The transformer form follows in D-19.

### D-7 Property backends: Clapeyron.jl first; 1-D safeguarded Newton on T with IFT derivatives — accepted
- **Chosen:**
  - v = V/m, s = S/m and the composition come from homogeneity;
  - solve s(T; v, x) = s with a bracketed, warm-started Newton;
  - Dual-number derivatives via the implicit-function theorem;
  - an (s, v, x)-flash for two-phase;
  - an isentropic (s, p) → v inversion for pumps;
  - Clapeyron.jl as the first backend, with a CoolProp wrapper later.
- **Alternatives:**
  - a 2-D Newton on (T, ρ);
  - AD through the iterations;
  - tabulated inverse maps;
  - forces as states;
  - CoolProp first (broader fluid set, but needs derivative and flash wrappers).
- **Why:** 1-D and monotone means robust and cheap. IFT gives exact Jacobians. Clapeyron's Helmholtz models are AD-friendly, which fits this scheme directly.

### D-8 Capacities, stability and c_p/c_v on the fixed-mass block — accepted
- **Chosen:** these work on the Hessian block that excludes `Species` charges. The full Hessian is never inverted.
- **Alternative:** per-unit-mass variables throughout (equivalent).
- **Why:** U is homogeneous of degree 1 once mass is a charge, so the full Hessian is singular (`01` §4.4).

### D-9 Solvers: Rodas5P default, FBDF/QNDF for large sparse problems; NonlinearSolve for steady state — accepted
- **Alternatives:**
  - explicit methods (unsuitable for stiff thermal networks);
  - implicit Euler only;
  - Sundials CVODE/IDA.
- **Why:**
  - Rosenbrock and BDF methods preserve linear invariants given 1ᵀJ = 0.
  - Both handle index-1 mass-matrix DAEs, so D-19 needs no solver change.
  - The sparsity pattern is known from the graph.

### D-10 Ledgers: non-feedback states, separable from the physics — accepted
- **Chosen:** Ṡ_gen per element (with parts split), plus reservoir, source, pump and boundary flows, are appended as u = [u_phys; q_ledger].
  - They are excluded from the integrator's error norm (`internalnorm` over u_phys only) and from steady-state solves (f_phys = 0 only; ledger rates reported).
  - Every transient run reports the energy, entropy and charge closure residuals.
- **Alternatives:**
  - post-hoc quadrature on saved solutions;
  - step callbacks that accumulate.
- **Why:** they integrate at the integrator's accuracy, never restrict steps, never make a steady state ill-posed, and the closure residuals become first-class outputs.

### D-11 Internal transformations — accepted
- **Chosen:**
  - **Reactions:** J_r = M_r(A/T) with M_r PSD; dm = νJ; dS = A·J/T; dU = 0; atom balance validated.
  - **Nuclear processes:** prescribed-rate sources, with nucleon number and charge validated; energy deposited as Q̇/T at the local T; booked as external.
- **Alternatives:**
  - decay with an affinity (wrong, being far from equilibrium);
  - enforcing elements or mass for nuclear processes (violated by T → ³He).
- **Why:** CLAUDE.md non-negotiable 1 and `07` §5.

### D-12 Coupling expressed by material composition — accepted
- **Chosen:** a `Composite` of parts on disjoint domain subsets, each part one coupling group. The Hessian is block-diagonal by part, and at most one part contains `Thermal`.
- **Alternatives:**
  - per-body coupling flags;
  - a dense U with zero detection by AD;
  - coupling declared on the domain.
- **Why:** Veinik's groups become structural, and A7 holds by construction.

### D-13 Flow drivers, pumps and boundaries — accepted
- **Chosen:**
  - Between bodies, `Stream` ṁ comes only from sign-preserving drivers (ṁ = κΔp, κ ≥ 0), which deposit nothing.
  - `IdealPump(u, d, ṁ_spec)` sets ṁ regardless of Δp, raising the fluid isentropically to p_d and logging Ẇ = ṁ(h(s_u, p_d) − h_u) as external work. For liquids that is ṁvΔp.
  - `FlowBoundary` requires `supply = :passive | :pumped`; unlabeled is an error.
  - A passive boundary with p_in < p_d is a runtime error.
  - `PassThroughOutlet` is external and has no generation.
- **Alternatives:**
  - prescribed ṁ anywhere with a runtime sign check;
  - logging negative inlet generation as a warning (revision 1's proposal);
  - pumps only as a motor converter plus losses (needs the momentum domain and slice 5);
  - ṁvΔp for all fluids (inexact for gases; the exact isentropic Δh reduces to it for liquids).
- **Why:** the reviewer asked for labeled boundaries and a work-supplied pump so pumped loops work now. Generation stays ≥ 0 structurally, because the pump delivers at p_d and its work is ledgered, not hidden. The motor-plus-converter pump remains the physical long-term model.

### D-14 Switching via gates on element outputs — accepted
- **Chosen:** g ∈ {0, 1} as a parameter scaling the `Transfer`, changed by callbacks.
- **Alternatives:**
  - rebuilding the system at each switch;
  - smooth gates as states.
- **Why:** g·M stays PSD and g·W stays antisymmetric, and no rebuild is needed.

### D-15 Units: Float64 SI internally; Unitful only at constructor boundaries — accepted
- **Alternatives:** Unitful throughout (slow, awkward with AD); no unit support.
- **Why:** non-negotiable 6, plus speed and AD compatibility.

### D-16 Domains as values with traits — accepted
- **Chosen:** a small struct per domain, with traits for symbols, units, `ncomponents` and `producible`.
- **Alternatives:** symbols only; a parametric `Domain{name}`.
- **Why:** adding a form of motion means a struct plus traits, with no solver change.

### D-17 Topology as a bipartite graph (Graphs.jl + MetaGraphsNext.jl), compiled to flat arrays — proposed
- **Chosen:**
  - Vertices are nodes (bodies, reservoirs) and elements. Edges are ports, carrying role, domains and direction as metadata.
  - The graph is used for building, validation, queries (connected components, sparsity, adjacency) and export.
  - `compile` lowers it to type-stable flat arrays, and the graph is never touched in the right-hand side.
- **Alternatives:**
  - elements as edges of a node graph (cannot hold parallel elements or ports beyond two without a multigraph or hypergraph);
  - a plain vector-of-elements registry (loses topology queries);
  - MTK's connection graph (excluded by D-2).
- **Why:** the reviewer recommends graph packages, specifically MetaGraphs. The bipartite form handles parallel elements and multi-port elements uniformly. Compiling to plain data keeps the hot loop fast and the representation inspectable, with no generated code. Package versions are to be checked when the slice-1 environment is created.

### D-18 Repository hygiene — accepted
- **Chosen:**
  - Keep the reviewer's `.gitignore` excluding `_reference/` (scanned copyrighted sources).
  - Add `reference/` (the planned snapshot of the old library, per `docs/EXISTING_LIBRARY.md`, which `_reference/veinik_source/README.md` also names) and `Manifest.toml` (this is a package).
- **Alternatives:** tracking the sources with git-excluded paths only via `.git/info/exclude` (not shared across clones).
- **Why:** the sources must never be committed, and the rule should travel with the repository.

### D-19 Transformer converters and zero-capacity non-thermal nodes via mass-matrix DAE — accepted (scheduled, slice 5)
- **Chosen:**
  - 𝕄 du/dt = f(u) with singular 𝕄. One algebraic row and a flux variable λ per constraint.
  - A transformer contributes F₁ = λ and F₂ = −rλ under P₁ = rP₂, so power is zero on the constraint.
  - This must be in place before any balance-of-plant cycle work.
- **Alternatives:**
  - penalty or stiff-spring approximations;
  - always inserting a kinetic DOF (D-6 only);
  - an index-reduction framework (excluded with MTK).
- **Why:** the reviewer requires transformer converters for cycle work. Rodas5P and FBDF already handle index-1 mass-matrix DAEs, so the core and solver stay unchanged.

### D-20 Slice 1 implementation: deviations from DESIGN.md — accepted (slice 1, 2026-10-14)
Each item lists what was built instead of the design text, and why.

1. **Tolerance scaling uses capacity, not charge magnitude.**
   - DESIGN §11 said "abstol per state from the initial charge magnitudes".
   - Built instead: abstol for each charge = `atol_T`/A_kk, i.e. the charge equivalent of `atol_T` kelvin at the initial state (`abstol_vector`).
   - Why: entropy has an arbitrary zero (it depends on T_ref), so |S| is not a meaningful scale. **DESIGN §11 should be read this way from now on.**
   - The relative tolerance still multiplies |S|. The default reltol = 1e-8 keeps that harmless, but it is reference-dependent; revisit when real-fluid reference states arrive.
2. **Conductor kernels: thermal row of kernel A only.**
   - DESIGN §15 put "the general kernel A/B code with identity tests" in slice 1.
   - The reviewer asked for a minimal `src/`, and no non-thermal domain exists yet to use non-thermal rows or kernel B. Both move to slice 4, where they are first used.
   - The bookkeeping shape is already the general one: heat deliveries Q_i, the remainder rule Σ Q_i = D, and S_gen = Σ Q_i/T_i.
3. **`Reservoir(name; T)` is thermal-only.**
   - Not the general `Reservoir(name, forces; material)`. That arrives with the second domain (slice 2) and with streams (slice 3).
4. **`steady_state` refuses isolated components.**
   - The constraint-row treatment of DESIGN §12.2 is not implemented. `steady_state` throws an `ArgumentError` naming the unanchored bodies. Isolated networks are equilibrated by time integration (A6).
5. **No `valid_state` / `isoutofdomain` yet.**
   - `ConstantCapacitySolid` has T = T_ref·exp(S/mc) > 0 for every S, so no state can be invalid. These arrive with the first material that has a restricted domain (the ideal gas, slice 2).
6. **No PreallocationTools.**
   - The node-temperature buffer is preallocated for Float64 evaluations, and allocated per call when the integrator differentiates with Dual numbers. This removes a dependency and avoids chunk-size handling.
   - Revisit if profiling shows allocation pressure.
7. **Sparse Jacobian prototype only above 50 physical states.**
   - Below that, dense AD is faster and simpler. The `sparse` keyword overrides this. D4 (300 cells) exercises the sparse path.
8. **Test-case adaptations** (no non-thermal domain yet):
   - **D1** is run as its thermal analogue: the entropy produced while relaxing to a reservoir is C[ln(T_r/T₀) + (T₀ − T_r)/T_r], independent of G (`03` §3). D1 proper and D2 come in slice 4.
   - **D3** is covered by D4 for the thermal domain only.
   - **D5** uses a slab, not a sphere, because the rod helper is Cartesian. Criterion: max |θ̄_distributed − θ_lumped| < 0.05 for Bi ≤ 0.1 (measured 0.006 at Bi = 0.05 and 0.012 at Bi = 0.1), and > 0.05 at Bi = 10 (measured 0.45).
   - **A2** is vacuous (no non-thermal charges). **A4** is trivially 1×1. **A5** is the `ThermalConductance` constructor rejecting G < 0 or NaN. **A7** waits for slice 2.
9. **`closure` measures from `sol.prob.u0`**, not from the first saved point, so it is correct when `saveat` omits t = 0. This was a bug caught during slice 1; D4 is the regression test.

### D-21 The earlier library is not consulted — accepted (2026-10-14)
- **Chosen:**
  - The author's earlier ThermalModeling/TSM library is not read or reviewed.
  - The review phase in `docs/EXISTING_LIBRARY.md` is dropped.
  - CLAUDE.md keeps only the background: it was a ModelingToolkit-based lumped-element code, and that dependency was the main lesson (D-2).
  - Specific pieces (correlations, property wrappers, validated cases for `08` §F) are requested from the author when needed.
- **Alternatives:** the phased review planned in `EXISTING_LIBRARY.md` (inventory, reuse list, do-not-carry list, mapping table).
- **Why:** the architecture is settled and shares nothing with the old code. The main architectural lesson is already recorded. Asking for a specific correlation or case costs less than a full review, and it avoids anchoring on the old structure.
- `docs/EXISTING_LIBRARY.md` stays in the repository as history. Its phase rules no longer apply.

---

## Book versus non-negotiables (CLAUDE.md "How to work" 4)

| Book | Non-negotiable | Resolution |
|---|---|---|
| Volume work dQ = −p dV with p as the force; sign flag on that domain (§2.24) | 3 | Force is −p = ∂U/∂V. No sign flags. The 1966 reciprocity sign flip disappears. |
| (135) capacity-weighted average applied to the thermal chain, conserving entropy on joining | 1 | Joining conserves energy and non-thermal charges, and S follows from the state (B2). |
| Lumped entropy conduction at finite ΔT; Onsager at one T | 1, 2 | Exact discrete kernels (D-4). |
| (78) μ = BCT; concentration as a driving force (163) | 5 | μ from U (ideal: μ₀ + RT ln C). Fick data converted to μ-based L. |
| (248)/(257) cross-coefficients from the static EoS | 4 | L_ij are material data; only symmetry and PSD are enforced. |
| Darcy, Poiseuille and barrier flow as volume-charge transport driven by p (§27, §49, (291), (316)) | 3, 5 | Mass transport (streams or mass conductors on μ). C4 and C6 coefficients are transformed to the mass basis, and the transformation is tested. |
| (217)/(494) irreversibility via ΔP/P₁ | — (needs an absolute zero) | Report Ṡ_gen and exergy destruction instead. |
| 1973 book: entropy conservation | 1 | Excluded. |
| kgf, kcal, A = 1/427, γ | 6 | SI only; ρ instead of γ. |
