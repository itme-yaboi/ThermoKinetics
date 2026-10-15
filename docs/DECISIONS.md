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

### D-22 Slice 2 implementation: deviations from DESIGN.md and D-20 — accepted (slice 2, 2026-10-14)
1. **Kernel B, `Split`/`Into` and non-thermal rows are built now**, not in slice 4 (D-20 §2). The volume conductor is the first non-thermal conductor, and DESIGN §5.3 gives it kernel B.
   - All three policies go through one bookkeeping path: charges moved port to port, heat Q_i to thermal targets, remainder rule.
   - Coupled conductors (thermal + non-thermal rows) are rejected at construction until slice 4.
2. **Kernel B ledger in a sign-exact form.**
   - The ledger records D/T_m (dissipation) + Σ G_i (T_m − T_i)²/(T_i T_m) (conduction), which is ≥ 0 in floating point.
   - The entropy written is Σ Q_i/T_i. The two agree to round-off, and entropy closure stays at ~1e-15.
   - Alternative: ledger = written entropy. Rejected: it can come out at −eps when T_a ≈ T_b, which would break the per-element A3 check for no physical reason.
3. **`ElasticVolume(C; V0, p0)` added** (not in DESIGN §2.3, now listed there).
   - It is an uncoupled volumetric store: Veinik's linear one-DOF body.
   - It is needed for A7 (a separable composite), for the §11 check that an uncoupled body is a pure accumulator, and for D1 proper.
   - Its docstring says it is not a gas.
4. **D1 proper is delivered now** in the volumetric domain: exponential relaxation, and ½K X₀² dissipated for any conductance. D-20 §8's thermal analogue stays. D2 (capacitors) still waits for `Electrical` in slice 4.
5. **Tolerances are per domain.**
   - `default_force_atol(domain)`: 1e-8 K, 1e-4 Pa, 1e-6 J/kg. Each is scaled through 1/A_kk.
   - `simulate`'s `atol_T` keyword is replaced by `atol_scale` (API change).
6. **Charge-closure threshold is 1e-10, not 1e-12.**
   - With reservoirs, Σ body charge + reservoir ledger drifted up to ~8e-12 relative.
   - That is round-off amplified by the conditioning of each Rosenbrock linear solve, not integrator error, which would be ~1e-8.
   - DESIGN §12.1 is annotated.
7. **Reservoir ledgers:** one slot per domain (the charge received; entropy for `Thermal`), then energy. Energy includes the P·J carried by charges as well as heat. `Reservoir(name, domain => force, …)` and `Reservoir(name; T, p)` replace the thermal-only constructor (closes D-20 §3).
8. **`valid_state` and `isoutofdomain` step rejection** are implemented (closes D-20 §5). `Body` also refuses invalid initial charges.
9. **Element groups are sorted by type name.**
   - The concrete `CompiledSystem` type, and so the compiled solver, no longer depends on insertion order.
   - First-solve compile latency for a new combination of element types is still several seconds; it dominates the ~75 s test run.
   - Alternatives for later: type-erased element storage (FunctionWrappers) or precompile workloads.
10. **No mass conductor yet.** `Species` exists, and A2 covers species mass, but nothing transports mass until slice 3 (streams) or slice 4 (mass conductors).

### D-23 Network plotting: Freivalds–Glagoļevs compact orthogonal layout — accepted (user, 2026-10-14); details proposed
- **Chosen (user):**
  - Basic network plotting is part of the intended functionality.
  - The layout implements the node-placement algorithm of Freivalds & Glagoļevs, *Graph Compact Orthogonal Layout Algorithm*. The local copy is `_reference/graphing/`, alongside Lengauer's book for routing and compaction.
  - Not implemented yet; specified in DESIGN §17 as plotting track P.
- **Proposed details** (DESIGN §17), to settle when the track is built:
  - lay out the full bipartite graph (elements as small boxes), with an option to collapse two-port elements into edges;
  - display metadata in a side-table on `ThermoSystem`, never read by `compile`;
  - the layout in core Julia, rendering through a Makie package extension;
  - edge routing by grid maze/A* with bend penalty;
  - normalization and the eq. (3) QP via a separation-constraint solver.
- **Alternatives:**
  - general-purpose layouts (spring or stress layouts from NetworkLayout.jl or GraphMakie) — not orthogonal, so poor as engineering schematics;
  - the topology–shape–metrics approach (Tamassia) — minimizes bends but gives poor area and edge length on practical graphs, per the paper's §1;
  - Graphviz `dot` with orthogonal splines — an external binary, with little control over node sizes and pinning;
  - plotting inside the core package — forces a plotting stack on every user.
- **Why:** orthogonal, compact layouts with prescribed box sizes match how thermal and process schematics are drawn. The algorithm handles mixed node sizes, which the bipartite element graph needs. Keeping layout separate from rendering makes it testable and keeps the core light.

### D-24 Track `_reference/` in git, except page images — accepted (user, 2026-10-14); supersedes part of D-18
- **Chosen:**
  - `.gitignore` no longer excludes `_reference/`. It excludes `_reference/veinik_source/page_images/`.
  - The graphing references (Freivalds–Glagoļevs; Lengauer) and the Veinik OCR text and notes become trackable, so the paths cited in DESIGN §17 and CLAUDE.md resolve in every clone.
  - D-18's other entries (`reference/`, `Manifest.toml`) stand.
- **Alternatives:**
  - keep `_reference/` ignored, so citations point at files only this machine has;
  - Git LFS for the large PDFs (Lengauer is ~20 MB);
  - track only a manifest of references (titles, sources), not the files.
- **Caution, raised with the user before any `_reference` file was staged:**
  - `_reference/` holds scanned or OCR'd copyrighted works (Veinik 1964; Lengauer 1990) whose own README says "keep out of git".
  - Committing them puts copies in every clone and in history permanently, including if the repository is ever pushed or shared.
  - The ~20 MB PDF also stays in history even if deleted later.
  - Staging those files waits for the user's explicit go-ahead.

### D-25 Which `_reference/` files are committed — accepted (user, 2026-10-14); narrows D-24
- **Chosen (option b):**
  - Commit the Freivalds–Glagoļevs paper, `_reference/veinik_source/README.md` and `notes_by_section_ch6-18.md`.
  - Leave the books local and git-ignored: Lengauer's textbook, the 1964 OCR texts, and the 1966/1973 Russian texts. `page_images/` stays ignored too (D-24).
- **Alternatives:** commit everything (D-24 as first written); Git LFS for the large PDFs.
- **Why:** copyrighted books should not be copied into every clone and into permanent history.
- **Consequence:** in a fresh clone, the Lengauer citations in DESIGN §17.1 and the Veinik source paths resolve only after the author copies the books in locally. The paper, the notes and `docs/veinik/` are always available.

### D-26 Slice 3 implementation: deviations from DESIGN.md — accepted (slice 3, 2026-10-14)
1. **Element groups held in a `Vector{Any}`, not a `Tuple`** (supersedes D-22 §9 and DESIGN §11 step 3; DESIGN updated).
   - Each group is still a concretely typed vector with a specialized loop (a function barrier).
   - `CompiledSystem` is no longer parametric, so every model shares one compiled right-hand side and one compiled solver.
   - Measured: a new model with already-seen element types went from a full solver recompile (~20–30 s) to ~0.01 s, and the test suite from ~10 min to the time recorded in STATUS.
   - Cost: one dynamic dispatch per group per call, negligible next to the AD Jacobian.
   - Alternative: FunctionWrappers or full type erasure per element — not needed.
2. **Passive-boundary admissibility is checked on accepted states, not inside the right-hand side.**
   - Solver trial points (Rosenbrock stages, Newton iterates) can momentarily put a node above the supply pressure.
   - The check (`passive_violations`) runs in `rates`, rejects steps through `isoutofdomain`, and is re-checked on the final steady state.
   - It raises the explicit error named in DESIGN §7.3.
3. **`steady_state` changes:**
   - It solves only the *active* charges (those some element writes); a rigid vessel's V is held fixed.
   - It checks anchoring per conserved quantity: energy, and each charge moved by an element.
   - Inadmissible trial states return NaN, and `TrustRegion` is the fallback after Newton.
   - Cold starts of strongly coupled flow networks (the 2 × N-cell heat exchanger) can still fail. The documented route is a transient warm start. `PseudoTransient` was tried and rejected: it failed too, and cost ~50 s.
4. **Stream production is computed as [D + J_s (T_up − T_dn)]/T_dn.**
   - It is ≥ 0 in exact arithmetic (concavity, flow down the pressure gradient) but not sign-exact in floating point, because the advected terms carry |s| with an arbitrary reference.
   - Tests allow −1e-12 absolute per element, and entropy closure ≤ 1e-9 (vs 1e-10 without streams).
5. **Single species only** for `IdealPump` and `FlowBoundary`. `isentropic_state` exists only for `IdealGas`; other materials throw. Streams themselves carry any species set, provided both ends hold the same one.
6. **Material reservoirs:** `Reservoir(name; T, p, material)` fixes μ(T, p) and s. `charges_at(::IdealGas)` accepts T and any two of m, p, V.
7. **Outlets and boundaries:**
   - `PassThroughOutlet` uses ṁ = κ max(p − p_out, 0), or a prescribed ṁ(t). There is no backflow, since an outlet cannot supply a state.
   - Boundaries and outlets attach to bodies only.
   - A turbulent (√Δp) driver is not built yet.
8. **`closure` charge scale:** the maximum over the run of the bodies' Σ|charge|, not its initial value. A tank filled from 1 Pa starts with 1e-7 kg, which made round-off look like 3e-8.
9. **Heat-exchanger application test:** checked as first-order convergence (error ratio N = 10 → 20 within 1.8–2.2, Richardson value within 2e-3 of ε–NTU), not a fixed tolerance. Measured errors: −0.042, −0.021, −0.011 at N = 10, 20, 40.
10. **Still deferred:** inertance-driven flows (accepted earlier), multi-species pumps and boundaries, real-fluid backends.

### D-27 Slice 4 implementation: deviations from DESIGN.md — accepted (slice 4, 2026-10-14)
1. **Coupled laws supply a factor, not a matrix.**
   - `coupled_factor(law, …) → (B, W)` at the current state, with J = (BBᵀ + W) x on the exact discrete forces x = (X_Q, ΔP c).
   - The ledger records ‖Bᵀx‖², which is sign-exact in floating point. The entropy written, Σ Q_i/T_i, equals it to round-off.
   - This is DESIGN §5.2's `Conductance(B, W)` in functional form.
2. **Coupled conductors default to `Split(½)` and accept only `Split`.**
   - With coupling, D can be negative (thermoelectric generation), so a third node is not allowed (DESIGN §5.1).
   - The internal node (kernel B) is not offered for coupled laws.
3. **Exactly calibrated coupled laws** (beyond DESIGN §5.2's "matrix at a state-dependent temperature"):
   - **`ThermoelectricConductance(R, K, α)`** reproduces J_e = (Δφ + αΔT)/R at any ΔT. Its PSD margin is exactly K T_a T_b, and its Peltier heat is αT̄ J_e.
     - A single lumped element then reproduces the classical TEG balance Q_h = KΔT + αT_h I − ½I²R exactly at steady state.
     - A discretized constant-α leg has exactly the Joule-only temperature profile, which is Kelvin's first relation with zero Thomson heat.
   - **`KnudsenBarrier(gas; g, K)`** has zero mass flux exactly at p_a/√T_a = p_b/√T_b. Near equilibrium it equals Veinik's (291)–(299) barrier transformed to the mass basis, with ξ = g v²/T.
4. **Basis transformations are full congruences** M = A L Aᵀ, not only diagonal S L S (DESIGN §5.2 updated). The (J_v, J_s) → (heat, mass) map needs A = [1/v 0; Ts/v T], which is lower-triangular. Congruence preserves symmetry, PSD and σ. Provided as `barrier_mass_basis`, `electrokinetic_mass_basis` and `thermoelectric_from_veinik`.
5. **`OnsagerConductance(rows, M; W, magnetic)`:**
   - Validated symmetric (1e-12 relative) and PSD (eigenvalues ≥ −1e-12·scale). Eigenvalues inside that validation tolerance are set to zero when B is formed; this is the only clamp, and it is round-off sized.
   - W must be exactly antisymmetric and requires `magnetic = true`.
   - Without a thermal row it is an engineering multi-row conductor, J = (M + W) ΔP.
6. **`LinearConductance(domain, L)` replaces the `VolumeConductance` struct.** `VolumeConductance`, `Resistor(R)` (R > 0) and `MassConductance(species, L)` are constructors. A mass conductor moves mass at μ without entropy (isothermal diffusion); non-isothermal mass transport with transport heat uses a coupled law (Knudsen, or `OnsagerConductance` with a thermal row).
7. **`steady_state` scales each residual row by its round-off floor**, Σ_j |∂r_k/∂x_j| s_j with s_j = capacity × max(|P_j|, 1).
   - It stops at a scaled residual of 1e-14, or accepts a stall at ≤ 1e-10.
   - Without this, a 1e6 W/K tie (realistic for an electrode on a sink plate) left an entropy-rate noise of ~eps·G·T above Newton's absolute tolerance, and the solve stalled at the true solution. Newton steps are unchanged by row scaling.
8. **Test scope:**
   - Kelvin's first relation is tested for constant α only (Thomson heat zero). A temperature-dependent α(T) is not supported yet.
   - The couple's open-circuit voltage is checked to 1e-5, because the electrodes sit ~0.2 mK above T_c; the per-leg relation (278) is checked exactly using measured temperatures.
9. **Electrical force is φ.** The electrochemical potential across material boundaries (contact potentials, 04) is not modeled yet: all electrical connections so far are within one conductor material.
10. **Mass-transfer application test (membrane permeation, Sieverts) is deferred.** The mass conductor and the Knudsen barrier are delivered; the DESIGN §0.1 row stays partly open.

### D-28 Plotting references complete; compaction and normalization solver chosen — accepted (user supplied the references, 2026-10-14); details proposed
- **References:** all the layout paper's references that the implementation needs are now in `_reference/graphing/`, catalogued in DESIGN §17.1 with what each is used for:
  - Freivalds–Kikusts 2001 ([10]);
  - IPSep-CoLa 2006;
  - Dwyer–Marriott–Stuckey 2005 ([7]);
  - Biedl–Madden–Tollis 2000 ([3], the journal version);
  - Lengauer ([14], local only).
- **Chosen:** the eq. (3) compaction and the normalization phase both use IPSep-CoLa's `solve_QPSC`: gradient projection, with projection by block merging and splitting on Lagrange multipliers. It is the same method as Freivalds–Kikusts' component shift/merge, so it is what the layout paper's authors used.
  - Normalization follows Freivalds–Kikusts: segment coordinates, obstacle-graph ordering, node drift + 10 × size.
  - Port assignment follows Biedl et al. §3.3.
  - Non-overlap constraints come from the Dwyer–Marriott–Stuckey scan line.
- **Alternatives:**
  - a generic QP solver (OSQP and similar) — another dependency, no warm start, slower per the IPSep-CoLa §3.3 experiments;
  - Freivalds–Kikusts' max-flow positive-cut formulation — superseded by the spanning-tree cut in the same paper;
  - force-scan overlap removal — lower quality per Dwyer–Marriott–Stuckey.
- **Why:** these are published, specific algorithms, so DESIGN §17.1 no longer relies on "standard textbook methods" for any phase. The separation-constraint structure keeps the solver small and dependency-free.
- **Supersedes:** the open solver choice in D-23 and DESIGN §17.1 (options (a)/(b)).

### D-30 All copyrighted PDFs are local only — accepted (user, 2026-10-14); supersedes the paper part of D-25
- **Chosen:**
  - Every PDF in `_reference/graphing/` is git-ignored (`_reference/graphing/*.pdf`), papers as well as books.
  - The Freivalds–Glagoļevs paper, committed under D-25, is removed from the tree (`git rm --cached`); the user keeps backups.
  - Its blob remains in history before this commit. History was not rewritten; there is no remote yet.
- **Still tracked:** `_reference/veinik_source/README.md` and `notes_by_section_ch6-18.md` (notes, not copies of works).
- **Consequence:** in a fresh clone, DESIGN §17.1 names every file exactly, but the PDFs must be copied in from the author's backups.
- **Alternatives:** commit the papers (D-25 as written); rewrite history to purge the earlier blob. That is possible if the repository is ever published; ask first.

### D-29 Slice 5 implementation: deviations from DESIGN.md — accepted (slice 5, 2026-10-14)
1. **Transformers are index 2, not index 1** (DESIGN §6.2 corrected).
   - Between two storing bodies the constraint P_a = r P_b involves only states, so λ needs one differentiation to appear.
   - Built in index-reduced form: λ's algebraic row is d/dt(P_a − r P_b) + (P_a − r P_b)/τ, with Hessian rows by AD and Baumgarte τ = 1e-4 s by default. It is formed in a second assembly pass, after all fluxes.
   - Power λ(P_a − r P_b) is zero on the manifold. The off-manifold mismatch is held at integrator tolerance and reported by `constraint_residual`.
   - Initial charges must satisfy the constraint.
   - Alternatives:
     - gyrator-form fluxes F_a = μ P_b, F_b = −μ P_a (power-exact but singular at P_b = 0, e.g. a gear at rest);
     - a stiff compliance (not ideal);
     - index reduction by a symbolic layer (excluded, D-2);
     - explicit elimination of λ inside the right-hand side (equivalent, but bypasses the DAE path the junctions need anyway).
2. **Junctions are index-1 algebraic nodes** (as designed). Their charge-target rows accumulate net inflow, so Kirchhoff's law needs no special code. DAE systems use a **dense Jacobian** (no sparsity prototype yet).
3. **Ideal converters `Converter(W, domains)`:** an exactly antisymmetric W across n ports on distinct nodes; thermal ports are forbidden.
   - Convenience constructors `Gyrator(r; a, b)` and `InductorBranch()`.
   - Inductors are bodies (`LinearInductor`, flux linkage Λ, force i), joined to circuits by the 3-port inductor-branch converter: dq_a = −i, dq_b = i, dΛ = φ_a − φ_b.
   - A DC motor is a 4-port converter (winding plus rotor).
   - Converters do not conserve their domains' totals, so `closure` skips converted domains, and `steady_state` treats converters as anchors.
4. **`Momentum` is one-axis.** The vector `Momentum{N}` of DESIGN §1 is deferred until a model needs it.
5. **Switching and time-varying reservoirs were pulled forward from "later"**, because E1/E2 need them.
   - `Switched(element; on = t -> Bool)` gates an element's whole application (Boolean only; transformers cannot be switched).
   - Reservoir forces may be functions of t: `domain => t -> P(t)`.
   - Discontinuity times are passed as `tstops`.
   - Accessor `transfer` on a switched element reports the inner element even while it is off.
6. **`assemble!` returns NaN for states outside the materials' domains.** Before, a trial stage (e.g. V < 0 inside a rejected Rosenbrock stage) threw a `DomainError` before `isoutofdomain` could reject the step.
7. **`closure` relative measures have floors of 1e-3 × |U₀| and 1e-3 × max(|S₀|, 1).** When almost nothing moves (a frictionless gas spring), round-off divided by round-off was meaningless; round-off in these sums scales with the magnitudes themselves.
8. **"Exactly isentropic" checks are to the ulp, not bit-exact.** A charge with an identically zero rate can still move by ~1 ulp per run, because the dense LU of each implicit stage mixes rows. Tests use 1e-12 to 1e-13 relative.
9. **E1 is shown by convergence**: the Carnot deficit and Ṡ_gen both fall about tenfold for a tenfold G. A strictly reversible bridge needs G → ∞. E2 (internal temperatures, (490)) and the Gouy–Stodola balance are checked exactly on the periodic cycle.
10. **Left limits at discontinuities.**
    - With `tstops`, time-varying reservoir forces and `Switched` gates are evaluated as left limits up to the next stop: the right-hand side uses min(t, prevfloat(next stop)).
    - A discrete callback moves to the next interval once the step landing on the stop is accepted.
    - Before this, the final Rosenbrock stage of a step ending on a stop saw the next interval's drive speed and gate state.
11. **Ledger error control: an open question for the user (it touches D-10).**
    - D-10 excludes ledgers from the integrator's error norm. When the physical states are easy to integrate but a ledger integrand is not, the ledgers become uncontrolled quadratures.
      - Example: a prescribed drive makes V linear in t while the adiabats keep S fixed, so each 100 s leg is one step.
      - Measured on the Carnot cycle: W off by 2–5 % (closure 2e-3 to 5e-3) without control; with control, W = Q_h − Q_c to 1e-8 (closure 4e-9).
    - Built: `simulate(…; ledger_control = true)` adds the ledgers to the error norm (they still never feed back, and steady solves still exclude them). The default stays `false` per D-10.
    - **Resolved by D-31:** `true` is now the default, with scaled ledger tolerances.
    - **Recommendation (adopted):** make `true` the default. Correctness and loss attribution matter more than speed, and `closure` exposes the problem but does not prevent it.
    - Alternative: post-hoc high-order quadrature of the ledgers from dense output.

### D-31 Ledgers are error-controlled by default, with scaled tolerances — accepted (user, 2026-10-14); supersedes the error-norm part of D-10 and resolves D-29 §11
- **Chosen:**
  - `simulate(…; ledger_control = true)` by default: ledgers are in the integrator's error norm. They still never feed back into the dynamics, and steady solves still exclude them.
  - Ledger absolute tolerances are scaled like the states they account for, per ledger kind recorded at compile time (`ledger_kind`):
    - entropy-type: the median body entropy tolerance;
    - energy-type: that × the median body temperature (the reservoirs' temperatures if there are no thermal bodies);
    - charge-type: the median body tolerance of that domain.
  - `atol_ledger = x` overrides. Guidance is in DESIGN §12.5 and the `simulate` docstring.
- **Measured** (on vs off, after warm-up; on uses the scaled tolerances):

  | Model | Steps (off → on) | Time | Energy closure (off → on) |
  |---|---|---|---|
  | slab, 300 cells | 143 → 143 | same | 9e-11 → 9e-11 |
  | rod, 25 cells | 157 → 162 | +1 % | 4e-12 → 3e-12 |
  | gas network | 205 → 210 | same | 6e-11 → 2e-11 |
  | heat exchanger | 212 → 221 | +9 % | 2e-14 → 2e-14 |
  | RLC (DAE) | 411 → 382 | −4 % | 4e-9 → 5e-9 |
  | Carnot drive | 254 → 568 | 2.3× | **2.4e-3 → 4e-9** |

  With a flat ledger tolerance of 1e-12 instead, the slab took 3.3× and the rod 1.4× for no accuracy gain.
- **Alternatives:**
  - keep ledgers out of the norm (D-10): fastest, but loss figures can be wrong by percent with no warning other than `closure`;
  - post-hoc high-order quadrature of ledgers from dense output: accurate, but every element's ledger integrand must be re-evaluated, and it does not compose with event handling;
  - a flat ledger tolerance: simple, but wastes steps while ledgers are near zero.
- **Why:** loss attribution is a consulting deliverable (CLAUDE.md). Correctness should be the default, and the measured overhead is small except where it buys accuracy.

### D-32 Mechanical and rotational systems: elastic forms of motion, spring branches, sources, gravity, monotone friction — accepted (user, 2026-10-15, with D-34; DESIGN §18; built in slice 6, D-35)
- **Asked for by the user:** basic mechanical systems (mass–spring–damper, pendulum) and rotational systems (gear trains, rotational dynamics) as intended uses (DESIGN §0.1).
- **Chosen:**
  - **Domains:** `Displacement` (x, F) and `Angle` (θ, M), from Veinik (5), (6), (12), (13).
  - **Materials:** `LinearSpring`, `TorsionSpring` and `PendulumPotential`, each a one-DOF accumulator.
  - **Converters:** the 3-port `SpringBranch` / `TorsionBranch`, mechanical twins of `InductorBranch`. Conversion between kinetic and elastic forms happens only between bodies.
  - **Friction:** viscous friction is the existing momentum conductor (kernel B). New monotone closed-form laws, `QuadraticDrag` and `CoulombFriction` (tanh-regularized), have D = ΔP f(ΔP) ≥ 0 by type.
  - **Gears:** the existing `Transformer`; clutches and brakes are a `Switched` conductor.
  - **Prescribed forces and torques:** `Source(domain => rate(t))`, a `SourceTransfer` with its work ledgered.
  - **Gravity:** `GravityField(g)`, a momentum source −Mg on each moving body, from the body's own inertial mass, with its work and ∫ w dt ledgered. The pendulum keeps gravity in its material instead.
  - **Tolerance scaling** needs a material-supplied charge scale where A_kk can be ≤ 0 (`PendulumPotential`).
- **Alternatives:**
  - (a) Springs as `ElasticVolume` through a unit-area gyrator. This works today, but puts newtons into the −p slot and metres into m³, which is misleading, and it is non-negotiable 3's volumetric force.
  - (b) Height as a charge with U ∋ Mgh (Veinik (11)). Potential energy would stay in the model and μ + gz would come automatically. It needs a converter between two domains of one body (ports must be on distinct nodes, D-29 §3), and it pays off only for hydrostatic heads. Revisit with balance-of-plant elevation.
  - (c) Gravity as a constant-force `Displacement` reservoir per body. Wrong once the mass varies (rockets, §19).
  - (d) A pendulum as one body with an internal converter. That breaks the rule that a body with no coupled part is a pure accumulator (DESIGN §2.2, Veinik §11).
  - (e) Coulomb friction with true sticking (switched or complementarity). Later.
  - (f) Euler's 3-D rigid-body equations. Needs vector angular momentum and a state-dependent antisymmetric operator inside one body, a new element kind (DESIGN §18.7). Not proposed.
- **Why:** each piece is a new (E, P) pair, a material, a converter with constant antisymmetric W, a conductance type that is sign-safe by type, or a ledgered source. The assembler, the kernels and the solver are unchanged. That is the architecture's own test (DESIGN §0.1).
- **To review:** the domain names; `GravityField` as a source rather than (b); whether mechanics goes before reactions (it is independent). *(Answered in D-34: the source, flagged; mechanics first, as slice 6.)*

### D-33 Rockets: feasible with moving variable-mass bodies, momentum-carrying streams and an ideal nozzle — accepted (user, 2026-10-15, D-47; DESIGN §19, slice 8)
- **Asked for by the user:** the rocket equation with a simplified engine (combustion, then expansion through a nozzle), with confirmation of whether the framework allows it.
- **Finding:** it does, after the reactions slice (now slice 7: reactions, plus an ideal-gas mixture with formation enthalpies) and the mechanics slice (now slice 6: gravity, drag). No new solver and no new `Transfer` kind.
- **Chosen:**
  - **`Moving(material; m_dry)`:** U = U_int + p²/(2M), with M from the body's own mass charges. This gives μ_k = μ_k,int − w²/2.
  - **Rigid links:** every body on the vehicle is `Moving`, held at one velocity by `Transformer(1)` on `Momentum` (index 2, D-29).
  - **Streams and pumps** carry momentum ṁ w_u between moving bodies. The generation gains ṁΔw²/(2T_d) ≥ 0.
  - **`Nozzle(A_t, ε; η)`:** a stream-kind outlet with choked mass flow, frozen isentropic exit state, v_e = √(2(h_c − h_e)), and a momentum bundle ṁ(w − v_e) − p_e A_e (stream thrust). Energy closes exactly per evaluation, the flow work at the moving exit plane included. Its generation is ṁ(s_e − s_c).
  - **Ambient pressure** uses the existing piston–atmosphere gyrator with area A_e.
  - **Test stand:** a w = 0 `Momentum` reservoir, whose momentum ledger rate is the thrust.
- **Alternatives:**
  - (a) One vehicle momentum with an inertia summing other bodies' masses. Non-local U (DESIGN §2.1).
  - (b) Stationary tanks feeding a moving chamber. Books a spurious ½ṁw² as dissipation.
  - (c) Thrust as a prescribed force source. No coupling to the engine, no energy closure, no loss attribution.
  - (d) Resolving the nozzle as a chain of bodies with inertance streams. Needs inertance (deferred), compressible supersonic lumped flow and shocks; not worth it for a simplified engine.
  - (e) A converter from the chemical form straight to momentum. Hides the thermal stage, so no I_sp physics and no loss map.
- **Why:** thrust then *follows from* momentum conservation in Veinik's own terms: the jet is an advected charge bundle, and the vehicle gains what the jet takes away. Energy closure and the loss map come for free, attributed to injection, combustion and the nozzle.
- **Known limits** (DESIGN §19.6):
  - a lumped chamber and a quasi-1-D nozzle;
  - constant c_p until the §2.4 Newton inversion exists (c_p(T), dissociation species);
  - a frozen nozzle, and an equilibrium chamber rather than combustion kinetics;
  - gaseous propellants until the real-fluid backend;
  - one axis.
- **To review:** whether slice 6 should support a sign-preserving nonlinear (mass-action) rate law as well as the linear affinity law. A rocket chamber does not need it; ignition and kinetics-limited combustion do. *(Answered in D-34: yes, in the reactions slice; future work only if it presents an issue.)*

### D-34 Review of revision 3: gravity as a source, special-purpose elements over rule changes, slice order — accepted (user, 2026-10-15)
- **Working principle (user):** do not sacrifice generality for a specific case. A few explicit, ledgered, special-purpose ("virtual", one-off) components that enable a specific case are preferred to changing the framework's rules to accommodate it. This is recorded in CLAUDE.md "How to work".
- **Gravity: `GravityField(g)`, a momentum source** (D-32 as proposed). This decision is **flagged: it may need to change**.
  - **Evaluated against the alternative,** height h as a charge with U ∋ Mgh (Veinik (11)).
    - The source fits the current framework more closely. It is an explicit external source with ledgered work, an element kind non-negotiable 1 already provides, and it needs no rule change.
    - Height as a charge is physically the same as dQ = F dx. But for a body of variable mass it must sit in the same body as the mass (M is that body's charge), and dh/dt = w then needs a converter with both ports on one node. That changes a framework rule (D-29 §3), and it puts a conversion inside a body whose p–h part is uncoupled in U, against the accumulator principle (DESIGN §2.2, Veinik §11).
    - For a constant mass, height-as-charge works today without any rule change: a separate weight body U = mgh joined by a converter. `PendulumPotential` is that pattern.
  - **Revisit when:**
    - hydrostatic or elevation heads matter in flow networks (balance-of-plant loops), where μ + gz is needed in streams;
    - a model needs altitude-dependent gravity;
    - the vehicle's height must be a state rather than a derived ledger.
    - Then weigh relaxing the converter rule (ports on one node with distinct domains; W stays antisymmetric) against a dedicated height element.
- **Mass-action rate law:** the reactions slice includes a sign-preserving nonlinear (mass-action) rate law alongside the linear affinity law. Only if it presents an issue does it move to future work, with the issue recorded.
- **Slice order:** mechanics (DESIGN §18) is built first, as slice 6, and merged after the user's approval before the reactions work starts. Reactions and nuclear sources become slice 7, rockets slice 8.
- **`docs/veinik/08_validation_cases.md`:** the user approved adding the mechanical validation cases there (new §G). The changed lines are reported with the slice.

### D-35 Slice 6 implementation: deviations from DESIGN §18 — accepted (slice 6, 2026-10-15)
1. **`stiffness_floor(mat)`** is the tolerance hook DESIGN §18.2 asked for.
   - Per charge, default 0; `PendulumPotential` uses mgℓ; a `Composite` concatenates its parts' floors.
   - `abstol_vector` and the steady-state row scaling (`_charge_scales`) use max(|A_kk|, floor).
   - Test: at θ₀ = π/2, where A ≈ 0, the angle's tolerance is 1e-8/(mgℓ).
2. **`inertial_mass(mat)`** is what `GravityField` reads: `InertialMass` m, or the sum over a `Composite`'s parts; `nothing` otherwise, which is a compile error. The rate is constant (−Mg) because every mass so far is a parameter. `Moving` (slice 8) will make it state-dependent.
3. **`GravityField` ledgers its work only.** The displacement ledger ∫ w dt is deferred to slice 8. For a constant mass the height change is −W/(Mg).
4. **`Source` is non-thermal only**, and acts on bodies only.
   - A thermal (heat) source is a separate element, not built yet; `Source(Thermal() => …)` is refused with that message.
   - Its ledgers are the charge received (negative: it supplies) and the external work supplied. The work sits with the pump work (`pump_work` lists it), and `source_ledger` reads both.
   - A source anchors its charge and energy for `steady_state`, like a boundary.
5. **`QuadraticDrag(c; domain)`, `CoulombFriction(F_c; w_ε, domain)`:** the domain defaults to `Momentum`. Coulomb friction is tanh-regularized, with no true sticking.
6. **Unchanged framework rules:** converter ports stay on distinct nodes. Gears, rack and pinion, clutches and brakes use the existing `Transformer` and `Switched`. No lossy transformer (load-dependent mesh efficiency) is built.
7. **Test layout:**
   - 08 §G G1, G2, G5 (translational) and G6–G8 (rotational) are application tests in `test_applications.jl`.
   - G3, G4, G9, construction and validation, the kernel identities, the source ledgers, a random mechanical-chain A1/A3 test, the pendulum stiffness floor, rack and pinion, and a DC motor with a gearbox are in `test_mechanics.jl`.
8. **`docs/veinik/08_validation_cases.md`** (user-approved, D-34): lines 60–75 were added, a blank line, then §G with the G1–G9 table. No existing line changed.

### D-36 Veinik section notes are local only; history purged before the first remote — accepted (user, 2026-10-15); supersedes the notes part of D-25 and the "still tracked" line of D-30
- **What happened:** the first push to the private GitHub remote (`itme-yaboi/ThermoKinetics`) carried `_reference/veinik_source/notes_by_section_ch6-18.md` (tracked since D-25) and, in older commits, the Freivalds–Glagoļevs PDF (untracked by D-30 but left in history). D-30 had asked for a purge before any publication; it was missed. The user deleted the GitHub repository.
- **Chosen:**
  - `notes_by_section_ch6-18.md` is git-ignored and local only, like the OCR texts. Only `_reference/veinik_source/README.md` stays tracked under `_reference/veinik_source/`.
  - History was rewritten (`git filter-branch --index-filter … --prune-empty`, all branches) to remove both files from every commit; the old refs, reflog and objects were purged. All commit hashes changed. A backup bundle of the pre-rewrite repository was taken first.
  - Before any push, the full list of tracked files and every path in history is shown to the user.
- **Alternatives:** keep the notes tracked (D-25 as written); delete them in a new commit only (leaves them in history); force-push without deleting the remote (old commits can stay reachable on GitHub).
- **Consequence:** in a fresh clone, the README's notes row and CLAUDE.md's pointer resolve only after the author copies the notes in locally. `docs/veinik/` is unaffected and remains the working reference.

### D-37 `docs/veinik/` is local only — accepted (user, 2026-10-15); extends D-36
- **Chosen:** `docs/veinik/` (the curated summaries of Veinik's books) is git-ignored and removed from the index (`git rm -r --cached`); the files stay on disk and remain the working reference CLAUDE.md points to. Purging it from earlier commits is a separate history rewrite, done only with the user's go-ahead.
- **Alternatives:** keep it tracked in a private repository (the summaries are close paraphrase of copyrighted books).
- **Consequence:** a fresh clone has no `docs/veinik/`; the author copies it in from a local backup. No code or test reads it.

### D-38 Package renamed `VeinikThermo` → `ThermoKinetics` — accepted (user, 2026-10-15)
- **Chosen:** the package, module and main file are `ThermoKinetics` (`Project.toml` name, `src/ThermoKinetics.jl`, `using ThermoKinetics`), matching the private GitHub repository `itme-yaboi/ThermoKinetics`. The UUID is unchanged. The local working directory keeps its old name.
- **Alternatives:** keep `VeinikThermo.jl`; `GeneralizedThermo.jl`, `ThermoNetworks.jl`, `EntropyLedger.jl`, `ChargeForce.jl`.
- **Why:** the user chose the repository name; Julia expects the repository `Name.jl` and the package `Name` to agree.

### D-39 Slice 7 scope: mixture, reactions, nuclear and heat sources — accepted (user, 2026-10-15, with D-40; DESIGN §8.1–8.2)
Detail of D-11 for the slice. No new element kind: reactions are `InternalTransfer`, nuclear and heat sources are `SourceTransfer` (DESIGN §4.2).
1. **Constant c_v per species, closed-form U.** Alternative: NASA polynomials now, which need the Newton inversion of §2.4. Deferred to the property backends; costs high-temperature accuracy in slice 8.
2. **Molar masses derived from atoms** through a small atomic-mass table, so reactions conserve mass to round-off. Alternative: user-given M with a tolerance check, which leaks mass at data precision (non-negotiable 1).
3. **Every species stays present (m_i > 0).** Alternative: a regularized zero-mass species (μ → −∞ makes it delicate). Future work.
4. **Mass action with the reverse rate implied:** ξ̇ = V k(T) Π c^order · (−expm1(−A/R̄T)). A·ξ̇ ≥ 0 holds exactly in floating point. Alternative: user-given forward and reverse constants (can contradict the thermodynamics); r_f − r_b computed separately (sign can fail near equilibrium in floating point).
5. **Nuclear energy convention: dU = E_dep ξ̇ exactly.** The species change and the deposited heat enter the entropy as (Q̇ + A ξ̇)/T, booked as boundary entropy, never as production. Alternative: species exchanged with the outside at the body's μ, like `Source` (energy Q̇ + Σμ dm external). Rejected: the reported nuclear energy would depend on the chemical datum.
6. **`HeatSource` is a separate element.** Alternative: `Source(Thermal() => Q̇)`. Rejected: `Source` prescribes a charge rate, and a heat source prescribes an energy rate.
7. **Validation cases H1–H8** go into `docs/veinik/08_validation_cases.md` §H only with the user's approval (D-34 precedent); until then they live in DESIGN §8.2.
- **To review:** items 1–7, and whether the tritium storage bed is the right nuclear application test.

### D-40 Review of the slice 7 scope — accepted (user, 2026-10-15)
1. **Constant c_v:** accepted for now.
2. **Trace amounts:** accepted, on the condition that a species approaching zero never makes a run fail. Made a requirement with its own test (H10).
3. **Nuclear energy convention (dU = E_dep ξ̇):** accepted.
4. **08 §H:** the user approved adding H1–H10 to `docs/veinik/08_validation_cases.md`.
5. **Application tests:** the tritium storage bed, and, at the user's request, tritium generation in irradiated water (H9).
   - H9 needs a liquid mixture, so slice 7 adds `IdealSolution`: a rigid condensed mixture with ideal mixing, closed-form U(S, m), no new element kind. Alternative: steam as an `IdealGasMixture` (wrong phase for a coolant).
   - Atom symbols carry a nucleon number for the nuclear checks; element symbols use their main isotope's, isotopes their own.

### D-41 Slice 7 implementation: deviations from DESIGN §8.1 — accepted (user, 2026-10-15, with D-42 and D-43) (slice 7); supersedes "every species stays present" in D-39 §3 and D-40 §2
1. **Trace floor in the mixing entropy** (the user's condition in D-40 §2: approaching zero must not fail).
   - **Problem, measured:** with exact ideal mixing (n ln n, μ → −∞ as n → 0), a burnt-out reactant whose equilibrium amount is far below round-off sits in round-off noise (≈ 2 eps of the body's mass) at the edge of the valid states. In the H4 calorimeter, 40 % of steps came out slightly negative and were rejected, at any step size, and the run hit MaxIters.
   - **Chosen:** below δ = `trace` × the body's total mass (default 1e-13, a material keyword), n ln n is continued by its second-order Taylor expansion at δ. U stays a convex, C², first-order homogeneous potential, μ is finite for every amount, and mass action uses the matching activity amounts, so A·ξ̇ ≥ 0 and the rate still vanishes at A = 0. The same calorimeter now runs 20 s in 276 steps.
   - **Consequences:** species may start at zero. A depleted species settles at a small finite value that can read slightly negative (about −1e-12 of the body's mass in H3/H4). Species above the floor are exactly ideal. A trace species of real interest below 1e-13 of its body's mass needs a smaller `trace`.
   - **Alternatives tried:** exact log with m > 0 (the stall above); tighter species tolerances (170,000 steps, more rejections: the noise is round-off, not tolerance); a log shift m + δ (moves the singularity, same stall); Rodas4P, whose stability function stays positive on decay (worse here). Not tried: log-mass states (a new state representation, a design change).
   - **Found on the way:** Rodas5P's stability function is negative for hλ ≳ 6 on linear decay. That is harmless with the floor, but worth knowing for any charge that relaxes to a bound.
2. **Closure scales.** Energy is measured against the energy moved, now including what is converted inside bodies (Σ|P·ΔE|, e.g. chemical to thermal); before, an adiabatic reactor had only the 1e-3·|U₀| floor and read integrator-level error (2e-9 of U) as 2e-6. A species is measured against the total mass of the bodies holding it; before, round-off from 1 kg of water read as 4e-9 drift on 1e-10 kg of HTO. Both make closure no stricter for existing tests.
3. **Nucleon convention:** `compile` checks Σ ν_i A_i = external.nucleons (what is absorbed from outside), not "+ … = 0" as first written in §8.1.
4. **Ledgers and accessors:** heat and nuclear sources are read with `boundary_ledger` (received: −U is the energy supplied). New `reaction_rates` and `element_totals`. `boundary_ledger` with no mass slots no longer fails on an empty sum.
5. **Steady state:** a reaction does not anchor its species (closed reacting groups conserve element totals, the unimplemented constraint-row case, so `steady_state` refuses them). A nuclear source anchors its species, like a `Source`; a heat source anchors nothing.
6. **Tests:** H10 runs H6 for 60 half-lives (not 20). At that depth the tritium amount is checked only for |n_T| ≤ 1e-10 n₀ (it is below the solver's absolute tolerance). 08 §H says 60 (user's file, edited with the §H approval).
- **To review:** item 1 (the floor and its default) and item 2. *(Answered: item 1 in D-42; item 2 in D-43; the rest accepted.)*

### D-42 The trace floor: when it matters and how to handle it — accepted (user, 2026-10-15); accepts D-41 §1
- **Chosen:** the default `trace = 1e-13` stays.
- **What it changes:** only species below δ = 1e-13 of their body's total mass. Above δ the mixture is exactly ideal, with no error at all. Below δ, the species' μ (its activity) is no longer ideal, and its mass can read slightly negative (about −1e-12 of the body's mass).
- **What it never changes:** conservation of mass, elements, energy and nucleons (all exact); the sign of entropy production; first-order nuclear rates, which use the true amount, not the activity.
- **When it matters:** when a species below 1e-13 of its body's mass has a driving force you care about:
  - its permeation or diffusion flux (a mass conductor runs on μ);
  - a reaction rate that depends on it (mass action uses the activity);
  - its equilibrium partial pressure or concentration.

  It does not matter for burnt-out reactants, for daughters early in a decay, or for inventories tracked only through nuclear sources.
- **How to handle it:**
  1. Lower `trace` on that body's material, but not below about 1e-15. The floor must stay above round-off noise (≈ 2 eps of the body's mass); below that, the stall of D-41 returns.
  2. If the species is below even that, model it in its own body (a separate trace inventory, or the stream or volume where it is not swamped), so it is not a 1e-15 fraction of a large body.
  3. When reporting, read a slightly negative amount as zero.
- **Revisit when** a fusion case needs accurate μ for a species below 1e-15 of its body (e.g. tritium in bulk coolant at those levels): that would need log-amount states, a design change.

### D-43 Closure scales accepted, with a per-species absolute report — accepted (user, 2026-10-15); accepts D-41 §2
- **Chosen:** D-41 §2 stands (energy against the energy moved including conversion inside bodies; a species against the mass of the bodies holding it). Because the second scale dilutes a real loss of a trace species, `closure` also returns:
  - `species`: per conserved species, `(; abs, rel)`: the absolute drift [kg], and the drift relative to that species' own largest body amount;
  - `elements`: when reactions act, per element `(; abs, rel)` [mol].
- **Why:** `charge_rel` stays a scale-aware summary for pass/fail, and an exact trace inventory (tritium accounting) is checked against itself, not against the body.
- **Alternatives:** revert to per-species scales (trace species read round-off as drift); leave the check to the user via `total_charge` (easy to forget).
- **How to use:** for tritium accounting, check `closure(cs, sol).species[:HTO].abs` against an absolute limit (it sits at about eps × the body's mass, e.g. ≤ 1e-15 kg in 1 kg of water).

### D-44 Low-mass-fraction species in reactions and nuclear processes: a modelling note — accepted (user, 2026-10-15); extends D-42
- **The note:** when a chemical reaction or a nuclear process involves a species whose transport or rate you care about at a very low mass fraction of its body, check it against the trace floor (1e-13 of the body's mass, D-41, D-42) before trusting the result.
  - **Exact regardless of the floor:** energy, entropy and every conservation check, and the heat carried with the species (μ_b dm_b stays consistent; its effect is ~1e-11 J for a kilogram-scale body).
  - **Inexact below the floor:** μ_b itself (100× below the floor, about 4.6 RT too high), so any rate driven by differences in μ_b: permeation and diffusion between bodies (about 25× too small for a 1e-15 vs 1e-14 pair), partition equilibria, and mass-action rates involving it. First-order nuclear rates are not affected (they use the true amount).
- **Where tritium lands (HTO in 1 kg of water):** 1 Ci ≈ 1e-7 of the mass (exact); 1 µCi ≈ 1e-13 (at the floor); drinking-water and discharge levels, ~10 kBq/L, ≈ 2e-13 (at the floor). Plant inventories, coolant activity, storage beds and breeding blankets are well above it.
- **What to do:** for discharge-level or environmental tritium, or any trace species whose transport matters near or below 1e-13, lower `trace` on that body to ~1e-15 (not lower: round-off), or model the dilute species in its own smaller body. Read slightly negative amounts as zero.
- **Alternatives:** lower the default now (cost: robustness margin over round-off shrinks for every model); leave it implicit in D-42 (easy to miss).

### D-45 HTML documentation with Documenter.jl, in `docs/` — accepted (user, 2026-10-15; layout left to Claude); slice 7a
- **Chosen:** the standard Julia layout inside the existing `docs/` folder: `docs/Project.toml` (Documenter plus the package, its own environment, so the package gains no dependency), `docs/make.jl`, pages in `docs/src/`. Build: `julia --project=docs docs/make.jl`; open `docs/build/index.html`. The build is strict (no `warnonly`), runs the getting-started examples, and fails if an export is missing from the reference (`checkdocs = :exports`).
- **Project records:** DESIGN, STATUS and DECISIONS are copied into `docs/src/project/` at build time (git-ignored), with local-file links made plain text. `docs/veinik/` and `_reference/` are never read or included.
- **Not deployed:** the repository is private, so there is no `deploydocs` and no `gh-pages` branch; `docs/build/` is git-ignored.
- **Alternatives:** a separate `human_docs/` folder (non-standard, tooling expects `docs/`); a docs branch (only useful for GitHub Pages deployment); hand-written HTML (no docstring extraction, drifts from the code).
- **Built with:** the Documenter setup, API pages and examples by a Sonnet subagent, reviewed and corrected (an example compared a production rate with an integrated value; misfiled API entries); the architecture page and the missing docstrings directly.

### D-46 The author's textbooks for examples and validation are local only — accepted (user, 2026-10-15)
- **Chosen:** `_reference/verification/` (textbooks the author has read) is git-ignored as a whole. They are used to set up example problems and validation cases, cited by book, section and problem number; their text, figures and tables are never copied into the repository. They are not a source for the framework's physics or code (`docs/veinik/` and DESIGN are). The list and the use of each book are in CLAUDE.md.
- **Alternatives:** track them (copyrighted, D-25/D-30/D-36 precedent); ignore per file type (the folder holds `.pdf` and `.PDF`; a folder rule also covers future additions).

### D-47 Review of DESIGN §19 and D-33 (rockets) — accepted (user, 2026-10-15)
- **Chosen:** DESIGN §19 and D-33 are accepted as written; slice 8 builds them.
- **Clarification (slice numbering):** D-33's "To review" line asks whether "slice 6" should support a mass-action rate law. When D-33 was written, the reactions slice was planned as slice 6; D-34 reordered the slices (mechanics became 6, reactions 7) and answered yes. Slice 7 built both `AffinityLaw` and `MassAction` (D-39–D-41). Slice 6 (mechanics) has no rate laws.

### D-48 Slice 8 implementation: deviations from DESIGN §19 — accepted (user, 2026-10-15; items 1 and 3 kept as future work) (slice 8)
1. **Loss attribution in an equilibrium chamber.** With a fast reaction the chamber sits at A ≈ 0, so the reaction element produces almost nothing (≈ 5e-10 W/K in the static fire) and the combustion irreversibility (≈ 2300 W/K) is booked on the injectors: their generation includes the jump from the feed's μ to the chamber's equilibrium μ. Total production still equals the jet's entropy minus the feeds', exactly. This is correct for a well-stirred lumped reactor; DESIGN §19.5 expected "combustion dominant", which holds only for a finite-rate or two-zone chamber (a frozen mixing zone feeding a reaction zone). Corrected in §19.5.
2. **Species datum consistency.** An `IdealGas` has its own datum (h = c_p T, s = 0 at the reference state), different from a species' formation data, so an `IdealGas(:H2)` tank feeding a mixture holding `Species(:H2)` would exchange inconsistent μ. `compile` now refuses a species held both by a material with species data and by one without. Rocket tanks are one-species `IdealGasMixture`s, and `isentropic_state` was added for them (closed form).
3. **Feeds into a richer mixture.** Pumps and flow boundaries may deliver a species set that is a subset of the receiving body's; streams still need equal sets (a stream can reverse, and would carry species the upstream body lacks). So the engine is pump- or boundary-fed; a pressure-fed (stream) feed needs a check valve, not built. A flow boundary accepts a one-species mixture as its fluid.
4. **The nozzle is always treated as choked.** No back-pressure check (the design allowed one); ambient pressure on the exit area is the existing gyrator, and a test stand reports the vacuum thrust ṁ v_e + p_e A_e.
5. **η < 1:** the exit pressure stays the isentropic p_e; T_e = T_es + (1 − η)(T_c − T_es), and the generation is ṁ c_p ln(T_e/T_es), a form that is ≥ 0 exactly (the isentropic part of s_e − s_c cancels analytically; computing it numerically could give −eps).
6. **Gravity** is a separate compiled kind (`CompiledGravity`) with the current mass of a `Moving` body (m_dry + fluid), and a third ledger, ∫ w dt (`source_ledger(...).x`), for every gravity field.
7. **`consistent_state`** now scales its tolerance by the algebraic residuals' size: a light chamber pushed by its thrust has transformer residuals ~1e5 m/s², above which a fixed 1e-12 sits in round-off (Newton stalled).
8. **The flight test is exact.** With adiabatic tanks and ideal pumps the delivered state depends only on the tanks' constant specific entropy, so the chamber stays steady and c_eff is constant: Beer & Johnston Sample Problem 14.9 (velocity, and its integral for the altitude) replaces the design's Tsiolkovsky bound.
- **Built with:** the kernel tests by a Sonnet subagent (reviewed; its `@test a ≈ b && c ≈ d` lines were split, since that form does not parse as two comparisons); the physics, the application tests and the hand calculation directly.
- **To review:** items 1 and 3 (whether a two-zone chamber or a pressure-fed feed is wanted). *(User, 2026-10-15: accepted as built; both kept as future work, to be addressed later: a finite-rate or two-zone chamber so the combustion loss is booked on the reaction, and a check valve for pressure-fed feeds.)*

### D-49 Slice 9: cycles E3 and E4 — accepted (user, 2026-10-15: E3/E4, then the Clapeyron backend, then plotting)
1. **`EntropicCharge(domain; C, B, T_ref)`**, U(S, e) = C T_ref exp((S + B e²/2)/C), so T = T_ref exp((S + B e²/2)/C) and P = B e T: the simplest convex potential with Veinik's P = B e T (det ∇²U = B T²/C > 0). P has an absolute zero at e = 0, which the force-ratio efficiency needs (`05`). Alternatives: a gas body on the volumetric domain (its force −p has no B e T form at constant volume); a non-potential law (would break P = ∂U/∂E).
2. **E3 is built from existing elements:** one body, switched charge taps to two charge reservoirs (on the electrical form of motion), and thermal contact. For the reversible limit the heat reservoir's temperature is prescribed along the ideal path (T = P/(B e) on the charge legs, constant e on the others), so switching is force-matched and exchange tracks the body; production then falls like 1/G. With fixed-temperature reservoirs it stays finite. η = 1 − P″/P′ holds exactly at the reservoir ledgers in every case: it is the charge leg's force ratio, not a heat-engine efficiency (`05`).
3. **E4 (Otto)** is an application test (power cycles row): a gyrator piston driven by a smooth prescribed velocity, a `HeatSource` at constant volume, a switched rejection conductor. Validated against Çengel 9e Example 9-3 (cold-air standard, η = 0.565 for r = 8).
4. **Bug fixed:** `force_buffer` typed its buffer from the state only, so a reservoir force varying smoothly in t (the Otto piston, the E3 temperature path) failed when a Rosenbrock method differentiated f in t. The buffer now promotes over t's type. Earlier tests used only piecewise-constant inputs, which never carried the dual part.
- **Built with:** the E4 test by a Sonnet subagent (reviewed; it passed unchanged once the fix was in); the material, E3 and the fix directly.

### D-50 The real-fluid backend: split, packaging and fluids — split, extension and fluids accepted (user, 2026-10-15); slice 10a scope proposed (DESIGN §2.5)
- **Split (user):** slice 10a is the single-phase backend with Brayton cycles; slice 10b adds the two-phase flash for pure fluids and Rankine cycles. Rankine boils and condenses its working fluid, which the single-phase backend cannot represent (§2.4 item 7). Alternative: one slice (larger, no checkpoint between the inversion and the flash).
- **Packaging (user, after explanation):** Clapeyron is a weak dependency and the backend a package extension in `ext/`, in the same repository. Measured: Clapeyron 0.6.29 precompiles in ~30 s and loads in 0.76 s; the extension keeps the core free of its dependency stack (CSV, databases). Alternative: a normal dependency (simpler; every user pays the load).
- **Fluids (user):** helium, water and CO₂. Helium uses GERG-2008 (Clapeyron bundles no helium multiparameter equation; GERG-2008 includes helium, and its density at 300 K, 1 MPa is 1.5970 kg/m³); CO₂ Span–Wagner (`SingleFluid("carbon dioxide")`); water IAPWS-95 (s = 392.8 J/(kg K) at 300 K, 1 MPa, matching the steam tables' 0.3930). Air (GERG-2008 N₂/O₂/Ar) for the textbook Brayton check. CoolProp stays a later wrapper (D-7).
- **Probe findings:** all three models differentiate with ForwardDiff (gradients and Hessians of `eos(model, V, T, n)`). Clapeyron's bulk property functions take (p, T), not (V, T): inside the material only the V-based Helmholtz energy is used.
- **To review:** the slice 10a scope in DESIGN §2.5 (material interface, guards, datum rule, tests).

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