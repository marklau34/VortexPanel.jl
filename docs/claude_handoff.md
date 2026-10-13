# VortexPanel.jl: development notes for the next Claude

These notes are for a Claude instance continuing work on this repo. They summarize how the project was built (one long session, 12–13 Oct 2026): how the user works, the tooling, what went wrong, and what is still open.

The technical record lives elsewhere, and these notes don't repeat it:
- `README.md`: every module, convention and validation table, stage by stage.
- `docs/free_wake_guide.md`: the free-wake settings and their sensitivities.
- `Plan Galerkin Linear-Vortex Panel Solver (Rigid Wake Phase).md` (repo root): the user's brief, which the README calls "the design doc".

Read the brief and the README's Status section first.

## 1. Where things stand

**The project.** VortexPanel.jl is a Galerkin linear-vortex panel method:
- nodal doublets μ on flat triangles, with a constant vortex sheet γ = n × ∇μ on each;
- prescribed sources σ = −V∞·n, and the no-penetration condition tested against the hat functions;
- a constraint matrix C, which absorbs mesh pathologies and the nodal Kutta condition;
- a bordered dense LU.

**Reference codes** (sibling checkouts under `C:\Users\MarkLau\GitHub\`):
- **MESA.jl** (a Morino code) is the baseline. The brief says to use it unmodified. Its reference case is `examples/wing/run_wing_surfpan.jl`: a NACA 0015 rectangular wing, AR 3.88, at α = 5°.
- **FLOWPanel.jl** was consulted for trailing-edge and wake handling.
- **PAGE** (`_tools_release/PAGE`) and **MeshLibrary.jl** generate the geometry. Both are pinned through `[sources]`; the README's Environments section says why.

**Branch** `galerkin-rigid-wake`. Commits, in order:

| commit | content |
|---|---|
| 824c0e4 | geometry layer, mesh abuse, DOF map, MESA baseline (steps 1–2) |
| fa68d88 | source-triangle kernel (step 3) |
| 2554c3e | Galerkin assembly, solve, post-processing; sphere (stage 3) |
| 2e329df | ellipsoids; wide-corner and contact quadrature (stage 4) |
| 16eb51f | rigid wake, lifting assembly, force routes; wing (stage 5) |
| f5adf48 | comparison with MESA, 2-D section referee, zero-tip Kutta option (stage 6) |
| 20233a0 | zero-tip Kutta by default; test-weight comparison (stage 8) |
| b49ccda | isoparametric quads; element comparison (stage 7) |
| cdcc189 | free wake with vortex-particle conversion (phase 2) |
| 0ee0eb2 | automatic cores and conversion distance, row merging, ParaView example |
| a39a76d | the user's own commit: free-wake guide, `notes.md`, `.vscode/` |

**Done:**
- All seven modules of the brief and its eight validation stages.
- Phase 2: the brief's first two deferred items, the particle wake and unsteady flow.

The last full test run passed 11,181 tests.

**User-confirmed defaults.** Don't change these without asking:
- `AssemblyOptions(test=:hat)`.
- `DOFMapOptions(tip_kutta=:zero)`, kept "for now". It doubles the lift's sensitivity to quad-diagonal flips next to the trailing edge (see §5).
- Split triangles as the default element. The isoparametric quads (`element=:quad`) cost 9–19× the assembly.
- Free wake: `max_rows = :auto`, `core = :auto` (overlap 1.5) and `merge_rows = 1`. When offered merging, the user said "don't change defaults"; `merge_rows = :auto` stays opt-in.

## 2. Working with the user

- **Stage-by-stage cadence.** The user drives with short prompts: "do the next step", "commit, then do the next stage", "commit and do stage 6".
  - Finish a stage completely: code, tests, validation script, README.
  - Then report the results and propose the next stage. Don't start it unasked.
  - Big phases (step 1, phase 2) went through plan mode. `AskUserQuestion` was used only for real design choices.
- **Commits.**
  - Commit only when asked.
  - End every message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - Leave the user's own files (`notes.md`, `.vscode/`) for them to commit.
  - The git author email is for authorship only.
- **Defaults are the user's decisions.** Present the options with numbers and let them choose. Examples:
  - "make tip_kutta=:zero the default";
  - "keep the hat functions. keep zero tip for now";
  - "dont change defaults".
- **The brief is binding where it is specific; MESA is the model where the user says so.** For phase 2, the user chose:
  - MESA's time stepping: an implicit trailing-edge panel, free rows advected by forward Euler, `maxwakelen`-style row counts, no stretching, coupling through velocity with no matrix rebuild;
  - the brief's conversion: local to each wake triangle.
  - When I wrote "conversion done MESA's way", the user caught it at once. Say exactly which parts follow MESA and which follow the brief.
- **Explain physics, not just numbers.** The user knows panel methods and MESA well, and asks why the Galerkin-specific things happen. Questions came up like:
  - "why can't the body carry filaments";
  - "what is the wake's tip filament";
  - "does the loading always go to zero at the tip";
  - "what is Symmetric Galerkin".
- **Make percentages concrete.**
  - Always say what a percentage is relative to, and give the absolute CL values next to it. The user asked "what are the CL values that got this table" and "is this the lift vs how far the wake has advected?".
  - Name the lift route: pressure, Kutta–Joukowski or impulse.
- **Runnable examples.** The user runs things themselves and visualizes in ParaView.
  - Give the exact command, including the environment: `julia -t auto --project=validation ...`. The user once asked "which julia env should i run in?".
  - Write VTK/PVD output with viewing tips in the script header.
- **Notes.** The user saves findings to `notes.md`. It currently lists candidate next steps (§6).

## 3. Tooling on this machine (Windows 11, Julia 1.12.6, 16 threads)

**Julia through the MCP tool:**
- Use the temporary session, then `Pkg.activate("C:/Users/MarkLau/GitHub/VortexPanel.jl/test")`, or `.../validation`. Passing the test env as `env_path` fails with a TestEnv error.
- Revise doesn't track `src/`. Call `julia_restart` after source edits, then re-activate the environment.
- Full suite: `Pkg.activate(<repo root>); Pkg.test(; julia_args=`--threads=16`)`. It takes about 2 minutes.
- `Pkg` rewrites `Project.toml` files and strips TOML comments. Explanations of the `[sources]` pins live in the README.

**Long validation runs:**
- Durations: stage 6 ~25 min, stage 7 ~2 h, free wake ~35 min.
- Run them section by section in the session: `include` the script and call its section functions.
- Or run them in the background through Bash. Julia fully buffers stdout redirected to a file, so follow progress through the CSVs each section writes, or pipe through `cat`.
- Keep other CPU work off the machine meanwhile, since the scripts report wall-clock times.
- To load a script without overwriting its outputs, use a throwaway module: `@eval module X include("...") end`.

**Shell quirks:**
- There is no Python. Parse JSON with PowerShell (`ConvertFrom-Json`), or work in Julia.
- Bash heredocs containing apostrophes fail; append prose with the Edit tool.
- `grep -o '.\{0,N\}'` on the 30 MB session transcript hangs. Find the byte offset with `grep -obF`, then cut with `tail -c`/`head -c`.
- Headless Edge hangs when rasterizing SVGs. Check the plots structurally instead.

**Outputs** go to `validation/output/` (untracked). The CSVs there are the source of truth for the README's numbers. Some are stale relative to later REPL sweeps: `freewake/sensitivity.csv` still has the old core labels.

## 4. Validation practice that worked

- **Clean and abused meshes at every stage** (`MeshTools.abuse`), as the brief requires.
  - Stage 3 was first validated on clean icospheres only. Those meshes are centrally symmetric, which made λ vanish to rounding whatever the quadrature. A wrong "compatible by construction" claim went into the README, and a T-junction quadrature defect stayed hidden.
  - Prefer generic meshes and a generic V∞ direction.
- **Use diagnostics that can't cancel.** λ is one number and can be small by symmetry. The per-unknown column-sum defect of A·C is the real quadrature check.
- **Stay away from triangulation symmetry lines.** The default quad split mirrors its diagonals at y = 0. That dents the circulation there, by up to 0.6% on wide strips, exactly where stages 6 and 8 measured the midspan section (found in stage 7). Cross-check with `split=:d13` or the quads.
- **Separate the method's error from triangulation artifacts.**
  - Most of the default split's skew sensitivity was diagonal flips next to the trailing edge, not skew.
  - `GeometryOptions(split=:pattern)` freezes the clean mesh's diagonals and isolates the skew itself. At L2 that is 5× less sensitive than the default split.
  - The quads have no diagonals to flip.
- **Use independent referees.**
  - Exact flows: the sphere, and Lamb's ellipsoids via Carlson's R_D.
  - A 2-D linear-vortex panel method (`validation/common/panel2d.jl`) for sections.
  - Adaptive cubature for kernels and pair blocks.
  - The referee settled the MESA comparison. The codes bracket the exact section: ours converges at first order from above, MESA's at about order 0.5 from below.
- **Extrapolate carefully.**
  - With the zero tip, the wing's lift has an observed order near 0 under uniform refinement, because the chordwise and tip errors cancel. Refine chord and span separately.
  - Apply Richardson only where the observed order is at least 0.3.
  - Extrapolate 2-D sections quadratically in 1/AR.
- **Time marching approaches steady state algebraically,** like (span / wake length)². Compare against closed rigid wakes of the same rows, extrapolate in 1/L², and don't march to convergence.
- **Check every number in a document against the outputs before handing it over.** A final pass on the free-wake guide found four slips: a rounding, a range, an unlabelled reference, and a prescribed-wake result described as a free-wake one.

## 5. Numerics that bit, and what fixed them

**Quadrature:**
- *Graded rules saturated near 7e-5 at fixed order.* hp grading fixed it: order +1 per level, with ratio 0.15 over 6 levels.
- *Wide corners.* The vertex rule's error reaches 4e-4 at 150°. Split such corners into pieces of at most 75°.
- *Contacts.* T-junctions, gaps and sliver edges need a contact tier.
  - Before it, the column defect was 2.3e-3.
  - When all three corners are within tolerance, fall back to the self rule. A 4-way split recursion made things worse.
- *The explicit right-hand side of the free wake must use the assembly's own tiers.* A fixed 7-point rule (the `external=` hook) put the lift 1.6% low when a quarter-panel trailing-edge panel brought wake vorticity close to the trailing edge.

**Signs and singularities:**
- *Filament orientation.* A doublet layer is bounded by a filament running clockwise about its normal. Reversed tip filaments gave CL 0.54 and nearly two-dimensional loading. Check new sign conventions against a physical loading shape.
- *Filaments on the body.* In a Galerkin row the swirl of a body filament is ~1/d normal to the panel. On a free edge this has no finite discrete equation. So:
  - tips must be capped (zip and open are for debugging);
  - the zero-tip Kutta condition removes the tip filament and the cap-seam jump.
  - MESA avoids the problem because its potential formulation sees doublets through their bounded solid angle.
- *Kernel near edges.* Form the log term without cancellation. MESA's `vel_srcpan` (AeroSingularityElements) zeroes the log term within 1e-8 of an edge, so it was not reused.
- *Border scaling.* Per-entry scaling of the border vector inflated cond₁ to 7e4. Scale both border vectors to a 2-norm equal to mean|diag|.

**Mesh abuse:**
- Jitter bounded by the shortest edge collapsed triangles to 176°. Bound it by each node's clearance instead.

**Quads (stage 7):**
- Newton projection started from the centre missed on trapezoids.
- Abandoning a step that left the parameter square returned an edge point 5 mm away, which gave Cp −33 at leading-edge slivers. The fix: start from a 5×5 grid and clamp the steps.
- Near-degenerate quads (corners at 5°/175°) still integrate to only about 1e-3. This is accepted and documented.

**Julia:**
- Submodules export clashing names: `sheet_velocity`, `wake_elements`, `trefftz_drag`, and `RigidWake` against MESA's.
  - Check exports across modules and MESA before naming anything.
  - Extending an imported type's constructor needs an explicit `import`.
- Two silent bugs in the validation scripts:
  - `3E2`, meant as 3·E2, parses as the literal 300.0;
  - CSV labels containing commas shift the columns.

**Free wake (phase 2):**
- *Particle cores.* Tie the core to the element (overlap × half the longest edge), not to the area: area-based cores leave spanwise gaps at small dt.
- *Conversion distance.* Count it in cores behind the trailing edge, not in rows. Particles created within about 2 cores bias the lift up (+1.9% with immediate conversion).
- *Total vorticity drifts* without stretching: 0.7% at L1 and 4% at L2 after 5 spans. The impulse force route is therefore unreliable in free mode.

## 6. Open items

From the brief:
- **Deferred items not started:**
  - body motion (the only hook is a constant `V_body` in the triangle path's sources);
  - the symmetric Galerkin (Nédélec) form;
  - fast summation with GMRES. Our cond₁ is about 13× MESA's, so GMRES needs a preconditioner. Direct particle sums dominate the free-wake cost: 1.8 s per step at L2 with 33k particles.
  - the pressure-type Kutta condition, which matters more now that flow is unsteady;
  - the one-particle-per-quad option.
- **Small leftovers:**
  - the leak diagnostic per quadrature tier, and the net flux through caps and stitched gaps;
  - a test for the `external=` velocity hook;
  - the optional Bempp-cl matrix cross-check, and checking the brief's prior-art references.

From phase 2:
- vortex stretching;
- diffusion and remeshing;
- quads in the time loop (`explicit_rhs` is triangles only);
- a free trailing-edge panel;
- ground effect.

From `notes.md` and stage 7:
- a quad split that doesn't flip under small perturbations; this is also what the zero-tip default decision waits on;
- faster quad assembly: most of the cost is in the touching pairs;
- splitting badly distorted quads into triangles.

## 7. Documentation conventions

The README and the guide share one style. Keep it when editing them:
- short sentences, with bold lead-ins on bullets;
- tables of numbers, each from a named script with its outputs under `validation/output/`;
- percentages with their reference and route stated;
- a finding that turns out wrong gets corrected in place, saying what changed (e.g. the stage-3 λ claim, the midspan dent).

The full session transcript is `C:\Users\MarkLau\.claude\projects\c--Users-MarkLau-GitHub-VortexPanel-jl\96ad2a58-eb70-46da-b583-96c6946b3faa.jsonl`. It is large, so search it by byte offset (§3).
