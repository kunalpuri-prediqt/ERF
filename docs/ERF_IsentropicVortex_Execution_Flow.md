# ERF execution flow (isentropic vortex, no acoustic substepping, hybrid 5/6 advection)

This note cross-references:
- `docs/ERF_Discretization_Explained.md` (conceptual discretization note), and
- the implementation in `Source/*` for the concrete execution path and formulas.

Scope assumed here:
- Compressible dycore (not anelastic),
- `erf.substepping_type = None` (no acoustic substeps),
- `Blended_5th6th` advection for dycore/scalars,
- Isentropic vortex problem initialization.

---

## 1) High-level call sequence for one coarse-level timestep

```mermaid
flowchart TD
    A[ERF::TimeStep] --> B[ERF::Advance(lev,time,dt)]
    B --> C[Swap old/new state pointers]
    C --> D[FillPatch + VelocityToMomentum]
    D --> E[advance_dycore(...)]

    E --> F[set MRI callbacks: slow_rhs_pre, no_substep, slow_rhs_post]
    F --> G[mri_integrator.advance(S_old,S_new,time,dt)]

    G --> H[RK stage 1: dtau = dt/3]
    G --> I[RK stage 2: dtau = dt/2]
    G --> J[RK stage 3: dtau = dt]

    H --> K[slow_rhs_pre(F_slow)]
    I --> K
    J --> K

    K --> L[no_substep_fun: update fast vars only]
    L --> M[slow_rhs_post: update slow scalars]

    M --> N[final S_new at t^{n+1}]
```

Implementation anchors:
- Entry and prep: `ERF::Advance` in `Source/TimeIntegration/ERF_Advance.cpp`.
- Dycore advance + MRI wiring: `Source/TimeIntegration/ERF_AdvanceDycore.cpp`.
- MRI stage logic (no-substepping branch): `Source/TimeIntegration/ERF_MRI.H`.
- No-substep update kernel: `Source/TimeIntegration/ERF_TI_no_substep_fun.H`.

---

## 2) How “no acoustic substepping” changes RK stage evolution

In `MRISplitIntegrator::advance` (compressible path):

- Stage 1 (`nrk=0`): if `no_substepping`, `nsubsteps = 1`, `dtau = dt/3`.
- Stage 2 (`nrk=1`): if `no_substepping`, `nsubsteps = 1`, `dtau = dt/2`.
- Stage 3 (`nrk=2`): if `no_substepping`, `nsubsteps = 1`, `dtau = dt`.

Then, instead of acoustic subcycling, it calls:

\[
\texttt{no\_substep}(S_{sum}, S_{old}, F_{slow}, t_{fp}, \Delta t_{stage}, nrk).
\]

### Fast variables that are updated in this branch

From `ERF_TI_no_substep_fun.H`:
- `scomp_fast = {0,0,0,0}`
- `ncomp_fast = {2,1,1,1}`

So only:
- conserved components 0..1: \(\rho\), \(\rho\theta\),
- and face momenta \(\rho u, \rho v, \rho w\)

are advanced by this no-substep update kernel.

In the non-moving-terrain case, the update is explicit Euler over stage step:
\[
S^{sum} = S^{old} + \Delta t_{stage}\,F_{slow}.
\]

(implemented component-wise for the fast components in the kernel).

---

## 3) Where hybrid 5/6 is selected in code

`Blended_5th6th` input string is converted to `AdvType::Upwind_5th` in `ERF_AdvStruct.H`, and blending is controlled by `*_upw_frac`.

So operationally, “hybrid 5/6” means:
- use the 5th-order upwind-capable stencil class (`UPWIND5`),
- with tunable upwind fraction (`upw_frac`) multiplying the dissipative term.

This is consistent with the conceptual description in `docs/ERF_Discretization_Explained.md`.

---

## 4) Direct formulas extracted from implementation

## 4.1 Isentropic vortex initialization formulas

From `Source/Prob/ERF_InitCustomPert.H`:
\[
\Omega(x,y) = \beta\exp\left(-\frac{r^2}{2\sigma^2}\right), \quad
r^2 = \frac{(x-x_c)^2 + (y-y_c)^2}{R^2}.
\]

From `Source/Prob/ERF_InitCustomPert_IsentropicVortex.H`:
\[
\Delta T = -\frac{\gamma-1}{2\sigma^2}\,\Omega^2,
\]
\[
\rho_{norm} = (1+\Delta T)^{1/(\gamma-1)},
\]
\[
\rho' = \rho_{norm}\rho_0 - \rho_{hse},
\]
\[
p = \frac{\rho_{norm}^{\Gamma}}{\Gamma}\,\rho_0 a_\infty^2,
\]
\[
\rho\theta = \rho_0\rho_{norm}\left(T\left(\frac{p_0}{p}\right)^{R_d/c_p}\right), \quad T=(1+\Delta T)T_\infty,
\]
\[
(\rho\theta)' = \rho\theta - (\rho\theta)_{hse}(p_{hse}).
\]

From `Source/Prob/ERF_InitCustomPertVels_IsentropicVortex.H`:
\[
u = \left(M_\infty\cos\alpha - \frac{y-y_c}{R}\Omega\right)a_\infty,
\]
\[
v = \left(M_\infty\sin\alpha + \frac{x-x_c}{R}\Omega\right)a_\infty.
\]

## 4.2 Hybrid 5/6 face interpolation used in advection

From `Source/Utils/ERF_Interpolation_UPW.H`, `UPWIND5::Evaluate`:

Define
\[
a_1=s+s_{-1},\ a_2=s_{+1}+s_{-2},\ a_3=s_{+2}+s_{-3},
\]
\[
d_1=s-s_{-1},\ d_2=s_{+1}-s_{-2},\ d_3=s_{+2}-s_{-3}.
\]
Then
\[
\phi_{face}=
\frac{37}{60}a_1-\frac{2}{15}a_2+\frac{1}{60}a_3
-\texttt{upw}\,\frac{1}{60}(d_3-5d_2+10d_1).
\]

In `InterpolateIn*`, `upw` is normalized to \(\pm 1\) from advecting momentum sign, then scaled by `m_upw_frac`:
\[
\texttt{upw} = \operatorname{sign}(u_{face})\cdot \texttt{upw\_frac}.
\]

So equivalently:
- `upw_frac = 0` \(\Rightarrow\) centered 6th part only,
- `upw_frac = 1` \(\Rightarrow\) full 5th upwind dissipation.

For reference, `CENTERED6::Evaluate` in same file is:
\[
\phi_{face}=\frac{37}{60}a_1-\frac{2}{15}a_2+\frac{1}{60}a_3.
\]

## 4.3 Scalar flux and divergence form actually used

From `Source/Advection/ERF_AdvectionSrcForScalars.H` + `ERF_AdvectionSrcForState.cpp`:

Face fluxes for a scalar-like component use
\[
F_x = (\overline{\rho u})_{face}\,\phi_{face},\quad
F_y = (\overline{\rho v})_{face}\,\phi_{face},\quad
F_z = (\overline{\rho w})_{face}\,\phi_{face},
\]
where \(\phi_{face}\) is from the interpolation above.

Tendency is assembled as
\[
\mathcal{A}_\phi
= -\frac{m_xm_y}{J}\left[
\frac{F_x(i+1)-F_x(i)}{\Delta x}+
\frac{F_y(j+1)-F_y(j)}{\Delta y}+
\frac{F_z(k+1)-F_z(k)}{\Delta z}
\right].
\]

This is exactly what is coded in `AdvectionSrcForScalars` before additional source-term contributions are applied.

## 4.4 Constant-\(\Delta z\) momentum advection form (representative simple case)

From `Source/Advection/ERF_AdvectionSrcForMom_ConstantDz.cpp`:

For each momentum component, code computes directional fluxes (e.g. `xflux_hi - xflux_lo`, etc.) and sets
\[
\mathcal{A}_{\rho u} = -\left[(\Delta_x F_x)+(\Delta_y F_y)+(\Delta_z F_z)\right],
\]
(and analogously for \(\rho v,\rho w\), with map-factor multipliers in x/y contributions).

When higher-order scheme is active, templated momentum advection uses `UPWIND5`/`CENTERED6` interpolation classes through `AdvectionSrcForMomVert_N<...>`.

---

## 5) End-to-end stage flow (simple isentropic-vortex run)

For each RK stage in no-substep mode:

1. `slow_rhs_pre(...)`
   - refresh primitives as needed,
   - builds scalar sources (`make_sources`), pressure-gradient (`make_gradp_pert`), momentum sources (`make_mom_sources`), and advection tendencies (via advection kernels).

2. `no_substep_fun(...)`
   - updates fast variables \((\rho,\rho\theta,\rho u,\rho v,\rho w)\) by stage \(\Delta t\).

3. `slow_rhs_post(...)`
   - updates remaining slow components (e.g., extra scalars/tracers), using stage-consistent state.

After stage 3, `S_new` is the advanced solution at \(t^{n+1}\).

---

## 6) Quick reconciliation with `ERF_Discretization_Explained.md`

The implementation matches the note’s structure:
- C-grid staggering and face-flux divergence form are directly reflected in scalar/momentum advection kernels.
- Hybrid 5/6 is implemented as 6th-order centered core plus an upwind dissipation term scaled by `upw_frac`.
- For no acoustic substepping, the MRI driver bypasses acoustic inner loops and calls a stagewise fast-variable update once per RK stage.

