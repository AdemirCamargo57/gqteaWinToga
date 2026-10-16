# Fukui Functions: Calculation, Multiwfn Workflow and the `fukuiFunctionAnalyze.py` Tool — Design

**Date:** 2026-10-16 · **Status:** design approved for planning · **Implementation plan:** [../plans/2026-10-16-fukui-function-analyze.md](../plans/2026-10-16-fukui-function-analyze.md)

Fukui functions show **where a molecule's electron density responds most when it gains or loses an electron**. They are used to rank sites by how susceptible each is to attack. This document treats them as **indicators of likely reactive sites, not predictions of reaction outcome**: sterics, solvent, transition-state energies, kinetics and thermodynamics are all outside the model (§5).

## 0. How claims are marked

| Mark | Meaning |
|---|---|
| **[M §x.y]** | Stated in the supplied manual, *Multiwfn_manual_2026.7.10.pdf*, section x.y. Page numbers are the manual's printed numbers. |
| **[RUN]** | Checked in this session by running **Multiwfn 2026.7.15** (Windows binary, already installed) and **Gaussian 09W** on the phenol example. Transcripts: §3.6. |
| **[WEB n]** | An online source in §8. The two sobereva.com articles were written by Multiwfn's author. |
| **[REC]** | A recommendation of this document: reasoned, but not stated by the manual or the sources. |
| **[NOT SUPPORTED]** | Not provided by the gQTEA tool. Where Multiwfn offers it, the Multiwfn route is given instead. |

The manual's menu numbers were checked against the running program: every option quoted in §3 appears exactly as printed by Multiwfn 2026.7.15 [RUN]. The manual is dated 2026.7.10 and the binary 2026.7.15, and the two agree on every menu used here.

---

## 1. Definitions and interpretation

### 1.1 The Fukui function

The Fukui function is the derivative of the electron density with respect to the number of electrons at fixed external potential v (the nuclei, for an isolated molecule): f(r) = [∂ρ(r)/∂N]_v. Source: Parr & Yang 1984 [WEB 1]; [M §4.5.4.1, p. 521]. The derivative is discontinuous at integer N, so three finite differences are used [M §4.5.4.1, p. 521]:

| Function | Finite difference | Frozen-orbital limit | Large where the site is susceptible to… |
|---|---|---|---|
| **f⁺(r)** | ρ_{N+1}(r) − ρ_N(r) | \|φ_LUMO\|² | **nucleophilic attack** (the molecule *accepts* electrons there) |
| **f⁻(r)** | ρ_N(r) − ρ_{N−1}(r) | \|φ_HOMO\|² | **electrophilic attack** (the molecule *donates* electrons from there) |
| **f⁰(r)** | [ρ_{N+1}(r) − ρ_{N−1}(r)]/2 = (f⁺ + f⁻)/2 | (\|φ_HOMO\|² + \|φ_LUMO\|²)/2 | **radical attack** (little net electron transfer) |

**Wording rule.** "Susceptible to electrophilic attack" means the site is *nucleophilic*: an incoming electrophile attacks it. The Multiwfn author uses the same convention [WEB 3, §0]. The labels "f⁻ ↔ electrophilic attack" and "f⁺ ↔ nucleophilic attack" describe the **attacking reagent**, not the site.

**Interpretation.**
- Sites with larger f values are favoured for the matching attack [M §4.5.4.1, p. 521].
- Fukui functions compare sites **within one molecule**. They are not size-consistent: f integrates to 1, so it thins out as the molecule grows. To compare different molecules, use local softness s(r) = S·f(r), or for the dual descriptor the local hyper-softness S²·Δf(r) [M §3.25.1, pp. 317–318]; [WEB 3, §2.7].
- **The finite-difference form is preferred over the frozen-orbital (HOMO/LUMO density) form.** The Multiwfn author strongly advises against the frozen-orbital form, because it neglects orbital relaxation and is no cheaper in practice [WEB 3, §2.1].
- **Reliability differs between the functions.** In the author's comparison studies, f⁻ predicts electrophilic-attack sites reasonably well, while f⁺ is *poor* at predicting nucleophilic-attack sites; the dual descriptor, or Hirshfeld charges, do better there [WEB 3, §2.1, §2.2, §4]. This is based on [WEB 9] and [WEB 10].

### 1.2 Dual descriptor

Δf(r) = f⁺(r) − f⁻(r) = ρ_{N+1} − 2ρ_N + ρ_{N−1}. Sources: Morell, Grand & Toro-Labbé 2005 [WEB 2]; [M §4.5.4.2, p. 523].

- Δf > 0: favourable for **nucleophilic** attack.
- Δf < 0: favourable for **electrophilic** attack.

One function shows both kinds of site. The manual recommends comparing sites by **how positive or how negative** Δf is, rather than reading the sign alone [M §4.5.4.2, p. 523]. A spin-density approximation, Δf ≈ ρ^s_{N+1} − ρ^s_{N−1}, is described in [M §4.5.4.2, p. 523]. Main function 22 always evaluates the exact three-density form [M §4.22.1, p. 904].

The **dual-descriptor potential** DDP(r) = V_{f⁺} − V_{f⁻} is smoother and "in principle more rigorous" than Δf, but costs one ESP evaluation per state [M §3.25.1, p. 320]; [M §4.A.4, p. 1064]. **[NOT SUPPORTED]** in gQTEA. Use Multiwfn main function 22, option 9 [M §3.25.2]; [RUN: listed].

### 1.3 Condensed (atom) indices

Condensing replaces the density by atomic populations p_A. With q_A = Z_A − p_A, the nuclear charges cancel [M §4.7.3, pp. 540–541]:

```
f⁺_A = q_A(N) − q_A(N+1)        nucleophilic attack
f⁻_A = q_A(N−1) − q_A(N)        electrophilic attack
f⁰_A = [q_A(N−1) − q_A(N+1)]/2  radical attack
Δf_A = f⁺_A − f⁻_A              (CDD in Multiwfn's output)
```

The condensed form goes back to Yang & Mortier 1986 [WEB 4]. **Hirshfeld charges are the recommended population scheme** [M §4.7.3, p. 541]; [M §4.A.4, p. 1063]. The comparison in [WEB 6] found Hirshfeld the best of 15 charge schemes. Multiwfn's main function 22 uses Hirshfeld automatically [M §3.25.2, p. 321]. Mulliken charges should not be used [M §4.A.4, p. 1063], and NPA is discouraged [WEB 5] [WEB 7].

**Real-space vs condensed.**

| | Real-space f(r) | Condensed f_A |
|---|---|---|
| What it is | A 3D function: isosurfaces, or a map on a molecular surface | One number per atom |
| Strength | Shows the *shape and region*: above or below a π face, a lone-pair lobe, a bond midpoint | Quantitative ranking between atoms; easy to tabulate |
| Weakness | Visual reading is "somewhat ambiguous and subjective" [M §4.5.4.2, p. 525] | Depends on the partition scheme (§5.3); cannot show directionality; can come out negative (§5.3) |
| Quantitative alternative | Surface extrema, and averages on local atomic surfaces [M §4.12.4] | — |

### 1.4 Global indices used alongside

These are from [M §3.25.1, p. 317], all at the N-state geometry:

- vertical IP: I = E(N−1) − E(N)
- vertical EA: A = E(N) − E(N+1)
- electronegativity χ = (I+A)/2; chemical potential μ = −χ
- hardness **η = I − A**. Multiwfn drops the conventional 1/2 [M §3.25.1, p. 317].
- softness S = 1/η
- electrophilicity ω = μ²/2η
- nucleophilicity N_Nu = E_HOMO − E_HOMO(TCE)

Multiwfn uses **E_HOMO(TCE) = −0.335198 Hartree**, computed at B3LYP/6-31G*. At any other level, recompute TCE and evaluate N_Nu by hand [M §3.25.2, p. 321]. All of these formulas were checked digit-for-digit against Multiwfn's own `CDFT.txt` for phenol [RUN].

---

## 2. Required quantum-chemical calculations

### 2.1 States and geometry

1. **Optimize the N-electron state**, i.e. the species as it reacts. Use one fixed geometry for all states [M §3.25.1, p. 317]; [M §3.25.2, p. 321].
2. **Confirm it is a minimum with a frequency calculation.** [REC]
3. **Run single points for N, N+1 and N−1 at that same geometry** with the **same method and basis set** [M §3.25.2, p. 321]. Do **not** optimize the ions: the derivative is taken at constant external potential (fixed nuclei) [M §4.5.4.1, p. 522].
4. Flexible molecules: analyse the lowest-energy conformer, but note that some reactions proceed through a higher-energy conformer [WEB 3, §0]. [REC: a conformer search before step 1.]

### 2.2 Charge and spin multiplicity

| N state | N | N+1 | N−1 | Source |
|---|---|---|---|---|
| Neutral, closed shell | 0 1 | −1 2 | +1 2 | Multiwfn defaults [M §4.22.1, p. 902]; [RUN: option 1 output] |
| Neutral radical | 0 2 | −1 1 | +1 1 | [WEB 5, §4] |
| Cation radical | +1 2 | 0 1 | +2 1 | [WEB 5, §4] |
| Anion radical | −1 2 | −2 1 | 0 1 | [WEB 5, §4] |
| Closed-shell cation | +1 1 | 0 2 | +2 2 | [WEB 5, §4] |
| Closed-shell anion | −1 1 | −2 2 | 0 2 | [WEB 5, §4] |
| Degenerate LUMO (p) / HOMO (q) | — | N+p, multiplicity p+1 | N−q, multiplicity q+1 | for a closed-shell N [M §3.25.4, p. 324] |

The rows above assume ordinary organic systems. For **transition-metal systems**, compute the lowest few spin states of N+1 and N−1 and use the lowest in energy [WEB 5, §4].

**Consistency checks the gQTEA tool runs:**
- identical atoms and coordinates;
- identical basis set;
- identical method once the R/U/RO prefix is stripped;
- N+1 holds more electrons than N and N−1 fewer.

A multiplicity change other than the expected ±1 is reported as a note, not an error (§4 of the plan).

### 2.3 Level of theory

- **Minimum level.** The manual calls B3LYP/6-31G\* "the lowest acceptable level for producing meaningful result" [M §4.5.4.1, p. 521]. It is also the default route written by Multiwfn's input generator: `#P B3LYP/6-31G* SCF=conver=7 out=wfn` [RUN].
- **Diffuse functions cut both ways.** They are needed for accurate *energies* of anions, and so for EA, χ and ω; the manual suggests at least 6-311+G\* [M §4.22.1, p. 903]. But with a self-interaction-prone functional such as B3LYP, they make the N+1 density over-diffuse (Rydberg-like), and **f⁺ and Δf lose chemical meaning** [M §4.22.3.1, p. 913].
- **The manual's recommendation for rigorous f⁺/Δf** is a long-range-corrected functional (ωB97XD) with a triple-ζ basis *without* diffuse functions, e.g. 6-311G\* [M §4.22.3.1, p. 913]. Even a triple-ζ basis without diffuse functions gives a meaningful dual descriptor [M §4.22.3.1, p. 913], citing *Chem. Phys. Lett.* 724, 29 (2019).
- **[REC] Pragmatic split:** take the real-space and condensed functions from the no-diffuse level. If accurate global indices matter, compute E(N±1) separately with diffuse functions, and state in the report which level each number comes from.
- **[REC] Solvent:** if the reaction occurs in solution, use the same implicit-solvent model for all three states. A solvated anion is more often bound than a gas-phase one. Whether this helps a given system is not covered by the manual.
- **Anion SCF convergence.** The N+1 state often converges with difficulty; add convergence aids [WEB 5, §3.1].

### 2.4 Files

| Program | Accepts | Notes |
|---|---|---|
| Multiwfn main function 22 | `.wfn`, `.wfx`, `.fch`/`.fchk`, `.mwfn` holding wavefunction **and** energy [M §3.25.2, p. 321] | It looks for `N.wfn`, `N+1.wfn`, `N-1.wfn` in the current folder. If any is missing it asks for a path, and `.fch` is accepted [M §3.25.2]; [RUN]. |
| Multiwfn, orbital-weighted Fukui (options 4–7) | `.mwfn`, `.fch`, `.molden`, `.gms` (needs virtual orbitals) [M §3.25.3, p. 323] | Closed-shell single determinant only; no diffuse functions [M §3.25.3, p. 323]. |
| gQTEA `fukuiFunctionAnalyze` | Gaussian `.fchk`/`.fch` (three files) | Reads `Total SCF Density`, `Total Energy` and orbital energies through the project's verified fchk reader. |

**Gaussian route.** Use `%chk=` plus `formchk name.chk name.fchk`, or add `out=wfn` as Multiwfn's generated inputs do [RUN]. ORCA users can let Multiwfn write ORCA inputs: option −2, then option 1 [M §3.25.2, pp. 321–322]; [RUN: option −2 listed].

---

## 3. Multiwfn workflow, verified

### 3.1 Main function 22: the recommended route

The manual "strongly suggest[s]" main function 22 for daily work over manual difference densities [M §4.5.4, pp. 520–521]; [M §4.7.3, p. 542]. The CDFT menu exactly as printed by Multiwfn 2026.7.15 [RUN]:

```
 -3 Set degree of FMO degeneracy, current: Nondegenerate
 -2 Choose the quantum chemistry program used in option 1, current: Gaussian
 -1 Toggle calculating electrophilicity index (w_cubic) and electrophilic descriptor (epsilon) by option 2
 0 Return
 1 Generate .wfn files for N, N+1, N-1 electrons states
 2 Calculate various quantitative indices
 3 Calculate grid data of Fukui function, dual descriptor and related functions
 4 Set delta in orbital-weighted (OW) calculation, current: 0.1000 a.u.
 5 Print current orbital weights used in orbital-weighted (OW) calculation
 6 Calculate condensed OW Fukui function and OW dual descriptor
 7 Calculate grid data of OW Fukui function and OW dual descriptor
 8 Calculate nucleophilic and electrophilic superdelocalizabilities
 9 Calculate grid data of Fukui potential and dual descriptor potential
 10 Calculate bond dual descriptor (BDD)
```

**Step by step:**

1. Load a file that holds the optimized N geometry. Any format with atoms works, e.g. `.xyz` [M §3.25.2, p. 320]. Then enter `22`.
2. **Option 1** writes `N.gjf`, `N+1.gjf`, `N-1.gjf`.
   - Press ENTER for the default B3LYP/6-31G\*, or type your own keywords.
   - Press ENTER for the default charge/multiplicity (0 1), (−1 2), (1 2), or give your own.
   - If `gaupath` is set in `settings.ini`, Multiwfn can run Gaussian itself [M §3.25.2, p. 321]; [M §4.22.1, p. 902]; [RUN].
   - For a mixed basis, put the basis definition in `basis.txt` [M §3.25.2, p. 321].
3. Run the three jobs, then place `N.wfn`, `N+1.wfn`, `N-1.wfn` in the folder, or give the paths when prompted [RUN].
4. **Option 2** writes `CDFT.txt` [M §3.25.2]; [RUN]. It contains:
   - Hirshfeld q(N), q(N+1), q(N−1);
   - condensed f⁻, f⁺, f⁰, CDD;
   - condensed local electrophilicity and nucleophilicity (e·eV);
   - condensed s⁻, s⁺, s⁰, s⁺/s⁻, s⁻/s⁺ and s⁽²⁾;
   - E(N), E(N±1), E_HOMO(N), E_HOMO(N±1), IP, EA, χ, μ, η, S, S², ω and N_Nu.

   Hirshfeld integration uses the Becke grid; `radpot`=75 and `sphpot`=434 are the defaults in `settings.ini` [RUN: printed "Radial grids: 75 Angular grids: 434"].
5. **Option 3** computes grids of ρ for the three states, then shows this grid-setup menu [RUN]:

   ```
   -10 Set extension distance of grid range for mode 1~4, current:  3.000 Bohr
   1 Low quality grid,    covering whole system, about 125000 points in total
   2 Medium quality grid, covering whole system, about 512000 points in total
   3 High quality grid,   covering whole system, about 1728000 points in total
   4 Input the number of points or grid spacing in X,Y,Z, covering whole system
   ... 5-8, 10, 11, 100 (other grid definitions)
   ```

   followed by this post-processing menu [RUN]:

   ```
   -1 Set the scale factor to various grid data, current: 1.000000
   0 Return
   1 Visualize isosurface of f+          5 Export grid data of f+ as f+.cub in current folder
   2 Visualize isosurface of f-          6 Export grid data of f- as f-.cub in current folder
   3 Visualize isosurface of f0          7 Export grid data of f0 as f0.cub in current folder
   4 Visualize isosurface of dual descriptor   8 Export grid data of dual descriptor as DD.cub in current folder
   ```

   Set the scale factor (option −1) to S to get s(r), to S² for the local hyper-softness, or to ω or N_Nu for the local indices [M §3.25.2, p. 321]; [M §4.22.1, pp. 904–905].
6. **Options for special cases:**
   - **−3**: degenerate frontier orbitals (§5.5).
   - **4–7**: orbital-weighted f_w and Δf_w, from one `.fch`/`.mwfn` [M §3.25.3].
   - **9**: Fukui potential and DDP [M §3.25.2]; [M §4.22.4].
   - **10**: bond dual descriptor, from Hirshfeld fuzzy bond orders [M §3.25.1, p. 318]; [M §4.22.5].
   - **−1**: needs an N−2 state, for ω_cubic and ε [M §3.25.2, p. 322].

### 3.2 Manual alternative: main function 5, custom operation [M §4.5.4]

The route for f⁻ with `phenol.wfn` loaded [M §4.5.4.1, p. 522]:

```
5 → 0 (custom operation) → 1 → -,phenol_N-1.wfn → 1 (electron density) → 2 (medium grid) → -1 (isosurface)
```

- **f⁰:** load N+1, subtract N−1, then option 6 "Divide all grid data by a factor", 2 [M §4.5.4.1, pp. 522–523].
- **Exact Δf:** load N+1; operations `-,N`, `-,N`, `+,N-1` [M §4.5.4.2, p. 524].
- **Spin-density Δf:** load N+1, `-,N-1`, then property 5 (spin density) [M §4.5.4.2, pp. 523–524].

### 3.3 Fukui function on a molecular surface: main function 12 [M §4.12.4, pp. 677–679]

**Use a ρ = 0.01 a.u. surface**, because f is "often too small" on the default 0.001 surface. Two Multiwfn processes are needed:

- **Process A:** `12 → 2 → 0` (mapped function from an external file) `→ 1 → 1 → 0.01 → 0`, which writes `surfptpos.txt`.
- **Process B:** `5 → 0 → 1 → -,N-1 → 1 → 100 → surfptpos.txt → t.txt`.
- **Back in A:** type `t.txt`. Extrema are listed.
  - Option **0** shows them.
  - Option **11** gives f averages on each atom's local surface.
  - Option **2** exports `surfanalysis.pdb`.

The gQTEA tool does the extrema part in one step (§4). It does **not** do the local-surface averages [NOT SUPPORTED].

### 3.4 Orbital contributions [M §4.200.13.1]

Main function 200, subfunction 13 decomposes an f⁻ cube into MO or NBO contributions. **[NOT SUPPORTED]** in gQTEA.

### 3.5 Visualization outside Multiwfn

- Exported cubes plot directly in VMD [M §4.22.1, p. 904].
- Atoms can be coloured by condensed values via a `.pqr` file [M §4.A.10]. The manual points to sobereva.com/425 for that procedure.

### 3.6 Session verification transcript [RUN]

- **Multiwfn 2026.7.15 on its own `examples/phenol*.wfn`.** Option 2 reproduced the manual's `CDFT.txt` [M §4.22.1, pp. 902–903] digit for digit; s columns differ in the 4th decimal at a few atoms.
- **Phenol at the same geometry, G09W `#p B3LYP/6-31G(d) scf=tight`,** charge/multiplicity 0 1, −1 2, 1 2:
  - Multiwfn option 2 on the three `.fchk` files reproduces the manual's condensed table to ±0.0001.
  - E(N) = −307.464870 Ha, against −307.464860 Ha for the G16-derived example file.
- **gQTEA-native f⁺, f⁻, f⁰ and Δf** (density from `gaussianWavefunction.electron_density`) **agree with Multiwfn's exported cubes** on the identical 64×75×28 grid to **≤ 8×10⁻⁶ a.u.** The largest value of |f| is 0.083.
- **Global indices** computed from the three fchk energies agree with `CDFT.txt` to every printed digit.
- **Scripted Multiwfn route.** It works with this stdin script: `N.fch` given on the command line, `settings.ini` copied into the working folder, then `22 / 2 / N.fch / N+1.fch / N-1.fch / 0 / q`. Setting the `Multiwfnpath` environment variable did **not** make this Windows build find `settings.ini`, although the manual describes that variable [M §2.1]. The tool therefore copies the file.

---

## 4. Visualization and interpretation procedure

1. **Read the global indices first.** IP, EA and η show how hard the molecule is to oxidize or reduce. A negative vertical EA, as for phenol (−2.21 eV at B3LYP/6-31G\*), means the N+1 state is unbound at this level. [REC] Treat f⁺ and Δf with extra caution then (§5.1).
2. **Look at isosurfaces of f⁻, f⁺, f⁰ and Δf.**
   - The manual's isovalues are 0.007 a.u. (phenol), 0.01 (propylene f⁰), 0.005 (benzene) and 0.0003 (orbital-weighted, C₆₀) [M §4.5.4, §4.22]. The gQTEA default is 0.005 [REC].
   - Lower the isovalue until the competing sites can just be told apart [M §4.5.4.2, p. 524].
   - Note *where* each lobe sits (π face, lone pair, σ-hole direction) as well as *which atom* it is on. Examples: the toroidal negative Δf_w around Cl in CH₃Cl, and the positive Δf_w at both ends of the C–Cl bond [M §4.22.2, p. 910].
3. **Map each function onto the ρ(N) = 0.01 a.u. surface and list its extrema** [M §4.12.4]. The largest surface maxima of f⁻ (or the most negative minima of Δf) mark the most exposed, most reactive regions for electrophilic attack. The same reasoning applies to f⁺ and maxima of Δf for nucleophilic attack.
4. **Compare condensed indices between plausible sites only.** Ignore atoms that cannot react: the manual excludes the ipso C4 of phenol [M §4.5.4.2, p. 524], and [WEB 3, §0] says to ignore sites where a reaction is known to be impossible. Rank by f⁻ (electrophilic attack), f⁺ (nucleophilic) or f⁰ (radical). Cross-check with CDD: most negative means electrophilic attack, most positive means nucleophilic attack.
5. **Agreement between descriptors raises confidence; disagreement means a transition-state calculation is needed.** If f⁻, Δf, the surface extrema and (optionally) ALIE or ESP all point to the same site, that site is a strong candidate. If they disagree (e.g. one says para, another meta), the reliable route is to locate transition states and compare free-energy barriers [WEB 3, §1].
6. **Worked reference: phenol, B3LYP/6-31G\*** [M §4.7.3, p. 542]; [RUN]. Among ring carbons:
   - f⁻ is largest at C1 (para, 0.1439), and also high at the oxygen (0.1460);
   - CDD is most negative at C1 (−0.0841);
   - the meta carbons C2 and C6 have the smallest f⁻ and the most positive CDD.

   The predicted ortho/para direction for electrophilic substitution agrees with the known behaviour of OH as an ortho/para director [M §4.5.4.1, p. 522]. The gQTEA test suite pins this ranking.

---

## 5. Methodological limitations

### 5.1 Electronic-structure method

- **Finite differences use one whole electron.** That is a large perturbation. The result depends on how well the method describes the ion, especially an unbound anion.
- **Self-interaction error** in B3LYP-like functionals over-delocalizes the extra electron of N+1 [M §4.22.3.1, p. 913]. Range-separated functionals reduce this.
- **The density used** here is the SCF (HF/DFT) density. For correlated methods, Gaussian may also store an MP2 or CC density, but the gQTEA tool always reads `Total SCF Density` and reports a note when another density is present.
- **[REC]** Check spin contamination in the unrestricted N±1 states (⟨S²⟩ in the Gaussian output). A heavily contaminated ion density contaminates f.

### 5.2 Basis set

See §2.3. Diffuse functions improve EA but can spoil f⁺ and Δf. Polarization functions are the minimum [M §4.5.4.1, p. 521].

### 5.3 Population scheme (condensed indices only)

- Condensed values depend on how space is divided among atoms. **Hirshfeld** is recommended [M §4.7.3]; [WEB 6]. Mulliken is the worst choice [M §4.A.4, p. 1063].
- **Negative condensed values** do occur with other schemes. Hirshfeld usually avoids them. Proposed causes include the partition itself and the large ΔN = 1 step [WEB 7] [WEB 8]. One view holds that a genuinely negative condensed Fukui function is chemically meaningful (reduction of the molecule oxidizing one centre) [WEB 11].
- **[REC]** Report negative values; do not rank on them.
- The real-space functions do not depend on any partition. The extrema in §4 step 3 sidestep this choice entirely.

### 5.4 Charged and open-shell systems

- Main function 22 works for charged and open-shell systems once the charges and multiplicities are set correctly [WEB 5, §4]; see the table in §2.2.
- **For an open-shell N state,** the Fukui functions above are *total-density* (spin-free) quantities. A spin-polarized CDFT framework (variations in N_α and N_β separately) exists, e.g. Fukui and dual-descriptor matrices in spin-polarized DFT [WEB 12]. **[NOT SUPPORTED]** in gQTEA and not part of main function 22's documented quantities.
- **Orbital-weighted forms and the degenerate f_Q form are closed-shell only** [M §3.25.3, p. 323]; [M §3.25.4, p. 324].
- **Multiply charged anions** (an N+1 state of an anion) are often unbound in the gas phase. [REC] Use an implicit solvent, or treat f⁺ and Δf as unreliable.

### 5.5 Degenerate or near-degenerate frontier orbitals

- The standard f⁻/f⁺ can break the molecular symmetry and mislead (benzene, C₆₀) [M §3.25.3, p. 322]; [M §4.22.3.1, p. 912].
- **Remedies.**
  - **N+p / N−q states** (f⁺_Q = [ρ_{N+p} − ρ_N]/p, f⁻_Q = [ρ_N − ρ_{N−q}]/q) with multiplicities p+1 and q+1. Multiwfn option −3 sets p and q [M §3.25.4, pp. 323–325]. Multiwfn proposes a degeneracy by counting orbitals within **0.01 eV** of the HOMO or LUMO [M §4.22.3.1, p. 911].
  - **Orbital-weighted f_w**, with δ = 0.1 Ha by default. It neglects orbital relaxation, so it should not be used when the frontier orbitals are clearly non-degenerate [M §3.25.3, pp. 322–323].
- The gQTEA tool detects degeneracy by the same 0.01 eV rule and reports it. It computes the real-space N+p/N−q functions when it is given such files. Condensed values in that case must come from Multiwfn option −3 [NOT SUPPORTED in the scripted engine].
- The manual's C₆₀ example creates files named `N-3.gjf` and `N+5.gjf` for LUMO degeneracy 3 and HOMO degeneracy 5 [M §4.22.3.2, p. 913], which looks inverted. **[REC]** Trust the program's own prompts over that sentence.

### 5.6 What Fukui functions cannot tell you

Collected in [WEB 3, §10]:
- steric hindrance;
- distortion of the reactants;
- non-covalent steering;
- special transition-state effects such as hyperconjugation;
- concerted processes;
- solvent effects;
- diffusion;
- reversibility, product stability, later steps and catalysis.

Fukui functions also do not rank *reactivity between molecules* without a global scale (S, ω, N_Nu). **They identify where a given molecule is most susceptible, not whether or how fast it reacts.**

---

## 6. Practical checklist

| # | Step | File(s) produced | Verified by |
|---|---|---|---|
| 1 | Optimize N; run frequencies [REC]; conformer search if flexible [REC] | `mol_opt.log`, `mol_opt.chk` | no imaginary frequencies |
| 2 | Choose the level (§2.3); note the diffuse vs no-diffuse trade-off | — | written into the report header |
| 3 | Single points N, N+1, N−1 at the N geometry, same route; charge/multiplicity from §2.2 | `mol_N.chk`, `mol_Np1.chk`, `mol_Nm1.chk` (+ `out=wfn` files if wanted) | `Normal termination` ×3; ⟨S²⟩ acceptable [REC] |
| 4 | `formchk` each `.chk` | `mol_N.fchk`, `mol_Np1.fchk`, `mol_Nm1.fchk` | gQTEA validation passes (atoms, geometry, basis, method, electron counts) |
| 5 | Check frontier-orbital degeneracy | — | gQTEA note, or Multiwfn option −3 listing |
| 6 | Multiwfn `22` → `2` (or gQTEA "Run Multiwfn") | `CDFT.txt` | f columns consistent with the q columns (gQTEA checks this to ±2×10⁻⁴) |
| 7 | Real-space grids | gQTEA: `mol_N_fplus.cub`, `_fminus.cub`, `_fzero.cub`, `_dual.cub`; Multiwfn: `f+.cub`, `f-.cub`, `f0.cub`, `DD.cub` | grid integrals ≈ 1 (f) and ≈ 0 (Δf) |
| 8 | Isosurfaces at 0.005–0.01 a.u.; maps on ρ = 0.01 a.u. with extrema | images; `mol_N_Fukui.txt` extrema tables | — |
| 9 | Rank plausible sites by f⁻ / f⁺ / f⁰ and CDD; compare with surface extrema | `mol_N_Fukui.txt` SITE RANKING | descriptors agree? |
| 10 | Write up with level, scheme (Hirshfeld), isovalues and caveats; cite Multiwfn and the CDFT chapter [M §3.25] | — | — |
| 11 | If descriptors disagree or the answer matters: transition-state barriers [WEB 3, §1] | — | — |

---

## 7. The gQTEA tool: `fukuiFunctionAnalyze.py` (summary; details in the plan)

**Decisions:**

1. **Input: three Gaussian `.fchk`/`.fch`** (N, N+1, N−1). This reuses `gaussianWavefunction`, which the ESP tool already verified against cubegen.
   - The real-space functions are computed natively. They match Multiwfn's cubes to ≤ 8×10⁻⁶ a.u. [RUN].
2. **[Superseded by §9 — kept for the record.]** **Condensed indices come from Multiwfn, never from a reimplementation.** Hirshfeld needs free-atom reference densities, which the project does not have; a home-made partition would silently disagree with the literature's values. Three sources:
   - **(a)** run Multiwfn as an external engine (scripted option 2, verified [RUN]);
   - **(b)** import an existing `CDFT.txt`;
   - **(c)** skip.

   As with the ESP tool's cubegen engine, a missing or failing Multiwfn **never aborts the run**: the real-space results are still produced, and a `NOTE:` is added.
3. **Outputs:**
   - `<stem>_Fukui.txt`: header; global indices; grid integrals; condensed table with s⁻, s⁺, s⁰; site ranking; surface extrema of f⁺, f⁻, f⁰ and Δf on ρ(N) = 0.01 a.u.; notes.
   - Four cube files, `<stem>_fplus.cub`, `_fminus.cub`, `_fzero.cub`, `_dual.cub` (switchable). Bohr grid, Gaussian cube layout; they open in VMD and in Multiwfn.
4. **3D view** reuses `ESPSurfaceViewer` through one backward-compatible addition, a `caption=` argument. Two display modes:
   - ±isovalue isosurfaces (blue +, red −);
   - the function mapped on ρ(N) = 0.01 a.u. with Multiwfn-rule extrema.
5. **Out of scope for gQTEA** [NOT SUPPORTED]; use Multiwfn instead:
   - Fukui potential and DDP (option 9);
   - orbital-weighted forms (4–7);
   - BDD (10);
   - ω_cubic and ε (−1);
   - condensed indices for N+p/N−q (−3, then 2);
   - local-surface averages (main function 12, option 11);
   - orbital decomposition (200 → 13);
   - spin-polarized CDFT.

---

## 9. Revision (2026-10-16, same day): no Multiwfn at run time, Numba, atom view

At the user's request the tool no longer uses Multiwfn. Decision 2 of §7 is reversed. The premise behind it ("Hirshfeld needs free-atom densities the project does not have") was removed by generating those densities. The evidence:

- **Free-atom densities** [RUN]. Multiwfn's built-in densities for Z ≤ 18 are UB3LYP/cc-pVQZ ground-state atoms, sphericalized [M §6.3, Appendix 3, pp. 1123–1124]. The same 18 atoms (H–Ar) were run in G09W and spherically averaged with the project's own fchk reader. They are stored as `proatomDensities.py` and rebuilt by `tools/make_proatom_densities.py`.
- **Hirshfeld integration** (`hirshfeldCharges.py`) uses a Becke grid of 75 radial × 434 Lebedev points per atom, the values Multiwfn printed [RUN], with populations normalised to N as in Multiwfn's "Final atomic charges".
- **Agreement with Multiwfn 2026.7.15** [RUN]:
  - max |Δq| = 1.4×10⁻⁴ e over 8 molecules spanning H, Li, B–F and Na–Cl;
  - ≤ 10⁻⁴ for phenol N, N+1 and N−1;
  - condensed f⁻ and f⁺ within 5×10⁻⁵ of `CDFT.txt`.
- **Scope.** Hirshfeld charges cover elements H–Ar. Heavier elements and ECP files give a NOTE; the real-space results are still produced. Condensed values for N+p/N−q states are now supported, with the same /p and /q formulas as §5.5.
- **Speed.** `densityNumba.py`, used automatically when Numba is installed and falling back to numpy, makes densities about 4× faster, with agreement ≤ 10⁻¹² of max ρ. Phenol at 0.25 bohr now takes about 8 s for everything, against about 11 s with Multiwfn.
- **3D view.** Atom style Line / Ball and stick ("CPK", as named in the project's viewer) / Space-filling (vdW), plus 1-based atom numbers, selectable on the Results tab.
- **Multiwfn's remaining role** is to produce the test references through the dev-only `make_reference.py` scripts. §3 still documents the Multiwfn workflow for users who want its extra quantities.

## 8. Sources

**Manual (supplied):** T. Lu, *Multiwfn Manual*, version 2026.7.10, sections cited inline. Section 3.25 asks users to cite T. Lu, Q. Chen, "Realization of Conceptual Density Functional Theory and Information-Theoretic Approach in Multiwfn Program", in *Conceptual Density Functional Theory*, Wiley-VCH (2022), pp. 631–647, [doi:10.1002/9783527829941.ch31](https://doi.org/10.1002/9783527829941.ch31). Multiwfn itself asks for T. Lu, F. Chen, *J. Comput. Chem.* 33, 580 (2012), and T. Lu, *J. Chem. Phys.* 161, 082503 (2024) [RUN: printed at start-up].

**Official site:** Multiwfn home page, <http://sobereva.com/multiwfn>. It was reachable only over plain HTTP in this session (the HTTPS handshake failed).

1. R. G. Parr, W. Yang, "Density functional approach to the frontier-electron theory of chemical reactivity", *J. Am. Chem. Soc.* 106, 4049 (1984), [doi:10.1021/ja00326a036](https://doi.org/10.1021/ja00326a036). Bibliographic record: [Duke Scholars](https://scholars.duke.edu/publication/649687).
2. C. Morell, A. Grand, A. Toro-Labbé, "New dual descriptor for chemical reactivity", *J. Phys. Chem. A* 109, 205 (2005), [doi:10.1021/jp046577a](https://doi.org/10.1021/jp046577a). Record: [repositorio UC](https://repositorio-az.uc.cl/dspace/items/1b25fe4a-c557-4364-8055-a34b5f1a80ef).
3. T. Lu (sobereva), "An overview of methods supported by Multiwfn for predicting reaction sites and reactivity of chemical systems" (in Chinese), 2026-04-01, <http://sobereva.com/767>. **Author of Multiwfn; read in full.**
4. W. Yang, W. J. Mortier, *J. Am. Chem. Soc.* 108, 5708 (1986), the condensed Fukui function. Record: [Duke FDS](https://fds.duke.edu/db/pratt/FIP/weitao.yang/publications/235145).
5. T. Lu (sobereva), "Using Multiwfn to easily calculate various quantities defined in conceptual density functional theory" (in Chinese), last update 2024-04-17, <http://sobereva.com/484>. **Author of Multiwfn; read in full.** Source of the charge/multiplicity rules for radicals and ions.
6. B. Wang, C. Rong, P. K. Chattaraj, S. Liu, "A comparative study to predict regioselectivity, electrophilicity and nucleophilicity with Fukui function and Hirshfeld charge", *Theor. Chem. Acc.* 138, 124 (2019), [doi:10.1007/s00214-019-2515-1](https://doi.org/10.1007/s00214-019-2515-1). Record: [repositorio U. Chile](https://repositorio.uchile.cl/handle/2250/154443) (abstract only).
7. P. Fuentealba, P. Pérez, R. Contreras, "On the condensed Fukui function", *J. Chem. Phys.* 113, 2544 (2000), [doi:10.1063/1.1305879](https://doi.org/10.1063/1.1305879). Record: [repositorio UC](https://repositorio-az.uc.cl/dspace/items/775cbc8c-f6ff-44cc-9c51-989b4dc703e0).
8. R. K. Roy et al., "On non-negativity of Fukui function indices" (1999) and part II (2000). Records: [BITS Pilani](https://dspace.bits-pilani.ac.in/items/5fa5351a-f87d-43c2-a72b-3a7740771cb9), [IAS repository](https://repository.ias.ac.in/25895).
9. "Comparative study on the methods for predicting the reactive site of nucleophilic reaction", *Sci. China Chem.* 58, 1845 (2015), [doi:10.1007/s11426-015-5494-7](https://doi.org/10.1007/s11426-015-5494-7). Cited by [M §4.A.4] and [WEB 3]; not opened.
10. T. Lu, F. Chen et al., *Acta Phys.-Chim. Sin.* 30, 628 (2014), [doi:10.3866/PKU.WHXB201401211](https://doi.org/10.3866/PKU.WHXB201401211). Cited by [M §4.A.4] and [WEB 3]; not opened.
11. P. W. Ayers et al., *Phys. Chem. Chem. Phys.* (2006), on the meaning of negative condensed Fukui functions: [RSC landing page](https://pubs.rsc.org/en/Content/ArticleLanding/2006/CP/B606167B). Abstract only.
12. "Fukui and dual-descriptor matrices within the framework of spin-polarized density functional theory", *Phys. Chem. Chem. Phys.* (2013): [RSC landing page](https://pubs.rsc.org/en/Content/ArticleLanding/2013/CP/C3CP50736J). Record: [CONICET](https://repositoriosdigitales.mincyt.gob.ar/vufind/Record/CONICETDig_25702d76b13239d1d0389397e0b0803e). Abstract only.
13. P. Geerlings, F. De Proft, W. Langenaeker, "Conceptual density functional theory", *Chem. Rev.* 103, 1793 (2003), [doi:10.1021/cr990029p](https://doi.org/10.1021/cr990029p) (DOI confirmed; general review); and Geerlings et al., "Conceptual DFT: status, prospects, issues", *Theor. Chem. Acc.* (2020), [HAL](https://hal.archives-ouvertes.fr/hal-02473087).

Where only an abstract or record was available ([WEB 6–8], [WEB 11–12]), only the claims stated there are relied on.
