# Chapter 03: Small-Molecule Ligand Preparation and Cheminformatics

---

### 1. What this part does
This chapter details how small-molecule drug candidates entered as 1D/2D SMILES strings are parsed, sanitized, desalted, and transformed into 3D energy-minimized conformations. It covers the physics of bond flexibility, how Meeko sets up rotational torsion trees in PDBQT format, and how drug-likeness descriptors (Lipinski, Veber, PAINS, Brenk) are calculated using RDKit.

---

### 2. The biology or chemistry behind it
Small molecules are flexible chemical entities whose biological activity depends on their ability to adopt a low-energy 3D shape that matches the electronic and steric contours of a protein's active site.

- **SMILES Limitations**: A SMILES string (e.g. `CC(=O)Oc1ccccc1C(=O)O`) is a 2D topological graph. It defines which atoms are connected to which, but contains zero information about 3D bond lengths, bond angles, or dihedral angles.
- **Conformational Search & Distance Geometry**: To generate a realistic 3D structure, computational chemists use Distance Geometry. RDKit implements the **ETKDGv3** algorithm (Experimental-Torsion Knowledge Distance Geometry), which combines experimental small-molecule crystal torsion preferences with distance bounds to construct chemically realistic ring conformations and stereocenters.
- **Torsional Degrees of Freedom**: Inside an active site, a drug molecule rarely binds in a rigid form. It rotates around single covalent bonds (rotatable bonds) while maintaining rigid planar aromatic rings and double bonds. AutoDock Vina models this flexibility by constructing a **torsion tree**: selecting a central "root" rigid cluster and defining rotating branches.

---

### 3. The logic step by step
In DockLab, ligand preparation is executed by `backend/app/services/ligand_service.py`:

1. **SMILES Parsing & Desalting (`parse_and_standardize_smiles`)**:
   - `Chem.MolFromSmiles(smiles)` validates chemical valency and atom types.
   - If the input contains multiple salt fragments (e.g., indinavir sulfate), `Chem.GetMolFrags` extracts the largest heavy-atom fragment to isolate the active pharmaceutical ingredient.
   - `Chem.SanitizeMol` verifies aromaticity, kekulization, and formal charges.
2. **Drug-Likeness & Structural Alert Screening (`compute_drug_likeness`)**:
   - Molecular Weight (MW), $\log P$ (Crippen partition coefficient), Hydrogen Bond Donors (HBD), and Hydrogen Bond Acceptors (HBA) are computed to evaluate **Lipinski's Rule of 5**.
   - Rotatable bonds and Polar Surface Area (TPSA) are evaluated against **Veber's Rules** for oral bioavailability.
   - **PAINS** (Pan-Assay Interference Compounds) and **Brenk** alert catalogs are screened using RDKit's `FilterCatalog` to flag reactive or problematic moieties.
   - `Chem.FindMolChiralCenters(mol, includeUnassigned=True)` detects undefined stereocenters.
3. **3D Embedding & Force Field Minimization (`prepare_ligand_from_smiles`)**:
   - Explicit hydrogens are added (`Chem.AddHs`).
   - `AllChem.EmbedMolecule(mol_3d, params)` runs ETKDGv3 using a pinned random seed for deterministic 3D conformer coordinates.
   - `AllChem.MMFFOptimizeMolecule(mol_3d)` minimizes steric strain using the Merck Molecular Force Field (MMFF94).
4. **Meeko PDBQT Parameterization**:
   - `MoleculePreparation()` detects rotatable bonds, macrocycles, and assigns Gasteiger partial charges.
   - `PDBQTWriterLegacy.write_string()` writes the AutoDock PDBQT format containing `ROOT`, `BRANCH`, and `ENDBRANCH` torsion records.

---

### 4. A worked example with real numbers from this project
Consider the approved protease inhibitor **darunavir**:
- **SMILES**: `CC(C)CN(C[C@H]([C@H](Cc1ccccc1)NC(=O)O[C@H]2CO[C@H]3[C@@H]2OCCO3)O)S(=O)(=O)c4ccc(cc4)N`
- **Computed Descriptors**:
  - Molecular Weight: **$547.67\text{ g/mol}$** (violates Lipinski rule $>500$)
  - $\log P$: **$2.68$** (passes $\le 5$)
  - H-Bond Donors: **$3$** (passes $\le 5$)
  - H-Bond Acceptors: **$8$** (passes $\le 10$)
  - Rotatable Bonds: **$11$** (violates Veber rule $>10$)
  - Topological Polar Surface Area (TPSA): **$140.4\text{ \AA}^2$**
  - Lipinski Violations: **$1$** (still classified as drug-like, as 1 violation is allowed)
  - Heavy Atom Count: **$38$ heavy atoms**
  - PAINS / Brenk Alerts: **$0$ alerts detected**
- **3D Conformer & PDBQT**:
  - Conformer embedded in **$0.08\text{ seconds}$**.
  - PDBQT contains **11 active torsions** and **75 total atoms** (including hydrogens).

---

### 5. A simple diagram
```
   [ 2D SMILES String ]
          │
          ▼ (RDKit Desalting & Sanitization)
   [ Desalted 2D Molecule ] ──► [ Drug-Likeness: Lipinski, Veber, PAINS ]
          │
          ▼ (ETKDGv3 3D Embedding + MMFF94 Minimization)
   [ Energy-Minimized 3D Conformer ]
          │
          ▼ (Meeko MoleculePreparation)
   [ PDBQT File with Torsion Tree ]
   ROOT
     ATOM ... (rigid core)
   BRANCH ... (rotatable bond)
     ATOM ... (flexible arm)
   ENDBRANCH
   TORSDOF 11
```

---

### 6. What can go wrong
1. **Unassigned Stereochemistry (Racemic Inputs)**: If a chemist enters a chiral drug with a flat bond (e.g. `CC(O)C` without `@` or `@@`), RDKit arbitrarily chooses $(R)$ or $(S)$. If the biologically active enantiomer is $(S)$, docking the random $(R)$ isomer will produce an inactive pose. DockLab flags unassigned stereocenters in the QC report.
2. **Embedding Highly Flexible Molecules**: Molecules with $>25$ rotatable bonds have huge conformational search spaces ($3^{25} \approx 8.4 \times 10^{11}$ possible states). Vina's Monte Carlo search will struggle to converge, leading to poor scores. DockLab emits a flexibility warning for $>20$ rotatable bonds.
3. **Improper Salt Separation**: If indinavir is docked with its counterion (e.g., indinavir $\cdot$ sulfate), Vina treats them as two disconnected molecules inside one box, causing docking to fail completely. DockLab automatically keeps the largest organic fragment.

---

### 7. How to say it in an interview
> *"In DockLab, ligand preparation transforms raw 2D SMILES into fully parameterized 3D docking candidates. We desalt the structure, compute Lipinski and Veber drug-likeness rules, screen for PAINS and Brenk structural alerts, and verify stereochemistry. We then generate realistic 3D coordinates using RDKit's ETKDGv3 distance geometry algorithm, minimize energy with the MMFF94 force field, and use Meeko to parameterize the flexible torsion tree into AutoDock PDBQT format."*

---

### 8. Check yourself
1. **What is the difference between a 2D SMILES string and a 3D conformer?**  
   *Answer*: A SMILES string is a topological graph specifying atom connectivity and valence. A 3D conformer defines the exact $(x, y, z)$ spatial Cartesian coordinates of every atom in space.
2. **Why does RDKit's ETKDGv3 algorithm produce better initial conformations than random distance geometry?**  
   *Answer*: ETKDGv3 incorporates experimental torsion angle preferences extracted from thousands of crystal structures in the Cambridge Structural Database (CSD), ensuring realistic ring conformations and rotamer states.
3. **What is a PAINS alert, and why is it important in drug discovery?**  
   *Answer*: PAINS (Pan-Assay Interference Compounds) are chemical motifs that produce false positives in biochemical assays through non-specific mechanisms like protein reactivity, fluorescence interference, or colloidal aggregation.
4. **Why do many approved HIV-1 protease inhibitors violate Lipinski's Rule of 5?**  
   *Answer*: Because they are peptidomimetics designed to span the large, symmetrical active site tunnel of a homodimeric enzyme, requiring higher molecular weights ($500 \dots 720\text{ Da}$) and more rotatable bonds ($>10$) than standard oral small molecules.

---

### Primary sources to read
1. Riniker, S., & Landrum, G. A. (2015). Better informed distance geometry: using experimental molecular geometries to improve the generation of 3D small-molecule conformations. *Journal of Chemical Information and Modeling*, 55(12), 2562-2574.  
   DOI: [10.1021/acs.jcim.5b00654](https://doi.org/10.1021/acs.jcim.5b00654)
2. Lipinski, C. A., et al. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Advanced Drug Delivery Reviews*, 46(1-3), 3-26.  
   DOI: [10.1016/s0169-409x(00)00129-0](https://doi.org/10.1016/s0169-409x(00)00129-0)
3. Baell, J. B., & Holloway, G. A. (2010). New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries and for their exclusion in promising hits. *Journal of Medicinal Chemistry*, 53(7), 2719-2740.  
   DOI: [10.1021/jm901137j](https://doi.org/10.1021/jm901137j)

