# Chapter 02: Protein Ingestion, Cleaning, and Receptor Preparation

---

### 1. What this part does
This chapter explains how experimental protein structures from the Protein Data Bank (PDB) are cleaned and transformed into computational docking receptors. Raw X-ray crystal structures contain missing side-chain atoms, crystallization salts, solvent waters, and lack hydrogen atoms; DockLab repairs these defects using PDBFixer and converts the receptor into the PDBQT format required by AutoDock Vina.

---

### 2. The biology or chemistry behind it
Proteins in physiological conditions are surrounded by water and ions at roughly pH 7.4. However, X-ray crystallography requires proteins to form ordered crystals, which is typically achieved using high concentrations of precipitating salts (such as ammonium sulfate $\text{(NH}_4)_2\text{SO}_4$, glycerol, polyethylene glycol (PEG), or sodium chloride). 

When crystallographers solve the electron density map:
1. **Missing Hydrogen Atoms**: X-rays scatter off electrons. Because a hydrogen atom has only a single electron, its scattering power is too weak to be resolved at standard crystallographic resolutions ($1.5\text{ \AA} \dots 2.5\text{ \AA}$). Consequently, standard PDB files contain zero hydrogen coordinates.
2. **Missing Heavy Atoms in Side Chains**: Highly flexible surface residues (such as lysine or arginine side chains) often fluctuate across multiple conformations, resulting in smeared electron density that crystallographers cannot model, leaving incomplete side chains.
3. **Protonation States at Physiological pH**: The electrostatic charge of amino acid side chains depends strictly on their $pK_a$ relative to the solution pH. At pH 7.4, aspartic acid ($\text{Asp}^-$) and glutamic acid ($\text{Glu}^-$) are negatively charged carboxylates; lysine ($\text{Lys}^+$) and arginine ($\text{Arg}^+$) are positively charged; and histidine can exist as neutral $\delta$-tautomer (HID), $\epsilon$-tautomer (HIE), or positively charged doubly-protonated (HIP).

Without adding hydrogens and assigning partial atomic charges, an AutoDock scoring function cannot calculate electrostatic or hydrogen-bonding energies.

---

### 3. The logic step by step
In DockLab, receptor preparation is handled by `backend/app/services/protein_service.py`:

1. **Structure Ingestion (`fetch_or_load_pdb`)**: The system fetches the `.pdb` coordinate file from RCSB PDB (`https://files.rcsb.org/download/{ID}.pdb`) or accepts a user upload, caching it locally.
2. **Structural Inspection (`inspect_structure`)**: Biopython parses the coordinate hierarchy. The service enumerates polypeptide chains (e.g. Chains A and B), identifies crystallographic waters, and filters out common crystallization additives (`COMMON_ADDITIVES = {'SO4', 'GOL', 'EDO', 'PEG', 'CL', 'NA', ...}`).
3. **Crystal Ligand Extraction (`extract_crystal_ligand`)**: Before any receptor modification occurs, the co-crystallized drug molecule (e.g., residue `MK1` in 1HSG) is extracted directly from the raw coordinates along with its `CONECT` records. Preserving these exact crystal coordinates is critical for redocking validation.
4. **Defect Repair with PDBFixer (`clean_and_prepare_receptor`)**:
   - `fixer.removeHeterogens(keepWater=False)`: Strips all additives and bulk waters (while structural waters are recorded for QC).
   - `fixer.findMissingAtoms()` and `fixer.addMissingAtoms()`: Rebuilds missing heavy atoms in existing residue side chains based on internal amino acid topology templates.
   - `fixer.addMissingHydrogens(7.4)`: Places hydrogens on all residues and assigns histidine tautomers corresponding to physiological pH 7.4.
5. **Rigid PDBQT Parameterization**: Open Babel (`obabel -ipdb ... -opdbqt -xr`) assigns Gasteiger partial atomic charges and AutoDock 4 atom types (such as `OA` for hydrogen bond acceptor oxygens, `SA` for sulfur, `A` for aromatic carbon, `HD` for polar hydrogen donor), producing the final rigid receptor grid input.

---

### 4. A worked example with real numbers from this project
When preparing **HIV-1 protease** (`examples/1hsg/1hsg.pdb`):

- **Raw Crystal Coordinates**:
  - Total observed atoms: **1,686 atoms**
  - Chains: **Chain A** (99 residues) and **Chain B** (99 residues)
  - Co-crystallized ligand: **`MK1`** (indinavir analog, 45 heavy atoms)
  - Crystallographic water molecules: **127 waters**
- **After PDBFixer Preparation (`1hsg_cleaned.pdb`)**:
  - All 127 waters and ligand `MK1` removed.
  - Incomplete heavy atoms checked: 0 missing heavy atoms in the core catalytic dyad.
  - Added explicit hydrogens at pH 7.4: Hydrogen count increased from 0 to **1,514 hydrogens**.
  - Total prepared receptor atoms: **3,028 atoms**.
- **Final PDBQT Receptor (`1hsg_cleaned.pdbqt`)**:
  - Output file size: **148,108 bytes**.
  - Successfully formatted with Gasteiger charges and AutoDock atom types.

---

### 5. A simple diagram
```
   [ Raw PDB File (1HSG.pdb) ]
          │
          ├──► Step 1: Save Crystal Ligand (MK1) ──► crystal_ligand.pdb
          │
          ▼
   [ PDBFixer Processing ]
   • Strip crystallization salts (SO4, GOL)
   • Strip bulk waters
   • Restore missing heavy atoms in side chains
   • Add explicit Hydrogens at pH 7.4
          │
          ▼
   [ 1hsg_cleaned.pdb ]
          │
          ▼ (Open Babel -xr / Meeko)
   [ 1hsg_cleaned.pdbqt ] (Rigid receptor with Gasteiger partial charges)
```

---

### 6. What can go wrong
1. **Modeler Adding Hydrogens at the Wrong pH**: Adding hydrogens at pH 2.0 would protonate aspartates to neutral carboxylic acids ($\text{Asp}-\text{COOH}$), completely destroying the catalytic dyad's negative charge and making indinavir unable to dock. DockLab locks protonation to pH 7.4.
2. **Missing Loops Near the Active Site**: Disordered peptide loops near the active site can leave structural voids. PDBFixer can model missing loops using geometry, but artificial loop modeling can introduce steric clashes. DockLab records missing gaps in `warnings` rather than guessing arbitrary loop coordinates.
3. **Loss of Co-Crystallized Ligand**: In poorly designed pipelines, cleaning removes the ligand before its coordinates are recorded, making redocking validation impossible. DockLab extracts the crystal ligand to disk as Step 1.

---

### 7. How to say it in an interview
> *"In DockLab, receptor preparation goes far beyond a simple file conversion. Raw PDB crystal structures lack hydrogen atoms, contain non-physiological crystallization additives, and may have missing heavy atoms in flexible side chains. We use Biopython to isolate the co-crystallized ligand in its true experimental coordinates, PDBFixer to restore missing heavy atoms and assign protonation states at pH 7.4, and Open Babel to parameterize Gasteiger charges into rigid PDBQT format, ensuring an accurate electrostatic active site."*

---

### 8. Check yourself
1. **Why do standard X-ray crystal structures in the PDB lack hydrogen atoms?**  
   *Answer*: Hydrogen atoms possess only a single electron and thus scatter X-rays too weakly to be detected at typical crystallographic resolutions ($>1.2\text{ \AA}$).
2. **Why is it essential to extract the crystal ligand before running PDBFixer?**  
   *Answer*: PDBFixer's `removeHeterogens()` strips all non-protein molecules; saving the ligand first preserves its exact experimental coordinates for redocking validation and search box centering.
3. **What is the protonation state of catalytic Asp25 in HIV-1 protease at neutral pH?**  
   *Answer*: In the catalytic dyad, one Asp25 residue is deprotonated (carboxylate anion, $\text{COO}^-$) and one is protonated ($\text{COOH}$), sharing a proton to coordinate the nucleophilic water molecule.
4. **What does the `-xr` flag in Open Babel signify during receptor conversion?**  
   *Answer*: It instructs Open Babel to treat the macromolecule as a **rigid receptor**, outputting atomic charges and atom types without creating internal flexible torsion trees.

---

### Primary sources to read
1. Eastman, P., et al. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Computational Biology*, 13(7), e1005659.  
   DOI: [10.1371/journal.pcbi.1005659](https://doi.org/10.1371/journal.pcbi.1005659)
2. Cock, P. J., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.  
   DOI: [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)
3. O'Boyle, N. M., et al. (2011). Open Babel: An open chemical toolbox. *Journal of Cheminformatics*, 3, 33.  
   DOI: [10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)

