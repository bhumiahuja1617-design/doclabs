# DockLab Scientific & Computational Glossary

This glossary defines computational chemistry, structural biology, and software terms used across DockLab. Every term introduced in the learning guides is archived here for interview preparation and conceptual clarity.

---

### Active Site (Binding Pocket)
A specialized cleft or cavity on the surface or interior of a protein lined with amino acid side chains designed to bind substrates, cofactors, or inhibitor molecules.

### Affinity ($\Delta G_{\text{bind}}$ / $K_d$ / $K_i$)
The thermodynamic strength of the non-covalent association between a protein and a ligand. Formally expressed as the change in Gibbs free energy:
$$\Delta G = -RT \ln(K_a) = RT \ln(K_d)$$
where a more negative $\Delta G$ indicates tighter, more spontaneous binding.

### AutoDock Vina
A widely used open-source molecular docking program developed by Dr. Oleg Trott at the Scripps Research Institute. It uses an iterated local search algorithm (Monte Carlo combined with BFGS quasi-Newton minimization) and an empirical scoring function.

### Automorphism (Symmetry-Aware RMSD)
A graph permutation that maps a chemical molecule onto itself while preserving all atom types and bonds. Essential for computing true RMSD when a ligand has symmetric functional groups (such as rotating a phenyl ring or swapping equivalent carboxylate oxygens).

### Decoy
A chemical molecule in virtual screening benchmarking that has similar physical properties (molecular weight, charge, rotatable bonds) to known active compounds, but is assumed or confirmed to have no biological binding activity. Used to calculate ROC-AUC and Enrichment Factor.

### Desolvation
The energy cost required to strip water molecules away from a protein's active site and from a ligand's polar groups so they can form direct protein-ligand contacts.

### Enrichment Factor (EF)
A virtual screening performance metric measuring how many times more actives are found in the top $X\%$ (e.g., EF1%, EF5%) of a ranked list compared to a random selection:
$$\text{EF}_{X\%} = \frac{\text{Actives in top } X\% / \text{Total molecules in top } X\%}{\text{Total actives in dataset} / \text{Total molecules in dataset}}$$

### ETKDG (Experimental-Torsion Knowledge Distance Geometry)
An algorithm in RDKit that generates realistic, low-energy 3D conformations of small molecules using distance geometry combined with torsional angle preferences extracted from small-molecule crystal databases (CSD).

### Exhaustiveness
In AutoDock Vina, an integer parameter controlling the depth of the global search. Roughly corresponds to the number of independent Monte Carlo search trajectories initialized. Default in DockLab is 8 (Balanced); 4 for Fast; 32 for Thorough.

### Fpocket
An open-source C program that detects pockets and cavities on protein surfaces using Voronoi tessellation and alpha spheres without requiring prior knowledge of the binding site.

### Heavy Atom
Any atom in a molecule that is not hydrogen (e.g., Carbon, Nitrogen, Oxygen, Sulfur, Phosphorus). Most docking RMSD calculations and ligand efficiency metrics consider heavy atoms exclusively.

### HETATM (Heteroatom)
A record type in a Protein Data Bank (PDB) file designating non-standard residues, such as small-molecule drugs, cofactors, metal ions, or solvent water molecules, as opposed to standard amino acid peptide chains (`ATOM` records).

### Homodimer
A protein quaternary complex composed of two identical polypeptide subunits associated non-covalently. For example, HIV-1 protease is a C2-symmetric homodimer formed by two identical 99-amino-acid chains.

### Key-Residue Recovery
A chemical validation metric defined in DockLab: the fraction of specific residue contacts observed in the experimental crystal structure that are successfully re-established by the top docked pose:
$$\text{Recovery} = \frac{|\text{Contacts}_{\text{crystal}} \cap \text{Contacts}_{\text{docked}}|}{|\text{Contacts}_{\text{crystal}}|}$$

### Lipinski's Rule of Five (Ro5)
An empirical guideline formulated by Christopher Lipinski stating that poor oral absorption is more likely when a molecule has:
- Molecular weight $> 500\text{ Da}$
- $\log P > 5$
- Hydrogen bond donors $> 5$
- Hydrogen bond acceptors $> 10$

### Meeko
The modern Python preparation tool developed by the Forli Lab at Scripps to convert RDKit molecules and proteins into AutoDock PDBQT format, properly parameterizing flexible bonds, partial charges, and macrocycles.

### PDBFixer
A Python utility built on OpenMM that repairs structural defects in experimental PDB files, such as missing heavy side-chain atoms, unmodeled loops, non-standard residues, and missing hydrogens.

### PDBQT Format
An extension of the PDB coordinate format used by AutoDock and AutoDock Vina. Adds partial atomic **Q** (charges) and AutoDock atom **T** (types), along with `ROOT`, `BRANCH`, and `ENDBRANCH` records defining the rotational flexibility tree.

### ProLIF (Protein-Ligand Interaction Fingerprints)
A Python library that analyzes non-covalent interactions (hydrogen bonds, hydrophobic contacts, $\pi$-stacking, ionic salt bridges) across molecular dynamics trajectories and docking poses.

### Receiver Operating Characteristic (ROC) & AUC
A plot of the True Positive Rate (Sensitivity) versus False Positive Rate ($1 - \text{Specificity}$) across all docking score ranking thresholds. The Area Under the Curve (AUC) ranges from 0.5 (random guessing) to 1.0 (ideal discrimination between actives and decoys).

### Redocking Validation
A gold-standard control experiment where a ligand co-crystallized with a protein is extracted, stripped of its coordinates, computationally redocked into the binding site, and evaluated against the original crystal pose. An RMSD $< 2.0\text{ \AA}$ is considered a successful redock.

### Reliability Check
DockLab's signature feature: re-docking top-ranked ligands across multiple independent random seeds ($N=5 \dots 10$) to measure score standard deviation ($\sigma$), pose clustering consistency ($>80\%$), and ranking stability, isolating stable poses from stochastic noise.

### RMSD (Root-Mean-Square Deviation)
A measure of the average distance between the atoms of two superimposed conformations:
$$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \|\mathbf{r}_{i, \text{docked}} - \mathbf{r}_{i, \text{crystal}}\|^2}$$

### Scissile Bond
The specific covalent peptide bond in a substrate protein that is targeted and cleaved by a protease enzyme.

### SMILES (Simplified Molecular Input Line Entry System)
A concise line notation for encoding molecular structures and stereochemistry using short ASCII strings (e.g., `CC(=O)Oc1ccccc1C(=O)O` for aspirin).

### Spyrmsd
A fast Python tool that calculates root-mean-square deviations between molecular conformations using graph theory to identify atomic symmetries and automorphisms.

