# Chapter 14: HIV-1 Protease and the Structural Mechanism of Drug Resistance

---

### 1. What this part does
This chapter explores HIV-1 protease, the retroviral enzyme central to the maturation of infectious HIV virions and one of modern medicine's greatest drug design triumphs. It explains how the enzyme functions as an aspartic protease homodimer, why approved protease inhibitors successfully block viral replication, and how specific amino acid mutations in the enzyme cause clinical drug resistance.

---

### 2. The biology or chemistry behind it
Human Immunodeficiency Virus Type 1 (HIV-1) produces its structural proteins and essential enzymes as long polyprotein precursor chains: **Gag** (group-specific antigen) and **Gag-Pol**. To convert an immature, non-infectious viral particle into a fully infectious, virulent particle, the viral **HIV-1 protease** must cleave these polyprotein precursors at nine specific peptide cleavage sites.

Structurally, HIV-1 protease is a **C2-symmetric homodimer**:
- It is composed of two identical polypeptide chains of 99 amino acids each (total molecular weight $\approx 21.5\text{ kDa}$).
- The active site lies at the dimer interface, where two catalytic aspartic acid residues—**Asp25** from Chain A and **Asp25** from Chain B—coordinate a catalytic water molecule to hydrolyze peptide scissile bonds.
- Two flexible peptide loops called **flaps** (residues 45–55 on each monomer) cover the active site. In the open conformation, substrates or inhibitors enter; when bound, the flaps close tightly down over the ligand, locking it inside a hydrophobic tunnel.

Because HIV-1 reverse transcriptase lacks proofreading activity, the virus replicates with an error rate of $\approx 3 \times 10^{-5}$ mutations per nucleotide per cycle. Under the selective pressure of protease inhibitor therapy, mutant variants emerge. Single amino acid substitutions alter the shape, volume, or electrostatics of the binding cleft, decreasing inhibitor affinity while retaining sufficient enzymatic activity to cleave viral polyproteins.

---

### 3. The logic step by step
In DockLab's case study (`docs/case_study/`, `benchmark/run_case_study.py`), we model these resistance mechanisms structurally:

1. **Receptor Variant Selection (`app/services/protein_service.py`)**:
   - The wild-type structure is loaded from crystal coordinates (e.g. `1HSG`, resolution $2.00\text{ \AA}$).
   - Clinically validated mutant crystal structures are retrieved from RCSB PDB (or modeled via verified in-silico mutation):
     - **V82A**: Mutation of valine to alanine at residue 82 removes two methyl groups ($-\text{CH}_3$), creating an active site cavity that weakens hydrophobic van der Waals contacts with bulky inhibitor groups (e.g., indinavir, ritonavir).
     - **I50V**: Isoleucine to valine mutation located in the flap tips; alters flap flexibility and contact geometry with darunavir and amprenavir.
     - **I84V**: Located in the substrate-binding pocket; alters steric packing against aromatic inhibitor moieties.
     - **L90M**: Leucine to methionine mutation located outside the immediate catalytic pocket (secondary/allosteric mutation) that alters dimer interface packing and indirectly reshapes the active site.
     - **D30N**: Aspartate to asparagine mutation that selectively destroys a crucial ionic/hydrogen-bonding interaction with nelfinavir.
2. **Ligand Docking Across Variants (`app/services/docking_service.py`)**:
   - Approved protease inhibitors (darunavir, saquinavir, indinavir, ritonavir, nelfinavir, lopinavir, atazanavir) are docked into both wild-type and mutant variants.
3. **Differential Metric Calculation (`app/services/analysis_service.py`)**:
   - $\Delta\text{Score} = \text{Score}_{\text{mutant}} - \text{Score}_{\text{WT}}$ (a positive value indicates weaker predicted binding).
   - Lost Contacts: specific residue contacts present in wild-type but absent in the mutant are tabulated.
   - Key-Residue Recovery: conservation of catalytic contacts (Asp25A, Asp25B, Gly27, Ile50) is quantified.

---

### 4. A worked example with real numbers from this project
Consider the impact of the **V82A** mutation on **indinavir**:

- **Wild-Type (`1HSG`)**:
  - Residue Val82 side chain consists of an isopropyl group ($-\text{CH(CH}_3)_2$).
  - In `1HSG`, Val82 forms 6 close hydrophobic contacts ($< 4.0\text{ \AA}$) with indinavir's benzyl and pyridyl rings.
  - Wild-type Vina docking score: **$-11.4\text{ kcal/mol}$**.
- **Mutant (`1QBR` or V82A model)**:
  - Ala82 has only a single methyl group ($-\text{CH}_3$), leaving a void of $\approx 35\text{ \AA}^3$ in the pocket.
  - Hydrophobic contacts at position 82 drop to 2.
  - Vina docking score: **$-10.3\text{ kcal/mol}$** ($\Delta\text{Score} = +1.1\text{ kcal/mol}$).
- **Biological Correlation**: In the Stanford HIV Drug Resistance Database, the V82A mutation produces a clinical resistance fold-change of $\approx 4\times \dots 10\times$ resistance to indinavir, consistent with the loss of van der Waals packing detected computationally.

---

### 5. A simple diagram
```
      [ Flap A (Ile50) ]                 [ Flap B (Ile50) ]
              \                                 /
               \                               /
       ──────────────────────────────────────────────
                [ Inhibitor (e.g. Darunavir) ]
       ──────────────────────────────────────────────
               /                               \
         Val82 (WT: bulky)            Val82 (WT: bulky)
        [Ala82 (Mut: void)]          [Ala82 (Mut: void)]
               \                               /
                [ Catalytic Dyad: Asp25A / Asp25B ]
```

---

### 6. What can go wrong
1. **Homodimer Asymmetry**: In real crystals, the two identical monomers are labeled Chain A and Chain B. A mutation might appear in one or both chains depending on crystallographic symmetry. DockLab verifies chain labeling and biological assemblies (REMARK 350).
2. **Missing Flap Density in Crystals**: When unliganded (apo), the flaps are highly dynamic and often exhibit poor electron density. Crystallographers sometimes omit these residues from the PDB file. DockLab's QC engine flags missing flap residues (residues 45–55) as critical failures before docking.
3. **Overpredicting Resistance from Static Structures**: Protein backbones relax in solution. A static docking grid cannot capture backbone conformational adaptation unless multi-conformer or induced-fit methods are used. DockLab explicitly notes this physical limitation.

---

### 7. How to say it in an interview
> *"HIV-1 protease is an essential C2-symmetric homodimeric enzyme required for viral particle maturation. Standard therapy uses competitive active-site inhibitors that mimic the tetrahedral transition state of peptide cleavage. However, clinical mutations such as V82A, I50V, and L90M alter hydrophobic packing and flap dynamics, reducing inhibitor binding affinity. In DockLab, we model both wild-type and mutant variants, correlating the computational loss of specific atomic contacts with published clinical resistance data from the Stanford HIV Drug Resistance Database."*

---

### 8. Check yourself
1. **What is the biological consequence if HIV-1 protease is completely inhibited by a drug?**  
   *Answer*: The viral polyproteins Gag and Gag-Pol cannot be cleaved, preventing maturation; the newly budded viral particles remain immature, non-infectious, and unable to infect new host T-cells.
2. **Which two catalytic residues form the enzymatic active site in the HIV-1 protease homodimer?**  
   *Answer*: Aspartate 25 from Chain A (Asp25A) and Aspartate 25 from Chain B (Asp25B).
3. **How does the V82A mutation physically cause resistance to indinavir?**  
   *Answer*: Mutating Valine (with two methyl side-chain branches) to Alanine (with only one) creates a steric cavity that reduces hydrophobic van der Waals contacts with the inhibitor, decreasing binding affinity.
4. **Why is the D30N mutation uniquely associated with resistance to nelfinavir rather than other protease inhibitors?**  
   *Answer*: Nelfinavir relies on a specific hydrogen bond / salt-bridge interaction with the carboxylate group of Asp30; mutating Asp30 to neutral Asparagine disrupts this specific interaction without severely crippling cleavage of natural viral substrates.

---

### Primary sources to read
1. Wlodawer, A., & Erickson, J. W. (1993). Protease inhibitors: An overview of the design of inhibitors of the aspartic protease of HIV-1. *Annual Review of Biochemistry*, 62(1), 543-585.  
   DOI: [10.1146/annurev.bi.62.070193.002551](https://doi.org/10.1146/annurev.bi.62.070193.002551)
2. Rhee, S. Y., et al. (2003). Human immunodeficiency virus reverse transcriptase and protease sequence database. *Nucleic Acids Research*, 31(1), 298-303.  
   Website & Database: [Stanford HIV Drug Resistance Database](https://hivdb.stanford.edu/)
3. Ali, A., et al. (2010). Molecular basis for drug resistance in HIV-1 protease. *Viruses*, 2(11), 2509-2535.  
   DOI: [10.3390/v2112509](https://doi.org/10.3390/v2112509)

