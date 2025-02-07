# ARHGAP5 Mutations in Parkinson's Disease: Whole Exome Sequencing and Molecular Dynamics Analysis

## Project Overview
This project aims to identify and analyze novel mutations associated with Parkinson's disease (PD) using whole exome sequencing (WES) data from the Parkinson's Progression Markers Initiative (PPMI). We focus on the **ARHGAP5** gene, which encodes a Rho GTPase Activating Protein, and investigate its potential role in PD pathogenesis through **bioinformatics analysis** and **molecular dynamics simulations**.

---

## Workflow

### 1. Data Acquisition and Preprocessing
- **Downloaded PPMI whole exome data:** `PPMI_SI_3051.raw.vcf`
- **Source:** [PPMI Exome Sequencing Methods](https://ida.loni.usc.edu/download/files/genetic/61e82b8a-d44d-4a34-b0d5-bf6ed82987ac/ppmi/PPMI_Methods_Exome_Sequencing_116_20150311.pdf)

### 2. Variant Filtering and Normalization
```bash
bcftools view PPMI_SI_3051.raw.vcf -o filtered.vcf -i 'QUAL>30 && INFO/DP>10 && FORMAT/DP>10'
bcftools norm -m -any -f hg19.fa filtered.vcf -o normalized.vcf
```

### 3. Variant Annotation
```bash
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20240917 humandb/
table_annovar.pl normalized.vcf humandb/ -buildver hg19 -out annotated -protocol refGene,clinvar_20240917,gnomad211_exome -operation g,f,f -nastring . -vcfinput
```

### 4. Mutation Analysis
- **Executed** `mutation.py` for initial mutation identification
- **Ran** `check.py` to compare with known PD-associated genes (no novel mutations found)
- **Executed** `Novel.py` to identify potential novel mutations


### 5. Molecular Dynamics Simulations
- **Executed** `wt_md_run.sh and mut_md_run.sh` with the following parameters:
  - 50,000 steps of **energy minimization**
  - **NVT equilibration** at **300K**
  - **Production run** (extendable as needed)
- **Analysis Focus:**
  - **Structural stability (RMSD)**
  - **Local flexibility (RMSF)**
  - **Protein aggregation (TANGO)**

---

## Dependencies
Ensure the following tools and libraries are installed:
- **bcftools** (Variant processing)
- **ANNOVAR** (Variant annotation)
- **Python** (for custom scripts: `mutation.py`, `check.py`, `Novel.py`)
- **GROMACS** (for molecular dynamics simulations)

## Usage
1. **Set up environment** and install required dependencies.
2. **Download** PPMI whole exome data as specified in the workflow.
3. **Run preprocessing, filtering, and annotation steps** using the provided commands.
4. **Execute** custom Python scripts (`mutation.py`, `check.py`, `Novel.py`) for mutation analysis.
5. **Perform molecular dynamics simulations** using `modified-md-script.sh`.
6. **Analyze the results** with a focus on **RMSD, TANGO, and RMSF**.

## Results
- **Identified novel stopgain mutations** in ARHGAP5 (**R480X and E481X**) that may be relevant to **Parkinson's disease pathogenesis**.
- **Molecular dynamics simulations** provided insights into the **structural and functional implications** of these mutations.

## Contact
For questions or collaborations, please contact **Ananaya Jain** at **ananayajain2024@gmail.com**.

## References
- **PPMI Website**
- **ANNOVAR Documentation**
- **GROMACS Manual**
