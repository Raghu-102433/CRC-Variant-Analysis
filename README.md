

# 🧬 Colorectal Cancer (CRC) Variant Analysis with WES

This project, titled **“CRC Variant Analysis with WES”** utilized whole-exome sequencing data from the NCBI BioProject **[PRJNA540552]** 
The analysis compared a **tumour sample (SRR8993342)** and its **matched normal (SRR8993344)** to identify both somatic and germline genomic variants contributing to colorectal carcinogenesis.  


---

## Somatic Variant Analysis  
> **Note:** Somatic variants were called using all three approaches — Tumor-only, Tumor vs Panel of Normals (PoN), and Tumor vs Normal. Among these, the **tumour vs normal** vcf was prioritized for interpretation, as it produced a realistic and high-confidence variant set by accurately removing germline and sample-specific artifacts, unlike the inflated variant counts observed in the other modes.

- **Key Mutated Genes:** Major somatic alterations were detected in **RECK**, **ADAMTS9**, **WNK2**, **SGK3**, and the large structural gene **SYNE1**, suggesting deregulation of *extracellular matrix remodeling, kinase signaling,* and *cytoskeletal organization* — hallmarks of **colorectal tumour progression and invasion**.  
- **Ontology / Pathway Summary:** Enrichr analysis (notably involving **SYNE1**, **MUSK**, **DPP9**, and **SLC27A4**) revealed enrichment in *muscle cell differentiation*, *neuromuscular synapse morphology*, and *protein autophosphorylation*, indicating **structural and signaling perturbations** influencing tumour cell adhesion and polarity rather than canonical APC–KRAS–TP53 pathways.  
- **Variant Classification (AMP):** One **Tier 3 Pathogenic (AMP)** nonsense mutation, *SYNE1 p.R31\** (*c.91C>T, Exon 4*), was identified with **high confidence**, while all other somatic variants were classified as **Tier 3 VUS (AMP)**.

---

## Germline Variant Analysis  

- **Shared Variants:** Comparison of tumour and matched normal identified 14 shared genes — *RBFOKL, CSMD1, RYR2, LRP1B, PKD1L2, WWOX, SORCS2, PCDH15, PCDH11, OBSCN, MUC16, MAP2K3, KALRN,* and *SYNE1* — representing a set of germline variants that may confer genomic background influences relevant to CRC susceptibility.  
- **Ontology Summary:** Functional enrichment (genes: **SYNE1, RYR2, MAP2K3, PCDH15**) connected these to *cardiac and skeletal muscle contraction*, *MAPK cascade*, and *calcium ion transport*. This reflects shared involvement in **cytoskeletal coordination and calcium-mediated signaling**, processes with known links to epithelial integrity and colorectal epithelial stress response.  
- **Variant Classification (ACMG):** All shared variants were classified as **VUS, Benign, or Likely Benign**, though predictive scores marked some as functionally **deleterious**:  
  - **MAP2K3:** *c.243G>A (p.Val81=, rs62057672)* and *c.281G>T (p.Arg94Leu, rs56067280)*  
  - **PCDH15:** *c.1304A>C (p.Asp435Ala, rs4935502)*  
  These variants likely represent *low-penetrance functional modifiers*, with possible roles in **MAPK-mediated stress signaling** and **cell adhesion maintenance**, both critical to colorectal epithelial homeostasis and barrier integrity.  

---

## Potential Actionables for the Case  

- **SYNE1 p.R31\*** (*c.91C>T, Exon 4*) was identified as a **pathogenic (PVS1–Very Strong)** null variant introducing a premature stop codon early in the transcript.  
  - **Mechanistic Implication:** The truncation is predicted to undergo **nonsense-mediated mRNA decay (NMD)**, leading to loss of SYNE1 protein function.  
  - **Clinical Context:** SYNE1 loss-of-function is a known disease mechanism, supported by >350 pathogenic null variants reported across 132 exons in ClinVar. The variant’s gnomAD observed/expected score (0.423) aligns with known LoF intolerance.
  - **Occurrence:** Variants in SYNE1 were detected in both somatic and germline analyseis pointing to a possible dual role in tumor biology and genetic background
  - **CRC Relevance:** Although SYNE1 pathogenic variants are primarily reported in neuromuscular and neurodegenerative contexts, the gene’s structural role in **nuclear-cytoskeletal coupling** suggests that truncating events could compromise nuclear integrity and mechanotransduction — processes recently linked to **CRC cell invasion and chromatin instability**.  
  - **Clinical Trial Reference:** Variants in SYNE1 are represented in non-CRC contexts (e.g., congenital muscle disease studies, NCT01403402; methotrexate clearance studies, NCT0219791), but **no direct therapeutic or prognostic evidence** is currently established in CRC.  

---

## How To Use:

System Requirements:
> All pipelines and scripts were tested with a **minimum of 16 GB RAM** and **12 CPU threads**. Users can adjust these settings based on their own system resources.

1️⃣ Environment Setup
Make sure Conda is installed, then create and activate the environment using provided yaml file
```bash
conda env create -f variant_env.yaml 
conda activate variant_env
```

2️⃣ Prepare Reference & Resource Files
All the scripts can be found in workflow_scripts/
Before running any pipeline, set up required references (genome, PoN, BQSR known sites, annotation resources, etc.)
```bash
./setup_requirements.sh
```

3️⃣ Run Pipelines (Variant Calling)

All variant-calling workflows are implemented as .sh scripts.
Use the relevant pipeline according to your use case.

```bash
# Tumor vs Normal 
./somatic_vavriantcalling_pipeline_tumour_vs_normal <tumour_id> <normal_id>
# → Detects true somatic variants by comparing tumor to its matched normal

# Tumor only
./somatic_variantcalling_pipeline_tumour_only.sh <tumour_id>
# → Used when no matched normal sample is available

# Tumor vs Panel of Normals (PoN)
./somatic_variantcalling_pipeline_tumour_vs_PoN.sh <tumour_id>
# → Filters technical artifacts using a pre-generated or custom PoN

# Germline Variant Calling
./germline_variantcalling_pipelne.sh <sample_id>
# → Identifies shared or inherited variants between paired samples

# Create Panel of Normals
./create_PoN.sh
#use the ouput PoN in somatic_variantcalling_pipeline_tumour_vs_PoN.sh
# → Aggregates multiple normal BAMs to create PoN resource for artifact filtering 
```
Each pipeline outputs an annotated VCF file for downstream analysis.

Vcf Annotation can also be done with VEP GUI tool for better analysis.

4️⃣ VCF Analysis & Visualization with vcfanalyser
```
chmod +x vcf_analyser.sh

./vcf_analyser.sh <input_vcf> <output_dir>
```

Example:

```
./vcfanalysis.sh sample_vcf.vcf analysisresults/
```

5️⃣ Downstream Functional & Clinical Interpretation
```bash
Annotated VCF
   ↓
Variant consequence and gene filtering
   ↓
Functional enrichment → Enrichr
   ↓
Clinical significance classification → Franklin by Genoox (query HGVS variant notations)
```

