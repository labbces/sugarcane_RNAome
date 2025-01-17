# Table Documentation

This document provides a detailed description of the columns present in the table sugarcane-panRNAome.tar.gz.

**md5sum:** c88e8b7cd3326f5253afc6cfff832de7  sugarcane-panRNAome.tar.gz

It encapsulates the main results derived from the master's thesis "Multi-genotype analyses of long-ncRNA in Sugarcane", authored by Felipe Vaz Peres. 

The data presented in this table highlights key findings from the study, focusing on the identification, classification, and functional characterization of non-coding RNAs and long non-coding RNAs (lncRNAs) across 48 sugarcane genotypes present in the Sugarcane Pan-Transcriptome.

## Columns Description

### 1. **Pan Category**
   - **Description**: This column categorizes the gene based on its presence across different species in the pan-transcriptome. It helps in understanding the conservation and distribution of the gene across the sugarcane pan-transcriptome. 
   - **Possible Values**: Hard-core, Soft-core, Accessory, Exclusive.

### 2. **Gene**
   - **Description**: This column lists the gene identifier. In the context of this project, a gene is defined as a cluster of transcripts grouped using **MMSeqs2** with the following parameters: `-s 5.7 --cov-mode 2 --cluster-mode 2 -c 0.8 --min-seq-id 0.8`.
   - **Possible Values**: Gene IDs (e.g., `OG0000000`).

### 3. **Transcript**
   - **Description**: This column lists the transcript identifier.
   - **Possible Values**: Transcript IDs (e.g., `B1_k25_TRINITY_DN12555_c1_g1_i10`).

### 4. **Gene Category**
   - **Description**: This column categorizes the gene based on the functional role of its grouped transcripts. If a gene consists exclusively of long non-coding RNAs, it is categorized as **lncRNA**. If the gene contains both protein-coding and non-coding transcripts, it is categorized as **protein and non-coding**. 
   - **Possible Values**: lncRNA, ncRNA, protein-coding, protein and non-coding.

### 5. **Transcript Size**
   - **Description**: This column indicates the length of the transcript in nucleotides (nt).
   - **Possible Values**: Integer values representing the number of nucleotides (e.g., `500`).

### 6. **Transcript Category**
   - **Description**: This column categorizes the transcript according to its potential to be translated into functional proteins or not. Transcripts classidied as **protein and non-coding** have conflicting evidence when using software tools to classify them as either coding or non-coding. As a result, this category was added to accomodate transcripts with dual evidence. 
   - **Possible Values**: lncRNA, ncRNA, protein-coding, protein and non-coding.

### 7. **Transcript Rfam Family**
   - **Description**: This column specifies the Rfam family to which the transcript belongs. 
   - **Possible Values**: Rfam family names or IDs (e.g., `RF00645`, `RF00201`).

### 8. **Transcript GO (protein)**
   - **Description**: This column lists the Gene Ontology (GO) terms associated with the protein product of the transcript. GO terms describe the biological process.
   - **Possible Values**: GO terms (e.g., `GO:0009741`).

### 9. **Gene GO**
   - **Description**: This column lists the Gene Ontology (GO) terms associated with all transcripts within the gene. Each term linked to the transcripts is comma separated. 
   - **Possible Values**: GO terms (e.g., `GO:0090227,GO:0048519,GO:0090229,GO:0080113,GO:0048523,GO:0009741,GO:2000030`).

### 10. **GO (Guilt by association)**
   - **Description**: This column lists GO terms inferred for the lncRNAs and ncRNAs genes (and also their transcripts) based on "guilt by association". This means that the gene is assigned GO terms based on its co-expression with protein-coding genes that have known functions. Each term linked to the gene is comma separated.
   - **Possible Values**: GO terms (e.g., `GO:0009733,GO:0009725,GO:0009719,GO:0010033,GO:0042221`).

## Example Row

| Pan Category | Gene     | Transcript       | Gene Category | Transcript Size | Transcript Category | Transcript Rfam Family | Transcript GO (protein) | Gene GO                | GO (Guilt by association) |
|--------------|----------|------------------|---------------|-----------------|---------------------|-------------------------|-------------------------|-------------------------|---------------------------|
| Accessory         | OG185435  | QBYN04-26041_k25_TRINITY_DN4395_c0_g1_i1  | lncRNA | 610.0            | lncRNA                | RF00206                |    |   | GO:0006952,GO:0050896,GO:0006950,GO:0016145,GO:0019759,GO:0019762,GO:0006388,GO:0050898,GO:0080028,GO:0000394,GO:0010043,GO:0009820,GO:0016560,GO:0044273,GO:0006605,GO:0016143,GO:0019757,GO:0019760,GO:1901658,GO:0006625,GO:0016558,GO:0072662,GO:0072663,GO:0008150,GO:0072594,GO:0015919,GO:0043574,GO:0033365,GO:0045039,GO:0007007,GO:0051204,GO:0090151,GO:0007006,GO:0044743,GO:0007031,GO:0010038,GO:0051205,GO:1901136,GO:0043207,GO:0051707,GO:0006626,GO:0065002,GO:0044419,GO:0009607,GO:0070585,GO:0072655,GO:0071806    |

## Notes

If you use this data, please cite: Aa
