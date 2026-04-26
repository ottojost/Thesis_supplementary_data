This repository contains the digital supplementary materials for my master's thesis


### 1. DNA_Data:
Contains the complete template sequences used throughout the experiments:
* **PCR_template.txt:** Contains the original target sequence in its plasmid form, alongside the specific forward and reverse primers used for the successful PCR amplification.
* **RPA_template.txt:** Contains the linear DNA sequence that resulted from the PCR amplification. This sequence was subsequently used as the template for all downstream RPA reactions.
* *Note on formatting:* Within these files, uppercase letters indicate the gBlock sequence representing the actual target template (HINTW gene). Lowercase letters denote the flanking vector DNA required to complete the plasmid backbone, which is not part of the template under investigation.

### 2. Scripts
Contains the custom computational workflows:
* **Primer_generating_pipeline.py:** A custom Python script developed to generate primer pairs suitable for RPA with similar forward and reverse primers, using the modified output of PrimedRPA as input.
