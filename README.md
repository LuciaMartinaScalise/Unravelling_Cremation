# Unravelling Cremation: A Multi-Analytical Approach to Bronze Age Cremated Human Remains

---

This repository accompanies **Lucia Martina Scalise’s PhD thesis in Biological Anthropology**, titled  
**_“Unravelling cremation: a multi-analytical approach to the study of Bronze Age cremated human remains”_**  
(**University of Cambridge, 2025**), supervised by **Dr. Emma Pomeroy**.

The repository contains data files, scripts, and outputs supporting **Chapters 4, 5, and 6** of the thesis.  
Materials related to **Chapter 3** are available separately on Zenodo:  
👉 [https://doi.org/10.5281/zenodo.15042290](https://doi.org/10.5281/zenodo.15042290)  
(as a scientific paper has already been formulated from this chapter and will be sumbitted for peer-review).

---

## 📘 Thesis Overview

The study integrates virtual imaging, physical excavation, osteological, taphonomic, and proteomic approaches to explore the **mortuary behaviour** and **biological profiles** of Bronze Age individuals with the aim to improve the quality and quantity of information that can be retrieved from the study of cremated remains.

**Temporal & cultural context:** Burials date to the **Middle–Recent Bronze Age (ca. 1450–1150 BCE)** and belong to the **Terramare culture (Northern Italy)**.

Analyses were conducted on:
- **12 urns from the necropolis of Vicofertile (Parma)**  
- **5 urns from the necropolis of Nogarole Rocca (Verona)**  

The research applies a **multi-proxy analytical framework** combining:

-  **Virtual analysis** of urns and their contents through **CT imaging**
-  **Physical micro-excavation** of cremated remains
-  **Assessment of biological profiles** (sex, age-at-death, pathology)
-  **Evaluation of heat-induced changes** (e.g. colour, fracture patterns, fragmentation)
-  **Spatial distribution** of bones inside urns to assess possible **intentionality**
-  **Proteomic analyses** (amelogenin & ZooMS) for **sex estimation and species identification**
-  **Statistical modelling** to integrate and test quantitative evidence

---

## 📂 Repository Structure

The repository follows a reproducible research structure to facilitate data access and code execution:

```bash
📦 Unravelling_Cremation
│
├── 📁 data/
│   ├── 📁 Chapter_4/
│   │   └── [Database_Vicofertile_Nogarole.csv] # Comprehensive database
│   │
│   ├── 📁 Chapter_5/
│   │   ├── [Database_Vicofertile_Nogarole.csv] # Comprehensive database
│   │   ├── [Phases_fract_stai_dim.csv]         # Subset with colour, fracture, staining, and dimension data
│   │   ├── [Database_Fragmentation.csv]        # Subset used to calculate Fragmentation Index
│   │   └── [phy_distribution.csv]              # Subset for bone distribution analysis inside urns
│   │
│   └── 📁 Chapter_6/
│       ├── 📁 Amelogenin_data/
│       │   ├── [pFind.protein]          # Peptide search results (pFind)
│       │   ├── [pFind.spectra]          # Spectra results (pFind)
│       │   └── README.md                # Link to raw amelogenin dataset on **Zenodo**
│       │
│       └── 📁 ZooMS_raw_data/
│           └── [2 txt files]            # Raw peptide fingerprint data from ZooMS analyses
│
├── 📁 figures/
│   └── 📁 Chapter_5/
│       ├── 📁 Distribution_remains/    # Figures for the "Distribution remains" paragraph
│       ├── 📁 Fragmentation_weight/    # Figures for the "Fragmentation & weight" paragraph
│       └── 📁 Thermal_alterations/     # Figures for the "Thermal alterations" paragraph
│
├── 📁 scripts/
│   └── 📁 Chapter_5/
│       ├── [Distribution_remains.R]   # R script for distribution analyses & figures
│       ├── [Fragmentation_weight.R]   # R script for fragmentation analyses & figures
│       └── [Thermal_alterations.R]    # R script for thermal alteration analyses & figures
│
└── 📄 README.md                       # Project overview and usage guide
```
---

### 💡 Notes

- All scripts use the `here()` package to maintain **relative paths** and ensure full **reproducibility**.  
- Figures and tables can be **reproduced directly** from the provided R scripts.  
- Raw files from amelogenin analysis are stored on **Zenodo** and linked in the README contained in the folder Data/Chapter6/Amelogenin_raw data.
- The order in which the analyses are performed and the figures generated in the R scripts corresponds to the order in which they are presented in the thesis.

## 🛠️ Software & Key Packages

The analysis was conducted in **R (≥ 4.3)** using the following core packages:

---
```
library(vcd)          # for visualising and analysing categorical data (e.g., mosaic plots)
library(tidyverse)    # for data manipulation and cleaning
library(ggplot2)      # for creating high-quality visualisations
library(here)         # for managing relative file paths and reproducibility
library(RColorBrewer) # for managing colour palettes in plots
```
---

## Chapters Overview

### **Chapter 3**
Presents a **virtual approach** to urned cremations using **CT-based qualitative and quantitative analysis** and **permutation testing**.  
The corresponding dataset and code are available on Zenodo:  
👉 [https://doi.org/10.5281/zenodo.15042290](https://doi.org/10.5281/zenodo.15042290)

### **Chapter 4**
Focuses on the **micro-excavation** and **osteological analysis** of cremated remains to assess the biological profile of the individuals

### **Chapter 5**
Explores **heat-induced changes** and **bone distribution** to infer **cremation practices** and **possible ritual intent**.  
Includes:
- Analysis of **colour variations**, **fracture morphology**, and calculation of the **fragmentation index**
- Application of **mosaic plots** to test associations among variables

### **Chapter 6**
Presents **proteomic analyses**:
- **Amelogenin-based sex estimation** from 11 burned tooth crowns  
- **ZooMS identification** on one burned and one unburned bone fragment  

---

## Acknowledgements

The author gratefully acknowledges:

- **Prof. Maria Pia Morigi** (University of Bologna) and her team for the collaboration on **CT scans**
- **Dr. Enrico Crema** (University of Cambridge) for his **guidance and statistical expertise**
- **Dr. Sara Silvestrini** (University of Bologna) for performing the **ZooMS analysis**
- **Dr. Miranda Evans** (University of Cambridge) for leading the **amelogenin analysis**
- **Dr. Federico Lugli** (University of Modena and Reggio Emilia) for the collaboration on the **amelogenin analysis**
- **Cambridge Trust, St John's College, and Department of Archaeology (University of Cambridge)** for the financial support

---

## 📜 Citation

If you use the data, code, or figures from this repository, please cite as:

> **Scalise, L. M. (2025). _Unravelling cremation: a multi-analytical approach to the study of Bronze Age cremated human remains_. PhD Thesis, University of Cambridge.**

---

## 📧 Contact

**Lucia Martina Scalise**  
PhD student in Biological Anthropology, St John's College, University of Cambridge  
📩 [lms217@cam.ac.uk](mailto:lms217@cam.ac.uk) or [luciamartina.scalise@gmail.com](mailto:luciamartina.scalise@gmail.com)

---


