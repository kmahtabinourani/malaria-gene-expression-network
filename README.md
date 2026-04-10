# 🧬 Malaria Gene Expression Network Analysis

Computational analysis of malaria single-cell gene expression data using gene regulatory network (GRN) modeling and systems biology approaches.

---

## 📌 Overview

This project presents a structured pipeline for analyzing malaria gene expression data derived from **single-cell RNA sequencing (scRNA-seq)** datasets from the **Malaria Cell Atlas (~5110 genes)**.

The goal is to identify:
- Gene regulatory structures
- Highly connected **hub genes**
- Functional **gene communities**
- Target-centered subnetworks for biological interpretation

---

## 🧠 Workflow Overview

The analysis follows a two-stage systems biology pipeline:

### 🔹 1. Gene Regulatory Network (GRN) Construction
- Preprocessing and filtering of gene expression data
- Normalization (CP10K + log transformation)
- Selection of highly variable genes
- Correlation-based network construction
- Identification of **hub genes**
- Detection of **gene communities (modules)**

### 🔹 2. Target-Centered Network Analysis
- Integration of structural (PDB-based) information
- Extraction of target genes
- Subnetwork construction around selected targets
- Identification of **neighboring genes**
- Highlighting potential **essential genes**

---

## 📂 Repository Structure
- code/ → MATLAB scripts for network analysis
- data/ → Input malaria gene expression dataset
- results/ → Generated figures and output tables

---


---

## 💻 Code Description

### 📄 01_malaria_grn_hub_and_community_analysis.m
This script performs:

- Gene filtering based on expression levels
- Data normalization (CP10K + log1p)
- Selection of highly variable genes
- Construction of a **correlation-based GRN**
- Identification of **hub genes (degree centrality)**
- Detection of **gene communities using hierarchical clustering**

---

### 📄 02_malaria_gene_network_and_essential_hub_analysis.m
This script performs:

- Structural processing of protein data (PDB files)
- Construction of protein interaction adjacency matrices
- Extraction of target genes
- Subnetwork analysis around target genes
- Identification of **directly connected neighboring genes**
- Highlighting of **essential genes within the network**

---

## 📊 Example Output

### 🧬 Gene Regulatory Network Visualization

Below is an example GRN constructed from malaria scRNA-seq data:

<br>

![Gene Regulatory Network](results/GRN%20malaria.png)

<br>

---

## 🔬 Data Source

- **Malaria Cell Atlas (scRNA-seq dataset)**
- ~5110 genes across single parasite transcriptomes

---

## ⚙️ Requirements

- MATLAB (recommended R2020a or later)
- Statistics and Machine Learning Toolbox

---

## 📈 Key Features

- Correlation-based GRN construction  
- Hub gene detection using network centrality  
- Community detection via hierarchical clustering  
- Integration of structural and expression-based analysis  
- Modular and reproducible pipeline  

---

## ✍️ Author

**Kiana Mahtabi Nourani**  
Master’s Student – Electrical & Electronic Engineering  
Eastern Mediterranean University  

---

## 🚀 Future Work

- Integration with pathway databases (Reactome / KEGG)
- Drug–gene interaction analysis
- Machine learning-based biomarker prediction
- Multi-omics data integration
