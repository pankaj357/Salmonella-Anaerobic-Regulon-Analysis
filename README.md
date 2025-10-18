# Salmonella Anaerobic Regulon Analysis

## 📖 Overview

This repository contains the computational workflow for analyzing anaerobic transcriptional regulation in *Salmonella enterica*. The code implements a comprehensive pipeline for comparative regulon prediction, from promoter sequence extraction to network analysis and classification.

## 🏗️ Project Structure
```
Salmonella-Anaerobic-Regulon-Analysis/
├── data/
│   ├── processed/
│   │   ├── *_promoters.fasta
│   │   └── training_genes/
│   │       └── *_training_genes.txt
│   ├── raw/
│   │   ├── all_genomes/
│   │   │   └── NC_*/ 
│   │   │       └── *_promoters.fasta/.csv, genomic.fna, genomic.gff
│   │   └── reference_genome/
│   │       ├── *.fna
│   │       └── *.gff
│   ├── results/
│   │   ├── classification_pvalue_1e-4/
│   │   │   └── *.csv
│   │   ├── classification_qvalue_0.05/
│   │   │   └── *.csv
│   │   ├── compiled_results_*/ 
│   │   │   └── *_hits.tsv
│   │   ├── fimo_results/
│   │   │   └── */fimo.tsv
│   │   ├── meme_outputs/
│   │   │   └── */meme.html
│   │   ├── motifs/
│   │   │   └── *.txt
│   │   └── threshould_selection_result/
│   │       └── final_validation_*.csv
│   └── threshould_selection_data/
│       ├── fimo_raw.tsv
│       └── gold.csv
├── figures/
│   ├── pvalue_0.01/
│   │   └── *.png
│   ├── qval_0.05/
│   │   └── *.png
│   └── Threshold_validation/
│       └── *.png
├── scripts/
│   └── *.py
├── README.md
└── requirements.txt

```

## 🔬 Methodology

### 1. Data Preparation
- **Promoter Extraction**: Scripts to extract promoter sequences from genomic data  
- **Training Sets**: Curated gene lists for motif discovery  
- **Sequence Processing**: Multi-strain promoter consolidation  

### 2. Motif Discovery
- **MEME Suite**: De novo motif discovery implementation  
- **Parameters**: Standard bioinformatics workflow for transcription factor motif identification  
- **Quality Control**: Statistical validation of discovered motifs  

### 3. Binding Site Prediction  
- **FIMO Analysis**: Genome-wide transcription factor binding site prediction  
- **Statistical Filtering**: Multiple testing correction and significance thresholds  
- **Cross-strain Integration**: Aggregation of predictions across multiple genomes  

### 4. Network Analysis
- **Regulon Characterization**: Computational identification of target genes  
- **Connectivity Metrics**: Network analysis using similarity measures  
- **Classification Algorithm**: Unsupervised learning for regulator categorization  

### 5. Visualization
- **Publication Figures**: Scripts for generating research visualizations  
- **Network Plots**: Regulatory network representations  
- **Statistical Graphics**: Data exploration and result presentation  

## 🚀 Installation & Usage

### Requirements
```bash
pip install -r requirements.txt
```

**Dependencies**
```
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.5.0
seaborn>=0.11.0
scikit-learn>=1.0.0
biopython>=1.79
adjustText>=0.8
```

### Basic Usage
```bash
# Run complete analysis pipeline
./scripts/run_analysis.sh

# Or run individual steps
cd scripts/01_data_preparation
python extract_promoters.py

cd ../02_threshold_selection
python perfect_thershould.py

cd ../03_results_compilation  
./compile_fimo_results.sh

cd ../04_analaysis_and_visualization
python glob_loc.py
```

## 📁 Script Descriptions

### Data Preparation (`scripts/01_data_preparation/`)
- `extract_promoters.py`: Extract promoter sequences from genome annotations  
- `extract_training_sequences.sh`: Prepare training sets for motif discovery  

### Threshold_Selection (`scripts/02_threshold_selection/`
- `perfect_thershould.py`: Script to find optimized threshold

### Results Compilation (`scripts/03_results_compilation/`)
- `compile_fimo_results.sh`: Process FIMO binding site predictions and filter significant hits 

### Network Analysis and Visualization (`scripts/04_analaysis_and_visualization/`)
 - `glob_loc.py`: Advanced network analysis and classification algorithms and plots generation

## 🔒 Data Availability
- Processed results are available in the `data/results/` directory  
- Training gene lists are provided for reproducibility  

## 🤝 Contributing
This repository is maintained for transparency and reproducibility. For questions about the methodology or code implementation, please open an issue.

## 📄 License
This project is licensed under the MIT License.

## 🙏 Acknowledgments
- MEME Suite for motif discovery tools  
- Biopython for bioinformatics utilities  
- Scikit-learn for machine learning algorithms  
