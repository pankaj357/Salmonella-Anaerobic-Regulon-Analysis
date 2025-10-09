# Salmonella Anaerobic Regulon Analysis

## 📖 Overview

This repository contains the computational workflow for analyzing anaerobic transcriptional regulation in *Salmonella enterica*. The code implements a comprehensive pipeline for comparative regulon prediction, from promoter sequence extraction to network analysis and classification.

## 🏗️ Project Structure
```
Salmonella-Anaerobic-Regulon-Analysis/
├── data/
│   ├── processed/
│   │   └── training_genes/ # Training gene sets for motif discovery
│   └── results/            # Analysis outputs (see .gitignore)
├── scripts/
│   ├── 01_data_preparation/    # Promoter extraction and training data
│   ├── 03_results_compilation/ # FIMO results processing
│   ├── 04_network_analysis/    # Regulon analysis and classification
│   └── 05_visualization/       # Figure generation
├── figures/                    # Generated visualizations
├── requirements.txt
└── README.md
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

cd ../03_results_compilation  
./compile_fimo_results.sh

cd ../04_network_analysis
python analyze_regulators.py
python glob_loc.py

cd ../05_visualization
python visul.py
```

## 📁 Script Descriptions

### Data Preparation (`scripts/01_data_preparation/`)
- `extract_promoters.py`: Extract promoter sequences from genome annotations  
- `extract_training_sequences.sh`: Prepare training sets for motif discovery  

### Results Compilation (`scripts/03_results_compilation/`)
- `compile_fimo_results.sh`: Process FIMO binding site predictions  
- `regultor_sumry.py`: Compile and filter significant hits  

### Network Analysis (`scripts/04_network_analysis/`)
- `analyze_regulators.py`: Basic regulon characterization and matrix creation  
- `glob_loc.py`: Advanced network analysis and classification algorithms  

### Visualization (`scripts/05_visualization/`)
- `visul.py`: Generate publication-quality figures and plots  

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
