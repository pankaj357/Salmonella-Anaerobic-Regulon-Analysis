# Salmonella Anaerobic Regulon Analysis

## 📖 Overview

This repository contains the computational workflow for analyzing anaerobic transcriptional regulation in *Salmonella enterica*. The code implements a comprehensive pipeline for comparative regulon prediction, from promoter sequence extraction to network analysis and classification.

## 🏗️ Project Structure
```
Salmonella-Anaerobic-Regulon-Analysis/
├── data/
│   ├── processed/
│   │   ├── AEZ45052.1_promoters.fasta
│   │   ├── all_promoters.fasta
│   │   ├── CDU88905.1_promoters.fasta
│   │   ├── EBX3951613.1_promoters.fasta
│   │   ├── EDA9785606.1_promoters.fasta
│   │   ├── EGT0473696_promoters.fasta
│   │   ├── ETA88628.1_promoters.fasta
│   │   └── training_genes/
│   │       ├── arcA_training_genes.txt
│   │       ├── dcuR_training_genes.txt
│   │       ├── fnr_training_genes.txt
│   │       ├── narL_training_genes.txt
│   │       ├── narP_training_genes.txt
│   │       ├── nsrR_training_genes.txt
│   │       └── ttrR_training_genes.txt
│   ├── raw/
│   │   ├── all_genomes/
│   │   │   ├── NC_003197.2/
│   │   │   │   ├── EGT0473696_promoters.csv
│   │   │   │   ├── EGT0473696_promoters.fasta
│   │   │   │   ├── GCF_000006945.2_ASM694v2_genomic.fna
│   │   │   │   └── genomic.gff
│   │   │   ├── NC_016832.1/
│   │   │   │   ├── AEZ45052.1_promoters.csv
│   │   │   │   ├── AEZ45052.1_promoters.fasta
│   │   │   │   └── ncbi_dataset-14/
│   │   │   │       └── [...]
│   │   │   └── [...]
│   │   └── reference_genome/
│   │       ├── GCF_000006945.2_ASM694v2_genomic.fna
│   │       └── genomic.gff
│   ├── results/
│   │   ├── classification_pvalue_1e-4/
│   │   │   ├── jaccard_similarity_matrix.csv
│   │   │   ├── overlap_count_matrix.csv
│   │   │   └── regulator_classification_comprehensive.csv
│   │   ├── classification_qvalue_0.05/
│   │   │   ├── jaccard_similarity_matrix.csv
│   │   │   ├── overlap_count_matrix.csv
│   │   │   └── regulator_classification_comprehensive.csv
│   │   ├── compiled_results_pvalue_1e-4/
│   │   │   ├── ArcA_hits.tsv
│   │   │   ├── DcuR_hits.tsv
│   │   │   ├── Fnr_hits.tsv
│   │   │   └── [...]
│   │   ├── compiled_results_qvalue_0.05/
│   │   │   ├── ArcA_hits.tsv
│   │   │   ├── DcuR_hits.tsv
│   │   │   └── [...]
│   │   ├── fimo_results/
│   │   │   ├── arcA/
│   │   │   │   ├── fimo.tsv
│   │   │   │   └── [...]
│   │   │   ├── dcuR/
│   │   │   │   └── [...]
│   │   │   └── [...]
│   │   ├── meme_outputs/
│   │   │   ├── arcA/
│   │   │   │   ├── meme.html
│   │   │   │   └── [...]
│   │   │   └── [...]
│   │   ├── motifs/
│   │   │   ├── ArcA.txt
│   │   │   ├── DcuR.txt
│   │   │   └── [...]
│   │   └── threshould_selection_result/
│   │       └── final_validation_pvalue_0.0001.csv
│   └── threshould_selection_data/
│       ├── fimo_raw.tsv
│       └── gold.csv
├── figures/
│   ├── pvalue_0.01/
│   │   ├── jaccard_heatmap.png
│   │   ├── overlap_heatmap.png
│   │   └── [...]
│   ├── qval_0.05/
│   │   └── [...]
│   └── Threshold_validation/
│       └── [...]
├── scripts/
│   └── [...]
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
- **Classification Algorithm**: percentile ranks of the Globalness Score 

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

### Threshold_Selection (`scripts/02_results_compilation/`)
- `perfect_thershould.py`: Selcect optimimum threshold for further classification 

### Results Compilation (`scripts/03_results_compilation/`)
- `compile_fimo_results.sh`: Process and Filter FIMO binding site predictions    

### Network Analysis and Visualization (`scripts/04_analaysis_and_visualization/`)
- `glob_loc.py`: Advanced network analysis, classification algorithms and visualization 


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
