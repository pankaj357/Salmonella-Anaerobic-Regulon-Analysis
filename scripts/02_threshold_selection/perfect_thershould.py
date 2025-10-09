import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_and_analyze_predictions(gold_standard_path, fimo_raw_path):
    """Load FIMO predictions and analyze using gene name mapping"""
    
    # Load files
    gold_standard = pd.read_csv(gold_standard_path)
    gold_standard.columns = gold_standard.columns.str.strip()
    
    fimo_raw = pd.read_csv(fimo_raw_path, sep='\t')
    fimo_raw.columns = fimo_raw.columns.str.strip()
    
    print("🔬 ANALYZING PREDICTION THRESHOLDS")
    print("="*50)
    print(f"Gold standard targets: {len(gold_standard)}")
    print(f"Raw FIMO hits: {len(fimo_raw)}")
    
    # Thresholds to test
    thresholds = [
        ('p < 1e-5', fimo_raw[fimo_raw['p-value'] < 1e-5]),
        ('p < 1e-4', fimo_raw[fimo_raw['p-value'] < 1e-4]),
        ('q < 0.05', fimo_raw[fimo_raw['q-value'] < 0.05]),
        ('q < 0.10', fimo_raw[fimo_raw['q-value'] < 0.10]),
        ('q < 0.20', fimo_raw[fimo_raw['q-value'] < 0.20])
    ]
    
    results = []
    
    # Create mapping dict: gene_name -> LT2_Locus_Tag
    gene_to_locus = dict(zip(gold_standard['Gene_Name'].str.lower(), gold_standard['LT2_Locus_Tag']))
    
    for threshold_name, filtered_data in thresholds:
        filtered_data = filtered_data.copy()
        
        # Map sequence_name to LT2_Locus_Tag
        def map_sequence(seq_name):
            seq_lower = str(seq_name).lower()
            for gene, locus in gene_to_locus.items():
                if gene in seq_lower:
                    return locus
            return None
        
        filtered_data['LT2_Locus_Tag'] = filtered_data['sequence_name'].apply(map_sequence)
        valid_predictions = filtered_data[filtered_data['LT2_Locus_Tag'].notna()]
        
        for regulator in gold_standard['Regulator'].unique():
            regulator_gs = gold_standard[gold_standard['Regulator'] == regulator]
            known_targets = set(regulator_gs['LT2_Locus_Tag'])
            predicted_targets = set(valid_predictions['LT2_Locus_Tag'])
            
            true_positives = known_targets.intersection(predicted_targets)
            recall = len(true_positives) / len(known_targets) if known_targets else 0
            precision = len(true_positives) / len(predicted_targets) if predicted_targets else 0
            
            results.append({
                'Threshold': threshold_name,
                'Regulator': regulator,
                'Known_Targets': len(known_targets),
                'True_Positives': len(true_positives),
                'Recall': recall,
                'Precision': precision,
                'Total_Predictions': len(predicted_targets),
                'Recovered_Genes': ', '.join(sorted(true_positives))
            })
    
    return pd.DataFrame(results)

def plot_threshold_comparison(results_df):
    """Plot recall and prediction counts across thresholds"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    pivot_recall = results_df.pivot_table(index='Threshold', columns='Regulator', values='Recall')
    pivot_recall.plot(kind='bar', ax=ax1, width=0.8)
    ax1.set_title('Recall by Threshold and Regulator')
    ax1.set_ylabel('Recall')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    
    total_preds = results_df.groupby('Threshold')['Total_Predictions'].mean()
    total_preds.plot(kind='bar', ax=ax2, color='lightblue', width=0.8)
    ax2.set_title('Average Predictions by Threshold')
    ax2.set_ylabel('Number of Predictions')
    ax2.tick_params(axis='x', rotation=45)
    
    for i, v in enumerate(total_preds):
        ax2.text(i, v + 0.5, f'{v:.0f}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('threshold_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

def recommend_best_threshold(results_df):
    """Select the threshold with best average recall and precision"""
    threshold_perf = results_df.groupby('Threshold').agg({'Recall':'mean', 'Precision':'mean'}).round(3)
    threshold_perf['Score'] = (threshold_perf['Recall'] + threshold_perf['Precision']) / 2
    best_threshold = threshold_perf['Score'].idxmax()
    print("\n🎯 Recommended threshold:", best_threshold)
    print(threshold_perf.loc[best_threshold])
    return best_threshold

def create_final_validation_with_pvalue(gold_standard_path, fimo_raw_path, pvalue_threshold=1e-4):
    """Final validation using best p-value threshold"""
    gold_standard = pd.read_csv(gold_standard_path)
    gold_standard.columns = gold_standard.columns.str.strip()

    fimo_raw = pd.read_csv(fimo_raw_path, sep='\t')
    fimo_raw.columns = fimo_raw.columns.str.strip()

    fimo_filtered = fimo_raw[fimo_raw['p-value'] < pvalue_threshold].copy()
    
    gene_to_locus = dict(zip(gold_standard['Gene_Name'].str.lower(), gold_standard['LT2_Locus_Tag']))
    
    def map_sequence(seq_name):
        seq_lower = str(seq_name).lower()
        for gene, locus in gene_to_locus.items():
            if gene in seq_lower:
                return locus
        return None
    
    fimo_filtered['LT2_Locus_Tag'] = fimo_filtered['sequence_name'].apply(map_sequence)
    
    results = []
    for regulator in gold_standard['Regulator'].unique():
        regulator_gs = gold_standard[gold_standard['Regulator'] == regulator]
        known_targets = set(regulator_gs['LT2_Locus_Tag'])
        predicted_targets = set(fimo_filtered['LT2_Locus_Tag'].dropna())
        
        true_positives = known_targets.intersection(predicted_targets)
        false_negatives = known_targets - predicted_targets
        recall = len(true_positives) / len(known_targets) if known_targets else 0
        precision = len(true_positives) / len(predicted_targets) if predicted_targets else 0
        
        results.append({
            'Regulator': regulator,
            'Known_Targets': len(known_targets),
            'True_Positives': len(true_positives),
            'Predicted_Regulon_Size': len(predicted_targets),
            'Recall': recall,
            'Precision': precision,
            'Recovered_Genes': ', '.join(sorted(true_positives)),
            'Missed_Genes': ', '.join(sorted(false_negatives))
        })
    
    final_df = pd.DataFrame(results)
    total_known = final_df['Known_Targets'].sum()
    total_recovered = final_df['True_Positives'].sum()
    overall_recall = total_recovered / total_known if total_known else 0
    print(f"\n✅ Overall Recall: {overall_recall:.3f} ({overall_recall*100:.1f}%)")
    
    return final_df

def plot_validation_results(validation_results):
    """Publication-ready figure"""
    plt.style.use('default')
    sns.set_palette("colorblind")
    plt.rcParams['font.size'] = 11
    plt.rcParams['font.family'] = 'Arial'
    
    regulators = sorted(validation_results['Regulator'].unique())
    validation_results = validation_results.set_index('Regulator').loc[regulators].reset_index()
    
    x_pos = np.arange(len(regulators))
    width = 0.35
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Recall & Precision
    bars1 = ax1.bar(x_pos - width/2, validation_results['Recall'], width, label='Recall', color='#2E86AB', alpha=0.8, edgecolor='black', linewidth=0.5)
    bars2 = ax1.bar(x_pos + width/2, validation_results['Precision'], width, label='Precision', color='#A23B72', alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Transcription Factor', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Performance Metric', fontsize=12, fontweight='bold')
    ax1.set_title('A. Validation Against Experimental Gold Standard', fontsize=14, fontweight='bold', pad=20)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(regulators, rotation=45, ha='right')
    ax1.legend()
    ax1.set_ylim(0, 1.2)
    ax1.grid(True, alpha=0.3, axis='y')
    for bar in bars1:
        ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.03, f'{bar.get_height():.2f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
    for bar in bars2:
        ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.03, f'{bar.get_height():.2f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
    
    # Plot 2: Recovered vs Missed
    recovered = validation_results['True_Positives']
    missed = validation_results['Known_Targets'] - validation_results['True_Positives']
    ax2.bar(regulators, recovered, label='Recovered', color='#4CAF50', alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.bar(regulators, missed, bottom=recovered, label='Missed', color='#F44336', alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Transcription Factor', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Targets', fontsize=12, fontweight='bold')
    ax2.set_title('B. Target Recovery Summary', fontsize=14, fontweight='bold', pad=20)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(regulators, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    for i, (rec, mis) in enumerate(zip(recovered, missed)):
        ax2.text(i, rec/2, f'{rec}', ha='center', va='center', fontweight='bold', color='white')
        if mis > 0:
            ax2.text(i, rec + mis/2, f'{mis}', ha='center', va='center', fontweight='bold', color='white')
    
    # Plot 3: Predicted regulon sizes
    ax3.bar(regulators, validation_results['Predicted_Regulon_Size'], color='#FF9800', alpha=0.8, edgecolor='black', linewidth=0.5)
    ax3.set_xlabel('Transcription Factor', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Predicted Regulon Size', fontsize=12, fontweight='bold')
    ax3.set_title('C. Predicted Regulon Sizes\n(Uniform p < 1e-4 threshold)', fontsize=14, fontweight='bold', pad=20)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(regulators, rotation=45, ha='right')
    ax3.grid(True, alpha=0.3, axis='y')
    for i, size in enumerate(validation_results['Predicted_Regulon_Size']):
        ax3.text(i, size + 0.5, f'{size}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 4: Overall metrics
    overall_metrics = ['Recall', 'Precision']
    overall_values = [validation_results['Recall'].mean(), validation_results['Precision'].mean()]
    colors = ['#2E86AB', '#A23B72']
    bars = ax4.bar(overall_metrics, overall_values, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax4.set_ylabel('Score', fontsize=12, fontweight='bold')
    ax4.set_title('D. Overall Performance Metrics', fontsize=14, fontweight='bold', pad=20)
    ax4.set_ylim(0, 1.0)
    ax4.grid(True, alpha=0.3, axis='y')
    for bar, value in zip(bars, overall_values):
        ax4.text(bar.get_x() + bar.get_width()/2., value + 0.02, f'{value:.3f}', ha='center', va='bottom', fontweight='bold')
    
    total_known = validation_results['Known_Targets'].sum()
    total_recovered = validation_results['True_Positives'].sum()
    overall_recall = total_recovered / total_known
    ax4.text(0.5, 0.8, f'Overall Recall: {overall_recall:.1%}\n({total_recovered}/{total_known} targets)',
             ha='center', va='center', transform=ax4.transAxes, fontsize=12,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7))
    
    plt.tight_layout()
    plt.savefig('Figure_4_Validation_Analysis.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('Figure_4_Validation_Analysis.pdf', bbox_inches='tight', facecolor='white')
    plt.show()


def main():
    GOLD_STANDARD_CSV = "gold.csv"
    FIMO_RAW_FILE = "fimo_raw.tsv"
    
    # Step 1: Analyze multiple thresholds
    results_df = load_and_analyze_predictions(GOLD_STANDARD_CSV, FIMO_RAW_FILE)
    
    # Step 2: Plot threshold comparison
    plot_threshold_comparison(results_df)
    
    # Step 3: Recommend best threshold
    best_threshold = recommend_best_threshold(results_df)
    if 'p <' in best_threshold:
        pval = float(best_threshold.split('p < ')[1])
    else:
        pval = 1e-4
    
    # Step 4: Create final validation table using best threshold
    final_results = create_final_validation_with_pvalue(GOLD_STANDARD_CSV, FIMO_RAW_FILE, pvalue_threshold=pval)
    final_results.to_csv(f"final_validation_pvalue_{pval}.csv", index=False)
    print(f"\n💾 Results saved to: final_validation_pvalue_{pval}.csv")
    
    # Step 5: Generate publication-ready figure
    plot_validation_results(final_results)


if __name__ == "__main__":
    main()
