# analyze_regulators_enhanced.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import warnings
import os

# Set style for better plots
plt.style.use('default')
sns.set_palette("husl")

# ----------------------------------------------------------------------------
# 1. Load binary regulatory matrix with error handling
# ----------------------------------------------------------------------------
def load_data(file_path):
    """Load and validate the regulatory binary matrix"""
    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Data file not found: {file_path}")
        
        binary_matrix = pd.read_csv(file_path, index_col=0)
        
        # Validate data
        if binary_matrix.empty:
            raise ValueError("Loaded matrix is empty")
        
        if not all(binary_matrix.isin([0, 1]).all()):
            warnings.warn("Matrix contains values other than 0 and 1. Converting to binary.")
            binary_matrix = binary_matrix.astype(bool).astype(int)
        
        print(f"✓ Successfully loaded binary matrix: {binary_matrix.shape}")
        print(f"✓ Regulators: {list(binary_matrix.columns)}")
        print(f"✓ Number of genes: {len(binary_matrix)}")
        
        return binary_matrix
        
    except Exception as e:
        print(f"✗ Error loading data: {e}")
        return None

# ----------------------------------------------------------------------------
# 2. Enhanced analysis functions
# ----------------------------------------------------------------------------
def calculate_regulon_sizes(binary_matrix):
    """Calculate regulon sizes and basic statistics"""
    regulon_sizes = binary_matrix.sum(axis=0)
    size_stats = {
        'mean': regulon_sizes.mean(),
        'std': regulon_sizes.std(),
        'min': regulon_sizes.min(),
        'max': regulon_sizes.max()
    }
    
    print("\n" + "="*50)
    print("REGULON SIZE ANALYSIS")
    print("="*50)
    print(regulon_sizes.sort_values(ascending=False).to_string())
    print(f"\nStatistics: Mean={size_stats['mean']:.1f}, Std={size_stats['std']:.1f}, "
          f"Range=[{size_stats['min']}-{size_stats['max']}]")
    
    return regulon_sizes, size_stats

def calculate_network_connectivity(binary_matrix):
    """Calculate Jaccard similarity-based connectivity with enhanced metrics"""
    regulator_names = binary_matrix.columns.tolist()
    n_regulators = len(regulator_names)
    
    # Initialize matrices
    jaccard_matrix = np.zeros((n_regulators, n_regulators))
    overlap_matrix = np.zeros((n_regulators, n_regulators))
    
    for i in range(n_regulators):
        for j in range(n_regulators):
            if i != j:
                A = binary_matrix.iloc[:, i]
                B = binary_matrix.iloc[:, j]
                
                intersection = np.logical_and(A, B).sum()
                union = np.logical_or(A, B).sum()
                
                # Jaccard similarity
                if union > 0:
                    jaccard_matrix[i, j] = intersection / union
                
                # Raw overlap count
                overlap_matrix[i, j] = intersection
    
    # Create dataframes
    jaccard_df = pd.DataFrame(jaccard_matrix, 
                             index=regulator_names, 
                             columns=regulator_names)
    
    overlap_df = pd.DataFrame(overlap_matrix,
                            index=regulator_names,
                            columns=regulator_names)
    
    # Calculate centrality metrics
    connectivity_centrality = jaccard_df.sum(axis=1)
    degree_centrality = (binary_matrix > 0).sum(axis=0)  # Same as regulon sizes
    betweenness = calculate_simple_betweenness(jaccard_df)
    
    print("\n" + "="*50)
    print("NETWORK CONNECTIVITY ANALYSIS")
    print("="*50)
    print("Jaccard-based connectivity:")
    print(connectivity_centrality.sort_values(ascending=False).to_string())
    
    return jaccard_df, overlap_df, connectivity_centrality, betweenness

def calculate_simple_betweenness(connectivity_df, threshold=0.01):
    """Simple betweenness centrality approximation"""
    n = len(connectivity_df)
    betweenness = pd.Series(0.0, index=connectivity_df.index)
    
    binary_connections = connectivity_df > threshold
    
    for i in range(n):
        for j in range(n):
            if i != j and binary_connections.iloc[i, j]:
                betweenness.iloc[i] += 1
                betweenness.iloc[j] += 1
    
    return betweenness / betweenness.max() if betweenness.max() > 0 else betweenness

def classify_regulators(metrics_df):
    """Enhanced classification with multiple metrics"""
    metrics_df['Size_Normalized'] = metrics_df['Regulon_Size'] / metrics_df['Regulon_Size'].max()
    metrics_df['Connectivity_Normalized'] = metrics_df['Network_Connectivity'] / metrics_df['Network_Connectivity'].max()
    
    metrics_df['Globalness_Score'] = (metrics_df['Size_Normalized'] + metrics_df['Connectivity_Normalized']) / 2
    
    percentiles = metrics_df['Globalness_Score'].rank(pct=True)
    
    def classify_regulator(pct_rank):
        if pct_rank >= 0.67:
            return 'Global'
        elif pct_rank <= 0.33:
            return 'Local'
        else:
            return 'Intermediate'
    
    metrics_df['Classification'] = percentiles.apply(classify_regulator)
    metrics_df['Percentile_Rank'] = percentiles
    metrics_df['Strength'] = pd.cut(metrics_df['Globalness_Score'], bins=[0,0.33,0.67,1.0],
                                   labels=['Weak','Moderate','Strong'])
    
    return metrics_df

# ----------------------------------------------------------------------------
# 3. Visualization functions (all separate figures)
# ----------------------------------------------------------------------------
from adjustText import adjust_text

def create_scatter_plot(metrics_df, save_path=None):
    """Scatter plot of normalized size vs connectivity with automatic label adjustment"""
    plt.figure(figsize=(8,6))
    color_map = {'Global': '#E74C3C', 'Intermediate': '#F39C12', 'Local': '#3498DB'}
    size_map = {'Global': 200, 'Intermediate': 150, 'Local': 100}
    
    texts = []
    
    # Plot points
    for classification in color_map.keys():
        mask = metrics_df['Classification'] == classification
        plt.scatter(metrics_df.loc[mask, 'Size_Normalized'], 
                    metrics_df.loc[mask, 'Connectivity_Normalized'],
                    s=size_map[classification], alpha=0.8, 
                    color=color_map[classification], label=classification,
                    edgecolor='white', linewidth=1)
        
        # Add labels to adjust later
        for i, row in metrics_df.loc[mask].iterrows():
            texts.append(plt.text(row['Size_Normalized'], row['Connectivity_Normalized'],
                                  row['Regulator'], fontsize=10, fontweight='bold'))

    # Automatically adjust labels to avoid overlaps
    adjust_text(texts, 
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
                expand_points=(1.2,1.2), expand_text=(1.2,1.2), force_text=0.5)
    
    plt.xlabel('Normalized Regulon Size', fontsize=12, fontweight='bold')
    plt.ylabel('Normalized Network Connectivity', fontsize=12, fontweight='bold')
    plt.title('Regulator Classification: Size vs Connectivity', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(alpha=0.3)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved scatter plot: {save_path}")
    plt.show()


def create_bar_chart(metrics_df, save_path=None):
    """Bar chart of Globalness Score"""
    plt.figure(figsize=(8,6))
    color_map = {'Global': '#E74C3C', 'Intermediate': '#F39C12', 'Local': '#3498DB'}
    colors = [color_map[cls] for cls in metrics_df['Classification']]
    y_pos = np.arange(len(metrics_df))
    
    plt.barh(y_pos, metrics_df['Globalness_Score'], color=colors, alpha=0.8)
    plt.yticks(y_pos, metrics_df['Regulator'], fontweight='bold')
    plt.xlabel('Globalness Score', fontsize=12, fontweight='bold')
    plt.title('Regulator Globalness Ranking', fontsize=14, fontweight='bold')
    plt.grid(axis='x', alpha=0.3)
    
    for i, v in enumerate(metrics_df['Globalness_Score']):
        plt.text(v + 0.01, i, f'{v:.3f}', va='center', fontweight='bold')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved bar chart: {save_path}")
    plt.show()

def create_regulon_size_barplot(metrics_df, save_path=None):
    """Publication-ready bar plot of regulon sizes colored by classification"""
    plt.figure(figsize=(10,6))
    sns.set_style("whitegrid")

    # Sort by size
    sorted_df = metrics_df.sort_values("Regulon_Size", ascending=False)

    # Define color map (consistent with scatter plot)
    color_map = {'Global': '#E74C3C', 'Intermediate': '#F39C12', 'Local': '#3498DB'}
    bar_colors = [color_map[cls] for cls in sorted_df["Classification"]]

    ax = sns.barplot(x=sorted_df["Regulator"], 
                     y=sorted_df["Regulon_Size"],
                     palette=bar_colors, alpha=0.9)

    # Annotate values on bars
    for i, v in enumerate(sorted_df["Regulon_Size"]):
        ax.text(i, v + max(sorted_df["Regulon_Size"])*0.01, str(v),
                ha='center', va='bottom', fontweight='bold', fontsize=9)

    plt.xlabel("Regulator", fontsize=12, fontweight='bold')
    plt.ylabel("Regulon Size (Number of Target Genes)", fontsize=12, fontweight='bold')
    plt.title("Regulon Sizes of Transcriptional Regulators", fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10)

    # Add legend
    handles = [plt.Rectangle((0,0),1,1, color=color_map[cls]) for cls in color_map]
    plt.legend(handles, color_map.keys(), title="Classification", fontsize=10)

    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved regulon size bar plot: {save_path}")
    plt.show()


def create_jaccard_heatmap(jaccard_df, save_path=None):
    plt.figure(figsize=(8,6))
    sns.heatmap(jaccard_df, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                square=True, cbar_kws={'label': 'Jaccard Similarity'})
    plt.title('Regulator Co-regulation Network (Jaccard Similarity)', fontweight='bold')
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved Jaccard heatmap: {save_path}")
    plt.show()

def create_overlap_heatmap(overlap_df, save_path=None):
    plt.figure(figsize=(8,6))
    sns.heatmap(overlap_df, annot=True, fmt='.0f', cmap='YlOrRd', square=True,
                cbar_kws={'label': 'Shared Target Count'})
    plt.title('Regulator Target Overlap (Shared Gene Count)', fontweight='bold')
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved overlap heatmap: {save_path}")
    plt.show()

def create_statistical_summary(metrics_df, regulon_sizes, connectivity_centrality):
    print("\n" + "="*50)
    print("STATISTICAL SUMMARY")
    print("="*50)
    
    corr_coef, p_value = pearsonr(regulon_sizes, connectivity_centrality)
    print(f"Correlation between regulon size and connectivity: {corr_coef:.3f} (p={p_value:.3f})")
    
    class_summary = metrics_df['Classification'].value_counts()
    print(f"\nClassification distribution:")
    for cls, count in class_summary.items():
        regulators = metrics_df[metrics_df['Classification'] == cls]['Regulator'].tolist()
        print(f"  {cls}: {count} regulators - {regulators}")
    
    top_global = metrics_df[metrics_df['Globalness_Score'] > 0.5]
    if not top_global.empty:
        print(f"\nTop global regulators (Score > 0.5):")
        for _, row in top_global.iterrows():
            print(f"  {row['Regulator']}: Score={row['Globalness_Score']:.3f}")

# ----------------------------------------------------------------------------
# 4. Main execution
# ----------------------------------------------------------------------------
def main():
    print("🧬 ENHANCED REGULATOR CLASSIFICATION ANALYSIS")
    print("="*60)
    
    binary_matrix = load_data('/Users/admin/Desktop/Salmonella-Anaerobic-Regulon-Analysis/data/results/compiled_results_pvalue_1e-4/regulatory_binary_matrix_p1e-4.csv')
    if binary_matrix is None:
        return
    
    os.makedirs('enhanced_results', exist_ok=True)
    
    regulon_sizes, size_stats = calculate_regulon_sizes(binary_matrix)
    jaccard_df, overlap_df, connectivity_centrality, betweenness = calculate_network_connectivity(binary_matrix)
    
    metrics_df = pd.DataFrame({
        'Regulator': binary_matrix.columns.tolist(),
        'Regulon_Size': regulon_sizes.values,
        'Network_Connectivity': connectivity_centrality.values,
        'Betweenness_Centrality': betweenness.values
    })
    
    metrics_df = classify_regulators(metrics_df)
    metrics_df = metrics_df.sort_values('Globalness_Score', ascending=False).reset_index(drop=True)
    
    print("\n" + "="*50)
    print("FINAL CLASSIFICATION RESULTS")
    print("="*50)
    display_columns = ['Regulator', 'Regulon_Size', 'Network_Connectivity', 
                      'Globalness_Score', 'Percentile_Rank', 'Classification']
    print(metrics_df[display_columns].round(4).to_string(index=False))
    
    # ------------------------------
    # Create visualizations
    # ------------------------------
    create_scatter_plot(metrics_df, 'enhanced_results/regulator_scatter.png')
    create_bar_chart(metrics_df, 'enhanced_results/regulator_bar.png')
    create_regulon_size_barplot(metrics_df, 'enhanced_results/regulon_size_barplot.png')  # NEW PLOT
    create_jaccard_heatmap(jaccard_df, 'enhanced_results/jaccard_heatmap.png')
    create_overlap_heatmap(overlap_df, 'enhanced_results/overlap_heatmap.png')
    
    create_statistical_summary(metrics_df, regulon_sizes, connectivity_centrality)
    
    output_files = {
        'classification_results': 'enhanced_results/regulator_classification_comprehensive.csv',
        'jaccard_matrix': 'enhanced_results/jaccard_similarity_matrix.csv',
        'overlap_matrix': 'enhanced_results/overlap_count_matrix.csv'
    }
    
    metrics_df.to_csv(output_files['classification_results'], index=False)
    jaccard_df.to_csv(output_files['jaccard_matrix'])
    overlap_df.to_csv(output_files['overlap_matrix'])
    
    print(f"\n✓ Analysis complete! Files saved:")
    for desc, path in output_files.items():
        print(f"  - {desc}: {path}")

if __name__ == "__main__":
    main()
