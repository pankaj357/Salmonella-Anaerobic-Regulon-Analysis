import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

def create_flowchart_pipeline():
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    # Colors
    colors = {
        'data': '#FF9F1C',
        'analysis': '#2EC4B6', 
        'discovery': '#E71D36',
        'validation': '#6A4C93',
        'results': '#1982C4'
    }
    
    # Define steps with consistent alignment
    steps = [
        # Phase 1: Data Preparation
        {"name": "Strain Selection &\nPromoter Extraction", 
         "pos": (1, 6), "size": (2, 1), "color": colors['data'],
         "details": "6 S. enterica strains\n300bp promoter regions"},
        
        {"name": "Gold Standard\nCompilation", 
         "pos": (1, 4), "size": (2, 1), "color": colors['data'],
         "details": "17 validated targets\nLiterature curation"},
         
        # Phase 2: Motif Discovery  
        {"name": "De Novo Motif\nDiscovery", 
         "pos": (4, 6), "size": (2, 1), "color": colors['discovery'],
         "details": "MEME suite\nTraining set guided"},
         
        {"name": "Binding Site\nPrediction", 
         "pos": (4, 4), "size": (2, 1), "color": colors['discovery'],
         "details": "FIMO scanning\np < 1e-5 initial"},
         
        # Phase 3: Threshold Optimization
        {"name": "Threshold\nOptimization", 
         "pos": (7, 5), "size": (2, 1), "color": colors['validation'],
         "details": "Gold standard validation\np < 1e-4 selected"},
         
        # Phase 4: Network Analysis
        {"name": "Regulon\nCharacterization", 
         "pos": (10, 6), "size": (2, 1), "color": colors['analysis'],
         "details": "Regulon size calculation\nTarget identification"},
         
        {"name": "Network\nConnectivity", 
         "pos": (10, 4), "size": (2, 1), "color": colors['analysis'],
         "details": "Jaccard similarity\nCo-regulation analysis"},
         
        # Phase 5: Classification & Results
        {"name": "Globalness Score\nClassification", 
         "pos": (13, 6), "size": (2, 1), "color": colors['results'],
         "details": "Rank-based hierarchy\n3-tier classification"},
         
        {"name": "Biological\nInterpretation", 
         "pos": (13, 4), "size": (2, 1), "color": colors['results'],
         "details": "Hierarchical organization\nNetwork insights"}
    ]
    
    # Draw steps
    for step in steps:
        x, y = step["pos"]
        width, height = step["size"]
        
        # Create rounded rectangle
        box = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0.2",
                           facecolor=step["color"], alpha=0.9,
                           edgecolor='black', linewidth=1.8)
        ax.add_patch(box)
        
        # Add text
        ax.text(x + width/2, y + height/2 + 0.2, step["name"], 
                ha='center', va='center', fontweight='bold', fontsize=10)
        ax.text(x + width/2, y + height/2 - 0.3, step["details"],
                ha='center', va='center', fontsize=8)
    
    # Draw arrows (start from box centers)
    arrows = [
        ((3, 6.5), (4, 6.5)), ((3, 4.5), (4, 4.5)),  # Data → Discovery
        ((6, 6.5), (7, 5.5)), ((6, 4.5), (7, 5.5)),  # Discovery → Validation  
        ((9, 5.5), (10, 6.5)), ((9, 5.5), (10, 4.5)), # Validation → Analysis
        ((12, 6.5), (13, 6.5)), ((12, 4.5), (13, 4.5)) # Analysis → Results
    ]
    
    for start, end in arrows:
        ax.annotate('', xy=end, xytext=start,
                   arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # Add phase labels
    phases = [
        ("Data\nPreparation", 2, 7.5),
        ("Motif\nDiscovery", 5, 7.5), 
        ("Threshold\nOptimization", 8, 7.5),
        ("Network\nAnalysis", 11, 7.5),
        ("Classification &\nResults", 14, 7.5)
    ]
    
    for phase, x, y in phases:
        ax.text(x, y, phase, ha='center', va='center', 
                fontweight='bold', fontsize=11, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray'))
    
    # Layout
    ax.set_xlim(0, 16)
    ax.set_ylim(3, 8.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Comprehensive Pipeline for Anaerobic Regulon Analysis', 
                 fontsize=16, fontweight='bold', pad=30)
    
    plt.tight_layout()
    plt.savefig('pipeline_overview_detailed.png', dpi=300, bbox_inches='tight')
    plt.show()

create_flowchart_pipeline()
