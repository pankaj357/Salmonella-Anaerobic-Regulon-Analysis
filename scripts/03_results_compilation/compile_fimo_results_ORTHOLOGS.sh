#!/bin/bash

# compile_fimo_results_ORTHOLOGS.sh
# Compiles FIMO hits and generates a regulatory binary matrix using ORTHOLOG GROUPS only
# Uses optimized p-value threshold of 1e-4

# Define regulators
regulators=("ArcA" "Fnr" "NarL" "NarP" "NsrR" "DcuR" "TtrR")

# Create output directory
mkdir -p compiled_results_orthologs

echo "=== Compiling significant FIMO hits (p < 1e-4) ==="
for reg in "${regulators[@]}"; do
    echo "Processing $reg..."
    folder_name=$(echo "$reg" | tr '[:upper:]' '[:lower:]')
    fimo_file="fimo_results/${folder_name}/fimo.tsv"
    
    if [[ -f "$fimo_file" ]]; then
        awk -F'\t' 'BEGIN {OFS="\t"} \
            NR==1 {print "Regulator",$0; next} \
            $8 < 0.0001 {print "'"$reg"'",$0}' "$fimo_file" > "compiled_results_orthologs/${reg}_hits.tsv"
        echo "  Saved to compiled_results_orthologs/${reg}_hits.tsv"
    else
        echo "  WARNING: $fimo_file not found. Skipping $reg."
    fi
done

echo "✅ All FIMO hits compiled with p < 1e-4 threshold."
echo "=== Generating ORTHOLOG regulatory_binary_matrix.csv ==="

python3 - << 'END_PYTHON'
import pandas as pd
import glob, os

compiled_dir = "compiled_results_orthologs"
regulator_files = glob.glob(os.path.join(compiled_dir, "*_hits.tsv"))

ortholog_sets = {}
all_orthologs = set()

def extract_ortholog_group(seq_name):
    parts = str(seq_name).split('|')
    if len(parts) >= 3:
        # Prefer gene symbol if it contains letters
        return parts[1] if any(c.isalpha() for c in parts[1]) else parts[0]
    elif len(parts) == 2:
        return parts[1] if any(c.isalpha() for c in parts[1]) else parts[0]
    else:
        return parts[0]

print("Loading FIMO results and grouping by orthologs...")
for file in regulator_files:
    reg = os.path.basename(file).split("_hits.tsv")[0]
    try:
        df = pd.read_csv(file, sep="\t")
        if 'sequence_name' not in df.columns:
            print(f"  WARNING: Column 'sequence_name' not found in {file}")
            continue

        df['ortholog_group'] = df['sequence_name'].dropna().apply(extract_ortholog_group)
        ortholog_groups = set(df['ortholog_group'].astype(str).unique())  # ensure strings
        ortholog_sets[reg] = ortholog_groups
        all_orthologs.update(ortholog_groups)

        sample_genes = df['sequence_name'].head(3).tolist()
        sample_orthologs = df['ortholog_group'].head(3).tolist()
        print(f"  {reg}: {len(ortholog_groups)} ortholog groups")
        print(f"    Example: {sample_genes[0]} → {sample_orthologs[0]}")

    except Exception as e:
        print(f"  ERROR reading {file}: {e}")

# Build binary matrix
all_orthologs = sorted(str(x) for x in all_orthologs)  # ensure all strings
matrix = pd.DataFrame(0, index=all_orthologs, columns=sorted(ortholog_sets.keys()))
for reg, orthologs in ortholog_sets.items():
    matrix.loc[list(orthologs), reg] = 1

matrix.index.name = "Ortholog_Group"
output_file = os.path.join(compiled_dir, "regulatory_binary_matrix_ORTHOLOGS_p1e-4.csv")
matrix.to_csv(output_file)

# Summary
print(f"\n📊 ORTHOLOG MATRIX SUMMARY:")
print(f"   Total ortholog groups: {len(matrix)}")
print(f"   Total regulatory interactions: {matrix.sum().sum()}")

print(f"\n🔍 INTERACTIONS PER REGULATOR:")
for regulator in sorted(matrix.columns):
    count = matrix[regulator].sum()
    percentage = (count / len(matrix)) * 100
    print(f"     {regulator}: {count} ortholog groups ({percentage:.1f}% of total)")

# Show multi-strain grouping examples
print(f"\n🔄 MULTI-STRAIN GENE GROUPING EXAMPLES:")
multi_strain_examples = []
for file in regulator_files[:2]:  # first 2 regulators
    df = pd.read_csv(file, sep="\t")
    df['ortholog_group'] = df['sequence_name'].dropna().apply(extract_ortholog_group)
    counts = df['ortholog_group'].value_counts()
    for gene, c in counts[counts>1].head(2).items():
        strains = df[df['ortholog_group']==gene]['sequence_name'].head(2).tolist()
        multi_strain_examples.append(f"    {gene}: appears in {c} strains (e.g., {strains[0].split('|')[0]})")

for example in multi_strain_examples[:4]:
    print(example)

print(f"\n✅ regulatory_binary_matrix_ORTHOLOGS_p1e-4.csv saved in {compiled_dir}")
print("   This matrix contains TRUE ortholog groups (no strain duplicates)")
END_PYTHON
