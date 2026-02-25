import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import glob
import os
import argparse

# Force non-interactive backend
matplotlib.use('Agg')

def get_clean_sample_name(filename):
    return filename.split('_all_results.tsv')[0]

def process_results(input_dir):
    all_data = []
    result_files = glob.glob(os.path.join(input_dir, "*_all_results.tsv"))
    
    if not result_files:
        print(f"Error: No '*_all_results.tsv' files found in {input_dir}")
        return None, []

    all_sample_ids = sorted([get_clean_sample_name(os.path.basename(f)) for f in result_files])
    print(f"Found {len(all_sample_ids)} total genomes.")

    for f in result_files:
        sample_id = get_clean_sample_name(os.path.basename(f))
        try:
            df = pd.read_csv(f, sep='\t')
            
            if df.empty or 'Gene_Protein' not in df.columns or 'Type' not in df.columns:
                continue
            
            # --- FIXED SECTION ---
            # Use .str.lower() instead of .lower() for Pandas Series
            df = df[df['Gene_Protein'].notna()]
            df = df[df['Gene_Protein'].astype(str).str.lower() != 'none']
            # ---------------------
            
            if df.empty:
                continue

            df['Sample_ID'] = sample_id
            
            def categorize(row):
                # Ensure the value is treated as a string before lowering
                t = str(row['Type']).lower()
                if 'acquired' in t: return 'Acquired Gene'
                elif 'alignment' in t or 'mutation' in t: return 'Mutation'
                elif 'amplicon' in t: return 'Amplicon'
                return 'Other'

            df['Category'] = df.apply(categorize, axis=1)
            all_data.append(df)
        except Exception as e:
            # This captures the filename so you know which one failed
            print(f"Warning: Could not process {f}: {e}")

    master_df = pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()
    return master_df, all_sample_ids

def generate_csv_summary(df, all_sample_ids, output_path):
    if df.empty:
        summary_matrix = pd.DataFrame(index=all_sample_ids)
        summary_matrix.index.name = "Sample_ID"
    else:
        # Create a combined label for the CSV header
        df['Gene_Label'] = df['Category'] + ": " + df['Gene_Protein']
        summary_matrix = df.pivot_table(
            index='Sample_ID', 
            columns='Gene_Label', 
            aggfunc='size', 
            fill_value=0
        )
        summary_matrix = summary_matrix.reindex(all_sample_ids, fill_value=0)

    summary_matrix.to_csv(output_path)
    print(f"Matrix saved to: {output_path}")

def generate_heatmap(df, all_sample_ids, category, output_path):
    cat_df = df[df['Category'] == category] if not df.empty else pd.DataFrame()

    if cat_df.empty:
        matrix = pd.DataFrame(index=all_sample_ids)
    else:
        matrix = cat_df.pivot_table(index='Sample_ID', columns='Gene_Protein', aggfunc='size', fill_value=0)
        matrix = matrix.reindex(all_sample_ids, fill_value=0)

    plt.figure(figsize=(max(12, matrix.shape[1]*0.6 + 2), max(6, len(all_sample_ids)*0.3 + 2)))
    cmap = {"Acquired Gene": "YlGnBu", "Mutation": "OrRd", "Amplicon": "Purples"}.get(category, "Greys")

    if matrix.shape[1] == 0:
        plt.text(0.5, 0.5, f'No {category} Detected', ha='center', va='center')
        plt.yticks(range(len(all_sample_ids)), all_sample_ids, fontsize=8)
    else:
        sns.heatmap(matrix, annot=(matrix.shape[1] < 30), fmt="d", cmap=cmap, 
                    cbar_kws={'label': 'Evidence Count'}, linewidths=.5)

    plt.title(f"Heatmap: {category} (N={len(all_sample_ids)})", fontsize=16)
    plt.xlabel("Determinant")
    plt.ylabel("Sample ID")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Aggregate AMR results.")
    parser.add_argument("-i", "--input", required=True, help="Input directory")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    df_master, all_sample_ids = process_results(args.input)

    if not all_sample_ids:
        print("No samples found.")
        return

    # CSV Summary (Includes all genomes)
    generate_csv_summary(df_master, all_sample_ids, os.path.join(args.output, "amr_matrix_summary.csv"))

    # Heatmaps (Includes all genomes)
    for cat in ['Acquired Gene', 'Mutation', 'Amplicon']:
        safe_name = cat.lower().replace(" ", "_")
        generate_heatmap(df_master, all_sample_ids, cat, os.path.join(args.output, f"{safe_name}_heatmap.png"))

    print(f"\nDone. Results in {args.output}")

if __name__ == "__main__":
    main()
