import pandas as pd
import os
from config_loader import SEQUENCES, GA_PARAMS, MOTIF_LENGTHS
from utils.folding import Fold
from core.ga_engine import GAEngine

A_MOTIF_LEN = MOTIF_LENGTHS["A_MOTIF_LEN"]
B_MOTIF_LEN = MOTIF_LENGTHS["B_MOTIF_LEN"]

def add_motif_columns_if_missing(df):
    """Ensure all motif coordinate columns exist in the DataFrame"""
    motif_cols = ["a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m"]
    if all(col in df.columns for col in motif_cols):
        return df
    
    # Calculate motifs based on existing coordinates if missing
    df["a1m"] = df["a2"] - A_MOTIF_LEN
    df["a2m"] = df["a2"]
    df["A1m"] = df["A2"] - A_MOTIF_LEN
    df["A2m"] = df["A2"]
    df["b1m"] = df["b2"] - B_MOTIF_LEN
    df["b2m"] = df["b2"]
    df["B1m"] = df["B2"] - B_MOTIF_LEN
    df["B2m"] = df["B2"]
    
    print("Added missing motif columns to DataFrame")
    return df

def main():
    # 1) Load initial data
    input_csv = "initial_candidates_ISEc11_ISPa11.csv"
    if not os.path.exists(input_csv):
        # Fallback for demonstration if file is in parent dir
        input_csv = os.path.join("..", input_csv)
        
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return

    C1 = pd.read_csv(input_csv)
    C1 = add_motif_columns_if_missing(C1)
    
    # Ensure integer types for coordinates
    int_cols = ["i", "a1", "a2", "A1", "A2", "b1", "b2", "B1", "B2",
                "a1m", "a2m", "A1m", "A2m", "b1m", "b2m", "B1m", "B2m"]
    for col in int_cols:
        if col in C1.columns:
            C1[col] = C1[col].astype(int)

    # 2) Pre-fold reference sequences
    print("Folding reference sequences A and B...")
    struct_A, _ = Fold(SEQUENCES["A"], mode="full")
    struct_B, _ = Fold(SEQUENCES["B"], mode="full")

    # 3) Initialize and run GA
    engine = GAEngine(C1, struct_A, struct_B)
    final_pop = engine.run()

    # 4) Save results
    output_path = "ga_results.csv"
    final_pop.to_csv(output_path, index=False)
    print(f"GA complete. Results saved to {output_path}")

    # Motif Verification for top 5
    print("\n" + "="*50)
    print("Motif Verification for Top 5 Candidates:")
    print("="*50)

    from utils.sequence import rev_comp
    for idx, row in final_pop.head(5).iterrows():
        seq = row["seq"]
        print(f"\nCandidate {idx} (Score: {row['Score']:.4f})")
        print(f"  A_motif_top : {seq[int(row['a1m']):int(row['a2m'])]}")
        print(f"  A_motif_bot : {seq[int(row['A1m']):int(row['A2m'])]} (rc match expected)")
        print(f"  B_motif_top : {seq[int(row['b1m']):int(row['b2m'])]}")
        print(f"  B_motif_bot : {seq[int(row['B1m']):int(row['B2m'])]} (rc match expected)")

if __name__ == "__main__":
    main()

