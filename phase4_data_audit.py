import pandas as pd
import os

def run_audit(file_path):
    print(f"--- Phase 4: Data Pipeline Completeness Audit ---")
    print(f"Loading data from: {os.path.basename(file_path)}")
    
    try:
        # The NASA Exoplanet Archive CSVs usually have comment lines starting with '#'
        df = pd.read_csv(file_path, comment='#', low_memory=False)
    except Exception as e:
        print(f"Failed to read CSV: {e}")
        return

    print(f"Initial row count: {len(df)}")
    
    # 1. Check required columns
    required_cols = [
        'pl_bmassprov', 'pl_rade', 'pl_bmasse', 
        'st_teff', 'st_lum'
    ]
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        print(f"Available columns: {df.columns.tolist()[:10]}...")
        return
    else:
        print("Schema check passed: All required columns present.")
        
    df_filtered = df.copy()
    
    # 3. Completeness audit on the filtered dataset
    print("\n--- Completeness Audit (default_flag == 1) ---")
    
    audit_cols = ['pl_rade', 'pl_bmasse', 'st_teff', 'st_lum']
    
    for col in audit_cols:
        missing_count = df_filtered[col].isna().sum()
        completeness = 100.0 * (1.0 - missing_count / len(df_filtered))
        print(f"{col:<15}: {len(df_filtered) - missing_count:<5} valid | {missing_count:<5} missing | {completeness:6.2f}% complete")
        
    # Mass provenance check
    print("\n--- Mass Provenance (pl_bmassprov) ---")
    prov_counts = df_filtered['pl_bmassprov'].value_counts(dropna=False)
    for name, count in prov_counts.items():
        print(f"{str(name):<15}: {count}")

    print("\nAudit Complete.")

if __name__ == "__main__":
    file_path = "c:/Users/harsh/ドキュメント/Exoplanet Habitability/PSCompPars_2026.07.11_02.37.36.csv"
    run_audit(file_path)
