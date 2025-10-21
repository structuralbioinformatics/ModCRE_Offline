import pandas as pd 
import matplotlib.pyplot as plt
import sys
file_path = sys.argv[1]
output_path = sys.argv[2]
column = sys.argv[3] # normal_s3dc_dd
# Read the CSV (handle commas and possible index column)
df = pd.read_csv(file_path) 
# In case the first column is unnamed (like an index), drop it
if df.columns[0] == '':
    df = df.drop(columns=df.columns[0])
 
 
# Ensure the relevant columns exist
if "Position" not in df.columns or column not in df.columns:
    print("❌ Required columns not found. Available columns:")
    print(df.columns.tolist())
else:
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(df["Position"], df[column], color="darkblue", linewidth=2)
    plt.xlabel("Position", fontsize=12)
    plt.ylabel(column, fontsize=12)
    plt.title(f"Profile: {column} across Positions", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path,dpi=300)
    print(f'plot saved to {output_path}')
    
