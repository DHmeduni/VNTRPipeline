#!/usr/bin/env python3
import os
import sys
import pandas as pd

def merge_xlsx_to_xlsx(search_dir, filename_pattern, output_xlsx):
    all_files = []
    for root, dirs, files in os.walk(search_dir):
        for file in files:
            if file == filename_pattern:
                all_files.append(os.path.join(root, file))

    if not all_files:
        print(f"No files named '{filename_pattern}' found in {search_dir}")
        return

    dfs = []
    for file in all_files:
        try:
            df = pd.read_excel(file)
            dfs.append(df)
            print(f"Read {file}")
        except Exception as e:
            print(f"Warning: Could not read {file}: {e}")

    if dfs:
        result = pd.concat(dfs, ignore_index=True)
        result.to_excel(output_xlsx, index=False)
        print(f"Saved merged Excel file to {output_xlsx}")
    else:
        print("No dataframes to concatenate.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <search_directory> <filename_to_search> <output_xlsx_path>")
        sys.exit(1)

    search_directory = sys.argv[1]
    filename = sys.argv[2]
    output_xlsx_path = sys.argv[3]

    merge_xlsx_to_xlsx(search_directory, filename, output_xlsx_path)

