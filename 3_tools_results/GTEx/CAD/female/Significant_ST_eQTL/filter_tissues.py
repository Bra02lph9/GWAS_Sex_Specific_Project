import os
import glob
import pandas as pd
import argparse

def filter_csv_files(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    csv_files = glob.glob(os.path.join(input_dir, "*.csv"))

    if not csv_files:
        print(f"No CSV files found in '{input_dir}'")
        return

    total_in, total_out, skipped = 0, 0, 0

    for filepath in csv_files:
        filename = os.path.basename(filepath)
        try:
            df = pd.read_csv(filepath)

            if "Tissue" not in df.columns:
                print(f"  [SKIP] {filename} — no 'Tissue' column found")
                skipped += 1
                continue

            mask = df["Tissue"].str.contains("heart|artery", case=False, na=False)
            filtered = df[mask]

            total_in  += len(df)
            total_out += len(filtered)

            out_path = os.path.join(output_dir, filename)
            filtered.to_csv(out_path, index=False)
            print(f"  [OK] {filename}: {len(df)} → {len(filtered)} rows kept")

        except Exception as e:
            print(f"  [ERROR] {filename}: {e}")
            skipped += 1

    print(f"\nDone. {len(csv_files) - skipped} file(s) processed.")
    print(f"Total rows: {total_in} → {total_out} kept "
          f"({total_in - total_out} removed).")
    if skipped:
        print(f"{skipped} file(s) skipped.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter CSVs to rows whose Tissue contains 'heart' or 'artery'."
    )
    parser.add_argument("input_dir",  help="Directory containing the source CSV files")
    parser.add_argument("output_dir", help="Directory where filtered CSVs will be written")
    args = parser.parse_args()

    filter_csv_files(args.input_dir, args.output_dir)
