import csv
import glob
import os

for csv_path in glob.glob('*.csv'):
    tsv_path = os.path.splitext(csv_path)[0] + '.tsv'
    with open(csv_path, 'r') as csv_file, open(tsv_path, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerows(csv.reader(csv_file))
    print(f"Converted: {csv_path} → {tsv_path}")