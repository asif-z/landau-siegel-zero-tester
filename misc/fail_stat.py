#generate a number of failed statistics from output_rank_*.csv

import os
import glob
import csv
import argparse
from collections import Counter

def process_files(indir, outfile, d):
    counts = Counter()

    def get_bucket_lower_bound(val):
        """Return lower bound of the bucket for a given value."""
        return (val // d) * d

    files = glob.glob(os.path.join(indir, "output_rank_*.csv"))
    print(f"Found {len(files)} files")

    for fpath in files:
        print(f"Processing {fpath} ...")
        with open(fpath, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) < 3:
                    continue  # skip malformed lines
                if(row[2]=="-1"):
                    lower_bound = get_bucket_lower_bound(abs(int(row[0])))
                    counts[lower_bound] += 1

    # Write results
    with open(outfile, "w") as f:
        f.write("bucket_lower_bound,count\n")
        for key in sorted(counts.keys()):
            f.write(f"{key},{counts[key]}\n")

    print(f"Results written to {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bucket 3rd column of CSVs into buckets of size d")
    parser.add_argument("--indir", required=True, help="Directory containing output_rank_*.csv")
    parser.add_argument("--outfile", required=True, help="Output summary file")
    parser.add_argument("--d", type=int, required=True, help="Bucket size")
    args = parser.parse_args()

    process_files(args.indir, args.outfile, args.d)
