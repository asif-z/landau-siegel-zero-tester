import os
import glob
import csv
import argparse
from collections import Counter

def process_files(indir, outfile, d, N):
    basket_counts = Counter()

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
                val = int(row[2])  # only 3rd column matters
                if val == -1:
                    basket_counts["-1"] += 1
                else:
                    lower_bound = get_bucket_lower_bound(val)
                    basket_counts[lower_bound] += 1

    # Write results
    with open(outfile, "w") as f:
        f.write("bucket_lower_bound,count\n")
        for key in sorted(k for k in basket_counts.keys() if k != "-1"):
            f.write(f"{key},{basket_counts[key]}\n")
        f.write(f"-1,{basket_counts.get('-1',0)}\n")

    print(f"Results written to {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bucket 3rd column of CSVs into buckets of size d")
    parser.add_argument("--indir", required=True, help="Directory containing output_rank_*.csv")
    parser.add_argument("--outfile", required=True, help="Output summary file")
    parser.add_argument("--d", type=int, required=True, help="Bucket size")
    parser.add_argument("--N", type=int, required=True, help="Maximum integer value")
    args = parser.parse_args()

    process_files(args.indir, args.outfile, args.d, args.N)
