#find the failed modulus from output_rank_*.csv

import os
import csv
import re
import glob
import csv

def search_fail_in_rank_files(directory='.'):

    out = open("output.txt","w")

    keyword = "fail"
    keyword_lower = keyword.lower()

    print(f"Searching for '{keyword}' in output_rank_*.csv files under: {directory}")
    print("-" * 60)

    files = glob.glob(os.path.join(directory, "output_rank_*.csv"))

    for fpath in files:
        print(f"Processing {fpath} ...")
        with open(fpath,"r") as f:
            reader = csv.reader(f)
            for row in reader:
                if(row[1].lower()==keyword_lower):
                    out.write(row[0]+"\n")


# Run the search
if __name__ == '__main__':
    search_dir = input("Enter directory to search (leave blank for current): ").strip()
    if not search_dir:
        search_dir = '.'
    search_fail_in_rank_files(search_dir)
    print("Done!")
