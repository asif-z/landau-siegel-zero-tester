#draw a histgram from a data source

import csv
import matplotlib.pyplot as plt

def draw_histogram(filename, d):
    buckets = []
    counts = []

    # Read CSV file, skip header
    with open(filename, "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header row
        for row in reader:
            if len(row) < 2:
                continue
            bucket_start = int(row[0])
            count = int(row[1])
            buckets.append(bucket_start)
            counts.append(count)

    # Draw histogram as bar chart
    plt.bar(buckets, counts, width=d, align="edge")

    plt.xlabel("Number of Primes Needed to Violate the Inequality")
    plt.ylabel("Count")
    plt.title("Distribution of Number of Primes Needed")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    plt.show()

if __name__ == "__main__":
    draw_histogram("run2.txt", d=5000)
