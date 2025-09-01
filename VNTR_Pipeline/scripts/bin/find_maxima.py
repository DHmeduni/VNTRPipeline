import argparse
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam



def extract_lengths(input_file):
    """
    Extract sequence lengths from a FASTQ or BAM file.
    """
    lengths = []
    if input_file.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        for record in SeqIO.parse(input_file, "fastq"):
            lengths.append(len(record.seq))
    elif input_file.endswith('.bam'):
        with pysam.AlignmentFile(input_file, "rb") as bam:
            for read in bam:
                lengths.append(read.query_length)
    else:
        raise ValueError("Unsupported file format. Please use FASTQ or BAM files.")
    return lengths

def generate_length_histogram(input_file, bin_size=10, output_file=None, x_min=None, x_max=None):
    """
    Generate a histogram of sequence lengths from a file.
    """
    lengths = extract_lengths(input_file)
    
    # Determine x-axis limits
    min_length = min(lengths) if x_min is None else max(int(x_min), min(lengths))
    max_length = max(lengths) if x_max is None else min(int(x_max), max(lengths))
    
    # Create histogram plot
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=range(min_length, max_length + bin_size, bin_size), edgecolor='black')
    plt.title(f"Sequence Length Histogram for {input_file}")
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Frequency")
    plt.xlim(min_length, max_length)
    
    return lengths, plt

import numpy as np

import numpy as np

import numpy as np

import numpy as np

def find_two_maxima_above_threshold(lengths, threshold=2500, min_separation=50,
                                    min_frequency=20, bin_size=10, percentage_min_frequency=0.01):
    """
    Find two maxima in a histogram of sequence lengths above a threshold.
    Rules:
      - Peaks are separated by at least `min_separation`
      - First maxima is the peak with highest frequency
      - If a candidate peak is shorter than first maxima, its frequency must be >= 100% of first maxima frequency
      - If a candidate peak is longer than first maxima, its frequency must be >= adaptive_min_freq
    """
    # Filter lengths above threshold
    filtered_lengths = [length for length in lengths if length > threshold]
    if not filtered_lengths:
        return []

    # Create histogram
    bins = np.arange(min(filtered_lengths), max(filtered_lengths) + bin_size, bin_size)
    histogram, bin_edges = np.histogram(filtered_lengths, bins=bins)

    # Adaptive minimum frequency
    max_peak_freq = np.max(histogram)
    adaptive_min_freq = max(min_frequency, percentage_min_frequency * max_peak_freq)

    # Find all peaks above adaptive minimum frequency
    peak_indices = np.where((histogram[1:-1] > histogram[:-2]) &
                            (histogram[1:-1] > histogram[2:]) &
                            (histogram[1:-1] >= adaptive_min_freq))[0] + 1

    if len(peak_indices) == 0:
        return []

    # First maxima: highest frequency
    first_idx = peak_indices[np.argmax(histogram[peak_indices])]
    first_freq = histogram[first_idx]
    first_peak = bin_edges[first_idx]

    maxima = [(first_freq, first_peak)]

    # Collect valid candidates for second maxima
    candidates = []
    for idx in peak_indices:
        if idx == first_idx:
            continue
        freq = histogram[idx]
        peak = bin_edges[idx]

        # Check separation
        if abs(peak - first_peak) < min_separation:
            continue

        # Apply rules based on relative peak position
        if peak < first_peak:  # shorter candidate
            if freq < 1 * first_freq:
                continue
        else:  # longer candidate
            if freq < adaptive_min_freq:
                continue

        candidates.append((freq, peak))

    # Pick best candidate by frequency if available
    if candidates:
        second_freq, second_peak = max(candidates, key=lambda x: x)
        maxima.append((second_freq, second_peak))

    # Return the x-values of maxima (sorted)
    return sorted([p for p in maxima])


def find_maxima_from_file(input_file, threshold=2500, separation=50, output_file=None, min_frequency=20, bin_size=10, x_min=None, x_max=None, percentage_min_frequency=0.01):
    lengths, _ = generate_length_histogram(input_file, bin_size, None, x_min, x_max)
    maxima = find_two_maxima_above_threshold(lengths, threshold, separation, min_frequency, bin_size, percentage_min_frequency)
    peaks = [peak for freq, peak in maxima]

    if output_file is not None:
        with open(output_file, 'w') as f:
            for peak in peaks:
                f.write(f"{peak}\n")

    return peaks


def main():
    parser = argparse.ArgumentParser(description="Generate a sequence length histogram and find two maxima above a threshold.")
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ or BAM file')
    parser.add_argument('-b', '--bin_size', type=int, default=10, help='Bin size for the histogram (default is 10)')
    parser.add_argument('-o', '--output', help='Output file to save the histogram (e.g., histogram.png)')
    parser.add_argument('--x_min', type=int, help="Minimum value for the x-axis (lower limit)")
    parser.add_argument('--x_max', type=int, help="Maximum value for the x-axis (upper limit)")
    parser.add_argument('-t', '--threshold', type=int, default=2500, help="Threshold for finding maxima (default is 2500)")
    parser.add_argument('-s', '--separation', type=int, default=50, help="Minimum separation between maxima (default is 100)")
    parser.add_argument('-m', '--min_frequency', type=int, default=20, help="Minimum frequency to be determined to be a maxima")
    
    args = parser.parse_args()
    
    # Generate histogram and extract sequence lengths
    lengths, plt = generate_length_histogram(args.input, args.bin_size, args.output, args.x_min, args.x_max)
    
    # Find two maxima above the threshold with a minimum separation
    maxima = find_two_maxima_above_threshold(lengths, args.threshold, args.separation, args.min_frequency, args.bin_size)
    
    if len(maxima) == 2:
        print(f"The two maxima above the threshold of {args.threshold} bp, separated by at least {args.separation} bp are: {maxima[0][1]:.0f} bp and {maxima[1][1]:.0f} bp")
        for i in range(len(maxima)):
            plt.axvline(maxima[i][1], linestyle='--', label=f"Maxima: {maxima[i][1]:.0f} bp")
        plt.legend()
        
        if args.output:
            plt.savefig(args.output)
            print(f"Histogram saved to {args.output}")
        else:
            plt.show()
            
    elif len(maxima) == 1:
        print(f"Only one maxima above the threshold of {args.threshold} bp, separated by at least {args.separation} bp is: {maxima[0][1]:.0f} bp")
        for i in range(len(maxima)):
            plt.axvline(maxima[i][1], linestyle='--', label=f"Maxima: {maxima[i][1]:.0f} bp")
        plt.legend()
        
        if args.output:
            plt.savefig(args.output)
            print(f"Histogram saved to {args.output}")
        else:
            plt.show()
        
        
    else:
        if args.output:
            print(f"Could not find one or two maxima above {args.threshold} bp separated by at least {args.separation} bp")
            plt.savefig(args.output)
            print(f"Histogram saved to {args.output}")
        else:
            print(f"Could not find two maxima above {args.threshold} bp separated by at least {args.separation} bp")
            plt.show()
        

if __name__ == "__main__":
    main()

