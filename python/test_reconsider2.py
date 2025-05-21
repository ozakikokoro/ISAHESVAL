import sys
import pysam

def calculate_length_excluding_clips(cigar_tuples):
    length = 0
    for operation, length_op in cigar_tuples:
        if operation in [0, 2, 8]:  # Skip soft and hard clips, only count matches (M=0), deletions (D=2), and mismatches (X=8)
            length += length_op
    return length

# Check if the script has been given an argument
if len(sys.argv) < 4:
    print("Usage: python calculate_read_lengths.py <bam_file_path> <output_file_path> <spvar_allelesnv_pair_id> <evidence_class>")
    sys.exit(1)  # Exit the script if no argument is given

bam_file_path = sys.argv[1]  # Use the first argument provided as the BAM file path
output_file_path = sys.argv[2]
spvar_allelesnv_pair_id = sys.argv[3]
evidence_class = sys.argv[4]
with pysam.AlignmentFile(bam_file_path, "rb") as bam, open(output_file_path, "a") as outfile:
    for read in bam.fetch():
        read_length_excl_clips = calculate_length_excluding_clips(read.cigar)
        print(f"Read {read.query_name} length excluding clips: {read_length_excl_clips}")
        outfile.write(f"{spvar_allelesnv_pair_id}\t{read.query_name}\t{read_length_excl_clips}\t{evidence_class}\n")
