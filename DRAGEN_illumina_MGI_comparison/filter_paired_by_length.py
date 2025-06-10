#!/usr/bin/env python3
import sys
import gzip

def filter_paired_fastq(r1_in, r2_in, r1_out, r2_out):
    """
    Read R1 and R2 in parallel; keep only those records where both
    seq_len == qual_len on R1 and seq_len == qual_len on R2.
    """
    with gzip.open(r1_in, 'rt') as f1, \
         gzip.open(r2_in, 'rt') as f2, \
         gzip.open(r1_out, 'wt') as o1, \
         gzip.open(r2_out, 'wt') as o2:

        while True:
            # Read 4 lines from R1
            h1 = f1.readline()
            if not h1:
                break
            s1    = f1.readline()
            plus1 = f1.readline()
            q1    = f1.readline()

            # Read 4 lines from R2
            h2 = f2.readline()
            if not h2:
                # If one file ends prematurely, stop
                break
            s2    = f2.readline()
            plus2 = f2.readline()
            q2    = f2.readline()

            # Strip trailing newline and compare lengths
            seq1 = s1.rstrip('\n')
            qual1 = q1.rstrip('\n')
            seq2 = s2.rstrip('\n')
            qual2 = q2.rstrip('\n')

            if len(seq1) == len(qual1) and len(seq2) == len(qual2):
                # Write both records if both pass the length check
                o1.write(h1)
                o1.write(s1)
                o1.write(plus1)
                o1.write(q1)

                o2.write(h2)
                o2.write(s2)
                o2.write(plus2)
                o2.write(q2)
            # else: drop both reads

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("Usage: {} R1_in.fastq.gz R2_in.fastq.gz R1_out.fastq.gz R2_out.fastq.gz\n"
                         .format(sys.argv[0]))
        sys.exit(1)

    r1_in, r2_in, r1_out, r2_out = sys.argv[1:5]
    filter_paired_fastq(r1_in, r2_in, r1_out, r2_out)

