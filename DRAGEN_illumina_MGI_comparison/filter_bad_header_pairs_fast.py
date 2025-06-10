#!/usr/bin/env python3
import sys
import gzip
from itertools import islice

def valid_header(line):
    # line is a bytes object ending in \n
    # It must start with '@' and have at least one non-whitespace after it
    if not line.startswith(b'@'):
        return False
    # strip newline and any spaces; must have something left
    return len(line[1:].strip()) > 0

def filter_pairs(r1_in, r2_in, r1_out, r2_out):
    with gzip.open(r1_in, 'rb') as f1, \
         gzip.open(r2_in, 'rb') as f2, \
         gzip.open(r1_out, 'wb') as o1, \
         gzip.open(r2_out, 'wb') as o2:

        while True:
            # Grab the next 4 lines (one FASTQ record) from each
            rec1 = list(islice(f1, 4))
            rec2 = list(islice(f2, 4))

            # If either file is exhausted, stop
            if len(rec1) < 4 or len(rec2) < 4:
                break

            # Check both headers
            if valid_header(rec1[0]) and valid_header(rec2[0]):
                # Write them out in binary
                o1.writelines(rec1)
                o2.writelines(rec2)
            # else: drop both records

if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.stderr.write(f"Usage: {sys.argv[0]} R1.fastq.gz R2.fastq.gz OUT_R1.fastq.gz OUT_R2.fastq.gz\n")
        sys.exit(1)

    r1_in, r2_in, r1_out, r2_out = sys.argv[1:]
    filter_pairs(r1_in, r2_in, r1_out, r2_out)
