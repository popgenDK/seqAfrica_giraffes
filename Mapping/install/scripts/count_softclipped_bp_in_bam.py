#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pysam
from tqdm import tqdm

BAM_CMATCH = 0  # M
BAM_CINS = 1  # I
BAM_CDEL = 2  # D
BAM_CREF_SKIP = 3  # N
BAM_CSOFT_CLIP = 4  # S
BAM_CHARD_CLIP = 5  # H
BAM_CPAD = 6  # P
BAM_CEQUAL = 7  # =
BAM_CDIFF = 8  # X
BAM_CBACK = 9  # B


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def count_cigar(cigar, selection=(BAM_CSOFT_CLIP, BAM_CHARD_CLIP)):
    count = 0
    for operation, length in cigar:
        if operation in selection:
            count += length
        else:
            # Clipping is trailing
            break

    return count


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("root", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    # Silence log-messages from HTSLIB
    pysam.set_verbosity(0)

    columns = ("n", "n5p", "n3p", "bp5p", "bp3p")
    counts = defaultdict(lambda: dict.fromkeys(columns, 0))
    for it in args.root.iterdir():
        if it.suffix.lower() == ".bam":
            print("reading", it, file=sys.stderr)
            sample, genome, _ = it.name.split(".", 2)
            junk = "Y" if ".junk." in it.name else "N"

            with pysam.AlignmentFile(it, threads=3) as handle:
                record_count = 0
                clipped_5p_reads = 0
                clipped_5p_bases = 0
                clipped_3p_reads = 0
                clipped_3p_bases = 0

                for record in tqdm(handle, unit=" records", unit_scale=True):
                    if record.is_unmapped:
                        continue

                    record_count += 1

                    if clip_5p := count_cigar(record.cigar):
                        clipped_5p_reads += 1
                        clipped_5p_bases += clip_5p

                    if clip_3p := count_cigar(reversed(record.cigar)):
                        clipped_3p_reads += 1
                        clipped_3p_bases += clip_3p

                sample_counts = counts[(sample, genome, junk)]
                sample_counts["n"] += record_count
                sample_counts["n5p"] += clipped_5p_reads
                sample_counts["bp5p"] += clipped_5p_bases
                sample_counts["n3p"] += clipped_3p_reads
                sample_counts["bp3p"] += clipped_3p_bases

    print("sample", "genome", "junk", *columns, sep="\t")
    for (sample, genome, junk), sample_counts in sorted(counts.items()):
        values = (sample_counts[key] for key in columns)
        print(sample, genome, junk, *values, sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
