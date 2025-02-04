#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import collections
import sys

from pathlib import Path

import pysam

from paleomix.common.timer import BAMTimer


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("files", nargs="+", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    print("File\tLength\tPctMapped\tCount\tFrac")
    for filepath in args.files:
        counts = collections.Counter()
        with pysam.AlignmentFile(filepath) as handle:
            for read in BAMTimer(handle):
                if not read.is_unmapped and read.query_name.startswith("M_"):
                    length = read.query_length
                    mapped_bases = sum(
                        count for op, count in read.cigartuples if op in (0, 7, 8)
                    )
                    pct = (100 * mapped_bases) // length

                    counts[(length, pct)] += 1

        totals = collections.Counter()
        for (length, _), count in counts.items():
            totals[length] += count

        for (length, mapped_bases), count in sorted(counts.items()):
            print(
                "{}\t{}\t{}\t{}\t{:.3f}".format(
                    filepath.name,
                    length,
                    mapped_bases,
                    count,
                    count / totals[length],
                )
            )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
