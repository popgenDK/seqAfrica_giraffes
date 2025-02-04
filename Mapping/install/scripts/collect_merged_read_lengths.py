#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import sys
from pathlib import Path


def find_settings(paths):
    paths = list(paths)
    while paths:
        path = paths.pop()

        if path.is_dir():
            paths.extend(path.iterdir())
        elif path.name.endswith(".settings"):
            yield path


def parse_settings(filepath):
    sections = {}

    with filepath.open("rt") as handle:
        section = None
        for line in handle:
            if line := line.strip():
                if line.startswith("["):
                    section = []
                    sections[line] = section
                elif section is not None:
                    section.append(line)

    statistics = {}
    for line in sections["[Trimming statistics]"]:
        key, value = line.split(":", 1)
        statistics[key] = value

    npairs = int(statistics["Total number of read pairs"])
    nmerged = int(statistics["Number of full-length collapsed pairs"]) + int(
        statistics["Number of truncated collapsed pairs"]
    )

    total_reads = 0
    total_bases = 0

    lengths = sections["[Length distribution]"]
    header = lengths[0].split("\t")
    for line in lengths[1:]:
        row = dict(zip(header, line.split("\t")))
        length = int(row["Length"])
        for key in ("Collapsed", "CollapsedTruncated"):
            n = int(row[key])
            total_reads += n
            total_bases += n * length

    return {
        "n_pairs": npairs,
        "n_merged": nmerged,
        "f_merged": nmerged / npairs,
        "average_len": total_bases / total_reads,
    }


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        description="Collects AdapterRemoval settings files (*.settings) and prints "
        "the fractin of merged reads and the average lengths of these",
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        metavar="path",
        help=".settings files or folders to search (recursively) for .settings files",
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    print("Path", "NPairs", "NMerged", "FracMerged", "MergedLength", sep="\t")
    for path in sorted(find_settings(args.paths)):
        row = parse_settings(path)
        print(
            path,
            row["n_pairs"],
            row["n_merged"],
            row["f_merged"],
            row["average_len"],
            sep="\t",
        )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
