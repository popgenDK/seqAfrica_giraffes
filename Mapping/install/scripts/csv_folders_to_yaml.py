#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import json
import os
import sys
from itertools import groupby
from pathlib import Path


def read_csv(filepath, sep, required_headers):
    rows = []
    with filepath.open("rt") as handle:
        header = None
        for idx, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split(sep)
            if header is None:
                header = fields
                continue

            if len(fields) != len(header):
                raise RuntimeError(f"line {filepath}:{idx} has wrong number of columns")

            row = dict(zip(header, fields))
            missing_headers = set(required_headers) - row.keys()
            if missing_headers:
                print("ERROR: Missing columns:", missing_headers, file=sys.stderr)
                sys.exit(1)

            rows.append(row)

    return rows


def find_fq_fpaths(folder):
    files = []
    for it in folder.iterdir():
        if it.is_file():
            if it.name.lower().endswith((".fastq.gz", ".fq.gz")):
                files.append(it)
        elif it.is_dir():
            files.extend(find_fq_fpaths(it))

    return files


def is_fastq_pair(fpath1, fpath2):
    fpath1 = str(fpath1)
    fpath2 = str(fpath2)

    if len(fpath1) != len(fpath2):
        return False
    elif fpath1 > fpath2:
        fpath1, fpath2 = fpath2, fpath1

    found_mate = False
    for a, b in zip(fpath1, fpath2):
        if a != b:
            if found_mate:
                return False
            elif (a, b) == ("1", "2"):
                found_mate = True

    return found_mate


def merge_fpath_pair(fpath1, fpath2):
    fpath1 = str(fpath1)
    fpath2 = str(fpath2)

    assert len(fpath1) == len(fpath2)
    return "".join((a if a == b else "{Pair}") for a, b in zip(fpath1, fpath2))


def merge_fpaths(files):
    result = []
    files = sorted(files)
    while len(files) > 1:
        fastq1 = files[-2]
        fastq2 = files[-1]

        if is_fastq_pair(fastq1, fastq2):
            result.append(merge_fpath_pair(fastq1, fastq2))

            files.pop()
            files.pop()
        else:
            result.append(files.pop())

    result.extend(files)

    return result


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--sep", default=",")
    parser.add_argument("--col-name", default="Sample_ID")
    parser.add_argument("--col-folder", default="data")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    files = read_csv(
        filepath=args.csv,
        sep=args.sep,
        required_headers=(args.col_name, args.col_folder),
    )

    files.sort(key=lambda it: (it[args.col_name], it[args.col_folder]))

    for name, rows in groupby(files, key=lambda it: it[args.col_name]):
        print("  {}:".format(json.dumps(name)))
        print('    "Library1":')

        for row in rows:
            for fpath in merge_fpaths(find_fq_fpaths(Path(row[args.col_folder]))):
                basename = []
                for field in os.path.basename(fpath).split("_"):
                    basename.append(field)
                    if len(field) > 1 and field.startswith("L") and field[1].isdigit():
                        break
                basename = "_".join(basename)

                print(
                    "      {}: {}".format(
                        json.dumps(basename),
                        json.dumps(str(fpath)),
                    ),
                )

        print()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
