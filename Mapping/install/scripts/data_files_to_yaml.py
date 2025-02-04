#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import re
import sys

from pathlib import Path


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("--root", default=Path("data"), type=Path)
    parser.add_argument("--regex-from", default="_1.fq.gz")
    parser.add_argument("--regex-to", default="_{Pair}.fq.gz")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    for sample_path in sorted(args.root.iterdir()):
        if not sample_path.is_dir():
            continue

        print("{}:".format(sample_path.name))
        print("  {}:".format(sample_path.name))
        print("    Library1:")
        for filepath in sorted(sample_path.iterdir()):
            if re.search(args.regex_from, filepath.name):
                name = re.sub(args.regex_from, "", filepath.name)
                tmpl = filepath.parent / re.sub(
                    args.regex_from, args.regex_to, filepath.name
                )

                print("      {}: {}".format(name, tmpl))

        print()
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
