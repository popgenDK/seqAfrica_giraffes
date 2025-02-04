#!/usr/bin/env python3
# -*- coding: utf8 -*-
# pyright: strict
import argparse
import copy
import sys
import warnings
from pathlib import Path
from typing import Any, List, IO, Optional
from dataclasses import dataclass

import ruamel.yaml
import ruamel.yaml.error

OPTIONS_KEY = "Options"
GENOMES_KEY = "Genomes"


@dataclass
class Fastq:
    batch: str
    group: str
    sample: str
    library: str
    fastq_1: str
    fastq_2: Optional[str] = None


def collect_samples(batch: str, data: object):
    assert isinstance(data, dict), data
    for group, samples in sorted(data.items()):
        assert isinstance(group, str)
        assert isinstance(samples, dict)

        if group in (OPTIONS_KEY, GENOMES_KEY):
            continue

        for sample, libraries in sorted(samples.items()):
            assert isinstance(sample, str)
            assert isinstance(libraries, dict)

            if sample == OPTIONS_KEY:
                continue

            for library, lanes in sorted(libraries.items()):
                assert isinstance(library, str)
                assert isinstance(lanes, dict)

                if library == OPTIONS_KEY:
                    continue

                for lane, path in sorted(lanes.items()):
                    assert isinstance(lane, str)
                    assert isinstance(path, (str, dict))

                    def create(path: str) -> Fastq:
                        if "{Pair}" in path:
                            fastq_1 = path.format(Pair=1)
                            fastq_2 = path.format(Pair=2)
                        else:
                            fastq_1 = path
                            fastq_2 = None

                        return Fastq(
                            batch=batch,
                            group=group,
                            sample=sample,
                            library=library,
                            fastq_1=fastq_1,
                            fastq_2=fastq_2,
                        )

                    if lane == OPTIONS_KEY:
                        continue
                    elif isinstance(path, dict):
                        for path in path.values():
                            yield create(path)
                    else:
                        yield create(path)


def safe_load(stream: IO[str]) -> Any:
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    # TODO: Switch to 1.2; requires handling of "yes"/"no" values for bools
    yaml.version = (1, 1)  # type: ignore

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ruamel.yaml.error.MantissaNoDotYAML1_1Warning)

        return yaml.load(stream)


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.add_argument("yaml", nargs="+", type=Path)

    return parser.parse_args(argv)


def main(argv: List[str]) -> int:
    args = parse_args(argv)

    print("batch", "group", "sample", "library", "fastq1", "fastq2", sep="\t")
    for filepath in args.yaml:
        with filepath.open() as handle:
            data = safe_load(handle)

        batch, *_ = filepath.name.split(".", 1)
        for it in collect_samples(batch, data):
            print(
                it.batch,
                it.group,
                it.sample,
                it.library,
                Path(it.fastq_1).resolve(),
                "." if it.fastq_2 is None else Path(it.fastq_2).resolve(),
                sep="\t",
            )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
