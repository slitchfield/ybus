#!/usr/bin/env python3
"""Summarizes samply profiles of profile_pf into stages and top functions, for
scripts/profile-pf.sh.

Usage: profile_summary.py <profile.json> [<profile.json> ...]
  Each profile needs its symbol sidecar (<name>.syms.json, from samply's
  --unstable-presymbolicate). Prints a Markdown section per profile: the samples, the average
  number of busy threads (samples over the run's wall time at the sampling rate, when the run's
  .log from profile_pf gives the wall time), each stage's share, and the functions with the most
  samples of their own and in total.

A sample's stage is decided by the innermost frame on its stack that matches one of STAGES, in
their order within a frame; samples matching none are "other".
"""

import bisect
import collections
import json
import re
import sys
from pathlib import Path

# (stage, patterns): a frame matches a stage if its demangled name contains any of the patterns
STAGES = [
    ("allocation and copying", ["malloc", "_libc_free", "cfree", "realloc", "brk",
                                "finish_grow", "grow_one", "memmove", "memcpy", "memset",
                                "_int_malloc", "_int_free"]),
    ("parsing", ["case::matpower::", "Case>::from_path", "Case>::from_bytes", "dec2flt",
                 "FromStr"]),
    ("ordering and symbolic analysis", ["approximate_minimum_degree", "colamd", "ordering::",
                                        "Adjacency", "analyze_lu", "analyze_cholesky",
                                        "symbolic", "elimination_tree", "column_counts"]),
    ("refactorization", ["refactor"]),
    ("numeric factorization", ["gilbert_peierls", "Lu>::factor", "Cholesky>::factor",
                               "cholesky::", "factor"]),
    ("triangular solves", ["Lu>::solve", "Cholesky>::solve", "solve_sparse", "triangular",
                           "lower_solve", "upper_solve"]),
    ("outage updates", ["outage", "woodbury", "Woodbury"]),
    ("assembly and refill", ["JacobianAssembly", "SupersetJacobian", "source_value", "refill",
                             "assemble", "from_triplets", "stamps", "Bdc>::new",
                             "Admittances>::new"]),
    ("model build", ["AcModel>::new", "DcModel>::new", "Topology", "Network>::new",
                     "bustypes"]),
    ("mismatch and flows", ["mismatch", "dot_kernel", "mul_vec", "vec_mul", "pfsoln",
                            "share_reactive", "hypot", "sbus"]),
    ("parallel runtime", ["rayon", "crossbeam"]),
]


def symbolizer(syms):
    strings = syms["string_table"]
    tables = {}
    for lib in syms["data"]:
        rows = sorted(lib["symbol_table"], key=lambda row: row["rva"])
        tables[lib["debug_name"]] = ([row["rva"] for row in rows], rows)

    def name(lib, address):
        if lib in tables:
            starts, rows = tables[lib]
            i = bisect.bisect_right(starts, address) - 1
            if i >= 0 and address < rows[i]["rva"] + rows[i]["size"]:
                return strings[rows[i]["symbol"]]
        return f"[{lib}]"

    return name


def stage_of(names):
    for frame in names:
        for stage, patterns in STAGES:
            if any(pattern in frame for pattern in patterns):
                return stage
    return "other"


def summarize(path):
    profile = json.loads(Path(path).read_text())
    syms_path = Path(str(path).removesuffix(".json") + ".syms.json")
    name = symbolizer(json.loads(syms_path.read_text()))
    interval_ms = profile["meta"]["interval"]
    own, total, stages = collections.Counter(), collections.Counter(), collections.Counter()
    samples = 0
    for thread in profile["threads"]:
        table = thread["samples"]
        if not table["length"]:
            continue
        frames, funcs = thread["frameTable"], thread["funcTable"]
        stacks, resources = thread["stackTable"], thread["resourceTable"]

        def frame_name(frame):
            resource = funcs["resource"][frames["func"][frame]]
            lib = "?"
            if resource is not None and resource >= 0:
                lib = profile["libs"][resources["lib"][resource]]["debugName"]
            return name(lib, frames["address"][frame])

        cache = {}
        for i in range(table["length"]):
            stack = table["stack"][i]
            if stack is None:
                continue
            if stack not in cache:
                names, s = [], stack
                while s is not None:
                    names.append(frame_name(stacks["frame"][s]))
                    s = stacks["prefix"][s]
                cache[stack] = (names, stage_of(names))
            names, stage = cache[stack]
            weight = table["weight"][i] if table["weight"] else 1
            samples += weight
            own[names[0]] += weight
            stages[stage] += weight
            for frame in set(names):
                total[frame] += weight

    title = Path(path).name.removesuffix(".json")
    lines = [f"### {title}", ""]
    log = Path(str(path).removesuffix(".json") + ".log")
    busy = ""
    if log.exists():
        match = re.search(r"(\d+) runs, ([\d.]+) ms per run", log.read_text())
        if match:
            wall_ms = int(match.group(1)) * float(match.group(2))
            busy = f", {samples * interval_ms / wall_ms:.1f} threads busy on average"
            lines[0] += f" ({match.group(1)} runs, {float(match.group(2)):.3f} ms per run)"
    lines += [f"{samples} samples at {1e3 / interval_ms:.0f} Hz{busy}.", ""]
    lines += ["| Stage | Share |", "|---|--:|"]
    for stage, count in stages.most_common():
        lines.append(f"| {stage} | {count / samples:.1%} |")
    lines += ["", "| Most samples of their own | Share |", "|---|--:|"]
    for frame, count in own.most_common(12):
        lines.append(f"| `{shorten(frame)}` | {count / samples:.1%} |")
    ours = [(f, c) for f, c in total.most_common()
            if ("powerflow::" in f or "sparse_mat::" in f) and c < samples]
    lines += ["", "| Most samples in total (powerflow and sparse_mat) | Share |", "|---|--:|"]
    for frame, count in ours[:12]:
        lines.append(f"| `{shorten(frame)}` | {count / samples:.1%} |")
    return "\n".join(lines) + "\n"


def strip_generics(text):
    depth, out = 0, []
    for ch in text:
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth -= 1
        elif depth == 0:
            out.append(ch)
    return "".join(out)


def shorten(frame):
    """A readable name: `<path::Type<T> as Trait>::method<U>` becomes `Type::method`, and other
    names lose their generic arguments and crate prefixes."""
    if frame.startswith("<"):
        depth = 0
        for j, ch in enumerate(frame):
            depth += {"<": 1, ">": -1}.get(ch, 0)
            if depth == 0:
                break
        inner, rest = frame[1:j], frame[j + 1:]
        inner = strip_generics(inner).split(" as ")[0]
        frame = inner.split("::")[-1] + strip_generics(rest)
    else:
        frame = strip_generics(frame)
    for prefix in ("sparse_mat::compressed_matrix::", "powerflow::"):
        frame = frame.replace(prefix, "")
    return frame.rstrip(":").replace("|", "\\|")[:100]


if __name__ == "__main__":
    for path in sys.argv[1:]:
        print(summarize(path))
