#!/usr/bin/env python3
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
from math import comb
from pathlib import Path
import json
import subprocess
import sys


def parse_factor_output(text: str) -> tuple[int, dict[int, int]]:
    line = text.strip().splitlines()[-1]
    left, right = line.split(":", 1)
    factors = Counter(int(token) for token in right.split())
    return int(left), dict(factors)


def is_sum_of_two_squares(factors: dict[int, int]) -> bool:
    return all(not (p % 4 == 3 and e % 2 == 1) for p, e in factors.items())


def parse_candidates(path: Path):
    rows = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 8 + 29:
            raise ValueError((line_number, len(fields)))
        score1m, score200, base, q1, r1, q2, r2 = map(int, fields[:7])
        n = int(fields[7])
        words = list(map(int, fields[8:]))
        rows.append({
            "line_number": line_number,
            "score1m": score1m,
            "score200": score200,
            "base_score": base,
            "q1": q1,
            "r1": r1,
            "q2": q2,
            "r2": r2,
            "n": n,
            "words": words,
        })
    return rows


B = [comb(2*k+1, k) for k in range(60)]
SHIFTS = []
PAIRS = []
for c in range(60):
    for d in range(c, 60):
        SHIFTS.append(B[c] + B[d])
        PAIRS.append((c, d))


def live_indices(row):
    result = []
    for word_index, original in enumerate(row["words"]):
        word = original
        while word:
            low = word & -word
            bit = low.bit_length() - 1
            word ^= low
            index = 64 * word_index + bit
            if index < len(SHIFTS) and SHIFTS[index] <= row["n"]:
                result.append(index)
    if len(result) != row["score1m"]:
        raise AssertionError((row["n"], len(result), row["score1m"]))
    return result


def factor_one(value: int, timeout_seconds: int):
    completed = subprocess.run(
        ["factor", str(value)],
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
        check=True,
    )
    parsed_value, factors = parse_factor_output(completed.stdout)
    if parsed_value != value:
        raise AssertionError((parsed_value, value))
    return factors


def verify_candidate(row, timeout_seconds=90):
    obstructed = []
    timeouts = []
    errors = []
    for index in live_indices(row):
        pair = PAIRS[index]
        remainder = row["n"] - SHIFTS[index]
        try:
            factors = factor_one(remainder, timeout_seconds)
        except subprocess.TimeoutExpired:
            timeouts.append({"index": index, "pair": pair, "remainder": remainder})
            continue
        except Exception as exc:
            errors.append({"index": index, "pair": pair, "remainder": remainder, "error": repr(exc)})
            continue

        entry = {
            "index": index,
            "pair": pair,
            "remainder": remainder,
            "factors": {str(p): e for p, e in factors.items()},
        }
        if is_sum_of_two_squares(factors):
            return {
                **{k: v for k, v in row.items() if k != "words"},
                "status": "REJECTED",
                "witness": entry,
                "obstructed_checked": len(obstructed),
                "timeouts_before_witness": timeouts,
                "errors_before_witness": errors,
            }
        obstructed.append(entry)

    if not timeouts and not errors:
        return {
            **{k: v for k, v in row.items() if k != "words"},
            "status": "FOUND_COUNTEREXAMPLE",
            "certificates": obstructed,
        }
    return {
        **{k: v for k, v in row.items() if k != "words"},
        "status": "UNRESOLVED",
        "obstructed": obstructed,
        "timeouts": timeouts,
        "errors": errors,
    }


def main():
    if len(sys.argv) < 3:
        raise SystemExit("usage: verify_two_coord.py CANDIDATES OUTPUT [workers=2] [timeout=90]")
    source = Path(sys.argv[1])
    output = Path(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    timeout_seconds = int(sys.argv[4]) if len(sys.argv) > 4 else 90

    rows = parse_candidates(source)
    results = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(verify_candidate, row, timeout_seconds): row for row in rows}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                result["status"],
                result["score1m"],
                result["n"],
                flush=True,
            )

    results.sort(key=lambda item: (item["score1m"], item["n"]))
    output.write_text(json.dumps(results, indent=2))

    counts = Counter(item["status"] for item in results)
    print("SUMMARY", dict(counts))
    found = [item for item in results if item["status"] == "FOUND_COUNTEREXAMPLE"]
    if found:
        Path("FOUND_COUNTEREXAMPLE.json").write_text(json.dumps(found, indent=2))
        print("FOUND_VALUES", [item["n"] for item in found])


if __name__ == "__main__":
    main()
