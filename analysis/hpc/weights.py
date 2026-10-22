#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
from pathlib import Path


FG = {
    "prs_cvd_mortality": "I9_CVD_HARD",
    "prs_cad": "I9_CHD",
    "prs_stroke": "I9_STR",
    "prs_heart_failure": "I9_HEARTFAIL",
    "prs_atrial_fibrillation": "I9_AF",
    "prs_type2_diabetes": "T2D",
    "prs_depression": "F5_DEPRESSIO",
    "prs_anxiety": "F5_ALLANXIOUS",
    "prs_dementia": "F5_DEMENTIA_INCLAVO",
}


def opener(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def prepare_finngen(raw, out, traits):
    for trait in traits:
        code = FG[trait]
        src = Path(raw) / f"finngen_R12_{code}.gz"
        dst = Path(out) / f"{trait}.txt"
        with opener(src) as fi, opener(dst, "wt") as fo:
            rd = csv.DictReader(fi, delimiter="\t")
            wr = csv.writer(fo, delimiter="\t", lineterminator="\n")
            wr.writerow(["SNP", "A1", "A2", "BETA", "P"])
            seen = set()
            for x in rd:
                snp = x["rsids"].split(",")[0]
                if not snp.startswith("rs") or snp in seen:
                    continue
                if len(x["alt"]) != 1 or len(x["ref"]) != 1:
                    continue
                seen.add(snp)
                wr.writerow([snp, x["alt"].upper(), x["ref"].upper(), x["beta"], x["pval"]])


def prepare_income(raw, out):
    src = Path(raw) / "income.txt"
    dst = Path(out) / "prs_income.txt"
    with opener(src) as fi, opener(dst, "wt") as fo:
        rd = csv.DictReader(fi, delimiter="\t")
        wr = csv.writer(fo, delimiter="\t", lineterminator="\n")
        wr.writerow(["SNP", "A1", "A2", "BETA", "P"])
        seen = set()
        for x in rd:
            snp = x["SNP"]
            if not snp.startswith("rs") or snp in seen:
                continue
            if len(x["A1"]) != 1 or len(x["A2"]) != 1:
                continue
            seen.add(snp)
            wr.writerow([snp, x["A1"].upper(), x["A2"].upper(), x["BETA"], x["P"]])


def snp_map(ref):
    path = Path(ref) / "snpinfo_1kg_hm3"
    ans = {}
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for x in rd:
            ans[(x["CHR"], x["BP"])] = (x["SNP"], x["A1"].upper(), x["A2"].upper())
    return ans


def write_bim(ref):
    src = Path(ref) / "snpinfo_1kg_hm3"
    dst = Path(ref).parent / "hm3.bim"
    with open(src) as fi, open(dst, "w") as fo:
        rd = csv.DictReader(fi, delimiter="\t")
        for x in rd:
            fo.write(f'{x["CHR"]}\t{x["SNP"]}\t0\t{x["BP"]}\t{x["A1"]}\t{x["A2"]}\n')


def prepare_longevity(raw, ref, out):
    src = Path(raw) / "longevity.gz"
    dst = Path(out) / "prs_all_cause_mortality.txt"
    lookup = snp_map(ref)
    with opener(src) as fi, opener(dst, "wt") as fo:
        rd = csv.DictReader(fi, delimiter="\t")
        wr = csv.writer(fo, delimiter="\t", lineterminator="\n")
        wr.writerow(["SNP", "A1", "A2", "BETA", "P"])
        seen = set()
        for x in rd:
            key = (x["Chr"], x["Position"])
            if key not in lookup:
                continue
            snp, r1, r2 = lookup[key]
            a1, a2 = x["EA"].upper(), x["NEA"].upper()
            if {a1, a2} != {r1, r2} or snp in seen:
                continue
            seen.add(snp)
            wr.writerow([snp, a1, a2, -float(x["Beta"]), x["P-value"]])


def prepare_education(raw, out):
    src = Path(raw) / "education.txt"
    dst = Path(out) / "prs_education.txt"
    with opener(src) as fi, opener(dst, "wt") as fo:
        rd = csv.DictReader(fi, delimiter="\t")
        wr = csv.writer(fo, delimiter="\t", lineterminator="\n")
        wr.writerow(["SNP", "A1", "BETA"])
        for x in rd:
            if x["MarkerName"].startswith("rs"):
                wr.writerow([x["MarkerName"], x["A1"].upper(), x["Beta"]])


def weights(src, out, traits):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for trait in traits:
        dst = out / f"{trait}.score"
        with open(dst, "w") as fo:
            fo.write(f"SNP\tA1\t{trait}\n")
            for chrom in range(1, 23):
                pat = Path(src).glob(f"{trait}*_chr{chrom}.txt")
                files = list(pat)
                if len(files) != 1:
                    raise RuntimeError(f"Expected one PRS-CS file for {trait} chr{chrom}, found {len(files)}")
                with open(files[0]) as fi:
                    for line in fi:
                        z = line.split()
                        if len(z) >= 6:
                            fo.write(f"{z[1]}\t{z[3]}\t{z[5]}\n")


def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers(dest="cmd", required=True)
    q = sp.add_parser("prepare")
    q.add_argument("--raw", required=True)
    q.add_argument("--ref", required=True)
    q.add_argument("--out", required=True)
    q = sp.add_parser("weights")
    q.add_argument("--src", required=True)
    q.add_argument("--out", required=True)
    a = p.parse_args()

    if a.cmd == "prepare":
        Path(a.out).mkdir(parents=True, exist_ok=True)
        write_bim(a.ref)
        prepare_finngen(a.raw, a.out, list(FG))
        prepare_income(a.raw, a.out)
        prepare_longevity(a.raw, a.ref, a.out)
        prepare_education(a.raw, a.out)
    else:
        traits = ["prs_income", "prs_all_cause_mortality", *FG]
        weights(a.src, a.out, traits)


if __name__ == "__main__":
    main()
