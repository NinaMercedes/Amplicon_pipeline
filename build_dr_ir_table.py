#!/usr/bin/env python3
"""
build_dr_ir_table.py

Joins a MalariaProfiler-style protein-level resistance marker table
(Gene=PF3D7_xxx, Mutation=p.AsnXXTyr, Info=key=value;key=value...)
against a per-position SNP annotation file (chrom, pos, ref, alt,
effect_string) to produce a genomic-coordinate DR/IR table usable by
amplicon_lib.find_special_snp_hits().

Usage:
    python build_dr_ir_table.py \
        --markers variants__4_.csv \
        --annotations annotations.tsv \
        --out dr_ir_snps.csv

Matching logic:
    - Parses the 3-letter-style p.AsnXXTyr mutation string into
      (wt_aa_1letter, codon_number, mut_aa_1letter).
    - Parses each annotations.tsv effect field, which looks like:
        missense|GENE|PF3D7_ID|protein_coding|strand|72C>72S|403612T>A
      or, for multi-nucleotide codon changes:
        missense|GENE|PF3D7_ID|protein_coding|strand|76K>76T|403625A>C,missense|...|403625A>C+403626A>T
      (the row may carry the field for MULTIPLE alt alleles comma-separated,
       and a single allele may itself span multiple genomic positions via '+').
    - Matches a marker to every annotation row whose gene ID and codon-level
      aa change matches, regardless of how many genomic positions that change
      spans.
    - Multi-position changes are emitted as a single combined row with a
      '+'-joined pos/ref/alt, flagged for manual VCF-representation check
      since how your variant caller represents a 2-base substitution across
      adjacent positions (two SNPs vs one MNP) needs confirming.

This is a best-effort automatic join -- ALWAYS spot check a handful of rows
against the literature/db manually before trusting it in production, and
confirm REF/ALT strings against a real VCF from a known-positive sample.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

PMUT_RE = re.compile(r"^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$")


def parse_p_mutation(mutation_str):
    """p.Asn86Tyr -> ('N', 86, 'Y'). Returns None if unparseable."""
    m = PMUT_RE.match(mutation_str.strip())
    if not m:
        return None
    wt3, codon, mut3 = m.groups()
    wt1 = AA3_TO_1.get(wt3)
    mut1 = AA3_TO_1.get(mut3)
    if wt1 is None or mut1 is None:
        return None
    return wt1, int(codon), mut1


def parse_info(info_str):
    """type=drug_resistance;drug=chloroquine;source=... -> dict"""
    out = {}
    for kv in info_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_effect_field(effect_str, gene_id, target_codon, target_wt, target_mut):
    """
    Parse one row's effect column. Returns a list of (chrom_positions, refs, alts)
    combined-change tuples that match the target gene+codon+aa change.
    The field may contain multiple comma-separated effect descriptors. Some
    descriptors are duplicated with a leading '*' (canonical-transcript marker
    in this annotation format) -- strip it and dedupe. Some descriptors are
    bare back-references like '@403625' pointing at another row that carries
    the full multi-position change -- skip just that token, not the whole field.
    """
    matches = []
    if effect_str in (".", ""):
        return matches

    seen = set()
    for descriptor in effect_str.split(","):
        descriptor = descriptor.strip()
        if descriptor.startswith("@"):
            continue  # partial back-reference token, not a full descriptor
        descriptor = descriptor.lstrip("*")  # canonical-transcript marker, not semantically different
        if descriptor in seen:
            continue
        seen.add(descriptor)

        parts = descriptor.split("|")
        if len(parts) < 7:
            continue
        consequence, gene_name, gid, biotype, strand, aa_change_field, nt_change_field = parts[:7]
        if gid != gene_id:
            continue
        aa_m = re.match(r"^(\d+)([A-Za-z\*]+)>(\d+)([A-Za-z\*]+)$", aa_change_field)
        if not aa_m:
            continue
        codon_from, aa_from, codon_to, aa_to = aa_m.groups()
        if int(codon_from) != target_codon:
            continue
        if aa_from != target_wt or aa_to != target_mut:
            continue
        nt_changes = nt_change_field.split("+")
        positions, refs, alts = [], [], []
        for nt_change in nt_changes:
            nt_m = re.match(r"^(\d+)([A-Za-z])>([A-Za-z])$", nt_change)
            if not nt_m:
                continue
            pos, ref, alt = nt_m.groups()
            positions.append(int(pos))
            refs.append(ref)
            alts.append(alt)
        if positions:
            matches.append((positions, refs, alts))
    return matches


def load_annotations(path):
    """
    Returns dict: (gene_id) -> list of (chrom, pos, ref, alt, effect_str) rows,
    for fast filtering, since the file is only ~2500 rows.
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, ref, alt, effect = fields[:5]
            rows.append((chrom, int(pos), ref, alt, effect))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", required=True, help="MalariaProfiler-style CSV (Gene,Mutation,Info)")
    ap.add_argument("--annotations", required=True, help="Per-position SNP annotation TSV")
    ap.add_argument("--out", required=True)
    ap.add_argument("--unmatched-out", default=None,
                     help="Optional: write markers that couldn't be auto-matched here for manual review")
    args = ap.parse_args()

    annot_rows = load_annotations(args.annotations)

    out_rows = []
    unmatched = []

    with open(args.markers) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = row["Gene"].replace("PF3D7_", "PF3D7_")  # normalize if needed
            parsed = parse_p_mutation(row["Mutation"])
            info = parse_info(row["Info"])
            if parsed is None:
                unmatched.append({**row, "reason": "unparseable mutation string"})
                continue
            wt1, codon, mut1 = parsed

            found_any = False
            for chrom, pos, ref, alt, effect in annot_rows:
                matches = parse_effect_field(effect, gene_id, codon, wt1, mut1)
                for positions, refs, alts in matches:
                    found_any = True
                    out_rows.append({
                        "chrom": chrom,
                        "pos": "+".join(str(p) for p in positions),
                        "ref": "".join(refs) if len(positions) > 1 else refs[0],
                        "alt": "".join(alts) if len(positions) > 1 else alts[0],
                        "gene": row["Gene"],
                        "codon": codon,
                        "aa_change": f"{wt1}{codon}{mut1}",
                        "drug": info.get("drug", ""),
                        "who_status": info.get("type", ""),
                        "reference": info.get("source", ""),
                        "note": "multi-position change -- confirm VCF representation (MNP vs adjacent SNPs)"
                                if len(positions) > 1 else "",
                    })

            if not found_any:
                unmatched.append({**row, "reason": "no matching annotation row found for this gene+codon+change"})

    fieldnames = ["chrom", "pos", "ref", "alt", "gene", "codon", "aa_change",
                  "drug", "who_status", "reference", "note"]
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in out_rows:
            writer.writerow(r)

    sys.stderr.write(f"Matched {len(out_rows)} genomic rows from {len(out_rows)} marker hits.\n")
    sys.stderr.write(f"Unmatched markers: {len(unmatched)}\n")

    if unmatched:
        unmatched_path = args.unmatched_out or (args.out + ".unmatched.csv")
        with open(unmatched_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["Gene", "Mutation", "Info", "reason"])
            writer.writeheader()
            for u in unmatched:
                writer.writerow(u)
        sys.stderr.write(f"Unmatched markers written to {unmatched_path} -- review these manually.\n")


if __name__ == "__main__":
    main()
