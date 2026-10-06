#!/usr/bin/env python
"""
amplicon_pipeline.py

Species-agnostic amplicon sequencing pipeline.

Modes:
  nanopore-plate    : pooled Nanopore plate fastq(s) + plate-layout csv + known
                       inner barcodes.csv -> demux -> ONT track
  nanopore-index    : single Nanopore fastq + known inner barcodes.csv
                       (id,forward,reverse) -> demux -> ONT track
  illumina-index    : Illumina R1/R2 + I1/I2 index csv -> demux -> Illumina track
  discover-barcodes : single Nanopore fastq + a CANDIDATE barcode panel csv
                       -> reports which barcodes actually have signal, so you
                       can build the barcodes.csv used by the other two
                       Nanopore modes when the inner barcode set isn't known.

Downstream of demux, both platform tracks converge on:
  map -> split per amplicon (BED) -> [optional downsample]
  -> coverage -> variant call -> merge -> filter -> split snps/indels
  -> bcftools csq annotate -> DR/IR special-SNP lookup -> markdown report

--ploidy is REQUIRED (haploid|diploid) for any mode that calls variants.
No species/ploidy default is assumed.

CHANGES IN THIS VERSION:
  - Each inner-demuxed sample's BAM is now split into one BAM per amplicon
    (using the BED file's name column), because a single inner barcode can
    contain multiple pooled amplicons. All downstream steps (coverage,
    variant calling, merge) now operate on these (sample, amplicon) units,
    not on samples alone.
  - samclip removed from mapping (see amplicon_lib.py docstring for why).
  - New --downsample-depth flag (default: None / off). When set, each
    (sample, amplicon) BAM is downsampled to approximately that mean depth
    before variant calling. Off by default -- no behaviour change unless
    explicitly requested.
"""

import sys
import argparse
import amplicon_lib as alib
import amplicon_demux as ademux


def add_common_downstream_args(p):
    p.add_argument('--ref', required=True, help='Reference fasta')
    p.add_argument('--gff', required=True, help='GFF file for bcftools csq annotation')
    p.add_argument('--bed', required=True,
                    help='BED file with amplicon target regions. Column 4 (name) '
                         'is REQUIRED and must be unique per amplicon -- it is used '
                         'as the amplicon ID when splitting pooled samples.')
    p.add_argument('--ploidy', required=True, choices=['haploid', 'diploid'],
                    help='Required: ploidy assumption for variant calling. No default.')
    p.add_argument('--min-depth', type=int, default=10, help='Minimum FMT/DP to keep a genotype')
    p.add_argument('--min-qual', type=int, default=30, help='Minimum QUAL to keep a variant')
    p.add_argument('--dr-ir-snps', type=str, default=None,
                    help='Optional drug-resistance marker table (gene,aa_change[,linked_group,...])')
    p.add_argument('--ir-snps', type=str, default=None,
                    help='Optional insecticide-resistance marker table, same format as --dr-ir-snps. '
                         'Both can be supplied together; hits from each are reported separately.')
    p.add_argument('--threads', type=int, default=4, help='Threads for trimming/mapping')
    p.add_argument('--downsample-depth', type=int, default=None,
                    help='Optional: downsample each (sample, amplicon) BAM to ~this mean '
                         'depth before variant calling. Default: None (no downsampling).')
    p.add_argument('--vaf-mixed-low', type=float, default=0.05,
                    help='VAF below this = wildtype-dominant/likely noise (default 0.05). '
                         'NOT validated against your data -- review before trusting labels.')
    p.add_argument('--vaf-mixed-high', type=float, default=0.95,
                    help='VAF above this = mutant-fixed (default 0.95). Same caveat as --vaf-mixed-low.')
    p.add_argument('--medaka-model', type=str, default='r1041_e82_400bps_sup_v5.0.0',
                    help='Medaka model to use (ONT mode, --variant-caller medaka)')
    p.add_argument('--variant-caller', choices=['medaka', 'clair3'], default='medaka',
                    help='ONT variant caller to use. Default: medaka. (Illumina track always uses GATK.)')
    p.add_argument('--clair3-model', type=str, default=None,
                    help='Path to Clair3 model directory matched to your basecaller model/chemistry. '
                         'Required when --variant-caller clair3.')


def _amplicon_units_for_sample(sample, bed):
    """Returns the list of (unit_id, amplicon) pairs for a given sample, e.g. [('s1_ampA','ampA'), ...]."""
    amplicons = alib.load_amplicons(bed)
    return [(f"{sample}_{amp}", amp) for amp in amplicons]


def split_and_prepare_units(samples, bed, args):
    """
    For each sample's {sample}.bam, split into per-amplicon BAMs and optionally
    downsample. Returns the flat list of unit IDs (sample_amplicon strings)
    to be used by every downstream step from here on.
    """
    units = []
    bed_regions = {}
    with open(bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end, amp = fields[0], fields[1], fields[2], fields[3]
            bed_regions[amp] = f"{chrom}:{int(start) + 1}-{end}"

    for s in samples:
        bam = f"{s}.bam"
        out_bams = alib.split_bam_by_amplicon(s, bam, bed)
        for amp, unit_bam in out_bams.items():
            unit_id = f"{s}_{amp}"
            final_bam = alib.downsample_bam(
                unit_bam,
                target_depth=args.downsample_depth,
                bed_region=bed_regions[amp],
            )
            # normalize naming so downstream steps can always find {unit}.bam
            if final_bam != f"{unit_id}.bam":
                alib.run_cmd(f"cp {final_bam} {unit_id}.bam", verbose=0)
                alib.run_cmd(f"samtools index {unit_id}.bam", verbose=0)
            units.append(unit_id)
    return units


def run_downstream(units, args, platform):
    """
    platform: 'ont' or 'illumina'. Assumes per-unit (sample x amplicon) bams
    already exist as {unit}.bam (produced by split_and_prepare_units()).
    """
    if platform == "ont" and args.variant_caller == "clair3" and not args.clair3_model:
        sys.stderr.write(
            "[error] --variant-caller clair3 requires --clair3-model "
            "(path to a model dir matched to your basecaller/chemistry)\n"
        )
        sys.exit(1)

    alib.prep_reference(args.ref)

    for u in units:
        alib.calc_coverage(u, f"{u}.bam", args.bed)
        if platform == "ont":
            alib.call_variants_ont(u, f"{u}.bam", args.ref, args.bed, args.ploidy,
                                    caller=args.variant_caller,
                                    medaka_model=args.medaka_model,
                                    clair3_model=args.clair3_model,
                                    threads=args.threads)
        else:
            vcf = alib.call_variants_illumina(u, f"{u}.bam", args.ref, args.bed, args.ploidy)
            if vcf != f"{u}.vcf.gz":
                alib.run_cmd(f"cp {vcf} {u}.vcf.gz", verbose=0)
                alib.run_cmd(f"tabix -f {u}.vcf.gz", verbose=0)

    combined = alib.merge_samples(units)
    filtered = alib.filter_variants(combined, args.ref, args.min_depth, args.min_qual)
    snps, indels = alib.split_snps_indels(filtered)
    snps_csq = alib.annotate_csq(snps, args.ref, args.gff)
    indels_csq = alib.annotate_csq(indels, args.ref, args.gff)

    snps_table = alib.export_variant_table(snps_csq, "combined.snps.txt")
    indels_table = alib.export_variant_table(indels_csq, "combined.indels.txt")

    special_hits = []
    for marker_table in [args.dr_ir_snps, args.ir_snps]:
        if marker_table:
            # Pass BOTH snp and indel csq VCFs -- find_special_snp_hits matches
            # on gene+aa_change regardless of whether the underlying variant is
            # a SNP, indel, or multi-codon combined change, so no separate
            # handling is needed for indels here.
            special_hits += alib.find_special_snp_hits(
                [snps_csq, indels_csq], marker_table,
                mixed_low=args.vaf_mixed_low, mixed_high=args.vaf_mixed_high)

    alib.load_amplicon_depth(units, args.bed)

    report = alib.render_report(units, args.bed, special_hits, snps_table, indels_table)
    sys.stderr.write(f"\nDone. Report written to {report}\n")


# ---------------------------------------------------------------------------
# Mode: illumina-index
# ---------------------------------------------------------------------------

def main_illumina_index(args):
    samples = alib.load_samples(args.index_file)
    counts = ademux.illumina_index_demux(
        args.read1, args.read2, args.index_file,
        max_mismatches=args.max_mismatches,
        search_flipped_index=args.search_flipped_index,
    )
    sys.stderr.write(f"Demux read counts: {counts}\n")

    for s in samples:
        t1, t2 = alib.trim_illumina(s, f"{s}_1.fastq.gz", f"{s}_2.fastq.gz", args.threads)
        alib.map_illumina(s, t1, t2, args.ref, args.threads)

    units = split_and_prepare_units(samples, args.bed, args)
    run_downstream(units, args, platform="illumina")


# ---------------------------------------------------------------------------
# Mode: nanopore-direct (Dorado barcode == sample, no pooling)
# ---------------------------------------------------------------------------

def main_nanopore_direct(args):
    """
    No inner demux step at all -- Dorado's outer barcode already equals the
    sample. Just maps the manifest's fastq filenames to sample IDs and runs
    straight to trim -> map -> amplicon split -> variant call.
    """
    import csv
    import shutil

    samples = []
    with open(args.manifest) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample, fastq_name = row["sample"], row["fastq"]
            src = f"{args.fastq_dir}/{fastq_name}"
            dst = f"{sample}.fastq"
            if fastq_name.endswith(".gz"):
                alib.run_cmd(f"zcat {src} > {dst}", target=dst)
            else:
                shutil.copyfile(src, dst)
            samples.append(sample)

    sys.stderr.write(f"nanopore-direct: {len(samples)} samples, no inner demux applied\n")

    for s in samples:
        trimmed = alib.trim_ont(s, f"{s}.fastq", args.threads)
        alib.map_ont(s, trimmed, args.ref, args.threads)

    units = split_and_prepare_units(samples, args.bed, args)
    run_downstream(units, args, platform="ont")


# ---------------------------------------------------------------------------
# Mode: nanopore-index
# ---------------------------------------------------------------------------

def main_nanopore_index(args):
    counts = ademux.nanopore_index_demux(
        args.fastq, args.barcodes,
        max_mismatch=args.max_mismatch, edge_size=args.edge_size,
        log_prefix=args.log_prefix,
    )
    sys.stderr.write(f"Demux read counts: {counts}\n")
    samples = [bid for bid in counts if bid != "unassigned"]

    for s in samples:
        trimmed = alib.trim_ont(s, f"{s}.fastq", args.threads)
        alib.map_ont(s, trimmed, args.ref, args.threads)

    units = split_and_prepare_units(samples, args.bed, args)
    run_downstream(units, args, platform="ont")


# ---------------------------------------------------------------------------
# Mode: nanopore-plate
# ---------------------------------------------------------------------------

def main_nanopore_plate(args):
    output_rows = ademux.nanopore_plate_demux(
        args.plate_layout, args.barcodes, args.fastq_dir,
        max_mismatch=args.max_mismatch, edge_size=args.edge_size,
    )
    samples = [r["sample"] for r in output_rows]

    for s in samples:
        trimmed = alib.trim_ont(s, f"{s}.fastq", args.threads)
        alib.map_ont(s, trimmed, args.ref, args.threads)

    units = split_and_prepare_units(samples, args.bed, args)
    run_downstream(units, args, platform="ont")


# ---------------------------------------------------------------------------
# Mode: discover-barcodes
# ---------------------------------------------------------------------------

def main_discover_barcodes(args):
    results = ademux.discover_barcodes(
        args.fastq, args.candidate_barcodes,
        max_mismatch=args.max_mismatch, edge_size=args.edge_size,
        min_reads=args.min_reads, log_prefix=args.log_prefix,
    )
    outfile = ademux.write_discovery_report(results, args.outfile)
    sys.stderr.write(f"\nBarcode discovery complete. Report written to {outfile}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(description="Species-agnostic amplicon sequencing pipeline")
    sub = parser.add_subparsers(dest="mode", required=True)

    p_ont_direct = sub.add_parser(
        "nanopore-direct",
        help="Dorado barcode == sample already (no pooling, no inner demux needed)"
    )
    p_ont_direct.add_argument(
        "--fastq-dir", required=True,
        help="Directory of already-demuxed per-sample fastqs, one per Dorado barcode"
    )
    p_ont_direct.add_argument(
        "--manifest", required=True,
        help="CSV: sample,fastq -- maps sample ID to the fastq filename (within --fastq-dir) for that barcode"
    )
    add_common_downstream_args(p_ont_direct)
    p_ont_direct.set_defaults(func=main_nanopore_direct)

    p_ont_idx = sub.add_parser("nanopore-index", help="Single Nanopore fastq, known inner barcodes")
    p_ont_idx.add_argument("--fastq", required=True)
    p_ont_idx.add_argument("--barcodes", required=True, help="CSV: id,forward,reverse")
    p_ont_idx.add_argument("--max-mismatch", type=int, default=1)
    p_ont_idx.add_argument("--edge-size", type=int, default=150)
    p_ont_idx.add_argument("--log-prefix", default=None)
    add_common_downstream_args(p_ont_idx)
    p_ont_idx.set_defaults(func=main_nanopore_index)

    p_ont_plate = sub.add_parser("nanopore-plate", help="Pooled Nanopore plate fastqs + plate layout")
    p_ont_plate.add_argument("--plate-layout", required=True, help="CSV: well,sample,plate")
    p_ont_plate.add_argument("--barcodes", required=True, help="CSV: well/id,forward,reverse")
    p_ont_plate.add_argument("--fastq-dir", required=True)
    p_ont_plate.add_argument("--max-mismatch", type=int, default=1)
    p_ont_plate.add_argument("--edge-size", type=int, default=150)
    add_common_downstream_args(p_ont_plate)
    p_ont_plate.set_defaults(func=main_nanopore_plate)

    p_ill = sub.add_parser("illumina-index", help="Illumina R1/R2 + I1/I2 index csv")
    p_ill.add_argument("--read1", required=True)
    p_ill.add_argument("--read2", required=True)
    p_ill.add_argument("--index-file", required=True, help="CSV: id,i1,i2")
    p_ill.add_argument("--max-mismatches", type=int, default=1)
    p_ill.add_argument("--search-flipped-index", action="store_true")
    add_common_downstream_args(p_ill)
    p_ill.set_defaults(func=main_illumina_index)

    p_disc = sub.add_parser("discover-barcodes", help="Find which inner barcodes have signal")
    p_disc.add_argument("--fastq", required=True)
    p_disc.add_argument("--candidate-barcodes", required=True, help="CSV: id,forward,reverse")
    p_disc.add_argument("--max-mismatch", type=int, default=1)
    p_disc.add_argument("--edge-size", type=int, default=150)
    p_disc.add_argument("--min-reads", type=int, default=10)
    p_disc.add_argument("--log-prefix", default=None)
    p_disc.add_argument("--outfile", default="barcode_discovery_report.txt")
    p_disc.set_defaults(func=main_discover_barcodes)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
