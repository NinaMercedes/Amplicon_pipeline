"""
amplicon_lib.py

Shared library functions for the species-agnostic amplicon pipeline.
Covers: trimming, mapping, per-amplicon BAM splitting, optional downsampling,
coverage, variant calling, filtering, annotation (bcftools csq), DR/IR
special-SNP lookup, and report rendering.

Design notes (updated):
- ONT track:      fastplong (trim) -> minimap2 map-ont -> sort/index
                   -> split per amplicon (BED) -> [optional downsample]
                   -> medaka_variant (per sample x amplicon)
- Illumina track: fastp (trim)     -> bwa mem          -> sort/index
                   -> split per amplicon (BED) -> [optional downsample]
                   -> GATK HaplotypeCaller (per sample x amplicon)

Changes from previous version:
  1. samclip removed from the mapping pipe entirely. Amplicon reads
     routinely soft-clip at primer/adapter boundaries -- that's expected,
     not a sign of misalignment -- so an edge-soft-clip filter tuned for
     shotgun/WGS data was dropping legitimate on-target reads, especially
     near amplicon ends. Region scoping to the BED file during variant
     calling does the job samclip was being asked to do here.
  2. New split_bam_by_amplicon(): splits a per-sample BAM into one BAM per
     (sample, amplicon) pair using the BED file's 4th column (name) as the
     amplicon ID. This is the fix for "one inner barcode can contain
     multiple pooled amplicons" -- everything downstream of mapping now
     operates on (sample, amplicon) units, not just sample.
  3. New optional downsample_bam(): per (sample, amplicon) unit, downsamples
     to a target read depth before variant calling, off by default.
"""

import sys
import os
import subprocess as sp
import csv
import random
from collections import defaultdict
import gzip


def run_cmd(cmd, verbose=1, target=None, terminate_on_error=True):
    """
    Run a shell command, with optional skip-if-target-exists and
    consistent error reporting.
    """
    if target and os.path.isfile(target):
        if verbose:
            sys.stderr.write(f"[skip] {target} exists, skipping: {cmd}\n")
        return
    if verbose:
        sys.stderr.write(f"[run] {cmd}\n")
    result = sp.run(cmd, shell=True)
    if result.returncode != 0:
        msg = f"[error] command failed (exit {result.returncode}): {cmd}\n"
        if terminate_on_error:
            sys.stderr.write(msg)
            sys.exit(result.returncode)
        else:
            sys.stderr.write(msg)
    return result.returncode


def load_samples(index_file):
    """Load sample IDs from an index/barcode CSV (expects an 'id' or 'sample' column)."""
    samples = []
    with open(index_file) as f:
        reader = csv.DictReader(f)
        id_col = "id" if "id" in reader.fieldnames else "sample"
        for row in reader:
            samples.append(row[id_col])
    return samples


def load_amplicons(bed):
    """
    Read amplicon names from a BED file's 4th column.
    Expects: chrom, start, end, name[, ...]
    Raises if any row is missing a name, since amplicon splitting depends on it.
    """
    amplicons = []
    with open(bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"BED file {bed} is missing a 4th (name) column on line: {line!r}. "
                    "Amplicon splitting requires a unique name per region."
                )
            amplicons.append(fields[3])
    if len(amplicons) != len(set(amplicons)):
        dupes = sorted({a for a in amplicons if amplicons.count(a) > 1})
        raise ValueError(f"BED file {bed} has duplicate amplicon names: {dupes}")
    return amplicons


# ---------------------------------------------------------------------------
# Reference prep
# ---------------------------------------------------------------------------

def prep_reference(ref):
    run_cmd(f"samtools faidx {ref}", target=f"{ref}.fai")


# ---------------------------------------------------------------------------
# Trimming
# ---------------------------------------------------------------------------

def trim_ont(sample, fastq_in, threads=4):
    """Adapter/quality trim a single ONT sample with fastplong."""
    out = f"{sample}.trimmed.fastq.gz"
    run_cmd(
        f"fastplong -i {fastq_in} -o {out} "
        f"--thread {threads} "
        f"-j {sample}.fastplong.json -h {sample}.fastplong.html",
        target=out
    )
    return out


def trim_illumina(sample, r1, r2, threads=4):
    """Adapter/quality trim a PE Illumina sample with fastp."""
    out1 = f"{sample}_1.trimmed.fastq.gz"
    out2 = f"{sample}_2.trimmed.fastq.gz"
    run_cmd(
        f"fastp -i {r1} -I {r2} -o {out1} -O {out2} "
        f"--thread {threads} "
        f"-j {sample}.fastp.json -h {sample}.fastp.html",
        target=out1
    )
    return out1, out2


# ---------------------------------------------------------------------------
# Mapping (samclip removed -- see module docstring point 1)
# ---------------------------------------------------------------------------

def map_ont(sample, fastq, ref, threads=10):
    bam = f"{sample}.bam"
    run_cmd(
        f"minimap2 -ax map-ont -t {threads} "
        f"-R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:nanopore' "
        f"{ref} {fastq} | samtools sort -o {bam} -",
        target=bam
    )
    run_cmd(f"samtools index {bam}", target=f"{bam}.bai")
    return bam


def map_illumina(sample, r1, r2, ref, threads=10):
    bam = f"{sample}.bam"
    run_cmd(
        f"bwa mem -t {threads} "
        f"-R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:Illumina' "
        f"{ref} {r1} {r2} | samtools sort -o {bam} -",
        target=bam
    )
    run_cmd(f"samtools index {bam}", target=f"{bam}.bai")
    return bam


# ---------------------------------------------------------------------------
# Per-amplicon BAM splitting (NEW -- fix for pooled amplicons within one
# inner barcode)
# ---------------------------------------------------------------------------

def split_bam_by_amplicon(sample, bam, bed):
    """
    Split a per-sample BAM into one BAM per amplicon, using each BED region's
    name (4th column) as the amplicon ID. Returns a dict: {amplicon_id: bam_path}.

    Reads overlapping multiple amplicon regions (e.g. tiled/overlapping schemes)
    will appear in more than one output BAM -- that's intentional, matching how
    region-restricted variant calling already treats overlapping BED entries.
    """
    out_bams = {}
    with open(bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end, amplicon = fields[0], fields[1], fields[2], fields[3]
            unit = f"{sample}_{amplicon}"
            out_bam = f"{unit}.bam"
            region = f"{chrom}:{int(start) + 1}-{end}"  # BED is 0-based, samtools region is 1-based
            run_cmd(
                f"samtools view -b -F 0x900 {bam} {region} "  # drop secondary/supplementary
                f"| samtools sort -o {out_bam} -",
                target=out_bam
            )
            run_cmd(f"samtools index {out_bam}", target=f"{out_bam}.bai")
            out_bams[amplicon] = out_bam
    return out_bams


# ---------------------------------------------------------------------------
# Optional downsampling (NEW -- off by default)
# ---------------------------------------------------------------------------

def downsample_bam(unit_bam, target_depth=None, bed_region=None, seed=42):
    """
    Optionally downsample a (sample x amplicon) BAM to a target mean depth
    over its amplicon region, using samtools view -s (fraction-based).

    target_depth: int or None. If None (default), this is a no-op and the
                  input BAM path is returned unchanged -- no downsampling
                  happens unless explicitly requested via --downsample-depth.
    bed_region:   "chrom:start-end" string used to measure current mean
                  depth via `samtools depth`. Required if target_depth is set.

    Uses a fixed seed by default for reproducibility across reruns.
    """
    if target_depth is None:
        return unit_bam

    if bed_region is None:
        raise ValueError("downsample_bam requires bed_region when target_depth is set")

    # Measure current mean depth over the amplicon region
    depth_cmd = f"samtools depth -a -r {bed_region} {unit_bam}"
    result = sp.run(depth_cmd, shell=True, capture_output=True, text=True)
    depths = [int(line.split("\t")[2]) for line in result.stdout.strip().split("\n") if line]
    current_depth = (sum(depths) / len(depths)) if depths else 0

    if current_depth <= target_depth:
        sys.stderr.write(
            f"[downsample] {unit_bam}: current depth {current_depth:.1f} <= "
            f"target {target_depth}, skipping downsample\n"
        )
        return unit_bam

    fraction = target_depth / current_depth
    out_bam = unit_bam.replace(".bam", f".ds{target_depth}.bam")
    run_cmd(
        f"samtools view -b -s {seed}.{int(fraction * 1000):03d} {unit_bam} > {out_bam}",
        target=out_bam
    )
    run_cmd(f"samtools index {out_bam}", target=f"{out_bam}.bai")
    sys.stderr.write(
        f"[downsample] {unit_bam}: {current_depth:.1f}x -> ~{target_depth}x "
        f"(fraction {fraction:.3f})\n"
    )
    return out_bam


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------

def calc_coverage(unit, bam, bed):
    """unit is now a (sample, amplicon) identifier string, e.g. 'sample1_ampA'."""
    run_cmd(
        f"mosdepth -x -b {bed} {unit} --thresholds 1,10,20,30 {bam}",
        target=f"{unit}.thresholds.bed.gz"
    )
    run_cmd(
        f"bedtools coverage -a {bed} -b {bam} -mean > {unit}_coverage_mean.txt",
        target=f"{unit}_coverage_mean.txt"
    )


def load_amplicon_depth(units, bed):
    """
    Build a per-unit, per-position depth dict restricted to the amplicon
    bed regions, from mosdepth's per-base.bed.gz output. 'units' are
    (sample, amplicon) identifiers.
    """
    bedlines = [l.strip().split() for l in open(bed)]

    def overlap_bedlines(a):
        overlaps = []
        for b in bedlines:
            if b[0] == a[0]:
                lo = max(int(a[1]), int(b[1]))
                hi = min(int(a[2]), int(b[2]))
                if hi > lo:
                    overlaps.append([b[0], lo, hi])
        return overlaps

    dp = defaultdict(dict)
    for u in units:
        path = f"{u}.per-base.bed.gz"
        if not os.path.isfile(path):
            continue
        for l in gzip.open(path):
            row = l.decode().strip().split()
            for ov in overlap_bedlines(row):
                for pos in range(ov[1], ov[2]):
                    dp[u][(row[0], pos)] = int(row[3])
    return dp


# ---------------------------------------------------------------------------
# Variant calling (now per sample x amplicon unit)
# ---------------------------------------------------------------------------

def call_variants_ont_medaka(unit, bam, ref, bed, medaka_model, ploidy):
    """
    medaka_variant for a single (sample, amplicon) BAM, restricted to the
    amplicon bed. ploidy: 'haploid' or 'diploid' (explicit, no default upstream).
    """
    outdir = f"{unit}_medaka"
    ploidy_flag = "--haploid" if ploidy == "haploid" else ""
    run_cmd(
        f"medaka_variant -i {bam} -f {ref} -o {outdir} "
        f"-m {medaka_model} -r {bed} {ploidy_flag}",
        target=f"{outdir}/round_1.vcf.gz"
    )
    vcf = f"{unit}.vcf.gz"
    run_cmd(f"cp {outdir}/round_1.vcf.gz {vcf}", target=vcf)
    run_cmd(f"tabix -f {vcf}")
    return vcf


def call_variants_ont_clair3(unit, bam, ref, bed, clair3_model, ploidy, threads=4):
    """
    Clair3 for a single (sample, amplicon) BAM, restricted to the amplicon bed.

    clair3_model: path to a Clair3 model directory matched to your basecaller
    model/flow cell chemistry (e.g. r1041_e82_400bps_sup_v500). Required --
    there is no sensible default, since the wrong model silently degrades
    accuracy rather than erroring.

    ploidy: 'haploid' or 'diploid'. Clair3 doesn't take an explicit ploidy
    flag the way medaka/GATK do. For 'haploid' we pass --no_phasing_for_fa
    and rely on the existing ploidy-aware downstream filter to handle any
    heterozygous calls -- this matches common amplicon/haploid-pathogen
    practice (e.g. Plasmodium) but should be validated against known truth
    samples before trusting silently.
    """
    outdir = f"{unit}_clair3"
    haploid_flags = "--no_phasing_for_fa" if ploidy == "haploid" else ""
    run_cmd(
        f"run_clair3.sh "
        f"--bam_fn={bam} --ref_fn={ref} --threads={threads} "
        f"--platform=ont --model_path={clair3_model} "
        f"--bed_fn={bed} --output={outdir} "
        f"--include_all_ctgs {haploid_flags}",
        target=f"{outdir}/merge_output.vcf.gz"
    )
    vcf = f"{unit}.vcf.gz"
    run_cmd(f"cp {outdir}/merge_output.vcf.gz {vcf}", target=vcf)
    run_cmd(f"tabix -f {vcf}")
    return vcf


def call_variants_ont(unit, bam, ref, bed, ploidy, caller="medaka",
                       medaka_model=None, clair3_model=None, threads=4):
    """
    Dispatches to the selected ONT variant caller.
    caller: 'medaka' (default) or 'clair3'.
    """
    if caller == "medaka":
        if not medaka_model:
            raise ValueError("medaka_model is required when --variant-caller medaka")
        return call_variants_ont_medaka(unit, bam, ref, bed, medaka_model, ploidy)
    elif caller == "clair3":
        if not clair3_model:
            raise ValueError("clair3_model is required when --variant-caller clair3")
        return call_variants_ont_clair3(unit, bam, ref, bed, clair3_model, ploidy, threads=threads)
    else:
        raise ValueError(f"Unknown ONT variant caller: {caller!r} (expected 'medaka' or 'clair3')")


def call_variants_illumina(unit, bam, ref, bed, ploidy):
    """
    GATK HaplotypeCaller for a single (sample, amplicon) BAM, restricted to
    the amplicon bed. ploidy: 'haploid' (-ploidy 1) or 'diploid' (-ploidy 2).
    """
    ploidy_n = 1 if ploidy == "haploid" else 2
    vcf = f"{unit}.gatk.vcf.gz"
    run_cmd(
        f"gatk HaplotypeCaller -R {ref} -L {bed} -I {bam} "
        f"-ploidy {ploidy_n} -O {vcf}",
        target=vcf
    )
    return vcf


# ---------------------------------------------------------------------------
# Merge / filter / annotate
# ---------------------------------------------------------------------------

def merge_samples(units, outfile="combined.vcf.gz"):
    """units are (sample, amplicon) identifiers; each gets its own VCF column."""
    with open("unit_vcf_list.txt", "w") as O:
        for u in units:
            O.write(f"{u}.vcf.gz\n")
    run_cmd(f"bcftools merge -l unit_vcf_list.txt -Oz -o {outfile}", target=outfile)
    run_cmd(f"tabix -f {outfile}")
    return outfile


def filter_variants(vcf_in, ref, min_depth=10, min_qual=30):
    """DP + QUAL filter, then normalize/split multiallelics."""
    tmp_vcf = "tmp.filtered.vcf.gz"
    run_cmd(
        f"bcftools filter -i 'FMT/DP>{min_depth}' -S . {vcf_in} -Oz -o {tmp_vcf}",
        target=tmp_vcf
    )
    qual_vcf = "tmp.qualfiltered.vcf.gz"
    run_cmd(
        f"bcftools view -i 'QUAL>{min_qual}' {tmp_vcf} -Oz -o {qual_vcf}",
        target=qual_vcf
    )
    norm_vcf = "combined.filtered.norm.vcf.gz"
    run_cmd(
        f"bcftools norm -f {ref} -m -any {qual_vcf} -Oz -o {norm_vcf}",
        target=norm_vcf
    )
    run_cmd(f"tabix -f {norm_vcf}")
    return norm_vcf


def split_snps_indels(vcf_in):
    snps = "combined.snps.vcf.gz"
    indels = "combined.indels.vcf.gz"
    run_cmd(f"bcftools view -v snps {vcf_in} -Oz -o {snps}", target=snps)
    run_cmd(f"bcftools view -v indels {vcf_in} -Oz -o {indels}", target=indels)
    run_cmd(f"tabix -f {snps}")
    run_cmd(f"tabix -f {indels}")
    return snps, indels


def annotate_csq(vcf_in, ref, gff):
    out = vcf_in.replace(".vcf.gz", ".csq.vcf.gz")
    run_cmd(
        f"bcftools csq -p a -f {ref} -g {gff} {vcf_in} -Oz -o {out}",
        target=out
    )
    run_cmd(f"tabix -f {out}")
    return out


def export_variant_table(vcf_in, outfile):
    run_cmd(
        f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t[%GT\\t%DP\\t%AD\\t]\\n' "
        f"{vcf_in} > {outfile}",
        target=outfile
    )
    return outfile


# ---------------------------------------------------------------------------
# DR/IR special-SNP lookup
# ---------------------------------------------------------------------------

def extract_csq_annotations(vcf_path):
    """
    Parse a bcftools-csq-annotated VCF and return one row per (sample, variant,
    consequence) triple:
      {sample, chrom, pos, ref, alt, gene, aa_change, gt, ad, dp, vaf}

    vaf (variant allele fraction = alt depth / total depth) matters
    specifically because Pf is called diploid here to capture clonal
    heterogeneity within a single blood-draw sample, not true biological
    diploidy. A 0/1 genotype therefore means "this sample's infection is a
    mixture of parasite clones, some with and some without this allele" --
    the GT alone doesn't say what fraction of the infection carries the
    resistance/IR allele, which is exactly the number that matters for
    interpreting a mixed-infection or partial-resistance result. VAF
    surfaces that fraction directly from the read counts.
    """
    query = (
        f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BCSQ\\t"
        f"[%SAMPLE=%GT,%AD;]\\n' {vcf_path}"
    )
    out = sp.run(query, shell=True, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stderr.write(f"[warn] extract_csq_annotations: bcftools query failed for {vcf_path}: {out.stderr}\n")
        return []

    rows = []
    for line in out.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        chrom, pos, ref, alt, bcsq, sample_gts = fields[:6]
        if bcsq in (".", ""):
            continue

        sample_info = {}  # sample -> (gt, ad_ref, ad_alt, dp, vaf)
        for tok in sample_gts.strip(";").split(";"):
            if "=" not in tok:
                continue
            s, rest = tok.split("=", 1)
            gt, _, ad_str = rest.partition(",")
            ad_ref, ad_alt, dp, vaf = "NA", "NA", "NA", "NA"
            ad_parts = ad_str.split(",") if ad_str else []
            if len(ad_parts) >= 2:
                try:
                    ad_ref, ad_alt = int(ad_parts[0]), int(ad_parts[1])
                    dp = ad_ref + ad_alt
                    vaf = round(ad_alt / dp, 3) if dp > 0 else 0.0
                except ValueError:
                    pass
            sample_info[s] = (gt, ad_ref, ad_alt, dp, vaf)

        for descriptor in bcsq.split(","):
            descriptor = descriptor.lstrip("*").strip()
            parts = descriptor.split("|")
            if len(parts) < 6:
                continue
            consequence, gene, transcript, biotype, strand, aa_change_field = parts[:6]
            aa_m = re.match(r"^(\d+)([A-Za-z\*]+)>(\d+)([A-Za-z\*]+)$", aa_change_field)
            if not aa_m:
                continue
            codon, aa_from, _, aa_to = aa_m.groups()
            aa_change = f"{aa_from}{codon}{aa_to}"
            for sample, (gt, ad_ref, ad_alt, dp, vaf) in sample_info.items():
                if gt in ("0/0", "0|0", ".", "./.", ".|."):
                    continue  # homozygous reference or missing -- not a hit
                rows.append({
                    "sample": sample, "chrom": chrom, "pos": int(pos),
                    "ref": ref, "alt": alt, "gene": gene,
                    "aa_change": aa_change, "gt": gt,
                    "ad_ref": ad_ref, "ad_alt": ad_alt, "dp": dp, "vaf": vaf,
                })
    return rows


def load_marker_table(marker_csv):
    """
    Loads a DR/IR marker table keyed by gene + aa_change (chrom/pos/ref/alt
    columns are optional and ignored for matching -- this is deliberate,
    since exact coordinates are often unconfirmed while the literature's
    gene+codon association is solid).

    Required columns: gene, aa_change
    Optional: linked_group (markers sharing a non-empty linked_group are only
    reported as a hit when EVERY aa_change in that group is present in the
    same sample -- this is how the crt 403618/403622 indel pair and the kdr
    L1014F+N1575Y super-kdr combo are both handled, with no special-casing
    needed for indels vs SNPs vs multi-codon changes).
    Any other columns (drug, insecticide_class, who_status, reference, note,
    marker_type, ...) are carried through into the report untouched.
    """
    markers = []
    with open(marker_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row.get("gene") or not row.get("aa_change"):
                continue
            markers.append(row)
    return markers


def classify_vaf(vaf, mixed_low=0.05, mixed_high=0.95):
    """
    Bands a variant allele fraction into a human-readable call classification
    for clonal-heterogeneity reporting (e.g. diploid-called P. falciparum,
    where heterozygous calls represent mixed-clone infections rather than
    true biological heterozygosity).

    Thresholds are a judgment call -- the defaults here are NOT validated
    against your data. Review/adjust mixed_low and mixed_high before trusting
    the label; the underlying VAF number is always reported alongside it so
    you can re-judge it yourself regardless of where these sit.
    """
    if vaf == "NA" or vaf is None:
        return "unknown (no AD)"
    if vaf < mixed_low:
        return "wildtype-dominant (check vs noise)"
    if vaf > mixed_high:
        return "mutant-fixed"
    return "mixed infection"


def find_special_snp_hits(units_or_vcf_paths, marker_csv, mixed_low=0.05, mixed_high=0.95):
    """
    Gene+aa_change based DR/IR marker detection. Works identically for DR
    (Pf) and IR (Anopheles) marker tables, and for SNPs, indels, and
    multi-codon "weird combination" changes alike, since matching happens
    on the bcftools-csq-derived aa_change string rather than on genomic
    coordinates.

    units_or_vcf_paths: list of per-unit VCF paths to scan (e.g. the
    {unit}.csq.vcf.gz files produced by annotate_csq() for both snps and
    indels -- pass both lists concatenated to catch everything).
    marker_csv: a DR or IR marker table (see load_marker_table docstring).
    mixed_low/mixed_high: VAF thresholds passed to classify_vaf() -- see that
    function's docstring, these are not validated defaults.

    Returns a list of hit dicts, one per (sample, linked_group-or-single-marker)
    that fully matched, with all marker metadata carried through, plus
    ad_ref/ad_alt/dp/vaf/call_class for clonal-heterogeneity interpretation.
    """
    markers = load_marker_table(marker_csv)

    # Index markers by (gene, aa_change) for fast lookup, and separately
    # group any markers sharing a non-empty linked_group.
    marker_index = defaultdict(list)
    linked_groups = defaultdict(list)
    for m in markers:
        marker_index[(m["gene"], m["aa_change"])].append(m)
        if m.get("linked_group"):
            linked_groups[m["linked_group"]].append(m)

    # Collect every (sample, gene, aa_change) observed across all VCFs.
    observed = defaultdict(set)  # sample -> set of (gene, aa_change)
    observed_rows = defaultdict(dict)  # (sample, gene, aa_change) -> annotation row

    for vcf_path in units_or_vcf_paths:
        for row in extract_csq_annotations(vcf_path):
            key = (row["gene"], row["aa_change"])
            observed[row["sample"]].add(key)
            observed_rows[(row["sample"], row["gene"], row["aa_change"])] = row

    hits = []
    reported = set()  # (sample, group_key) already emitted, avoid duplicates

    for sample, keys_present in observed.items():
        # Standalone (non-linked) markers: report each one directly.
        for (gene, aa_change) in keys_present:
            for m in marker_index.get((gene, aa_change), []):
                if m.get("linked_group"):
                    continue  # handled in the linked-group pass below
                report_key = (sample, gene, aa_change)
                if report_key in reported:
                    continue
                reported.add(report_key)
                row = observed_rows[(sample, gene, aa_change)]
                hits.append({
                    **m, "sample": sample, "chrom": row["chrom"],
                    "pos": row["pos"], "gt": row["gt"],
                    "ad_ref": row["ad_ref"], "ad_alt": row["ad_alt"], "dp": row["dp"],
                    "vaf": row["vaf"],
                    "call_class": classify_vaf(row["vaf"], mixed_low, mixed_high),
                })

        # Linked groups: only report if EVERY member's (gene, aa_change) is
        # present in this sample. This is what makes the crt indel pair and
        # the kdr+super-kdr combo work without special-casing indels at all --
        # it's the same logic regardless of how many genomic positions or
        # nucleotides each member's underlying variant actually spans.
        for group_name, group_markers in linked_groups.items():
            group_keys = [(m["gene"], m["aa_change"]) for m in group_markers]
            if all(k in keys_present for k in group_keys):
                report_key = (sample, "GROUP", group_name)
                if report_key in reported:
                    continue
                reported.add(report_key)
                rows = [observed_rows[(sample, g, a)] for (g, a) in group_keys]
                vafs = [r["vaf"] for r in rows if r["vaf"] != "NA"]
                # For a genuinely linked indel pair, every component's VAF
                # should be near-identical (same underlying reads/haplotype) --
                # report the min as a conservative estimate, and flag if the
                # components disagree, since that would mean they're NOT
                # actually as tightly linked in this sample as expected.
                min_vaf = min(vafs) if vafs else "NA"
                max_vaf = max(vafs) if vafs else "NA"
                vaf_note = ""
                if len(vafs) > 1 and (max_vaf - min_vaf) > 0.10:
                    vaf_note = (f" [WARNING: linked-group component VAFs disagree by "
                                f"{max_vaf - min_vaf:.2f} -- check if these are really linked in this sample]")
                hits.append({
                    **group_markers[0],
                    "sample": sample,
                    "chrom": rows[0]["chrom"],
                    "pos": "+".join(str(r["pos"]) for r in rows),
                    "gt": "+".join(r["gt"] for r in rows),
                    "aa_change": "+".join(m["aa_change"] for m in group_markers),
                    "ad_ref": "+".join(str(r["ad_ref"]) for r in rows),
                    "ad_alt": "+".join(str(r["ad_alt"]) for r in rows),
                    "dp": "+".join(str(r["dp"]) for r in rows),
                    "vaf": min_vaf,
                    "call_class": classify_vaf(min_vaf, mixed_low, mixed_high) if min_vaf != "NA" else "unknown (no AD)",
                    "note": (group_markers[0].get("note", "") +
                             f" [linked group '{group_name}': all {len(group_markers)} "
                             f"components confirmed present in this sample]" + vaf_note),
                })

    return hits


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------

def render_report(units, bed, special_hits, snps_table_file, indels_table_file, outfile="report.md"):
    lines = ["# Amplicon Pipeline Report\n"]
    lines.append(f"\n## Units processed (sample x amplicon): {len(units)}\n")
    for u in units:
        lines.append(f"- {u}\n")

    lines.append("\n## Variant counts\n")
    for label, f in [("SNPs", snps_table_file), ("Indels", indels_table_file)]:
        n = sum(1 for _ in open(f)) if os.path.isfile(f) else 0
        lines.append(f"- **{label}:** {n} calls\n")

    lines.append("\n## Drug / insecticide resistance marker hits\n")
    if special_hits:
        lines.append("| Sample | Gene | AA change | Chrom | Pos | GT | AD (ref,alt) | DP | VAF | Call | Drug/Insecticide | Status | Note |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for h in special_hits:
            drug = h.get("drug") or h.get("insecticide_class") or ""
            status = h.get("who_status") or h.get("marker_type") or ""
            vaf = h.get("vaf", "NA")
            vaf_str = f"{vaf:.1%}" if isinstance(vaf, float) else str(vaf)
            lines.append(
                f"| {h.get('sample','')} | {h.get('gene','')} | {h.get('aa_change','')} | "
                f"{h.get('chrom','')} | {h.get('pos','')} | {h.get('gt','')} | "
                f"{h.get('ad_ref','')},{h.get('ad_alt','')} | {h.get('dp','')} | {vaf_str} | "
                f"{h.get('call_class','')} | {drug} | {status} | {h.get('note','')} |"
            )
        lines.append(
            "\n_VAF (variant allele fraction) and the Call column matter specifically because "
            "Pf is being called diploid here to capture clonal heterogeneity within a single "
            "infection, not true biological diploidy. A 'mixed infection' call means some but "
            "not all parasite clones in this sample carry the marker -- the GT column alone "
            "cannot tell you that. Call thresholds (wildtype-dominant / mixed / mutant-fixed) "
            "are NOT validated against your data -- review and adjust in `classify_vaf()` "
            "before trusting the labels; the raw VAF number is always shown alongside so you "
            "can re-judge it yourself._\n"
        )
    else:
        lines.append("_No DR/IR marker hits found in this run._\n")

    with open(outfile, "w") as O:
        O.write("\n".join(lines) + "\n")
    return outfile
