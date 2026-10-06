# Amplicon Pipeline — Final Usage Guide

Species-agnostic amplicon sequencing pipeline (Nanopore + Illumina), with
per-amplicon splitting for pooled samples, selectable ONT variant callers,
and gene+amino-acid-change DR/IR marker detection that works for SNPs,
indels, and multi-codon "weird combination" changes alike.

## 0. Setup

```bash
conda env create -f environment.yml
conda activate amplicon-pipeline
```

Flagged, unverified: `fastplong` package availability in conda hasn't been
checked (`mamba search fastplong` first). `medaka` and `clair3` both have
heavy/specific dependency pins -- if `conda env create` fails or resolves
very slowly, consider separate environments for each, called via subprocess.

## 1. Files you need

| File | Used by | Format |
|---|---|---|
| Reference FASTA (`ref.fasta`) | all modes | indexed automatically |
| GFF3 (`ref.gff3`) | all modes | same coordinate system as `ref.fasta` |
| Amplicon BED (`amplicons.bed`) | all modes | **column 4 (name) required & unique** -- used as the amplicon ID for splitting pooled samples |
| Manifest CSV | `nanopore-direct` | `sample,fastq` |
| Inner barcodes CSV | `nanopore-index` | `id,forward,reverse` |
| Plate layout CSV | `nanopore-plate` | `well,sample,plate` |
| Plate barcodes CSV | `nanopore-plate` | `well/id,forward,reverse` |
| Illumina index CSV | `illumina-index` | `id,i1,i2` |
| Candidate barcode panel CSV | `discover-barcodes` | `id,forward,reverse` |
| DR marker table | optional, all modes | `gene,aa_change,linked_group,...` -- see below |
| IR marker table | optional, all modes | same schema as DR table |

## 2. Choosing a Nanopore mode

| Your situation | Mode |
|---|---|
| Dorado barcode == sample, nothing pooled inside it | `nanopore-direct` |
| One barcode contains multiple pooled samples, inner barcodes known | `nanopore-index` |
| Multiple plates, same reusable barcode set across plates | `nanopore-plate` |
| Pooled, inner barcodes unknown | `discover-barcodes` first, then `nanopore-index`/`nanopore-plate` |

The deciding factor for `nanopore-index` vs `nanopore-plate` isn't sample
count -- it's whether barcode identity *alone* determines sample identity
(`nanopore-index`), or only does so combined with which plate/fastq it came
from, because the same barcode is being reused across multiple plates
(`nanopore-plate`).

## 3. Example input files

### `amplicons.bed`
```
Pf3D7_05_v3   429146  430547  kelch13
Pf3D7_07_v3   403177  403832  crt
Pf3D7_14_v3   1953500 1954500 mdr1
```

### `manifest.csv` (nanopore-direct)
```csv
sample,fastq
patient001,barcode01.fastq.gz
patient002,barcode02.fastq.gz
```

### `barcodes.csv` (nanopore-index)
```csv
id,forward,reverse
sample01,AAGAAAGTTGTCGGTGTCTTTGTG,CTGCCAAGGCACACAGGGTC
sample02,GAGAGGACAAAGGTTTCAACGC,CTGCCAAGGCACACAGGGTC
```

### `plate_layout.csv` (nanopore-plate)
```csv
well,sample,plate
A01,patient001,plate1
A02,patient002,plate1
A01,patient097,plate2
```

### `illumina_index.csv` (illumina-index)
```csv
id,i1,i2
sample01,AACGTGAT,GCCTAAGT
```

### DR/IR marker table (`dr_markers.csv` / `ir_markers.csv`)

**No genomic coordinates required.** Matching happens on gene + amino-acid
change, extracted live from the `bcftools csq` annotation already present
in your VCFs -- this is what lets the same logic handle SNPs, indels, and
multi-codon combined changes (like the PfCRT 72-76 haplotype) without any
special-casing.

```csv
gene,aa_change,linked_group,drug,who_status,reference,note
PF3D7_0709000,K76T,,chloroquine,drug_resistance,https://www.who.int/...,
PF3D7_0709000,M74I,crt_72_76_indel,chloroquine,drug_resistance,...,
PF3D7_0709000,N75E,crt_72_76_indel,chloroquine,drug_resistance,...,
```

- `gene`: must match the gene name bcftools csq uses in your annotation.
- `aa_change`: plain `{wt}{codon}{mut}` string, e.g. `K76T`, `L1014F`.
- `linked_group`: markers sharing a non-empty group name are only reported
  together, when every member is present in the same sample. This is how
  the PfCRT indel pair (`M74I`+`N75E`) and the kdr super-kdr combo
  (`L1014F`+`N1575Y`) are both handled.
- Other columns (drug, insecticide_class, who_status, reference, note,
  marker_type) carry straight through into the report.

`dr_markers.csv` and `ir_markers.csv` (provided) are starting points. Two
items in the Pf table still need your decision:
- `p.Asn75Glu` only exists as part of the linked crt indel, never standalone -- already handled via `linked_group`.
- `A437G` direction conflict (the annotation source's reference allele looked reversed vs literature convention) -- confirm via `samtools faidx ref.fasta Pf3D7_08_v3:549684-549686` before trusting that row.

The Anopheles IR table has confirmed genomic coordinates (AgamP4) for most
VGSC, Rdl, Ace1, and GSTe2 markers, sourced from a user-supplied genomic
table that **resolved two earlier conflicts**: Rdl wild-type is confirmed
Ala (A296G, not the alternative "Cys-wildtype" claim from another source),
and Ace1 G280S (native numbering) is confirmed as the same mutation as the
commonly-cited G119S (Torpedo/Musca numbering) -- don't add both as
separate rows. One conflict was NOT resolved and is kept as two distinct
rows rather than merged: GSTe2 codon 119 has a genomic-table-confirmed
L119V allele and a separately-sourced, functionally-validated-but-
coordinate-unconfirmed L119F allele -- these cannot both come from the
same single-nucleotide substitution, so they're real, different markers.
Cyp4j5, Coeae1d, and the Cyp6p4 structural marker remain unconfirmed/out of
scope (see inline notes in the CSV).

`build_dr_ir_table.py` is a reusable utility for building marker tables from
a MalariaProfiler-style export plus a per-position annotation file -- it's
not consumed by the pipeline itself, since matching moved to gene+aa_change.

## 4. Running each mode

```bash
# Dorado barcode == sample (DR Pf data)
python amplicon_pipeline.py nanopore-direct \
  --fastq-dir ./dorado_demux_output/ --manifest manifest.csv \
  --ref ref.fasta --gff ref.gff3 --bed amplicons.bed \
  --ploidy diploid --dr-ir-snps dr_markers.csv

# Pooled, single fastq, known inner barcodes
python amplicon_pipeline.py nanopore-index \
  --fastq run1.fastq.gz --barcodes barcodes.csv \
  --ref ref.fasta --gff ref.gff3 --bed amplicons.bed \
  --ploidy diploid --dr-ir-snps dr_markers.csv

# Pooled plates (Anisa anopheles)
python amplicon_pipeline.py nanopore-plate \
  --plate-layout plate_layout.csv --barcodes barcodes.csv \
  --fastq-dir ./plate_fastqs/ \
  --ref ref.fasta --gff ref.gff3 --bed amplicons.bed \
  --ploidy diploid --ir-snps ir_markers.csv

# Illumina, indexed (COI Pf data)
python amplicon_pipeline.py illumina-index \
  --read1 run_R1.fastq.gz --read2 run_R2.fastq.gz \
  --index-file illumina_index.csv \
  --ref ref.fasta --gff ref.gff3 --bed amplicons.bed \
  --ploidy diploid --dr-ir-snps dr_markers.csv

# Discover unknown inner barcodes
python amplicon_pipeline.py discover-barcodes \
  --fastq run1.fastq.gz --candidate-barcodes candidate_barcodes.csv \
  --min-reads 20 --outfile barcode_discovery_report.txt
```

`--ploidy diploid` is used for Pf deliberately, to capture clonal
heterogeneity within a single infection (mixed parasite clones), not
because Pf is biologically diploid. See Section 6 for how to read the
resulting heterozygous calls.

## 5. Optional flags

| Flag | Default | Effect |
|---|---|---|
| `--variant-caller` | `medaka` | ONT only: `medaka` or `clair3` |
| `--clair3-model` | `None` | Required if `--variant-caller clair3`. Must match your basecaller model/chemistry exactly. |
| `--downsample-depth N` | `None` (off) | Downsample each (sample, amplicon) BAM to ~Nx mean depth before variant calling |
| `--dr-ir-snps` | `None` | DR marker table path |
| `--ir-snps` | `None` | IR marker table path. Can be combined with `--dr-ir-snps`. |
| `--vaf-mixed-low` | `0.05` | Below this VAF, classified "wildtype-dominant (check vs noise)". Not validated -- review before trusting. |
| `--vaf-mixed-high` | `0.95` | Above this VAF, classified "mutant-fixed". Same caveat. |
| `--min-depth` / `--min-qual` | `10` / `30` | Genotype/variant QC filters |
| `--threads` | `4` | |

## 6. Reading the DR/IR section of `report.md`

Because Pf is called diploid for clonal-heterogeneity reasons, a `0/1`
genotype means *some but not all parasite clones in this infection carry the
marker* -- not noise, not true heterozygosity. The report gives you the
numbers to interpret that, not just a yes/no:

| Column | Meaning |
|---|---|
| GT | Raw genotype call |
| AD (ref,alt) | Read counts backing the call |
| DP | Total depth at this position |
| VAF | Alt allele fraction -- what fraction of the infection carries it |
| Call | wildtype-dominant / mixed infection / mutant-fixed, from VAF vs thresholds |

For linked-group hits, VAF is reported as the minimum across components,
and the report flags a warning if components' VAFs disagree by more than
0.10 -- a truly linked indel pair should show near-identical VAF since it's
the same underlying reads; disagreement is worth a manual look.

VAF thresholds are not validated against your data -- treat 5%/95% as a
starting point. Known mono-clonal controls should cluster near 0% or 100%;
if they don't, revisit the thresholds or sequencing depth.

## 7. Outputs

- `{sample}_{amplicon}.bam` / `.vcf.gz` -- per-unit alignments and calls
- `combined.vcf.gz` / `combined.filtered.norm.vcf.gz` -- merged multi-sample VCF
- `combined.snps.txt`, `combined.indels.txt` -- flat tables (counts only; marker matching reads annotated VCFs directly)
- `report.md` -- units, variant counts, DR/IR marker hits with VAF/Call

## 8. Known gaps before trusting real results

- `amplicon_demux.py` internals reconstructed from discussion -- diff against your actual implementation if one exists.
- `nanopore_plate_demux` plate-to-fastq filename matching untested against real layouts.
- `illumina_index_demux` is still a stub (no real I1/I2 logic yet).
- Clair3 haploid flags unvalidated against truth data.
- `extract_csq_annotations`'s assumed bcftools csq field format should be confirmed against one real annotated VCF from your pipeline before trusting marker matching end-to-end.
- `A437G` direction conflict in `dr_markers.csv` (Section 3) -- still open.
- `ir_markers.csv`: Rdl wild-type and Ace1 native/Torpedo numbering
  conflicts are now RESOLVED via confirmed AgamP4 coordinates (Section 3).
  GSTe2 codon 119 (L119V vs L119F) remains an open, deliberately-unmerged
  discrepancy -- kept as two separate marker rows.
- `Cyp4j5`, `Coeae1d`, `Cyp6p4` rows in `ir_markers.csv` lack confirmed
  coordinates; `Coeae1d`'s aa_change encoding is still ambiguous;
  `Cyp6p4-236M-TE-Cyp6aap-Dup1` is a structural variant out of scope for
  point-mutation matching.
- VAF thresholds unvalidated (Section 6).
- Numbering-scheme consistency between `ir_markers.csv` and your actual
  bcftools csq annotation output is still unverified for VGSC/Rdl/Ace1/GSTe2
  -- a mismatch would cause silent false negatives (Section 3). The L995S
  and A296S alt-allele nucleotides were inferred, not confirmed -- verify
  against your own VCF.
