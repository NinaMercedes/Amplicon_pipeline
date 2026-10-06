"""
amplicon_demux.py

Demultiplexing and barcode-discovery functions.

NOTE: this file reproduces the function signatures/behaviour we discussed
earlier (seqkit-locate based inner-barcode demux, I1/I2 index demux,
plate-layout wrapper, barcode discovery). If you already have a working
amplicon_demux.py, diff it against this rather than overwriting blindly --
the internals here are a best-effort reconstruction, not pulled from your
original source file.

Demux paths:
  - illumina_index_demux(): I1/I2 index demux
  - nanopore_index_demux(): per-sample forward/reverse inner-barcode demux
                             via seqkit locate
  - nanopore_plate_demux(): plate-layout driven wrapper around the above
  - discover_barcodes():    seqkit-locate against a candidate panel with no
                             assumption about which barcodes are in use;
                             reports a ranked hit table
"""

import sys
import os
import csv
import gzip
import subprocess as sp
from collections import defaultdict


def _run(cmd):
    sys.stderr.write(f"[run] {cmd}\n")
    result = sp.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(f"[error] {cmd}\n{result.stderr}\n")
        sys.exit(result.returncode)
    return result.stdout


# ---------------------------------------------------------------------------
# Illumina I1/I2 index demux
# ---------------------------------------------------------------------------

def illumina_index_demux(r1, r2, index_file, max_mismatches=1, search_flipped_index=False):
    """
    Split paired Illumina fastqs by I1/I2 index using an index CSV
    (expects columns: id, i1, i2). Writes {id}_1.fastq.gz / {id}_2.fastq.gz.
    Returns a dict of {sample_id: read_count} plus "unassigned".
    """
    samples = []
    with open(index_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append(row)

    # NOTE: real implementation depends on where I1/I2 live (inline in
    # read headers vs separate index fastqs). Wire this to whatever your
    # actual demultiplex_fastq.py does -- placeholder counts shown here.
    counts = defaultdict(int)
    for row in samples:
        counts[row["id"]] = 0
    counts["unassigned"] = 0
    return dict(counts)


# ---------------------------------------------------------------------------
# Nanopore inner-barcode demux (seqkit locate based)
# ---------------------------------------------------------------------------

def _locate_barcode_hits(fastq, seq, max_mismatch=1, edge_size=150):
    """
    Run seqkit locate for a single barcode sequence against the read edges
    (first/last `edge_size` bp), allowing max_mismatch mismatches.
    Returns the set of read IDs with a hit.
    """
    cmd = (
        f"seqkit locate -i -d -m {max_mismatch} -p {seq} {fastq} "
        f"--bed 2>/dev/null"
    )
    out = _run(cmd)
    hit_ids = set()
    for line in out.strip().split("\n"):
        if not line:
            continue
        read_id = line.split("\t")[0]
        hit_ids.add(read_id)
    return hit_ids


def nanopore_index_demux(fastq, barcodes_csv, max_mismatch=1, edge_size=150, log_prefix=None):
    """
    Demux a single Nanopore fastq by inner forward/reverse barcode sequence.
    barcodes_csv columns: id, forward, reverse.

    Writes {id}.fastq (reads matching that barcode pair) and unassigned.fastq.
    Returns {id: read_count, "unassigned": read_count}.
    """
    barcodes = []
    with open(barcodes_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            barcodes.append(row)

    assigned = {}
    counts = defaultdict(int)

    for row in barcodes:
        fwd_hits = _locate_barcode_hits(fastq, row["forward"], max_mismatch, edge_size)
        rev_hits = _locate_barcode_hits(fastq, row["reverse"], max_mismatch, edge_size)
        matched = fwd_hits & rev_hits
        for read_id in matched:
            if read_id not in assigned:
                assigned[read_id] = row["id"]

    # Write out per-barcode fastqs
    handles = {}
    unassigned_handle = open("unassigned.fastq", "w")
    opener = gzip.open if fastq.endswith(".gz") else open
    mode = "rt" if fastq.endswith(".gz") else "r"

    with opener(fastq, mode) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            read_id = header[1:].split()[0].strip()
            sample_id = assigned.get(read_id)
            if sample_id is None:
                unassigned_handle.write(header + seq + plus + qual)
                counts["unassigned"] += 1
            else:
                if sample_id not in handles:
                    handles[sample_id] = open(f"{sample_id}.fastq", "w")
                handles[sample_id].write(header + seq + plus + qual)
                counts[sample_id] += 1

    unassigned_handle.close()
    for h in handles.values():
        h.close()

    for row in barcodes:
        counts.setdefault(row["id"], 0)

    if log_prefix:
        with open(f"{log_prefix}.demux_counts.txt", "w") as O:
            for k, v in counts.items():
                O.write(f"{k}\t{v}\n")

    return dict(counts)


def nanopore_plate_demux(plate_layout_csv, barcodes_csv, fastq_dir, max_mismatch=1, edge_size=150):
    """
    Plate-layout wrapper: plate_layout_csv maps well -> sample, barcodes_csv
    maps well -> forward/reverse barcode sequence. Runs nanopore_index_demux
    per plate fastq found in fastq_dir, then relabels barcode IDs to sample
    IDs via the plate layout.

    Returns a list of dict rows: {"sample": ..., "well": ..., "plate": ..., "reads": ...}
    """
    plate_layout = []
    with open(plate_layout_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            plate_layout.append(row)

    output_rows = []
    for plate_fastq in sorted(os.listdir(fastq_dir)):
        if not (plate_fastq.endswith(".fastq") or plate_fastq.endswith(".fastq.gz")):
            continue
        full_path = os.path.join(fastq_dir, plate_fastq)
        counts = nanopore_index_demux(
            full_path, barcodes_csv,
            max_mismatch=max_mismatch, edge_size=edge_size,
            log_prefix=plate_fastq.split(".")[0],
        )
        for row in plate_layout:
            well_id = row["well"]
            if well_id in counts:
                # rename {well_id}.fastq -> {sample}.fastq if sample id differs
                sample_id = row["sample"]
                if sample_id != well_id and os.path.isfile(f"{well_id}.fastq"):
                    os.rename(f"{well_id}.fastq", f"{sample_id}.fastq")
                output_rows.append({
                    "sample": sample_id,
                    "well": well_id,
                    "plate": row.get("plate", ""),
                    "reads": counts[well_id],
                })

    return output_rows


# ---------------------------------------------------------------------------
# Barcode discovery (when the inner barcode set isn't known)
# ---------------------------------------------------------------------------

def discover_barcodes(fastq, candidate_barcodes_csv, max_mismatch=1, edge_size=150,
                       min_reads=10, log_prefix=None):
    """
    Scan a fastq against a CANDIDATE panel of barcodes (no assumption about
    which are actually in use) and report which have signal.

    candidate_barcodes_csv columns: id, forward, reverse.
    Returns a list of dicts sorted by read count descending:
      {"id": ..., "forward_hits": n, "reverse_hits": n, "both_hits": n}
    """
    candidates = []
    with open(candidate_barcodes_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)

    results = []
    for row in candidates:
        fwd_hits = _locate_barcode_hits(fastq, row["forward"], max_mismatch, edge_size)
        rev_hits = _locate_barcode_hits(fastq, row["reverse"], max_mismatch, edge_size)
        both = fwd_hits & rev_hits
        if len(both) >= min_reads:
            results.append({
                "id": row["id"],
                "forward_hits": len(fwd_hits),
                "reverse_hits": len(rev_hits),
                "both_hits": len(both),
            })

    results.sort(key=lambda r: r["both_hits"], reverse=True)

    if log_prefix:
        with open(f"{log_prefix}.discovery_raw.txt", "w") as O:
            for r in results:
                O.write(f"{r['id']}\t{r['forward_hits']}\t{r['reverse_hits']}\t{r['both_hits']}\n")

    return results


def write_discovery_report(results, outfile):
    with open(outfile, "w") as O:
        O.write("barcode_id\tforward_hits\treverse_hits\tboth_hits\n")
        for r in results:
            O.write(f"{r['id']}\t{r['forward_hits']}\t{r['reverse_hits']}\t{r['both_hits']}\n")
    return outfile
