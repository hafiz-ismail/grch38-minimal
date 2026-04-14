#!/usr/bin/env bash
# subset_annotations.sh
# ---------------------
# Downloads the Ensembl GRCh38 release 113 GTF and GFF3 annotation files
# and subsets them to 3 marker genes:
#
#   Gene   | Chromosome | Strand | Biotype        | Ensembl ID
#   -------|------------|--------|----------------|------------------
#   ACTB   | 7          | +      | protein_coding | ENSG00000075624
#   GAPDH  | 12         | -      | protein_coding | ENSG00000111640
#   RPPH1  | 14         | .      | misc_RNA       | ENSG00000264726
#
# For the GFF3, two additional biological_region features are retained.
#
# Usage:  bash subset_annotations.sh
# Requirements: bash, wget (or curl), zcat/gzip

set -euo pipefail

# ── configuration ──────────────────────────────────────────────
ENSEMBL_RELEASE=113
BASE_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/homo_sapiens"
GTF_GZ="Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz"
GFF3_GZ="Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz"
GTF_URL="${BASE_URL}/gtf/homo_sapiens/${GTF_GZ}"
GFF3_URL="${BASE_URL}/gff3/homo_sapiens/${GFF3_GZ}"

GENE_IDS="ENSG00000075624|ENSG00000111640|ENSG00000264726"
GENE_NAMES="ACTB|GAPDH|RPPH1"

OUT_GTF="Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.subset.gtf"
OUT_GFF3="Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.subset.gff3"

# ── helper: prefer wget, fall back to curl ──────────────────────
download() {
    local url="$1" dest="$2"
    if command -v wget &>/dev/null; then
        wget -q -O "$dest" "$url"
    else
        curl -fsSL -o "$dest" "$url"
    fi
}

# ── Step 1: download ────────────────────────────────────────────
echo "[1/4] Downloading GTF ..."
download "$GTF_URL" "$GTF_GZ"

echo "[2/4] Downloading GFF3 ..."
download "$GFF3_URL" "$GFF3_GZ"

# ── Step 2: subset GTF ─────────────────────────────────────────
# Strategy:
#   • keep all header lines (starting with #)
#   • keep any data line whose attribute field contains one of the
#     three target gene IDs – this captures gene, transcript, exon,
#     CDS, start_codon, stop_codon and UTR rows in one pass.
echo "[3/4] Subsetting GTF ..."
{
    # Header block
    zcat "$GTF_GZ" | grep '^#'

    # Data lines: match any of the three gene IDs in column 9
    zcat "$GTF_GZ" | grep -v '^#' \
        | grep -E "gene_id \"(${GENE_IDS})\""
} > "$OUT_GTF"

echo "  → written: $OUT_GTF  ($(grep -vc '^#' "$OUT_GTF") data lines)"

# ── Step 3: subset GFF3 ────────────────────────────────────────
# Strategy:
#   • keep all header lines
#   • keep chromosome/region lines for the 3 relevant chromosomes (7, 12, 14)
#   • keep up to 2 biological_region features (non-gene/non-transcript rows)
#   • for each target gene, keep:
#       – the gene feature  (ID=gene:<ID>)
#       – transcript/mRNA   (Parent=gene:<ID>)
#       – exon/CDS/UTR      (Parent=transcript:<transcript_ID>)
#
# GFF3 parent–child relationships require a two-pass approach:
#   pass 1 – collect gene IDs and their child transcript IDs
#   pass 2 – emit lines whose ID or Parent matches collected IDs
echo "[4/4] Subsetting GFF3 ..."
python3 - "$GFF3_GZ" "$OUT_GFF3" "$GENE_IDS" <<'PYEOF'
import sys, gzip, re

gff3_gz, out_path, gene_ids_str = sys.argv[1], sys.argv[2], sys.argv[3]
target_gene_ids = set(gene_ids_str.split("|"))

opener = gzip.open if gff3_gz.endswith(".gz") else open

header_lines   = []
chr_lines      = []   # chromosome / region features
bioreg_lines   = []   # biological_region features (collect all, keep 2)
gene_lines     = {}   # gene_id  -> line
txn_lines      = {}   # transcript_id -> line
child_lines    = []   # exon / CDS / UTR – need parent lookup

gene_to_txns   = {}   # gene_id -> [transcript_id, ...]
txn_ids        = set()

# ── pass 1: collect all relevant lines ──────────────────────────
with opener(gff3_gz, "rt") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            header_lines.append(line)
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        feat_type = cols[2]
        attrs     = cols[8]

        if feat_type in ("chromosome", "region"):
            # keep only the chromosomes we care about
            if cols[0] in ("7", "12", "14"):
                chr_lines.append(line)
            continue

        if feat_type == "biological_region":
            bioreg_lines.append(line)
            continue

        # Extract ID and Parent
        id_m     = re.search(r'(?:^|;)ID=([^;]+)', attrs)
        parent_m = re.search(r'(?:^|;)Parent=([^;]+)', attrs)
        feat_id  = id_m.group(1)     if id_m     else ""
        parent   = parent_m.group(1) if parent_m else ""

        # Gene features
        if feat_type == "gene":
            gid = feat_id.replace("gene:", "")
            if gid in target_gene_ids:
                gene_lines[gid] = line
                gene_to_txns.setdefault(gid, [])
            continue

        # Transcript / mRNA features
        if feat_type in ("mRNA", "transcript", "ncRNA", "lnc_RNA",
                          "pseudogenic_transcript", "misc_RNA"):
            parent_gid = parent.replace("gene:", "")
            if parent_gid in target_gene_ids:
                tid = feat_id.replace("transcript:", "")
                txn_lines[tid] = line
                txn_ids.add(tid)
                gene_to_txns.setdefault(parent_gid, []).append(tid)
            continue

        # Child features (exon, CDS, UTR, etc.)
        parent_tid = parent.replace("transcript:", "")
        child_lines.append((parent_tid, line))

# ── pass 2: write output ─────────────────────────────────────────
with open(out_path, "w") as out:
    for h in header_lines:
        out.write(h + "\n")

    # chromosome lines
    for l in chr_lines:
        out.write(l + "\n")

    # up to 2 biological_region lines
    for l in bioreg_lines[:2]:
        out.write(l + "\n")

    # gene → transcript → children
    for gid in sorted(target_gene_ids):
        if gid not in gene_lines:
            print(f"  WARNING: gene {gid} not found in GFF3", file=sys.stderr)
            continue
        out.write(gene_lines[gid] + "\n")
        for tid in gene_to_txns.get(gid, []):
            out.write(txn_lines[tid] + "\n")
            for parent_tid, child_line in child_lines:
                if parent_tid == tid:
                    out.write(child_line + "\n")

print(f"  → written: {out_path}")
PYEOF

echo ""
echo "Done.  Output files:"
echo "  GTF  : $OUT_GTF"
echo "  GFF3 : $OUT_GFF3"
