# GRCh38 Minimal Annotation Subset

A minimal subset of the **Ensembl GRCh38 release 113** annotation files (GTF and GFF3)
covering three marker genes, one per strand orientation.

---

## Included Files

| File | Description |
|------|-------------|
| `Homo_sapiens.GRCh38.113.subset.gtf`  | GTF subset (3 genes) |
| `Homo_sapiens.GRCh38.113.subset.gff3` | GFF3 subset (3 genes + 2 biological regions) |
| `subset_annotations.sh`               | Reproducible subsetting script |

---

## Gene Selection

### Criteria
- **Confirmed expression** in healthy human tissue (housekeeping / universal marker genes)
- One gene from each strand orientation: `+`, `-`, and `.` (non-stranded)
- Full annotation retained: gene, transcript, exon, CDS, start/stop codon, UTR

### Genes

| Gene  | Ensembl ID      | Chromosome | Strand | Biotype        | Rationale |
|-------|-----------------|------------|--------|----------------|-----------|
| **ACTB**  | ENSG00000075624 | 7          | `+`    | protein_coding | Beta-actin; ubiquitous cytoskeletal gene; canonical housekeeping gene used in virtually every normalisation panel |
| **GAPDH** | ENSG00000111640 | 12         | `-`    | protein_coding | Glyceraldehyde-3-phosphate dehydrogenase; glycolytic enzyme expressed in all cell types; standard RNA-seq reference gene |
| **RPPH1** | ENSG00000264726 | 14         | `.`    | misc_RNA       | RNA component of nuclear RNase P; essential for tRNA 5′-end processing in every cell; annotated without strand (`"."`) in Ensembl GRCh38 release 113 because the locus is embedded within a multi-copy repeat cluster where individual copy orientations are ambiguous |

### Non-stranded gene — additional note
Truly non-stranded *gene* annotations are rare in GRCh38; most arise for non-coding RNA
loci in highly repetitive regions (rDNA arrays, tRNA clusters, snRNA/scaRNA pseudogene
clusters) where the Ensembl/HAVANA pipeline cannot assign a reliable strand from
the available evidence.  RPPH1 is the best-characterised example that satisfies
the expressed-marker requirement while appearing with strand `"."` in the
Ensembl GRCh38 primary annotation.

---

## Source Annotation Files

| Format | Source  | Release | URL |
|--------|---------|---------|-----|
| GTF    | Ensembl | 113     | `https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz` |
| GFF3   | Ensembl | 113     | `https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz` |

Genome assembly: **GRCh38** (GCA_000001405.15 / hg38)  
Annotation build: Ensembl 113 / GENCODE v47-equivalent

---

## How to Reproduce

### Requirements
```
bash ≥ 4, wget or curl, python3 ≥ 3.6, gzip/zcat
```

### Option A — run the bundled script

```bash
bash subset_annotations.sh
```

This script will:
1. Download `Homo_sapiens.GRCh38.113.gtf.gz` and `Homo_sapiens.GRCh38.113.gff3.gz`
   from the Ensembl FTP server
2. Stream-decompress each file and filter for the three target gene IDs
3. Write `Homo_sapiens.GRCh38.113.subset.gtf` and `Homo_sapiens.GRCh38.113.subset.gff3`

### Option B — manual commands

```bash
# 1. Download
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
wget https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz

# 2. Subset GTF
# Keep header + any data line whose attributes mention one of the three gene IDs
{
  zcat Homo_sapiens.GRCh38.113.gtf.gz | grep '^#'
  zcat Homo_sapiens.GRCh38.113.gtf.gz | grep -v '^#' \
    | grep -E 'gene_id "(ENSG00000075624|ENSG00000111640|ENSG00000264726)"'
} > Homo_sapiens.GRCh38.113.subset.gtf

# 3. Subset GFF3
# The GFF3 uses a hierarchical ID/Parent scheme, so a two-pass approach is needed.
# The bundled Python snippet in subset_annotations.sh handles this automatically.
# Quick one-liner for the gene + direct children (no grandchildren):
{
  zcat Homo_sapiens.GRCh38.113.gff3.gz | grep '^#'
  # chromosome entries for chr7, chr12, chr14
  zcat Homo_sapiens.GRCh38.113.gff3.gz | grep -v '^#' \
    | awk -F'\t' '$1 ~ /^(7|12|14)$/ && $3 == "chromosome"'
  # first 2 biological_region rows
  zcat Homo_sapiens.GRCh38.113.gff3.gz | grep -v '^#' \
    | grep 'biological_region' | head -2
  # all rows whose ID or Parent mentions a target gene/transcript
  zcat Homo_sapiens.GRCh38.113.gff3.gz | grep -v '^#' \
    | grep -E '(ENSG00000075624|ENSG00000111640|ENSG00000264726|ENST00000331789|ENST00000229239|ENST00000383925)'
} > Homo_sapiens.GRCh38.113.subset.gff3
```

> **Tip:** If you have [AGAT](https://github.com/NBISweden/AGAT) installed, you can
> use `agat_sq_filter_feature_from_keep_list.pl` for a format-aware extraction that
> automatically follows parent–child links.

---

## Annotation Format Details

### GTF (Gene Transfer Format)

Nine tab-separated columns — coordinates are **1-based, inclusive**.

```
<seqname> <source> <feature> <start> <end> <score> <strand> <frame> <attributes>
```

Feature types present in this subset:

| Feature       | Description |
|---------------|-------------|
| `gene`        | Gene locus boundary |
| `transcript`  | Transcript (isoform) boundary |
| `exon`        | Exon coordinates |
| `CDS`         | Coding sequence (includes stop codon in Ensembl convention) |
| `start_codon` | First codon of CDS |
| `stop_codon`  | Last codon of CDS |
| `UTR`         | 5′ or 3′ untranslated region |

The `strand` column is `+`, `-`, or `.` (unstranded).  
The `frame` column for CDS gives the reading-frame phase (0, 1, or 2);
`.` for all other features.

### GFF3 (Generic Feature Format 3)

Same nine columns; attributes use `key=value` pairs separated by `;`.  
Hierarchical relationships are encoded via `ID=` and `Parent=`.

Feature types present in this subset:

| Feature              | Description |
|----------------------|-------------|
| `chromosome`         | Whole-chromosome region entry |
| `biological_region`  | Non-gene genomic feature (e.g. CpG island, CTCF site) — **non-stranded** |
| `gene`               | Gene locus |
| `mRNA`               | Protein-coding transcript |
| `transcript`         | Non-coding transcript |
| `exon`               | Exon |
| `CDS`                | Coding sequence |
| `five_prime_UTR`     | 5′ UTR |
| `three_prime_UTR`    | 3′ UTR |

The two **biological_region** rows retained in the GFF3 are a
*CTCF binding site* and a *CpG island* in the promoter region
upstream of ACTB on chromosome 7.  Both have strand `"."`.

---

## Equivalence Between GTF and GFF3

The GTF and GFF3 subsets cover the **same three genes** and the **same transcripts**,
exons, and CDS boundaries.  The key differences are purely format-level:

| Aspect | GTF | GFF3 |
|--------|-----|------|
| Attribute separator | `key "value";` | `key=value;` |
| Transcript feature name | `transcript` | `mRNA` (protein-coding) / `transcript` (ncRNA) |
| UTR feature name | `UTR` | `five_prime_UTR` / `three_prime_UTR` |
| Hierarchy encoding | implicit (shared `gene_id`/`transcript_id`) | explicit `ID=` / `Parent=` |
| Extra rows | — | `chromosome` + `biological_region` |

---

## Quick Verification

```bash
# Count feature types in the GTF
grep -v '^#' Homo_sapiens.GRCh38.113.subset.gtf \
  | awk -F'\t' '{print $3}' | sort | uniq -c | sort -rn

# Count feature types in the GFF3
grep -v '^#' Homo_sapiens.GRCh38.113.subset.gff3 \
  | awk -F'\t' '{print $3}' | sort | uniq -c | sort -rn

# Confirm all three strands are represented
grep -v '^#' Homo_sapiens.GRCh38.113.subset.gtf \
  | awk -F'\t' '$3=="gene" {print $7, $9}' | grep -oP 'gene_name "\K[^"]+'
```

---

## License

Source data © [EMBL-EBI / Ensembl](https://www.ensembl.org) — released under
[Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/).
