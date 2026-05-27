#!/bin/bash
# =============================================================================
# CL1s (Class 1 Integron) Analysis Pipeline for " Employing Class 1 Integron as
# a Universal Indicator of Antimicrobial Resistance Risk in Livestock Farming "
# =============================================================================
# This script performs comprehensive analysis of IntI1-positive contigs,
# including antimicrobial resistance (ARG), metal resistance (MRG),
# mobile genetic elements (MGE), integron structures,
# plasmid prediction, taxonomic classification, and phylogenetic analysis.
#
# Software versions:
#   abricate v1.0.1 | integron_finder v2.0.5 | diamond v2.1.9.163
#   MMseqs2 v14.7e284 | cd-hit v4.8.1 | muscle v3.8.31 | trimal v1.4
#   IQ-TREE v2.2.6 | PlasClass v0.1.1 | PLASMe v1.1 | PlasmidHunter v1.4
#   MOB-Typer v3.1.9 | GTDB-Tk v2.4.0 (R220) | blast+ v2.15.0+
#   ARGs-OAP v3.2.2 | Trimmomatic v2.3
#
# Databases:
#   ResFinder v2.4.0 | BacMet (31 Jan 2025) | MGEs DB (31 Jan 2025)
#   TnCentral (31 Jan 2025) | IMG/VR (1 Jan 2025) | GTDB R220
# =============================================================================

set -euo pipefail

# =============================================================================
# 0. Configuration — EDIT THESE PATHS
# =============================================================================

# Input data
CONTIGS="path/to/all_intI1_contigs.fa"       # IntI1-positive contigs
FLANK="path/to/all_intI1_flank.fa"           # Flanking sequences
CONTIGS_1000="path/to/intI1_1000_contigs.fa" # Contigs >= 1000 bp

# Output directory
OUTDIR="path/to/output"
mkdir -p "$OUTDIR"

# Database paths
BACMET_DB="path/to/BacMet2_EXP_database.fasta"
TNCENTRAL_DB="path/to/tncentral_fa/tn"
IMG_DB="path/to/IMG_VR/img"

# Threads
THREADS=8

# =============================================================================
# 1. Antimicrobial Resistance Genes (ARGs) — abricate v1.0.1 / ResFinder v2.4.0
# =============================================================================
echo "=== 1. ARG detection ==="

abricate --db resfinder "$CONTIGS" --threads $THREADS \
    --minid 80 --mincov 80 > "$OUTDIR/resfinder.out"

# =============================================================================
# 2. Integron Identification — integron_finder v2.0.5
# =============================================================================
echo "=== 2. Integron identification ==="

integron_finder --func-annot "$CONTIGS" --cpu $THREADS --local-max

# =============================================================================
# 3. Metal Resistance Genes (MRGs) — diamond v2.1.9.163 / BacMet (31 Jan 2025)
# =============================================================================
echo "=== 3. Metal resistance gene detection ==="

diamond blastx \
    --db "$BACMET_DB" \
    -q "$CONTIGS" \
    -o "$OUTDIR/bac.out" \
    --sensitive \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen \
               qstart qend sstart send evalue bitscore qlen slen \
    -e 1e-3 -p $THREADS

# =============================================================================
# 4. Mobile Genetic Elements (MGEs) — abricate v1.0.1 / MGEs DB (31 Jan 2025)
# =============================================================================
echo "=== 4. MGE detection (abricate) ==="

abricate --db MGEs "$CONTIGS" --threads $THREADS \
    --minid 80 --mincov 80 > "$OUTDIR/mges.out"

# =============================================================================
# 5. Transposon Detection — blast+ v2.15.0+ / TnCentral (31 Jan 2025)
# =============================================================================
echo "=== 5. Transposon detection (TnCentral) ==="

blastn -perc_identity 60 -culling_limit 1 -evalue 1e-5 \
    -query "$CONTIGS" -db "$TNCENTRAL_DB" \
    -num_threads $THREADS \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen \
              qstart qend sstart send evalue bitscore qlen slen' \
    -out "$OUTDIR/tn.blast.out"

# =============================================================================
# 6. Taxonomic Classification — MMseqs2 v14.7e284 / GTDB R220
# =============================================================================
echo "=== 6. Taxonomy (MMseqs2) ==="

mmseqs easy-taxonomy "$CONTIGS" \
    "/path/to/GTDB_database" \
    "$OUTDIR/mmseq" "$OUTDIR/mmseq_tmp" \
    --split-memory-limit 40G --threads $THREADS --tax-lineage 1

# =============================================================================
# 7. Plasmid Prediction — Multiple Tools
# =============================================================================
echo "=== 7. Plasmid prediction ==="
mkdir -p "$OUTDIR/plasmid"

# ---------------------------------------------------------------------------
# 7a. PlasClass v0.1.1 (score >= 0.5)
# ---------------------------------------------------------------------------
echo "  -> PlasClass"
python /path/to/PlasClass/classify_fasta.py \
    -f "$CONTIGS_1000" -o "$OUTDIR/plasmid/plasclass.out" -p $THREADS

# ---------------------------------------------------------------------------
# 7b. PLASMe v1.1
# ---------------------------------------------------------------------------
echo "  -> PLASMe"
python /path/to/PLASMe/PLASMe.py \
    "$CONTIGS_1000" "$OUTDIR/plasmid/plasme.fa" -t $THREADS \
    -d /path/to/plasme_database --temp "$OUTDIR/plasmid/temp"

# ---------------------------------------------------------------------------
# 7c. PlasmidHunter v1.4 (score = 1.0)
# ---------------------------------------------------------------------------
echo "  -> PlasmidHunter"
plasmidhunter -i "$CONTIGS_1000" -c $THREADS -o "$OUTDIR/plasmid/plasmidhunter"

# ---------------------------------------------------------------------------
# 7d. IMG/VR plasmid BLAST — blast+ v2.15.0+ / IMG/VR (1 Jan 2025)
# ---------------------------------------------------------------------------
echo "  -> IMG/VR BLAST"
mkdir -p "$OUTDIR/plasmid/img"
blastn -perc_identity 80 -culling_limit 1 -evalue 1e-5 \
    -query "$CONTIGS_1000" -db "$IMG_DB" \
    -num_threads $THREADS \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen \
              qstart qend sstart send evalue bitscore qlen slen' \
    -out "$OUTDIR/plasmid/img/img.out"

# ---------------------------------------------------------------------------
# 7e. MOB-typer v3.1.9 for plasmid mobility
# ---------------------------------------------------------------------------
echo "  -> MOB-typer"
mob_typer --multi -n $THREADS \
    --infile "$OUTDIR/plasmid/plasmid_contigs.fa" \
    --out_file "$OUTDIR/plasmid/mobtyper_results.txt"

# =============================================================================
# 8. Phylogenetic Analysis — cd-hit v4.8.1 / muscle v3.8.31 / trimal v1.4 /
#                           IQ-TREE v2.2.6 / ggtree v3.0.4 (visualization)
# =============================================================================
echo "=== 8. Phylogenetic analysis ==="

# 8a. cd-hit clustering at 98% identity
# cd-hit-est -i "$CONTIGS" -o cluster.fa -aS 0.98 -c 0.98 -G 0 -g 0 -T $THREADS -M 40000

# 8b. Multiple sequence alignment
muscle -in "path/to/intI1_sequences.fa" -out "$OUTDIR/phylogeny/intI1.afa"

# 8c. Trim alignment (gap threshold 0.1)
/path/to/trimal -in "$OUTDIR/phylogeny/intI1.afa" \
    -out "$OUTDIR/phylogeny/intI1.trim.afa" -gt 0.1

# 8d. Maximum-likelihood tree (1,000 ultrafast bootstrap replicates)
iqtree -s "$OUTDIR/phylogeny/intI1.trim.afa" \
    -m MFP -alrt 1000 -B 1000 -T $THREADS \
    --prefix "$OUTDIR/phylogeny/intI1_tree"

# =============================================================================
# 9. Genome Classification — GTDB-Tk v2.4.0 (R220)
# =============================================================================
echo "=== 9. GTDB-tk classification ==="

gtdbtk classify_wf \
    --genome_dir "path/to/genomes" \
    --out_dir "$OUTDIR/gtdb" --extension fa --cpus $THREADS \
    --scratch_dir "$OUTDIR/gtdb/SCRATCH" --skip_ani_screen

# =============================================================================
# Summary
# =============================================================================
echo "=== Analysis complete ==="
echo "Output directory: $OUTDIR"
echo ""
echo "Results:"
echo "  - ARGs:              $OUTDIR/resfinder.out"
echo "  - Integrons:         $OUTDIR/*.integrons"
echo "  - Metal resistance:  $OUTDIR/bac.out"
echo "  - MGEs:              $OUTDIR/mges.out"
echo "  - Transposons:       $OUTDIR/tn.blast.out"
echo "  - Taxonomy:          $OUTDIR/taxonomy.txt"
echo "  - Phylogeny:         $OUTDIR/phylogeny/intI1_tree"
echo "  - GTDB-tk:           $OUTDIR/gtdb/"
