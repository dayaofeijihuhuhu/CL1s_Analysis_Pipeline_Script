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
