#!/usr/bin/env python3

"""
Synk Synteny Plotting Pipeline

Synk runs compleasm for one main/reference genome and any number of comparison
genomes, then uses the shared complete BUSCO hits to generate one pairwise
synteny result for each main-vs-comparison species pair. It can also keep the
original pairwise mode when compleasm/BUSCO tables, karyotypes, and replacement
maps have already been prepared.

WRAPPER MODE EXAMPLE:
---------------------
python synk.py \
  --main_name pfas \
  --main_assembly pfas.fa \
  --compare tiliqua=tiliqua.fa \
  --compare egernia=egernia.fa \
  --lineage sauropsida \
  --threads 16 \
  --outdir synk_output \
  --plot

WRAPPER INPUTS:
---------------
--main_name              Name for the main/reference species
--main_assembly          Main/reference species genome FASTA
--compare NAME=FASTA     Comparison species genome FASTA; repeat for multiple species
--lineage                BUSCO lineage to pass to compleasm
--autolineage            Let compleasm choose the lineage instead of using --lineage
--threads                Threads passed to compleasm (default: 8)
--compleasm              compleasm command or path (default: compleasm)
--compleasm_library      Optional compleasm lineage library path passed with -L
--reuse_compleasm        Reuse existing compleasm full_table outputs when found
--min_contig_length      Minimum FASTA sequence length to include in auto-karyotypes
--outdir                 Directory for all compleasm, karyotype, and pairwise outputs

OPTIONAL CURATED KARYOTYPES:
----------------------------
--main_karyotype                 Main karyotype file
--main_replacement               Main chromosome replacement map
--compare_karyotype NAME=FILE    Comparison karyotype; repeat as needed
--compare_replacement NAME=FILE  Comparison replacement map; repeat as needed

If curated karyotypes are not supplied, Synk creates karyotype and replacement
files from FASTA sequence lengths. Contigs are sorted longest to shortest and
renumbered 1..N for RIdeogram compatibility.

PAIRWISE MODE EXAMPLE:
----------------------
python synk.py \
  --karyotype1 pfas_karyotype.txt \
  --karyotype2 tiliqua_karyotype.txt \
  --busco1 pfas_full_table.tsv \
  --busco2 tiliqua_full_table.tsv \
  --rep1 replacement_pfas.txt \
  --rep2 replacement_tiliqua.txt \
  --outdir pfas_vs_tiliqua \
  --plot

OUTPUT STRUCTURE:
-----------------
Wrapper mode:
- compleasm/                         : one compleasm run per species
- karyotypes/                        : auto-generated karyotypes and replacement maps
- pairwise/<main>_vs_<comparison>/   : Synk outputs for each species pair

Each pairwise output folder contains:
- chr_color_map.txt                  : Map of main-species chromosomes to colors
- color_replace.txt                  : Map of numeric chromosome order to colors
- merged_busco.txt                   : Shared complete BUSCO hits
- final_synteny.txt                  : Synteny links colored by main-species chromosome
- dual_karyotype.txt                 : Combined karyotype file with label and fill info
- chromosome.png/.svg                : RIdeogram synteny plot when --plot is used

Author: Jon J. Hoffman

Copyright 2025 Jon J. Hoffman

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


# Python dependencies for Synk processing
# pip install pandas matplotlib

# External tools for full wrapper/plotting mode
# compleasm must be available on PATH, or passed with --compleasm
# R, Rscript, RIdeogram, and rsvg are required when --plot is used

# Or with conda for the Python environment
# conda create -n synteny_env python=3.9 pandas matplotlib
# conda activate synteny_env

import argparse
import os
import sys
import csv
import gzip
import shutil
import shlex
from collections import OrderedDict
import subprocess

# === VISUALS === #
ERROR_ART = """
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⣿⣷⠀⠀⠀⠀⣸⣶⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣼⣿⣿⡞⣿⣷⣮⣻⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣿⣿⣿⣾⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⡝⢿⣿⣿⣿⣿⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣤⡀⠀⠀⠀⠀⠀⠀⠻⣿⣿⣿⠸⣸⣻⣏⣿⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⣿⣿⡿⡀⠀⠀⠀⠀⠀⣾⡞⡝⣿⢿⣿⣿⣿⣿⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠩⣾⣿⣶⢦⣤⣀⠸⠻⢭⣥⡻⣧⠀⡙⠛⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⣄⢠⣴⣾⣿⣿⣿⣏⣶⣾⡽⣿⣷⣟⣿⣿⣿⣻⣷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⣀⣀⣀⠀⠀⠀⠸⣿⡿⠘⠻⢿⣿⣿⠟⠛⠿⠿⠃⢍⣿⣿⢸⣿⣿⣿⡽⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⣰⣟⠛⠛⢿⣿⣦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣿⣜⢿⣿⡿⡷⡿⣼⣶⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⢰⣿⠃⠀⠀⠀⠈⢿⣿⣧⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣷⣯⣾⣿⡀⠀⠙⠻⢿⣶⣄⠀⠀⠀⠀⠀⠀⠀
⢸⣿⠀⠀⠀⠀⠀⠀⢻⣿⣷⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣧⡀⠀⠀⠀⠙⢿⣧⡀⠀⠀⠀⠀⠀
⢸⣿⡀⠀⠀⠀⠀⠀⠀⢻⣿⣿⣦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣬⣽⣿⣿⢟⣛⣳⠀⠀⠀⠀⠀⠹⣿⣆⠀⠀⠀⠀
⠀⣿⣇⠀⠀⠀⠀⠀⠀⠈⣿⣿⣿⣷⡀⠀⠀⠀⠀⠀⠀⠀⣴⣿⣿⣿⣿⣷⢻⣾⣿⣿⣷⡽⣄⠀⠀⢀⣾⣿⣷⣄⠀⠀
⠀⠘⣿⣆⠀⠀⠀⠀⠀⠀⠘⣿⣿⣿⣿⣷⣄⡀⠀⠀⢀⣾⣿⣿⣿⣿⣿⣿⡇⣿⣿⣿⣿⣿⢹⣦⠀⢸⣇⠀⠹⣏⢧⡀
⠀⠀⠹⣿⣷⡀⠀⠀⠀⠀⠀⠘⣿⣿⣿⣿⣿⣿⣿⡆⣿⣿⣿⣿⣿⣿⣿⣿⣧⣿⣿⣿⣿⣿⢸⣿⡄⠈⠛⠀⣶⠟⠼⠇
⠀⠀⠀⠹⣿⣿⣷⣤⡀⠀⠀⠀⠘⢿⣿⣿⣿⣿⣿⢸⣿⣿⣿⣿⣿⣿⣿⡿⣼⣿⣿⣿⣿⡿⣾⣿⠁⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠙⣿⣿⣿⣿⣶⣄⠀⠀⠈⠻⣿⣿⣿⣿⢸⣿⣿⣿⣿⣿⣿⡿⣱⣿⣿⣿⣿⢟⣼⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠈⢻⣿⣿⣿⣿⣧⡀⠀⠀⠈⠻⢿⣿⢸⣿⣿⣿⡿⢟⣫⣾⣿⣿⠿⣛⣵⣿⡿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠙⠿⣿⣿⣿⡇⠀⠀⠀⠀⠀⢈⣾⣿⡟⠙⠚⠛⠛⠋⠉⠀⠘⣿⣿⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠛⠁⠀⠀⠀⠀⢀⣾⣿⡟⠀⠀⠀⠀⠀⠀⠀⠀⢰⣿⡿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣾⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣿⡿⡏⠀⠀⠀⠀⠀⠀⠀⠀⢠⣾⣯⢻⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣿⣧⡀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠋⠘⠻⣿⣿⣷⣶⣒⣒⢢⡄⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣿⣿⣿⡿⣏⣃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠻⠿⠿⠟⠈⠁⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⣿⡿⠿⠿⠿⣿⣿⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
"""

SUCCESS_ART = """
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⡀⠀⠄⣀⡀⡀⠤⠐⣢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⠠⢤⠔⠈⠀⠀⠀⠀⠀⠀⠀⠁⠀⣾⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢛⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢿⠏⠀⣀⣀⡀⠀⠀⠀⠀⢀⠀⡔⢻⣦⠀⢃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⣀⡀⠀⠀⠀⠀⣀⡀⠀⠀⠀⠀⠀⠀⠸⠀⠰⠛⣇⣹⡜⡄⠀⠀⠸⢠⣿⣿⣀⠇⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⢀⠔⡉⠤⠐⠒⠒⠒⠂⠠⠬⣁⠒⠠⢄⡀⠀⢠⠀⠐⠠⠿⢿⡇⠀⠀⠀⠀⠈⠛⠉⠀⡀⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⡐⡡⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠐⠂⠄⡁⠒⠱⢤⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠤⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠤⠄⠒⣒⣀⣴
⠰⠰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠐⠢⠌⣉⣶⣶⣦⣄⣀⣠⠔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠠⠂⢁⣠⣴⣾⣿⣿⡟⠁
⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠀⠉⠛⠿⠿⠛⠆⠤⠄⠀⣀⣀⣀⣀⣀⡀⠔⢁⣤⣾⣿⣿⣿⡿⠟⠁⠀⠀
⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠁⢀⠀⠀⠄⠀⡀⠀⠁⠐⠒⠂⠠⠤⢤⣤⣤⣶⣾⣿⠿⠿⠛⠋⠁⠀⠀⠀⠀⠀
⢇⢃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠁⠀⠀⠈⡄⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠘⡈⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠁⠀⠀⠉⢂⠀⢱⡀⡰⢠⡃⣸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠐⢄⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⠀⠀⠀⠀⢸⠀⠀⠈⠀⠀⠉⡉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠑⢄⡁⠂⠤⢀⣀⠀⠀⠀⠀⣀⣀⣼⣿⠀⠀⠀⠀⣮⣤⣶⣶⣦⠋⢀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠈⠑⠂⠤⠄⠀⠀⠠⠾⠿⠿⢻⣿⠀⠀⢀⣴⣿⣿⣿⡿⠁⡀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣶⡿⢃⠤⠐⠁⠀⢰⣾⠟⡀⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⡏⠀⠆⠀⠀⠀⠀⢸⡟⠀⢡⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⠀⢰⠀⠀⠀⠀⠀⠀⣇⠀⠈⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡟⠀⠸⠀⠀⠀⠀⠀⠀⢹⠀⠀⢁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⡇⠀⡀⠀⠀⠀⠀⠀⠀⢸⣇⠀⠘⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⣿⡇⠀⡇⠀⠀⠀⠀⠀⠀⢸⣿⡀⠀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⡿⠁⣆⠁⠀⠀⠀⠀⠀⠀⠀⣿⣷⢠⡸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠉⠁⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
"""

def fail(message):
    print(ERROR_ART)
    print(f"\n\U0001F4A5 {message}")
    sys.exit(1)

def succeed(message):
    print(SUCCESS_ART)
    print(f"\n\u2705 {message}")

def log(message, icon="\U0001F527"):
    print(f"{icon} {message}")

def ensure_outdir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        log(f"Created output directory: {path}", icon="\U0001F4C1")

def require_pandas():
    try:
        import pandas as pd
        return pd
    except ImportError:
        fail("Missing Python dependency: pandas. Install with `pip install pandas` or conda.")

def require_matplotlib():
    try:
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        return cm, mcolors
    except ImportError:
        fail("Missing Python dependency: matplotlib. Install with `pip install matplotlib` or conda.")

def slugify_name(name):
    safe = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in name.strip())
    return safe.strip("._") or "species"

def copy_to_workdir(path, outdir, label):
    dest = os.path.join(outdir, f"{label}_{os.path.basename(path)}")
    if os.path.abspath(path) != os.path.abspath(dest):
        shutil.copyfile(path, dest)
    return dest

# === COLOR MAP ===
ROMAN_NUMERAL_MAP = {
    'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8,
    'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15,
    'XVI': 16, 'XVII': 17, 'XVIII': 18, 'XIX': 19, 'XX': 20
}

def try_parse_chr(chr_label):
    try:
        return int(chr_label)
    except ValueError:
        return ROMAN_NUMERAL_MAP.get(chr_label.upper(), chr_label)

#generate color gradient
def generate_color_gradient(num_colors, cmap_name='plasma'):
    cm, mcolors = require_matplotlib()
    cmap = cm.get_cmap(cmap_name, num_colors)
    return [mcolors.to_hex(cmap(i)).replace('#', '') for i in range(num_colors)]

#Generate Color map
def create_chr_color_map(karyo_file, output_file, replace_output_file, cmap_name='plasma'):
    try:
        pd = require_pandas()
        df = pd.read_csv(karyo_file, sep="\t")
        if 'Chr' not in df.columns:
            fail("Karyotype file must contain a 'Chr' column.")
        chromosomes = df['Chr'].astype(str).unique()
        chromosomes_sorted = sorted(chromosomes, key=try_parse_chr)
        colors = generate_color_gradient(len(chromosomes_sorted), cmap_name)

        with open(output_file, 'w') as out, open(replace_output_file, 'w') as rep:
            for i, (chr_label, color) in enumerate(zip(chromosomes_sorted, colors)):
                out.write(f"{chr_label}\t{color}\n")
                rep.write(f"{i+1}\t{color}\n")
        log(f"Wrote chromosome-color map to: {output_file}", icon="\U0001F3A8")
    except Exception as e:
        fail(f"Failed to generate color map: {e}")

def filter_non_integer_chrs(filepath):
    try:
        cleaned = []
        with open(filepath, "r") as f:
            for i, line in enumerate(f):
                parts = line.strip().split("\t")
                if i == 0 or (parts[0].isdigit() and parts[3].isdigit()):
                    cleaned.append(line)
                else:
                    log(f"Removing line {i+1} with non-integer Chr: {line.strip()}", icon="\U0001F9F9")
        with open(filepath, "w") as f:
            f.writelines(cleaned)
        log("Removed rows with non-integer chromosome labels.", icon="\u2705")
    except Exception as e:
        fail(f"Failed to filter non-integer chromosome labels: {e}")

#combine karyotypes
def augment_karyotype(karyo_file, color_map, fill_all=None, size=None, color=None):
    try:
        pd = require_pandas()
        df = pd.read_csv(karyo_file, sep="\t")

        if fill_all:
            df['fill'] = fill_all
        else:
            chr_color_map = pd.read_csv(color_map, sep="\t", header=None, names=["Chr", "fill"])
            df = df.merge(chr_color_map, on="Chr", how="left")

            # Safely resolve fill columns
            fill_x = df['fill_x'] if 'fill_x' in df.columns else pd.Series([None] * len(df))
            fill_y = df['fill_y'] if 'fill_y' in df.columns else pd.Series([None] * len(df))
            fill   = df['fill']    if 'fill' in df.columns    else pd.Series([None] * len(df))

            df['fill'] = fill_x.combine_first(fill_y).combine_first(fill)
            df.drop(columns=[col for col in ['fill_x', 'fill_y'] if col in df.columns], inplace=True)

        df['size'] = size
        df['color'] = color

        # Enforce final column order expected by RIdeogram
        expected_order = ['Chr', 'Start', 'End', 'fill', 'species', 'size', 'color']
        df = df[[col for col in expected_order if col in df.columns]]

        df.to_csv(karyo_file, sep="\t", index=False)
        log(f"🧰 Augmented karyotype file: {karyo_file}", icon="🧰")

    except Exception as e:
        fail(f"Error augmenting {karyo_file}: {e}")

# === BUSCO MERGE ===
def merge_busco(busco1_path, busco2_path, output_path):
    try:
        def resolve_columns(fieldnames, path):
            normalized = {name.strip().lower(): name for name in fieldnames or []}
            aliases = {
                "busco_id": ["# busco id", "busco id", "busco_id", "#busco id"],
                "status": ["status"],
                "sequence": ["contig", "sequence", "scaffold"],
                "start": ["gene start", "start"],
                "end": ["gene end", "end"],
            }
            resolved = {}
            missing = []
            for key, names in aliases.items():
                match = next((normalized[name] for name in names if name in normalized), None)
                if match:
                    resolved[key] = match
                else:
                    missing.append(names[0])
            if missing:
                fail(
                    f"{path} is not a compatible BUSCO-format full table. "
                    f"Missing columns: {', '.join(missing)}"
                )
            return resolved

        def read_busco(path, label):
            data = {}
            with open(path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                columns = resolve_columns(reader.fieldnames, path)
                for row in reader:
                    if row.get(columns["status"], "").strip() != "Complete":
                        continue
                    bid = row[columns["busco_id"]].strip().lower()
                    data.setdefault(bid, {})[label] = row
            return data, columns

        log(f"Using BUSCO table for species 1: {busco1_path}", icon="\U0001F9EC")
        log(f"Using BUSCO table for species 2: {busco2_path}", icon="\U0001F9EC")
        b1, columns1 = read_busco(busco1_path, "1")
        b2, columns2 = read_busco(busco2_path, "2")
        shared = [bid for bid in b1 if bid in b2]
        print(f"🧬 Shared complete BUSCOs: {len(shared)}")

        with open(output_path, "w") as out:
            out.write("Species_1\tStart_1\tEnd_1\tSpecies_2\tStart_2\tEnd_2\tfill\n")
            for bid in shared:
                try:
                    row1 = b1[bid]["1"]
                    row2 = b2[bid]["2"]
                    chr1 = row1.get(columns1["sequence"])
                    chr2 = row2.get(columns2["sequence"])
                    start1 = int(float(row1[columns1["start"]]))
                    end1 = int(float(row1[columns1["end"]]))
                    start2 = int(float(row2[columns2["start"]]))
                    end2 = int(float(row2[columns2["end"]]))
                    if not chr1 or not chr2:
                        continue
                    out.write(f"{chr1}\t{start1}\t{end1}\t{chr2}\t{start2}\t{end2}\tplaceholder\n")
                except Exception as e:
                    log(f"Skipping BUSCO {bid}: {e}", icon="\u26A0")
    except Exception as e:
        fail(f"Failed to merge BUSCO files: {e}")


def apply_chr_replacements(synteny_file, rep1_path, rep2_path):
    try:
        def load_map(path):
            m = {}
            with open(path) as f:
                for line in f:
                    old, new = line.strip().split("\t")
                    m[old] = new
            return m

        rep1 = load_map(rep1_path)
        rep2 = load_map(rep2_path)

        out_lines = []
        with open(synteny_file) as f:
            header = f.readline().strip()
            out_lines.append(header)
            for i, line in enumerate(f, 2):
                parts = line.strip().split("\t")
                if len(parts) != 7:
                    continue
                parts[0] = rep1.get(parts[0], parts[0])
                parts[3] = rep2.get(parts[3], parts[3])
                out_lines.append("\t".join(parts))

        with open(synteny_file, "w") as f:
            for line in out_lines:
                f.write(line + "\n")
        log("Chromosome names replaced using provided maps.", icon="\U0001F504")
    except Exception as e:
        fail(f"Failed to apply chromosome name replacements: {e}")

def replace_fill(input_path, color_map_path, output_path):
    try:
        color_dict = {}
        with open(color_map_path) as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        color_dict[parts[0]] = parts[1]

        with open(input_path) as infile, open(output_path, "w") as outfile:
            for i, line in enumerate(infile):
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if i == 0:
                    outfile.write(line)
                    continue
                if len(parts) != 7:
                    continue
                key = parts[0]
                fill = color_dict.get(key)
                if fill:
                    parts[-1] = fill
                    outfile.write("\t".join(parts) + "\n")
    except Exception as e:
        fail(f"Failed to replace fill values: {e}")

def count_data_rows(path):
    try:
        with open(path) as f:
            return sum(1 for i, line in enumerate(f) if i > 0 and line.strip())
    except Exception as e:
        fail(f"Failed to count rows in {path}: {e}")

def write_and_run_rscript(karyotype_path, synteny_path, working_dir, script_path):
    rscript = shutil.which("Rscript")
    if not rscript:
        fail("Rscript was not found on PATH. Install/load R, or rerun without --plot.")

    r_script = f"""
if (!requireNamespace("RIdeogram", quietly=TRUE)) {{
  if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
  devtools::install_github("TickingClock1992/RIdeogram")
}}
library(RIdeogram)
if (!requireNamespace("rsvg", quietly=TRUE)) install.packages("rsvg")
library(rsvg)

setwd("{working_dir}")
print(getwd())
print(list.files())

kary <- read.table("{os.path.basename(karyotype_path)}", header=TRUE)
synt <- read.table("final_synteny.txt", header=TRUE, colClasses = c("numeric", "integer", "integer", "numeric", "integer", "integer", "character"))

ideogram(karyotype=kary, synteny=synt)
rsvg_png("chromosome.svg", "chromosome.png", width=1000)
"""
    with open(script_path, "w") as f:
        f.write(r_script)

    result = subprocess.run([rscript, script_path], capture_output=True)
    if result.returncode != 0 or not os.path.exists(os.path.join(working_dir, "chromosome.png")):
        fail(f"R plotting failed: {result.stderr.decode()}")
    else:
        log("RIdeogram plot created successfully!", icon="\U0001F3A8")

def run_pairwise_pipeline(karyotype1, karyotype2, busco1, busco2, rep1, rep2, outdir,
                          cmap="plasma", plot=False, rscript_path="plot_ideogram.R",
                          karyo_size="12", karyo_color="black"):
    pd = require_pandas()
    ensure_outdir(outdir)

    karyotype1_work = copy_to_workdir(karyotype1, outdir, "species1_karyotype")
    karyotype2_work = copy_to_workdir(karyotype2, outdir, "species2_karyotype")

    chr_color_map = os.path.join(outdir, "chr_color_map.txt")
    color_replace_map = os.path.join(outdir, "color_replace.txt")
    merged_busco = os.path.join(outdir, "merged_busco.txt")
    final_synteny = os.path.join(outdir, "final_synteny.txt")

    create_chr_color_map(karyotype1_work, chr_color_map, color_replace_map, cmap_name=cmap)
    augment_karyotype(karyotype1_work, chr_color_map, size=karyo_size, color=karyo_color)
    augment_karyotype(karyotype2_work, chr_color_map, fill_all="cccccc", size=karyo_size, color=karyo_color)

    df1 = pd.read_csv(karyotype1_work, sep="\t")
    df2 = pd.read_csv(karyotype2_work, sep="\t")
    combined = pd.concat([df1, df2], ignore_index=True)
    combined_karyo_path = os.path.join(outdir, "dual_karyotype.txt")
    combined.to_csv(combined_karyo_path, sep="\t", index=False)
    log(f"Wrote combined karyotype: {combined_karyo_path}", icon="\U0001F91D")

    merge_busco(busco1, busco2, merged_busco)
    apply_chr_replacements(merged_busco, rep1, rep2)
    filter_non_integer_chrs(merged_busco)
    replace_fill(merged_busco, color_replace_map, final_synteny)

    synteny_rows = count_data_rows(final_synteny)
    if synteny_rows == 0:
        log(
            f"No synteny rows were written to {final_synteny}. "
            "Check that both species used the same BUSCO lineage/table format.",
            icon="\u26A0"
        )
        if plot:
            log("Skipping RIdeogram plot because there are no synteny links.", icon="\u26A0")
        return final_synteny

    if plot:
        write_and_run_rscript(
            karyotype_path=combined_karyo_path,
            synteny_path=final_synteny,
            working_dir=outdir,
            script_path=os.path.join(outdir, rscript_path)
        )

    return final_synteny

def open_fasta(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def fasta_lengths(path):
    lengths = OrderedDict()
    current = None
    with open_fasta(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in lengths:
                    fail(f"Duplicate FASTA sequence id '{current}' in {path}")
                lengths[current] = 0
            elif current is None:
                fail(f"FASTA sequence data found before first header in {path}")
            else:
                lengths[current] += len(line)
    if not lengths:
        fail(f"No FASTA sequences found in {path}")
    return lengths

def make_karyotype_and_replacement(assembly_path, species_name, outdir, min_contig_length=0):
    ensure_outdir(outdir)
    lengths = fasta_lengths(assembly_path)
    rows = [(contig, length) for contig, length in lengths.items() if length >= min_contig_length]
    if not rows:
        fail(f"No contigs in {assembly_path} passed --min_contig_length {min_contig_length}")

    rows.sort(key=lambda item: item[1], reverse=True)
    karyotype_path = os.path.join(outdir, f"{slugify_name(species_name)}_karyotype.txt")
    replacement_path = os.path.join(outdir, f"{slugify_name(species_name)}_replacement.txt")

    with open(karyotype_path, "w") as karyo, open(replacement_path, "w") as rep:
        karyo.write("Chr\tStart\tEnd\tspecies\n")
        for i, (contig, length) in enumerate(rows, 1):
            karyo.write(f"{i}\t0\t{length}\t{species_name}\n")
            rep.write(f"{contig}\t{i}\n")

    log(f"Generated karyotype and replacement map for {species_name}", icon="\U0001F9EC")
    return karyotype_path, replacement_path

def parse_named_path(value, option_name):
    if "=" not in value:
        fail(f"{option_name} must use NAME=PATH format: {value}")
    name, path = value.split("=", 1)
    name = name.strip()
    path = path.strip()
    if not name or not path:
        fail(f"{option_name} must use NAME=PATH format: {value}")
    return name, path

def run_command(command, cwd=None):
    log("Running: " + " ".join(command), icon="\U0001F680")
    result = subprocess.run(command, cwd=cwd)
    if result.returncode != 0:
        fail(f"Command failed with exit code {result.returncode}: {' '.join(command)}")

def run_compleasm_for_species(species_name, assembly_path, outdir, args):
    species_dir = os.path.join(outdir, slugify_name(species_name))
    ensure_outdir(species_dir)

    command = shlex.split(args.compleasm) + ["run", "-a", assembly_path, "-o", species_dir, "-t", str(args.threads)]
    if args.autolineage:
        command.append("--autolineage")
    else:
        command.extend(["-l", args.lineage])
    if args.compleasm_library:
        command.extend(["-L", args.compleasm_library])
    if args.odb:
        command.extend(["--odb", args.odb])
    if args.miniprot_execute_path:
        command.extend(["--miniprot_execute_path", args.miniprot_execute_path])
    if args.hmmsearch_execute_path:
        command.extend(["--hmmsearch_execute_path", args.hmmsearch_execute_path])
    if args.sepp_execute_path:
        command.extend(["--sepp_execute_path", args.sepp_execute_path])
    if args.retrocopy:
        command.append("--retrocopy")

    if args.reuse_compleasm:
        existing = find_compleasm_busco_table(species_dir)
        if existing:
            log(f"Reusing existing compleasm output for {species_name}: {existing}", icon="\U0001F501")
            return existing

    run_command(command)
    return require_compleasm_busco_table(species_dir)

def find_compleasm_busco_table(species_dir):
    priority_names = [
        "full_table_busco_format.tsv",
        "full_table_BUSCO.tsv",
        "full_table.tsv",
        "full_table.txt",
    ]
    candidates = []
    for root, _, files in os.walk(species_dir):
        for filename in files:
            lower = filename.lower()
            if filename in priority_names or ("full_table" in lower and lower.endswith((".tsv", ".txt"))):
                path = os.path.join(root, filename)
                try:
                    priority = priority_names.index(filename)
                except ValueError:
                    priority = len(priority_names)
                candidates.append((priority, len(path), path))
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][2]

def require_compleasm_busco_table(species_dir):
    table = find_compleasm_busco_table(species_dir)
    if not table:
        fail(f"Could not find a compleasm full_table file under {species_dir}")
    log(f"Selected compleasm BUSCO table: {table}", icon="\U0001F4CB")
    return table

def build_parser():
    parser = argparse.ArgumentParser(
        description="Synk synteny plotting from compleasm/BUSCO outputs, with optional one-shot compleasm wrapper mode."
    )

    legacy = parser.add_argument_group("pairwise Synk inputs")
    legacy.add_argument("--karyotype1")
    legacy.add_argument("--karyotype2")
    legacy.add_argument("--busco1")
    legacy.add_argument("--busco2")
    legacy.add_argument("--rep1")
    legacy.add_argument("--rep2")

    wrapper = parser.add_argument_group("one-shot compleasm wrapper inputs")
    wrapper.add_argument("--main_name")
    wrapper.add_argument("--main_assembly")
    wrapper.add_argument("--main_karyotype")
    wrapper.add_argument("--main_replacement")
    wrapper.add_argument("--compare", action="append", default=[], metavar="NAME=ASSEMBLY",
                         help="Comparison species assembly. Repeat for multiple species.")
    wrapper.add_argument("--compare_karyotype", action="append", default=[], metavar="NAME=KARYOTYPE")
    wrapper.add_argument("--compare_replacement", action="append", default=[], metavar="NAME=REPLACEMENT")
    wrapper.add_argument("--lineage", help="BUSCO lineage for compleasm, e.g. sauropsida or eukaryota.")
    wrapper.add_argument("--autolineage", action="store_true")
    wrapper.add_argument("--threads", type=int, default=8)
    wrapper.add_argument("--compleasm", default="compleasm")
    wrapper.add_argument("--compleasm_library")
    wrapper.add_argument("--odb")
    wrapper.add_argument("--miniprot_execute_path")
    wrapper.add_argument("--hmmsearch_execute_path")
    wrapper.add_argument("--sepp_execute_path")
    wrapper.add_argument("--retrocopy", action="store_true")
    wrapper.add_argument("--reuse_compleasm", action="store_true")
    wrapper.add_argument("--min_contig_length", type=int, default=0)

    shared = parser.add_argument_group("shared outputs and plotting")
    shared.add_argument("--outdir", required=True)
    shared.add_argument("--cmap", default="plasma")
    shared.add_argument("--plot", action="store_true")
    shared.add_argument("--rscript_path", default="plot_ideogram.R")
    shared.add_argument("--karyo_size", default="12")
    shared.add_argument("--karyo_color", default="black")
    return parser

def is_wrapper_mode(args):
    return bool(args.main_name or args.main_assembly or args.compare)

def validate_pairwise_args(args):
    missing = [
        name for name in ("karyotype1", "karyotype2", "busco1", "busco2", "rep1", "rep2")
        if not getattr(args, name)
    ]
    if missing:
        fail("Missing required pairwise arguments: " + ", ".join(f"--{name}" for name in missing))

def validate_wrapper_args(args):
    missing = []
    if not args.main_name:
        missing.append("--main_name")
    if not args.main_assembly:
        missing.append("--main_assembly")
    if not args.compare:
        missing.append("--compare")
    if not args.autolineage and not args.lineage:
        missing.append("--lineage or --autolineage")
    if missing:
        fail("Missing required wrapper arguments: " + ", ".join(missing))

def run_wrapper(args):
    validate_wrapper_args(args)
    ensure_outdir(args.outdir)

    compare_assemblies = OrderedDict(parse_named_path(item, "--compare") for item in args.compare)
    compare_karyotypes = dict(parse_named_path(item, "--compare_karyotype") for item in args.compare_karyotype)
    compare_replacements = dict(parse_named_path(item, "--compare_replacement") for item in args.compare_replacement)

    karyo_dir = os.path.join(args.outdir, "karyotypes")
    compleasm_dir = os.path.join(args.outdir, "compleasm")
    pairwise_dir = os.path.join(args.outdir, "pairwise")
    ensure_outdir(karyo_dir)
    ensure_outdir(compleasm_dir)
    ensure_outdir(pairwise_dir)

    main_busco = run_compleasm_for_species(args.main_name, args.main_assembly, compleasm_dir, args)
    if args.main_karyotype:
        if not args.main_replacement:
            fail("--main_replacement is required when --main_karyotype is supplied")
        main_karyotype, main_replacement = args.main_karyotype, args.main_replacement
    else:
        main_karyotype, main_replacement = make_karyotype_and_replacement(
            args.main_assembly, args.main_name, karyo_dir, args.min_contig_length
        )

    outputs = []
    for compare_name, assembly_path in compare_assemblies.items():
        compare_busco = run_compleasm_for_species(compare_name, assembly_path, compleasm_dir, args)
        if compare_name in compare_karyotypes:
            if compare_name not in compare_replacements:
                fail(f"--compare_replacement {compare_name}=PATH is required with --compare_karyotype {compare_name}=PATH")
            compare_karyotype = compare_karyotypes[compare_name]
            compare_replacement = compare_replacements[compare_name]
        else:
            compare_karyotype, compare_replacement = make_karyotype_and_replacement(
                assembly_path, compare_name, karyo_dir, args.min_contig_length
            )

        outdir = os.path.join(pairwise_dir, f"{slugify_name(args.main_name)}_vs_{slugify_name(compare_name)}")
        final_synteny = run_pairwise_pipeline(
            main_karyotype, compare_karyotype, main_busco, compare_busco,
            main_replacement, compare_replacement, outdir,
            cmap=args.cmap, plot=args.plot, rscript_path=args.rscript_path,
            karyo_size=args.karyo_size, karyo_color=args.karyo_color
        )
        outputs.append(final_synteny)

    return outputs

def main():
    parser = build_parser()
    args = parser.parse_args()

    if is_wrapper_mode(args):
        outputs = run_wrapper(args)
        succeed("Wrapper pipeline complete! Final synteny files:\n" + "\n".join(outputs))
    else:
        validate_pairwise_args(args)
        final_synteny = run_pairwise_pipeline(
            args.karyotype1, args.karyotype2, args.busco1, args.busco2, args.rep1, args.rep2,
            args.outdir, cmap=args.cmap, plot=args.plot, rscript_path=args.rscript_path,
            karyo_size=args.karyo_size, karyo_color=args.karyo_color
        )
        succeed(f"Pipeline complete! Final synteny file: {final_synteny}")

if __name__ == "__main__":
    main()
