# Synk Dependencies

This file lists the programs needed to run Synk in wrapper mode and pairwise mode.

## Required For All Synk Runs

Synk itself is a Python script. The Python environment needs:

```bash
python >= 3.7
pandas
matplotlib
```

## Recommended Conda Environment

The easiest setup is to create the environment from the included file:

```bash
conda env create -f environment.yml
conda activate synk
```

Then install the RIdeogram R package:

```bash
Rscript -e 'remotes::install_github("TickingClock1992/RIdeogram")'
```

Check the full environment:

```bash
python -c "import pandas, matplotlib; print('Python dependencies OK')"
which compleasm
which miniprot
which hmmsearch
which Rscript
Rscript -e 'library(RIdeogram); library(rsvg); cat("R plotting dependencies OK\n")'
```

Install with conda:

```bash
conda create -n synk -c conda-forge python=3.10 pandas matplotlib
conda activate synk
```

Or add the packages to an existing environment:

```bash
conda install -c conda-forge pandas matplotlib
```

Check that Python can see the packages:

```bash
python -c "import pandas, matplotlib; print('Python dependencies OK')"
```

## Required For Wrapper Mode

Wrapper mode runs compleasm for each genome before making Synk pairwise outputs.

Required programs:

```bash
compleasm
miniprot
hmmsearch
```

Compleasm usually finds `miniprot` and `hmmsearch` from the same conda environment. Install with:

```bash
conda create -n compleasm -c conda-forge -c bioconda python=3.10 compleasm pandas matplotlib
conda activate compleasm
```

If your cluster has trouble with newer BUSCO/ODB filenames, use compleasm `0.2.8` or newer.

Check the install:

```bash
which compleasm
which miniprot
which hmmsearch
compleasm --version
```

## BUSCO Lineage Data

Synk passes the lineage to compleasm with:

```bash
--lineage sauropsida
```

Compleasm will download the lineage data automatically unless you provide a local library path.

To use a local compleasm library:

```bash
python synk.py \
  --compleasm compleasm \
  --compleasm_library /path/to/compleasm/library \
  --lineage sauropsida \
  ...
```

If you already ran compleasm once, rerun Synk with:

```bash
--reuse_compleasm
```

Synk will look for compleasm full-table outputs and prefers:

```text
full_table_busco_format.tsv
```

## Optional Plotting Dependencies

Synk can make `final_synteny.txt` without R. R is only required when you use:

```bash
--plot
```

Plotting requires:

```bash
R
Rscript
RIdeogram
XML
rsvg
grImport2
devtools
```

Install R packages from conda where possible:

```bash
conda install -c conda-forge r-base r-xml r-rsvg r-grimport2 r-devtools r-remotes
```

Then install RIdeogram:

```bash
Rscript -e 'if (!requireNamespace("RIdeogram", quietly=TRUE)) remotes::install_github("TickingClock1992/RIdeogram")'
```

If RIdeogram fails with missing `XML`, `rsvg`, or `grImport2`, install those dependencies with conda rather than letting R compile them:

```bash
conda install -c conda-forge r-xml r-rsvg r-grimport2
Rscript -e 'remotes::install_github("TickingClock1992/RIdeogram", upgrade="never")'
```

Check plotting dependencies:

```bash
which Rscript
Rscript -e 'library(RIdeogram); library(rsvg); cat("R plotting dependencies OK\n")'
```

If `Rscript` is missing, run Synk without `--plot`. The synteny text files will still be created.

## Recommended SLURM Environment

Example SLURM header:

```sh
#!/bin/sh
#SBATCH --job-name=synk
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=200gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.edu

source ~/.bashrc
conda activate compleasm
```

Example Synk wrapper run:

```sh
python synk.py \
  --main_name pfas \
  --main_assembly /path/to/main.fasta \
  --compare emoia=/path/to/emoia.fasta \
  --lineage sauropsida \
  --threads ${SLURM_TASKS_PER_NODE:-8} \
  --reuse_compleasm \
  --outdir /path/to/synk_output
```

Add `--plot` only if `Rscript`, `RIdeogram`, and `rsvg` are available.

## Quick Troubleshooting

If compleasm fails while downloading BUSCO data, update compleasm:

```bash
conda install -c conda-forge -c bioconda compleasm
```

If compleasm/miniprot fails with `failed to open/build the index`, check the FASTA:

```bash
ls -lh genome.fasta
file genome.fasta
head -n 5 genome.fasta
grep -c '^>' genome.fasta
```

If Synk reports `Shared complete BUSCOs: 0`, check that both species used the same lineage and that Synk selected `full_table_busco_format.tsv`:

```bash
find synk_output/compleasm -name '*full_table*' -print
head -n 3 synk_output/compleasm/*/*/full_table_busco_format.tsv
```

If plotting fails with `Rscript was not found on PATH`, either install/load R or rerun without `--plot`.
