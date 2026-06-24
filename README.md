
# Synk - Rapid synteny plotting tool using [Compleasm](https://github.com/huangnengCSU/compleasm) outputs
![synk.pdf](https://github.com/user-attachments/files/27210094/synk.pdf)

Synk is a tool to plot synteny using BUSCO genes. With BUSCO genes, synteny is plotted rapidly and avoids the issue of mapping paralogs, sinch BUSCO genes are conserved, single copy orthologs. This makes Synk an effective tool for analysing synteny between distantly related species. 

## Installation
Synk is a python script. Download it!
```
git clone https://github.com/jomhoff/Synk.git
cd Synk
```
Dependencies:
  Pandas;
  Matplotlib;
  Rideogram

See [DEPENDENCIES.md](DEPENDENCIES.md) for full Python, compleasm, miniprot, hmmsearch, R, RIdeogram, and SLURM setup notes.

Recommended conda setup:
```
conda env create -f environment.yml
conda activate synk
Rscript -e 'remotes::install_github("TickingClock1992/RIdeogram")'
```
  
## Use

### One-shot compleasm wrapper

Synk can now run compleasm first and then make every pairwise plot from one main species to any number of comparison species.

Example with automatic karyotypes generated from FASTA sequence lengths:
```
python synk.py \
  --main_name pfas \
  --main_assembly pfas.fa \
  --compare tiliqua=tiliqua.fa \
  --compare egernia=egernia.fa \
  --lineage sauropsida \
  --threads 16 \
  --outdir synk_output \
  --plot
```

Use `--compare NAME=ASSEMBLY` once for each comparison species. Synk will write:
```
synk_output/
  compleasm/        # one compleasm run per species
  karyotypes/       # auto-generated karyotype and replacement files
  pairwise/         # one Synk result folder per main-vs-comparison pair
```

If compleasm has already been run, add `--reuse_compleasm` and Synk will look for an existing `full_table` file under each species output folder before rerunning compleasm.

You can also provide curated karyotypes and replacement maps. This is useful when you want chromosome labels or ordering to follow a published karyotype instead of FASTA length order:
```
python synk.py \
  --main_name pfas \
  --main_assembly pfas.fa \
  --main_karyotype pfas_karyotype.txt \
  --main_replacement pfas_replacement.txt \
  --compare tiliqua=tiliqua.fa \
  --compare_karyotype tiliqua=tiliqua_karyotype.txt \
  --compare_replacement tiliqua=tiliqua_replacement.txt \
  --lineage sauropsida \
  --outdir synk_output \
  --plot
```

Wrapper inputs:
```
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
```

Optional curated karyotype inputs:
```
--main_karyotype                 Main karyotype file
--main_replacement               Main chromosome replacement map
--compare_karyotype NAME=FILE    Comparison karyotype; repeat as needed
--compare_replacement NAME=FILE  Comparison replacement map; repeat as needed
```

### Pairwise mode

Usage example:
```
python synk.py \
  --karyotype1 pfas_karyotype.txt \
  --karyotype2 tiliqua_karyotype.txt \
  --busco1 pfas_busco_format.tsv \
  --busco2 tiliqua_busco_format.tsv \
  --rep1 replacement_pfas.txt \
  --rep2 replacement_tiliqua.txt \
  --outdir output_directory \
  --cmap viridis \
  --plot
```

Required inputs:
```
--karyotype1       First species karyotype file (must have columns: Chr, Start, End, Species)
--karyotype2       Second species karyotype file (same format as above)
--busco1           BUSCO full table from species 1 (full_table_BUSCO.tsv)
--busco2           BUSCO full table from species 2
--rep1             Replacement map for species 1 chromosome names (ie. chr1	1)
--rep2             Replacement map for species 2 chromosome names
--outdir           Directory for output files
```

Optional inputs:
```
--cmap             Colormap for chromosome gradient (default: plasma; other: viridis, magma)
--plot             Include this flag to generate an ideogram PNG/SVG with RIdeogram
--rscript_path     Path to output R script (default: plot_ideogram.R)
--karyo_size       Label text size (default: 5)
--karyo_color      Label text color (default: black)
```

Output Files:
```
- chr_color_map.txt      : Map of chromosome names to colors
- color_replace.txt      : Map of numeric chromosome order to colors
- merged_busco.txt       : Merged BUSCO hits shared between both species
- final_synteny.txt      : Synteny links colored by species 1 chromosome
- dual_karyotype.txt     : Combined karyotype file with label and fill info
- chromosome.png/.svg    : RIdeogram synteny plot (if --plot is used)
```

## Result
<img width="1209" alt="pfas-tiliqua_synteny" src="https://github.com/user-attachments/assets/2d7e9f15-501e-4991-a400-17588c6a1cca" />
This is the result (with labels edited slightly in illustrator)

## Citation
If you use Synk in your work, please cite: 

Hoffman, J.J., Burbrink, F.T., Pyron, R.A., Raxworthy, C.J., 2025. Telomere-to-telomere reference genome of the common five-lined skink, Plestiodon fasciatus (Squamata: Scincidae). https://doi.org/10.1101/2025.07.03.663019
