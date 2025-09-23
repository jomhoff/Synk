# Synk - Rapid synteny plotting tool using [Compleasm](https://github.com/huangnengCSU/compleasm) outputs
![Synk logo](https://github.com/user-attachments/assets/e3af80d2-1236-4655-88d8-c4b27855fba0)

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
  
## Use

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

