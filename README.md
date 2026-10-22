# MOVE

Code for *Migration distance and direction show distinct associations with socioeconomic status and health*.

## Structure

```text
analysis/
  00_common/ ... 07_output/  Numbered analysis modules
  reference/        Variable dictionary and lookup tables
  rap/              UKB-RAP launchers
  hpc/              Summary-statistic and atlas workflows
figures/
  01_Main/          Figures 1-6
  02_Extended_Data/ Extended Data Figures 1-10
  03_Supplementary/ Supplementary Figures 1-10
  style.R           Shared plotting functions
```

Analysis folders follow data preparation, genetics, associations, regional models, omics, enrichment and output. `analysis/modules.tsv` lists the commands. Package versions are in the two `versions.tsv` files.

## Inputs

Participant data are not included. `analysis/reference/fields.tsv` lists the required files, variables, field IDs and definitions. Prepare these inputs using approved data access. Feature and reference lists are in the same folder.

Use RAP for UKB participant-level analyses, GWAS and genotype scoring, HPC for summary-statistic and atlas analyses, and an approved local environment for UKHLS/HRS and plotting.

## Analysis

Run from this code directory:

```bash
export MOVE_ROOT=/path/to/workspace
export MOVE_CODE="$PWD/analysis"
bash analysis/run.sh -h
```

Defaults: inputs in `$MOVE_ROOT/data/phe` and `data/ukb`, results in `results/`, source workbooks in `source/`. Override with `MOVE_PHE`, `MOVE_UKB`, `MOVE_RESULTS` or `MOVE_SOURCE`.

Examples, once the required inputs are ready:

```bash
bash analysis/rap/run.sh assoc       # UKB associations
bash analysis/run.sh external       # UKHLS/HRS associations
bash analysis/run.sh tables         # Supplementary Tables
bash analysis/run.sh export Fig2    # Figure source data
```

Run each stage after its inputs are complete. RAP/HPC launchers contain the required tool and input settings. Do not run jobs that write the same result file simultaneously.

## Figures

Each figure script reads its matching `source/` workbook and assets, without participant data or analysis scripts. Use Arial for matching fonts.

```bash
Rscript figures/01_Main/Figure_02.R
Rscript figures/02_Extended_Data/Extended_Data_Figure_08.R
Rscript figures/03_Supplementary/Supplementary_Figure_08.R
```

Outputs go to `output/`. Optional arguments specify another workbook and output folder:

```bash
Rscript figures/01_Main/Figure_02.R /path/to/source.xlsx /path/to/output
```

Keep assets beside the workbooks. Fig1a is editable artwork. ED2 point maps use approved raster layers, with precise coordinates retained in RAP.

Code is released under the MIT License. Third-party data retain their original licences.

