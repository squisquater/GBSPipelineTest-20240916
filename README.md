# Population Genetics Pipeline

## Overview

This Snakemake workflow automates the process of aligning genomic reads, identifying SNPs, performing population structure analysis, and conducting landscape genetics analysis. The workflow is modular and can be run in parts or as a complete pipeline.

*IMPORTANT: This pipeline does rely on the prior setup of conda or mamba as well as a snakemake slurm profile (which allows snakemake to submit jobs to the different computing nodes on an hpc)*

## Key Features

- **Read Alignment and Quality Control**: Align paired-end reads to a reference genome, add read group information, and generate quality control metrics.
- **SNP Calling**: Use the Stacks `ref_map` and `populations` pipelines to call SNPs.
- **Population Genetics**: Perform PCA, MDS, and ADMIXTURE analyses to estimate population structure.
- **Landscape Genetics**: Generate pairwise FST and Nei's D matrices, calculate centroids, and test isolation-by-distance using Mantel tests.
- **Modular Workflow**: Individual snakefiles allow for running different sections independently, but the entire pipeline can also be run from the main snakefile.

## Directory Structure

```plaintext
.
├── 00.snakefile_main
├── README.md
├── accessory_scripts
│   ├── GENHETv3.1.R
│   ├── GenDist.R
│   ├── admixture_metrics.R
│   ├── allelic_balance.py
│   ├── append_population.py
│   ├── calculate_ind_centroids.R
│   ├── calculate_pop_centroids.R
│   ├── compare_replicates.R
│   ├── export_coords.py
│   ├── extract_sample_coords.py
│   ├── feems_spatial_graph.py
│   ├── fit_feems.py
│   ├── generate_concave_hull.py
│   ├── generate_popmap.py
│   ├── genhet_vs_depth.R
│   ├── isolation_by_distance.R
│   ├── mds_plot.R
│   ├── merge_centroids_with_fam.py
│   ├── pca_plot.R
│   ├── ped2genhet.py
│   ├── plot_admixture_barplots.R
│   ├── plot_genhet_barplots.R
│   ├── plot_genhet_vs_depth.R
│   ├── plot_polygon.py
│   ├── rename_chromosomes.py
│   └── summarize_admixture_qvalues.R
├── config
│   ├── 01.snakefile_alignPE.yml
│   ├── 02.snakefile_stacks.yml
│   ├── 04.snakefile_admixture.yml
│   ├── 05.snakefile_genhet.yml
│   ├── 06.snakefile_landgen.yml
│   └── 07.snakefile_feems.yml
├── envs
│   ├── GBS.yml
│   ├── LandGen.yml
│   ├── alphashape.yml
│   ├── bcftools.yml
│   ├── cowplot.yml
│   ├── feems.yml
│   └── snakemake.yml
├── input_files
│   ├── MasterDB_test.txt
│   ├── hexgrid10km.shp
│   ├── hexgrid10km.shx
│   ├── hexgrid5km.shp
│   └── hexgrid5km.shx
└── modules
    ├── admixture
    │   └── 04.snakefile_admixture
    ├── alignPE
    │   └── 01.snakefile_alignPE
    ├── feems
    │   └── 07.snakefile_FEEMS
    ├── genhet
    │   └── 05.snakefile_genhet
    ├── landgen
    │   └── 06.snakefile_landgen
    └── stacks
```
## Getting Started

### 1. Clone the Repository

First, clone the repository from GitHub to your local machine:
```
git clone https://github.com/squisquater/KitFoxGBS-PipelineTest20240916.git

cd KitFoxGBS-PipelineTest20240916
```

### 2. Set Up Conda Environment

Once the repository is cloned, create and activate the Snakemake environment (using conda or mamba):
```
mamba env create -f envs/snakemake.yml

mamba activate snakemake
```
The pipeline will automatically create other necessary Conda environments as it runs, but they are defined in the envs/ directory (e.g., GBS.yml, LandGen.yml).

### 3. Customize Input Files

Place your input files in the input_files/ directory. Key files include:

- **Master Database:** See [MasterDB_test.txt]() as an example but your master database file must contain at least 3 columns ('Individual.ID', 'Library.ID', 'and 'Region.ID') in order to run the alignment and stacks pipelines. If you want to run the LandGen pipeline you will also need a 'Lat' and 'Long' column. 
- **Shapefiles:** Hexagonal grid shapefiles for landscape genetics pipeline (hexgrid5km.shp, hexgrid10km.shp). If your project spans outside California or uses a different region, replace these shapefiles with your own. (NOTE: At the moment feems isn't working in this pipeline anyways so these files aren't super useful)

### 4. Modify Configuration Files

Each module of the pipeline has a corresponding configuration file in the config/ directory. These files define paths to input files, output directories, and specific parameters for each analysis.

## Steps to Modify:

Open the relevant .yml file in the config/ directory.
Modify the paths to your input files, set parameters like depth thresholds, and update other analysis options as needed.

## Key Config Files:

* 01.snakefile_alignPE.yml: Defines paths for read alignment, reference genome, and sample metadata.
* 02.snakefile_stacks.yml: Defines settings for SNP calling, population mapping, and filtering.
* 04.snakefile_admixture.yml: Defines the K-values range for ADMIXTURE analysis.
* 05.snakefile_genhet.yml: Defines parameters for running genetic diversity adn inbreeding stats.
* 06.snakefile_landgen.yml: Defines settings for landscape genetics, including IBD analysis.


Make sure all file paths and parameters match your data and desired analysis.

### 5. Run the Pipeline

After customizing the input files and configuration files, you can run the entire workflow or individual modules.

**Run the entire workflow:**

```
## Dry Run
snakemake -s 00.snakefile_main --profile slurm -n -r 

## Actual Run
snakemake -s 00.snakefile_main --profile slurm --jobs 30
```

**Run individual modules:**

To run a specific analysis (e.g., read alignment, SNP calling), use the corresponding snakefile:

```
## Dry Run
snakemake -s modules/alignPE/01.snakefile_alignPE  --profile slurm -n -r

## Actual Run
snakemake -s modules/alignPE/01.snakefile_alignPE  --profile slurm --jobs 30
```
## Output Files

After running the workflow, the output will include:

- BAM files and alignment metrics.
- VCF files containing SNPs.
- PLINK files for population genetics analysis.
- ADMIXTURE results including Q-value summaries and bar plots.
- Genetic diversity and inbreeding across populations/regions.
- Distance matrices (FST and Nei’s D) and heatmaps for landscape genetics.
- Isolation-by-distance plots from the Mantel test.

Output directories and files are defined in the config files.

## Troubleshooting

1. Ensure that input file paths are correctly specified in the configuration files.
2. Check Snakemake logs for error messages if the workflow fails.

## Acknowledgements

This pipeline integrates multiple tools and custom scripts. If you use this workflow in your research, please cite the following:

- Snakemake: https://snakemake.readthedocs.io
- ADMIXTURE: https://github.com/chrchang/admixture
- Stacks: https://catchenlab.life.illinois.edu/stacks/
- PLINK: https://www.cog-genomics.org/plink/
- BCFtools: http://samtools.github.io/bcftools/bcftools.html
