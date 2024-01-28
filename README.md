# NDUFAF6 DMS analysis pipelines

## Download raw data
1. Download raw sequence data from `https://www.ncbi.nlm.nih.gov/bioproject/1007392` and make sure the files are in the `raw-data/` directory.
2. The following 30 files should be in the `raw-data/` directory:
    - `1_ko1a_p0.R1.fastq.gz`
    - `1_ko1a_p0.R2.fastq.gz`
    - `2_ko1a_p3.R1.fastq.gz`
    - `2_ko1a_p3.R2.fastq.gz`
    - `3_ko1a_p6.R1.fastq.gz`
    - `3_ko1a_p6.R2.fastq.gz`
    - `4_ko1b_p0.R1.fastq.gz`
    - `4_ko1b_p0.R2.fastq.gz`
    - `5_ko1b_p3.R1.fastq.gz`
    - `5_ko1b_p3.R2.fastq.gz`
    - `6_ko1b_p6.R1.fastq.gz`
    - `6_ko1b_p6.R2.fastq.gz`
    - `7_ko2a_p0.R1.fastq.gz`
    - `7_ko2a_p0.R2.fastq.gz`
    - `8_ko2a_p3.R1.fastq.gz`
    - `8_ko2a_p3.R2.fastq.gz`
    - `9_ko2a_p6.R1.fastq.gz`
    - `9_ko2a_p6.R2.fastq.gz`
    - `10_ko2b_p0.R1.fastq.gz`
    - `10_ko2b_p0.R2.fastq.gz`
    - `11_ko2b_p3.R1.fastq.gz`
    - `11_ko2b_p3.R2.fastq.gz`
    - `12_ko2b_p6.R1.fastq.gz`
    - `12_ko2b_p6.R2.fastq.gz`
    - `13_ko2c_p0.R1.fastq.gz`
    - `13_ko2c_p0.R2.fastq.gz`
    - `14_ko2c_p3.R1.fastq.gz`
    - `14_ko2c_p3.R2.fastq.gz`
    - `15_ko2c_p6.R1.fastq.gz`
    - `15_ko2c_p6.R2.fastq.gz`

## Read processing and variant counting:
1. Install `Bowtie2` v2.4.4 and `FLASH` v1.2.11.
1. Create Conda environment using provided `environment.yml` with `conda env create -n af6-dms -f environment.yml`.
2. Activate Conda environment with `conda activate af6-dms`
3. Navigate to `ngs-data/` directory.
4. Execute Snakemake pipeline using `Snakemake` (set number of cores with the `--cores` flag).

## Compile variant counts into variant count table
1. Run the `compile_counts.sh` script with `sh compile_counts.sh` to compile variant and wild-type read counts into variant count table needed for DiMSum analysis.
    - `combine_variant_counts.py` produces the `variant_counts.txt` file.
    - `wt_counts.py` produces the `indexed_wt_counts.txt` file.
    - `dimsum_preprocess.py` produces the `dimsum/variant_count_table.txt` file.

## Run DiMSum analysis
1. Create and activate new Conda environment for the DiMSum package according to https://github.com/lehner-lab/DiMSum.
2. Navigate to the `dimsum/` directory.
3. Run the `run_dimsum.sh` script with `sh run_dimsum.sh`.
4. Activate the `af6-dms` environment again with `conda activate af6-dms`.
5. Run `dimsum_postprocess.py` with `python dimsum_postprocess.py` to produce the final `dms_fitness.txt` file in the `data` directory.

## Gaussian Mixture Model analysis
The Jupyter notebook file `gmm_analysis.ipynb` contains the code used to generate the Gaussian Mixture Model and subsequent calculations of probability of abnormal function. Running this notebook will generate the `dms_fitness_gmm.txt` file, which corresponds to Supplemental Table S1.
