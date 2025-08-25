**A reproducible pipeline to infer quartet gene trees and infer the distribution of possible topologies**.
**v1.0**
[TOC]

## Introduction

**nf-GTs** is a bioinformatics pipeline that takes a reference genome, variant calls and mask files, generate consensus sequences and infers the distribution of quartet trees.

The pipeline is built using [Nextflow](https://www.nextflow.io) (DSL1), a workflow tool to run tasks across multiple computing infrastructures in a very portable manner. It currently creates a custom Conda environment from scratch and therefore does not require any pre-compiled software packages (besides Nextflow itself). Future development may include containerized versions as well to further enhance reproducibility. If possible, e.g. when running on a HPC cluster, the pipeline will process alignments, infer trees, etc. in parallel and all batch submission jobs are handled internally through the Nextflow workflow.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/) (version >= 19.04) 
2. Install [`Conda`](https://conda.io/miniconda.html) (version >= 4.10)
3. Download the pipeline and adjust the config file. Here you need to specify the location of the input data, filtering settings etc.

    ```bash
    ./nf-GTs/nextflow.config
    ```
4. Create a nextflow config profile that matches your cluster set-up ( [`profile`]( https://www.nextflow.io/docs/latest/config.html#config-profiles) and start running your own analysis!

    ```bash
    nextflow run ./nf-GTs/main.nf -profile mfn
    ```

6. Once your run has completed successfully, clean up the intermediate files.

    ```bash
    nextflow clean -f -k
    ```

If a run stalls for a given reason, inspect the error message, adjust and you can rerun using the 'resume' flag.

    ```bash
    nextflow run ./nf-GTs/main.nf -profile mfn -resume
    ```


**N.B.** Please make sure that:
* The scripts in `./bin/xxx.py` can be executed by any user. You can check the user permissions via:

    ```bash
    ls -l ./bin/
    ```
And if need be changed via:

    ```bash
    cd ./bin/
    chmod 755 *
    ```

* There is **sufficient storage capacity** on the location where you will execute the pipeline from (and where the nextflow work directory will be stored). On the MfN cluster, running on `/home/` will easily lead to problems and all pipelines will need to be executed from a project folder on main storage.


## Pipeline Summary

### Default Steps

By default the pipeline currently performs the following:

* Call consensus sequences for each chromosome/scaffold and each individual (`samtools`, `bcftools`)
* Merge individual consensus sequences into multiple sequence alignments per chromosome/scaffold 
* Window based subsetting of chromosomes/scaffolds and filtering (`bin/01_window_subset_MSA.py`)
* Phylogenetic inference using maximum-likelihood, for each window (`IQtree2`)
* Topology frequency estimation and visualization for each comparison: indiv_P1 x indiv_P2 x indiv_P3 x indiv_OUT (`bin/02_calc_topo_freqs.py`)
* Summarise quartet frequencies across all comparisons (`bin/03_topo_comp_summary.py`)


### Input

The pipeline can start with consensus calling (each individual, each scaffold/chromosome) or with consensus sequences themselves. If starting with consensus calling, the pipeline requires a `variant call set` and a `mask set` for each chromosome and stored in a single folder per individual:

variant call set = `/inputdir/[indiv]/[chromo]_vars_filt_indels.vcf.gz`
mask set = `/inputdir/[indiv]/[chromo]_mask_cov_het.vcf`

- A text file with all the chromosomes to be included, e.g.:

chr1
chr2
chr3
etc.

- A tab-delim file with information regarding the length of chromosomes and their relative start and end positions (for plotting purposes)

| chromo | start |     end    |    mid   | start_rel |   end_rel  |  mid_rel  | 
| :----: | :----:| :--------: | :------: | :--------:| :--------: | :-------: | 
| chr1   | 0     |  119122530 | 59561265 | 0         |  119122530 | 59561265  | 
| chr1A  | 0     |  74173823  | 37086912 | 119122530 |  193296353 | 156209442 | 
| chr2   | 0     |  152213986 | 76106993 | 193296353 |  345510339 | 269403346 | 


Run specific options can  be specified in the `./nextflow.config` script. All other values for programs are set at default. Boolean values can be specified as `true` or `false` and path names need to be absolute. If preferred, the same options listed in the config file can also be directly modified by using a parameter flag when initiating the nextflow run. Example given:

    ```bash
    nextflow run nf-phylo/phylo.nf -profile mfn --indivs_file /path/to/indivs.txt --chromos_file /path/to/chromos.txt --outdir /path/to/results --consensus_calling true
    ```

**NOTE**:
- NF-GTs only does a 'pseudo-alignment', meaning that it assumes that the consensus sequences for each chromosome/scaffold is the same length across all individuals. Indel variation should therefore have been removed from the variant call set or not incorporated into the existing sequences.
- Partially motivated for the required above, variant and mask sets need to be called on the SAME reference assembly across all individuals.
- Sex chromosomes should be treated differently due to their distint inheritance history. Such chromosomes should either be excluded from the (chromo) input list or assigned in the `sex_chromos` parameter.
- See environment.yml for all the required software AND also some suggestion on how to bypass the incompatibility between iqtree2 and ete3 in a single conda environment


### Output
The pipeline returns both tsv summary files and plots

