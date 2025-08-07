[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16759028.svg)](https://doi.org/10.5281/zenodo.16759028)

# Installation

In order to get a working environment, we recommend to clone the git repository :

``` bash
git clone https://github.com/gcorre/GNT_GuideSeq
```

Your folder architecture should look similar to :

``` bash
GNT_GuideSeq
├── 00-pipeline/
├── 01-envs/
├── 02-resources/
├── 03-misc/
├── test_dataset/
```

## Conda environments

Install the miniconda environment manager :

``` bash
## https://www.anaconda.com/docs/getting-started/miniconda/install

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash ~/Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc

conda update -n base -c conda-forge conda

conda install mamba -c conda-forge # faster packages manager
```

Install the workflow manager snakemake in the base environment

``` bash
mamba install snakemake openpyxl # worflow manager 
```

Create the environment in your favorite path (-p) from the `environment.yaml` file:

``` bash
cd GNT_GuideSeq/

mamba env create -f 01-envs/environment.yml
```

You should now have an environment named 'guideseq' containing all the required programs when running:

``` bash
#list environments
mamba env list 

# programs version in guideseq environment
mamba list -n guideseq
```

## Reference genomes

The pipeline uses the bowtie2 program to align reads on the reference genome (@langmead2018).

For efficient genome alignment, you can use a pre-built index for Bowtie2. These indices can be downloaded from the [Bowtie2 website](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Using a pre-built index saves computational resources and time, as it eliminates the need to build the index from scratch.

### Build index manualy

If you decide to build the index, we recommend to **not** use haplotypes, scaffolds and unplaced chromosomes to avoid unnecessary multihits alignments. Instead, use the "primary Assembly" version from ensembl or gencode for example. Manually remove unwanted chromosomes if necessary using the seqkit program:

``` bash
conda activate guideseq # load the 'guideseq' environment to access tools

# keep only chromosomes, dicard scaffolds
seqkit grep  -r -p '^[chr]?[0-9XYM]+' Homo_sapiens.GRCh38.dna.primary_assembly.fa > your_clean_fasta
```

Then use the bowtie2 command line:

``` bash
bowtie2-build -@ {threads} {your_clean_fasta} {index_prefix}

conda activate           # return to "base" environment
```

The index prefix path will be used in the configuration file.

## Annotation

### Genes

Off-Target site annotation is performed from a GTF file that can be downloaded from any source (ensembl, gencode ...). The GTF file will be processed to keep only gene and exon features for annotation.

Both the GTF and fasta reference files should be downloaded from the same source for compatibility (especially in chromosome naming).

### Oncogenes

Off-targets sites can be annotated with an oncogene list if provided in the configuration file, otherwise, `NA` will be added.

We provide an example in the `02-ressources/` folder. This file is derived from the [oncoKB cancer gene list](https://www.oncokb.org/cancer-genes). A script is also available in `03-misc` to convert file downloaded from oncoKB website to compatible file format.

User can provide its own annotation file but it must contain the following columns separated by tabs:

-   Ensembl **transcript** ID (ie ENST00000318560, without version)

-   Onco_annotation (collapsed strings, ex : "oncogene \| TSG")

## Prediction tool

The pipeline uses the SWOffinder program (@yaish2024) to predict gRNA off-targets on the fly and annotate OT accordingly as predicted or not.

-   Install from the github repository : <https://github.com/OrensteinLab/SWOffinder>
-   Set the path to SWOffinder in the config file ([Prediction])

# Running the pipeline

In order to run the analysis, 4 elements are mandatory:

1.  The conda environment (see above for installation)
2.  The Sample Data Sheet
3.  The configuration file
4.  The input un-demultiplexed reads in fastq.gz format (R1,R2,I1,I2). *If working with libraries that were already demultiplexed, please read section .....*

## Prepare samples data sheet (SDS) {#prepare-samples-data-sheet-sds}

The sample data-sheet (SDS) is a simple file provided as delited ( ; ) or xlsx format that contains information about each sample to process in the run. An example is proposed in ``` ./test_dataset/sampleInfo.csv and``./test_dataset/sampleInfo.xlsx ```.

| sampleName | CellType | Genome | gRNA_name | gRNA_sequence | orientation | Cas | PAM_sequence | PAM_side | Cut_Offset | type | index1 | index2 |
|------|------|------|------|------|------|------|------|------|------|------|------|------|
| VEGFA_s1_K562_pos | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | positive | Cas9 | NGG | 3 | -4 | guideseq | AGGCAGAA | CTAAGCCT |
| VEGFA_s1_K562_neg | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | negative | Cas9 | NGG | 3 | -4 | guideseq | TCCTGAGC | CTAAGCCT |

: Mandatory columns are:

-   **`sampleName`** : Sample name to use.
    -   Will be use to name output files and folders.
    -   Samples with the same `sampleName` will be merged before processing ([Merging samples])
    -   Samples name [**must not include**]{.underline} '-' in their name as this symbol is used in the pipeline for a special purpose. Instead, use "\_" or any other separator of your choice.
-   **`Genome`**: Uses to define which reference genome to use for each sample.
    -   This **must** be one of the values present in config file **`genome`** key (see [Reference genome]).
-   **`gRNA_name`** : name of the gRNA used.
-   **`gRNA_sequence`**: Sequence of the gRNA [**without**]{.underline} the PAM sequence.
-   **`orientation`** : which PCR orientation was chosen ["positive", "negative", "mix" if they share the same indexes]
    -   (this info is only used for metadata purposes, both PCR orientations are automatically processes by the pipeline).
-   **`Cas`**: name of the Cas used
-   **`PAM_sequence`**: Sequence of the PAM
-   `PAM_side`: [5 or 3] indicating if the PAM is 5' or 3'
-   **`Cut_Offset`**: Distance from gRNA end where the cut occurs [default : -4]
-   **`type`**: Type of experiment ["guideseq", "iguideseq"] or other. This value will define which sequence to trim ([Reads adapter & ODN trimming sequences](#reads-adapter-odn-trimming-sequences)). Type of library in the SDS file must match one identical type in the configuration file.
-   **`index1`**: Sequence of index 1 for demultiplexing
-   **`index2`**: Sequence of index 2 for demultiplexing

Additional columns can be added for metadata annotation purpose.

SDS format is automatically validated using the `snakemake.utils - validate` function.

### Merging samples

In certain scenarios, it may be beneficial to merge different samples from a library during the reads processing stage. This can be particularly useful when dealing with Multiple Replicates or +/-PCR orientations.

To achieve this, you can use the same sample name for multiple rows in the Sample Description Sheet (SDS). Samples that share the same name will be merged during the reads processing, provided they meet the following criteria:

1.  **`Reference Genome`**: The samples must use the same reference genome.

2.  **`gRNA`**: The samples must have the same gRNA, both in terms of sequence and name.

3.  **`Cas Protein`**: The samples must use the same Cas protein, with identical PAM (Protospacer Adjacent Motif) and Offset values.

4.  **`Protocole`**: The samples must use the same ODN (Oligodeoxynucleotide), defined as guideseq or iguideseq.

If any of the above conditions are not met, an error will be raised, and the pipeline will be stopped. This ensures that only compatible samples are merged, maintaining the integrity of the data processing workflow.

Using the test SDS above, if you want to merge both positive and negative libraries, give the same sample name to both rows. As they use the same genome, gRNA, Cas & method, they will be aggregated in a single library.

| sampleName | CellType | Genome | gRNA_name | gRNA_sequence | orientation | Cas | PAM_sequence | PAM_side | Cut_Offset | type | index1 | index2 |
|------|------|------|------|------|------|------|------|------|------|------|------|------|
| VEGFA_s1_K562 | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | positive | Cas9 | NGG | 3 | -4 | guideseq | AGGCAGAA | CTAAGCCT |
| VEGFA_s1_K562 | K562 | GRCh38 | VEGFAs1 | GGGTGGGGGGAGTTTGCTCC | negative | Cas9 | NGG | 3 | -4 | guideseq | TCCTGAGC | CTAAGCCT |

: UMI and reads counts will still be breakdown to each PCR orientation in the final result table.

### Working with already demultiplexed libraries

To be documented

## Prepare the configuration file

The configuration file is a yaml formatted file with key-value dictionary used to fine-tune the pipeline behavior. Settings will apply to [**all samples of the run**]{.underline}.

An example of configuration file is proposed in `./test_dataset/guideSeq_GNT.yml`

[config file must be present in the working directory when starting the pipeline.]{.underline}

### Metadata

Author and affiliation will be printed on the final report.

``` yaml
author: "Guillaume CORRE, PhD"
affiliation: "Therapeutic Gene Editing - GENETHON, INSERM U951, Paris-Saclay University, Evry, France"
contact: "gcorre@genethon.fr"
version: 'V1.0'
```

### Remove intermediate files

Whether to remove all intermediate files after run completion.

``` yaml
clean_intermediates_files: 'TRUE'
```

### Path to Sample Data Sheet and read files

SDS path is relative to current folder (from where the pipeline is started). Absolute path can be used.

``` yaml
## Library informations
sampleInfo_path: "sampleInfo.csv"
read_path: "."  # path to reads if not current folder
R1: "Undetermined_S0_L001_R1_001.fastq.gz"
R2: "Undetermined_S0_L001_R2_001.fastq.gz"
I1: "Undetermined_S0_L001_I1_001.fastq.gz"
I2: "Undetermined_S0_L001_I2_001.fastq.gz"
```

### Reference genome

Genome name is user-defined but must be referenced using the exact same name in the Sample Data Sheet.

'Index' corresponds to the prefix used when bowtie2 index was built ([Build index manualy]).

``` yaml
## path to references
genome:
  GRCh37:
    fasta: "/PATH_TO_REFERENCE/GRCh37/Sequences/GRCh37.primary_assembly.genome.fa"
    index: "/PATH_TO_REFERENCE/GRCh37/Indexes/Bowtie2/GRCh37.primary_assembly.genome" 
    annotation: "/PATH_TO_REFERENCE/GRCh37/Annotations/gencode.v19.annotation.gtf.gz"
  GRCh38:
    fasta: "/PATH_TO_REFERENCE/GRCh38/Sequences/GRCh38.primary_assembly.genome.fa"
    index: "/PATH_TO_REFERENCE/GRch38/Indexes/Bowtie2/GRCh38.primary_assembly.genome" 
    annotation: "/PATH_TO_REFERENCE/GRCh38/Annotations/gencode.v46.annotation.gtf.gz"
    oncogene_list: "/PATH_TO_REFERENCE/GRCh38/Annotations/OncoList_OncoKB_GRCh38_2025-07-04.tsv"
  GRCm39:
    fasta: "/PATH_TO_REFERENCE/GRCm39/Sequences/GRCm39.primary_assembly.genome.fa"
    index: "/PATH_TO_REFERENCE/GRCm39/Indexes/Bowtie2/GRCm39.primary_assembly.genome" 
    annotation: "/PATH_TO_REFERENCE/GRCm39/Annotations/gencode.vM36.annotation.gtf.gz"
    oncogene_list: ""
```

### Reads filtering

After adaptor & ODN triming and before alignment to the reference genome, reads that are too short are discarded. This filter is applied if **any** of the mates size is below this threshold.

``` yaml
################################################
minLength: 25 ## Minimal read length after trimming, before alignment
################################################
```

### Alignment to reference genome

Define which aligner to use and the range of fragment size.

``` yaml
## Alignement 
################################################
aligner: "bowtie2"   ## Aligner to use (bowtie2 only for now)
minFragLength: 100         # Minimal fragment length after alignment
maxFragLength: 1500        # Maximal fragment length after alignment 
rescue_R2: "TRUE" # keep R2 from unaligned pairs or pairs with too short R1 reads and align them as single-end reads.
################################################
```

### Insertion sites calling

After alignment on the genome, R2 reads are collapse to single insertion points and aggregated if they cluster in a distance smaller than `ISbinWindow` defined here. Filters can be applied to exclude UMI with few reads (`minReadsPerUMI`) or insertion sites with few UMIs (`minUMIPerIS`).

Here, you can also define if you want to tolerate bulges in the alignment between the gRNA and gDNA.

``` yaml
################################################
## Off targets calling
tolerate_bulges: "TRUE"           # whether to include gaps in the gRNA alignment (this will change the gap penalty during SW pairwise alignment)
max_edits_crRNA: 6                # filter clusters with less or equal than n edits in the crRNA sequence (edits = substitutions + INDELs)
ISbinWindow: 100                  # insertion sites closer than 'ISbinWindow' will be clustered together
minReadsPerUMI: 4                 # Min number of reads per UMIs (>=)
minUMIPerIS: 4                      # Min number of UMI per IS (>=)
slopSize: 50                      # window size (bp) around IS (both directions) to identify gRNA sequence (ie 50bp = -50bp to +50bp)
min_predicted_distance: 100       # distance between cut site and predicted cut site to consider as predicted
################################################
```

### UMI correction

Due to potential sequencing errors, additionnal UMIs may be detected and a correction step is required.

For each cluster of Off Targets, a similarity matrix between all UMIs detected is calculated and similar UMI collapsed together if the editing distance is smaller than `UMI_hamming_distance` . The Adjacency method described in UMI-tools (@Smith2017) is used by default (see <https://umi-tools.readthedocs.io/en/latest/the_methods.html> for details).

``` yaml
################################################
# post alignment
minMAPQ: 20                      # Min MAPQ score to keep alignments 
UMI_hamming_distance: 1         # min distance to cluster UMI using network-based deduplication, use [0] to keep raw UMIs
UMI_deduplication: "Adjacency"  # method to correct UMI (cluster or Adjacency)
UMI_pattern: "NNWNNWNN"  
UMI_filter: "TRUE"              # If TRUE, remove UMIs that do no match the expected pattern [FALSE or TRUE]
################################################
```

### Reporting

``` yaml
################################################
# reporting
max_clusters: 100                 # max number of cluster alignments to report
minUMI_alignments_figure: 1       # filter clusters with more than n UMI in the report alignment figure (set to 0 to keep all clusters -> can be slow)
################################################
```

### Prediction

``` yaml
# Prediction
################################################
SWoffFinder:
  path: "/opt/SWOffinder" ## Path to SWoffinder on your server (downloaded from https://github.com/OrensteinLab/SWOffinder)
  maxE: 6                 # Max edits allowed (integer).
  maxM: 6                 # Max mismatches allowed without bulges (integer).
  maxMB: 6                # Max mismatches allowed with bulges (integer).
  maxB: 3                 # Max bulges allowed (integer).
  window_size: 100
################################################
```

### Reads adapter & ODN trimming sequences {#reads-adapter-odn-trimming-sequences}

Indicate which sequence will be trimmed from R1 & R2 reads ends depending on the PCR orientation and ODN used.

``` yaml
################################################
# Sequences for the trimming steps

guideseq:
  positive:
    R1_trailing: "GTTTAATTGAGTTGTCATATGT"
    R2_leading: "ACATATGACAACTCAATTAAAC"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
  negative:
    R1_trailing: "ATACCGTTATTAACATATGACAACTCAA"
    R2_leading: "TTGAGTTGTCATATGTTAATAACGGTAT"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"



iguideseq:
  positive:
    R2_leading: "ACATATGACAACTCAATTAAACGCGAGC"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
    R1_trailing: "GCTCGCGTTTAATTGAGTTGTCATATGT"
  negative:
    R1_trailing: "TCGCGTATACCGTTATTAACATATGACAACTCAA"
    R2_leading: "TTGAGTTGTCATATGTTAATAACGGTATACGCGA"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
    
    
    
olitagseq:
  positive:
    R1_trailing: "GGGGTTTAATTGAGTTGTCATATGTT"
    R2_leading: "AACATATGACAACTCAATTAAACCCC"
    R2_trailing: "TCCGCTCCCTCG"
  negative:
    R1_trailing: "CCCATACCGTTATTAACATATGAC"
    R2_leading: "GTCATATGTTAATAACGGTATGGG"
    R2_trailing: "TCCGCTCCCTCG"

tagseq:
  positive:
    R1_trailing: "TGCGATAACACGCATTTCGCATAA"
    R2_leading: "CTTATGCGAAATGCGTGTTATCGCA"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
  negative:
    R1_trailing: "ATCTCTGAGCCTTATGCGAAATGC"
    R2_leading: "CGCATTTCGCATAAGGCTCAGAGAT"
    R2_trailing: "AGATCGGAAGAGCGTCGTGT"
    
################################################
```

## Folder structure

In order to start a run:

-   create a new directory in the installation folder and `cd` into.

-   Then :

    -   add the sample data sheet ([Prepare samples data sheet (SDS)](#prepare-samples-data-sheet-sds))

    -   add the configuration file ([Prepare the configuration file])

    -   add the sequencing `undeterminded_R1/R2/I1/I2.fastq.gz` files (undemultiplexed) or symblic link to the raw data.

> Input fastq files should respect the following structure from original paper :
>
> -   R1 : contains fragment sequence starting in gDNA
>
> -   R2: Starts with ODN sequence followed by gDNA sequence and potential trailing adaptor sequence
>
> -   i1 : Contains barcode 1 (usually 8 nucleotides)
>
> -   i2 : Contains barcode 2 and UMI sequence (usually 8 + 8 nucleotides)

Your folder architecture should look similar to :

``` bash
Path/to/guideseq/
├── 00-pipeline/
├── 01-envs/
├── 02-resources/
├── 03-misc/
├── test_dataset/
├── My/folder/
    ├── guideSeq_GNT.yml
    ├── sampleInfo.csv
    ├── undertermined_R1.fastq.gz
    ├── undertermined_R2.fastq.gz
    ├── undertermined_I1.fastq.gz
    └── undertermined_I2.fastq.gz
```

From inside your analysis folder, run the command below after adjusting for number of CPU (-j) :

``` bash
snakemake -s ../00-pipeline/guideseq_GNT.smk \
  -j 24 \
  -k \
  --use-conda \ 
  --report-after-run --report runtime_report.html  \
  -n
    
# remove the -n argument to start the pipeline
```

Each rule takes a maximum of 6 threads (alignment, trimming) to speed up the data processing. Set `-j` as a multiple of 6 to process 2 (12threads) , 3 (18 threads) ... samples in parallel.

``` bash
## usefull snakemake arguments 
# see https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options

--quiet
--notemp # keep all intermediate files
--rerun-trigger {code,input,mtime,params,software-env} 
--rerun-incomplete

 --filegraph | sed -n '/digraph/,$p' | dot -Tpdf > dag_files.pdf #DAG representation of workflow
 --rulegraph | sed -n '/digraph/,$p' | dot -Tpdf > dag_rules.pdf 
```

# Output

If everything goes well, the pipeline should end successfully :

``` latex
Workflow finished, no error
 _____________
< You rock baby !! >
 -------------
        \   ^__^
         \  (♥♥)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```

Otherwise, an error will be raised and the origine of the problem reported by snakemake:

``` latex
< Houston, we have a problem >
 ----------------------------
       \   \_______
 v__v   \  \   O   )
 (xx)      ||----w |
 (__)      ||     ||  \/\
```

Upon pipeline completion & if `'clean_intermediates_files: 'FALSE'`, your folder should now look like (test dataset provided in ./test) :

``` latex
my_folder_name/
├── 00-demultiplexing
├── 01-trimming
├── 02-filtering
│   ├── VEGFA_s1_K562_neg_R1.UMI.ODN.trimmed.filtered.fastq.gz
│   ├── VEGFA_s1_K562_neg_R2.UMI.ODN.trimmed.filtered.fastq.gz
│   ├── VEGFA_s1_K562_pos_R1.UMI.ODN.trimmed.filtered.fastq.gz
│   └── VEGFA_s1_K562_pos_R2.UMI.ODN.trimmed.filtered.fastq.gz
├── 03-align
│   ├── VEGFA_s1_K562_neg_R1.UMI.ODN.trimmed.unmapped.fastq.gz
│   ├── VEGFA_s1_K562_neg_R2rescued.UMI.ODN.trimmed.filtered.sorted.bam
│   ├── VEGFA_s1_K562_neg_R2rescued.UMI.ODN.trimmed.filtered.sorted.bam.bai
│   ├── VEGFA_s1_K562_neg_R2.UMI.ODN.trimmed.unmapped.fastq.gz
│   ├── VEGFA_s1_K562_neg.UMI.ODN.trimmed.filtered.sorted.filtered.bam
│   ├── VEGFA_s1_K562_neg.UMI.ODN.trimmed.filtered.sorted.filtered.bam.bai
│   ├── VEGFA_s1_K562_pos_R1.UMI.ODN.trimmed.unmapped.fastq.gz
│   ├── VEGFA_s1_K562_pos_R2rescued.UMI.ODN.trimmed.filtered.sorted.bam
│   ├── VEGFA_s1_K562_pos_R2rescued.UMI.ODN.trimmed.filtered.sorted.bam.bai
│   ├── VEGFA_s1_K562_pos_R2.UMI.ODN.trimmed.unmapped.fastq.gz
│   ├── VEGFA_s1_K562_pos.UMI.ODN.trimmed.filtered.sorted.filtered.bam
│   └── VEGFA_s1_K562_pos.UMI.ODN.trimmed.filtered.sorted.filtered.bam.bai
├── 04-IScalling
│   ├── VEGFA_s1_K562_neg.cluster_slop.bed
│   ├── VEGFA_s1_K562_neg.cluster_slop.fa
│   ├── VEGFA_s1_K562_neg_R2rescued.reads_per_UMI_per_IS.bed
│   ├── VEGFA_s1_K562_neg.reads_per_UMI_per_IS.bed
│   ├── VEGFA_s1_K562_neg.reads_per_UMI_per_IS_corrected.bed
│   ├── VEGFA_s1_K562_neg.UMIs_per_IS_in_Cluster.bed
│   ├── VEGFA_s1_K562_pos.cluster_slop.bed
│   ├── VEGFA_s1_K562_pos.cluster_slop.fa
│   ├── VEGFA_s1_K562_pos_R2rescued.reads_per_UMI_per_IS.bed
│   ├── VEGFA_s1_K562_pos.reads_per_UMI_per_IS.bed
│   ├── VEGFA_s1_K562_pos.reads_per_UMI_per_IS_corrected.bed
│   └── VEGFA_s1_K562_pos.UMIs_per_IS_in_Cluster.bed
├── 05-Report
│   ├── VEGFA_s1_K562_neg.rdata
│   ├── VEGFA_s1_K562_neg.stat
│   ├── VEGFA_s1_K562_pos.rdata
│   └── VEGFA_s1_K562_pos.stat
├── 06-offPredict
│   └── GRCh38_GGGTGGGGGGAGTTTGCTCC_NGG_3.csv
├── guideSeq_GNT.yml
├── logs
│   ├── demultiplexing_R1.log
│   ├── demultiplexing_R2.log
│   ├── VEGFA_s1_K562_neg.filter.log
│   ├── VEGFA_s1_K562_neg.odn.log
│   ├── VEGFA_s1_K562_neg_R2rescued.UMI.ODN.trimmed.filtered.align.log
│   ├── VEGFA_s1_K562_neg.trailing.log
│   ├── VEGFA_s1_K562_neg.UMI.log
│   ├── VEGFA_s1_K562_neg.UMI.ODN.trimmed.filtered.align.log
│   ├── VEGFA_s1_K562_pos.filter.log
│   ├── VEGFA_s1_K562_pos.odn.log
│   ├── VEGFA_s1_K562_pos_R2rescued.UMI.ODN.trimmed.filtered.align.log
│   ├── VEGFA_s1_K562_pos.trailing.log
│   ├── VEGFA_s1_K562_pos.UMI.log
│   └── VEGFA_s1_K562_pos.UMI.ODN.trimmed.filtered.align.log
├── results
│   ├── report-files
│   │   ├── VEGFA_s1_K562_neg_offtargets_dynamic_files
│   │   │   ├── .... truncated
│   │   ├── VEGFA_s1_K562_neg_offtargets_dynamic.html
│   │   ├── VEGFA_s1_K562_neg_offtargets.html
│   │   ├── VEGFA_s1_K562_pos_offtargets_dynamic_files
│   │   │   ├── .... truncated
│   │   ├── VEGFA_s1_K562_pos_offtargets_dynamic.html
│   │   └── VEGFA_s1_K562_pos_offtargets.html
│   ├── report.rdata
│   ├── test_report.html
│   ├── VEGFA_s1_K562_neg_summary.xlsx
│   └── VEGFA_s1_K562_pos_summary.xlsx
├── sampleInfo.csv
├── test_datasheet.csv
└── test_datasheet.xlsx
```

## Report

A general report is generated. It summarizes all main QC and key features obtained from the run using graphical representations and tables.

An example is available in the `./test_dataset/results` folder.

## Off-targets files

For each sample, an Excel file containing all the OT information is created. This file includes numerous columns, some of which are particularly noteworthy. It lists all detected positions, even those without any gRNA matches. For each position, the file reports the total number of UMIs and reads, providing the same information for both positive and negative PCRs. Insertion sites are annotated to genes and oncogenes if they are defined in the configuration file. For OT sites with gRNA matches, a summary of indel and mismatch positions is also provided.

Here is the data formatted into a markdown table:

| Index | Header value | Description |
|------------------------|------------------------|------------------------|
| 1 | Chromosome | Chromosome ID where cluster is located |
| 2 | Start_cluster | Genomic coordinates of the cluster |
| 3 | end_cluster |  |
| 4 | ClusterID | ID of cluster |
| 5 | MedianMAPQ_cluster | Median of MAPQ score for alignments in the cluster [0, 42] |
| 6 | N_IS_cluster | Number of unique cut positions in the cluster [1, +inf] |
| 7 | N_orientations_cluster | Number of ODN integration orientation in the cluster [1, 2] |
| 8 | N_orientation_PCR | Number of PCR orientation detected in the cluster [1, 2] |
| 9 | N_UMI_cluster | Sum of UMIs for all sites in the cluster [1, +inf] |
| 10 | N_reads_cluster | Sum of reads for all sites in the cluster [1, +Inf] |
| 11 - 20 |  | Same as 5-10, broken down by PCR orientation |
| 21 | Sequence_window | Nucleotide sequence of the cluster |
| 22 | grna_orientation | Strand the gRNA sequence is matched [Watson, crick] |
| 23 | Seq_gDNA | Sequence of the gDNA matching the gRNA |
| 24 | Seq_gRNA | Sequence of the gRNA |
| 25 | Alignment | Pairwise alignment representation between gRNA and gDNA |
| 26-27 | Gibbs hybridization | Efficacy and dG_0 |
| 28 | GC_content | GC% of the gDNA sequence matched |
| 29 | Alignment_score | Alignment score for pairwise alignment of gRNA & gDNA |
| 30 | Identity_pct | \% identity between gRNA & gDNA sequences |
| 31 | N_edits | Total number of Edit in the alignment (mismatches + indels) |
| 32 | mismatches_position_gRNA | Position of mismatches in the gRNA |
| 33 | \- |  |
| 34 | N_indels | Number of indels in the pairwise alignment |
| 35 | Insertion_length | Total length of insertions in the pairwise alignment |
| 36 | Start_insertion | Start, end and width of each insertion in the pairwise alignment |
| 37 | End_insertion |  |
| 38 | Width_insertion |  |
| 39 - 43 |  | Start, end and width of each deletion in the pairwise alignment |
| 44 | Pam_gDNA | PAM sequence identified in gDNA |
| 45 | Pam_gRNA | Theoretical PAM sequence (sample datasheet) |
| 46 | Pam_side | Which side is the PAM located (sample datasheet) [5, 3] |
| 47 | PAM_indel_count | Number of indels in the PAM sequence |
| 48 | PAM_indel_pos | Position of indels in the PAM sequence |
| 49 | Cut_gRNA_alignment | Theoretical cut site from pairwise alignment and offset (sample datasheet) |
| 50 | Start_IS | Genomic Position, abundance and relative abundance of the most abundant cut site in the cluster |
| 51 | Count_UMI |  |
| 52 | UMI_proportion |  |
| 53 | Cut_modal_position | Cut site with the higher UMI counts in the cluster |
| 54 | gRNA | gRNA sequence (sample datasheet) |
| 55 | gRNA_name | gRNA name (sample datasheet) |
| 56 | Gene_ensemblID | Gene ID, symbol, type and position (Intron, Exon) for intragenic cut sites |
| 57 | Symbol |  |
| 58 | Gene_type |  |
| 59 | Position |  |
| 60 | Onco_annotation | Annotation with user defined oncogene |
| 61 | Best | Is cluster the best candidate relative to number of edits |
| 62 | Relative_abundance | Relative abundance of cluster among all clusters |
| 63 | Rank | Rank of cluster (by decreasing relative abundance) |
| 64 | Predicted | Was position predicted? |
| 65 | bulge | If indels are present, do they form gRNA or gDNA bulge, or both |

# Pipeline step-by-step explanations

## Demultiplexing

**Tool** : cutadapt

Undetermined fastq files are demultiplexed to `sampleName` fastq files according to barcodes present in the sample data sheet.

-   i1 and i2 fastq files are concatenated to a single fastq file (i3)

-   R1 and R2 fastq are demultiplexed according to i1+i2 sequence present in i3 fastq files.

## ODN trimming:

**Tool** : cutadapt

The leading ODN sequence is remove from R2 reads according to the `method` defined in the sample data sheet for each sample.

Reads that **do not start** with the ODN sequence are discarded.

PCR orientation is automatically detected and added to read name for future processing:

``` yaml
@M02111:194:000000000-LT722:1:1101:15810:6788_positive
```

## UMI extraction

**Tool** : cutadapt

UMI is extracted from the I3 read generate before based on UMI pattern length defined in the configuration file.

UMI is added to reads name:

``` yaml
@M02111:194:000000000-LT722:1:1101:15810:6788_positive_GCTGTAGG
```

## Adaptor trimming:

**Tool** : cutadapt

The adaptor trailing sequences are trimmed from R1 and R2 reads respectively if present.

## Read filtering:

**Tool** : cutadapt

After trimming, only read pairs with **both** mates longer than a `minLength` defined in the configuration file are selected for alignment.

If the `rescue_R2` variable is set to `TRUE`, then R2 reads that are longer that `minLength` is a discrarded pair are rescued and processed as single-end reads.

## Genome alignment:

**Tools**: Bowtie2

Reference sequences are specified in the configuration file. Genome index is build if it does not already exist.

Trimmed reads pairs that passed the filtering steps are then aligned on the `reference genome` specified in the sample data sheet for each sample using the `aligner` specified in configuration file.

### Multi-hits management:

Multihits are reads with multiple possible equally good alignment positions in the genome. We choose to keep only a single random alignment for each read instead of reporting all possible positions.

Multihit alignments with low MAPQ score are discarded except if the Alignment score (AS tag) is equal to the score of the second best alignment score (XS tag) .

On the long run, if multi read arise from the same cuting site, they will distribute randomly to all sites (explain more).

## Cutting site calling:

**Tool** : bedtools and awk

Following alignment on the reference genome, nuclease cutting sites are extracted from the start position of R2 read alignment.

Reads that align at the same cutting site with the same UMI are aggregated together keeping track of total reads count.

## UMI correction:

**Tool** : R script

UMI sequences are corrected for potential sequencing error. UMI with less than n edits (hamming distance defined in the configuration file) are clustered together.

## Cut site clustering:

**Tool**: bedtools

Cut sites than fall in the same window of `ISbinWindow` defined in the configuration file are clustered together.

## gRNA match:

**Tool**: R script

For each cluster of cutting sites, the gRNA sequence defined in the sample data sheet for each sample is looked up in a window of +/- `slopSize` bp using the Swith-Watterman algorithm.

Gap tolerance can be accepted to detect bulges in gDNA or gRNA if the `tolerate_bulges` variable is set.

## Annotation of clusters:

**Tool** : R script

Clusters of cutting sites are annotated using the gtf file specified in the configuration file for each organism.

A first step prepare the annotation file to extract only gene and exons features.

A second step annotate clusters to gene and position relative to those genes (exon/intron). Multiple annotations may be present for each cluster and are all reported.

## Off target prediction:

**Tool** : SWOffinder

For each gRNA sequence, PAM sequence and each organism specified in the sample data sheet, a prediction of OTS is realized using the SWOffinder tool with edits and mismatches tolerance defined in the configuration file.

## Reporting:

A report is generated with main tables and graphical representations to better understand pipeline results.

# Troubleshooting

Potential issues :

Write permissions : in case of error, check that the user has read/write permission to the different reference/annotation files and working directory

# References
