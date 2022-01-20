# MRSDexplorer
Toolkit for generation of Minimum Required Sequencing Depth (MRSD) calculations from RNA-seq datasets.

## Table of Contents

- [Overview](#overview)
- [Generation of MRSD input files](#inputfiles)
	- [Junction counts](#counts)
	- [Sequencing depths](#depths)
	- [Genes/Transcripts of interest](#genes)
- [Running MRSD](#run)
- [Output file formats](#output)

## <a name="overview"></a>Overview
MRSD (Minimum Required Sequencing Depth) is a statistical framework to predict the depth of RNA sequencing required to achieve a user-specified level of coverage of a transcript(s) of interest. MRSD score calculations and utility have been described in detail in our recent manuscript (https://doi.org/10.1016/j.ajhg.2021.12.014).

These MRSD scores are derived from RNA-seq data available through the Genotype-Tissue Expression project (GTEx, https://www.gtexportal.org/home/). GTEx samples are sequenced using an Illumina 75 bp paired-end read poly-A enrichment workflow. As such, the pre-computed MRSDs are most accurate for assessing transcript coverage with similar RNA-seq workflows. For alternative sequencing methodologies, parameter combinations, or tissues, users are invited and encouraged to generate their own MRSD scores using the scripts provided in the `bin/` subdirectory of the above GitHub toolkit, named **MRSDexplorer**.

## <a name="input_files"></a>Generation of MRSD input files
A minimum of three files are required for each run of MRSD calculations. 

For each set of RNA-seq data to be analysed, the three input files are: 

### <a name="counts"></a>1. Splice junction read counts
To generate transcriptome-wide junction counts, we make use of a pipeline developed by Cummings et al. (2017; https://doi.org/10.1126/scitranslmed.aal5209). The computational steps required to produce a count file from FASTQs, including split-read alignment, are described in detail at https://macarthurlab.org/2017/05/31/improving-genetic-diagnosis-in-mendelian-disease-with-transcriptome-sequencing-a-walk-through/, and the corresponding scripts are available at https://github.com/berylc/MendelianRNA-seq.

The resulting junction count files (with suffix *.csf.splicing.txt*) must be merged and converted for compatibility with MRSDexplorer, which can be done using the `process_cummings_output.py` script:

**Usage:**
```sh
python process_cummings_output.py output_list output_prefix
```
Where `output_prefix` is the prefix of the output file (and must uniquely identify this dataset from any that have been previously generated) and `output_list` is a one-per-line list of paths to the junction count files to be merged (a single-line file is also acceptable):
```sh
$> head example_files/example.junc-count-files.txt

/Users/Documents/SomeDirectory/All.sample1.csf.splicing.txt
/Users/Documents/SomeDirectory/All.sample2.csf.splicing.txt
/Users/Documents/SomeDirectory/All.sample3.csf.splicing.txt
/Users/Documents/SomeDirectory/All.sample4.csf.splicing.txt
/Users/Documents/SomeDirectory/All.sample5.csf.splicing.txt
…
```
This will generate a reformatted junction count file named *output_prefix.junc-counts.txt*. We recommend gzipping this file to preserve disk space. Junction counts generated through alternative approaches may be used for downstream analysis provided they are converted to an identical format to the *.junc-counts.txt* file, as shown here:
```sh
$> head example_files/example.junc-counts.txt

SJID	sample1	sample2	sample3	sample4	…	sample9
chr1-11995-899764	0	0	0	0	…	0
chr4-53525843-53527059	11	8	8	4	…	12
chr8-11700403-11701321	109	130	211	214	…	240
chr19-39370136-39692624	0	0	1	0	…	0
```

### <a name="depths"></a>2. Sample sequencing depths
To calculate MRSD, a text file must be provided detailing the sequencing depths of the samples listed in the *.junc-counts.txt* file (when using the STAR split-read aligner, these values can be found in the *.Log.final.out* file for the run). The file must be named *output_prefix.seq-depths.txt* - sharing the prefix with the corresponding *.junc-counts.txt* file - and must have four columns of the following format:
```sh
$> head example_files/example.seq-depths.txt

#sample_id	total_input	uniquely_mapping	multimapping
sample1	43854892	41543875	-
sample2	20394758	19343265	-
…	…	…	…
sample9	39432729	38295432	-
```
The columns must be ordered as above. However, the MRSD calculation will only use one read type (of the user’s choice), and columns deemed unnecessary for analysis may be filled with any non-numerical value; they must be filled with something to pass input file validation.

Importantly, the IDs in the sample_id column must precisely match those listed in the *output_prefix.junc-counts.txt* file. Sample IDs identified in only one of the two files will be excluded from analysis.

### <a name="genes"></a>3. Genes or transcripts of interest
MRSDexplorer allows users to predict the required sequencing depth for coverage of splice junctions in both whole genes and individual transcripts, and a file must be provided showing the desired gene(s) or transcript(s) for analysis. A single transcript model for each gene in the GENCODE human v19 annotation has been assigned according to a hierarchical approach, as detailed in our accompanying manuscript.

Where analysis of genes, rather than transcripts, is desired, users can provide a text file containing a list of genes of interest, one per line. This file can have up to two columns, containing the GENCODE v19 gene symbol and/or the Ensembl gene ID.
```sh
$> head example_files/example.genes.txt

RPL30P11    ENSG00000244573
QTRTD1 ENSG00000151576
C5orf64    ENSG00000178722
CASS4  ENSG00000087589
…	…
```
If the MRSDs of individual transcripts are desired, intronic coordinates for these transcripts can be provided in a BED format:
```sh
$> head example_files/example.introns.bed

#chr    start   end     id
chr3    99536887        99730547        ENST00000496116.1_CMSS1
chr3    99730602        99742163        ENST00000496116.1_CMSS1
chr3    99742366        99744272        ENST00000496116.1_CMSS1
chr3    99744402        99758619        ENST00000496116.1_CMSS1
chr20   57226462        57227129        ENST00000496117.1_STX16
chr20   57227143        57234679        ENST00000496117.1_STX16
…	…	…	…
```
The supplied coordinates should correspond to the first and final bases of the intron (with respect to chromosomal coordinates), and the ID column should contain a uniquely identifying ID for the transcript, which will be displayed in the output file(s). We have provided a script, `transcript_bed_generator.py`, which allows users to generate such BED files from a supplied GTF file and list of transcript IDs of interest.

**Usage:**
```sh
python transcript_bed_generator.py transcripts annotation_gtf output_prefix

Where `transcripts` is a one-per-line list of transcript IDs, `annotation_gtf` is the GTF and `output_prefix` is the prefix of the resulting file name. Running this script yields a BED file with the suffix *.introns.bed* that can be used for downstream score calculation.
```
## <a name="run"></a>Running MRSD
Prior to running the final MRSD generation step, the *.junc-counts.txt* and *.seq-depths.txt* files should be moved to the `files/` subdirectory of the MRSDexplorer directory, the default location in which they are searched for. The `MRSDexplorer.py` script can then be used to generate MRSD scores. A typical MRSDexplorer command may look something like this:
```sh
python MRSDexplorer.py --tissues all --number_reads 10 --sj_prop 0.8 --mrsd_param 0.9 --output_prefix example transcripts.introns.bed
```
The only required argument is the name of the gene list file or transcript BED, however, multiple optional parameters allow customisation of the MRSD query, including:

`--tissues = a comma-delimited list of prefixes for the desired analysis datasets. MRSDexplorer will search in the files/ subdirectory (or other directory specified by the --splice_dir option) for both a .seq-depths.txt and .junc-counts.txt file with the listed prefixes and include them in analysis if both are present (default = “all”, analysing all matching pairs of .seq-depths.txt and .junc-counts.txt files identified in the --splice_dir subdirectory).`

`--number_reads = desired number of reads to cover proportion of splice junctions specified in --sj_prop (default = 8)`

`--sj_prop = desired proportion of transcript splice junctions to be covered by read count specified in --number_reads (default = 0.75, i.e. 75%)`

`--mrsd_param = percentile MRSD value among control samples to be returned as overall MRSD value (default = 0.95, i.e. 95%)`

`--read_type = selected type of read to form the basis of the MRSD calculation (options = total, unique or multimap; default = unique, i.e. uniquely mapping)`

Once `MRSDexplorer.py` has loaded the sequencing depths and junctions, and evaluated the MRSDs for the selected genes/transcripts, one or more files will be output into the working directory.

##<a name="output"></a>Output file formats
The primary output file of `MRSDexplorer.py` is a text file with the suffix *.results.mrsd.txt*. This file begins with a set of summary statistics detailing the number of datasets and genes/transcripts successfully surveyed, followed by a list of MRSD scores (in millions of reads) for each of the successfully mapped gene/transcript identifiers, stratified by control dataset.
```sh
$> head example.results.mrsd.txt

INPUT PARAMETERS
-------------

Selected dataset(s): example
Selected read type: unique
Desired coverage of splice junctions: 8
Proportion of splice junctions to reach this coverage per gene: 0.75
Confidence level: 0.95


MINIMUM REQUIRED SEQUENCING DEPTH
---------------------

Where '-' is given as a required read depth, median coverage of splice junctions for the given gene was 0 reads.

Number of genes successfully cross-referenced: 3322

id1	id2	MRSD(example)
AAAS	ENSG00000094914	20.76
AAGAB	ENSG00000103591	72.58
AARS	ENSG00000090861	71.32
AARS2	ENSG00000124608	-
AASS	ENSG00000008311	-
ABAT	ENSG00000183044	-
ABCA1	ENSG00000165029	-
ABCA12	ENSG00000144452	-
ABCA3	ENSG00000167972	-
ABCA4	ENSG00000198691	-
ABCB11	ENSG00000073734	-
ABCB4	ENSG00000005471	-
ABCB6	ENSG00000115657	-
ABCB7	ENSG00000131269	155.42
ABCC2	ENSG00000023839	-
ABCC6	ENSG00000091262	-
ABCC8	ENSG00000006071	-
...
```
If the user supplies the `--output_junctions` option when running the `MRSDexplorer.py` script, a second output file with suffix *.junc-coverage.txt* will be generated, listing the per-million read coverage of each junction for each sample in the control datasets:
```sh
$> head example.junc-coverage.txt 

#id      gene    tissue  chr     start   end     junc_no per_M_coverage
1       TUBB    example 6       30688340        30690314        1       37.51470823370052
2       TUBB    example 6       30688340        30690314        1       40.513965057139494
3       TUBB    example 6       30688340        30690314        1       45.10835016830831
4       TUBB    example 6       30688340        30690314        1       34.01621412157552
5       TUBB    example 6       30688340        30690314        1       33.83243377436169
6       TUBB    example 6       30688340        30690314        1       42.46621481799056
7       TUBB    example 6       30688340        30690314        1       45.870132602176476
8       TUBB    example 6       30688340        30690314        1       58.94871008770753
9       TUBB    example 6       30688340        30690314        1       37.367312127719416
10      TUBB    example 6       30690422        30690695        2       33.48681324439794
11      TUBB    example 6       30690422        30690695        2       39.289858900442304
12      TUBB    example 6       30690422        30690695        2       41.95609341442062
13      TUBB    example 6       30690422        30690695        2       35.06171795846258
14      TUBB    example 6       30690422        30690695        2       29.052667332037863
15      TUBB    example 6       30690422        30690695        2       39.30192446745443
16      TUBB    example 6       30690422        30690695        2       45.90391031396011
17      TUBB    example 6       30690422        30690695        2       54.83690306072924
18      TUBB    example 6       30690422        30690695        2       34.48974133680593
19      TUBB    example 6       30690805        30691117        3       33.68425907720689
20      TUBB    example 6       30690805        30691117        3       39.763706444970246
21      TUBB    example 6       30690805        30691117        3       39.114622537676794
22      TUBB    example 6       30690805        30691117        3       34.12823238981342
23      TUBB    example 6       30690805        30691117        3       29.55142556949774
24      TUBB    example 6       30690805        30691117        3       40.16491274487337
25      TUBB    example 6       30690805        30691117        3       47.08613022638734
26      TUBB    example 6       30690805        30691117        3       54.60633444239401
27      TUBB    example 6       30690805        30691117        3       35.517445190703604
...
```
The *.junc-coverage.txt* files can be very large, particularly when evaluating MRSDs for larger datasets, and so `--output_junctions` should be used with caution.
