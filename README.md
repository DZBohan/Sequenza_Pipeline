# Sequenza CNV Pipeline

#### Bohan Zhang September 13th

## Index

- [1. Introduction of the Pipeline](#1)

	- [1.1 Preparation step](#1.1)

	- [1.2 CNV analysis step](#1.2)

	- [1.3 Gene level CNV analysis](#1.3)

- [2. Files Preparation](#2)

	- [2.1 Project level files](#2.1)

	- [2.2 Sample level files](#2.2)

	- [2.3 Geneinfo.txt](#2.3)

	- [2.4 BamIdsUniq.txt](#2.4)

- [3. Environment Setting](#3)

	- [3.1 Installation of sequenza-utils](#3.1)

	- [3.2 R environment setting](#3.2)

	- [3.3 HPC modules loading](#3.3)

- [4. Introduction of Scripts](#4)

	- [4.1 Preparation step](#4.1)

	- [4.2 Main step](#4.2)

- [5. Introduction of Results](#5)

	- [5.1 CNV plots](#5.1)

	- [5.2 Segment level CNV table](#5.2)

	- [5.3 Gene level CNV table](#5.3)

## <h2 id="1">1. Introduction of the Pipeline</h2>

The pipeline uses three tools and two R packages.

* sequenza-utils (need installation)
* samtools (HPC module)
* tabix (HPC module)
* sequenza (R package)
* CNTools (R package)

Here are the official web pages introducing the [Sequenza CNV analysis steps](https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html#plots-and-results) and the [CNTools package](https://bioconductor.org/packages/release/bioc/html/CNTools.html).

The pipeline is theoretically divided into three steps, the preparation step, the CNV analysis step, and the gene level CNV analysis step, but practically, the pipeline includes two scripts, a project level preparation script and a main script containing both the CNV analysis step, and the gene level CNV analysis step.

Now, I am going to have a brief introduction to the three steps first.

### <h2 id="1.1">1.1 Preparation step</h2>

Here is a flowchart of the preparation step. The wave shape stands for files, and the rectangle stands for tools.

<div align=center><img width =20% src ="https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/sequenza_prepare.png?raw=true"/></div>

This step aims to convert the reference file into a GC wiggle track file, which is needed in the following step, by using a tool called `sequenza-utils`.

### <h2 id="1.2">1.2 CNV analysis step</h2>

Here is also a flowchart of the CNV analysis step.

![Sequenza](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/sequenza.png?raw=true)

For the `sequenza-utils`, if your samples are cram files, you must first convert them into bam files using `samtools`. However, this pipeline can recognize cram files and convert them into bams automaticly, so you do not need to do this step yourself if your files are crams.

Thus, bam files and cram files are the same in this pipeline, so I will always say bam files to stand for both bam and cram files in the following documentation.

Before the `Sequenza` analysis, you need to process your bams and wiggle track file to generate a .seqz file, and then post-process the original .seqz file by binning it.

Sequenza analysis is processed in R. The results tables and plots will be introduced in chapter 6.

### <h2 id="1.3">1.3 Gene level CNV analysis</h2>

Sequenza provides the CNV analysis on the segment level. However, if you want to know the CNV of a specific gene, you need to introduce another tool called `CNTools`. With the help of genomic information, CNTools can convert the output table of Sequenza pipeline from the segment level to the gene level.

## <h2 id="2">2. Files Preparation</h2>

### <h2 id="2.1">2.1 Project level files</h2>

These six files are the same among different samples in one project.

* reference.fa (exist)
* reference.fai (along with reference.fa)
* bamIdsUniq.txt (need to be generated yourself)
* geneinfo.txt (can be downloaded)

### <h2 id="2.2">2.2 Sample level files</h2>

These four files are sample-specific.

* tumor.cram or bam (exist)
* tumor.crai or bai (along with tumor.carm or bam)
* normal.cram or bam (exist)
* normal.crai or bai (along with normal.carm or bam)

### <h2 id="2.3">2.3 Geneinfo.txt</h2>

Geneinfo.txt is the table of information on each gene. It includes genes' locations, names, and ids. Usually, we use GRCh38 as the geneinfo file, and you can directly download the GRCh38 geneinfo file from my GitHub.

![Geneinfo ](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/geneinfo.png?raw=true)

Sometimes the chromosome format ("chrN" or "N") of the geneinfo is different from the bam files, but you do not need to modify the geneinfo file since I have a chromosome format judgment step in the CNTools R script.

If you want to use your geneinfo file, you need to remove the duplicate data and data, not on chromosomes 1-22 or X and Y.

### <h2 id="2.4">2.4 BamIdsUniq.txt</h2>

Using this pipeline, you are required to create a bamIdsUniq.txt file. The file should have two columns. The first column is all the tumor bam file names, and the second column is all the normal bam file names. Use commas to separate two columns. Here is an example of what the file looks like.

![BamIdsUniq](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/bamidsuniq.png?raw=true)

## <h2 id="3">3. Environment Setting</h2>

### <h2 id="3.1">3.1 Installation of sequenza-utils</h2>

Use this command to install `sequenza-utils` tool. You also need Python for the installation.

```
pip install sequenza-utils
```

### <h2 id="3.2">3.2 R environment setting</h2>

In this pipeline, I used R with version 4.1.2. You must install two packages, Sequenza and CNTools, in the r/4.1.2.

```
setRepositories(graphics = FALSE, ind = 1:6)
install.packages("sequenza")
```

```
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("CNTools")
```

### <h2 id="3.3">3.3 HPC module loading</h2>

In this pipeline, I load seven modules on the HPC. They are available on the USC CARC, so if you are using another HPC, please check whether these modules are available.

* module load gcc/9.2.0
* module load intel/19.0.4
* module load samtools/18.0.4
* module load htslib/1.10.2 (for tabix)
* module load gcc/11.2.0
* module load openblas/0.3.18
* module load r/4.1.2

You don't need to load the modules yourself since I have written them inside the scripts. If some of the modules are unavailable on the HPC you are using, try installing the programs yourself.

## <h2 id="4">4. Introduction of Scripts</h2>

This Sequenza CNV analysis pipeline includes two Slurm scripts, sequenza_prepare.slurm and sequenza.slurm. Both scripts have config files, config_sequenza.txt and config_sequenza_prepare.txt, for inputting the variables or files. In addition, two R scripts, sequenza.R and cntools_gatk.R, should be used in the second Slurm script. You can download these six files on GitHub. You should locate these files in the same directory when using the pipeline.

### <h2 id="4.1">4.1 Preparation step</h2>

The first step also called the preparation step, is for the whole set of samples since they share the same reference.

This step aims to generate an intermedia file, `filename.gc50Base.wig.gz`, for the next step.

For inputting the variables and files, you do not need to modify the script, `sequenza_prepare.slurm`, but write them into the config file `config_sequenza_prepare.txt`. Now, let us have a look at the config file.

![Config_pre](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/config_pre.png?raw=true)

There are two variables in this config file, `output_path` and `refer`. `output_path` is the path of a directory to store the output file `filename.gc50Base.wig.gz`. `refer` is the absolute path of the .fa file. One important thing is that you should put the .fai in the same directory as the .fa file.

After completing the config file, you can submit the Slurm script, `sequenza_prepare.slurm`, to run the preparation step.

```
sbatch --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16GB --time=10:00:00 sequenza_prepare.slurm
```

You can modify the memory, time, and other parameters, but they are usually enough in this step.

Now, you should have a new file, `filename.gc50Base.wig.gz` in the target directory and be able to go to the next step.

### <h2 id="4.2">4.2 Main step</h2>

You can do the main step with the newly generated intermedia file. This step contains a Slurm script `sequenza.slurm`, a config file `config_sequenza.txt`, and two R scripts `sequenza.R` and `cntools_sequenza.R`.

For inputting the variables and files, you need to add them into the config file `config_sequenza.txt`. Now, let us have a look at this one.

![Config](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/config.png?raw=true)

`BAMDIR` is the absolute path of the directory storing bam files. You need to put the set of bam files (bai files as well) in one directory.

`FILE` is the absolute path of the bamIdsUniq.txt file you create in chapter 2.4.

`refer` is the absolute path of the .fa file. One important thing is that you should put the .fai in the same directory as the .fa file.

`refer_wig` is the absolute path of the intermedia file generated in chapter 4.1.

`output_path` is the absolute path of the directory for storing the Sequenza pipeline outputs. Finally, there should be a set of directories named by the tumor bam file names in this directory. Inside individual directories, you can find the output files of each sample.

`geneinfo_path` is the absolute path of the geneinfo file mentioned in chapter 2.3. You can download the GRCh38 version of the geneinfo file through this GitHub.

After completing the config file, you can submit the Slurm script `sequenza.slurm`, to run the main step.

```
sbatch --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16GB --time=12:00:00 --array=[1-2]%2 sequenza.slurm
```

The parameter of the array depends on the number of samples you want to do. In other word, the total number of the array is supposed to be the same as the number of rows in the `bamIdsUniq.txt` files. The details of `bamIdsUniq.txt` can be found in chapter 2.4. For example, if you have 50 samples requiring CNV analysis, and you want to run five jobs each time, you are supposed to set the array parameter as `--array=[1-50]%5`.

It usually takes 6-8 hours to run the pipeline once, meaning if you set the array as `--array=[1-50]%5`, it will take 60-80 hours to run all 50 samples. However, for some larger bams (maybe larger than 50G), you need to modify the command on the parameter time to a longer one. I recommend setting 24 hours when the bams are large in case the process shuts down.

After running the main step, you will get a set of directories, which are the same number as the samples, inside the output directory you set, and the name of each directory is supposed to be the same as the tumor bam's filename.

![Output1](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/output1.png?raw=true)

Inside each directory, there should be 12 files (eight `.pdf` files, five `.txt` files, and two .RData files) and a directory called cntools_result.

![Output2](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/output2.png?raw=true)

This main step should have generated several intermedia files, but you cannot see them in the final results since I added a clean-up step in the script. Here I also list all intermedia files generated by this step in case you are interested in them.

```
tumor_bam_name_t.bam
normal_bam_name_g.bam
tumor_bam_name.seqz.gz
tumor_bam_name.seqz.gz.tbi
tumor_bam_name_small.seqz.gz
tumor_bam_name_small.seqz.gz.tbi
```

## <h2 id="5">5. Introduction of Results</h2>

The results of the Sequenza CNV Pipeline can be divided into three categories, CNV plots, segment level CNV table, and gene level CNV table. First, let us have a look at the CNV plots.

### <h2 id="5.1">5.1 CNV plots</h2>

There are eight plot files in Sequenza CNV Pipeline's outputs. Here I will introduce two of them. The first one is `tumor_bam_name_CP_contours.pdf`.

![CP_contours](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/CP_contours.png?raw=true)

Here is an example of CP_contours plot. The y-axis `Cellularity` indicates the percentage of tumor cells to all cells (tumor + normal) in the tumor sample. The x-axis `Ploidy` is the average ploidy of the cells in the tumor sample. `LPP` stands for Log-posterior probability of the segment. The three points on the plots stand for three possibilities of the sample, and the algorithm selects the first point (Cellularity=0.9, Ploidy-2) as the correct one. This plot is the ploidy and cellularity foundation of the next one.

However, the algorithm cannot always choose the correct one. Here is an example of an algorithm choosing the wrong point.

![CP_contours_W](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/CP_contours_w.png?raw=true)

In this case, the algorithm selected the point with the ploidy of 4, which is supposed to be 2. Thus, this plot can help us know if the algorithm selects the right situation. In addition, knowing the percentage of tumor and normal cells is also essential.

The second polt file is `tumor_bam_name_genome_view`. This one has three plots inside.

![Gene_view1](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/gene_view1.png?raw=true)

This plot individually tells the copy number of both alleles (blue and red).

![Gene_view2](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/gene_view2.png?raw=true)

This one tells the copy number of the two alleles.

![Gene_view3](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/gene_view3.png?raw=true)

The last plot is the one we are supposed to focus on. There are two subplots in this one. 

The plot above tells the allele frequency or allele fraction. Allele fraction can be used to infer whether a mutation at a specific locus is somatic or germline. It is calculated as the proportion of sequence reads mutated at a locus divided by the total coverage of the locus.

An allele fraction of roughly 50% for diploid means that the mutation in this sample at the particular locus is a heterozygous germline mutation, i.e., recombination from the parents' germ cells.

An allele fraction of 100% or 0 means that the allele is homozygous at the locus, but this cannot determine whether it is a somatic or germline mutation. The allele fraction of the reference genome is also 0. If the allele fraction is not 50%, 100%, or 0, then there are potential errors other than germline mutations, i.e., somatic mutations, at that locus for that sample.

Allele fraction's representation of somatic variation requires consideration of three additional issues. First, tumor samples often contain contamination by normal cells with only germline mutations, so the correspondence between allele fraction and somatic mutations must also include contamination by the normal cells in the calculation.

Second, changes in allele copy number or chromosome ploidy can intrinsically affect the allele fraction of heterozygous mutations. However, the allele fraction of homozygotes is independent of copy number or ploidy and is always 0 or 100%.

Third, the allele fraction for germline mutations is not 50% because some significant germline mutations, such as large insertions or deletions (INDEL), can alter allele copy numbers, even the ploidy. Therefore, the allele fraction should be used to determine tumors' somatic mutations and match normal samples.

The plot below tells the depth ratio or copy ratio, meaning the copy number of the tumor sample divided by the copy number of the normal sample. For example, 3/2 means a single copy gain, and 1/2 means a single copy gain. 

The copy ratio presentation is better than copy number since the algorithm sometimes would pick the wrong ploidy.

### <h2 id="5.2">5.2 Segment level CNV table</h2>

The output table we need to focus on is `tumor_bam_name_segemts.txt`. Here is an example of this table.

![Segments](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/segments.png?raw=true)

Here is the explanation of each column. 

![Explain](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/explain.png?raw=true)

As I mentioned in chapter 5.1, Bf (minor allel fraction) and depth.ratio (copy ratio) is supposed to be the two features we need to focus on. Along with the position information (chromosome, start.pos, and end.pos), you can know the segment level CNV of a sample.

### <h2 id="5.3">5.3 Gene level CNV table</h2>

After running the Sequenza CNV pipeline, there is supposed to be a directory called cntools_result in the output directory, and you can find a file called `tumor_bam_filename.cntools.sequenza.txt`.

![Gene_level](https://github.com/DZBohan/Sequenza_Pipeline/blob/main/images/genn_level.png?raw=true)

In this file, you can get the genenames, geneids, genes' position, and genes' log2 copy ratio.

