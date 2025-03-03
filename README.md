# rowan-imputation-tests

Repository for manually running the test files for the initial imputation pipeline.

## Pipeline overview

### Step 1: Synthetic downsampling

Our test file contains 36 individuals and XX SNPs. Troy sent over a file with SNP50 positions in the format of CHR, SNP_NAME, BLANK, POSITION. To perform the synthetic downsampling, I'm using `bcftools view`. When subsampling, `bcftools view` expects the `--regions-file` argument to be supplied with this file, in the format of of CHROM POS (see help documentation if you'd like to include start and end positions for a different file format). So first, we need to remove the two middle columns that disrupt `bcftools view` from reading the file.

#### Modify downsampling coordinates file

We'll do so by: 

```
awk 'BEGIN{FS=OFS="\t"}{print $1,$4}' GGP_Bovine_50K_C.map > snp50_coordinates.txt
```

What this does: `BEGIN` tells `awk` to execute before reading the file. Specifies that we have a tab-delimited file `(OFS="\t")` and then asks it to print the two columns we want to keep, which were 1 and 4 (`{print $1,$4}`). Then, we tell it the input file name and specify the output file name. You can use `less snp50_coordinates.txt` to check to make sure it worked as intended. The `FS=OBS` might be redundant but I haven't tested that yet.

The next problem is that this file contains some NAs. We will need to remove any rows that contain NAs, because these aren't meaningful to `bcftools view` either. 

```
awk -F "\t" '{if($1 != "NA") { print }}' snp50_coordinates.txt > snp50_coordinates_filtered.txt
```

I then checked the number of lines between the files (`wc -l snp50_coordinates_filtered.txt'`). I checked the number of NAs in the original file (`grep 'NA' snp50_coordinates.txt | awk '{print $1}' | wc -l`) and made sure the number of lines no longer in the filtered file made sense. By comparing the before, after, and number of NAs, I can confirm only rows with NA values were removed. We can also make sure something funky didn't happen by checking `tail snp50_coordinates_filtered.txt`, where if the columns were now of uneven length, we could check.

#### Index the files

`bcftools view` expects that the files be indexed first - so we'll do `bcftools index`. In this, I just specify where I want my output file to go and give it 8 threads to work with (which is probably overkill with how fast it ran it).

```
singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools index \
--threads 8 \
--o samples/Testing_animals_seq.chr25.vcf.gz.csi \
samples/Testing_animals_seq.chr25.vcf.gz
```

#### Perform downsampling

We're going to now downsample with `bcftools view`. In this, I specify the `--regions-file`, which is our coordinates that we want to downsample to, the number of threads to use, and what format I want the output file as. In this case, `-Ob` means output the file as a gzipped BCF file. (Note: I chose this option because, although our files are small for this test run, in the future we will be working with much larger files and gzipped BCF files will be the most storage savvy). I then put the input file and specified where I wanted the output.

```
singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools view \
--regions-file samples/snp50_coordinates_filtered.txt \
--threads 8 \
-Ob \
samples/Testing_animals_seq.chr25.vcf.gz >> samples/snp50_testing_animals_seq.chr25.bcf.gz
```

Now that we've performed the subsampling, it would be good for us to just perform a sanity check and make sure we did it correctly. We're first going to print the stats for the original files, and then for the subsampled files.

The original files:

```
singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools stats \
samples/Testing_animals_seq.chr25.vcf.gz >> stats/Testing_animals_seq.chr25.stats.out
```

From this, we can see that the original file contains 36 samples with 312,811 SNPs. 

The subsampled files:

```
singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools stats \
samples/snp50_testing_animals_seq.chr25.bcf.gz >> stats/snp50_testing_animals_seq.chr25.stats.out
```

Once we've done this, we see that our new file contains 786 SNPs. This makes sense because we're only looking at CHR 25 and we've taken the subset from the SNP50 coordinates. 

### Step 2: Phase samples

Please note that our reference in this case is already phased, so it doesn't need to be phased. We need to phase our sample. Our default in the pipeline will be to phase to the reference, with the option to phase to a pedigree file. We're going to perform phasing using `SHAPEIT5`, and in this workflow, we're going to phase to the reference. In order to phase the samples, we need to index both the reference and our subset samples. 

Index the reference...

```
singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools index \
--threads 2 \
--o databases/MU_HD_only.chr25.bcf.csi \
databases/MU_HD_only.chr25.bcf
```

Index the samples...

```
singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools index \
--threads 2 \
--o samples/snp50_testing_animals_seq.chr25.bcf.gz.csi \
samples/snp50_testing_animals_seq.chr25.bcf.gz
```

Now we can phase the samples to the reference. As a side note, I think it may not be possible for `SHAPEIT5` to zip the files back when they output. Once this is tested, it may be ideal to combine this with a `gzip` procedure for downstream sample size. 

```
singularity exec https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--hb60d31d_0 SHAPEIT5_phase_common \
--input samples/snp50_testing_animals_seq.chr25.bcf.gz \
--reference databases/MU_HD_only.chr25.bcf \
--thread 24 \
--region 25 \
--output phasing/phased_snp50_testing_animals_seq.chr25.bcf
```

### Step 3: Impute samples

We're going to perform imputation with `IMPUTE5`. The imputation steps will have multiple parts. As a brief overview of what will happen in this section, we will first create imputation chunks from the target and reference panel. While this isn't exactly necessary for our test samples, we'll be working with much larger samples in the future and this will be helpful. Then, we'll convert the file into XCF file format which woks a bit faster for `IMPUTE5`. We will impute each chunk of data  and then ligate the chunks together to make one file per chromosome. Ideally, this should speed up the imputation process.

#### Chunking step

For this step, we have a couple of parts. We're going to specify the reference with `--h`, the test samples with `--g`, the chromosome with `--r`, the name of the log file with `--l`, and the output chunked coordinates with `--o`. This will output a TXT file that contains 7 columns. The first column is numbered 0-n number of chunks, the second column is the chromosome, the third column and fourth column are the chromosome and range (i.e., 25:31475-5742958), the sixth column is ???, and the seventh column is ???. 

```
imp5Chunker_v1.2.0_static \
--h databases/MU_HD_only.chr25.bcf \
--g phasing/phased_snp50_testing_animals_seq.chr25.bcf \
--r 25 \
--l stats/phased_chunking.log \
--o imputation/phased_snp50_testing_animals_seq_chunked_coords.txt
```

#### Convert reference to XCF

The `IMPUTE5` documentation recommends that you create a XCF file for each chromosome in order to perform imputation on different chunks of the chromosome. You can also optionally convert the target panel to the XCF file format, but check for different output options.

For this, we will specify `--i` for the input file in VCF/BCF format with associated index file, `--r` will be the region containing the whole imputation region (generally a chromosome number), `--o` is the output file name, `-O` is the output encoding (sh, sparse haplotypes for ref panel; bh, binary haplotypes for SNP array phased dat;, or bcf, converts XCF back to BCF). We then also have the `-m` argument, which is the sparse MAF threshold to use for the sparse haplotype reference panel file. The `IMPUTE5` documentation states: "For reference panels used with `IMPUTE5`, please always use `-m 0.03125`, representing the optimal value in SNP array imputation settings." Please note that this value for `--m` is not the default value (the default is 0.001), so it must be included to specify.

```
xcftools_static view -i databases/MU_HD_only.chr25.bcf \
-o databases/MU_HD_only.chr25_xcf.bcf \
-O sh \
-r 25 \
-T8 \
-m 0.03125
```

