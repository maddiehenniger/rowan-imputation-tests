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

Please note that our reference in this case is already phased, so it doesn't need to be phased. We need to phase our sample. Our default in the pipeline will be to phase to the reference, with the option to phase to a pedigree file. We're going to perform phasing using `SHAPEIT5`. 