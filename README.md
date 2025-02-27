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

I then checked the number of lines between the files (`wc -l snp50_coordinates_filtered.txt'). I checked the number of NAs in the original file (`grep 'NA' snp50_coordinates.txt | awk '{print $1}' | wc -l`) and made sure the number of lines no longer in the filtered file made sense. By comparing the before, after, and number of NAs, I can confirm only rows with NA values were removed. We can also make sure something funky didn't happen by checking `tail snp50_coordinates_filtered.txt`, where if the columns were now of uneven length, we could check.