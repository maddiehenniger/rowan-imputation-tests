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

For this step, we have a couple of parts. We're going to specify the reference with `--h`, the test samples with `--g`, the chromosome with `--r`, the name of the log file with `--l`, and the output chunked coordinates with `--o`. This will output a TXT file that contains 7 columns. The first column is numbered 0-n number of chunks, the second column is the chromosome, the third column is the buffered region, the fourth column is the imputation region, the fifth column is the length, the sixth column is the number of target arkers, and the seventh column is the number of reference markers.

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

For this, we will specify `--input` for the input file in VCF/BCF format with associated index file, `--region` will be the region containing the whole imputation region (generally a chromosome number), `--output` is the output file name, `--format` is the output encoding (sh, sparse haplotypes for ref panel; bh, binary haplotypes for SNP array phased dat;, or bcf, converts XCF back to BCF). We then also have the `--maf` (`-m`) argument, which is the sparse MAF threshold to use for the sparse haplotype reference panel file. The `IMPUTE5` documentation states: "For reference panels used with `IMPUTE5`, please always use `-m 0.03125`, representing the optimal value in SNP array imputation settings." Please note that this value for `--m` is not the default value (the default is 0.001), so it must be included to specify. We also specify the number of threads with `--thread` and where to output the log file with `--log`. 

```
xcftools_static view \
--input databases/MU_HD_only.chr25.bcf \
--output databases/MU_HD_only.chr25_xcf.bcf \
--format sh \
--region 25 \
--thread 8 \
--maf 0.03125 \
--log stats/xcf_log.out
```

#### Run IMPUTE5

Now we can run the actual imputation part. In order to impute using the chunked data, we're going to have to add in a few extra arguments. From the documentation: "Each chunk of imputed data is expanded by a buffer region - this helps prevent imputation quality from deteriorating near the edges of the region. Markers in the buffer region will help the inference but do not appear in the output files, and larger buffers can improve accuracy." From our chunked file we generated earlier, we have the buffered region (column 3). We also have our imputation region (column 4). IN this test run, I'm not supplying `--m`, but this would be the fine-scale recombination map for the region to be analyzed, and if not defined, a constant recombination rate is assumed.

We're going to use a while loop to loop through the lines in our coordinates file. To do this, we're going to have our first line be `while IFS= read -r line; do`. Including `IFS=` prevents field splitting on the line by line basis and then the `read -r line` tells the while loop to read line by line. Each line in our file contains the imputation and buffering regions, so this is handy. The fourth column is going to be our region, so we'll first define that we want to, by line, print column 4's value. Then, we'll also do the same to store the buffer region. I'll use the region info to name the output and log files. Then, we define the reference database (formatted as XCF) with `--h`, the input phased samples with `--g`, the stored region info from our chunked text file as our new variable, `${region}`, with `--r`, the stored buffer region from our chunked text file as our new variable, `${buffer}`, with `--buffer-region`, the output file as our stored variable, `${out_file}`, with `--o`, and our log file as our stored variable, `${log_file}`, with `--l`. We specify `done` at the end to terminate the while loop, and feed in our input chunked coordinates TXT file that we generated earlier to tell the while loop what to read. Note that you can perform multi-threading also, but the authors of `IMPUTE5` state that parallelization by chunk is more efficient. We will have to test this later to see, since our run time was ~30s with the size of the data we currently have.

```
while IFS= read -r line; do
    chr=$(echo "$line" | awk '{print $2}')
    region=$(echo "$line" | awk '{print $4}')
    buffer=$(echo "$line" | awk '{print $3}')
    count=$(echo "$line" | awk '{print $1}')
    out_file="imputation/imputed_phased_snp50_testing_animals_seq.${chr}_${count}.bcf"
    log_file="stats/imputated_chunks_testing_animals_seq.${chr}_${count}_${region}.out"
    ../../z_rowan_imp_pipeline_tests/3_imputation/impute5/impute5_v1.2.0/impute5_v1.2.0_static \
    --h databases/MU_HD_only.chr25_xcf.bcf \
    --g phasing/phased_snp50_testing_animals_seq.chr25.bcf \
    --r ${region} \
    --buffer-region ${buffer} \
    --o ${out_file} \
    --l ${log_file}
done < imputation/phased_snp50_testing_animals_seq_chunked_coords.txt
```

#### Ligation step

We now need to join the several imputed files we generated. According the the documentation: "The simplest way to ligate imputated chunks back is using `bcftools concat` providing the list of files in the right order." The thing about this bit is that it must "be in the right order". I originally had saved all the file names with the regions, but this makes generating a list "in the right order" kinda annoying. So, instead I went back and changed it to print the chromosome (because I anticipate more chromosomes in our future) and then, the first column is a count column, so it prints that also, which we can either inspect the log file or look back on our chunked coordinates TXT file if we need to reference the specific region. This essentially should print everything in order, since the chunked coordinates file is organized "in order". Note that this may be entirely messed up if there are >10 numbers, need to test later.

We can then run `bcftools concat` to put together all of these from a TXT file of the list of chunked, imputed samples. First, we'll make the list of everything in order.

```
ls imputation/*bcf >> imputation/imputed_phased_snp50_file_names.txt
```

We can look at this file to validate that it's "in order" (for example, using `less`; note: you press `q` to exit `less`).

Now to run `bcftools concat`. 

```
singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0 bcftools concat \
-n \
-f imputation/imputed_phased_snp50_file_names.txt \
-Ob \
-o imputation/ligated/ligated_imputed_phased_snp50_testing_animals_seq.chr25.bcf
```







