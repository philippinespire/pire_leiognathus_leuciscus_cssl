# Lle Data Processing Log

Log to track progress through capture bioinformatics pipeline for the Albatross and Contemporary *Leiognathus leuciscus* samples from Hamilo Cove.
  * **NOTE:** *Leiognathus leuciscus* & *Equulites leuciscus* are the same species. Furthermore, subsequent species identification work has revealed that many of these individuals are actually *Equulites laterofenestra* and NOT *Leiognathus leuciscus*.

---

## 0. Rename files for dDocent HPC

Raw data in `/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/raw_fq_capture` (check *Leiognathus leuciscus* channel on Slack). Starting analyses in `/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus`.

Used decode file from Sharon Magnuson & Chris Bird.

```bash
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/raw_fq_capture

salloc
bash

#check got back sequencing data for all individuals in decode file
ls | wc -l #256 files (2 additional files for README & decode.tsv = XX/2 = XX individuals (R&F)
wc -l Lle_CaptureLibraries_SequenceNameDecode_fixed.tsv #129 lines (1 additional line for header = XX individuals), checks out

#run renameFQGZ.bash first to make sure new names make sense
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ.bash Lle_CaptureLibraries_SequenceNameDecode_fixed.tsv

#run renameFQGZ.bash again to actually rename files
#need to say "yes" 2X
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ.bash Lle_CaptureLibraries_SequenceNameDecode_fixed.tsv rename
```

---

## 1.  Check data quality with fastqc

Ran [`Multi_FASTQC.sh`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/Multi_FASTQC.sh).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/raw_fq_capture

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/raw_fq_capture" "fq.gz"
```

Potential issues:  
  * % duplication - 
    * Alb: 79.20%, Contemp: 54.65%
  * GC content - 
    * Alb: 46%, Contemp: 47%
  * number of reads - 
    * Alb: ~19 mil, Contemp: ~6 mil

---

## 2. 1st fastp

Ran [`runFASTP_1st_trim.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runFASTP_1st_trim.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFASTP_1st_trim.sbatch raw_fq_capture fq_fp1
```

Potential issues:  
  * % duplication - 
    * Alb: 73.35%, Contemp: 49.06%
  * GC content -
    * Alb: 43.88%, Contemp: 46.16%
  * passing filter - 
    * Alb: 96.52%, Contemp: 95.92%
  * % adapter - 
    * Alb: 85.02%, Contemp: 49.25%
  * number of reads - 
    * Alb: ~31 mil, Contemp: ~8 mil

---

## 3. Clumpify

Ran [`runCLUMPIFY_r1r2_array.bash`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runCLUMPIFY_r1r2_array.bash).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runCLUMPIFY_r1r2_array.bash fq_fp1 fq_fp1_clmp /scratch/mmalabag 10
```

Ran [`checkClumpify_EG.R`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/checkClumpify_EG.R) to see if any failed.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

salloc
module load container_env mapdamage2

#had to install tidyverse package first
crun R
install.packages("tidyverse") #said yes when prompted, when finished, exited & didn't save env

crun R < /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/checkClumpify_EG.R --no-save
#all files ran successfully
```

Ran [`runMULTIQC.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runMULTIQC.sbatch)  to get the MultiQC ouput

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "fq_fp1_clmp" "fqc_clmp_report"  "fq.gz"
```

---

## 4. 2nd fastp

Ran [`runFASTP_2_cssl.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runFASTP_2_cssl.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFASTP_2_cssl.sbatch fq_fp1_clmp fq_fp1_clmp_fp2
```

Potential issues:  
  * % duplication - 
    * Alb: 24.17%, Contemp: 12.42%
  * GC content - 
    *  Alb: 45%, Contemp: 46.51%
  * passing filter - 
    * Alb: 97.88%, Contemp: 99.00%
  * % adapter - 
    * Alb: 2.74%, Contemp: 0.88%
  * number of reads - 
    * Alb: ~7 mil, Contemp: ~4.5 mil

---

## 5. Run fastq_screen

Ran [`runFQSCRN_6.bash`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runFQSCRN_6.bash).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 20
```
Some files didn't run (either due to storage problems or general errors) so I had to rerun them individually.

```sh
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_089.clmp.fp2_r2.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_092.clmp.fp2_r1.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_092.clmp.fp2_r2.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_093.clmp.fp2_r2.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_094.clmp.fp2_r1.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_094.clmp.fp2_r2.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_095.clmp.fp2_r2.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_096.clmp.fp2_r1.fq.gz
bash ../pire_fq_gz_procesing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Lle-CNas_096.clmp.fp2_r2.fq.gz
```

Checked that all files were successfully completed.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

#checked that all 5 output files from fastqc screen were created for each file (should be 256 for each = 128 R1 & 128 R2)
ls fq_fp1_clmp_fp2_fqscrn/*tagged.fastq.gz | wc -l #256
ls fq_fp1_clmp_fp2_fqscrn/*tagged_filter.fastq.gz | wc -l #256 
ls fq_fp1_clmp_fp2_fqscrn/*screen.txt | wc -l #256
ls fq_fp1_clmp_fp2_fqscrn/*screen.png | wc -l #256
ls fq_fp1_clmp_fp2_fqscrn/*screen.html | wc -l #256

#checked all out files for any errors
grep 'error' slurm-fqscrn.*out #nothing
grep 'No reads in' slurm-fqscrn.*out #nothing
```

Everything looks good, no errors/missing files.

Ran [`runMultiQC.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runMULTIQC.sbatch) separately.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runMULTIQC.sbatch fq_fp1_clmp_fp2_fqscrn fastqc_screen_report
```

Potential issues:

  * one hit, one genome, no ID - 
    * Alb: 85%, Contemp: 90%
  * no one hit, one genome to any potential contaminators (bacteria, virus, human, etc) - 
    * Alb: 4%, Contemp: 2%

---

## 6. Re-pair fastq_screen paired end files

Ran [`runREPAIR.sbatch`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/runREPAIR.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_procesing/runREPAIR.sbatch fq_fp1_clmp_fp2_fqscrn fq_fp1_clmp_fp2_fqscrn_repaired 40
```

Once finished, ran [`Multi_FASTQC.sh`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/Multi_FASTQC.sh) to assess quality.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh ".../leiognathus_leuciscus/fq_fp1_clmp_fp2_fqscrn_repaired" "fqc_rprd_report" "fq.gz"

# check to be sure the job is running
watch squeue -u mmalabag
```

Potential issues:  
  * % duplication - 
    * Alb: 22.33%, Contemp: 13.17%
  * GC content - 
    * Alb: 44%, Contemp: 45%
  * number of reads - 
    * Alb: ~3.5 mil, Contemp: ~2 mil

---

## 7. Calculate the percent of reads lost in each step

Executed [`read_calculator_cssl.sh`](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/read_calculator_cssl.sh).

```
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/read_calculator.sh "/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus" "raw_fq_capture"
```

Generated the [percent_read_loss](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/leiognathus_leuciscus/preprocess_read_change/readLoss_table.tsv) and [percent_reads_remaining](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/leiognathus_leuciscus/preprocess_read_change/readsRemaining_table.tsv) tables.

Reads lost:

  * fastp1 dropped 3.93% of the reads
  * 58.46% of reads were duplicates and were dropped by Clumpify
  * fastp2 dropped 1.28% of the reads after deduplication

Reads remaining:

  * Total reads remaining: 34.56%

---

## 8. Set up mapping dir and get reference genome

Made mapping directory and moved `*fq.gz` files over.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus
mkdir mkBAM

mv fq_fp1_clmp_fp2_fqscrn_repaired/*fq.gz mkBAM
```

Pulled latest changes from dDocentHPC repo & copied `config.6.cssl` over.

```sh
cd /pire_cssl_data_procesing/scripts/dDocentHPC
git pull

cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkBAM

cp ../../scripts/dDocentHPC/configs/config.6.cssl .
```

Found the best genome by running `wrangleData.R`, sorted tibble by busco single copy complete, quast n50, and filtered by species in Rstudio. The best genome to map to for *Leiognathus leuciscus* is: `<scaffolds.fasta>` in `</home/e1garcia/shotgun_PIRE/lle_spades/scaffolds.faasta`. Copied this to `mkBAM`. 
  * The best genome had to to first be extracted from `out_Lle-C_3NR_R1R2ORPH_contam_noisolate_covcutoff-off.tar.gz` and was placed in main `lle_spades` directory for easy access.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkBAM

cp /home/e1garcia/shotgun_PIRE/lle_spades/scaffolds.fasta .

#the destination reference fasta should be named as follows: reference.<assembly type>.<unique assembly info>.fasta
#<assembly type> is `ssl` for denovo assembled shotgun library or `rad` for denovo assembled rad library
#this naming is a little messy, but it makes the ref 100% tracable back to the source
#it is critical not to use `_` in name of reference for compatibility with ddocent and freebayes

mv scaffolds.fasta ./reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.fasta
```

Updated the config file with the ref genome info.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkBAM

nano config.6.cssl
```

Inserted `<assembly type>` into the `Cutoff1` variable and `<unique assembly info>` into the `Cutoff2` variable.

```
----------mkREF: Settings for de novo assembly of the reference genome--------------------------------------------
PE                                     Type of reads for assembly (PE, SE, OL, RPE)         PE=ddRAD & ezRAD pairedend, non-overlapping reads; SE=singleend reads; OL=ddRAD & ezRAD overlapping reads, miseq; RPE=oregonRAD, restriction site + random shear
ssl                                    Cutoff1 (integer)                                    Use unique reads that have at least this much coverage for making the reference genome
Lle-C-3NR-R1R2ORPH-contam-noIsolate    Cutoff2 (integer)                                    Use unique reads that occur in at least this many individuals for making the reference genome
0.05                                   rainbow merge -r <percentile> (decimal 0-1)          Percentile-based minimum number of seqs to assemble in a precluster
0.95                                   rainbow merge -R <percentile> (decimal 0-1)          Percentile-based maximum number of seqs to assemble in a precluster
------------------------------------------------------------------------------------------------------------------
```

---

## 9. Map reads to reference

Ran [`dDocentHPC.sbatch`](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/dDocentHPC.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkBAM

sbatch ../../../dDocentHPC/dDocentHPC.sbatch mkBAM config.6.cssl
```

---

## 10. Filter BAM files

Ran [`dDocentHPC.sbatch`](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/dDocentHPC.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkBAM

sbatch ../../../dDocentHPC/dDocentHPC.sbatch fltrBAM config.6.cssl
```

---

## 11. Generate mapping stats for capture targets

Ran `getBAITcvg.sbatch` to calculate the breadth and depth of coverage for the targeted bait regions:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch ../scripts/getBAITcvg.sbatch ./mkBAM /home/e1garcia/shotgun_PIRE/pire_probe_sets/07_Leiognathus_leuciscus/Leiognathus_leuciscus_Chosen_baits.singleLine.bed
```

Ran `mappedReadStats.sbatch` to calculate the number of reads in each filtered `.bam` file, along with their mean length and depth:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

sbatch ../../../pire_fq_gz_processing/mappedReadStats.sbatch . coverageMappedReads
```

---

## 12. Run mapDamage

Ran mapDamage:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkBAM

sbatch ../../scripts/runMAPDMG.2.sbatch "Lle-*RG.bam" reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.fasta
```

Cleaning up directories:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkBAM

mkdir mapDamage_output
mv results*-RG/ mapDamage_output

cd ..
mkdir mapdamageBAM

cd mapDamageBAM
mv ../mkBAM/mapDamage_output/results*/*bam .
```

Renamed the rescaled .bam files so that dDocent will recognize them (made them end in *-rescaled-RG.bam).

---

## 13. Run mkVCF on BAM files

Copied and renamed reference faasta and `config.6.cssl` to `mapDamageBAM`:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mapDamageBAM

cp ../mkBAM/reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.fasta ./reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.fasta
cp ../mkBAM/config.6.cssl ./config.6.cssl.mkVCF.rescale
```

Edited `config.6.cssl.mkVCF.rescale` so that the file names match and the settings are as desired:

```
----------mkREF: Settings for de novo assembly of the reference genome--------------------------------------------
PE                                     Type of reads for assembly (PE, SE, OL, RPE)         PE=ddRAD & ezRAD pairedend, non-overlapping reads; SE=singleend reads; OL=ddRAD & ezRAD overlapping reads, miseq; RPE=oregonRAD, restriction site + random shear
ssl                                    Cutoff1 (integer)                                    Use unique reads that have at least this much coverage for making the reference genome
Lle-C-3NR-R1R2ORPH-contam-noIsolate    Cutoff2 (integer)                                    Use unique reads that occur in at least this many individuals for making the reference genome
0.05                                   rainbow merge -r <percentile> (decimal 0-1)          Percentile-based minimum number of seqs to assemble in a precluster
0.95                                   rainbow merge -R <percentile> (decimal 0-1)          Percentile-based maximum number of seqs to assemble in a precluster
------------------------------------------------------------------------------------------------------------------
```

Ran mkVCF:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mapDamageBAM

sbatch ../../../dDocentHPC/dDocentHPC.sbatch mkVCF config.6.cssl.mkVCF.rescale
```

---

## 14. Filter VCF File

Made a filtering directory and copied config file over:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus
mkdir filterVCF
cd filterVCF

cp ../../scripts/fltrVCF/config_files/config.fltr.ind.cssl .
```

Updated config file with correct paths:

```
fltrVCF Settings, run fltrVCF -h for description of settings
        # Paths assume you are in `filterVCF dir` when running fltrVCF, change as necessary
        fltrVCF -f 01 02 03 04 14 07 05 16 15 06 11 09 10 04 13 05 16 07                                   # order to run filters in
        fltrVCF -c ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled                                        # cutoffs, ie ref description
        fltrVCF -b ../mapDamageBAM                                                                         # path to *.bam files
        fltrVCF -R ../../scripts/fltrVCF/scripts                                                           # path to fltrVCF R scripts
        fltrVCF -d ../mapDamageBAM/mapped.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.bed             # bed file used in genotyping
        fltrVCF -v ../mapDamageBAM/TotalRawSNPs.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.vcf       # vcf file to filter
        fltrVCF -g ../mapDamageBAM/reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.fasta        # reference genome
        fltrVCF -p ../mapDamageBAM/popmap.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled                 # popmap file
        fltrVCF -w ../../scripts/fltrVCF/filter_hwe_by_pop_HPC.pl                                          # path to HWE filter script
        fltrVCF -r ../../scripts/rad_haplotyper/rad_haplotyper.pl                                          # path to rad_haplotyper script
        fltrVCF -o Lle.A                                                                                   # prefix on output files, use to track settings
        fltrVCF -t 40                                                                                      # number of threads [1]

        01 vcftools --min-alleles       2               #Remove sites with less alleles [2]
        01 vcftools --max-alleles       2               #Remove sites with more alleles [2]
        02 vcftools --remove-indels                     #Remove sites with indels.  Not adjustable
        03 vcftools --minQ              100             #Remove sites with lower QUAL [20]
        04 vcftools --min-meanDP        5:15            #Remove sites with lower mean depth [15]
        05 vcftools --max-missing       0.55:0.6        #Remove sites with at least 1 - value missing data (1 = no missing data) [0.5]
        06 vcffilter AB min             0.375           #Remove sites with equal or lower allele balance [0.2]
        06 vcffilter AB max             0.625           #Remove sites with equal or lower allele balance [0.8]
        06 vcffilter AB nohet           0               #Keep sites with AB=0. Not adjustable
        07 vcffilter AC min             0               #Remove sites with equal or lower MINOR allele count [3]
        09 vcffilter MQM/MQMR min       0.25            #Remove sites where the difference in the ratio of mean mapping quality between REF and ALT alleles is greater than this proportion from 1. Ex: 0 means the mapping quality must be equal between REF and ALTERNATE. Smaller numbers are more stringent. Keep sites where the following is true: 1-X < MQM/MQMR < 1/(1-X) [0.1]
        10 vcffilter PAIRED                             #Remove sites where one of the alleles is only supported by reads that are not properly paired (see SAM format specification). Not adjustable
        11 vcffilter QUAL/DP min        0.2             #Remove sites where the ratio of QUAL to DP is deemed to be too low. [0.25]
        13 vcftools --max-meanDP        400             #Remove sites with higher mean depth [250]
        14 vcftools --minDP             5               #Code genotypes with lesser depth of coverage as NA [5]
        15 vcftools --maf               0               #Remove sites with lesser minor allele frequency.  Adjust based upon sample size. [0.005]
        15 vcftools --max-maf           1               #Remove sites with greater minor allele frequency.  Adjust based upon sample size. [0.995]
        16 vcftools --missing-indv      0.6:0.5         #Remove individuals with more missing data. [0.5]
```

Did not adjust the filter settings (left them as the default).

Ran [`fltrVCF.sbatch`](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/fltrVCF.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/filterVCF

sbatch ../../scripts/fltrVCF.sbatch config.fltr.ind.cssl
```

---

## 15. Check for cryptic species

Made a `pop_structure` directory and copied filtered VCF file there.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

mkdir pop_structure
cd pop_structure

cp ../filterVCF/Lle.A.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.Fltr07.18.vcf .
```

There were too many "_" in sample ID names for PLINK. Fixed with the following code:

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/pop_structure

module load bcftools
bash
export SINGULARITY_BIND=/home/e1garcia

#created list of sample names
crun bcftools query -l Lle.A.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.Fltr07.18.vcf > sample_names.txt

#modified sample_names.txt so that the only _ in the sample name was between the population and the individual number (ex: Lle-AHam_001)

#renamed
crun bcftools reheader --samples sample_names.txt -o Lle.A.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.Fltr07.18.vcf

exit
```

Ran PCA w/PLINK.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/pop_structure

module load container_env python3
bash
export SINGULARITY_BIND=/home/e1garcia

#VCF file has split chromosome, so running PCA from bed file
crun.python3 -p ~/.conda/envs/popgen plink --vcf Lle.A.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.Fltr07.18.vcf
 --alow-extra-chr --make-bed --out PIRE.Lle.Ham.rescaled.preHWE
crun.python3 -p ~/.conda/envs/popgen plink --pca --allow-extra-chr --bfile PIRE.Lle.Ham.rescaled.preHWE --out PIRE.Lle.Ham.rescaled.preHWE

exit
```

Made input files for ADMIXTURE with PLINK.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/pop_structure

bash
export SINGULARITY_BIND=/home/e1garcia

#bed and bim files already made (for PCA)
awk '{$1=0;print $0}' PIRE.Lle.Ham.rescaled.preHWE.bim > PIRE.Lle.Ham.rescaled.preHWE.bim.tmp
mv PIRE.Lle.Ham.rescaled.preHWE.bim.tmp PIRE.Lle.Ham.rescaled.preHWE.bim

exit
```

Ran ADMIXTURE (K = 1-5).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/pop_structure

module load container_env python3
bash
export SINGULARITY_BIND=/home/e1garcia

crun.python3 -p ~/.conda/envs/popgen admixture PIRE.Lle.Ham.rescaled.preHWE.bed 1 --cv > PIRE.Lle.Ham.rescaled.preHWE.log1.out
crun.python3 -p ~/.conda/envs/popgen admixture PIRE.Lle.Ham.rescaled.preHWE.bed 2 --cv > PIRE.Lle.Ham.rescaled.preHWE.log2.out
crun.python3 -p ~/.conda/envs/popgen admixture PIRE.Lle.Ham.rescaled.preHWE.bed 3 --cv > PIRE.Lle.Ham.rescaled.preHWE.log3.out
crun.python3 -p ~/.conda/envs/popgen admixture PIRE.Lle.Ham.rescaled.preHWE.bed 4 --cv > PIRE.Lle.Ham.rescaled.preHWE.log4.out
crun.python3 -p ~/.conda/envs/popgen admixture PIRE.Lle.Ham.rescaled.preHWE.bed 5 --cv > PIRE.Lle.Ham.rescaled.preHWE.log5.out

exit
```

Copied `*.eigenval`, `*.eigenvec` & `*.Q` files to local computer. Ran `pire_cssl_data_processing/scripts/popgen_analyses/pop_structure.R` on local computer to visualize PCA & ADMIXTURE results (figures in `/home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/pop_structure`).

---

## 16. Filter VCF file for HWE

PCA & ADMIXTURE showed cryptic structure. AHam (Albatross individuals) were roughly split 50/50 between two demes (identified as demes "A" and "B" for now). CNas (contemporary individuals) were almost entirely assigned to deme "B" (with one exception).
  * Species IDS are unknown at this point, however Roy did barcode one contemporary individual from this site (and some from other sites). CNas_054 barcoded to *Equulites laterofenestra*, an identification confirmed by Kent & Jem morphologically as well (they also looked at a few other contemporary individuals). Thus, we are tentatively considering deme "B" to be *Equulites laterofenestra* for now (although it is coded in the `popmap` file as "Lle-CNas-B" or "Lle-AHam-B").
    * The other individuals (from other sites, NOT Hamilo Cove) Roy barcoded came out as *Equulites leuciscus* (aka *Leiognathus leuciscus*). We tentatively expect non-Ela individuals here (those assigned to deme "A") to also be truly *Leiognathus leuciscus* as well, but cannot be sure.

List of AHam & CNas individuals in deme "A" (20 in all):
  * **Albatross (19 total):** AHam_001, AHam_002, AHam_003, AHam_004, AHam_011, AHam_012, AHam_014, AHam_015, AHam_016, AHam_017, AHam_018, AHam_023, AHam_024, AHam_025, AHam_026, AHam_027, AHam_028, AHam_029 & AHam_030
  * **Contemporary (1 total):** CNas_043

Adjusted popmap file to reflect new structure.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/filterVCF

cp ../mapDamageBAM/popmap.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled ./popmap.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.HWEsplit

#added -A or -B to end of pop assignment (second column) to assign individual to either group A or group B.
```

Made a copy of the `config.fltr.ind.cssl` file called `config.fltr.ind.cssl.HWE` with correct file paths, extensions, and filters.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/filterVCF

cp config.fltr.ind.cssl ./config.fltr.ind.cssl.HWE
```

```sh
fltrVCF Settings, run fltrVCF -h for description of settings
        # Paths assume you are in `filterVCF dir` when running fltrVCF, change as necessary
        fltrVCF -f 18 17                                                                                   # order to run filters in
        fltrVCF -c ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled                                        # cutoffs, ie ref description
        fltrVCF -b ../mapDamageBAM                                                                         # path to *.bam files
        fltrVCF -R ../../scripts/fltrVCF/scripts                                                           # path to fltrVCF R scripts
        fltrVCF -d ../mapDamageBAM/mapped.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.bed             # bed file used in genotyping
        fltrVCF -v Lle.A.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.Fltr07.18.vcf                             # vcf file to filter
        fltrVCF -g ../mapDamageBAM/reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.fasta        # reference genome
        fltrVCF -p ../mapDamageBAM/popmap.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.HWEsplit        # popmap file
        fltrVCF -w ../../scripts/fltrVCF/filter_hwe_by_pop_HPC.pl                                          # path to HWE filter script
        fltrVCF -r ../../scripts/rad_haplotyper/rad_haplotyper.pl                                          # path to rad_haplotyper script
        fltrVCF -o Lle.A                                                                                   # prefix on output files, use to track settings
        fltrVCF -t 40                                                                                      # number of threads [1]

        17 vcftools --missing-sites     0.5             #Remove sites with more data missing in a pop sample. [0.5]
        18 filter_hwe_by_pop_HPC        0.001           #Remove sites with <p in test for HWE by pop sample. Adjust based upon sample size [0.001]
```

Did not change the filter settings.

Ran [`fltrVCF.sbatch`](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/fltrVCF.sbatch).

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/filterVCF

sbatch ../../scripts/fltrVCF.sbatch config.fltr.ind.cssl.HWE 
```

---

## 17. Make a VCF file with monomorphic loci

Created a `mkVCF_monomorphic` dir to make an "all sites" VCF (with monomorphic loci included) and moved necessary files over.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus

mkdir mKVCF_monomorphic

ln mapDamage/*bam* mkVCF_monomorphic
cp mapDamgeBAM/*fasta mkVCF_monomorphic
cp mapDamageBAM/config.6.cssl.mkVCF.rescale mkVCF_monomorphic/config.6.cssl.rescale.monomorphic
```

Changed the config file so that the last setting (monomorphic) is set to yes.

```
yes      freebayes    --report-monomorphic (no|yes)          Report even loci which appear to be monomorphic, and report allconsidered alleles,
```

Genotyped with `dDocentHPC.sbatch`.

```sh
cd /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/leiognathus_leuciscus/mkVCF_monomorphic

sbatch ../../../dDocentHPC/dDocentHPC.sbatch mkVCF config.6.cssl.rescale.monomorphic 
```

---

## 18. Filter VCF with monomorphic loci

*Everything from this step down was run in `/archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl`.*

Set-up filtering for monomorphic sites only.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic

cp /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/config.fltr.ind.cssl.mono .
```

Updated the `config.fltr.ind.cssl.mono` file with file paths and file extensions based on your species. The VCF path should point to the "all sites" VCF file you just made. **The settings for filters 04, 14, 05, 16, 13 & 17 should match the settings used when filtering the original VCF file.**

```
fltrVCF Settings, run fltrVCF -h for description of settings
        # Paths assume you are in `filterVCF dir` when running fltrVCF, change as necessary
        fltrVCF -f 01 02 04 14 05 16 04 13 05 16 17                                                                  # order to run filters in
        fltrVCF -c ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled                                                  # cutoffs, ie ref description
        fltrVCF -b ../mapDamageBAM                                                                                   # path to *.bam files
        fltrVCF -R /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/fltrVCF/scripts                     # path to fltrVCF R scripts
        fltrVCF -d ../mapDamageBAM/mapped.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.bed                       # bed file used in genotyping
        fltrVCF -v TotalRawSNPs.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.vcf                                 # vcf file to filter
        fltrVCF -g reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.fasta                                  # reference genome
        fltrVCF -p ../filterVCF/popmap.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.HWEsplit                     # popmap file
        fltrVCF -w /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/fltrVCF/filter_hwe_by_pop_HPC.pl    # path to HWE filter script
        fltrVCF -r /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/rad_haplotyper/rad_haplotyper.pl    # path to rad_haplotyper script
        fltrVCF -o lle.mono                                                                                          # prefix on output files, use to track settings
        fltrVCF -t 40                                                                                                # number of threads [1]

        01 vcftools --min-alleles       1               #Remove sites with less alleles [2]
        01 vcftools --max-alleles       1               #Remove sites with more alleles [2]
        02 vcftools --remove-indels                     #Remove sites with indels.  Not adjustable
        04 vcftools --min-meanDP        5:15            #Remove sites with lower mean depth [15]
        05 vcftools --max-missing       0.55:0.6        #Remove sites with at least 1 - value missing data (1 = no missing data) [0.5]
        13 vcftools --max-meanDP        400             #Remove sites with higher mean depth [250]
        14 vcftools --minDP             5               #Code genotypes with lesser depth of coverage as NA [5]
        16 vcftools --missing-indv      0.6:0.5         #Remove individuals with more missing data. [0.5]
        17 vcftools --missing-sites     0.5             #Remove sites with more data missing in a pop sample. [0.5]
```

Ran [`fltrVCF.sbatch`](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/fltrVCF.sbatch) for monomorphic sites.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic

#before running, make sure the config file is updated with file paths and file extensions based on your species
#VCF file should be the VCF file made after the "make monomorphic VCF" step
#settings for filters 04, 14, 05, 16, 13 & 17 should match the settings used when filtering the original VCF file (step 10)
sbatch /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/fltrVCF.sbatch config.fltr.ind.cssl.mono
```

Set-up filtering for polymorphic sites only.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic

mkdir polymorphic_filter
cd polymorphic_filter

cp /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/config.fltr.ind.cssl.poly .
```

Updated the `config.fltr.ind.cssl.poly` file with file paths and file extensions based on your species. The VCF path should point to the "all sites" VCF file you just made AND the HWEsplit popmap you made if you had any cryptic population structure. **The settings for all your filters should match the settings used when filtering the original VCF file.**

```
fltrVCF Settings, run fltrVCF -h for description of settings
        # Paths assume you are in `filterVCF dir` when running fltrVCF, change as necessary
        fltrVCF -f 01 02 03 04 14 07 05 16 15 06 11 09 10 04 13 05 16 07 18 17                                       # order to run filters in
        fltrVCF -c ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled                                                  # cutoffs, ie ref description
        fltrVCF -b ../../mapDamageBAM                                                                                # path to *.bam files
        fltrVCF -R /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/fltrVCF/scripts                     # path to fltrVCF R scripts
        fltrVCF -d ../../mapDamageBAM/mapped.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.bed                    # bed file used in genotyping
        fltrVCF -v ../TotalRawSNPs.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.vcf                              # vcf file to filter
        fltrVCF -g ../reference.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.fasta                               # reference genome
        fltrVCF -p ../../filterVCF/popmap.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.HWEsplit                  # popmap file
        fltrVCF -w /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/fltrVCF/filter_hwe_by_pop_HPC.pl    # path to HWE filter script
        fltrVCF -r /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/rad_haplotyper/rad_haplotyper.pl    # path to rad_haplotyper script
        fltrVCF -o lle.poly                                                                                          # prefix on output files, use to track settings
        fltrVCF -t 40                                                                                                # number of threads [1]

        01 vcftools --min-alleles       2               #Remove sites with less alleles [2]
        01 vcftools --max-alleles       2               #Remove sites with more alleles [2]
        02 vcftools --remove-indels                     #Remove sites with indels.  Not adjustable
        03 vcftools --minQ              100             #Remove sites with lower QUAL [20]
        04 vcftools --min-meanDP        5:15            #Remove sites with lower mean depth [15]
        05 vcftools --max-missing       0.55:0.6        #Remove sites with at least 1 - value missing data (1 = no missing data) [0.5]
        06 vcffilter AB min             0.375           #Remove sites with equal or lower allele balance [0.2]
        06 vcffilter AB max             0.625           #Remove sites with equal or lower allele balance [0.8]
        06 vcffilter AB nohet           0               #Keep sites with AB=0. Not adjustable
        07 vcffilter AC min             0               #Remove sites with equal or lower MINOR allele count [3]
        09 vcffilter MQM/MQMR min       0.25            #Remove sites where the difference in the ratio of mean mapping quality between REF and ALT alleles is greater than this proportion from 1. Ex: 0 means the mapping quality must be equal between REF and ALTERNATE. Smaller numbers are more stringent. Keep sites where the following is true: 1-X < MQM/MQMR < 1/(1-X) [0.1]
        10 vcffilter PAIRED                             #Remove sites where one of the alleles is only supported by reads that are not properly paired (see SAM format specification). Not adjustable
        11 vcffilter QUAL/DP min        0.2             #Remove sites where the ratio of QUAL to DP is deemed to be too low. [0.25]
        13 vcftools --max-meanDP        400             #Remove sites with higher mean depth [250]
        14 vcftools --minDP             5               #Code genotypes with lesser depth of coverage as NA [5]
        15 vcftools --maf               0               #Remove sites with lesser minor allele frequency.  Adjust based upon sample size. [0.005]
        15 vcftools --max-maf           1               #Remove sites with greater minor allele frequency.  Adjust based upon sample size. [0.995]
        16 vcftools --missing-indv      0.6:0.5         #Remove individuals with more missing data. [0.5]
        17 vcftools --missing-sites     0.5             #Remove sites with more data missing in a pop sample. [0.5]
        18 filter_hwe_by_pop_HPC        0.001           #Remove sites with <p in test for HWE by pop sample. Adjust based upon sample size [0.001]
  ```

Ran [`fltrVCF.sbatch`](https://github.com/philippinespire/pire_cssl_data_processing/blob/main/scripts/fltrVCF.sbatch) for polymorphic sites.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic/polymorphic_filter

sbatch /home/e1garcia/shotgun_PIRE/pire_cssl_data_processing/scripts/fltrVCF.sbatch config.fltr.ind.cssl.poly 
```

---

## 19. Merge monomorphic & polymorphic VCF files

Checked the *filtered* monomorphic & polymorphic VCF files to make sure that filtering removed the same individuals.
* **mono.VCF filtering removed:** AHam_008, AHam_013, CNas_063
* **poly.VCF filtering removed:** AHam_008, AHam_013, CNas_063

Same individuals removed from both, so okay to go ahead.

Created `indv_missing.txt` in `mkVCF_monomorphic` directory. This is a list of all the individuals removed from either  file (total of XX for *spp*). Used this list to make sure number of individuals matched in both filtered VCFs.

Sorted each VCF file.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic

module load container_env bcftools
bash

#sort monomorphic
crun vcf-sort lle.mono.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.11.vcf > lle.mono.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.11.recode.sorted.vcf

#sort polymorphic
cd polymorphic_filter

crun vcf-sort lle.poly.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.20.recode.vcf > lle.poly.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.20.recode.sorted.vcf

exit
```

Zipped each VCF file.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic

module load container_env samtools
bash

#zip monomorphic
crun bgzip -c lle.mono.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.11.recode.sorted.vcf > lle.mono.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.11.recode.sorted.vcf.gz

#zip polymorphic
cd polymorphic_filter

crun bgzip -c lle.poly.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.20.recode.sorted.vcf > lle.poly.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.20.recode.sorted.vcf.gz

exit
```

Indexed each VCF file.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic

module load container_env samtools
bash

#index monomorphic
crun tabix lle.mono.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.11.recode.sorted.vcf.gz

#index polymorphic
cd polymorphic_filter
crun tabix lle.poly.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate.Fltr17.20.recode.sorted.vcf.gz
```

Merged files.

```sh
cd /archive/carpenterlab/pire/pire_leiognathus_leuciscus_cssl/mkVCF_monomorphic

module unload samtools
module load container_env bcftools
bash

mv ppolymorphic_filter *Fltr17.20*sorted.vcf.gz* .

crun bcftools concat --allow-overlaps lle.mono.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.11.recode.sorted.vcf.gz lle.poly.ssl.Lle-C-3NR-R1R2ORPH-contam-noIsolate-rescaled.Fltr17.20.recode.sorted.vcf.gz -O z -o lle.all.recode.sorted.vcf.gz

exit

module load samtools
bash 
crun tabix lle.all.recode.sorted.vcf.gz #index all sites VCF for downstream analyses

exit
```
