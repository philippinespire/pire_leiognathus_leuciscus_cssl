#### Code used for NeEstimator - species = Lle!

cp /projects/f_mlp195/brendanr/vcf2genepop.pl ./

### convert vcf to genepop
perl vcf2genepop.pl vcf=lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.vcf pops=Lle-CNas,Lle-AHam > lle.genepop

### create file associating loci with "chromosomes"
tail -n +2 lle.genepop | head -1 | tr " " "\n" | awk -F "," {'print $1'} > lle_locfile
tail -n +2 lle.genepop | head -1 | tr " " "\n" | awk -F "_" {'print $2'} > lle_chromfile
paste lle_chromfile lle_locfile > lle_chromlocfile

### NeEst wants shorter sample names
cp lle.genepop lle.rename.genepop
sed -i 's/-LlA.*repr//g' lle.rename.genepop
sed -i 's/-LlC.*repr//g' lle.rename.genepop

### run LD estimates per population
(echo "1"; echo "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/" ; echo "lle.rename.genepop" ; echo "2" ; echo "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/" ; echo "lle_ldout" ; echo "2" ; echo "0.05 0.01"; echo "0") > lle_neest_info
(echo "0 0" ; echo "0" ; echo "0" ; echo "0" ; echo "1" ; echo "1" ; echo "0" ; echo "0" ; echo "0" ; echo "2 lle_chromlocfile") > lle_neest_option
/projects/f_mlp195/brendanr/Ne2-1L i:lle_neest_info o:lle_neest_option

### run temporal estimates
(echo -e "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/lle.rename.genepop" ; echo -e "8\n" ; echo -e "2\n" ; echo -e "0\n" ; echo -e "110\n") | /projects/f_mlp195/brendanr/Ne2-1L

### try with thinning first - no big difference in estimates
vcftools --vcf lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.vcf --recode --thin 100000
perl vcf2genepop.pl vcf=out.recode.vcf pops=Lle-CNas,Lle-AHam > lle.thin.genepop
tail -n +2 lle.thin.genepop | head -1 | tr " " "\n" | awk -F "," {'print $1'} > lle_thin_locfile
tail -n +2 lle.thin.genepop | head -1 | tr " " "\n" | awk -F "_" {'print $2'} > lle_thin_chromfile
paste lle_thin_chromfile lle_thin_locfile > lle_thin_chromlocfile
cp lle.thin.genepop lle.thin.rename.genepop
sed -i 's/-LlA.*repr//g' lle.thin.rename.genepop
sed -i 's/-LlC.*repr//g' lle.thin.rename.genepop

(echo "1"; echo "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/" ; echo "lle.thin.rename.genepop" ; echo "2" ; echo "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/" ; echo "lle_thin_ldout" ; echo "2" ; echo "0.05 0.01"; echo "0") > lle_thin_neest_info
(echo "0 0" ; echo "0" ; echo "0" ; echo "0" ; echo "1" ; echo "1" ; echo "0" ; echo "0" ; echo "0" ; echo "2 lle_thin_chromlocfile") > lle_thin_neest_option
/projects/f_mlp195/brendanr/Ne2-1L i:lle_thin_neest_info o:lle_thin_neest_option

### does order of populations matter for temporal estimate? switched order and it did not change estimates

head -n 3 lle.rename.genepop > lle.switch.genepop
tail -n 94 lle.rename.genepop >> lle.switch.genepop
head -n 30 lle.rename.genepop | tail -n 28 >> lle.switch.genepop
(echo -e "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/lle.switch.genepop" ; echo -e "8\n" ; echo -e "2\n" ; echo -e "0\n" ; echo -e "110\n") | /projects/f_mlp195/brendanr/Ne2-1L

