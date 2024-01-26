#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH --time 10:00:00
#SBATCH --mem-per-cpu 4000
#SBATCH --job-name NeEst_lle
#SBATCH --output NeEst_lle.txt

cd /home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/

perl vcf2genepop.pl vcf=lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.vcf pops=Lle-CNas,Lle-AHam > lle.genepop

tail -n +2 lle.genepop | head -1 | tr " " "\n" | awk -F "," {'print $1'} > lle_locfile

tail -n +2 lle.genepop | head -1 | tr " " "\n" | awk -F "_" {'print $2'} > lle_chromfile

paste lle_chromfile lle_locfile > lle_chromlocfile

cp lle.genepop lle.rename.genepop

sed -i 's/-LlA.*repr//g' lle.rename.genepop

sed -i 's/-LlC.*repr//g' lle.rename.genepop

(echo "1"; echo "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/" ; echo "lle.rename.genepop" ; echo "2" ; echo "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/" ; echo "lle_ldout" ; echo "2" ; echo "0.05 0.01"; echo "0") > lle_neest_info
(echo "0 0" ; echo "0" ; echo "0" ; echo "0" ; echo "1" ; echo "1" ; echo "0" ; echo "0" ; echo "0" ; echo "2 lle_chromlocfile") > lle_neest_option
/projects/f_mlp195/brendanr/Ne2-1L i:lle_neest_info o:lle_neest_option
(echo -e "/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/NeEst/lle.rename.genepop" ; echo -e "8\n" ; echo -e "5\n" ; echo -e "110\n" ; echo -e "0\n") | Ne2-1L

