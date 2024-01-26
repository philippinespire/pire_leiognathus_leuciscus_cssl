###running on Rutgers Amarel cluster

###in bash - activate conda momi envelope, generate ind2pop files associating individuals with populations

conda activate momi-env

grep "#CHROM" lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.vcf | cut -d$'\t' -f10- | tr '\t' '\n' > ind.txt 

cat ind.txt | sed 's/Lle-//g' | sed 's/_.*//g' > pop.txt 

paste ind.txt pop.txt > ind2pop.txt 

###enter python3, do allele counts

import momi

vcffile="lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.vcf.gz"
ind2popfile="ind2pop.txt"
outfile="lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.allelecounts" 

with open(ind2popfile) as f: 
	ind2popdict = dict([l.split() for l in f]) 

counts=momi.SnpAlleleCounts.read_vcf(vcf_file=vcffile,ind2pop=ind2popdict,ancestral_alleles=False) 

counts.dump(outfile)

###exit python, zip allele counts file, extract SFS

gzip lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.allelecounts

python -m momi.extract_sfs /home/br450/pire_cssl_data_processing/leiognathus_leuciscus/momi2/lle.sfs.gz 100 /home/br450/pire_cssl_data_processing/leiognathus_leuciscus/momi2/lle.D.ssl.Lle-C-3NR-R1R2ORPH-contam-noisolate-off.Fltr17.1.recode.allelecounts.gz

###enter python, import sfs and examine to see if import worked/get % missing data

import momi

sfsfile="/home/br450/pire_cssl_data_processing/leiognathus_leuciscus/momi2/lle.sfs.gz"
sfs = momi.Sfs.load(sfsfile)

print("populations", sfs.populations)
print("percent missing data per population", sfs.p_missing)
print("length", sfs.length)

### some preliminary models for inference

#set model for inference - constant size, contemp samples only
NeConstant=1e4
model_inf_constant_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_constant_contemp.set_data(sfs, length=102300)

#set parameter to infer - inferring constant pop size using contemporary data only
model_inf_constant_contemp.add_size_param("n_constant")
model_inf_constant_contemp.add_leaf("CNas",N="n_constant")
model_inf_constant_contemp.optimize(method="TNC")

#            fun: 0.17717824918254976
#             jac: array([8.39148395e-10])
#   kl_divergence: 0.17717824918254976
#  log_likelihood: -127941.78580176886
#         message: 'Converged (|f_n-f_(n-1)| ~= 0)'
#            nfev: 13
#             nit: 4
#      parameters: ParamsDict({'n_constant': 426446.15603132})
#          status: 1
#         success: True
#               x: array([12.96324139])

#set model for inference - constant size, contemp samples only
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_constant_temporal.set_data(sfs, length=102300)
#set parameter to infer - inferring constant pop size using temporal data
model_inf_constant_temporal.add_size_param("n_constant")
model_inf_constant_temporal.add_leaf("CNas",N="n_constant")
model_inf_constant_temporal.add_leaf("AHam",N="n_constant",t=109)
model_inf_constant_temporal.move_lineages("CNas","AHam",t=110)
model_inf_constant_temporal.optimize(method="TNC")

#             fun: 3.725472089935075
#             jac: array([0.40803155])
#   kl_divergence: 3.725472089935075
#  log_likelihood: -344750.53193931637
#         message: 'Converged (|x_n-x_(n-1)| ~= 0)'
#            nfev: 96
#             nit: 4
#      parameters: ParamsDict({'n_constant': 10677.553419719381})
#          status: 2
#         success: True
#               x: array([9.27589901])

######

#model with recent size change, contemp samples only

model_inf_change_contemp =  momi.DemographicModel(N_e=1e4, gen_time=1, muts_per_gen=2.5e-8)
model_inf_change_contemp.set_data(sfs,length=102300)
model_inf_change_contemp.add_size_param("n_alb")
model_inf_change_contemp.add_size_param("n_bot")
model_inf_change_contemp.add_time_param("t_bot",upper=1e2)
model_inf_change_contemp.add_leaf("CNas",N="n_bot")
model_inf_change_contemp.set_size("CNas", N="n_alb", t="t_bot")
model_inf_change_contemp.optimize(method="TNC")

#             fun: 0.1573366711703606
#             jac: array([-8.01577831e-06, -7.65389766e-06,  6.85985379e-05])
#   kl_divergence: 0.1573366711703606
#  log_likelihood: -127239.53283118345
#         message: 'Converged (|f_n-f_(n-1)| ~= 0)'
#            nfev: 14
#             nit: 6
#      parameters: ParamsDict({'n_alb': 429178.6631071992, 'n_bot': 3042.421645958911, 't_bot': 34.650840655949125})
#          status: 1
#         success: True
#               x: array([12.96962858,  8.02040907, -0.63442259])


#model with recent size change, temporal samples

model_inf_change_temporal =  momi.DemographicModel(N_e=1e4, gen_time=1, muts_per_gen=2.5e-8)
model_inf_change_temporal.set_data(sfs,length=102300)
model_inf_change_temporal.add_size_param("n_alb")
model_inf_change_temporal.add_size_param("n_bot")
model_inf_change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_change_temporal.add_leaf("CNas",N="n_bot")
model_inf_change_temporal.set_size("CNas", N="n_alb", t="t_bot")
model_inf_change_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_change_temporal.move_lineages("CNas","AHam",t=110)
model_inf_change_temporal.optimize(method="TNC")

#             fun: 3.540431263458293
#             jac: array([ 0.44079185,  0.10002892, -0.02129006])
#   kl_divergence: 3.540431263458293
#  log_likelihood: -337244.16577245924
#         message: 'Converged (|x_n-x_(n-1)| ~= 0)'
#            nfev: 55
#             nit: 2
#      parameters: ParamsDict({'n_alb': 7559.447087473259, 'n_bot': 5712.851132367439, 't_bot': 12.038419075830616})
#          status: 2
#         success: True
#               x: array([ 8.93055333,  8.6504735 , -1.98879701])
