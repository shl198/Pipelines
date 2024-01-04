
path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/NIH_UDP/p02_UDP11794_fam'

# path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/schizo'


# #----------------- 1. transfer vep to table
# import glob, os
# from natsort import natsorted
# py_file = '/home/lis262/Code/Variant_Analysis/Annotate_VCF/m01_vcf2tab.py'

# vep_path = path + '/f03_vep_chr_vcf'

# vep_files = natsorted(glob.glob(vep_path + '/*.vcf.gz'))
# tab_path = path + '/f04_vep_chr_tab'

# if not os.path.exists(tab_path):
#     os.mkdir(tab_path)

# for vep_file in vep_files:
#     chrom = vep_file.split('.')[-4]
#     tab_file = tab_path + '/' + vep_file.split('/')[-1][:-6] + 'tsv.gz'
#     cmd = ('bsub -q short -R "select[hname!=amrndhl1372]" -J {job} -o {log} \"python {py} -i {vep} -o {tab} \"').format(log=tab_path+'/log.'+chrom, py=py_file,vep=vep_file,tab=tab_file,job=chrom)
#     os.system(cmd)



# #--------------- 2. add gnomad v3 NFE frequency
# import glob, os
# from natsort import natsorted


# tab_path = path + '/f04_vep_chr_tab'
# out_path = path + '/f05_add_gnomadV3_tab'
# if not os.path.exists(out_path):
#     os.mkdir(out_path)

# py_file = '/home/lis262/Code/Variant_Analysis/Annotate_VCF/m02_add_gnomad_v3.py'
# gnomad_v3_path = '/hpc/grid/wip_drm_targetsciences/projects/gnomAD/gnomad_v3'

# tab_files = natsorted(glob.glob(tab_path + '/*.tsv.gz'))
# for tab_fn in tab_files:
#     chrom = tab_fn.split('.')[-4]
#     out_file = out_path + '/' + tab_fn.split('/')[-1]
#     log = out_path + '/log.' + chrom
#     cmd = 'bsub -R "select[hname!=amrndhl1372]" -o {log} -J {chrom} -q short \"python {py} -i {tab} -o {out} -g {gnomad}\"'.format(log=log,chrom=chrom,py=py_file,tab=tab_fn,out=out_file,gnomad=gnomad_v3_path)
#     os.system(cmd)



# #-------------- 3. split coding, regu and noncoding
# import pandas as pd
# import os, glob

# py = '/home/lis262/Code/Variant_Analysis/Annotate_VCF/m03_split_tab.py'

# tab_path = path + '/f05_add_gnomadV3_tab'
# code_path = path + '/f06_vep_chr_code'
# regu_path = path + '/f06_vep_chr_regu'
# ncode_path = path + '/f06_vep_chr_ncode'
# if not os.path.exists(code_path):
# 	os.mkdir(code_path)
# 	os.mkdir(regu_path)
# 	os.mkdir(ncode_path)

# tab_files = glob.glob(tab_path + '/*tsv.gz')
# for tab in tab_files:
# 	chrom = tab.split('.')[-4]
# 	log = code_path + '/log.' + chrom
# 	c = code_path + '/' + tab.split('/')[-1]
# 	r = regu_path + '/' + tab.split('/')[-1]
# 	n = ncode_path + '/' + tab.split('/')[-1]
# 	cmd = 'bsub -R "select[hname!=amrndhl1372]" -o {log} -J {chrom} -q short \"python {py} -i {tab} -c {c} -r {r} -n {n}\"'.format(log=log,chrom=chrom,py=py,tab=tab,c=c,r=r,n=n)
# 	os.system(cmd)
# 	# print(cmd)



# #--------------------------------------------------------
# # 4. transfer tab to sql
# import glob,os
# from natsort import natsorted
# py_file = '/home/lis262/Code/Variant_Analysis/Annotate_VCF/m04_tab2sql.py'

# tab_path = path + '/f06_vep_chr_code'
# tab_files = natsorted(glob.glob(tab_path + '/*.tsv.gz'))

# sql_path = path + '/f07_vep_chr_sql'
# varcards_path = '/hpc/grid/wip_drm_targetsciences/projects/VEP_plugin/Varcards_GRCh38'
# varcards_files = natsorted(glob.glob(varcards_path + '/*.gz'))

# if not os.path.exists(sql_path):
#     os.mkdir(sql_path)

# for tab_file in tab_files:
#     chrom = tab_file.split('.')[-4]
#     try:
#         varcards_file = [i for i in varcards_files if '.'+chrom+'.' in i][0]
#         print()
#     except:
#         print('doesnot find varcard files')
#         continue
#     sql_file = sql_path + '/' + tab_file.split('/')[-1][:-6] + 'sql'
#     cmd = ('bsub -M 21457280 -R "select[hname!=amrndhl1372]" -q short -J {job} -o {log} \"python {py} -t {tab} -s {sql} -d {db}\"').format(py=py_file,
#                     tab = tab_file,sql=sql_file,db=varcards_file,
#                     log=sql_path+'/log.'+chrom,job=chrom)
#     os.system(cmd)



# #--------------------------------------------------------------------
# # 5. OVERLAP VCF TABLE WITH VARCARDS TABLE
# import glob,os
# from natsort import natsorted


# sql_path = path + '/f07_vep_chr_sql'
# out_path = path + '/f08_vep_chr_varcards_tab'

# if not os.path.exists(out_path):
#     os.mkdir(out_path)

# sql_files = natsorted(glob.glob(sql_path + '/*.sql'))

# py_file = '/home/lis262/Code/Variant_Analysis/Annotate_VCF/m05_vcf_sql_overlap.py'

# for sql in sql_files:
#     chrom = sql.split('.')[-3]
#     log = out_path + '/log.' + sql.split('.')[1]
#     cmd = 'bsub -q short -R "select[hname!=amrndhl1372]" -o {log} -J {job} \"python {py} -i {sql} -o {out_dir}\"'.format(
#         py=py_file,log=log,job=chrom,sql=sql,out_dir=out_path)
#     os.system(cmd)



#-----------------------------------------------------------------------
# 6.add regBase to regulatory variants
import glob,os

regBase_path = '/hpc/grid/wip_drm_targetsciences/projects/VEP_plugin/regData/GRCh38'
py = '/home/lis262/Code/Variant_Analysis/Annotate_VCF/m06_add_regBase_Score.py'
regu_path = path + '/f06_vep_chr_regu'
regu_files = glob.glob(regu_path + '/*.tsv.gz')

out_path = path + '/f09_regu_anno'
if not os.path.exists(out_path):
	os.mkdir(out_path)

for regu in regu_files:
	chrom = regu.split('.')[-4]
	log = out_path + '/log.' + chrom
	cmd = ('bsub -q medium -R "select[hname!=amrndhl1372]" -o {log} -J {job} \"python {py} -i {tab} -o {out} -r {ref}\"').format(log=log, job=chrom,py=py,tab=regu,out=out_path,ref=regBase_path)
	os.system(cmd)
