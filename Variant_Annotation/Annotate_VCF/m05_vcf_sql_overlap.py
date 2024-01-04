from sqlalchemy import create_engine
import sqlite3,gzip
import argparse

parser = argparse.ArgumentParser(description='overlap vcf with database')


parser.add_argument('-o','--out',action='store',dest='out_dir',help='output directory')
parser.add_argument('-i','--sql',action='store',dest='sql',help='input sql file')

args = parser.parse_args()


sql_file = args.sql
out_dir = args.out_dir
# sql_file = '/lustre/scratch/lis262/NIH_UDP/p01_FSGS_cohort/f05_vep_chr_sql/fsgs.chrY.vep.sql'
# out_dir = '/lustre/scratch/lis262/NIH_UDP/p01_FSGS_cohort/f06_vep_chr_varcards_tab'

out_overlap_file = out_dir + '/' + sql_file.split('/')[-1][:-3] + 'tsv.gz'

conn = sqlite3.connect(sql_file)
cur = conn.cursor()

# get sample columns
study_name = sql_file.split('/')[-1].split('.')[0]
cursor = conn.execute('select * from {t} limit 2'.format(t=study_name))
fsgs_col = list(map(lambda x: x[0], cursor.description))
# get varcard columns
cursor = conn.execute('select * from varcards limit 2')
varcards_col = list(map(lambda x: x[0], cursor.description))
varcards_col = [i for i in varcards_col if i not in \
               ['Chr','Start','Ref','Alt','Ratio_of_tools (ReVe,gt,0.7)']]
cursor.close()

# define output columns
new_fsgs_col = [study_name + '.' + i for i in fsgs_col]
new_varcards_col = []
for c in varcards_col:
    if '.' in c or '-' in c or '+' in c:
        new_varcards_col.append('\"' + c + '\"')
    else:
        new_varcards_col.append('varcards.' + c)
all_columns = new_fsgs_col + new_varcards_col
out_header = [i.split('.')[1] for i in new_fsgs_col] + varcards_col

# output to file
tmp_dir = out_dir
cursor = conn.cursor()
cmd = "PRAGMA temp_store_directory = '{tmp}'".format(tmp=tmp_dir)
cursor.execute(cmd)
cmd = ('SELECT {c} FROM {s} \
        LEFT JOIN varcards \
        ON {s}.chr = varcards.Chr AND \
        {s}.pos = varcards.Start AND \
        {s}.ref = varcards.Ref AND \
        {s}.alt = varcards.Alt').format(
        c=','.join(all_columns),s=study_name)
cursor.execute(cmd)

rows = cursor.fetchall()
with gzip.open(out_overlap_file,'wb') as out:
    header_line = '\t'.join(out_header)+'\n'
    out.write(header_line.encode('utf-8'))
    for row in rows:
        record = [str(i) for i in row]
        record = ['-' if i == 'None' else i for i in record]
        outline = '\t'.join(record) + '\n'
        out.write(outline.encode('utf-8'))