'''
This file converts the vcf from table format to sql format
'''



from sqlalchemy import create_engine
import sqlalchemy
import pandas as pd
from collections import defaultdict
import gzip,os
import argparse
import numpy as np


parser = argparse.ArgumentParser(description='transfer vcf table format to sql')

parser.add_argument('-t','--table',action='store',dest='table',help='input vcf in table format')
parser.add_argument('-d','--db',action='store',dest='database',help='input database file')
parser.add_argument('-s','--sql',action='store',dest='sql',help='output sql file')

args = parser.parse_args()


# tab_file = '/lustre/scratch/lis262/NIH_UDP/p01_FSGS_cohort/f04_vep_chr_tab/fsgs.chrY.vep.tsv.gz'
# varcards_file = '/hpc/grid/wip_drm_targetsciences/projects/Varcards_GRCh38/hg38.chrY.extreme.xls.gz'
# sql_file = '/lustre/scratch/lis262/NIH_UDP/p01_FSGS_cohort/f05_vep_chr_sql/fsgs.chrY.vep.sql'
tab_file = args.table
varcards_file = args.database
sql_file = args.sql

if os.path.exists(sql_file):
    os.remove(sql_file)

def get_dtype(columns):
    '''
    this function change the data type of some columns in sql
    '''
    type_dict = defaultdict()
    for k in columns:
        if k.endswith('_AF') or k.endswith('_score'):
            dtype = sqlalchemy.FLOAT
        elif k in ['pos']:
            dtype = sqlalchemy.INTEGER
        else:
            dtype = sqlalchemy.TEXT
        type_dict[k] = dtype
    return type_dict


#--------------- this part dump variants from study into sql
# get column names
with gzip.open(tab_file, 'rt') as f:
    tab_cols = f.readline().strip().split('\t')
type_dict = get_dtype(tab_cols)
# build {column:dtype} dictionary
study_name = tab_file.split('/')[-1].split('.')[0]
table = study_name
primary_key = ['chr','pos','ref','alt']
df = pd.read_csv(tab_file, sep='\t', header=0,index_col=primary_key,low_memory=False)
df = df.replace('-',np.nan)

# connect to sql file
db_url = 'sqlite:///' + sql_file
engine = create_engine(db_url)

# write the case to datbase
columns = df.columns.tolist()
columns = ['_'.join(i.split('-')) for i in columns]
df.columns = columns
df.to_sql(table, engine, index_label = primary_key,
                    if_exists='replace',dtype=type_dict)



def get_varcards_dtype(columns, funct_columns):
    '''
    this function change the data type of some columns in sql
    '''
    type_dict = defaultdict()
    for k in funct_columns:
        type_dict[k+'_score'] = sqlalchemy.FLOAT
        type_dict[k+'_pred'] = sqlalchemy.TEXT
    
    for k in columns:
        if k in funct_columns:
            continue
        if k in ['Kaviar'] or 'ExAC' in k or 'gnomAD' in k  \
             or 'HRC' in k or k.endswith('_score'):
            dtype = sqlalchemy.FLOAT
        elif k in ['Start']:
            dtype = sqlalchemy.INTEGER
        else:
            dtype = sqlalchemy.TEXT
        type_dict[k] = dtype
    return type_dict



def split_score_pred(df, columns):
    '''
    the varcards prediction is in the format of score:pred
    this function split them and delete the origional merged one
    '''
    for column in columns:
        if column in df.columns:
            score = column + '_score'
            pred = column + '_pred'
            df[score] = df[column].map(lambda x: x.split(':')[0] if ':' in x else'-')
            df[pred] = df[column].map(lambda x: x.split(':')[1] if ':' in x else '-')
            del df[column]
    return df


#--------------- this part dump variants from varcard db into sql
# read varcard database
with gzip.open(varcards_file,'rt') as in_f:
    varcards_col = in_f.readline().strip().split('\t')
# drop some columns
drop_cols = ['End','GeneFullName.refGene',
            'GeneFunction.refGene','MGI.refGene'] + \
            [i for i in varcards_col if i.startswith('1000g')]
use_cols = [i for i in varcards_col if i not in drop_cols]
use_idxs = [varcards_col.index(i) for i in use_cols]
v_primary_key = ['Chr','Start','Ref','Alt']

# get functional prediciton columns
funct_start_idx = use_cols.index('SIFT')
funct_end_idx = use_cols.index('ReVe')
funct_columns = use_cols[funct_start_idx:funct_end_idx+1]

db_type_dict = get_varcards_dtype(use_cols, funct_columns)

# clean the table
connection = engine.raw_connection()
cursor = connection.cursor()
command = "DROP TABLE IF EXISTS varcards;"
cursor.execute(command)
connection.commit()
cursor.close()
# loop to output 
parser = pd.read_csv(varcards_file, sep='\t',header=0,
                     index_col=v_primary_key,usecols=use_cols,
                     encoding='iso-8859-1',chunksize=100000)


for v_df in parser:
    v_df = split_score_pred(v_df, funct_columns)
    v_df = v_df.replace('-',np.nan)
    v_df = v_df.replace('.',np.nan)
    for k,v in db_type_dict.items():
        if v == sqlalchemy.FLOAT:
            v_df[k] = pd.to_numeric(v_df[k],errors='coerce')
    v_df.to_sql('varcards',engine,index_label=v_primary_key,
                       if_exists='append',dtype=db_type_dict)