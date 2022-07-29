import tabix
import gzip
import argparse


parser = argparse.ArgumentParser(description='add regBase scores')

parser.add_argument('-i','--input',action='store',dest='tab',help='input tble file')
parser.add_argument('-o','--out',action='store',dest='out_dir',help='output directory')
parser.add_argument('-r','--reg',action='store',dest='reg_dir',help='regBase directory')

args = parser.parse_args()


tab_file = args.tab
regBase_path = args.reg_dir
out_path = args.out_dir

regu_col_fn = regBase_path + '/columns.txt'
with open(regu_col_fn) as c:
    cols = c.readlines()[5:]
    cols = [c.strip().split('\t')[1] for c in cols]


out_file = out_path + '/' + tab_file.split('/')[-1]
chrom = tab_file.split('.')[-4]
reg_ref = regBase_path + '/hg38.{c}.sort.gz'.format(c=chrom)

tb = tabix.open(reg_ref)
with gzip.open(tab_file,'rt') as in_tab, gzip.open(out_file,'wb') as out:
    header = in_tab.readline().strip().split('\t')
    out_header = header + cols
    outline = '\t'.join(out_header) + '\n'
    out.write(outline.encode('utf-8'))
    for line in in_tab:
        item = line.strip().split('\t')
        pos = item[1]
        ref = item[2]
        alt = item[3]
        try:
            records = tb.query(chrom[3:], int(pos)-1, int(pos))
        except:
            continue
        for r in records:
            if r != '':
                if r[3] == ref and r[4] == alt:
                    outline = '\t'.join(item+r[5:]) + '\n'
                    out.write(outline.encode('utf-8'))