import sys
import pandas as pd

from_g = sys.argv[1]
to_g = sys.argv[2]
pre = sys.argv[3]

raw_file_path = f'{pre}/from{from_g}to{to_g}_raw.tsv'

df = pd.read_csv(raw_file_path, sep='\t')
df = df.melt(id_vars=['gene_name'], var_name='evd', value_name='target')


def functmp(x):
    if '@' in x:
        result = []
        tmp1 = x.split(',')
        for i in tmp1:
            g = i.split('@')[0]
            if g not in result:
                result.append(g)
        return result
    else:
        return [item.strip() for item in x.split(',')]


df['target'] = df['target'].map(functmp)
df = df.explode('target')
df = df.drop(df[df.target == 'None'].index).reset_index(drop=True)


result_df = df.groupby(['gene_name', 'target'], as_index=False).agg(
    {'evd': lambda x: ','.join(x)})

result_df.to_csv(f'{pre}/from{from_g}to{to_g}_final.tsv', index=False)
