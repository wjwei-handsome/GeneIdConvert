import pandas as pd
from utils import joint_valid_from_to
# import logging
# logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
#                     level=logging.DEBUG)


# TODO   filter e 1e-10

def GetReciprocalBlastFile(from_g: str, to_g: str):
    if joint_valid_from_to(from_g, to_g):
        query_from_target_to_blast_file = f'00_blast/{from_g}.{to_g}.blast'
        query_to_target_from_blast_file = f'00_blast/{to_g}.{from_g}.blast'
    else:
        exit('genome not valid')
    return query_from_target_to_blast_file, query_to_target_from_blast_file


def read_blast_file(from_genome: str, to_genome: str, e_filter: float = 1e-10) -> pd.DataFrame:

    query_from_target_to_blast_file, query_to_target_from_blast_file = GetReciprocalBlastFile(
        from_genome, to_genome)
    fwd_data = pd.read_csv(
        query_from_target_to_blast_file, sep='\t', header=None)
    rev_data = pd.read_csv(
        query_to_target_from_blast_file, sep='\t', header=None)
    column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                    'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    fwd_data.columns = column_names
    rev_data.columns = column_names
    fwd_data = fwd_data.drop_duplicates(subset=['qseqid'], keep='first')
    rev_data = rev_data.drop_duplicates(subset=['qseqid'], keep='first')
    rbh = pd.merge(fwd_data,
                   rev_data[['qseqid', 'sseqid']],
                   left_on='sseqid',
                   right_on='qseqid',
                   how='outer')
    rbh = rbh.loc[rbh.qseqid_x == rbh.sseqid_y]
    rbh = rbh.groupby(['qseqid_x', 'sseqid_x']).max()

    return rbh
