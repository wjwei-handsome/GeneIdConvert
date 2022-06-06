from utils import joint_valid_from_to
import pandas as pd
# tmp =tmp.groupby(['from_genome_x', 'to_genome_x']).max()


def GetReciprocalSyntenyFile(from_g: str, to_g: str):
    if joint_valid_from_to(from_g, to_g):
        query_from_target_to_synteny_file = f'01_synteny/{from_g}.{to_g}.anchor'
        query_to_target_from_synteny_file = f'01_synteny/{to_g}.{from_g}.anchor'
    else:
        exit('genome not valid')
    return query_from_target_to_synteny_file, query_to_target_from_synteny_file


def read_synteny_file(from_genome: str, to_genome: str) -> pd.DataFrame:

    query_from_target_to_synteny_file, query_to_target_from_synteny_file = GetReciprocalSyntenyFile(
        from_genome, to_genome)
    exclude_fwd = [i for i, line in enumerate(
        open(query_from_target_to_synteny_file)) if line.startswith('#')]
    exclude_rev = [i for i, line in enumerate(
        open(query_to_target_from_synteny_file)) if line.startswith('#')]

    fwd_data = pd.read_csv(
        query_from_target_to_synteny_file, sep='\t', header=None, skiprows=exclude_fwd, usecols=[0, 1])
    rev_data = pd.read_csv(
        query_to_target_from_synteny_file, sep='\t', header=None, skiprows=exclude_rev, usecols=[0, 1])
    colnames = ['from_genome', 'to_genome']

    fwd_data.columns = colnames
    rev_data.columns = colnames

    merged_df = pd.merge(fwd_data, rev_data, left_on='from_genome',
                         right_on='to_genome', how='outer')
    merged_df = merged_df.loc[merged_df.to_genome_x == merged_df.from_genome_y]

    return merged_df
