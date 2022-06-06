import pandas as pd
from utils import joint_valid_from_to


def read_orthogroup_file(orthogroup_file_path: str, from_genome: str, to_genome: str) -> pd.DataFrame:
    if joint_valid_from_to(from_genome, to_genome):
        used_cols = ['Orthogroup', from_genome, to_genome]
        orthogroup_df = pd.read_csv(
            orthogroup_file_path, sep='\t', usecols=used_cols)[used_cols]
        orthogroup_df.columns = ['Orthogroup', 'from_genome', 'to_genome']
        orthogroup_df.fillna('None', inplace=True)
        return orthogroup_df
    else:
        exit('genome not valid')
