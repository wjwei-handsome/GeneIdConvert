#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :main.py
@说明        :
@时间        :2022/04/15 22:18:16
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com

requirement files:
1.{genome1}_to_REF-{genome2}.chain.gz
2.query__Zm_{genome1}__target__Zm_{genome2}.anchor
3.query__Zm_{genome1}__target__Zm_{genome2}.blast
4.Zm-{genome2}.bed : ""

'''


import sys
import chain
import blast
import orthogroup
import synteny
import logging
import argparse
import utils
import pybedtools
import pandas as pd
from tqdm import tqdm
logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG)


def GetArgs():
    """recceive the argument

    Returns:
        args: args
    """

    parser = argparse.ArgumentParser(
        description=__doc__, prog='main.py')

    input_group = parser.add_argument_group(title='input options')
    input_group.add_argument('-g', action='store', dest='in_gene', type=str,
                             help='gene name')  # FIXME will delete
    input_group.add_argument('-l', action='store', dest='in_file', type=str,
                             help='gene list file name')
    parser.add_argument('-f', action='store', dest='from_genome', type=str, required=True,
                        help='from this genome')
    parser.add_argument('-t', action='store', dest='to_genome', type=str, required=True,
                        help='convert to this genome')
    parser.add_argument('-o', action='store', dest='out_prefix',
                        type=str, required=True, help='output files prefix')
    return parser.parse_args()


def retrieve_orthogroup_result(orthogroup_df, gene_name: str) -> str:
    result = []
    filter_df = orthogroup_df.loc[orthogroup_df['from_genome'].str.contains(
        gene_name)]
    if filter_df.empty:
        return 'None'
    else:
        for index, row in filter_df.iterrows():
            result.append(row['to_genome'])
        result_str = ','.join(result)
        return result_str


def retrieve_synteny_result(synteny_df, gene_name: str) -> str:
    result = []
    filter_df = synteny_df.loc[synteny_df['from_genome_x'] == gene_name]
    if filter_df.empty:
        return 'None'
    else:
        for index, row in filter_df.iterrows():
            result.append(row['to_genome_x'])
        result_str = ','.join(result)
        return result_str


def retrieve_rbh_result(rbh_df, gene_name: str) -> str:
    result = []
    filter_df = rbh_df.query(f'qseqid_x=="{gene_name}"')
    if filter_df.empty:
        return 'None'
    else:
        for index, row in filter_df.iterrows():
            result.append(row['qseqid_y'])
        result_str = ','.join(result)
        return result_str


def create_bed_string_from_matches(matches: list) -> str:
    '''chr1 111 222 jj jj +\nchr1 111 222 jj jj +\n'''
    raw_result_list = []
    for match in matches:
        formatted_match = list(match[:3])
        formatted_match.append('junk_name')
        formatted_match.append('junk_score')
        formatted_match.append(match[3])
        single_str = ' '.join(map(str, formatted_match))
        raw_result_list.append(single_str)
    result_str = '\n'.join(raw_result_list)
    return result_str


def retrieve_chain_result(chain_result: tuple, gene_name: str) -> str:
    from_gene_chr = from_gene_info_dict[gene_name][0]
    from_gene_start = int(from_gene_info_dict[gene_name][1])
    from_gene_end = int(from_gene_info_dict[gene_name][2])
    from_gene_strand = from_gene_info_dict[gene_name][5]
    from_bx_ivl = chain_result[0]
    # ref to crossmap bed usage
    # NOTE flexable chr name
    matches = utils.map_coordinates(
        from_bx_ivl, from_gene_chr, from_gene_start, from_gene_end, from_gene_strand)
    if matches:
        bedtools_string = create_bed_string_from_matches(matches)
        bedtools_ivl = pybedtools.bedtool.BedTool(
            bedtools_string, from_string=True)
        merged_bedtools_ivl = bedtools_ivl.sort().merge(d=10)
        target_bed = pybedtools.BedTool(target_bed_file_path)
        intersect_ivl = target_bed.intersect(
            merged_bedtools_ivl, nonamecheck=True)
        if intersect_ivl:
            result = []
            intersect_ivl_lines = intersect_ivl.__str__().splitlines()
            for line in intersect_ivl_lines:
                line_chr, line_start, line_end, line_junk, line_name, line_strand = line.split(
                    '\t')
                result.append(
                    f'{line_name}@({line_chr}:{line_start}-{line_end}-{line_strand})')
            str_result = ','.join(result)
            return str_result  # NOTE private var __str__
        else:
            result = []
            for ivl in merged_bedtools_ivl:
                result.append(
                    f'({ivl.chrom}:{ivl.start}-{ivl.end})')
            str_result = ','.join(result)
            return str_result
    else:
        return 'None'


if __name__ == '__main__':
    args = GetArgs()
    if (args.in_gene and args.in_file) or (not args.in_gene and not args.in_file):
        sys.exit('error! file and gene only one!')
    # print(args.in_gene)

    # get gene position and name from bed file
    if utils.joint_valid_from_to(args.from_genome, args.to_genome):
        from_bed_file_path = f'99_genome_gff/{args.from_genome}.bed'
        # to_bed_file_path = f'Zm-{args.to_genome}.bed'
        chain_file_path = f'03_wga_chain/{args.from_genome}.{args.to_genome}.chain.gz'
        # chain_file_path = 'B73v4_to_REF-B73v5.chain.gz'
        target_bed_file_path = f'99_genome_gff/{args.to_genome}.bed'
    else:
        sys.exit('genome not valid')
    from_gene_info_dict = utils.get_gene_info_dict(from_bed_file_path)
    # to_genes_info_dict = utils.get_gene_info_dict(to_bed_file_path)

    # get chain result dict
    # chain_tmp_path = 'B73v442_to_REF-B73v5.chain.gz'
    chain_result = chain.read_chain_file(chain_file_path)

    # get rbh result Dataframe
    logging.info("Read the blast file")
    rbh_result_df = blast.read_blast_file(
        args.from_genome, args.to_genome, e_filter=1e-10)

    # get Orthogroup result
    logging.info("Read the orthogroup file")
    orthogroup_df = orthogroup.read_orthogroup_file(
        '02_ortholog/Orthogroups.tsv', args.from_genome, args.to_genome)

    # get synteny result
    logging.info("Read the synteny file")
    synteny_df = synteny.read_synteny_file(args.from_genome, args.to_genome)

    def get_single_gene_result(gene_name: str):
        # result_pd = pd.DataFrame(colnames=['gene_name', 'crossmap', 'rbh', 'orthogroup'])
        rbh_result = retrieve_rbh_result(rbh_result_df, gene_name)
        orthogroup_result = retrieve_orthogroup_result(
            orthogroup_df, gene_name)
        crossmap_result = retrieve_chain_result(chain_result, gene_name)
        return rbh_result, orthogroup_result, crossmap_result

    def get_multi_gene_result(gene_name_list: list) -> pd.DataFrame:
        result_pd = pd.DataFrame(
            columns=['gene_name', 'crossmap', 'rbh', 'orthogroup'])
        for gene_name in tqdm(gene_name_list):
            rbh_result = retrieve_rbh_result(rbh_result_df, gene_name)
            orthogroup_result = retrieve_orthogroup_result(
                orthogroup_df, gene_name)
            crossmap_result = retrieve_chain_result(chain_result, gene_name)
            synteny_result = retrieve_synteny_result(synteny_df, gene_name)
            raw_result_dict = {'gene_name': gene_name,
                               'rbh': rbh_result,
                               'orthogroup': orthogroup_result,
                               'synteny': synteny_result,
                               'crossmap': crossmap_result}
            result_pd = result_pd.append([raw_result_dict], ignore_index=True)
        return result_pd

    def write_result(gene_list_file, outpre):
        with open(gene_list_file, 'r') as f:
            gene_list = f.read().splitlines()
        result_pd = get_multi_gene_result(gene_list)
        with open(f'{outpre}_raw.tsv', 'w') as f:
            logging.info(f"Write the raw result to {outpre}_raw.tsv")
            result_pd.to_csv(f, sep='\t', index=False)  # write to raw file
        logging.info(f"Process the raw result and write to {outpre}_final.tsv")
        process_df = result_pd.melt(
            id_vars=['gene_name'], var_name='evd', value_name='target')

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
        process_df['target'] = process_df['target'].map(functmp)
        process_df = process_df.explode('target')
        process_df = process_df.drop(
            process_df[process_df.target == 'None'].index).reset_index(drop=True)

        result_df = process_df.groupby(['gene_name', 'target'], as_index=False).agg({
            'evd': lambda x: ','.join(x)})

        result_df.to_csv(f'{outpre}_fianl.tsv', sep='\t', index=False)

    if args.in_file:
        write_result(args.in_file, args.out_prefix)
    elif args.in_gene:
        print(get_single_gene_result(args.in_gene))
