def joint_valid_from_to(from_g: str, to_g: str) -> bool:
    VALID_GENOMES_NAMES = ['Zd-TEO07',
                           'Zh-TEO10',
                           'Zl-TEO02',
                           'Zm-B73v5',
                           'Zm-SK',
                           'Zn-TEO11',
                           'Zv-TEO01',
                           'Zx-TEO09']  # define valid genome name
    if from_g in VALID_GENOMES_NAMES and to_g in VALID_GENOMES_NAMES:
        return True
    else:
        return False


def get_gene_info_dict(bed_file_path: str) -> dict:
    gene_info_dict = {}
    with open(bed_file_path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            gene_info_dict[line[4]] = line[:]
    return gene_info_dict


def intersectBed(lst1, lst2):
    '''
    Return intersection of two bed regions.

    Parameters
    ----------
    lst1 : list
            The 1st genomic region. List of chrom, start, end.
            Example: ['chr1',10, 100]

    lst2 : list
             The 2nd genomic region. List of chrom, start, end.
             Example: ['chr1',50, 120]

    Examples
    --------
    >>> intersectBed(['chr1',10, 100],['chr1',50, 120])
    ('chr1', 50, 100)
    >>> intersectBed(['chr1',10, 100],['chr1',20, 30])
    ('chr1', 20, 30)

    '''
    (chr1, st1, end1) = lst1
    (chr2, st2, end2) = lst2
    if int(st1) > int(end1) or int(st2) > int(end2):
        raise Exception("Start cannot be larger than end")
    if chr1 != chr2:
        return None
    if int(st1) > int(end2) or int(end1) < int(st2):
        return None
    return (chr1, max(st1, st2), min(end1, end2))


def map_coordinates(mapping, q_chr, q_start, q_end, q_strand='+'):
    '''
    Map coordinates from source (i.e. original) assembly to target (i.e. new) assembly.

    Parameters
    ----------
    mapping : dict
            Dictionary with source chrom name as key, IntervalTree object as value.

    q_chr : str
            Chromosome ID of query interval

    q_start : int
            Start position of query interval.

    q_end : int
            End position of query interval.

    q_strand : str
            Strand of query interval.

    '''

    matches = []
    complement = {'+': '-', '-': '+'}

    if q_chr in mapping:
        targets = mapping[q_chr].find(q_start, q_end)
    elif q_chr.replace('chr', '') in mapping:
        targets = mapping[q_chr.replace('chr', '')].find(q_start, q_end)
    elif ('chr' + q_chr) in mapping:
        targets = mapping['chr' + q_chr].find(q_start, q_end)
    else:
        return None
    if len(targets) == 0:
        return None
    elif len(targets) == 1:
        s_start = targets[0].start
        s_end = targets[0].end
        t_chrom = targets[0].value[0]
        # t_chrom = update_chromID(q_chr, t_chrom, chr_style = chrom_style)
        t_start = targets[0].value[1]
        t_end = targets[0].value[2]
        t_strand = targets[0].value[3]

        (chr, real_start, real_end) = intersectBed(
            (q_chr, q_start, q_end), (q_chr, s_start, s_end))
        l_offset = abs(real_start - s_start)
        #r_offset = real_end - s_end
        size = abs(real_end - real_start)

        # matches.append((chr, real_start, real_end, q_strand))
        if t_strand == '+':
            i_start = t_start + l_offset
            if q_strand == '+':
                matches.append((t_chrom, i_start, i_start + size, t_strand))
            else:
                matches.append((t_chrom, i_start, i_start +
                               size, complement[t_strand]))
        elif t_strand == '-':
            i_start = t_end - l_offset - size
            if q_strand == '+':
                matches.append((t_chrom, i_start,	i_start + size, t_strand))
            else:
                matches.append((t_chrom, i_start,	i_start +
                               size, complement[t_strand]))
        else:
            raise Exception(
                "Unknown strand: %s. Can only be '+' or '-'." % q_strand)

    elif len(targets) > 1:
        for t in targets:
            s_start = t.start
            s_end = t.end
            t_chrom = t.value[0]
            # t_chrom = update_chromID(q_chr, t_chrom, chr_style=chrom_style)
            t_start = t.value[1]
            t_end = t.value[2]
            t_strand = t.value[3]

            (chr, real_start, real_end) = intersectBed(
                (q_chr, q_start, q_end), (q_chr, s_start, s_end))

            l_offset = abs(real_start - s_start)
            #r_offset = abs(real_end - s_end)
            size = abs(real_end - real_start)
            # matches.append((chr, real_start, real_end, q_strand))
            if t_strand == '+':
                i_start = t_start + l_offset
                if q_strand == '+':
                    matches.append(
                        (t_chrom, i_start, i_start + size, t_strand))
                else:
                    matches.append(
                        (t_chrom, i_start, i_start + size, complement[t_strand]))
            elif t_strand == '-':
                i_start = t_end - l_offset - size
                if q_strand == '+':
                    matches.append(
                        (t_chrom, i_start,	i_start + size, t_strand))
                else:
                    matches.append(
                        (t_chrom, i_start,	i_start + size, complement[t_strand]))
            else:
                raise Exception(
                    "Unknown strand: %s. Can only be '+' or '-'." % q_strand)

    # if print_match:
    #     print(matches)
        # input: 'chr1',246974830,247024835
        # output: [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' ), ('chr1', 247024833, 247024835, '+'), ('chr1', 249058210, 249058212,'+')]
        # [('chr1', 246974830, 246974833), ('chr1', 248908207, 248908210)]

    return matches
