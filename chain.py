import ireader
import logging
from intersection import Interval, Intersecter
logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG)


def read_chain_file(chain_file, print_table=False):
    '''
    Read chain file.

    Parameters
    ----------
    chain_file : file
            Chain format file. Input chain_file could be either plain text, compressed file
            (".gz",".Z", ".z", ".bz", ".bz2", ".bzip2"), or a URL pointing to the chain file
            ("http://","https://", "ftp://"). If url was used, chain file must be plain text.

    print_table : bool, optional
            Print mappings in human readable table.

    Returns
    -------
    maps : dict
            Dictionary with source chrom name as key, IntervalTree object as value. An
            IntervalTree contains many intervals. An interval is a start and end position
            and a value. eg. Interval(11, 12, strand="-", value = "abc")

    target_chromSize : dict
            Chromosome sizes of target genome

    source_chromSize : dict
            Chromosome sizes of source genome
    '''

    logging.info("Read the chain file \"%s\" " % chain_file)
    maps = {}
    target_chromSize = {}
    source_chromSize = {}
    if print_table:
        blocks = []

    for line in ireader.reader(chain_file):
        # Example: chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
        if not line.strip():
            continue
        line = line.strip()
        if line.startswith(('#', ' ')):
            continue
        fields = line.split()

        if fields[0] == 'chain' and len(fields) in [12, 13]:
            # score = int(fields[1])		  # Alignment score
            source_name = fields[2]		  # E.g. chrY
            source_size = int(fields[3])  # Full length of the chromosome
            source_strand = fields[4]	  # Must be +
            if source_strand != '+':
                raise Exception(
                    "Source strand in a chain file must be +. (%s)" % line)
            source_start = int(fields[5])  # Start of source region
            # source_end = int(fields[6])	  # End of source region

            target_name = fields[7]		  # E.g. chr5
            target_size = int(fields[8])  # Full length of the chromosome
            target_strand = fields[9]	  # + or -
            target_start = int(fields[10])
            #target_end = int(fields[11])
            target_chromSize[target_name] = target_size
            source_chromSize[source_name] = source_size

            if target_strand not in ['+', '-']:
                raise Exception("Target strand must be - or +. (%s)" % line)
            #chain_id = None if len(fields) == 12 else fields[12]
            if source_name not in maps:
                maps[source_name] = Intersecter()

            sfrom, tfrom = source_start, target_start

        # Now read the alignment chain from the file and store it as a list (source_from, source_to) -> (target_from, target_to)
        elif fields[0] != 'chain' and len(fields) == 3:
            size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
            if print_table:
                if target_strand == '+':
                    blocks.append((source_name, sfrom, sfrom+size, source_strand,
                                  target_name, tfrom, tfrom+size, target_strand))
                elif target_strand == '-':
                    blocks.append((source_name, sfrom, sfrom+size, source_strand, target_name,
                                  target_size - (tfrom+size), target_size - tfrom, target_strand))

            if target_strand == '+':
                maps[source_name].add_interval(
                    Interval(sfrom, sfrom+size, (target_name, tfrom, tfrom+size, target_strand)))
            elif target_strand == '-':
                maps[source_name].add_interval(Interval(
                    sfrom, sfrom+size, (target_name, target_size - (tfrom+size), target_size - tfrom, target_strand)))

            sfrom += size + sgap
            tfrom += size + tgap

        elif fields[0] != 'chain' and len(fields) == 1:
            size = int(fields[0])
            if print_table:
                if target_strand == '+':
                    blocks.append((source_name, sfrom, sfrom+size, source_strand,
                                  target_name, tfrom, tfrom+size, target_strand))
                elif target_strand == '-':
                    blocks.append((source_name, sfrom, sfrom+size, source_strand, target_name,
                                  target_size - (tfrom+size), target_size - tfrom, target_strand))

            if target_strand == '+':
                maps[source_name].add_interval(
                    Interval(sfrom, sfrom+size, (target_name, tfrom, tfrom+size, target_strand)))
            elif target_strand == '-':
                maps[source_name].add_interval(Interval(
                    sfrom, sfrom+size, (target_name, target_size - (tfrom+size), target_size - tfrom, target_strand)))
        else:
            raise Exception("Invalid chain format. (%s)" % line)
    # if (sfrom + size) != source_end  or (tfrom + size) != target_end:
    #	 raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)

    if print_table:
        for i in blocks:
            print('\t'.join([str(n) for n in i]))

    return (maps, target_chromSize, source_chromSize)