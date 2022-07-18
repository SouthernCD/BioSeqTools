import argparse
from bioseqtools.pipelines import *


def args_init():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='SeqParser',
    )

    subparsers = parser.add_subparsers(
        title='subcommands', dest="subcommand_name")

    # argparse for FastaFormatClean
    parser_a = subparsers.add_parser('FastaFormatClean',
                                     help='clean fasta format to a good way\n',
                                     description='')

    parser_a.add_argument('raw_fasta_file', type=str, help='raw fasta file')
    parser_a.add_argument('output_fasta_file', type=str,
                          help='output fasta file')
    parser_a.add_argument("-t", "--seq_type", help="what type sequences you have? (default:prot)", default="prot",
                          choices=['prot', 'nucl'])
    parser_a.add_argument("-d", "--degenerate_sites_allowed",
                          help="allow degenerate sites", action='store_true')
    parser_a.add_argument("-g", "--gap_allowed",
                          help="allow gap sites", action='store_true')
    parser_a.add_argument("-s", "--translation_stop_allowed",
                          help="allow translation stop sites", action='store_true')
    parser_a.add_argument("-w", "--wrap", help="all lines of text be shorter than warp characters in length", type=int,
                          default=75)
    parser_a.add_argument("-r", "--replace_flag",
                          help="if find a unknown or not allowed site, should we replace to any site? (default: no and skip sequence)",
                          action='store_true')
    parser_a.add_argument("-n", "--full_name_flag",
                          help="use full sequence name in output", action='store_true')
    parser_a.add_argument('-l', "--log_file", type=str,
                          help='path for log file (default:None)', default=None)

    # argparse for FastaStats
    parser_a = subparsers.add_parser('FastaStats',
                                     help='Stats fasta file like N50 and so on\n',
                                     description='')

    parser_a.add_argument('fasta_file', type=str, help='fasta file')
    parser_a.add_argument('-o', "--output_file", type=str, help='output file')
    parser_a.add_argument(
        "-old", "--old_way", help="use old way from xuyuxing", action='store_true')

    # argparse for FastaToSQL
    parser_a = subparsers.add_parser('FastaToSQL',
                                     help='store a fasta file into a sqlite database\n',
                                     description='store a fasta file into a sqlite database')

    parser_a.add_argument('fasta_file', type=str, help='a fasta file')
    parser_a.add_argument('sqlite_file', type=str,
                          help='output sqlite database')
    parser_a.add_argument('-l', "--log_file", type=str,
                          help='path for log file (default:None)', default=None)
    parser_a.add_argument("-gzip", "--gzip_flag",
                          help="if fasta is gzip file", action='store_true')

    # argparse for FastaByID
    parser_a = subparsers.add_parser('FastaByID',
                                     help='Extract fasta records by IDs from a database fasta file\n',
                                     description='Extract fasta record by ID from a database fasta file')

    parser_a.add_argument('db_fasta_file', type=str,
                          help='a database fasta file')
    parser_a.add_argument('-i', "--ID_list", type=str,
                          help='ID list like ID1,ID2,ID3')
    parser_a.add_argument('-f', '--ID_file', type=str,
                          help='file with ID in a column')
    parser_a.add_argument('-o', "--output_file", type=str, help='output file')
    parser_a.add_argument('-l', "--log_file", type=str,
                          help='path for log file (default:None)', default=None)
    parser_a.add_argument(
        "-old", "--old_way", help="use old way from xuyuxing", action='store_true')
    parser_a.add_argument("-v", "--invert_flag",
                          help="use old way from xuyuxing", action='store_true')

    # argparse for FastaByGroup
    parser_a = subparsers.add_parser('FastaByGroup',
                                     help='Extract many subset fasta records by groups from a database fasta file\n',
                                     description='')

    parser_a.add_argument("-n", "--name_flag",
                          help="if group file frist column is not name of group, auto make group name",
                          action='store_false')
    parser_a.add_argument('group_file', type=str,
                          help='a csv file which different rows represent different groups and the first column should be the group name')
    parser_a.add_argument('db_fasta_file', type=str,
                          help='a database fasta file')
    parser_a.add_argument('out_dir', type=str, help='output dir')

    # argparse for FastaRename
    parser_a = subparsers.add_parser('FastaRename',
                                     help='Extract fasta records by IDs and rename them\n',
                                     description='Extract fasta records by IDs and rename them')

    parser_a.add_argument('db_fasta_file', type=str,
                          help='a database fasta file')
    parser_a.add_argument('-f', '--ID_map_file', type=str,
                          help='file with ID in first column and new ID in a second column')
    parser_a.add_argument('-o', "--output_file", type=str, help='output file')
    parser_a.add_argument('-p', "--seq_prefix", type=str,
                          help='rename seq by a given prefix')
    parser_a.add_argument(
        "-old", "--old_way", help="use old way from xuyuxing", action='store_true')

    # argparse for SeqMerge
    parser_a = subparsers.add_parser('SeqMerge',
                                     help='merge give sequences into a single sequence\n',
                                     description='merge give sequences into a single sequence')

    parser_a.add_argument('ID_merge_file', type=str,
                          help='file with ID in first column and merged new ID in a second column')
    parser_a.add_argument('db_fasta_file', type=str,
                          help='a database fasta file')
    parser_a.add_argument('output_file', type=str, help='output file')

    # argparse for GeneStructure
    parser_a = subparsers.add_parser('GeneStructure',
                                     help='Extract gene sequence by they structure\n')

    parser_a.add_argument('genome_fasta_file', type=str,
                          help='reference genome')
    parser_a.add_argument('gff_file', type=str,
                          help='reference genome gff file')
    parser_a.add_argument('gene_id', type=str, help='gene id which you want')
    parser_a.add_argument('-o', "--output_file", type=str, help='output file')

    # argparse for SeqLength
    parser_a = subparsers.add_parser('SeqLength',
                                     help='Get the length of each sequence in fasta file\n',
                                     description='Get the length of each sequence in fasta file')

    parser_a.add_argument('seq_file', type=str,
                          help='a fasta or fastq file as input')
    parser_a.add_argument(
        '-s', "--stat", help='make a histogram for seq length', action='store_true')
    parser_a.add_argument('-fq', "--fq_flag",
                          help='input is fastq file', action='store_true')
    parser_a.add_argument('-p', "--picture", type=str,
                          help='make a histogram picture for seq length picture file path need here', default="")
    parser_a.add_argument('-o', "--output_file", type=str, help='output file')

    # argparse for FilterByLen
    parser_a = subparsers.add_parser('FilterByLen',
                                     help='Filter a fasta file by length of sequence\n',
                                     description='Filter a fasta file by length of sequence')

    parser_a.add_argument('seq_file', type=str, help='a fasta as input')
    parser_a.add_argument('output_file', type=str, help='output file')
    parser_a.add_argument('condition_string', type=str,
                          help='condition string, such as: "seq>100", "seq<100", "seq>=100", "seq<=100", "seq<100 and seq>50", "seq<100 or seq>200"')

    # argparse for SubFasta
    parser_a = subparsers.add_parser('SubFasta',
                                     help='Extract local sequence by ID, start, end, and strand site from a database fasta file\n',
                                     description='Extract local sequence by ID, start, end, and strand site from a database fasta file')

    parser_a.add_argument('db_fasta_file', type=str,
                          help='a database fasta file')
    parser_a.add_argument('-s', "--search_string", type=str,
                          help='search_string like Chr1:1-100:+')
    parser_a.add_argument('-f', "--input_file", type=str,
                          help='csv file [subseqname,refseqname,start,end,strand]')
    parser_a.add_argument('-b', "--input_a_bed_file", action='store_true')
    parser_a.add_argument(
        "-r", "--RNA", help="change T to U", action='store_true')
    parser_a.add_argument('-o', "--output_file", type=str, help='output file')

    # argparse for GCproportion
    parser_a = subparsers.add_parser('GCproportion',
                                     help='get GC proportion for a fasta file\n')

    parser_a.add_argument('fasta_file', type=str, help='a database fasta file')
    parser_a.add_argument('output_file', type=str, help='output file')

    # argparse for ReComp
    parser_a = subparsers.add_parser('ReComp',
                                     help='reverse and complement a sequence\n')

    parser_a.add_argument('input_fasta_file', type=str, help='path for input')
    parser_a.add_argument('output_fasta_file', type=str,
                          help='path for output')
    parser_a.add_argument('-r', "--only_reverse", action='store_true')
    parser_a.add_argument('-c', "--only_complement", action='store_true')

    # argparse for SeqEntropy
    parser_a = subparsers.add_parser('SeqEntropy',
                                     help='get seq entropy value\n')

    parser_a.add_argument('input_fasta', type=str, help='path for input')
    parser_a.add_argument('output_file', type=str, help='path for output')
    parser_a.add_argument('-w', "--word_size", type=int, default=3)

    # argparse for SplitFasta
    parser_a = subparsers.add_parser('SplitFasta',
                                     help='Split a big fasta file into some smaller file\n',
                                     description='Split a big fasta file into some smaller file')

    parser_a.add_argument('fasta_file', type=str,
                          help='a big fasta file as input')
    parser_a.add_argument('-r', "--records", type=int, default=None,
                          help='put -r sequences per output fasta file')
    parser_a.add_argument('-c', "--contig_model",
                          help="use contig model, store one contig in one file and named with contig",
                          action='store_true')
    parser_a.add_argument('-s', "--size", type=int,
                          default=None, help='put size per output fasta file')
    parser_a.add_argument('-o', "--output_prefix", type=str, default="new_",
                          help='prefix for output, default is "new_"')

    # argparse for FragmentGenome
    parser_a = subparsers.add_parser('FragmentGenome',
                                     help='cut genome to fragment')
    parser_a.add_argument(
        "genome_file", help="Path of genome file with fasta format", type=str)
    parser_a.add_argument(
        "output_file", help="Path of cutted genome file", type=str)
    parser_a.add_argument("step", help="step of genome cutting", type=int)
    parser_a.add_argument(
        "length", help="length of sequence of cutted genome", type=int)
    parser_a.add_argument("-f", "--shift_start", help="fragment start at which site, (default as 0, meaning not shift)",
                          default=0, type=int)
    parser_a.add_argument("-s", "--Consider_scaffold", help="If consider scaffold in the genome (default as True)",
                          default=True)
    parser_a.add_argument("-e", "--entropy_threshold", help="filter simple repeat by seq entropy (default as 3)",
                          default=3, type=int)

    # argparse for ReGetContig
    parser_a = subparsers.add_parser('ReGetContig',
                                     help='cut scaffold to contig')
    parser_a.add_argument(
        "scaff_file", help="Path of scaffold version to contig again", type=str)
    parser_a.add_argument(
        "output_gff3", help="Path of cutted contig gff file", type=str)
    parser_a.add_argument(
        "output_fasta", help="Path of cutted contig fasta file", type=str)
    parser_a.add_argument(
        "-t", "--threads", help="threads number (default as 5)", default=5, type=int)

    # argparse for AlnStat
    parser_a = subparsers.add_parser('AlnStat', help='Alignment statistics tool\n',
                                     description='Alignment statistics tool')

    parser_a.add_argument('aln_file', type=str,
                          help='Alignments in fasta file')
    parser_a.add_argument('output_file', type=str, help='output file')
    parser_a.add_argument('-q', "--query", type=str,
                          help='query seq for marker record, default is marking in all colunm')

    # argparse for FastqSize
    parser_a = subparsers.add_parser('FastqSize',
                                     help='count base and reads number in a fastq file')

    parser_a.add_argument('-1', "--fastq_1", type=str, help='left end fastq')
    parser_a.add_argument('-2', "--fastq_2", type=str,
                          help='right end fastq (option)')
    parser_a.add_argument("-gzip", "--gzip_flag",
                          help="if fastq is gzip file", action='store_true')

    # argparse for PreFastqByID
    parser_a = subparsers.add_parser('PreFastqByID',
                                     help='build a sqlite3 database for fastq records\n',
                                     description='build a sqlite3 database for fastq records')

    parser_a.add_argument('-1', "--fastq_1", type=str, help='left end fastq')
    parser_a.add_argument('-2', "--fastq_2", type=str,
                          help='right end fastq (option)')
    parser_a.add_argument('-o', "--output_file", type=str,
                          help='name of output file')
    parser_a.add_argument("-gzip", "--gzip_flag",
                          help="if fastq is gzip file", action='store_true')

    # argparse for FastqByID
    parser_a = subparsers.add_parser('FastqByID',
                                     help='Extract fastq records by IDs from a database fastq file\n',
                                     description='Extract fastq record by ID from a database fastq file')

    parser_a.add_argument('db_file', type=str,
                          help='database file made by PreFastqByID')
    parser_a.add_argument('-i', "--ID_list", type=str,
                          help='ID list like ID1,ID2,ID3')
    parser_a.add_argument('-f', '--ID_file', type=str,
                          help='file with ID in a column')
    parser_a.add_argument('-o', "--output_prefix",
                          type=str, help='prefix of output file')

    # argparse for PreExtractNr
    parser_a = subparsers.add_parser('PreExtractNr',
                                     help='convert nr fasta file to a sqlite3 database for extract record\n',
                                     description='convert nr fasta file to a sqlite3 database for extract record')

    parser_a.add_argument('nr_fasta_file', type=str,
                          help='NCBI nr database in the fasta format')
    parser_a.add_argument('sql3_db_file', type=str,
                          help='sqlite3 which output by there command')

    # argparse for ExtractNr
    parser_a = subparsers.add_parser('ExtractNr',
                                     help='Extract local sequence by ID, start, end, and strand site from a database fasta file\n',
                                     description='Extract local sequence by ID, start, end, and strand site from a database fasta file')
    parser_a.add_argument('db_fasta_file', type=str,
                          help='a database fasta file')
    parser_a.add_argument('-i', "--ID_list", type=str,
                          help='ID list like ID1,ID2,ID3')
    parser_a.add_argument('-f', '--ID_file', type=str,
                          help='file with ID in a column')
    parser_a.add_argument("-l", "--long_name",
                          help="use full sequence name", action='store_false')
    parser_a.add_argument('-o', "--output_file", type=str, help='output file')

    args = parser.parse_args()

    return args


def main():

    args = args_init()
    args_dict = vars(args)

    # FastaByID
    if args_dict["subcommand_name"] == "FastaByID":
        FastaByID_main(args)

    # FastaByGroup
    elif args_dict["subcommand_name"] == "FastaByGroup":
        FastaByGroup_main(args)

    # subfasta
    elif args_dict["subcommand_name"] == "SubFasta":
        subfasta_main(args)

    # PreExtractNr
    elif args_dict["subcommand_name"] == "PreExtractNr":
        PreExtractNr_main(args)

    # ExtractNr
    elif args_dict["subcommand_name"] == "ExtractNr":
        ExtractNr_main(args)

    # SeqLength
    elif args_dict["subcommand_name"] == "SeqLength":
        SeqLength_main(args)

    # SplitFasta
    elif args_dict["subcommand_name"] == "SplitFasta":
        SplitFasta_main(args)

    # AlnStat
    elif args_dict["subcommand_name"] == "AlnStat":
        AlnStat_main(args)

    # FastqByID
    elif args_dict["subcommand_name"] == "FastqByID":
        FastqByID_main(args)

    # PreFastqByID
    elif args_dict["subcommand_name"] == "PreFastqByID":
        PreFastqByID_main(args)

    # FastaRename
    elif args_dict["subcommand_name"] == "FastaRename":
        FastaRename_main(args)

    # FilterByLen
    elif args_dict["subcommand_name"] == "FilterByLen":
        FilterByLen_main(args)

    # FastaToSQL
    elif args_dict["subcommand_name"] == "FastaToSQL":
        FastaToSQL_main(args)

    # FastaFormatClean
    elif args_dict["subcommand_name"] == "FastaFormatClean":
        FastaFormatClean_main(args)

    # ReComp
    elif args_dict["subcommand_name"] == "ReComp":
        ReComp_main(args)

    # FastqSize
    elif args_dict["subcommand_name"] == "FastqSize":
        FastqSize_main(args)

    elif args_dict["subcommand_name"] == "FastaStats":
        FastaStats_main(args)

    # ReGetContig
    elif args_dict["subcommand_name"] == "ReGetContig":
        ReGetContig_main(args)

    # GeneStructure
    elif args_dict["subcommand_name"] == "GeneStructure":
        GeneStructure_main(args)

    # SeqMerge
    elif args_dict["subcommand_name"] == "SeqMerge":
        SeqMerge_main(args)

    # GCproportion
    elif args_dict["subcommand_name"] == "GCproportion":
        GCproportion_main(args)

    # SeqEntropy
    elif args_dict["subcommand_name"] == "SeqEntropy":
        SeqEntropy_main(args)

    # FragmentGenome
    elif args_dict["subcommand_name"] == "FragmentGenome":
        FragmentGenome_main(args)

if __name__ == '__main__':
    main()
