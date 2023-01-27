from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from collections import Counter
from toolbiox.lib.common.os import mkdir, multiprocess_running
from toolbiox.lib.xuyuxing.base.common_command import log_print
from toolbiox.lib.common.util import logging_init
from toolbiox.api.xuyuxing.file_parser.fileIO import tsv_file_parse, tsv_file_dict_parse
from toolbiox.lib.common.math.interval import merge_intervals
from toolbiox.lib.xuyuxing.math.set_operating import uniqify
from toolbiox.api.xuyuxing.plot.histogram import int_barplot
from toolbiox.lib.common.genome.seq_base import read_fasta_by_faidx, BioSeq
from pyfaidx import Fasta
import toolbiox.lib.common.sqlite_command as sc
import toolbiox.lib.common.genome.seq_base as sb
import re
import sqlite3
import sys
import time
import os


def FastaByID(ID_list, fasta_file, output_file, log_file, invert_flag=False):
    logger = logging_init("FastaByID", log_file)

    logger.info("Step1: parsing query ID list")
    if output_file is not None:
        F1 = open(output_file, 'w')
    ID_list_short = [re.search('^(\S+)', i).group(1) for i in ID_list]
    logger.info("Step1: finished")

    logger.info("Step2: loading fasta file")
    fasta_dict = Fasta(fasta_file)
    logger.info("Step2: finished")

    num = 0

    if invert_flag is False:
        for query_id in ID_list_short:
            if query_id in fasta_dict:
                num = num + 1
                seq_record = fasta_dict[query_id]
                full_name = seq_record.long_name
                # print(full_name)
                if output_file is not None:
                    F1.write(">%s\n" % full_name)
                    for line in seq_record:
                        F1.write(str(line) + "\n")
                else:
                    print(">%s" % full_name)
                    for line in seq_record:
                        print(line)
    elif invert_flag is True:
        ID_list_short = set(ID_list_short)
        for query_id in fasta_dict.keys():
            if query_id not in ID_list_short:
                num = num + 1
                seq_record = fasta_dict[query_id]
                full_name = seq_record.long_name
                if output_file is not None:
                    F1.write(">%s\n" % full_name)
                    for line in seq_record:
                        F1.write(str(line) + "\n")
                else:
                    print(">%s" % full_name)
                    for line in seq_record:
                        print(line)

    if output_file is not None:
        F1.close()
    logger.info("Found query sequences")
    logger.critical("%d records found in %d queries" %
                    (num, len(ID_list_short)))
    del logger.handlers[:]


def FastaByID_xyx(ID_list, db_fasta_file, output_file, log_file):
    logger = logging_init("FastaByID", log_file)
    logger.info("Loading fasta file")
    if output_file is not None:
        F1 = open(output_file, 'w')
    ID_list_short = [re.search('^(\S+)', i).group(1) for i in ID_list]
    output_dict = []
    for record in sb.read_fasta_big(db_fasta_file, upper=False, log_file=log_file):
        if record.seqname_short() in ID_list_short or record.seqname in ID_list:
            output_dict.append(record)
    output_dict = uniqify(output_dict)
    logger.info("Loaded fasta file")
    logger.info("Finding query sequences")
    for record in output_dict:
        record.wrap()
        if output_file is not None:
            F1.write(">%s\n%s\n" % (record.seqname_short(), record.seq))
        else:
            print(">%s\n%s" % (record.seqname_short(), record.seq))
    if output_file is not None:
        F1.close()
    logger.info("Found query sequences")
    logger.critical("%d records found in %d queries" %
                    (len(output_dict), len(ID_list)))


def FastaByGroup(group_file, db_fasta_file, out_dir, group_name_flag=True):
    if group_name_flag is True:
        group_list = tsv_file_parse(group_file, key_col=1, seq=",")
    else:
        group_list = tsv_file_parse(group_file, seq=",", prefix="Group_")
    group_list_hash = {}
    for group_name in group_list:
        for record_name in group_list[group_name]:
            if record_name in group_list_hash:
                group_list_hash[record_name].append(group_name)
            else:
                group_list_hash[record_name] = [group_name]
            group_list_hash[record_name] = list(
                set(group_list_hash[record_name]))

    group_list_record = {}
    for record in sb.read_fasta_big(db_fasta_file, upper=False):
        record_name = record.seqname_short()
        if record_name in group_list_hash:
            for group_name in group_list_hash[record_name]:
                if group_name in group_list_record:
                    group_list_record[group_name].append(record)
                else:
                    group_list_record[group_name] = [record]
                group_list_record[group_name] = uniqify(
                    group_list_record[group_name])

    mkdir(out_dir)

    for group_name in group_list_record:
        with open(out_dir + "/" + group_name + ".fasta", 'w') as f:
            for record in group_list_record[group_name]:
                f.write(">%s\n%s\n" % (record.seqname_short(), record.seq))


def subfasta(db_fasta_file, search_string, input_file, RNA_flag, output_file, bed_file):
    out_dir = {}
    if search_string is not None:
        refseq, start, end, strand = re.search(
            r'^(\S+):(\d+)-(\d+):(\S)$', search_string).groups()
        out = sb.sub_fasta(db_fasta_file, refseq, int(
            start), int(end), strand, RNA_flag)
        out_dir[search_string] = out

    if input_file is not None:
        if not bed_file:
            sub_table = tsv_file_parse(
                input_file, key_col=1, fields="2,3,4,5", seq=",")
        else:
            sub_table = tsv_file_parse(
                input_file, key_col=4, fields="1,2,3,6", seq="\t")
        print(sub_table)
        if RNA_flag:
            for i in sub_table.keys():
                sub_table[i] = sub_table[i] + (True,)
            out_dir = sb.sub_fasta_many(db_fasta_file, sub_table)
        else:
            for i in sub_table.keys():
                sub_table[i] = sub_table[i] + (False,)
            out_dir = sb.sub_fasta_many(db_fasta_file, sub_table)

    if output_file is not None:
        F1 = open(output_file, 'w')

    for record_name in out_dir:
        record = out_dir[record_name]
        record.wrap()
        if output_file is not None:
            F1.write(">%s\n%s" % (record.seqname_short(), record.seq))
        else:
            print(">%s\n%s" % (record.seqname_short(), record.seq))
    if output_file is not None:
        F1.close()


def split_fasta(input_file, output_dir, outpre, contig_model=False, unit_num=None, file_size=None):
    mkdir(output_dir, True)

    record_dict = sb.read_fasta_by_faidx(input_file)
    if contig_model:
        for i in record_dict:
            record_tmp = record_dict[i]
            with open(output_dir + "/" + record_tmp.seqname_short() + ".fa", 'w') as f:
                f.write(">%s\n%s" % (record_tmp.seqname, record_tmp.seq))
    elif unit_num:
        all_used = 0
        num_used = 0
        file_num = 0
        F1 = open(output_dir + "/%s%d.fa" % (outpre, file_num), 'w')

        for record_id in record_dict:
            record = record_dict[record_id]
            num_used = num_used + 1
            F1.write(">%s\n%s\n" % (record.seqname, record.seq))
            all_used += 1
            if num_used >= unit_num:
                F1.close()
                file_num = file_num + 1
                if all_used < len(record_dict):
                    F1 = open(output_dir + "/%s%d.fa" %
                              (outpre, file_num), 'w')
                num_used = 0

        try:
            F1.close()
        except:
            pass

    elif file_size:
        sum_size = 0
        file_num = 0
        F1 = open(output_dir + "/%s%d.fa" % (outpre, file_num), 'w')

        for record_id in record_dict:
            record = record_dict[record_id]
            if sum_size == 0 and record.len() >= file_size:
                new_file_flag = True
                add_to_new_file = False
                F1.write(">%s\n%s\n" % (record.seqname, record.seq))
                sum_size += record.len()
                # print(sum_size)
            elif sum_size + record.len() >= file_size:
                new_file_flag = True
                add_to_new_file = True
                sum_size += record.len()
                # print(sum_size)
            else:
                new_file_flag = False
                add_to_new_file = False
                F1.write(">%s\n%s\n" % (record.seqname, record.seq))
                sum_size += record.len()

            if new_file_flag:
                F1.close()
                file_num = file_num + 1
                F1 = open(output_dir + "/%s%d.fa" % (outpre, file_num), 'w')
                sum_size = 0
                if add_to_new_file:
                    F1.write(">%s\n%s\n" % (record.seqname, record.seq))
                    sum_size += record.len()

        F1.close()


def FastqByID(ID_list, db_file, output_prefix):
    def output_result(db_file, ID_list_tmp, paired_flag, output_prefix, frist_write):
        if paired_flag:
            record_tuple_list_f = sc.sqlite_select_by_a_key(db_file, "record", "fseqname_short",
                                                            tuple(ID_list_tmp))
            record_tuple_list_r = sc.sqlite_select_by_a_key(db_file, "record", "rseqname_short",
                                                            tuple(ID_list_tmp))
            record_tuple_list = list(
                set(record_tuple_list_f + record_tuple_list_r))
        else:
            record_tuple_list = sc.sqlite_select_by_a_key(db_file, "record", "fseqname_short",
                                                          tuple(ID_list_tmp))
            record_tuple_list = list(set(record_tuple_list))

        if output_prefix is not None:
            if paired_flag:
                if frist_write:
                    F1 = open(output_prefix + "_1.fq", 'w')
                    F2 = open(output_prefix + "_2.fq", 'w')
                else:
                    F1 = open(output_prefix + "_1.fq", 'a')
                    F2 = open(output_prefix + "_2.fq", 'a')
            else:
                if frist_write:
                    F1 = open(output_prefix + ".fq", 'w')
                else:
                    F1 = open(output_prefix + ".fq", 'a')

        for record in record_tuple_list:
            if output_prefix is not None:
                if paired_flag:
                    record1_seqname_short, record1_seqname, record1_seqs, record1_quality, record2_seqname_short, record2_seqname, record2_seqs, record2_quality = record
                    F1.write(
                        "@%s\n%s\n+%s\n%s\n" % (
                            record1_seqname, record1_seqs, record1_seqname, record1_quality))
                    F2.write(
                        "@%s\n%s\n+%s\n%s\n" % (
                            record2_seqname, record2_seqs, record2_seqname, record2_quality))
                else:
                    record1_seqname_short, record1_seqname, record1_seqs, record1_quality = record
                    F1.write(
                        "@%s\n%s\n+%s\n%s\n" % (
                            record1_seqname, record1_seqs, record1_seqname, record1_quality))
            else:
                if paired_flag:
                    record1_seqname_short, record1_seqname, record1_seqs, record1_quality, record2_seqname_short, record2_seqname, record2_seqs, record2_quality = record
                    print("F: @%s\nF: %s\nF: +%s\nF: %s\n" % (
                        record1_seqname, record1_seqs, record1_seqname, record1_quality))
                    print("R: @%s\nR: %s\nR: +%s\nR: %s\n" % (
                        record2_seqname, record2_seqs, record2_seqname, record2_quality))
                else:
                    record1_seqname_short, record1_seqname, record1_seqs, record1_quality = record
                    print("@%s\n%s\n+%s\n%s\n" % (
                        record1_seqname, record1_seqs, record1_seqname, record1_quality))

        if output_prefix is not None:
            if paired_flag:
                F1.close()
                F2.close()
            else:
                F1.close()

        return len(record_tuple_list)

    log_print("Begin: Extract fastq files from sqlite3 db")
    start_time = time.time()

    conn = sqlite3.connect(db_file)
    content = conn.execute("SELECT * FROM info").fetchall()[0]
    conn.close()

    if content[0] == 'paired_end':
        paired_flag = True
    else:
        paired_flag = False

    num = 0
    frist_write = True
    ID_list_tmp = []
    output_num = 0
    for ID in ID_list:
        ID_list_tmp.append(ID)
        num = num + 1
        if num % 10000 == 0:
            output_num = output_num + output_result(db_file, ID_list_tmp, paired_flag, output_prefix,
                                                    frist_write)
            frist_write = False
            ID_list_tmp = []

        round_time = time.time()
        if round_time - start_time > 10:
            log_print("\tparsed: %d" % (num))
            start_time = round_time

    if len(ID_list_tmp) > 0:
        output_num = output_num + \
            output_result(db_file, ID_list_tmp, paired_flag,
                          output_prefix, frist_write)
        ID_list_tmp = []
        log_print("\tparsed: %d" % (num))

    log_print("End: Extract fastq files from sqlite3 db")

    return output_num, len(ID_list)


def store_fastq_record_into_sqlite(fq1, fq2, output_file, gzip_flag):
    log_print("Begin: convert fastq files to sqlite3 db")
    start_time = time.time()

    if fq2 is None:
        table_columns_dict = {
            "info": ["end_type", "ffq_file_name"],
            "record": ["fseqname_short", "fseqname", "fseqs", "fquality"]
        }
        sc.init_sql_db_many_table(output_file, table_columns_dict)
        sc.insert_one_record_to_sql_table(
            ("single_end", fq1), table_columns_dict["info"], output_file, "info")
    else:
        table_columns_dict = {
            "info": ["end_type", "ffq_file_name", "rfq_file_name"],
            "record": ["fseqname_short", "fseqname", "fseqs", "fquality", "rseqname_short", "rseqname", "rseqs",
                       "rquality"]
        }
        sc.init_sql_db_many_table(output_file, table_columns_dict)
        sc.insert_one_record_to_sql_table(("paired_end", fq1, fq2), table_columns_dict["info"], output_file,
                                          "info")

    num = 0
    record_tmp_dict = []
    for records in sb.read_fastq_big(fq1, fq2, gzip_flag):
        if fq2 is not None:
            record1 = records[0]
            record2 = records[1]
            record_tmp_dict.append(
                (
                    record1.seqname_short(), record1.seqname, record1.seq, record1.quality,
                    record2.seqname_short(),
                    record2.seqname, record2.seq, record2.quality))
        else:
            record1 = records[0]
            record_tmp_dict.append(
                (record1.seqname_short(), record1.seqname, record1.seq, record1.quality))

        num = num + 1

        if num % 10000 == 0:
            sc.sqlite_write(record_tmp_dict, output_file,
                            "record", table_columns_dict["record"])
            record_tmp_dict = []

        round_time = time.time()
        if round_time - start_time > 10:
            log_print("\tparsed: %d" % (num))
            start_time = round_time

    if len(record_tmp_dict) > 0:
        sc.sqlite_write(record_tmp_dict, output_file,
                        "record", table_columns_dict["record"])
        record_tmp_dict = []
        log_print("\tparsed: %d" % (num))

    conn = sqlite3.connect(output_file)
    if fq2 is not None:
        conn.execute(
            "CREATE UNIQUE INDEX fseqname_index on record (fseqname_short)")
        conn.execute(
            "CREATE UNIQUE INDEX rseqname_index on record (rseqname_short)")
    else:
        conn.execute(
            "CREATE UNIQUE INDEX fseqname_index on record (fseqname_short)")
    conn.close()

    log_print("End: convert fastq files to sqlite3 db")


def FastaRename(ID_dict, db_fasta_file, output_file):
    """
    ID_dict =
    db_fasta_file =
    output_file:
    """

    if output_file is not None:
        F1 = open(output_file, 'w')

    output_dict = []
    for record in sb.read_fasta_big(db_fasta_file, upper=False):
        if record.seqname_short() in ID_dict:
            record.newname = ID_dict[record.seqname_short()][0]
            output_dict.append(record)
        elif record.seqname in ID_dict:
            record.newname = ID_dict[record.seqname][0]
            output_dict.append(record)

    output_dict = uniqify(output_dict)

    for record in output_dict:
        record.wrap()
        if output_file is not None:
            F1.write(">%s\n%s" % (record.newname, record.seq))
        else:
            print(">%s\n%s" % (record.newname, record.seq))
    if output_file is not None:
        F1.close()

    return len(output_dict), len(ID_dict)


def FastaRename_faidx(ID_dict, db_fasta_file, output_file):
    # if output_file is not None:
    #     F1 = open(output_file, 'w')

    db_file_dict = sb.read_fasta_by_faidx(db_fasta_file)
    num = 0
    for old_name in ID_dict:
        new_name = ID_dict[old_name][0]
        if old_name in db_file_dict:
            num = num + 1
            record = db_file_dict[old_name]
            record.seqname = new_name
            # record.wrap()
            # print("wrap done")
            if output_file is not None:
                record.write_to_file(output_file)
                # F1.write(">%s\n%s" % (new_name, record.seq))
            else:
                print(">%s\n%s" % (new_name, record.seq))
        else:
            print('failed to find %s' % old_name)
    #
    # if output_file is not None:
    #     F1.close()

    return num, len(ID_dict)


# main

# FastaByID
def FastaByID_main(args):

    db_fasta_file = args.db_fasta_file
    output_file = args.output_file
    log_file = args.log_file
    old_way = args.old_way

    ID_file = args.ID_file
    ID_list = args.ID_list
    if ID_file:
        ID_list = tsv_file_parse(ID_file, key_col=1)
        ID_list = ID_list.keys()
    elif ID_list:
        ID_list = ID_list.split(",")

    if not old_way:
        FastaByID(ID_list, db_fasta_file, output_file,
                  log_file, args.invert_flag)
    else:
        FastaByID_xyx(ID_list, db_fasta_file, output_file, log_file)

# FastaByGroup


def FastaByGroup_main(args):
    group_file = args.group_file
    db_fasta_file = args.db_fasta_file
    out_dir = args.out_dir
    group_name_flag = args.name_flag
    FastaByGroup(group_file, db_fasta_file, out_dir, group_name_flag)

# subfasta


def subfasta_main(args):
    db_fasta_file = args.db_fasta_file
    search_string = args.search_string
    input_file = args.input_file
    RNA_flag = args.RNA
    output_file = args.output_file

    if input_file is None and search_string is None:
        raise ValueError("input_file or search_string need at least one")

    subfasta(db_fasta_file, search_string, input_file,
             RNA_flag, output_file, args.input_a_bed_file)

# PreExtractNr


def PreExtractNr_main(args):
    nr_fasta_file = args.nr_fasta_file
    sql3_db_file = args.sql3_db_file
    sb.PreExtractNr(nr_fasta_file, sql3_db_file)

# ExtractNr


def ExtractNr_main(args):

    db_fasta_file = args.db_fasta_file
    output_file = args.output_file
    long_name = args.long_name

    ID_file = args.ID_file
    ID_list = args.ID_list
    if ID_file:
        ID_list = tsv_file_parse(ID_file, key_col=1)
        ID_list = ID_list.keys()
    elif ID_list:
        ID_list = ID_list.split(",")

    found_num, want_num = sb.ExtractNr(db_fasta_file, tuple(
        ID_list), output_file, short_name=long_name)
    print("%d records found in %d queries" % (found_num, want_num))


# SeqLength


def SeqLength_main(args):
    seq_file = args.seq_file
    output_file = args.output_file
    stat = args.stat
    picture = args.picture
    fq_flag = args.fq_flag

    if output_file is not None:
        F1 = open(output_file, 'w')

    seq_len_list = []
    if not fq_flag:
        seq_file_gen = sb.read_fasta_big(seq_file, upper=False)
    else:
        seq_file_gen = sb.read_fastq_big(seq_file)

    for record in seq_file_gen:
        if fq_flag:
            record = record[0]
        if stat is False:
            if output_file is not None:
                F1.write("%s\t%d\n" %
                         (record.seqname_short(), len(record.seq)))
            else:
                print("%s\t%s" % (record.seqname_short(), len(record.seq)))
        else:
            seq_len_list.append(len(record.seq))

    if stat:
        if not picture == "":
            int_barplot(seq_len_list, picture)
        else:
            from collections import Counter

            Counter_data = Counter(seq_len_list)
            for i in range(min(Counter_data.keys()), max(Counter_data.keys()) + 1):
                if output_file is not None:
                    F1.write("%d\t%d\n" % (i, Counter_data[i]))
                else:
                    print("%d\t%d" % (i, Counter_data[i]))

    if output_file is not None:
        F1.close()

# SplitFasta


def SplitFasta_main(args):
    # split_fasta(input_file, unit_num, output_dir, outpre, args.contig_model)

    split_fasta(args.fasta_file, os.getcwd(), args.output_prefix, contig_model=args.contig_model,
                unit_num=args.records, file_size=args.size)

# AlnStat


def AlnStat_main(args):
    aln_file = args.aln_file
    query = args.query
    output_file = args.output_file

    record_dict, seqname_list = sb.read_fasta(aln_file)
    col_num = len(record_dict[seqname_list[0]].seq)
    seq_num = len(seqname_list)
    count_dict = {}
    for i in range(0, col_num):
        count_dict[i] = dict(Counter([record_dict[seqname].seq[i]
                                      for seqname in record_dict]))

    if query is None:
        with open(output_file, "w") as f:
            for site in count_dict:
                highest = sorted(
                    count_dict[site], key=lambda x: count_dict[site][x], reverse=True)[0]
                if not highest == "-":
                    ratio_tmp = float(
                        count_dict[site][highest]) / float(seq_num)
                else:
                    ratio_tmp = 0.0
                f.write("%d\t%s\n" % (site, ratio_tmp))
    else:
        with open(output_file, "w") as f:
            num = 0
            for site in count_dict:
                if record_dict[query].seq[site] == "-":
                    continue
                num = num + 1
                highest = sorted(
                    count_dict[site], key=lambda x: count_dict[site][x], reverse=True)[0]
                if not highest == "-":
                    ratio_tmp = float(
                        count_dict[site][highest]) / float(seq_num)
                else:
                    ratio_tmp = 0.0
                f.write("%d\t%s\n" % (num, ratio_tmp))

# FastqByID


def FastqByID_main(args):
    db_file = args.db_file
    ID_list = args.ID_list
    ID_file = args.ID_file
    output_prefix = args.output_prefix

    if ID_file:
        ID_list = tsv_file_parse(ID_file, key_col=1)
        ID_list = ID_list.keys()
    elif ID_list:
        ID_list = ID_list.split(",")

    ID_list_used = []
    for i in ID_list:
        ID_list_used.append(re.search('^(\S+)', i).group(1))

    found_num, want_num = FastqByID(ID_list_used, db_file, output_prefix)
    print("%d records found in %d queries" % (found_num, want_num))

# PreFastqByID


def PreFastqByID_main(args):
    fastq_1 = args.fastq_1
    fastq_2 = args.fastq_2
    output_file = args.output_file
    gzip_flag = args.gzip_flag

    store_fastq_record_into_sqlite(fastq_1, fastq_2, output_file, gzip_flag)

# FastaRename


def FastaRename_main(args):
    ID_map_file = args.ID_map_file
    db_fasta_file = args.db_fasta_file
    output_file = args.output_file

    if args.seq_prefix is not None:
        num = 0
        rename_map = []
        with open(args.output_file, 'w') as f:
            for record in sb.read_fasta_big(db_fasta_file, upper=False):
                record.wrap()
                new_name = args.seq_prefix + str(num)
                f.write(">%s\n%s" % (new_name, record.seq))
                rename_map.append((record.seqname, new_name))
                num = num + 1
        with open(args.output_file + ".rename.map", 'w') as f:
            for i in rename_map:
                f.write("%s\t%s\n" % (i[1], i[0]))

    else:
        # ID_map_file = "/lustre/home/xuyuxing/Work/miRNA/Cau_two_assembly/remove/remove/Cucumis_Cau/Cau+Cucumis.del.map"
        ID_dict = tsv_file_parse(ID_map_file, key_col=1, fields="2")
        if args.old_way:
            found_num, want_num = FastaRename(
                ID_dict, db_fasta_file, output_file)
        else:
            found_num, want_num = FastaRename_faidx(
                ID_dict, db_fasta_file, output_file)
        print("%d records found in %d queries" % (found_num, want_num))

# FilterByLen


def FilterByLen_main(args):
    seq_file = args.seq_file
    output_file = args.output_file
    condition_string = args.condition_string

    with open(output_file, "w") as f:
        for record in sb.read_fasta_big(seq_file, upper=False):
            seq = record.seqs_length()
            record.wrap()
            # print(seq,condition_string,eval(condition_string))
            if eval(condition_string):
                f.write(">%s\n%s\n" % (record.seqname, record.seq))

# FastaToSQL


def FastaToSQL_main(args):
    fasta_file = args.fasta_file
    sqlite_file = args.sqlite_file
    log_file = args.log_file
    gzip_flag = args.gzip_flag

    sb.read_fasta_to_sqlite(fasta_file, sqlite_file, gzip_flag, log_file)

# FastaFormatClean


def FastaFormatClean_main(args):
    raw_fasta_file = args.raw_fasta_file
    output_fasta_file = args.output_fasta_file

    seq_type = args.seq_type
    if seq_type == 'nucl':
        seq_type = 0b0001
    elif seq_type == 'prot':
        seq_type = 0b0000

    dsa_flag = args.degenerate_sites_allowed
    if dsa_flag:
        dsa_flag = 0b0010

    ga_flag = args.gap_allowed
    if ga_flag:
        ga_flag = 0b0100

    tsa_flag = args.translation_stop_allowed
    if tsa_flag:
        tsa_flag = 0b1000

    wrap = args.wrap
    replace_flag = args.replace_flag
    full_name_flag = args.full_name_flag
    log_file = args.log_file

    bit_field = seq_type + dsa_flag + ga_flag + tsa_flag

    logger = logging_init("FastaFormatClean", log_file)

    logger.info("Start: parsing fasta file")

    if bit_field == 0b0111:
        all_allowed = sb.FastaRecord.nucl + sb.FastaRecord.any_nucl + \
            sb.FastaRecord.nucl_degenerate + sb.FastaRecord.gap_site
        safe_remove = []
        safe_replace = sb.FastaRecord.any_nucl
    elif bit_field == 0b0011:
        all_allowed = sb.FastaRecord.nucl + \
            sb.FastaRecord.any_nucl + sb.FastaRecord.nucl_degenerate
        safe_remove = sb.FastaRecord.gap_site
        safe_replace = sb.FastaRecord.any_nucl
    elif bit_field == 0b0101:
        all_allowed = sb.FastaRecord.nucl + \
            sb.FastaRecord.any_nucl + sb.FastaRecord.gap_site
        safe_remove = []
        safe_replace = sb.FastaRecord.any_nucl
    elif bit_field == 0b0001:
        all_allowed = sb.FastaRecord.nucl + sb.FastaRecord.any_nucl
        safe_remove = sb.FastaRecord.gap_site
        safe_replace = sb.FastaRecord.any_nucl
    elif bit_field == 0b0000:
        all_allowed = sb.FastaRecord.prot + sb.FastaRecord.any_prot
        safe_remove = sb.FastaRecord.gap_site + sb.FastaRecord.stop_site
        safe_replace = sb.FastaRecord.any_prot
    elif bit_field == 0b0100:
        all_allowed = sb.FastaRecord.prot + \
            sb.FastaRecord.any_prot + sb.FastaRecord.gap_site
        safe_remove = sb.FastaRecord.stop_site
        safe_replace = sb.FastaRecord.any_prot
    elif bit_field == 0b1100:
        all_allowed = sb.FastaRecord.prot + sb.FastaRecord.any_prot + \
            sb.FastaRecord.gap_site + sb.FastaRecord.stop_site
        safe_remove = []
        safe_replace = sb.FastaRecord.any_prot
    elif bit_field == 0b1000:
        all_allowed = sb.FastaRecord.prot + \
            sb.FastaRecord.any_prot + sb.FastaRecord.stop_site
        safe_remove = sb.FastaRecord.gap_site
        safe_replace = sb.FastaRecord.any_prot
    elif bit_field == 0b0010:
        all_allowed = sb.FastaRecord.prot + \
            sb.FastaRecord.any_prot
        safe_remove = sb.FastaRecord.gap_site
        safe_replace = sb.FastaRecord.stop_site
    else:
        raise ValueError("Bad input flag")

    with open(output_fasta_file, 'w') as f:
        for record in sb.read_fasta_big(raw_fasta_file, upper=True, log_file=log_file):
            logger.info("parsing %s with code %d" %
                        (record.seqname_short(), bit_field))
            record.no_wrap()
            new_seq = record.seq
            char_site = set(new_seq)
            bad_flag = False
            for i in char_site:
                if i in all_allowed:
                    continue
                elif i in safe_remove:
                    if i == "*":
                        new_seq = re.sub("\*", "", new_seq)
                    else:
                        new_seq = re.sub(i, "", new_seq)
                else:
                    if replace_flag:
                        if i == "*":
                            new_seq = re.sub("\*", "", new_seq)
                        else:
                            new_seq = re.sub(i, safe_replace[0], new_seq)
                    else:
                        logger.error("%s have bad site %s" %
                                     (record.seqname_short(), i))
                        bad_flag = True

            if bad_flag is True:
                continue

            record.seq = new_seq
            record.wrap(wrap)

            if full_name_flag:
                new_name = record.seqname
            else:
                new_name = record.seqname_short()

            f.write(">%s\n%s" % (new_name, record.seq))

# ReComp


def ReComp_main(args):
    with open(args.output_fasta_file, 'w') as f:
        for record in sb.read_fasta_big(args.input_fasta_file, upper=True):
            record.no_wrap()
            new_seq = record.seq

            if args.only_reverse:
                new_seq = sb.reverse_seq(new_seq)
                new_name = record.seqname + "_reverse"
            elif args.only_complement:
                new_seq = sb.complement_seq(new_seq)
                new_name = record.seqname + "_complement"
            else:
                new_seq = sb.reverse_complement(new_seq)
                new_name = record.seqname + "_reverse_complement"
            record.seq = new_seq
            record.wrap(60)

            f.write(">%s\n%s" % (new_name, record.seq))

# FastqSize


def FastqSize_main(args):
    count_read1 = 0
    count_read1_base = 0
    count_read2 = 0
    count_read2_base = 0
    for records in sb.read_fastq_big(args.fastq_1, args.fastq_2, args.gzip_flag):
        if args.fastq_2 is not None:
            record1 = records[0]
            record2 = records[1]
            count_read1 = count_read1 + 1
            count_read2 = count_read2 + 1
            count_read1_base = count_read1_base + len(record1.seq)
            count_read2_base = count_read2_base + len(record2.seq)
        else:
            record1 = records[0]
            count_read1 = count_read1 + 1
            count_read1_base = count_read1_base + len(record1.seq)

    if args.fastq_2 is not None:
        print("The fastq is paired:\n%s\t%s\n\n\tfastq_1\tfastq_2\tsum" %
              (args.fastq_1, args.fastq_2))
        print("base_num\t%d\t%d\t%d" % (count_read1_base,
                                        count_read2_base, count_read1_base + count_read2_base))
        print("count_num\t%d\t%d\t%d" %
              (count_read1, count_read2, count_read1 + count_read2))
    else:
        print("The fastq is single:\n%s\n\n\tfastq" % args.fastq_1)
        print("base_num\t%d" % (count_read1_base))
        print("count_num\t%d" % (count_read1))

# FastaStats


def FastaStats_main(args):
    """
    class abc():
        pass

    args = abc()
    args.old_way = False
    args.fasta_file = '/lustre/home/xuyuxing/Database/Cuscuta/Cau/genomev1.1/Cuscuta.genome.v1.1.fasta'


    """
    from toolbiox.lib.xuyuxing.math.stats import n50_stats
    import toolbiox.lib.common.genome.seq_base as sb
    from pyfaidx import Fasta

    if args.old_way is False:
        len_list = (len(i) for i in Fasta(args.fasta_file))
    else:
        len_list = (
            record.seqs_length for record in sb.read_fasta_big(args.fasta_file))

    percent_stats, index_stats_dict, count_sum, sum_len, min_len, max_len = n50_stats(
        len_list)

    output_string = "%s\n%s\t%s\t%s\t%s\n" % (
        args.fasta_file, "Nx", "Size", "Count_Sum", "Base_Sum")
    for i in index_stats_dict:
        output_string = output_string + "N%d\t%d\t%d\t%d\n" % (
            i, index_stats_dict[i]['length'], index_stats_dict[i]['count'], index_stats_dict[i]['base'])

    output_string = output_string + "\nN-perc\tSize\n"

    args.index_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for i in range(0, len(percent_stats)):
        output_string = output_string + \
            "%d\t%.2f\n" % (args.index_list[i], percent_stats[i])

    output_string = output_string + \
        "\nTotal Length/Count\t%d\t%d\n" % (sum_len, count_sum)
    output_string = output_string + "Maximum Length\t%d\n" % max_len
    output_string = output_string + "Minimum Length\t%d\n" % min_len

    print(output_string)

    if not args.output_file is None:
        with open(args.output_file, 'w') as f:
            f.write(output_string)

# FragmentGenome


def FragmentGenome_main(args):
    genome_file = args.genome_file
    output_file = args.output_file
    step = args.step
    length = args.length
    Consider_scaffold = args.Consider_scaffold
    shift_start = args.shift_start

    with open(output_file, 'w') as f:
        for i in sb.read_fasta_big(genome_file):
            record_tmp = i
            seq_name = record_tmp.seqname_short()
            sequence = record_tmp.seq
            # print(shift_start)
            for start, sub_seq in sb.scaffold_cutting(sequence, step=step, length=length,
                                                      Consider_scaffold=Consider_scaffold, shift_start=shift_start):
                # print(start)
                end = len(sub_seq) + start - 1
                sub_seq_name = seq_name + "_" + str(start) + "-" + str(end)

                if sb.sequence_entropy(sub_seq, 3) > args.entropy_threshold:

                    f.write(">%s\n%s\n" % (sub_seq_name, sub_seq))

# ReGetContig


def contig_range(SeqRecord_input):
    sequence = str(SeqRecord_input.seq)

    num = 1
    n_site = []
    for i in sequence:
        if i != "N":
            n_site.append((num, num))
        num = num + 1
    contig_range_output = merge_intervals(n_site, True)

    SeqRecord_input.features = []

    num = 0
    for contig_range_tmp in contig_range_output:
        qualifiers = {
            "source": "ReGetContig",
            "ID": "%s_%d" % (SeqRecord_input.id, num)
        }
        top_feature = SeqFeature(FeatureLocation(contig_range_tmp[0] - 1, contig_range_tmp[1]),
                                 type="contig", qualifiers=qualifiers)
        SeqRecord_input.features.append(top_feature)
        num = num + 1

    return SeqRecord_input


def ReGetContig_main(args):
    """
    class abc():
        pass

    args=abc()

    args.scaff_file = "/lustre/home/xuyuxing/Database/Cuscuta/Csp/Csp.fa"
    args.threads = 5
    args.output_gff3 = "/lustre/home/xuyuxing/Work/Csp/HiC/contig.gff3"
    args.output_fasta = "/lustre/home/xuyuxing/Work/Csp/HiC/contig.fasta"
    """

    scaff_dir = Fasta(args.scaff_file)

    SeqRecord_list = []

    for scaff in scaff_dir:
        scaff_seq = Seq(str(scaff))
        SeqRecord_tmp = SeqRecord(scaff_seq, scaff.name)
        SeqRecord_list.append((SeqRecord_tmp,))

    cmd_result = multiprocess_running(
        contig_range, SeqRecord_list, args.threads, silence=True)

    with open(args.output_gff3, 'w') as f:
        for i in cmd_result:
            GFF.write([cmd_result[i]['output']], f)

    with open(args.output_fasta, 'w') as f:
        for i in cmd_result:
            SeqRecord_tmp = cmd_result[i]['output']

            record_feature = SeqRecord_tmp.features
            for contig_feature in record_feature:
                start = int(contig_feature.location.start)
                end = int(contig_feature.location.end)
                ID = contig_feature.qualifiers['ID'][0]
                contig_seq = SeqRecord_tmp.seq[start:end]

                contig_record = SeqRecord(contig_seq, id=ID,
                                          description="%s:%d-%d" % (SeqRecord_tmp.id, start + 1, end))
                f.write(contig_record.format("fasta"))
                # SeqIO.write([contig_record], args.output_fasta, "fasta")

# GeneStructure


def GeneStructure_main(args):
    """
    class abc():
        pass

    args = abc()
    args.gff_file = "/lustre/home/xuyuxing/Database/Plant_genome/raw_data/Gastrodia_elata/GWHAAEX00000000.gff"
    args.genome_fasta_file = "/lustre/home/xuyuxing/Database/Plant_genome/raw_data/Gastrodia_elata/GWHAAEX00000000.genome.fasta"
    args.gene_id = "evm.TU.scaffold_0.420"
    args.output_file = None
    """

    if args.output_file is None:
        file = sys.stdout
    else:
        file = open(args.output_file, 'w')

    genome_dict = read_fasta_by_faidx(args.genome_fasta_file)

    in_file = args.gff_file

    with open(in_file, 'r') as in_handle:
        get_flag = False
        for rec in GFF.parse(in_handle):
            for gene in rec.features:
                if gene.id == args.gene_id:
                    get_flag = True
                    break
            if get_flag:
                break

    contig_record = genome_dict[rec.id]
    sub_seq_tmp = contig_record.sub(1, 100, '+', RNA=False)

    if get_flag:
        if gene.strand == 1:
            strand = "+"
        else:
            strand = "-"

        file.write("# Find gene: %s, on contig %s, strand %s\n" %
                   (gene.id, rec.id, strand))
        file.write("# There are %d mRNA in gene\n\n" % len(gene.sub_features))

        num = 1
        for mRNA in gene.sub_features:
            file.write("mRNA %d (ID: %s): \n\n" % (num, mRNA.id))
            file.write("total seq: \n")
            start = int(mRNA.location.start) + 1
            end = int(mRNA.location.end)
            seq_tmp = contig_record.sub(start, end, strand, RNA=False)
            file.write(">%s:%d-%d:%s\n" % (rec.id, start, end, strand))
            file.write("%s\n\n" % seq_tmp)

            # upstream_range 3k
            file.write("upstream range (3k): \n")
            if strand == "+":
                start = int(mRNA.location.start) + 1 - 3000
                end = int(mRNA.location.start) + 1 - 1
            else:
                end = int(mRNA.location.end) + 3000
                start = int(mRNA.location.end) + 1
            seq_tmp = contig_record.sub(start, end, strand, RNA=False)
            file.write(">%s:%d-%d:%s\n" % (rec.id, start, end, strand))
            file.write("%s\n\n" % seq_tmp)

            # mRNA sub feature
            if strand == "+":
                print_sub_feature = sorted([i for i in mRNA.sub_features if i.type != 'exon'],
                                           key=lambda x: int(x.location.start) + 1, reverse=False)
            else:
                print_sub_feature = sorted([i for i in mRNA.sub_features if i.type != 'exon'],
                                           key=lambda x: int(x.location.start) + 1, reverse=True)

            for index in range(0, len(print_sub_feature)):
                feature_tmp = print_sub_feature[index]

                file.write("%s\t(ID: %s): \n" %
                           (feature_tmp.type, feature_tmp.id))
                start = int(feature_tmp.location.start) + 1
                end = int(feature_tmp.location.end)
                seq_tmp = contig_record.sub(start, end, strand, RNA=False)
                file.write(">%s:%d-%d:%s\n" % (rec.id, start, end, strand))
                file.write("%s\n\n" % seq_tmp)

                if index != len(print_sub_feature) - 1:

                    if strand == "+":
                        start = int(feature_tmp.location.end) + 1
                        feature_tmp_next = print_sub_feature[index + 1]
                        end = int(feature_tmp_next.location.start) + 1 - 1
                    else:
                        # start = int(feature_tmp.location.start) + 1 - 1
                        # feature_tmp_next = print_sub_feature[index + 1]
                        # end = int(feature_tmp_next.location.end) + 1
                        end = int(feature_tmp.location.start) + 1 - 1
                        feature_tmp_next = print_sub_feature[index + 1]
                        start = int(feature_tmp_next.location.end) + 1

                    if abs(start - end) > 1:
                        # print(start, end)
                        seq_tmp = contig_record.sub(
                            start, end, strand, RNA=False)
                        file.write("intron: \n")
                        file.write(">%s:%d-%d:%s\n" %
                                   (rec.id, start, end, strand))
                        file.write("%s\n\n" % seq_tmp)

            # downstream_range 3k
            file.write("downstream range (3k): \n")
            if strand == "+":
                start = int(mRNA.location.end) + 1
                end = int(mRNA.location.end) + 3001
            else:
                start = int(mRNA.location.start) + 1 - 3000
                end = int(mRNA.location.start) + 1 - 1
            seq_tmp = contig_record.sub(start, end, strand, RNA=False)
            file.write(">%s:%d-%d:%s\n" % (rec.id, start, end, strand))
            file.write("%s\n\n" % seq_tmp)

            num += 1

    else:
        file.write("# Can not find gene: %s\n" % gene.id)

# SeqMerge


def SeqMerge_main(args):
    """
    class abc():
        pass

    args = abc()

    args.ID_merge_file = "/lustre/home/xuyuxing/Database/Gel/genome/assembly/hic/adjust/merge.txt"
    args.db_fasta_file = "/lustre/home/xuyuxing/Database/Gel/genome/assembly/hic/adjust/Gel.genome.v1.0.FINAL.fasta"
    args.output_file = "/lustre/home/xuyuxing/Database/Gel/genome/assembly/hic/adjust/merged.fasta"
    """

    file_info = tsv_file_dict_parse(
        args.ID_merge_file, fieldnames=['old_id', 'new_id'])
    ID_dict = {}
    for i in file_info:
        old_id, new_id = file_info[i]['old_id'], file_info[i]['new_id']
        if new_id not in ID_dict:
            ID_dict[new_id] = []
        ID_dict[new_id].append(old_id)

    fasta_dict = read_fasta_by_faidx(args.db_fasta_file)

    merge_used_id = []
    for new_id in ID_dict:
        merge_seq = ""
        for old_id in ID_dict[new_id]:
            merge_seq += fasta_dict[old_id].seq
            merge_used_id.append(old_id)
        new_contig_seq = BioSeq(seq=merge_seq, seqname=new_id)
        ID_dict[new_id] = new_contig_seq
        # new_contig_seq.write_to_file(args.masked_ref_genome)

    for id in fasta_dict:
        if id not in merge_used_id:
            fasta_record = fasta_dict[id]
            fasta_record.write_to_file(args.output_file)

    for id in ID_dict:
        fasta_record = ID_dict[id]
        fasta_record.write_to_file(args.output_file)

# GCproportion


def GCproportion_main(args):
    """
    class abc():
        pass

    args = abc()

    args.fasta_file = "/lustre/home/xuyuxing/Database/Gel/genome/assembly/contig_filter/Gel.genome.v1.0.final.rename.fasta"
    args.output_file = "/lustre/home/xuyuxing/Database/Gel/genome/assembly/contig_filter/Gel.genome.v1.0.final.rename.GC.txt"
    """
    seq_dict = read_fasta_by_faidx(args.fasta_file)

    with open(args.output_file, 'w') as f:
        for contig_id in seq_dict:
            tmp_seq = str(seq_dict[contig_id].seq)
            GC_num = tmp_seq.count('C') + tmp_seq.count('G') + \
                tmp_seq.count('c') + tmp_seq.count('g')
            N_num = tmp_seq.count('N') + tmp_seq.count('n')
            GC_pro = GC_num / (len(tmp_seq) - N_num)
            f.write("%s\t%.5f\n" % (contig_id, GC_pro))

# SeqEntropy


def SeqEntropy_main(args):
    """
    class abc():
        pass

    args = abc()

    args.input_fasta = "//lustre/home/xuyuxing/Database/Plant_genome/clean_data/Phelipanche_aegyptiaca/no_clean/T99112N0.gene_model.protein.fasta"
    args.output_file = "/lustre/home/xuyuxing/Database/Plant_genome/clean_data/Phelipanche_aegyptiaca/no_clean/T99112N0.gene_model.entropy.txt"
    args.word_size = 3
    """

    record_dict = read_fasta_by_faidx(args.input_fasta)

    with open(args.output_file, 'w') as f:
        for gene_id in record_dict:
            entropy = sb.sequence_entropy(
                record_dict[gene_id].seq, args.word_size)
            f.write("%s\t%f\n" % (gene_id, entropy))
