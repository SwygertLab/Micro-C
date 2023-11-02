import sys

pairs_file = sys.argv[1]
pairs_file_name = pairs_file[0:pairs_file.find(".")]
noIN_pairs_file_name = pairs_file_name + "_noIN.pairs"
noIN_pairs_file = open(noIN_pairs_file_name, 'x')  # file opening for noIN
IN_read_file_name = pairs_file_name + '_IN_reads.pairs'
IN_reads_file = open(IN_read_file_name, 'x')  # file opening for only IN reads
SAME_reads_file_name = pairs_file_name + '_SAME_reads.pairs'
SAME_reads_file = open(SAME_reads_file_name, 'x')  # file opening for only SAME reads
OUT_reads_file_name = pairs_file_name + '_OUT_reads.pairs'
OUT_reads_file = open(OUT_reads_file_name, 'x')  # file opening for only OUT reads
log_file_name = pairs_file_name + '_filter.log'
log_file = open(log_file_name, 'x')
# python filter_orientations.py all_r1_mergedall_r2_merged_output.pairs


def read_file(pairs_file2):
    pair_obj_ls = []
    dup = 0
    unk_chr = 0
    IN_reads_count = 0
    OUT_reads_count = 0
    SAME_reads_count = 0
    noIN_reads_count = 0
    included_reads = 0
    filtered_reads = 0
    read_pairs_file = open(pairs_file2, 'r')
    pairs_line = read_pairs_file.readline()

    while pairs_line[0] == '#':  # includes the header info in each .pairs file orientation
        noIN_pairs_file.write(pairs_line)
        IN_reads_file.write(pairs_line)
        OUT_reads_file.write(pairs_line)
        SAME_reads_file.write(pairs_line)
        pairs_line = read_pairs_file.readline()
    while pairs_line != '':  # splits the information in each line to more accurately compare then determines orientation then writes to appropriate file location
        pairs_arg_list = pairs_line.split(' ')
        pairs_arg_list = list(filter(None, pairs_arg_list))
        pairs_arg_list = pairs_arg_list[0].split('\t')
        temp = pairs_arg_list[len(pairs_arg_list) - 1].split('\n')
        pairs_arg_list.pop()
        pairs_arg_list.append(temp[0])
        pair_obj = Pair(pairs_arg_list[0], pairs_arg_list[1], pairs_arg_list[2], pairs_arg_list[3], pairs_arg_list[4],
                        pairs_arg_list[5], pairs_arg_list[6], pairs_arg_list[7])
        if not pair_obj_ls.__contains__(pair_obj):
            pair_obj_ls.append(pair_obj)
            if pair_obj.if_unk_chr():
                filtered_reads = filtered_reads + 1
                unk_chr = unk_chr + 1
            elif pair_obj.if_rDNA():
                filtered_reads = filtered_reads + 1
            elif pair_obj.determine_orientation() == 'SAME':
                noIN_pairs_file.write(pairs_arg_list[0] + '\t' + pairs_arg_list[1] + '\t' + pairs_arg_list[2] + '\t' +
                                      pairs_arg_list[3] + '\t' + pairs_arg_list[4] + '\t' + pairs_arg_list[5] + '\t' +
                                      pairs_arg_list[6] + '\t' + pairs_arg_list[7] + '\n')
                noIN_reads_count = noIN_reads_count + 1
                SAME_reads_file.write(pairs_arg_list[0] + '\t' + pairs_arg_list[1] + '\t' + pairs_arg_list[2] + '\t' +
                                      pairs_arg_list[3] + '\t' + pairs_arg_list[4] + '\t' + pairs_arg_list[5] + '\t' +
                                      pairs_arg_list[6] + '\t' + pairs_arg_list[7] + '\n')
                SAME_reads_count = SAME_reads_count + 1
            elif pair_obj.determine_orientation() == 'OUT':
                noIN_pairs_file.write(pairs_arg_list[0] + '\t' + pairs_arg_list[1] + '\t' + pairs_arg_list[2] + '\t' +
                                      pairs_arg_list[3] + '\t' + pairs_arg_list[4] + '\t' + pairs_arg_list[5] + '\t' +
                                      pairs_arg_list[6] + '\t' + pairs_arg_list[7] + '\n')
                noIN_reads_count = noIN_reads_count + 1
                OUT_reads_file.write(pairs_arg_list[0] + '\t' + pairs_arg_list[1] + '\t' + pairs_arg_list[2] + '\t' +
                                     pairs_arg_list[3] + '\t' + pairs_arg_list[4] + '\t' + pairs_arg_list[5] + '\t' +
                                     pairs_arg_list[6] + '\t' + pairs_arg_list[7] + '\n')
                OUT_reads_count = OUT_reads_count + 1
            elif pair_obj.determine_orientation() == 'IN':
                IN_reads_file.write(pairs_arg_list[0] + '\t' + pairs_arg_list[1] + '\t' + pairs_arg_list[2] + '\t' +
                                    pairs_arg_list[3] + '\t' + pairs_arg_list[4] + '\t' + pairs_arg_list[5] + '\t' +
                                    pairs_arg_list[6] + '\t' + pairs_arg_list[7] + '\n')
                IN_reads_count = IN_reads_count + 1
            included_reads = included_reads + 1
        else:
            dup = dup + 1
            filtered_reads = filtered_reads + 1
        # Writing out to log file filtering info.
        log_file.write("The number of duplicates filtered out is " + dup + ".")
        log_file.write("The number of reads filtered out to unknown chromosomes is  " + unk_chr + ".")
        log_file.write("The total number of IN reads is " + IN_reads_count + ".")
        log_file.write("The total number of OUT reads is " + OUT_reads_count + ".")
        log_file.write("The total number of SAME reads is " + SAME_reads_count + ".")
        log_file.write("The total number of noIN reads is " + noIN_reads_count + ".")
        log_file.write("The total number of filtered reads is" + filtered_reads + ".")
        log_file.write("The total number of included reads is " + included_reads + ".")
        total_reads = filtered_reads + included_reads
        log_file.write("The total number of filtered and included reads is " + total_reads + ".")
        pairs_line = read_pairs_file.readline()


class Pair:
    def __init__(self, read_id, chrom1, pos1, chrom_mate, pos_mate, strand1, strand2, pair_type):
        self.read_id = read_id
        self.chrom1 = chrom1
        self.pos1 = pos1
        self.chrom_mate = chrom_mate
        self.pos_mate = pos_mate
        self.strand1 = strand1
        self.strand2 = strand2
        self.pair_type = pair_type

    def if_unk_chr(self):
        if self.chrom1 == 'chrM':
            return True

    def if_rDNA(self):
        if self.chrom1 == 'chrXII' or self.chrom_mate == 'chrXII':
            if 451573 < self.pos1 < 470000 or 451573 < self.pos_mate < 470000:
                return True

    def determine_orientation(self):
        pos1 = int(self.pos1)
        pos2 = int(self.pos_mate)

        if self.chrom1 != self.chrom_mate:
            return 'OUT'
        # note that I have intentionally left out the instance where one of the reads is in this rDNA area but the mate is not on the same chromosome, so you can still potentially get reads from this region
        if self.strand1 == self.strand2:
            return 'SAME'
        if pos1 < pos2:
            if self.strand1 == '+':
                return 'IN'
            else:
                return 'OUT'
        else:
            if self.strand1 == '+':
                return 'OUT'
            else:
                return 'IN'


read_file(pairs_file)
