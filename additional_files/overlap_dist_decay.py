import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

pairs_file1 = sys.argv[1]
pairs_file2 = sys.argv[2]
reads_in_pairs_file1 = float(sys.argv[3])
reads_in_pairs_file2 = float(sys.argv[4])
pairs_file_name1 = pairs_file1[0:pairs_file1.find('.')]
pairs_file_name2 = pairs_file2[0:pairs_file2.find('.')]
contacts = np.zeros((1531935,), dtype=float)
distance = np.arange(0, 15319350, 10)
contact_probability = np.zeros((1531935,), dtype=float)
dist_orientation_decay = np.vstack((distance, contact_probability)).T
contacts2 = np.zeros((1531935,), dtype=float)
distance2 = np.arange(0, 15319350, 10)
contact_probability2 = np.zeros((1531935,), dtype=float)
dist_orientation_decay2 = np.vstack((distance, contact_probability)).T


def create_distance_decay_plot_list(file_to_read, num_of_reads, file_to_read2, num_of_reads2):
    """
    :param num_of_reads2:
    :param file_to_read2:
    :param num_of_reads:
    :param file_to_read:
    :return:
    this method will parse through a .pairs file and fill in arrays with counts for each distance into respective bins
    then normalize the counts across all counts
    """
    read_pairs_file = open(file_to_read, 'r')
    pairs_line = read_pairs_file.readline()
    while pairs_line[0] == '#':  # will skip over header so that it is not read
        pairs_line = read_pairs_file.readline()
    while pairs_line != '':
        pairs_arg_list = pairs_line.split(' ')
        pairs_arg_list = list(filter(None, pairs_arg_list))
        pairs_arg_list = pairs_arg_list[0].split('\t')
        temp = pairs_arg_list[len(pairs_arg_list) - 1].split('\n')
        pairs_arg_list.pop()
        pairs_arg_list.append(temp[0])
        pair_obj = Pair(pairs_arg_list[0], pairs_arg_list[1], pairs_arg_list[2], pairs_arg_list[3], pairs_arg_list[4],
                        pairs_arg_list[5], pairs_arg_list[6], pairs_arg_list[7])
        if pairs_arg_list[1] == pairs_arg_list[3]:
            contacts[int(pair_obj.calculate_distance() / 10)] = contacts[int(pair_obj.calculate_distance() / 10)] + 1
        pairs_line = read_pairs_file.readline()
    dist_orientation_decay[:, 1] += (contacts / num_of_reads)

    read_pairs_file2 = open(file_to_read2, 'r')
    pairs_line2 = read_pairs_file2.readline()
    while pairs_line2[0] == '#':  # will skip over header so that it is not read
        pairs_line2 = read_pairs_file2.readline()
    while pairs_line2 != '':
        pairs_arg_list2 = pairs_line2.split(' ')
        pairs_arg_list2 = list(filter(None, pairs_arg_list2))
        pairs_arg_list2 = pairs_arg_list2[0].split('\t')
        temp2 = pairs_arg_list2[len(pairs_arg_list2) - 1].split('\n')
        pairs_arg_list2.pop()
        pairs_arg_list2.append(temp2[0])
        pair_obj2 = Pair(pairs_arg_list2[0], pairs_arg_list2[1], pairs_arg_list2[2], pairs_arg_list2[3],
                         pairs_arg_list2[4],
                         pairs_arg_list2[5], pairs_arg_list2[6], pairs_arg_list2[7])
        if pairs_arg_list2[1] == pairs_arg_list2[3]:
            contacts2[int(pair_obj2.calculate_distance() / 10)] = contacts2[int(pair_obj2.calculate_distance() / 10)] + 1
        pairs_line2 = read_pairs_file2.readline()
    dist_orientation_decay2[:, 1] += (contacts2 / num_of_reads2)

    contacts.tofile(pairs_file_name1 + '_contacts_in_bin.txt', sep=' ', format='%s')
    contacts2.tofile(pairs_file_name2 + '_contacts_in_bin.txt', sep=' ', format='%s')
    create_distance_decay_plot(dist_orientation_decay, dist_orientation_decay2)


def create_distance_decay_plot(lis, lis2):
    """
    :param lis2:
    :param lis:
    :return:
    this method creates matplotlib graphs for passed in lists
    """
    plot_name = pairs_file_name1 + pairs_file_name2 + 'ori_decay_.png'
    df = pd.DataFrame({
        'x_axis': [lis[i][0] for i in range(0, 201, 2)],
        'y_axis': [lis[i][1] for i in range(1, 202, 2)]
    })

    df2 = pd.DataFrame({
        'x_axis': [lis2[i][0] for i in range(0, 201, 2)],
        'y_axis': [lis2[i][1] for i in range(1, 202, 2)]
    })

    plt.plot('x_axis', 'y_axis', data=df, linestyle='-', marker='o', label='Quiescent_data')
    plt.plot('x_axis', 'y_axis', data=df2, linestyle='-', marker='o', label='Log_data')
    plt.xlim(0, 1000)
    plt.ylim(0.0008, 0.0122)
    plt.xlabel('bins')
    plt.ylabel('relative counts')
    plt.title('SAME reads Q + Log')
    plt.legend(fontsize='10')
    plt.savefig(plot_name)


class Pair:
    def __init__(self, read_id, chrom1, pos1, chrom_2, pos_2, strand1, strand2, pair_type):
        self.read_id = read_id
        self.chrom1 = chrom1
        self.pos1 = pos1
        self.chrom2 = chrom_2
        self.pos2 = pos_2
        self.strand1 = strand1
        self.strand2 = strand2
        self.pair_type = pair_type

    def calculate_distance(self):
        return abs(int(self.pos1) - int(self.pos2))


create_distance_decay_plot_list(pairs_file1, reads_in_pairs_file1, pairs_file2, reads_in_pairs_file2)
