import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

pairs_file = sys.argv[1]
reads_in_pairs_file = float(sys.argv[2])  # the script needs to know the total amount of reads in the respective orientation
distance_graphed = sys.argv[3]  # the distance you want to get graphed
# note that python stops 1 before the number so if you want to stop at 200 the input needs to be 201

#in order to set graph limits change False to True and then add the 4 limit numbers
set_limit = sys.argv[4]
if set_limit is "True"
    x_limit = [sys.argv[5], sys.argv[6]]  # example values [0, 1000]
    y_limit = [sys.argv[7], sys.argv[8]]  # example values [0, 0.025]


pairs_file_name = pairs_file[0:pairs_file.find('.')]
contacts = np.zeros((1531935,),
                    dtype=float)  # the number of bins is supposed to be the summation of all the chromosomes base
# pair count to the next highest 10 and then take out the 0
distance = np.arange(0, 15319350, 10)
contact_probability = np.zeros((1531935,), dtype=float)
dist_orientation_decay = np.vstack((distance, contact_probability)).T


def create_distance_decay_plot_list(file_to_read):
    """
    :param file_to_read:
    :return:
    this method will parse through a .pairs file and a graph for that .pairs file
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
    dist_orientation_decay[:, 1] += (contacts / reads_in_pairs_file)
    contacts.tofile(pairs_file_name + '_contacts_in_bin.txt', sep=' ', format='%s')
    create_distance_decay_plot(dist_orientation_decay, 'ori_decay_', set_limit)


def create_distance_decay_plot(lis, type_of_lis, set_limit="False"):
    plot_name = pairs_file_name + type_of_lis + '.png'
    plot_name2 = pairs_file_name + type_of_lis + '.svg'
    df = pd.DataFrame({
        'x_axis': [lis[i][0] for i in range(0, int(distance_graphed), 2)],
        'y_axis': [lis[i][1] for i in range(1, int(distance_graphed)+1, 2)]
    })
    plt.plot('x_axis', 'y_axis', data=df, linestyle='-', marker='o', label=type_of_lis)
    if set_limit is "True"
        plt.xlim(x_limit[0], x_limit[1])
        plt.ylim(y_limit[0], y_limit[1])
    plt.xlabel('base pair')
    plt.ylabel('relative counts')
    plt.title(pairs_file_name + 'distance decay plot')
    plt.legend(fontsize='10')
    plt.savefig(plot_name)
    plt.plot_name(plot_name2)
    plt.close()


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


create_distance_decay_plot_list(pairs_file)
