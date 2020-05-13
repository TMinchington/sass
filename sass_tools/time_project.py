def load_backed_data(data_path, outfile):
    import numpy as np
    data_dic = {}
    pos_dic = {}
    outfile = open(outfile, 'w')
    with open(data_path) as o_data:

        first_line = True

        for line in o_data:
            # print(line)
            split_line = line.strip().split('\t')

            if first_line:

                tdex = split_line.index('time')
                ndex = split_line.index('nuc')
                ndexX = split_line.index('nucx')
                ndexY = split_line.index('nucy')
                chandex = split_line.index('channel')
                intent_type = split_line.index('variable')
                valuedex = split_line.index('value')

                first_line = False

                continue


            nuc = split_line[ndex]
            variable = split_line[intent_type]
            channel = split_line[chandex]
            # print(split_line[valuedex])
            value = float(split_line[valuedex])

            try:

                data_dic[nuc][variable][channel].append(value)

            except KeyError:

                try:

                    data_dic[nuc][variable][channel] = [value]

                except KeyError:

                    try:

                        data_dic[nuc][variable] = {channel: [value]}

                    except KeyError:

                        data_dic[nuc] = {variable: {channel: [value]}}

            try:

                pos_dic[nuc]['x'].append(float(split_line[ndexX]))
                pos_dic[nuc]['y'].append(float(split_line[ndexY]))

            except KeyError:

                pos_dic[nuc] = {'x': [float(split_line[ndexX])], 'y': [float(split_line[ndexY])]}

    outfile.write('nuc\tmean_x\tmean_y\tvariable\tchannel\tt_mean\tt_median\tt_max\tt_sum\n')
    for nuc in data_dic:

        for variable in data_dic[nuc]:

            for chan in data_dic[nuc][variable]:

                out_ls = [nuc, np.mean(pos_dic[nuc]['x']), np.mean(pos_dic[nuc]['y']), variable,
                          chan, np.mean(data_dic[nuc][variable][chan]), np.median(data_dic[nuc][variable][chan]),
                          np.max(data_dic[nuc][variable][chan]), np.sum(data_dic[nuc][variable][chan])]

                out_str = '\t'.join([str(x) for x in out_ls])+'\n'

                outfile.write(out_str)


if __name__ == "__main__":

    import argparse
    import os

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser()

    parser.add_argument("exp_dir")

    args = parser.parse_args()

    data_path = args.exp_dir + '/time_data/backed_down.txt'
    data_path_out = args.exp_dir + '/time_data/time_projected.txt'

    load_backed_data(data_path, data_path_out)