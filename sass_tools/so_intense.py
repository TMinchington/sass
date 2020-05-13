"""

Python script takes the output file from elipse and then adds in the intensities

"""

def assign_intensities(exp_dir, spot_dic, pos_file):

    import os

    out_path = exp_dir + '/time_data/position_data-intense.txt'

    outfile = open(out_path, 'w')

    with open(pos_file) as o_pos:

        first_line = True

        for line in o_pos:

            split_line = line.strip().split('\t')

            if first_line:

                first_line = False
                outfile.write(line.replace('\n', '\tvariable\tchannel\tvalue\n'))
                tdex = split_line.index('time')
                sdex = split_line.index('spot')

                continue

            t = float(split_line[tdex])
            spot = split_line[sdex]
            # print(set(spot_dic))
            for variable in spot_dic[t][spot]:

                for chan in spot_dic[t][spot][variable]:

                    value = spot_dic[t][spot][variable][chan]
                    # print(t, spot, channel, variable, value)
                    outfile.write(line.replace('\n', f'\t{variable}\t{chan}\t{value}\n'))


def get_spot_info(exp_dir, spot_files):

    spot_dic = {}

    for file in spot_files['_Intensity_']:
        print(file)
        with open(file) as o_file:

            data_found = False

            for line in o_file:
                # print(line)
                split_line = line.strip().split(',')
                # print(end='\r')

                if 'ID' in line and not data_found:

                    data_found = True

                    tdex = split_line.index('Time')
                    sdex = split_line.index('ID')
                    chandex = split_line.index('Channel')
                    intent_type = split_line[0]

                    continue

                elif not data_found:

                    continue

                t = float(split_line[tdex])
                spot = split_line[sdex]
                channel = split_line[chandex]
                value = split_line[0]

                try:

                    spot_dic[t][spot][intent_type][channel] = value

                except KeyError:

                    try:

                        spot_dic[t][spot][intent_type] = {channel: value}

                    except KeyError:

                        try:

                            spot_dic[t][spot] = {intent_type: {channel: value}}

                        except KeyError:

                            spot_dic[t] = {spot: {intent_type: {channel: value}}}

    return spot_dic



if __name__ == "__main__":


    import argparse
    import os
    from spot_the_dffierence import get_files, get_experiments_folder

    parser = argparse.ArgumentParser()

    parser.add_argument("exp_dir")

    args = parser.parse_args()

    folder_dic = get_experiments_folder(args.exp_dir)
    spot_files = get_files(folder_dic['spot'])

    pos_file = args.exp_dir+'/time_data/position_data.txt'

    if not os.path.isfile(pos_file):

        exit('Position file not found')

    spotted_dic = get_spot_info(args.exp_dir, spot_files)

    assign_intensities(args.exp_dir, spotted_dic, pos_file)