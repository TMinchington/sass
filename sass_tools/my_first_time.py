def get_first_times(data_path):
    """
    Cycles through files and assigns first 
    time points
    """
    file_ls = os.listdir(data_path)

    bd = data_path + '/backed_down.txt'
    bd_hc = data_path + '/backed_down-hc.txt'


    files_ls = [bd, bd_hc]

    bd_first_dic = {}

    times_dic = {}

    for file in files_ls:

        with open(file) as o_file:

            first_line = True
            second_line = True

            for line in o_file:

                if first_line:

                    first_line = False
                    bd_first_dic[file] = {}

                    continue

                (time, nuc, num_spots, nucx, nucy, nucz,
                 spot, spotx, spoty, spotz, variable, channel, value,
                 background, corrected_val, dis_mid, above_line) = line.strip().split('\t')

                if second_line:

                    second_line = False

                    # frame_interval = float(time)
                
                times_dic[time] = 0

                # if corrected value is less than 0 then the nucleus will be lost

                if nuc not in bd_first_dic[file] and float(corrected_val) > 0:

                    bd_first_dic[file][nuc] = float(time)

                elif nuc in bd_first_dic[file] and float(corrected_val) > 0:

                    if float(time) < bd_first_dic[file][nuc]:

                        bd_first_dic[file][nuc] = float(time)

                else:

                    continue

    times_ls = sorted([float(t) for  t in list(times_dic)])
    frame_interval = times_ls[1] - times_ls[0]    

    return bd_first_dic, frame_interval


def load_corrected(data_path):

    import os

    time_path = data_path+'/key_time.txt'

    if not os.path.isfile(time_path):



        check_time = input('No key time file detected. Would you like to stop and create this file now? (y/n)')

        if check_time == 'y':

            time_file = open(time_path, 'w')
            time_file.close()
            exit(f'key_time file created at {time_path}')

        else:

            return []

    correct_ls = []

    with open(time_path) as time_file:

        for line in time_file:

            split_line = line.strip().split('\t')
            if len(split_line) < 3:

                continue
            split_line = [float(x) for x in split_line]
            if len(split_line) != 3:

                continue

            correct_ls.append(tuple(split_line))

    return correct_ls


def correct_times_and_mark_first(first_dic, data_path, frame_interval):

    """
    Uses the dictionary of first timepoints (first_dic) from get_first_times to annotate the first time points in the bd
    and bd_hc.

    flipped allows the code to determine if the flipper script has been run on the data so it can include the
    files or not.

    data path is the path to the data so that it can be accessed


    :param first_dic:
    :param data_path:
    :param data_path:
    :return:
    """

    bd = data_path + '/backed_down.txt'
    bd_hc = data_path + '/backed_down-hc.txt'

    files_ls = [bd, bd_hc]

    temp_files = []
    correct_ls = load_corrected(data_path)
    # print(correct_ls)
    # exit()
    for file in files_ls:

        with open(file) as o_file:
            first_line = True
            temp_out = file.replace('.', '___temp___.')
            temp_files.append(temp_out)
            o_temp_out = open(temp_out, 'w')

            for line in o_file:

                if first_line:
                    first_line = False
                    o_temp_out.write(line.replace('\n', '\ttime_shift\tfirst_time\n'))
                    continue

                (time, nuc, num_spots, nucx, nucy, nucz,
                 spot, spotx, spoty, spotz, variable, channel, value,
                 background, corrected_val, dis_mid, above_line) = line.strip().split('\t')

                time = float(time)
                nucx = float(nucx)
                
                try:

                    nuc_first = first_dic[file][nuc]

                except KeyError:

                    continue

                first_nuc = 0

                if nuc_first == time:

                    first_nuc = 1

                corrected_time = time

                for x in correct_ls:
                    x1, x2, shift = x

                    if x1 <= nucx < x2:
                        corrected_time += shift*frame_interval
                        # if corrected_time <0:
                        #     print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
                        #     print(corrected_time)
                        #     print(line)
                        #     print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
                o_temp_out.write('\t'.join([str(i) for i in [time, nuc, num_spots, nucx, nucy, nucz,
                 spot, spotx, spoty, spotz, variable, channel, value,
                 background, corrected_val, dis_mid, above_line, corrected_time, first_nuc]])+'\n')

def tidy_up():

    bd = data_path + '/backed_down.txt'
    bd_hc = data_path + '/backed_down-hc.txt'

    file_ls = [bd, bd_hc]

    for file in file_ls:

        temp_file = file.replace('.', '___temp___.')
        os.remove(file)
        os.rename(temp_file, file)

if __name__ == "__main__":

    """
    
    Looks through backed down and other files to determine the first timepoint for each nucleus where expression is determined.
    
    """

    import argparse
    import os
    import shutil

    parser = argparse.ArgumentParser()

    parser.add_argument("experiment_dir")

    args = parser.parse_args()

    data_path = args.experiment_dir + '/time_data'

    if os.path.isfile(args.experiment_dir+'/key_time.txt'):

        shutil.move(args.experiment_dir+'/key_time.txt', data_path + '/key_time.txt')

    first_dic, frame_interval = get_first_times(data_path)

    correct_times_and_mark_first(first_dic, data_path, frame_interval)

    tidy_up()