"""
Import and calculate background

"""

def plot_spots(back_dic):

    import matplotlib.pyplot as plt

    plot_dic = {}

    for t in back_dic:

        for spot in back_dic[t]:

            for variable in back_dic[t][spot]:

                for chan in back_dic[t][spot][variable]:

                    value = back_dic[t][spot][variable][chan]

                    try:

                        plot_dic[variable][chan][t].append(value)

                    except KeyError:

                        try:

                            plot_dic[variable][chan][t] = [value]

                        except KeyError:

                            try:

                                plot_dic[variable][chan] = {t : [value]}

                            except KeyError:

                                    plot_dic[variable] = {chan: {t: [value]}}

    background_eq_dic = {}

    for variable in plot_dic:
        background_eq_dic[variable] = {}
        for chan in plot_dic[variable]:

            t_ls = []
            v_ls = []

            for t in plot_dic[variable][chan]:

                for val in plot_dic[variable][chan][t]:

                    t_ls.append(float(t))
                    v_ls.append(float(val))

            background_eq_dic[variable][chan] = get_best_fit(t_ls, v_ls, 1)
            # plt.scatter(t_ls, v_ls)
            # plt.title(f'{variable} {chan}')
            # plt.show()
            # plt.close()

    return background_eq_dic

def get_best_fit(x_ls, y_ls, curvey):

    import numpy as np

    f = np.poly1d(np.polyfit(x_ls, y_ls, curvey))
    # x = range(0, int(max(x_ls)))
    x_new = np.linspace(min(x_ls), max(x_ls), 100)

    y_new = f(x_new)

    for i in x_new:

        print(f(i))

    # import matplotlib.pyplot as plt
    #
    # plt.plot(x_new, y_new)
    # plt.scatter(x_ls, y_ls,  s=50, marker='x')
    # plt.show()
    # plt.close()

    return f

def remove_background(intense_path, back_eq_dic, t_interval, nuc_dic, mid_line):
    import os
    split_dic = {}
    head_line = ''
    time_ls = []
    chan_ls = []
    intent_type_ls = []

    outfile = open(os.path.split(intense_path)[0]+'/backed_down.txt', 'w')
    points_file = open(os.path.split(intense_path)[0] + '/mid_points2.txt', 'w')
    with open(intense_path) as o_intense:

        first_line = True

        for line in o_intense:

            split_line = line.strip().split('\t')

            if first_line:

                outfile.write(line.replace('\n', '\tbackground\tcorrected_val\tdis_mid\tabove_line\n'))
                # print(line)
                tdex = split_line.index('time')
                ndex = split_line.index('nuc')
                ndexX = split_line.index('nucx')
                ndexY = split_line.index('nucy')
                chandex = split_line.index('channel')
                intent_type = split_line.index('variable')
                value = split_line.index('value')
                first_line = False
                continue


            nuc = split_line[ndex]
            t = int(float(split_line[tdex]))
            # print(time_ls)
            time_ls.append(t)
            intent_type_ls.append(split_line[intent_type])
            chan_ls.append(split_line[chandex])
            try:

                split_dic[nuc][t].append(split_line)

            except KeyError:

                try:

                    split_dic[nuc][t] = [split_line]

                except KeyError:

                    split_dic[nuc] = {t: [split_line]}


    min_t = min(time_ls)
    max_t = max(time_ls)
    print(min_t, max_t)
    chan_ls = set(chan_ls)
    # time_ls = set(time_ls)
    intent_type_ls = set(intent_type_ls)

    for nuc in split_dic:

        for t in range(1, max_t+1):
            print(t, end='\r')
            t_out = t_interval*float(t)/60

            if t in split_dic[nuc]:

                for line in split_dic[nuc][t]:
                    line[0] = t_out
                    val = line[value]
                    chan = line[chandex]
                    intent = line[intent_type]
                    backG = back_eq_dic[intent][chan](t)
                    back_correct = float(val) - backG

                    if back_correct < 0:

                        back_correct = 0
                    nucx, nucy, nucz = nuc_dic[t][nuc]
                    dis_mid, above_line = get_distance_mid(mid_line, float(nucx), float(nucy), points_file)

                    out_ls = line + [backG, back_correct, dis_mid, above_line]
                    out_str = '\t'.join([str(x) for x in out_ls])+'\n'
                    outfile.write(out_str)
            else:

                for intent in intent_type_ls:

                    for chan in chan_ls:

                        backG = back_eq_dic[intent][chan](t)

                        try:

                            nucx, nucy, nucz = nuc_dic[t][nuc]

                        except KeyError:

                            continue

                        dis_mid, above_line = get_distance_mid(mid_line, float(nucx), float(nucy), points_file)

                        out_ls = [t_out, nuc, 0, nucx, nucy, nucz, 'NA', 'NA', 'NA', 'NA', intent, chan, 0, backG, 0, dis_mid, above_line]
                        out_str = '\t'.join([str(x) for x in out_ls]) + '\n'
                        outfile.write(out_str)

                continue


def get_mid(intense_path):

    import numpy as np
    import matplotlib.pyplot as plt
    import os

    curvey = 2
    x_ls = []
    y_ls = []

    image_file = os.path.split(intense_path)[0]+'mid_line_image.png'

    with open(intense_path) as o_intense:

        first_line = True

        for line in o_intense:

            split_line = line.strip().split('\t')

            if first_line:

                first_line = False

                x_dex = split_line.index('spotx')
                y_dex = split_line.index('spoty')

                continue

            x_ls.append(float(split_line[x_dex]))
            y_ls.append(float(split_line[y_dex]))

    f = np.poly1d(np.polyfit(x_ls, y_ls, curvey))
    # x = range(0, int(max(x_ls)))
    x_new = np.linspace(min(x_ls)-10, max(x_ls)+10, 300)
    y_new = f(x_new)

    plt.scatter(x_ls, y_ls, c='black')
    plt.plot(x_new, y_new, c='red')

    plt.savefig(image_file)
    plt.close()

    return (x_new, y_new, f)


def get_distance_mid(mid_line, x, y, points_file):

    x_mid, y_mid, f = mid_line

    dis_ls = []

    for i in range(0, len(x_mid)):

        mid_dis = ((x_mid[i] - x)**2) + ((y_mid[i] - y)**2)

        dis_ls.append(mid_dis)

    if y < f(x):

        trans = -1

    elif y > f(x):

        trans = 1

    else:

        trans = 0

    points_file.write(f'{x}\t{y}\t{f(x)}\t{min(dis_ls)*trans}\t{trans}\n')

    return min(dis_ls) * trans, trans


def OFFget_distance_mid(mid_line, nucx, nucy, points_file):

    x_ls, y_ls, f = mid_line
    dis_ls = []
    for i in range(0, len(x_ls)):

        x = x_ls[i]
        y = y_ls[i]

        temp_dis = (nucx-x)**2 + (nucy-y)**2

        dis_ls.append(temp_dis)

    mindex = dis_ls.index(min(dis_ls))

    x = x_ls[mindex]
    y = y_ls[mindex]
    mindis = min(dis_ls)

    if nucx < y:

        trans = -1

    else:

        trans = 1

    del dis_ls

    mindis = mindis**.5
    mindis = mindis*trans

    return mindis, trans


def return_mid_line(mid_line, path):

    with open(path+'/mid_line.txt', 'w') as o_mid:

        o_mid.write(f'x\ty\n')

        for x in range(0, len(mid_line[0])):

            o_mid.write(f'{mid_line[0][x]}\t{mid_line[1][x]}\n')



def plot_data(csv, bins_ls, bin_numbers):

    import matplotlib.pyplot as plt


    # print(len(bins_ls))

    counts = []

    for x in range(0, len(bins_ls) -1):

        counts.append(len(csv[(csv['min_dis'] > bins_ls[x]) & (csv['min_dis'] < bins_ls[x+1]) & (csv['spots'] == 1)]))

    counts2 = []
    plt.bar(x=bin_numbers, height=counts, color='orange')

    for x in range(0, len(bins_ls) - 1):

        counts2.append(len(csv[(csv['min_dis'] > bins_ls[x]) & (csv['min_dis'] < bins_ls[x + 1]) & (csv['spots'] >= 2)]))

    plt.bar(x=bin_numbers, height=counts, color='red')
    plt.show()
    plt.close()


if __name__ == "__main__":

    import argparse
    import os
    from spot_the_dffierence import get_experiments_folder, get_files, get_spot_positions, \
         get_nuc_positions, get_distance_filter, check_nuc_spots
    from so_intense import get_spot_info
    from pprint import pprint as pp
    import multiprocessing as mp

    parser = argparse.ArgumentParser()

    parser.add_argument("exp_dir")

    parser.add_argument("--i", default=20, help="time interval in seconds between points. Default: 20 seconds", type=float)

    args = parser.parse_args()

    intense_path = args.exp_dir+'/time_data/position_data-intense.txt'

    if not os.path.isfile(intense_path):

        exit('intensity file is missing please run previous scripts before running take it back again.')


    folders_dic = get_experiments_folder(args.exp_dir)

    pp(folders_dic)

    background_files = get_files(folders_dic['backG'])
    nuc_files = get_files(folders_dic['nuc'])

    nuc_dic = get_nuc_positions(nuc_files)

    back_dic = get_spot_info(args.exp_dir, background_files)

    back_eq_lines = plot_spots(back_dic)
    pp(nuc_dic)
    mid_line = get_mid(intense_path)
    remove_background(intense_path, back_eq_lines, args.i, nuc_dic, mid_line)