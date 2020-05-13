def import_axis(exp_dir, spot_files):

    from spot_the_dffierence import build_index_dic

    files_to_get = ['A', 'B', 'C']

    axis_dic = {}

    for x in spot_files['Ellipsoid_Axis']:

        for filex in files_to_get:

            # axis_dic[filex] = {}

            if 'Ellipsoid_Axis_'+filex+'.' in x:
                print(x)
                with open(x) as openX:

                    data_found = False

                    for line in openX:

                        split_line = line.strip().split(',')

                        if split_line[0].startswith('Ellipsoid Axis '+filex+' X') and not data_found:

                            data_found = True

                            column_index = build_index_dic(split_line)
                            print('coloumns indexed')
                            continue

                        elif not data_found:

                            continue

                        nuc_id = split_line[column_index['TrackID']]

                        nuc_time = float(split_line[column_index['Time']])

                        pos_x = float(split_line[column_index['Ellipsoid Axis '+filex+' X']])
                        pos_y = float(split_line[column_index['Ellipsoid Axis '+filex+' Y']])
                        pos_z = float(split_line[column_index['Ellipsoid Axis '+filex+' Z']])

                        try:

                            axis_dic[filex][nuc_time][nuc_id] = (pos_x, pos_y, pos_z)

                        except KeyError:

                            try:

                                axis_dic[filex][nuc_time] = {nuc_id: (pos_x, pos_y, pos_z)}

                            except KeyError:

                                axis_dic[filex] = {nuc_time: {nuc_id: (pos_x, pos_y, pos_z)}}

    return axis_dic


def plot_axis_at_time(axis_dic, t, nuc_dic):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    print(list(axis_dic))

    for axis in axis_dic:

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for nuc in axis_dic[axis][t]:

            nuc_pos = nuc_dic[t][nuc]

            axis_info = axis_dic[axis][t][nuc]

            fx_ls, fy_ls, fz_ls = make_further_ls(-2, 12, 200, nuc_pos[0], nuc_pos[1], nuc_pos[2], axis_info)

            xls = [nuc_pos[0]]+fx_ls
            yls = [nuc_pos[1]]+fy_ls
            zls = [nuc_pos[2]]+fz_ls

            ax.scatter(xls, yls, zls)

        plt.show()
        plt.close()


def make_further_ls(start, stop, num_pos, x, y, z, axis_info):

    import numpy as np

    fx_ls = []
    fy_ls = []
    fz_ls = []

    pos_to_plot = np.linspace(start, stop, num_pos)

    for ni in pos_to_plot:

        fx_ls.append(x+axis_info[0]*ni)
        fy_ls.append(y+axis_info[1]*ni)
        fz_ls.append(z+axis_info[2]*ni)

    return fx_ls, fy_ls, fz_ls


def make_positions_dic(axis_dic, nuc_dic):

    print("Generating all possible positions, this may take a short while depending on your computer...")

    from pprint import pprint as pp

    pos_dic = {}
    # tup_ls = []

    for t in axis_dic['C']:

        for nuc in axis_dic['C'][t]:

            nuc_pos = nuc_dic[t][nuc]

            axis_info = axis_dic['C'][t][nuc]

            fx_ls, fy_ls, fz_ls = make_further_ls(0, 25, 40, nuc_pos[0], nuc_pos[1], nuc_pos[2], axis_info)

            for i in range(0, len(fx_ls)):

                tup_key = (fx_ls[i], fy_ls[i], fz_ls[i])
                # tup_ls.append((fx_ls[i], fy_ls[i], fz_ls[i]))
                # print(type(tup_key))
                try:

                    pos_dic[t][tup_key].append(nuc)

                except KeyError:

                    try:

                        pos_dic[t][tup_key] = [nuc]

                    except KeyError:

                        pos_dic[t] = {tup_key: [nuc]}

    print("Generation finished")

    for x in pos_dic:

        print(x, len(pos_dic[x]))

    return pos_dic

def assign_spots_new(spots_dic, nuc_dic, pos_dic):
    from spot_the_dffierence import get_distance_filter
    nuc_spot_dic = {}
    counter2 = 0
    lentimes = len(spots_dic)

    for time in spots_dic:
        counter2 += 1
        nuc_spot_dic[time] = {}

        dis_filter = get_distance_filter(nuc_dic, time)*3

        counter = 0
        spotslen = len(spots_dic[time])

        for spot in spots_dic[time]:

            print(str(counter2/lentimes*100)[:4]+'\t'+str(counter/spotslen*100)[:4]+'\t\t\t', end='\r')

            counter += 1

            spotx, spoty, spotz = spots_dic[time][spot]

            dis_dic = {}

            for pos in pos_dic[time]:

                nucx, nucy, nucz = pos

                dis = ((nucx - spotx)**2) + ((nucy - spoty)**2) + ((nucz - spotz)**2)
                # dis = ((nucx - spotx)**2) + ((nucy - spoty)**2)     # without z

                dis_dic[dis] = pos

            min_dis = min(list(dis_dic))

            if min_dis > dis_filter:

                continue

            # print(min_dis, dis_dic[min_dis])
            # print(pos_dic[dis_dic[min_dis]])

            try:

                nuc_spot_dic[time][pos_dic[time][dis_dic[min_dis]][0]].append(spot)

            except KeyError:

                nuc_spot_dic[time][pos_dic[time][dis_dic[min_dis]][0]] = [spot]

    return nuc_spot_dic


def assign_spots_new_multi(spots_dic, nuc_dic, pos_dic, t_list, nuc_spots):

    import multiprocessing as mp

    counter2 = 0
    lentimes = len(spots_dic)

    for time in t_list:

        counter2 += 1
        # nuc_spots[time] = {}

        dis_filter = get_distance_filter(nuc_dic, time)*3

        counter = 0
        spotslen = len(spots_dic[time])

        for spot in spots_dic[time]:

            # print(str(counter2/lentimes*100)+'\t'+str(counter/spotslen*100)+'\t\t\t', end='\r')
            counter+=1

            spotx, spoty, spotz = spots_dic[time][spot]

            dis_dic = {}

            for pos in pos_dic[time]:

                nucx, nucy, nucz = pos

                dis = ((nucx - spotx)**2) + ((nucy - spoty)**2) + ((nucz - spotz)**2)
                # dis = ((nucx - spotx)**2) + ((nucy - spoty)**2)     # without z

                dis_dic[dis] = pos

            min_dis = min(list(dis_dic))

            if min_dis > dis_filter:

                continue

            # print(min_dis, dis_dic[min_dis])
            # print(pos_dic[dis_dic[min_dis]])

            try:

                nuc_spots[time][pos_dic[time][dis_dic[min_dis]][0]].append(spot)

            except KeyError:

                nuc_spots[time][pos_dic[time][dis_dic[min_dis]][0]] = [spot]

    # return nuc_spots

def make_times_dic(pos_dic, cores):

    times_dic = {}

    counter = 1

    for time in pos_dic:

        try:

            times_dic[counter].append(time)

        except KeyError:

            times_dic[counter] = [time]

        counter += 1

        if counter > cores:

            counter = 1

    return times_dic


def write_out_spots_dic(nuc_spots, nuc_dic, spot_dic, exp_dir):

    import os

    out_dir = exp_dir + '/time_data'

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    outfile = out_dir + '/position_data.txt'

    with open(outfile, 'w') as o_out:

        o_out.write('time\tnuc\tnum_spots\tnucx\tnucy\tnucz\tspot\tspotx\tspoty\tspotz\n')

        for t in nuc_spots:

            for nuc in nuc_spots[t]:

                num_spots = len(nuc_spots[t][nuc])

                nucx, nucy, nucz = nuc_dic[t][nuc]

                for spot in nuc_spots[t][nuc]:

                    spotx, spoty, spotz = spot_dic[t][spot]

                    out_ls = [str(x) for x in [t, nuc, num_spots, nucx, nucy, nucz, spot, spotx, spoty, spotz]]

                    outstr = '\t'.join(out_ls)+'\n'

                    o_out.write(outstr)


if __name__ == "__main__":

    import argparse
    from spot_the_dffierence import get_experiments_folder, get_files, get_spot_positions, \
        build_index_dic, get_nuc_positions, get_distance_filter, check_nuc_spots
    from pprint import pprint as pp
    import multiprocessing as mp

    parser = argparse.ArgumentParser()

    parser.add_argument("exp_dir")


    parser.add_argument("--i", default=20, help="time interval in seconds between points. Default: 20 seconds", type=float)

    args = parser.parse_args()

    folders_dic = get_experiments_folder(args.exp_dir)

    pp(folders_dic)

    spot_files = get_files(folders_dic['spot'])
    nuc_files = get_files(folders_dic['nuc'])

    spot_pos = get_spot_positions(spot_files)
    nuc_pos = get_nuc_positions(nuc_files)
    pp(spot_files)

    nuc_axis_dic = import_axis(folders_dic['nuc'], nuc_files)

    # plot_axis_at_time(nuc_axis_dic, 90, nuc_pos)

    pos_dic = make_positions_dic(nuc_axis_dic, nuc_pos)

    nuc_spots = assign_spots_new(spot_pos, nuc_pos, pos_dic)

    # times_dic = make_times_dic(spot_pos, 7)
    #
    # pp(times_dic)
    #
    # jobs = []
    # manager = mp.Manager()
    # nuc_spots = manager.dict()

    # for time in nuc_pos:
    #
    #     for nuc in nuc_pos[time]:
    #
    #         nuc_spots[time] = {nuc:  []}
    #
    # print(nuc_spots)
    #
    # for tt in times_dic:
    #     t_list = times_dic[tt]
    #     jobs.append(mp.Process(target=assign_spots_new_multi,
    #                            args=[spot_pos, nuc_pos, pos_dic, t_list, nuc_spots]))
    #
    # for j in jobs:
    #
    #     j.start()
    #
    # for i in jobs:
    #
    #     i.join()
    #
    # print(nuc_spots)
    check_nuc_spots(nuc_spots)

    write_out_spots_dic(nuc_spots, exp_dir=args.exp_dir, spot_dic=spot_pos, nuc_dic=nuc_pos)
