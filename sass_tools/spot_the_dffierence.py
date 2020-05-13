"""

Better version of all_in_good_time.

Import and assign spots to a nucleous over time.

"""

def get_experiments_folder(exp_dir):

    import os


    folder_dic = {}


    for x in os.listdir(exp_dir):
        # print(x, 'background' in x.lower())
        if 'nucleus' in x.lower() or 'cells' in x.lower():

            folder_dic['nuc'] = exp_dir+'/'+x

        elif 'spot' in x.lower():

            folder_dic['spot'] = exp_dir + '/' + x

        elif 'background' in x.lower():

            folder_dic['backG'] = exp_dir + '/' + x

        else:

            continue

    return folder_dic

def get_files(item_dir):

    import os

    files_to_get = ["_Position.", "_Intensity_", "Ellipsoid_Axis"]

    files_ls = [x for x in os.listdir(item_dir) if ".csv" in x]

    files_dic = {}

    for get_file in files_to_get:

        for file in files_ls:

            if get_file in file:
                print(get_file)
                try:

                    files_dic[get_file].append(item_dir+'/'+file)

                except KeyError:

                    files_dic[get_file] = [item_dir+'/'+file]

    import pprint

    pprint.pprint(files_dic)


    return files_dic


def build_index_dic(ls):

    dic = {}

    for x in range(0, len(ls)):

        dic[ls[x]] = x

    return dic


def get_nuc_positions(nuc_files_dic):

    from pprint import pprint
    nuc_dic = {}

    if len(nuc_files_dic['_Position.']) != 1:

        for x in nuc_files_dic['_Position.']:

            print('--->', x)

            if '_Track_Position' in x:
                print('++++>', x)

                nuc_files_dic['_Position.'].remove(x)

        if len(nuc_files_dic['_Position.']) != 1:

            exit("More than one positions file ***detected!!")

    with open(nuc_files_dic['_Position.'][0]) as o_pos:

        data_found = False

        for line in o_pos:

            split_line = line.strip().split(',')

            if split_line[0].startswith('Position X') and not data_found:

                data_found = True

                column_index = build_index_dic(split_line)

                continue

            elif not data_found:

                continue

            nuc_id = split_line[column_index['TrackID']]
            nuc_time = float(split_line[column_index['Time']])
            pos_x = float(split_line[column_index['Position X']])
            pos_y = float(split_line[column_index['Position Y']])
            pos_z = float(split_line[column_index['Position Z']])

            try:

                nuc_dic[nuc_time][nuc_id] = (pos_x, pos_y, pos_z)

            except KeyError:

                nuc_dic[nuc_time] = {nuc_id: (pos_x, pos_y, pos_z)}

    pprint(nuc_dic)
    return nuc_dic


def get_spot_positions(spot_files_dic):

    from pprint import pprint
    spot_dic = {}

    if len(spot_files_dic['_Position.']) != 1:

        for x in spot_files_dic['_Position.']:

            print('--->', x)

            if '_Track_Position' in x:
                print('++++>', x)

                spot_files_dic['_Position.'].remove(x)

        if len(spot_files_dic['_Position.']) != 1:

            exit("More than one positions file ***detected!!")

    with open(spot_files_dic['_Position.'][0]) as o_pos:

        data_found = False

        for line in o_pos:

            split_line = line.strip().split(',')

            if split_line[0].startswith('Position X') and not data_found:

                data_found = True

                column_index = build_index_dic(split_line)

                continue

            elif not data_found:

                continue

            spot_id = split_line[column_index['ID']]
            spot_time = float(split_line[column_index['Time']])
            pos_x = float(split_line[column_index['Position X']])
            pos_y = float(split_line[column_index['Position Y']])
            pos_z = float(split_line[column_index['Position Z']])

            try:

                spot_dic[spot_time][spot_id] = (pos_x, pos_y, pos_z)

            except KeyError:

                spot_dic[spot_time] = {spot_id: (pos_x, pos_y, pos_z)}

    # pprint(spot_dic)
    return spot_dic

def assign_spots(spots_dic, nuc_dic):

    nuc_spot_dic = {}

    for time in spots_dic:

        nuc_spot_dic[time] = {}
        dis_filter = get_distance_filter(nuc_dic, time)
        for spot in spots_dic[time]:

            spotx, spoty, spotz = spots_dic[time][spot]

            dis_dic = {}

            for nuc in nuc_dic[time]:

                nucx, nucy, nucz = nuc_dic[time][nuc]

                dis = ((nucx - spotx)**2) + ((nucy - spoty)**2) + ((nucz - spotz)**2)
                # dis = ((nucx - spotx)**2) + ((nucy - spoty)**2)     # without z

                dis_dic[dis] = nuc

            min_dis = min(list(dis_dic))

            if min_dis > dis_filter:

                continue

            try:

                nuc_spot_dic[time][dis_dic[min_dis]].append(spot)

            except KeyError:

                nuc_spot_dic[time][dis_dic[min_dis]] = [spot]

    from pprint import pprint

    pprint(nuc_spot_dic)

    return nuc_spot_dic

def get_distance_filter(nuc_dic, t):

    from numpy import mean

    min_dis_ls = []

    for nuc1 in nuc_dic[t]:

        nuc1x, nuc1y, nuc1z = nuc_dic[t][nuc1]

        dis_dic = {}

        for nuc in nuc_dic[t]:

            if nuc1 == nuc:

                continue

            nucx, nucy, nucz = nuc_dic[t][nuc]

            dis = ((nucx - nuc1x) ** 2) + ((nucy - nuc1y) ** 2) + ((nucz - nuc1z) ** 2)
            # dis = ((nucx - nuc1x) ** 2) + ((nucy - nuc1y) ** 2)     # without Z

            dis_dic[dis] = nuc

        min_dis = min(list(dis_dic))

        min_dis_ls.append(min_dis)

    return mean(min_dis_ls)

def get_dis(p1, p2):

    p1x, p1y, p1z = p1
    p2x, p2y, p2z = p2

    dis = ((p1x - p2x) ** 2) + ((p1y - p2y) ** 2) + ((p1z - p2z) ** 2)
    # dis = ((p1x - p2x) ** 2) + ((p1y - p2y) ** 2)   # without z

    return dis

def minimise_difference(nuc, spot, nuc_dic, spot_dic, t):

    exit()

def fix_assignment(nuc_spots, spot_dic, nuc_dic):
    counter = 0
    for t in nuc_spots:

        dis_filter = get_distance_filter(nuc_dic, t)

        for nuc in nuc_spots[t]:

            if len(nuc_spots[t][nuc]) != 1:
                dis_dic = {}

                print_with_closest_nucs(nuc_spots, t, nuc, nuc_dic, spot_dic)
                counter += 1

                if counter > 10:

                    exit()
                for spot in nuc_spots[t][nuc]:

                    spotx, spoty, spotz = spot_dic[t][spot]
                    dis_dic[spot] = {}

                    for nuc in nuc_dic[t]:

                        nucx, nucy, nucz = nuc_dic[t][nuc]

                        dis = ((nucx - spotx)**2) + ((nucy - spoty)**2) + ((nucz - spotz)**2)
                        # dis = ((nucx - spotx) ** 2) + ((nucy - spoty) ** 2)  # without z

                        dis_dic[spot][dis] = nuc

                    min_dis_ls = sorted(list(dis_dic[spot]))
                    filtered_dis = [i for i in min_dis_ls if i <= dis_filter*1.2]

                    for i in filtered_dis:

                        print(spot, nuc, i, dis_dic[spot][i], dis_filter)

                # exit()

def print_with_closest_nucs(nuc_spots, t, nuc, nuc_dic, spots_dic):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    dis_dic = {}

    for nuc1 in nuc_dic[t]:

        nuc1x, nuc1y, nuc1z = nuc_dic[t][nuc1]

        nucx, nucy, nucz = nuc_dic[t][nuc]

        dis = ((nucx - nuc1x) ** 2) + ((nucy - nuc1y) ** 2) + ((nucz - nuc1z) ** 2)
        # dis = ((nucx - nuc1x) ** 2) + ((nucy - nuc1y) ** 2)     # without Z

        dis_dic[dis] = nuc1

    dis_ls = sorted(list(dis_dic))[:11]

    nucs_x = []
    nucs_y = []
    nucs_z = []

    for i in range(0, len(dis_ls)):

        if nuc == dis_dic[dis_ls[i]]:
            continue
        if i > 10:

            break

        print(nuc_dic[t][dis_dic[dis_ls[i]]], len(nuc_dic[t]))
        nucs_x.append(nuc_dic[t][dis_dic[dis_ls[i]]][0])
        nucs_y.append(nuc_dic[t][dis_dic[dis_ls[i]]][1])
        nucs_z.append(nuc_dic[t][dis_dic[dis_ls[i]]][2])


    spots_x = []
    spots_y = []
    spots_z = []

    for spot in nuc_spots[t][nuc]:

        spots_x.append(spots_dic[t][spot][0])
        spots_y.append(spots_dic[t][spot][1])
        spots_z.append(spots_dic[t][spot][2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(spots_x, spots_y, spots_z)
    ax.scatter(nucs_x, nucs_y, nucs_z)
    ax.scatter(nuc_dic[t][nuc][0], nuc_dic[t][nuc][1], nuc_dic[t][nuc][2], c='black')

    # plt.show()

    # exit()


def check_nuc_spots(nuc_spots, path):

    import matplotlib.pyplot as plt

    x_ls = []

    counts_dic = {}
    total_count = 0

    perc_1 = []
    perc_2 = []
    perc_3_plus = []

    tim_ls = []

    for t in nuc_spots:
        count1 = 0
        count2 = 0
        count3plus = 0
        total2 = 0
        for nuc in nuc_spots[t]:

            total_count += 1
            total2 += 1
            x_ls.append(len(nuc_spots[t][nuc]))

            if len(nuc_spots[t][nuc]) == 1:

                count1 += 1

            elif len(nuc_spots[t][nuc]) == 2:

                count2 += 1

            elif len(nuc_spots[t][nuc]) > 2:

                count3plus += 1

            try:

                counts_dic[len(nuc_spots[t][nuc])] += 1

            except KeyError:

                counts_dic[len(nuc_spots[t][nuc])] = 1
        
        if total2 != 0:

            perc_1.append(count1/total2*100)
            perc_2.append(count2/total2*100)
            perc_3_plus.append(count3plus/total2*100)
            tim_ls.append(t)

        else:
            
            perc_1.append(0)
            perc_2.append(0)
            perc_3_plus.append(0)
            tim_ls.append(t)
            

    # plt.hist(x_ls)
    #
    # plt.show()
    #
    # plt.close()
    #
    plt.plot(tim_ls, perc_1)
    plt.plot(tim_ls, perc_2)
    plt.plot(tim_ls, perc_3_plus)
    plt.axhline(y=95)
    # plt.show()
    plt.savefig(path)
    plt.close()

    for x in counts_dic:

        print(x, counts_dic[x], counts_dic[x]/total_count*100)

def plot_nucs(nuc_dic, spot_dic, t):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    x_ls = []
    y_ls = []
    z_ls = []

    sx_ls = []
    sy_ls = []
    sz_ls = []

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for nuc in nuc_dic[t]:

        x_ls.append(nuc_dic[t][nuc][0])
        y_ls.append(nuc_dic[t][nuc][1])
        z_ls.append(nuc_dic[t][nuc][2])

    for spot in spot_dic[t]:

        sx_ls.append(spot_dic[t][spot][0])
        sy_ls.append(spot_dic[t][spot][1])
        sz_ls.append(spot_dic[t][spot][2])

    ax.scatter(x_ls, y_ls, z_ls, c=z_ls)
    ax.scatter(sx_ls, sy_ls, sz_ls, c=sz_ls, marker='x')
    plt.show()
    plt.close()

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("exp_dir")


    parser.add_argument("--i", default=20, help="time interval in seconds between points. Default: 20 seconds", type=float)

    args = parser.parse_args()

    folders_dic = get_experiments_folder(args.exp_dir)

    nuc_files = get_files(folders_dic['nuc'])
    spot_files = get_files(folders_dic['spot'])

    nuc_pos = get_nuc_positions(nuc_files)

    spot_pos = get_spot_positions(spot_files)

    # nuc_spots = assign_spots(spot_pos, nuc_pos)

    # fix_assignment(nuc_spots, spot_pos, nuc_pos)
    # check_nuc_spots(nuc_spots)

    for x in range(90, 100):

        plot_nucs(nuc_pos, spot_pos, x)