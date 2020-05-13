def get_imaris_files(folder_path, folder_list, file_list):
    """
    Searches through the embyro folder (folder_path) and gets all the relevant files to run the data analysis.

    Folder names are in the folder_list varaible and the files will be passed as the file_list variable.

    Data is returned as a dictionary with the top level folders as keys containing the te relevant internal files.
    """

    import os

    file_dic = {}

    for folder in os.listdir(folder_path):
        for x in folder_list:
            if x in folder.lower():
                file_dic[x] = {'path': folder, 'files': []}
                files_in_folder = os.listdir(os.path.join(folder_path, folder))

                for file in files_in_folder:
                    for y in file_list:
                        if y.lower() in file.lower():
                            file_dic[x]['files'].append(file)

    return file_dic

def get_positions(positions_file):

    """
    ======
    inputs
    ======

    positions_file: posiiton file from imaris containing the x, y and z coordinates of an imaris object

    ======
    output
    ======

    pos_dic: dictionary of positions
    
    """

    from pprint import pprint

    pos_dic = {}

    with open(positions_file) as o_pos:

        data_found = False

        for line in o_pos:

            split_line = line.strip().split(',')

            if split_line[0].startswith('Position X') and not data_found:

                data_found = True

                column_index = build_index_dic(split_line)

                continue

            elif not data_found:

                continue

            pos_id = split_line[column_index['ID']]
            pos_time = float(split_line[column_index['Time']])
            pos_x = float(split_line[column_index['Position X']])
            pos_y = float(split_line[column_index['Position Y']])
            pos_z = float(split_line[column_index['Position Z']])

            try:

                pos_dic[pos_time][pos_id] = (pos_x, pos_y, pos_z)

            except KeyError:

                pos_dic[pos_time] = {pos_id: (pos_x, pos_y, pos_z)}
  
    return pos_dic

def build_index_dic(ls):

    dic = {}

    for x in range(0, len(ls)):

        dic[ls[x]] = x

    return dic

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

def assign_nuc_axis_dic(axis_dic, nuc_dic):

    """
    Uses the axis of the nucleus to generate many points along the axis for each nucleus.

    """

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


def import_axis(elipsoid_files, folder_path):
    import os
    files_to_get = ['A', 'B', 'C']

    axis_dic = {}

    for x in elipsoid_files:

        for filex in files_to_get:

            if 'Ellipsoid_Axis_'+filex+'.' in x:
        
                with open(os.path.join(folder_path, x)) as openX:

                    data_found = False

                    for line in openX:

                        split_line = line.strip().split(',')

                        if split_line[0].startswith('Ellipsoid Axis '+filex+' X') and not data_found:

                            data_found = True

                            column_index = build_index_dic(split_line)
                            # pp(column_index)
                            continue

                        elif not data_found:

                            continue

                        nuc_id = split_line[column_index['ID']]

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

def assign_spots_new(spots_dic, nuc_dic, pos_dic):
    
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


def check_nuc_spots(nuc_spots, path):

    import matplotlib.pyplot as plt

    pp(nuc_spots)


def write_out_spots_dic(nuc_spots, nuc_dic, spot_dic, exp_dir):

    import os

    out_dir = exp_dir + '/time_data'

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    outfile = out_dir + '/position_data.txt'

    with open(outfile, 'w') as o_out:

        o_out.write('time\tnuc\tnum_spots\tnucx\tnucy\tnucz\tspot\tspotx\tspoty\tspotz\n')

        for t in nuc_dic:

            for nuc in nuc_dic[t]:
                try:
                    num_spots = len(nuc_spots[t][nuc])

                    nucx, nucy, nucz = nuc_dic[t][nuc]

                    for spot in nuc_spots[t][nuc]:

                        spotx, spoty, spotz = spot_dic[t][spot]

                        out_ls = [str(x) for x in [t, nuc, num_spots, nucx, nucy, nucz, spot, spotx, spoty, spotz]]

                        outstr = '\t'.join(out_ls)+'\n'

                        o_out.write(outstr)

                except KeyError:

                    num_spots = 0

                    nucx, nucy, nucz = nuc_dic[t][nuc]
                    spot = 'na'
                    spotx, spoty, spotz = 0, 0, 0

                    out_ls = [str(x) for x in [t, nuc, num_spots, nucx, nucy, nucz, spot, spotx, spoty, spotz]]

                    outstr = '\t'.join(out_ls)+'\n'

                    o_out.write(outstr)

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
            if spot == 'na':
                outfile.write(line.replace('\n', f'\tna\tna\tna\n'))
                continue
            # print(type(spot), spot)
            for variable in spot_dic[t][spot]:

                for chan in spot_dic[t][spot][variable]:

                    value = spot_dic[t][spot][variable][chan]
                    # print(t, spot, channel, variable, value)
                    outfile.write(line.replace('\n', f'\t{variable}\t{chan}\t{value}\n'))


def get_spot_info(exp_dir, spot_intensity_files):
    import os
    spot_dic = {}

    for file in spot_intensity_files:
        
        with open(os.path.join(exp_dir, file)) as o_file:
            data_found = False

            for line in o_file:
                split_line = line.strip().split(',')
              

                if 'ID' in line and not data_found:
                    data_found = True
                    tdex = split_line.index('Time')
                    sdex = split_line.index('ID')
                    print(file)
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

def get_mid(intense_path):

    import numpy as np
    import matplotlib.pyplot as plt
    import os

    curvey = 2
    x_ls = []
    y_ls = []
    x2 = []
    y2 = []
    x0 = []
    y0 = []
    x1 = []
    y1 = []
    image_file = os.path.split(intense_path)[0]+'/mid_line_image.png'

    with open(intense_path) as o_intense:

        first_line = True

        for line in o_intense:

            split_line = line.strip().split('\t')

            if first_line:

                first_line = False

                x_dex = split_line.index('spotx')
                y_dex = split_line.index('spoty')
                num_spots = split_line.index('num_spots')

                nucx_dex = split_line.index('nucx')
                nucy_dex = split_line.index('nucy')

                continue
            if float(split_line[num_spots]) == 0:

                x0.append(float(split_line[nucx_dex]))
                y0.append(float(split_line[nucy_dex]))
                continue

            elif float(split_line[num_spots]) == 1:

                x1.append(float(split_line[nucx_dex]))
                y1.append(float(split_line[nucy_dex]))
                continue
            
            elif float(split_line[num_spots]) > 1:

                x2.append(float(split_line[nucx_dex]))
                y2.append(float(split_line[nucy_dex]))

            x_ls.append(float(split_line[x_dex]))
            y_ls.append(float(split_line[y_dex]))

    f = np.poly1d(np.polyfit(x_ls, y_ls, curvey))
    # x = range(0, int(max(x_ls)))
    x_new = np.linspace(min(x_ls)-10, max(x_ls)+10, 300)
    y_new = f(x_new)

    plt.scatter(x0, y0, c='black')
    plt.scatter(x1, y1, c='blue')
    plt.scatter(x2, y2, c='green')
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

    return min(dis_ls)**.5 * trans, trans

def return_mid_line(mid_line, path):

    with open(path+'/mid_line.txt', 'w') as o_mid:

        o_mid.write(f'x\ty\n')

        for x in range(0, len(mid_line[0])):

            o_mid.write(f'{mid_line[0][x]}\t{mid_line[1][x]}\n')

def distance_to_midline(intense_path, mid_line):

    import os
    split_dic = {}
    head_line = ''
    time_ls = []
    chan_ls = []
    intent_type_ls = []

    outfile = open(os.path.split(intense_path)[0]+'/all_data_with_distance.txt', 'w')
    points_file = open(os.path.split(intense_path)[0] + '/mid_points2.txt', 'w')
    with open(intense_path) as o_intense:

        first_line = True

        for line in o_intense:

            split_line = line.strip().split('\t')

            if first_line:

                outfile.write(line.replace('\n', '\tdis_mid\tabove_line\n'))
                # print(line)
                # tdex = split_line.index('time')
                # ndex = split_line.index('nuc')
                ndexX = split_line.index('nucx')
                ndexY = split_line.index('nucy')
                # chandex = split_line.index('channel')
                # intent_type = split_line.index('variable')
                # value = split_line.index('value')
                first_line = False
                continue

            nucx, nucy = split_line[ndexX], split_line[ndexY]

            dis_mid, above_line = get_distance_mid(mid_line, float(nucx), float(nucy), points_file)

            outfile.write(line.replace('\n', f'\t{dis_mid}\t{above_line}\n'))


def get_inter_nuclaer_distance(nuc_dic):

    """
    nuc_dic: dictionary of nuclear positions

    loops throughall nuclei and finds the ddistance to it's closest neighbour, the mean minimal distance is then used as an estiamte of cell size.

    returns mean minimum nuclear distances

    """

    from numpy import mean

    min_dis_ls = []
    for nuc in nuc_dic:
        dis_ls = []
        # print(nuc_dic[nuc])
        for nuc2 in nuc_dic:
            if nuc == nuc2:
                # print(nuc, nuc2)
                continue
            
            nucx1, nucy1, nucz1 = nuc_dic[nuc]
            nucx2, nucy2, nucz2 = nuc_dic[nuc2]

            distance_calculated = (nucx1-nucx2)**2 + (nucy1-nucy2)**2 + (nucz1-nucz2)**2
            # print(distance_calculated)
            dis_ls.append(distance_calculated)
        # print(dis_ls)
        min_dis_ls.append(min(dis_ls)**.5)

    return mean(min_dis_ls)

def get_bins(intense_tsv_pd_obj, bin_size):

    """
    Returns distance bins from the mid_line
    """

    from numpy import max, min
    import pandas as pd

    max_abs_dis = abs(max(intense_tsv_pd_obj['dis_mid']))
    min_abs_dis = abs(min(intense_tsv_pd_obj['dis_mid']))

    target_dis = max_abs_dis

    if max_abs_dis < min_abs_dis:

        target_dis = min_abs_dis

    counter = 0
    counter2 = 0
    counter3 = 0

    bins_ls = []
    bin_numbers = []

    while counter <= target_dis:

        bins_ls.append(counter)
        bins_ls.append(counter2)

        bin_numbers.append(counter3)
        bin_numbers.append(counter3*-1)

        counter += bin_size
        counter2 -= bin_size
        counter3 += 1

    bins_ls = sorted(list(set(bins_ls)))
    bins_numbers = sorted(list(set(bin_numbers)))
    bins_numbers = [x for x in bins_numbers if x != 0]

    # print(bins_ls)
    # print(bins_numbers)
    return bins_ls, bins_numbers

def plot_data(intense_tsv_pd_obj, bins_ls, bin_numbers, bin_size, bin_factor, data_path):

    import matplotlib.pyplot as plt


    # print(len(bins_ls))

    counts = []

    bin_axis = []

    bins_out = open(data_path+'/bins.txt', 'w')
    bins_out.write(f'bin\tcount\tnum_spots\n')
    for x in range(0, len(bins_ls)-1):
        bin_axis.append(f'{round(bin_ls[x]/bin_size*bin_factor)} : {round(bin_ls[x+1]/bin_size*bin_factor)}')

    for x in range(0, len(bins_ls) -1):
        # print(x)
        count_out = len(set(intense_tsv_pd_obj[(intense_tsv_pd_obj['dis_mid'] > bins_ls[x]) & (intense_tsv_pd_obj['dis_mid'] < bins_ls[x + 1]) & (intense_tsv_pd_obj['num_spots'] == 1)]['nuc']))
        counts.append(count_out)
        bins_out.write(f'{bin_axis[x]}\t{count_out}\t1\n')
    counts2 = []
    plt.bar(x=bin_axis, height=counts, color='orange')

    for x in range(0, len(bins_ls) - 1):
        count_out = len(set(intense_tsv_pd_obj[(intense_tsv_pd_obj['dis_mid'] > bins_ls[x]) & (intense_tsv_pd_obj['dis_mid'] < bins_ls[x + 1]) & (intense_tsv_pd_obj['num_spots'] >= 2)]['nuc']))
        counts2.append(count_out)
        bins_out.write(f'{bin_axis[x]}\t{count_out}\t2\n')
    print(bin_numbers)
    
    print(bins_ls)
    plt.bar(x=bin_axis, height=counts2, color='red')
    plt.xlabel('Distance to midline (Cell widths)')
    plt.ylabel('Number of Nuclei')
    # plt.show()
    plt.xticks(rotation=90)
    plt.savefig(data_path+'/binned.png')
    plt.close()
    bins_out.close()

def filter_spots(spot_dic, filter_dis):
    from pprint import pprint as pp
    import matplotlib.pyplot as plt
    # pp(spot_dic)
    dis_ls = []
    spots_to_remove = []
    for time in spot_dic:
        for spot1 in spot_dic[time]:
            # print(spot1, spot_dic[time][spot1])
            x1, y1, z1 = spot_dic[time][spot1]
            temp_dis_list = []
            for spot2 in spot_dic[time]:

                x2, y2, z2 = spot_dic[time][spot2]

                if spot1 == spot2:
                    continue
                
                distance_temp = ((x1 - x2)**2) + ((y1 - y2)**2) + ((z1 - z2)**2)
                if distance_temp < filter_dis**2:
                    spots_to_remove.append(spot2)
                    print(spot1, spot_dic[time][spot1], spot2, spot_dic[time][spot2])
                temp_dis_list.append(distance_temp)
            dis_ls.append(min(temp_dis_list)**.5)
    
        for spot in set(spots_to_remove):
            del spot_dic[time][spot]

    return spot_dic

def reassign_spots(pos_dic, nuc_spots, nuc_with_axis):
    from numpy import mean
    # pprint.pprint(nuc_spots)
    spots_dic = pos_dic['spots']
    for t in nuc_spots:
        nuc_counts = {}
        for nuc in nuc_spots[t]:
            try:
                nuc_counts[len(nuc_spots[t][nuc])].append(nuc)
            
            except KeyError:
                nuc_counts[len(nuc_spots[t][nuc])]= [nuc]

        for x in [i for i in nuc_counts if i > 2]:
            for nuc in nuc_counts[x]:
                spot_pos = {}
                for spot in nuc_spots[t][nuc]:
                
                    spotx, spoty, spotz = spots_dic[t][spot]

                    dis_dic = {}
                    # print(list(pos_dic))
                    for pos in nuc_with_axis[t]:

                        nucx, nucy, nucz = pos

                        dis = ((nucx - spotx)**2) + ((nucy - spoty)**2) + ((nucz - spotz)**2)
                        # dis = ((nucx - spotx)**2) + ((nucy - spoty)**2)     # without z

                        dis_dic[dis] = pos

                    temp_ls = []

                    sorted_dis = sorted(list(dis_dic))
                    counter = 0 
                    for xi in sorted_dis:
                        if counter == 20:
                            break
                        try:
                            spot_pos[spot][nuc_with_axis[t][dis_dic[xi]][0]].append(xi**.5)

                        except KeyError:
                            try:
                                spot_pos[spot][nuc_with_axis[t][dis_dic[xi]][0]] = [xi**.5]
                            except KeyError:
                                spot_pos[spot] = {nuc_with_axis[t][dis_dic[xi]][0]: [xi**.5]}
                        # temp_ls.append((nuc_with_axis[t][dis_dic[xi]], xi**.5))
                        counter += 1

                     
                    for spot in spot_pos:
                        print('---- spot ------')
                        for nuc in spot_pos[spot]:
                            
                            print(nuc, mean(spot_pos[spot][nuc]), spot_pos[spot][nuc])

	
    exit()
if __name__ == "__main__":

    import argparse
    import os
    from pprint import pprint as pp
    import pandas as pd
  
    parser = argparse.ArgumentParser()
    parser.add_argument('folder_path')
    parser.add_argument('-b', type=float, default=3, help='bin_size in cell widths')
    parser.add_argument('-ex', type=float, default=1, help='signal channel of exonic spots')
    parser.add_argument('-in', type=float, default=2, help='signal channel of intronic spots')
    parser.add_argument('-ti', type=float, default=20, help='time interval')
    parser.add_argument('-sf', type=float, default=.5, help='spot distance filter')
    args = parser.parse_args()
    
    folder_list = ['spots', 'nucl']
    file_list = ["_Position.", "_Intensity_", "Ellipsoid_Axis"]

    files_dic = get_imaris_files(args.folder_path, folder_list, file_list)
    pp(files_dic)
    positions_dic = {}

    for folder in folder_list:

        folder_pos = [x for x in files_dic[folder]['files'] if "_Position." in x]
        
        if len(folder_pos) != 1:

            exit(f'wrong number of positions files detected: {folder} || {folder_pos}')

        positions_dic[folder] = get_positions(os.path.join(args.folder_path, files_dic[folder]['path'], folder_pos[0]))
    len1 = len(positions_dic['spots'][1])
    positions_dic['spots'] = filter_spots(positions_dic['spots'].copy(), args.sf)
    len2 = len(positions_dic['spots'][1])
    print(len1, len2, len1-len2)
    
    nuc_axis_dic = import_axis([x for x in files_dic[folder]['files'] if "Ellipsoid_Axis" in x], os.path.join(args.folder_path, files_dic[folder]['path']))

    nuc_with_axis = assign_nuc_axis_dic(nuc_axis_dic, positions_dic['nucl'])

    pp(nuc_with_axis)

    nuc_spots = assign_spots_new(positions_dic['spots'], positions_dic['nucl'], nuc_with_axis)
    # nuc_spots = reassign_spots(positions_dic, nuc_spots, nuc_with_axis)
    # exit()
    check_nuc_spots(nuc_spots, args.folder_path)

    write_out_spots_dic(nuc_spots, exp_dir=args.folder_path, spot_dic=positions_dic['spots'], nuc_dic=positions_dic['nucl'])

    pos_file = args.folder_path+'/time_data/position_data.txt'

    spot_file_path = os.path.join(args.folder_path, files_dic['spots']['path'])

    spotted_dic = get_spot_info(spot_file_path, [x for x in files_dic['spots']['files']if "_Intensity_" in x])

    assign_intensities(args.folder_path, spotted_dic, pos_file)

    intense_path = args.folder_path+ '/time_data/position_data-intense.txt'

    mid_line = get_mid(intense_path)
    return_mid_line(mid_line, args.folder_path)

    distance_to_midline(intense_path, mid_line)

    inter_nuc_dis = get_inter_nuclaer_distance(positions_dic['nucl'][1.0])
    bin_size = inter_nuc_dis*args.b
    pandas_frame = os.path.split(intense_path)[0]+'/all_data_with_distance.txt'
    bin_ls, bin_numbers = get_bins(pd.read_table(pandas_frame, sep='\t', header=0), bin_size)
    plot_data(pd.read_table(pandas_frame, sep='\t', header=0), bin_ls, bin_numbers, bin_size, args.b, args.folder_path)
