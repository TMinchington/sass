"""
SpotMe_VD2 with improved midline detection.

Bug fix:

-  counts from anterior position instead of random. (This may fail in some circumstances though runs on all test flies)
"""


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


def build_index_dic(ls):

    dic = {}

    for x in range(0, len(ls)):

        dic[ls[x]] = x

    return dic

def straight_mid(points, min_x, max_x, max_y):

    from numpy import ones,vstack, linspace
    from numpy.linalg import lstsq
    # points = [(1,5),(3,4)]
    x_coords, y_coords = zip(*points)
    A = vstack([x_coords,ones(len(x_coords))]).T
    m, c = lstsq(A, y_coords, rcond=None)[0]
    print("Line Solution is y = {m}x + {c}".format(m=m,c=c))

    new_x = linspace(min_x, max_x, 300)
    new_y = [x*m+c for x in new_x]

    return new_x, new_y

    
def get_number_of_neighbours(nucDic, interNucDis):

    import matplotlib.pyplot as plt
    from pprint import pprint as pp
    neighbourDic = {}

    for nuc1 in nucDic:
        for nuc2 in nucDic:
            if nuc1 == nuc2:
                continue

            dis2Nuc = get_distance(nucDic[nuc1], nucDic[nuc2])
            
            if dis2Nuc <= 2 * interNucDis:
                # print(dis2Nuc, interNucDis)

                try:
                    neighbourDic[nuc1][nuc2] = 0

                except KeyError:
                    neighbourDic[nuc1] = {nuc2: 0}
    
    pp(neighbourDic)
    count_ls = []
    count_ls2 = []
    xls = []
    yls = []

    for neighbour in neighbourDic:
        neighbourNum = len(neighbourDic[neighbour])
        count_ls.append(neighbourNum)
        if neighbourNum >4:
            count_ls2.append(1)

        else:
            count_ls2.append(5)
        xls.append(nucDic[neighbour][0])
        yls.append(nucDic[neighbour][1])

    # plt.scatter(range(0, len(count_ls)), count_ls)
    plt.hist(count_ls)
    plt.show()
    plt.close()

    plt.scatter(xls, yls, c=count_ls)
    plt.show()
    plt.close()

    plt.scatter(xls, yls, c=count_ls2)
    plt.show()


    return neighbourDic


def find_edge_cells(nucPosDic):
    from numpy import mean, std, var, percentile
    import numpy as np
    interNucDis = get_inter_nuclaer_distance(nucPosDic)
    neighbourDic = {}

    for nuc1 in nucPosDic:
        for nuc2 in  nucPosDic:
            if nuc1 == nuc2:
                continue
            
            dis2Nuc = get_distance(nucPosDic[nuc1], nucPosDic[nuc2])
            
            if dis2Nuc < interNucDis*6:
                # print(dis2Nuc, interNucDis)
                try:
                    neighbourDic[nuc1][nuc2] = dis2Nuc

                except KeyError:
                    neighbourDic[nuc1] = {nuc2: dis2Nuc}


    score_ls = []
    xls = []
    yls = []

    for nuc in nucPosDic:
        temp_score = []
        for neighbour in neighbourDic[nuc]:
            positionScore = getVectorScore(nucPosDic[nuc], nucPosDic[neighbour], neighbourDic[nuc][neighbour])
            temp_score.append(positionScore)

        score_ls.append(var(temp_score))
        xls.append(nucPosDic[nuc][0])
        yls.append(nucPosDic[nuc][1])

    import matplotlib.pyplot as plt
    score2_ls = []

    upperPerc = percentile(score_ls, 95)
    lowerPerc = percentile(score_ls, 5)
    minX, maxX, minY, maxY = min(xls), max(xls), min(yls), max(xls)
    interNucDis*=1.2

    xls2 = []
    yls2 = []

    for score, x, y in zip(score_ls, xls, yls):
        if minX+interNucDis > x or x > maxX-interNucDis or minY+interNucDis > y or y > maxY - interNucDis:
            # print()
        
            score2_ls.append(3)
        elif score > upperPerc or score < lowerPerc:
            score2_ls.append(8)
            xls2.append(x)
            yls2.append(y)
        else:
            score2_ls.append(1)

    # print(score_ls)
    # score2

    from numpy import ones,vstack, linspace
    from numpy.linalg import lstsq

    A = vstack([xls2,ones(len(yls2))]).T
    m, c = lstsq(A, yls2, rcond=None)[0]
    print("Line Solution is y = {m}x + {c}".format(m=m,c=c))

    new_x = linspace(minX, maxX, 300)
    new_y = [x*m+c for x in new_x]

    plt.scatter(xls, yls, c='black')
    plt.scatter(xls, yls, c=score2_ls, s=0.8)
    # plt.scatter(xls, yls, c=score_ls)
    plt.Line2D(new_x, new_y)
    plt.show()
    plt.close()

    plt.hist(score_ls)
    plt.show()
    plt.close()


def getVectorScore(nuc1, nuc2, distanceNeighbour):

    from math import atan2

    x1, y1, _ = nuc1
    x2, y2, _ = nuc2

    xd = x1-x2
    yd = y1-y2

    vectorScore = atan2(xd, yd)/distanceNeighbour

    return vectorScore


def get_distance(nuc1pos, nuc2pos):

    x1, y1, z1 = nuc1pos
    x2, y2, z2 = nuc2pos

    distance = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**.5

    return distance


def easy_edges(nucPosDic):
    import pandas as pd
    import numpy as np

    x_ls = []
    y_ls = []
    nuc_ls = []
    for nuc in nucPosDic:
        x, y, _ = nucPosDic[nuc]
        x_ls.append(x)
        y_ls.append(y)
        nuc_ls.append(nuc)

    nuc_frame = pd.DataFrame({'nuc':nuc_ls, 'x':x_ls, 'y':y_ls})
    nuc_frame['xbins'] = pd.cut(nuc_frame.x, 20)
    nuc_frame['ybins'] = pd.cut(nuc_frame.y, 20)

    minXnuc = []
    maxXnuc = []
    minYnuc = []
    maxYnuc = []

    uniqueYbins = set(nuc_frame.ybins)
    uniqueXbins = set(nuc_frame.xbins)

    for yBin in uniqueYbins:
        tempFrame = nuc_frame.loc[nuc_frame.ybins == yBin]

        minX = min(tempFrame.x)

        nuc_ls = list(tempFrame.loc[nuc_frame.x == minX].nuc)
        # print(nuc_ls)
        
        if len(nuc_ls) != 1:
            error_str = '\n'+ "\n".join(nuc_ls)
            exit(f'Too many nucs :{error_str}')
        
        minXnuc.append(nuc_ls[0])

        maxX = max(tempFrame.x)

        nuc_ls = list(tempFrame.loc[nuc_frame.x == maxX].nuc)

        if len(nuc_ls) != 1:
            error_str = '\n'+ "\n".join(nuc_ls)
            exit(f'Too many nucs :{error_str}')
        
        maxXnuc.append(nuc_ls[0])

    for xBin in uniqueXbins:
        tempFrame = nuc_frame.loc[nuc_frame.xbins == xBin]

        minY = min(tempFrame.y)

        nuc_ls = list(tempFrame.loc[nuc_frame.y == minY].nuc)
        # print(nuc_ls)
        
        if len(nuc_ls) != 1:
            error_str = '\n'+ "\n".join(nuc_ls)
            exit(f'Too many nucs :{error_str}')
        
        minYnuc.append(nuc_ls[0])

        maxY = max(tempFrame.y)

        nuc_ls = list(tempFrame.loc[nuc_frame.y == maxY].nuc)

        if len(nuc_ls) != 1:
            error_str = '\n'+ "\n".join(nuc_ls)
            exit(f'Too many nucs :{error_str}')
        
        maxYnuc.append(nuc_ls[0])


    edge_nuncs = set(minYnuc+maxYnuc+maxXnuc+minXnuc)

    minX, maxX, minY, maxY = min(nuc_frame.x), max(nuc_frame.x), min(nuc_frame.y), max(nuc_frame.y)
    interNucDis = get_inter_nuclaer_distance(nucPosDic)

    edge_nucs2 = []

    for nuc in edge_nuncs:
        x, y, _ = nucPosDic[nuc]
        
        if minX+interNucDis*1.5 > x or x > maxX-interNucDis*1.5 or minY+interNucDis*1.5 > y or y > maxY - interNucDis*1.5:
            continue
        
        
        # print(temp_dis_ls)
        # if min(temp_dis_ls) > 1.5*interNucDis:
        #     continue
            
        edge_nucs2.append(nuc)

    for nuc in set(edge_nucs2):
        # print()
        temp_dis_dic = {}
        for nuc2 in edge_nucs2:
                
            if nuc == nuc2:
                continue

            x, y, _ = nucPosDic[nuc]
            x2, y2, _ = nucPosDic[nuc2]
            temp_dis = ((x-x2)**2 + (y-y2)**2)**.5

            
            temp_dis_dic[temp_dis] = (nuc, x2, y2)

        # print('>>>', min(temp_dis_dic), interNucDis, x, y, '|', nuc, temp_dis_dic[min(temp_dis_dic)])
        if min(temp_dis_dic) > 4*interNucDis:
            # print('>', nuc, temp_dis_dic[min(temp_dis_dic)])
            edge_nucs2.pop(edge_nucs2.index(nuc))
    

    mid_cloud, anteriorPosition = get_mid_cloud(nucPosDic, edge_nucs2, interNucDis)

    mid_points = [(x, y) for x, y in zip(mid_cloud[0], mid_cloud[1])]

    midline = straight_mid(mid_points, minX, maxX, maxY)
    # print(midline)
    plot_nucs(nucPosDic, edge_nucs2, mid_cloud, midline)
    x_new, y_new = midline
    f = np.poly1d(np.polyfit(x_new, y_new, 1))


    return (x_new, y_new, f), interNucDis, anteriorPosition


def get_mid_cloud(nucPosDic, edge_nucs, interNucDis):

    from numpy import mean
    cloudx = []
    cloudy = []

    cloudx2 = []
    cloudy2 = []

    for nuc1 in edge_nucs:
        x1, y1, _ = nucPosDic[nuc1]
        for nuc2 in edge_nucs:
            if nuc1 == nuc2:
                continue
        
            x2, y2, _ = nucPosDic[nuc2]
            mx = (x1+x2)/2
            my = (y1+y2)/2

            # print(f'{x1, y1} | {x2, y2} | {mx, my}')

            cloudx.append(mx)
            cloudy.append(my)

    pos_remove = []

    for clX, clY in zip(cloudx, cloudy):
        outofrange =True
        for nuc in edge_nucs:
            checkx, checky, _ = nucPosDic[nuc]

            if ((clX-checkx)**2 + (clY-checky)**2)**.5 < interNucDis*5:
                # print(clX, clY, '\t', checkx, checky)
                outofrange = False
                continue

        if outofrange:
            cloudx2.append(clX)
            cloudy2.append(clY)
    # print('XC', cloudx2, cloudy2)


    bigPair = get_most_distant(nucPosDic)

    anteriorPosition = find_anterior(bigPair, interNucDis, nucPosDic)

    closest_x, closest_y = get_10_closest(anteriorPosition, nucPosDic)

    return (cloudx2+closest_x*3, cloudy2+closest_y*3), anteriorPosition


def distance_to_midline(intense_path, mid_line, inter_nuc_dis, anteriorPosition):

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

                outfile.write(line.replace('\n', '\tdis_mid\tabove_line\tdis_along_line\tdis_along_mid_in_nuclei\n'))
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
            mid_start = get_mid_start(mid_line, anteriorPosition)
            dis_mid, above_line, dis_along_mid = get_distance_mid(mid_line, float(nucx), float(nucy), points_file, mid_start)

            outfile.write(line.replace('\n', f'\t{dis_mid}\t{above_line}\t{dis_along_mid}\t{dis_along_mid/inter_nuc_dis}\n'))

def get_mid_start(mid_line, anteriorPosition):

    aX, aY = anteriorPosition
    xList, yList, _ = mid_line
    
    dist_dic = {}

    for n, (x, y) in enumerate(zip(xList, yList)):

        dis = ((aX - x)**2 + (aY - y)**2)**.5

        dist_dic[dis] = n

    return dist_dic[min(dist_dic)]

    

def get_better_midlines(test_path):

    import pandas as pd
    from pprint import pprint as pp

    spotData = pd.read_csv(test_path, sep='\t', header=0)

    nucDic = {}

    for nuc, x, y, z in zip(spotData.nuc, spotData.nucx, spotData.nucy, spotData.nucz):
        nucDic[nuc] = (x, y, z)

    pp(nucDic)

    # interNucDis = get_inter_nuclaer_distance(nucDic)

    # get_number_of_neighbours(nucDic, interNucDis)
    # find_edge_cells(nucDic)
    midline, inter_nuc_dis, anteriorPosition  = easy_edges(nucDic)

    return midline, inter_nuc_dis, anteriorPosition


def find_anterior(big_pair, interNucDis, nucPosDic):

    bigPair1, bigPair2 = big_pair

    bigCount1 = 0
    bigCount2 = 0

    x1, y1, _ = nucPosDic[bigPair1]

    for nuc in nucPosDic:
        if nuc == bigPair1:
            continue

        x2, y2, _ = nucPosDic[nuc]

        if  abs(x1-x2) < 1 * interNucDis or abs(y1-y2) < 1 * interNucDis:
            bigCount1 += 1

    x1, y1, _ = nucPosDic[bigPair2]

    for nuc in nucPosDic:
        if nuc == bigPair2:
            continue

        x2, y2, _ = nucPosDic[nuc]

        if  abs(x1-x2) < 1 * interNucDis or abs(y1-y2) < 1 * interNucDis:
            bigCount2 += 1

    if bigCount1 > bigCount2:
        x1, y1, _ = nucPosDic[bigPair2]
        print(f'Anterior is {x1}, {y1}')
        return (x1, y1)

    else:
        x1, y1, _ = nucPosDic[bigPair1]
        print(f'Anterior is {x1}, {y1}')
        return (x1, y1)


def get_10_closest(source_nuc, nucs):
    from numpy import mean
    dis_dic = {}

    for nuc2 in nucs:

        x1, y1 = source_nuc
        x2, y2, _ = nucs[nuc2]

        distance = (x1-x2)**2 + (y1-y2)**2

        dis_dic[distance] = nuc2

    sorted_dis = sorted(list(dis_dic))

    cx, cy = [], []

    for x in sorted_dis[:11]:

        cx.append(nucs[dis_dic[x]][0])
        cy.append(nucs[dis_dic[x]][1])

    return cx, cy


def get_most_distant(nuc_dic2):

    import matplotlib.pyplot as plt
    biggest = 0
    big_pair = ('nothing to', 'see here')

    for nuc1 in nuc_dic2:
        for nuc2 in nuc_dic2:
            if nuc1 == nuc2:
                continue
            x1, y1, _ = nuc_dic2[nuc1]
            x2, y2, _ = nuc_dic2[nuc2]

            distance = (x1-x2)**2 + (y1-y2)**2

            if distance > biggest:
                biggest = distance
                big_pair = (nuc1, nuc2)

    biggest = biggest**.5
    print(biggest, big_pair)

    
    return big_pair
    

def plot_nucs(nucDic, nucs_to_colour, plotcross, midline):

    import matplotlib.pyplot as plt

    xls, yls, cxls, cyls = [],[],[],[]

    for nuc in nucDic:
        xls.append(nucDic[nuc][0])
        yls.append(nucDic[nuc][1])

    for nuc in nucs_to_colour:
        cxls.append(nucDic[nuc][0])
        cyls.append(nucDic[nuc][1])
    plt.figure(figsize=(8,7))
    plt.scatter(xls, yls, c='grey')
    plt.scatter(cxls, cyls, c='yellow', s=0.5)

    midx, midy = plotcross
    plt.scatter(midx, midy, c='orange', marker='x')
    plt.plot(midline[0], midline[1])
    plt.savefig(args.folder_path + '/newMid.png', dpi=300)
    
    # plt.show()
    plt.close()


def get_distance_mid(mid_line, x, y, points_file, mid_start):

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

    pos_along_mid = dis_ls.index(min(dis_ls))

    x_match = x_mid[pos_along_mid]
    y_match = y_mid[pos_along_mid]

    dis_along_mid = ((x_mid[mid_start] - x_match)**2) + ((y_mid[mid_start] - y_match)**2)

    return min(dis_ls)**.5 * trans, trans, dis_along_mid**.5


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

def get_imaris_files(folder_path, folder_list, file_list):
    """
    Searches through the embyro folder (folder_path) and gets all the relevant files to run the data analysis.

    Folder names are in the folder_list varaible and the files will be passed as the file_list variable.

    Data is returned as a dictionary with the top level folders as keys containing the te relevant internal files.
    """

    import os
    from pprint import pprint as pp

    file_dic = {}

    spots_added = False

    for folder in os.listdir(folder_path):
        for x in folder_list:
            if x in folder.lower():
                if 'nuc' not in x:
                    spots_added = True

                file_dic[x] = {'path': folder, 'files': []}
                files_in_folder = os.listdir(os.path.join(folder_path, folder))

                for file in files_in_folder:
                    for y in file_list:
                        if y.lower() in file.lower():
                            file_dic[x]['files'].append(file)

    if not spots_added:
        print('no spots', spots_added)
        for folder in os.listdir(folder_path):
            if 'spot' in folder.lower():
                file_dic['spot'] = {'path': folder, 'files': []}

                files_in_folder = os.listdir(os.path.join(folder_path, folder))

                for file in files_in_folder:
                    
                    if 'spot' in file.lower():
                        file_dic['spot']['files'].append(file)

    # pp(file_dic)
    # exit()
    return file_dic


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


def write_out_spots_dic(o_out, nuc_spots, nuc_dic, spot_dic, spot_type):

    import os

    o_out.write('time\tnuc\tnum_spots\tnucx\tnucy\tnucz\tspot\tspot_type\tspotx\tspoty\tspotz\n')

    for t in nuc_dic:

        for nuc in nuc_dic[t]:
            try:
                num_spots = len(nuc_spots[t][nuc])

                nucx, nucy, nucz = nuc_dic[t][nuc]

                for spot in nuc_spots[t][nuc]:

                    spotx, spoty, spotz = spot_dic[t][spot]

                    out_ls = [str(x) for x in [t, nuc, num_spots, nucx, nucy, nucz, spot, spot_type, spotx, spoty, spotz]]

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
                outfile.write(line.replace('\n', f'\tna\tna\tna\tna\n'))
                continue
            # print(type(spot), spot)
            print(spot)
            for variable in spot_dic[t][spot]:

                for chan in spot_dic[t][spot][variable]:

                    value = spot_dic[t][spot][variable][chan]
                    # print(t, spot, channel, variable, value)
                    outfile.write(line.replace('\n', f'\t{variable}\t{chan}\t{value}\n'))


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


if __name__ == "__main__":

    import argparse
    import os
    from pprint import pprint as pp
    import pandas as pd
    # import multiprocessing as mp

  
    parser = argparse.ArgumentParser()
    parser.add_argument('folder_path', help="folder path where the spots and nuc files should be")
    parser.add_argument('-b', type=float, default=3, help='bin_size in cell widths')
    parser.add_argument('-ti', type=float, default=20, help='time interval')
    parser.add_argument('-sf', type=float, default=0, help='spot distance filter')
    args = parser.parse_args()
    

    # for debugging
    # args.folder_path = test_path


    folder_list = ['nucl', 'exonic', 'intronic']
    file_list = ["_Position.", "_Intensity_", "Ellipsoid_Axis"]

    files_dic = get_imaris_files(args.folder_path, folder_list, file_list)
    pp(files_dic)
    positions_dic = {}
    
    for folder in files_dic:

        folder_pos = [x for x in files_dic[folder]['files'] if "_Position." in x]
        
        if len(folder_pos) != 1:

            exit(f'wrong number of positions files detected: {folder} || {folder_pos}')

        positions_dic[folder] = get_positions(os.path.join(args.folder_path, files_dic[folder]['path'], folder_pos[0]))

    # len1 = len(positions_dic['spots'][1])

    if args.sf != 0:

        if 'spot' in positions_dic:
            positions_dic['spot'] = filter_spots(positions_dic['spot'].copy(), args.sf)

        else:
            # print(list(positions_dic))
            positions_dic['exonic'] = filter_spots(positions_dic['exonic'].copy(), args.sf)
            positions_dic['intronic'] = filter_spots(positions_dic['intronic'].copy(), args.sf)

    # len2 = len(positions_dic['spots'][1])
    # print(len1, len2, len1-len2)
    
    nuc_axis_dic = import_axis([x for x in files_dic['nucl']['files'] if "Ellipsoid_Axis" in x], os.path.join(args.folder_path, files_dic['nucl']['path']))
    print(nuc_axis_dic)
    pp(files_dic)
    
    nuc_with_axis = assign_nuc_axis_dic(nuc_axis_dic, positions_dic['nucl'])

    # pp(nuc_with_axis)
    
    if 'spot' in positions_dic:
        nuc_spots = assign_spots_new(positions_dic['spot'], positions_dic['nucl'], nuc_with_axis)

    else:
        nuc_spots_ex = assign_spots_new(positions_dic['exonic'], positions_dic['nucl'], nuc_with_axis)
        nuc_spots_in = assign_spots_new(positions_dic['intronic'], positions_dic['nucl'], nuc_with_axis)

    # nuc_spots = reassign_spots(positions_dic, nuc_spots, nuc_with_axis)
    # exit()
    
    out_dir = args.folder_path + '/time_data'

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    outfile = out_dir + '/position_data.txt'

    with open(outfile, 'w') as o_out:

        if 'spot' in positions_dic:
            write_out_spots_dic(o_out, nuc_spots, spot_dic=positions_dic['spot'], nuc_dic=positions_dic['nucl'], spot_type='unknown')
        
        else:
            write_out_spots_dic(o_out, nuc_spots_ex, spot_dic=positions_dic['exonic'], nuc_dic=positions_dic['nucl'], spot_type='exonic')
            write_out_spots_dic(o_out, nuc_spots_in, spot_dic=positions_dic['intronic'], nuc_dic=positions_dic['nucl'], spot_type='intronic')

    pos_file = args.folder_path+'/time_data/position_data.txt'

    if 'spot' in positions_dic:

        spot_file_path = os.path.join(args.folder_path, files_dic['spot']['path'])

        spotted_dic = get_spot_info(spot_file_path, [x for x in files_dic['spot']['files']if "_Intensity_" in x])
        pp(spotted_dic)
        assign_intensities(args.folder_path, spotted_dic, pos_file)
    
    intense_path = args.folder_path+ '/time_data/position_data-intense.txt'

    # mid_line = get_mid(intense_path)
    # exit()
    # return_mid_line(mid_line, args.folder_path)


    mid_line, inter_nuc_dis, anteriorPosition = get_better_midlines(intense_path)
    distance_to_midline(intense_path, mid_line, inter_nuc_dis, anteriorPosition)

    
    bin_size = inter_nuc_dis*args.b
    pandas_frame = os.path.split(intense_path)[0]+'/all_data_with_distance.txt'
    bin_ls, bin_numbers = get_bins(pd.read_table(pandas_frame, sep='\t', header=0), bin_size)
    plot_data(pd.read_table(pandas_frame, sep='\t', header=0), bin_ls, bin_numbers, bin_size, args.b, args.folder_path)

    # correct offset

    dataLoad = pd.read_table(pandas_frame, sep='\t', header=0)
    print(dataLoad.tail())
    print(min(dataLoad.dis_along_line))
    # dataLoad.dis_along_line =  [x- min(dataLoad.dis_along_line) for x in list(dataLoad.dis_along_line)]
    dataLoad.dis_along_mid_in_nuclei = dataLoad.dis_along_mid_in_nuclei - min(dataLoad.dis_along_mid_in_nuclei)

    pd.DataFrame.to_csv(dataLoad,pandas_frame, sep='\t', index=None)