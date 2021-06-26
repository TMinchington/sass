def get_mid_soph(pos_frame, out_path):
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import os

    image_file = out_path+'/mid_line_image.png'

    # overexpand position list by spot count
    curvey = 2
    x_ls = []
    y_ls = []

    for x, y, count in zip(pos_frame['Cell Position X'], pos_frame['Cell Position Y'], pos_frame['expressing']):
        for n in range(0, int(count)):
            x_ls.append(x)
            y_ls.append(y)
    
    f = np.poly1d(np.polyfit(x_ls, y_ls, curvey))
    # x = range(0, int(max(x_ls)))
    x_new = np.linspace(min(x_ls)-10, max(x_ls)+10, 300)
    y_new = f(x_new)

    return (x_new, y_new, f)

def get_bins(series, bin_size):

    import numpy as np

    max_s = max(series)
    min_s = min(series)

    return np.linspace(min_s, max_s, len(range(int(min_s-bin_size), int(max_s+bin_size)))/bin_size)

def example_plots(df, out_path, bin_size):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    [print(x) for x in set(df.variable) if 'distance' in x.lower()]
    # exit()
    variablesNeeded = ['Cell Position X', 'Cell Position Y', 'Cell Intensity Sum of Vesicles', 'Cell Intensity Mean of Vesicles', 'Cell Number Of Vesicles', 'Cell Volume', 'disMid', 'distanceInCells']

    plottable_data = df.loc[df.variable.isin(variablesNeeded)]
    plottable_data.replace(to_replace="-nan(ind)", value=np.nan, inplace=True)
    plottable_data.replace(to_replace=" Vesicles hnt", value="Vesicles hnt", inplace=True)
    plottable_data.replace(to_replace=" Vesicles ush", value="Vesicles ush", inplace=True)

    del df
    del plottable_data['CellID']

    plottable_data.value = [float(x) for x in plottable_data.value]

    cell_data = plottable_data.loc[plottable_data.variable.isin(['Cell Position X', 'Cell Position Y', 'Cell Volume', 'disMid', 'distanceInCells'])]
    print(cell_data)
    del cell_data['Channel']
    del cell_data['VesicleType']

    cell_data = cell_data.set_index(['ID', 'variable']).unstack(level=-1)
    cell_data.columns = cell_data.columns.droplevel([0])
    cell_data = cell_data.rename_axis([None], axis=1).reset_index()

    print(cell_data)

    vesicle_data = plottable_data.loc[plottable_data.variable.isin(['Cell Intensity Sum of Vesicles', 'Cell Intensity Mean of Vesicles', 'Cell Number Of Vesicles'])].copy()

    vesicle_data = vesicle_data.set_index(['ID', 'VesicleType', 'Channel', 'variable']).unstack(level=-1)
    vesicle_data.columns = vesicle_data.columns.droplevel([0])
    vesicle_data = vesicle_data.rename_axis([None], axis=1).reset_index()
    
    print(vesicle_data)

    vesicle_data.Channel = [float(x) for x in vesicle_data.Channel]

    plotting_table = pd.merge(cell_data, vesicle_data, on="ID", how="outer", validate="one_to_many")

    # add bins
    # bin_size = 2
    # plotting_table['binnedCells'] = pd.qcut(plotting_table['distanceInCells'], q=((max(plotting_table['distanceInCells'])-min(plotting_table['distanceInCells']))/bin_size))
    myBins = get_bins(plotting_table.distanceInCells, bin_size)

    plotting_table['binnedCells'] = pd.cut(plotting_table['distanceInCells'], bins=myBins, labels=myBins[1:])
    
    print(plotting_table)
    print(set(plotting_table.VesicleType))
    g = sns.FacetGrid(plotting_table.loc[plotting_table.Channel.isin([1.0,2.0])], row="Channel", hue='VesicleType', palette='mako', margin_titles=True)
    g.map(sns.scatterplot, "disMid", "Cell Intensity Mean of Vesicles", alpha=.5)
    g.add_legend()

    plt.savefig(f'{out_dir}/VesicleMeanIntensity.png', transparent=True)
    # print(type(plotting_table['Cell Number Of Vesicles'][0]))
    # print(plotting_table.loc[pd.notna(plotting_table['Cell Number Of Vesicles'])])
    plt.close()
    sns.relplot(x="disMid", y="Cell Number Of Vesicles", hue="VesicleType",
            sizes=(40, 400), alpha=.5, palette="mako",
            height=6, data=plotting_table)
    plt.savefig(f'{out_dir}/NumberOfVeiscles.png', transparent=True)
    plt.close()


    sns.lineplot(x="binnedCells", y="Cell Number Of Vesicles", hue="VesicleType", err_style="band",
             ci='sd',
             estimator="mean", palette="rocket",
             data=plotting_table)

    # plt.show()
    
    plt.savefig(f'{out_dir}/bins.png', transparent=True)
    plt.close()
    # plt.show()
    # plt.savefig("D:\\Dropbox (The University of Manchester)\\0-PhD\\9-Other-Lab_Members\\1-CAROL\\Live imaging blahhhhhh\\Sophie's dataset\\Embryo1_Statistics\\time_data\\bins.png")
   
    pd.DataFrame.to_csv(plotting_table, f"{out_dir}/easy_format_data.tsv", sep='\t', index=None)

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


def distance_to_midline(pos_frame, mid_line):

    import os
    import seaborn as sns
    import matplotlib.pyplot as plt

    # outfile = open(os.path.split(out_dir)[0]+'/all_data_with_distance.txt', 'w')
    points_file = open(out_dir + '/mid_points2.txt', 'w')
    dis_ls = []
    for x, y in zip(pos_frame['Cell Position X'], pos_frame['Cell Position Y']):

        dis_mid, above_line = get_distance_mid(mid_line, x, y, points_file)
        dis_ls.append(dis_mid)
    
    pos_frame['dis_mid'] = dis_ls
    pos_frame['abs_dis_mid'] = abs(pos_frame['dis_mid'])
    # print(pos_frame)
    
    sns.scatterplot(data=pos_frame, x="Cell Position X", y="Cell Position Y", s=80, hue="dis_mid", legend=None, palette="rocket")
    sns.lineplot(x=mid_line[0], y=mid_line[1])
    plt.savefig(f"{out_dir}/dis_mid.png")
    plt.close()

    del pos_frame['abs_dis_mid']


def get_inter_cell_distance(cell_dic):

    """
    cell_dic: dictionary of nuclear positions

    loops throughall nuclei and finds the ddistance to it's closest neighbour, the mean minimal distance is then used as an estiamte of cell size.

    returns mean minimum nuclear distances

    """

    from numpy import mean

    min_dis_ls = []
    for cell in cell_dic:
        dis_ls = []
        # print(cell_dic[cell])
        for cell2 in cell_dic:
            if cell == cell2:
                # print(cell, cell2)
                continue
            
            cellx1, celly1, cellz1 = cell_dic[cell]
            cellx2, celly2, cellz2 = cell_dic[cell2]

            distance_calculated = (cellx1-cellx2)**2 + (celly1-celly2)**2 + (cellz1-cellz2)**2
            # print(distance_calculated)
            dis_ls.append(distance_calculated)
        # print(dis_ls)
        min_dis_ls.append(min(dis_ls)**.5)

    return mean(min_dis_ls)


def gather_all_cell_data(folder_path):

    import os
    import pandas as pd

    folder_ls = list(set([ x for x in os.listdir(folder_path) if 'cell' in x.lower() and '.csv' in x] + [x for x in os.listdir(folder_path) if 'vesicle' in x.lower() and '.csv' in x]))

    # print(folder_ls)
    [print(x) for x in folder_ls if 'Vesicle' in x]
    vesicle_position_data = 'bob'
    cell_positions = 'empty_str'
    all_data = 'str'
    for filex in folder_ls:
        filex_path = os.path.join(folder_path, filex)
        
        data_temp = pd.read_csv(filex_path, header=2)
        

        if 'x' in data_temp.columns[0].lower() and 'y' in data_temp.columns[1].lower() and 'z' in data_temp.columns[2].lower():
            ids_needed = ['ID', 'Time', 'Channel', 'CellID', 'VesicleType']
            missing_ids = [x for x in ids_needed if x not in data_temp.columns]
            data_temp2 = data_temp.melt(value_vars=list(data_temp.columns[:3]), id_vars=[x for x in data_temp.columns if x in ids_needed])
            # print(data_temp.columns)
            # print(data_temp.head())
            # print(data_temp2.head())

            for missID in missing_ids:
                data_temp2[missID] = "NAN"

            if data_temp.columns[0] == 'Cell Position X':
                cell_positions = data_temp

            if data_temp2.variable[0].startswith('Vesicles'):
                # print(data_temp.columns)
                if type(vesicle_position_data) == str:
                    vesicle_position_data = data_temp2
                    # print('>', list(set(list(vesicle_position_data.variable))))
                else:

                    vesicle_position_data = pd.concat([vesicle_position_data, data_temp2], axis=0)
                    # print(list(set(list(vesicle_position_data.variable))))

            # reshape variables 

            # print(data_temp2.head())

            data_temp2 = data_temp2[['ID', 'variable', 'value', 'Channel', 'CellID', 'VesicleType']]

            if type(all_data) == str:
                all_data = data_temp2
                
                # print('>', list(set(list(vesicle_position_data.variable))))
            else:

                all_data = pd.concat([all_data, data_temp2], axis=0)
                continue

        data_temp['variable'] = data_temp.columns[0]
        headList = list(data_temp.columns)
        headList[0] = 'value'
        data_temp.columns = headList
        # print(data_temp.head())
        for missingHead in ['CellID', 'Channel', 'VesicleType']:
            if missingHead not in headList:
                data_temp[missingHead] = 'NAN'

        data_temp = data_temp[['ID', 'variable', 'value', 'Channel', 'CellID', 'VesicleType']]

        if type(all_data) == str:
            all_data = data_temp
            
            # print('>', list(set(list(vesicle_position_data.variable))))
        else:

            all_data = pd.concat([all_data, data_temp], axis=0)
            
            # print(list(set(list(vesicle_position_data.variable))))

    positive_cells = list(set(list(vesicle_position_data.CellID)))
    # print(all_data)
    return positive_cells, all_data, cell_positions


def show_positive_cells(count_dic, cell_positions, out_dir):
    import seaborn as sns
    # import 
    sns.set()
    # print(count_dic)
    expressing = []

    for cell in cell_positions['ID']:
        try:
            expressing.append(count_dic[cell])
            # print(cell, count_dic[cell])
        except KeyError:
            expressing.append(0)


    cell_positions['expressing'] = expressing
    
    # print(cell_positions)
    # sns_plot = sns.scatterplot(data=cell_positions, x="Cell Position X", y="Cell Position Y", hue="expressing")
    sns_plot = sns.scatterplot(data=cell_positions, x="Cell Position X", y="Cell Position Y", size=100, hue="expressing", legend=None, palette="mako")
    fig = sns_plot.get_figure()
    fig.savefig(f"{out_dir}/spots_total.png")

def get_number_of_vesicles(all_data):
    # print(all_data.head())
    vesicle_number_vars = list(set([x for x in all_data.variable if 'Cell Number Of Vesicles' in x]))
    # print(vesicle_number_vars)
    
    count_frame = all_data.loc[all_data.variable.isin(vesicle_number_vars)]

    # print(count_frame)
    # exit()
    count_dic = {}

    for cell, count in zip(count_frame.ID, count_frame.value):
        
        try:
            try:
                # print(type(cell))
                count_dic[cell] += float(count)
            except ValueError:
                continue

        except KeyError:
            count_dic[cell] = float(count)
    import pprint
    # pprint.pprint(count_dic)
    # exit()
    return count_dic

def make_compatible_cell_dic(df):
    """
    Take a dataframe containing x,y,z coordinates and an ID
    Returns a dictionary in the format dic[ID] = (x, y, z)
    """

    cell_dic = {}

    for cell_id, x, y, z in zip(df['ID'], df['Cell Position X'], df['Cell Position Y'], df['Cell Position Z']):
        cell_dic[cell_id] = (x, y, z)

    return cell_dic

if __name__ == "__main__":

    import argparse
    import os
    from pprint import pprint as pp
    import pandas as pd
  
    parser = argparse.ArgumentParser()
    parser.add_argument('folder_path')
    parser.add_argument('-b', type=float, default=3, help='bin_size in cell widths')
    parser.add_argument('-ti', type=float, default=20, help='time interval')
    # parser.add_argument('-sf', type=float, default=.5, help='spot distance filter')
    args = parser.parse_args()
    # args.folder_path = "D:\\Dropbox (The University of Manchester)\\0-PhD\\9-Other-Lab_Members\\1-CAROL\\Live imaging blahhhhhh\\Sophie's dataset\\Embryo1_Statistics"
    positive_cells, all_data, cell_positions = gather_all_cell_data(args.folder_path)
    
    count_dic = get_number_of_vesicles(all_data)
    # for cell in count_dic:
        # print(cell, count_dic[cell])

    # print(count_dic.keys())
    # print(8 in count_dic.keys())
        
    # exit()
    out_dir = args.folder_path + '/time_data'

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    show_positive_cells(count_dic, cell_positions, out_dir)
    
    pos_file = args.folder_path+'/time_data/position_data.txt'

    mid_line = get_mid_soph(cell_positions, out_dir)
    
    return_mid_line(mid_line, out_dir)

    distance_to_midline(cell_positions, mid_line)
    cell_dic = make_compatible_cell_dic(cell_positions)
    inter_nuc_dis = get_inter_cell_distance(cell_dic)
    bin_size = inter_nuc_dis*args.b
    cell_positions['distanceInCells'] = cell_positions['dis_mid']/inter_nuc_dis

    cell_out = cell_positions[['Cell Position X', 'Cell Position Y', 'Cell Position Z', 'Time', 'ID', 'expressing', 'dis_mid', 'distanceInCells']].copy()
    cell_out.columns = ['CellX', 'CellY', 'CellZ', 'Time', 'ID', 'totalSpots', 'disMid', 'distanceInCells']
    cell_out['absDisMid'] = abs(cell_positions['dis_mid'])

    pd.DataFrame.to_csv(cell_out, f"{out_dir}/cell_positions.tsv", sep='\t', index=None)
    cell_for_all = cell_out[['Time', 'ID', 'totalSpots', 'disMid', 'distanceInCells']]
    cell_for_all = cell_for_all.melt(value_vars=['totalSpots', 'disMid', 'distanceInCells'], id_vars=['Time', 'ID'])

    ids_needed = ['ID', 'variable', 'value', 'Channel', 'CellID', 'VesicleType']


    for head in ids_needed:
        if head not in cell_for_all.columns:
            cell_for_all[head] = 'NAN'

    cell_for_all = cell_for_all[['ID', 'variable', 'value', 'Channel', 'CellID', 'VesicleType']]
    all_data = pd.concat([all_data, cell_for_all], axis=0)
    del cell_for_all
    del cell_positions
    # print(all_data.columns)
    # [print(x) for x in list(set(list(all_data.variable))) if 'Position' in x]
    pd.DataFrame.to_csv(all_data, f"{out_dir}/all_data.tsv", sep='\t', index=None)
    
    example_plots(all_data, out_dir, args.b)
    # pandas_frame = os.path.split(intense_path)[0]+'/all_data_with_distance.txt'
    # bin_ls, bin_numbers = get_bins(pd.read_table(pandas_frame, sep='\t', header=0), bin_size)
    # plot_data(pd.read_table(pandas_frame, sep='\t', header=0), bin_ls, bin_numbers, bin_size, args.b, args.folder_path)

    del all_data

    exit()