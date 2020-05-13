"""

Takes the backed_down file and just keeps cells that are in frame for the whole movie from the onset of expression.

Also makes a plot for the expression domain.

"""

def get_time_points(backed_down):

    t_dic = {}

    firsT = 1_000_000
    lasT = 0

    with open(backed_down) as obd:

        first_line = True

        for line in obd:

            if first_line:

                head_line = line
                first_line = False

                continue

            (time, nuc, num_spots, nucx, nucy,
             nucz, spot, spotx, spoty, spotz, variable,
             channel, value, background, corrected_val, dis_mid, above_line) = line.strip().split('\t')

            time = float(time)

            try:

                t_dic[nuc][time] = 0

            except KeyError:

                t_dic[nuc] = {time:0}

            if float(corrected_val) != 0 and time < firsT:

                firsT = time

            if time > lasT:

                lasT = time

    with open(backed_down) as obd:

        hotel_cal = open(backed_down.replace('.txt', '-hc.txt'), 'w')
        hotel_cal.write(head_line)
        first_line = True

        for line in obd:

            if first_line:

                head_line = line
                first_line = False

                continue

            (time, nuc, num_spots, nucx, nucy,
             nucz, spot, spotx, spoty, spotz, variable,
             channel, value, background, corrected_val, dis_mid, above_line) = line.strip().split('\t')

            if firsT in t_dic[nuc] and lasT in t_dic[nuc]:

                hotel_cal.write(line)

            else:

                continue

    return t_dic



def import_mid(mid_line):

    midx, midy = [], []

    with open(mid_line) as o_mid:
        first_line = True

        for line in o_mid:

            if first_line:

                first_line = False

                continue

            x, y = line.strip().split('\t')

            midx.append(float(x))
            midy.append(float(y))
    print(midx)

    return midx, midy



def get_expression_domain(t_dic, position_file, outpath, mid):

    import matplotlib.pyplot as plt

    o_out = open(outpath, 'w')
    o_out.write('t\tx\ty\tdomain\n')

    # x_ls, y_ls, domain = [], [], []
    plot_dic = {}
    with open(position_file) as opo:

        in_data = False

        for line in opo:
            print('>', line)
            if not in_data and line.startswith('Position X'):
                print(line)
                in_data = True

                continue

            elif not in_data:

                continue

            (x, y, z, unit, cat, collection, time, nuc, ID, un) = line.strip().split(',')
            time = float(time)

            x = float(x)
            y=float(y)

            try:

                plot_dic[time]['x'].append(x)

            except KeyError:

                plot_dic[time] = {'x': [x], 'y': [] , 'domain': []}

            plot_dic[time]['y'].append(y)

            if nuc in t_dic:

                o_out.write(f'{time}\t{x}\t{y}\tin\n')

                plot_dic[time]['domain'].append('green')

            else:

                o_out.write(f'{time}\t{x}\t{y}\tout\n')
                plot_dic[time]['domain'].append('blue')


    t_plot = max(plot_dic)

    plt.scatter(x=plot_dic[t_plot]['x'], y=plot_dic[t_plot]['y'], c=plot_dic[t_plot]['domain'])
    plt.scatter(x=mid[0], y=mid[1])
    plt.savefig(args.exp_dir+'/time_data/mid_line_expression.png')
    plt.close()

if __name__ == "__main__":

    import argparse
    import os
    from spot_the_dffierence import get_experiments_folder
    parser = argparse.ArgumentParser()

    parser.add_argument('exp_dir')

    args = parser.parse_args()
    print('Running')
    t_nuc = get_time_points(args.exp_dir+'/time_data/backed_down.txt')

    folders_dic = get_experiments_folder(args.exp_dir)
    print(folders_dic['nuc'])

    cell_file = [x for x in os.listdir(folders_dic['nuc']) if '_position.' in x.lower() and 'Track' not in x][0]

    get_expression_domain(t_nuc, folders_dic['nuc']+'/'+cell_file, args.exp_dir+'/time_data/expression_domain.txt', import_mid(args.exp_dir+'/mid_line.txt'))