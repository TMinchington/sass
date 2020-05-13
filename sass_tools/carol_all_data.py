"""

Runs all scripts in order to run Carols live analysis

"""

if __name__ == "__main__":

    import argparse
    import os
    import subprocess
    from spot_the_dffierence import get_experiments_folder, get_files, get_spot_positions, \
        get_nuc_positions, check_nuc_spots
    from pprint import pprint as pp

    from so_intense import assign_intensities, get_spot_info

    from elipse import import_axis, make_positions_dic, assign_spots_new, write_out_spots_dic

    from take_it_back_again import plot_spots, get_mid, remove_background, return_mid_line

    parser = argparse.ArgumentParser()

    parser.add_argument("exp_dir")

    parser.add_argument("--i", default=20, help="time interval in seconds between points. Default: 20 seconds",
                        type=float)

    args = parser.parse_args()

    """---------------------------- Elipse ------------------------------------------------------------------"""

    print('Running Elipse')
    #
    folders_dic = get_experiments_folder(args.exp_dir)

    pp(folders_dic)

    spot_files = get_files(folders_dic['spot'])
    nuc_files = get_files(folders_dic['nuc'])

    spot_pos = get_spot_positions(spot_files)
    nuc_pos = get_nuc_positions(nuc_files)
    pp(spot_files)

    nuc_axis_dic = import_axis(folders_dic['nuc'], nuc_files)

    pos_dic = make_positions_dic(nuc_axis_dic, nuc_pos)

    nuc_spots = assign_spots_new(spot_pos, nuc_pos, pos_dic)


    check_nuc_spots(nuc_spots, args.exp_dir+'/assignemnets.png')

    write_out_spots_dic(nuc_spots, exp_dir=args.exp_dir, spot_dic=spot_pos, nuc_dic=nuc_pos)


    """---------------------------- so intense ----------------------------------------------------------------------"""

    folder_dic = get_experiments_folder(args.exp_dir)
    spot_files = get_files(folder_dic['spot'])

    pos_file = args.exp_dir + '/time_data/position_data.txt'

    if not os.path.isfile(pos_file):
        exit('Position file not found')

    spotted_dic = get_spot_info(args.exp_dir, spot_files)

    assign_intensities(args.exp_dir, spotted_dic, pos_file)


    """---------------------------- take it back again --------------------------------------------------------------"""

    print('Running take it back again')

    intense_path = args.exp_dir + '/time_data/position_data-intense.txt'

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
    return_mid_line(mid_line, args.exp_dir)
    remove_background(intense_path, back_eq_lines, args.i, nuc_dic, mid_line)

    subprocess.run(['python', './sass_tools/hotel_california.py', f'{args.exp_dir}'])
    subprocess.run(['python', './sass_tools/my_first_time.py', f'{args.exp_dir}'])
    subprocess.run(['python', './sass_tools/time_project.py', f'{args.exp_dir}'])
    






