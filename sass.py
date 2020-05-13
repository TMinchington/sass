"""

This is the master controller for SASS (Simple Assignment of Spots to Surfaces)
This script is responsible for handelling all other modules each will be documented seperately
Use of this script and all components requires aknowldgement to the author Thomas G Minchington.

For any support please contact thomas.minchington@manchester.ac.uk

usage: sass.py [-h] [--i I] fly_folder

positional arguments:
  fly_folder  Folder containing folders of flies, each sub folder should
              contain background folder, spot folder and a cell/nucleus folder
              

optional arguments:
  -h, --help  show this help message and exit
  --i I       time interval in seconds between points. Default: 20 seconds

"""

if __name__ == "__main__":

    import subprocess
    import multiprocessing as mp
    import argparse
    import os

    print("---------------------------------------------------------\n\nWelcome to SASS (Simple Assignemnt of Spots to Surfaces)\n\n Modulation of promoter occupancy dictates the transcriptional response to graded BMP signalling levels in the Drosophila embry. \nCaroline Hoppe, Jonathan Bowles,  Thomas G. Minchington, Catherine Sutcliffe,  Priyanka Upadhyai,  Magnus Rattray,  Hilary L. Ashe \ndoi: https://doi.org/10.1101/837179 \n\nUse of this script and all components requires aknowldgement to the author Thomas G Minchington and citation of the above paper.\n---------------------------------------------------------\n\n")

    parser = argparse.ArgumentParser()
    parser.add_argument('fly_folder', help='Folder containing folders of flies, each sub folder should contain background folder, spot folder and a cell/nucleous folder')
    parser.add_argument("--i", default=20, help="time interval in seconds between points. Default: 20 seconds",
                        type=float)
    args = parser.parse_args()
    
    file_ls = [x for x in os.listdir(args.fly_folder) if '.' not in x]

    print('Running files:')

    for file in file_ls:

        print(file)

        subprocess.run(['python', './sass_tools/carol_all_data.py', args.fly_folder+'/'+file, '--i', f'{args.i}'])

