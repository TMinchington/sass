"""

This is the master controller for SASS (Spot Assignment)
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

    print("Spot Assignement to Surface Script\n_________________________________________\n\n")

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

