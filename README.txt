---------------------------------------------------------

Welcome to SASS (Simple Assignemnt of Spots to Surfaces)

---------------------------------------------------------


This script executes a set of subscripts used to assign spots detected by imaris to their closest surface in Drosophila embryos. 
It will also return a midline for an expression domain and the distances of cells from that expression domain. 

This is currently only tested on expression domains running along the AP axis as used in: 

 Modulation of promoter occupancy dictates the transcriptional response to graded BMP signalling levels in the Drosophila embryo
 Caroline Hoppe, Jonathan Bowles,  Thomas G. Minchington, Catherine Sutcliffe,  Priyanka Upadhyai,  Magnus Rattray,  Hilary L. Ashe
doi: https://doi.org/10.1101/837179

Use of this script and all components requires aknowldgement to the author Thomas G Minchington and citation of the above paper.

For any support please contact thomas.minchington@manchester.ac.uk or use github.

--------------
Installation
--------------

1. Download and install python3. This can be from either python of Anaconda.
2. create a virtual environment:
    a) anaconda: conda create --name sass_env
    b) python3: python -m venv sass_env

3. enter virtual environment:
    a) anaconda: activate sass_env/bin/activate
    b) python3: source sass_env/bin/activate

4. install requirements: pip install -r requirements.txt (probably not required on anaconda)

Then you should be good to go.



This repository can be downloaded and as long as all the scripts remain in their relative locations it should work fine. 

Usage is outlined below. The script runs from the command line but, is very easy to use. On MAC/Unix based machines you should use terminal and on windows command prompt.

You can navigate to the directory containing the folder easily by typing cd into the terminal window and then dragging the folder containing the sass.py script into terminal/cmd window.

If you press enter on the keyboard that should take you to the correct folder.

You can then run the script by typing (do not include the square brackets): python sass.py [drag in the folder containing the surface folder, background (spots) folder and spot data (folder)] -i [the time interval between frames]

usage: sass.py [-h] [--i I] fly_folder

positional arguments:
  fly_folder  Folder containing folders of flies, each sub folder should
              contain background folder, spot folder and a cell/nucleus folder
              

optional arguments:
  -h, --help  show this help message and exit
  --i I       time interval in seconds between points. Default: 20 seconds
