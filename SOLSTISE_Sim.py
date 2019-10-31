# SOLSTISE Simulation Code Front End
from Sim_File import sim
from Event_Builder import BuildEvts
from PipeMaker import makepipe
import glob
import os
import re

print("")
print("########   ########  ##       ########  ########  ########  ########  ######## ")
print("##         ##    ##  ##       ##           ##        ##     ##        ##          ")
print("##         ##    ##  ##       ##           ##        ##     ##        ##          ")
print("########   ##    ##  ##       ########     ##        ##     ########  ########    ")
print("      ##   ##    ##  ##             ##     ##        ##           ##  ##          ")
print("      ##   ##    ##  ##             ##     ##        ##           ##  ##          ")
print("########   ########  #######  ########     ##     ########  ########  ########   ")
print("")
print("            Solenoid & Supersonic Target In Structure Experiments")
print("")
print("                          Particle Simulation Code")
input("\n\n\nTo continue, press ENTER")

list_files = glob.glob('*evts.txt')

if len(list_files) == 0:
    print("It appears that no event files exist, so we'll create one now!")
    filein = BuildEvts()
else:
    latest_file = max(list_files, key=os.path.getctime)
    print("\nThe most recently created input file is: " + latest_file)
    yn = input("\nWould you like to use this file? [Y/N] ")

    if yn == "N" or yn == "n":
        list_file2 = []
        numuscore = 0
        while len(list_file2) == 0 or numuscore != 4:
            filein = input("\nEnter the name of the new input file (.dat or .txt): ")
            list_file2 = glob.glob(filein)
            numuscore = filein.count("_")
            if len(list_file2) == 0 or numuscore != 4:
                print("\nERROR: Incorrect file syntax or the file does not exist...")
    else:
        filein = latest_file

print("\nThe file to be used is: " + filein)

# Get the reaction out of the file name
# the locations of the underscores in the file name are found here and put into an array
uslocs = [m.start() for m in re.finditer('_', filein)]

reac = filein[:uslocs[0]] + "(" + filein[(uslocs[0]+1):uslocs[1]] + "," + filein[(uslocs[1]+1):uslocs[2]] + ")"

print("The reaction is: " + reac)

# Get the beam energy from the file name:

ebeam = int(filein[(uslocs[2]+1):uslocs[3]])

print("The beam energy is: " + str(ebeam))

pipeyn = input("\nWould you like to use an ISO 160 pipe for the gas return? [Y/N] ")

if pipeyn == "N" or pipeyn == "n":

    list_pipe_files = glob.glob('PipeOut*.txt')

    if len(list_pipe_files) == 0:
        print("\nNo pipe setup file exists in your directory, so we'll make one now.")

        pipedeffile = makepipe()
    else:
        latest_pipe_file = max(list_pipe_files, key=os.path.getctime)
        print("The most recently created pipe setup file is: " + latest_pipe_file)
        pfileyn = input("Would you like to use this file? [Y/N] ")
        if pfileyn == "N" or pfileyn == "n":
            list_file_check = []
            while len(list_file_check) == 0:
                pipedeffile = input("Enter the name of the pipe input file: ")
                list_file_check = glob.glob(pipedeffile)
                if len(list_file_check) == 0:
                    print("ERROR: The file does not exist...")
        else:
            pipedeffile = latest_pipe_file

    file = open(pipedeffile, "r")

    lines = file.readlines()

    sim(float(lines[0]), float(lines[1]), float(lines[2]), float(lines[3]), float(lines[4]), ebeam, filein, reac)

else:
    sim(1, 1, 1, 1, 1, ebeam, filein, reac)
