# SOLSTISE Simulation Code Front End
from Sim_File import sim
from Event_Builder import BuildEvts
from PipeMaker import makepipe
import glob
import os
import re
import math

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
print("                     Particle Shadowing Simulation Code")
input("\n\n\nTo continue, press ENTER")

list_files = glob.glob('*evts*.txt')

if len(list_files) == 0:
    print("It appears that no event files exist, so we'll create one now!")
    filein = BuildEvts()
else:
    latest_file = max(list_files, key=os.path.getctime)
    print("\nThe most recently created input file is: " + latest_file)
    yn = input("\nWould you like to use this file? [Y/N] ")

    if yn == "N" or yn == "n":
        newevtfileyn = input("Would you like to generate a new event simulation file? [Y/N] ")

        if newevtfileyn == "N" or newevtfileyn == "n":
            print("\n")
            for i in range(len(list_files)):
                print(str(i+1) + ") " + list_files[i])

            filenum = 1000000
            while filenum > len(list_files):
                filenum = int(input("\nChoose a number from the list, or enter 0 to manually type the file name: "))

                if len(list_files) >= filenum > 0:
                    filein = list_files[filenum-1]
                elif filenum > len(list_files):
                    print("ERROR: Number entered is greater than the number of simulation files...")
                else:
                    list_file2 = []
                    numuscore = 0
                    while len(list_file2) == 0 or numuscore != 4:
                        filein = input("\nEnter the name of an existing input file (.dat or .txt): ")
                        list_file2 = glob.glob(filein)
                        numuscore = filein.count("_")
                        if len(list_file2) == 0 or numuscore != 4:
                            print("\nERROR: Incorrect file name syntax or the file does not exist...")
        else:
            filein = BuildEvts()
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

print("The beam energy is: " + str(ebeam) + " MeV")

pipeyn = input("\nWould you like to use a custom shaped pipe for the gas return? [Y/N] ")

if pipeyn == "N" or pipeyn == "n":

    hors = 0

    while hors != 1 and hors != 2:

        hors = int(input("\nWill the reaction be occurring in HELIOS (Enter 1) or SOLARIS (Enter 2)? "))

        if hors != 1 and hors != 2:
            print("\nERROR: Unknown entry...")

    if hors == 1:
        r1 = 0.92 / 2
    else:
        r1 = 0.9 / 2

    piper = 0.148082/2
    pipecenter = 0.148082/2-r1
    phi1 = (180 + 90)*math.pi/180 - math.asin(piper/pipecenter)
    phi2 = 2*math.pi - (phi1 - math.pi)

    sim(r1, piper, pipecenter, math.pi, 2*math.pi, ebeam, filein, reac)

else:
    list_pipe_files = glob.glob('PipeOut*.txt')

    if len(list_pipe_files) == 0:
        print("\nNo pipe setup file exists in your directory, so we'll make one now.")

        pipedeffile = makepipe()
    else:
        latest_pipe_file = max(list_pipe_files, key=os.path.getctime)
        print("The most recently created pipe setup file is: " + latest_pipe_file)
        pfileyn = input("Would you like to use this file? [Y/N] ")
        if pfileyn == "N" or pfileyn == "n":
            pfileyn2 = input("Would you like to generate a new custom pipe input file? [Y/N] ")
            if pfileyn2 == "N" or pfileyn2 == "n":
                print("\n")
                for i in range(len(list_pipe_files)):
                    print(str(i + 1) + ") " + list_pipe_files[i])

                filenum = 1000000
                while filenum > len(list_pipe_files):
                    filenum = int(input("\nChoose a number from the list, or enter 0 to manually type the file name: "))

                    if len(list_pipe_files) >= filenum > 0:
                        pipedeffile = list_pipe_files[filenum - 1]
                    elif filenum > len(list_pipe_files):
                        print("ERROR: Number entered is greater than the number of simulation files...")
                    else:
                        list_file_check = []
                        while len(list_file_check) == 0:
                            pipedeffile = input("Enter the name of the pipe input file: ")
                            list_file_check = glob.glob(pipedeffile)
                            if len(list_file_check) == 0:
                                print("ERROR: The file does not exist...")
            else:
                pipedeffile = makepipe()
        else:
            pipedeffile = latest_pipe_file

    file = open(pipedeffile, "r")

    lines = file.readlines()

    sim(float(lines[0]), float(lines[1]), float(lines[2]), float(lines[3]), float(lines[4]), ebeam, filein, reac)
