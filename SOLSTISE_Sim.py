# SOLSTISE Simulation Code Front End
from Sim_File_pd import sim_pd
from Event_Builder import BuildEvts
from PipeMaker import makepipe
from ConeMaker import makecone
from NozzleMaker import makenozz
import numpy as np
import glob
import os
import re
import math
import fnmatch

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
                        numuscopre = filein.count("_")
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

# need to load the file into a numpy array to do genfromtxt to determine whether or not the reaction is going to be
# measured in normal or inverse kinematics, then ask which way the return pipe is facing. If the pipe is facing the
# opposite direction, we'll just skip over the pipe definition and make the radius of the pipe tiny so it doesn't
# interfere with the reaction products...

datas = np.genfromtxt(filein)
if datas[:, 0].mean() > 90:
    invkin = 1
else:
    invkin = 0

if not fnmatch.fnmatch(filein, '*eloss_s*'):
    try:
        coneopt = int(input("\nWould you like to use the default SOLSTISE cone (0) or a custom cone (1)?: "))
    except ValueError:
        coneopt = 0
else:
    coneopt = 0

if coneopt == 0:
    conetxt = "SOLSTISE_cone_3_2-6in.txt"

else:
    list_cones = glob.glob('*cone*.txt')
    latest_cone = max(list_cones, key=os.path.getctime)
    print("\nThe most recently created cone input file is: " + latest_cone)
    yn = input("\nWould you like to use this file? [Y/N] ")

    if yn == 'n' or yn == 'N':
        newconeyn = input("Would you like to create a new cone input file? [Y/N] ")

        if newconeyn == 'n' or newconeyn == 'N':
            print("\n")
            for i in range(len(list_cones)):
                print(str(i + 1) + ") " + list_cones[i])

            filenum = 1000000
            while filenum > len(list_cones):
                filenum = int(input("\nChoose a number from the list, or enter 0 to manually type the file name: "))

                if len(list_cones) >= filenum > 0:
                    conetxt = list_cones[filenum - 1]
                elif filenum > len(list_cones):
                    print("ERROR: Number entered is greater than the number of cone input files...")
                else:
                    list_file_check = []
                    while len(list_file_check) == 0:
                        conetxt = input("Enter the name of the cone input file: ")
                        list_file_check = glob.glob(conetxt)
                        if len(list_file_check) == 0:
                            print("ERROR: The file does not exist...")
        else:
            conetxt = makecone()
    else:
        conetxt = latest_cone

if not fnmatch.fnmatch(filein, '*eloss_s*'):
    try:
        nozzopt = int(input("\nWould you like to use the default SOLSTISE nozzle (0) or a custom-shaped nozzle (1)?: "))
    except ValueError:
        nozzopt = 0
else:
    nozzopt = 0

# We'll put the nozzle questions here, since there are only ~3 and it isn't worth a new function or code.
if nozzopt == 0:
    nozztxt = 'SOLSTISE_nozz_default.txt'
else:
    list_nozz = glob.glob('*nozz*.txt')
    latest_nozz = max(list_nozz, key=os.path.getctime)
    print("\nThe most recently created nozzle input file is: " + latest_nozz)
    yn = input("\nWould you like to use this file? [Y/N] ")

    if yn == 'n' or yn == 'N':
        newnozzyn = input("Would you like to create a new nozzle input file? [Y/N] ")

        if newnozzyn == 'n' or newnozzyn == 'N':
            print("\n")
            for i in range(len(list_nozz)):
                print(str(i + 1) + ") " + list_nozz[i])

            filenum = 1000000
            while filenum > len(list_nozz):
                filenum = int(input("\nChoose a number from the list, or enter 0 to manually type the file name: "))

                if len(list_nozz) >= filenum > 0:
                    nozztxt = list_nozz[filenum - 1]
                elif filenum > len(list_nozz):
                    print("ERROR: Number entered is greater than the number of nozzle input files...")
                else:
                    list_file_check = []
                    while len(list_file_check) == 0:
                        nozztxt = input("Enter the name of the nozzle input file: ")
                        list_file_check = glob.glob(nozztxt)
                        if len(list_file_check) == 0:
                            print("ERROR: The file does not exist...")
        else:
            nozztxt = makenozz()
    else:
        nozztxt = latest_nozz

# Remove this question if the user is using a solid target:
if not fnmatch.fnmatch(filein, '*eloss_s*'):
    try:
        pipefb = int(input("\nIs the pipe for the gas return in the downstream (0) or upstream (1) half of the "
                           "magnet?: "))
    except ValueError:
        pipefb = 0
else:
    pipefb = 0

if pipefb != invkin:
    pipeyn = 'n'
else:
    pipeyn = input("\nWould you like to use a custom shaped pipe for the gas return? [Y/N] ")

if pipeyn == "N" or pipeyn == "n":

    hors = 0

    while hors != 1 and hors != 2:
        try:
            hors = int(input("\nWill the reaction be occurring in HELIOS (Enter 1) or SOLARIS (Enter 2)? "))
        except ValueError:
            hors = 1

        if hors != 1 and hors != 2:
            print("\nERROR: Unknown entry...")

    if hors == 1:
        r1 = 0.92 / 2
    else:
        r1 = 0.9 / 2

    if pipefb != invkin:
        piper = 0.001
        pipecenter = -r1
        phi1 = 3/2*math.pi
        phi2 = 3/2*math.pi
    else:
        piper = 0.152908/2
        pipecenter = 0.152908/2-r1
        phi1 = (180 + 90)*math.pi/180 - math.asin(piper/pipecenter)
        phi2 = 2*math.pi - (phi1 - math.pi)

    sim_pd(r1, piper, pipecenter, math.pi, 2*math.pi, ebeam, filein, reac, conetxt, nozztxt)

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

    sim_pd(float(lines[0]), float(lines[1]), float(lines[2]), float(lines[3]), float(lines[4]),
           ebeam, filein, reac, conetxt, nozztxt)
