# SOLSTISE Simulation Code Front End
from Sim_File import sim

print("")
print("########   ########  ##       ########  ########  ########  ########  ######## ")
print("##         ##    ##  ##       ##           ##        ##     ##        ##          ")
print("##         ##    ##  ##       ##           ##        ##     ##        ##          ")
print("########   ##    ##  ##       ########     ##        ##     ########  ########    ")
print("      ##   ##    ##  ##             ##     ##        ##           ##  ##          ")
print("      ##   ##    ##  ##             ##     ##        ##           ##  ##          ")
print("########   ########  #######  ########     ##     ########  ########  ########   ")
print("            Solenoid & Supersonic Target In Structure Experiments")
print("")
print("                          Particle Simulation Code")
input("\n\n\nTo continue, press ENTER")
print("\nThe default input file is: eventsout.txt")
yn = input("\nWould you like to use this file? (Y/N) ")

filein = "eventsout.txt"

# if the user wants to enter a new file, get the file name here.
# it must end in .dat or .txt to be recognized.

ebeam = 168

if yn == "N" or yn == "n":
    filein = "buff"

    print("\nThe default reaction is d(28Si,p)")
    yn2 = input("\nWould you like to input a new reaction? (Y/N)")

    if yn2 == "N" or yn2 == "n":
        while filein[-4:] != ".dat" and filein[-4:] != ".txt":
            filein = input("\nEnter the name of the new input file (.dat or .txt): ")
            if filein[-4:] != ".dat" and filein[-4:] != ".txt":
                print("\nERROR: Incorrect file extension...")
    else:
        ebeam = BuildEvts()

print("\nThe file to be used is: " + filein)

pipeyn = input("Would you like to use an ISO 160 pipe for the gas return? [Y/N] ")

if pipeyn == "N" or pipeyn == "n":

    print("If you haven't already, please run PipeMaker.py to generate a pipe input file.")

    pfileyn = input("Would you like to use the default pipe file PipeOut_3_260.txt? [Y/N] ")

    if pfileyn == "N" or pfileyn == "n":
        pipedeffile = input("Enter the name of the pipe input file: ")

    else:
        pipedeffile = "PipeOut_3_260.txt"

    file = open(pipedeffile, "r")

    lines = file.readlines()
else:
    sim(1,1,1,1,1,ebeam,filein)