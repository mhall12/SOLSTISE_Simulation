import math
import random
from massreader import readmass

def BuildEvts():

    reac = input("Enter a reaction of the form d(17F,p): ")

    kinemat = 3
    while kinemat > 2 or kinemat ==0:
        kinemat = int(input("Will the reaction be measured in normal (1) or inverse kinematics (2)? "))
        if kinemat > 2 or kinemat == 0:
            print("ERROR: Please choose 1 for normal kinematics or 2 for inverse kinematics...")


    masses = readmass(reac)

    utoMeV = 931.4941

    for i in range(4):
        masses[i] = masses[i] * utoMeV
    #    print(masses[i])

    beamenergy = float(input("Enter a beam energy in MeV: "))

    levnum = int(input("Enter the number of energy levels to be populated in the recoil (or 0 for "
                       "all energies): "))

    energyend = 0

    if levnum == 0:
        print("All energies will be simulated.")
        energyend = int(input("Enter the highest excitation energy (integer MeV) you would like to simulate: "))

    while levnum < 0:
        print("ERROR: The number of levels cannot be less than 0.")
        levnum = int(input("Enter the number of energy levels to be populated in the recoil: "))

    levels = []

    if levnum > 0:
        for i in range(levnum):
            levels.append(float(input("Enter the energy of a level in MeV: ")))

    numevents = int(input("Input the number of events you would like to generate: "))

    reac = reac.replace("(", "_")
    reac = reac.replace(",", "_")
    reac = reac.replace(")", "_")

    if levnum > 0:
        outfilename = reac + str(int(beamenergy)) + "_evts.txt"
    else:
        outfilename = reac + str(int(beamenergy)) + "_evts_allE.txt"

    file = open(outfilename, "w+")

    thmindeg = 0
    dead = 0

    for i in range(numevents):

        if levnum > 0:
            randlevnum = random.randrange(0, levnum, 1)
            excitation_energy = levels[randlevnum]
        else:
            excitation_energy = random.randrange(0, energyend*100, 1)
            excitation_energy = excitation_energy/100

        qval = (masses[0] + masses[1] - masses[2] - masses[3]) - excitation_energy

        #print(qval)
        #print(excitation_energy)
        if (masses[2] > masses[0]):
            thmin = math.acos(-1*math.sqrt(-((masses[3] + masses[2]) * (masses[3] * qval + (masses[3] - masses[1]) *
                                                      beamenergy))/(masses[1] * masses[2] * beamenergy)))
            thmindeg = thmin * 180 / math.pi
        else:
            thmindeg = 90

        #print(thmindeg)

        if kinemat == 1:
            theta = random.random() * 90
        elif kinemat == 2:
            theta = random.random() * (180-thmindeg) + thmindeg

        #print(theta)
        trad = theta * math.pi / 180

        pm = random.randrange(-1, 2, 2)

        while True:
            try:
                Eejec2 = ((math.sqrt(masses[1] * masses[2] * beamenergy) * math.cos(trad) +
                            math.sqrt(masses[1] * masses[2] * beamenergy * math.cos(trad)**2 +
                                    (masses[3] + masses[2]) * (masses[3] * qval + (masses[3] - masses[1]) *
                                                              beamenergy))) / (masses[3] + masses[2]))**2
                break
            except ValueError:
                dead = dead+1
                Eejec2 = 0
                Theta = 90
                break

        Eejec2 = random.gauss(Eejec2, .01)
        theta = random.gauss(theta, .1)


        file.write(str(theta) + '\t' + str(Eejec2) + '\n')

    print(dead)
    file.close()
    return outfilename



#BuildEvts()