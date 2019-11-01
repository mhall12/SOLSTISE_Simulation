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

    beamenergy = float(input("Enter a beam energy in MeV: "))

    levnum = int(input("Enter the number of energy levels to be populated in the recoil: "))

    while levnum < 1:
        print("ERROR: The number of levels cannot be less than 1.")
        levnum = int(input("Enter the number of energy levels to be populated in the recoil: "))

    levels = []

    for i in range(levnum):
        levels.append(float(input("Enter the energy of a level: ")))

    numevents = int(input("Input the number of events you would like to generate: "))

    reac = reac.replace("(", "_")
    reac = reac.replace(",", "_")
    reac = reac.replace(")", "_")

    outfilename = reac + str(int(beamenergy)) + "_evts.txt"

    file = open(outfilename, "w+")

    for i in range(numevents):

        randlevnum = random.randrange(0, levnum, 1)

        qval = (masses[0] + masses[1] - masses[2] - masses[3]) - levels[randlevnum]

        # print(qval)

        if kinemat == 1:
            theta  = random.random() * 90
        elif kinemat == 2:
            theta = random.random() * 90 + 90

        trad = theta * math.pi / 180

        pm = random.randrange(-1, 2, 2)

        Eejec2 = ((math.sqrt(masses[1] * masses[2] * beamenergy) * math.cos(trad) +
                    math.sqrt(masses[1] * masses[2] * beamenergy * math.cos(trad)**2 +
                            (masses[3] + masses[2]) * (masses[3] * qval + (masses[3] - masses[1]) *
                                                      beamenergy))) / (masses[3] + masses[2]))**2

        Eejec2 = random.gauss(Eejec2, .01)
        theta = random.gauss(theta, .1)

        file.write(str(theta) + '\t' + str(Eejec2) + '\n')

    return outfilename


    file.close()


#BuildEvts()