import math
import random
from massreader import readmass

def BuildEvts():

    masses = readmass()

    utoMeV = 931.4941

    for i in range(4):
        masses[i] = masses[i] * utoMeV

    beamenergy = float(input("Enter a beam energy in MeV: "))

    levnum = int(input("Enter the number of energy levels to be populated in the recoil: "))

    levels = []

    for i in range(levnum):
        levels.append(float(input("Enter the energy of a level: ")))

    numevents = int(input("Input the number of events you would like to generate: "))

    file = open("eventsout.txt", "w+")

    for i in range(numevents):

        randlevnum = random.randrange(0, levnum, 1)

        qval = (masses[0] + masses[1] - masses[2] - masses[3]) - levels[randlevnum]

        # print(qval)

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

        return beamenergy


    file.close()


BuildEvts()