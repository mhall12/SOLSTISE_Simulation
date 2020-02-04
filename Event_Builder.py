import math
import random
from massreader import readmass

def BuildEvts():

    # Have the user enter a reaction of the form Target(Beam, Ejectile). The recoil is calculated.
    reac = input("Enter a reaction of the form d(17F,p): ")

    # Let the user say whether or not the reaction is to be measured in inverse of normal kin.
    kinemat = 3
    while kinemat > 2 or kinemat ==0:
        kinemat = int(input("Will the reaction be measured in normal (1) or inverse kinematics (2)? "))
        if kinemat > 2 or kinemat == 0:
            print("ERROR: Please choose 1 for normal kinematics or 2 for inverse kinematics...")

    # Go into readmass to get the masses of the four particles from the entered reaction.
    masses = readmass(reac)

    utoMeV = 931.4941

    # Now that we have the particles masses we can get their mass in MeV:
    for i in range(4):
        masses[i] = masses[i] * utoMeV

    beamenergy = float(input("Enter a beam energy in MeV: "))

    # Entering 0 makes the program randomly generate excitation energies within the specified range.
    levnum = int(input("Enter the number of energy levels to be populated in the recoil (or 0 for "
                       "all energies): "))

    energyend = 0

    # Let the user specify the max excitation energy.
    if levnum == 0:
        print("All energies will be simulated.")
        energyend = int(input("Enter the highest excitation energy (integer MeV) you would like to simulate: "))

    # Or, have the user specify the number of levels they want simulated:
    while levnum < 0:
        print("ERROR: The number of levels cannot be less than 0.")
        levnum = int(input("Enter the number of energy levels to be populated in the recoil: "))

    levels = []

    # and their energies:
    if levnum > 0:
        for i in range(levnum):
            levels.append(float(input("Enter the energy of a level in MeV: ")))

    numevents = int(input("Input the number of events you would like to generate: "))

    # Take the reaction and replace all the parentheses with underscores so we can use it in the text file name
    reac = reac.replace("(", "_")
    reac = reac.replace(",", "_")
    reac = reac.replace(")", "_")

    if levnum > 0:
        outfilename = reac + str(int(beamenergy)) + "_evts.txt"
    else:
        outfilename = reac + str(int(beamenergy)) + "_evts_allE.txt"

    file = open(outfilename, "w+")

    # dead is for debugging and can be ignored for now.
    thmindeg = 0
    dead = 0

    for i in range(numevents):
        # randomly get a level nnumber if the levels were specified and then its excitation E
        if levnum > 0:
            randlevnum = random.randrange(0, levnum, 1)
            excitation_energy = levels[randlevnum]
        else:
            # otherwise randomly choose an excitation energy within the specified range.
            excitation_energy = random.randrange(0, energyend*100, 1)
            excitation_energy = excitation_energy/100

        # calculate the Q-Ex
        qval = (masses[0] + masses[1] - masses[2] - masses[3]) - excitation_energy

        # Here I'm trying to calculate the min theta that the particles can actually come out at from the math
        if masses[2] > masses[0]:
            thmin = math.acos(-1*math.sqrt(-((masses[3] + masses[2]) * (masses[3] * qval + (masses[3] - masses[1]) *
                                                                        beamenergy))/(masses[1] * masses[2] * beamenergy)))
            thmindeg = thmin * 180 / math.pi
        else:
            thmindeg = 90

        #print(thmindeg)

        # get a random angle in inverse or normal kinematics.
        if kinemat == 1:
            theta = random.random() * 90
        elif kinemat == 2:
            theta = random.random() * (180-thmindeg) + thmindeg

        # Convert the random degree angle to radians
        trad = theta * math.pi / 180

        # This while actually doesn't need to be here and can be removed in the future. At one point it
        # was used for debugging.
        while True:
            try:
                # calculate the ejectile energy here based on the angle and q-value.
                Eejec2 = ((math.sqrt(masses[1] * masses[2] * beamenergy) * math.cos(trad) +
                           math.sqrt(masses[1] * masses[2] * beamenergy * math.cos(trad)**2 +
                                    (masses[3] + masses[2]) * (masses[3] * qval + (masses[3] - masses[1]) *
                                                               beamenergy))) / (masses[3] + masses[2]))**2
                break
            except ValueError:
                Eejec2 = 0
                Theta = 90
                break

        # Smear the results a little to make it look like there is energy loss.
        Eejec2 = random.gauss(Eejec2, .01)
        theta = random.gauss(theta, .1)

        # Write the results to a text file so we can read it in later.
        file.write(str(theta) + '\t' + str(Eejec2) + '\n')

    # Close the file and return the file name
    file.close()
    return outfilename


if __name__ == "__main__":
    # We can also run this code separately by running >>python3 Event_Builder.py
    BuildEvts()
