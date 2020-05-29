import math
import random
from massreader import readmass
from stopyt import desorb
import pickle
import pandas as pd
import numpy as np

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
    masses, ztarget, atarget, zejectile, aejectile, zbeam, abeam = readmass(reac)

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

    # We now have an energy loss code, so we need to see if the user wants artificial
    # smearing if the allE option isn't used.
    elossopt = 2
    if levnum > 0:
        elossopt = int(input("If you want to calculate energy loss in the code, "
                             "enter (1), otherwise enter (0): "))

    # Here we need to calculate the energy loss of the beam in the magnet and target. Moving that section of code
    # from SOLSTISE_Sim.py:
    if elossopt == 1:
        # Have them define the target...
        targetparms = []
        # need to set up the absorber data here, which should be easy because we know the target:

        if ztarget == 1:
            gs = int(input("\nIs the target a gas (1) or solid (2)?: "))
            if gs == 2:
                if atarget == 1:
                    print("\nA CH2 target will be assumed.")
                    zt = [6, 1]
                    at = [12, 1]
                    num = [1, 2]
                elif atarget == 2:
                    print("\nA CD2 target will be assumed.")
                    zt = [6, 1]
                    at = [12, 2]
                    num = [1, 2]
                density = [0.94]
                thickness = [float(input("\nEnter the thickness in mg/cm^2: "))]
                # Assume the beam interacts at the center of the target:
                thickness[0] = thickness[0] / 2
                jetpress = [0]
                champress = [0]
                jetrad = [0]
                chamdist = [0]
                gas = [False]
        else:
            gs = 1

        if gs == 1:
                zt = [ztarget]
                at = [atarget]
                num = [int(input("Enter the number of atoms in the molecule: "))]
                density = [0]
                thickness = [0]
                jetpress = [float(input("\nEnter the pressure in the jet in Torr: "))]
                jetrad = [float(input("\nEnter the radius of the jet in mm: "))]
                champress = [float(input("\nEnter the ambient pressure in the magnet in Torr: "))]
                chamdist = [float(input("\nEnter the distance from the end of the magnet to the jet in cm "
                                        "\n(distance to the center of HELIOS (SOLARIS) is ~117 cm (~136 cm)): "))]
                # Convert jetrad into cm
                jetrad[0] = jetrad[0] / 10
                gas = [True]

        targetparms = [zt, at, num, density, thickness, jetpress, jetrad, champress, gas]
        dfout = pd.DataFrame()

        zbeam = np.array([zbeam])
        abeam = np.array([abeam])
        beame = np.array([beamenergy])
        beamei = np.array([beamenergy])
        angstrag = 0
        beamstrag = 0
        # If the absorber is a gas, we need the for loop to run twice. If it's a solid, we need it to run once.
        if gas[0]:
            rng = 2
        else:
            rng = 1
        for i in range(rng):
            if i == 0:
                prs = champress[0]
                leng = chamdist[0]
            if i == 1:
                prs = jetpress[0]
                leng = jetrad[0]

            dfout = desorb(zbeam, abeam, beame, zt, at, num, gas[0], density[0], thickness[0], prs, leng, beamei)
            beame = dfout['Energy_i'] - dfout['DeltaE_tot']
            # We have two layers, so get the total sum of the energy and angular straggling:
            beamstrag = dfout['E_strag_FWHM'][0]
            angstrag = dfout['AngleStrag'][0] + angstrag

        print(beame)
        print(beamstrag)
        print(angstrag)
        targetparms.append([beame[0]])

    else:
        print("\nThe energies and angles of the particles will be artificially smeared.")

    numevents = int(input("\nInput the number of events you would like to generate: "))

    # Take the reaction and replace all the parentheses with underscores so we can use it in the text file name
    reac = reac.replace("(", "_")
    reac = reac.replace(",", "_")
    reac = reac.replace(")", "_")

    if levnum > 0 and elossopt == 0:
        outfilename = reac + str(int(beamenergy)) + "_evts.txt"
    elif levnum > 0 and elossopt == 1:
        outfilename = reac + str(int(beamenergy)) + "_evts_eloss.txt"
    else:
        outfilename = reac + str(int(beamenergy)) + "_evts_allE.txt"

    file = open(outfilename, "w+")

    # Need to create a pickle of targetparms so we can open it again later.
    if elossopt == 1:
        pklname = reac + str(int(beamenergy)) + "_tgt.pkl"
        with open(pklname, 'wb') as f:
            pickle.dump(targetparms, f)

    # dead is for debugging and can be ignored for now.
    thmindeg = 0
    dead = 0

    for i in range(numevents):

        # if the energy loss is calculated, we want the beam energy to be smeared by the energy straggling:
        if elossopt == 1:
            beamenergy = random.gauss(beame[0], beamstrag)

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
        # or if the user has put elossopt = 0, we don't want the gaussian smearing.
        if elossopt == 0:
            Eejec2 = random.gauss(Eejec2, .01)
            theta = random.gauss(theta, .1)
        else:
            # If the energyloss is calculated, we need to smear the theta of the final particles by the angular
            # straggling.
            theta = random.gauss(theta, angstrag)

        # Write the results to a text file so we can read it in later.
        file.write(str(theta) + '\t' + str(Eejec2) + '\n')

    # Close the file and return the file name
    file.close()
    return outfilename


if __name__ == "__main__":
    # We can also run this code separately by running >>python3 Event_Builder.py
    BuildEvts()
