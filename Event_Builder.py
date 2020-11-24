import math
import random
from massreader import readmass
from stopyt import desorb
import pickle
import pandas as pd
import numpy as np
import os


def BuildEvts():

    if not os.path.exists("Event_Files"):
        print("An Event_Files directory was not found, so one was created for you!")
        os.mkdir("Event_Files")

    evtdir = "Event_Files"

    # Have the user enter a reaction of the form Target(Beam, Ejectile). The recoil is calculated.
    reac = input("Enter a reaction of the form d(17F,p): ")

    # Let the user say whether or not the reaction is to be measured in inverse of normal kin.
    kinemat = 3
    while kinemat > 2 or kinemat == 0:
        try:
            kinemat = int(input("Will the reaction be measured upstream (1) or downstream (2) of the target? "))
        except ValueError:
            kinemat = 1
        if kinemat > 2 or kinemat == 0:
            print("ERROR: Please choose 1 for upstream or 2 for downstream...")

    # Go into readmass to get the masses of the four particles from the entered reaction.
    masses, ztarget, atarget, zejectile, aejectile, zbeam, abeam = readmass(reac)

    utoMeV = 931.4941

    # Now that we have the particles masses we can get their mass in MeV:
    for i in range(4):
        masses[i] = masses[i] * utoMeV

    beamenergy = float(input("Enter a beam energy in MeV (or enter 0 for MeV/u): "))
    if beamenergy == 0:
        beamenergyu = float(input("Enter a beam energy in MeV/u: "))
        beamenergy = beamenergyu * abeam

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

    numevents = int(input("\nInput the number of events you would like to generate: "))

    beamenp = np.zeros(numevents) + beamenergy
    angstrag = np.array([0])

    # We now have an energy loss code, so we need to see if the user wants artificial
    # smearing if the allE option isn't used.
    # Elossopt also initialized for the allE option.
    elossopt = 0
    if levnum > 0:
        try:
            elossopt = int(input("If you want to calculate energy loss in the code, "
                                 "enter (1), otherwise enter (0): "))
        except ValueError:
            elossopt = 1

    # Here we need to calculate the energy loss of the beam in the magnet and target. Moving that section of code
    # from SOLSTISE_Sim.py:
    if elossopt == 1:
        # Have them define the target...
        targetparms = []
        # need to set up the absorber data here, which should be easy because we know the target:

        if ztarget == 1:
            gs = int(input("\nIs the target a gas (1) or solid (2)?: "))
            if gs == 2:
                # depth is a random number between 0 and 1, which when multiplied by the thickness of the absorber
                # gives the position in the target that the reaction occurs.
                depth = np.random.rand(numevents)
                gasorsolid = "s"
                if atarget == 1:
                    print("A CH2 target will be assumed.")
                    zt = [6, 1]
                    at = [12, 1]
                    num = [1, 2]
                elif atarget == 2:
                    print("A CD2 target will be assumed.")
                    zt = [6, 1]
                    at = [12, 2]
                    num = [1, 2]
                density = [0.94]
                thickness = [float(input("Enter the thickness in mg/cm^2: "))]
                # Not assuming center of target anymore, so multiply the thickness by the depth.
                thkin = thickness[0] * depth
                jetpress = [0]
                champress = [0]
                jetrad = [0]
                chamdist = [0]
                gas = [False]
        else:
            gs = 1
            print("\nA gas target will be assumed...")

        if gs == 1:
            gasorsolid = "g"
            zt = [ztarget]
            at = [atarget]
            if ztarget == 1:
                num = [2]
            else:
                num = [int(input("Enter the number of atoms in the molecule: "))]

            # Density and thickness are 0 for gasses.
            density = [0]
            thickness = [0]
            thkin = [0]
            jetpress = [float(input("Enter the pressure in the jet in Torr (typically ~400 Torr): "))]
            jetrad = [float(input("Enter the radius of the jet in mm (typically 1.0 to 1.75 mm): "))]
            champress = [float(input("Enter the ambient pressure in the magnet in Torr (typically less than 1.0 Torr)"
                                     ": "))]
            chamdist = [float(input("Enter the distance from the end of the magnet to the jet in cm \n"
                                    "(distance to the center of HELIOS (SOLARIS) is ~117 cm (~136 cm)): "))]
            # Convert jetrad into cm
            jetrad[0] = jetrad[0] / 10
            jetradin = np.zeros(numevents) + jetrad[0]
            # Reaction takes place
            # in a guss position on the jet diameter, centered on the jet radius
            # jetradin will be the gauss position, but we also must "fix" it to remove negative values and values
            # greater than the jet diameter. Do it with masks.
            jetradin = np.random.normal(jetradin, (jetrad[0] * 2) / 2.355)
            # save the interaction position, because we'll need that in the text file. Currently in cm.
            jetradin_2 = jetradin
            ltomask = jetradin < 0
            gtdmask = jetradin > jetrad[0] * 2
            jetradin[ltomask] = 0.0001
            jetradin[gtdmask] = jetrad[0] * 2
            gas = [True]

        targetparms = [zt, at, num, density, thickness, jetpress, jetrad, champress, gas]
        dfout = pd.DataFrame()

        zbeam = np.zeros(numevents) + zbeam
        abeam = np.zeros(numevents) + abeam
        beame = beamenp
        beamei = beamenp
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
                leng = jetradin

            dfout = desorb(zbeam, abeam, beame, zt, at, num, gas[0], density[0], thkin, prs, leng, beamei)
            beame = dfout['Energy_i'].to_numpy() - dfout['DeltaE_tot'].to_numpy()
            # We have two layers, so get the total sum of the energy and angular straggling:
            beamstrag = dfout['E_strag_FWHM'].to_numpy() + beamstrag
            angstrag = dfout['AngleStrag'].to_numpy() + angstrag

            # print(beame)
            # print(np.average(beame))

    else:
        if levnum > 0:
            print("\nThe energies and angles of the particles will be artificially smeared.")

    if elossopt == 1:
        beame2 = np.random.normal(beame, beamstrag)
        targetparms.append([np.average(beame)])
    else:
        beame2 = beamenp

    # Take the reaction and replace all the parentheses with underscores so we can use it in the text file name
    reac = reac.replace("(", "_")
    reac = reac.replace(",", "_")
    reac = reac.replace(")", "_")

    if levnum > 0 and elossopt == 0:
        outfilename = "./" + evtdir + "/" + reac + str(int(beamenergy)) + "_evts_artsm.txt"
    elif levnum > 0 and elossopt == 1:
        outfilename = "./" + evtdir + "/" + reac + str(int(beamenergy)) + "_evts_eloss_" + gasorsolid + ".txt"
    else:
        outfilename = "./" + evtdir + "/" + reac + str(int(beamenergy)) + "_evts_allE.txt"

    file = open(outfilename, "w+")

    # dead is for debugging and can be ignored for now.
    thmindeg = 0
    dead = 0

    # Append 3 seed events on to beame2 for Ex calculation later at the beginning of the array
    avgbeame = np.average(beame2)
    for i in range(3):
        beame2 = np.insert(beame2, 0, avgbeame)

    numevents = numevents + 3

    # The beam energy is now a numpy array, so no longer need the for loop here...
    # But we have to assign some level numbers:

    if levnum > 0:
        # Assigns a random level number to each of the events.
        nplevnums = np.random.randint(levnum, size=numevents)

        # Set the level numbers for the first 3 events used for the excitation energy spectrum:
        nplevnums[0] = 0
        nplevnums[1] = 0
        if levnum > 2:
            nplevnums[2] = levnum - 2
        else:
            nplevnums[2] = levnum - 1

        exenergy = np.zeros(len(nplevnums))
        for i in range(levnum):
            indmask = nplevnums == i
            exenergy[indmask] = levels[i]

    if levnum == 0:
        # In this case, generate a number between 0 and 1, and multiply it by the energyend (highest desired Ex).
        exenergy = np.random.rand(numevents)

        exenergy = exenergy * energyend

    qval = (masses[0] + masses[1] - masses[2] - masses[3]) - exenergy

    # Here I'm trying to calculate the min theta that the particles can actually come out at from the math
    if masses[2] > masses[0]:
        thmin = np.arccos(-1 * np.sqrt(-((masses[3] + masses[2]) *
                                       (masses[3] * qval + (masses[3] - masses[1]) * beame2)) /
                                     (masses[1] * masses[2] * beame2)))
        thmindeg = thmin * 180 / np.pi
    else:
        thmindeg = np.zeros(numevents) + 90

    # Normal kinematics, random angle between 0 and 90
    if kinemat == 2:
        theta = np.random.rand(numevents) * 90

        # Hardcode angles in for the first 3 events
        theta[0] = 75
        theta[1] = 25
        theta[2] = 45
    # Inverse kinematics, random angle between 90 and 180
    if kinemat == 1:
        theta = np.random.rand(numevents) * (180-thmindeg) + thmindeg
        # Hardcode angles in for the first 3 events
        theta[0] = 95
        theta[1] = 115
        theta[2] = 115

    # Also need to hard code in the angle straggling, which is used in a normal dist later
    for i in range(3):
        angstrag = np.insert(angstrag, 0, 0, axis=0)

    trad = theta * np.pi / 180

    # calculate the ejectile energy here based on the angle and q-value.
    eejec2 = ((np.sqrt(masses[1] * masses[2] * beame2) * np.cos(trad) +
               np.sqrt(masses[1] * masses[2] * beame2 * np.cos(trad) ** 2 +
                         (masses[3] + masses[2]) * (masses[3] * qval + (masses[3] - masses[1]) *
                                                    beame2))) / (masses[3] + masses[2])) ** 2

    # Smear the results a little to make it look like there is energy loss.
    # or if the user has put elossopt = 0, we don't want the gaussian smearing.
    if elossopt == 0:
        eejec2 = np.random.normal(eejec2, .01)
        theta = np.random.normal(theta, .1)
    else:
        # If the energyloss is calculated, we need to smear the theta of the final particles by the angular
        # straggling.
        theta = np.random.normal(theta, angstrag)

    if elossopt == 0:
        results = np.array([theta, eejec2])
    else:
        # Need to create a pickle of targetparms so we can open it again later.
        targetparms.append([exenergy[0]])
        targetparms.append([exenergy[2]])
        pklname = "./" + evtdir + "/" + reac + str(int(beamenergy)) + "_tgt_" + gasorsolid + ".pkl"
        with open(pklname, 'wb') as f:
            pickle.dump(targetparms, f)

        if gas[0]:
            # jetradin_2 comes out in cm
            for i in range(3):
                jetradin_2 = np.insert(jetradin_2, 0, jetrad[0])
            results = np.array([theta, eejec2, jetradin_2])
        else:
            # thkin comes out as mg/cm^2
            for i in range(3):
                thkin = np.insert(thkin, 0, thickness[0] / 2)
            results = np.array([theta, eejec2, thkin])

    # Write the numpy results array into a text file.
    np.savetxt(outfilename, results.T, delimiter='\t')

    return outfilename[14:]


if __name__ == "__main__":
    # We can also run this code separately by running >>python3 Event_Builder.py
    BuildEvts()
