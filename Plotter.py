import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle
import pandas as pd
import scipy.ndimage as ndimage
from mpl_toolkits import mplot3d
import glob
import os

def plot(pklin):

    # Read the pickled DataFrame File here.
    df = pd.read_pickle(pklin)

    print("\n******************************************************\n"
          "The Reacton is: " + df['Reaction'][0] +
          "\nBeam Energy: " + str(df['Beam Energy'][0]) + " MeV" +
          "\nMagnetic Field: " + str(df['Magnetic Field'][0]) + " Tesla" +
          "\nReaction Distance from Nozzle: " + str(df['Reaction Distance from Nozzle'][0]) + ' inches' +
          "\nNozzle-Cone Distance: " + str(df['Nozzle-Cone Distance'][0]) + ' inches'
          "\nMagnet Bore Radius: " + str(df['Bore Radius'][0]) + ' m' +
          "\nCustom Pipe? " + str(df['Custom Pipe?'][0]) +
          "\nPipe Radius: " + str(df["Pipe Radius"][0]) + " m"
          "\nPipe Left Edge Angle: " + str(df['Pipe Left Edge Angle'][0]) + " Degrees" +
          "\nPipe Right Edge Angle: " + str(df['Pipe Right Edge Angle'][0]) + " Degrees" +
          "\nReceiver Cone Opening Diameter: " + str(df['Cone Opening Diameter'][0]) + " inches" +
          "\nReceiver Cone Height: " + str(df["Cone Height"][0]) + " inches")

    # Create the detector array here:
    detarr = [df["Det2"], df["Det1"], df["Det3"], df["Det4"]]

    # Determine whether or not the particles are coming out at backward or forward angles for plotting purposes:
    if df["Theta_Deg"].mean() > 90:
        invkin = True
    else:
        invkin = False

    Reds = cm.get_cmap('Reds', 256)
    newcolors = Reds(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 0])
    red = np.array([1, 0, 0, 1])
    newcolors[:1, :] = white
    newcolors[1:, :] = red
    newcmpRed = ListedColormap(newcolors)

    Blues = cm.get_cmap('Blues', 256)
    newcolors = Blues(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 0])
    blue = np.array([0, .5, 1, 1])
    newcolors[:1, :] = white
    newcolors[1:, :] = blue
    newcmpBlue = ListedColormap(newcolors)

    Greens = cm.get_cmap('Greens', 256)
    newcolors = Greens(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 0])
    green = np.array([0, 1, 0, 1])
    newcolors[:1, :] = white
    newcolors[1:, :] = green
    newcmpGreen = ListedColormap(newcolors)

    Greys = cm.get_cmap('Greys', 256)
    newcolors = Greys(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 0])
    black = np.array([0, 0, 0, 1])
    newcolors[:1, :] = white
    newcolors[1:, :] = black
    newcmpBlack = ListedColormap(newcolors)

    BuGn_r = cm.get_cmap('BuGn_r', 256)
    newcolors = BuGn_r(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 0])
    newcolors[:10, :] = white
    BuGn_new = ListedColormap(newcolors)

    blk = newcmpBlack(1)
    grn = newcmpGreen(1)
    blu = newcmpBlue(1)
    red = newcmpRed(1)

    switch = 0

    # Histogram order to get the quadrants right.
    ho = [1, 0, 2, 3]

    # Set the max and min values for the various histogram axis parameters
    if invkin:
        zmax = 0
        zmin = df['zpos_final'].min()
    else:
        zmin = 0
        zmax = df['zpos_final'].max()

    emax = df['Energy'].max()
    # emin should always just be 0

    thmax = df['Theta_Deg'].max()
    thmin = df['Theta_Deg'].min()

    print("\n\nChoose from the list below to plot histograms from the generated data.\n"
          "The four histograms represent the four quadrants of a fictional cylindrical detector\n"
          "looking down the beam axis.")

    while switch == 0:
        plt.ion()
        # plt.pause(0.0001)
        while True:
            try:
                plotnum = int(input("\n1: Unblocked particles in all 4 detectors (2D).\n"
                                    "2: Blocked particles in all 4 detectors (2D).\n"
                                    "3: Unblocked and blocked particles in all 4 detectors (2D).\n"
                                    "4: Total blocked counts vs angle (1D).\n"
                                    "5: Blocked counts vs z in all 4 detectors (1D)\n"
                                    "6: Blocked counts vs z in all 4 detectors stacked (1D)\n"
                                    "7: Blocked counts vs lab angle in all 4 detectors stacked (1D)\n"
                                    "8: Unblocked particles Energy vs ejected lab angle in all 4 detectors (2D)\n"
                                    "9: Unblocked and blocked particles Energy vs ejected lab angle in all 4 detectors "
                                    "(2D)\n"
                                    "10: Unblocked particles Ex from detected energy and position (1D)\n"
                                    "11: Unblocked and blocked particles Ex from detected energy and position (1D)\n"
                                    "0: End\n"
                                    "Enter a Histogram Number: "))
                break
            except:
                print("\n*****Enter an integer number from the list!*****\n")

        if plotnum > 0:
            fig = plt.figure()

            plt.rc('axes', labelsize=15)
            plt.rc('xtick', labelsize=15)
            plt.rc('ytick', labelsize=15)
            plt.rcParams["font.family"] = "STIXGeneral"

            for i in range(4):
                if plotnum == 1 or plotnum == 3:
                    plt.subplot(2, 2, i + 1)
                    plt.hist2d(df['zpos_final'][detarr[i] & df["Unblocked"]],
                               df['Energy'][detarr[i] & df["Unblocked"]], bins=(750, 750),
                               range=[[zmin, zmax], [0, emax]], cmap=newcmpBlack)
                    plt.xlabel('z(m)')
                    plt.ylabel('Energy (MeV)')

                if plotnum == 2 or plotnum == 3:
                    plt.subplot(2, 2, i + 1)
                    plt.hist2d(df['zpos_final'][detarr[i] & df["Blocked_Cone"]],
                               df['Energy'][detarr[i] & df["Blocked_Cone"]], bins=(750, 750),
                               range=[[zmin, zmax], [0, emax]], cmap=newcmpGreen)
                    plt.hist2d(df['zpos_final'][detarr[i] & df["Blocked_Pipe"]],
                               df['Energy'][detarr[i] & df["Blocked_Pipe"]], bins=(750, 750),
                               range=[[zmin, zmax], [0, emax]], cmap=newcmpRed)
                    plt.hist2d(df['zpos_final'][detarr[i] & df["Blocked_Nozzle"]],
                               df['Energy'][detarr[i] & df["Blocked_Nozzle"]], bins=(750, 750),
                               range=[[zmin, zmax], [0, emax]], cmap=newcmpBlue)
                    plt.xlabel('z(m)')
                    plt.ylabel('Energy (MeV)')

                if plotnum == 4 and i == 0:
                    plt.hist(df['Theta_Deg'][df['Blocked_Cone'] | df['Blocked_Pipe'] | df['Blocked_Nozzle']])
                    plt.xlabel('Lab Angle (deg)')
                    plt.ylabel('Counts')

                if plotnum == 5:
                    plt.subplot(2, 2, i + 1)
                    plt.hist(df['zpos_final'][df['Blocked_Cone'] & detarr[i]], bins=375, range=[-0.5, 0],
                             color=grn, alpha=1)
                    plt.hist(df['zpos_final'][df['Blocked_Pipe'] & detarr[i]], bins=375, range=[-0.5, 0],
                             color=red, alpha=0.7)
                    plt.hist(df['zpos_final'][df['Blocked_Nozzle'] & detarr[i]], bins=375, range=[-0.5, 0],
                             color=blu, alpha=0.6)
                    plt.xlabel('z(m)')
                    plt.ylabel('Counts')

                if plotnum == 6:
                    plt.subplot(2, 2, i + 1)
                    plt.hist((df['zpos_final'][df['Blocked_Cone'] & detarr[i]],
                              df['zpos_final'][df['Blocked_Pipe'] & detarr[i]],
                              df['zpos_final'][df['Blocked_Nozzle'] & detarr[i]]),
                             bins=375, range=[-0.5, 0], color=(grn, red, blu), stacked=True)
                    plt.xlabel('z(m)')
                    plt.ylabel('Counts')

                if plotnum == 7:
                    plt.subplot(2, 2, i + 1)
                    plt.hist((df['Theta_Deg'][df['Blocked_Cone'] & detarr[i]],
                              df['Theta_Deg'][df['Blocked_Pipe'] & detarr[i]],
                              df['Theta_Deg'][df['Blocked_Nozzle'] & detarr[i]]),
                             bins=60, range=[90, 120], color=(grn, red, blu), stacked=True)
                    plt.xlabel('Lab Angle (Deg)')
                    plt.ylabel('Counts')

                if plotnum == 8 or plotnum == 9:
                    plt.subplot(2, 2, i + 1)
                    plt.hist2d(df['Theta_Deg'][df['Unblocked'] & detarr[i]],
                               df['Energy'][df['Unblocked'] & detarr[i]], bins=(750, 750), range=[[90, 180], [0, 11]],
                               cmap=newcmpBlack)
                    if plotnum == 9:
                        plt.hist2d(df['Theta_Deg'][df['Blocked_Cone'] & detarr[i]],
                                   df['Energy'][detarr[i] & df["Blocked_Cone"]], bins=(750, 750),
                                   range=[[90, 180], [0, 11]], cmap=newcmpGreen)
                        plt.hist2d(df['Theta_Deg'][df['Blocked_Pipe'] & detarr[i]],
                                   df['Energy'][detarr[i] & df["Blocked_Pipe"]], bins=(750, 750),
                                   range=[[90, 180], [0, 11]], cmap=newcmpRed)
                        plt.hist2d(df['Theta_Deg'][df['Blocked_Nozzle'] & detarr[i]],
                                   df['Energy'][detarr[i] & df["Blocked_Nozzle"]], bins=(750, 750),
                                   range=[[90, 180], [0, 11]], cmap=newcmpBlue)
                    plt.xlabel('Lab Angle (Deg)')
                    plt.ylabel('Energy (MeV)')

                if plotnum == 10 and i == 0:
                    plt.hist(df['Ex_Reconstructed'][df["Unblocked"]], bins=750, range=[0, 10])
                    plt.xlabel('Excitation Energy (MeV)')
                    plt.ylabel('Counts')

                if plotnum == 11:
                    plt.subplot(2, 2, i + 1)
                    plt.hist((df['Ex_Reconstructed'][df["Unblocked"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Cone"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Pipe"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Nozzle"] & detarr[i]]),
                             bins=1000, range=[0, 8], color=(blk, grn, red, blu), stacked=True)
                    plt.xlabel('Excitation Energy (MeV)')
                    plt.ylabel('Counts')

                if 11 < plotnum < 16 and i == 0:
                    # Make Energy vs Theta contour plot here. Theta goes from 90 to 180 and we'll use bins every 5
                    # degrees. So,
                    binstheta = np.zeros(91)
                    for j in range(91):
                        binstheta[j] = j * 1 + 90
                    # We'll also make the energy bins as well:
                    binse = np.zeros(151)
                    for j in range(151):
                        binse[j] = j * .0667

                    unblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['AllPossible']],
                                                                df['Energy'][df['AllPossible']],
                                                                bins=(binstheta, binse))

                    blockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Cone'] |
                                                                              df['Blocked_Pipe'] |
                                                                              df['Blocked_Nozzle']],
                                                              df['Energy'][df['Blocked_Cone'] |
                                                                           df['Blocked_Pipe'] |
                                                                           df['Blocked_Nozzle']],
                                                              bins=(binstheta, binse))

                    coneblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Cone']],
                                                                  df['Energy'][df['Blocked_Cone']],
                                                                  bins=(binstheta, binse))
                    pipeblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Pipe']],
                                                                  df['Energy'][df['Blocked_Pipe']],
                                                                  bins=(binstheta, binse))
                    nozzleblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Nozzle']],
                                                                    df['Energy'][df['Blocked_Nozzle']],
                                                                    bins=(binstheta, binse))

                    #ratio = np.divide(blockedevt, unblockedevt, out=np.zeros_like(blockedevt), where=blockedevt != 0)
                    ratio2 = np.divide((unblockedevt - blockedevt), unblockedevt, out=np.zeros_like(blockedevt),
                                       where=unblockedevt != 0)

                    ratiocone = np.divide((unblockedevt - coneblockedevt), unblockedevt,
                                          out=np.zeros_like(coneblockedevt), where=unblockedevt != 0)

                    ratiopipe = np.divide((unblockedevt - pipeblockedevt), unblockedevt,
                                          out=np.zeros_like(pipeblockedevt), where=unblockedevt != 0)

                    rationozzle = np.divide((unblockedevt - nozzleblockedevt), unblockedevt,
                                            out=np.zeros_like(nozzleblockedevt), where=unblockedevt != 0)

                    ratio2 = ratio2.T

                    ratiocone = ratiocone.T
                    ratiopipe = ratiopipe.T
                    rationozzle = rationozzle.T

                    tbins2 = np.zeros(90)
                    ebins2 = np.zeros(150)

                    for k in range(150):
                        if k < 90:
                            tbins2[k] = (tbins[k] + tbins[k + 1]) / 2
                        ebins2[k] = (ebins[k] + ebins[k + 1]) / 2

                    X, Y = np.meshgrid(tbins2, ebins2)
                    ratio2_blurr = ndimage.gaussian_filter(ratio2, sigma=1.5, order=0)
                    ratiocone_blurr = ndimage.gaussian_filter(ratiocone, sigma=1.5, order=0)
                    ratiopipe_blurr = ndimage.gaussian_filter(ratiopipe, sigma=1.5, order=0)
                    rationozzle_blurr = ndimage.gaussian_filter(rationozzle, sigma=1.5, order=0)

                    if plotnum == 12:

                        cf = plt.contourf(X, Y, ratio2_blurr, [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85,
                                                               .9, .95, 1], cmap='YlGnBu', linewidths=1.2)
                        cs = plt.contour(X, Y, ratio2_blurr, [.5, .6, .7, .8, .9], colors='k', linewidths=1.2)
                        manual_locations = [(92, 5), (94, 5), (97, 5), (104, 5), (110, 5)]
                        labels = plt.clabel(cs, inline=1, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected", fontdict={'fontsize': 16})
                    if plotnum == 13:

                        cf = plt.contourf(X, Y, ratiocone_blurr, [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85,
                                                               .9, .95, 1], cmap='Greens', linewidths=1.2)
                        cs = plt.contour(X, Y, ratiocone_blurr, [.5, .6, .7, .8, .9], colors='k', linewidths=1.2)
                        manual_locations = [(94, 5), (96, 5), (100, 5), (106, 5)]
                        labels = plt.clabel(cs, inline=1, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected (Not Blocked by the Cone)",
                                  fontdict={'fontsize': 16})

                    if plotnum == 14:

                        cf = plt.contourf(X, Y, ratiopipe_blurr,
                                          [.74, .78, .82, .86, .90, .94, .98, 1], cmap='Reds', linewidths=1.2)
                        cs = plt.contour(X, Y, ratiopipe_blurr, [.78, .86, .94], colors='k', linewidths=1.2)
                        manual_locations = [(98, 7.5), (95, 7.8), (92, 9)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-20, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected (Not Blocked by the Pipe)",
                                  fontdict={'fontsize': 16})

                    if plotnum == 15:

                        cf = plt.contourf(X, Y, rationozzle_blurr, [.80, .82, .84, .86, .88, .90, .92, .94,
                                                                    .96, .98, 1], cmap='Blues', linewidths=1.2)
                        cs = plt.contour(X, Y, rationozzle_blurr, [.88, .92, .96], colors='k', linewidths=1.2)
                        manual_locations = [(100, 7.5), (110, 6), (115, 4)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-10, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected (Not Blocked by the Nozzle)",
                              fontdict={'fontsize': 16})


                # Handle the legend here for each plot that needs it.
                if i == 0 and (plotnum == 1 or plotnum == 2 or plotnum == 3 or plotnum == 5 or
                               plotnum == 6 or plotnum == 7 or plotnum == 9 or plotnum == 11):
                    handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in [blk, grn, red, blu]]
                    labels = ["Unblocked", "Cone", "Pipe", "Nozzle"]
                    if plotnum == 11:
                        plt.legend(handles, labels, bbox_to_anchor=(1.0, .9), ncol=4)
                    else:
                        plt.legend(handles, labels, bbox_to_anchor=(1.4, 1.1), ncol=4)

        elif plotnum == 0:
            switch = 1

    input("\nPress ENTER to end.")


if __name__ == "__main__":
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
    print("                     Particle Shadowing Plotting Code")
    input("\n\n\nTo continue, press ENTER")

    list_pickles = glob.glob('*.pkl')

    if len(list_pickles) == 0:
        print("\nIt appears that no simulated DataFrame exists, run SOLSTISE_Sim.py first.")
    else:
        latest_file = max(list_pickles, key=os.path.getctime)
        print("\nThe most recently created simulated DataFrame file is: " + latest_file)
        yn = input("\nWould you like to use this file? [Y/N] ")
        if yn == "N" or yn == "n":
            if len(list_pickles) == 1:
                print("\nYou only have one DataFrame file. To create another, run SOLSTISE_Sim.py")
            if len(list_pickles) > 1:
                print("\n")
                for i in range(len(list_pickles)):
                    print(str(i + 1) + ") " + list_pickles[i])

                filenum = 1000000
                while filenum > len(list_pickles):
                    filenum = int(input("\nChoose a number from the list, or enter 0 to manually type the file name: "))
                    if len(list_pickles) >= filenum > 0:
                        filein = list_pickles[filenum-1]
                    elif filenum > len(list_pickles):
                        print("ERROR: Number entered is greater than the number of DataFrame files...")
        else:
            filein = latest_file

    print("\nThe file to be used is: " + filein)
    plot(filein)
