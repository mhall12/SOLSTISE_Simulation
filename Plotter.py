import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle
from massreader import readmass
import pandas as pd
import scipy.ndimage as ndimage
from mpl_toolkits import mplot3d
import glob
import os

def plot(pklin):

    # Read the pickled DataFrame File here.
    df = pd.read_pickle(pklin)

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
                                    "5: Blocked counts vs z  in all 4 detectors (1D)\n"
                                    "6: Blocked counts vs z  in all 4 detectors stacked (1D)\n"
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

                if (plotnum == 13 or plotnum == 12) and i == 0:
                    # Make Energy vs Theta contour plot here. Theta goes from 90 to 180 and we'll use bins every 5
                    # degrees. So,
                    binstheta = np.zeros(91)
                    for j in range(91):
                        binstheta[j] = j * 1 + 90
                    # We'll also make the energy bins as well:
                    binse = np.zeros(201)
                    for j in range(201):
                        binse[j] = j * .05

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

                    ratio = np.divide(blockedevt, unblockedevt, out=np.zeros_like(blockedevt), where=blockedevt != 0)
                    ratio2 = np.divide((unblockedevt - blockedevt), unblockedevt, out=np.zeros_like(blockedevt),
                                       where=unblockedevt != 0)

                    rtest = np.divide(unblockedevt, unblockedevt, out=np.zeros_like(unblockedevt),
                                      where=unblockedevt != 0)
                    # ratio = (unblockedevt-blockedevt)/unblockedevt
                    ratiocone = np.divide((unblockedevt - coneblockedevt), unblockedevt,
                                          out=np.zeros_like(coneblockedevt), where=unblockedevt != 0)

                    ratio = ratio.T
                    ratio2 = ratio2.T
                    rtest = rtest.T
                    ratiocone = ratiocone.T

                    tbins2 = np.zeros(90)
                    ebins2 = np.zeros(200)

                    for k in range(200):
                        if k < 90:
                            tbins2[k] = (tbins[k] + tbins[k + 1]) / 2
                        ebins2[k] = (ebins[k] + ebins[k + 1]) / 2

                    X, Y = np.meshgrid(tbins2, ebins2)
                    ratio2_blurr = ndimage.gaussian_filter(ratio2, sigma=1.2, order=0)
                    # X2, Y2 = np.meshgrid(tbins, ebins)

                    if plotnum == 12:
                        # fig, ax = plt.subplots()
                        cf = plt.contourf(X, Y, ratio2, 20, cmap='BuGn', linewidths=1.2)
                        cs = plt.contour(X, Y, ratio2, 10, colors='k', linewidths=1.2)
                        manual_locations = [(92, 5), (94, 5), (97, 5), (105, 5), (108, 5)]
                        plt.clabel(cs, inline=1, fontsize=14, manual=manual_locations)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                    if plotnum == 13:
                        # fig, ax = plt.subplots()
                        cf = plt.contourf(X, Y, ratiocone, 20, cmap='RdPu')
                        cs = plt.contour(X, Y, ratiocone, 10, colors='k')
                        manual_locations = [(92, 5), (94, 5), (97, 5), (105, 5)]
                        plt.clabel(cs, inline=1, fontsize=14, manual=manual_locations)
                        cbar = fig.colorbar(cf)
                        # ax.colobar()
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')

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
