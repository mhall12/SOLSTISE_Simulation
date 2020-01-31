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

    print("\n********************************Simulation Parameters*********************************\n"
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
          "\nReceiver Cone Height: " + str(df["Cone Height"][0]) + " inches"
          "\n**************************************************************************************\n")

    # Create the detector array here:
    detarr = [df["AllPossible"], df["Det2"], df["Det1"], df["Det3"], df["Det4"]]

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
          "looking down the beam axis.\n")

    print("\n1) Energy vs z: Unblocked particles in all 4 detectors (2D).\n"
          "2) Energy vs z: Blocked particles in all 4 detectors (2D).\n"
          "3) Energy vs z: Unblocked and blocked particles in all 4 detectors (2D).\n"
          "4) Counts vs Lab Angle: Total blocked counts vs angle (1D).\n"
          "5) Counts vs z: Blocked counts vs z in all 4 detectors (1D)\n"
          "6) Counts vs z: Blocked counts vs z in all 4 detectors stacked (1D)\n"
          "7) Counts vs Lab Angle: Blocked counts vs lab angle in all 4 detectors stacked (1D)\n"
          "8) Energy vs Lab Angle: Unblocked particles Energy vs ejected lab angle in all 4 detectors (2D)\n"
          "9) Energy vs Lab Angle: Unblocked and blocked particles Energy vs ejected lab "
          "angle in all 4 detectors (2D)\n"
          "10) Counts vs Ex: Unblocked particles Ex from detected energy and position (1D)\n"
          "11) Counts vs Ex: Unblocked and blocked particles Ex from detected energy and position (1D)\n"
          "12) Energy vs Lab Angle: Contour plot of detected particles. \n"
          "13) Energy vs Lab Angle: Contour plot of particles not blocked by the cone. \n"
          "14) Energy vs Lab Angle: Contour plot of particles not blocked by the pipe. \n"
          "15) Energy vs Lab Angle: Contour plot of particles not blocked by the nozzle. \n"
          "16) Energy vs Lab Angle: Contour plot of detected particles in all 4 detectors. \n"
          "17) Energy vs Lab Angle: Contour plot of particles not blocked by the cone in all 4 detectors.. \n"
          "18) Energy vs Lab Angle: Contour plot of particles not blocked by the pipe in all 4 detectors.. \n"
          "19) Energy vs Lab Angle: Contour plot of particles not blocked by the nozzle in all 4 detectors.. \n"
          "20) Energy vs z: Contour plot of detected particles. \n"
          "21) Energy vs z: Contour plot of particles not blocked by the cone. \n"
          "22) Energy vs z: Contour plot of particles not blocked by the pipe. \n"
          "23) Energy vs z: Contour plot of particles not blocked by the nozzle. \n"
          "24) Energy vs z: Contour plot of detected particles in all 4 detectors.. \n"
          "25) Energy vs z: Contour plot of particles not blocked by the cone in all 4 detectors. \n"
          "26) Energy vs z: Contour plot of particles not blocked by the pipe in all 4 detectors. \n"
          "27) Energy vs z: Contour plot of particles not blocked by the nozzle in all 4 detectors. \n"
          "0) End\n\n")

    while switch == 0:
        plt.ion()
        # plt.pause(0.0001)
        while True:
            try:
                plotnum = int(input("Enter a Histogram Number: "))
                break
            except:
                print("\n*****Enter an integer number from the list!*****\n")

        if plotnum > 0:
            fig = plt.figure()

            plt.rc('axes', labelsize=15)
            plt.rc('xtick', labelsize=15)
            plt.rc('ytick', labelsize=15)
            plt.rcParams["font.family"] = "STIXGeneral"

            for i in range(5):
                if (plotnum == 1 or plotnum == 3) and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist2d(df['zpos_final'][detarr[i] & df["Unblocked"]],
                               df['Energy'][detarr[i] & df["Unblocked"]], bins=(750, 750),
                               range=[[zmin, zmax], [0, emax]], cmap=newcmpBlack)
                    plt.xlabel('z(m)')
                    plt.ylabel('Energy (MeV)')

                if (plotnum == 2 or plotnum == 3) and i > 0:
                    plt.subplot(2, 2, i)
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

                if plotnum == 5 and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist(df['zpos_final'][df['Blocked_Cone'] & detarr[i]], bins=375, range=[-0.5, 0],
                             color=grn, alpha=1)
                    plt.hist(df['zpos_final'][df['Blocked_Pipe'] & detarr[i]], bins=375, range=[-0.5, 0],
                             color=red, alpha=0.7)
                    plt.hist(df['zpos_final'][df['Blocked_Nozzle'] & detarr[i]], bins=375, range=[-0.5, 0],
                             color=blu, alpha=0.6)
                    plt.xlabel('z(m)')
                    plt.ylabel('Counts')

                if plotnum == 6 and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist((df['zpos_final'][df['Blocked_Cone'] & detarr[i]],
                              df['zpos_final'][df['Blocked_Pipe'] & detarr[i]],
                              df['zpos_final'][df['Blocked_Nozzle'] & detarr[i]]),
                             bins=375, range=[-0.5, 0], color=(grn, red, blu), stacked=True)
                    plt.xlabel('z(m)')
                    plt.ylabel('Counts')

                if plotnum == 7 and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist((df['Theta_Deg'][df['Blocked_Cone'] & detarr[i]],
                              df['Theta_Deg'][df['Blocked_Pipe'] & detarr[i]],
                              df['Theta_Deg'][df['Blocked_Nozzle'] & detarr[i]]),
                             bins=60, range=[90, 120], color=(grn, red, blu), stacked=True)
                    plt.xlabel('Lab Angle (Deg)')
                    plt.ylabel('Counts')

                if (plotnum == 8 or plotnum == 9) and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist2d(df['Theta_Deg'][df['Unblocked'] & detarr[i]],
                               df['Energy'][df['Unblocked'] & detarr[i]], bins=(750, 750), range=[[90, 180], [0, 11]],
                               cmap=newcmpBlack)
                    if plotnum == 9 and i > 0:
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

                if plotnum == 11 and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist((df['Ex_Reconstructed'][df["Unblocked"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Cone"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Pipe"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Nozzle"] & detarr[i]]),
                             bins=1000, range=[0, 8], color=(blk, grn, red, blu), stacked=True)
                    plt.xlabel('Excitation Energy (MeV)')
                    plt.ylabel('Counts')

                if 11 < plotnum < 35:
                    # Make Energy vs Theta contour plot here. Theta goes from 90 to 180 and we'll use bins every 5
                    # degrees. So,
                    binstheta = np.zeros(91)
                    for j in range(91):
                        binstheta[j] = j * 1 + 90
                    # We'll also make the energy bins as well:
                    binse = np.zeros(151)
                    for j in range(151):
                        binse[j] = j * .0667
                    binsz = np.zeros(101)
                    # We'll also make the z bins here:
                    for j in range(101):
                        binsz[j] = (100-j) * -.0085

                    # Can't have detarr[i] alone. It has to include AllPossible because the Det masks do not contain
                    # the mask that detemines whether or not the particle hit the magnet bore.

                    unblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][detarr[i] & df["AllPossible"]],
                                                                df['Energy'][detarr[i] & df["AllPossible"]],
                                                                bins=(binstheta, binse))

                    unblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][detarr[i] & df["AllPossible"]],
                                                                df['Energy'][detarr[i] & df["AllPossible"]],
                                                                bins=(binsz, binse))

                    blockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][(df['Blocked_Cone'] |
                                                                              df['Blocked_Pipe'] |
                                                                              df['Blocked_Nozzle']) &
                                                                              detarr[i]],
                                                              df['Energy'][(df['Blocked_Cone'] |
                                                                           df['Blocked_Pipe'] |
                                                                           df['Blocked_Nozzle']) &
                                                                           detarr[i]],
                                                              bins=(binstheta, binse))

                    blockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][(df['Blocked_Cone'] |
                                                                              df['Blocked_Pipe'] |
                                                                              df['Blocked_Nozzle']) &
                                                                              detarr[i]],
                                                              df['Energy'][(df['Blocked_Cone'] |
                                                                           df['Blocked_Pipe'] |
                                                                           df['Blocked_Nozzle']) &
                                                                           detarr[i]],
                                                              bins=(binsz, binse))

                    coneblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Cone'] & detarr[i]],
                                                                  df['Energy'][df['Blocked_Cone'] & detarr[i]],
                                                                  bins=(binstheta, binse))
                    pipeblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Pipe'] & detarr[i]],
                                                                  df['Energy'][df['Blocked_Pipe'] & detarr[i]],
                                                                  bins=(binstheta, binse))
                    nozzleblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Nozzle'] & detarr[i]],
                                                                    df['Energy'][df['Blocked_Nozzle'] & detarr[i]],
                                                                    bins=(binstheta, binse))

                    coneblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][df['Blocked_Cone'] & detarr[i]],
                                                                  df['Energy'][df['Blocked_Cone'] & detarr[i]],
                                                                  bins=(binsz, binse))
                    pipeblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][df['Blocked_Pipe'] & detarr[i]],
                                                                  df['Energy'][df['Blocked_Pipe'] & detarr[i]],
                                                                  bins=(binsz, binse))
                    nozzleblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][df['Blocked_Nozzle'] & detarr[i]],
                                                                    df['Energy'][df['Blocked_Nozzle'] & detarr[i]],
                                                                    bins=(binsz, binse))

                    #ratio = np.divide(blockedevt, unblockedevt, out=np.zeros_like(blockedevt), where=blockedevt != 0)
                    ratioevt = np.divide((unblockedevt - blockedevt), unblockedevt, out=np.zeros_like(blockedevt),
                                         where=unblockedevt != 0)

                    ratioevz = np.divide((unblockedevz - blockedevz), unblockedevz, out=np.zeros_like(blockedevz),
                                         where=unblockedevz != 0)

                    ratioconeevt = np.divide((unblockedevt - coneblockedevt), unblockedevt,
                                             out=np.zeros_like(coneblockedevt), where=unblockedevt != 0)

                    ratiopipeevt = np.divide((unblockedevt - pipeblockedevt), unblockedevt,
                                             out=np.zeros_like(pipeblockedevt), where=unblockedevt != 0)

                    rationozzleevt = np.divide((unblockedevt - nozzleblockedevt), unblockedevt,
                                               out=np.zeros_like(nozzleblockedevt), where=unblockedevt != 0)

                    ratioconeevz = np.divide((unblockedevz - coneblockedevz), unblockedevz,
                                             out=np.zeros_like(coneblockedevz), where=unblockedevz != 0)

                    ratiopipeevz = np.divide((unblockedevz - pipeblockedevz), unblockedevz,
                                             out=np.zeros_like(pipeblockedevz), where=unblockedevz != 0)

                    rationozzleevz = np.divide((unblockedevz - nozzleblockedevz), unblockedevz,
                                               out=np.zeros_like(nozzleblockedevz), where=unblockedevz != 0)

                    ratioevt = ratioevt.T
                    ratioevz = ratioevz.T

                    ratioconeevt = ratioconeevt.T
                    ratiopipeevt = ratiopipeevt.T
                    rationozzleevt = rationozzleevt.T

                    ratioconeevz = ratioconeevz.T
                    ratiopipeevz = ratiopipeevz.T
                    rationozzleevz = rationozzleevz.T

                    tbins2 = np.zeros(90)
                    ebins2 = np.zeros(150)
                    zbins2 = np.zeros(100)

                    for k in range(150):
                        if k < 90:
                            tbins2[k] = (tbins[k] + tbins[k + 1]) / 2
                        if k < 100:
                            zbins2[k] = (zbins[k] + zbins[k + 1]) / 2
                        ebins2[k] = (ebins[k] + ebins[k + 1]) / 2

                    xevt, yevt = np.meshgrid(tbins2, ebins2)
                    xevz, yevz = np.meshgrid(zbins2, ebins2)

                    ratioevt_blurr = ndimage.gaussian_filter(ratioevt, sigma=1.5, order=0)
                    ratioevz_blurr = ndimage.gaussian_filter(ratioevz, sigma=1.5, order=0)

                    ratioconeevt_blurr = ndimage.gaussian_filter(ratioconeevt, sigma=1.5, order=0)
                    ratiopipeevt_blurr = ndimage.gaussian_filter(ratiopipeevt, sigma=1.5, order=0)
                    rationozzleevt_blurr = ndimage.gaussian_filter(rationozzleevt, sigma=1.5, order=0)

                    ratioconeevz_blurr = ndimage.gaussian_filter(ratioconeevz, sigma=1.6, order=0)
                    ratiopipeevz_blurr = ndimage.gaussian_filter(ratiopipeevz, sigma=1.6, order=0)
                    rationozzleevz_blurr = ndimage.gaussian_filter(rationozzleevz, sigma=1.6, order=0)

                    if plotnum == 12 and i == 0:

                        cf = plt.contourf(xevt, yevt, ratioevt_blurr, [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8,
                                                                       .85, .9, .95, 1], cmap='YlGnBu')
                        cs = plt.contour(xevt, yevt, ratioevt_blurr, [.5, .6, .7, .8, .9], colors='k', linewidths=1.2)
                        manual_locations = [(92, 5), (94, 5), (97, 5), (104, 5), (110, 5)]
                        labels = plt.clabel(cs, inline=1, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected", fontdict={'fontsize': 16})
                    if plotnum == 13 and i == 0:

                        cf = plt.contourf(xevt, yevt, ratioconeevt_blurr, [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75,
                                                                           .8, .85, .9, .95, 1], cmapbb='Greens')
                        cs = plt.contour(xevt, yevt, ratioconeevt_blurr, [.5, .6, .7, .8, .9],
                                         colors='k', linewidths=1.2)
                        manual_locations = [(94, 5), (96, 5), (100, 5), (106, 5)]
                        labels = plt.clabel(cs, inline=1, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected (Not Blocked by the Cone)",
                                  fontdict={'fontsize': 16})

                    if plotnum == 14 and i == 0:

                        cf = plt.contourf(xevt, yevt, ratiopipeevt_blurr,
                                          [.65, .7, .75, .8, .85, .9, .95, 1], cmap='Reds')
                        cs = plt.contour(xevt, yevt, ratiopipeevt_blurr, [.7, .8, .9], colors='k', linewidths=1.2)
                        manual_locations = [(112, 5), (109, 6.2), (92, 9)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-20, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected (Not Blocked by the Pipe)",
                                  fontdict={'fontsize': 16})

                    if plotnum == 15 and i == 0:

                        cf = plt.contourf(xevt, yevt, rationozzleevt_blurr, [.80, .82, .84, .86, .88, .90, .92, .94,
                                                                             .96, .98, 1], cmap='Blues')
                        cs = plt.contour(xevt, yevt, rationozzleevt_blurr, [.88, .92, .96], colors='k', linewidths=1.2)
                        manual_locations = [(100, 7.5), (110, 6), (115, 4)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-10, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-90)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected (Not Blocked by the Nozzle)",
                                  fontdict={'fontsize': 16})

                    if plotnum == 16 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevt, yevt, ratioevt_blurr,
                                          [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85,
                                           .9, .95, 1], cmap='YlGnBu')
                        cs = plt.contour(xevt, yevt, ratioevt_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55,
                                                                      .6, .65, .7, .75, .8, .85,
                                           .9, .95, 1], colors='k', linewidths=.2)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles Detected', fontsize=18)

                    if plotnum == 17 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevt, yevt, ratioconeevt_blurr,
                                          [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85,
                                           .9, .95, 1], cmap='Greens')
                        cs = plt.contour(xevt, yevt, ratioconeevt_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5,
                                                                          .55, .6, .65, .7, .75, .8, .85,
                                           .9, .95, 1], colors='k', linewidths=.3)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles That Could Be Detected (Not Blocked by the Cone)',
                                     fontsize=18)

                    if plotnum == 18 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevt, yevt, ratiopipeevt_blurr,
                                          [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85,
                                           .9, .95, 1], cmap='Reds')
                        cs = plt.contour(xevt, yevt, ratiopipeevt_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5,
                                                                          .55, .6, .65, .7, .75, .8, .85,
                                           .9, .95, 1], colors='k', linewidths=.2)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles That Could Be Detected (Not Blocked by the Pipe)',
                                     fontsize=18)

                    if plotnum == 19 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevt, yevt, rationozzleevt_blurr,
                                          [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85,
                                           .9, .95, 1], cmap='Blues')
                        cs = plt.contour(xevt, yevt, rationozzleevt_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                            .5, .55, .6, .65, .7, .75, .8, .85,
                                                                            .9, .95, 1], colors='k', linewidths=.2)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles That Could Be Detected (Not Blocked by the Nozzle)',
                                     fontsize=18)

                    if plotnum == 20 and i == 0:

                        cf = plt.contourf(xevz, yevz, ratioevz_blurr, [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8,
                                                                       .85, .9, .95, 1], cmap='YlGnBu')
                        cs = plt.contour(xevz, yevz, ratioevz_blurr, [.5, .6, .7, .8, .9], colors='k', linewidths=1.2)
                        manual_locations = [(-.35, 5), (-.28, 5), (-.15, 5), (-.11, 5.5), (-.09, 5.5)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-15, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-70)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected", fontdict={'fontsize': 16})

                    if plotnum == 21 and i == 0:

                        cf = plt.contourf(xevz, yevz, ratioconeevz_blurr, [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75,
                                                                           .8, .85, .9, .95, 1], cmap='Greens')
                        cs = plt.contour(xevz, yevz, ratioconeevz_blurr, [.5, .6, .7, .8, .9], colors='k',
                                         linewidths=1.2)
                        manual_locations = [(-.28, 5), (-.15, 5), (-.11, 5.5), (-.09, 5.5)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-15, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-70)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected", fontdict={'fontsize': 16})

                    if plotnum == 22 and i == 0:

                        cf = plt.contourf(xevz, yevz, ratiopipeevz_blurr, [.74, .78, .82, .86, .90, .94, .98, 1],
                                          cmap='Reds')
                        cs = plt.contour(xevz, yevz, ratiopipeevz_blurr, [.78, .86, .94], colors='k',
                                         linewidths=1.2)
                        manual_locations = [(-.35, 6.2), (-.22, 7.3), (-.1, 8.3)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-10, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(0)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected", fontdict={'fontsize': 16})

                    if plotnum == 23 and i == 0:

                        cf = plt.contourf(xevz, yevz, rationozzleevz_blurr, [.80, .82, .84, .86, .88, .90, .92, .94,
                                                                             .96, .98, 1], cmap='Blues')
                        cs = plt.contour(xevz, yevz, rationozzleevz_blurr, [.88, .92, .96],
                                         colors='k', linewidths=1.2)
                        manual_locations = [(-.4, 4.5), (-.35, 5), (-.2, 5)]
                        labels = plt.clabel(cs, inline=1, inline_spacing=-15, fontsize=14, manual=manual_locations)
                        for l in labels:
                            l.set_rotation(-65)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.title("Percentage of Particles Detected", fontdict={'fontsize': 16})

                    if plotnum == 24 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevz, yevz, ratioevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                       .5, .55, .6, .65, .7, .75, .8, .85,
                                                                       .9, .95, 1], cmap='YlGnBu')
                        cs = plt.contour(xevz, yevz, ratioevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                      .5, .55, .6, .65, .7, .75, .8, .85,
                                                                      .9, .95, 1], colors='k', linewidths=.2)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles Detected', fontsize=18)

                    if plotnum == 25 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevz, yevz, ratioconeevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                           .5, .55, .6, .65, .7, .75, .8, .85,
                                                                           .9, .95, 1], cmap='Greens')
                        cs = plt.contour(xevz, yevz, ratioconeevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                          .5, .55, .6, .65, .7, .75, .8, .85,
                                                                          .9, .95, 1], colors='k', linewidths=.2)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles That Could Be Detected (Not Blocked by the Cone)',
                                     fontsize=18)

                    if plotnum == 26 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevz, yevz, ratiopipeevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                           .5, .55, .6, .65, .7, .75, .8, .85,
                                                                           .9, .95, 1], cmap='Reds')
                        cs = plt.contour(xevz, yevz, ratiopipeevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                          .5, .55, .6, .65, .7, .75, .8, .85,
                                                                          .9, .95, 1], colors='k', linewidths=.2)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles That Could Be Detected (Not Blocked by the Pipe)',
                                     fontsize=18)

                    if plotnum == 27 and i > 0:
                        plt.subplot(2, 2, i)
                        cf = plt.contourf(xevz, yevz, rationozzleevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                             .5, .55, .6, .65, .7, .75, .8, .85,
                                                                             .9, .95, 1], cmap='Blues')
                        cs = plt.contour(xevz, yevz, rationozzleevz_blurr, [.05, .1, .15, .2, .25, .3, .35, .4, .45,
                                                                            .5, .55, .6, .65, .7, .75, .8, .85,
                                                                            .9, .95, 1], colors='k', linewidths=.2)
                        cbar = fig.colorbar(cf)
                        plt.xlabel('z (m)')
                        plt.ylabel('Energy (MeV)')
                        plt.suptitle('Percentage of Particles That Could Be Detected (Not Blocked by the Nozzle)',
                                     fontsize=18)

                    if plotnum > 28:
                        # Make a new dataframe to contain the angle information for cuts
                        ratiot = pd.DataFrame()
                        ratiot['thap'] = df['Theta_Deg'][detarr[i] & df["AllPossible"]]
                        ratiot['thunb'] = df['Theta_Deg'][detarr[i] & df["Unblocked"]]

                        ratioe = pd.DataFrame()
                        ratioe['eap'] = df['Energy'][detarr[i] & df["AllPossible"]]
                        ratioe['eunb'] = df['Energy'][detarr[i] & df["Unblocked"]]

                        ratioz = pd.DataFrame()
                        ratioz['zap'] = df['zpos_final'][detarr[i] & df["AllPossible"]]
                        ratioz['zunb'] = df['zpos_final'][detarr[i] & df["Unblocked"]]

                        cuttheta = pd.cut(ratiot['thap'], binstheta)
                        cuttheta_unblocked = pd.cut(ratiot['thunb'], binstheta)

                        cute = pd.cut(ratioe['eap'], binse)
                        cute_unblocked = pd.cut(ratioe['eunb'], binse)

                        cutz = pd.cut(ratioz['zap'], binsz)
                        cutz_unblocked = pd.cut(ratioz['zunb'], binsz)

                        divt = ratiot.groupby(cuttheta_unblocked)['thunb'].agg('count') / \
                               ratiot.groupby(cuttheta)['thap'].agg('count')
                        divt = divt.tolist()

                        dive = ratioe.groupby(cute_unblocked)['eunb'].agg('count') / \
                               ratioe.groupby(cute)['eap'].agg('count')
                        dive = dive.tolist()

                        divz = ratioz.groupby(cutz_unblocked)['zunb'].agg('count') / \
                               ratioz.groupby(cutz)['zap'].agg('count')
                        divz = divz.tolist()

                        if plotnum == 28:
                            plt.plot(tbins2, divt, marker='o')
                            plt.legend(['Total', 'Detector 2', 'Detector 1', 'Detector 3', 'Detector 4'],
                                       loc='lower right', fontsize=16)
                            plt.xlabel('Lab Angle (Deg)')
                            plt.ylabel('Fraction of Particles Blocked')

                        if plotnum == 29:
                            plt.plot(ebins2, dive, marker='o')
                            plt.legend(['Total', 'Detector 2', 'Detector 1', 'Detector 3', 'Detector 4'],
                                       loc='lower left', fontsize=16)
                            plt.xlabel('Energy (MeV)')
                            plt.ylabel('Fraction of Particles Blocked')

                        if plotnum == 30:
                            plt.plot(zbins2, divz, marker='o')
                            plt.legend(['Total', 'Detector 2', 'Detector 1', 'Detector 3', 'Detector 4'],
                                       loc='lower left', fontsize=16)
                            plt.xlabel('z (m)')
                            plt.ylabel('Fraction of Particles Blocked')

                # Handle the legend here for each plot that needs it.
                if i == 1 and (plotnum == 2 or plotnum == 3 or plotnum == 5 or
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
    print("########   ########  ##                 ########  ##        ########  ########  ")
    print("##         ##    ##  ##                 ##    ##  ##        ##    ##     ##     ")
    print("##         ##    ##  ##                 ##    ##  ##        ##    ##     ##     ")
    print("########   ##    ##  ##                 ########  ##        ##    ##     ##     ")
    print("      ##   ##    ##  ##                 ##        ##        ##    ##     ##     ")
    print("      ##   ##    ##  ##                 ##        ##        ##    ##     ##     ")
    print("########   ########  #######            ##        ########  ########     ##    ")
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
                while filenum > len(list_pickles) or filenum == 0:
                    filenum = int(input("\nChoose a number from the list: "))
                    if filenum == 0 or len(list_pickles) < filenum:
                        print("ERROR: Number entered does not correspond to a DataFrame file...")
                    else:
                        filein = list_pickles[filenum-1]
        else:
            filein = latest_file

    print("\nThe file to be used is: " + filein)
    plot(filein)
