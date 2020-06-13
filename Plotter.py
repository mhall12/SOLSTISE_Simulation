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
import fnmatch

def plot(pklin):

    # In the future I might hide the contour plots if the user doesn't have an allE pickle.

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

    # Make the color maps for the histograms here. Matplotlib does not have a standard Red/White etc color map.
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

    # Get the individual colors here we'll use for the legends
    blk = newcmpBlack(1)
    grn = newcmpGreen(1)
    blu = newcmpBlue(1)
    red = newcmpRed(1)

    # Switch is used in a while loop to see if we should continue running the progran. once it is !=0 the program closes
    switch = 0

    # Histogram order to get the quadrants right.
    ho = [0, 2, 1, 3, 4]

    # Set the max and min values for the various histogram axis parameters
    if invkin:
        zmax = 0
        zmin = df['zpos_final'].min()
    else:
        zmin = 0
        zmax = df['zpos_final'].max()

    emax = 2 * df['Energy'].mean() + 2

    exmax = 2 * df['Ex_Reconstructed'].mean() + 0.5

    # emin should always just be 0

    thmax = df['Theta_Deg'].max()
    thmin = df['Theta_Deg'].min()

    print("\n\nChoose from the list below to plot histograms from the generated data.\n"
          "The four histograms represent the four quadrants of a fictional cylindrical detector\n"
          "looking down the beam axis.\n")
    if not fnmatch.fnmatch(pklin, '*eloss_s*'):
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
              "11) Counts vs Ex: Unblocked and blocked particles Ex from detected energy and position in all 4 "
              "detectors (1D)\n"
              "12) Counts vs Ex: Unblocked particles Ex from detected energy and position in all 4 detectors (1D)\n")
    if fnmatch.fnmatch(pklin, '*allE*'):
        print("13) Energy vs Lab Angle: Contour plot of detected particles. \n"
              "14) Energy vs Lab Angle: Contour plot of particles not blocked by the cone. \n"
              "15) Energy vs Lab Angle: Contour plot of particles not blocked by the pipe. \n"
              "16) Energy vs Lab Angle: Contour plot of particles not blocked by the nozzle. \n"
              "17) Energy vs Lab Angle: Contour plot of detected particles in all 4 detectors. \n"
              "18) Energy vs Lab Angle: Contour plot of particles not blocked by the cone in all 4 detectors. \n"
              "19) Energy vs Lab Angle: Contour plot of particles not blocked by the pipe in all 4 detectors. \n"
              "20) Energy vs Lab Angle: Contour plot of particles not blocked by the nozzle in all 4 detectors. \n"
              "21) Energy vs z: Contour plot of detected particles. \n"
              "22) Energy vs z: Contour plot of particles not blocked by the cone. \n"
              "23) Energy vs z: Contour plot of particles not blocked by the pipe. \n"
              "24) Energy vs z: Contour plot of particles not blocked by the nozzle. \n"
              "25) Energy vs z: Contour plot of detected particles in all 4 detectors.. \n"
              "26) Energy vs z: Contour plot of particles not blocked by the cone in all 4 detectors. \n"
              "27) Energy vs z: Contour plot of particles not blocked by the pipe in all 4 detectors. \n"
              "28) Energy vs z: Contour plot of particles not blocked by the nozzle in all 4 detectors. \n")
    if not fnmatch.fnmatch(pklin, '*eloss_s*'):
        print("29) Fraction of particles blocked vs lab angle. \n"
              "30) Fraction of particles blocked vs energy. \n"
              "31) Fraction of particles blocked vs z position. \n")
    if fnmatch.fnmatch(pklin, '*eloss_s*'):
        print("33) Energy vs z: Unblocked particles in all 4 detectors (2D). \n"
              "34) Energy vs z: Blocked particles in all 4 detectors (2D). \n"
              "35) Counts vs Ex: Unblocked particles Ex from detected energy and position (1D) \n"
              "36) Counts vs Ex: Unblocked and blocked particles Ex from detected energy and position in all 4 "
              "detectors (1D)\n"
              "37) Counts vs Ex: Unblocked particles Ex from detected energy and position in all 4 detectors (1D)\n")
    print("0) End\n\n")

    while switch == 0:
        plt.ion()

        # try except ensures the program does not segfault if the user does not enter a number (i.e. presses enter too
        # many times or something.)
        while True:
            try:
                plotnum = int(input("Enter a Histogram Number: "))
                break
            except:
                print("\n*****Enter an integer number from the list!*****\n")

        if plotnum > 0:
            # initialize the figure here, it might not be necessary.
            fig = plt.figure()

            # Set the font sizes here and the font used below. The title font sizes are handled when the hist is made.
            plt.rc('axes', labelsize=18)
            plt.rc('xtick', labelsize=18)
            plt.rc('ytick', labelsize=18)
            plt.rcParams["font.family"] = "STIXGeneral"

            for i in range(5):
                # The first histograms are Energy vs Z 2D histograms. 1) unblocked particles only, 2) blocked only,
                # 3) Unblocked and blocked together on the same hist.
                if (plotnum == 1 or plotnum == 3) and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist2d(df['zpos_final'][detarr[i] & df["Unblocked"]],
                               df['Energy'][detarr[i] & df["Unblocked"]], bins=(750, 750),
                               range=[[zmin, zmax], [0, emax]], cmap=newcmpBlack)
                    plt.xlabel('z(m)')
                    plt.ylabel('Energy (MeV)')

                    print(df['Phi'][detarr[i]] * 180 / np.pi)

                if (plotnum == 2 or plotnum == 3) and i > 0:

                    # This section giving a KeyError when the DataFrame size is 1, so I'll try except it.
                    plt.subplot(2, 2, i)
                    try:
                        plt.hist2d(df['zpos_final'][detarr[i] & df["Blocked_Cone"]],
                                   df['Energy'][detarr[i] & df["Blocked_Cone"]], bins=(750, 750),
                                   range=[[zmin, zmax], [0, emax]], cmap=newcmpGreen)
                        if df["Pipe Radius"][0] > 0.01:
                            plt.hist2d(df['zpos_final'][detarr[i] & df["Blocked_Pipe"]],
                                       df['Energy'][detarr[i] & df["Blocked_Pipe"]], bins=(750, 750),
                                       range=[[zmin, zmax], [0, emax]], cmap=newcmpRed)
                        plt.hist2d(df['zpos_final'][detarr[i] & df["Blocked_Nozzle"]],
                                   df['Energy'][detarr[i] & df["Blocked_Nozzle"]], bins=(750, 750),
                                   range=[[zmin, zmax], [0, emax]], cmap=newcmpBlue)
                        plt.xlabel('z(m)')
                        plt.ylabel('Energy (MeV)')
                    except KeyError:
                        print("ERROR")

                # This hist is just a 1D number of blocked counts vs theta. A better percentage version broken up by
                # detector is made in hist 28, this was just a first pass.
                if plotnum == 4 and i == 0:
                    plt.hist(df['Theta_Deg'][df['Blocked_Cone'] | df['Blocked_Pipe'] | df['Blocked_Nozzle']])
                    plt.xlabel('Lab Angle (deg)')
                    plt.ylabel('Counts')

                # 5 is raw counts vs Z in all four detectors, where 6 is the same but a stacked hist of blocked
                # particles, broken up by blocking type and detector. Different version
                # broken up by detector only in plot 30
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

                # 7 is the same as 6 but for theta instead of z.
                if plotnum == 7 and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist((df['Theta_Deg'][df['Blocked_Cone'] & detarr[i]],
                              df['Theta_Deg'][df['Blocked_Pipe'] & detarr[i]],
                              df['Theta_Deg'][df['Blocked_Nozzle'] & detarr[i]]),
                             bins=60, range=[90, 120], color=(grn, red, blu), stacked=True)
                    plt.xlabel('Lab Angle (Deg)')
                    plt.ylabel('Counts')

                # 8 is a 2D histogram of E vs theta for the unblocked particles only, 9 is the same but with blocked
                # particles as well, broken up by blocking type.
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

                # 10 is the Excitation Energy reconstructed from the "detected" energy and z position
                if plotnum == 10 and i == 0:
                    plt.hist(df['Ex_Reconstructed'][df["Unblocked"]], bins=750, range=[-.2, exmax])
                    plt.xlabel('Excitation Energy (MeV)')
                    plt.ylabel('Counts')

                # 11 is the same as 10 but also showing blocked particles.
                if plotnum == 11 and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist((df['Ex_Reconstructed'][df["Unblocked"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Cone"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Pipe"] & detarr[i]],
                              df['Ex_Reconstructed'][df["Blocked_Nozzle"] & detarr[i]]),
                             bins=1000, range=[0, exmax], color=(blk, grn, red, blu), stacked=True)
                    plt.xlabel('Excitation Energy (MeV)')
                    plt.ylabel('Counts')

                if plotnum == 12 and i > 0:
                    plt.subplot(2, 2, i)
                    plt.hist(df['Ex_Reconstructed'][df["Unblocked"] & detarr[i]], bins=750, range=[-.2, exmax])
                    plt.xlabel('Excitation Energy (MeV)')
                    plt.ylabel('Counts')

                # 14-29 are contour plots of blocked particles. These should only be used with "allE" simulated files
                # because they don't really make sense with specific excited states populated.
                if 12 < plotnum < 33:

                    # Make Energy vs Theta contour plot here. Theta goes from 90 to 180 and we'll use bins every 1
                    # degree.
                    binstheta = np.zeros(91)
                    for j in range(91):
                        if invkin:
                            binstheta[j] = j * 1 + 90
                        else:
                            binstheta[j] = j * 1
                    # We'll also make the energy bins as well:
                    numebins = 150
                    binse = np.zeros(numebins + 1)
                    for j in range(numebins + 1):
                        binse[j] = j * (emax / numebins)

                    numzbins = 100
                    binsz = np.zeros(numzbins + 1)
                    # We'll also make the z bins here:
                    for j in range(numzbins + 1):
                        if invkin:
                            binsz[j] = (100-j) * zmin / numzbins
                        else:
                            binsz[j] = j * zmax / numzbins

                    numphibins = 45
                    binsphi = np.zeros(numphibins + 1)
                    for j in range(numphibins + 1):
                        binsphi[j] = j * (360 / numphibins)

                    # Can't have detarr[i] alone. It has to include AllPossible because the Det masks do not contain
                    # the mask that detemines whether or not the particle hit the magnet bore. So, with detarr[0] you
                    # technically have AllPossible & AllPossible but it's fine because it's always true.

                    unblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][detarr[i] & df["AllPossible"]],
                                                                df['Energy'][detarr[i] & df["AllPossible"]],
                                                                bins=(binstheta, binse))
                    unblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][detarr[i] & df["AllPossible"]],
                                                                df['Energy'][detarr[i] & df["AllPossible"]],
                                                                bins=(binsz, binse))
                    df['Phi_Deg'] = df['Phi'] * 180 / np.pi

                    unblockedevphi, pbins, ebins = np.histogram2d(df['Phi_Deg'][detarr[i] & df["AllPossible"]],
                                                                  df['Energy'][detarr[i] & df["AllPossible"]],
                                                                  bins=(binsphi, binse))

                    if plotnum < 29:
                        # The following lines bin the particles into either energy vs theta or energy vs z 2d histograms
                        # Ex: unbloxkedevt is an array that has given each particle a theta and E bin, and tbins and
                        # ebins are the corresponding bin edges.

                        # a note, unblocked is actually "AllPossible" NOT "Unblocked" like the mask...

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

                        # Here we get the ratios that we want to use to make the contour plots. Instead of unblocked minus
                        # blocked over unblocked we could have used the unblocked mask and just done "detected/unblocked"
                        # where "detected" comes from the unblocked mask.

                        # The ratios are also split up into separate ratios for the particles only blocked by the cone, etc.

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

                        # To actually plot the ratios into histograms we have to transpose the binned arrays:

                        ratioevt = ratioevt.T
                        ratioevz = ratioevz.T

                        ratioconeevt = ratioconeevt.T
                        ratiopipeevt = ratiopipeevt.T
                        rationozzleevt = rationozzleevt.T

                        ratioconeevz = ratioconeevz.T
                        ratiopipeevz = ratiopipeevz.T
                        rationozzleevz = rationozzleevz.T

                        # As mentioned, tbins, ebins, and zbins are the bin edges. Here initialize a new array:

                        tbins2 = np.zeros(90)
                        ebins2 = np.zeros(150)
                        zbins2 = np.zeros(100)

                        # And here get the bin centers by taking the averages of three bin edges.

                        for k in range(150):
                            if k < 90:
                                tbins2[k] = (tbins[k] + tbins[k + 1]) / 2
                            if k < 100:
                                zbins2[k] = (zbins[k] + zbins[k + 1]) / 2
                            ebins2[k] = (ebins[k] + ebins[k + 1]) / 2

                        # To plot them we need to make a mesh grid of the bins:

                        xevt, yevt = np.meshgrid(tbins2, ebins2)
                        xevz, yevz = np.meshgrid(zbins2, ebins2)

                        # Now, since we have so many bins, the contour plots won't look nice. If we put too few bins the
                        # contours also don't look great. So, the solution is to use a lot of bins and use these gaussian
                        # filters on the ratios to smooth them out.

                        ratioevt_blurr = ndimage.gaussian_filter(ratioevt, sigma=1.5, order=0)
                        ratioevz_blurr = ndimage.gaussian_filter(ratioevz, sigma=1.5, order=0)

                        ratioconeevt_blurr = ndimage.gaussian_filter(ratioconeevt, sigma=1.5, order=0)
                        ratiopipeevt_blurr = ndimage.gaussian_filter(ratiopipeevt, sigma=1.5, order=0)
                        rationozzleevt_blurr = ndimage.gaussian_filter(rationozzleevt, sigma=1.5, order=0)

                        ratioconeevz_blurr = ndimage.gaussian_filter(ratioconeevz, sigma=1.6, order=0)
                        ratiopipeevz_blurr = ndimage.gaussian_filter(ratiopipeevz, sigma=1.6, order=0)
                        rationozzleevz_blurr = ndimage.gaussian_filter(rationozzleevz, sigma=1.6, order=0)

                        # 13 is the contour plot of percentage of detected particles.
                        if plotnum == 13 and i == 0:

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

                        # 14 is a contour plot of percetage of particles not blocked by the cone.
                        if plotnum == 14 and i == 0:

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
                            plt.title("Percentage of Particles Not Blocked by the Cone",
                                      fontdict={'fontsize': 16})

                        # 15 is the percentage of particles not blocked by the pipe
                        if plotnum == 15 and i == 0:

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
                            plt.title("Percentage of Particles Not Blocked by the Pipe",
                                      fontdict={'fontsize': 16})

                        # 16 is the percentage of particles not blocked by the nozzle
                        if plotnum == 16 and i == 0:

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

                        # 17 is the Energy vs angle split into the four detector quadrants, same as 12
                        if plotnum == 17 and i > 0:
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

                        # 18 is the Energy vs angle not blocked by the cone split into the four detector quadrants,
                        # same as 13
                        if plotnum == 18 and i > 0:
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
                            plt.suptitle('Percentage of Particles Not Blocked by the Cone',
                                         fontsize=18)

                        # 19 is the Energy vs angle not blocked by the pipe split into the four detector quadrants,
                        # same as 14
                        if plotnum == 19 and i > 0:
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
                            title18 = "Percentage of Particles Not Blocked by the Pipe B = " + \
                                      str(df['Magnetic Field'][0]) + " T"
                            plt.suptitle(title18, fontsize=18)

                        # 20 is the Energy vs angle not blocked by the nozzle split into the four detector quadrants,
                        # same as 16
                        if plotnum == 20 and i > 0:
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

                        # The contour plot cycle repeats here but is instead made with Energy vs z position.

                        if plotnum == 21 and i == 0:

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

                        if plotnum == 22 and i == 0:

                            cf = plt.contourf(xevz, yevz, ratioconeevz_blurr, [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75,
                                                                               .8, .85, .9, .95, 1], cmap='Greens')
                            cs = plt.contour(xevz, yevz, ratioconeevz_blurr, [.5, .6, .7, .8, .9], colors='k',
                                             linewidths=1.2)
                            manual_locations = [(-.28, 5), (-.15, 5), (-.11, 5.5), (-.09, 5.5)]
                            labels = plt.clabel(cs, inline=1, inline_spacing=-15, fontsize=16, manual=manual_locations)
                            for l in labels:
                                l.set_rotation(-70)
                            cbar = fig.colorbar(cf)
                            plt.xlabel('z (m)')
                            plt.ylabel('Energy (MeV)')
                            title21 = df['Reaction'][0] + " Percentage of Particles Not Blocked by the Cone, B = " + \
                                      str(df['Magnetic Field'][0]) + " T"
                            plt.suptitle(title21, fontsize=18)

                        if plotnum == 23 and i == 0:

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

                        if plotnum == 24 and i == 0:

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

                        if plotnum == 25 and i > 0:
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

                        if plotnum == 26 and i > 0:
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
                            title25 = df['Reaction'][0] + " Percentage of Particles Not Blocked by the Cone, B = " + \
                                      str(df['Magnetic Field'][0]) + " T"
                            plt.suptitle(title25, fontsize=18)

                        if plotnum == 27 and i > 0:
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
                            title26 = df['Reaction'][0] + " Percentage of Particles Not Blocked by the Pipe B = " + \
                                      str(df['Magnetic Field'][0]) + " T"
                            plt.suptitle(title26, fontsize=18)

                        if plotnum == 28 and i > 0:
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
                            title27 = df['Reaction'][0] + " Percentage of Particles Not Blocked by the Nozzle B = " + \
                                      str(df['Magnetic Field'][0]) + " T"
                            plt.suptitle(title27, fontsize=18)

                        # Here we do something a little different, and instead of contour plots we do ratio plots,
                        # still using the bins that were defined above.

                if plotnum > 28 and plotnum < 33:
                    # Unfortunately I couldn't get the cuts to work by putting in the masks, so we have to create
                    # some dummy dataframes to contain them so we can cut them later. As mentioned before, here
                    # I actually use Unblocked and AllPossible since I'm not trying to break them up by cone, nozzle
                    # etc.

                    # Make a new dataframe to contain the angle information for cuts
                    ratiot = pd.DataFrame()
                    ratiot['thap'] = df['Theta_Deg'][detarr[ho[i]] & df["AllPossible"]]
                    ratiot['thunb'] = df['Theta_Deg'][detarr[ho[i]] & df["Unblocked"]]
                    # Make a new dataframe to contain the E information for cuts
                    ratioe = pd.DataFrame()
                    ratioe['eap'] = df['Energy'][detarr[ho[i]] & df["AllPossible"]]
                    ratioe['eunb'] = df['Energy'][detarr[ho[i]] & df["Unblocked"]]
                    # Make a new dataframe to contain the z information for cuts
                    ratioz = pd.DataFrame()
                    ratioz['zap'] = df['zpos_final'][detarr[ho[i]] & df["AllPossible"]]
                    ratioz['zunb'] = df['zpos_final'][detarr[ho[i]] & df["Unblocked"]]
                    # Make a new dataframe to contain the phi information for cuts
                    ratiop = pd.DataFrame()
                    ratiop['pap'] = df['Phi_Deg'][detarr[ho[i]] & df["AllPossible"]]
                    ratiop['punb'] = df['Phi_Deg'][detarr[ho[i]] & df["Unblocked"]]

                    # Now we can cut the data based on the bins that we defined before. This effectively makes
                    # a pandas series that shows the bin each particle falls into
                    cuttheta = pd.cut(ratiot['thap'], binstheta)
                    cuttheta_unblocked = pd.cut(ratiot['thunb'], binstheta)

                    cute = pd.cut(ratioe['eap'], binse)
                    cute_unblocked = pd.cut(ratioe['eunb'], binse)

                    cutz = pd.cut(ratioz['zap'], binsz)
                    cutz_unblocked = pd.cut(ratioz['zunb'], binsz)

                    cutp = pd.cut(ratiop['pap'], binsphi)
                    cutp_unblocked = pd.cut(ratiop['punb'], binsphi)

                    # Now we can make the ratios here. We have to group the cuts by the bins and aggregate the data
                    # End up making a pandas series of bin, counts.
                    divt = ratiot.groupby(cuttheta_unblocked)['thunb'].agg('count') / \
                           ratiot.groupby(cuttheta)['thap'].agg('count')
                    divt = divt.tolist()

                    dive = ratioe.groupby(cute_unblocked)['eunb'].agg('count') / \
                           ratioe.groupby(cute)['eap'].agg('count')
                    dive = dive.tolist()

                    divz = ratioz.groupby(cutz_unblocked)['zunb'].agg('count') / \
                           ratioz.groupby(cutz)['zap'].agg('count')
                    divz = divz.tolist()

                    unbcount = ratiop.groupby(cutp_unblocked)['punb'].agg('count')

                    totcount = ratiop.groupby(cutp)['pap'].agg('count')

                    zeromask = totcount > (np.amax(unbcount) / 25)

                    unbcount = unbcount[zeromask]
                    totcount = totcount[zeromask]

                    divp = unbcount / totcount

                    # As mentioned, tbins, ebins, and zbins are the bin edges. Here initialize a new array:

                    tbins2 = np.zeros(90)
                    ebins2 = np.zeros(150)
                    zbins2 = np.zeros(100)
                    pbins2 = np.zeros(numphibins)

                    # And here get the bin centers by taking the averages of the bin edges.

                    for k in range(150):
                        if k < 90:
                            tbins2[k] = (tbins[k] + tbins[k + 1]) / 2
                        if k < numphibins:
                            pbins2[k] = (pbins[k] + pbins[k + 1]) / 2
                        if k < 100:
                            zbins2[k] = (zbins[k] + zbins[k + 1]) / 2
                        ebins2[k] = (ebins[k] + ebins[k + 1]) / 2

                    pbins2 = np.array(pbins2)
                    pbins2 = pbins2[zeromask]

                    # The next 3 hists are the ratio plots broken up by detector on the same plot:
                    if plotnum == 29:
                        plt.plot(tbins2, divt, marker='o')
                        plt.legend(['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4'],
                                   loc='lower right', fontsize=16)
                        plt.xlabel('Lab Angle (Deg)')
                        plt.ylabel('Fraction of Particles Detected')

                    if plotnum == 30:
                        plt.plot(ebins2, dive, marker='o')
                        plt.legend(['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4'],
                                   loc='lower left', fontsize=16)
                        plt.xlabel('Energy (MeV)')
                        plt.ylabel('Fraction of Particles Detected')

                    if plotnum == 31:
                        plt.plot(zbins2, divz, marker='o')
                        plt.legend(['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4'],
                                   loc='lower left', fontsize=16)
                        plt.xlabel('z (m)')
                        plt.ylabel('Fraction of Particles Detected')

                    if plotnum == 32:

                        stdcolors = ['Blue', 'Orange', 'limegreen', 'Red', 'Purple']
                        labels32 = ['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4']
                        if i == 0:
                            plt.rc('axes', labelsize=16)
                            plt.rc('xtick', labelsize=16)
                            plt.rc('ytick', labelsize=16)
                            ax32 = plt.subplot(111, projection='polar')
                        pbins2 = pbins2 * np.pi / 180
                        if i != 3:
                            if i > 0:
                                ax32.plot(pbins2, divp, marker='o', color=stdcolors[i], label=labels32[i], markersize=8, MarkerEdgeColor='Black', alpha=0.7)
                            if i == 0:
                                ax32.plot(pbins2, divp, marker='o', color=stdcolors[i], label=labels32[i], markersize=14, MarkerEdgeColor='Black', alpha=0.7)
                        if i == 3:
                            pbins2gt0 = pbins2[pbins2 < np.pi]
                            divpgt0 = divp[pbins2 < np.pi]

                            pbins2lt0 = pbins2[pbins2 > np.pi]
                            divplt0 = divp[pbins2 > np.pi]

                            ax32.plot(pbins2gt0, divpgt0, marker='o', color=stdcolors[i], label=labels32[i], markersize=8, MarkerEdgeColor='Black', alpha=0.7)
                            ax32.plot(pbins2lt0, divplt0, marker='o', color=stdcolors[i], markersize=8, MarkerEdgeColor='Black', alpha=0.7)

                        #ax32.set_rorigin(.2)
                        ax32.grid(True)
                        plt.legend(fontsize=16, bbox_to_anchor=(1.0, .9))
                        plt.xlabel('Initial Phi Angle (Deg)')
                        plt.ylabel('Fraction of Particles Detected', rotation=0, size=14, labelpad=-370)

                # This next section is for solid targets only:
                if plotnum > 32:

                    if (plotnum == 33 or plotnum == 34) and i > 0:
                        plt.subplot(2, 2, i)
                        plt.hist2d(df['zpos_final'][detarr[i] & df["UnblockedSolidTarg"]],
                                    df['Energy'][detarr[i] & df["UnblockedSolidTarg"]], bins=(750, 750),
                                    range=[[zmin, zmax], [0, emax]], cmap=newcmpBlack)
                        plt.xlabel('z(m)')
                        plt.ylabel('Energy (MeV)')

                        if plotnum == 34:
                            plt.hist2d(df['zpos_final'][detarr[i] & ~df["UnblockedSolidTarg"]],
                                       df['Energy'][detarr[i] & ~df["UnblockedSolidTarg"]], bins=(750, 750),
                                       range=[[zmin, zmax], [0, emax]], cmap=newcmpRed)
                            plt.xlabel('z(m)')
                            plt.ylabel('Energy (MeV)')

                    if (plotnum == 35) and i > 0:
                        plt.subplot(2, 2, i)
                        plt.hist2d(df['zpos_final'][detarr[i] & ~df["UnblockedSolidTarg"]],
                                    df['Energy'][detarr[i] & ~df["UnblockedSolidTarg"]], bins=(750, 750),
                                    range=[[zmin, zmax], [0, emax]], cmap=newcmpBlack)
                        plt.xlabel('z(m)')
                        plt.ylabel('Energy (MeV)')

                    if plotnum == 36 and i == 0:
                        plt.hist(df['Ex_Reconstructed'][df["UnblockedSolidTarg"]], bins=750, range=[-0.2, exmax])
                        plt.xlabel('Excitation Energy (MeV)')
                        plt.ylabel('Counts')

                        # 35 is the same as 11 but also showing blocked particles.
                    if plotnum == 37 and i > 0:
                        plt.subplot(2, 2, i)
                        plt.hist((df['Ex_Reconstructed'][df["UnblockedSolidTarg"] & detarr[i]],
                                  df['Ex_Reconstructed'][~df["UnblockedSolidTarg"] & detarr[i]]),
                                 bins=1000, range=[-0.2, exmax], color=(blk, red), stacked=True)
                        plt.xlabel('Excitation Energy (MeV)')
                        plt.ylabel('Counts')

                    if plotnum == 38 and i > 0:
                        plt.subplot(2, 2, i)
                        plt.hist(df['Ex_Reconstructed'][df["UnblockedSolidTarg"] & detarr[i]], bins=750,
                                 range=[-0.2, exmax])
                        plt.xlabel('Excitation Energy (MeV)')
                        plt.ylabel('Counts')

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

    # Get a list of all the pkl files that are already created:
    list_pickles = glob.glob('*evt*.pkl')

    # If there aren't any pickles, prompt the user to make one.
    if len(list_pickles) == 0:
        print("\nIt appears that no simulated DataFrame exists, run SOLSTISE_Sim.py first.")
    else:
        # If there are pickles, get the pickle that was modified last:
        latest_file = max(list_pickles, key=os.path.getctime)
        print("\nThe most recently created simulated DataFrame file is: " + latest_file)
        yn = input("\nWould you like to use this file? [Y/N] ")
        if yn == "N" or yn == "n":
            if len(list_pickles) == 1:
                # If the user wants to use a different pickle file but there is only one, make the user rum the sim.
                print("\nYou only have one DataFrame file. To create another, run SOLSTISE_Sim.py")
            # if there is more than one, list them and give them a number so the user can choose.
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
