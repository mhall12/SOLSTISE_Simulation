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
import sys


class Color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


def plot(pklin, pkloverlay):

    # Read the pickled DataFrame File here.
    df = pd.read_pickle(pklin)

    df_main = df

    detarr = [df["AllPossible"], df["Det2"], df["Det1"], df["Det3"], df["Det4"]]

    if fnmatch.fnmatch(pkloverlay, "*evts*"):
        df_over = pd.read_pickle(pkloverlay)

    print(pklin)

    print("\n********************************Simulation Parameters*********************************\n"
          "The Reacton is: " + df['Reaction'][0] +
          "\nBeam Energy: " + str(df['Beam Energy'][0]) + " MeV" +
          "\nMagnetic Field: " + str(df['Magnetic Field'][0]) + " Tesla" +
          "\nReaction Distance from Nozzle: " + str(df['Reaction Distance from Nozzle'][0]) + ' inches' +
          "\nNozzle-Cone Distance: " + str(df['Nozzle-Cone Distance'][0]) + ' inches'
          "\nReceiver Cone Opening Diameter: " + str(df['Cone Opening Diameter'][0]) + " inches" +
          "\nReceiver Cone Height: " + str(df["Cone Height"][0]) + " inches"
          "\nMagnet Bore Radius: " + str(df['Bore Radius'][0]) + ' m' +
          "\nCustom Pipe? " + str(df['Custom Pipe?'][0]) +
          "\nPipe Radius: " + str(df["Pipe Radius"][0]) + " m"
          "\nPipe Left Edge Angle: " + str(df['Pipe Left Edge Angle'][0]) + " Degrees" +
          "\nPipe Right Edge Angle: " + str(df['Pipe Right Edge Angle'][0]) + " Degrees" +
          "\nReceiver Cone File Used: " + df['Cone File'][0] +
          "\nNozzle File Used: " + df['Nozzle File'][0] +
          "\nSolid Target Thickness: " + str(df['Solid Thickness'][0]) + " mg/cm^2" +
          "\nJet Pressure: " + str(df['Jet Pressure'][0]) + " Torr" +
          "\nJet Radius: " + str(df["Jet Radius"][0]) + " mm" +
          "\nChamber Pressure: " + str(df['Chamber Pressure'][0]) + " Torr" +
          "\n**************************************************************************************\n")

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

    # Switch is used in a while loop to see if we should continue running the program. once it is !=0 the program closes
    switch = 0

    # Histogram order to get the quadrants right.
    ho = [0, 2, 1, 3, 4]

    # Set the max and min values for the various histogram axis parameters
    if invkin:
        zmax = 0
        zmin = df['zpos_final'].min()
        thmax =180
        thmin = 90
    else:
        zmin = 0
        zmax = df['zpos_final'].max()
        thmax = 90
        thmin = 0

    # Get the ax energy from the data frame. This assumes no shadowing, so we can also see where the bore shadows.
    emax = df['Max Energy'][0]

    exmax = 2 * df['Ex_Reconstructed'].mean() + 0.5

    # emin should always just be 0

    thmax = df['Theta_Deg'].max()
    thmin = df['Theta_Deg'].min()

    # Change quadrant of CM angle to be in line with HELIOS spreadsheet
    df['CM_Deg'] = -1 * df["Theta_CM"] * 180 / np.pi + 180
    cmmin = df['CM_Deg'].min()
    cmmax = df['CM_Deg'].max()

    if df['CM_Deg'].mean() < 90:
        cmflag = True
    else:
        cmflag = False

    if fnmatch.fnmatch(pkloverlay, "*evts*"):
        df_over['CM_Deg'] = -1 * df_over["Theta_CM"] * 180 / np.pi + 180

    overlaybool = False
    detposbool = False
    overlayonoff = Color.RED + "OFF" + Color.END
    detposonoff = Color.RED + "OFF" + Color.END

    lastentry = ''

    while switch == 0:

        print('\n' + Color.BOLD + "Standard Plots:\n" + Color.END + Color.CYAN + Color.UNDERLINE +
              "A:" + Color.END + Color.CYAN + " 1) Energy vs z-position\n"
                                              "   2) Lab Angle vs Energy\n"
                                              "   3) Lab Angle vs z-position\n"
                                              "   4) Center-of-Mass Angle vs z-position\n"
                                              "   5) Excitation Energy Spectrum (Set C=0)" + Color.END)
        print(Color.RED + Color.UNDERLINE +
              "\tB:" + Color.END + Color.RED + " 0) Sum Spectrum\n"
                                               "\t   1) Broken down by detector" + Color.END)
        print(Color.GREEN + Color.UNDERLINE + "\t\tC:" + Color.END + Color.GREEN + " 0) 2D Histogram(s)" + Color.END)
        if not fnmatch.fnmatch(pklin, '*allE*') or overlaybool:
            print(
                Color.UNDERLINE + Color.YELLOW + "\t\t\tD:" + Color.END + Color.YELLOW + " 0) Unblocked Particles "
                                                                                         "Only\n"
                                                                                         "\t\t\t   1) Unblocked + "
                                                                                         "Shadowed Particles\n"
                                                                                         "\t\t\t   2) Shadowed "
                                                                                         "Particles Only" + Color.END)
        if fnmatch.fnmatch(pklin, '*allE*') and not overlaybool:
            print(Color.GREEN + "\t\t   1) Shadowing Contour Plots (no D needed)\n"
                                "\t\t   2) Cone Shadowing Contour Plots (no D needed)\n"
                                "\t\t   3) Pipe Shadowing Contour Plots (no D needed)\n"
                                "\t\t   4) Nozzle Shadowing Contour Plots (no D needed)" + Color.END)
        print(Color.BOLD + "\nSpecial Plots:\n" + Color.END + Color.CYAN + Color.UNDERLINE +
              "A:" + Color.END + Color.CYAN + " 6) Phi vs Lab Angle Polar Plots (no C/D needed)" + Color.END)
        print(Color.RED + Color.UNDERLINE + "\tB:" + Color.END + Color.RED + " 0) Ratio Plot")
        print("\t   1) Contour Plot" + Color.END)
        print(Color.CYAN + "   7) 1D Ratio Plots (no C/D needed)" + Color.END)
        print(Color.RED + Color.UNDERLINE + "\tB:" + Color.END + Color.RED + " 0) Fraction Blocked vs Lab Angle\n" +
              "\t   1) Fraction Blocked vs Energy\n"
              "\t   2) Fraction Blocked vs z-position\n" + Color.END)

        print(Color.BOLD + "Other Options:" + Color.END)
        print(Color.PURPLE + "100: Toggle ON z-axis detector positions.")
        print("101: Toggle OFF z-axis detector positions.")
        if fnmatch.fnmatch(pkloverlay, "*evts*"):
            print("\n200: Toggle ON overlay histograms.")
            print("201: Toggle OFF overlay histograms.")
            print("300: Auto Overlay.")
        print("400: Print Data to CSV.")
        print("0) End\n\n" + Color.END)

        print("\nThe syntax for choosing plots is as follows: " + Color.CYAN + "A" + Color.END +
              "." + Color.RED + "B" + Color.END + "." + Color.GREEN + "C" + Color.END + "." + Color.YELLOW + "D" +
              Color.END + " (Example, 2.1.0.3).\n")

        plt.ion()

        if overlaybool:
            overlayonoff = Color.GREEN + "ON" + Color.END
        else:
            overlayonoff = Color.RED + "OFF" + Color.END

        # try except ensures the program does not segfault if the user does not enter a number (i.e. presses enter too
        # many times or something.)
        while True:
            try:
                print("\nDetector Positions: ", detposonoff)
                if fnmatch.fnmatch(pkloverlay, "*evts*"):
                    print("Overlay: ", overlayonoff)
                plotnum = str(input("Entry: "))
                if 100 > int(plotnum[0]) > 7:
                    raise ValueError
                break
            except (ValueError, IndexError):
                print("\n*****Incorrect Syntax!*****\n")

        currentry = plotnum

        #if currentry == '300':
        #    autooverlaybool = True
        #else:
        #    autooverlaybool = False

        if plotnum == '300':
            if lastentry == '':
                print("ERROR: You must make a contour plot before using the Auto Overlay...")
                lastentry = '1000'
            plotnum = lastentry
            plotnum = plotnum[:4] + '0' + plotnum[5:]
            plotnum = plotnum[:6] + '0'

        # Grab A, B, C, D here. If the user didn't specify a C/D, set them to 0.

        if plotnum.count(".") == 3:
            aa = plotnum[0]
            bb = plotnum[2]
            cc = plotnum[4]
            dd = plotnum[6]
        elif plotnum.count(".") == 2:
            aa = plotnum[0]
            bb = plotnum[2]
            cc = plotnum[4]
            dd = "9"
            if cc == '0':
                dd = "0"
        elif plotnum.count(".") == 1:
            aa = plotnum[0]
            bb = plotnum[2]
            cc = "9"
            dd = "9"
        else:
            aa = plotnum
            bb = "9"
            cc = "9"
            dd = "9"

        plotnum = aa + '.' + bb + '.' + cc + '.' + dd

        sol = fnmatch.fnmatch(pklin, '*_s*')

        if aa != "0" and int(aa) < 400:
            # Set the font sizes here and the font used below. The title font sizes are handled when the hist is made.
            plt.rc('axes', labelsize=18)
            plt.rc('xtick', labelsize=18)
            plt.rc('ytick', labelsize=18)
            plt.rcParams["font.family"] = "STIXGeneral"

            # 1) On when aa=200 or the current entry is 300
            # 2) Off when aa=201 or the current entry isn't 30

            if currentry == '300':
                overlaybool = True
                df = df_over
                if detposbool:
                    detarr = [df["AllPossible"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det2"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det1"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det3"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det4"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6'])]
                else:
                    detarr = [df["AllPossible"], df["Det2"], df["Det1"], df["Det3"], df["Det4"]]

            if aa == "200" and fnmatch.fnmatch(pkloverlay, "*evts*"):
                overlaybool = True
                df = df_over
                if detposbool:
                    detarr = [df["AllPossible"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det2"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det1"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det3"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det4"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6'])]
                else:
                    detarr = [df["AllPossible"], df["Det2"], df["Det1"], df["Det3"], df["Det4"]]
                overlayonoff = Color.GREEN + "ON" + Color.END

            if aa == "201" and fnmatch.fnmatch(pkloverlay, "*evts*"):
                overlaybool = False
                df = df_main
                if detposbool:
                    detarr = [df["AllPossible"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det2"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det1"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det3"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det4"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6'])]
                else:
                    detarr = [df["AllPossible"], df["Det2"], df["Det1"], df["Det3"], df["Det4"]]
                overlayonoff = Color.RED + "OFF" + Color.END

            # initialize the figure here, it might not be necessary.
            # All of these are 1D histograms or single 2D histograms
            #if (plotnum == 4 or plotnum == 8 or 10 <= plotnum <= 11 or plotnum == 14 or 17 <= plotnum <= 20 or
            #    25 <= plotnum <= 29 or 33 <= plotnum <= 36 or
            #    42 <= plotnum <= 44 or plotnum == 41) and not overlaybool:

            if ((int(aa) < 6 and bb == "0") or 100 > int(aa) > 5) and not overlaybool:
                fig, axs = plt.subplots()
            # All these are 1D or 2D histograms broken down by detector
            if (int(aa) < 6 and bb == "1") and not overlaybool:
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
                axi = np.array([0, ax1, ax2, ax3, ax4])

            if aa == "100":
                for j in range(5):
                    detarr = [df["AllPossible"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det2"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det1"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det3"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6']),
                              df["Det4"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                              df['Detz5'] | df['Detz6'])]
                detposonoff = Color.GREEN + "ON" + Color.END
                detposbool = True

            if aa == "101":
                for j in range(5):
                    detarr = [df["AllPossible"],
                              df["Det2"],
                              df["Det1"],
                              df["Det3"],
                              df["Det4"]]
                detposonoff = Color.RED + "OFF" + Color.END
                detposbool = False

            # *****************************************************************************************************

            for i in range(5):

                if not sol:
                    # These set the E vs z histograms
                    if int(aa) == 1 and int(cc) == 0:
                        xmin, xmax, ymin, ymax, xlabel, ylabel = zmin, zmax, 0, emax, 'z (m)', 'Energy (MeV)'
                        binsmin, binsmax = 550, 550

                        if 0 <= int(dd) <= 1:
                            dfx1, dfy1 = df['zpos_final'][detarr[i] & df["Unblocked"]], \
                                         df['Energy'][detarr[i] & df["Unblocked"]]

                        if 1 <= int(dd) <= 2:
                            dfx2, dfy2 = df['zpos_final'][detarr[i] & df["Blocked_Cone"]], \
                                         df['Energy'][detarr[i] & df["Blocked_Cone"]]
                            dfx3, dfy3 = df['zpos_final'][detarr[i] & df["Blocked_Pipe"]], \
                                         df['Energy'][detarr[i] & df["Blocked_Pipe"]]
                            dfx4, dfy4 = df['zpos_final'][detarr[i] & df["Blocked_Nozzle"]], \
                                         df['Energy'][detarr[i] & df["Blocked_Nozzle"]]

                    # These set the E vs Theta
                    if int(aa) == 2 and int(cc) == 0:
                        xmin, xmax, ymin, ymax, xlabel, ylabel = thmin, thmax, 0, emax, 'Theta (Deg)', \
                                                                 'Energy (MeV)'
                        binsmin, binsmax = 750, 750

                        if 0 <= int(dd) <= 1:
                            dfx1, dfy1 = df['Theta_Deg'][df["Unblocked"] & detarr[i]], \
                                         df['Energy'][detarr[i] & df["Unblocked"]]

                        if 1 <= int(dd) <= 2:
                            dfx2, dfy2 = df['Theta_Deg'][df["Blocked_Cone"] & detarr[i]], \
                                         df['Energy'][detarr[i] & df["Blocked_Cone"]]
                            dfx3, dfy3 = df['Theta_Deg'][df["Blocked_Pipe"] & detarr[i]], \
                                         df['Energy'][detarr[i] & df["Blocked_Pipe"]]
                            dfx4, dfy4 = df['Theta_Deg'][df["Blocked_Nozzle"] & detarr[i]], \
                                         df['Energy'][detarr[i] & df["Blocked_Nozzle"]]

                    # These set the Theta vs z hists
                    if int(aa) == 3 and int(cc) == 0:
                        xmin, xmax, ymin, ymax, xlabel, ylabel = zmin, zmax, thmin, thmax, 'z (m)', 'Theta (Deg)'
                        binsmin, binsmax = 550, 550

                        if 0 <= int(dd) <= 1:
                            dfx1, dfy1 = df['zpos_final'][df["Unblocked"] & detarr[i]], \
                                         df['Theta_Deg'][df["Unblocked"] & detarr[i]]

                        if 1 <= int(dd) <= 2:
                            dfx2, dfy2 = df['zpos_final'][df["Blocked_Cone"] & detarr[i]], \
                                         df['Theta_Deg'][df["Blocked_Cone"] & detarr[i]]
                            dfx3, dfy3 = df['zpos_final'][df["Blocked_Pipe"] & detarr[i]], \
                                         df['Theta_Deg'][df["Blocked_Pipe"] & detarr[i]]
                            dfx4, dfy4 = df['zpos_final'][df["Blocked_Nozzle"] & detarr[i]], \
                                         df['Theta_Deg'][df["Blocked_Nozzle"] & detarr[i]]

                    # These set the CM angle vs z hists
                    if int(aa) == 4 and int(cc) == 0:
                        xmin, xmax, ymin, ymax, xlabel, ylabel = zmin, zmax, cmmin, cmmax, 'z (m)', 'CM Angle (Deg)'
                        binsmin, binsmax = 550, 550

                        if 0 <= int(dd) <= 1:
                            dfx1, dfy1 = df['zpos_final'][df["Unblocked"] & detarr[i]], \
                                         df['CM_Deg'][df["Unblocked"] & detarr[i]]

                        if 1 <= int(dd) <= 2:
                            dfx2, dfy2 = df['zpos_final'][df["Blocked_Cone"] & detarr[i]], \
                                         df['CM_Deg'][df["Blocked_Cone"] & detarr[i]]
                            dfx3, dfy3 = df['zpos_final'][df["Blocked_Pipe"] & detarr[i]], \
                                         df['CM_Deg'][df["Blocked_Pipe"] & detarr[i]]
                            dfx4, dfy4 = df['zpos_final'][df["Blocked_Nozzle"] & detarr[i]], \
                                         df['CM_Deg'][df["Blocked_Nozzle"] & detarr[i]]

                # Determine here if the axes are single or multiple, and what the bool itest needs to be.
                if int(bb) == 1:
                    itest = i > 0
                    axes = axi[i]
                if int(bb) == 0:
                    itest = i == 0
                    axes = axs

                if int(aa) < 5 and int(cc) == 0 and 0 <= int(dd) <= 1 and itest:

                    axes.hist2d(dfx1, dfy1, bins=(binsmin, binsmax), range=[[xmin, xmax], [ymin, ymax]],
                                  cmap=newcmpBlack, zorder=1)
                    axes.set_xlabel(xlabel)
                    axes.set_ylabel(ylabel)

                if int(aa) < 5 and int(cc) == 0 and 1 <= int(dd) <= 2 and itest:

                    # This section giving a KeyError when the DataFrame size is 1, so I'll try except it.
                    try:
                        axes.hist2d(dfx2, dfy2, bins=(binsmin, binsmax), range=[[xmin, xmax], [ymin, ymax]],
                                    cmap=newcmpGreen)
                        axes.hist2d(dfx3, dfy3, bins=(binsmin, binsmax), range=[[xmin, xmax], [ymin, ymax]],
                                    cmap=newcmpRed)
                        axes.hist2d(dfx4, dfy4, bins=(binsmin, binsmax), range=[[xmin, xmax], [ymin, ymax]],
                                    cmap=newcmpBlue)
                        axes.set_xlabel(xlabel)
                        axes.set_ylabel(ylabel)
                    except KeyError:
                        print("ERROR")

                if int(aa) == 5 and int(cc) == 0:
                    dfx1 = df['Ex_Reconstructed'][df["Unblocked"] & detarr[i]]
                    dfx2 = df['Ex_Reconstructed'][df["Blocked_Cone"] & detarr[i]]
                    dfx3 = df['Ex_Reconstructed'][df["Blocked_Pipe"] & detarr[i]]
                    dfx4 = df['Ex_Reconstructed'][df["Blocked_Nozzle"] & detarr[i]]

                if int(aa) == 5 and itest:
                    axes.set_xlabel('Excitation Energy (MeV)')
                    axes.set_ylabel('Counts')
                    if int(dd) == 0:
                        axes.hist(dfx1, bins=750, range=[-.2, exmax])
                    if int(dd) == 1:
                        axes.hist((dfx1, dfx2, dfx3, dfx4), bins=750, range=[-.2, exmax],  color=(blk, grn, red, blu),
                                  stacked=True)
                    if int(dd) == 2:
                        axes.hist((dfx2, dfx3, dfx4), bins=750, range=[-.2, exmax],  color=(grn, red, blu),
                                  stacked=True)

            # *****************************************************************************************************

                # The next section are contour plots of blocked particles. These should only be used with "allE"
                # simulated files because they don't really make sense with specific excited states populated.
                if int(cc) > 0 and int(aa) < 7:

                    # Make Energy vs Theta contour plot here. Theta goes from 90 to 180 and we'll use bins every 1
                    # degree.
                    binstheta = np.zeros(91)
                    for j in range(91):
                        if invkin:
                            binstheta[j] = j * 1 + 90
                        else:
                            binstheta[j] = j * 1

                    binscm = np.zeros(91)
                    for j in range(91):
                        if cmflag:
                            binscm[j] = j * 1
                        else:
                            binscm[j] = j * 1 + 90

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
                    unblockedcvz, zbins, cbins = np.histogram2d(df['zpos_final'][detarr[i] & df["AllPossible"]],
                                                                df['CM_Deg'][detarr[i] & df["AllPossible"]],
                                                                bins=(binsz, binscm))
                    unblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][detarr[i] & df["AllPossible"]],
                                                                df['Energy'][detarr[i] & df["AllPossible"]],
                                                                bins=(binsz, binse))
                    unblockedtvz, zbins, tbins = np.histogram2d(df['zpos_final'][detarr[i] & df["AllPossible"]],
                                                                df['Theta_Deg'][detarr[i] & df["AllPossible"]],
                                                                bins=(binsz, binstheta))
                    df['Phi_Deg'] = df['Phi'] * 180 / np.pi

                    unblockedevphi, pbins, ebins = np.histogram2d(df['Phi_Deg'][detarr[i] & df["AllPossible"]],
                                                                  df['Energy'][detarr[i] & df["AllPossible"]],
                                                                  bins=(binsphi, binse))

                    unblockedtvphi, pbins, tbins = np.histogram2d(df['Phi_Deg'][detarr[i] & df["AllPossible"]],
                                                                  df['Theta_Deg'][detarr[i] & df["AllPossible"]],
                                                                  bins=(binsphi, binstheta))

                    if int(dd) > 1:
                        # The following lines bin the particles into either energy vs theta or energy vs z 2d histograms
                        # Ex: unbloxkedevt is an array that has given each particle a theta and E bin, and tbins and
                        # ebins are the corresponding bin edges.

                        # a note, unblocked is actually "AllPossible" NOT "Unblocked" like the mask...

                        try:
                            blockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][(df['Blocked_Cone'] |
                                                                                      df['Blocked_Pipe'] |
                                                                                      df['Blocked_Nozzle']) &
                                                                                      detarr[i]],
                                                                      df['Energy'][(df['Blocked_Cone'] |
                                                                                   df['Blocked_Pipe'] |
                                                                                   df['Blocked_Nozzle']) &
                                                                                   detarr[i]],
                                                                      bins=(binstheta, binse))

                            blockedcvz, zbins, cbins = np.histogram2d(df['zpos_final'][(df['Blocked_Cone'] |
                                                                                      df['Blocked_Pipe'] |
                                                                                      df['Blocked_Nozzle']) &
                                                                                      detarr[i]],
                                                                      df['CM_Deg'][(df['Blocked_Cone'] |
                                                                                   df['Blocked_Pipe'] |
                                                                                   df['Blocked_Nozzle']) &
                                                                                   detarr[i]],
                                                                      bins=(binsz, binscm))

                            blockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][(df['Blocked_Cone'] |
                                                                                      df['Blocked_Pipe'] |
                                                                                      df['Blocked_Nozzle']) &
                                                                                      detarr[i]],
                                                                      df['Energy'][(df['Blocked_Cone'] |
                                                                                   df['Blocked_Pipe'] |
                                                                                   df['Blocked_Nozzle']) &
                                                                                   detarr[i]],
                                                                      bins=(binsz, binse))

                            blockedtvz, zbins, tbins = np.histogram2d(df['zpos_final'][(df['Blocked_Cone'] |
                                                                                      df['Blocked_Pipe'] |
                                                                                      df['Blocked_Nozzle']) &
                                                                                      detarr[i]],
                                                                      df['Theta_Deg'][(df['Blocked_Cone'] |
                                                                                   df['Blocked_Pipe'] |
                                                                                   df['Blocked_Nozzle']) &
                                                                                   detarr[i]],
                                                                      bins=(binsz, binstheta))

                            blockedtvphi, pbins, tbins = np.histogram2d(df['Phi_Deg'][(df['Blocked_Cone'] |
                                                                                      df['Blocked_Pipe'] |
                                                                                      df['Blocked_Nozzle']) &
                                                                                      detarr[i]],
                                                                        df['Theta_Deg'][(df['Blocked_Cone'] |
                                                                                        df['Blocked_Pipe'] |
                                                                                        df['Blocked_Nozzle']) &
                                                                                        detarr[i]],
                                                                        bins=(binsphi, binstheta))

                            coneblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Cone'] &
                                                                                          detarr[i]],
                                                                          df['Energy'][df['Blocked_Cone'] & detarr[i]],
                                                                          bins=(binstheta, binse))
                            pipeblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Pipe'] &
                                                                                          detarr[i]],
                                                                          df['Energy'][df['Blocked_Pipe'] & detarr[i]],
                                                                          bins=(binstheta, binse))
                            nozzleblockedevt, tbins, ebins = np.histogram2d(df['Theta_Deg'][df['Blocked_Nozzle'] &
                                                                                            detarr[i]],
                                                                            df['Energy'][df['Blocked_Nozzle'] &
                                                                                         detarr[i]],
                                                                            bins=(binstheta, binse))

                            coneblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][df['Blocked_Cone'] &
                                                                                           detarr[i]],
                                                                          df['Energy'][df['Blocked_Cone'] & detarr[i]],
                                                                          bins=(binsz, binse))
                            pipeblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][df['Blocked_Pipe'] &
                                                                                           detarr[i]],
                                                                          df['Energy'][df['Blocked_Pipe'] & detarr[i]],
                                                                          bins=(binsz, binse))
                            nozzleblockedevz, zbins, ebins = np.histogram2d(df['zpos_final'][df['Blocked_Nozzle'] &
                                                                                             detarr[i]],
                                                                            df['Energy'][df['Blocked_Nozzle'] &
                                                                                         detarr[i]],
                                                                            bins=(binsz, binse))

                            coneblockedcvz, zbins, cbins = np.histogram2d(df['zpos_final'][df['Blocked_Cone'] &
                                                                                           detarr[i]],
                                                                          df['CM_Deg'][df['Blocked_Cone'] & detarr[i]],
                                                                          bins=(binsz, binscm))
                            pipeblockedcvz, zbins, cbins = np.histogram2d(df['zpos_final'][df['Blocked_Pipe'] &
                                                                                           detarr[i]],
                                                                          df['CM_Deg'][df['Blocked_Pipe'] & detarr[i]],
                                                                          bins=(binsz, binscm))
                            nozzleblockedcvz, zbins, cbins = np.histogram2d(df['zpos_final'][df['Blocked_Nozzle'] &
                                                                                             detarr[i]],
                                                                            df['CM_Deg'][df['Blocked_Nozzle'] &
                                                                                         detarr[i]],
                                                                            bins=(binsz, binscm))

                            coneblockedtvz, zbins, tbins = np.histogram2d(df['zpos_final'][df['Blocked_Cone'] &
                                                                                           detarr[i]],
                                                                          df['Theta_Deg'][df['Blocked_Cone'] & detarr[i]],
                                                                          bins=(binsz, binstheta))
                            pipeblockedtvz, zbins, tbins = np.histogram2d(df['zpos_final'][df['Blocked_Pipe'] &
                                                                                           detarr[i]],
                                                                          df['Theta_Deg'][df['Blocked_Pipe'] & detarr[i]],
                                                                          bins=(binsz, binstheta))
                            nozzleblockedtvz, zbins, tbins = np.histogram2d(df['zpos_final'][df['Blocked_Nozzle'] &
                                                                                             detarr[i]],
                                                                            df['Theta_Deg'][df['Blocked_Nozzle'] &
                                                                                         detarr[i]],
                                                                            bins=(binsz, binstheta))
                        except KeyError:
                            print("KeyError suppressed.")

                        # Here we get the ratios that we want to use to make the contour plots. Instead of unblocked
                        # minus
                        # blocked over unblocked we could have used the unblocked mask and just done
                        # "detected/unblocked where "detected" comes from the unblocked mask.

                        # The ratios are also split up into separate ratios for the particles only
                        # blocked by the cone, etc.

                        # ratio = np.divide(blockedevt, unblockedevt, out=np.zeros_like(blockedevt),
                        # where=blockedevt != 0)
                        ratioevt = np.divide((unblockedevt - blockedevt), unblockedevt, out=np.zeros_like(blockedevt),
                                             where=unblockedevt != 0)

                        ratioevz = np.divide((unblockedevz - blockedevz), unblockedevz, out=np.zeros_like(blockedevz),
                                             where=unblockedevz != 0)

                        ratiocvz = np.divide((unblockedcvz - blockedcvz), unblockedcvz, out=np.zeros_like(blockedcvz),
                                             where=unblockedcvz != 0)

                        ratiotvz = np.divide((unblockedtvz - blockedtvz), unblockedtvz, out=np.zeros_like(blockedtvz),
                                             where=unblockedtvz != 0)

                        ratiotvp = np.divide((unblockedtvphi - blockedtvphi), unblockedtvphi,
                                             out=np.zeros_like(blockedtvphi), where=unblockedtvphi != 0)

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

                        ratioconecvz = np.divide((unblockedcvz - coneblockedcvz), unblockedcvz,
                                                 out=np.zeros_like(coneblockedcvz), where=unblockedcvz != 0)

                        ratiopipecvz = np.divide((unblockedcvz - pipeblockedcvz), unblockedcvz,
                                                 out=np.zeros_like(pipeblockedcvz), where=unblockedcvz != 0)

                        rationozzlecvz = np.divide((unblockedcvz - nozzleblockedcvz), unblockedcvz,
                                                   out=np.zeros_like(nozzleblockedcvz), where=unblockedcvz != 0)

                        ratioconetvz = np.divide((unblockedtvz - coneblockedtvz), unblockedtvz,
                                                 out=np.zeros_like(coneblockedtvz), where=unblockedtvz != 0)

                        ratiopipetvz = np.divide((unblockedtvz - pipeblockedtvz), unblockedtvz,
                                                 out=np.zeros_like(pipeblockedtvz), where=unblockedtvz != 0)

                        rationozzletvz = np.divide((unblockedtvz - nozzleblockedtvz), unblockedtvz,
                                                   out=np.zeros_like(nozzleblockedtvz), where=unblockedtvz != 0)

                        # To actually plot the ratios into histograms we have to transpose the binned arrays:

                        ratioevt = ratioevt.T
                        ratioevz = ratioevz.T
                        ratiotvp = ratiotvp
                        ratiocvz = ratiocvz.T
                        ratiotvz = ratiotvz.T

                        ratioconeevt = ratioconeevt.T
                        ratiopipeevt = ratiopipeevt.T
                        rationozzleevt = rationozzleevt.T

                        ratioconeevz = ratioconeevz.T
                        ratiopipeevz = ratiopipeevz.T
                        rationozzleevz = rationozzleevz.T

                        ratioconecvz = ratioconecvz.T
                        ratiopipecvz = ratiopipecvz.T
                        rationozzlecvz = rationozzlecvz.T

                        ratioconetvz = ratioconetvz.T
                        ratiopipetvz = ratiopipetvz.T
                        rationozzletvz = rationozzletvz.T

                        # As mentioned, tbins, ebins, and zbins are the bin edges. Here initialize a new array:

                        tbins2 = np.zeros(90)
                        cbins2 = np.zeros(90)
                        ebins2 = np.zeros(150)
                        zbins2 = np.zeros(100)
                        pbins2 = np.zeros(numphibins)

                        # And here get the bin centers by taking the averages of three bin edges.

                        for k in range(180):
                            if k < 90:
                                tbins2[k] = (tbins[k] + tbins[k + 1]) / 2
                                cbins2[k] = (cbins[k] + cbins[k + 1]) / 2
                            if k < numphibins:
                                pbins2[k] = (pbins[k] + pbins[k + 1]) / 2
                            if k < 100:
                                zbins2[k] = (zbins[k] + zbins[k + 1]) / 2
                            if k < 150:
                                ebins2[k] = (ebins[k] + ebins[k + 1]) / 2

                        # To plot them we need to make a mesh grid of the bins:

                        xevt, yevt = np.meshgrid(tbins2, ebins2)
                        xevz, yevz = np.meshgrid(zbins2, ebins2)
                        xcvz, ycvz = np.meshgrid(zbins2, cbins2)
                        xtvz, ytvz = np.meshgrid(zbins2, tbins2)
                        xtvp, ytvp = np.meshgrid(tbins2, pbins2 * np.pi/180)

                        # Now, since we have so many bins, the contour plots won't look nice. If we put too few bins the
                        # contours also don't look great. So, the solution is to use a lot of bins
                        # and use these gaussian filters on the ratios to smooth them out.

                        ratioevt_blurr = ndimage.gaussian_filter(ratioevt, sigma=1.5, order=0)
                        ratioevz_blurr = ndimage.gaussian_filter(ratioevz, sigma=1.5, order=0)
                        ratiocvz_blurr = ndimage.gaussian_filter(ratiocvz, sigma=1.5, order=0)
                        ratiotvz_blurr = ndimage.gaussian_filter(ratiotvz, sigma=1.5, order=0)
                        ratiotvp_blurr = ndimage.gaussian_filter(ratiotvp, sigma=1, order=0)

                        ratioconeevt_blurr = ndimage.gaussian_filter(ratioconeevt, sigma=1.5, order=0)
                        ratiopipeevt_blurr = ndimage.gaussian_filter(ratiopipeevt, sigma=1.5, order=0)
                        rationozzleevt_blurr = ndimage.gaussian_filter(rationozzleevt, sigma=1.5, order=0)

                        ratioconeevz_blurr = ndimage.gaussian_filter(ratioconeevz, sigma=1.6, order=0)
                        ratiopipeevz_blurr = ndimage.gaussian_filter(ratiopipeevz, sigma=1.6, order=0)
                        rationozzleevz_blurr = ndimage.gaussian_filter(rationozzleevz, sigma=1.6, order=0)

                        ratioconecvz_blurr = ndimage.gaussian_filter(ratioconecvz, sigma=1.5, order=0)
                        ratiopipecvz_blurr = ndimage.gaussian_filter(ratiopipecvz, sigma=1.5, order=0)
                        rationozzlecvz_blurr = ndimage.gaussian_filter(rationozzlecvz, sigma=1.5, order=0)

                        ratioconetvz_blurr = ndimage.gaussian_filter(ratioconetvz, sigma=1.5, order=0)
                        ratiopipetvz_blurr = ndimage.gaussian_filter(ratiopipetvz, sigma=1.5, order=0)
                        rationozzletvz_blurr = ndimage.gaussian_filter(rationozzletvz, sigma=1.5, order=0)

                        # 13 is the contour plot of percentage of detected particles.
                        if int(aa) < 5 and int(dd) == 9:
                            # Energy vs z contour plot parameters:
                            if int(aa) == 1:
                                xdata, ydata = xevz, yevz
                                if int(cc) == 1:
                                    ratiodata = ratioevz_blurr
                                    cmp = 'GnBu'
                                    title = df['Reaction'][0] + " Fraction of Particles Detected, B = " + \
                                            str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 2:
                                    ratiodata = ratioconeevz_blurr
                                    cmp = 'Greens'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by " \
                                                                "the Cone, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 3:
                                    ratiodata = ratiopipeevz_blurr
                                    cmp = 'Reds'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Pipe, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 4:
                                    ratiodata = rationozzleevz_blurr
                                    cmp = 'Blues'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Nozzle, B = " + str(df['Magnetic Field'][0]) + " T"
                                if int(bb) == 0:
                                    if 1 <= int(cc) <= 2:
                                        conts = [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1]
                                    if 3 <= int(cc) <= 4:
                                        conts = [.5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1]
                                else:
                                    conts = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8,
                                             .85, .9, .95, 1]

                                color = 'k'
                                lnwdth = 0.2
                                xlabel, ylabel = 'z (m)', 'Energy (MeV)'

                            # Energy vs theta lab contour plot parameters:
                            if int(aa) == 2:
                                xdata, ydata = xevt, yevt
                                if int(cc) == 1:
                                    ratiodata = ratioevt_blurr
                                    cmp = 'YlGnBu'
                                    title = df['Reaction'][0] + " Fraction of Particles Detected, B = " + \
                                            str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 2:
                                    ratiodata = ratioconeevt_blurr
                                    cmp = 'Greens'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by " \
                                                                "the Cone, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 3:
                                    ratiodata = ratiopipeevt_blurr
                                    cmp = 'Reds'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Pipe, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 4:
                                    ratiodata = rationozzleevt_blurr
                                    cmp = 'Blues'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Nozzle, B = " + str(df['Magnetic Field'][0]) + " T"
                                if int(bb) == 0:
                                    if 1 <= int(cc) <= 2:
                                        conts = [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1]
                                    if 3 <= int(cc) <= 4:
                                        conts = [.65, .7, .75, .8, .85, .9, .95, 1]
                                else:
                                    conts = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8,
                                             .85, .9, .95, 1]
                                color = 'k'
                                lnwdth = 0.2
                                xlabel, ylabel = 'Lab Angle (Deg)', 'Energy (MeV)'

                            # Lab angle vs z contour plot parameters:
                            if int(aa) == 3:
                                xdata, ydata = xtvz, ytvz
                                if int(cc) == 1:
                                    ratiodata = ratiotvz_blurr
                                    cmp = 'YlGnBu'
                                    title = df['Reaction'][0] + " Fraction of Particles Detected, B = " + \
                                            str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 2:
                                    ratiodata = ratioconetvz_blurr
                                    cmp = 'Greens'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by " \
                                                                "the Cone, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 3:
                                    ratiodata = ratiopipetvz_blurr
                                    cmp = 'Reds'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Pipe, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 4:
                                    ratiodata = rationozzletvz_blurr
                                    cmp = 'Blues'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Nozzle, B = " + str(df['Magnetic Field'][0]) + " T"
                                if int(bb) == 0:
                                    if 1 <= int(cc) <= 2:
                                        conts = [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1]
                                    if 3 <= int(cc) <= 4:
                                        conts = [.65, .7, .75, .8, .85, .9, .95, 1]
                                else:
                                    conts = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8,
                                             .85, .9, .95, 1]
                                color = 'k'
                                lnwdth = 0.2
                                xlabel, ylabel = 'z (m)', 'Lab Angle (Deg)'

                            # CM Angle vs z
                            if int(aa) == 4:
                                xdata, ydata = xcvz, ycvz
                                if int(cc) == 1:
                                    ratiodata = ratiocvz_blurr
                                    cmp = 'YlGnBu'
                                    title = df['Reaction'][0] + " Fraction of Particles Detected, B = " + \
                                            str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 2:
                                    ratiodata = ratioconecvz_blurr
                                    cmp = 'Greens'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by " \
                                                                "the Cone, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 3:
                                    ratiodata = ratiopipecvz_blurr
                                    cmp = 'Reds'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Pipe, B = " + str(df['Magnetic Field'][0]) + " T"
                                elif int(cc) == 4:
                                    ratiodata = rationozzlecvz_blurr
                                    cmp = 'Blues'
                                    title = df['Reaction'][0] + " Fraction of Particles Not Blocked by the " \
                                                                "Nozzle, B = " + str(df['Magnetic Field'][0]) + " T"
                                if int(bb) == 0:
                                    if 1 <= int(cc) <= 2:
                                        conts = [.25, .3, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1]
                                    if 3 <= int(cc) <= 4:
                                        conts = [.65, .7, .75, .8, .85, .9, .95, 1]
                                else:
                                    conts = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8,
                                             .85, .9, .95, 1]
                                color = 'k'
                                lnwdth = 0.2
                                xlabel, ylabel = 'z (m)', 'CM Angle (Deg)'

                        if int(aa) < 5 and int(dd) == 9 and itest:
                            cf = axes.contourf(xdata, ydata, ratiodata, conts, cmap=cmp)
                            cs = axes.contour(xdata, ydata, ratiodata, conts, colors=color, linewidths=lnwdth)

                            cbar = fig.colorbar(cf, ax=axes)
                            axes.set_xlabel(xlabel)
                            axes.set_ylabel(ylabel)
                            plt.suptitle(title, fontsize=18)

                        # Phi vs Lab Angle contour
                        if plotnum == "6.1.9.9" and i == 0:
                            plt.rc('axes', labelsize=16)
                            plt.rc('xtick', labelsize=16)
                            plt.rc('ytick', labelsize=16)
                            fig2, ax = plt.subplots(subplot_kw=dict(projection='polar'))
                            cf = ax.contourf(ytvp, xtvp, ratiotvp_blurr, cmap='hot')
                            ax.contour(ytvp, xtvp, ratiotvp_blurr, colors='k', linewidths=.2)
                            plt.xlabel('Initial Phi Angle (Deg)')
                            plt.ylabel('Lab Angle (Deg)', rotation=0, size=14, labelpad=-370)
                            cbar = fig2.colorbar(cf)
                            title29 = df['Reaction'][0] + " Percentage of Particles Detected, B = " + \
                                      str(df['Magnetic Field'][0]) + " T"
                            plt.suptitle(title29, fontsize=18)

                        # Here we do something a little different, and instead of contour plots we do ratio plots,
                        # still using the bins that were defined above.

                if int(aa) == 7:
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
                    if plotnum == "7.1.9.9":
                        axs.plot(tbins2, divt, marker='o')
                        axs.legend(['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4'],
                                   loc='lower right', fontsize=16)
                        axs.set_xlabel('Lab Angle (Deg)')
                        axs.set_ylabel('Fraction of Particles Detected')

                    if plotnum == "7.2.9.9":
                        axs.plot(ebins2, dive, marker='o')
                        axs.legend(['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4'],
                                   loc='lower left', fontsize=16)
                        axs.set_xlabel('Energy (MeV)')
                        axs.set_ylabel('Fraction of Particles Detected')

                    if plotnum == "7.3.9.9":
                        axs.plot(zbins2, divz, marker='o')
                        axs.legend(['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4'],
                                   loc='lower left', fontsize=16)
                        axs.set_xlabel('z (m)')
                        axs.set_ylabel('Fraction of Particles Detected')

                    if plotnum == "7.4.9.9":

                        stdcolors = ['Blue', 'Orange', 'limegreen', 'Red', 'Purple']
                        labels32 = ['Total', 'Detector 1', 'Detector 2', 'Detector 3', 'Detector 4']
                        if i == 0:
                            plt.rc('axes', labelsize=16)
                            plt.rc('xtick', labelsize=16)
                            plt.rc('ytick', labelsize=16)
                            ax32 = plt.subplot(111, projection='polar')
                        pbins2 = pbins2 * np.pi / 180
                        if i !=2 and i != 3:
                            if i > 0:
                                ax32.plot(pbins2, divp, marker='o', color=stdcolors[i], label=labels32[i],
                                          markersize=8, MarkerEdgeColor='Black', alpha=0.7)
                            if i == 0:
                                ax32.plot(pbins2, divp, marker='o', color=stdcolors[i], label=labels32[i],
                                          markersize=14, MarkerEdgeColor='Black', alpha=0.7)
                        if i == 2 or i == 3:
                            pbins2gt0 = pbins2[pbins2 < np.pi]
                            divpgt0 = divp[pbins2 < np.pi]

                            pbins2lt0 = pbins2[pbins2 > np.pi]
                            divplt0 = divp[pbins2 > np.pi]

                            ax32.plot(pbins2gt0, divpgt0, marker='o', color=stdcolors[i], label=labels32[i],
                                      markersize=8, MarkerEdgeColor='Black', alpha=0.7)
                            ax32.plot(pbins2lt0, divplt0, marker='o', color=stdcolors[i], markersize=8,
                                      MarkerEdgeColor='Black', alpha=0.7)

                        ax32.grid(True)
                        plt.legend(fontsize=16, bbox_to_anchor=(1.0, .9))
                        plt.xlabel('Initial Phi Angle (Deg)')
                        plt.ylabel('Fraction of Particles Detected', rotation=0, size=14, labelpad=-370)

                # This next section is for solid targets only:
                if sol:

                    if (plotnum == "1.1.0.0" or plotnum == "1.1.0.1") and i > 0:

                        axi[i].hist2d(df['zpos_final'][detarr[i] & df["UnblockedSolidTarg"]],
                                    df['Energy'][detarr[i] & df["UnblockedSolidTarg"]], bins=(750, 750),
                                    range=[[zmin, zmax], [0, emax]], cmap=newcmpBlack)
                        axi[i].set_xlabel('z(m)')
                        axi[i].set_ylabel('Energy (MeV)')

                        if plotnum == "1.1.0.1":
                            axi[i].hist2d(df['zpos_final'][detarr[i] & ~df["UnblockedSolidTarg"]],
                                       df['Energy'][detarr[i] & ~df["UnblockedSolidTarg"]], bins=(750, 750),
                                       range=[[zmin, zmax], [0, emax]], cmap=newcmpRed)
                            axi[i].set_xlabel('z(m)')
                            axi[i].set_ylabel('Energy (MeV)')

                    if (plotnum == "1.1.0.2") and i > 0:

                        axi[i].hist2d(df['zpos_final'][detarr[i] & ~df["UnblockedSolidTarg"]],
                                    df['Energy'][detarr[i] & ~df["UnblockedSolidTarg"]], bins=(750, 750),
                                    range=[[zmin, zmax], [0, emax]], cmap=newcmpBlack)
                        axi[i].set_xlabel('z(m)')
                        axi[i].set_ylabel('Energy (MeV)')

                    if plotnum == "5.0.0.0" and i == 0:
                        axs.hist(df['Ex_Reconstructed'][df["UnblockedSolidTarg"] & detarr[i]], bins=750,
                                 range=[-0.2, exmax])
                        axs.set_xlabel('Excitation Energy (MeV)')
                        axs.set_ylabel('Counts')

                        # 35 is the same as 11 but also showing blocked particles.
                    if plotnum == "5.1.0.1" and i > 0:

                        axi[i].hist((df['Ex_Reconstructed'][df["UnblockedSolidTarg"] & detarr[i]],
                                  df['Ex_Reconstructed'][~df["UnblockedSolidTarg"] & detarr[i]]),
                                 bins=1000, range=[-0.2, exmax], color=(blk, red), stacked=True)
                        axi[i].set_xlabel('Excitation Energy (MeV)')
                        axi[i].set_ylabel('Counts')

                    if plotnum == "5.1.0.0" and i > 0:

                        axi[i].hist(df['Ex_Reconstructed'][df["UnblockedSolidTarg"] & detarr[i]], bins=750,
                                 range=[-0.2, exmax])
                        axi[i].set_xlabel('Excitation Energy (MeV)')
                        axi[i].set_ylabel('Counts')

                # Handle the legend here for each plot that needs it.
                if i == 1 and int(dd) == 1:
                    handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in [blk, grn, red, blu]]
                    labels = ["Unblocked", "Cone", "Pipe", "Nozzle"]
                    if plotnum == "11":
                        plt.legend(handles, labels, bbox_to_anchor=(1.0, .9), ncol=4)
                    else:
                        plt.legend(handles, labels, bbox_to_anchor=(1.4, 1.1), ncol=4)

        elif aa == "0":
            switch = 1
        elif aa == "400":
            # Here, we'll export the pandas data frame to a text file for use with Tableau or another visualization tool
            # We need to turn the masks into strings instead.
            df["Detector Selection"] = ""
            df["z-axis Detectors"] = ""
            df["Shadowing"] = ""
            detznames = ["z1", "z2", "z3", "z4", "z5", "z6"]
            detnames = ["Detector 2", "Detector 1", "Detector 3", "Detector 4"]
            detzarr = [df["Detz1"], df["Detz2"], df["Detz3"], df["Detz4"], df["Detz5"], df["Detz6"]]
            blkmasks = [df["Unblocked"], df["Blocked_Cone"], df["Blocked_Pipe"], df["Blocked_Nozzle"]]
            unb = ["Unblocked", "Cone", "Pipe", "Nozzle"]
            for j in range(4):
                df["Detector Selection"] = np.where(detarr[j+1], detnames[j], df["Detector Selection"])
                df["Shadowing"] = np.where(blkmasks[j], unb[j], df["Shadowing"])
            for k in range(6):
                df["z-axis Detectors"] = np.where(detzarr[k], detznames[k], df["z-axis Detectors"])

            empmask1 = df["z-axis Detectors"] == ""
            empmask2 = df["Shadowing"] == ""
            df["z-axis Detectors"] = np.where(empmask1, "Not Detected", df["z-axis Detectors"])
            df["Shadowing"] = np.where(empmask2, "Magnet Bore", df["Shadowing"])

            csvfile = pklin[:-4] + ".csv"

            print("Saving DataFrame to csv...")
            df.to_csv(csvfile, index=True, header=True)
            print("csv file created in ./Output_Files!")

        if aa != '100' and aa != '101':
            lastentry = plotnum

        if currentry == '300':
            print("SWITCHBACKDF")
            overlaybool = False
            df = df_main
            if detposbool:
                detarr = [df["AllPossible"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                                               df['Detz5'] | df['Detz6']),
                          df["Det2"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                                        df['Detz5'] | df['Detz6']),
                          df["Det1"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                                        df['Detz5'] | df['Detz6']),
                          df["Det3"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                                        df['Detz5'] | df['Detz6']),
                          df["Det4"] & (df['Detz1'] | df['Detz2'] | df['Detz3'] | df['Detz4'] |
                                        df['Detz5'] | df['Detz6'])]
            else:
                detarr = [df["AllPossible"], df["Det2"], df["Det1"], df["Det3"], df["Det4"]]

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

    outdir = "./Output_Files/"

    # Get a list of all the pkl files that are already created:
    list_pickles = glob.glob(outdir + '*evt*.pkl')
    filein2 = ""

    # If there aren't any pickles, prompt the user to make one.
    if len(list_pickles) == 0:
        print("\nIt appears that no simulated DataFrame exists, run SOLSTISE_Sim.py first.")
    else:
        # If there are pickles, get the pickle that was modified last:
        latest_file = max(list_pickles, key=os.path.getctime)
        latest_file = latest_file[15:]
        print("\nThe most recently created simulated DataFrame file is: " + latest_file)
        yn = input("\nWould you like to use this file? [Y/N]: ")
        if yn == "N" or yn == "n":
            if len(list_pickles) == 1:
                # If the user wants to use a different pickle file but there is only one, make the user rum the sim.
                print("\nYou only have one DataFrame file. To create another, run SOLSTISE_Sim.py")
                sys.exit()
            # if there is more than one, list them and give them a number so the user can choose.
            if len(list_pickles) > 1:
                print("\n")
                for i in range(len(list_pickles)):
                    print(str(i + 1) + ") " + list_pickles[i][15:])

                filenum = 1000000
                while filenum > len(list_pickles) or filenum == 0:
                    filenum = int(input("\nChoose a number from the list: "))
                    if filenum == 0 or len(list_pickles) < filenum:
                        print("ERROR: Number entered does not correspond to a DataFrame file...")
                    else:
                        filein = list_pickles[filenum-1][15:]
        else:
            filein = latest_file

        if fnmatch.fnmatch(filein, "*allE*"):
            yn2 = input("\nWould you like to load a file to overlay on top of the contour histograms? [Y/N]: ")

            if yn2 == "N" or yn2 == "n":
                print("")
            else:
                list_overlays = glob.glob(outdir + filein[:-12] + '*.pkl')
                print("\nFiles that can be used: ")
                for i in range(len(list_overlays)):
                    if not fnmatch.fnmatch(list_overlays[i], "*allE*"):
                        print(str(i + 1) + ") " + list_overlays[i][15:])

                filenum2 = int(input("\nChoose a number from the list: "))
                if filenum2 == 0 or len(list_overlays) < filenum2:
                    print("ERROR: Number entered does not correspond to a DataFrame file...")
                else:
                    filein2 = list_overlays[filenum2 - 1][15:]

        print("\nThe file to be used is: " + filein)
        plot(outdir + filein, outdir + filein2)
