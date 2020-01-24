import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle
from massreader import readmass
import pandas as pd
from mpl_toolkits import mplot3d


def sim_pd(rbore, rblock, cheight, phi1block, phi2block, ebeam, filein, reac):

    # The z axis points in beam direction, the x-axis points to the left, and the y-axis points down

    masses = readmass(reac)

    # reaction of form t(b,e)R
    utoMeV = 931.4941

    mt = masses[0]
    mb = masses[1]
    me = masses[2]
    mr = masses[3]

    qvalnoex = (mt + mb - me - mr)*utoMeV

    #ebeam = 168 # MeV, for d(28Si,p) it is 6 MeV/u

    rblock = np.array(rblock, dtype=np.float64)
    cheight = np.array(cheight, dtype=np.float64)
    phi1block = np.array(phi1block, dtype=np.float64)
    phi2block = np.array(phi2block, dtype=np.float64)

    mevtoj = 1.6021766e-13
    c = 2.998e8

    # finite radius of the detector array
    r0 = 0.011

    #tcm is a constant
    tcm = mt/(mb+mt)*ebeam

    amutokg = 1.66053907e-27  # amu to kg conversion
    B = 1.915  # teslas
    q = 1.6e-19  # 1 elemental charge in coulombs


    # Generates a pandas data frame of shape (xxx,2) whose columns are theta angle and energy.
    df = pd.read_csv(filein, sep="\t", header=None, low_memory = False)
    df.columns = ["Theta_Deg", "Energy"]

    # Determine whether or not the particles are coming out at backward or forward angles:
    if df["Theta_Deg"].mean() > 90:
        invkin = True
    else:
        invkin = False

    # Convert the angle to radians and add it to the dataframe
    df['Theta_Rad'] = df['Theta_Deg'] * np.pi/180

    # Cyclotron frequency and period.
    omega = (q * B) / (me * amutokg)
    tcyc = (2 * np.pi) / omega

    # Velocity of the ejectile in the lab frame (m/s) added to dataframe
    df['vel_ejec'] = np.sqrt((2 * df['Energy'] * mevtoj) / (me * amutokg))
    # Velocities parallel to the z-axis and perpendicular to the z-axis.
    df['vel_perp'] = df['vel_ejec'] * np.sin(df['Theta_Rad'])
    df['vel_par'] = df['vel_ejec'] * np.cos(df['Theta_Rad'])

    # Calculates Q-Ex from the energy and angle of the ejected particle (Krane 11.10)
    df['Q-Ex'] = df['Energy'] * (1 + me / mr) - ebeam * (1 - mb / mr) - 2 * np.sqrt(
        mb * me / mr**2 * df['Energy'] * ebeam) * np.cos(df['Theta_Rad'])

    # The total energy of the reaction, Tcm + Q-value - Ex
    df['E_tot'] = tcm + df['Q-Ex']

    # The velocity of the ejectile in the CM frame.
    df['v0'] = np.sqrt(2 * mr * utoMeV * df['E_tot'] / (me * utoMeV * (me * utoMeV + mr * utoMeV))) * c

    # The velocity of the CM frame, it's a constant and dependent on the beam energy and species.
    vcm = np.sqrt(2 * ebeam / (mb * utoMeV)) * (mb * utoMeV/(mb * utoMeV+mt * utoMeV)) * c

    # Calculates the CM angle of the ejectile and recoil (they're the same). arccos takes values from -1 to 1 and won't
    # break if it gets a value outside of that range.
    # However, we'll get rid of all the NaN entries before the arccos line with a mask.
    df['cosarg'] = (df['vel_ejec']**2 - df['v0']**2 - vcm**2)/(2 * df['v0'] * vcm)
    maskarg = (df['cosarg'] >= -1) & (df['cosarg'] <= 1)

    df = df[maskarg]

    df = df.reset_index(drop=True)
    #print(df)

    df['Theta_CM'] = np.arccos((df['vel_ejec']**2 - df['v0']**2 - vcm**2)/(2 * df['v0'] * vcm))

    # Reduced cyclotron frequency because of the finite size of the detector array.
    df['t_reduced'] = tcyc - r0/(df['v0']*np.sin(df['Theta_CM']))

    # makes a phi array the same size as the theta array, random number 0 to 1
    #np.random.seed = 18
    phi = np.random.rand(len(df))
    # then multiply the phi array by 2pi to get a real phi value and put it into the dataframe
    df['Phi'] = phi * 2 * np.pi

    #print(df['Phi'])

    # creates a mask the same shape as the energy array
    # maskmaster is the mask that keeps track of the overall mask in the loop
    maskmaster = df['Energy'] > 0
    # maskmasters for the pipe and cone defined here.
    maskmaster_pipe = df['Energy'] > 0
    maskmaster_cone = df['Energy'] > 0
    maskmaster_nozzle = df['Energy'] > 0
    # maskrbore takes care of the mask for the bore radius because some particles will hit that
    maskrbore = df['Energy'] > 0
    # maskcone initialized like maskrbore
    masknozzle = df['Energy'] > 0
    # this initializes phic, which tracks the current position of the phi particle.
    phic = df['Phi']

    # ***************************************************************************************
    # Parameters of the nozzle, cone, and pipe get entered here:

    nozzleconedistin = 3.83  # dist between nozzle and cone in inches
    reacdistbelownozzle = 0.09843  # dist below nozzle the reaction happens in inches
    conedia = 2.6  # cone outer diameter in inches
    coneheight = 3.02  # cone height in inches as measured from the top of the ISO base.

    # Height above the cone that the reaction occursfrom massreader import readmass
    reacheight = ((nozzleconedistin - reacdistbelownozzle) * 2.54) / 100
    rcone = ((conedia / 2) * 2.54) / 100

    # ISO base is 5.12 inches outer diameter, below is converted to meters.
    baseheight = reacheight + coneheight * 2.54 / 100
    rISObase = (5.12 / 2) * 2.54 / 100

    # distance from the reaction that the cone side equation starts (this equation starts at ~5.5, now 3.83).
    # sideheight = (5.5 - reacdistbelownozzle) * 2.54 / 100
    sideheight = (3.83 - reacdistbelownozzle) * 2.54 / 100

    # Polynomial coefficients for the current best cone: Top of cone is 3.83 in away from nozzle

    poly3 = [0.0611, -0.8077, 3.5255, -3.7867]

    rconeside = lambda y: (poly3[0] * (y + reacdistbelownozzle) ** 3 + poly3[1] * (y + reacdistbelownozzle) ** 2 +
                           poly3[2] * (y + reacdistbelownozzle) + poly3[3]) * 2.54 / 100

    # parameters for nozzle shadowing here:
    # cone tip dist is the extra vertical height added onto the nozzle "cone" shape if it extended out to a point
    # it was determined using the 22 degree slope of the nozzle side and the nozzle opening radius 0.0455 inches
    nozzletipdist = 0.112616
    nozzleang = 22 * np.pi / 180  # the nozzle slope angle, which is 22 degrees (converted to rads)

    rnozzle = lambda y: ((y + np.abs(reacdistbelownozzle - nozzletipdist)) * np.tan(nozzleang)) * 2.54 / 100

    # function determines the r coordinates of the 2nd circle that makes up the gas pipe.
    rpipe = lambda ph: cheight*np.sin(ph) + np.sqrt(cheight**2*np.sin(ph)**2 - cheight**2 + rblock**2)

    # *********************************************************************************************

    dummy = df['Energy']

    # Simulating events status bar for the for loop
    print("Simulating Events...")
    statbar = "[                              ]"

    # Splits the flight time into 300 segments for tracking purposes to see whether or not the particle is blocked.
    for i in range(300):

        if i % 10 == 0:
            statbar = statbar.replace(" ", "=", 1)
            print(statbar, end='\r', flush=True)

        t = df['t_reduced']/300 * (i+1)
        xpos = (-(df['vel_perp']/omega)*np.cos((omega*t)+df['Phi']))+((df['vel_perp']/omega)*np.cos(df['Phi']))
        ypos = ((df['vel_perp']/omega)*np.sin(omega*t+df['Phi']))-df['vel_perp']/omega*np.sin(df['Phi'])
        zpos = df['vel_par']*t

        #print(ypos)

        # r is the radial position of the particle
        r = np.sqrt(xpos**2 + ypos**2)

        # phic is the phi current position of the particle.
        phic = np.arctan2(ypos, xpos) + np.pi

        # rxzplane determines the radial position in the xz-plane
        rxzplane = np.sqrt(xpos**2 + zpos**2)

        # rpipe determines the r position of the 2nd circle boundary
        # so if the particle radius is greater than that, it gets blocked
        maskrpipe = (r > rpipe(phic))

        # maskphipipe is the mask that determines whether or not the particle is within the phi boundaries of the pipe
        # if maskphipipe and maskrpipe are true, then the particle is blocked by the pipe
        maskphipipe = (phic > phi1block) & (phic < phi2block)

        # maskcone determines if the particle is within the opening of the cone
        # Since if statements on the arrays are so slow, we'll break up the cone mask into three: tube (top), sides,
        # and ISO base. In the future, the base could be gotten rid of depending on final geometry.
        # masktop = (rxzplane < rcone) & (ypos > reacheight) & (ypos < sideheight)
        # With the new cone, masktop does not need to be used anymore since the long straight neck at the top of
        # the cone has been removed.
        masksides = (rxzplane < rconeside(ypos * 100 / 2.54)) & (ypos > (sideheight)) & (ypos < baseheight)
        masktest = (rxzplane < rconeside(ypos * 100 / 2.54)) & (ypos > (sideheight))
        maskbase = (rxzplane < rISObase) & (ypos > baseheight)

        #print(rconeside(ypos * 100 / 2.54),"\n")
        #print(dummy[masksides])

        #print(rconeside(ypos * 100 / 2.54), rxzplane)
        #print(ypos[masktest])

        maskcone = masksides | maskbase

        # we want only particles that come out at backward angles
        if invkin:
            maskz = zpos < 0
        else:
            maskz = zpos > 0

        # masknozzle determines if the the particle hits the nozzle.
        masknozzle = masknozzle*((rxzplane < rnozzle(-1 * ypos * 100 / 2.54)) & ((-1 * ypos) > reacdistbelownozzle *
                                                                                 2.54 / 100))

        #print(rxzplane[masksides] * 100 / 2.54)
        #print(maskr.shape, phir.shape, zpos[np.invert(maskz)].shape)

        #maskmaster = maskmaster*np.invert(maskrpipe & maskphipipe)#np.invert(maskcone)*np.invert(masky)#*np.invert(maskz)
        maskmaster = maskmaster*np.invert(maskcone)*np.invert(maskrpipe & maskphipipe)*np.invert(masknozzle)
        maskmaster_cone = maskmaster_cone*np.invert(maskcone)
        maskmaster_pipe = maskmaster_pipe*np.invert((maskrpipe & maskphipipe)*np.invert(maskcone)*np.invert(masknozzle))
        maskmaster_nozzle = maskmaster_nozzle * np.invert(masknozzle)
        #maskmaster = maskmaster*np.invert(maskcone | masknozzle | maskrpipe | maskphipipe)
        maskrbore = maskrbore & (r < rbore)

    # Adds the final phi position to the dataframe
    df['zpos_final'] = zpos
    # Adds the final phi position to the dataframe
    df['Phi_final'] = phic

    # splits the particles up onto four quadrants of a fictional detector
    df["Det1"] = (phic > 0) & (phic < np.pi/2)
    df["Det2"] = (phic > np.pi/2) & (phic < np.pi)
    df["Det3"] = (phic > np.pi) & (phic < 3*np.pi/2)
    df["Det4"] = (phic > 3*np.pi/2) & (phic < 2*np.pi)

    detarr = [df["Det2"], df["Det1"], df["Det3"], df["Det4"]]

    # Each element of df_det contains the info for the particles that hit that detector quadrant.
    # i.e. element 0 = quad 1, 1 = quad 2 etc..

    # Add boolean masks to the dataframe to pick out the specific particles we want.
    # All possible rejects all the particles that go in the wrong half of the magnet or hit the bore.
    df["AllPossible"] = maskz & maskrbore
    # All unshadowed particles
    df["Unblocked"] = maskmaster & maskz & maskrbore
    # All particles blocked by the cone
    df["Blocked_Cone"] = np.invert(maskmaster_cone) & maskz & maskrbore
    # All particles blocked by the nozzle
    df["Blocked_Nozzle"] = np.invert(maskmaster_nozzle) & maskz & maskrbore
    # All particles blocked by the pipe
    df["Blocked_Pipe"] = np.invert(maskmaster_pipe) & maskz & maskrbore

    # We'll also reconstruct the Q-value spectrum from the "detected" particles.
    # First we have to calculate the CM Energy:
    df['EnergyCM'] = df['Energy'] + .5 * me * utoMeV * (vcm / c) ** 2 - me * utoMeV * (vcm / c ** 2) / tcyc * \
                     df['zpos_final']
    df['Ex_Reconstructed'] = tcm + qvalnoex - df['EnergyCM'] * (me + mr) / mr

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
    #emin should always just be 0

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
                    binstheta = np.zeros(19)
                    for j in range(19):
                        binstheta[j] = j*5+90
                    # We'll also make the energy bins as well:
                    binse = np.zeros(26)
                    for j in range(26):
                        binse[j] = j*.4

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

                    #ratio = np.divide(blockedevt, unblockedevt, out=np.zeros_like(blockedevt), where=unblockedevt != 0)
                    ratio = np.divide((unblockedevt-blockedevt), unblockedevt, out=np.zeros_like(blockedevt), where=unblockedevt != 0)
                    #ratio = (unblockedevt-blockedevt)/unblockedevt
                    ratio = ratio.T

                    tbins2 = np.zeros(18)
                    ebins2 = np.zeros(25)

                    for k in range(25):
                        if k < 18:
                            tbins2[k] = (tbins[k] + tbins[k + 1]) / 2
                        ebins2[k] = (ebins[k] + ebins[k + 1]) / 2

                    X, Y = np.meshgrid(tbins2, ebins2)
                    X2, Y2 = np.meshgrid(tbins, ebins)

                    if plotnum == 12:
                        ax = plt.contourf(X, Y, ratio, 10, cmap='BuGn')
                        #plt.clabel(ax, inline=-3, fontsize=10)
                        #ax.colobar()
                    if plotnum == 13:
                        plt.pcolormesh(X2, Y2, ratio)

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





