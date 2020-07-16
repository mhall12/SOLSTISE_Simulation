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
import os
import glob
from Plotter import plot
import fnmatch
from stopyt import desorb
import pickle
import warnings
from random import randrange


def sim_pd(rbore, rblock, cheight, phi1block, phi2block, ebeam, filein, reac, conetxt, nozztxt,
           bfield, beamdiamm, jetrin, detzi):
    # suppress warnings that occur in the code calculations:
    warnings.filterwarnings("ignore")

    evtdir = "./Event_Files/"
    outdir = "./Output_Files/"
    geodir = "./Geometry_Files/"

    # Open the targetparms pkl:
    # targetparms now contains all the information needed for the desorb calculation.
    if fnmatch.fnmatch(filein, '*eloss*'):
        pklname = filein[:-16] + 'tgt_' + filein[-5] + '.pkl'

        with open(evtdir + pklname, 'rb') as f:
            targetparms = pickle.load(f)

        # The list has to contain other lists, so many of the parameters will require [][0]
        ebeam = targetparms[9][0]

        # Put the targetparms into useful variables:
        ztarg = targetparms[0]
        atarg = targetparms[1]
        numtarg = targetparms[2]
        density = targetparms[3][0]
        thickness = targetparms[4][0]
        jetpress = targetparms[5][0]
        jetrad = targetparms[6][0]
        champress = targetparms[7][0]
        gas = targetparms[8][0]
        if len(targetparms) == 12:
            ex1 = targetparms[10][0]
            ex2 = targetparms[11][0]
        else:
            ex1 = 1
            ex2 = 2
            print("***You're using an old file, so the Excitation Energy calculation will not work correctly!***")

        elossbool = True
    else:
        elossbool = False
        thickness = 'N/A'
        jetpress = 'N/A'
        champress = 'N/A'
        jetrad = 'N/A'

    if rblock < rbore:
        phi1block = 3/2*np.pi - np.arctan(rblock/(-1*cheight))
        phi2block = 3/2*np.pi + np.arctan(rblock/(-1*cheight))

    # Control which way the reaction products bend in the field based on the solenoid.
    # HELIOS particles bend CCW looking downstream at target, Ben thinks they go the same way in SOLARIS, but it can be
    # changed here easily if need be. cwsign = -1 is CCW, cwsign = 1 is CW.:
    # HELIOS
    if rbore == 0.46:
        cwsign = -1
    # SOLARIS
    if rbore == 0.45:
        cwsign = -1
    else:
        cwsign = -1

    # The z axis points in beam direction, the x-axis points to beam right, and the y-axis points up

    # Get the various parameters from readmass. In this case, we don't actually want z/atarg to be defined here
    # because they'll be wrong if it's a solid target
    masses, ztarg_buff, atarg_buff, zeject, aeject, zbeam, abeam = readmass(reac)


    # reaction of form t(b,e)R
    utoMeV = 931.4941

    mt = masses[0]
    mb = masses[1]
    me = masses[2]
    mr = masses[3]

    mevtoj = 1.6021766e-13
    c = 2.998e8

    qvalnoex = (mt + mb - me - mr)*utoMeV

    #ebeam = 168 # MeV, for d(28Si,p) it is 6 MeV/u

    rblock = np.array(rblock, dtype=np.float64)
    cheight = np.array(cheight, dtype=np.float64)
    phi1block = np.array(phi1block, dtype=np.float64)
    phi2block = np.array(phi2block, dtype=np.float64)

    # finite radius of the detector array
    r0 = 0.011

    #tcm is a constant
    tcm = mt/(mb+mt)*ebeam

    amutokg = 1.66053907e-27  # amu to kg conversion
    B = bfield
   # B = 1.915  # teslas
    q = 1.6e-19  # 1 elemental charge in coulombs

    # Generates a pandas data frame of shape (xxx,2) whose columns are theta angle and energy.
    df = pd.read_csv(evtdir + filein, sep="\t", header=None, low_memory=False)
    if not fnmatch.fnmatch(filein, '*eloss*'):
        df.columns = ["Theta_Deg", "Energy"]
    else:
        if fnmatch.fnmatch(filein, '*_g*'):
            # If it's an eloss text file, we have an additional column,
            # which is thickness for a solid or z offset for a jet.
            df.columns = ["Theta_Deg", "Energy", "z_Offset"]
        elif fnmatch.fnmatch(filein, '*_s*'):
            df.columns = ["Theta_Deg", "Energy", "Tgt_Thick"]

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
    phi = np.random.rand(len(df))
    # then multiply the phi array by 2pi to get a real phi value and put it into the dataframe
    df['Phi'] = phi * 2 * np.pi
    # Phi for debugging
    # df['Phi'] = np.zeros_like(df['Energy']) + np.pi

    # creates a mask the same shape as the energy array

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
    # ***************************************************************************************
    # ***************************************************************************************
    # ***************************************************************************************

    # Parameters of the nozzle, cone, and pipe get entered here:

    # Have to open the cone text file:
    conefile = open(geodir + conetxt, "r")

    # coneparms in order will have conenozzdist (in), coneheight (in), iso160height (in), conedia (in), polyorder,
    # ^3 coeff, ^2 coeff, ^1 coeff, ^0 coeff
    coneparms = conefile.readlines()

    nozzfile = open(geodir + nozztxt, "r")

    # nozzparms in order will have the reaction distance (mm), nozzle diameter (in), nozzle angle (deg), cone length (m)
    # cylinder radius (m), cylinder height (m), holder radius (m), and holder height (m).
    nozzparms = nozzfile.readlines()

    nozzleconedistin = float(coneparms[0])  # dist between nozzle and cone in inches
    reacdistbelownozzle = float(nozzparms[0]) / 10 / 2.54  # dist below nozzle the reaction happens in inches
    nozzconelen = float(nozzparms[3])  # nozzle cone length in m
    nozzcylrad = float(nozzparms[4]) / 2  # nozzle cylinder radius in m
    nozzcylh = float(nozzparms[5])  # nozzle cylinder height in m
    nozzholdrad = float(nozzparms[6]) / 2  # nozzle holder radius in m
    nozzholdh = float(nozzparms[7])  # nozzle holder height in m

    # Distance from the bottom of the nozzle holder to the bottom of the box
    # Distance from the bottom of the box to the top of the box
    nozzboxstart = float(nozzparms[8])
    nozzboxh = float(nozzparms[9])

    nozzholdcylstart = reacdistbelownozzle * 2.54 / 100 + nozzconelen + nozzcylh
    nozzholdcylend = reacdistbelownozzle * 2.54 / 100 + nozzconelen + nozzcylh + nozzholdh
    nozzholdboxstart = reacdistbelownozzle * 2.54 / 100 + nozzconelen + nozzcylh + nozzboxstart
    nozzholdboxend = reacdistbelownozzle * 2.54 / 100 + nozzconelen + nozzcylh + nozzboxstart + nozzboxh

    nozzboxhalfl = float(nozzparms[10]) / 2
    nozzboxhalfw = float(nozzparms[11]) / 2

    conedia = float(coneparms[3])  # cone outer diameter in inches
    coneheight = float(coneparms[1])  # cone height in inches as measured from the top of the ISO-100 cylinder.

    # Height above the cone that the reaction occurs
    reacheight = ((nozzleconedistin - reacdistbelownozzle) * 2.54) / 100

    conefitdist = reacheight + coneheight * 2.54 / 100

    # In the new SOLSTISE designs, the ISO base isn't really present, so we'll make this the size of the bottom
    # of the cone instead.
    basedia = 4.2  # inches
    riso100cyl = (basedia / 2) * 2.54 / 100

    # distance from the reaction that the cone side equation starts (this equation starts at ~5.5, now 3.83).
    # sideheight = (5.5 - reacdistbelownozzle) * 2.54 / 100
    sideheight = (nozzleconedistin - reacdistbelownozzle) * 2.54 / 100

    # Polynomial coefficients for the current best cone: Top of cone is 3.83 in away from nozzle

    # poly3 = [0.0611, -0.8077, 3.5255, -3.7867]
    poly3 = [float(coneparms[5]), float(coneparms[6]), float(coneparms[7]), float(coneparms[8])]

    rconeside = lambda y: (poly3[0] * (y + reacdistbelownozzle - nozzleconedistin) ** 3 +
                           poly3[1] * (y + reacdistbelownozzle - nozzleconedistin) ** 2 +
                           poly3[2] * (y + reacdistbelownozzle - nozzleconedistin) + poly3[3]) * 2.54 / 100

    # parameters for nozzle shadowing here:
    nozzang = float(nozzparms[2]) * np.pi / 180  # the nozzle slope angle, which is 22 degrees (converted to rads)
    nozzdia = float(nozzparms[1])

    # feed it inches, gives back meters. Calc the nozzle radius at a specified y position height.
    rnozzle = lambda y: (nozzdia / 2 + (y - reacdistbelownozzle) * np.tan(nozzang)) * 2.54 / 100

    # function determines the r coordinates of the 2nd circle that makes up the gas pipe.
    # Need pipepm to multiply rpipe by a minus sign if we are using an ISO 160 Pipe:
    if rblock < rbore:
        pipepm = -1
    else:
        pipepm = 1

    rpipe = lambda ph: cheight*np.sin(ph) + pipepm * np.sqrt(cheight**2*np.sin(ph)**2 - cheight**2 + rblock**2)

    # Distance from the reaction to the top of the vertical ISO160 Pipe in m
    iso160pipeheight = (float(coneparms[2]) - reacdistbelownozzle) * 2.54 / 100

    # ***************************************************************************************
    # ***************************************************************************************
    # ***************************************************************************************
    # ***************************************************************************************

    # dummy = df['Energy']

    # The last position of the particle, to determine how far the particle travelled in the last time step
    xlast = np.zeros_like(phic)
    ylast = np.zeros_like(phic)
    zlast = np.zeros_like(phic)

    if elossbool:
        if gas:
            # Set the jet radius here, gets set from targetparms. Redefinition could be removed:
            jetr = jetrad / 100

            # zoff is now the z offset in m.
            zoff = df['z_Offset'].to_numpy() / 100.0 - jetr

            jetroff = zoff
            jetroff = np.where(zoff < -jetr, -jetr, jetroff)
            jetroff = np.where(zoff > jetr, jetr, jetroff)
        else:
            zoff = np.zeros_like(df['Energy'].to_numpy())

    else:
        zoff = np.random.normal(np.zeros_like(df['Energy'].to_numpy()), ((jetrin * 2)/1000)/2.355)

    beamspotdia = beamdiamm / 1000  # beamdiamm mm FWHM
    xoff = np.random.normal(np.zeros_like(df['Energy'].to_numpy()), beamspotdia / 2.355)
    yoff = np.random.normal(np.zeros_like(df['Energy'].to_numpy()), beamspotdia / 2.355)

    # For the ejectile energy loss, the energy loss will be split into two layers (if gas target).
    # The first layer will be the jetlength, defined below which is the approximate straight-line distance the
    # ejectile will traverse the jet, and the total distance (disttravl) that the particle will traverse in its orbit.
    # I make the assumption that the energy loss and angular spread will be tiny and don't take it into account here.
    disttravl = np.zeros_like(phic)

    helpfultip = ["You probably don\'t have enough time to make a cup of coffee…",
                  "If you're getting up for a bathroom break, you\'d better be quick!",
                  "Jeopardy jingle: do do do do do do do, do do do do doot! do do do do do",
                  "The Tennessee jumping snake can grow up to 5 feet in length and jump up to 7 feet high!",
                  "Close your eyes and pretend that you\'re not waiting for this code to run.",
                  "If your office has windows, look outside! If your office doesn\'t, pretend to.",
                  "The coffee cart beckons... later.",
                  "Honk if this code is taking too long to run!",
                  "Hold on to your butts...",
                  "An actual tip: You can still draw each histogram in Plotter, even if it isn't listed "
                  "(though it might not make sense...).",
                  "An actual tip: Each section of the code can be run individually Ex: >>python3 Plotter.py",
                  "An actual tip: If the code keeps crashing, you may need to upgrade your version of pandas or numpy."]

    htipnum = randrange(12)

    i_final = np.zeros_like(df['Energy'])

    # Simulating events status bar for the for loop
    print("Simulating Events...")
    statbar = "[                                   ]"

    vperp = df['vel_perp'].to_numpy()
    phii = df['Phi'].to_numpy()
    vpar = df['vel_par'].to_numpy()
    tred = df['t_reduced'].to_numpy()

    # Splits the flight time into 300 segments for tracking purposes to see whether or not the particle is blocked.
    for i in range(350):

        if i == 100:
            print(helpfultip[htipnum], end='\n')

        if i % 10 == 0:
            statbar = statbar.replace(" ", "=", 1)
            print(statbar, end='\r', flush=True)

        # Loop goes to 350 because we want to make sure all of the particles have hit the detector. With the beam spot
        # size offset, it doesn't necessarily happen at t_reduced (though most particles do).
        t = tred/300 * (i+1)

        # If we aren't doing the energy loss I want the code to function like normal.

        xpos = ((vperp/omega)*np.sin((omega*t) - cwsign * phii)) - ((vperp/omega)*np.sin(-1 * cwsign * phii)) + xoff
        ypos = (cwsign * (vperp/omega)*np.cos(omega*t - cwsign * phii)) - \
               cwsign * vperp/omega*np.cos(-1 * cwsign * phii) + yoff
        zpos = vpar*t + zoff

        # For debugging:
        # print(str(xpos[0]) + " " + str(ypos[0]))

        if (i+1) % 10 == 0 and elossbool and i < 300:
            disttravl = np.sqrt((xlast - xpos) ** 2 + (ylast - ypos) ** 2 + (zlast - zpos) ** 2) + disttravl

        # r is the radial position of the particle
        r = np.sqrt(xpos**2 + ypos**2)

        # Grab the iteration where the radius is less than the detector radius (i.e. it hits the detector)
        if i > 295:
            i_final = np.where((r < r0) & (i_final == 0), i, i_final)

        # phic is the phi current position of the particle. Before I had phic = np.arctan2(ypos, xpos) + np.pi which is
        # wrong. We actually need to add 360 deg if y<0 and add 0 if y > 0
        phic = np.where(ypos < 0, np.arctan2(ypos, xpos) + 2 * np.pi, np.arctan2(ypos, xpos))

        # rxzplane determines the radial position in the xz-plane
        rxzplane = np.sqrt(xpos**2 + zpos**2)

        # rpipe determines the r position of the 2nd circle boundary
        # so if the particle radius is greater than that, it gets blocked
        maskrpipe = (r > rpipe(phic))

        # Shadowing from the vertical portion of the ISO160 pipe. Radius is 3.25 inches.
        masktoppipe = (rxzplane < 0.08255) & (-1*ypos > iso160pipeheight)

        # maskphipipe is the mask that determines whether or not the particle is within the phi boundaries of the pipe
        # if maskphipipe and maskrpipe are true, then the particle is blocked by the pipe
        maskphipipe = (phic > phi1block) & (phic < phi2block)

        maskpipe = (maskrpipe & maskphipipe) | masktoppipe

        # maskcone determines if the particle is within the opening of the cone
        # Since if statements on the arrays are so slow, we'll break up the cone mask into three: tube (top), sides,
        # and ISO base. In the future, the base could be gotten rid of depending on final geometry.
        # masktop = (rxzplane < rcone) & (ypos > reacheight) & (ypos < sideheight)
        # With the new cone, masktop does not need to be used anymore since the long straight neck at the top of
        # the cone has been removed.
        masksides = (rxzplane < rconeside(-1 * ypos * 100 / 2.54)) & ((-1 * ypos) > sideheight) & \
                    ((-1 * ypos) < conefitdist)
        # masktest = (rxzplane < rconeside(-1 * ypos * 100 / 2.54)) & ((-1 * ypos) > sideheight)
        maskbase = (rxzplane < riso100cyl) & ((-1 * ypos) > conefitdist) & ((-1 * ypos) < iso160pipeheight)

        maskcone = masksides | maskbase

        # masknozzle determines if the the particle hits the nozzle.
        # masknozzle = masknozzle & ((rxzplane < rnozzle(-1 * ypos * 100 / 2.54)) & ((-1 * ypos) > reacdistbelownozzle *
        #                                                                         2.54 / 100))
        # if the ypos is between the nozzle cone exaust and top of nozzle cone:
        masknozzlecone = (rxzplane < rnozzle(ypos * 100 / 2.54)) & \
                         (reacdistbelownozzle < (ypos * 100 / 2.54)) & \
                         ((ypos * 100 / 2.54) < (reacdistbelownozzle + nozzconelen * 100 / 2.54))

        masknozzlecyl = ((ypos * 100 / 2.54) > (reacdistbelownozzle + nozzconelen * 100 / 2.54)) & \
                        (rxzplane < nozzcylrad) & \
                        (ypos < (reacdistbelownozzle * 2.54 / 100 + nozzconelen + nozzcylh))

        masknozzleholdercyl = (ypos > nozzholdcylstart) & (rxzplane < nozzholdrad) & (ypos < nozzholdcylend)

        masknozzleholderbox = (ypos > nozzholdboxstart) & (ypos < nozzholdboxend) & (np.abs(xpos) < nozzboxhalfl) & \
                              (np.abs(zpos) < nozzboxhalfw)

        masknozzle = masknozzlecone | masknozzlecyl | masknozzleholdercyl | masknozzleholderbox

        #print(np.where((ypos > 0.0025) & (ypos < 0.0127), rnozzle(ypos*100/2.54), np.zeros_like(ypos)))

        # Combine all the different masks into maskmasters down here. Because of this section,
        # doing maskcone = maskcone & ... is unnecessary above because it is done here.
        # Here they get inverted so particles that hit are False, and particles that don't are True
        maskmaster_cone = maskmaster_cone & np.invert(maskcone)
        maskmaster_pipe = maskmaster_pipe & np.invert(maskpipe)
        maskmaster_nozzle = maskmaster_nozzle & np.invert(masknozzle)

        # Here, I'll attempt to remove similar events from the masks (i.e. ones that would hit the cone and pipe, etc)
        # Basically, if it already hit one of the other two blocking elements, we want to update the maskmaster to
        # prevent overlap.
        maskmaster_cone = np.where(np.invert(maskmaster_nozzle) | np.invert(maskmaster_pipe), True, maskmaster_cone)
        maskmaster_nozzle = np.where(np.invert(maskmaster_cone) | np.invert(maskmaster_pipe), True, maskmaster_nozzle)
        maskmaster_pipe = np.where(np.invert(maskmaster_nozzle) | np.invert(maskmaster_cone), True, maskmaster_pipe)

        #maskmaster = maskmaster*np.invert(maskcone | masknozzle | maskrpipe | maskphipipe)
        maskrbore = maskrbore & (r < rbore)

        if (i+1) % 10 == 0 and elossbool and i < 300:
            xlast = xpos
            ylast = ypos
            zlast = zpos

    # Calculate the final z and phi of the particles based on the iteration where they hit the detector array
    t_final = tred/300 * (i_final + 1)
    xpos = ((vperp / omega) * np.sin((omega * t_final) - cwsign * phii)) - (
                (vperp / omega) * np.sin(-1 * cwsign * phii)) + xoff
    ypos = (cwsign * (vperp / omega) * np.cos(omega * t_final - cwsign * phii)) - \
           cwsign * vperp / omega * np.cos(-1 * cwsign * phii) + yoff
    phic = np.where(ypos < 0, np.arctan2(ypos, xpos) + 2 * np.pi, np.arctan2(ypos, xpos))
    zpos = vpar*t_final + zoff

    # Move maskmaster from line 367 and maskz out of the loop
    # we want only particles that come out at backward angles for inverse, and forward angles for normal kin
    if invkin:
        maskz = zpos < 0
    else:
        maskz = zpos > 0

    maskmaster = maskmaster_cone & maskmaster_pipe & maskmaster_nozzle

    print("\n")

    # Grab the max energy before the eloss section because sometimes we get one really large value that messes up the
    # histogram ranges in plotter...
    maxe = df['Energy'].max()

    if elossbool:
        # Need to set up the projectile data here that goes into desorb:
        zp = np.zeros_like(phic) + zeject
        ap = np.zeros_like(phic) + aeject
        proj_e = df['Energy'].to_numpy()

        emaxinit = df['Energy'].max()

        # Make an empty data frame to store the output:
        df_elossout = pd.DataFrame()
        proj_ein = proj_e

        if gas:
            # The distance that the particle must traverse to get out of the jet is set here.
            jetlength = np.abs(
                jetr / np.sqrt(np.sin(df['Theta_Rad'].to_numpy()) ** 2 * np.cos(df['Phi'].to_numpy()) ** 2 +
                               np.cos(df['Theta_Rad'].to_numpy()) ** 2) + jetroff /
                np.sqrt(np.sin(df['Theta_Rad'].to_numpy()) ** 2 * np.cos(df['Phi'].to_numpy()) ** 2 +
                        np.cos(df['Theta_Rad'].to_numpy()) ** 2))

            # convert jetlength to cm:
            jetlength = jetlength * 100.0

            chamlength = disttravl * 100.0 - jetlength

            for j in range(2):
                # j = 0 corresponds to the jet thickness, and 1 corresponds to the chamber
                if j == 0:
                    df_elossout = desorb(zp, ap, proj_ein, ztarg, atarg, numtarg, gas, 0, 0, jetpress, jetlength,
                                         proj_e)

                    proj_ein = df_elossout['Energy_i'].to_numpy() - df_elossout['DeltaE_tot'].to_numpy()
                if j == 1:
                    df_elossout = desorb(zp, ap, proj_ein, ztarg, atarg, numtarg, gas, 0, 0, champress, chamlength,
                                         proj_e)

                    proj_e = df_elossout['Energy_i'].to_numpy() - df_elossout['DeltaE_tot'].to_numpy()
                    estragtot = df_elossout['E_strag_FWHM'].to_numpy()
        if not gas:
            print("\nYou're using a solid target, so the energy loss calculation is going to take a minute or two...")

            if invkin:
                # If in inverse kinematics, we want the target thickness to be the thickness traversed by the beam,
                # which is what we get from event builder.
                indthickness = df['Tgt_Thick'].to_numpy() / np.sin(df['Theta_Rad'].to_numpy() - np.pi/2)
            else:
                # If in normal kinematics, we want to subtract the thickness seen by the beam from the target thickness
                # So, if the beam sees 0.95 mg/cm^2 of a 1 mg/cm^2 target, the light ejectile will see 0.05 mg/cm^2
                indthickness = (thickness - df['Tgt_Thick'].to_numpy()) / np.sin(df['Theta_Rad'].to_numpy())

            # We don't need a for loop here because there's only one layer for the protons to lose energy
            df_elossout = desorb(zp, ap, proj_ein, ztarg, atarg, numtarg, gas, density, indthickness, 0, 0, proj_e)
            proj_e = df_elossout['Energy_i'].to_numpy() - df_elossout['DeltaE_tot'].to_numpy()
            estragtot = df_elossout['E_strag_FWHM'].to_numpy()

        emask = emaxinit > proj_e
        df['Energy'] = np.random.normal(proj_e, estragtot)

        # Detector energy resolution assumed to be 25 keV, divide by 2.355 to get sigma
        df['Energy'] = np.random.normal(df['Energy'], 0.025 / 2.355)

        # If we have a solid target, the thickness in mm is way too small to make a difference, so we ignore it.

        # Detector position resolution depends on particle energy, so it'll be more difficult to put in
        e0to2 = df['Energy'] < 2
        e2to4 = (df['Energy'] > 2) & (df['Energy'] < 4)
        e4to6 = (df['Energy'] > 4) & (df['Energy'] < 6)
        egt6 = df['Energy'] > 6

        # Mask the energies to include the detector position resolution:
        zpos = np.where(e0to2, np.random.normal(zpos, 0.00117 / 2.355), zpos)
        zpos = np.where(e2to4, np.random.normal(zpos, 0.00085 / 2.355), zpos)
        zpos = np.where(e4to6, np.random.normal(zpos, 0.000532 / 2.355), zpos)
        zpos = np.where(egt6, np.random.normal(zpos, 0.0004 / 2.355), zpos)


    # Adds the final phi position to the dataframe
    df['zpos_final'] = zpos
    # Adds the final phi position to the dataframe
    df['Phi_final'] = phic

    # splits the particles up onto four quadrants of a fictional detector
    df["Det1"] = (phic > 0) & (phic < np.pi/2)
    df["Det2"] = (phic > np.pi/2) & (phic < np.pi)
    df["Det3"] = (phic > np.pi) & (phic < 3*np.pi/2)
    df["Det4"] = (phic > 3*np.pi/2) & (phic < 2*np.pi)

    # mask for individual detectors here. Need det_zpos_i which gives the start of the detector array.
    # Dets are 5 mm wide with 1 cm gaps
    det_zpos_i = detzi
    df['Detz1'] = (np.abs(zpos) > np.abs(det_zpos_i)) & (np.abs(zpos) < (np.abs(det_zpos_i) + 0.05))
    df['Detz2'] = (np.abs(zpos) > (np.abs(det_zpos_i) + 0.06)) & (np.abs(zpos) < (np.abs(det_zpos_i) + 0.11))
    df['Detz3'] = (np.abs(zpos) > (np.abs(det_zpos_i) + 0.12)) & (np.abs(zpos) < (np.abs(det_zpos_i) + 0.17))
    df['Detz4'] = (np.abs(zpos) > (np.abs(det_zpos_i) + 0.18)) & (np.abs(zpos) < (np.abs(det_zpos_i) + 0.23))
    df['Detz5'] = (np.abs(zpos) > (np.abs(det_zpos_i) + 0.24)) & (np.abs(zpos) < (np.abs(det_zpos_i) + 0.29))
    df['Detz6'] = (np.abs(zpos) > (np.abs(det_zpos_i) + 0.3)) & (np.abs(zpos) < (np.abs(det_zpos_i) + 0.35))

    if invkin:
        masktheta = df["Theta_Deg"] > 95
    else:
        masktheta = df["Theta_Deg"] < 85

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

    df["UnblockedSolidTarg"] = maskz & maskrbore & masktheta

    # We'll also reconstruct the Q-value spectrum from the "detected" particles.
    # First we have to calculate the CM Energy:
    if not elossbool:
        df['EnergyCM'] = df['Energy'] + .5 * me * utoMeV * (vcm / c) ** 2 - me * utoMeV * (vcm / c) / t_final * \
                         df['zpos_final'] / c
        df['Ex_Reconstructed'] = tcm + qvalnoex - df['EnergyCM'] * (me + mr) / mr
    else:
        # Calculate the Ex spectrum differently for energy loss.
        # The first 3 events are hard coded in to be "average" events
        # Rotate the E vs z spectrum:
        slope = (df['Energy'][1] - df['Energy'][0]) / (zpos[1] - zpos[0])

        icept = df['Energy'][0] - slope * zpos[0]

        exi = df['Energy'] - slope * zpos

        exslope = (ex1 - ex2) / (exi[0] - exi[2])
        exicept = ex1 - exslope * exi[0]

        df['Ex_Reconstructed'] = exslope * exi + exicept

    if rblock < rbore:
        custpipe = False
    else:
        custpipe = True

    dictparams = {
        "Reaction": reac,
        "Beam Energy": ebeam,
        "Magnetic Field": B,
        "Reaction Distance from Nozzle": reacdistbelownozzle,
        "Nozzle-Cone Distance": nozzleconedistin,
        "Bore Radius": rbore,
        "Custom Pipe?": custpipe,
        "Pipe Radius": rblock,
        "Pipe Left Edge Angle": int(phi1block*180/np.pi),
        "Pipe Right Edge Angle": int(phi2block * 180 / np.pi),
        "Cone Opening Diameter": conedia,
        "Cone Height": coneheight,
        "Calculated Energy Loss?": elossbool,
        "Cone File": conetxt,
        "Nozzle File": nozztxt,
        "Solid Thickness": thickness,
        "Jet Pressure": jetpress,
        "Chamber Pressure": champress,
        "Jet Radius": jetrad * 10,
        "Max Energy": maxe
        }

    dfparams = pd.DataFrame([dictparams])

    df_all = pd.concat([df, dfparams], axis=1)

    writefilebase = filein[:-4]

    writefile = writefilebase + str(1)
    lnum = 1

    pklstring = writefilebase + '*.pkl'

    list_pkls = glob.glob(outdir + pklstring)

    # The point of the following section is to change writefile if a file with the same data already exists
    # so we don't overwrite something we don't want to.

    if len(list_pkls) > 0:
        latest_pkl = max(list_pkls, key=os.path.getctime)
        latest_pkl = latest_pkl[15:]
        if not latest_pkl[-6].isdigit():
            lnum = int(latest_pkl[-5])
        else:
            # Gets lnum if the filename is in double digits.
            lnum = int(latest_pkl[-6:-4])

        fileyn = input("\n\nAn event file using this data already exists. Would you like to overwrite the file? [Y/N] ")
        if fileyn == "n" or fileyn == "N":
            print("\nA number will be appended onto the end of the file name.")
            while os.path.exists(outdir + writefile + ".pkl"):
                lnum = lnum + 1
                writefile = writefilebase + str(lnum)
        else:
            # Latest pkl already contains lnum
            writefile = latest_pkl[:-4]

    df_all.to_pickle(outdir + writefile + ".pkl")

    print("\nOutput file named " + writefile + ".pkl was generated.")

    plotyn = input("\nWould you like to plot the simulated data? [Y/N] ")
    if plotyn == "N" or plotyn == "n":
        print("\nNow exiting. You can plot the simulated data any time by running Plotter.py.")
    else:
        plot(outdir + writefile + ".pkl", "")

 #   input("\nPress ENTER to end.")





