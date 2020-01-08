import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle
from massreader import readmass
from mpl_toolkits import mplot3d


def sim(rbore, rblock, cheight, phi1block, phi2block, ebeam, filein, reac):

    # The z axis points in beam direction, the x-axis points to the left, and the y-axis points down

    masses = readmass(reac)

    # reaction of form t(b,e)R
    utoMeV = 931.4941

    mt = masses[0]
    mb = masses[1]
    me = masses[2]
    mr = masses[3]

    #mt = 2.0135532
    #mb = 27.97692653246
    #me = 1.0072765
    #mr = 28.9764947

    #ebeam = 168 # MeV, for d(28Si,p) it is 6 MeV/u

    rblock = np.array(rblock, dtype=np.float64)
    cheight = np.array(cheight, dtype=np.float64)
    phi1block = np.array(phi1block, dtype=np.float64)
    phi2block = np.array(phi2block, dtype=np.float64)

    mevtoj = 1.6021766e-13
    c = 2.998e8

    # finite radius of the detector array
    r0 = 0.011

    tcm = mt/(mb+mt)*ebeam

    mass = 1.6726219e-27  # proton mass in MeV/c^2
    B = 1.915  # teslas
    q = 1.6e-19  # 1 elemental charge in coulombs


    # Generates a numpy array of shape (xxx,2) whose columns are theta angle and energy.
    data = np.genfromtxt(filein)

    # read in the energy and angle from the file (all rows, 1st column), (all rows, 0th column)
    energy = data[:, 1]
    theta = data[:, 0]

    # Convert the angle to radians
    theta = theta * np.pi/180

    # Cyclotron frequency and period.
    omega = (q*B)/mass
    tcyc = (2*np.pi)/omega

    # Velocity of the ejectile in the lab frame (m/s)
    vel = np.sqrt((2 * energy * mevtoj) / mass)
    # Velocities parallel to the z-axis and perpendicular to the z-axis.
    vper = vel * np.sin(theta)
    vpar = vel * np.cos(theta)

    # Calculates Q-Ex from the energy and angle of the ejected particle (Krane 11.10)
    qval = energy * (1 + me / mr) - ebeam * (1 - mb / mr) - 2 * np.sqrt(
        mb * me / mr**2 * energy * ebeam) * np.cos(theta)

    # The total energy of the reaction, Tcm + Q-value - Ex
    etot = tcm + qval

    # The velocity of the ejectile in the lab frame.
    v0 = np.sqrt(2 * mr * utoMeV * etot / (me * utoMeV * (me * utoMeV + mr * utoMeV))) * c

    # The velocity of the CM frame, it's a constant and dependent on the beam energy and species.
    vcm = np.sqrt(2 * ebeam / (mb * utoMeV)) * (mb * utoMeV/(mb * utoMeV+mt * utoMeV)) * c

    # Calculates the CM angle of the ejectile and recoil. arccos takes values from -1 to 1 and won't break if it
    # gets a value outside of that range.
    thetacm = np.arccos((vel**2 - v0**2 - vcm**2)/(2 * v0 * vcm))

    # Could mask all of the arrays to get rid of "bad" values if needed, but it probably isn't worth it since
    # it doesn't actually break the code.
    cosarg = (vel ** 2 - v0 ** 2 - vcm ** 2) / (2 * v0 * vcm)
    maskarg = cosarg >= -1

    treduced = tcyc - r0/(v0*np.sin(thetacm))

    energy = energy[maskarg]
    theta = theta[maskarg]
    treduced = treduced[maskarg]
    vper = vper[maskarg]
    vpar = vpar[maskarg]

    # makes a phi array the same size as the theta array, random number 0 to 1
    phi = np.random.rand(theta.size)
    # then multiply the phi array by 2pi to get a real phi value
    phi = phi * 2 * np.pi

    # debugging
    # print(energy.shape)
    # print(phi.shape)

    # creates a mask the same shape as the energy array
    # maskmaster is the mask that keeps track of the overall mask in the loop
    maskmaster = energy > 0
    # maskmasters for the pipe and cone defined here.
    maskmaster_pipe = energy > 0
    maskmaster_cone = energy > 0
    maskmaster_nozzle = energy > 0
    # maskrbore takes care of the mask for the bore radius because some particles will hit that
    maskrbore = energy > 0
    # maskcone initialized like maskrbore
    masknozzle = energy > 0
    # this just initializes phic, but it might be unnecessary.
    phic = energy

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

    #print(sideheight)
    #print(baseheight)

    # polynomial coefficients

    #poly3 = [0.0349, -0.4655, 1.9711, -1.4778] These were for the tall cone

    # Polynomial coefficients for the current best cone: Top of cone is 3.83 in away from nozzle

    poly3 = [0.0611, -0.8077, 3.5255, -3.7867]

    rconeside = lambda y: (poly3[0] * (y + reacdistbelownozzle)**3 + poly3[1] * (y + reacdistbelownozzle)**2 +
                           poly3[2] * (y + reacdistbelownozzle) + poly3[3]) * 2.54 / 100

    # parameters for nozzle shadowing here:
    # cone tip dist is the extra vertical height added onto the nozzle "cone" shape if it extended out to a point
    # it was determined using the 22 degree slope of the nozzle side and the nozzle opening radius 0.0455 inches
    nozzletipdist = 0.112616
    nozzleang = 22 * np.pi / 180  # the nozzle slope angle, which is 22 degrees (converted to rads)

    rnozzle = lambda y: ((y+np.abs(reacdistbelownozzle - nozzletipdist))*np.tan(nozzleang)) * 2.54 / 100

    #This worked fine, but it was SLOW! Also, using the lambda function above produces the same amount of shadowing
    #So at some point the actual geometry doesn't matter.
    #def rnozzle(y):
    #    if y < 0.6081:
    #        return ((y+np.abs(reacdistbelownozzle - nozzletipdist))*np.tan(nozzleang)) * 2.54 / 100
    #    else:
    #        return 0.25 * 2.54 / 100

    #vrnozzle = np.vectorize(rnozzle)

    #reacheight = 0.0356 #m, y direction height above the cone
    #rcone = 0.01524 #m cone radius

    # function determines the r coordinates of the 2nd circle that makes up the gas pipe.
    rpipe = lambda ph: cheight*np.sin(ph) + np.sqrt(cheight**2*np.sin(ph)**2 - cheight**2 + rblock**2)

    # Simulating events status bar for the for loop
    print("Simulating Events...")
    statbar = "[                              ]"

    # Splits the flight time into 300 segments for tracking purposes to see whether or not the particle is blocked.
    for i in range(300):

        if i % 10 == 0:
            statbar = statbar.replace(" ", "=", 1)
            print(statbar, end='\r', flush=True)

        t = treduced/300 * (i+1)
        xpos = (-(vper/omega)*np.cos((omega*t)+phi))+((vper/omega)*np.cos(phi))
        ypos = ((vper/omega)*np.sin(omega*t+phi))-vper/omega*np.sin(phi)
        zpos = vpar*t

        # r is the radial position of the particle
        r = np.sqrt(xpos**2 + ypos**2)
        # phic is the phi position of the particle.
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
        # and ISO base. In the future, the base could be gotten rid of.
        # masktop = (rxzplane < rcone) & (ypos > reacheight) & (ypos < sideheight)
        # With the new cone, masktop does not need to be used anymore since the long straight neck at the top of
        # the cone has been removed.
        masksides = (rxzplane < rconeside(ypos * 100 / 2.54)) & (ypos > (sideheight)) & (ypos < baseheight)
        maskbase = (rxzplane < rISObase) & (ypos > baseheight)

        maskcone = masksides | maskbase

        # we want only particles that come out at backward angles
        maskz = zpos < 0

        # masknozzle determines if the the particle hits the nozzle.
        masknozzle = masknozzle*((rxzplane < rnozzle(-1 * ypos * 100 / 2.54)) & ((-1 * ypos) > reacdistbelownozzle *
                                                                                 2.54 / 100))

        #print(rxzplane[masksides] * 100 / 2.54)
        #print(maskr.shape, phir.shape, zpos[np.invert(maskz)].shape)

        #maskmaster = maskmaster*np.invert(maskrpipe & maskphipipe)#np.invert(maskcone)*np.invert(masky)#*np.invert(maskz)
        maskmaster = maskmaster*np.invert(maskcone)*np.invert(maskrpipe & maskphipipe)*np.invert(masknozzle)
        maskmaster_cone = maskmaster_cone*np.invert(maskcone)
        maskmaster_pipe = maskmaster_pipe*np.invert(maskrpipe & maskphipipe)
        maskmaster_nozzle = maskmaster_nozzle * np.invert(masknozzle)
        #maskmaster = maskmaster*np.invert(maskcone | masknozzle | maskrpipe | maskphipipe)
        maskrbore = maskrbore & (r < rbore)

    # splits the particles up onto four quadrants of a fictional detector
    maskdet1 = (phic > 0) & (phic < np.pi/2)
    maskdet2 = (phic > np.pi/2) & (phic < np.pi)
    maskdet3 = (phic > np.pi) & (phic < 3*np.pi/2)
    maskdet4 = (phic > 3*np.pi/2) & (phic < 2*np.pi)
    # inverting maskmaster will pick out the blocked particles
    maskmasterdet = [maskmaster & maskdet1 & maskz & maskrbore, maskmaster & maskdet2 & maskz & maskrbore,
                     maskmaster & maskdet3 & maskz & maskrbore, maskmaster & maskdet4 & maskz & maskrbore]
    maskmasterdet_cone = [np.invert(maskmaster_cone) & maskdet1 & maskz & maskrbore,
                          np.invert(maskmaster_cone) & maskdet2 & maskz & maskrbore,
                          np.invert(maskmaster_cone) & maskdet3 & maskz & maskrbore,
                          np.invert(maskmaster_cone) & maskdet4 & maskz & maskrbore]
    maskmasterdet_pipe = [np.invert(maskmaster_pipe) & maskdet1 & maskz & maskrbore,
                          np.invert(maskmaster_pipe) & maskdet2 & maskz & maskrbore,
                          np.invert(maskmaster_pipe) & maskdet3 & maskz & maskrbore,
                          np.invert(maskmaster_pipe) & maskdet4 & maskz & maskrbore]
    maskmasterdet_nozzle = [np.invert(maskmaster_nozzle) & maskdet1 & maskz & maskrbore,
                          np.invert(maskmaster_nozzle) & maskdet2 & maskz & maskrbore,
                          np.invert(maskmaster_nozzle) & maskdet3 & maskz & maskrbore,
                          np.invert(maskmaster_nozzle) & maskdet4 & maskz & maskrbore]

    invmaskmasterdet1 = np.invert(maskmaster) & maskdet1 & maskz & maskrbore
    invmaskmasterdet2 = np.invert(maskmaster) & maskdet2 & maskz & maskrbore
    invmaskmasterdet3 = np.invert(maskmaster) & maskdet3 & maskz & maskrbore
    invmaskmasterdet4 = np.invert(maskmaster) & maskdet4 & maskz & maskrbore


    #print("I get here")
    thetablocked = theta[np.invert(maskmaster) & maskz & maskrbore]
    thetablocked = thetablocked*180/np.pi

    #print(thetablocked.shape)
    f2 = plt.figure(2)
    plt.rc('axes', labelsize=15)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    plt.hist(thetablocked)
    plt.xlabel('Lab Angle (deg)')
    plt.ylabel('Counts')
    f2.show()

    input("ENTER")

    #
    energyarr = [energy[maskmasterdet[0]], energy[maskmasterdet[1]], energy[maskmasterdet[2]], energy[maskmasterdet[3]]]
    thetaarr = [theta[maskmasterdet[0]], theta[maskmasterdet[1]], theta[maskmasterdet[2]], theta[maskmasterdet[3]]]
    zposarr = [zpos[maskmasterdet[0]], zpos[maskmasterdet[1]], zpos[maskmasterdet[2]], zpos[maskmasterdet[3]]]

    energyarr_blockedcone = [energy[maskmasterdet_cone[0]], energy[maskmasterdet_cone[1]],
                             energy[maskmasterdet_cone[2]], energy[maskmasterdet_cone[3]]]
    zposarr_blockedcone = [zpos[maskmasterdet_cone[0]], zpos[maskmasterdet_cone[1]],
                           zpos[maskmasterdet_cone[2]], zpos[maskmasterdet_cone[3]]]

    energyarr_blockedpipe = [energy[maskmasterdet_pipe[0]], energy[maskmasterdet_pipe[1]],
                             energy[maskmasterdet_pipe[2]], energy[maskmasterdet_pipe[3]]]
    zposarr_blockedpipe = [zpos[maskmasterdet_pipe[0]], zpos[maskmasterdet_pipe[1]],
                           zpos[maskmasterdet_pipe[2]], zpos[maskmasterdet_pipe[3]]]

    energyarr_blockednozzle = [energy[maskmasterdet_nozzle[0]], energy[maskmasterdet_nozzle[1]],
                             energy[maskmasterdet_nozzle[2]], energy[maskmasterdet_nozzle[3]]]
    zposarr_blockednozzle = [zpos[maskmasterdet_nozzle[0]], zpos[maskmasterdet_nozzle[1]],
                           zpos[maskmasterdet_nozzle[2]], zpos[maskmasterdet_nozzle[3]]]

    invenergyarr = [energy[invmaskmasterdet1], energy[invmaskmasterdet2], energy[invmaskmasterdet3], energy[invmaskmasterdet4]]
    invzposarr = [zpos[invmaskmasterdet1], zpos[invmaskmasterdet2], zpos[invmaskmasterdet3], zpos[invmaskmasterdet4]]

    #print(energy.shape)
    #print(theta.shape)

    #plt.hist2d(zpos, energy, bins=(1000,1000))
    #fig = plt.figure()
    #ax = plt.axes(projection='3d')
    #ax.scatter3D(xpos[:][0], ypos[:][0], zpos[:][0])

    print(invenergyarr[0])
    print(invenergyarr[0].shape)

    print(invzposarr[0])
    print(invzposarr[0].shape)


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


    fig = plt.figure()

    #for i in range(4):
    plt.rc('axes', labelsize=15)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)


    plt.subplot(2, 2, 1)
    p1 = plt.hist2d(zposarr[1], energyarr[1], bins=(750, 750), range=[[-0.8, 0], [0, 11]], cmap=newcmpBlue)
    p2 = plt.hist2d(zposarr_blockedcone[1], energyarr_blockedcone[1], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpGreen)
    p3 = plt.hist2d(zposarr_blockedpipe[1], energyarr_blockedpipe[1], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpRed)
    p4 = plt.hist2d(zposarr_blockednozzle[1], energyarr_blockednozzle[1], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpBlack)
    plt.xlabel('z(m)')
    plt.ylabel('Energy (MeV)')

    handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in [blu, grn, red, blk]]
    labels = ["Unblocked", "Cone", "Pipe", "Nozzle"]
    plt.legend(handles, labels, bbox_to_anchor=(1.4, 1.1), ncol=4)

    plt.subplot(2, 2, 2)
    plt.hist2d(zposarr[0], energyarr[0], bins=(750, 750), range=[[-0.8,0],[0,11]], cmap=newcmpBlue)
    plt.hist2d(zposarr_blockedcone[0], energyarr_blockedcone[0], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpGreen)
    plt.hist2d(zposarr_blockedpipe[0], energyarr_blockedpipe[0], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpRed)
    plt.hist2d(zposarr_blockednozzle[0], energyarr_blockednozzle[0], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpBlack)
    plt.xlabel('z(m)')
    plt.ylabel('Energy (MeV)')

    plt.subplot(2, 2, 3)
    plt.hist2d(zposarr[2], energyarr[2], bins=(750, 750), range=[[-0.8,0],[0,11]], cmap=newcmpBlue)
    plt.hist2d(zposarr_blockedcone[2], energyarr_blockedcone[2], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpGreen)
    plt.hist2d(zposarr_blockedpipe[2], energyarr_blockedpipe[2], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpRed)
    plt.hist2d(zposarr_blockednozzle[2], energyarr_blockednozzle[2], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpBlack)
    #plt.xlim(-.8, 0)
    #plt.ylim(0, 11)
    plt.xlabel('z(m)')
    plt.ylabel('Energy (MeV)')

    plt.subplot(2, 2, 4)
    plt.hist2d(zposarr[3], energyarr[3], bins=(750, 750), range=[[-0.8,0],[0,11]], cmap=newcmpBlue)
    plt.hist2d(zposarr_blockedcone[3], energyarr_blockedcone[3], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpGreen)
    plt.hist2d(zposarr_blockedpipe[3], energyarr_blockedpipe[3], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpRed)
    plt.hist2d(zposarr_blockednozzle[3], energyarr_blockednozzle[3], bins=(750, 750),
               range=[[-0.8, 0], [0, 11]], cmap=newcmpBlack)
    plt.xlabel('z(m)')
    plt.ylabel('Energy (MeV)')

    plt.show()


    #plt.hist2d(zpos, energy, bins=(500, 500), cmap=plt.cm.gist_earth_r)
    #plt.hist2d(energy, zpos, bins=(500, 500), cmap=plt.cm.gist_earth_r)
    #plt.show()

    #print(zpos)
    #print(ypos)
    #print(xpos)

    input("\nPress ENTER to end.")





