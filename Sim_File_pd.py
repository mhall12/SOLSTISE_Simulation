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

    print(mt,'\n')
    print(mb,'\n')
    print(me,"\n")
    print(mr,"\n")

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
    df = pd.read_csv(filein, sep="\t", header=None)
    df.columns = ["Theta_Deg", "Energy"]

    # Convert the angle to radians and add it to the dataframe
    df['Theta_Rad'] = df['Theta_Deg'] * np.pi/180

    # Cyclotron frequency and period.
    omega = (q * B) / (me*amutokg)
    tcyc = (2 * np.pi) / omega

    # Velocity of the ejectile in the lab frame (m/s) added to dataframe
    df['vel_ejec'] = np.sqrt((2 * df['Energy'] * mevtoj) / (me*amutokg))
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
    df['cosarg'] = (df['vel_ejec']**2 - df['v0']**2 - vcm**2)/(2 * df['v0'] * vcm)
    df['Theta_CM'] = np.arccos((df['vel_ejec']**2 - df['v0']**2 - vcm**2)/(2 * df['v0'] * vcm))

    print(df.iloc[:,-3:])
    print(df)

    input("\nPress ENTER to end.")





