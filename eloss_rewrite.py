import numpy as np
import pandas as pd


def desorb(z_projectile, a_projectile, energy, z_absorber, a_absorber, numa_absorber, isgas, density, thickness, pressure, length):

    # desorb functions close to what the old desorb does. However, it is rewritten to be readable...
    # nsand is the number of elements in the absorber.
    # z_, a_, en_ are the Z, A and energy of the projectile (the thing that is losing energy!)
    # isgas is a boolean True/False to determine whether or not the absorber is a gas or solid.
    # density is g/cm^2 and thickness is mg/cm^2.

    # z_absorber and a_absorber and numa_absorber are going to be arrays. So for example, if we had CD2,
    # z_absorber = [[6, 1]], a_absorber = [[12, 2]], numa_absorber = [[1,2]]. These are 2D arrays. Z, A are columns and
    # rows are the layers.
    # Density and thickness are 1D arrays where each element corresponds to a layer.

    # Input length, pressure only for gas. Input density, thickness only for non gasses

    # Here we'll set the partial thicknesses assuming it's a solid absorber
    thick_frac = []
    den_frac = []
    z_ind = []
    a_ind = []
    length_frac = []
    numa_ind = []
    gas = []
    for i in range(len(z_absorber)):
        a_tot = 0
        if not isgas[i]:
            for j in range(len(z_absorber[i])):
                a_tot = a_tot + a_absorber[i][j] * numa_absorber[i][j]
            for j in range(len(z_absorber[i])):
                a_frac = a_absorber[i][j] * numa_absorber[i][j] / a_tot
                thick_frac.append(thickness[i] * a_frac)
                den_frac.append(density[i] * a_frac)
                length_frac.append(thickness[i] / density[i] * 1000.0)
                z_ind.append(z_absorber[i][j])
                a_ind.append(a_absorber[i][j])
                numa_ind.append(numa_absorber[i][j])
                gas.append(False)
        if isgas[i]:
            for j in range(len(z_absorber[i])):
                a_tot = a_tot + a_absorber[i][j] * numa_absorber[i][j]
                thick_frac.append((pressure[i] / 760) * (length[i] / 22.4) * a_absorber[i] * numa_absorber[i])
                den_frac.append(thick_frac[i] / length[i])
                length_frac.append(length[i])
                gas.append(True)

    # We can set up the DataFrame down here:
    df_abs = pd.DataFrame()

    df_abs['Z_absorber'] = z_ind
    df_abs['A_absorber'] = a_ind
    df_abs['Number_Atoms'] = numa_ind
    df_abs['Partial_Density'] = den_frac
    df_abs['Partial_Thickness'] = thick_frac
    df_abs['Gas?'] = gas

    df_proj = pd.DataFrame()

    df_proj['Z_proj'] = z_projectile
    df_proj['A_proj'] = a_projectile
    df_proj['Energy_i'] = energy

    #print(df_proj)
    #print(df_abs)

    df_fin = eloss(df_abs, df_proj)

    return df_fin


def eloss(absorberdf, projectiledf):
    # initial number of integrations that's hard coded into the original desorb code.
    k = 0
    projectiledf['Energy_curr'] = projectiledf['Energy_i']
    projectiledf_dead = pd.DataFrame()
    dednext = np.zeros_like(projectiledf['Energy_i'])
    ded1st = np.zeros_like(dednext)
    eps = 0.0001
    valsdf = pd.DataFrame()
    zeros = np.zeros(len(projectiledf.index))
    valsdf['dednext'] = pd.Series(zeros)
    valsdf['ded1st'] = pd.Series(zeros)
    valsdf['ddd'] = pd.Series(zeros)
    valsdf['ddr'] = pd.Series(zeros)
    valsdf['dds'] = pd.Series(zeros)

    valsdf['k'] = pd.Series(zeros)
    valsdf['num_int'] = pd.Series(zeros) + 2
    num_integrations = valsdf['num_int'].max()
    k_min = valsdf['k'].min()

    while k_min < num_integrations:
        valsdf['k'] = valsdf['k'] + 1
        k_min = valsdf['k'].min()
        j_end = len(absorberdf.index)
        j = 0

        #absorberdf['FX'] = absorberdf['Partial_Thickness'] / num_integrations
        while j < j_end:
            # Loop through each row of the absorber dataframe by putting each row into df_currabsorber and using that
            df_currabsorber = absorberdf.iloc[j]
            valsdf['FX'] = df_currabsorber['Partial_Thickness'] / valsdf['num_int']
            # here we loop through each element, thickness etc to calculate how the energy changes after going through
            # each partial layer.
            projectiledf['Velocity'] = np.sqrt(2.13e-3 * projectiledf['Energy_curr'] / projectiledf['A_proj'])

            # Projectile velocity will change every loop step depending on the energy.
            # need to change Energy_curr in the dataframe.
            #print("Enter")
            #print(projectiledf)
            projectiledf['deltaE'] = dedx(df_currabsorber, projectiledf)

            # sign is always -1, so just hard code it in
            projectiledf['Energy_curr'] = projectiledf['Energy_curr'] + \
                                          projectiledf['deltaE'] * -1 * valsdf['FX']

            projectiledf['Egt0'] = projectiledf['Energy_curr'] > 0
            #print(projectiledf['Egt0'])

            #print(projectiledf)

            # Now we want to remove any values where the energy is 0 so the energy loss doesn't get calculated again

            #print(projectiledf_dead)

            # Collect all the particles that have 0 energy here:
            projectiledf_dead = projectiledf_dead.append(projectiledf[~projectiledf['Egt0']], sort=True)
            valsdf = valsdf[projectiledf['Egt0']]
            projectiledf = projectiledf[projectiledf['Egt0']]

            #print(projectiledf_dead)

            if k <= 2:
                valsdf['dednext'] = valsdf['dednext'] + projectiledf['deltaE'] * valsdf['FX']

            j = j + 1

            #if k < 50:
             #   print(k)

        k1mask = valsdf['k'] == 1

        valsdf.loc[k1mask, 'ded1st'] = valsdf['dednext']
        valsdf.loc[k1mask, 'dednext'] = pd.Series(zeros)

        k2mask = valsdf['k'] == 2
        valsdf.loc[k2mask, 'ddd'] = valsdf['ded1st'] - valsdf['dednext']

        dddmask = valsdf['ddd'] < 0
        k2dddmask = k2mask & dddmask

        valsdf.loc[k2dddmask, 'ddd'] = valsdf['ddd'] * -1

        valsdf.loc[k2mask, 'dds'] = valsdf['ded1st'] + valsdf['dednext']
        valsdf.loc[k2mask, 'ddr'] = valsdf['ddd'] / valsdf['dds']

        ddr_mask = valsdf['ddr'] > eps
        k2ddrmask = k2mask & ddr_mask

        if len(valsdf[k2ddrmask]) > 0:
            #projectiledf_dead = projectiledf_dead.append(projectiledf[np.invert(k2ddrmask)])
            #projectiledf = projectiledf[k2ddrmask]

            valsdf.loc[k2ddrmask, 'num_int'] = valsdf['num_int'] * 2
            j = -1
            valsdf.loc[k2ddrmask, 'k'] = 0 * valsdf['k']
            valsdf.loc[k2ddrmask, 'dednext'] = pd.Series(zeros)
            projectiledf.loc[k2ddrmask, 'Energy_curr'] = projectiledf['Energy_i']

        #print(valsdf['k'], valsdf['num_int'])

        num_integrations = valsdf['num_int'].max()
        k_min = valsdf['k'].min()

        projectiledf['keqmax_mask'] = valsdf['k'] == valsdf['num_int']
        valsdf = valsdf[~projectiledf['keqmax_mask']]
        projectiledf_dead = projectiledf_dead.append(projectiledf[projectiledf['keqmax_mask']], sort=True)
        projectiledf = projectiledf[~projectiledf['keqmax_mask']]
        #print(projectiledf, valsdf['k'], valsdf['num_int'])
        if projectiledf.empty:
            break


       # if k == 2:
       #     valsdf['ddd'] = valsdf['ded1st'] - valsdf['dednext']

       #     dddmask = valsdf['ddd'] < 0
        #    valsdf.loc[dddmask, 'ddd'] = valsdf['ddd'] * -1

        #    valsdf['dds'] = valsdf['ded1st'] + valsdf['dednext']
        #    valsdf['ddr'] = valsdf['ddd'] / valsdf['dds']



         #   ddr_mask = valsdf['ddr'] > eps

         #   if len(valsdf[ddr_mask]) > 0:
         #       projectiledf = projectiledf[ddr_mask]
         #       projectiledf_dead = projectiledf_dead.append(projectiledf[np.invert(ddr_mask)])

         #       num_integrations = num_integrations * 2
         #       j = -1
         #       k = 0
         #       valsdf['dednext'] = pd.Series(zeros)
         #       projectiledf['Energy_curr'] = projectiledf['Energy_i']

                #print(projectiledf['Energy_curr'])

    #print(projectiledf_dead)

    # Remake the original dataframe to include all of the rows that were removed previously.
    projectiledf = (projectiledf.append(projectiledf_dead, sort=True)).sort_index(ascending=True)
    elt0mask = projectiledf['Energy_curr'] < 0
    projectiledf.loc[elt0mask, 'Energy_curr'] = 0

    projectiledf['DeltaE_tot'] = projectiledf['Energy_i'] - projectiledf['Energy_curr']

    return projectiledf

def dedx(df_currabsorber, projectiledf):
    # Function calculates the differential energy loss dE/dX in solid targets using a semiempirical formula deduced
    # from experimental work

    # This program is modified for gas absorbers
    # H(Z2) is the sum of 5 guassian functions
    # A1 is the Mass Number - Projectile
    # Z2 is the Atomic number of the absorber
    # A1 is the Mass number of the absorber
    # RHO is the density of the absorber in grams/cm**3 (meanless if gas absorber)
    # ENER is the energy of the projectile in MeV
    # v is the velocity of the projectile in MeV/(mg/cm**2)
    # Z1 is atomic number - projectile


    calcdf = pd.DataFrame()

    rho = 0

    # Set the density here. Not sure how it gets the density if it's a gas, but we'll see later maybe
    if not df_currabsorber['Gas?']:
        rho = df_currabsorber['Partial_Density']
    elif df_currabsorber['Gas?']:
        rho = 1

    calcdf['xi'] = projectiledf['Velocity']**2 / df_currabsorber['Z_absorber']

    # Absorber function
    # G(XI) = Y(EXP) - Y(Theory) is deduced from experimental energy loss measurements
    fy = 0
    # fy is function y
    if not df_currabsorber['Gas?']:
        fy = 54721.0 * (1.0 + 5.15e-2 * np.sqrt(df_currabsorber['A_absorber'] / rho) -
                        np.exp(-0.23 * df_currabsorber['Z_absorber']))
    elif df_currabsorber['Gas?']:
        fy = 54721.0 * (1.35 - np.exp(df_currabsorber['Z_absorber'] * (-0.13 + 0.0014 * df_currabsorber['Z_absorber'])))

    # G(XI) is the derivation of a gaussian with variable height H(Z2)

    g1 = 0
    if df_currabsorber['Z_absorber'] <= 26.0:
        g1 = 19.84 * np.exp(-0.17 * (df_currabsorber['Z_absorber'] - 4.25) * (df_currabsorber['Z_absorber'] - 4.25))
    elif df_currabsorber['Z_absorber'] > 26.0:
        g1 = 0.000001

    g2 = 0
    if df_currabsorber['Z_absorber'] < 38.0:
        g2 = 17.12 * np.exp(-0.12 * (df_currabsorber['Z_absorber'] - 11.63) * (df_currabsorber['Z_absorber'] - 11.63))
    elif df_currabsorber['Z_absorber'] > 38.0:
        g2 = 0.0000001

    g3 = 7.95 * np.exp(-0.015 * (df_currabsorber['Z_absorber'] - 30.2) * (df_currabsorber['Z_absorber'] - 30.2))
    g4 = 5.84 * np.exp(-0.022 * (df_currabsorber['Z_absorber'] - 48.63) * (df_currabsorber['Z_absorber'] - 48.63))
    g5 = 7.27 * np.exp(-0.005 * (df_currabsorber['Z_absorber'] - 73.06) * (df_currabsorber['Z_absorber'] - 73.06))
    hz2 = (9.0 - (g1 + g2 + g3 + g4 + g5)) * 1.32e-5

    z2zwd = np.cbrt(df_currabsorber['Z_absorber']) * np.cbrt(df_currabsorber['Z_absorber'])

    # Multiplication factors of G(XI)
    fg = 1.2e-4 * df_currabsorber['Z_absorber'] * df_currabsorber['Z_absorber'] + \
         (2.49e-2 * df_currabsorber['A_absorber'] / rho)

    if df_currabsorber['Gas?']:
        fg = 1.3 / (1.0 + np.exp(3.0 - (df_currabsorber['Z_absorber'] / 5.0)))

    alefg = np.log(2.7e-5 / fg)

    # Calculation of G(XI)

    calcdf['gxi'] = calcdf['xi'] * 0.0

    xi_mask = (1.0e-9 <= calcdf['xi']) & (calcdf['xi'] <= 5.0e-4)

    calcdf_masked = calcdf[xi_mask]
    calcdf = calcdf[~xi_mask]

    #print(calcdf,calcdf_masked)

    if xi_mask[xi_mask].size > 0:
        sqxi = np.sqrt(calcdf_masked['xi'])
        c = 2.0 / df_currabsorber['Z_absorber'] * (sqxi / (1.0 + 1.0e4 * sqxi))

        if df_currabsorber['Gas?']:
            c = c / 2.0

        fg0 = 1.0 / (1.0 + (calcdf_masked['xi'] * 10000.0) * (calcdf_masked['xi'] * 10000.0) *
                     (calcdf_masked['xi'] * 10000.0))
        al = np.log(calcdf_masked['xi']) - alefg
        calcdf_masked['gxi'] = (c - hz2 * al * np.exp(-0.32 * al * al)) * fg0

    calcdf = (calcdf.append(calcdf_masked)).sort_index(ascending=True)

    # Calculation of Y(XI)
    calcdf['y'] = 3.3e-4 * np.log(1.0 + (calcdf['xi'] * fy)) + calcdf['gxi']

    # Energy loss of heavy ions
    # Effective charge

    calcdf['vv0'] = projectiledf['Velocity'] * 137.0

    calcdf['fv'] = calcdf['xi'] * 0.0 + 1.0
    #print(calcdf['fv'])

    velmask = projectiledf['Velocity'] >= 0.62

    calcdf.loc[velmask, 'fv'] = 1.0 - np.exp(-calcdf['vv0'])

    calcdf['az1'] = np.log(1.035 - 0.4 * np.exp(-0.16 * projectiledf['Z_proj']))

    calcdf['qq'] = projectiledf['Velocity'] / projectiledf['Z_proj']**0.509

    calcdf['ghi'] = projectiledf['Z_proj']

    calcdf['vz1'] = (-116.79 - 3350.4 * calcdf['qq']) * calcdf['qq']

    vz1mask = calcdf['vz1'] > -85.2

    calcdf.loc[vz1mask, 'ghi'] = projectiledf['Z_proj'] * (1.0 - np.exp(calcdf['vz1']))

    zprojmask = projectiledf['Z_proj'] > 2.0
    calcdf.loc[zprojmask, 'ghi'] = projectiledf['Z_proj'] * (1.0 - np.exp(calcdf['fv'] * calcdf['az1'] - 0.879 *
                                                                          (calcdf['vv0'] /
                                                                           projectiledf['Z_proj']**0.65)))

    # Effective charge of protons and aphas
    # Electronic energy loss DEDXHI

    calcdf['dedxhi'] = calcdf['ghi']**2 * df_currabsorber['Z_absorber'] * calcdf['y'] / \
                       (df_currabsorber['A_absorber'] * projectiledf['Velocity']**2)

    # nuclear energy loss DEDXNU
    calcdf['za'] = np.sqrt(np.cbrt(projectiledf['Z_proj'])**2 + z2zwd)

    calcdf['eps'] = 3.25e4 * df_currabsorber['A_absorber'] * projectiledf['Energy_curr'] / \
                    (projectiledf['Z_proj'] * df_currabsorber['Z_absorber'] *
                     (projectiledf['A_proj'] + df_currabsorber['A_absorber']) * calcdf['za'])

    calcdf['sigman'] = 1.7 * np.sqrt(calcdf['eps']) * np.log(calcdf['eps'] + 2.1718282) / \
                       (1.0 + 6.8 * calcdf['eps'] + 3.4 * np.sqrt(calcdf['eps']) ** 3)

    calcdf['dedxnu'] = calcdf['sigman'] * 5.105 * projectiledf['Z_proj'] * df_currabsorber['Z_absorber'] * \
                       projectiledf['A_proj'] / (calcdf['za'] * df_currabsorber['A_absorber'] *
                                                 (projectiledf['A_proj'] + df_currabsorber['A_absorber']))

    dedx_tot = calcdf['dedxnu'] + calcdf['dedxhi']

    #print("IN DEDX WE HAVE")
    #print(calcdf)
    #print("DEDX END")

    return dedx_tot

if __name__ == "__main__":

    # only work with tritons going through Alumimum now so we can see how this works...
    z_absorber = [[13]]
    a_absorber = [[27]]
    numa_absorber = [[1]]
    pressure = [0]
    length = [0]
    density = [2.7]
    thick = [137.16]
    isgas = [False]

    proj_z = [1]
    proj_a = [3]
    proj_ei = [13]

    df_f = desorb(proj_z, proj_a, proj_ei, z_absorber, a_absorber, numa_absorber, isgas, density, thick, pressure, length)

    print(df_f)
