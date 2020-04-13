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


def eloss(absorberdf, projectiledf):
    # initial number of integrations that's hard coded into the original desorb code.
    num_integrations = 2
    num_elements = 1  # keep this just one for now, can expand later

    k = 0
    projectiledf['Energy_curr'] = projectiledf['Energy_i']

    while k < num_integrations:
        k = k + 1
        j_end = len(absorberdf.index)
        j = 0
        absorberdf['FX'] = absorberdf['Thickness'] / num_integrations
        while j < j_end:
            df_currabsorber = absorberdf.iloc[j]
            # here we loop through each element, thickness etc to calculate how the energy changes after going through
            # each partial layer.
            projectiledf['Velocity'] = np.sqrt(2.13e-3 * projectiledf['Energy_curr'] / projectiledf['A_proj'])
            # Projectile velocity will change every loop step depending on the energy.
            # need to change Energy_curr in the dataframe.
            projectiledf = dedx(df_currabsorber, projectiledf)


    return 0

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

    rho = 0

    # Set the density here. Not sure how it gets the density if it's a gas, but we'll see later maybe
    if not df_currabsorber['Gas?']:
        rho = df_currabsorber['Partial_Density']
    elif df_currabsorber['Gas?']:
        rho = 1

    xi = projectiledf['Velocity']**2 / df_currabsorber['Z_absorber']

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

    gxi = 0

    if 1.0e-9 <= xi <= 5.0e-4:
        sqxi = np.sqrt(xi)
        c = 2.0 / df_currabsorber['Z_absorber'] * (sqxi / (1.0 + 1.0e4 * sqxi))

        if df_currabsorber['Gas?']:
            c = c / 2.0

        fg0 = 1.0 / (1.0 + (xi * 10000.0) * (xi * 10000.0) * (xi * 10000.0))
        al = np.log(xi) - alefg
        gxi = (c - hz2 * al * np.exp(-0.32 * al * al)) * fg0


    # Calculation of Y(XI)
    y = 3.3e-4 * np.log(1.0 + (xi * fy)) + gxi

    # Energy loss of heavy ions
    # Effective charge

    vv0 = projectiledf['Velocity'] * 137.0
    fv = 1.0

    if projectiledf['Velocity'] >= 0.62:
        fv = 1.0 - np.exp(-vv0)

    az1 = np.log(1.035 - 0.4 * np.exp(-0.16 * projectiledf['Z_proj']))

    qq = projectiledf['Velocity'] / projectiledf['Z_proj']**0.509
    
    ghi = projectiledf['Z_proj']

    vz1 = (-116.79 - 3350.4 * qq) * qq

    return delta_E

if __name__ == "__main__":
    desorb(1, 3, 25, 13, 27, 1, False, 2.7, 137.16)
