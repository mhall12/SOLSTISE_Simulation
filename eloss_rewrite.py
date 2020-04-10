import numpy as np
import pandas as pd


def desorb(z_projectile, a_projectile, energy, z_absorber, a_absorber, numa_absorber, isgas, density, thickness):

    # desorb functions close to what the old desorb does. However, it is rewritten to be readable...
    # nsand is the number of elements in the absorber.
    # z_, a_, en_ are the Z, A and energy of the projectile (the thing that is losing energy!)
    # isgas is a boolean True/False to determine whether or not the absorber is a gas or solid.
    # density is g/cm^2 and thickness is mg/cm^2.

    # z_absorber and a_absorber and numa_absorber are going to be arrays. So for example, if we had CD2,
    # z_absorber = [[6, 1]], a_absorber = [[12, 2]], numa_absorber = [[1,2]]. These are 2D arrays. Z, A are columns and
    # rows are the layers.
    # Density and thickness are 1D arrays where each element corresponds to a layer.

    # Here we'll set the partial thicknesses assuming it's a solid absorber
    thick_frac = []
    den_frac = []
    z_ind = []
    a_ind = []
    numa_ind = []
    for i in range(len(z_absorber)):
        a_tot = 0
        for j in range(len(z_absorber[i])):
            a_tot = a_tot + a_absorber[i][j] * numa_absorber[i][j]
        for j in range(len(z_absorber[i])):
            a_frac = a_absorber[i][j] * numa_absorber[i][j] / a_tot
            thick_frac.append(thickness[i] * a_frac)
            den_frac.append(density[i] * a_frac)
            z_ind.append(z_absorber[i][j])
            a_ind.append(a_absorber[i][j])
            numa_ind.append(numa_absorber[i][j])

    # We can set up the DataFrame down here:
    df_abs = pd.DataFrame()

    df_abs['Z_absorber'] = z_ind
    df_abs['A_absorber'] = a_ind
    df_abs['Number_Atoms'] = numa_ind
    df_abs['Partial_Density'] = den_frac
    df_abs['Partial_Thickness'] = thick_frac


    projectiledf = pd.DataFrame()

    projectiledf['Z_proj'] = z_projectile
    projectiledf['A_proj'] = a_projectile
    projectiledf['Energy_i'] = energy

    # Right now only accept non gasses...
    if not isgas:
        absorberdf['Z_Abs'] = z_absorber
        absorberdf['A_Abs'] = a_absorber
        absorberdf['NumbA_Abs'] = numa_absorber
        absorberdf['Density'] = density
        absorberdf['Thickness'] = thickness
        absorberdf['Length'] = absorberdf['Thickness'] / absorberdf['Density'] * 1000.0
        # Thickness is in mg/cm^2 and Density is in g/cm^3, so need to multiply by 1000. The original code also takes
        # the length of the absorber and converts it to inches but I don't know why that'd be necessary...

    # calculate partial densities and thicknesses:

    absorberdf['Partial Density'] =

    # Eventually define a function here that will take the composite absorber and calc partial densities and thicknesses
    # also eventually define a function to set up the gas part...


def eloss(z_projectile, a_projectile, energy, absorberdf, projectiledf):
    # initial number of integrations that's hard coded into the original desorb code.
    num_integrations = 2
    num_elements = 1  # keep this just one for now, can expand later

    k = 0
    while k < num_integrations:
        k = k + 1
        j_end = num_elements
        j = 0
        while j < j_end:
            absorberdf['FX'] = absorberdf['Thickness'] / num_integrations
            projectiledf['Velocity'] = np.sqrt(2.13e-3 * energy / a_projectile)

    return 0

def dedx(absorberdf, projectiledf):
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

    xi = projectiledf['Velocity']**2 / absorberdf['Z_Abs']




if __name__ == "__main__":
    desorb(1, 3, 25, 13, 27, 1, False, 2.7, 137.16)
