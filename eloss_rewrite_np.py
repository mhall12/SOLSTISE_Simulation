import numpy as np
import pandas as pd
import random


class ab():
    thick_frac = np.empty(0)
    den_frac = np.empty(0)
    # length_frac = np.empty(0)
    isgas = np.empty(0)
    a = np.empty(0)
    z = np.empty(0)
    numa = np.empty(0)


def desorb(z_projectile, a_projectile, energy, z_absorber, a_absorber, numa_absorber, isgas,
           density, thickness, pressure, length):

    # desorb functions close to what the old desorb does. However, it is rewritten to use pandas and numpy.
    # the old functions SETABG and SETABS are no longer used, and are replaced by the for loops below.
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

    for i in range(len(z_absorber)):
        a_tot = 0
        if not isgas[i]:
            for j in range(len(z_absorber[i])):
                a_tot = a_tot + a_absorber[i][j] * numa_absorber[i][j]
            for j in range(len(z_absorber[i])):
                a_frac = (a_absorber[i][j] * numa_absorber[i][j]) / a_tot
                ab.thick_frac = np.append(ab.thick_frac, thickness[i] * a_frac)
                ab.den_frac = np.append(ab.den_frac, density[i] * a_frac)
                # ab.length_frac = np.append(ab.length_frac, thickness[i] / density[i] * 1000.0)
                ab.z = np.append(ab.z, z_absorber[i][j])
                ab.a = np.append(ab.a, a_absorber[i][j])
                ab.numa = np.append(ab.numa, numa_absorber[i][j])
                ab.isgas = np.append(ab.isgas, False)
        if isgas[i]:
            for j in range(len(z_absorber[i])):
                a_tot = a_tot + a_absorber[i][j] * numa_absorber[i][j]
                ab.thick_frac = np.append(ab.thick_frac, (pressure[i] / 760) * (length[i] / 22.4) * a_absorber[i][j] *
                                          numa_absorber[i][j])
                ab.den_frac = np.append(ab.den_frac, ab.thick_frac[j] / length[i])
                ab.z = np.append(ab.z, z_absorber[i][j])
                ab.a = np.append(ab.a, a_absorber[i][j])
                ab.numa = np.append(ab.numa, numa_absorber[i][j])
                ab.isgas = np.append(ab.isgas, True)

    index = np.arange(len(energy))

    df_fin = eloss(z_projectile, a_projectile, energy, index)
    de_mask = df_fin['DeltaE_tot'] > 0.05

    df_fin.loc[de_mask, 'E_strag_FWHM'] = 0.03520 * df_fin['DeltaE_tot']**(-0.69964) * df_fin['DeltaE_tot']
    df_fin.loc[~de_mask, 'E_strag_FWHM'] = 0.0207 * df_fin['DeltaE_tot'] ** (-0.583) * df_fin['DeltaE_tot']

    df_fin['dE_wStrag'] = np.random.normal(df_fin['DeltaE_tot'], df_fin['E_strag_FWHM'])

    allelostmask = df_fin['Energy_i'] == df_fin['DeltaE_tot']

    df_fin['E_fin'] = df_fin['Energy_i'] - df_fin['dE_wStrag']

    df_fin.loc[allelostmask, 'E_fin'] = 0

    negeloss_mask = df_fin['E_fin'] > df_fin['Energy_i']
    df_fin.loc[negeloss_mask, 'E_fin'] = df_fin['Energy_i']

    return df_fin


def eloss(z_projectile, a_projectile, energy, index):
    # initial number of integrations that's hard coded into the original desorb code.
    k = 0
    e_init = energy
    e_curr = energy
    a_proj = a_projectile
    z_proj = z_projectile
    indx = index

    e_init_dead = np.empty(0)
    e_curr_dead = np.empty(0)
    a_dead = np.empty(0)
    z_dead = np.empty(0)
    indx_dead = np.empty(0)

    dednext = np.zeros_like(energy)
    ded1st = np.zeros_like(dednext)
    eps = 0.0001

    ddd = np.zeros_like(energy)
    ddr = np.zeros_like(energy)
    dds = np.ones_like(energy)

    k = np.zeros_like(energy)
    num_int = np.zeros_like(energy) + 2
    num_integrations = np.amax(num_int)
    k_min = np.amin(k)

    while k_min < num_integrations:
        k = k + 1
        k_min = np.amin(k)
        j_end = len(ab.z)
        j = 0

        #absorberdf['FX'] = absorberdf['Partial_Thickness'] / num_integrations
        while j < j_end:
            # Loop through each row of the absorber dataframe by putting each row into df_currabsorber and using that

            fx = ab.thick_frac[j] / num_int
            # here we loop through each element, thickness etc to calculate how the energy changes after going through
            # each partial layer.
            vel = np.sqrt(2.13e-3 * e_curr / a_proj)

            # Projectile velocity will change every loop step depending on the energy.
            # need to change Energy_curr in the dataframe.
            #print("Enter")
            #print(projectiledf)
            delta_e = dedx(z_proj, a_proj, e_curr, vel, j)

            # sign is always -1, so just hard code it in
            e_curr = e_curr + delta_e * -1 * fx

            egt0mask = e_curr > 0

            #print(projectiledf['Egt0'])

            #print(projectiledf)

            # Now we want to remove any values where the energy is 0 so the energy loss doesn't get calculated again

            #print(projectiledf_dead)

            # Collect all the particles that have 0 energy here:
            z_dead = np.append(z_dead, z_proj[np.invert(egt0mask)])
            a_dead = np.append(a_dead, a_proj[np.invert(egt0mask)])
            e_init_dead = np.append(e_init_dead, e_init[np.invert(egt0mask)])
            e_curr_dead = np.append(e_curr_dead, e_curr[np.invert(egt0mask)])
            indx_dead = np.append(indx_dead, indx[np.invert(egt0mask)])

            k = k[egt0mask]
            fx = fx[egt0mask]
            num_int = num_int[egt0mask]
            dednext = dednext[egt0mask]
            ded1st = ded1st[egt0mask]
            ddd = ddd[egt0mask]
            ddr = ddr[egt0mask]
            dds = dds[egt0mask]

            z_proj = z_proj[egt0mask]
            a_proj = a_proj[egt0mask]
            e_init = e_init[egt0mask]
            e_curr = e_curr[egt0mask]
            vel = vel[egt0mask]
            delta_e = delta_e[egt0mask]
            indx = indx[egt0mask]

            #print(projectiledf_dead)

            klt2mask = k <= 2
            dednext = np.where(np.invert(klt2mask), dednext, dednext + delta_e * fx)

            j = j + 1

            #if k < 50:
             #   print(k)

        k1mask = k == 1

        ded1st = np.where(np.invert(k1mask), ded1st, dednext)
        dednext = np.where(np.invert(k1mask), dednext, np.zeros_like(dednext))

        k2mask = k == 2
        ddd = np.where(np.invert(k2mask), ddd, ded1st - dednext)

        dddmask = ddd < 0
        k2dddmask = k2mask & dddmask

        ddd = np.where(np.invert(k2dddmask), ddd, ddd * -1.0)

        dds = np.where(np.invert(k2mask), dds, ded1st + dednext)

        ddr = np.where(np.invert(k2mask), ddr, ddd / dds)

        ddr_mask = ddr > eps
        k2ddrmask = k2mask & ddr_mask

        if len(ddr[k2ddrmask]) > 0:
            num_int = np.where(np.invert(k2ddrmask), num_int, num_int * 2)

            j = -1
            k = np.where(np.invert(k2ddrmask), k, k * 0.0)
            dednext = np.where(np.invert(k2ddrmask), dednext, np.zeros_like(dednext))

            e_curr = np.where(k2ddrmask, e_init, e_curr)

        if e_init.size == 0:
            break

        num_integrations = np.amax(num_int)
        k_min = np.amin(k)
        keqmax_mask = k == num_int
        k = k[np.invert(keqmax_mask)]
        num_int = num_int[np.invert(keqmax_mask)]
        dednext = dednext[np.invert(keqmax_mask)]

        ded1st = ded1st[np.invert(keqmax_mask)]
        ddd = ddd[np.invert(keqmax_mask)]
        ddr = ddr[np.invert(keqmax_mask)]
        dds = dds[np.invert(keqmax_mask)]

        z_dead = np.append(z_dead, z_proj[keqmax_mask])
        a_dead = np.append(a_dead, a_proj[keqmax_mask])
        e_init_dead = np.append(e_init_dead, e_init[keqmax_mask])
        e_curr_dead = np.append(e_curr_dead, e_curr[keqmax_mask])
        indx_dead = np.append(indx_dead, indx[keqmax_mask])

        z_proj = z_proj[np.invert(keqmax_mask)]
        a_proj = a_proj[np.invert(keqmax_mask)]
        e_init = e_init[np.invert(keqmax_mask)]
        e_curr = e_curr[np.invert(keqmax_mask)]
        indx = indx[np.invert(keqmax_mask)]



    # Remake the original dataframe to include all of the rows that were removed previously.
    z_proj = np.append(z_proj, z_dead)
    a_proj = np.append(a_proj, a_dead)
    e_init = np.append(e_init, e_init_dead)
    e_curr = np.append(e_curr, e_curr_dead)
    indx = np.append(indx, indx_dead)

    sortindx = np.argsort(indx)

    z_proj = z_proj[sortindx]
    a_proj = a_proj[sortindx]
    e_init = e_init[sortindx]
    e_curr = e_curr[sortindx]

    elt0mask = e_curr < 0
    #print(e_curr, )
    e_curr = np.where(elt0mask, np.zeros_like(e_curr), e_curr)

    projectiledf = pd.DataFrame()

    projectiledf['Energy_i'] = e_init
    projectiledf['Z_proj'] = z_proj
    projectiledf['A_proj'] = a_proj

    projectiledf['DeltaE_tot'] = e_init - e_curr

    return projectiledf

def dedx(z_proj_in, a_proj_in, e_curr_in, vel_in, j):
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

    z_abs = ab.z[j]
    a_abs = ab.a[j]
    isg = ab.isgas[j]
    part_den = ab.den_frac[j]

    z_proj = z_proj_in
    a_proj = a_proj_in
    en = e_curr_in
    vel = vel_in

    rho = 0

    # Set the density here. Not sure how it gets the density if it's a gas, but we'll see later maybe
    if not isg:
        rho = part_den
    elif isg:
        rho = 1

    xi = vel**2 / z_abs


    # Absorber function
    # G(XI) = Y(EXP) - Y(Theory) is deduced from experimental energy loss measurements
    fy = 0
    # fy is function y
    if not isg:
        fy = 54721.0 * (1.0 + 5.15e-2 * np.sqrt(a_abs / rho) -
                        np.exp(-0.23 * z_abs))
    elif isg:
        fy = 54721.0 * (1.35 - np.exp(z_abs * (-0.13 + 0.0014 * z_abs)))

    # G(XI) is the derivation of a gaussian with variable height H(Z2)

    g1 = 0
    if z_abs <= 26.0:
        g1 = 19.84 * np.exp(-0.17 * (z_abs - 4.25) * (z_abs - 4.25))
    elif z_abs > 26.0:
        g1 = 0.000001

    g2 = 0
    if z_abs < 38.0:
        g2 = 17.12 * np.exp(-0.12 * (z_abs - 11.63)**2)
    elif z_abs > 38.0:
        g2 = 0.0000001

    g3 = 7.95 * np.exp(-0.015 * (z_abs - 30.2) * (z_abs - 30.2))
    g4 = 5.84 * np.exp(-0.022 * (z_abs - 48.63) * (z_abs - 48.63))
    g5 = 7.27 * np.exp(-0.005 * (z_abs - 73.06) * (z_abs - 73.06))
    hz2 = (9.0 - (g1 + g2 + g3 + g4 + g5)) * 1.32e-5

    z2zwd = np.cbrt(z_abs) * np.cbrt(z_abs)

    # Multiplication factors of G(XI)
    fg = 1.2e-4 * z_abs * z_abs + \
         (2.49e-2 * a_abs / rho)

    if isg:
        fg = 1.3 / (1.0 + np.exp(3.0 - (z_abs / 5.0)))

    alefg = np.log(2.7e-5 / fg)

    # Calculation of G(XI)

    gxi = xi * 0.0

    xi_mask = (1.0e-9 <= xi) & (xi <= 5.0e-4)

    #print(calcdf,calcdf_masked)

    if xi_mask[xi_mask].size > 0:
        sqxi = np.sqrt(xi)

        c = 2.0 / z_abs * (sqxi / (1.0 + 1.0e4 * sqxi))

        if isg:
            c = c / 2.0

        fg0 = 1.0 / (1.0 + (xi * 10000.0) * (xi * 10000.0) *
                     (xi * 10000.0))
        al = np.log(xi) - alefg

        gxi = np.where(np.invert(xi_mask), gxi, (c - hz2 * al * np.exp(-0.32 * al * al)) * fg0)


    # Calculation of Y(XI)
    y = 3.3e-4 * np.log(1.0 + (xi * fy)) + gxi

    # Energy loss of heavy ions
    # Effective charge

    vv0 = vel * 137.0

    fv = xi * 0.0 + 1.0
    #print(calcdf['fv'])

    velmask = vel <= 0.62

    fv = np.where(np.invert(velmask), fv, 1.0 - np.exp(-vv0))

    az1 = np.log(1.035 - 0.4 * np.exp(-0.16 * z_proj))

    qq = vel / z_proj**0.509

    ghi = z_proj

    vz1 = (-116.79 - 3350.4 * qq) * qq

    vz1mask = vz1 > -85.2

    ghi = np.where(np.invert(vz1mask), ghi, z_proj * (1.0 - np.exp(vz1)))

    zprojmask = z_proj > 2.0

    ghi = np.where(np.invert(zprojmask), ghi, z_proj * (1.0 - np.exp(fv * az1 - 0.879 * (vv0 / z_proj**0.65))))


    # Effective charge of protons and aphas
    # Electronic energy loss DEDXHI

    dedxhi = ghi**2 * z_abs * y / \
                       (a_abs * vel**2)

    # nuclear energy loss DEDXNU
    za = np.sqrt(np.cbrt(z_proj)**2 + z2zwd)

    eps = 3.25e4 * a_abs * en / \
                    (z_proj * z_abs *
                     (a_proj + a_abs) * za)

    sigman = 1.7 * np.sqrt(eps) * np.log(eps + 2.1718282) / \
                       (1.0 + 6.8 * eps + 3.4 * np.sqrt(eps) ** 3)

    dedxnu = sigman * 5.105 * z_proj * z_abs * \
                       a_proj / (za * a_abs *
                                                 (a_proj + a_abs))

    dedx_tot = dedxnu + dedxhi

    return dedx_tot

if __name__ == "__main__":

    z_absorber = [[6,8],[6,1],[13]]
    a_absorber = [[12,16],[12,1],[26]]
    numa_absorber = [[1,2],[1,2],[1]]
    pressure = [100,0,0]
    length = [2,0,0]
    density = [0,.9, 2.7]
    thick = [0,1, 130]
    isgas = [True, False, False]

    proj_z = np.array([1,1,1,1,1,1,1,1])
    proj_a = np.array([3,3,3,3,3,3,3,3])
    proj_ei = np.array([4,10,20,50,15,7,33,30])

    df_f = desorb(proj_z, proj_a, proj_ei, z_absorber, a_absorber, numa_absorber, isgas, density, thick, pressure, length)

    print(df_f)
