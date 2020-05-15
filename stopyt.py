import numpy as np
import pandas as pd
import random
import re
import warnings


class ab():
    # The ab class defines the numpy arrays used for the absorber in the energy loss calculations
    # thick_frac is the fractional thickness of each isotope based on the layer composition and total layer thickness
    # in units of mg/cm^2
    thick_frac = np.empty(0)
    # den_frac is the fractional density (g/cm^3)
    den_frac = np.empty(0)
    # isgas is a boolean array that specifies if the layer is a solid or gas
    isgas = np.empty(0)
    # a is the atomic mass A of an isotope in the layer
    a = np.empty(0)
    # a is the number of Z of an isotope in the layer
    z = np.empty(0)
    # numa is the number of atoms contained in the compound of each isotope
    numa = np.empty(0)


def desorb(z_projectile, a_projectile, energy, z_absorber, a_absorber, numa_absorber, isgas,
           density, thickness, pressure, length):
    # desorb functions close to what the old desorb does. However, it is rewritten to use pandas and numpy.
    # the old functions SETABG and SETABS are no longer used, and are replaced by the for loops below.

    # Generally, it takes the input parameters for the absorber and makes numpy arrays for them, sends those
    # params and the projectile data to the energy loss code, then outputs the results.

    # Resets the absorber between desorb calls, so if the user running the program changes the absorber or runs it again
    # the program doesn't save the old absorber data and just append the new data onto the old.
    ab.thick_frac = np.empty(0)
    ab.den_frac = np.empty(0)
    # length_frac = np.empty(0)
    ab.isgas = np.empty(0)
    ab.a = np.empty(0)
    ab.z = np.empty(0)
    ab.numa = np.empty(0)

    # z_absorber and a_absorber and numa_absorber are regular python lists. So for example, if we had CD2,
    # z_absorber = [[6, 1]], a_absorber = [[12, 2]], numa_absorber = [[1,2]]. These are 2D arrays. Z, A are columns and
    # rows are the layers.
    # Density and thickness are 1D arrays where each element corresponds to a layer.

    # Input length, pressure only for gas. Input density, thickness only for non gasses

    # Here we'll set the partial thicknesses assuming it's a solid absorber

    a_tot = 0
    if not isgas:
        # if the layer isn't a gas, sum the mass numbers of the isotopes in the layer first
        for j in range(len(z_absorber)):
            a_tot = a_tot + a_absorber[j] * numa_absorber[j]
        for j in range(len(z_absorber)):
            # the sum is then used to calculate the fracitional thicknesses and densities for each element
            a_frac = (a_absorber[j] * numa_absorber[j]) / a_tot
            ab.thick_frac = np.append(ab.thick_frac, thickness * a_frac)
            ab.den_frac = np.append(ab.den_frac, density * a_frac)
            # the z, a, numa, and isgas bool from each layer element are appended to the numpy arrays
            ab.z = np.append(ab.z, z_absorber[j])
            ab.a = np.append(ab.a, a_absorber[j])
            ab.numa = np.append(ab.numa, numa_absorber[j])
            ab.isgas = np.append(ab.isgas, False)

    if isgas:
        # we have to do something different if it's a gas, though it's similar.
        for j in range(len(z_absorber)):
            a_tot = a_tot + a_absorber[j] * numa_absorber[j]
            # the fractional thickness for the gas layer is calculated here. The pressure is in torr, so the 760
            # converts that to atms, and the 22.4 is the molar gas volume at STP in L.
            ab.thick_frac = np.append(ab.thick_frac, (pressure / 760) * (length / 22.4) * a_absorber[j] *
                                          numa_absorber[j])
            ab.den_frac = np.append(ab.den_frac, ab.thick_frac[j] / length)
            # the rest of this section proceeds the same as the first. Actually, some of this is redundant and could
            # be reorganized in the future. i.e. move to one for loop
            ab.z = np.append(ab.z, z_absorber[j])
            ab.a = np.append(ab.a, a_absorber[j])
            ab.numa = np.append(ab.numa, numa_absorber[j])
            ab.isgas = np.append(ab.isgas, True)

    # index is the original position of each of the stopping particles. Keeping this consistent is important for the
    # simulation especially.
    index = np.arange(len(energy))

    # eloss runs the energy loss code and outputs a dataframe with the final energies and projectile params
    df_fin = eloss(z_projectile, a_projectile, energy, index)

    # Here, the energy straggling is estimated. The equation used changes based on the energy loss. The mask picks out
    # the particles that lost a low amount of energy.
    de_mask = df_fin['DeltaE_tot'] > 0.05

    # These equations were found based on LISE calculations. The straggling is close to those calculated values but
    # are just estimates.
    df_fin.loc[de_mask, 'E_strag_FWHM'] = 0.03520 * df_fin['DeltaE_tot']**(-0.69964) * df_fin['DeltaE_tot']
    df_fin.loc[~de_mask, 'E_strag_FWHM'] = 0.0207 * df_fin['DeltaE_tot'] ** (-0.583) * df_fin['DeltaE_tot']

    nanmask = np.isnan(df_fin['DeltaE_tot'])
    nanmask2 = np.isnan(df_fin['E_strag_FWHM'])
    df_fin.loc[nanmask, 'DeltaE_tot'] = 0.0
    df_fin.loc[nanmask2, 'E_strag_FWHM'] = 0.0

    # The final dE is calculated for each particle as a normal distribution centered on the calculated energy whose
    # FWHM is the straggling energy calculated above

    df_fin['dE_wStrag'] = np.random.normal(df_fin['DeltaE_tot'], df_fin['E_strag_FWHM'])

    # Need to make the final energy 0 if all the energy is lost. With the normal distribution, the energy loss could
    # actually be negative, so we need to fix that here as well...
    allelostmask = df_fin['Energy_i'] == df_fin['DeltaE_tot']

    df_fin['E_fin'] = df_fin['Energy_i'] - df_fin['dE_wStrag']

    df_fin.loc[allelostmask, 'E_fin'] = 0

    negeloss_mask = df_fin['E_fin'] > df_fin['Energy_i']
    df_fin.loc[negeloss_mask, 'E_fin'] = df_fin['Energy_i']

    # Return the final dataframe.
    return df_fin


def eloss(z_projectile, a_projectile, energy, index):
    # initial number of integrations that's hard coded into the original desorb code.
    k = 0
    # intialize the values for the projectile. We could just use the four inputs though. The absorber data isn't
    # passed through because it's a class.
    e_init = energy
    e_curr = energy
    a_proj = a_projectile
    z_proj = z_projectile
    indx = index

    # the _dead numpy arrays collect the particles that have hit the max number of integrations or have lost all their
    # energy.
    e_init_dead = np.empty(0)
    e_curr_dead = np.empty(0)
    a_dead = np.empty(0)
    z_dead = np.empty(0)
    indx_dead = np.empty(0)

    # initialize the parameters used for calculations with zero arrays of the correct length.
    dednext = np.zeros_like(energy)
    ded1st = np.zeros_like(dednext)
    eps = 0.0001

    ddd = np.zeros_like(energy)
    ddr = np.zeros_like(energy)
    dds = np.ones_like(energy)

    # intial number of integrations is 2. We need to keep track of the integration number and max integration number for
    # each particle separately, so they each get their own arrays.
    k = np.zeros_like(energy)
    num_int = np.zeros_like(energy) + 2
    # we need to run the loop through the maximum number of integrations so we need to keep track of the minimum k.
    num_integrations = np.amax(num_int)
    k_min = np.amin(k)

    # the while loop runs as long as the minimum k has not hit the max number of integrations.
    while k_min < num_integrations:
        # increase k for each particle
        k = k + 1
        # get the new kmin
        k_min = np.amin(k)
        j_end = len(ab.z)
        j = 0

        while j < j_end:

            # Loop through each row of the absorber array (jth element).

            # fx depends on the absorber element we're using.
            fx = ab.thick_frac[j] / num_int
            # here we loop through each element, thickness etc to calculate how the energy changes after going through
            # each partial layer.

            # Calculate the velocity here:
            vel = np.sqrt(2.13e-3 * e_curr / a_proj)

            # dedx calculates the dE energy loss and returns it.
            delta_e = dedx(z_proj, a_proj, e_curr, vel, j)

            # sign is always -1, so just hard code it in from the original desorb code.
            e_curr = e_curr + delta_e * -1 * fx

            # make a mask for the particles whose current energy is greater than 0. We'll use this to move particles
            # to the _dead arrays.
            egt0mask = e_curr > 0

            # Now we want to remove any values where the energy is 0 so the energy loss doesn't get calculated again
            # Collect all the particles that have 0 energy here:
            z_dead = np.append(z_dead, z_proj[np.invert(egt0mask)])
            a_dead = np.append(a_dead, a_proj[np.invert(egt0mask)])
            e_init_dead = np.append(e_init_dead, e_init[np.invert(egt0mask)])
            e_curr_dead = np.append(e_curr_dead, e_curr[np.invert(egt0mask)])
            indx_dead = np.append(indx_dead, indx[np.invert(egt0mask)])

            # We also need to remove those elements from the calculation parameters because we don't need them anymore.
            k = k[egt0mask]
            fx = fx[egt0mask]
            num_int = num_int[egt0mask]
            dednext = dednext[egt0mask]
            ded1st = ded1st[egt0mask]
            ddd = ddd[egt0mask]
            ddr = ddr[egt0mask]
            dds = dds[egt0mask]

            # also remove those E=0 particles from the current projectile arrays so we don't try to get an energy loss
            # for them.
            z_proj = z_proj[egt0mask]
            a_proj = a_proj[egt0mask]
            e_init = e_init[egt0mask]
            e_curr = e_curr[egt0mask]
            vel = vel[egt0mask]
            delta_e = delta_e[egt0mask]
            indx = indx[egt0mask]

            # if k <=2, we need to calculate dednext
            klt2mask = k <= 2
            dednext = np.where(np.invert(klt2mask), dednext, dednext + delta_e * fx)

            # Now add one to j so we can loop through the next absorber element.

            j = j + 1



        # need a mask for k = 1, then change ded1st and dednext if it is 1.
        k1mask = k == 1

        # np.where(bool, value if true, value if false): if we wanted to, we could not invert the mask and swap the
        # last two values...
        ded1st = np.where(np.invert(k1mask), ded1st, dednext)
        dednext = np.where(np.invert(k1mask), dednext, np.zeros_like(dednext))

        # need an additional mask for k = 2:
        k2mask = k == 2
        ddd = np.where(np.invert(k2mask), ddd, ded1st - dednext)

        dddmask = ddd < 0
        k2dddmask = k2mask & dddmask

        ddd = np.where(np.invert(k2dddmask), ddd, ddd * -1.0)

        dds = np.where(np.invert(k2mask), dds, ded1st + dednext)

        ddr = np.where(np.invert(k2mask), ddr, ddd / dds)

        ddr_mask = ddr > eps
        k2ddrmask = k2mask & ddr_mask

        # Need to increase the number of integrations when ddr > eps:
        if len(ddr[k2ddrmask]) > 0:
            num_int = np.where(np.invert(k2ddrmask), num_int, num_int * 2)

            j = -1
            k = np.where(np.invert(k2ddrmask), k, k * 0.0)
            dednext = np.where(np.invert(k2ddrmask), dednext, np.zeros_like(dednext))

            e_curr = np.where(k2ddrmask, e_init, e_curr)

        # Need to end the while loop if all the particles are in the _dead arrays
        if e_init.size == 0:
            break

        # Similar to the E=0 condition above, we need to remove particles that have hit their max number of
        # integrations and add them to the _dead arrays.
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

    # Here we need to add back in the particles that were removed to the original arrays. This step maybe isn't
    # necessary, since everything should be contained in the _dead array anyway...
    z_proj = np.append(z_proj, z_dead)
    a_proj = np.append(a_proj, a_dead)
    e_init = np.append(e_init, e_init_dead)
    e_curr = np.append(e_curr, e_curr_dead)
    indx = np.append(indx, indx_dead)

    # Sort the particles based on their index to ensure they are in the same location they started and don't get
    # jumbled up:
    sortindx = np.argsort(indx)

    z_proj = z_proj[sortindx]
    a_proj = a_proj[sortindx]
    e_init = e_init[sortindx]
    e_curr = e_curr[sortindx]

    # If the particle loses all of its energy, e_curr will be negative, so we make it 0 here.
    enanmask = np.isnan(e_curr)
    e_curr = np.where(enanmask, np.zeros_like(e_curr), e_curr)

    elt0mask = e_curr < 0

    e_curr = np.where(elt0mask, np.zeros_like(e_curr), e_curr)

    # Add all the particles into a dataframe for easy transport out of the function:
    projectiledf = pd.DataFrame()

    projectiledf['Energy_i'] = e_init
    projectiledf['Z_proj'] = z_proj
    projectiledf['A_proj'] = a_proj

    projectiledf['DeltaE_tot'] = e_init - e_curr

    return projectiledf


def dedx(z_proj, a_proj, en, vel, j):
    # Function calculates the differential energy loss dE/dX in solid targets using a semiempirical formula deduced
    # from experimental work

    # Set the A, Z, etc for the absorber from the array in the class. Setting them into their own variable isn't
    # really necessary but makes the code a little cleaner.
    z_abs = ab.z[j]
    a_abs = ab.a[j]
    isg = ab.isgas[j]
    part_den = ab.den_frac[j]

    # initialize the density to 0
    rho = 0

    # Set the density here. If the stopping layer is a gas, the density is just set to 1. If it's a gas, rho isn't
    # actually used below.
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
    fg = 1.2e-4 * z_abs * z_abs + (2.49e-2 * a_abs / rho)

    if isg:
        fg = 1.3 / (1.0 + np.exp(3.0 - (z_abs / 5.0)))

    alefg = np.log(2.7e-5 / fg)

    # Calculation of G(XI)

    gxi = xi * 0.0

    xi_mask = (1.0e-9 <= xi) & (xi <= 5.0e-4)

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

    eps = 3.25e4 * a_abs * en / (z_proj * z_abs * (a_proj + a_abs) * za)

    sigman = 1.7 * np.sqrt(eps) * np.log(eps + 2.1718282) / (1.0 + 6.8 * eps + 3.4 * np.sqrt(eps) ** 3)

    dedxnu = sigman * 5.105 * z_proj * z_abs * a_proj / (za * a_abs * (a_proj + a_abs))

    dedx_tot = dedxnu + dedxhi

    # Returns the dedx for that integration step
    return dedx_tot


if __name__ == "__main__":

    print("")
    print("########  ########  ########   ########  ##      ##  ########  ")
    print("##           ##     ##    ##   ##    ##   ##    ##      ##     ")
    print("##           ##     ##    ##   ##    ##    ##  ##       ##     ")
    print("########     ##     ##    ##   ########     ####        ##    ")
    print("      ##     ##     ##    ##   ##            ##         ##    ")
    print("      ##     ##     ##    ##   ##            ##         ##     ")
    print("########     ##     ########   ##           ####        ##    ")
    print("")


    # Here is the front end for stopyt, the python version of stopit. It functions in much the same way and can be
    # called from the terminal using >>python3 stopyt.py

    # isotopetable.txt contains elemental symbols, Z, and A for the most abundant isotope. This table allows the user
    # to enter compounds like CH2, for example, without specifying A.
    inFile = "isotopetable.txt"
    isotable = np.genfromtxt(inFile, delimiter='\t', dtype='unicode')

    # Split the data from the text file into individual numpy arrays:
    symb = isotable[:, 0]
    zarr = isotable[:, 1]
    aarr = isotable[:, 2]

    # Initialize option to some random high number that won't be used:
    option = 100

    # Initialize the layer data and projectile data:
    individuallayers = []
    individual_proj = []

    # While loop to keep the program running:
    while option != 0:

        print("\n***********************************************************\n"
              "1) Define the stopping particles. \n"
              "2) Define the absorber layers. \n"
              "3) Run the energy loss code. \n"
              "4) Find the layer thickness for a final projectile energy. \n"
              "0) Exit."
              "\n***********************************************************")

        # Need the try except to make sure the program doesn't crash if the user just presses enter.
        option = 100
        try:
            option = int(input("\nInput: "))
        except ValueError:
            print("\n*****Enter an integer number from the list!*****\n")

        print("\n")

        if option == 1:
            # When defining the projectile(s) the parameters for those particles must be reset first
            proj_zin = []
            proj_ain = []
            proj_ei = []

            proj_z = []
            proj_a = []
            proj_symblist = []

            # Here we need the user to specify the A of the isotope, because it makes the coding a little easier
            proj = input("\nEnter the isotopes of the particles to undergo the energy loss calculations, "
                         "\nseparated by spaces (i.e. 3He 18F 16O): ")
            projnum = proj.count(" ") + 1
            # split the isotopes by the spaces and put them into an array:
            individual_proj = proj.split()

            # This for loop does the handling splitting the A and symbol from the isotopes and getting the Z from the
            # text file.
            for i in range(len(individual_proj)):
                out = 0
                while out == 0:
                    try:
                        projsplit = re.split('(\d+)', individual_proj[i])

                    # The first element of the split list is an empty string "", the second is the A, and the 3rd is the
                    # symbol, which we then can use to mask the numpy array made by the text file and get the Z

                        proj_ain.append(int(projsplit[1]))
                        symbmask = symb == projsplit[2]
                        proj_zin.append(int(zarr[symbmask][0]))
                        out = 1
                    except IndexError:
                        print("The " + individual_proj[i] + " isotope was not found!")
                        oldin = individual_proj[i]
                        individual_proj[i] = input("Re-enter the isotope: ")

                # We let the user define the energies here. One of the neat things about stopyt is that you can enter
                # any number of particles and energies. The program also allows the user to define the initial, final
                # Energy and step size.
                out  = 0
                while out == 0:
                    try:
                        eninput = input("\nEnter the energies of the " + individual_proj[i] +
                                        " in MeV, OR enter an energy range and step size, \nseparated "
                                        "by spaces (i.e. 1 4 0.5): ")

                        # Split the energies into a list like we have split the isotopes before:
                        einsplit = eninput.split()

                        for c in range(len(einsplit)):
                            checknum = float(einsplit[c])
                        out = checknum
                    except ValueError:
                        print("\nERROR: The number you entered is not valid!")

                # If the user input 3 energies, we need to decide if they're 3 unique energies or a range....
                if len(einsplit) == 3:
                    # is it divisible? if not, enter each energy as an individual number. If the difference between the
                    # first two energies is divisible by the 3rd, it's the step size.
                    if (float(einsplit[1]) - float(einsplit[0])) % float(einsplit[2]) == 0 and \
                            (float(einsplit[1]) > float(einsplit[0])):
                        ensteps = int((float(einsplit[1]) - float(einsplit[0])) / float(einsplit[2])) + 1
                    else:
                        ensteps = len(einsplit)
                else:
                    ensteps = len(einsplit)

                # Now append each energy to the array individually, and if the user entered a range and steps, calculate
                # each energy the user wanted to input here:
                for j in range(ensteps):

                    if len(einsplit) == 3:
                        if (float(einsplit[1]) - float(einsplit[0])) % float(einsplit[2]) == 0 and \
                                (float(einsplit[1]) > float(einsplit[0])):
                            proj_ei.append(float(einsplit[0]) + j*float(einsplit[2]))
                        else:
                            proj_ei.append(float(einsplit[j]))
                    else:
                        proj_ei.append(float(einsplit[j]))

                    # We have to append the Z, A and element symbol to the respective lists for each energy
                    # that was specified:
                    proj_z.append(proj_zin[i])
                    proj_a.append(proj_ain[i])
                    proj_symblist.append(individual_proj[i])

            # We also need to convert these to numpy arrays. In the simulation code we'll have to make numpy arrays for
            # each of these for when we feed them into the energyloss code.
            proj_z = np.array(proj_z)
            proj_a = np.array(proj_a)
            proj_ei = np.array(proj_ei)
            # The symbol list is used for the output and contains the mass number and symbol.
            proj_symblistnp = np.array(proj_symblist)

        if option == 2:
            # Define the absorbers here: we need to first reset all of the parameters in case this gets called twice
            numa_absorber = []
            ele_absorber = []
            a_absorber = []
            z_absorber = []
            prs = []
            den = []
            thk = []
            leng = []
            # Like the projectile definition, we want the user to specify the materials with a space. In this case,
            # since we're defining materials, in general you'll just want the most abundant isotope, so A is taken
            # from the text file for compounds like CO2.
            layerdata = input("Enter the material in each layer, separated by spaces (i.e. CO2 3He Si): ")
            numlayers = layerdata.count(" ") + 1

            individuallayers = layerdata.split()
            # print(individuallayers)

            for i in range(len(individuallayers)):
                if individuallayers[i] == 'mylar':
                    individuallayers[i] = 'C10H8O4'
                if individuallayers[i] == 'P10':
                    individuallayers[i] = 'C2H8Ar90'

            for i in individuallayers:
                # also need to reset these arrays if it gets run twice
                ele = []
                numele = []
                aind = []
                zind = []
                isgas = []

                # It gets a little complicated because I want the user to be able to enter a specific isotope if they
                # need to. So, we need to check whether or not the first character is a number:
                if i[0].isdigit():
                    out = 0
                    while out == 0:
                        try:
                            # if it's a number, we want to split the number and the symbol like we do for the projectile
                            isotope = re.split('(\d+)', i)
                            a_absorber.append([int(isotope[1])])
                            ele_absorber.append([isotope[2]])
                            # Normally you'd want to specify the number of atoms in the molecule, but that's only for
                            # calculating the partial densities and thicknesses, and if there's only one atom like 3He
                            # we don't actually need to do that...
                            numa_absorber.append([1])
                            symbmask = symb == isotope[2]
                            # Get Z using the symbol mask.
                            z_absorber.append([int(zarr[symbmask][0])])
                            out = 1
                        except IndexError:
                            print("\nThe isotope was not found!")
                            i = input("Re-enter the " + i + " isotope: ")
                else:
                    # If the first char isn't a number it must be a compound or standard isotope like Si:
                    buff = re.split('(\d+)', i)
                    for j in range(len(buff)):
                        if buff[j] != '' and not buff[j].isdigit():
                            # if the element of the split isn't a digit, we need to grab the element symbols. This is
                            # complicated by the fact that we could have some complicated compound like CaFC, where
                            # there are no numbers, and so we need a way to split that up by upper case letter followed
                            # by a lower case letter, so we'll get ['Ca', 'F', 'C']:
                            ind = re.findall('([A-Z][a-z]*)', buff[j])
                            for k in range(len(ind)):
                                out = 0
                                while out == 0:
                                    try:
                                        # Then we get the Z and A like normal
                                        symbmask = symb == ind[k]
                                        ele.append(ind[k])
                                        zind.append(int(zarr[symbmask][0]))
                                        aind.append(int(aarr[symbmask][0]))
                                        out = 1
                                    except IndexError:
                                        print("\nERROR: The " + ind[k] + " isotope in the " + i +
                                              " layer was not found!")
                                        ind[k] = input("Re-enter this isotope: ")

                                # We need to also append the number of elements onto the numele array, if they're
                                # present, otherwise we'll just append one
                                if k == len(ind) - 1 and len(buff) > 1 and j < len(buff) - 1:
                                    numele.append(int(buff[j + 1]))
                                else:
                                    numele.append(1)
                    a_absorber.append(aind)
                    z_absorber.append(zind)
                    numa_absorber.append(numele)
                    ele_absorber.append(ele)

            gslist = []

            for i in range(numlayers):
                # Define standardized media gas/solid, else indicate for each layer
                if individuallayers[i] == "CH2" or individuallayers[i] == "Si" or individuallayers[i] == "CD2" or \
                        individuallayers[i] == "C" or individuallayers[i] == "Al" or individuallayers[i] == "C10H8O4":
                    isgas.append(False)
                elif individuallayers[i] == "CO2" or individuallayers[i] == "C4H10" or individuallayers[i] == "CF4" or \
                        individuallayers[i] == 'C2H8Ar90':
                    isgas.append(True)
                else:
                    gs = input("Indicate whether the " + individuallayers[i] + " layer is a gas (g) or solid (s): ")
                    if gs == 'g':
                        gs = True
                    else:
                        gs = False

                    isgas.append(gs)

                # Ask the user for the pressure and length of the gas absorber:
                if isgas[i]:
                    out = 0
                    while out == 0:
                        try:
                            prs.append(float(input("For the " + individuallayers[i] +
                                                   " layer, enter the pressure in Torr: ")))

                            leng.append(float(input("For the " + individuallayers[i] +
                                                    " layer, enter the length in cm: ")))
                            out = 1
                        except ValueError:
                            print("ERROR: The number you entered was not valid!")

                    # Density and thickness are 0 for the gas Thickness is calculated from the pressure and length
                    den.append(0)
                    thk.append(0)
                # Append the densities for the standardized media here:
                if not isgas[i]:
                    if individuallayers[i] == "CH2" or individuallayers[i] == "CD2":
                        den.append(0.94)
                    elif individuallayers[i] == "C":
                        den.append(2.267)
                    elif individuallayers[i] == "Al":
                        den.append(2.7)
                    elif individuallayers[i] == "C10H8O4":
                        den.append(1.38)
                    # Silicon is a little different because we want to define the thickness in micons because that's
                    # how we measure detector thicknesses...
                    elif individuallayers[i] == "Si":
                        den.append(2.33)
                        thsi = float(input("For the " +
                                           individuallayers[i] + " (" +
                                           str(i+1) + ") layer, enter the thickness in microns: "))

                        # thickness in micron x density x 0.1 to convert um*g/cm^3 to mg/cm^2
                        thk.append(thsi * 2.33 * 0.1)

                    # if it isn't standardized, we want to specifically ask for the density
                    else:
                        den.append(float(input("For the " + individuallayers[i] +
                                               " layer, enter the density in g/cm^3: ")))

                    # We also need to ask for the thickness, but this was already done above for the silicon
                    if individuallayers[i] != "Si":
                        thk.append(float(input("For the " + individuallayers[i] +
                                               " layer, enter the thickness in mg/cm^2: ")))

                    # Append 0 here for the pressure and length because we don't need it.
                    prs.append(0)
                    leng.append(0)

        if option == 3:
            # If the user wants to run the energyloss code, they first need to run option 1 and 2 because it won't work
            # otherwise.
            if len(individual_proj) == 0:
                print("Please define the stopping particles using Option 1 first.")
            if len(individuallayers) == 0:
                print("Please define the stopping layers using Option 2 first.")
            if len(individual_proj) > 0 and len(individuallayers) > 0:
                proj_energycurr = proj_ei
                de_tot_layer = []
                for i in range(len(z_absorber)):
                    warnings.filterwarnings("ignore")
                    df_f = desorb(proj_z, proj_a, proj_energycurr, z_absorber[i], a_absorber[i], numa_absorber[i], isgas[i],
                              den[i], thk[i], prs[i], leng[i])
                    proj_energycurr = df_f['Energy_i'] - df_f['DeltaE_tot']
                    proj_energycurr = df_f['Energy_i'] - df_f['DeltaE_tot']
                    de_tot_layer.append(df_f['DeltaE_tot'])

                # Define two new dataframes to be the final outputs:
                df_out = pd.DataFrame()
                absorberdf_out = pd.DataFrame()
                absorberdf_out['Layer'] = individuallayers
                # We'll use list comprehension to assign values to the different absorber parameters or write N/A
                # if the parameters don't make sense for the specific state.
                absorberdf_out['State'] = ['Gas' if i else 'Solid' for i in isgas]
                absorberdf_out['Density (g/cm^3)'] = [den[i] if not isgas[i] else 'N/A' for i in range(len(isgas))]
                absorberdf_out['Thickness (mg/cm^2)'] = [thk[i] if not isgas[i] else 'N/A' for i in range(len(isgas))]
                absorberdf_out['Pressure (Torr)'] = [prs[i] if isgas[i] else 'N/A' for i in range(len(isgas))]
                absorberdf_out['Length (cm)'] = [leng[i] if isgas[i] else 'N/A' for i in range(len(isgas))]

                print("Absorber Data: ")
                print(absorberdf_out)

                # We also want to put the particle results into a nice dataframe, so that's done here:
                df_out['Isotope'] = proj_symblistnp
                df_out['Initial Energy (MeV)'] = proj_ei
                df_out['Energy lost (MeV)'] = proj_ei - proj_energycurr
                df_out['Final Energy (MeV)'] = proj_ei - df_out['Energy lost (MeV)']
                df_out['Estimated Straggling (MeV)'] = df_f['E_strag_FWHM']
                print("\nProjectile Data: ")
                print(df_out)

                if len(individuallayers) > 1:
                    layeropt = int(input("\nEnter 1 to print the energy lost in each layer, or 0 to continue: "))
                    if layeropt == 1:

                        df_out_indlayer = pd.DataFrame()
                        df_out_indlayer['Isotope'] = proj_symblistnp
                        df_out_indlayer['Initial Energy (MeV)'] = proj_ei
                        for i in range(len(individuallayers)):
                            stringgy = str(i+1) + ": " + individuallayers[i] + " dE (MeV)"
                            df_out_indlayer[stringgy] = de_tot_layer[i]
                        print(df_out_indlayer)
                else:
                    input("\nPress ENTER to continue: ")

        if option == 4:
            # Option 4 finds the thickness/length/pressure needed in a specific layer to get some final energy for
            # a stopping particle.
            proj_ei4 = []
            proj_z4 = []
            proj_a4 = []

            # These are the upper bounds used for the different absorber parameters and can be changed if need be.
            solidupper = 1000  # mg/cm2
            gaslengthup = 1000  # cm
            gaspressup = 2000  # Torr

            # Need the layers defined before we can run this part of the code.
            if len(individuallayers) == 0:
                print("Please define the stopping layers using Option 2 first.")
            else:
                # We can only do this for one projectile at once. It'd be pretty slow to do more.
                proj = input("\nEnter the isotope to use in the thickness estimation (i.e. 12C): ")
                # Split the A and symbol like above
                projsplit = re.split('(\d+)', proj)

                # The 4 at the end is just to distinguish these for option 4. Splitting and finding Z works the same
                proj_a4.append(int(projsplit[1]))
                symbmask = symb == projsplit[2]
                proj_z4.append(int(zarr[symbmask][0]))

                # The user must define the initial energy and the final energy they want the particle to have
                proj_ei4.append(float(input("\nEnter the initial energy in MeV: ")))
                efin = float(input("\nEnter the final energy in MeV: "))

                proj_z4 = np.array(proj_z4)
                proj_a4 = np.array(proj_a4)
                proj_ei4 = np.array(proj_ei4)

                # The calculations we do rely on the ratio of the calculated energy and the final energy that the
                # user inputs.
                efin_calc = 0
                eratio = efin_calc / efin

                # If there are multiple layers we need the user to say which layer parameters they want to vary:
                if len(individuallayers) > 1:
                    layer = int(input("\nWhich layer would you like to vary?: ")) - 1
                else:
                    layer = 0

                # The initial value of err is the upper limit of the parameter that the user wants to vary. err then
                # gets added (subtracted) to (from) the parameter until the calculation converges.
                if isgas[layer]:
                    # If the layer is a gas the user can vary the pressure or length. For solid layers, you can only
                    # vary the thickness and not the density.
                    lenprsq = int(input("\nWould you like to vary the gas length (1) or gas pressure (2)?: "))
                    if lenprsq == 1:
                        err = gaslengthup
                    else:
                        err = gaspressup
                else:
                    err = solidupper
                    lenprsq = 29

                # The broken boolean allows the program to print the results if the calculation converges.
                broken = False

                # We want the final energy to be found to within 1%, so the while loop ends when the energy ratio
                # reaches that threshold.
                while eratio < .99 or eratio > 1.01:
                    # calculate the energy:
                    df_f = desorb(proj_z4, proj_a4, proj_ei4, z_absorber, a_absorber, numa_absorber, isgas,
                                  den, thk, prs, leng)
                    efin_calc = proj_ei4[0] - df_f['DeltaE_tot'][0]
                    eratio = efin_calc / efin

                    # if the calculated energy is bigger than the desired energy, we want to add err to the parameter
                    if eratio > 1:
                        if lenprsq == 29:
                            thk[layer] = thk[layer] + err
                        if lenprsq == 1:
                            leng[layer] = leng[layer] + err
                        else:
                            prs[layer] = prs[layer] + err
                    # if the calculated energy is smaller, we want to subtract err
                    if eratio < 1:
                        if lenprsq == 29:
                            thk[layer] = thk[layer] - err
                            if thk[layer] < 0:
                                thk[layer] = 0.01
                        if lenprsq == 1:
                            leng[layer] = leng[layer] - err
                            if leng[layer] < 0:
                                leng[layer] = 0.01
                        else:
                            prs[layer] = prs[layer] - err
                            if prs[layer] < 0:
                                prs[layer] = 0.01

                    # every pass through the loop we divide err by 2
                    err = err / 2

                    # if err gets small enough, we want to stop the program
                    if err < 0.01:
                        print("No results found for the parameters specified...")
                        broken = True
                        break

                # if broken is False still (err > 0.01), we can print the final result:
                if not broken:
                    print("\nThis is only an estimate, and the desired final energy was found to within ~1%.\n")
                    if lenprsq == 29:
                        print("The thickness you seek is: " + str(round(thk[layer], 2)) + " mg/cm^2...")
                    elif lenprsq == 1:
                        print("The gas length you seek is: " + str(round(leng[layer], 2)) + " cm...")
                    else:
                        print("The pressure you seek is: " + str(round(prs[layer], 2)) + " Torr...")

                input("Press ENTER to continue.")




