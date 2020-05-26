import numpy as np


class Common:

    # Each of these in the original desorb code have 19 entries. Define them blank for now
    aISG, aINN = np.zeros(19), np.zeros(19)
    aDEN, aTHK, aPRS, aXLN, aARDEN = np.zeros(19), np.zeros(19), np.zeros(19), np.zeros(19), np.zeros(19)
    loste = np.zeros(19)

    # Each of thse are 19x4
    aZNUMB, aANUMB, aELNUM, aCONCN, aPRDEN, aPRTHK = np.zeros((19, 4)), np.zeros((19, 4)), np.zeros((19, 4)), \
                                                     np.zeros((19, 4)), np.zeros((19, 4)), np.zeros((19, 4))

    def desorb(self, ianz, zp, ap, ep):

        # each of these have 4 entries
        aZNUMBW, aANUMBW, aELNUMW, aCONCNW = np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4)
        aPRDENW, aPRTHKW = np.zeros(4), np.zeros(4)

        # these have 19
        aE, aDE, aZmem = np.zeros(19), np.zeros(19), np.zeros(19)
        aTOUT, aTOUTE = np.zeros(19), np.zeros(19)

        # first is 50 x 500 x 2, second is 50
        sEptable, aEmintabz = [[[]]], []
        INW, ISTORE, INS, iArrayValue, aIANZV = 0, 0, 0, 0, []

        io1, IO2, IO3, io0, iopt, ianzi, ianzide, ISTAT, XNS, DEI = 9, 11, 12, 1, 1, 2, 1, 0, 2.0, 0.0

        # the mass table is to be used only for iopt = 5, 6 use atomic masses to average for isotipic composition.
        # taken from Formulas Facts and Constants, H.J.Fischbeck K.H.Fischbeck.Springer - Verlag 1987 2nd ed, pages
        # 164 - 183.

        amass = [1.01, 4.00, 6.94, 9.01, 10.81, 12.01, 14.01, 16.00, 19.00,
                 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07, 35.45, 39.95,
                 39.10, 40.08, 44.96, 47.88, 50.94, 52.00, 54.94, 55.85, 58.93,
                 58.69, 63.55, 65.39, 69.72, 72.59, 74.92, 78.96, 79.90, 83.80,
                 85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 98., 101.07, 102.91,
                 106.42, 107.87, 112.41, 114.82, 118.71, 121.75, 127.60, 126.90,
                 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 147.,
                 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93,
                 173.04]

        # For X > 70 you are in trouble!

        # IOPT = 1 - SUPPLY ENERGY OF PARTICLE ENTERING
        #            THE ABSORBER ARRAY AND GET LOSS AND
        #            RANGES
        # IOPT = 2 - SUPPLY TARGET, PROECTILE AND EJECTILE
        #           INFO. AND THEN GO THROUGH ABSORBER
        #           SANDWICH
        # IOPT = 3 - CALCULATE ENERGY DEPOSITS IN DETECTOR
        #           DETECTOR DIMENSIONS ARE STANDARD AND
        #           THE VARIABLE -'IDET' - CHOOSES BETWEEN
        #           VARIETY OF AVAILABLE DETECTORS
        # IOPT = 4 - FINDS MAXIMUM ENERGY THAT CAN BE STOPPED IN
        #           IANZ ELEMENTS OF THE SANDWICH FOR GIVEN
        #           ZP, AP.
        #           WHEN CALCULATION IS FINISHED, THE PROGRAM READS
        #           IN NEW VALUES OF ZP, AP AND RESTARTS. TO END
        #           THE PROGRAM, GIVE ZP < 0.
        #           IN ORDER TO HELP THE SPEED OF THE PROGRAM,
        #           GIVE THE PARTICLE'S "Z" IN increasing  ORDER.
        # IOPT = 5 - STORES ARRAYS OF Edet AS A FUNCTION OF INCIDENT
        #           ENERGY AND THE PARTICLE'S ID (Z,A)
        #           ARRAY LABELED  eptable(Z-Zth,Einc,ipunch)
        #           ipunch = 1 stopped,  = 2 punched through
        #           Einc = E(incident)/detable
        #           Zth  = lowest Z considered - 1

        # ianz is the number of elements in absorber sandqich - the particle deoposits all its energy in these layers
        # ianzi is the index of last layer in which the energy of the particle is not recorded
        # this unrecorded energy is used in the DST production in two modes: when makint DST taoe from data:
        # Since the detector records only deposited energy the output tables are used to correct for this deficiency
        # and for a given charge and mass extrapolate the measured energy to the target exit energy when making DST
        # tape from model calcs: The lost energy is a source of broadening of the energy spectra due to straggling
        # this smearing is estimated and superimposed on the calc spectra

        # ianzide is the element number of the dE calc
        # zth is the threshold Z for incident energy tble calc (iopt=5,6)
        # zmax is the maximum Z for table calc
        # detable is the energy step size used for array storage
        # emin is the starting incident energy for the table
        # emax is the maximum incident energy for the table calc
        # EP is ignored when iopt = 5 or 6

        for ISN in range(ianz):
            if Common.aISG[ISN] != 1:
                aTOUT[ISN] = Common.aTHK[ISN] / (Common.aDEN[ISN] * 1000.0)
                aTOUTE[ISN] = aTOUT[ISN] / 2.54

        for i in range(ianz):
            INNW = int(Common.aINN[i])
            for j in range(INNW):
                aANUMBW[j] = Common.aANUMB[i][j]
                aZNUMBW[j] = Common.aZNUMB[i][j]
                aELNUMW[j] = Common.aELNUM[i][j]
                aCONCNW[j] = Common.aCONCN[i][j]

            DENW = Common.aDEN[i]
            XLNW = Common.aXLN[i]
            PRSW = Common.aPRS[i]
            THKW = Common.aTHK[i]

            if Common.aISG[i] == 1:
                INNW, THKW, DENW, PRSW, XLNW = Common.setabg(1, INNW, aANUMBW, aZNUMBW, aELNUMW, aPRTHKW, aCONCNW,
                                                             THKW, aPRDENW, DENW, PRSW, XLNW)
                Common.aDEN[i] = 0.0
                Common.aTHK[i] = 0.0

                for j in range(INNW):
                    Common.aPRDEN[i][j] = aPRDENW[j]
                    Common.aPRTHK[i][j] = aPRTHKW[j]
                    Common.aDEN[i] = Common.aDEN[i] + Common.aPRDEN[i][j]
                    Common.aTHK[i] = Common.aTHK[i] + Common.aPRTHK[i][j]
            elif Common.aISG[i] != 1:
                INNW, THKW, DENW = Common.setabs(1, INNW, aANUMBW, aZNUMBW, aELNUMW, aPRTHKW, THKW,
                                                             aPRDENW, DENW)
                for j in range(INNW):
                    Common.aPRDEN[i][j] = aPRDENW[j]
                    Common.aPRTHK[i][j] = aPRTHKW[j]


        zp = zp + 0.5
        izp = int(zp)
        zp = zp-0.5
        # EI energy in
        EI = ep
        XUPDN = -1.0
        EPS = 0.0001
        I1STPASS = 1

        ipunch = 2

        # skipping this part where it writes to the file

        # XNS initial number of the intervals for integration of dE
        # DEI energy out, energy in (<0 for energy loss)
        # E[i] Energy left after ith element (EP=DE[0]-DE[1]...)
        # if particle stopped in detector this is equal to energy lost in remaining layers

        i = 0
        i, ISTAT, XUPDN, XNS, EPS, ap, zp, EI, DEI = Common.ads(1, i, XUPDN, XNS, EPS, ap, zp, EI, DEI, ISTAT)
        EIOLD = EI
        aDE[i] = DEI
        aE[i] = EI+DEI
        EI = aE[i]
        XNS = XNS + 0.1
        INS = int(XNS)

        Common.loste[i] = -1.0 * aDE[i]

        iArrayValue = ianz + 1
        Common.loste[iArrayValue] = aE[ianz]

        ISTORE = i

        if (EI < EPS):
            ipunch = 1


    def setabs(self, pINW, A, Z, AN, T, pTH, D, pDN):
        # Function for setting up composite absorber data
        # partial densities and thicknesses

        INW = pINW
        DN = pDN
        TH = pTH

        AW = 0

        for i in range(INW):
            AW = AW + (A[i]*AN[i])
        for i in range(INW):
            AN[i] = A[i] * AN[i] / AW
            T[i] = TH * AN[i]
            D[i] = DN * AN[i]

        pINW = INW
        pDN = DN
        pTH = TH

        return pINW, pTH, pDN

    def setabg(self, pINW, A, Z, AN, CN, T, pTH, D, pDN, pPR, pXL):
        # subroutine for setting up composite absorber data for gaseous layers

        TH = pTH
        DN = pDN
        PR = pPR
        XL = pXL
        INW = pINW

        P = PR / 760  # I think PR gets converted to atm here?
        X = XL / 22.4
        AWW = 0.0
        AW = 0.0

        for i in range(INW):
            AW = AW + (A[i] * AN[i])
            AWW = AWW + (A[i] * AN[i] * CN[i])
            T[i] = P * X * A[i] * AN[i] * CN[i]
            D[i] = T[i] / XL

        pTH = TH
        pDN = DN
        pPR = PR
        pXL = XL
        pINW = INW

        return pINW, pTH, pDN, pPR, pXL

    def ads(self, pI1, psign, pXN1, pEPS, pA, pZ, pE, pDEE, pISTAT):
        # Subroutine for energy loss calculations
        # call DEDX for stopping power calculations

        I1 = pI1
        ISTAT = pISTAT
        sign = psign
        XN1 = pXN1
        EPS = pEPS
        A = pA
        Z = pZ
        E = pE
        DEE = pDEE
        DE = 0.0
        DEX = 0.0
        # NI number of integrations for energy loss
        EH = E
        XN1 = XN1 + 0.001
        N1 = int(XN1)
        DEDNEXT = 0.0

        k = 0
        while k < N1:
            k = k + 1
            #print(N1)
            J1 = int(Common.aINN[I1])
            ISGW = Common.aISG[I1]
            I = I1
            j = 0
            while j < J1:
                AX = Common.aANUMB[I][j]
                ZX = Common.aZNUMB[I][j]
                FX = Common.aPRTHK[I][j] / XN1
                #print("FX " + str(FX))
                DENST = Common.aPRDEN[I][j]
                VH = np.sqrt(2.13e-3 * EH / A)

                Z, A, ZX, AX, DENST, EH, VH, ISGW, DEX, DE = Common.dedx(1, Z, A, ZX, AX, DENST, EH, VH, ISGW,
                                                                             DEX, DE)
                EH = EH + DE * sign * FX
                #print(DE)
                #if k == 1:
                    #print(EH, DE, sign, FX)
                if EH <= 0.0:
                    if k <= 2:

                        N1 = N1 * 2
                        XN1 = float(N1)
                        j = -1
                        k = 0
                        DEDNEXT = 0.0
                        EH = E

                    ISTAT = -1
                    break
                if k <= 2:
                    DEDNEXT = DEDNEXT + DE * FX

                j = j + 1
                #print(k)
            if k == 1:
                DED1ST = DEDNEXT
                DEDNEXT = 0.0
            if k == 2:
                DDD = DED1ST - DEDNEXT
                #print(DED1ST, DEDNEXT)
                if DDD < 0.0:
                    DDD = -DDD

                DDS = DED1ST + DEDNEXT
                DDR = DDD / DDS

                if DDR > EPS:
                    print("DDR " + str(DDR))
                    print("EPS " + str(EPS))
                    #print("derp")
                    N1 = N1 * 2
                    print("yellllo")
                    XN1 = float(N1)
                    j = -1
                    k = 0
                    DEDNEXT = 0.0
                    EH = E
                    print(EH)



        ISTAT = 0
        DEE = EH - E
        #print(DEE)

        pI1 = I1
        pISTAT = ISTAT
        psign = sign
        pXN1 = XN1
        pEPS = EPS
        pA = A
        pZ = Z
        pE = E
        pDEE = DEE

        return pI1, pISTAT, psign, pXN1, pEPS, pA, pZ, pE, pDEE

    def dedx(self, pZ1, pA1, pZ2, pA2, pRHO, pENR, pV, pIFG, pDEDXHI, pDEDXTO):
        Z1 = pZ1
        A1 = pA1
        Z2 = pZ2
        A2 = pA2
        RHO = pRHO
        ENER = pENR
        V = pV
        DEDXHI = pDEDXHI
        DEDXTO = pDEDXTO
        IFG = pIFG

        A2SAV = 0.0
        Z2SAV = 0.0

        # Program calculates the differential energy loss dE/dX in solid targets using a semiempirical formula deduced
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

        if (IFG == 1):
            RHO = 1

        XI = V * V / Z2

        # Absorber function
        # G(XI) = Y(EXP) - Y(Theory) is deduced from experimental energy loss measurements

        if A2 != A2SAV and Z2 != Z2SAV:
            A2SAV = A2
            Z2SAV = Z2

            # FY is function Y
            FY = 54721.0 * (1.0 + 5.15e-2 * np.sqrt(A2 / RHO) - np.exp(-0.23 * Z2))

            if IFG == 1:
                FY = 54721.0 * (1.35 - np.exp(Z2*(-0.13 + 0.0014 * Z2)))

        # G(XI) is the derivation of a gaussian with variable height H(Z2)
        if Z2 <= 26.0:
            G1 = 19.84 * np.exp(-0.17 * (Z2 - 4.25)*(Z2 - 4.25))
        elif Z2 > 26.0:
            G1 = 0.000001

        if Z2 < 38.0:
            G2 = 17.12 * np.exp(-0.12 * (Z2 - 11.63) * (Z2 - 11.63))
        elif Z2 > 38.0:
            G2 = 0.0000001

        G3 = 7.95 * np.exp(-0.015 * (Z2 - 30.2) * (Z2 - 30.2))
        G4 = 5.84 * np.exp(-0.022 * (Z2 - 48.63) * (Z2 - 48.63))
        G5 = 7.27 * np.exp(-0.005 * (Z2 - 73.06) * (Z2 - 73.06))
        HZ2 = (9.0 - (G1 + G2 + G3 + G4 + G5)) * 1.32e-5

        Z2ZWD = np.cbrt(Z2) * np.cbrt(Z2)

        # Multiplication factors of G(XI)

        FG = 1.2e-4 * Z2 * Z2 + (2.49e-2 * A2 / RHO)

        if IFG == 1:
            FG = 1.3 / (1.0 + np.exp(3.0 - (Z2 / 5.0)))

        ALEFG = np.log((2.7e-5) / FG)

        # Calculation of G(XI)

        GXI = 0.0

        if XI >= 1.0e-9 and XI <= 5.0e-4:
            SQXI = np.sqrt(XI)
            C = (2.0 / Z2) * (SQXI / (1.0 + 1.0e4 * SQXI))

            if (IFG == 1):
                C = C / 2.0

            FG0 = 1.0 / (1.0 + (XI * 10000.0) * (XI * 10000.0) * (XI * 10000.0))
            AL = np.log(XI) - ALEFG
            GXI = (C - HZ2 * AL * np.exp(-0.32 * AL * AL)) * FG0

        # Calculation of Y(XI)
        Y = 3.3e-4 * np.log(1.0 + (XI * FY)) + GXI

        # Energy loss of heavy ions
        # Effective charge

        VV0 = V * 137.0
        FV = 1.0

        if (V >= 0.62):
            FV = 1.0 - np.exp(-VV0)

        AZ1 = np.log(1.035 - 0.4 * np.exp(-0.16 * Z1))

        QQ = V / Z1**0.509
        GHI = Z1
        VZ1 = (-116.79 - 3350.4 * QQ) * QQ
        if VZ1 > -85.2:
            GHI = Z1 * (1.0 - np.exp(VZ1))
        if Z1 > 2.0:
            GHI = Z1 * (1.0 - np.exp(FV * AZ1 - 0.879 * (VV0 / Z1**0.65)))

        # Effective charge of protons and aphas
        # Electronic energy loss DEDXHI
        DEDXHI = GHI * GHI * Z2 * Y / (A2 * V * V)
        # nuclear energy loss DEDXNU
        ZA = np.sqrt(np.cbrt(Z1) * np.cbrt(Z1) + Z2ZWD)

        EPS = 3.25e4 * A2 * ENER / (Z1 * Z2 * (A1 + A2) * ZA)
        SIGMAN = 1.7 * np.sqrt(EPS) * np.log(EPS + 2.1718282) / (1.0 + 6.8 * EPS + 3.4 * np.sqrt(EPS)**3)

        DEDXNU = SIGMAN * 5.105 * Z1 * Z2 * A1 / (ZA * A2 * (A1 + A2))

        # Total energy loss
        DEDXTO = DEDXHI + DEDXNU


        pZ1 = Z1
        pA1 = A1
        pZ2 = Z2
        pA2 = A2
        pRHO = RHO
        pENER = ENER
        pV = V
        pDEDXHI = DEDXHI
        pDEDXTO = DEDXTO
        pIFG = IFG

        return pZ1, pA1, pZ2, pA2, pRHO, pENER, pV, pIFG, pDEDXHI, pDEDXTO

def defprta(energy, thickness, density, zstop, astop, zproj, aproj):
    Target  = Common()

    Target.aPRS[0] = 10.0
    Target.aXLN[0] = 10.0
    Target.aCONCN[0] = 1.0

    Target.aISG[0] = 0
    Target.aINN[0] = 1
    Target.aDEN[0] = density  #this is the density of the stopping medium
    Target.aZNUMB[0][0] = zstop  #this is the Z of the stopping medium
    Target.aANUMB[0][0] = astop  # this is the A of the stopping medium
    Target.aELNUM[0][0] = 1
    Target.loste[0] = 0
    Target.aTHK[0] = thickness  # thickness is g/cm2

    Target.desorb(1, zproj, aproj, energy) # zproj and aproj are the projectile A and A


    elost = Target.loste[0]

    print(Target.loste[0])

    return elost


if __name__ == "__main__":

    # zp = input("Z of projectile ")
    # ap = input("A of projectile ")
    # ep = input("energy of projectile ")
    # zt = input("Z of target")
    # at = input("A of target")
    # dt = input("Density of target in g/cm^3")
    # tt = input("Thickness of target in g/cm^2")

    de = defprta(14.2, 137.16, 2.7, 13, 27, 1, 3)



