import re
import numpy as np

def readmass():

    # reac = input("Enter a reaction of the form d(17F,p): ")
    reac = 'd(28Si,p)'

    # firstpos grabs the position of the ( so we can get the target
    firstpos = reac.find('(')

    # target mass and symbol grabbed here
    target = reac[:firstpos]
    # splittarget initialized. It's used later to split the mass from the symbol
    splittarget = []

    # if statements on case the user typed p, d, etc.
    # The else condition splits the mass from the symbol into an array
    # with an empty element in the front (not sure why)
    if target == 'p':
        splittarget = ['', '1', 'H']
    elif target == 'd':
        splittarget = ['', '2', 'H']
    elif target == 't':
        splittarget = ['', '3', 'H']
    elif target == 'a':
        splittarget = ['', '4', 'He']
    else:
        splittarget = re.split('(\d+)',target)

    # beam happens the same as the target:

    secondpos = reac.find(',')

    beam = reac[(firstpos+1):secondpos]

    splitbeam = []

    if beam == 'p':
        splitbeam = ['', '1', 'H']
    elif beam == 'd':
        splitbeam = ['', '2', 'H']
    elif beam == 't':
        splitbeam = ['', '3', 'H']
    elif beam == 'a':
        splitbeam = ['', '4', 'He']
    else:
        splitbeam = re.split('(\d+)', beam)

    # ejectile happens the same as the last two

    thirdpos = reac.find(')')

    ejectile = reac[(secondpos+1):thirdpos]

    splitejectile = []

    if ejectile == 'p':
        splitejectile = ['', '1', 'H']
    elif ejectile == 'd':
        splitejectile = ['', '2', 'H']
    elif ejectile == 't':
        splitejectile = ['', '3', 'H']
    elif ejectile == 'a':
        splitejectile = ['', '4', 'He']
    else:
        splitejectile = re.split('(\d+)', ejectile)

    # recoil is a bit different, since the user does not need to specify it. Instead, it has to
    # be calculated from the inputs that were given.

    # File name is specified here for the masses. It has the structure:
    # Z A Symbol Mass_MeV Mass_u
    infile = 'masses.txt'

    # Generate a numpy array from the mass file. The dtype is needed because otherwise the symbols try to
    # get read in as numbers and fill the arrays with NaN
    data = np.genfromtxt(infile, delimiter='\t', dtype = 'unicode')

    # Grab each column and put it into its own array to make it a little easier
    z = data[:, 0]
    a = data[:, 1]
    symb = data[:, 2]
    massMeV = data[:, 3]
    massu = data[:, 4]

    # These generate a mask for the target, beam, and ejectile that are True on the row where the
    # specific isotope is located, and false on all the other lines. The masks will be used to grab
    # the masses later. The string.lower() command takes the isotopic symbol and makes it lower case
    # so it can be matched to the one in the file. So, we're matching the mass and the symbol.
    masktarget = (splittarget[2].lower() == symb) & (a == splittarget[1])
    maskbeam = (splitbeam[2].lower() == symb) & (a == splitbeam[1])
    maskejectile = (splitejectile[2].lower() == symb) & (a == splitejectile[1])

    # The Z for the target, beam, and ejectile are found here by masking the z array. It will now only have one
    # element, which is an integer.
    ztarget = int(z[masktarget])
    zbeam = int(z[maskbeam])
    zejectile = int(z[maskejectile])

    # The proton number of the recoil is found here.
    zrecoil = ztarget + zbeam - zejectile
    # The mass A of the recoil is found here.
    arecoil = int(splittarget[1]) + int(splitbeam[1]) - int(splitejectile[1])

    # A recoil mask is made using the Z and A this time instead of symbol and A.
    maskrecoil = (zrecoil == z.astype(np.int)) & (arecoil == a.astype(np.int))

    # Now, the masses can be returned for the four particles by applying the mask to the mass arrays. Easy peasy.
    masses = [float(massu[masktarget][0]), float(massu[maskbeam][0]), float(massu[maskejectile][0]),
              float(massu[maskrecoil][0])]

    return masses


readmass()