import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def makenozz():

    # Code for defining nozzle parameters. The information from the text file is then
    # used in the main simulation code to define the nozzle for shadowing..

    hp = input("To open a figure that will help explain some of the required dimensions, enter H. Otherwise, press "
               "ENTER to continue.")

    if hp == "H" or hp == "h":
        plt.ion()
        img = mpimg.imread("./Jupyter_Pics/NozzleFig.png")
        imgplot = plt.imshow(img)
        plt.show()
        plt.ioff()

    # Distance below the nozzle that the reaction occurs (output in mm):
    while True:
        try:
            reacdist = float(input("Enter the distance that the reactions occur below the nozzle in mm (Enter 0 for "
                                   "inches): "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if reacdist == 0:
        while True:
            try:
                reacdistin = float(input("Enter the distance in inches: "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")
        reacdist = reacdistin * 2.54 * 10

    # Nozzle opening diameter in inches
    while True:
        try:
            nozzdia = float(input("Enter the nozzle exhaust opening diameter in mm (Enter 0 for in): "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    if nozzdia > 0:
        nozzdia = nozzdia / 10 / 2.54
    elif nozzdia == 0:
        while True:
            try:
                nozzdia = float(input("Enter the diameter in inches: "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")

    # Angle the outside of the nozzle makes from vertical.
    print("The outside of the nozzle is defined from the exhaust and up, "
          "using a cone shape and then a cylindrical shape.")

    while True:
        try:
            nozzang = float(input("Enter the angle in degrees that the nozzle cone makes with the vertical axis: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    # Output this one in m so we don't need to convert again.
    while True:
        try:
            nozzconelen = float(input("Enter the distance from the cone exhaust to the "
                                    "cylindrical neck portion in mm (Enter 0 for inches): "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    if nozzconelen == 0:
        while True:
            try:
                nozzconelen = float(input("Enter the distance in inches: "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")

        nozzconelen = nozzconelen * 2.54 / 100
    elif nozzconelen > 0:
        nozzconelen = nozzconelen / 1000

    nozzcyldia = 2 * nozzconelen * math.tan(nozzang * math.pi/180) + nozzdia * 2.54 / 100
    print("The nozzle cylinder diameter was calculated to be: ", nozzcyldia * 1000, " mm.")

    try:
        chyn = input("Would you like to change the calculated cylinder diameter? [Y/N] ")
    except ValueError:
        chyn = "N"

    if chyn == "Y" or chyn == "y":
        while True:
            try:
                nozzcyldia = float(
                    input("Enter the diameter of the nozzle cylinder (above the cone) in mm (Enter 0 for inches): "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")
        if nozzcyldia == 0:
            while True:
                try:
                    nozzcyldia = float(input("Enter the diameter in inches: "))
                    break
                except ValueError:
                    print("ERROR: Your entry is invalid. Please try again. \n")
            nozzcyldia = nozzcyldia * 2.54 / 100
        elif nozzcyldia > 0:
            nozzcyldia = nozzcyldia / 1000

    # Nozzle cylinder height is in m
    while True:
        try:
            nozzcylh = float(input("Enter the distance from the bottom of the nozzle cylinder to the bottom of the "
                                   "nozzle holder in mm (Enter 0 for inches): "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if nozzcylh == 0:
        while True:
            try:
                nozzcylh = float(input("Enter the distance in inches: "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")
        nozzcylh = nozzcylh * 2.54 / 100
    elif nozzcylh > 0:
        nozzcylh = nozzcylh / 1000

    print("Now, we'll define the nozzle holder. The holder is made of a cylinder and a box portion.")

    # Diameter of the nozzle holder is also in m
    while True:
        try:
            nozzholderdia = float(input("Enter the diameter of the nozzle holder cylinder in mm (default is 19.05 mm) "
                                        "or enter 0 for inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if nozzholderdia > 0:
        nozzholderdia = nozzholderdia / 1000
    elif nozzholderdia == 0:
        while True:
            try:
                nozzholderdia = float(input("Enter the diameter in inches (default is 0.75 in): "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")
        nozzholderdia = nozzholderdia * 2.54 / 100

    # default holder height is 1.15 in, should be in m

    while True:
        try:
            nozzholderheight = float(input("Enter the height of the nozzle holder cylinder in mm (default is 29.2 mm) "
                                           "or enter 0 for inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if nozzholderheight == 0:
        while True:
            try:
                nozzholderheight = float(input("Enter the height in inches (default is 1.15 in): "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")
        nozzholderheight = nozzholderheight * 2.54 / 100
    elif nozzholderheight > 0:
        nozzholderheight = nozzholderheight / 1000

    while True:
        try:
            nozzboxstart = float(input("Enter the distance between the bottom of the holder cylinder \n"
                                       "and the bottom of the holder box in mm (default is 5.08 mm) or enter 0 for "
                                       "inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if nozzboxstart == 0:
        while True:
            try:
                nozzboxstart = float(input("Enter the distance in inches (default is 0.2 in): "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")
        nozzboxstart = nozzboxstart * 2.54 / 100
    elif nozzboxstart > 0:
        nozzboxstart = nozzboxstart / 1000

    while True:
        try:
            nozzboxh = float(input("Enter the height of the box in mm (default is 19.05 mm) or enter 0 for inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if nozzboxh == 0:
        while True:
            try:
                nozzboxh = float(input("Enter the distance in inches (default is 0.75 in): "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")

        nozzboxh = nozzboxh * 2.54 / 100
    elif nozzboxh > 0:
        nozzboxh = nozzboxh / 1000

    while True:
        try:
            nozzboxl = float(input("Enter the (side-to-side) width of the box in mm (default is 43.142 mm) or "
                                   "enter 0 for inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if nozzboxl == 0:
        while True:
            try:
                nozzboxl = float(input("Enter the length in inches (default is 1.7 in): "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")
        nozzboxl = nozzboxl * 2.54 / 100
    elif nozzboxl > 0:
        nozzboxl = nozzboxl / 1000

    while True:
        try:
            nozzboxw = float(input("Enter the (back-to-front) length of the box in mm (default is 14.199 mm) or "
                                   "enter 0 for inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")
    if nozzboxw == 0:
        while True:
            try:
                nozzboxw = float(input("Enter the distance in inches (default is 0.559 in): "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")

        nozzboxw = nozzboxw * 2.54 / 100
    elif nozzboxw > 0:
        nozzboxw = nozzboxw / 1000

    params = [reacdist, nozzdia, nozzang, nozzconelen, nozzcyldia, nozzcylh, nozzholderdia, nozzholderheight,
              nozzboxstart, nozzboxh, nozzboxl, nozzboxw]

    geodir = "./Geometry_Files/"
    fname = "cust_nozz_" + str(int(nozzang)) + "deg.txt"

    file = open(geodir + fname, "w+")

    for i in params:
        file.write(str(i) + "\n")

    print("A custon nozzle file was created with the name " + fname + "!")

    return fname


if __name__ == "__main__":
    makenozz()