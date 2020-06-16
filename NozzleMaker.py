def makenozz():

    # Code for defining nozzle parameters. The information from the text file is then
    # used in the main simulation code to define the nozzle for shadowing..

    # Distance below the nozzle that the reaction occurs (output in mm):
    reacdist = float(input("Enter the distance that the reactions occur below the nozzle in mm (Enter 0 for inches): "))
    if reacdist == 0:
        reacdistin = float(input("Enter the distance in inches: "))
        reacdist = reacdistin * 2.54 * 10

    # Nozzle opening diameter in inches
    nozzdia = float(input("Enter the nozzle exhaust opening diameter in mm (Enter 0 for in): "))
    if nozzdia > 0:
        nozzdia = nozzdia / 10 / 2.54
    elif nozzdia == 0:
        nozzdia = float(input("Enter the diameter in inches: "))

    # Angle the outside of the nozzle makes from vertical.
    print("The outside of the nozzle is defined from the exhaust and up, "
          "using a cone shape and then a cylindrical shape.")
    nozzang = float(input("Enter the angle in degrees that the outside of the nozzle makes with the vertical axis: "))

    # Output this one in m so we don't need to convert again.
    nozzconelen = float(input("Enter the distance from the cone exhaust to the "
                              "cylindrical neck portion in mm (Enter 0 for inches): "))
    if nozzconelen == 0:
        nozzconelen = float(input("Enter the distance in inches: "))
        nozzconelen = nozzconelen * 2.54 / 100
    elif nozzconelen > 0:
        nozzconelen = nozzconelen / 1000

    # Cylinder dia for the rest of the nozzle. Also in meters.
    nozzcyldia = float(input("Enter the diameter of the nozzle cylinder (above the cone) in mm (Enter 0 for inches): "))
    if nozzcyldia == 0:
        nozzcyldia = float(input("Enter the diameter in inches: "))
        nozzcyldia = nozzcyldia * 2.54 / 100
    elif nozzcyldia > 0:
        nozzcyldia = nozzcyldia / 1000

    # Nozzle cylinder height is in m
    nozzcylh = float(input("Enter the distance from the bottom of the nozzle cylinder to the bottom of the nozzle "
                           "holder in mm (Enter 0 for inches): "))
    if nozzcylh == 0:
        nozzcylh = float(input("Enter the distance in inches: "))
        nozzcylh = nozzcylh * 2.54 / 100
    elif nozzcylh > 0:
        nozzcylh = nozzcylh / 1000

    print("Now, we'll define the nozzle holder. The holder is made of a cylinder and a box portion.")

    # Radius of the nozzle holder is also in m
    nozzholderdia = float(input("Enter the diameter of the nozzle holder cylinder in mm (default is 19.05 mm) "
                                "or enter 0 for inches: "))
    if nozzholderdia > 0:
        nozzholderdia = nozzholderdia / 1000
    elif nozzholderdia == 0:
        nozzholderdia = float(input("Enter the diameter in inches (default is 0.75 in): "))
        nozzholderdia = nozzholderdia * 2.54 / 100

    # default holder height is 1.15 in, should be in m

    nozzholderheight = float(input("Enter the height of the nozzle holder cylinder in mm (default is 29.2 mm) or "
                                   "enter 0 for inches: "))
    if nozzholderheight == 0:
        nozzholderheight = float(input("Enter the height in inches (default is 1.15 in): "))
        nozzholderheight = nozzholderheight * 2.54 / 100
    elif nozzholderheight > 0:
        nozzholderheight = nozzholderheight / 1000

    nozzboxstart = float(input("Enter the distance between the bottom of the holder cylinder \n"
                               "and the bottom of the holder box in mm (default is 5.08 mm) or enter 0 for inches: "))
    if nozzboxstart == 0:
        nozzboxstart = float(input("Enter the distance in inches (default is 0.2 in): "))
        nozzboxstart = nozzboxstart * 2.54 / 100
    elif nozzboxstart > 0:
        nozzboxstart = nozzboxstart / 1000

    nozzboxh = float(input("Enter the height of the box in mm (default is 19.05 mm) or enter 0 for inches: "))
    if nozzboxh == 0:
        nozzboxh = float(input("Enter the distance in inches (default is 0.75 in): "))
        nozzboxh = nozzboxh * 2.54 / 100
    elif nozzboxh > 0:
        nozzboxh = nozzboxh / 1000

    nozzboxl = float(input("Enter the (side-to-side) length of the box in mm (default is 43.142 mm) or "
                           "enter 0 for inches: "))
    if nozzboxl == 0:
        nozzboxl = float(input("Enter the length in inches (default is 1.7 in): "))
        nozzboxl = nozzboxl * 2.54 / 100
    elif nozzboxl > 0:
        nozzboxl = nozzboxl / 1000

    nozzboxw = float(input("Enter the (back-to-front) width of the box in mm (default is 14.199 mm) or "
                           "enter 0 for inches: "))
    if nozzboxw == 0:
        nozzboxw = float(input("Enter the distance in inches (default is 0.559 in): "))
        nozzboxw = nozzboxw * 2.54 / 100
    elif nozzboxw > 0:
        nozzboxw = nozzboxw / 1000

    params = [reacdist, nozzdia, nozzang, nozzconelen, nozzcyldia, nozzcylh, nozzholderdia, nozzholderheight,
              nozzboxstart, nozzboxh, nozzboxl, nozzboxw]

    fname = "cust_nozz_" + str(int(nozzang)) + "deg.txt"

    file = open(fname, "w+")

    for i in params:
        file.write(str(i) + "\n")

    print("A custon nozzle file was created with the name " + fname + "!")

    return fname


if __name__ == "__main__":
    makenozz()