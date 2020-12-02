from CircleAreaCalc import area_calc, new_height_calc
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def makepipe():
    # main will ask the user questions and call the relevant functions to track the particles from VIKAR

    print("\nNow the SOLSTISE pipe will be defined. It is bounded by two circles of radius r1 and r2,")
    print("as well as two angles \u03B8\u2081 and \u03B8\u2082.")

    try:
        hyn = input("For additional explanation and a helpful figure, enter H, otherwise"
                    " press ENTER to continue.")
    except:
        hyn = "N"

    if hyn == "H" or hyn == "h":
        plt.ion()
        img = mpimg.imread("./Jupyter_Pics/NewPipePic.png")
        imgplot = plt.imshow(img)
        plt.show()
        plt.ioff()

        print("\n Here, you will make a pipe bounded by a circle of radius r2 (the larger the radius, the flatter the "
              "pipe will be on the top. 3-5 meters is a good start), the magnet bore, and two lines at "
              "angles \u03B8\u2081 and \u03B8\u2082 (see figure for reference). The pipe is constrained such that the "
              "cross sectional area is equal to that of an ISO-160 pipe. This is achieved by calculating the height "
              "that the second circle (r2) must be placed, and will be recalculated if the angles of the two lines "
              "are changed. Note that the created pipe cannot be taller than the ISO-160 pipe and if your chosen "
              "angle is too large, you will have to choose another. \n")

    hors = 0

    while hors != 1 and hors != 2 and hors != 333:
        try:
            hors = int(input("\nWill the reaction be occurring in HELIOS (Enter 1) or SOLARIS (Enter 2)? "))
        except ValueError:
            hors = 0

        if hors != 1 and hors != 2 and hors != 333:
            print("\nERROR: Unknown entry...")

    # 333 is for the SOLSTISE chamber/debugging. Shhhhh it's a secret!
    # hors just chooses the magnet bore radius.
    if hors == 1:
        r1 = 0.92 / 2
    elif hors == 333:
        r1 = 0.376174  # 14.81 inches
    else:
        r1 = 0.9 / 2

    plt.ion()
    circle1 = plt.Circle((0, 0), r1, ec='black', color='white')

    # ISO 160 pipe dia in m
    pipedia = 5.83 * 2.54 / 100
    pipecenter = pipedia / 2 - r1
    circlepipe = plt.Circle((0, pipecenter), pipedia / 2, ec='black', fill=False)

    fig, ax = plt.subplots()
    ax.add_artist(circle1)
    ax.add_artist(circlepipe)
    ax.axis('equal')
    ax.set_xlim((-.5, .5))
    ax.set_ylim((-.6, .4))
    # fig.show()

    # matplotlib spits out some warnings and errors that can be ignored, but print the radius and require an input
    # because otherwise the error goes at the end of the next input.
    input("The radius is " + str(r1) + " m. Press ENTER to continue.")
    while True:
        try:
            r2 = float(input("\nEnter the radius of the 2nd circle in meters (or enter 0 for inches): "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    if r2 == 0:
        while True:
            try:
                r2in = float(input("\nEnter the radius of the 2nd circle in inches: "))
                break
            except ValueError:
                print("ERROR: Your entry is invalid. Please try again. \n")

        r2 = r2in * 2.54 / 100

    ang1, h = area_calc(r1, r2)

    # print(h)

    circle2 = plt.Circle((0, h), r2, ec='black', fill=False)
    plt.gcf().gca().add_artist(circle2)
    # fig.show()

    ang1deg = ang1 * 180 / math.pi

    if ang1deg > 270:
        ang1deg = 360 - (ang1deg - 180)

    ang2deg = 360 - (ang1deg - 180)

    print("The circle intersection points are", int(ang1deg), "degrees and", int(ang2deg), "degrees")

    try:
        yn2 = input("\nWould you like to use new angles for the pipe geometry? (Y/N) ")
    except ValueError:
        yn2 = "Y"

    # for now accept only one angle, but later we'll do multiple angles.
    # anglist = []
    pipethickcheck = .15
    iso160pipedia = .148082
    c2height = 0.0
    angin = 0

    if yn2 == "N" or yn2 == "n":
        # anglist.append(ang1)
        angin = ang1deg

        anginrad = angin * math.pi / 180
        angin2 = 360 - (angin - 180)
        angin2rad = angin2 * math.pi / 180

        c2height = h

    else:
        circle2.set_visible(False)
        # It has to be less than 270 because that's the bottom of the circle. Actually, 270 is also dumb because
        # the pipe would block the entire path then.
        print("\nNote: the new angle must be greater than", int(ang1deg), "and less than 270 degrees.")

        while (angin < ang1deg or angin > 270) or pipethickcheck > iso160pipedia:
            angin = float(input("\nInput the new angle (in degrees): "))

            if angin < ang1deg:
                print("Error: The angle must be greater than", int(ang1deg), "and less than 270 degrees.")
            elif angin > ang1deg:
                c2height = new_height_calc(r1, r2, angin)
                pipethickcheck = r1 - r2 + c2height
                if pipethickcheck > iso160pipedia:
                    print("Error: The calculated pipe height is", pipethickcheck,
                          "m and must be less than", iso160pipedia, "m.")
        # anglist.append(angin)

        anginrad = angin * math.pi / 180
        angin2 = 360 - (angin - 180)
        angin2rad = angin2 * math.pi / 180

        circle3 = plt.Circle((0, c2height), r2, ec='black', fill=False)
        plt.gcf().gca().add_artist(circle3)
        x1, y1 = [0, r1 * math.cos(anginrad)], [0, r1 * math.sin(anginrad)]
        x2, y2 = [0, r1 * math.cos(angin2rad)], [0, r1 * math.sin(angin2rad)]
        plt.plot(x1, y1, x2, y2, color='black')
        input("Press ENTER")

    plt.ioff()

    fname = "PipeOut_" + str(int(r2)) + "_" + str(int(angin)) + ".txt"

    params = [r1, r2, c2height, anginrad, angin2rad]

    # input("\nPress ENTER to continue.")

    geodir = "./Geometry_Files/"

    file = open(geodir + fname, "w+")

    # Build the pipe setup file here
    for i in range(5):
        file.write(str(params[i]) + "\n")

    print("An output file named " + fname + " was created in " + geodir + " !")

    file.close()

    # final pipe geometry is defined by anglist[], r1 ,r2, and the output of new_height_calc, which is the
    # height of the 2nd circle above the x-axis (c2height)

    # The equation of an offset circle in cylindrical coordinates rho = d*sin(th) + Sqrt(d^2 sin(th)^2 - d^2 + r^2)
    # for a normal circle, rho = r

    # read in the energy and angle from the text file into a numpy array:

    # sim(r1, r2, c2height, anginrad, angin2rad, vikar_in)

    return fname


if __name__ == "__main__":
    makepipe()



