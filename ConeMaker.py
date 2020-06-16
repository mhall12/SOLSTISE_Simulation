import pickle


def makecone():

    # Code for defining cone parameters and fitting parameters. The information from the text file is then
    # used in the main simulation code to define the cone.

    # Distance from the nozzle to the cone in Inches
    conenozzdist = float(input("Enter the distance from the bottom of the nozzle to the top of the cone in inches: "))

    coneheight = float(input("Enter the cone height in inches, from the top of the ISO-100 cylinder"
                             " to the top of the cone: "))

    iso160height = float(input("Enter the distance from the bottom of the nozzle to the top of the "
                               "ISO-160 base in inches: "))

    # Outer diameter of the top of the cone
    conedia = float(input("Enter the outer diameter of the top of the cone in inches: "))

    print("The side of the cone is described by a polynomial (up to order 3) of the form: radius = a_i * y^i, "
          "\n where y and radius are in inches. The y-distance is defined as 0 at the top of the cone and"
          "\nincreases going down the cone.")

    # et polyorder to some high value so that the user gets stuck in the while loop until they enter 3 or less.
    polyorder = 4
    while polyorder > 3 or polyorder < 0:
        polyorder = int(input("Please enter the order of the polynomial you'd like to use to "
                              "describe the sides of the cone (max 3): "))

        if polyorder > 3:
            print("\nERROR: The maximum order of the polynomial is 3... ")

    # Convert conedia into a string and replace the . with a - to put it in the pickle file name
    conediastring = str(conedia)
    conediastring = conediastring.replace('.', '-')

    # Add all the input parameters to a new list called coneparms:
    coneparms = [conenozzdist, coneheight, iso160height, conedia, polyorder]

    # In the for loop, add 0's for the coefficients not being used and let the user
    # input coefficients for the ones that are.
    for i in range(4):
        expon = 3 - i

        if polyorder < expon:
            coneparms.append(0)
        else:
            coneparms.append(float(input("Enter the coefficient on the y^" + str(expon) +
                                         " term (include the sign!): ")))

    fname = "cust_cone_" + str(polyorder) + "_" + conediastring + "in.txt"

    file = open(fname, "w+")

    for i in coneparms:
        file.write(str(i) + "\n")

    return fname


if __name__ == "__main__":

    makecone()
