import pickle


def makecone():

    # Code for defining cone parameters and fitting parameters. The information from the text file is then
    # used in the main simulation code to define the cone.

    # Distance from the nozzle to the cone in Inches
    conenozzdist = float(input("Enter the distance from the bottom of the nozzle to the top of the cone in inches: "))

    # Outer diameter of the top of the cone
    conedia = float(input("Enter the outer diameter of the top of the cone in inches: "))

    # Distance from the top of the ISO 100 base to the top of the cone
    coneheight = float(input("Enter the cone height in inches, from the top of the ISO base to the top of the cone: "))

    print("The side of the cone is described by a polynomial (up to order 3) of the form radius = a * height + b, "
          "\n where height and radius are in inches. ")

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
    coneparms = [conenozzdist]
    coneparms.append(conedia)
    coneparms.append(coneheight)

    coneparms.append(polyorder)
    # In the for loop, add 0's for the coefficients not being used and let the user
    # input coefficients for the ones that are.
    for i in range(4):
        expon = 3 - i

        if polyorder < expon:
            coneparms.append(0)
        else:
            coneparms.append(float(input("Enter the coefficient on the x^" + str(expon) + " term (inlude the sign!): ")))

    fname = "cust_cone_" + str(polyorder) + "_" + conediastring + "in.pkl"

    # Write coneparms to a pickle.
    with open(fname, 'wb') as f:
        pickle.dump(coneparms, f)

    return fname

if __name__ == "__main__":
    makecone()