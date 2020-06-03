import pickle


def makecone():

    # Code for defining cone parameters and fitting parameters. The information from the text file is then
    # used in the main simulation code to define the cone.

    conenozzdist = float(input("Enter the distance from the bottom of the nozzle to the top of the cone in inches: "))

    conedia = float(input("Enter the outer diameter of the top of the cone in inches: "))

    coneheight = float(input("Enter the cone height in inches, from the top of the ISO base to the top of the cone: "))

    print("The side of the cone is described by a polynomial (up to order 3) of the form radius = a * height + b, "
          "\n where height and radius are in inches. ")

    polyorder = int(input("Please enter the order of the polynomial you'd like to use to "
                          "describe the sides of the cone: "))

    conediastring = str(conedia)
    conediastring = conediastring.replace('.', '-')

    coneparms = []
    coneparms.append(conenozzdist)
    coneparms.append(conedia)
    coneparms.append(coneheight)

    coneparms.append(polyorder)
    for i in range(polyorder+1):
        expon = polyorder - i

        coneparms.append(float(input("Enter the coefficient on the x^" + str(expon) + " term (inlude the sign!): ")))

    fname = "cust_cone_" + str(polyorder) + "_" + conediastring + "in.pkl"

    with open(fname, 'wb') as f:
        pickle.dump(coneparms, f)

    return fname

if __name__ == "__main__":
    makecone()