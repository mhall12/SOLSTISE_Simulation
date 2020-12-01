import pickle
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def makecone():

    # Code for defining cone parameters and fitting parameters. The information from the text file is then
    # used in the main simulation code to define the cone.

    print("Here, you will define some of the setup parameters, including the shape of the receiver cone. \n")

    # Distance from the nozzle to the cone in Inches
    while True:
        try:
            conenozzdist = float(input("\nEnter the distance from the bottom of the nozzle to the top of the cone "
                                       "in inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    while True:
        try:
            coneheight = float(input("\nEnter the cone height in inches, from the top of the ISO-100 cylinder"
                                     " to the top of the cone: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    while True:
        try:
            iso160height = float(input("\nEnter the distance from the bottom of the nozzle to the top of the "
                                       "ISO-160 base in inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    # Outer diameter of the top of the cone
    while True:
        try:
            conedia = float(input("\nEnter the outer diameter of the top of the cone in inches: "))
            break
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

    print("\nThe side of the cone is described by a polynomial (up to order 3) of the form: radius = a_i * y^i, "
          "\n where y and radius are in inches. The y-distance is defined as 0 at the top of the cone and"
          "\nincreases going down the cone. \n")

    hp = input("Confused? For a helpful image, enter H now. Otherwise, press ENTER to continue: ")

    if hp == "H" or hp == "h":
        plt.ion()
        img = mpimg.imread("./Jupyter_Pics/ConeDiagram.png")
        imgplot = plt.imshow(img)
        plt.show()
        plt.ioff()

        print("\n\n")
        print("The image demonstrates how the fitting should be done. As mentioned, \n"
              "you can use a polynomial up to order 3 here. \n\n")

    # set polyorder to some high value so that the user gets stuck in the while loop until they enter 3 or less.
    polyorder = 4
    while polyorder > 3 or polyorder < 1:
        try:
            polyorder = int(input("Please enter the order of the polynomial you'd like to use to "
                                  "describe the sides of the cone (max 3): "))
        except ValueError:
            print("ERROR: Your entry is invalid. Please try again. \n")

        if polyorder > 3:
            print("\nERROR: The maximum order of the polynomial is 3... ")
        if polyorder < 1:
            print("\nERROR: The minimum order of the polynomial is 1... ")

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
            while True:
                try:
                    coneparms.append(float(input("Enter the coefficient on the y^" + str(expon) +
                                                 " term (include the sign!): ")))
                    break
                except ValueError:
                    print("\n ERROR: The coefficient you entered is not valid.")

    geodir = "./Geometry_Files/"
    fname = "cust_cone_" + str(polyorder) + "_" + conediastring + "in.txt"

    file = open(geodir + fname, "w+")

    for i in coneparms:
        file.write(str(i) + "\n")

    print("A file named ", fname, " was saved to ", geodir)

    return fname


if __name__ == "__main__":

    makecone()
