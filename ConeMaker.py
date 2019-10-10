def main():

    f = open("ConeParams.txt", "w")

    for i in range(100):
        for j in range(100):
            for k in range(100):
                dist = (i+1) / 100 * 10
                dia = (j+1) / 100 * 10
                h = (k+1) / 100 * 10

                ddm = dist * dia
                ddmch = ddm / h
                rat = dist/dia

                if (0.5 < rat < 2.5) and (6 < ddm < 10) and (2.4 < ddmch < 3.8):
                    f.write(str(dist) + " " + str(dia) + " " + str(h) + " " + str(ddm) + " " + str(ddmch) + "\n")




main()