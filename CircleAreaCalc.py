from sympy import nsolve, symbols, acos, sqrt, pi, re, sin, integrate, cos, log, N


def area_calc(r1, r2):

    # convert r1 and r2 from the main program in meters to inches for the calculation
    r1 = r1*100/2.54
    r2 = r2*100/2.54

    # ISO 160 pipe dia in inches
    pipedia = 5.83

    a1 = pi*r1**2
    asm = pi*((pipedia/2)**2)

    # d is the distance between the center points of the two circles
    # d1 is the distance from the origin to the intersection point
    def d1(d):
        return (r1**2-r2**2+d**2)/(2*d)

    # d2 is the distance from the intersection point to the center of circle 2
    def d2(d):
        return d-d1(d)

    # t is defined but not used (used in main program as pipethickcheck) used for debugging at one point
    def t(d):
        return r1 - r2 + d

    # a(d) is the area of the crescent shape
    def a(d):
        return a1 - r1**2 * acos(d1(d)/r1) + d1(d)*sqrt(
            r1**2 - d1(d)**2) - r2**2 * acos(d2(d)/r2) + d2(d) * sqrt(r2**2-d2(d)**2)

    # Solve for the height of the second circle above the x-axis assuming the are is the same
    x = symbols('x', real=True)

    h = re(nsolve(a(x)-asm, x, 10))

    # print(t(h))
    # Solve for the intersection angles assuming h is the height of the 2nd circle above the x-axis
    ang = symbols('ang')

    th = nsolve(r1 - h*sin(ang) - sqrt(h**2*sin(ang)**2 - h**2 + r2**2), ang, 2)

    hmeter = h*2.54/100
    # return this angle and height to the main program
    return th, hmeter


def new_height_calc(r1,r2,th1):
    # th1 and th2 in degrees

    # convert r1 and r2 from the main program in meters to inches for the calculation
    r1 = r1*100/2.54
    r2 = r2*100/2.54

    # ISO 160 pipe dia in inches
    pipedia = 5.83

    # cross sectional area of the ISO 160 pipe
    asm = pi * ((pipedia / 2) ** 2)

    # angles converted to radians for the calculation
    th1rad = th1*pi/180
    th2rad = 2*pi - (th1rad - pi)


    def a2(d):

        # returns the area of a crescent shape that's cut off on the ends. (Integral done in Mathematica)
        return re((r1 ** 2 * th2rad - r2 ** 2 * th2rad - r1 ** 2 * th1rad + r2 ** 2 * th1rad - (d * (-(sqrt(2) * sqrt(
                -d ** 2) * cos(th2rad) * sqrt(-d ** 2 + 2 * r2 ** 2 - d ** 2 * cos(2 * th2rad))) + sqrt(2) * sqrt(
                    -d ** 2) * cos(th1rad) * sqrt(-d ** 2 + 2 * r2 ** 2 - d ** 2 * cos(2 * th1rad)) + 2 * r2 ** 2 * (
                        -log(sqrt(2) * sqrt(-d ** 2) * cos(th2rad) + sqrt(-d ** 2 + 2 * r2 ** 2 - d ** 2 * cos(
                                2 * th2rad))) + log(sqrt(2) * sqrt(-d ** 2) * cos(th1rad) + sqrt(
                                    -d ** 2 + 2 * r2 ** 2 - d ** 2 * cos(2 * th1rad)))))) / (2. * sqrt(
                                        -d ** 2)) + d ** 2 * cos(th2rad) * sin(th2rad) - d ** 2 * cos(th1rad) * sin(
                                                th1rad)) / 2.)

    # Solve the equation above assuming the area is the same as the ISO 160 pipe cross sectional area.
    x = symbols('x', real=True)
    s2 = nsolve(a2(x) - asm, x, 14)

    # convert the answer from inches to meters. s2 is the height above the x-axis of the 2nd circle
    s2meters = s2*2.54/100

    return s2meters

