import math

class Particle:
    def __init__(self, amass, protons, energy, theta, phi):
        self.amass = amass
        self.protons = protons
        self.energy = energy
        self.theta = theta
        self.phi = phi
        self.vlab = 2*energy/(amass*931.5)
        self.vpar = vlab * math.cos(theta)
        self.vperp = 0

    def calc_vlab(self):
        self.vlab = 2*self.energy/(self.amass*931.5)

    def calc_vpar(self):
        self.vpar = self.vlab * math.cos(self.theta)


