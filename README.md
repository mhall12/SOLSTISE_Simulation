# SOLSTISE_Simulation
SOLSTISE Simulation Code

The SOLSTISE Simulation code takes output from Vikar and propogates the particles through the magnetic field of a solenoid (HELIOS or SOLARIS). SOLSTISE is a gas jet target currently in development at Oak Ridge National Laboratory. 

The code is written in Python 3, so it (obviously) is a prerequisite to run the code. Most of the packages necessary to run the code are included in Anaconda, so it is recommended to just install Anaconda on your machine from https://docs.anaconda.com/anaconda/install/ or using the method of your choice.

If you don't want to install Anaconda for whatever reason (not to worry, these snakes don't jump), the required packages are: matplotlib, sympy, and numpy. You may also need X11 installed and running for the matplotlib plots. 

Required Files:

SOLSTISE_Sim.py

CircleAreaCalc.py

Sim_File.py

massreader.py

masses.txt

Code Summary:

The code can be run by typing:

>>python3 SOLSTISE_Sim.py

The code currently (10/10/19) reads in a text file (.txt or .dat) that is an output from Vikar (two columns, lab angle and energy). It then starts by asking the user questions about the geometry of the solenoid and the pipe at the bottom of the solenoid. The area of the pipe is calculated to be equivalent to an ISO 160. 

Then, the code reads the text file and propogates the particles through the magnetic field. The final plot that appears represents the particles detected in the four quadrants of a cylindrical detector along the beam axis. 

Questions? Email: mhall12@alumni.nd.edu
