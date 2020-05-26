# SOLSTISE_Simulation
SOLSTISE Simulation Code

The SOLSTISE Simulation code simulates reactions within a magnetic spectrometer (HELIOS or SOLARIS) and propogates the particles through the magnetic field of the solenoid. SOLSTISE is a gas jet target currently in development at Oak Ridge National Laboratory. 

The code is written in Python 3, so it (obviously) is a prerequisite to run the code. Most of the packages necessary to run the code are included in Anaconda, so it is recommended to just install Anaconda on your machine from https://docs.anaconda.com/anaconda/install/ or using the method of your choice. 

If you don't want to install Anaconda for whatever reason (not to worry, these snakes don't jump), the required packages are: matplotlib, sympy, and numpy, and pandas. You may also need X11 installed and running for the matplotlib plots. 

Required Files:

SOLSTISE_Sim.py

CircleAreaCalc.py

Sim_File_pd.py

EventBuilder.py

Plotter.py

stopyt.py

massreader.py

masses.txt

isotopetable.txt

Code Summary:

The code can be run by typing:

>>python3 SOLSTISE_Sim.py

Additionally, individual codes (EventBuilder, Plotter.py, stopyt.py) can be run in "standalone" mode if necessary in the same fashion. stopyt.py has its own front-end that can calculate energy loss like the old stopit fortran code. 

The code currently (5/26/20) reads in a text file (.txt or .dat) that is an output from EventBuilder.py (two columns, lab angle and energy). Within the EventBuilder code, the user can specify whether or not they want to produce results with artificial energy loss/smearing, calculated energy loss (using stopyt), or an "all energies" option, which can be used if the user does not want to specify specific populated energy levels. Using the "all energies" option allows for the creation of blocked particle contour plots in Plotter.py. The EventBuilder code will produce a .txt file with the particle information and a .pkl file that has the target information (if using the stopyt energy loss option).

***
The code starts by asking the user questions about the input files and geometry of the solenoid and the pipe at the bottom of the solenoid. The area of the pipe is calculated to be equivalent to an ISO 160. Additionally, the code now will ask which way the pipe is facing (downstream or upstream). If a custom pipe is created, a PipeOut file will be saved. 

Then, the code reads the text file and propogates the particles through the magnetic field, taking into account shadowing by the jet nozzle, receiver cone, and pipe (if it is defined in the same half of the solenoid as the detectors). The "detectors" are currently split into four quadrants.

The code finishes by saving a *.pkl file of the final data frame, which can be run in Plotter.py at any time. Plotter can then be used to draw various histograms (ublocked particles, shadowed particles, etc). Additionally, contour plots of the shadowed particles can be drawn if the "all energies" option was used in EventBuilder. 

Questions? Email: mhall12@alumni.nd.edu
