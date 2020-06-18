# SOLSTISE_Simulation
SOLSTISE Simulation Code

The SOLSTISE Simulation code simulates reactions within a magnetic spectrometer (HELIOS or SOLARIS) and propogates the particles through the magnetic field of the solenoid. SOLSTISE is a gas jet target currently in development at Oak Ridge National Laboratory. 

The code is written in Python 3, so it (obviously) is a prerequisite to run the code. Most of the packages necessary to run the code are included in Anaconda, so it is recommended to just install Anaconda on your machine from https://docs.anaconda.com/anaconda/install/ or using the method of your choice. 

If you don't want to install Anaconda for whatever reason (not to worry, these snakes don't jump), *probably* the only packages you'd need to install are: matplotlib, sympy, and numpy, and pandas. You may also need X11 installed and running for the matplotlib plots. 

Required Files:

SOLSTISE_Sim.py

CircleAreaCalc.py

Sim_File_pd.py

EventBuilder.py

Plotter.py

PipeMaker.py

ConeMaker.py

NozzleMaker.py

stopyt.py

massreader.py

masses.txt

isotopetable.txt

Code Summary:

The code can be run by typing:

>>python3 SOLSTISE_Sim.py

Additionally, individual codes (EventBuilder, Plotter.py, stopyt.py) can be run in "standalone" mode if necessary in the same fashion. stopyt.py has its own front-end that can calculate energy loss like the old stopit fortran code. 

The code currently (6/18/20) on first run when downloaded from Github will make 3 new directories, it will then copy a couple of hidden text files to the Geometry_Files folder if they aren't there already. Then, since no event files have been made, it will immediately run EventBuilder. EventBuilder has 3 different modes of operation. The user can specify whether or not they want to produce results with artificial energy loss/smearing (artsm.txt), calculated energy loss (using stopyt) (eloss.txt), or an "all energies" option (allE.txt), which can be used if the user does not want to specify specific populated energy levels. Using the "all energies" option allows for the creation of blocked particle contour plots in Plotter.py. The EventBuilder code will produce a .txt file with the particle information and a .pkl file that has the target information in the Event_Files folder (if using the stopyt energy loss option). Be aware that if you try to make a file with the same reaction and energy as one that has been previously made, it will be overwritten (if it has the same extension).

It then reads in the text file (.txt or .dat) that is output from EventBuilder.py (two columns, lab angle and energy when no energy loss is used. Three columns with energy loss, lab angle, energy, thickness/z-position). 

***

The code then will ask if you want to make files for new Cones, Nozzles, and Pipes. The default cone and nozzle files are contained in the Geometry_Files folder, and can be readily used. If you want to make a new cone or nozzle, you can follow the instructions in the code to do so. In addition, the area of the pipe is calculated to be equivalent to an ISO 160 if no custom pipe is specified. The code will ask which way the pipe is facing (downstream or upstream). If a custom pipe is created, a PipeOut file will be saved to the Geometry_Files folder.

Then, the code reads the text file and propogates the particles through the magnetic field, taking into account shadowing by the jet nozzle, receiver cone, and pipe (if the pipe was specified as being in the half of the solenoid without detectors, there is still some shadowing from the top portion of the pipe. 

The "detectors" are currently split into four quadrants when looking downstream at the target. In this configuration, the x-axis points to the right, the y-axis points up, and the z-axis points downstream. So, Detector 1 is defined as being 0 < Phi < 90, Detector 2 is 90 < Phi < 180, detector 3 is 180 < Phi < 270, and detector 4 is 270 < Phi < 0.

The code finishes by saving a .pkl file of the final data frame, which can be run in Plotter.py at any time. This file is not overwritten if you run the code again (if you don't want it to be). Plotter can then be used to draw various histograms (unblocked particles, shadowed particles, etc). Additionally, contour plots of the shadowed particles can be drawn if the "all energies" option was used in EventBuilder. Special histograms can be drawn if a solid target is used. 

Questions? Email: mhall12@alumni.nd.edu
