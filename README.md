# SOLSTISE_Simulation
SOLSTISE Simulation Code

The SOLSTISE Simulation code simulates reactions within a magnetic spectrometer (HELIOS or SOLARIS) and propogates the particles through the magnetic field of the solenoid. SOLSTISE is a gas jet target currently in development at Oak Ridge National Laboratory. 

The code is written in Python 3, so it (obviously) is a prerequisite to run the code. Most of the packages necessary to run the code are included in Anaconda, so it is recommended to just install Anaconda on your machine from https://docs.anaconda.com/anaconda/install/ or using the method of your choice. 

**Note** You need AT LEAST pandas version 1.0.3 and numpy version 1.15.1 to run the code, otherwise it will probably crash. You can check the versions via:

>>python3

>>import pandas as pd

>>import numpy as np

>>pd.\_\_version\_\_

>>np.\_\_version\_\_

If your version numbers are below the ones I used to write the code, you can update anaconda by typing:
>>conda update

If you don't want to install Anaconda for whatever reason (not to worry, these snakes don't jump), *probably* the only packages you'd need to install are: matplotlib, sympy, and numpy, and pandas. Everything else should come in the python standard library. You may also need X11 installed and running for the matplotlib plots. 

Code Summary:

For a longer summary and more explanation, check out the Jupyter notebook. It has much of the functionality of the full code and gives some useful explanation for first time users.

The code can be run in a terminal via:

>>python3 SOLSTISE_Sim.py

Additionally, individual codes (EventBuilder, Plotter.py, stopyt.py, etc) can be run in "standalone" mode if necessary in the same fashion. stopyt.py has its own front-end that can calculate energy loss like the old stopit fortran code. 

***

On first run, the code will make 3 new directories, it will then copy a couple of hidden text files to the Geometry_Files folder if they aren't there already. Then, since no event files have been made, it will immediately run EventBuilder. 

**Note** If the user is new to the SOLSTISE code, they can follow the example answers (... Ex) ##: ) at the end of each question to generate their first SOLSTISE input/output files. Then, they can plot the example data in Plotter to make sure the code is working as intended.  

EventBuilder has 3 different modes of operation. The user can specify whether or not they want to produce results with artificial energy loss/smearing (artsm.txt), calculated energy loss (using stopyt) (eloss.txt), or an "all energies" option (allE.txt), which can be used if the user does not want to specify specific populated energy levels. Using the "all energies" option allows for the creation of blocked particle contour plots in Plotter.py. The EventBuilder code will produce a .txt file with the particle information and a .pkl file that has the target information in the Event_Files folder (if using the stopyt energy loss option). Be aware that if you try to make a file with the same reaction and energy as one that has been previously made, it will be overwritten (if it has the same extension).

It then reads in the text file that is output from EventBuilder.py (two columns, lab angle and energy when no energy loss is used. Three columns with energy loss, lab angle, energy, thickness/z-position). 

***

The code then will ask if you want to make files for new Cones, Nozzles, and Pipes. The default cone and nozzle files are contained in the Geometry_Files folder, and can be readily used. If you want to make a new cone or nozzle, you can follow the instructions in the code to do so. In addition, the area of the pipe is calculated to be equivalent to an ISO 160 if no custom pipe is specified. The code will ask which way the pipe is facing (downstream or upstream). If a custom pipe is created, a PipeOut file will be saved to the Geometry_Files folder.

Then, the code reads the text file and propogates the particles through the magnetic field, taking into account shadowing by the jet nozzle, receiver cone, and pipe (if the pipe was specified as being in the half of the solenoid without detectors, there is still some shadowing from the top portion of the pipe. 

The "detectors" are currently split into four quadrants when looking downstream at the target. In this configuration, the x-axis points to the right, the y-axis points up, and the z-axis points downstream. So, Detector 1 is defined as being 0 < Phi < 90, Detector 2 is 90 < Phi < 180, detector 3 is 180 < Phi < 270, and detector 4 is 270 < Phi < 0.

The code finishes by saving a .pkl file of the final data frame, which can be run in Plotter.py at any time. This file is not overwritten if you run the code again (if you don't want it to be). Plotter can then be used to draw various histograms (unblocked particles, shadowed particles, etc). Additionally, contour plots of the shadowed particles can be drawn if the "all energies" option was used in EventBuilder. Special histograms can be drawn if a solid target is used. 

Questions? Email: mhall12@alumni.nd.edu

## Example output plots:

### E vs z contour plot in all four detector quadrants with gas jet target simulation overlay:
![image](https://user-images.githubusercontent.com/13751793/101208487-e0f22980-3637-11eb-9b66-1aeb7a275d35.png)
***
### E vs angle contour with overlay showing individual detector coverage:
![image](https://user-images.githubusercontent.com/13751793/101209935-46dfb080-363a-11eb-90a8-d45233f320d8.png)
***
### Polar plots showing fraction of particles detected vs initial phi angle:
![image](https://user-images.githubusercontent.com/13751793/101208830-5c53db00-3638-11eb-9782-bc7772248802.png)
![image](https://user-images.githubusercontent.com/13751793/101209060-bd7bae80-3638-11eb-9b44-0184a06b62e3.png)
***
### Fraction of particles detected vs initial lab angle:
![image](https://user-images.githubusercontent.com/13751793/101209193-fe73c300-3638-11eb-8d19-03f6229cf08d.png)


