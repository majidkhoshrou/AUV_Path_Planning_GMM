README for the Simulink Simualtions "Cooperation_Strategy_I/II/III"
12.02.2014
For Questions please e-mail to: and-hh@gmx.de

-----------STARTING A SiMULATION------------
Before starting a simulation, it is necessary to run the simulation_init.m script.
After this, the simulation can be started by clicking "run" in Simulink.

-----------The simulation_init.m----------------------
The simulation_init.m contains all parameters of the simulation such as intial states, noise gains, control gains, pathdata, etc.
An important property is the "TrsmtLaw" vector. It contains the order of the broadcasting vehicle/beacon.
Example: TrsmtLaw = [4  1  4  2  4  3] means that the beacon (4) is broadcasting between each vehicle (1,2,3).
With TrsmtLaw = [1 2 3]; the beacon is not active.

----------Displaying Results--------------
The Positions of the vehicles are stored in Position1, Position2, Position3. [x-pos, y-pos, psi]
The Results from the estimator are stored in Sensors1, Sensors2, Sensors3. [Vcx, Vcy, x-est, y-est, psi-est]
The script Plotting_Actual_vs_Estimated_Position.m extracts two simple plots and the RMS