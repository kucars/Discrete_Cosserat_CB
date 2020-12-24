# Composite-Body Algorithm of Discrete Cosserat Method

This MATLAB code implements the composite-body algorithm of the Discrete Cosserat method for soft manipulator dynamics, as described in

F. Renda, F. Boyer, J. Dias and L. Seneviratne, "Discrete Cosserat Approach for Multisection Soft Manipulator Dynamics," in IEEE Transactions on Robotics, vol. 34, no. 6, pp. 1518-1533, Dec. 2018, doi: 10.1109/TRO.2018.2868815.

Insert your soft manipulator parameters in:
```
piecewise_driver
```
and run it to simulate the dynamics.
The differential equations are implemented in:
```
piecewise_derivatives
```
The state vector is composed by one 6D vector of the section strain for each section of the manipulator. The integration quantities such as position, orientation and velocity of each cross section of the manipulator are calculated in:
```
piecewise_observables
```

# Cite

If you use this code for scientific publication purpose, please cite:
```
@ARTICLE{8500341,  author={F. {Renda} and F. {Boyer} and J. {Dias} and L. {Seneviratne}},  
journal={IEEE Transactions on Robotics},   
title={Discrete Cosserat Approach for Multisection Soft Manipulator Dynamics},   
year={2018},  
volume={34},  
number={6},  
pages={1518-1533},  
doi={10.1109/TRO.2018.2868815}}
```
# Author

The author of the algorithm and MATLAB implementation is:

Federico Renda (federico.renda@ku.ac.ae)
Department of Mechanical Engineering, Khalifa University Robotics Institute,
Khalifa University of Science and Technology, Abu Dhabi.
