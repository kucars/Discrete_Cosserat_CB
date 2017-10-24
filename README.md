# Composite-Body Algorithm of Discrete Cosserat Method for Soft Manipulators Dynamics

This MATLAB code implement the composite-body algorithm of the Discrete Cosserat method for soft manipulator dynamics, as described in

F. Renda, F. Boyer, J. Dias and L. Seneviratne. Discrete Cosserat Approach for Multi-Section Soft Robots Dynamics. arXiv:1702.03660 [cs.RO] (https://arxiv.org/abs/1702.03660).

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

