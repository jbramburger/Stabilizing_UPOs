This repository contains the MATLAB files to reproduce the data and figures from [**Data-driven stabilization of periodic orbits**](https://arxiv.org/abs/2010.13896) by Jason J. Bramburger, J. Nathan Kutz, and Steven L. Brunton (IEEE Access, 2021). Computations use the publicly available SINDy architecture found at https://faculty.washington.edu/kutz/page26/ and should be stored in a folder entitled 'Util'. 

The scripts associated to this repository are as follows:

- Henon_control.m: Stabilizing cyclic orbits of the Henon mapping. Produces the figures presented in Section 4.1. Imports the control matrices from Henon_control_matrices.mat.

- Isolated_discovery.m: Mapping discovery for the isolated periodic orbit in Section 4.2.

- Isolated_control.m: Stabilization proceedure of the isolated periodic orbit using the mapping found in Isolated_discovery.m. Produces the figures in Section 4.2.

- Rossler_discovery.m: Parameter-dependent mapping discovery for the Rossler system. Produces the mappings presented in Section 4.3.

- Rossler_control.m: Simulates the controlled orbits of the Rossler system. Produces the figures in Section 4.3.

- Sprott_discovery.m: Parameter-dependent mapping discovery for the Sprott system. Produces the mappings used for the results in Section 4.4.

- Sprott_control.m: Simulates the controlled orbits of the Sprott system. Produces the figures in Section 4.4.

- Satellite_discovery.m: Discovers a mapping near the Lagrange points in an Earth-moon restricted three-body problem. Produces the linearizations presented in Section 5.

- Satellite_control.m: Simulates the controlled orbits of satellites made to sit near the Lagrange points in an Earth-moon restricted three-body problem. Produces the figures in Section 5.

A video abstract associated to this code and the corresponding paper is available at: https://www.youtube.com/watch?v=72MeRaNki8E
