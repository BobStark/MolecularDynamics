-This repository contains code for a molecular dynamics simulation.

-The system in question is a bar composed of a system of particles. The
 bar is subjected to a given strain rate of 10^9/s. The goal of the code
 is to simulate how the particles interact and how the overall system fails.


 The method for calculating the potential energy is the Lennard-Jones
 potential energy method. This method is relatively simple and accurate.
 
 To calculate forces, the definition of potential energy is used. Specifically,
 in one dimension, F = -dV/dr. Newton's laws are then integrated to find the 
 velocities and displacements.