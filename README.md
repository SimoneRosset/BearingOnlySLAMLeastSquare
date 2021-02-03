# BearingOnlySLAMLeastSquare

**AUTHOR:** 

*Simone Rossetti, Rome, January 2021*

[Istitutional](mailto:rossetti.1900592@studenti.uniroma1.it)

[Private](mailto:simone.rossetti@live.com)

**NB: the use of this material is free as long as the content is not altered and the name of the author and the link to this repo are always cited.**


This script is implemented in Octave, install it and run octave-cli.


The dataset is composed by a g2o file which contain poses and bearing observations. The file contain also odometry edges that are use to construct the initial guess for the problem.

Steps:
  - Parse the whole dataset and initialize the landmarks (Linear Triangulation) using two bearing observations coming from the most orthogonal poses wrt the relative landmark.
  - Setup a LS optimization that involves all the poses and landmarks initialized
  - In the bearing edges the ID of the relative landmark is reported, it is used to identify it for both the triangulation and the global optimization
         
Expected output :
  - Robot trajectory and Map


**With noise on measurements:**

Optimization minimizes the global error but still the optimization is subject to local minima due to intrinsic structure of the problem. Open this repo in a terminal and run:

"""
octave-cli # to open octave

BearingOnlySLAMLeastSquare # in octave-cli

"""

<p float="center">
  <img src="/images/noise_chi.png" width="45%%" title=" "/ style="margin-right:10%"> <img src="/images/map_chi.png" width="45%%" title=" "/> 
</p>


**Zero noise on measurements:**

Optimization starting from the ground truth with zero noise on measurements produces always zero error, this proves the correctness of the implementation. Open this repo in a terminal and run:

"""
octave-cli 

test_no_noise

"""

<p float="center">
  <img src="/images/no_noise_chi.png" width="45%%" title=" "/ style="margin-right:10%"> <img src="/images/no_map_chi.png" width="45%%" title=" "/> 
</p>

