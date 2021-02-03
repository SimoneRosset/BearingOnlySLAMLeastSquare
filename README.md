# BearingOnlySLAMLeastSquare

**AUTHOR:** 

*Simone Rossetti, Rome, January 2021*

[Istitutional](mailto:rossetti.1900592@studenti.uniroma1.it)

[Private](mailto:simone.rossetti@live.com)

**NB: the use of this material is free as long as the content is not altered and the name of the author and the link to this repo are always cited.**

<p float="center">
  <img src="/IMAGES/uni_rrobot_u10.png" width="70%%" title=" "/>
</p>

Project Bearing only SLAM in Octave

The dataset is composed by a g2o file which contain poses and bearing observations. The file contain also odometry edges that are use to construct the initial guess for the problem.

Steps:
    - Parse the whole dataset and initialize the landmarks (Linear Triangulation) by using at least two bearing observations with the proper parallax
    - Setup a LS optimization that involves all the poses and landmarks that you have initialized
    - In the bearing edges the ID of the landmarkd is reported, use it to identify them for both the triangulation and the global optimization
         
Expected output :
  - Robot trajectory and Map

  **Manipulator:**
<p  float="center">
  <img src="/IMAGES/6R.png" width="50%%" title=" "/>
</p>
