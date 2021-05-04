# Singularity-Avoidance of Task-Redundant Robots in Pointing Tasks: On Nullspace Projection and Cardan Angles as Orientation Coordinates

Authors: Moritz Schappler, and Tobias Ortmaier

This code is supplemental material to the following submission to the ICINCO Conference 2021:

```bib
@InProceedings{SchapplerOrt2021,
  author    = {Schappler, Moritz and and Ortmaier, Tobias},
  booktitle = {ICINCO},
  title     = {Singularity-Avoidance of Task-Redundant Robots in Pointing Tasks: On Nullspace Projection and Cardan Angles as Orientation Coordinates},
  year      = {2021},
}
```


## Contents and Useage

This repository contains Matlab scripts and Latex code to reproduce results and statements of the paper. Beforehand, the steps from the [prerequesites section](#prerequisites) have to be performed.

### Replicate Results

Each section of the paper can be reproduced by a Matlab script.

* Sec. 3 (serial-link robot modeling regarding nullspace motion in pointing task)
  * `proof/SerRob_nullspace_proof.m`
* Sec. 4 (parallel robot modeling regarding nullspace motion in pointing task)
  * `proof/ParRob_nullspace_proof.m`
* Sec. 5 (aspects of control loop design for nullspace motion generation)
  * `statement_stability_I2_P_controlled.m` (one statement regarding stability of P controller in the linear case)
* Sec. 6 (simulation results for the nullspace motion)
  * `case_study/ParRob_nullspace_static_pose.m` (this also creates the figures of the section).

## Prerequisites <a name="prerequisites">

For the Matlab scripts to work, the following steps have to be performed:

1. Install Matlab (tested with R2021a)
2. Set up the Matlab mex compiler, if not already done
3. Download Matlab dependencies:
    * [Toolbox for geometry functions](https://github.com/SchapplM/robotics-dep-geometry) (tested with rev. 82c3215; 2021-02-05)
    * [Toolbox for trajectory functions](https://github.com/SchapplM/robotics-dep-traj) (tested with rev. 41ceefd; 2020-02-15)
    * [Collection of external dependencies](https://github.com/SchapplM/robotics-dep-ext) (tested with rev. ce08b62; 2021-03-08)
    * [Toolbox for various Matlab functions](https://github.com/SchapplM/matlab_toolbox) (tested with rev. 31298d1; 2021-03-30)
    * [Robotics toolbox](https://github.com/SchapplM/robotics-toolbox) (tested with rev. 91b01b4; 2021-05-03)
    * [Serial robot model database](https://github.com/SchapplM/robsynth-serroblib) (tested with rev. f608d0c; 2021-04-09)
    * [Parallel robot model database](https://github.com/SchapplM/robsynth-parroblib) (tested with rev. 8cb58ed; 2021-04-08)
    * Update the repositories to the given versions. If the master branch is on an older version, change of branch is required.
    * Short command to download all repositories (run in Git Bash in Windows or in Linux terminal):
      ```Bash
      git clone git@github.com:SchapplM/robotics-dep-geometry.git -b dev
      git clone git@github.com:SchapplM/robotics-dep-traj.git -b master
      git clone git@github.com:SchapplM/robotics-dep-ext.git -b dev
      git clone git@github.com:SchapplM/matlab_toolbox.git -b dev
      git clone git@github.com:SchapplM/robotics-toolbox.git -b taskred2
      git clone git@github.com:SchapplM/robsynth-serroblib.git -b dev
      git clone git@github.com:SchapplM/robsynth-parroblib.git -b dev
      ```
4. Run the path initialization scripts (`..._path_init.m`) of all downloaded software repos in Matlab. Put the run commands in the startup.m file of Matlab if used frequently. The command for 

