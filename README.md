sedInterFoam
=======

This repository provides the sedInterFoam solver.

Status
------

The sedInterFoam solver is in development.

Pull requests are encouraged!

Features
--------
A three-dimensional two-phase flow solver with resolution of a free surface, sedInterFoam, has been developed for sediment transport applications. The solver is extended from sedFoam, itself extended from twoPhaseEulerFoam available in the 2.1.0 release of the open-source CFD (computational fluid dynamics) toolbox OpenFOAM.

Installation
------------

```bash
cd $WM_PROJECT_USER_DIR
git clone --recurse-submodules https://github.com/mathieu7an/sedInterFoam sedInterFoam
cd sedInterFoam
./Allwclean
./Allwmake
```
