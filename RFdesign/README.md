# RF design

Pulse design program is wrapped inside `designScript.m`.

The program assumes [MIRT](https://github.com/JeffFessler/mirt) already installed as a dependency.

Other dependencies are included under `./dep`, run `./setup.m` to adding them into matlab path.

```md
makeIVpulse(iv, objSupportMask, fmap , [fovx, fovy, fovz]);
   div              [nx ny nz]    binary mask
   objSupportMask   [nx ny nz]    binary mask
   b0Map            [nx ny nz]    Hz
   fov              [1 3]         cm
   offset           [1 3]         cm
```

Cite [this](https://ieeexplore.ieee.org/document/7268765):

```bib
@article{Sun2016Joint,
author = {Sun, Hao and Fessler, Jeffrey A. and Noll, Douglas C. and Nielsen, Jon-Fredrik},
doi = {10.1109/TMI.2015.2478880},
issn = {0278-0062},
journal = {IEEE Transactions on Medical Imaging},
keywords = {Joint design,K-space trajectory design,MRI,RF pulse design,Tailored excitation},
month = {feb},
number = {2},
pages = {468--479},
pmid = {26390450},
title = {{Joint Design of Excitation k-Space Trajectory and RF Pulse for Small-Tip 3D Tailored Excitation in MRI}},
url = {http://ieeexplore.ieee.org/document/7268765/},
volume = {35},
year = {2016}
}
```
