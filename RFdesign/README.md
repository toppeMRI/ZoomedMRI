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
