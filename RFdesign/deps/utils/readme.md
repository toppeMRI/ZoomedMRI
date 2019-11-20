# General Purpose Utils Repo

https://gitlab.com/fMRI/utils/readme.md

The repo is dedicated to share and keep track of utilities that are non-specific
to projects, e.g. mri-fft-like functions, etc.\
So wheels are not reinvented overtime, and function calls are consistent.

It would be nice of you to have your own `dev` branch for add-in new utils, and
test them before merging into `master`.

## Configs

At first time, navigate your matlab to this repo, and run,\
`>> setup_fmri_utils()`\
to add all subfolders in this repo to matlab path.
If you also want to save the paths, run\
`>> setup_fmri_utils(true)`\
After setting up, at any where,\
`>> fmri_utilsPath;`\
will navigate your matlab to the root directory of this repo.\
These two functions are instances of functions in `./meta`.
