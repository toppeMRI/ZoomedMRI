# Center finders for Nd-array

Two short functions for finding the indices/subscripts of the center element of a nd-array.
Implemented to provide indices/subscripts consistency across functions.

Usage:
For an nd-array, `X`,

```matlab
Nd = size(X);
cInd = ctrInd(Nd);
cSub = ctrSub(Nd);
```
