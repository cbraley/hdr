hdr
===

C++ code for finding camera response functions(CTFs) and creating HDR images.  Implements the 1997 paper "Recovering High Dynamic Range Radiance Maps from Photographs".


Uses the Eigen library for linear algebra.


Current TODO list:
    -Make code fail gracefully when run with images that don't exist.
    -Include second executable for creating HDR images(this code here just creates the response curve)
    -Include a few example images w/ response curves and info on how to use this
    -Stratified sampling of the image; not pure random
    -Account for pixel bloom
