Lesion Segmentation Applictations for Slicer4
=============================================

Overview
--------

A collection of CLI modules for Slicer4 to enable users to train lesion segmentation models on their own data and then perform segmentation. A module to train lesion segmentation models and another to use those models to segment lesions are included. 

Dependencies
------------
Currently [Slicer4](http://www.slicer.org/slicerWiki/index.php/Documentation/4.0/Developers/Build_Instructions) is required.  

Checkout, Configure, and Build
------------------------------
    $ git clone git@github.com:msscully/LesionSegmentation.git
    $ mkdir LesionSegmentation-build
    $ cd LesionSegmentation
    $ cmake -DSlicer_DIR:PATH=/path/to/Slicer-Superbuild-Debug/Slicer-build ../LesionSegmentation 
    $ make

Slicer Integration
------------------
Once this project has been built add the build directory to Slicer4's modules directory.  This can be done by going to "View"->"Application Settings"->"Module settings" and add the path to the LesionSemenation-build directory and restart Slicer.

Copyright
--------
(The 3D Slicer License)

Copyright © 2012 Mark Scully

This software is distributed under the 3D Slicer license, a BSD-style open source license that is compatible with the Open Source Definition by [The Open Source Initiative](http://opensource.org/) and contains no restrictions on use of the software. Please read the full [3D Slicer License Agreement](http://www.slicer.org/pages/LicenseText) before downloading or using this code.
