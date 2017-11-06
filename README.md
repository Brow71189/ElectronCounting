ElectronCounting
============

This python package implements software electron counting for low-dose STEM images.

Installation
============

After downloading and extracting the package open a terminal in the package folder and run
```
python3 setup.py install --user
```
or
```
sudo python3 setup.py install
```
if you want to install it system-wide.

You should then be able to call
```
from ElectronCounting import c_electron_counting
```
within your python installation.

Requirements
============

* numpy
* scipy

Usage
=====

The module `c_electron_counting` contains two functions for image processing:

`electron_counting(image, baseline=0.0, countlevel=0.01, peakwidth=1, only_integrate=False):`
--------------------------------------------------------------------------------------------

This function is used to integrate the single electron peaks and (if enabled) perform electron counting in the raw image.

### Parameters:

* __image__: Input image as floating point numpy array.
* __baseline__ (optional): Zero offset that will be subtracted from all pixel intensities (usually very close to zero, so you can omit it most of the time). Defaults to 0.
* __countlevel__ (optional): The typical integral of the single electron impact peaks. A good way to find this parameter is to bin the raw image in x-direction by factor larger than the typical peakwidth and look for the single electron peak in the histogram of the binned image. The position of this peak is the countlevel. A binning factor of 10 is sufficient in most cases. Note that __countlevel__ is important for the program to work correctly so you should always pass it even though it is optional. Defaults to 0.1.
* __peakwidth__ (optional): The typical width of the single electron impact peaks. This parameter can be calculated from the slope of a horizontal line profile with logarithmic y-axis of the power spectrum of the raw image. This slope then has to be multiplied with the half of the image width and, most of the times, divided by 2 Pi. Whether the factor 2 Pi is necessary depends on the normalization in the respective Fourier transform implementation you are using. However, this is the most widely used definition. Note that __peakwidth__ is important for the program to work correctly so you should always pass it even though it is optional. Defaults to 1.
* __only_integrate__ (optional): If set to `True`, the actual electron counting part will be omitted and the result will only contain the integrated electron signals. Defaults to `False`.

### Returns:

* __result__: The processed image as 32-bit floating point numpy array. Note that the data type will always be float, even when electron counting is enabled. This is to have one common data type for both use cases of this function. In electron counting mode the result can simply be casted to integer.

`find_most_likely_counts(image, baseline=0.0, countlevel=0.1):`
--------------------------------------------------------------

This function is used for transforming a raw image into electron counts. If no integration of the single electron peaks is necessary (e.g. because the pixel dwell time was above the critical value) this is the right function to use.

### Parameters:

* __image__: Input image as floating point numpy array.
* __baseline__ (optional): Zero offset that will be subtracted from all pixel intensities (usually very close to zero, so you can omit it most of the time). Defaults to 0.
very close to zero, so you can omit it most of the time). Defaults to 0.
* __countlevel__ (optional): The typical integral of the single electron impact peaks. A good way to find this parameter is to bin the raw image in x-direction by factor larger than the typical peakwidth and look for the single electron peak in the histogram of the binned image. The position of this peak is the countlevel. A binning factor of 10 is sufficient in most cases. Note that __countlevel__ is important for the program to work correctly so you should always pass it even though it is optional. Defaults to 0.1.

### Returns:

* __result__: The processed image as 32-bit floating point numpy array. Note that the data type will always be float, even though the electron counts are integer numbers. This is to have one common data type for all functions in this module.

Command line script
===================

There is a command line script included in the package for fast and easy application of the above module.

It is a wrapper script for the two functions `electron_counting` and `find_most_likely_counts` from the `c_electron_counting` module. You can call `python3 electron_counting.py -h` to see the possible command line arguments. For a more detailed explanation of all parameters have a look at the above section of this `README` or at the documentation of the `c_electron_counting` module.

