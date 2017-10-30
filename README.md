AnalyzeMaxima
============

This python package finds the local maxima in an image. If there is noise present in the image, make sure to blur it before
passing it to "local maxima", otherwise almost every pixel will be marked as a maximum. You can also use a high value
for "noise tolerance", but usually using a blurred image and a lower value for "noise tolerance" gives better results.

The algorithm is a python implementation of the ImageJ plugin "MaximumFinder"².

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

If the package is intended to be used with Nion Swift¹ plugins make sure that you run "setup.py" with EXACTLY the same python interpreter that is used for running Nion Swift. Otherwise it will most likely not work. For a standard installation on Linux as described on Nion's website run
```
~/miniconda3/bin/python3 setup.py install
```
in a terminal that you opened in the project folder.


¹ www.nion.com/swift

³https://github.com/imagej/imagej1/blob/master/ij/plugin/filter/MaximumFinder.java
