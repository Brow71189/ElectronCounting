#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program implements software electron counting for low-dose STEM images.

It is a wrapper script for the two functions `electron_counting` and `find_most_likely_counts` from the
`c_electron_counting` module.
You can call `python3 electron_counting.py -h` to see the possible command line arguments.
For a more detailed explanation of all parameters have a look at the `README` in this repository or at the
documentation of the `c_electron_counting` module.
"""


import argparse
from scipy.misc import imread, imsave
from ElectronCounting import c_electron_counting
import os


def run(image_path, mode, baseline, countlevel, peakwidth):
    image = imread(image_path, mode='F')

    if baseline is None:
        baseline = 0.0
    if peakwidth is None and mode in ['integrate', 'both']:
        raise ValueError('"peakwidth" must be set for modes "integrate" and "both".')

    if mode == 'count':
        result = c_electron_counting.find_most_likely_counts(image, baseline, countlevel)
    elif mode == 'integrate':
        result = c_electron_counting.electron_counting(image, baseline, countlevel, peakwidth, only_integrate=True)
    else:
        result = c_electron_counting.electron_counting(image, baseline, countlevel, peakwidth, only_integrate=False)

    dirname = os.path.dirname(image_path)
    filename = os.path.basename(image_path)
    savepath = os.path.join(dirname, 'ecount_' + filename)

    imsave(savepath, result)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('image_path', type=str, help='Path to the image that is to be processed')
    parser.add_argument('-c', '--countlevel', type=float, required=True, help='Typical area under the Lorentzians. ' +
                                                                              '(required)')
    parser.add_argument('-b', '--baseline', type=float, help='Baseline is subtracted from all pixel intensities ' +
                                                             'before integrating or counting electrons.')
    parser.add_argument('-p', '--peakwidth', type=float, help='Typical width of the Lorentzians in pixels (HWHM). ' +
                                                              'Only needed for modes "integrate" and "both".')
    parser.add_argument('-m', '--mode', type=str, choices=['integrate', 'count', 'both'],
                        help='What processing steps should we execute on the image. Can be one of "integrate", ' +
                             '"count" or "both". Defaults to "both".')
    parser.parse_args()
    run(parser.image_path, parser.mode, parser.baseline, parser.countlevel, parser.peakwidth)

if __name__ == '__main__':
    main()