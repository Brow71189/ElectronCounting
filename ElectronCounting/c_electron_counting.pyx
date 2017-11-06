#!python
#cython: boundscheck=False
#cython: cdivision=True
# -*- coding: utf-8 -*-
"""
This python module implements software electron counting for low-dose STEM images. See the documentation of the
individual functions below for more detailed information.
"""

import numpy
from libc.math cimport round, atan, log, exp, pow, sqrt
from libc.stdlib cimport rand, RAND_MAX

DEF factorials_table_length = 20
DEF integration_range_gamma = 2.0
DEF pi = 3.14159

cdef long factorials_table[factorials_table_length]

def electron_counting(image, baseline=0.0, countlevel=0.1, peakwidth=1, only_integrate=False):
    """
    electron_counting(image, baseline=0.0, countlevel=0.1, peakwidth=1, only_integrate=False):

    This function is used to integrate the single electron peaks and, if enabled, perform electron counting in the
    raw image.

    Parameters
    ----------

    image : ndarray
        Input image as floating point numpy array.
    baseline : float, optional
        Zero offset that will be subtracted from all pixel intensities (usually very close to zero, so you can omit it
        most of the time). Defaults to 0.
    countlevel : float, optional
        The typical integral of the single electron impact peaks. A good way to find this parameter is to bin the raw
        image in x-direction by factor larger than the typical peakwidth and look for the single electron peak in the
        histogram of the binned image. The position of this peak is the countlevel. A binning factor of 10 is
        sufficient in most cases. Note that `countlevel` is important for the program to work correctly so you
        should always pass it even though it is optional. Defaults to 0.1.
    peakwidth : float, optional
        The typical width of the single electron impact peaks. This parameter can be calculated from the slope of a
        horizontal line profile with logarithmic y-axis of the power spectrum of the raw image. This slope then has
        to be multiplied with the half of the image width and, most of the times, divided by 2 Pi. Whether the factor
        2 Pi is necessary depends on the normalization in the respective Fourier transform implementation you are
        using. However, this is the most widely used definition. Note that `peakwidth` is important for the program
        to work correctly so you should always pass it even though it is optional. Defaults to 1.
    only_integrate : bool, optional
        If set to `True`, the actual electron counting part will be omitted and the result will only contain the
        integrated electron signals. Defaults to `False`.

    Returns
    -------

    result : ndarray
        The processed image as 32-bit floating point numpy array. Note that the data type will always be float, even
        when electron counting is enabled. This is to have one common data type for both use cases of this function.
        In electron counting mode the result can simply be casted to integer.

    """

    cdef float[:, :] c_image = image
    result = numpy.zeros(image.shape, dtype=numpy.float32)
    cdef float[:, :] c_result = result
    cdef float c_baseline = baseline
    cdef float c_countlevel = countlevel
    cdef float c_peakwidth = peakwidth
    cdef int c_height = image.shape[0]
    cdef int c_width = image.shape[1]
    cdef bint c_only_integrate = only_integrate
    c_electron_counting(c_image, c_result, c_height, c_width, c_baseline, c_countlevel, c_peakwidth, c_only_integrate)
    return result

def find_most_likely_counts(image, baseline=0.0, countlevel=0.1):

    """
    find_most_likely_counts(image, baseline=0.0, countlevel=0.1):

    This function is used for transforming a raw image into electron counts. If no integration of the single electron
    peaks is necessary (e.g. because the pixel dwell time was above the critical value) this is the right function to
    use.

    Parameters
    ----------

    image : Input image as floating point numpy array.
    baseline : float, optional
        Zero offset that will be subtracted from all pixel intensities (usually very close to zero, so you can omit
        it most of the time). Defaults to 0.
    countlevel : float, optional
        The typical integral of the single electron impact peaks. A good way to find this parameter is to bin the raw
        image in x-direction by factor larger than the typical peakwidth and look for the single electron peak in the
        histogram of the binned image. The position of this peak is the countlevel. A binning factor of 10 is
        sufficient in most cases. Note that `countlevel` is important for the program to work correctly so you
        should always pass it even though it is optional. Defaults to 0.1.

    Returns
    -------

    result : ndarray
        The processed image as 32-bit floating point numpy array. Note that the data type will always be float, even
        though the electron counts are integer numbers. This is to have one common data type for all functions in this
        module.

    """

    cdef float[:, :] c_image = image
    result = numpy.zeros(image.shape, dtype=numpy.float32)
    cdef float[:, :] c_result = result
    cdef float c_baseline = baseline
    cdef float c_countlevel = countlevel
    cdef int c_height = image.shape[0]
    cdef int c_width = image.shape[1]
    c_find_most_likely_counts_image(image, c_result, c_height, c_width, c_baseline, c_countlevel)
    return result

cdef void c_electron_counting(float[:, :] image, float[:, :] result, int height, int width, float baseline,
                              float countlevel, float peakwidth, bint only_integrate):
    global factorials_table
    cdef int k, i
    cdef float integral
    cdef int index
    cdef int randnum
    cdef float threshold
    cdef float counts_divisor
    cdef float peak_height
    cdef float peak_max
    cdef int peaklength

    create_factorials_table(factorials_table_length)
    # All following calculations are done with the height of the Lorentzians normalized to 1
    threshold = 1/(1.0+pow(integration_range_gamma, 2)) # height of our peak at +- integration_range_gamma
    # We assume that the peaks have been fully caputered when measuring their area. We therefore use a range of
    # 10 gamma for finding the peak height from the integral.
    peak_height = countlevel/peakwidth/(atan(10) - atan(-10))
    # Integration range going from -integration_range_gamma to +integration_range_gamma
    peaklength = <int>(round(2.0*integration_range_gamma*peakwidth))
    # correct integral for baseline
    counts_divisor = countlevel - baseline*(<float>peaklength)

    for k in range(height):
        integral = 0
        index = -1
        i = 0
        peak_max = peak_height
        while i < width:
            if index != -1 and (image[k, i] - baseline <= peak_max*threshold or i - index >= peaklength or i == width-1):
                if only_integrate:
                    result[k, index] = integral
                else:
                    result[k, index] = c_find_most_likely_counts(integral/counts_divisor)
                integral = 0
                if image[k, i] - baseline > peak_max*threshold:
                    randnum = <int>(<float>rand() / <float>RAND_MAX * <float>peaklength)
                    i = index + 1 + randnum
                    index = i
                    integral = image[k, i]
                else:
                    index = -1
                peak_max = peak_height

            elif index != -1:
                integral += image[k, i] - baseline
                if image[k, i] > peak_max:
                    peak_max = image[k, i]

            if index == -1 and image[k, i] - baseline > peak_max*threshold:
                index = i
                integral = image[k, i]

            i += 1

cdef void c_find_most_likely_counts_image(float[:, :] image, float[:, :] result, int height, int width,
                                          float baseline, float countlevel):
    cdef int k, i
    create_factorials_table(factorials_table_length)
    for k in range(height):
        for i in range(width):
            result[k, i] = c_find_most_likely_counts((image[k ,i]-baseline)/countlevel)

cdef float c_count_probability(int counts, float x):
    return ((1/sqrt(2.0*pi*(<float>counts)**2)*exp(-1.0*pow(x - <float>counts, 2.0)/(2.0*(<float>counts)**2))) *
            (pow((<float>counts), <float>counts)*exp(-1.0*(<float>counts))/(<float>factorials_table[counts])))

cdef float c_find_most_likely_counts(float x):
    if round(x+0.25) == 0.0:
        return 0.0
    cdef int counter
    cdef float max_probability, current_probability
    max_probability = 0
    counter = 0
    while counter < factorials_table_length:
        current_probability = c_count_probability(counter, x)
        if current_probability < max_probability:
            return counter - 1
        max_probability = current_probability
        counter += 1
    return factorials_table_length

cdef void create_factorials_table(int length):
    cdef int counter
    factorials_table[0] = 1
    counter = 1
    while counter < length:
        factorials_table[counter] = factorials_table[counter-1] * <long>counter
        counter += 1