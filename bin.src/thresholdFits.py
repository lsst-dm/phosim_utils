#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012, 2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import os
import argparse
from glob import glob
import astropy.io.fits as fits


def expand_files(filename_list=None):
    """Expand a set of filenames and globs, returning a list of filenames

    Parameters
    ------------
    filename_list:  list, optional
        A list of files and glob patterns

    Returns:
    -----------
    expandedList : list
        An expanded list of files that exist and files that fit the glob pattern
    """
    expanded_list = []
    for globPattern in filename_list:
        files = glob(globPattern)

        if not files:
            print("%s doesn't match any file" % globPattern)
            continue

        expanded_list.extend(files)

    return expanded_list


def get_size(filepath=None):
    """
    Convenience function to return the filesize.

    Parameters:
    -----------
    filepath : str, optional
        Path to the file to check

    Returns:
    --------
    size : str
        Human-readable filesize (Kilobytes or Megabytes)
    """
    size_bytes = os.path.getsize(filepath)
    if size_bytes < (1024 * 1024):
        size = str(round(size_bytes / 1024)) + "K"
    else:
        size = str(round((size_bytes / (1024 * 1024)), 1)) + "M"
    return size


def threshold_file(
    fits_file=None,
    threshold=None,
    not_binary=False,
    subtract_threshold=False,
    out_dir=None,
):
    """
    Threshold a FITS file, by changing the value of pixels for
    each .data element. By default, pixels below the threshold
    are set to 0, and those above to 1. The headers are kept
    unchanged.

    Parameters:
    -----------
    fits_file : str, optional
        Path to the FITS file to threshold
    threshold : float, optional
        A value of threshold to apply
    not_binary: bool, optional
        If set to True, do non-binary thresholding.
        In that mode by default we keep the original value of
        pixels above the threshold
    subtract_threshold: bool, optional
        If set to True, set pixels below threshold to 0,
        from those above, subtract the threshold
    out_dir : str, optional
        Path to output directory
    """
    print("\nThresholding %s" % fits_file)

    hdul = fits.open(fits_file)
    for i in range(len(hdul)):
        raw_data = hdul[i].data[:]
        thresh_data = raw_data.copy()

        if threshold is None:
            gain = hdul[0].header["GAIN"]
            seed = hdul[0].header["SEED"]
            threshold = seed + 3 * gain
            print("SEED=%.2f, GAIN=%f" % (seed, gain))
            print("threshold = %.2f + 3*%f = %f" % (seed, gain, threshold))

        # select data below threshold
        below_thresh = raw_data < threshold

        # the data below threshold is always set to 0
        thresh_data[below_thresh] = 0

        # the data above the threshold:
        if not_binary:

            # subtract the threshold value if needed
            if subtract_threshold:
                print("subtracting threshold")
                thresh_data[~below_thresh] = raw_data[~below_thresh] - threshold

        else:
            # binary mode: set to 1
            print("setting to 1 ")
            thresh_data[~below_thresh] = 1

        # replace the raw data with thresholded data
        hdul[i].data[:] = thresh_data

    # save the fits file
    out_file = os.path.join(out_dir, fits_file)
    hdul.writeto(out_file, overwrite=True)

    # close the fits file
    hdul.close()

    # print the file size reduction
    old_size = get_size(fits_file)
    new_size = get_size(out_file)
    print("The filesize has been reduced from %s to %s" % (old_size, new_size))


def main(
    filename=None,
    threshold=None,
    not_binary=False,
    subtract_threshold=False,
    out_dir=None,
):
    # choose binary mode
    if not_binary is False:
        print(
            "\nMaking a binary image: \
        px < thresh are set to 0,  \
        px > thresh are set to 1."
        )

    # choose non-binary mode: select what to do with pixels above the
    # threshold
    elif not_binary is True:
        print("\nMaking a not-binary image:")
        if subtract_threshold is True:
            print("subtract-threshold mode, i.e. ")
            print(
                "px < thresh are set to  0, \
            \npx > thresh are set to px - thresh "
            )

        else:
            print("By default, keeping original value above the threshold, i.e.")
            print(
                " all px < thresh are set to 0, all px > thresh are kept at original value"
            )

    # in both binary and non-binary mode
    if threshold is None:
        print(
            '\nThreshold is not provided; using by default \
        threshold =  "SEED" + 3 * "GAIN" from the header .'
        )

    fits_files = expand_files(filename)
    print("Files to threshold:", fits_files)

    # make the output directory
    if out_dir is None:
        print("Saving thresholded files in ")
        out_dir = os.path.join(os.getcwd(), "thresholded")
        print(out_dir)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # threshold each fits file
    for fits_file in fits_files:
        threshold_file(
            fits_file,
            threshold=threshold,
            not_binary=not_binary,
            subtract_threshold=subtract_threshold,
            out_dir=out_dir,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Threshold FITS image files to reduce file size"
    )
    parser.add_argument(
        "filename",
        nargs="+",
        help="file(s) to threshold. Accept explicit filenames, \
    eg. a.fits b.fits c.fits, or a glob pattern, eg *.fits",
    )

    parser.add_argument("--threshold", "-t", nargs=1, type=float)
    parser.add_argument(
        "--not_binary",
        "-nb",
        action="store_true",
        help="By default we make a binary image (by default pixels \
                        below threshold set to 0, and above threshold set to 1). Set this flag to do\
                        otherwise. If --notBinary flag is used, then by default we keep the\
                        original value of pixels above the threshold. To subtract the threshold\
                        value, set also --subtractThreshold flag.",
        default=False,
    )

    parser.add_argument(
        "--subtract_threshold", "-sub", action="store_true", default=False
    )
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory. By default, make /thresholded/ \
                        directory inside the working directory.",
    )
    args = parser.parse_args()

    main(
        filename=args.filename,
        threshold=args.threshold,
        not_binary=args.not_binary,
        subtract_threshold=args.subtract_threshold,
        out_dir=args.out_dir,
    )
