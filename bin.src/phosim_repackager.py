#!/usr/bin/env python
"""
Script to repackage phosim amplifier files into single sensor
multi-extension FITS files, with an HDU per amplifier.
"""
import argparse
import lsst.phosim.utils as phosim_utils

parser = argparse.ArgumentParser(
    description="Repackager for phosim files (default is amp files)")
parser.add_argument('visit_dir', type=str, help="visit directory")
parser.add_argument('--eimage', default=False, action='store_true',
                    help='repackage the eimage')
parser.add_argument('--out_dir', type=str, default=None,
                    help="output directory")
parser.add_argument('--inst', default='lsst', type=str,
                    help='Instrument to use: comcam or lsst. (default: lsst). \
                    It corresponds to the prefix of PhoSim amplifier files (usually lsst or comcam).\
                    Note: corner sensors are included in the lsst instrument.')
parser.add_argument('--verbose', default=False, action='store_true',
                    help='print time to process the data each sensor')
parser.add_argument('--image_type', type=str, default='skyexp',
                    help="image type, passed as 'IMGTYPE': skyexp, dark,\
                    flat, or bias (default: skyexp)")
args = parser.parse_args()

repackager = phosim_utils.PhoSimRepackager(instName=args.inst, image_type=args.image_type)

if (args.eimage):
    repackager.process_visit_eimage(args.visit_dir, out_dir=args.out_dir,
                                    instName=args.inst, verbose=args.verbose)
else:
    repackager.process_visit(args.visit_dir, out_dir=args.out_dir,
                             instName=args.inst, verbose=args.verbose)
