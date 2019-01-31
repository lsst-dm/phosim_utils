#!/usr/bin/env python
"""
Script to repackage phosim amplifier files into single sensor
multi-extension FITS files, with an HDU per amplifier.
"""
import argparse
import lsst.phosim.utils as phosim_utils

parser = argparse.ArgumentParser(description="Repackager for phosim amp files")
parser.add_argument('visit_dir', type=str, help="visit directory")
parser.add_argument('--out_dir', type=str, default=None,
                    help="output directory")
parser.add_argument('--prefix', default='lsst', type=str,
                    help='prefix of PhoSim amplifier files (usually lsst or comcam)')
parser.add_argument('--verbose', default=False, action='store_true',
                    help='print time to process the data each sensor')
args = parser.parse_args()

repackager = phosim_utils.PhoSimRepackager()
repackager.process_visit(args.visit_dir, out_dir=args.out_dir,
                         prefix=args.prefix, verbose=args.verbose)
