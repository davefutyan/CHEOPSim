#! /usr/bin/env python
import argparse
import numpy.random

import spotlc

# Useful function to underline text
def u(s):
    return '\033[4m'+s+'\033[0m'

# Initialise argument parser
parser = argparse.ArgumentParser(description='Produce light curves of '
                                 'spotted stars to use with CHEOPSsim.')

# Add arguments
parser.add_argument('prot', type=float,
                    help='Rotational period of the simulated star')
parser.add_argument('sptype', type=str,
                    help='Spectral type of the simulated star')

parser.add_argument('-t', '--ttime', type=float, metavar=u('float'),
                    default=2880,
                    help='Total length of simulation in minutes')

parser.add_argument('-c', '--cadence', type=float, metavar=u('float'),
                    default=120,
                    help='Resolution in minutes (defaults to 120 minutes)')

parser.add_argument('-n', '--ncurves', type=int, metavar=u('int'),
                    default=1,
                    help='Number of light curves to produce simultaneously')

parser.add_argument('--filetemplate', type=str, metavar=u('str'),
                    default='lc',
                    help='Template used to construct light curve files')

parser.add_argument('--nowrite', action='store_false',
                    help='Flag to block writing of output file')

parser.add_argument('--maxn', '-m', type=int, metavar=u('int'),
                    default=5000,
                     help='Maximum number of points in simulation. Above this'
                     ' number, a linear interpolation is used.')

parser.add_argument('--force', action='store_true',
                    help='Flag to force sampling to input value. Overrides '
                    'maxpoints')

parser.add_argument('--seed', type=int, metavar=u('int'), default=None,
                    help='Seed for random number generator. Passed to '
                    'numpy.random.seed')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')
#

# Parse arguments
args = parser.parse_args()

print('Simulating light curve with parameters')
print('Rotational period: {:.2f} days'.format(args.prot))
print('Spectral type: {}'.format(args.sptype))
print('########')
numpy.random.seed(args.seed)
time, flux = spotlc.spot_lc(args.prot, args.sptype, ttotal_minutes=args.ttime,
                            dt=args.cadence, ncurves=args.ncurves,
                            maxpoints=args.maxn, forcedt=args.force,
                            outtemplate=args.filetemplate, save=args.nowrite)
                            

