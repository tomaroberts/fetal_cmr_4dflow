import argparse

## parse arguments
parser = argparse.ArgumentParser(description='Create 4D flow cine MRI visualisation of whole fetal heart.')
parser.add_argument('fcmrDir', metavar='fcmrDir', type=str, nargs='+',
                    help='path to fcmr directory')
parser.add_argument('cineVtkFilename', metavar='cineVtk', type=str, nargs='+',
                    help='name of cine*.vtk files - nb: omit numbers and .vtk')
parser.add_argument('vel_volVtkFilename', metavar='vel_volVtk', type=str, nargs='+',
                    help='name of vel_vol*.vtk files - nb: omit numbers and .vtk')
parser.add_argument('numFrames', metavar='numFrames', type=int,
                    help='sum the integers (default: find the max)')

args = parser.parse_args()
print(args.fcmrDir)
print(args.cineVtkFilename)
print(args.vel_volVtkFilename)
print(args.numFrames)