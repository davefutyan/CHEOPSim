import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.style
matplotlib.style.use('classic')
import pylab as pl
import astropy.io.fits as pf
import sys,os,getopt
import numpy as np
import matplotlib.animation as animation
from matplotlib.colors import LogNorm

def main(argv):

    nimages = 10000
    
    try:
        opts, args = getopt.getopt(argv,"hi:n:")
    except getopt.GetoptError:
        print('animation.py -i <inputfile> -n <nimages>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('animation.py -i <inputfile> -n <nimages>')
            sys.exit()
        elif opt == "-i":
            inputdir = arg
        elif opt == "-n":
            nimages = arg

    if nimages == 10000:
        print('Generating animation from stacked images in',inputdir)
    else:
        print('Generating animation from the first',nimages,'stacked images in',inputdir)

    fig = pl.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    #ax.axis('off')

    first = True
    
    for file in os.listdir(inputdir):
        if file.find("SCI_RAW_SubArray")>-1:
            print(file)

            if first:
                data = pf.getdata(inputdir+"/"+file)
                im = ax.imshow(data[0],origin='lower',interpolation='nearest',norm=LogNorm())
                cbar = fig.colorbar(im,shrink=0.78)
                first = False
            else:
                data1 = pf.getdata(inputdir+"/"+file)
                data = np.concatenate((data,data1))

    def update_image(n) :
        im.set_data(data[n])
        return im

    nimages = min(int(nimages),len(data))
    ani = animation.FuncAnimation(fig,update_image,nimages,interval=200)

    ani.save('animation.mp4',bitrate=4000)

if __name__ == "__main__":
   main(sys.argv[1:])
