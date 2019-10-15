#!/usr/bin/env python
import numpy as np
from astropy.io import fits
import multiprocessing
from astropy.wcs import wcs
from astropy.coordinates import SkyCoord
from subprocess import call
import glob
import os
import sys

################# inputs ####################
## Bands to be extracted
band=['g','r','i','z']

## Filename of the object catalog (see next line)
fname="input_cat_example.txt"

## Galaxy_ID, Ra, Dec, Field_ID
## e.g. 3068640062 62.741055 -45.914261 DES0412-4540
gid,gra,gdec,fldname=np.loadtxt(fname,usecols=(0,1,2,3),unpack=1,dtype="string");

## Path to the input directory where original FITS are located
inpdirpath="/data/des80.a/data/kuropat/Y3A2"

## Half image size in pixels
half_im=20

## No. of processors on which to run e.g. "cat /proc/cpuinfo" to check the num of processors you have on a linux machine
Nproc=25

## Set flag=1 if you'd like to extract and  add exposure times in the output FITS files
flag=1

#############################################

if(flag):
    exparr=np.zeros(len(band))


def worker(num,nproc):
#   np.random.seed(num*10+29824);
    iimin=np.int(gid.size*num/nproc);
    iimax=np.int(gid.size*(num+1)/nproc);
    if(num==nproc-1):
        iimax=gid.size;
    for ii in range(iimin,iimax):

        if(not os.path.isdir("outfits/%s"%(fldname[ii]))):
            call("mkdir -p outfits/%s"%(fldname[ii]),shell=True)


        fitspath=glob.glob("%s/%s/p*/coadd/%s_*_g.fits.fz"%(inpdirpath,fldname[ii],fldname[ii]))
        foutg="outfits/%s/gal_%s_%s_%s.fits"%(fldname[ii],fldname[ii],gid[ii],band[0]) 
        print("fitspath", fitspath, foutg)

        ## If the g-band FITS does not exist, then create FITS files for all bands
        if(not os.path.isfile(foutg)):
            for kk in range(len(band)):
                print("Running band",kk)
                fitsname="%s_%s.fits.fz"%(fitspath[0][:-10],band[kk])
                hdulist = fits.open(fitsname)
                if(flag):
                    exparr[kk]=hdulist[1].header['EXPTIME']
                if(kk==0):
                    axis1=hdulist[1].header['NAXIS1']
                    axis2=hdulist[1].header['NAXIS1']
                    world = wcs.WCS(hdulist[1].header)
                    ## Extract the x,y pixel given the RA, Dec 
                    c_before = SkyCoord(gra[ii], gdec[ii], "icrs", unit="deg")
                    x_pix, y_pix = c_before.to_pixel(world, mode="all")
               
                    if(0.<x_pix<axis1 and 0.<y_pix<axis2):
                        print("### ",gid[ii], gra[ii], gdec[ii], fldname[ii], x_pix, y_pix, "%s_%s"%(fldname[ii],gid[ii]))
                        sys.stdout.flush();
                    else: ## Move to the next object if the current object is too close to the border of a tile 
                        print(gid[ii], gra[ii], gdec[ii], "offimage")
                        continue

                ## Slicing the FITS given the window size
                scidata=hdulist[1].data
                xl=np.int(x_pix)-half_im 
                xh=np.int(x_pix)+half_im
                yl=np.int(y_pix)-half_im 
                yh=np.int(y_pix)+half_im
             
                gal_im=scidata[yl:yh+1,xl:xh+1] 

                ## Setting output FITS name 
                fout="outfits/%s/gal_%s_%s_%s.fits"%(fldname[ii],fldname[ii],gid[ii],band[kk]) 
                print("writing", fout)
                hdu = fits.PrimaryHDU(gal_im);
                hdulist = fits.HDUList([hdu]);
                if(flag):
                    hdr = hdu.header
                    hdr.append(('EXPOSURE', '%d'%(exparr[kk])), end=True)
                hdulist.writeto(fout)

        else:  ## Do not make a FITS cutout if it exists already
            print("%s fits file exists already"%(foutg))
            continue

## Run this code faster by specifying Nproc (no. of processors)
jobs=[];
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p)
    p.start()
print( "#####################################################################")
