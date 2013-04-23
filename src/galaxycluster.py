import sys
import pyfits
import numpy as np

import stats

import pylab as pl

G= 1.3273e26 # cm^3/s/M_*
cm2Mpc=3.086e24 # cm/Mpc

def z2v(z, zc): 
    """Convert the redshift to km/s relative to the cluster center"""
    return 2.99792458e5*(z-zc)/(1+zc)

def calc_mass(sigma, h=0.7):
    """Calculates the mass of the clusters in Solar Masses"""
    #from finn et al. 2005
    return 1.2e15*(sigma/1000.0)**3/h


def calc_r200(sigma, z, h=0.7, o_b=0.3, o_l=0.7, o_r=0.0):
    """Calculates the r200 radius for a cluster in Mpc"""
    #from finn et al. 2005
    return 1.73*(sigma/1000.0)/h/(o_l+o_b*(1+z)**3+o_r*(1+z)**2)**0.5


def velocity_dispersion(z, niter=5, crad=None):
   """Use the biweight estimator to determine the velocity 
      dispersion for a distribution

   """
   zc=stats.biweight_location(z)
   v=z2v(z, zc)
   vc=stats.biweight_location(v)
   dv=stats.biweight_midvariance(v)
   if crad.any():
      r200=calc_r200(dv, zc)
      mask=crad<r200
   else:
      mask=vc*0.0+1.0
   for i in range(niter):
       vc=stats.biweight_location(v[mask], M=vc)
       dv=stats.biweight_midvariance(v[mask], M=vc)

   return vc, dv


def shifty_gapper(r, z, zc, vlimit=10000, ngap=30, glimit=1000):
   """Determine cluster membersip according to the shifty
      gapper technique

      The shifty gapper technique includes the following steps:
      1.  Remove outliers outside of +- vlimit
      2.  Locate the Ngap targets with radius close to r_i
      3.  Within this sample of N targets, identify sources with
          |v_pec| < |v_pec_i|
      4.  Meaure the largest gap within this subset of targets
      5.  If v_gap is larger than glimit, reject the source

      Parameters
      -----------
      vlimit: float
         Initial limit in velocity offset from the cluster center
      ngap: int
         Number of sources to use in in estimating the gap

      glimit:  float
         Maximum gap size for rejecting objects
 
      Parameters
      -----------
      incluster: ndarray
          Returns a boolean array where objects in the cluster have
          a value of True
   """
 
   #convert to the velocity scale
   v = z2v(z,zc)

   #limit the source to only sources within the vlimit
   vmask = abs(v) < vlimit

   nobj=len(r)
   incluster=np.zeros(nobj, dtype=bool)


   if nobj<ngap:
      raise Exception('Number of sources is less thant number of gap sources')

   for i in range(nobj):
     if abs(v[i])<vlimit:
       #find the ngap closest sources
       r_j=abs(r[vmask]-r[i]).argsort()
       vg=v[vmask][r_j[0:ngap]]

       #find the sources with |v_pec| < |v_pec_i|
       mask=abs(vg)<=abs(v[i])
       if mask.sum()>1:
          vg=vg[mask]
          #now sort these sources and find the biggest gap
          vg.sort()
          if np.diff(vg).max()<glimit: incluster[i]=True
       else:
          incluster[i]=True


   return incluster

def calcds(i, x, y, v, vc, sc, nobj=11):
    """Calculate the dressler-schectman test for source i using the following 
       equation:

       d**2=(11/s**2)*[(ave(v_local)-ave(v))**2+(sigma_local-sigma)**2]

      Parameters
      -----------
      i: int 
         Index of source to calculate DS value
      x: ndarray
         X-positions of the sources.  Can be in pixel or spherical coordinates
      y: ndarray
         Y-positions of the sources.  Can be in pixel or spherical coordinates
      v: ndarray
         Velocities  of the sources.  Can be in pixel or spherical coordinates
      vc: float
         Mean velocity of the cluster
      sc: float
         Velocity dispersion of the cluster
      nobj: int
         Number of nearby neighbors to average over

      Output
      -----------
      ds: float 
          Returns a float with the ds value 

    """
    sc=float(sc)
    val=nobj/sc**2
    #set the central value for the sources
    xc=x[i]
    yc=y[i]
 
    #find the distance to all of its neighbors
    dist=((xc-x)**2+(yc-y)**2)**0.5
    j=dist.argsort()

    #calculate the mean velocity and dispersion for the sub-group
    vl=varr[j[0:nobj]].mean()
    sl=varr[j[0:nobj]].std()
 
    #return the ds value
    return (val * ((vl-vc)**2+(sl-sc)**2))**0.5

