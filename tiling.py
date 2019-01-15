#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 13:53:13 2018

conda create --name nbodykit-env python=3
source activate nbodykit-env

@author: burtin
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import os,sys
import astropy.io.fits as fits
import astropy.coordinates as coord
import astropy.units as u
import mangle
#import pandas

#import nbodykit.transform as nk_tsf
#import nbodykit.cosmology.cosmology as nk_cosmo

import nbodykit.algorithms.fibercollisions as fc
#import dask.array as da

#myCLPT    = RSDmodels.Model('WMAP7')
power=1


tracer='ELG'
if tracer=='QSO':
    dirEZ=os.environ['CLU']+'/EZ-mock/v1.8/'
    fileEZ=dirEZ+'zevoEZmock_QSO_v1.8_veto_ngc_0001.dat'
elif tracer=='ELG':
    dirEZ=os.environ['CLU']+'/EZmock/ELG_v1.1/NGC/'
    fileEZ=dirEZ+'EZmock_eboss_ELG_v1.1_NGC_0001.dat'
else:
    sys.exit()
catLoad=np.loadtxt(fileEZ)

chunk='25'
geomFile='geometry-eboss'+chunk+'.ply'
tilesFile='tiles-eboss'+chunk+'.par'

def inChunk(geomFile,cat):
    mng=mangle.Mangle(geomFile)
    object_mask = mng.get_polyids(cat[:,0],cat[:,1])
#    print(object_mask)
    return cat[(object_mask>0)]

catEZ=inChunk(geomFile,catLoad)
    
print('len catEZ',len(catEZ))
##cat=coord.SkyCoord(catEZ[:,0]*u.degree,catEZ[:,1]*u.degree)
#ida,idb,sep,dist=coord.search_around_sky(catEZ,cat,1.5*u.degree)
#

def getPlates(plateFile,cat):
    plChunk=np.loadtxt(plateFile,skiprows=41,usecols=(3,4))
    catPlates=coord.SkyCoord(plChunk[:,0]*u.degree,plChunk[:,1]*u.degree)
    catObj=coord.SkyCoord(cat[:,0]*u.degree,cat[:,1]*u.degree)
    objIndex,plateIndex,sep,dist=coord.search_around_sky(catObj,catPlates,1.49*u.degree)
#    print('objIndex:',objIndex,'plateIndex:',plateIndex)
#    print('lencat',len(catObj),'objIndex',len(objIndex),'platejIndex',len(plateIndex))
    objId,counts=np.unique(objIndex,return_counts=True)
#    print(len(catObj),len(catPlates))
#    print(len(objId),objId,counts)
    val,countsC=np.unique(counts,return_counts=True)
    print("objects belonging to",val,"plates, fibers:",countsC)
    return counts

def assign(cat,nPlates):
    '''
    find a decollided set of input catalog cat and discards object belonging to n <= nPlates
    output :
        - catColl : catalog of objects that did not get a fiber for the next round
        - catDisc : catalog of objects that could never get a fiber
    '''
    print('------------------------------------------------------------------------')
    print('input catalog:')
    counts=getPlates(tilesFile,cat)
    coll    = fc.FiberCollisions(cat[:,0],cat[:,1])
    coll.run()
    noFiber    = np.array(coll.labels['Collided'])
    catNoFiber = cat[(noFiber==1)]
    print('collided set:')
    counts     = getPlates(tilesFile,catNoFiber)
    catColl    = catNoFiber[counts>nPlates]
    discarded  = counts<=nPlates
    catDisc    = catNoFiber[discarded]
#    plotGeometry(geomFile)
#    plt.scatter(cat[:,0]    ,cat[:,1] ,s=50,color='blue')
#    plt.show()
    print('len(catcoll):',len(catColl),'discarded',np.sum(discarded))
    return catColl,catDisc

#cat=coord.SkyCoord(catEZ[:,0]*u.degree,catEZ[:,1]*u.degree)
#ida,idb,sep,dist=cat.search_around_sky(cat,62.*u.arcsecond)

def plotGeometry(geomFile):
    mng=mangle.Mangle(geomFile)
    polys = mng.graphics()
#mng.drawpolys()
    for poly in polys:
        plt.plot(poly['ra'],poly ['dec'],lw=1,color='black')#,color=poly['weight']
#
#print(np.shape(ida),np.shape(idb))

#print('parent catalog:')
#counts=getPlates(tilesFile,catEZ)

#coll = fc.FiberCollisions(catEZ[:,0],catEZ[:,1])
#coll.run()
#collisionGroup = np.array(coll.labels['Label'])
#noFiber        = np.array(coll.labels['Collided'])
#NeighborID     = np.array(coll.labels['NeighborID'])
#
#catCollision= catEZ[(collisionGroup>0)]
#catFiber    = catEZ[(noFiber==0)]
#onePlate    = catEZ[counts==1]
#twoPlates   = catEZ[counts==2]
#threePlates = catEZ[counts==3]

# second round, only object which did not get a fiber and are in n>1 plate are considered
#print('parent catalog:')
catColl ,catDisc  = assign(catEZ,   1)
catColl2,catDisc2 = assign(catColl, 2)
catColl3,catDisc3 = assign(catColl2,3)
#sys.exit()

plt.subplot(121)
plotGeometry(geomFile)
plt.scatter(catEZ[:,0]    ,catEZ[:,1] ,s=50,color='pink')
plt.scatter(catColl[:,0]    ,catColl[:,1] ,s=50,color='blue')
plt.scatter(catColl2[:,0]   ,catColl2[:,1],s=20,color='red')
plt.scatter(catColl3[:,0]   ,catColl3[:,1],s=100,color='green')

plt.subplot(122)
plotGeometry(geomFile)
plt.scatter(catDisc[:,0]    ,catDisc[:,1] ,s=50,color='blue')
plt.scatter(catDisc2[:,0]   ,catDisc2[:,1],s=50,color='red')
plt.scatter(catDisc3[:,0]   ,catDisc3[:,1],s=50,color='green')
plt.show()

#plt.scatter(catEZ[:,0]       ,catEZ[:,1],s=50,color='blue')
#plt.scatter(onePlate[:,0]    ,onePlate[:,1],s=50,color='cyan')
#plt.scatter(twoPlates[:,0]   ,twoPlates[:,1],s=50,color='magenta')
#plt.scatter(catCollision[:,0],catCollision[:,1],s=20,color='red')
#plt.scatter(catFiber[:,0]    ,catFiber[:,1],s=5,color='yellow')
#plotGeometry(geomFile)
#plt.show()
#

sys.exit()


fc=62.*u.arcsecond
print(fc)

nCollisionGroup=0

#collisionGroup=np.zeros(len(catEZ))

collidedTargetsTreated=[]
collidedTargetsGroup=[]

cond = sep>0.*u.arcsec

def getGroup(obj,targetList,groupList):
    for i in np.arange(len(targetList)):
        if targetList[i]==obj:
            return groupList[i]
    return None

for i in np.arange(len(catEZ)):
    if cond[i]==True:    
#        print('--------------------------------')
#        print(i,ida[i],idb[i],sep[i],dist[i],cond[i])
        isInCollided=np.isin(collidedTargetsTreated,ida[i])
#        print(isInCollided)
#        print(collidedTargetsTreated[isInCollided])
#        
#        cut=(isInCollided==True)
#        print(collidedTargetsGroup[cut])
#        print(isInCollided)
        if np.sum(isInCollided):
            if np.isin(idb[i],collidedTargetsTreated):
                pass
            else:
                collidedTargetsTreated.append(idb[i])
                collidedTargetsGroup.append(nCollisionGroup)

        else:
            if np.isin(idb[i],collidedTargetsTreated):
                collidedTargetsGroup.append(getGroup(idb[i],collidedTargetsTreated,collidedTargetsGroup)) 
            else:
                nCollisionGroup+=1
                collidedTargetsTreated.append(ida[i])
                collidedTargetsTreated.append(idb[i])
                collidedTargetsGroup.append(nCollisionGroup)                
                collidedTargetsGroup.append(nCollisionGroup)
            
#        print(collidedTargetsTreated)  
            
#        print(collidedTargetsTreated)  

print('groups',collidedTargetsGroup)  
a,counts=np.unique(collidedTargetsGroup,return_counts=True)
print(counts)
b,coll=np.unique(counts,return_counts=True)
print(b,coll)
#collided=np.where(sep.to(u.arcsec)<62.*u.arcsec)
#print(idx[collided])
    
print(np.sum(cond),len(catEZ))


cap='NGC'

#catD   =fits.open(os.environ['CLU']+'/qso_catalogs/4-year/eBOSS_QSO_clustering_'+cap+'_test.dat.fits',memmap=True)[1].data
#catR   =fits.open(os.environ['CLU']+'/qso_catalogs/4-year/eBOSS_QSO_clustering_'+cap+'_test.ran.fits',memmap=True)[1].data
#
#cosmo=nk_cosmo.Cosmology()
#
#size=len(catD)
##size=10
#cart=nk_tsf.SkyToCartesian(da.from_array(catD['RA'],chunks=(size)),da.from_array(catD['DEC'],chunks=(size)),da.from_array(catD['Z'],chunks=(size)),cosmo)
#
#np_cart=np.array(cart)
#print np_cart
