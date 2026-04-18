# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
import numpy
import scipy.interpolate as SI
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
from scipy.spatial import KDTree as KDTree
        
class ipoint:
    """
    Interpolated object to rapidly get values
      Note: grids should be ghosted in the +i,j+ directions
    """
    def __init__(self,xpt,ypt,xgrid,ygrid,xrng=None,yrng=None):
        

        self.value = 0.0
        self.onProc = False
        
        self.x = xpt
        self.y = ypt

        # Bounding box for determining proc location
        #if not xrng:
        #    onX = ( (self.x > xgrid.min() ) and  (self.x <= xgrid.max() ) )
        #else:
        #    onX = ( (self.x > xrng[0] ) and  (self.x <= xrng[1] ) )

        #if not yrng:
        #    onY = ( (self.y > ygrid.min() ) and  (self.y <= ygrid.max() ) )
        #else:
        #    onY = ( (self.y > yrng[0] ) and  (self.y <= yrng[1] ) )
        
        
        #if onX and onY:
        #    self.onProc = True
        #else:
        #    return
            
        self.ii1 = -1
        self.iin = -1

        self.jj1 = -1
        self.jjn = -1

        [self.nx,self.ny] = xgrid.shape
               
        self.weight  = []

        x_y = numpy.dstack([xgrid.ravel(),ygrid.ravel()])[0]

        pts = numpy.zeros( (2,1) )
        pts[0] = self.x*1.0
        pts[1] = self.y*1.0
        points = list(pts.transpose())

        mytree = KDTree(x_y)
        dist, indexes = mytree.query(points)

        ii = int(indexes[0]/self.ny)
        jj = indexes[0]%self.ny

        # 
        nonX = ( ((ii-1) < 0) or ( (ii+1) > (self.nx-1) ) )
        nonY = ( ((jj-1) < 0) or ( (jj+1) > (self.ny-1) ) )

        self.onProc = True
        if nonX or nonY:
            self.onProc = False
            return
        
        x0 = numpy.abs(xgrid[ii,jj] - self.x)
        y0 = numpy.abs(ygrid[ii,jj] - self.y)

        self.ii0 = ii
        self.jj0 = jj

        self.xlocal = numpy.array([ [xgrid[ii-1,jj-1],xgrid[ii,jj-1],xgrid[ii+1,jj-1]],
                                    [xgrid[ii-1,jj  ],xgrid[ii,jj  ],xgrid[ii+1,jj  ]],
                                    [xgrid[ii-1,jj+1],xgrid[ii,jj+1],xgrid[ii+1,jj+1]] ])

        self.ylocal = numpy.array([ [ygrid[ii-1,jj-1],ygrid[ii,jj-1],ygrid[ii+1,jj-1]],
                                    [ygrid[ii-1,jj  ],ygrid[ii,jj  ],ygrid[ii+1,jj  ]],
                                    [ygrid[ii-1,jj+1],ygrid[ii,jj+1],ygrid[ii+1,jj+1]] ])


        myWeights = []
        myIndices = []
        for io in range( self.xlocal.shape[0] ):
            for jo in range( self.xlocal.shape[1] ):
                test = numpy.zeros( [3,3] )
                test[io,jo] = 1.0

                pts = [ numpy.array([self.x]), numpy.array([self.y]) ]
                dc = getPoints(self.xlocal,self.ylocal,test,pts,'linear')
                if ( dc != 0.0 ):
                    myWeights.append( dc )
                    myIndices.append( [ ii-1+io  , jj-1+jo ] )
                

        self.indices = myIndices
        self.weights = myWeights
        #self.xlocal = self.xlocal[1:,1:]
        #self.ylocal = self.ylocal[1:,1:]
        #self.xlocal = self.xlocal[:-1,:-1]
        #self.ylocal = self.ylocal[:-1,:-1]
        

        self.getWeights()
        

            
            
        

    def interp(self,vals,imethod='linear'):


        if (not self.onProc):
            return 0.0
        
        ii = self.ii0
        jj = self.jj0
        vlocal = numpy.array([ [vals[ii-1,jj-1],vals[ii,jj-1],vals[ii+1,jj-1]],
                             [vals[ii-1,jj  ],vals[ii,jj  ],vals[ii+1,jj  ]],
                             [vals[ii-1,jj+1],vals[ii,jj+1],vals[ii+1,jj+1]] ])

        #vlocal = vlocal[1:,1:]
        #vlocal = vlocal[:-1,:-1]
        
        pts = [ numpy.array([self.x]), numpy.array([self.y]) ]
        dc = getPoints(self.xlocal,self.ylocal,vlocal,pts,imethod)

        return dc




    def interpFast(self,vals):


        if (not self.onProc):
            return 0.0
        
        ii = self.ii0
        jj = self.jj0
        vlocal = numpy.array([ [vals[ii-1,jj-1],vals[ii,jj-1],vals[ii+1,jj-1]],
                             [vals[ii-1,jj  ],vals[ii,jj  ],vals[ii+1,jj  ]],
                             [vals[ii-1,jj+1],vals[ii,jj+1],vals[ii+1,jj+1]] ])

        #vali = interpolateFast( vlocal.flatten(), self.vtx,self.SIweight)
        #def interpolateFast(values, vtx, wts):
        
        vals =  numpy.einsum('nj,nj->n', numpy.take(vlocal.flatten(), self.vtx), self.SIweight)
        
        return vals[0]
    
    def getWeights(self):
        
        X = self.xlocal
        Y = self.ylocal
        
        pts = [ numpy.array([self.x]), numpy.array([self.y]) ]
        
        xy=numpy.zeros([X.shape[0]*X.shape[1],2])
        xy[:,0]=X.flatten()
        xy[:,1]=Y.flatten()
        
        uv=numpy.zeros([1,2])
        uv[0,0]=pts[0].item()
        uv[0,1]=pts[1].item()
        
        vtx, wts = interp_weights(xy, uv)

        self.SIweight = wts
        self.vtx = vtx
    
def getPoints(x,y,f,pts,imethod):

    xx = x.ravel()
    yy = y.ravel()
    ff = f.ravel()
    xpts = pts[0].ravel()
    ypts = pts[1].ravel()
    
    dc = SI.griddata(numpy.array([xx,yy]).T,ff,numpy.array([xpts,ypts]).T,method=imethod)
        
    return dc[0]


def interp_weights(xy, uv,d=2):
    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = numpy.take(tri.simplices, simplex, axis=0)
    temp = numpy.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = numpy.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, numpy.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

