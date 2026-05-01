################################################################################
# Copyright (c) 2023, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
import numpy
import sys
try:
    # Public release of "leospy" only works for python3
    sys.path.append("/usr/apps/leos/latest_v8/lib/")
    import leospy
    
except:
    pass

from .pyrandaPackage import pyrandaPackage



class pyrandaLEOS(pyrandaPackage):
    """
    LEOS package api
    """
    def __init__(self,pysim):

        PackageName = "LEOS"
        pyrandaPackage.__init__(self,PackageName,pysim)
        self.pysim = pysim                

        self.nx = pysim.PyMPI.ax
        self.ny = pysim.PyMPI.ay
        self.nz = pysim.PyMPI.az
        self.shape = (self.nx, self.ny, self.nz)
        self.size = numpy.prod( self.shape )

        
        # Materials dictionary
        self.materials = {}
                
        
    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        #sMap['LEOS.pressure(']    = "self.packages['LEOS'].pressure("
        #sMap['LEOS.energy(']      = "self.packages['LEOS'].energy("
        #sMap['LEOS.soundSpeed(']  = "self.packages['LEOS'].soundSpeed("

        #sMap['LEOS.inv.pressure(']    = "self.packages['LEOS'].ipressure("
        #sMap['LEOS.inv.energy(']      = "self.packages['LEOS'].ienergy("
        #sMap['LEOS.inv.soundSpeed(']  = "self.packages['LEOS'].isoundSpeed("

        sMap['LEOS.pressure(']    = "self.packages['LEOS'].leos_lookup(False,'Pt',"
        sMap['LEOS.energy(']      = "self.packages['LEOS'].leos_lookup(False,'Et',"
        sMap['LEOS.soundSpeed(']  = "self.packages['LEOS'].leos_lookup(False,'Cs',"

        sMap['LEOS.inv.pressure(']  = "self.packages['LEOS'].leos_lookup(True,'Pt',"
        sMap['LEOS.inv.energy(']    = "self.packages['LEOS'].leos_lookup(True,'Et',"
        sMap['LEOS.inv.soundSpeed(']= "self.packages['LEOS'].leos_lookup(True,'Cs',"
               
        sMap['LEOS.temperature('] = "self.packages['LEOS'].temperature("
        
        self.sMap = sMap


    def leos_lookup(self,inverse,field,matID,rho,val):

        ieos = self.materials[matID]
        if not inverse:
            out = ieos[field].eval( rho.flatten(), val.flatten() )
        else:
            out = ieos[field].ieval( rho.flatten(), val.flatten() )

        return out.reshape(self.shape)
            
        
    def energy(self,matID,rho,Temp):

        ieos = self.materials[matID]        
        E = ieos['Et'].eval( rho.flatten() , Temp.flatten() )
        
        return E.reshape(self.shape)
    
    
    def pressure(self,matID,rho,Temp):

        ieos = self.materials[matID]        
        P = ieos['Pt'].eval( rho.flatten() , Temp.flatten() )
        
        return P.reshape(self.shape)

    def soundSpeed(self,matID,rho,Temp):

        ieos = self.materials[matID]        
        Cs = ieos['Cs'].eval( rho.flatten() , Temp.flatten() )
        
        return Cs.reshape(self.shape)   

    # Temperature is really a Et inverse lookup
    def temperature(self,matID,rho,ie):

        ieos = self.materials[matID]
        T = ieos['Et'].ieval( rho.flatten(), ie.flatten() )

        return T.reshape(self.shape)
    
        

    def addMaterial(self,EOS_number,function_names=['Pt','Et','Cs'],matID=None,dxdy=False):

        # Assign a name/ID for this material
        if not matID:
            matID = EOS_number
            
        self.materials[matID] = self.__setupLEOS__(EOS_number,function_names,dxdy=dxdy)
        


    def __setupLEOS__(self,EOS_File,fnames,dxdy=False):

        matnum = EOS_File

        # Set options
        lopts = leospy.LEOS_LookupOptions()
        if not dxdy:
            lopts.calculateFunction(True).calculateDFDX(False).calculateDFDY(False)
        else :
            lopts.calculateFunction(True).calculateDFDX(True).calculateDFDY(True)

        fopts = leospy.LEOS_FunctionOptions()
        fopts.interpolation(leospy.BICUBIC)    # Bi-cubic
        fopts.units(leospy.LEOS_UNITS_CGS)     # CGS units

        infile= "leos"
        eos = {}

        for func in fnames:

            lfunc = leospy.getFunction(infile, matnum, func, fopts)
            lfunc.initialize()

            # Save to dictionary
            eos[func] = leosFunc(lfunc,lopts,self.size)


        # Also get rho0 and T0 of this table
        db = leospy.getDatabase("leos")
        mat = db.getMaterial(matnum)
        eos['rho0'] = float( mat.getMetadata('rho0')[1])
        eos['T0']   = float( mat.getMetadata('t0')[1])    

        # Also compute the "fnames"_0 at the reference state
        rho0 = eos['rho0']
        T0   = eos['T0']
        for func in fnames:
            eos[func+"_0"] = eos[func].eval( rho0, T0 , N=1 )
        
        
        return eos

        



    
class leosFunc():
    def __init__(self,eosFunc,lopts,Npts):
        self.eosFunc = eosFunc
        self.lopts   = lopts
        self.name    = self.eosFunc.getName()
        self.Npts    = int(Npts)
        
    def eval(self,rho,T,Var=None,N=None,dxdy=False):

        # Allow for single pt. calls
        if not N:
            N = self.Npts

        # Force floats to arrays
        if type(rho) == type(0.0):
            rho = numpy.array( [rho] )

        if type(T) == type(0.0):
            T = numpy.array( [T] )        
        
        if hasattr(self.eosFunc, "eval_iDDDDDo"):
            x = leospy.VectorOfDoubles( rho )
            y = leospy.VectorOfDoubles( T   )
            if (type(Var) == type( numpy.array([]))) :
                f = leospy.VectorOfDoubles( Var   )
            else:
                f = leospy.VectorOfDoubles( rho*0.0   )
            dx = leospy.VectorOfDoubles( rho*0.0 )
            dy = leospy.VectorOfDoubles( rho*0.0 )
            
            self.eosFunc.eval_iDDDDDo( N, x, y, f, dx, dy, self.lopts )
            f = numpy.array( f )
            dx = numpy.array(dx)
            dy = numpy.array(dy)
        else:
            f, dx, dy = self.eosFunc.eval(int(N), rho.tolist(), T.tolist(), self.lopts)
            f = numpy.array(f)
            dx = numpy.array(dx)
            dy = numpy.array(dy)

        if not dxdy:
            return f

        else:
            return [ f, dx, dy ]

    

    def ieval(self,rho,E,Tg=None,N=None):

        if not N:
            N = self.Npts

        # Force floats to arrays
        if type(rho) == type(0.0):
            rho = numpy.array( [rho] )

        if type(E) == type(0.0):
            E = numpy.array( [E] )

        if hasattr(self.eosFunc, "inverseEval_iDDDDDo"):
            x = leospy.VectorOfDoubles( rho )
            E = leospy.VectorOfDoubles( E   )
            if (type(Tg) == type( numpy.array([]))) :
                T = leospy.VectorOfDoubles( Tg   )
            else:
                T = leospy.VectorOfDoubles( rho*0.0   )
            dx = leospy.VectorOfDoubles( rho*0.0 )
            dy = leospy.VectorOfDoubles( rho*0.0 )
            
            self.eosFunc.inverseEval_iDDDDDo( N, x, E, T, dx, dy, self.lopts )
            return numpy.array( T )

        T, dx, dy = self.eosFunc.inverseEval(int(N), rho.tolist(), E.tolist(), self.lopts)
        return numpy.array(T)       
        

        
