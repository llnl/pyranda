################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
from __future__ import print_function
from mpi4py import MPI
import numpy 
from numpy import nan
import re
import sys,os
import time,random
import glob
import inspect
import scipy

from .pyrandaMPI   import pyrandaMPI
from .pyrandaVar   import pyrandaVar
from .pyrandaEq    import pyrandaEq
from .pyrandaMesh  import pyrandaMesh, defaultMeshOptions
from .pyrandaIO    import pyrandaIO
from .pyrandaPlot  import pyrandaPlot
from .pyrandaUtils import *
from .pyrandaTex   import pyrandaTex
from .pyrandaElliptical import pyrandaPoisson
                                              
class pyrandaSim:

    def __init__(self,name,meshOptions,silent=False):

        self.name = name
        self.silent = silent
        self.mesh = pyrandaMesh()

        # Parse shorthand mesh description if string
        if type(meshOptions) == type(''):
            self.mesh.makeMeshStr( meshOptions )
            meshOptions = self.mesh.options

        defMeshOptions = defaultMeshOptions()
        defMeshOptions.update( meshOptions )

        meshOptions =      defMeshOptions
        self.meshOptions = defMeshOptions
        self.mesh.options =defMeshOptions 

        if 'coordsys' in meshOptions:
            self.mesh.coordsys = meshOptions['coordsys']
        else:            
            raise ValueError('No suitable mesh type specified.')

        nx = meshOptions['nn'][0]
        ny = meshOptions['nn'][1]
        nz = meshOptions['nn'][2]

        dx = (meshOptions['xn'][0]-meshOptions['x1'][0])/max(nx-1,1)
        dy = (meshOptions['xn'][1]-meshOptions['x1'][1])/max(ny-1,1)
        dz = (meshOptions['xn'][2]-meshOptions['x1'][2])/max(nz-1,1)
        periodic = meshOptions['periodic']                    

        self.dx = dx
        self.dy = dy
        self.dz = dz
                
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.npts = nx*ny*nz

        # Initialize some lists/dictionaries
        self.equations = []
        self.conserved = []
        self.variables = {}
        self.initial_conditions = []
        self.nPDE = 0
        self.nALG = 0
        self.nELP = 0

        # Some solver settings
        #self.linearSolver = "scipy"  # User can specify "scipy" or "petsc"
        
        self.mesh.options = meshOptions
        self.mesh.dims  = meshOptions['dim']

        if ("comm" in meshOptions):
            my_comm = meshOptions['comm']
        else:
            my_comm = None
            
        # Setup mpi (and parcop)
        self.PyMPI = pyrandaMPI( self.mesh, comm=my_comm )
        self.mesh.PyMPI = self.PyMPI

        # Create the mesh
        self.mesh.makeMesh()
        
        # IO setup
        self.PyIO = pyrandaIO( self.name, self.PyMPI)

        # Plotting setup
        self.plot = pyrandaPlot( self )
        
        # Grab some scalars off of parcop.mesh
        self.zero = self.PyMPI.emptyScalar()

        # Package info
        self.packages = {}
        self.packagesRestart = []

        # User functions
        self.userDefined = {}
        
        # Compute sMap
        self.get_sMap()

        # Time
        self.time = 0.0
        self.cycle = 0
        self.vizDumpHistory = []
        self.deltat = 0.0

        
        # Print startup message
        self.iprint( code() , 4000 )
        self.iprint( version() , 1000)
        self.iprint( icopyright(), 1000 )

        
    def iprint(self,sprnt,wpm=0):
        if self.PyMPI.master and (not self.silent):

            if wpm == 0:
                print(sprnt)
                sys.stdout.flush()
            else:            
                typing_speed = wpm                
                for l in sprnt:
                    sys.stdout.write(l)
                    sys.stdout.flush()
                    time.sleep(random.random()*10.0/typing_speed)
                print('')
    
        
    def addPackage(self,package):
        
        self.packages[package.name] = package
        self.packagesRestart.append( package.__module__ )
        package.get_sMap()
        self.sMap.update( package.sMap )

    def addUserDefinedFunction(self,name,equation):
        if name in self.sMap:
            self.iprint("Error: Cant add user-defined equation '%s', it already exists" % name)
        self.userDefined[name] = equation
        imap = { "%s(" % name: "self.userDefined['%s'](self," % name  }
        self.sMap.update( imap )
    
        
    def allocate(self):
        """
        Loop over vars and allocate the data
        """
        for v in self.variables:
            self.variables[v].__allocate__(self.PyMPI)
        
    
    def var(self,name):
        if name in self.variables.keys():
            return self.variables[name]
        else:
            raise ValueError('Error: variable name: %s not found in database' % name)

    def x(self):
        return self.variables['meshx']

    def y(self):
        return self.variables['meshy']

    def z(self):
        return self.variables['meshz']
        
    def eval(self,expression):
        return eval(fortran3d(expression,self.sMap))

    def parse(self,expression,parseDict=None):

        if parseDict:
            expression = self.updateStringFromDictionary( expression, parseDict, "parse")
        
        exec(fortran3d(expression,self.sMap))
    
    def addVar(self,name,kind=None,rank='scalar'):
        if name not in self.variables:
            var = pyrandaVar(name,kind,rank)
            self.variables[name] = var

    def addEqu(self,equation):

        peq = pyrandaEq(equation,self.sMap,self.PyMPI) 
        self.equations.append( peq )
        if peq.kind == 'PDE':
            self.nPDE += 1
            self.conserved.append( peq.LHS[0] )
            if len( peq.LHS ) > 1:
                self.iprint('Warning... only single return values for PDEs allowed')
                exit()
        elif peq.kind == 'ALG':
            self.nALG += 1
        elif peq.kind == 'ELP':
            self.nELP += 1
        else:
            raise ValueError('No suitable equation type specified for %s' % peq.eqstr)

    def updateStringFromDictionary(self,string,dictionary,handle):

        # Apply a dictionary to the string
        for item in dictionary:
            if item in string:
                string = string.replace( item , str( dictionary[item] ) )
            else:
                self.iprint("Notice: Given %s dictionary value: '%s' does not exists in the %s string." % (handle,item,handle) )

        return string        

        
    def EOM(self,eom,eomDict=False):
        """
        Higher level wrapper to make equations of motion from a single string
        Will add eqautions and variables as needed.
        """
        
        # Apply a dictionary to the string
        if eomDict:
            eom = self.updateStringFromDictionary( eom, eomDict, "EOM")
        
        self.eom = eom
        
        # Split up the equation lines
        eom_lines = filter(None,eom.split('\n'))
        eom_lines = [el for el in eom_lines if el.strip()[0] != '#']  # Comments work
        var_names = []

        #### Get unique set of variables - scalars ####
        for eq in eom_lines:
            evars = findVar( eq,'scalar')
            var_names += evars
        var_names = list(set(var_names))
        for evar in var_names:
            self.addVar(evar,kind='conserved')   # Todo: classify variables

        #### Get unique set of variables - vectors ####
        var_names = []
        for eq in eom_lines:
            evars = findVar( eq,'vector')
            var_names += evars
        var_names = list(set(var_names))
        for evar in var_names:
            self.addVar(evar,rank='vector',kind='conserved')   
                                
        for variable in self.variables:
            self.iprint('Adding variables: %s' %  variable )

        #### Add equations of motion ####
        self.iprint('Adding equations of motion: ')

        eqN = 1
        for eq in eom_lines:
            self.addEqu( eq )
            self.iprint( '(%s)   %s' % (eqN,eq) )
            eqN += 1
            

        # Set up the variables in memory
        self.allocate()

        # Mesh vars point to mesh.coords
        self.variables['meshx'] = self.mesh.coords[0]
        self.variables['meshy'] = self.mesh.coords[1]
        self.variables['meshz'] = self.mesh.coords[2]

        
            
    def checkForNan(self,names=[]):

        nans = False
        svars = ''
        if not names:
            names = self.variables
        for ivar in names:
            try:
                myvar = self.variables[ivar]
                if numpy.isnan(myvar.data).any():
                    nans = True
                    svars += ivar + ' '
            except:
                self.iprint("%s is not a variable" % ivar)

        return svars

    def setIC(self,ics,icDict=False,addOnly=False,verbose=False):
        """
        Evaluate the initial conditions and then update variables
        """

        # Apply a dictionary to the string
        if icDict:
            ics = self.updateStringFromDictionary( ics, icDict, "ICs")
            
        self.ics = ics
        
        # Split up the equation lines
        ic_lines = filter(None,ics.split('\n'))
        ic_lines = [el.replace(' ','') for el in ic_lines ]  # Comments work
        ic_lines = filter(None,ic_lines)
        ic_lines = [el for el in ic_lines if el.strip()[0] != '#']  # Comments work

        var_names = []

        #### Get unique set of variables ####
        for eq in ic_lines:
            evars = findVar( eq,'scalar')
            var_names += evars

        var_names = list(set(var_names))

        for evar in var_names:
            self.addVar(evar,kind='conserved')   # Todo: classify variables
        self.allocate()


        for ic in ic_lines:
            self.initial_conditions.append( ic )
        
        # For addOnly option (used for restarts) dont evaluate or update
        if addOnly:
            return        
        
        # Actually compute the Initial Conditions
        for ic in ic_lines:
            ic_mod = ic #+ '+self.emptyScalar()'

            if verbose:
                self.iprint( ic )
                
            try:
                exec(fortran3d(ic_mod,self.sMap))
            except Exception as e:
                self.iprint("Error: cant parse following string from ICs")
                self.iprint(ic_mod)
                self.iprint("Error message: %s" % e )
                exit()
        
            
        # Check for nans here
        snans = self.checkForNan( var_names )
        if snans:
            self.iprint("Found some nans in inits: %s" % snans)
            exit()

        self.updateVars(verbose=verbose)
        snans = self.checkForNan()
        if snans:
            self.iprint("Found some nans in Update after init: %s" % snans)
            exit()


        
    def emptyScalar(self,val=0.0):
        return self.PyMPI.emptyScalar() + val
        
    def updateFlux(self):
        #
        ncons = self.nPDE
        shape = self.mesh.shape[:]
        shape.append( ncons )        
        # Only ddt() terms are extracted
        flux = {} #numpy.asfortranarray( numpy.zeros( shape ) )
        #ieq = 0
        for eqo in self.equations:
            eq = eqo.eqstr
            if ( eqo.kind == 'PDE' ):
                lhs = eqo.LHS[0]
                Srhs = eq.split('=')[1]  # This is a string to evaluate
                flux[lhs] = eqo.RHS(self)

        return flux

    def updateVars(self,verbose=False):
        #
        # Update the equations
        #
        for eq in self.equations:
            if not eq.active:
                continue

            if ( eq.kind == 'ALG'):
            
                if verbose:
                    self.iprint( eq.eqstr )

                rhs = eq.RHS(self)

                # If no LHS , then on to the next eq
                if not eq.LHS:
                    continue
                
                if eq.rank == 1:
                    self.variables[eq.LHS[0]].data = rhs
                else:
                    for ii in range(len(rhs)):
                        self.variables[eq.LHS[ii]].data = rhs[ii]

            # Solve an elliptical PDE
            #if ( eq.kind == 'ELP'):

                # Setup solver
                #if not eq.solver:
                #    eq.solver = pyrandaPoisson(None,self,solver_type=self.linearSolver)
                
                # \Delta \phi = rhs
                #rhs = eq.RHS(self)
                #self.variables[eq.LHS[0]].data = eq.solver.solve(rhs)

        
        
    def write(self,wVars=[]):
        """ 
        Write viz file 
        """
        procInt = 6
        visInt  = 7
        if not wVars:
            wVars = self.conserved

        dumpFile = 'vis' + str(self.cycle).zfill(visInt)
        
        if self.PyMPI.master == 1:
            try:
                os.mkdir(os.path.join(self.PyIO.rootname, dumpFile))
            except:
                pass

        self.PyMPI.comm.barrier()   # Wait for directory to be made
        rank = self.PyMPI.comm.rank
        dumpFile = os.path.join( self.PyIO.rootname,
                                 dumpFile,
                                 'proc-%s.%s' % (str(rank).zfill(procInt),str(self.cycle).zfill(visInt)))


        suff = 'vtk'
        self.PyIO.makeDumpVTK(self.mesh,self.variables,wVars,dumpFile,self.cycle, self.time)

        self.vizDumpHistory.append( [self.cycle, self.time] )

        # Write .visit file
        if self.PyMPI.master == 1:
            vid = open( os.path.join(self.PyIO.rootname, 'pyranda.visit' ) , 'w')
            vid.write("!NBLOCKS %s \n" % self.PyMPI.comm.Get_size() )
            for vdump in self.vizDumpHistory:
                iv = vdump[0]
                vid.write("!TIME %s \n" % vdump[1] )  # DOES NOT WORK, EITHER.
                for p in range(self.PyMPI.comm.Get_size()):
                    vid.write("%s\n" % os.path.join('vis' + str(iv).zfill(visInt),
                                                    'proc-%s.%s.%s' % (str(p).zfill(procInt),str(iv).zfill(visInt),suff ) ) )
            vid.close()
                    
            


    def writeRestart(self,ivars=None,suffix=None,wVars=[]):
        """
        -writeRestart-
        Description - Main driver to write the entire pyrandaSim state
        in parallel for later use at restart.
        """
        if not wVars:
            wVars = self.variables.keys()

        # Use cycle number if no suffix is given
        restartInt = 7
        if not suffix:
            suffix = '_' + str(self.cycle).zfill(restartInt)

        # Prep directory
        dumpDir = os.path.join(self.PyIO.rootname, "restart" + suffix)
        if self.PyMPI.master == 1:
            try:
                os.mkdir( dumpDir )
            except:
                pass

        self.PyMPI.comm.Barrier()
            
        ## Persistent data ##
        serial_data = {}

        supported_types = []
        supported_types.append(type(''))
        supported_types.append(type(0))
        supported_types.append(type(1.0))
        supported_types.append(type([]))
        supported_types.append(type({}))
        supported_types.append(type(None))
        supported_types.append(type(numpy.float64))
        supported_types.append(type(True))
        
        # Check if ivars dictionary is passed
        if ivars:
            # Sanitize "ivars" (dont allow any
            serial_data["local_vars"] = {}
            serial_data['local_vars']['numpyArrays'] = {}
            
            # For ivars that are numpy arrays, save them individually
            for vv in ivars:
                itype = type(ivars[vv])

                isType = False
                for typ in supported_types:
                    isType = isType or (typ == itype)
                
                if itype == type(numpy.ones(1)):
                    filename = os.path.join(dumpDir,"%s" % vv)
                    serial_data['local_vars']['numpyArrays'][vv] = filename
                    if self.PyMPI.master:
                        numpy.save( filename, ivars[vv] )
                elif isType:
                    serial_data["local_vars"][vv] = ivars[vv]
                #else:
                #    self.iprint("Warning: The variable %s is type %s is not supported by restarts." % (vv,itype))


                
        # Original domain-decomp
        serial_data['mesh'] = self.meshOptions.copy()
        serial_data['EOM']  = self.eom
        serial_data['ICs']  = self.ics
        serial_data['decomp'] = [ self.PyMPI.px, self.PyMPI.py, self.PyMPI.pz ]
        serial_data['procMap'] = self.PyMPI.procMap
        serial_data['packages'] = self.packagesRestart
        serial_data['time'] = self.time
        serial_data['deltat'] = self.deltat
        serial_data['cycle']  = self.cycle
        serial_data['vizDumpHistory'] = self.vizDumpHistory
        # Added by D. Lavacot 04/08/2022: save sMap to restart file so that user added functions are also saved
        serial_data['sMap'] = self.sMap  
        
        # Serialize the mesh function (if it exists)
        if 'function' in serial_data['mesh']:
            if serial_data['mesh']['function']:
                serial_data['mesh']['function'] = None

        if 'comm' in serial_data['mesh']:
            serial_data["mesh"]["comm"] =  "user specified"
                
                
        # Variable map
        serial_data['vars'] = {}
        cnt = 0
        for ivar in sorted(wVars):   # BJO: was getting proc dependent ordering here, force it to be uniform
            serial_data['vars'][ivar] = cnt
            cnt += 1

        PO = numpy.set_printoptions()
        numpy.set_printoptions(threshold=sys.maxsize)
        numpy.set_printoptions(linewidth=sys.maxsize)
        if self.PyMPI.master == 1:
            fid = open(os.path.join(dumpDir,'serial.dat'),'w')
            fid.write(str(serial_data))
        numpy.set_printoptions(PO)
            
            
        # Parallel
        # Variables
        DATA = self.PyMPI.emptyVector(len(wVars))
        for ivar in sorted(wVars):   # BJO: was getting proc dependent ordering here, force it to be uniform
            DATA[:,:,:,serial_data['vars'][ivar]] = self.variables[ivar].data

            
        # Write this big thing
        rank = self.PyMPI.comm.rank
        procFile = open(os.path.join(dumpDir,'proc-%s.bin' % str(rank).zfill(5)),'w')
        DATA.tofile(procFile)
        procFile.close()
        

    def writeGrid(self):
        """ 
        Write grid file 
        """
        wlen = 3
        shape = list(self.mesh.shape)
        shape.append( wlen )
        iodata = numpy.zeros( shape )
        for i in range( 3 ):
            iodata[:,:,:,i] = self.mesh.coords[i].data

        dumpName = 'grid' 
        self.PyIO.makeDump(iodata,dumpName)

        
                                        
    def ddx(self,val):
        if self.nx <= 1:
            return 0.0
        return self.PyMPI.der.ddx( val )

    def ddy(self,val):
        if self.ny <= 1:
            return 0.0
        return self.PyMPI.der.ddy( val )

    def ddz(self,val):
        if self.nz <= 1:
            return 0.0
        return self.PyMPI.der.ddz( val )
        
    def interp_z2fx(self,val):
        if self.nx <= 1:
            return val
        return self.PyMPI.der.interp_z2fx( val )

    def interp_z2fy(self,val):
        if self.ny <= 1:
            return val
        return self.PyMPI.der.interp_z2fy( val )

    def interp_z2fz(self,val):
        if self.nz <= 1:
            return val
        return self.PyMPI.der.interp_z2fz( val )

    def dd4x(self,val):
        return self.PyMPI.der.dd4x( val )

    def dd4y(self,val):
        return self.PyMPI.der.dd4y( val )

    def dd4z(self,val):
        return self.PyMPI.der.dd4z( val )

    def dd8x(self,val):
        return self.PyMPI.der.dd8x( val )

    def dd8y(self,val):
        return self.PyMPI.der.dd8y( val )

    def dd8z(self,val):
        return self.PyMPI.der.dd8z( val )

    def div(self,f1,f2=None,f3=None):

        if (type(f2) == type(None) and type(f3) == type(None)):
            if self.nx > 1:
                return self.PyMPI.der.div(f1,self.zero,self.zero)
            if self.ny > 1:
                return self.PyMPI.der.div(self.zero,f1,self.zero)
            if self.nz > 1:
                return self.PyMPI.der.div(self.zero,self.zero,f1)
            # If we're actually 0D
            if ((self.nx <= 1) and (self.ny <= 1) and (self.nx <= 1)):
                return 0.0

        elif type(f3) == type(None):
            if (self.nx > 1) and (self.ny > 1):
                return self.PyMPI.der.div(f1,f2,self.zero)
            if (self.nz > 1) and (self.ny > 1):
                return self.PyMPI.der.div(self.zero,f1,f2)
            if (self.nz > 1) and (self.nx > 1):
                return self.PyMPI.der.div(f1,self.zero,f2)
            # If we're actually 1D
            if self.nx > 1:
                return self.PyMPI.der.div(f1,self.zero,self.zero)
            if self.ny > 1:
                return self.PyMPI.der.div(self.zero,f1,self.zero)
            if self.nz > 1:
                return self.PyMPI.der.div(self.zero,self.zero,f1)
            # If we're actually 0D
            if ((self.nx <= 1) and (self.ny <= 1) and (self.nx <= 1)):
                return 0.0
        else:            
            return self.PyMPI.der.div(f1,f2,f3)

    def divT(self,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz):
        return self.PyMPI.der.divT( fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz )
        
    def grad(self,val):
        #return [self.ddx(val),self.ddy(val),self.ddz(val)]
        return self.PyMPI.der.grad( val )    

    def laplacian(self,val):
        return self.PyMPI.der.laplacian( val )

    def ring(self,val):
        return self.PyMPI.der.ring( val )

    def ringV(self,vx,vy,vz):
        return self.PyMPI.der.ringV( vx,vy,vz )

    def filterx(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_x(val, f_tilde)
        return f_tilde

    def filtery(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_y(val, f_tilde)
        return f_tilde

    def filterz(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_z(val, f_tilde)
        return f_tilde 


    def filter(self,val):
        return self.PyMPI.fil.filter(val)

    def filterOld(self,val):
        if self.nx > 1:
            f1 = self.filterx(val)
        else:
            f1 = val
        if self.ny > 1:
            f2 = self.filtery(f1)
        else:
            f2 = f1
        if self.nz > 1:
            f1 = self.filterz(f2)
        else:
            f1 = f2
        return f1

    def getVar(self,vname):
        return self.PyMPI.getVar(vname)    
                                 
    def gfilterx(self,val):
        return self.PyMPI.gfil.filterDir(val, 1)

    def gfiltery(self,val):
        return self.PyMPI.gfil.filterDir(val, 2)

    def gfilterz(self,val):
        return self.PyMPI.gfil.filterDir(val, 3)

    def gfilter(self,val):
        return self.PyMPI.gfil.filter(val)
    
    def gfilterOld(self,val):
        if self.nx > 1:
            f1 = self.gfilterx(val)
        else:
            f1 = val
        if self.ny > 1:            
            f2 = self.gfiltery(f1)
        else:
            f2 = f1
        if self.nz > 1:
            f1 = self.gfilterz(f2)
        else:
            f1 = f2
        return f1

    def step(self,nsteps=1):
        for tt in range(nsteps):
            self.rk4(self.time,self.deltat)
    
    
    def rk4(self,time,dt):

        Ark = [0.0]*5
        Ark[0] = 0.0;
        Ark[1] = -6234157559845./12983515589748.;
        Ark[2] = -6194124222391./4410992767914.;
        Ark[3] = -31623096876824./15682348800105.;
        Ark[4] = -12251185447671./11596622555746.;

        Brk = [0.0]*5
        Brk[0] = 494393426753./4806282396855.;
        Brk[1] = 4047970641027./5463924506627.;
        Brk[2] = 9795748752853./13190207949281.;
        Brk[3] = 4009051133189./8539092990294.;
        Brk[4] = 1348533437543./7166442652324.;

        eta = [0.0]*5
        eta[0] = 494393426753./4806282396855.;
        eta[1] = 4702696611523./9636871101405.;
        eta[2] = 3614488396635./5249666457482.;
        eta[3] = 9766892798963./10823461281321.;
        eta[4] = 1.0;

        #	Initialize some intermediate arrays
        ncons = self.nPDE
        shape = tuple(self.mesh.shape)

        tmp1 = {}
        tmp2 = {}
        PHI  = {}
        for U in self.conserved: 
            tmp1[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            tmp2[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            PHI[U]  = numpy.asfortranarray( numpy.zeros(  shape ) )
        
        # Get primative flow variables
        #self.updateVars()
        time_i = time
        self.deltat = dt
        for ii in range(5):
            #    ii
            FLUX = self.updateFlux()
            for U in self.conserved:
                tmp1[U] =  Ark[ii]*PHI[U]
                PHI[U]  =  dt*FLUX[U] + tmp1[U]
                tmp2[U] =  Brk[ii]*PHI[U]
                self.variables[U].data =  self.variables[U].data + tmp2[U]
            time = time_i + eta[ii]*dt
            self.time = time
            self.updateVars()

        self.cycle += 1
        return time

    def get_sMap(self):
        sMap = {}
        

        # Simple find/replace mappings
        sMap['div(' ] = 'self.div('
        sMap['divT(' ] = 'self.divT('
        sMap['ddx(' ] = 'self.ddx('
        sMap['ddy(' ] = 'self.ddy('
        sMap['ddz(' ] = 'self.ddz('
        sMap['iz2fx(' ] = 'self.interp_z2fx('
        sMap['iz2fy(' ] = 'self.interp_z2fy('
        sMap['iz2fz(' ] = 'self.interp_z2fz('
        sMap['fbar('] = 'self.filter('
        sMap['gbar('] = 'self.gfilter('
        sMap['gbarx('] = 'self.gfilterx('
        sMap['gbary('] = 'self.gfiltery('
        sMap['gbarz('] = 'self.gfilterz('
        sMap['grad('] = 'self.grad('
        sMap['simtime'] = 'self.time'
        sMap['deltat'] = 'self.deltat'
        sMap['lap(' ] = 'self.laplacian('
        sMap['ring(' ] = 'self.ring('
        sMap['ringV(' ] = 'self.ringV('
        sMap['dd4x(' ] = 'self.dd4x('
        sMap['dd4y(' ] = 'self.dd4y('
        sMap['dd4z(' ] = 'self.dd4z('
        sMap['dd8x(' ] = 'self.dd8x('
        sMap['dd8y(' ] = 'self.dd8y('
        sMap['dd8z(' ] = 'self.dd8z('
        sMap['sum(' ] = 'self.PyMPI.sum3D('
        sMap['mean('] = '1.0/float(self.npts) * self.PyMPI.sum3D('
        sMap['mean3D('] = 'self.PyMPI.mean3D('
        sMap['sign(' ] = 'numpy.sign('
        sMap['dot(' ] = 'numpy.dot('
        sMap['abs(' ] = 'numpy.abs('
        sMap['sqrt(' ] = 'numpy.sqrt('
        sMap['sin('] = 'numpy.sin('
        sMap['cos('] = 'numpy.cos('
        sMap['tanh('] = 'numpy.tanh('
        sMap['exp('] = 'numpy.exp(' 
        sMap['where('] = 'numpy.where('
        sMap['max('] = 'self.PyMPI.max3D('
        sMap['min('] = 'self.PyMPI.min3D('
        sMap['3d('] = 'self.emptyScalar('
        sMap['pi'] = 'numpy.pi'
        sMap['meshVar('] = 'self.PyMPI.getVar('
        sMap['gridLen'] = 'self.mesh.GridLen'
        sMap['random3D()'] = 'numpy.random.random( self.mesh.coords[0].data.shape) '
        sMap['erf('] = "scipy.special.erf("
        
        sMap['meshx']   = 'self.mesh.coords[0].data'
        sMap['meshy']   = 'self.mesh.coords[1].data'
        sMap['meshz']   = 'self.mesh.coords[2].data'

        sMap['meshi']   = 'self.mesh.indices[0].data'
        sMap['meshj']   = 'self.mesh.indices[1].data'
        sMap['meshk']   = 'self.mesh.indices[2].data'

        
        self.sMap = sMap
        
    def euler(self,time,dt):

        #	Initialize some intermediate arrays
        ncons = self.nPDE
        shape = tuple(self.mesh.shape)

        # Get primative flow variables
        time_i = time
        FLUX = self.updateFlux()
        for U in self.conserved:
            self.variables[U].data += FLUX[U]*dt  

        time = time_i + dt
        self.time = time
        self.updateVars()
        return time

    def rk4_step(self,dt,ii,PHI):

        Ark = [0.0]*5
        Ark[0] = 0.0;
        Ark[1] = -6234157559845./12983515589748.;
        Ark[2] = -6194124222391./4410992767914.;
        Ark[3] = -31623096876824./15682348800105.;
        Ark[4] = -12251185447671./11596622555746.;

        Brk = [0.0]*5
        Brk[0] = 494393426753./4806282396855.;
        Brk[1] = 4047970641027./5463924506627.;
        Brk[2] = 9795748752853./13190207949281.;
        Brk[3] = 4009051133189./8539092990294.;
        Brk[4] = 1348533437543./7166442652324.;

        eta = [0.0]*5
        eta[0] = 494393426753./4806282396855.;
        eta[1] = 4702696611523./9636871101405.;
        eta[2] = 3614488396635./5249666457482.;
        eta[3] = 9766892798963./10823461281321.;
        eta[4] = 1.0;

        #	Initialize some intermediate arrays
        ncons = self.nPDE
        shape = tuple(self.mesh.shape)

        tmp1 = {}
        tmp2 = {}

        if not PHI:
            PHI  = {}
            for U in self.conserved: 
                PHI[U]  = numpy.asfortranarray( numpy.zeros(  shape ) )
        
        for U in self.conserved: 
            tmp1[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            tmp2[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
        
        # Get primative flow variables

        #for ii in range(5):
            #    ii
        FLUX = self.updateFlux()
        for U in self.conserved:
            tmp1[U] =  Ark[ii]*PHI[U]
            PHI[U]  =  dt*FLUX[U] + tmp1[U]
            tmp2[U] =  Brk[ii]*PHI[U]
            self.variables[U].data += tmp2[U]
        self.updateVars()        
        return PHI

    

    def get_tMap(self):
        tMap = {}
        
        # Simple find/replace mappings
        tMap['ddx(#arg#)' ] = r'\frac{\partial #arg#}{\partial x}'
        tMap['ddy(#arg#)' ] = r'\frac{\partial #arg#}{\partial y}'
        tMap['ddz(#arg#)' ] = r'\frac{\partial #arg#}{\partial z}'
        tMap['ddt(#arg#)' ] = r'\frac{\partial #arg#}{\partial t}'

        tMap['meshx'] = 'x'
        tMap['meshy'] = 'y'
        tMap['meshz'] = 'z'

        tMap[':pi:'] = r'\pi '

        #tMap['exp(#arg#)'] = r'\exp( #arg# )'
        #tMap['sqrt(#arg#)'] = r'\sqrt( #arg# )'
        tMap['exp('] = r'\exp('
        tMap['sqrt('] = r'\sqrt('

        tMap['*' ] = ''
        tMap['='] = r'&='
        

        return tMap


    def setupLatex(self):

        # Latex object only available on master
        if self.PyMPI.master:
            self.latex = pyrandaTex(self)
        else:
            self.latex = None
            


def pyrandaRestart(rootname,suffix=None,comm=None):
    from numpy import array,int32
    import numpy as npy
    npy.set_printoptions(threshold=npy.inf)
    """
    Non-member function; return a valid pyrandaSim object
    with data from restart
    """
    # Use suffix else, use largest file
    if not suffix:
        searchDumps = os.path.join(rootname, "restart*" )
        dumps = sorted(glob.glob( searchDumps ))
        dump = dumps[-1]
    else:
        dump = os.path.join(rootname, "restart_%s" % suffix )

    if not os.path.isdir( dump ):
        self.iprint("Error: Cant read restart file %s" % dump)
        return None

    # Get serial data
    fid = open(os.path.join(dump,'serial.dat'))
    dstr = fid.readline()

    # Load local data into current state
    serial_data = eval( dstr )

    message = ''
    

    # Unpack any local vars that have been serialized to pickel
    my_local_vars = {}
    if "local_vars" in serial_data:
        for vv in serial_data["local_vars"]:
            exec("%s = %s" % (vv,serial_data["local_vars"][vv]) )
            message += "    variable %s added \n" % vv
        my_local_vars.update( serial_data['local_vars'] )

        for na in serial_data['local_vars']['numpyArrays']:
            filename = "%s%s" % (serial_data['local_vars']['numpyArrays'][na],'.npy')
            data = numpy.load(filename)
            my_local_vars.update( {na:data} )
            exec("%s = numpy.load('%s')" % (na,filename))
            message += "    numpy array %s added \n" % na
            
        serial_data.pop("local_vars",None)

    # Add user supplied comm to mesh dictionary
    if comm:
        serial_data["mesh"]["comm"] = comm

    if "comm" in serial_data["mesh"]:
        if (type(serial_data["mesh"]["comm"]) != type(MPI.COMM_WORLD) ):
            print("Error: valid comm is not present")
            exit()
        
        
    pysim = pyrandaSim(rootname,serial_data['mesh'])

    # Load packages
    for pack in serial_data['packages']:            
        ipack = pack.split('.')[1] 
        exec("from pyranda import %s" % ipack ) # Edited by D. Lavacot 04/08/2022: change to "from pyranda import" 
        pk = eval("%s(pysim)" % ipack ) 
        #exec("import %s" % ipack )
        #pk = eval("%s.%s(pysim)" % (ipack,ipack) )
        pysim.addPackage( pk )
        message += "    package %s added \n" % ipack

    # EOM and IC's
    # Added by D. Lavacot 04/08/2022: Retrieve user-added functions through sMap
    pysim.sMap = serial_data['sMap']  
    pysim.EOM( serial_data['EOM'] )
    pysim.setIC( serial_data['ICs'], addOnly=True )

    # Simulation state data
    pysim.time = serial_data['time']
    pysim.deltat = serial_data['deltat']
    pysim.cycle  = serial_data['cycle']  
    
    # Loop over processor dumps of restart data:
    # Loop through variables
    procs = serial_data['decomp']
    procMap = serial_data['procMap']

    # Vizualization time history
    pysim.vizDumpHistory = serial_data['vizDumpHistory']
    
    readChunk(pysim,procs,procMap,dump,serial_data)

    # Recompute the mesh metrics
    pysim.mesh.makeMesh(xr=pysim.x().data,
                        yr=pysim.y().data,
                        zr=pysim.z().data)

    # Re-init packages when there is grid dependence
    for pack in pysim.packages:
        pysim.packages[pack].__init__(pysim)
    
    # Print restart summary
    pysim.iprint("Successful read of restart file: %s " % dump )
    pysim.iprint(message)

    # Remove the array meta data
    my_local_vars.pop("numpyArrays",None)
    
    return [ pysim, my_local_vars ]

        
        
def readChunk(pysim,procs,procMap,dump,serial_data):
    """
    Same as readData but only reads in global range of data given by
    irange.
    """

    Rx = [0]*2
    Ry = [0]*2
    Rz = [0]*2

    # This procs extents
    Rx[0] = pysim.PyMPI.chunk_3d_lo[0]    #irange[0]
    Rx[1] = pysim.PyMPI.chunk_3d_hi[0]+1  #irange[1]
    Ry[0] = pysim.PyMPI.chunk_3d_lo[1]    #irange[2]
    Ry[1] = pysim.PyMPI.chunk_3d_hi[1]+1  #irange[3]
    Rz[0] = pysim.PyMPI.chunk_3d_lo[2]    #irange[4]
    Rz[1] = pysim.PyMPI.chunk_3d_hi[2]+1  #irange[5]

    # Restart deomain info
    nprocs = procs[0]*procs[1]*procs[2]
    nprocs = int(nprocs)
    ax = int(pysim.nx / procs[0])
    ay = int(pysim.ny / procs[1])
    az = int(pysim.nz / procs[2])
    nshape = (ax,ay,az,len(serial_data['vars']))
    
    for iproc in range(nprocs):

        g1 = procMap['%s-g1' % iproc] 
        gn = procMap['%s-gn' % iproc] + 1

        # Shift left point if node data
        iff = 0;jff = 0;kff = 0;

        c1 = (Rx[1] in range(g1[0],gn[0]) )
        c2 = (Rx[0] in range(g1[0],gn[0]) )
        c3 = ( (g1[0] and gn[0]) in range(Rx[0],Rx[1]+1) )
        CX = c1 or c2 or c3

        c1 = (Ry[1] in range(g1[1],gn[1]) )
        c2 = (Ry[0] in range(g1[1],gn[1]) )
        c3 = ( (g1[1] and gn[1]) in range(Ry[0],Ry[1]+1) )
        CY = c1 or c2 or c3

        c1 = (Rz[1] in range(g1[2],gn[2]) )
        c2 = (Rz[0] in range(g1[2],gn[2]) )
        c3 = ( (g1[2] and gn[2]) in range(Rz[0],Rz[1]+1) )
        CZ = c1 or c2 or c3

        if ( CX and CY and CZ ):

            Li1 = numpy.max( (0 , Rx[0] - g1[0] ) ) + iff
            Lif = numpy.min( (Rx[1] , gn[0] ) ) - g1[0] + iff
            Ki1 = numpy.max( (Rx[0] , g1[0]) ) - Rx[0]
            Kif = Ki1 + (Lif-Li1)

            Lj1 = numpy.max( (0 , Ry[0] - g1[1] ) ) + jff
            Ljf = numpy.min( (Ry[1] , gn[1] ) ) - g1[1] + jff
            Kj1 = numpy.max( (Ry[0] , g1[1]) ) - Ry[0]
            Kjf = Kj1 + (Ljf-Lj1)

            Lk1 = numpy.max( (0 , Rz[0] - g1[2] ) ) + kff
            Lkf = numpy.min( (Rz[1] , gn[2] ) ) - g1[2] + kff
            Kk1 = numpy.max( (Rz[0] , g1[2]) ) - Rz[0]
            Kkf = Kk1 + (Lkf-Lk1)

            #pdata = self.readDataProc(time,iproc,variable)
            fid = open(os.path.join(dump,"proc-%s.bin" % str(iproc).zfill(5)))
            DATA = numpy.reshape(numpy.fromfile( fid ),nshape,order='C')
                       
            #for ii in range(len(variable)):
            for var in serial_data['vars']:
                #vdata[ii][Ki1:Kif,Kj1:Kjf,Kk1:Kkf] = pdata[ii][Li1:Lif,Lj1:Ljf,Lk1:Lkf]
                if isinstance(pysim.variables[var].data,numpy.ndarray):
                    pysim.variables[var].data[Ki1:Kif,
                                              Kj1:Kjf,
                                              Kk1:Kkf] = DATA[Li1:Lif,
                                                              Lj1:Ljf,
                                                              Lk1:Lkf,
                                                              serial_data['vars'][var]]
                else:
                    pysim.variables[var].data = DATA[0,0,0,
                                                     serial_data['vars'][var]]
