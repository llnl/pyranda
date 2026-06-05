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
from .pyrandaUtils import fortran3d, findVar
import numpy  # noqa: F401 (used implicity)


class pyrandaEq:
    def __init__(self, eqstr, sMap, pympi):
        """
        Read in eqstr and extract the LHS variable
        and the kind of equation (PDE or ALG)
        """
        self.eqstr = eqstr
        self.kind = "ALG"
        self.active = True
        self.rank = 1
        self.sRHS = ""
        if "=" in eqstr:
            self.LHS = findVar(
                eqstr.split("=")[0], "scalar", unique=False
            )  # Return values are assumed to all be unique.  Unique to false ensure order is preserved
            self.rank = len(self.LHS)
        else:
            self.LHS = None  # No return functions

        # Make a lambda for this equation
        if self.LHS:
            Srhs = fortran3d(self.eqstr.split("=")[1], sMap)
        else:
            Srhs = fortran3d(self.eqstr, sMap)
        self.sRHS = Srhs
        self.RHS = eval("lambda self: " + Srhs)

        # Check to see if this is conserved PDE
        if "ddt(" in eqstr:
            self.kind = "PDE"

        # Check to see if this is an elliptic equation
        # if ( 'Delta(' in eqstr ):
        #    self.kind = 'ELP'
        #    self.solver = None
