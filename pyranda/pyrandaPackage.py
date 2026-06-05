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


class pyrandaPackage:
    """
    Case physics package module for adding new physics packages to pyranda
    """

    def __init__(self, name, pyranda):

        self.name = name
        self.pyranda = pyranda

    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        self.sMap = sMap
