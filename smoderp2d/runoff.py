# @package smoderp2d.runoff loop of the modul
#
#  The computing area is determined  as well as the boundary cells.
#
#  \e vypocet probiha v zadanem casovem kroku, pripade je cas kracen podle \b "Couranotva kriteria":
#    - vystupy jsou rozdelieny do \b zakladnich a \b doplnkovych, podle zvoleneh typu vypoctu
#    - \b zakladni
#        - @return \b h0 maximalni vyska haladiny plosneho odtoku
#


#!/usr/bin/python
# -*- coding: latin-1 -*-
# Surface and subsurface rutine
# Created by Petr Kavka, FCE, CTU Prague, 2015

__author__ = "edlman"
__date__ = "$29.12.2015 18:31:25$"

# INITIAL SETTINGS:
#
# importing system moduls
# import math
import numpy as np
import time
import os
import sys
#from   smoderp2d.classes_main_arrays import *
#from   smoderp2d.tools.resolve_partial_computing import *

# importing classes


from smoderp2d.arrs.General import Globals
from smoderp2d.arrs.Solve import ImplicitSolver
from smoderp2d.iter_crit import IterCrit
from smoderp2d.providers.logger import Logger


# def init_classes():

# return LS, courant, delta_t


class Runoff():

    def __init__(self, provider):
        gl = Globals()

        self.iter_crit = IterCrit()
        self.iter_crit.set_time_step(Globals.maxdt)

        self.LS = ImplicitSolver(self.iter_crit)
        self.provider = provider

    def run(self):
        import smoderp2d.flow_algorithm.D8 as D8_

        gl = Globals()

        self.LS.solveStep(self.iter_crit)
                
        self.provider.progress(self.iter_crit.dt, 0, self.LS.total_time)
        self.LS.total_time += self.iter_crit.dt
        self.iter_crit.check_time_step()

        ok = 1
        while (self.LS.total_time <= gl.end_time):
            self.iter_crit.reset()
            if ok == 1: 
                self.LS.hold = self.LS.hnew.copy()
            ok = self.LS.solveStep(self.iter_crit)
            if ok == 0 :
                self.provider.progress(self.iter_crit.dt, self.iter_crit.iter_, self.LS.total_time)
                self.iter_crit.check_time_step02()
                self.provider.message('repeating time step')
                self.provider.message("-" * 40)
            if ok == 1:
                self.provider.progress(self.iter_crit.dt, self.iter_crit.iter_, self.LS.total_time)
                self.LS.total_time += self.iter_crit.dt
                self.iter_crit.check_time_step()
        
        self.provider.message('total time  [s]' + str(time.time()-self.provider.startTime))
        #Logger.debug('hcrit natvrdo')
        
        return 0


















