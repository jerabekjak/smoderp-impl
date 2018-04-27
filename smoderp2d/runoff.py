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


from smoderp2d.main_classes.General import Globals
from smoderp2d.main_classes.Solve import ImplicitSolver

from smoderp2d.courant import Courant


def init_classes():

    gl = Globals()
    delta_t = gl.get_max_dt()

    courant = Courant()
    courant.set_time_step(delta_t)

    LS = ImplicitSolver()

    return LS, courant, delta_t


class Runoff():

    def run(self):
        import smoderp2d.flow_algorithm.D8 as D8_

        # taky se vyresi vztypbi soubory nacteni dat
        # vse se hodi do ogjektu Globals as Gl
        gl = Globals()
        LS, courant, delta_t = init_classes()
        
        print gl.end_time
        t1 = time.time()
        print LS.total_time+delta_t,
        LS.solveStep(delta_t)
        
        
        
        while (LS.total_time <= gl.end_time):
            #print LS.total_time/60, gl.end_time
            LS.hold = LS.hnew.copy()
            print LS.total_time+delta_t, 
            LS.solveStep(delta_t)

        return 0




