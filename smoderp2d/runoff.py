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



from smoderp2d.courant import Courant


#def init_classes():


    #return LS, courant, delta_t


class Runoff():
    
    def __init__(self,provider):
        gl = Globals()
        self.delta_t = Globals.maxdt

        self.courant = Courant()
        self.courant.set_time_step(self.delta_t)

        self.LS = ImplicitSolver()
        


    def run(self):
        import smoderp2d.flow_algorithm.D8 as D8_

        # taky se vyresi vztypbi soubory nacteni dat
        # vse se hodi do ogjektu Globals as Gl
        gl = Globals()
        #LS, courant, delta_t = init_classes()
        
        t1 = time.time()
        print self.LS.total_time+self.courant.dt,
        self.LS.solveStep(self.courant)
        
        self.courant.reset()
        
        
        
        while (self.LS.total_time <= gl.end_time):
            #print LS.total_time/60, gl.end_time
            self.LS.hold = self.LS.hnew.copy()
            print self.LS.total_time+self.courant.dt, 
            self.LS.solveStep(self.courant)
            
        print 'cas vypoctu', time.time() - t1
        return 0




