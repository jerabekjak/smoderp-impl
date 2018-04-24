# @package smoderp2d.src.runoff loop of the modul
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
#from   smoderp2d.src.classes_main_arrays import *
#from   smoderp2d.src.tools.resolve_partial_computing import *

# importing classes


from smoderp2d.src.main_classes.General import Globals
#from smoderp2d.src.main_classes.Vegetation    import Vegetation
#from smoderp2d.src.main_classes.Surface       import Surface
from smoderp2d.src.main_classes.Solve import ImplicitSolver


#from smoderp2d.src.main_classes.Subsurface    import Subsurface
#from smoderp2d.src.main_classes.CumulativeMax import Cumulative
#from smoderp2d.src.time_step                  import TimeStep

#import smoderp2d.src.constants                 as constants
from smoderp2d.src.courant import Courant
#import smoderp2d.src.tools.tools               as tools
#import smoderp2d.src.io_functions.post_proc    as post_proc
#import smoderp2d.src.io_functions.prt          as prt
#import smoderp2d.src.io_functions.progress_bar as progress_bar
from smoderp2d.src.tools.tools import comp_type
#from   smoderp2d.src.tools.times_prt       import TimesPrt
#from smoderp2d.src.tools.tools             import get_argv


def init_classes():

    import time
    gl = Globals()
    delta_t = gl.get_max_dt()

    #isRill, subflow, stream, diffuse, = comp_type()

    courant = Courant()
    courant.set_time_step(delta_t)

    IS = ImplicitSolver()

    # raw_input('...')
    return IS, courant, delta_t


class Runoff():

    def run(self):
        import smoderp2d.src.flow_algorithm.D8 as D8_

        # taky se vyresi vztypbi soubory nacteni dat
        # vse se hodi do ogjektu Globals as Gl

        A, courant, delta_t = init_classes()

        t1 = time.time()
        print A.total_time+delta_t,
        A.solveStep(delta_t)
        #print 'IS.total_time' , IS.total_time

        while (A.total_time/60 <= 20):
            A.hold = A.hnew.copy()
            print A.total_time+delta_t,
            A.solveStep(delta_t)
            #print 'IS.total_time' , IS.total_time

        #print "total time: ", time.time() - t1

        return 0
