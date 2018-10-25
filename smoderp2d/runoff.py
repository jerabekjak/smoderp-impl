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
# INITIAL SETTINGS:
#
# importing system moduls
# import math
import numpy as np
import time
import os
import sys

from smoderp2d.core.General import Globals


class Runoff():

    def __init__(self, provider):
        from smoderp2d.core.solve import ImplicitSolver
        from smoderp2d.iter_crit import IterCrit

        gl = Globals()

        # init iteration and time step controll
        self.iter_crit = IterCrit()

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
            if ok == 0:
                self.provider.progress(
                    self.iter_crit.dt,
                    self.iter_crit.iter_,
                    self.LS.total_time)
                self.iter_crit.check_time_step02()
                self.provider.message('repeating time step')
                self.provider.message("-" * 40)
            if ok == 1:
                self.provider.progress(
                    self.iter_crit.dt,
                    self.iter_crit.iter_,
                    self.LS.total_time)
                self.LS.total_time += self.iter_crit.dt
                self.iter_crit.check_time_step()

        self.provider.message(
            'total time  [s]' + str(time.time() - self.provider.startTime))

        return 0
