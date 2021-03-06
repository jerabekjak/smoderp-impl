#!/usr/bin/python
# -*- coding: latin-1 -*-
# SMODERP 2D
# Created by Jan Zajicek, FCE, CTU Prague, 2012-2013

import numpy as np
import sys
import smoderp2d.io_functions.prt as prt
from smoderp2d.core.general import Globals


# definice erroru  na urovni modulu
#
class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class NonCumulativeRainData(Error):
    """Exception raised bad rainfall record assignment.

    Attributes:
        msg  -- explanation of the error
    """

    def __init__(self):
        self.msg = 'Error: Rainfall record has to be cumulative'

    def __str__(self):
        return repr(self.msg)


def load_precipitation(fh):
    y2 = 0
    try:
        fh = open(fh, "r")
        x = []
        for line in fh.readlines():
            z = line.split()
            if len(z) == 0:
                continue
            elif z[0].find('#') >= 0:
                continue
            else:
                if len(z) == 0:
                    continue
                else:
                    y0 = float(z[0]) * 60.0  # prevod na vteriny
                    y1 = float(z[1]) / 1000.0  # prevod na metry
                    if y1 < y2:
                        raise NonCumulativeRainData()
                    y2 = y1
                    mv = y0, y1
                    x.append(mv)
        fh.close

        # Values ordered by time ascending
        dtype = [('cas', float), ('value', float)]
        val = np.array(x, dtype=dtype)
        x = np.sort(val, order='cas')
        # Test if time time is more than once the same
        state = 0
        k = 1
        itera = len(x)  # iter is needed in main loop
        for k in range(itera):
            if x[k][0] == x[k - 1][0] and itera != 1:
                state = 1
                y = np.delete(x, k, 0)

        if state == 0:
            x = x
        else:
            x = y
        # Amount of rainfall in individual intervals
        if len(x) == 0:
            sr = 0
        else:
            sr = np.zeros([itera, 2], float)
            for i in range(itera):
                if i == 0:
                    sr_int = x[i][1] / x[i][0]
                    sr[i][0] = x[i][0]
                    sr[i][1] = sr_int

                else:
                    sr_int = (x[i][1] - x[i - 1][1]) / (x[i][0] - x[i - 1][0])
                    sr[i][0] = x[i][0]
                    sr[i][1] = sr_int

        return sr, itera

    except IOError:
        prt.message("The file does not exist!")
    except BaseException:
        prt.message("Unexpected error:", sys.exc_info()[0])
        raise


class Rainfall():

    def __init__(self):

        self.tz = 0
        self.tz_save = 0
        self.veg = Globals.get_mat_nan().copy()
        self.veg.fill(int(1))
        self.veg_save = self.veg.copy()
        self.sum_interception = Globals.get_mat_nan().copy()
        self.sum_interception.fill(int(0))
        self.sum_interception_save = self.sum_interception.copy()

    def timestepRainfall(self, total_time, delta_t):
        """Function returns a rainfall amount for current time step
        if two or more rainfall records belongs to one time step
        the function integrates the rainfall amount.

        :param float total_time: total time of the computation
        :param float delta_t   : time step size

        :return float: potential precipitation
        """

        iterace = Globals.itera
        sr = Globals.sr

        z = self.tz
        # skontroluje jestli neni mimo srazkovy zaznam
        if z > (iterace - 1):
            rainfall = 0
        else:
            # skontroluje jestli casovy krok, ktery prave resi, je stale vramci
            # srazkoveho zaznamu z

            if sr[z][0] >= (total_time + delta_t):
                rainfall = sr[z][1] * delta_t
            # kdyz je mimo tak
            else:
                # dopocita zbytek ze zaznamu z, ktery je mezi total_time a
                # total_time + delta_t
                rainfall = sr[z][1] * (sr[z][0] - total_time)
                # skoci do dalsiho zaznamu
                z += 1
                # koukne jestli ten uz neni mimo
                if z > (iterace - 1):
                    rainfall += 0
                else:
                    # pokud je total_time + delta_t stale dal nez konec posunuteho zaznamu
                    # vezme celou delku zaznamu a tuto srazku pricte
                    while (sr[z][0] <= (total_time + delta_t)):
                        rainfall += sr[z][1] * (sr[z][0] - sr[z - 1][0])
                        z += 1
                        if z > (iterace - 1):
                            break
                    # nakonec pricte to co je v poslednim zaznamu kde je total_time + delta_t pred konce zaznamu
                    # nebo pricte nulu pokud uz tam zadny zaznam neni
                    if z > (iterace - 1):
                        rainfall += 0
                    else:
                        rainfall += sr[z][1] * \
                            (total_time + delta_t - sr[z - 1][0])

                self.tz = z

        return rainfall

    def current_rain(self, i, j, potential_rain):
        """ Reduces potetial rainfall by interception """
        ppl = Globals.get_ppl(i, j)
        pi = Globals.get_pi(i, j)
        if self.veg[i][j] != int(5):
            interc = ppl * potential_rain  # interception is konstant
            # jj nemelo by to byt interc = (1-rain_ppl) * rainfallm
            #                             -------------

            self.sum_interception[i][j] += interc  # sum of intercepcion
            NS = potential_rain - interc  # netto rainfallm
            # jj nemela by byt srazka 0 dokun neni naplnena intercepcni zona?
            #

            # if potentional interception is overthrown by intercepcion sum,
            # then the rainfall is effetive
            if self.sum_interception[i][j] >= pi:
                self.veg[i][j] = int(5)
        else:
            NS = potential_rain

        return NS

    def save_rainfall_vars(self):
        """ Store variables in case for iterations. """
        self.veg_save = self.veg.copy()
        self.sum_interception_save = self.sum_interception.copy()
        self.tz_save = self.tz

    def load_rainfall_vars(self):
        """ Restore variables in case for iterations. """
        self.veg = self.veg_save.copy()
        self.sum_interception = self.sum_interception_save.copy()
        self.tz = self.tz_save
