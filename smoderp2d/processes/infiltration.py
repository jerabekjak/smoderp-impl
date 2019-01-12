from smoderp2d.core.general import Globals

import numpy as np
import math


# combinatIndex muze byt tady jako globalni
# primenna, main loop bude pro infiltraci volat
# stejnou funkci, ale az taky bude jina globalni
# promenna nastavena. mene ifu v main loop
combinatIndex = []


class Infiltration (object):

    def __init__(self):
        """ Copy combinat index from Globals """
        self.combinatIndex = Globals.get_combinatIndex()

    def potential_infiltration(self, total_time, dt):
        """ Calculate potential infiltration (pi) based on Philips equation

        The pi is stored in combinatIndex

        :param float total_time : total time of the computation
        :param float dt         : length od time step
        """
        for iii in self.combinatIndex:
            index = iii[0]
            k = iii[1]
            s = iii[2]
            iii[3] = self.phlilip(
                k, s, dt, total_time, Globals.get_NoDataValue())

    def philip_infiltration(self, i, j):
        """ Find height of infiltrated water from the combinatIndex. """
        soil = Globals.get_mat_inf_index(i, j)
        for z in self.combinatIndex:
            if soil == z[0]:
                return z[3]

    def phlilip(self, k, s, deltaT, totalT, NoDataValue):
        """ Calculates infiltratin based on the Philips equation """
        if k and s == NoDataValue:
            infiltration = NoDataValue
        else:
            infiltration = (0.5 * s / math.sqrt(totalT + deltaT) + k) * deltaT

        return infiltration
