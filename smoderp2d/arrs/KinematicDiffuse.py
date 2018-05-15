

from smoderp2d.arrs.General import Globals as Gl
from smoderp2d.arrs.Flow import *


import smoderp2d.io_functions.prt as prt
import smoderp2d.flow_algorithm.flow_direction as flow_direction


class Kinematic(Mfda if Gl.mfda == True else D8):

    def __init__(self):
        prt.message("\tKinematic approach")
        super(Kinematic, self).__init__()

    def new_inflows(self):
        pass

    def update_H(self):
        pass


class Diffuse(Mfda if Gl.mfda == True else D8):

    def __init__(self):
        prt.message("\tDiffuse approach")
        if (Globals.r == None or Globals.r == None):
            exit("Global variables are not assigned")
        r = Gl.r
        c = Gl.c

        self.H = np.zeros([r, c], float)

    def new_inflows(self):
        fd = flow_direction.flow_direction(
            self.H, self.rr, self.rc, self.br, self.bc, self.pixel_area)
        self.update_inflows(fd)

    def update_H(self):

        arr = self.arr

        for i in self.rr:
            for j in self.rc[i]:
                arr.H[i][j] = arr.h[i][j] + arr.z[i][j]

        for i in self.rr:
            for j in self.rc[i]:
                arr.slope = self.slope_(i, j)