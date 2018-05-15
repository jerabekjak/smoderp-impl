# @package smoderp2d.courant defines iteration criterion


import math
import smoderp2d.constants as constants
from smoderp2d.tools.tools import comp_type
from smoderp2d.tools.tools import get_argv
import smoderp2d.io_functions.prt as prt


from smoderp2d.arrs.General import Globals as Gl


class IterCrit():

    # constructor
    #
    def __init__(self):
        self.max_delta_t = Gl.maxdt
        self.min_delta_t = 0.1
        self.max_iter = 10
        self.min_iter = 4
        self.pre_rill_count = 0

    # Store the original guess time step
    #
    def set_time_step(self, dt):
        self.dt = dt

    # Resets the self.cour_most and self.cour_speed after each time stop computation is successfully completed
    #
    def reset(self):
        self.pre_rill_count = 0

    def check_time(self, iter_):
        if iter_ >= self.max_iter:
            self.dt = max(0.1, self.dt * 0.7)
        if iter_ <= self.min_iter:
            self.dt = min(self.max_delta_t, self.dt * 1.3)

    def rill_check(self, rc):  # rc - rill count
        if rc > self.pre_rill_count:
            self.dt = max(0.1, self.dt * 0.3)
            self.pre_rill_count = rc
