# @package smoderp2d.courant defines iteration criterion


import math
import smoderp2d.constants as constants
from smoderp2d.tools.tools import comp_type
from smoderp2d.tools.tools import get_argv
import smoderp2d.io_functions.prt as prt
from smoderp2d.exceptions import MaxIterationExceeded

from smoderp2d.arrs.General import Globals as Gl


class IterCrit():

    # constructor
    #
    def __init__(self):
        self.iter_ = 0
        self.crit_iter_ = 50
        self.max_delta_t = Gl.maxdt
        self.min_delta_t = 0.1
        self.max_iter = 10
        self.min_iter = 4
        self.pre_rill_count = 0

    # Store the original guess time step
    #
    def set_time_step(self, dt):
        self.dt = dt

    def reset(self):
        self.pre_rill_count = 0
        self.iter_ = 0
        
    def iter_up(self):
        self.iter_ += 1
        
    def crit_iter_check(self,total_time):
        if (self.iter_ > self.crit_iter_):
            raise MaxIterationExceeded(self.crit_iter_, total_time)
        
    def check_time(self):
        if self.iter_ >= self.max_iter:
            self.dt = max(0.1, self.dt * 0.7)
        if self.iter_ <= self.min_iter:
            self.dt = min(self.max_delta_t, self.dt * 1.3)

    def rill_check(self, rc):  # rc - rill count
        if rc > self.pre_rill_count:
            self.dt = max(0.1, self.dt * 0.3)
            self.pre_rill_count = rc
