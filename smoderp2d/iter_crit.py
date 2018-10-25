# @package smoderp2d.courant defines iteration criterion


import math
import smoderp2d.constants as constants
from smoderp2d.tools.tools import comp_type
from smoderp2d.tools.tools import get_argv
import smoderp2d.io_functions.prt as prt
from smoderp2d.exceptions import MaxIterationExceeded

from smoderp2d.core.General import Globals as Gl


class IterCrit():

    # constructor
    #
    def __init__(self):
        self.iter_ = 0
        self.crit_iter_ = 7
        self.max_delta_t = Gl.maxdt
        self.dt = Gl.maxdt
        self.min_delta_t = 0.1
        self.max_iter = 7
        self.min_iter = 4
        self.pre_rill_count = 0

    # Reset the iteration count
    #
    def reset(self):
        self.pre_rill_count = 0
        self.iter_ = 0
    
    # Increase interataion count
    # 
    def iter_up(self):
        self.iter_ += 1
    
    # Check max iteration reached
    def crit_iter_check(self,total_time):
        if (self.iter_ >= self.crit_iter_):
            return True
            #raise MaxIterationExceeded(self.crit_iter_, total_time)
    
    # Adjust time step
    #
    def check_time_step(self):
        if self.iter_ >= self.max_iter:
            self.dt = max(0.1, self.dt * 0.7)
        if self.iter_ <= self.min_iter:
            self.dt = min(self.max_delta_t, self.dt * 1.3)
            
    def check_time_step02(self):
        self.dt = max(0.1, self.dt * 0.3)
        #if self.iter_ >= self.max_iter:
        #if self.iter_ <= self.min_iter:
            #self.dt = min(self.max_delta_t, self.dt * 1.3)
    
    # Tranck rill emergency
    # 
    def rill_check(self, rc):  # rc - rill count
        if rc > self.pre_rill_count:
            self.dt = max(0.1, self.dt * 0.3)
            self.pre_rill_count = rc
