# @package smoderp2d.courant defines iteration criterion

from smoderp2d.exceptions import MaxIterationExceeded
from smoderp2d.core.general import Globals


class IterCrit():

    # constructor
    #
    def __init__(self):
        self.iter_ = 0
        self.max_delta_t = Globals.maxdt
        self.dt = Globals.maxdt
        self.min_delta_t = 0.1
        self.max_iter = 5
        self.crit_iter = 3
        self.min_iter = 1
        self.pre_rill_count = 0
        self.multipiler = 0.8

    # Reset the iteration count
    #
    def reset(self):
        self.pre_rill_count = 0
        self.iter_ = 0

    # Increase interataion count
    #
    def iter_up(self):
        self.iter_ += 1

    # Check max iteration reached in this case the time is repeted
    #
    def max_iter_check(self, total_time):
        if (self.iter_ >= self.max_iter):
            return True
            #raise MaxIterationExceeded(self.crit_iter_, total_time)

    # Adjust time step
    #
    def check_time_step(self):
        if self.iter_ >= self.crit_iter:
            self.dt = max(0.1, self.dt * 0.7)
        if self.iter_ <= self.min_iter:
            self.dt = min(self.max_delta_t, self.dt * (1/self.multipiler))

    # Adjust time step
    #
    def check_time_step02(self):
        self.dt = max(0.1, self.dt * 0.5)

    # Tranck rill emergency
    #
    def rill_check(self, rc):  # rc - rill count
        if rc > self.pre_rill_count:
            self.dt = max(0.1, self.dt * 0.3)
            self.pre_rill_count = rc
