class Error(Exception):

    """Base class for exceptions in this module."""
    pass


class MaxIterationExceeded(Error):

    """Exception raised if number of iteration exceed max iteration criterion
    Attributes:
        maxIter  -- max n of iteration
    """

    def __init__(self, mi, t):
        self.msg = 'Maximum of iterations (maxIter = ' + str(
            mi) + ') was exceeded of at time [s]: ' + str(t) + '.'

    def __str__(self):
        return repr(self.msg)
    
    