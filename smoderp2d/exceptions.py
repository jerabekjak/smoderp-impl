class SmoderpError(Exception):
    pass


class MaxIterationExceeded(SmoderpError):
    """Number of iteration exceed max iteration criterion.
    """

    def __init__(self, mi, t):
        self.msg = 'Maximum of iterations (maxIter = {}) was exceeded of at time [s]: {}'.format(
            mi, t
        )

    def __str__(self):
        return repr(self.msg)



class NegativeVolumeInRIll(SmoderpError):
    """Number of iteration exceed max iteration criterion.
    """
    def __init__(self, V):
        self.msg = 'Negative volume goes to rill (V_to_rill = {}).'.format(
            V
        )

    def __str__(self):
        return repr(self.msg)
