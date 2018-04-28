import math


def sheet_flow(a,b,h):

    return math.pow(h, b) * a


def sheet_flowb_(a,b,h):

    return math.pow(h, b-1) * a
