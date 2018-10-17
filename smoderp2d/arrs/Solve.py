from smoderp2d.arrs.General import Globals
from smoderp2d.tools.tools import make_ASC_raster
from smoderp2d.tools.tools import make_sur_raster
from smoderp2d.arrs.Surface import Surface
import smoderp2d.processes.rainfall as rain_f
import smoderp2d.processes.infiltration as infilt
import smoderp2d.flow_algorithm.D8 as D8_
from smoderp2d.io_functions import hydrographs as wf

from smoderp2d.providers.logger import Logger
from scipy.sparse import csr_matrix


from ctypes import CDLL, POINTER, c_int, c_float

import numpy as np
import os
import sys
import time


def init_getIJel():
    """Documentation for a function init_getIJel()

    popis:
    funkce ma vytvorit pole a vektory na zjistovani
    odpovidajiciho radku a sloupce rasteru podle pozive v lin soustave
    a odpovidajiciho elementu v soustave podle sloupce a rarku v rasteru

    return nEl      - pocet elementu v soustave
    return getEl   - i j pole, obsahuje pozici v soustave
    return getElIN - i j pole, pozive vtekajicich elementu
    return getIJ   - vektor, vrati odpovidajici i j pro dany element v soustave
    """

    gl = Globals()
    inflows = D8_.new_inflows(gl.mat_fd)

    r = gl.get_rows()
    c = gl.get_cols()
    rr = gl.get_rrows()
    rc = gl.get_rcols()
    br = gl.get_bor_rows()
    bc = gl.get_bor_cols()

    getEl = np.zeros([r, c], int)
    getIJ = []
    getElIN = []

    nEl = int(-1)

    # prvni cyklus priradi
    # k i j bunky jeji poradi ve vektoru
    # a k elementu vektoru i j bunky
    for i in rr:
        for j in rc[i]:
            nEl += int(1)
            getIJ.append([i, j])
            getEl[i][j] = nEl

    for i in br:
        for j in bc[i]:
            nEl += int(1)
            getIJ.append([i, j])
            getEl[i][j] = nEl

    # druhy cyklus priradi
    # k elementu list elementu v inflows
    for i in rr:
        for j in rc[i]:
            n = len(inflows[i][j])
            tmp = [-int(99)] * 8

            for z in range(n):
                ax = inflows[i][j][z][0]
                bx = inflows[i][j][z][1]
                iax = i + ax
                jbx = j + bx
                tmp[z] = getEl[iax][jbx]

            getElIN.append(tmp)

    for i in br:
        for j in bc[i]:
            n = len(inflows[i][j])
            tmp = [-int(99)] * 8

            for z in range(n):
                ax = inflows[i][j][z][0]
                bx = inflows[i][j][z][1]
                iax = i + ax
                jbx = j + bx
                tmp[z] = getEl[iax][jbx]

            getElIN.append(tmp)

    # toto je pokazde stejne pro kinematickou vlnu
    indptr = [0]
    # toto je pokazde stejne pro kinematickou vlnu
    indices = []

    for iel in range(nEl):
        indices.append(iel)
        for inel in getElIN[iel]:
            if inel >= 0:
                indices.append(inel)
        indptr.append(len(indices))

    # convert co ctype F ordered array
    getIJ = np.array(getIJ, dtype=c_int, order='F')
    getElIN = np.array(getElIN, dtype=c_int, order='F')

    return nEl, getEl, getElIN, getIJ, indices, indptr


class ImplicitSolver:

    def __init__(self, iter_crit):

        gl = Globals()
        r = gl.get_rows()
        c = gl.get_cols()

        # interval srazky
        self.tz = 0
        # aktualni cas programu
        self.total_time = 0

        self.nEl, self.getEl, self.getElIN, self.getIJ, self.indices, self.indptr = init_getIJel()

        #self.A = np.zeros([self.nEl, self.nEl], float)
        self.b = np.zeros([self.nEl], dtype=c_float)
        self.hnew = np.ones([self.nEl], dtype=c_float)
        self.hold = np.zeros([self.nEl], dtype=c_float)

        # variable counts rills
        self.rill_count = 0

        self.sur = Surface()

        # opens files for storing hydrographs
        if gl.points and gl.points != "#":
            self.hydrographs = wf.Hydrographs()
            arcgis = gl.arcgis
            if not arcgis:
                with open(os.path.join(gl.outdir, 'points.txt'), 'w') as fd:
                    for i in range(len(gl.array_points)):
                        fd.write('{} {} {} {}'.format(
                            gl.array_points[i][0], gl.array_points[i][3],
                            gl.array_points[i][4], os.linesep
                        ))
        else:
            self.hydrographs = wf.HydrographsPass()

        for i in Globals.rr:
            for j in Globals.rc[i]:
                self.hydrographs.write_hydrographs_record(
                    i,
                    j,
                    iter_crit,
                    self,
                    first=True
                )
                
        gl.mat_inf_index = np.asarray(gl.mat_inf_index, dtype=c_int, order='F')
        gl.mat_aa = np.asarray(gl.mat_aa, dtype=c_float, order='F')
        gl.mat_b = np.asarray(gl.mat_b, dtype=c_float, order='F')
        gl.mat_hcrit = np.asarray(gl.mat_hcrit, dtype=c_float, order='F')
        gl.mat_n = np.asarray(gl.mat_n, dtype=c_float, order='F')
        gl.mat_slope = np.asarray(gl.mat_slope, dtype=c_float, order='F')
        gl.mat_efect_vrst = np.asarray(
            gl.mat_efect_vrst, dtype=c_float, order='F')

        self.fortran = CDLL('smoderp2d/f/fill_a_mat.so')

        # TODO to oby se asi nemelo delat pokazde
        self.fortran.fill_a_mat.argtypes = [POINTER(c_int),    # nel
                                       POINTER(c_int),    # sizes
                                       POINTER(c_int),    # getIJ
                                       POINTER(c_int),    # getElIN
                                       POINTER(c_float),  # data
                                       POINTER(c_float),  # hnew
                                       POINTER(c_float),  # hold
                                       POINTER(c_float),  # mat_aa
                                       POINTER(c_float),  # mat_b
                                       POINTER(c_float),  # mat_hcrit
                                       POINTER(c_float),  # mat_effect_vrst
                                       POINTER(c_int),    # mat_inf_index
                                       POINTER(c_float),  # mat_n
                                       POINTER(c_float),  # mat_slope
                                       POINTER(c_float),  # combinatIndex
                                       POINTER(c_float),  # dx
                                       POINTER(c_float)]  # dt

    def fillAmat(self, dt):

        gl = Globals()
        from smoderp2d.processes.surface import sheet_flowb_
        # if rills are calculated
        if gl.isRill:
            from smoderp2d.processes.rill import rill
        else:
            from smoderp2d.processes.rill import rill_pass
            rill = rill_pass

        # potential precipitation
        PS, self.tz = rain_f.timestepRainfall(self.total_time, dt, self.tz)

        # prepares infiltration table
        for iii in gl.get_combinatIndex():
            index = iii[0]
            k = iii[1]
            s = iii[2]
            iii[3] = infilt.phlilip(
                k, s, dt, self.total_time, gl.get_NoDataValue())

        # creates empty list for data
        data = np.ones((len(self.indices)), dtype=c_float)

        n_data = len(data)
        n_mat = gl.r
        m_mat = gl.c
        n_combinatIndex = len(gl.combinatIndex)
        m_combinatIndex = len(gl.combinatIndex[0])

        # combinatIndex must be numpy array
        combinatIndex = np.array(gl.combinatIndex, dtype=c_float, order='F')

        # do funkce to musi jit numpy array
        sizes = np.array(
            [n_data, n_mat, m_mat, n_combinatIndex, m_combinatIndex], dtype=c_int)
     

        self.fortran.fill_a_mat(c_int(self.nEl),
                           sizes.ctypes.data_as(POINTER(c_int)),
                           self.getIJ.ctypes.data_as(POINTER(c_int)),
                           self.getElIN.ctypes.data_as(POINTER(c_int)),
                           data.ctypes.data_as(POINTER(c_float)),
                           self.hnew.ctypes.data_as(POINTER(c_float)),
                           self.hold.ctypes.data_as(POINTER(c_float)),
                           gl.mat_aa.ctypes.data_as(POINTER(c_float)),
                           gl.mat_b.ctypes.data_as(POINTER(c_float)),
                           gl.mat_hcrit.ctypes.data_as(POINTER(c_float)),
                           gl.mat_efect_vrst.ctypes.data_as(POINTER(c_float)),
                           gl.mat_inf_index.ctypes.data_as(POINTER(c_int)),
                           gl.mat_n.ctypes.data_as(POINTER(c_float)),
                           gl.mat_slope.ctypes.data_as(POINTER(c_float)),
                           combinatIndex.ctypes.data_as(POINTER(c_float)),
                           c_float(gl.dx),
                           c_float(dt))

        for iel in range(self.nEl):
            i = self.getIJ[iel][0]
            j = self.getIJ[iel][1]

            # infiltration
            inf = infilt.philip_infiltration(
                gl.get_mat_inf_index(i, j), gl.get_combinatIndex())
            if inf >= self.hold[iel]:
                inf = self.hold[iel]

            self.b[iel] = self.hold[iel] / dt + PS / dt - inf / dt

        t3 = time.time()

        self.A = csr_matrix((data, self.indices, self.indptr),
                            shape=(self.nEl, self.nEl), dtype=float)
        t4 = time.time()

    def solveStep(self, iter_crit):
        from scipy.sparse.linalg import spsolve

        hewp = self.hnew.copy()

        err = 1

        while (err > 0.000001):

            iter_crit.iter_up()
            self.fillAmat(iter_crit.dt)
            hewp = self.hnew.copy()
            self.hnew = spsolve(self.A, self.b)

            if (iter_crit.crit_iter_check(self.total_time)):
                return 0

            err = np.sum(hewp - self.hnew)**2.0

        # write hydrograph record
        for i in Globals.rr:
            for j in Globals.rc[i]:
                self.hydrographs.write_hydrographs_record(
                    i,
                    j,
                    iter_crit,
                    self
                )

        return 1
