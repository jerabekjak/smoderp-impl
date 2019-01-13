from smoderp2d.core.general import Globals
from smoderp2d.tools.tools import make_ASC_raster
from smoderp2d.tools.tools import make_sur_raster


import smoderp2d.flow_algorithm.D8 as D8_
from smoderp2d.io_functions import hydrographs as wf

from smoderp2d.providers.logger import Logger
from scipy.sparse import csr_matrix

import numpy as np
import os
import sys
import time


def prt_arr(A):
    """ Print csr_matrix into file and on the screen """
    np.savetxt('arr.txt', A.toarray(), fmt='%1.4e', delimiter='\t')
    print ('Printing array...')
    Logger.info(A.toarray())


def init_getIJel():
    """Documentation for a function init_getIJel()

    popis:
    funkce ma vytvorit pole a vektory na zjistovani
    odpovidajiciho radku a sloupce rasteru podle pozive v lin soustave
    a odpovidajiciho elementu v soustave podle sloupce a rarku v rasteru

    return nEl     - pocet elementu v soustave
    return getEl   - i j pole, obsahuje pozici v soustave
    return getElIN - i j pole, pozive vtekajicich elementu
    return getIJ   - vektor, vrati odpovidajici i j pro dany element v soustave
    """

    gl = Globals()
    inflows = D8_.new_inflows(Globals.mat_fd)

    r = Globals.get_rows()
    c = Globals.get_cols()
    rr = Globals.get_rrows()
    rc = Globals.get_rcols()
    br = Globals.get_bor_rows()
    bc = Globals.get_bor_cols()

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

    # druhy cyklus priradi
    # k elementu list elementu v inflows
    for i in rr:
        for j in rc[i]:
            n = len(inflows[i][j])
            tmp = [-int(99)] * n

            for z in range(n):
                ax = inflows[i][j][z][0]
                bx = inflows[i][j][z][1]
                iax = i + ax
                jbx = j + bx
                if Globals.mat_boundary[iax][jbx] < 0.0:
                    tmp[z] = int(-99)
                else:
                    tmp[z] = getEl[iax][jbx]

            getElIN.append(tmp)

    return nEl, getEl, getElIN, getIJ


class ImplicitSolver:

    def __init__(self, iter_crit):

        gl = Globals()
        r = Globals.get_rows()
        c = Globals.get_cols()

        self.total_time = 0

        self.nEl, self.getEl, self.getElIN, self.getIJ = init_getIJel()

        self.b = np.zeros([self.nEl+1], float)
        self.hnew = np.ones([self.nEl+1], float)
        self.hold = np.zeros([self.nEl+1], float)

        # variable counts rills
        self.rill_count = 0

        # iteration criterion
        self.err = 1e-7

        from smoderp2d.processes.surface import sheet_flowb_
        self.sheet_flow = sheet_flowb_

        # if rills are calculated
        if Globals.isRill:
            from smoderp2d.processes.rill import rill
            self.rill_flow = rill
        else:
            from smoderp2d.processes.rill import rill_pass
            self.rill_flow = rill_pass

        from smoderp2d.processes.rainfall import Rainfall
        self.rainfall = Rainfall()

        from smoderp2d.processes.sur_retention import SurRetention
        self.sur_ret = SurRetention()

        from smoderp2d.processes.infiltration import Infiltration
        self.infiltration = Infiltration()

        # opens files for storing hydrographs
        if Globals.points and Globals.points != "#":
            self.hydrographs = wf.Hydrographs()
            arcgis = Globals.arcgis
            if not arcgis:
                with open(os.path.join(Globals.outdir, 'points.txt'), 'w') as fd:
                    for i in range(len(Globals.array_points)):
                        fd.write('{} {} {} {}'.format(
                            Globals.array_points[i][0], Globals.array_points[i][3],
                            Globals.array_points[i][4], os.linesep
                        ))
        else:
            self.hydrographs = wf.HydrographsPass()

    def fillAmat(self, dt):

        # potential precipitation
        PS = self.rainfall.timestepRainfall(self.total_time, dt)

        # calculate current infiltration
        self.infiltration.potential_infiltration(self.total_time, dt)

        # creates empty list for data
        data = []
        self.rill_count = 0

        # toto je pokazde jine
        indptr = [0]
        # toto je pokazde jine
        indices = []

        # takto puvodne ve funkci init_getIJel
        # for iel in range(nEl):
        # indices.append(iel)
        # for inel in getElIN[iel]:
        # indices.append(inel)
        # indptr.append(len(indices))

        for iel in range(self.nEl+1):
            i = self.getIJ[iel][0]
            j = self.getIJ[iel][1]

            # overland outflow
            if self.hnew[iel] > 0:
                hcrit = Globals.get_hcrit(i, j)
                a = Globals.get_aa(i, j)
                b = Globals.get_b(i, j)
                hsheet = min(hcrit, self.hnew[iel])
                hrill = max(0, self.hnew[iel] - hcrit)
                sf = self.sheet_flow(a, b, hsheet)
                rf = 0
                if (hrill > 0):
                    self.rill_count += 1
                    t1 = time.time()
                    rf = self.rill_flow(
                        i, j, hrill, dt) / self.hnew[iel]

                data.append(
                    1. + dt * Globals.dx * (sf) / Globals.pixel_area + dt * rf / Globals.pixel_area)
            else:
                data.append((1.))

            indices.append(iel)

            for inel in self.getElIN[iel]:
                if inel >= 0:  # skip inflows from bc
                    if self.hnew[inel] > 0:
                        ii = self.getIJ[inel][0]
                        jj = self.getIJ[inel][1]
                        hcrit = Globals.get_hcrit(ii, jj)
                        a = Globals.get_aa(ii, jj)
                        b = Globals.get_b(ii, jj)
                        hsheet = min(hcrit, self.hnew[inel])
                        hrill = max(0, self.hnew[inel] - hcrit)
                        sf = self.sheet_flow(a, b, hsheet)
                        rf = 0
                        if (hrill > 0):
                            rf = self.rill_flow(
                                ii, jj, hrill, dt) / self.hnew[iel]

                        data.append(- dt * Globals.dx * (sf) / Globals.pixel_area -
                                    dt * rf / Globals.pixel_area)
                    else:
                        data.append(0)

                    indices.append(inel)

            indptr.append(len(indices))

            # effective precipitation
            ES = self.rainfall.current_rain(i, j, PS)

            inf = self.infiltration.philip_infiltration(i, j)
            if inf >= ES:
                inf = ES

            #self.b[iel] = self.hold[iel] / dt + ES / dt - inf / dt
            self.b[iel] = self.hold[iel] + ES - inf
        self.A = csr_matrix((data, indices, indptr),
                            shape=(self.nEl+1, self.nEl+1), dtype=float)

    def solveStep(self, iter_crit):
        from scipy.sparse.linalg import spsolve

        hewp = self.hnew.copy()
        self.rainfall.save_rainfall_vars()

        err = 1

        while (err > self.err):
            self.rainfall.load_rainfall_vars()
            iter_crit.iter_up()
            self.fillAmat(iter_crit.dt)
            hewp = self.hnew.copy()
            self.hnew = spsolve(self.A, self.b)
            if (iter_crit.max_iter_check(self.total_time)):
                self.hnew = self.hold.copy()
                return 0

            err = np.sum(hewp - self.hnew)**2.0

        self.hnew = self.sur_ret.reduce_h(self.hnew, self.nEl, self.getIJ)
        # write hydrograph record
        for i in Globals.rr:
            for j in Globals.rc[i]:
                self.hydrographs.write_hydrographs_record(
                    i,
                    j,
                    iter_crit,
                    self
                )

        make_sur_raster(self, Globals.outdir, self.total_time)
        return 1
