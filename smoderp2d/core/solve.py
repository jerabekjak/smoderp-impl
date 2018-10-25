from smoderp2d.core.general import Globals
from smoderp2d.tools.tools import make_ASC_raster
from smoderp2d.tools.tools import make_sur_raster
import smoderp2d.processes.rainfall as rain_f
import smoderp2d.processes.infiltration as infilt
import smoderp2d.flow_algorithm.D8 as D8_
from smoderp2d.io_functions import hydrographs as wf

from smoderp2d.providers.logger import Logger
from scipy.sparse import csr_matrix

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
            tmp = [-int(99)] * n

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
            tmp = [-int(99)] * n

            for z in range(n):
                ax = inflows[i][j][z][0]
                bx = inflows[i][j][z][1]
                iax = i + ax
                jbx = j + bx
                tmp[z] = getEl[iax][jbx]

            getElIN.append(tmp)

    # toto je pokazde stejne
    indptr = [0]
    # toto je pokazde stejne
    indices = []

    for iel in range(nEl):
        indices.append(iel)
        for inel in getElIN[iel]:
            indices.append(inel)
        indptr.append(len(indices))

    return nEl, getEl, getElIN, getIJ, indices, indptr


class ImplicitSolver:

    def __init__(self, iter_crit):

        gl = Globals()
        r = gl.get_rows()
        c = gl.get_cols()

        # aktualni cas programu
        self.total_time = 0

        self.nEl, self.getEl, self.getElIN, self.getIJ, self.indices, self.indptr = init_getIJel()

        self.b = np.zeros([self.nEl], float)
        self.hnew = np.ones([self.nEl], float)
        self.hold = np.zeros([self.nEl], float)

        # variable counts rills
        self.rill_count = 0
        
        # iteration criterion
        self.err = 1e-5

        from smoderp2d.processes.surface import sheet_flowb_
        self.sheet_flow = sheet_flowb_

        # if rills are calculated
        if gl.isRill:
            from smoderp2d.processes.rill import rill
            self.rill_flow = rill
        else:
            from smoderp2d.processes.rill import rill_pass
            self.rill_flow = rill_pass
        
        from smoderp2d.processes.rainfall import Rainfall
        self.rainfall = Rainfall()
        
        
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

    def fillAmat(self, dt):
        gl = Globals()

        # potential precipitation
        PS = self.rainfall.timestepRainfall(self.total_time, dt)

        # prepares infiltration table
        for iii in gl.get_combinatIndex():
            index = iii[0]
            k = iii[1]
            s = iii[2]
            iii[3] = infilt.phlilip(
                k, s, dt, self.total_time, gl.get_NoDataValue())

        # creates empty list for data
        data = []
        self.rill_count = 0

        for iel in range(self.nEl):
            i = self.getIJ[iel][0]
            j = self.getIJ[iel][1]


            # overland outflow
            if self.hnew[iel] > 0:
                hcrit = gl.get_hcrit(i, j)
                a = gl.get_aa(i, j)
                b = gl.get_b(i, j)
                hsheet = min(hcrit, self.hnew[iel])
                hrill = max(0, self.hnew[iel] - hcrit)
                sf = self.sheet_flow(a, b, hsheet)
                rf = 0
                if (hrill > 0):
                    self.rill_count += 1
                    t1 = time.time()
                    rf = self.rill_flow(
                        i, j, hrill, dt) / self.hnew[iel]

                else:
                    pass

                data.append(
                    (1. / dt + gl.dx * (sf) / gl.pixel_area) + rf / gl.pixel_area)
            else:
                data.append((1. / dt))
            #print max(gl.dx * (sf) / gl.pixel_area*dt/gl.dx, rf/ gl.pixel_area * dt / gl.dx)
            # TODO to by meli byt jiny acka a becka
            # pokud to vteka z jineho lu nebo pudy
            for inel in self.getElIN[iel]:
                if self.hnew[inel] > 0:
                    i = self.getIJ[inel][0]
                    j = self.getIJ[inel][1]
                    hcrit = gl.get_hcrit(i, j)
                    a = gl.get_aa(i, j)
                    b = gl.get_b(i, j)
                    hsheet = min(hcrit, self.hnew[inel])
                    hrill = max(0, self.hnew[inel] - hcrit)
                    sf = self.sheet_flow(a, b, hsheet)
                    rf = 0
                    if (hrill > 0):
                        rf = self.rill_flow(
                            i, j, hrill, dt) / self.hnew[iel]

                    data.append(-gl.dx * (sf) / gl.pixel_area -
                                rf / gl.pixel_area)
                else:
                    data.append(0)



            # infiltration
            inf = infilt.philip_infiltration(
                gl.get_mat_inf_index(i, j), gl.get_combinatIndex())
            if inf >= self.hold[iel]:
                inf = self.hold[iel]
                
            # effective precipitation
            ES = self.rainfall.current_rain(i,j,PS)
                
            self.b[iel] = self.hold[iel] / dt + ES / dt - inf / dt

        self.A = csr_matrix((data, self.indices, self.indptr),
                            shape=(self.nEl, self.nEl), dtype=float)

    def solveStep(self, iter_crit):
        from scipy.sparse.linalg import spsolve

        hewp = self.hnew.copy()

        err = 1

        while (err > self.err):

            iter_crit.iter_up()
            t1 = time.time()
            self.fillAmat(iter_crit.dt)
            t2 = time.time()
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

        make_sur_raster(self, Globals.outdir, self.total_time)
        return 1
