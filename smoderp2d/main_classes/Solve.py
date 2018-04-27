from smoderp2d.main_classes.General import Globals
from smoderp2d.exceptions import MaxIterationExceeded
import smoderp2d.processes.rainfall as rain_f
import smoderp2d.processes.infiltration as infilt
import smoderp2d.flow_algorithm.D8 as D8_

from scipy.sparse import csr_matrix
import numpy as np
import os
import sys


def make_ASC_raster(name_, numpy_arr):

    gl = Globals()

    rr = gl.rr
    rc = gl.rc
    br = gl.br
    bc = gl.bc
    nrows = gl.r
    ncols = gl.c

    tmpStr = str(numpy_arr.dtype)[0:3]
    if tmpStr == 'int':
        noData = gl.NoDataInt
    else:
        noData = gl.NoDataValue

    tmp = np.copy(numpy_arr)
    tmp.fill(noData)

    f = open(name_, 'w')
    f.write("ncols " + str(ncols) + '\n')
    f.write("nrows " + str(nrows) + '\n')
    f.write("xllcorner " + str(gl.xllcorner) + '\n')
    f.write("yllcorner " + str(gl.yllcorner) + '\n')
    f.write("cellsize " + str(gl.dx) + '\n')
    f.write("nodata_value " + str(noData) + '\n')

    for i in rr:
        for j in rc[i]:
            tmp[i][j] = numpy_arr[i][j]

    # for i in br:
        # for j in bc[i]:
            # tmp[i][j] = numpy_arr[i][j]

    for i in range(nrows):
        line = ""
        for j in range(ncols):
            line += str(tmp[i][j]) + "\t"
        line += '\n'
        f.write(line)


def make_sur_raster(IS, output, t):
    gl = Globals()

    rrows = gl.rr
    rcols = gl.rc

    r = gl.r
    c = gl.c

    arr = np.zeros([r, c], float)

    for i in rrows:
        for j in rcols[i]:
            el = IS.IJtoEl_a[i][j]
            arr[i][j] = IS.hnew[el]

    outName = output + os.sep + str(int(t)).zfill(10) + 'h' + ".asc"
    make_ASC_raster(outName, arr)


def init_getIJel():
    """Documentation for a function init_getIJel()

    popis:
    funkce ma vytvorit pole a vektory na zjistovani
    odpovidajiciho radku a sloupce rasteru podle pozive v lin soustave
    a odpovidajiciho elementu v soustave podle sloupce a rarku v rasteru

    return el      - pocet elementu v soustave
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

    el = int(-1)
    # prvni cyklus priradi
    # k i j bunky jeji poradi ve vektoru
    # a k elementu vektoru i j bunky
    for i in rr:
        for j in rc[i]:
            el += int(1)
            getIJ.append([i, j])
            getEl[i][j] = el

    for i in br:
        for j in bc[i]:
            el += int(1)
            getIJ.append([i, j])
            getEl[i][j] = el

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

    for iel in range(el):
        indices.append(iel)
        for inel in getElIN[iel]:
            indices.append(inel)
        indptr.append(len(indices))

    return el, getEl, getElIN, getIJ, indices, indptr


class ImplicitSolver:

    def __init__(self):

        gl = Globals()
        r = gl.get_rows()
        c = gl.get_cols()

        # interval srazky
        self.tz = 0
        # aktualni cas programu
        self.total_time = 0

        self.nEl, self.IJtoEl_a, self.ELinEL_l, self.ELtoIJ, self.indices, self.indptr = init_getIJel()

        self.A = np.zeros([self.nEl, self.nEl], float)
        self.b = np.zeros([self.nEl], float)
        self.hnew = np.ones([self.nEl], float)
        self.hold = np.zeros([self.nEl], float)

    def fillAmat(self, dt):

        gl = Globals()

        PS, self.tz = rain_f.timestepRainfall(self.total_time, dt, self.tz)

        for iii in gl.get_combinatIndex():
            index = iii[0]
            k = iii[1]
            s = iii[2]
            iii[3] = infilt.phlilip(
                k, s, dt, self.total_time, gl.get_NoDataValue())

        data = []
        
        for iel in range(self.nEl):

            i = self.ELtoIJ[iel][0]
            j = self.ELtoIJ[iel][1]
            hcrit = gl.get_hcrit(i, j)
            a = gl.get_mat_aa(i, j)
            b = gl.get_mat_b(i, j)

            # infiltration
            inf = infilt.philip_infiltration(
                gl.get_mat_inf_index(i, j), gl.get_combinatIndex())
            if inf >= self.hold[iel]:
                inf = self.hold[iel]
            
            
            # overlad outflow    
            if self.hnew[iel] > 0:
                hsheet = gl.hcrit[i][j]
                hrill  = gl.hcrit[i][j]
                data.append(
                    (1. / dt + gl.dx * (a * self.hnew[iel]**(b - 1)) / gl.pixel_area))
            else:
                data.append((1. / dt))
                
            
            for inel in self.ELinEL_l[iel]:
                if self.hnew[inel] > 0:
                    data.append(-(gl.dx * a *
                                  (self.hnew[inel]**(b - 1)) / gl.pixel_area))
                else:
                    data.append(0)
            
            #print self.hold[iel] / dt, PS / dt, inf / dt
            self.b[iel] = self.hold[iel] / dt + PS / dt - inf / dt

        self.A = csr_matrix((data, self.indices, self.indptr),
                            shape=(self.nEl, self.nEl), dtype=float)


    def solveStep(self, dt):
        from scipy.sparse.linalg import spsolve

        iter_ = 0
        maxIter = 10
        hewp = self.hnew.copy()
        hewp.fill(0.0)
        
        while (abs(np.sum((hewp - self.hnew))) > 0.000001):
            iter_ += 1
            #print iter_
            self.fillAmat(dt)
            hewp = self.hnew.copy()
            self.hnew = spsolve(self.A, self.b)
            #raw_input()
            if (iter_ > maxIter):
                raise  MaxIterationExceeded(maxIter, self.total_time)


        #print self.hnew[20]
        self.total_time += dt
        print self.hnew[10], self.hnew[20]
        make_sur_raster(self, 'out', self.total_time)
