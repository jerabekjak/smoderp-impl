import sys
import numpy as np
import os
from smoderp2d.tools.tools import comp_type
import smoderp2d.io_functions.prt as prt
from smoderp2d.tools.tools import get_argv
import smoderp2d.constants as constants


from smoderp2d.core.general import Globals as Gl

# extraout = get_argv(constants.PARAMETER_EXTRA_OUTPUT)

# rill, subflow, stream, diffuse = comp_type()

rill = Gl.isRill
subflow = Gl.subflow
stream = Gl.isStream
extraout = Gl.extraOut


class Hydrographs:

    def __init__(self):

        gl = Gl()

        points = gl.get_array_points()
        ipi = points.shape[0]
        jpj = 5
        point_int = [[0] * jpj for i in range(ipi)]

        outdir = gl.get_outdir()
        rr = gl.get_rrows()
        rc = gl.get_rcols()
        pixel_area = gl.get_pixel_area()

        self.inSurface = []
        self.inStream = []

        for ip in range(ipi):
            for jp in [0, 1, 2]:
                point_int[ip][jp] = int(points[ip][jp])

        for ip in range(ipi):
            for jp in [3, 4]:
                point_int[ip][jp] = points[ip][jp]

        # for ttt in point_int:
            # print ttt

        # tento cylkus meze budy, ktere jsou
        # v i,j cylku o jednu vedle rrows a rcols
        outsideDomain = False
        del_ = []
        for ip in range(ipi):
            l = point_int[ip][1]
            m = point_int[ip][2]
            for ipp in rr:
                if l == ipp:
                    for jpp in rc[ipp]:
                        if m == jpp:
                            outsideDomain = True
            if not(outsideDomain):
                del_.append(ip)
            outsideDomain = False
        point_int = [i for j, i in enumerate(point_int) if j not in del_]
        ipi -= len(del_)

        # for ttt in point_int:
            # print ttt

        # for ip in range(ipi):
            # l = point_int[ip][1]
            # m = point_int[ip][2]
            # print mat_tok_usek[ip][jp]

        counter = 0

        # mat_tok_usek is alway presented if stream == True
        # if (mat_tok_usek != None) and (stream == True):
        if (stream):
            for ip in range(ipi):
                l = point_int[ip][1]
                m = point_int[ip][2]

                if gl.get_mat_tok_reach(l, m) >= 1000:
                    self.inStream.append(counter)
                    counter += 1
                else:
                    self.inSurface.append(counter)
                    counter += 1
        else:
            self.inSurface = [i for i in range(ipi)]

        self.inStream.append(-99)
        self.inSurface.append(-99)

        self.n = ipi
        self.point_int = point_int
        self.subflow = subflow
        self.rill = rill
        self.stream = stream
        self.pixel_area = pixel_area
        # print self.point_int
        # raw_input()

        iStream = 0
        iSurface = 0

        self.header = []

        for i in range(self.n):

            if i == self.inStream[iStream]:

                header = '# Hydrograph at the point with coordinates: ' + \
                    str(self.point_int[i][3]) + ' ' + str(
                        self.point_int[i][4]) + '\n'
                header += '# A pixel size is [m2]:\n'
                header += '# ' + str(self.pixel_area) + '\n'

                if not(extraout):
                    header += '# time[s];deltaTime[s];rainfall[m];reachWaterLevel[m];reachFlow[m3/s];reachVolRunoff[m3]\n'
                else:
                    header += '# time[s];deltaTime[s];Rainfall[m];Waterlevel[m];V_runoff[m3];Q[m3/s];V_from_field[m3];V_rests_in_stream[m3]\n'
                self.header.append(header)
                iStream += 1

            elif i == self.inSurface[iSurface]:
                header = '# Hydrograph at the point with coordinates: ' + \
                    str(self.point_int[i][3]) + ' ' + str(
                        self.point_int[i][4]) + '\n'
                header += '# A pixel size is [m2]:\n'
                header += '# ' + str(self.pixel_area) + '\n'
                header += '# time[s];hsheet[m];hrill[m];htot[m];dt[s];n-iter[-]'
                header += '\n'
                iSurface += 1
                self.header.append(header)

        self.files = []
        for i in range(self.n):
            name_ = outdir + os.sep + 'point' + \
                str(self.point_int[i][0]).zfill(3) + '.dat'
            file_ = open(name_, 'w')
            file_.writelines(self.header[i])
            self.files.append(file_)

        del self.inStream[-1]
        del self.inSurface[-1]

        prt.message("Hydrographs files has been created...")

    def write_hydrographs_record(
            self, i, j, iter_crit, LS, first=False, isStream = False,  sep=';'):
        gl = Gl()
        
        if (first) :
            for ip in self.inSurface:
                l = self.point_int[ip][1]
                m = self.point_int[ip][2]
                if i == l and j == m:
                    hcrit = gl.get_hcrit(i,j)
                    #Logger.debug('hcrit natvrdo')
                    #hcrit = 0.01 
                    hsheet = min(hcrit,LS.hold[LS.getEl[i][j]])
                    hrill  = max(0,    LS.hold[LS.getEl[i][j]]-hcrit)
                    line = str(LS.total_time + iter_crit.dt) + sep
                    line += str(hsheet) + sep
                    line += str(hrill) + sep
                    line += str(LS.hold[LS.getEl[i][j]]) + sep
                    line += str(iter_crit.dt) + sep
                    line += str(iter_crit.iter_)
                    line += '\n'
                    self.files[ip].writelines(line)    
        else :
            for ip in self.inSurface:
                l = self.point_int[ip][1]
                m = self.point_int[ip][2]
                if i == l and j == m:
                    hcrit = gl.get_hcrit(i,j)
                    #Logger.debug('hcrit natvrdo')
                    #hcrit = 0.01
                    hsheet = min(hcrit,LS.hnew[LS.getEl[i][j]])
                    hrill  = max(0,    LS.hnew[LS.getEl[i][j]]-hcrit)
                    line = str(LS.total_time + iter_crit.dt) + sep
                    line += str(hsheet) + sep
                    line += str(hrill) + sep
                    line += str(LS.hnew[LS.getEl[i][j]]) + sep
                    line += str(iter_crit.dt) + sep
                    line += str(iter_crit.iter_)
                    line += '\n'
                    self.files[ip].writelines(line)   
            
        
  


        
        
        
        
        
        

        ##ratio = fc.ratio
        ##total_time = fc.total_time + dt
        ##iter_ = fc.iter_
        #currRain = 0    

        ##courantMost = courant.cour_most
        ##courantRill = courant. cour_most_rill
        #inStream = False

        #if inStream:
            #for ip in self.inStream:
                #l = self.point_int[ip][1]
                #m = self.point_int[ip][2]
                #line = str(total_time) + sep
                #line += str(dt) + sep
                #line += str(currRain) + sep
                #line += surface.return_stream_str_vals(l, m, sep, dt, extraout)
                #line += '\n'
                #self.files[ip].writelines(line)

        #else:
            #for ip in self.inSurface:
                #l = self.point_int[ip][1]
                #m = self.point_int[ip][2]
                #if i == l and j == m:
                    #line = str(total_time) + sep
                    #line += str(dt) + sep
                    #line += str(currRain) + sep
                    #linebil = surface.return_str_vals(l, m, sep, dt, extraout)
                    #line += linebil[0] + sep
                    #line += str(linebil[1])  # + sep
                    ## line += subsurface.return_str_vals(l,m,sep,dt) + sep   #
                    ## prozatim
                    #if extraout:
                        #line += sep + str(surface.arr[l][m].V_to_rill) + sep
                        #line += str(ratio) + sep
                        #line += str(courantMost) + sep
                        #line += str(courantRill) + sep
                        #line += str(iter_)

                    #line += '\n'
                    #self.files[ip].writelines(line)

        
        
    # def write_hydrographs_usek(self,dt,total_time,surface,currRain,sep=';'):
        # line = str(total_time) + sep
        # line += str(dt) + sep
        # line += str(currRain) + sep
        # line += surface.return_stream_str_vals(0,0,sep,dt)
        # line += '\n'
        # self.tokusek.writelines(line)

    def closeHydrographs(self):
        for i in range(self.n):
            self.files[i].close()


class HydrographsPass:

    def write_hydrographs_record(
            self, i, j, fc, courant, dt, surface, subsurface, currRain, inStream=False, sep=';'):
        pass

    def closeHydrographs(self):
        pass
