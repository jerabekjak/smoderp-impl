import os
import sys
import numpy as np
from smoderp2d.core.general import Globals


# Class with the debug mark
#
#  method self.mark with given name return konsole message with generically increasing integer value


class DebugMark:
    n = 0

    def mark(self, name_, info=''):
        DebugMark.n += 1


# Class with the file name generation
#
#  method self.gen with given name returns the name plus increasing int number 000*
class FileNameGen:
    def __init__(self):
        self.n = 0

    def gen(self, name_, info=''):
        filename = name_+str(self.n).zfill(4)
        self.n += 1
        return filename


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
            el = IS.getEl[i][j]
            arr[i][j] = IS.hnew[el]

    outName = output + os.sep + str(int(t)).zfill(10) + 'h' + ".asc"
    make_ASC_raster(outName, arr)





#def make_sub_raster(subArr, G, t, output):
    #rrows = G.rr
    #rcols = G.rc
    #arr = np.zeros(np.shape(subArr), float)
    #for i in rrows:
        #for j in rcols[i]:
            #arr[i][j] = subArr[i][j].h

    #outName = output+os.sep+'prubeh'+os.sep+str(int(t)).zfill(10)+'hsub'+".asc"
    #make_ASC_raster(outName, arr, G)
## Creates a ascii raster from the numpy array
##


#def make_ASC_raster(name_, numpy_arr, G):
    #rr = G.rr
    #rc = G.rc
    #br = G.br
    #bc = G.bc
    #nrows = G.r
    #ncols = G.c

    #tmpStr = str(numpy_arr.dtype)[0:3]
    #if tmpStr == 'int':
        #noData = G.NoDataInt
    #else:
        #noData = G.NoDataValue

    #tmp = np.copy(numpy_arr)
    #tmp.fill(noData)

    #f = open(name_, 'w')
    #f.write("ncols " + str(ncols) + '\n')
    #f.write("nrows " + str(nrows) + '\n')
    #f.write("xllcorner " + str(G.xllcorner) + '\n')
    #f.write("yllcorner " + str(G.yllcorner) + '\n')
    #f.write("cellsize " + str(G.dx) + '\n')
    #f.write("nodata_value " + str(noData) + '\n')

    #for i in rr:
        #for j in rc[i]:
            #tmp[i][j] = numpy_arr[i][j]

    ## for i in br:
        ## for j in bc[i]:
            ##tmp[i][j] = numpy_arr[i][j]

    #for i in range(nrows):
        #line = ""
        #for j in range(ncols):
            #line += str(tmp[i][j]) + "\t"
        #line += '\n'
        #f.write(line)


# Returns boolean information about the components of the computation
#
#  Return 4 true/values for rill, subflow, stream, diffuse presence/non-presence.\n
#  Optionally string parameter co_ specify the process user ask for rill | subflow | stream | diffuse
#

# Returns boolean information about the components of the computation
#
#  Return 4 true/values for rill, subflow, stream, diffuse presence/non-presence.\n
#  Optionally string parameter co_ specify the process user ask for rill | subflow | stream | diffuse
#
def comp_type(tc, co_=""):
    """
    string_type_of_coputing = get_argv(constants.PARAMETER_TYPE_COMPUTING)
    string_type_of_coputing = string_type_of_coputing.lower().replace(' ','').replace(',','')  #jj .lower().replace(' ','').replace(',','') udela ze vsecho v tom stringu maly pismena, replace vyhodi mezery a carky

    if string_type_of_coputing == "onlyshallowsurface":
        type_of_computing = 0
    elif string_type_of_coputing == "shallowandrillsurface":
        type_of_computing = 1
    elif string_type_of_coputing == "diffuseshallowsurface":
        type_of_computing = 2
    elif string_type_of_coputing == "shallowrillstreamsurface":
        type_of_computing = 3
    elif string_type_of_coputing == "surfaceandsubsurfaceflow":
        type_of_computing = 4
    elif string_type_of_coputing == "surfaceandsubsurfacestreamflow":
        type_of_computing = 5
    """

    diffuse = False
    subflow = False
    stream = False
    rill = False
    only_surface = False

    if int(tc) == 1:
        rill = True

    if int(tc) == 3:
        stream = True
        rill = True

    if int(tc) == 4:
        subflow = True
        rill = True

    if int(tc) == 5:
        stream = True
        subflow = True
        rill = True
    if int(tc) == 0:
        only_surface = True
    if (co_ == "rill"):
        return rill
    if (co_ == "subflow"):
        return subflow
    if (co_ == "stream"):
        return stream
    if (co_ == "diffuse"):
        return diffuse
    if (co_ == "surface"):
        return only_surface

    return rill, subflow, stream, diffuse


def int_comp_type(intco_):

    if intco_ == 0:
        return 'onlyshallowsurface'
    elif intco_ == 1:
        return 'shallowandrillsurface'
    elif intco_ == 2:
        return 'diffuseshallowsurface'
    elif intco_ == 3:
        return 'shallowrillstreamsurface'
    elif intco_ == 4:
        return 'surfaceandsubsurfaceflow'
    elif intco_ == 5:
        return 'surfaceandsubsurfacestreamflow'
    else:
        print 'error in data_preparation, PARAMETER_TYPE_COMPUTING error'

# Returns input parameter from the sys.argv
#


def get_argv(id_):
    iid_ = id_
    #print '!!!Bacha tools.get_argv jsou plus 1!!!';
    iid_ += 1
    return sys.argv[iid_]


def set_argv(id_, value):
    iid_ = id_
    #print '!!!Bacha tools.get_argv jsou plus 1!!!';
    iid_ += 1
    print sys.argv[iid_]
    sys.argv[iid_] = value


def prt_sys_argv():
    for item in sys.argv:
        print item


# Transfers string 'true'/'false' in sys.argv to a logical True/False
#
def logical_argv(id_):
    iid_ = id_
    #print '!!!Bacha tools.get_argv jsou plus 1!!!';
    iid_ += 1
    if type(sys.argv[iid_]) == str:
        if sys.argv[iid_].lower().strip() == 'true':
            sys.argv[iid_] = True
        elif sys.argv[iid_].lower().strip() == 'false':
            sys.argv[iid_] = False
        else:
            sys.exit('Logical parameter are not assign correctly...')
    if type(sys.argv[iid_]) == int:
        sys.exit('Logical parameter are not assign correctly...')
    if type(sys.argv[iid_]) == float:
        sys.exit('Logical parameter are not assign correctly...')
