from smoderp2d.core.general import Globals

class SurRetention (object):
    
    def __init__(self):
        """ Create array for accumulating the surface runoff. """
        self.sum_sur_retention = Globals.get_mat_reten().copy()
        
    def reduce_h(self,h,n,getij):
        """ reduce the solution by surface retention
        
        :param array h     : total water level
        :param integer n   : number of element 
        :param array getij : i j position in mask array of n'th element
        
        
        :return array h: reduced solution
        """
        for iel in range(n):
            i = getij[iel][0]
            j = getij[iel][1]
            
            #ret is negative calue
            ret = self.sum_sur_retention[i][j]
            
            # if retention is filled
            if (ret >= 0) :
                return h
            if (h[iel] <= -ret) :
                ret += h[iel]
                h[iel] = 0
            else : 
                h[iel] = min(h[iel],h[iel]+ret)
                ret += h[iel]
                
            self.sum_sur_retention[i][j] = ret
            
        return h
        
        