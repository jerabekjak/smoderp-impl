# USAGE: ! gfortran -O3 -funroll-loops -ffast-math -floop-strip-mine -shared -fPIC -o fill_a_mat.so fill_a_mat.f90

subroutine fill_a_mat(nEl, sizes, getIJ, getElIN, data, hnew, hold, &
                      mat_aa, mat_b, mat_hcrit, mat_eff_vrst, mat_inf_index, &
                      mat_n, mat_slope, &
                      combinat_intex, &
                      dx, dt) bind(c, name='fill_a_mat')
    use iso_c_binding, only: c_float, c_int
    implicit none
    
    integer(c_int), intent(in) :: nEl
    integer(c_int), dimension(5), intent(in) :: sizes
    integer(c_int), dimension(nEl+1,2), intent(in) :: getIJ
    integer(c_int), dimension(nEl+1,8), intent(in) :: getElIN
    real(c_float), dimension(1:sizes(1)), intent(out) :: data
    real(c_float), dimension(nEl), intent(in) :: hnew, hold
    real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_aa
    real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_b
    real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_hcrit
    real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_n
    real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_slope
    real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_eff_vrst
    integer(c_int), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_inf_index
    real(c_float), dimension(1:sizes(4),1:sizes(5)), intent(in) :: combinat_intex
    
    real(c_float), intent(in) :: dx, dt

    integer :: iel, inel,  i, j, n
    integer :: idata
    real(c_float)    :: inf
    real ::  pixel_area
    
    real(c_float) :: sf, rf, a, b, hcrit, hrill, hsheet
    
    idata = 1
    data(idata) = 1.
    pixel_area = dx**2.
    
    do iel = 1, nEl
        i = getIJ(iel,1)
        j = getIJ(iel,2)
        
        inf = philips_infiltration(mat_inf_index(i,j),combinat_intex)
        
        if (hnew(iel) > 0.0_c_float) then
        
            hcrit = mat_hcrit(i,j)
            a     = mat_aa(i,j)
            b     = mat_b(i,j)
            hsheet = min(hcrit, hnew(iel))
            hrill  = max(0.0, hnew(iel) - hcrit)
            
            sf = a*hsheet**b
            rf = 0.0_c_float
            if (hrill > 0.0_c_float) then
                rf = rill_flow(hrill, mat_eff_vrst(i,j),dx**2., mat_n(i,j), mat_slope(i,j))
            end if 
            
            data(idata) = 1. / dt + dx * (sf) / pixel_area + rf / pixel_area
            idata = idata + 1
            
        else
            data(idata) = 1. / dt
            idata = idata + 1
        end if    
        
        do inel = getElIN(iel,1), getElIN(iel,8)
            if (inel > 0) then
                i = getIJ(inel,1)
                j = getIJ(inel,2)
            end if
        end do

end do
    
  
  
 contains
 
    function philips_infiltration(soil, ci) result(inf)
        use iso_c_binding, only: c_float, c_int
        implicit none
        
        integer(c_int), intent(in) :: soil
        real(c_float), dimension(:,:), intent(in) :: ci
        real(c_float):: inf

        integer :: i, n
        
        n = ubound(ci,2)
        
        do i=1,n
            if ( ci(i,1) == soil ) then
                inf = ci(i,4)
                exit
            end if 
        end do
        
    end function philips_infiltration
    
    function rill_flow(hrill, ef_cont,pa, n, slope) result(rf)
        use iso_c_binding, only: c_float, c_int
        implicit none

        real, intent(in) :: hrill
        real, intent(in) :: ef_cont
        real, intent(in) :: pa
        real, intent(in) :: n
        real, intent(in) :: slope
        real             :: rf
        
        real :: b = 0
        real :: v_to_rill
        real :: v
        real :: r_rill 
        real :: hloc
        
        
        v_to_rill = hrill*pa
        
        call update_hb(v_to_rill, ef_cont, b, hloc)
        
        r_rill = (hloc*b)/(b + 2.0*hloc)
        v = r_rill**(2./3.) * 1./n * (slope/100.)**0.5  
        
        rf = hloc*b*v
        
        
    end function rill_flow
    
    
    subroutine update_hb(v_to_rill, ef_cont, b, h) 
        use iso_c_binding, only: c_float, c_int
        implicit none
        
        real, intent(in) :: v_to_rill
        real, intent(in) :: ef_cont
        real, intent(inout) :: b
        real, intent(out) :: h
        
        real :: rr, l, newb
      
        rr = 0.7
        
        l = ef_cont
        
        if (v_to_rill < 0) then
            print *, 'error in update_hb'
            error stop
        end if 
        
        newb = sqrt(v_to_rill/(rr*l))
        
        if (v_to_rill > 0) then
            b = newb
            h = v_to_rill/(b*l)
        else
            h = v_to_rill/(b*l)
        end if 

        
        
        
        
        
    end subroutine update_hb
    
    
end subroutine fill_a_mat