! gfortran -O3 -funroll-loops -ffast-math -floop-strip-mine -shared -fPIC -o fill_a_mat.so fill_a_mat.f90








subroutine fill_a_mat(nEl, sizes, getIJ, data, hnew, hold, &
                      mat_aa, mat_inf_index, &
                      combinat_intex) bind(c, name='fill_a_mat')
    use iso_c_binding, only: c_float, c_int
    implicit none
    
    integer(c_int), intent(in) :: nEl
    integer(c_int), dimension(5), intent(in) :: sizes
    integer(c_int), dimension(nEl+1,2), intent(in) :: getIJ
    real(c_float), dimension(sizes(1)), intent(out) :: data
    real(c_float), dimension(nEl), intent(in) :: hnew, hold
    real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_aa
    integer(c_int), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_inf_index
    real(c_float), dimension(1:sizes(4),1:sizes(5)), intent(in) :: combinat_intex

    integer :: iel, i, j, n
    real(c_float)    :: inf
    
    real(c_type) :: sf, rf
    
    
    do iel = 1, nEl
        i = getIJ(iel,1)
        j = getIJ(iel,2)
        
        inf = philips_infiltration(mat_inf_index(i,j),combinat_intex)
        
        if (hnew(iel) > 0.0_c_float) then
            hcrit = mat_hcrit(i,j)
            a     = mat_aa(i,j)
            b     = mat_b(i,j)
            hsheet = min(hcrit, hnew(iel))
            hrill  = max(0, hnew(iel) - hcrit)
            
            sf = ! funkce 
            rf = 0.0_c_float
            
            if (hrill > 0.0_c_float) then
                rf
            end if 
            
        end if
        
        
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

  

end subroutine fill_a_mat