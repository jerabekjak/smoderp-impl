! gfortran -O3 -funroll-loops -ffast-math -floop-strip-mine -shared -fPIC -o fill_a_mat.so fill_a_mat.f90



subroutine fill_a_mat(nEl, sizes, getIJ, data, hnew, mat_aa) bind(c, name='fill_a_mat')
  use iso_c_binding, only: c_float, c_int
  implicit none
  
  integer(c_int), intent(in) :: nEl
  integer(c_int), dimension(3), intent(in) :: sizes
  integer(c_int), dimension(nEl+1,2), intent(in) :: getIJ
  real(c_float), dimension(sizes(1)), intent(out) :: data
  real(c_float), dimension(nEl), intent(in) :: hnew
  real(c_float), dimension(1:sizes(2),1:sizes(3)), intent(in) :: mat_aa

  integer :: iel, i, j, n
  
  do iel = 1, 10
    i = getIJ(iel,1)
    j = getIJ(iel,2)
    print *, i,j 
    
    
  end do
  
  
  
  
  

end subroutine fill_a_mat