! This subroutine uses the trilinear interpolation to interpolate the non-local part of the pseudopotential.
! Algorithm from Wikipedia


subroutine  trilinear_obf_nc(kptns,ik,sh,c)
  USE kinds, ONLY : DP 
  USE constants, ONLY : pi   ! DEBUG
  USE io_global, ONLY : stdout ! DEBUG
  USE obf_diag_objects
  !
  IMPLICIT NONE
  !
  TYPE(kpoints) :: kptns
  TYPE(shirley) :: sh
  INTEGER, INTENT(in) :: ik ! dense local k-grid index
  COMPLEX(kind=DP), INTENT(out) :: c(sh%nkb,sh%npol,sh%ntot_e)  ! Trilinear interpolated projector (non-collinear)
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: c000, c100, c010, c110, c001, c101, c011, c111 
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: c00, c01, c10, c11, c0, c1
  REAL(kind=DP) ::  delta(3)

  allocate(c000(sh%nkb,sh%npol,sh%ntot_e),c100(sh%nkb,sh%npol,sh%ntot_e),&
  & c010(sh%nkb,sh%npol,sh%ntot_e),c110(sh%nkb,sh%npol,sh%ntot_e), &
  & c001(sh%nkb,sh%npol,sh%ntot_e), c101(sh%nkb,sh%npol,sh%ntot_e), &
  & c011(sh%nkb,sh%npol,sh%ntot_e), c111(sh%nkb,sh%npol,sh%ntot_e) )

  allocate(c00(sh%nkb,sh%npol,sh%ntot_e),c01(sh%nkb,sh%npol,sh%ntot_e),&
  & c10(sh%nkb,sh%npol,sh%ntot_e),c11(sh%nkb,sh%npol,sh%ntot_e), &
  & c0(sh%nkb,sh%npol,sh%ntot_e), c1(sh%nkb,sh%npol,sh%ntot_e))
  
  c000(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(1,ik))
  c100(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(2,ik))
  c010(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(3,ik))
  c110(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(4,ik))
  c001(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(5,ik))
  c101(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(6,ik))
  c011(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(7,ik))
  c111(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(8,ik))
     
  delta(1:3) = kptns%pos_cube(1:3,ik)
  
  c00 = c000*(1.0 - delta(1)) + c100*delta(1)
  c01 = c001*(1.0 - delta(1)) + c101*delta(1)
  c10 = c010*(1.0 - delta(1)) + c110*delta(1)
  c11 = c011*(1.0 - delta(1)) + c111*delta(1)

  c0 = c00*(1.0 - delta(2)) + c10*delta(2)
  c1 = c01*(1.0 - delta(2)) + c11*delta(2)
  
  c = c0*(1.0 - delta(3)) + c1*delta(3)

  deallocate(c000, c100, c010, c110, c001, c101, c011, c111, c00, c01, c10, c11, c0, c1)

end subroutine trilinear_obf_nc

subroutine  trilinear_obfc(kptns,ik,sh,c)
  USE kinds, ONLY : DP 
  USE obf_diag_objects

  implicit none

  TYPE(kpoints) :: kptns
  TYPE(shirley) :: sh
  INTEGER, INTENT(in) :: ik ! dense local k-grid index
  COMPLEX(kind=DP), INTENT(out) :: c(sh%nkb,sh%ntot_e)  ! Trilinear interpolated projector (collinear)
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: c000, c100, c010, c110, c001, c101, c011, c111 
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: c00, c01, c10, c11, c0, c1
  REAL(kind=DP) ::  delta(3)

  allocate(c000(sh%nkb,sh%ntot_e),c100(sh%nkb,sh%ntot_e),&
  & c010(sh%nkb,sh%ntot_e),c110(sh%nkb,sh%ntot_e), &
  & c001(sh%nkb,sh%ntot_e), c101(sh%nkb,sh%ntot_e), &
  & c011(sh%nkb,sh%ntot_e), c111(sh%nkb,sh%ntot_e) )

  allocate(c00(sh%nkb,sh%ntot_e),c01(sh%nkb,sh%ntot_e),&
  & c10(sh%nkb,sh%ntot_e),c11(sh%nkb,sh%ntot_e), &
  & c0(sh%nkb,sh%ntot_e), c1(sh%nkb,sh%ntot_e))
 
  c000(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(1,ik))
  c100(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(2,ik))
  c010(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(3,ik))
  c110(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(4,ik))
  c001(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(5,ik))
  c101(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(6,ik))
  c011(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(7,ik))
  c111(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(8,ik))
     
  delta(1:3) = kptns%pos_cube(1:3,ik)
  
  c00 = c000*(1.0 - delta(1)) + c100*delta(1)
  c01 = c001*(1.0 - delta(1)) + c101*delta(1)
  c10 = c010*(1.0 - delta(1)) + c110*delta(1)
  c11 = c011*(1.0 - delta(1)) + c111*delta(1)

  c0 = c00*(1.0 - delta(2)) + c10*delta(2)
  c1 = c01*(1.0 - delta(2)) + c11*delta(2)
  
  c = c0*(1.0 - delta(3)) + c1*delta(3)

  deallocate(c000, c100, c010, c110, c001, c101, c011, c111, c00, c01, c10, c11, c0, c1)

end subroutine trilinear_obfc

