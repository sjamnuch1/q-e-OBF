!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE gk_sort_ngw_lmax( k, ngm, g, ecut, ngk, igk, gk,ngw_lmax )
   !----------------------------------------------------------------------------
   !! Sorts k+g in order of increasing magnitude, up to ecut.
   !
   !! NB: this version should yield the same ordering for different ecut
   !!     and the same ordering in all machines AS LONG AS INPUT DATA
   !!     IS EXACTLY THE SAME
   !! Find igk k+G up until ngw_lmax. We no longer care about ngk 
   !
   USE kinds,      ONLY: DP
   USE constants,  ONLY: eps8
   !USE wvfct,      ONLY: npwx              !SJ
   USE input_basis,         ONLY : npw_max !SJ
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN)  :: k(3)
   !! the k point
   INTEGER,  INTENT(IN)  :: ngm
   !! the number of g vectors
   REAL(DP), INTENT(IN)  :: g(3,ngm)
   !! the coordinates of G vectors
   REAL(DP), INTENT(IN)  :: ecut
   !! the cut-off energy
   INTEGER,  INTENT(IN)  :: ngw_lmax
   !! the number of g vectors
   INTEGER,  INTENT(OUT) :: ngk
   !! the number of k+G vectors inside the "ecut sphere"
   INTEGER,  INTENT(OUT) :: igk(ngw_lmax)
   !! the correspondence k+G <-> G
   REAL(DP), INTENT(OUT) :: gk(ngw_lmax)
   !! the moduli of k+G
   !
   !  ... local variables
   !
   INTEGER :: ng   ! counter on   G vectors
   INTEGER :: nk   ! counter on k+G vectors
 REAL(DP) :: q   ! |k+G|^2
   REAL(DP) :: q2x ! upper bound for |G|
   !
   ! ... first we count the number of k+G vectors inside the cut-off sphere
   !
   q2x = ( SQRT( SUM(k(:)**2) ) + SQRT( ecut ) )**2
   !
   !WRITE(6,*) 'In gk_sort() : q(k)',k(:)
   ngk = 0
   igk(:) = 0
   gk (:) = 0.0_DP
   !
   DO ng = 1, ngw_lmax
      q = SUM( ( k(:) + g(:,ng) )**2 )
      IF ( q <= eps8 ) q = 0.0_DP
      !
      ! ... here if |k+G|^2 <= Ecut
      !
      !
      ngk=ngk+1
      gk(ngk) = q
         !
         ! set the initial value of index array
       igk(ngk) = ng
   ENDDO
   !
   IF ( ng > ngm ) &
      CALL infomsg( 'gk_sort', 'unexpected exit from do-loop' )
   !
   ! ... order vector gk keeping initial position in index
   !
   CALL hpsort_eps( ngk, gk, igk, eps8 )
   !
   ! ... now order true |k+G|
   !
   DO nk = 1, ngk
      gk(nk) = SUM( (k(:) + g(:,igk(nk)) )**2 )
   ENDDO
   !
END SUBROUTINE gk_sort_ngw_lmax
