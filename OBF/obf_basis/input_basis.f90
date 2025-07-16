 MODULE input_basis

  USE kinds, ONLY : DP
  USE becmod, ONLY : bec_type,deallocate_bec_type
  USE klist, ONLY : nks

  SAVE

  INTEGER :: calc_mode=0  ! calculation mode: 0 = Gram-schmidt, 1 = Smatrix 
  INTEGER :: trace_mode=0 ! Method for keeping OBF in Smatrix calc_mode
                          ! 0 use trace difference nbnd(ntot) - trace
                          ! 1 use current trace (Overlap matrix trace truncation method similar to Prendergast)
  LOGICAL :: restart_obf =.false. !Will try to read OBF from previous run
  CHARACTER (len=256) :: obf_dir !for reading obf from previous 

  INTEGER :: num_val !valence bands to be included starting from HOMO 
  INTEGER :: num_cond !conduction bands to be included starting from LUMO
  LOGICAL :: no_gs = .false. !Whether to do GS at the end after s_basis  
  REAL(kind=DP) :: ecutobf=0.d0 ! New energy cutoff for OBF for restarting new run
  !INTEGER :: npw_obf=0 ! New total number of planewave from prevous run
  REAL(kind=DP) :: s_bands=1.d-9 !threshold for Shirley's algorithm
  REAL(kind=DP) :: gsthres=1.d-4 !threshold for GS routine after Smatrix. GS is used to make sure that only orthonormal basis is
                           ! added to current OBF.
  
  COMPLEX(kind=DP), POINTER :: wfc_e(:,:),obf_e(:,:)!common basis for wfcs at every k point
  INTEGER :: ntot_e !their total number of final basis
  INTEGER :: npw_max !max number of differ G for that  MPI task considering all k-points

  COMPLEX(kind=DP), POINTER :: vkb_max(:,:) !projectors for US peudos in commom G ordering
  TYPE(bec_type), POINTER ::  bec_e(:) !<beta_k|wfc_e>


  !COMPLEX(kind=DP), POINTER :: prod_e(:,:)!common basis for wfc products
  !INTEGER :: nprod_e!their number
  !REAL(kind=DP) :: s_product!treshold for products

  !LOGICAL :: l_truncated_coulomb=.false.   !if true truncate the Coulomb potential
  !REAL(kind=DP) :: truncation_radius    !in Bohr

  !INTEGER :: nkpoints(3)  !k-points grid
  !LOGICAL :: nonlocal_commutator=.true.
  !INTEGER :: interp_npw = -1    ! number of plane waves used for the Shirley interpolation
   
  INTEGER :: npw_gamma !New number of pw for OBF after expanding ecut. 
  INTEGER :: numpw !dimension of polarizability basis of GWW, if == 0 do not call routines, it adds automatically +1 for extended systems (from pw4gww)

  !LOGICAL :: l_debug  !options for using plane waves instead of KS states for debugging purposes only
  !INTEGER :: n_debug  !number of G waves along each cartesian direction for debug 

  !INTEGER :: w_type  !approximation used for W_c: 0 from GWW extrapolated, 1 diagonal model function for screening
  !REAL(kind=DP) :: epsm  !parameter eps_m for diagonal model of epsilon (dielectric constant, eps_{\infty})
  !REAL(kind=DP) :: lambdam  !parameter lambda for diagonal model of epsilon

  !INTEGER :: n_shrink!for using only a fraction of the k_points mesh 

  !SJ
  !LOGICAL :: kpoint_manual =.false. !for using kpoint file
  !LOGICAL :: lshrink = .false. !for using only input nbnd for output shirley wfc.
  !LOGICAL :: eigval_as_thres =.false. ! If .true. then the smallest eigenvalue from diagonalized Smatrix that
                                      ! is still within the trace tolerance is used as the GS tolerance
                                      ! in the GS routine after instead of specified gsthres
  !CHARACTER (len=256) :: kpoint_filename !for kpoint filename to read

CONTAINS

  subroutine deallocate_simple
    USE uspp, ONLY : okvan

    implicit none
    
    INTEGER :: i
    
    if(okvan) then
       do i=1,nks
          call deallocate_bec_type(bec_e(i))
       enddo
    endif
    if(associated(wfc_e)) deallocate(wfc_e)
    if(associated(bec_e)) deallocate(bec_e)

    return
  end subroutine deallocate_simple

  subroutine allocate_simple

    implicit none

    allocate(bec_e(nks))

  end  subroutine allocate_simple

END MODULE input_basis
