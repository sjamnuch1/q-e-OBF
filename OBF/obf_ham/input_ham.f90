 MODULE input_ham

  USE kinds, ONLY : DP
  USE becmod, ONLY : bec_type,deallocate_bec_type
  USE klist, ONLY : nks

  SAVE

  !INTEGER :: calc_mode  ! calculation mode: 0 = BSE calculation, 1 = IP calculation

  INTEGER :: num_val !valence bands to be included starting from HOMO 
  INTEGER :: num_cond !conduction bands to be included starting from LUMO  
  REAL(kind=DP) :: s_bands=0.01 !threshold for Shirley's algorithm (Not really used. Need to trim this out later)
  
  !COMPLEX(kind=DP), POINTER :: wfc_e(:,:)!common basis for wfcs at every k point
  INTEGER :: ntot_e !their total number from OBF part. 
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

  !INTEGER :: numpw !dimension of polarizability basis of GWW, if == 0 do not call routines, it adds automatically +1 for extended systems (from pw4gww)

  !LOGICAL :: l_debug  !options for using plane waves instead of KS states for debugging purposes only
  !INTEGER :: n_debug  !number of G waves along each cartesian direction for debug 

  !INTEGER :: w_type  !approximation used for W_c: 0 from GWW extrapolated, 1 diagonal model function for screening
  !REAL(kind=DP) :: epsm  !parameter eps_m for diagonal model of epsilon (dielectric constant, eps_{\infty})
  !REAL(kind=DP) :: lambdam  !parameter lambda for diagonal model of epsilon

  !INTEGER :: n_shrink!for using only a fraction of the k_points mesh 

  !SJ
  !LOGICAL :: kpoint_manual =.false. !for using kpoint file
  LOGICAL :: hamilk_byfile = .false. ! If true, the Hamiltonian_k files are written into multiple
                                     ! individual files instead of one large file.
  !LOGICAL :: lshrink = .false. !for using only input nbnd for output shirley wfc.
  CHARACTER (len=256) :: kpoint_style = 'crystal' !Style for kpoint. 
                                               !'crystal' expects number of kpoint. 
                                               !Followed by list of kpoints and weight.
                                               !'automatic' expects nk1 nk2 nk3 k1 k2 k3. Same as QE format.
                                               ! nk => number of kpoint. k => offset.
  LOGICAL :: kpoint_wrap  = .false. ! If true will add additional kpoint at 1.0 making it nk+1 kpoints total.
                                    ! Useful when building periodic OBF.
  CHARACTER (len=256) :: kpoint_filename = 'KPOINTS' !for kpoint filename to read
  INTEGER :: nk(3) !Mesh for kpoint 
  REAL(DP) :: kshift(3) !and offset

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
   ! if(associated(wfc_e)) deallocate(wfc_e)
    if(associated(bec_e)) deallocate(bec_e)

    return
  end subroutine deallocate_simple

  subroutine allocate_simple

    implicit none

    allocate(bec_e(nks))

  end  subroutine allocate_simple

END MODULE input_ham
