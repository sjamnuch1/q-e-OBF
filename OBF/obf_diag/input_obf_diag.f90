MODULE input_obf_diag
!this module provides input file routines 

  USE kinds, ONLY: DP

   TYPE input_options_obf_diag
      CHARACTER(len=256) :: prefix = 'prefix'!prefix to designate the files same as in PW        
      CHARACTER(len=256) :: outdir = './'!outdir to designate the files same as in PW   
      CHARACTER(len=256) :: unkdir = 'Out_unk' !Save directory for interpolated wfc files
      INTEGER :: h_level = 2 ! terms included in the k-hamiltonian (for debug) 0=kinetic, 1=kinetic+local >1=kinetic+local+non-local
      INTEGER :: interp_grid(3)    !  k-grid used for the calculation of the optical spectra (Shirley interpolation)
      REAL(DP) :: kshift(3)       !  k shift offset when building interpolated kmesh.
    !  LOGICAL :: nonlocal_commutator=.true.  ! if true it includes the non-local part of [r,H] in the calculation of the optical spectra
      LOGICAL :: nonlocal_interpolation = .false.   ! if true it interpolates the non-local part of the psp  
                                                    !based on the interp_grid and kshift (WARNING experimental feature)
                                                   
      LOGICAL :: tetrahedron_method = .false.    ! Use tetrahedron method to calculate Drude plasma frequency, Fermi energy and DOS
      REAL(kind=DP) :: fermi_degauss=0.02205   ! degauss (in Ry) for Fermi level calculation
      REAL(kind=DP) :: fermi_energy = -1   ! input Fermi energy (in Ry) (if not given in input it is calculated using the Shirley interpolated bands)
      INTEGER  :: fermi_ngauss=-99         ! ngauss for Fermi level calculation (n=-99 --> Fermi-Dirac ; n=-1 --> cold-smearing ; n>=0 --> Methfessel-Paxton)
      REAL(kind=DP) ::  delta_energy_dos = 0.000735  ! Energy spacing for DOS calculation (in Ry) (default is 0.01 eV)
      !REAL(kind=DP) :: drude_degauss   ! degauss (in Ry) for Drude plasma frequency calculation
      INTEGER  :: drude_ngauss = -99        ! ngauss for Drude plasma frequency calculation (n=-99 --> Fermi-Dirac) (WARNING: only Fermi-Dirac smearing is implemented now)
      REAL(kind=DP) :: elec_temp = 0.0018375      ! electronic temperature (in Ry) for the occupation of the electronic states. Default is room temperature
      REAL(kind=DP) ::  inter_broadening=0.0142   ! broadening used for the optical spectra (interband part)
      !INTEGER :: nw                 ! Number of points in the energy interval [wmin,wmax]
      !SJ 
      LOGICAL :: hamilk_byfile = .false. ! If true Hamiltonian_k files are individually writte. 
                                         ! Each proc would then try to read corresponding hamil files 
                                         ! According to the stored global kpoint index.
      LOGICAL :: lcollect_bnk = .false. !for saving diagonalized Hamiltonian matrix.
      LOGICAL :: lshrink = .false. !for using nbnd only diagonalization.
      INTEGER :: nband_save = 0 !for lshrink to save up 1 -> nband_save 
                                ! If == 0 then the saved OBF will be the size of DFT bands.
      LOGICAL :: build_unk =.false. !for building unk from bnk and wfc from simple.x
                                      !each block contains nkpoints(3) entry
      CHARACTER (len=256) :: kpoint_style ='crystal' !for kpoint style. either 'crystal' or 'automatic'
                                                     ! if crystal then it expects in kpoint_filename 
                                                     !to have a number of kpoints and list of the kpoint in crystl coord
                                                     ! if automatic then it follows same QE format nk1 nk2 nk3 k1 k2 k3
                                                     ! Note that k can be any real number and the offset is k/nk.
                                                     ! In QE it's k/(2*nk). 
      LOGICAL :: kpoint_wrap =.false. ! If true, the Hamiltonian is built with [0,1] BZ. Having NK+1 point.
                                      ! Preserve PBC better for OBF.
      CHARACTER (len=256) :: kpoint_filename ='KPOINTS' !for kpoint filename to read
      CHARACTER (len=256) :: interp_kpoint_filename ='KPOINTS_INTERP' !for kpoint filename to read for interpolate grid.
      CHARACTER (len=256) :: interp_style ='automatic' ! Whether to use a list of kpoint or uniform mesh.
                                                       ! crystal will read the list while automatic will use
                                                       ! interp_grid(3) and kshift(3) to build the mesh.      
   END TYPE input_options_obf_diag

   CONTAINS

     SUBROUTINE  read_input_obf_diag( obfdiag_in )
       USE io_global,            ONLY : stdout, ionode, ionode_id
       USE mp,                   ONLY : mp_bcast
       USE mp_world,             ONLY : world_comm
       USE io_files,             ONLY : tmp_dir, prefix
       

       implicit none

       CHARACTER(LEN=256), EXTERNAL :: trimcheck
       TYPE(input_options_obf_diag) :: obfdiag_in          !in output the input parameters

       NAMELIST/inputobfdiag/obfdiag_in

	CHARACTER(LEN=256) :: outdir

       if(ionode) then
          read(*,NML=inputobfdiag)
          outdir = trimcheck(obfdiag_in%outdir)
          tmp_dir = outdir
          prefix = trim(obfdiag_in%prefix)
       endif

       call mp_bcast( outdir,ionode_id, world_comm )
       call mp_bcast( obfdiag_in%unkdir,ionode_id, world_comm )
       call mp_bcast( tmp_dir,ionode_id, world_comm )
       call mp_bcast( prefix,ionode_id, world_comm )
       call mp_bcast( obfdiag_in%prefix,ionode_id,world_comm)
       call mp_bcast( obfdiag_in%interp_grid, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%kshift, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%h_level, ionode_id, world_comm)
       !call mp_bcast( obfdiag_in%nonlocal_commutator, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%nonlocal_interpolation, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%fermi_degauss, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%fermi_energy, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%fermi_ngauss, ionode_id, world_comm)
       !call mp_bcast( obfdiag_in%drude_degauss, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%drude_ngauss, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%elec_temp, ionode_id, world_comm)
       !call mp_bcast( obfdiag_in%wmin, ionode_id, world_comm)
       !call mp_bcast( obfdiag_in%wmax, ionode_id, world_comm)
       !call mp_bcast( obfdiag_in%nw, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%inter_broadening, ionode_id, world_comm)
       !call mp_bcast( obfdiag_in%intra_broadening, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%tetrahedron_method, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%delta_energy_dos, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%kpoint_style, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%kpoint_filename, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%interp_style, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%interp_kpoint_filename, ionode_id, world_comm)      
       call mp_bcast( obfdiag_in%kpoint_wrap, ionode_id, world_comm) 
       call mp_bcast( obfdiag_in%hamilk_byfile, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%lshrink, ionode_id, world_comm)
       call mp_bcast( obfdiag_in%nband_save,ionode_id,world_comm)
       call mp_bcast( obfdiag_in%lcollect_bnk,ionode_id, world_comm)
       call mp_bcast( obfdiag_in%build_unk,ionode_id, world_comm)
       return

     END SUBROUTINE read_input_obf_diag

END MODULE input_obf_diag

