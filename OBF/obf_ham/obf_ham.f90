!-----------------------------------------------------------------------
program obf_ham
  !!-----------------------------------------------------------------------
  !! 
  !! input:  namelist "&inputsimple", with variables
  !!   prefix       prefix of input files saved by program pwscf
  !!   outdir       temporary directory where files resides
  !!    
  use io_files,  ONLY : prefix, tmp_dir
  use io_files,  ONLY : psfile, pseudo_dir
  use io_global, ONLY : stdout, ionode, ionode_id
  USE control_flags,        ONLY : io_level, lsimple
  USE control_flags,        ONLY : lmanual_shirley_grid
  USE mp_global, ONLY: mp_startup
  USE mp_pools, ONLY : kunit
  use mp_world, ONLY: mpime, world_comm, nproc
  USE environment,   ONLY: environment_start
  USE mp, ONLY : mp_bcast
  use ldaU, ONLY : lda_plus_u
  use scf, only : vrs, vltot, v, kedtau
  USE fft_base,             ONLY : dfftp
  use pwcom, only :  nspin, nbnd
  use uspp, ONLY : okvan
  use realus, ONLY : generate_qpointlist
  USE wannier_gw, ONLY : num_nbndv 
  USE gvect, ONLY : ngm
  USE gvecs, ONLY : doublegrid
  USE constants, ONLY : BOHR_RADIUS_SI
  USE cell_base,        ONLY : alat
  USE input_ham, ONLY : num_val,  num_cond, s_bands, deallocate_simple,&
                       &allocate_simple,kpoint_filename, &
                       &npw_max , ntot_e, &
                       &kpoint_style,kpoint_wrap, &
                       &hamilk_byfile
  USE pwcom, ONLY : igk_k, qk, nksinterp_shirley
  !SJ
  !USE wfc_basis, ONLY : wfc_basis_hack
   !USE read_kpoint
  !USE write_interp
  !
  IMPLICIT NONE
  CHARACTER(len=9) :: code = 'HAM'
  INTEGER :: ios, kunittmp
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(len=200) :: pp_file
  LOGICAL :: uspp_spsi, ascii, single_file, raw
  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=256) :: prefix_nc 
  INTEGER :: nglobals
  LOGICAL :: l_dipols
  CHARACTER(LEN=256) :: outdir
  
  NAMELIST /inputham/ prefix,outdir,num_nbndv,num_val,num_cond,s_bands, &
       & kpoint_filename,&
       & kpoint_wrap,kpoint_style,&
       & hamilk_byfile !,ntot_e
  !
  CALL mp_startup ( )
  CALL environment_start ( code )
  !
  CALL start_clock('obf_ham')
  !
  !
  prefix='export'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  IF ( ionode ) THEN
     CALL input_from_file ( )  !This open the input file and assign read in unit to 5
     READ(5,inputham,IOSTAT=ios) !Read stuffs that matches NAMELIST
     IF (ios /= 0) CALL errore ('OBF', 'reading inputham namelist', ABS(ios) )
  ENDIF
  !
  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( outdir, ionode_id, world_comm  )
  CALL mp_bcast( tmp_dir, ionode_id , world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm  )
  CALL mp_bcast( num_nbndv, ionode_id , world_comm )
  CALL mp_bcast(num_val, ionode_id, world_comm )
  CALL mp_bcast(num_cond, ionode_id, world_comm )
  CALL mp_bcast(s_bands, ionode_id, world_comm )
  !CALL mp_bcast( lshrink, ionode_id, world_comm  )
  CALL mp_bcast( hamilk_byfile, ionode_id, world_comm)

  !variable for hamiltonian building
  !CALL mp_bcast(nonlocal_commutator,  ionode_id, world_comm)
  !CALL mp_bcast( kpoint_manual, ionode_id, world_comm  ) !Always manual from file now
  CALL mp_bcast( kpoint_filename, ionode_id , world_comm ) !SJ
  CALL mp_bcast( kpoint_style, ionode_id , world_comm ) !SJ
  CALL mp_bcast( kpoint_wrap, ionode_id , world_comm ) !SJ
  !CALL mp_bcast( ntot_e, ionode_id , world_comm ) !SJ
  !
  !
  CALL read_file !Read file from scf calculation. 
  !!SPIN STUFF. newd() is called during read_file
  !
  ! Read in k-point from file
  WRITE(stdout,*) 'Reading in from file', kpoint_filename
  IF (kpoint_style == 'crystal') THEN
  CALL read_kpoint()
  ELSEIF (kpoint_style == 'automatic') THEN
  CALL kpoint_mesh()
  ELSE 
  WRITE(stdout,*) 'Incorrect kpoint style. Expects crystal or automatic. Exiting'
  STOP
  ENDIF



  CALL mp_bcast( qk, ionode_id, world_comm  ) !SJ
  CALL mp_bcast( nksinterp_shirley, ionode_id, world_comm  ) !SJ
  
  
  CALL openfile_school
  !
#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif
  !
  pp_file= ' '
  uspp_spsi = .FALSE.
  ascii = .FALSE.
  single_file = .FALSE.
  raw = .FALSE.
  ntot_e = nbnd
  !
  !CALL read_export(pp_file,kunittmp,uspp_spsi, ascii, single_file, raw)
  !
  CALL summary()
  !
  CALL print_ks_energies()
  !
  IF (lda_plus_u) THEN
    CALL init_ns()
  ENDIF
  !
  WRITE(stdout,*) 'nbnd from data xml',nbnd
  WRITE(stdout,*) 'ntot_e',ntot_e
  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  !
  IF ( okvan) CALL generate_qpointlist()
  !
  !CALL allocate_simple
  !
  !CALL wfc_basis_hack
  
  !WRITE(stdout,*)     'wfc_basis completed'
  !
  IF ( hamilk_byfile ) THEN
          CALL khamiltonian_byfile
  ELSE
          CALL khamiltonian
  ENDIF

  IF ( hamilk_byfile) THEN
          CALL khamiltonian_byfile_spin
  ELSE
          CALL khamiltonian_spin
  ENDIF
  
  !WRITE(stdout,*) 'npw_max in each proc'
  !WRITE(*,*) 'mpime ',mpime, 'npw_max ', npw_max
  !PRINT *, 'mpime: ', mpime, 'npw_max ', npw_max


  !CALL write_shirley_ob() 
  !CALL mp_bcast(npw_max, ionode_id, world_comm)
  !CALL deallocate_simple !Will deallocate bec_e, wfc_e in wfc_basis()
  !
  CALL stop_clock('obf_ham')
  CALL print_clock('obf_ham')
  !
  STOP
  CALL stop_pp
  STOP
  !
END PROGRAM obf_ham



