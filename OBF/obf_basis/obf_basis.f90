!-----------------------------------------------------------------------
program obf_basis
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
  use pwcom, only :  nspin
  use uspp, ONLY : okvan
  use realus, ONLY : generate_qpointlist
  USE wannier_gw, ONLY : num_nbndv 
  USE gvect, ONLY : ngm
  USE gvecs, ONLY : doublegrid
  USE constants, ONLY : BOHR_RADIUS_SI
  USE cell_base,        ONLY : alat
  USE input_basis, ONLY : num_val,  num_cond, s_bands, deallocate_simple,&
                       &allocate_simple, gsthres, &
                       !&lshrink,
                       &npw_max, calc_mode, trace_mode, &
                       &restart_obf, obf_dir, ecutobf, &
                       &no_gs
  USE pwcom, ONLY : igk_k, qk, nksinterp_shirley
  !SJ
  USE wfc_basis, ONLY : wfc_basis_hack
  USE s_basis, ONLY : s_basis_obf
  USE f_basis, ONLY : f_basis_obf
  !USE write_interp
  !
  IMPLICIT NONE
  CHARACTER(len=9) :: code = 'OBF'
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
  
  NAMELIST /inputobf/ prefix,outdir,num_nbndv,num_val,num_cond,s_bands, &
       & calc_mode, gsthres, &
       & trace_mode, restart_obf, obf_dir, ecutobf, &
       & no_gs
  !
  CALL mp_startup ( )
  CALL environment_start ( code )
  !
  CALL start_clock('obf_basis')
  !
  !SJ lsimple flag to turn on write_collected_wf_simple in punch -> pw_restart_new 
  !lsimple =.true.
  !
  !n_shrink=1
  prefix='export'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  IF ( ionode ) THEN
     CALL input_from_file ( )  !This open the input file and assign read in unit to 5
     READ(5,inputobf,IOSTAT=ios) !Read stuffs that matches NAMELIST
     IF (ios /= 0) CALL errore ('SIMPLE', 'reading inputobf namelist', ABS(ios) )
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
  CALL mp_bcast(gsthres, ionode_id, world_comm )
  !CALL mp_bcast(eigval_as_thres,ionode_id,world_comm)
  !CALL mp_bcast( lshrink, ionode_id, world_comm  ) !SJ
  CALL mp_bcast( calc_mode, ionode_id , world_comm )
  CALL mp_bcast(trace_mode, ionode_id , world_comm ) 
  CALL mp_bcast( restart_obf, ionode_id , world_comm )
  CALL mp_bcast(obf_dir, ionode_id , world_comm )
  CALL mp_bcast( no_gs, ionode_id, world_comm  )
  !
  !
  CALL read_file !Read file from scf calculation. 
  !!SPIN STUFF. newd() is called during read_file

  !WRITE(stdout,*) 'Reading in from file', kpoint_filename
  !CALL read_kpoint()
  !CALL mp_bcast( qk, ionode_id, world_comm  ) !SJ
  !CALL mp_bcast( nksinterp_shirley, ionode_id, world_comm  ) !SJ
  !
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
  !
  CALL read_export(pp_file,kunittmp,uspp_spsi, ascii, single_file, raw)
  !
  CALL summary()
  !
  CALL print_ks_energies()
  !
  !IF (lda_plus_u) THEN
  !  CALL init_ns()
  !ENDIF
  !
  !CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  !
  IF ( okvan) CALL generate_qpointlist()
  !
  CALL allocate_simple

  !Read in OBF from previous run
   IF ( restart_obf) THEN
           IF (ionode) THEN
               CALL read_obf
           ENDIF
               !Now we have obf(igwx,nbnd) in memory.
               !Need to distribute among procs
           CALL split_obf
           !STOP
           !ENDIF
   ENDIF
  !
  !
  !
  IF ( calc_mode == 0 ) THEN
          WRITE(stdout,*) 'Using Gram-Schmidt to build OBF'
   CALL wfc_basis_hack
  ELSEIF ( calc_mode == 1) THEN
          WRITE(stdout,*) 'Using diagonalized Smatrix to build OBF'
   CALL s_basis_obf
  ELSEIF ( calc_mode == 2) THEN
          WRITE(stdout,*) 'Using full Smatrix to build OBF'
   CALL f_basis_obf
  ELSE
          WRITE(stdout,*) 'calc_mode is', calc_mode
          WRITE(stdout,*) 'calc_mode 0 for Gram-Schmidt 1 for Smatrix'
          WRITE(stdout,*) 'Incorrect calc_mode. Exiting'
          STOP
  ENDIF 
  
  WRITE(stdout,*)     'wfc_basis completed'
  !
  
  !WRITE(stdout,*) 'npw_max in each proc'
  !WRITE(*,*) 'mpime ',mpime, 'npw_max ', npw_max
  !PRINT *, 'mpime: ', mpime, 'npw_max ', npw_max


  CALL write_shirley_ob() 
  !CALL mp_bcast(npw_max, ionode_id, world_comm)
  CALL deallocate_simple !Will deallocate bec_e, wfc_e in wfc_basis()
  !
  CALL stop_clock('obf_basis')
  CALL print_clock('optimal_basis')
  CALL print_clock('obf_basis')
  !
  !SJ
  WRITE(stdout,*) 'io_level',io_level
  WRITE(stdout,*) 'lsimple', lsimple
  !CALL punch('all')
  !SJ
  STOP
  CALL stop_pp
  STOP
  !
END PROGRAM obf_basis



