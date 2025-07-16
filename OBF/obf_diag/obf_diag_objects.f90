MODULE obf_diag_objects
  ! this module describes the most important objects
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  !
  TYPE energies
     INTEGER :: nk!total number of k points
     INTEGER :: nk_loc!local number of k points
     INTEGER :: ik_first!first local k point
     INTEGER :: ik_last!last local k point  
     INTEGER :: num_bands   ! number of states
     REAL(kind=DP), DIMENSION(:,:), POINTER :: energy ! energies (num_bands,nk_loc)
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: energy_der ! derivatives of the energy  (3,num_bands,nk_loc)
  END TYPE
  !
  TYPE kpoints
     INTEGER :: nkgrid(3)         !  k-grid from obf_ham / for interpolation (interp_grid) 
     REAL(kind=DP), DIMENSION(:,:), POINTER :: qk ! kpoints coordinates (3,nkgrid(1)*nkgrid(2)*nkgrid(3))
     REAL(kind=DP) :: alat        ! lattice paramater 
     REAL(kind=DP) :: bg(3,3)     ! reciprocal basis vectors
     INTEGER :: nk_loc            ! local number of k points
     INTEGER :: ik_first          ! first local k point
     INTEGER :: ik_last           ! last local k point  
     INTEGER :: nk_smooth_loc     ! number of local k-points
     INTEGER :: nk                ! total number of k-points
     !INTEGER :: ncoarse(3)        ! k-grid from obf_ham
     REAL(kind=DP), DIMENSION(:,:), POINTER  ::  pos_cube   !   (3,nk_loc) 
     INTEGER, DIMENSION(:,:), POINTER :: coord_cube   ! (8,nk_loc)
     REAL(DP) :: ks(3)            ! k-grid offset from obf_ham
  END TYPE
  !
  TYPE shirley
      INTEGER :: ntot_e ! number of optimal basis states SJ:new nbnd after simple.x
      LOGICAL :: noncolin ! if it is non-collinear
      INTEGER :: nat ! number of atoms
      INTEGER :: ntyp ! number of types of atoms
      INTEGER :: nhm  ! max number of different beta functions per atom
      INTEGER :: nspin  ! number of spins (1=no spin, 2=LSDA)
      INTEGER :: nkb  ! total number of beta functions
      INTEGER :: npol  !
      INTEGER :: nks ! number of k-points of the smooth grid
      INTEGER :: nk_smooth_loc ! number of local k-points
      INTEGER :: num_val , num_cond  ! number of interpolated valence bands, number of interpolated conduction bands
      INTEGER :: num_nbndv(2)       ! total number of occupied bands
      INTEGER :: num_bands       ! total number of considered bands (num_bands = num_val + num_cond)
      !LOGICAL :: nonlocal_commutator  !
      INTEGER :: nkpoints
      !INTEGER, DIMENSION(3) :: nkpoints     ! smooth k-points grid on which H(k) is calculated
      INTEGER :: nkp_tot                    !SJ Total number of kpoint. nksinterp_shirley from simple_hack
      INTEGER, DIMENSION(:), POINTER :: ityp ! (nat)
      INTEGER, DIMENSION(:), POINTER :: nh   ! (ntyp)
      INTEGER, DIMENSION(:), POINTER :: ofsbeta ! (nat)
      REAL(kind=8) :: alat, nelec, omega        ! lattice paramater, number of electrons, volume prim. cell
      REAL(kind=8) :: bg(3,3)     ! reciprocal basis vectors
      REAL(kind=8) :: at(3,3)     ! direct basis vectors
      REAL(kind=DP) :: s_bands     ! threshold for the construction of the optimal basis
      COMPLEX(kind=DP), DIMENSION(:,:), POINTER ::  h0 , Vloc ! k-dependent Hamiltonian  (ntot_e,ntot_e)
      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: h1 ! (ntot_e,ntot_e,3)
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: deeq_nc  ! (nhm,nhm,nat,nspin)
      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: deeqc !(nhm,nhm,nat)
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: beck_nc ! (nkb,npol,ntot_e,nk)
      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: beckc ! (nkb,ntot_e,nk)
      !COMPLEX(kind=DP),DIMENSION(:,:,:,:), POINTER :: commut ! (ntot_e,ntot_e,3,nk)
  END TYPE
  !
  TYPE eigen
    INTEGER :: num_bands       ! number of states
    REAL(kind=DP) :: q(3)        ! k-point
    REAL(kind=DP), DIMENSION(:), POINTER  :: energy     ! eigenenergies (num_bands)
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER  :: wave_func      ! eigenfunctions (ntot_e,num_bands)
  END TYPE

   CONTAINS

   subroutine initialize_energies(element)
      implicit none
      TYPE(energies) :: element
      nullify(element%energy)
      nullify(element%energy_der)
      return
    end subroutine initialize_energies
    
    subroutine deallocate_energies(element)
      implicit none
      TYPE(energies) :: element
      if(associated(element%energy)) deallocate(element%energy)
      nullify(element%energy)
      if(associated(element%energy_der)) deallocate(element%energy_der)
      nullify(element%energy_der)
      return
    end subroutine deallocate_energies

 
   subroutine initialize_kpoints(element)
      implicit none
      TYPE(kpoints) :: element
      nullify(element%qk)
      nullify(element%pos_cube)
      nullify(element%coord_cube)
      return
    end subroutine initialize_kpoints
    
    SUBROUTINE deallocate_kpoints(element)
      implicit none
      TYPE(kpoints) :: element
      if(associated(element%qk)) deallocate(element%qk)
      nullify(element%qk)
      if(associated(element%pos_cube)) deallocate(element%pos_cube)
      nullify(element%pos_cube)
      if(associated(element%coord_cube)) deallocate(element%coord_cube)
      nullify(element%coord_cube)
      return
    END SUBROUTINE deallocate_kpoints

    SUBROUTINE initialize_shirley(element)
      implicit none
      TYPE(shirley) :: element
      nullify(element%ityp)
      nullify(element%nh)
      nullify(element%ofsbeta)
      nullify(element%h0)
      nullify(element%h1)
      nullify(element%Vloc)
      nullify(element%deeq_nc)
      nullify(element%deeqc)
      nullify(element%beck_nc)
      nullify(element%beckc)
      !nullify(element%commut)
      return
    END SUBROUTINE initialize_shirley
    
    SUBROUTINE deallocate_shirley(element)
      implicit none
      TYPE(shirley) :: element
      if(associated(element%ityp)) deallocate(element%ityp)
      nullify(element%ityp)
      if(associated(element%nh)) deallocate(element%nh)
      nullify(element%nh)
      if(associated(element%ofsbeta)) deallocate(element%ofsbeta)
      nullify(element%ofsbeta)
      if(associated(element%h0)) deallocate(element%h0)
      nullify(element%h0)
      if(associated(element%h1)) deallocate(element%h1)
      nullify(element%h1)
      if(associated(element%Vloc)) deallocate(element%Vloc)
      nullify(element%Vloc)
      if(associated(element%deeq_nc)) deallocate(element%deeq_nc)
      nullify(element%deeq_nc)
      if(associated(element%deeqc)) deallocate(element%deeqc)
      nullify(element%deeqc)
      if(associated(element%beck_nc)) deallocate(element%beck_nc)
      nullify(element%beck_nc)
      if(associated(element%beckc)) deallocate(element%beckc)
      nullify(element%beckc)
      !if(associated(element%commut)) deallocate(element%commut)
      !nullify(element%commut)
      return
    END SUBROUTINE deallocate_shirley

   SUBROUTINE initialize_eigen(element)
      implicit none
      TYPE(eigen) :: element
      nullify(element%energy)
      nullify(element%wave_func)
      return
    END SUBROUTINE initialize_eigen
    
    SUBROUTINE deallocate_eigen(element)
      implicit none
      TYPE(eigen) :: element
      if(associated(element%energy)) deallocate(element%energy)
      nullify(element%energy)
      if(associated(element%wave_func)) deallocate(element%wave_func)
      nullify(element%wave_func)
      return
    END SUBROUTINE deallocate_eigen



    SUBROUTINE read_shirley(obfdiag_in,sh) !SJ Todo change the input to (simpleip_in,sh,up/down) then point to read file 
      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout

      implicit none

      TYPE(input_options_obf_diag) :: obfdiag_in
      TYPE(shirley) :: sh
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun, idir, nk
      write(stdout,*)'obf_diag: opening file "hamiltonian"'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian', status='old',form='unformatted')
         write(stdout,*)'File opened'
         read(iun) sh%ntot_e
         write(stdout,*) 'Basis/nbnd from obf_basis.x', sh%ntot_e
      endif   
      call mp_bcast(sh%ntot_e,ionode_id,world_comm)
      allocate(sh%h0(sh%ntot_e,sh%ntot_e), sh%h1(sh%ntot_e,sh%ntot_e,3), sh%Vloc(sh%ntot_e,sh%ntot_e))

      if (ionode) read(iun) sh%h0(1:sh%ntot_e,1:sh%ntot_e) 

      do idir = 1,3
        if (ionode) read(iun) sh%h1(1:sh%ntot_e,1:sh%ntot_e,idir) 
      enddo

      if (ionode) read(iun) sh%Vloc(1:sh%ntot_e,1:sh%ntot_e) 
      !if (ionode) write(stdout,*) 'Vloc(1,1)',sh%Vloc(1,1)

      call mp_bcast(sh%h0,ionode_id,world_comm)
      call mp_bcast(sh%h1,ionode_id,world_comm)
      call mp_bcast(sh%Vloc,ionode_id,world_comm)
      
      if(ionode) then
        read(iun) sh%noncolin
        read(iun) sh%nat 
        read(iun) sh%ntyp 
        read(iun) sh%nhm
        read(iun) sh%nspin
        read(iun) sh%nkb 
        read(iun) sh%npol
        write(stdout,*)  "sh%noncolin", sh%noncolin
        write(stdout,*)  "sh%nat", sh%nat
        write(stdout,*)  "sh%ntyp", sh%ntyp
        write(stdout,*)  "sh%nhm", sh%nhm
        write(stdout,*)  "sh%nspin", sh%nspin
        write(stdout,*)  "sh%nkb", sh%nkb
        write(stdout,*)  "sh%npol", sh%npol
      endif
      call mp_bcast(sh%noncolin,ionode_id,world_comm)
      call mp_bcast(sh%nat,ionode_id,world_comm)
      call mp_bcast(sh%ntyp,ionode_id,world_comm)
      call mp_bcast(sh%nhm,ionode_id,world_comm)
      call mp_bcast(sh%nspin,ionode_id,world_comm)
      call mp_bcast(sh%nkb,ionode_id,world_comm)
      call mp_bcast(sh%npol,ionode_id,world_comm)
      

      allocate(sh%ityp(sh%nat), sh%nh(sh%ntyp), sh%ofsbeta(sh%nat))
      if(ionode) then
        read(iun) sh%ityp(1:sh%nat)
        read(iun) sh%nh(1:sh%ntyp)
        read(iun) sh%ofsbeta(1:sh%nat)
        read(iun) sh%nkpoints
        read(iun) sh%nkp_tot
      endif
      call mp_bcast(sh%ityp,ionode_id,world_comm)
      call mp_bcast(sh%nh,ionode_id,world_comm)
      call mp_bcast(sh%ofsbeta,ionode_id,world_comm)
      call mp_bcast(sh%nkpoints,ionode_id,world_comm)
      call mp_bcast(sh%nkp_tot,ionode_id,world_comm)
      
      !nk = (sh%nkpoints(1))*(sh%nkpoints(2))*(sh%nkpoints(3))

      allocate(sh%deeqc(sh%nhm,sh%nhm,sh%nat), sh%deeq_nc(sh%nhm,sh%nhm,sh%nat,sh%nspin))
      if (sh%noncolin) then
         if (ionode) read(iun) sh%deeq_nc(1:sh%nhm,1:sh%nhm,1:sh%nat,1:sh%nspin)
      else
         if (ionode) read(iun) sh%deeqc(1:sh%nhm,1:sh%nhm,1:sh%nat)
      endif
      
      if (sh%noncolin) then
         call mp_bcast(sh%deeq_nc,ionode_id,world_comm)
      else
         call mp_bcast(sh%deeqc,ionode_id,world_comm)
      endif
            
      if(ionode) then
         read(iun) sh%alat
         read(iun) sh%bg(1:3,1:3)
         WRITE(stdout,*) 'bg'
         WRITE(stdout,*) sh%bg(1:3,1:3)
         read(iun) sh%at(1:3,1:3)
         read(iun) sh%nelec
         read(iun) sh%omega
         read(iun) sh%num_val
         read(iun) sh%num_cond
         read(iun) sh%num_nbndv
         !read(iun) sh%nonlocal_commutator
         read(iun) sh%s_bands
      endif
      !sh%num_bands = sh%num_val + sh%num_cond
      sh%num_bands = sh%ntot_e !SJ
      if (obfdiag_in%lshrink) THEN  !SJ
              if (obfdiag_in%nband_save  < 1 ) THEN
              sh%num_bands = sh%num_val + sh%num_cond !SJ
              else
              write(stdout,*) 'Using lshrink, will save the OBF to size of',obfdiag_in%nband_save
              sh%num_bands = obfdiag_in%nband_save
              endif
              if ( sh%num_bands > sh%ntot_e) then
                      write(stdout,*) 'Requested number of bands is larger than the number of bands in OBF'
                      write(stdout,*) 'Adjusting it to maxmimum size of OBF'
                      sh%num_bands = sh%ntot_e
              endif
              
      endif !SJ
      call mp_bcast(sh%alat,ionode_id,world_comm)
      call mp_bcast(sh%bg,ionode_id,world_comm)
      call mp_bcast(sh%at,ionode_id,world_comm)
      call mp_bcast(sh%nelec,ionode_id,world_comm)
      call mp_bcast(sh%omega,ionode_id,world_comm)
      call mp_bcast(sh%num_val,ionode_id,world_comm)
      call mp_bcast(sh%num_cond,ionode_id,world_comm)
      call mp_bcast(sh%num_nbndv,ionode_id,world_comm)
      call mp_bcast(sh%num_bands,ionode_id,world_comm) 
      !call mp_bcast(sh%nonlocal_commutator,ionode_id,world_comm)
      call mp_bcast(sh%s_bands,ionode_id,world_comm)
      !
      if(ionode) then
         close(iun)
         write(stdout,*)'File closed'
      endif
      !
      !if ( .not. simpleip_in%nonlocal_interpolation .and. simpleip_in%kpoint_manual ) then
      !  write(stdout,*) 'Using no nonlocal interpolation and Manual kpoints' 
      !elseif (.not. simpleip_in%nonlocal_interpolation .and. ( sh%nkpoints(1) /= simpleip_in%interp_grid(1) .or. &
      !& sh%nkpoints(2) /= simpleip_in%interp_grid(2) .or. sh%nkpoints(3) /= simpleip_in%interp_grid(3) ) ) then
      !  call errore('SIMPLE_IP', 'W/o trilinear interpolation, k-grids from simple and simple_ip have to be equal',1)
      !elseif ( simpleip_in%nonlocal_interpolation .and. ( 2*sh%nkpoints(1) /= simpleip_in%interp_grid(1) .or. &
      !& 2*sh%nkpoints(2) /= simpleip_in%interp_grid(2) .or. 2*sh%nkpoints(3) /= simpleip_in%interp_grid(3) ) ) then
      !  call errore('SIMPLE_IP', 'W/ trilinear interpolation, simple_ip k-grid has to be the double of the simple k-grid',1)
      !endif

      !if( (sh%num_val .ne. sh%num_nbndv(1)) ) then
      !  call errore('OBF_DIAG', 'All the occupied bands must be included for Shirley interpolation (num_val=num_nbndv)',1)
      !endif
      !
    end subroutine read_shirley

        SUBROUTINE read_shirley_spin(obfdiag_in,sh)
      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout

      implicit none

      TYPE(input_options_obf_diag) :: obfdiag_in
      TYPE(shirley) :: sh
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun, idir, nk
      write(stdout,*)'obf_diag: opening file "hamiltonian"'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian_down', status='old',form='unformatted')
         write(stdout,*)'File opened'
         read(iun) sh%ntot_e
         write(stdout,*) 'Basis/nbnd from obf_basis.x', sh%ntot_e
      endif
      call mp_bcast(sh%ntot_e,ionode_id,world_comm)
      allocate(sh%h0(sh%ntot_e,sh%ntot_e), sh%h1(sh%ntot_e,sh%ntot_e,3), sh%Vloc(sh%ntot_e,sh%ntot_e))

      if (ionode) read(iun) sh%h0(1:sh%ntot_e,1:sh%ntot_e)

      do idir = 1,3
        if (ionode) read(iun) sh%h1(1:sh%ntot_e,1:sh%ntot_e,idir)
      enddo

      if (ionode) read(iun) sh%Vloc(1:sh%ntot_e,1:sh%ntot_e)

      call mp_bcast(sh%h0,ionode_id,world_comm)
      call mp_bcast(sh%h1,ionode_id,world_comm)
      call mp_bcast(sh%Vloc,ionode_id,world_comm)

      if(ionode) then
        read(iun) sh%noncolin
        read(iun) sh%nat
        read(iun) sh%ntyp
        read(iun) sh%nhm
        read(iun) sh%nspin
        read(iun) sh%nkb
        read(iun) sh%npol
        write(stdout,*)  "sh%noncolin", sh%noncolin
        write(stdout,*)  "sh%nat", sh%nat
        write(stdout,*)  "sh%ntyp", sh%ntyp
        write(stdout,*)  "sh%nhm", sh%nhm
        write(stdout,*)  "sh%nspin", sh%nspin
        write(stdout,*)  "sh%nkb", sh%nkb
        write(stdout,*)  "sh%npol", sh%npol
      endif
      call mp_bcast(sh%noncolin,ionode_id,world_comm)
      call mp_bcast(sh%nat,ionode_id,world_comm)
      call mp_bcast(sh%ntyp,ionode_id,world_comm)
      call mp_bcast(sh%nhm,ionode_id,world_comm)
      call mp_bcast(sh%nspin,ionode_id,world_comm)
      call mp_bcast(sh%nkb,ionode_id,world_comm)
      call mp_bcast(sh%npol,ionode_id,world_comm)


      allocate(sh%ityp(sh%nat), sh%nh(sh%ntyp), sh%ofsbeta(sh%nat))
      if(ionode) then
        read(iun) sh%ityp(1:sh%nat)
        read(iun) sh%nh(1:sh%ntyp)
        read(iun) sh%ofsbeta(1:sh%nat)
        read(iun) sh%nkpoints
        read(iun) sh%nkp_tot
      endif
      call mp_bcast(sh%ityp,ionode_id,world_comm)
      call mp_bcast(sh%nh,ionode_id,world_comm)
      call mp_bcast(sh%ofsbeta,ionode_id,world_comm)
      call mp_bcast(sh%nkpoints,ionode_id,world_comm)
      call mp_bcast(sh%nkp_tot,ionode_id,world_comm)

      !nk = (sh%nkpoints(1))*(sh%nkpoints(2))*(sh%nkpoints(3))

      allocate(sh%deeqc(sh%nhm,sh%nhm,sh%nat), sh%deeq_nc(sh%nhm,sh%nhm,sh%nat,sh%nspin))
      if (sh%noncolin) then
         if (ionode) read(iun) sh%deeq_nc(1:sh%nhm,1:sh%nhm,1:sh%nat,1:sh%nspin)
      else
         if (ionode) read(iun) sh%deeqc(1:sh%nhm,1:sh%nhm,1:sh%nat)
      endif

      if (sh%noncolin) then
         call mp_bcast(sh%deeq_nc,ionode_id,world_comm)
      else
         call mp_bcast(sh%deeqc,ionode_id,world_comm)
      endif

      if(ionode) then
         read(iun) sh%alat
         read(iun) sh%bg(1:3,1:3)
         WRITE(stdout,*) 'bg'
         WRITE(stdout,*) sh%bg(1:3,1:3)
         read(iun) sh%at(1:3,1:3)
         read(iun) sh%nelec
         read(iun) sh%omega
         read(iun) sh%num_val
         read(iun) sh%num_cond
         read(iun) sh%num_nbndv
         !read(iun) sh%nonlocal_commutator
         read(iun) sh%s_bands
      endif
      sh%num_bands = sh%ntot_e !SJ
      if (obfdiag_in%lshrink) THEN  !SJ
              if (obfdiag_in%nband_save  < 1) THEN
              sh%num_bands = sh%num_val + sh%num_cond !SJ
              else
              write(stdout,*) 'Using lshrink, will save the OBF to size of',obfdiag_in%nband_save
              sh%num_bands = obfdiag_in%nband_save
              endif
              if ( sh%num_bands > sh%ntot_e) then
                      write(stdout,*) 'Requested number of bands is larger than the number of bands in OBF'
                      write(stdout,*) 'Adjusting it to maxmimum size of OBF'
                      sh%num_bands = sh%ntot_e
              endif

      endif !SJ
      call mp_bcast(sh%alat,ionode_id,world_comm)
      call mp_bcast(sh%bg,ionode_id,world_comm)
      call mp_bcast(sh%at,ionode_id,world_comm)
      call mp_bcast(sh%nelec,ionode_id,world_comm)
      call mp_bcast(sh%omega,ionode_id,world_comm)
      call mp_bcast(sh%num_val,ionode_id,world_comm)
      call mp_bcast(sh%num_cond,ionode_id,world_comm)
      call mp_bcast(sh%num_nbndv,ionode_id,world_comm)
      call mp_bcast(sh%num_bands,ionode_id,world_comm)
      !call mp_bcast(sh%nonlocal_commutator,ionode_id,world_comm)
      call mp_bcast(sh%s_bands,ionode_id,world_comm)
      !
      if(ionode) then
         close(iun)
         write(stdout,*)'File closed'
      endif
      !
      !if ( .not. simpleip_in%nonlocal_interpolation .and. simpleip_in%kpoint_manual ) then
      !  write(stdout,*) 'Using no nonlocal interpolation and Manual kpoints'
      !elseif (.not. simpleip_in%nonlocal_interpolation .and. ( sh%nkpoints(1) /= simpleip_in%interp_grid(1) .or. &
      !& sh%nkpoints(2) /= simpleip_in%interp_grid(2) .or. sh%nkpoints(3) /= simpleip_in%interp_grid(3) ) ) then
      !  call errore('SIMPLE_IP', 'W/o trilinear interpolation, k-grids from simple and simple_ip have to be equal',1)
      !elseif ( simpleip_in%nonlocal_interpolation .and. ( 2*sh%nkpoints(1) /= simpleip_in%interp_grid(1) .or. &
      !1& 2*sh%nkpoints(2) /= simpleip_in%interp_grid(2) .or. 2*sh%nkpoints(3) /= simpleip_in%interp_grid(3) ) ) then
      !  call errore('SIMPLE_IP', 'W/ trilinear interpolation, simple_ip k-grid has to be the double of the simple k-grid',1)
      !endif

      !if( (sh%num_val .ne. sh%num_nbndv(1)) ) then
      !  call errore('OBF_DIAG', 'All the occupied bands must be included for Shirley interpolation (num_val=num_nbndv)',1)
      !endif
      !
    end subroutine read_shirley_spin

    !SJ
      SUBROUTINE read_ham_kgrid(obfdiag_in,kgrid, sh)
      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp_world,  ONLY : mpime, nproc
      USE io_global, ONLY : stdout
      !
      implicit none
      !
      TYPE(input_options_obf_diag) :: obfdiag_in
      TYPE(kpoints) :: kgrid
      TYPE(shirley) :: sh
      INTEGER ::  i, j, k, ii,jj,kk,ll
      INTEGER :: l_blk, ierr
      INTEGER :: iunkpt=999
      REAL, ALLOCATABLE    :: wk(:)
      !
      ! WARNING: this must be the same as in create_energies
      !kgrid%nkgrid(1:3) = simpleip_in%interp_grid(1:3) !SJ No need for grid. We want to just diagonalize
      !kgrid%nk = (kgrid%nkgrid(1))*(kgrid%nkgrid(2))*(kgrid%nkgrid(3)) !This should be read in from input file.
      

      !Manual grid will read list from file.
      IF (obfdiag_in%kpoint_style == 'crystal') THEN

      OPEN( UNIT = iunkpt, FILE = trim(obfdiag_in%kpoint_filename), FORM='FORMATTED', &
      STATUS = 'OLD', IOSTAT = ierr )
      read(iunkpt,*) kgrid%nk
      WRITE(stdout,*) 'Total number of kpoints:', kgrid%nk
      ALLOCATE(kgrid%qk(3,kgrid%nk),wk(kgrid%nk))
      DO i=1,kgrid%nk
        read(iunkpt,*) kgrid%qk(1,i),kgrid%qk(2,i),kgrid%qk(3,i),wk(i)
      ENDDO
      WRITE(stdout,*) 'Done reading in kpoint from file: ' , trim(obfdiag_in%kpoint_filename)
      DO i=1,kgrid%nk
        WRITE(stdout,*) 'xk', kgrid%qk(:,i),wk(i)
      ENDDO
      DEALLOCATE(wk)
      kgrid%alat = sh%alat
      kgrid%bg(1:3,1:3) = sh%bg(1:3,1:3)

      l_blk=kgrid%nk/nproc !FORTRAN always rounds down
      if(l_blk*nproc<kgrid%nk) l_blk=l_blk+1 

      if(l_blk*mpime+1 <= kgrid%nk) then
         kgrid%ik_first=l_blk*mpime+1
         kgrid%ik_last=kgrid%ik_first+l_blk-1
         if(kgrid%ik_last>kgrid%nk) kgrid%ik_last=kgrid%nk
         kgrid%nk_loc=kgrid%ik_last-kgrid%ik_first+1
      else
         kgrid%nk_loc=0
         kgrid%ik_first=0
         kgrid%ik_last=-1
      endif

      !If mesh is built from NK1 NK2 NK3 k1 k2 k3
      ELSEIF ( obfdiag_in%kpoint_style == 'automatic') THEN
      !Firs read the kpoint file
      OPEN( UNIT = iunkpt, FILE = trim(obfdiag_in%kpoint_filename), FORM='FORMATTED', &
      STATUS = 'OLD', IOSTAT = ierr )
      read(iunkpt,*) kgrid%nkgrid(1),kgrid%nkgrid(2),kgrid%nkgrid(3),kgrid%ks(1),kgrid%ks(2),kgrid%ks(3)
      !kgrid%nk=kgrid%ncoarse(1)*kgrid%ncoarse(2)*kgrid%ncoarse(3)
      !WRITE(stdout,*) 'Total number of kpoints:', kgrid%nk
      !ALLOCATE(kgrid%qk(3,kgrid%nk),wk(kgrid%nk))
      !wk=1/kgrid%nk
              !If wrap then we will have kpoint at 1 too. Therefore NK+1 points
              IF ( obfdiag_in%kpoint_wrap ) THEN

                  kgrid%nk=(kgrid%nkgrid(1)+1)*(kgrid%nkgrid(2)+1)*(kgrid%nkgrid(3)+1)
                  WRITE(stdout,*) 'Total number of kpoints:', kgrid%nk
                  WRITE(1000+mpime,*) 'kgrid%nk in read_ham',kgrid%nk
                  ALLOCATE(kgrid%qk(3,kgrid%nk),wk(kgrid%nk))
                  wk=1.d0/kgrid%nk
                  DO ii=0,kgrid%nkgrid(1)
                  DO jj=0,kgrid%nkgrid(2)
                  DO kk=0,kgrid%nkgrid(3)
                  ll = kk + jj*(kgrid%nkgrid(3)+1) + ii*(kgrid%nkgrid(3)+1)*(kgrid%nkgrid(2)+1)+1
                  kgrid%qk(1,ll)=dble(ii)/kgrid%nkgrid(1) !+ dble(kgrid%ks(1))/kgrid%ncoarse(1)
                  kgrid%qk(2,ll)=dble(jj)/kgrid%nkgrid(2) !+ dble(kgrid%ks(2))/kgrid%ncoarse(2)
                  kgrid%qk(3,ll)=dble(kk)/kgrid%nkgrid(3) !+ dble(kgrid%ks(3))/kgrid%ncoarse(3)
                  ENDDO
                  ENDDO
                  ENDDO
              ELSE !Unwrap have Nk point from [0,1). Also allow offset.
                  kgrid%nk=(kgrid%nkgrid(1))*(kgrid%nkgrid(2))*(kgrid%nkgrid(3))
                  WRITE(stdout,*) 'Total number of kpoints:', kgrid%nk
                  ALLOCATE(kgrid%qk(3,kgrid%nk),wk(kgrid%nk))
                  wk=1.d0/kgrid%nk
                  DO ii=1,kgrid%nkgrid(1)
                  DO jj=1,kgrid%nkgrid(2)
                  DO kk=1,kgrid%nkgrid(3)
                  ll = kk-1 + (jj-1)*(kgrid%nkgrid(3)) + (ii-1)*(kgrid%nkgrid(3))*(kgrid%nkgrid(2))+1
                  kgrid%qk(1,ll)=dble(ii-1)/kgrid%nkgrid(1) + dble(kgrid%ks(1))/kgrid%nkgrid(1)
                  kgrid%qk(2,ll)=dble(jj-1)/kgrid%nkgrid(2) + dble(kgrid%ks(2))/kgrid%nkgrid(2)
                  kgrid%qk(3,ll)=dble(kk-1)/kgrid%nkgrid(3) + dble(kgrid%ks(3))/kgrid%nkgrid(3)
                  ENDDO
                  ENDDO
                  ENDDO
              ENDIF

      kgrid%alat = sh%alat
      kgrid%bg(1:3,1:3) = sh%bg(1:3,1:3)

      l_blk=kgrid%nk/nproc !FORTRAN always rounds down
      if(l_blk*nproc<kgrid%nk) l_blk=l_blk+1

      if(l_blk*mpime+1 <= kgrid%nk) then
         kgrid%ik_first=l_blk*mpime+1
         kgrid%ik_last=kgrid%ik_first+l_blk-1
         if(kgrid%ik_last>kgrid%nk) kgrid%ik_last=kgrid%nk
         kgrid%nk_loc=kgrid%ik_last-kgrid%ik_first+1
      else
         kgrid%nk_loc=0
         kgrid%ik_first=0
         kgrid%ik_last=-1
      endif

      ELSE 
              WRITE(stdout,*) 'Wrong kpoint_style. Exiting'
              STOP
      ENDIF


      !Convert crystal kpoint to cart
      CALL cryst_to_cart( kgrid%nk, kgrid%qk, kgrid%bg, 1 )

      CLOSE(iunkpt)

      !
    END SUBROUTINE read_ham_kgrid
    !SJ

    SUBROUTINE interpolated_kgrid_create(obfdiag_in,kgrid,kgrid_fine,sh)
      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp_world,  ONLY : mpime, nproc
    implicit none

    TYPE(input_options_obf_diag) :: obfdiag_in
    TYPE(kpoints) :: kgrid,kgrid_fine
    TYPE(shirley) :: sh
    INTEGER  :: i,j,k,ii,jj,kk,ll,ierr
    INTEGER  :: iunkpt=999
    INTEGER  :: i_blk,ntot,l_blk
    REAL(DP) :: ks(3) 
    REAL(DP), ALLOCATABLE :: wk(:) 
    INTEGER  :: nk(3)

    WRITE(stdout,*) 'Using interpolation'

    !Real in list if the style is crystal
    IF (obfdiag_in%interp_style == 'crystal') THEN
      OPEN( UNIT = iunkpt, FILE = trim(obfdiag_in%interp_kpoint_filename), FORM='FORMATTED', &
      STATUS = 'OLD', IOSTAT = ierr )
      read(iunkpt,*) kgrid_fine%nk
      WRITE(stdout,*) 'Total number of kpoints:', kgrid_fine%nk
      ALLOCATE(kgrid_fine%qk(3,kgrid_fine%nk),wk(kgrid_fine%nk))
      DO i=1,kgrid_fine%nk
        read(iunkpt,*) kgrid_fine%qk(1,i),kgrid_fine%qk(2,i),kgrid_fine%qk(3,i),wk(i)
      ENDDO
      WRITE(stdout,*) 'Done reading in kpoint from file: ' , trim(obfdiag_in%interp_kpoint_filename)
      DO i=1,kgrid_fine%nk
        WRITE(stdout,*) 'xk', kgrid_fine%qk(:,i),wk(i)
      ENDDO
      DEALLOCATE(wk)
      kgrid_fine%alat = sh%alat
      kgrid_fine%bg(1:3,1:3) = sh%bg(1:3,1:3)

!If mesh is built from NK1 NK2 NK3 k1 k2 k3
      ELSEIF ( obfdiag_in%interp_style == 'automatic') THEN
      !Firs read the kpoint file
      OPEN( UNIT = iunkpt, FILE = trim(obfdiag_in%interp_kpoint_filename), FORM='FORMATTED', &
      STATUS = 'OLD', IOSTAT = ierr )
      read(iunkpt,*) kgrid_fine%nkgrid(1),kgrid_fine%nkgrid(2),kgrid_fine%nkgrid(3),kgrid_fine%ks(1),&
              &kgrid_fine%ks(2),kgrid_fine%ks(3)
      !No need to wrap to 1 anymore. Now we are only interested in builing u_nk in the first BZ [0,1).

      kgrid_fine%nk=kgrid_fine%nkgrid(1)*kgrid_fine%nkgrid(2)*kgrid_fine%nkgrid(3)
      kgrid_fine%bg(1:3,1:3) = sh%bg(1:3,1:3)
      kgrid_fine%alat = sh%alat
      WRITE(stdout,*) 'Total number of kpoints:', kgrid_fine%nk
      ALLOCATE(kgrid_fine%qk(3,kgrid_fine%nk),wk(kgrid_fine%nk))
                  wk=1.d0/kgrid_fine%nk
                  DO ii=1,kgrid_fine%nkgrid(1)
                  DO jj=1,kgrid_fine%nkgrid(2)
                  DO kk=1,kgrid_fine%nkgrid(3)
                  ll = kk-1 + (jj-1)*(kgrid_fine%nkgrid(3)) + (ii-1)*(kgrid_fine%nkgrid(3))*(kgrid_fine%nkgrid(2))+1
                  kgrid_fine%qk(1,ll)=dble(ii-1)/kgrid_fine%nkgrid(1) + dble(kgrid_fine%ks(1))/kgrid_fine%nkgrid(1)
                  kgrid_fine%qk(2,ll)=dble(jj-1)/kgrid_fine%nkgrid(2) + dble(kgrid_fine%ks(2))/kgrid_fine%nkgrid(2)
                  kgrid_fine%qk(3,ll)=dble(kk-1)/kgrid_fine%nkgrid(3) + dble(kgrid_fine%ks(3))/kgrid_fine%nkgrid(3)
                  WRITE(stdout,*) 'xk', kgrid_fine%qk(:,ll),wk(ll)
                  ENDDO
                  ENDDO
                  ENDDO

      ELSE
              WRITE(stdout,*) 'Wrong kpoint_style. Exiting'
              STOP
      ENDIF

      !Equally divide kpoint among each proc
      l_blk=kgrid_fine%nk/nproc !FORTRAN always rounds down
      if(l_blk*nproc<kgrid_fine%nk) l_blk=l_blk+1

      if(l_blk*mpime+1 <= kgrid_fine%nk) then
         kgrid_fine%ik_first=l_blk*mpime+1
         kgrid_fine%ik_last=kgrid_fine%ik_first+l_blk-1
         if(kgrid_fine%ik_last>kgrid_fine%nk) kgrid_fine%ik_last=kgrid_fine%nk
         kgrid_fine%nk_loc=kgrid_fine%ik_last-kgrid_fine%ik_first+1
      else
         kgrid_fine%nk_loc=0
         kgrid_fine%ik_first=0
         kgrid_fine%ik_last=-1
      endif 


    CALL cryst_to_cart( kgrid_fine%nk, kgrid_fine%qk, kgrid%bg, 1 ) !Conver from cryst to cart
    close(iunkpt)   

    END SUBROUTINE interpolated_kgrid_create

        SUBROUTINE create_energies_flat(sh,kgrid,ene)
      USE mp_world,  ONLY : mpime, nproc
      USE io_global, ONLY : stdout
      implicit none

      TYPE(shirley) :: sh
      TYPE(kpoints) :: kgrid
      TYPE(energies) :: ene
      INTEGER :: l_blk
      !SJ
      WRITE(stdout,*) 'Calling create_energies_flat'
      ene%nk = kgrid%nk
      WRITE(stdout,*) 'NK total from kgrid%nk', ene%nk
      ene%num_bands = sh%num_bands

      l_blk=ene%nk/nproc
      if(l_blk*nproc<ene%nk) l_blk=l_blk+1

      if(l_blk*mpime+1 <= ene%nk) then !mpime = index of proc starting from 0 -> nproc-1
         ene%ik_first=l_blk*mpime+1
         ene%ik_last=ene%ik_first+l_blk-1
         if(ene%ik_last>ene%nk) ene%ik_last=ene%nk
         ene%nk_loc=ene%ik_last-ene%ik_first+1
      else
         ene%nk_loc=0
         ene%ik_first=0
         ene%ik_last=-1
      endif
       
      allocate(ene%energy(ene%num_bands,ene%nk_loc), ene%energy_der(3,ene%num_bands,ene%nk_loc))

    END SUBROUTINE create_energies_flat

!SJ
    SUBROUTINE read_shirley_k_flat(obfdiag_in,sh,ene)
      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp,                   ONLY : mp_bcast, mp_barrier
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      !
      IMPLICIT NONE
      !
      TYPE(input_options_obf_diag) :: obfdiag_in
      TYPE(shirley) :: sh
      TYPE(energies) :: ene
      !
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun, ll, ii , jj , kk
      COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: sum_beckc_tmp
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: sum_beck_nc_tmp 
      !
      write(stdout,*)'obf_diag: opening file "hamiltonian_k" flat'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian_k', status='old',form='unformatted')
         write(stdout,*)'File opened'
      endif
      !

     
      if (sh%noncolin) then
        allocate(sum_beck_nc_tmp(sh%nkb,sh%npol,sh%ntot_e,sh%nkp_tot)) !From sh%nkpoints(3) flatening to nkp_tot
        allocate(sh%beck_nc(sh%nkb,sh%npol,sh%ntot_e,ene%nk_loc))
      else
        allocate(sum_beckc_tmp(sh%nkb,sh%ntot_e,sh%nkp_tot))
        allocate(sh%beckc(sh%nkb,sh%ntot_e,ene%nk_loc))
      endif
      !
      !if (sh%nonlocal_commutator) then
      !  allocate(sum_commut_tmp(sh%ntot_e,sh%ntot_e,3,sh%nkp_tot))
      !  allocate(sh%commut(sh%ntot_e,sh%ntot_e,3,ene%nk_loc))
      !endif
      !
      write(stdout,*) 'Total number of k-point (from Hamilotnian_k): ' , sh%nkp_tot
      if ( ene%nk /= sh%nkp_tot) then
              write(stdout,*) 'Numbers of k-point from the input and Hamiltonian_k are different'
              write(stdout,*) 'sh%nkp_tot',sh%nkp_tot
              write(stdout,*) 'ene%nk', ene%nk
              write(stdout,*) 'Stopping'
              STOP
      endif

      do ii=1,sh%nkp_tot
        !
        write(stdout,*) 'k-point index: ' , ii
        if (sh%noncolin) then
          if (ionode) read(iun) sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e,ii) !Read from blcok nkpoints(3) to 
        else                                                                       ! 1 at a time
          if (ionode) read(iun) sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e,ii)
        endif
        !
        !if (sh%nonlocal_commutator) then
        !  if (ionode) read(iun) sum_commut_tmp(1:sh%ntot_e,1:sh%ntot_e,1:3,ii)  ! read commutator matrix
        !endif
        !
        if (sh%noncolin) then
          call mp_bcast(sum_beck_nc_tmp,ionode_id,world_comm)
        else
          call mp_bcast(sum_beckc_tmp,ionode_id,world_comm)
        endif
        !
        !if (sh%nonlocal_commutator) then
        !  call mp_bcast(sum_commut_tmp,ionode_id,world_comm)
        !endif
        !
        !do kk=1,sh%nkpoints(3)  No longer this a loop since it's one by one kpoint now.
          !
          !ll = kk + (jj - 1)*sh%nkpoints(3) + (ii - 1)*sh%nkpoints(2)*sh%nkpoints(3)   ! global kindex
          !
          if (ene%ik_first <= ii .and. ii <= ene%ik_last) then !global index is now ll -> ii
            !
            if (sh%noncolin) then
              sh%beck_nc(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e,ii - ene%ik_first + 1) = &
              & sum_beck_nc_tmp(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e, ii) !SJ kk -> ii
            else
              sh%beckc(1:sh%nkb, 1:sh%ntot_e, ii - ene%ik_first + 1) = &
              & sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e, ii) !SJ kk -> ii
            endif
            !
            !if (sh%nonlocal_commutator) then
            !  sh%commut(1:sh%ntot_e, 1:sh%ntot_e, 1:3, ii - ene%ik_first + 1) = &
            !  & sum_commut_tmp(1:sh%ntot_e, 1:sh%ntot_e, 1:3, ii)
            !endif
          endif
          !
        !enddo
        !
        !enddo
      enddo
      !
      if(ionode) then
        close(iun)
        write(stdout,*)'File closed'
      endif
      !
      if (sh%noncolin) then
        deallocate(sum_beck_nc_tmp)
      else
        deallocate(sum_beckc_tmp)
      endif
      !
      !if (sh%nonlocal_commutator) then
      !  deallocate(sum_commut_tmp)
      !endif
      !
    END SUBROUTINE read_shirley_k_flat


    SUBROUTINE read_shirley_k_byfile(obfdiag_in,sh,ene)
      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp,                   ONLY : mp_bcast, mp_barrier
      USE mp_world,             ONLY : world_comm,mpime
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      !
      IMPLICIT NONE
      !
      TYPE(input_options_obf_diag) :: obfdiag_in
      TYPE(shirley) :: sh
      TYPE(energies) :: ene
      !
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun, ll, ii , jj , kk,kindex
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: sum_beckc_tmp
      COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: sum_beck_nc_tmp
      !
      write(stdout,*) 'obf_diag: opening file "hamiltonian_k" by file'
      write(stdout,*) 'obf_diag: expects hamiltonian_k to be written into multiple individual files'
      !if(ionode) then
      !   iun = find_free_unit()
      !   open( unit=iun, file=trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian_k', status='old',form='unformatted')
      !   write(stdout,*)'File opened'
      !endif
      !


      if (sh%noncolin) then
        allocate(sum_beck_nc_tmp(sh%nkb,sh%npol,sh%ntot_e)) !From sh%nkpoints(3) flatening to nkp_tot
        allocate(sh%beck_nc(sh%nkb,sh%npol,sh%ntot_e,ene%nk_loc))
      else
        allocate(sum_beckc_tmp(sh%nkb,sh%ntot_e))
        allocate(sh%beckc(sh%nkb,sh%ntot_e,ene%nk_loc))
      endif
      !
      !if (sh%nonlocal_commutator) then
      !  allocate(sum_commut_tmp(sh%ntot_e,sh%ntot_e,3,sh%nkp_tot))
      !  allocate(sh%commut(sh%ntot_e,sh%ntot_e,3,ene%nk_loc))
      !endif
      !
      write(stdout,*) 'Total number of k-point (from Hamilotnian_k): ' , sh%nkp_tot

      if ( ene%nk /= sh%nkp_tot) then
              write(stdout,*) 'Numbers of k-point from the input and Hamiltonian_k are different'
              write(stdout,*) 'sh%nkp_tot',sh%nkp_tot
              write(stdout,*) 'ene%nk', ene%nk
              write(stdout,*) 'Stopping'
              STOP
      endif

      !write(3000+mpime,*) 'ene%nk_loc',ene%nk_loc

      do ii=1,ene%nk_loc
        !
        kindex = ii + ene%ik_first - 1
        write(stdout,*) 'k-point index: ' , kindex
        !write(3000+mpime,*) 'k-point index: ', kindex
        !write(mpime+100,*) 'k-point index: ' , kindex

        iun = find_free_unit()

        !write(3000+mpime,*) 'filename ',trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian_k'//TRIM(int_to_char(kindex))
        !write(3000+mpime,*) 'iun',iun
        !flush(3000+mpime)
        open( unit=iun, file=trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian_k'//TRIM(int_to_char(kindex)),&
        status='old',form='unformatted')
        write(stdout,*)'File opened: Hamiltonian_k at ik', kindex

         if (sh%noncolin) then
                read(iun) sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e) !Read the block file.
         else                                                                          !No more k-point index.
                read(iun) sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e)
         endif
        !
            !
         if (sh%noncolin) then
              sh%beck_nc(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e,ii ) = &
              & sum_beck_nc_tmp(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e) 
         else
              sh%beckc(1:sh%nkb, 1:sh%ntot_e, ii ) = &
              & sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e) 
         endif
            !
          !
         close(iun)
         write(stdout,*)'File closed: Hamiltonian_k at ik', kindex
        !enddo
        !
        !enddo
      enddo
      !
      !if(ionode) then
       ! close(iun)
       ! write(stdout,*)'File closed'
      !endif
      !
      if (sh%noncolin) then
        deallocate(sum_beck_nc_tmp)
      else
        deallocate(sum_beckc_tmp)
      endif
      !
      !if (sh%nonlocal_commutator) then
      !  deallocate(sum_commut_tmp)
      !endif
      !
    END SUBROUTINE read_shirley_k_byfile


    SUBROUTINE read_shirley_k_interp(obfdiag_in,sh,ene,k,k_interp)

      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp,                   ONLY : mp_bcast, mp_barrier
      USE mp_world,             ONLY : world_comm,mpime
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      !
      IMPLICIT NONE
      !
      TYPE(input_options_obf_diag) :: obfdiag_in
      TYPE(shirley) :: sh
      TYPE(energies) :: ene
      TYPE(kpoints) :: k,k_interp
      !
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: ierr,iunkpt
      INTEGER :: i,iun, ll, ii , jj , kk,nleft(3),nright(3),ik,ncoarse(3),counter
      INTEGER,DIMENSION(:,:), ALLOCATABLE :: coord_cube_global
      INTEGER, DIMENSION(:), ALLOCATABLE  :: diff_kpoints
      REAL(DP) :: d,q(3),dq(3)
      COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: sum_beckc_tmp
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: sum_beck_nc_tmp !, sum_commut_tmp

      !We will try to build the coordinates of the kinterp. Each will point to global index of the Vnl
      !from the k-dependent Hamiltonian. 
      !This will be collected along with the respective postiion along interpolated path. 
      !All these will be used to build Vnl from trilinear interpolation
      !For distance and corresponding position

      !First check for Hamiltonian kpoint_style. If we are not using uniform mesh to build H.
      !This would not work.
      IF (obfdiag_in%kpoint_style == 'crystal') THEN
              WRITE(stdout,*) 'Need uniform mesh for Vnl interpolation.'
              WRITE(stdout,*) 'Built the Hamiltonian with kpoint_style = automatic'
              STOP
      ENDIF

      !Allocate for interpolation stuff.
      allocate(k_interp%pos_cube(3,k_interp%nk_loc),k_interp%coord_cube(8,k_interp%nk_loc))
      allocate(coord_cube_global(8,k_interp%nk_loc),diff_kpoints(k%nk))

      !Convert back to cryst.
      WRITE(stdout,*) 'sh%at',sh%at
      WRITE(stdout,*) 'sh%bg',sh%bg
      WRITE(stdout,*) 'k_interp%bg',k_interp%bg
      WRITE(stdout,*) 'k%bg', k%bg
      CALL cryst_to_cart( k_interp%nk, k_interp%qk, sh%at, -1 )

      !Determine number of kpoint in each direction.
      !Will be use for global indexing. Rounding part will be a bit manual.
      !If we do not wrap to 1.00 then we need to do interpolation using point (i-1)/NK.
      !This means we will jump back by 1 so nleft/nright will be the same.
      !IF we have 1.00. Then it will be go up to maximum of nkgrid(i)+1.
      IF (obfdiag_in%kpoint_wrap) THEN
              ncoarse(1:3) = k%nkgrid(1:3)+1
      ELSE
              ncoarse(1:3) = k%nkgrid(1:3)
      ENDIF

      !Coarse grid size
      do i=1,3
       dq(i) = dble(1)/k%nkgrid(i)
      enddo

      WRITE(stdout,*) 'ncoarse',ncoarse
      WRITE(stdout,*) 'dq',dq
      WRITE(stdout,*) 'k_interp%nk_loc',k_interp%nk_loc
      diff_kpoints=0

      do ik=1,k_interp%nk_loc
          !First get the q
          do i=1,3
          q(1:3)=k_interp%qk(1:3,ik+k_interp%ik_first -1) 
          enddo

          !Relative position to known k-point in obf_ham Hamiltonian mesh.
          nleft=0
          nright=0
          do i=1,3
             !Find relative position of the index of the left point.
             nleft(i) = int( q(i)*k%nkgrid(i) ) !Fortran will round down.
                                                ! Ex if we have nkgrid = 4. 
                                                !This gives 0,0.25,0.50,0.75 (and possibly 1.00)
                                                !Suppose q = 0.05 -> nleft = int(0.2) = 0
             nright(i) = nleft(i) + 1
             !Rounding is dependent on kpoint_wrap. Whether we have NK or NK+1 point.
             !For NK points, we will use previous k-point as there is no 1.00.
             !For NK+1 points, we will do correction such that the index is not out of bound just in case.
             IF (obfdiag_in%kpoint_wrap) THEN
                     IF ( nright(i) > k%nkgrid(i) + 1 ) nright(i)=k%nkgrid(i)+1
             ELSE
                     IF ( nright(i) >= k%nkgrid(i) ) nright(i) = k%nkgrid(i)-1
             ENDIF
          enddo

          WRITE(stdout,*) 'q',q
          WRITE(stdout,*) 'nleft',nleft
          WRITE(stdout,*) 'nright', nright

          !Now find the global index of each point
          !Loop z -> block z * index_y -> slab zxy * index_x
          coord_cube_global(1,ik) = nleft(3) + 1  + nleft(2) * ncoarse(3) + nleft(1) * ncoarse(3) * ncoarse(2)
          coord_cube_global(2,ik) = nleft(3) + 1  + nleft(2) * ncoarse(3) + nright(1) * ncoarse(3) * ncoarse(2)
          coord_cube_global(3,ik) = nleft(3) + 1  + nright(2) * ncoarse(3) + nleft(1) * ncoarse(3) * ncoarse(2)
          coord_cube_global(4,ik) = nleft(3) + 1  + nright(2) * ncoarse(3) + nright(1) * ncoarse(3) * ncoarse(2)
          coord_cube_global(5,ik) = nright(3) + 1  + nleft(2) * ncoarse(3) + nleft(1) * ncoarse(3) * ncoarse(2)
          coord_cube_global(6,ik) = nright(3) + 1  + nleft(2) * ncoarse(3) + nright(1) * ncoarse(3) * ncoarse(2)
          coord_cube_global(7,ik) = nright(3) + 1  + nright(2) * ncoarse(3) + nleft(1) * ncoarse(3) * ncoarse(2)
          coord_cube_global(8,ik) = nright(3) + 1  + nright(2) * ncoarse(3) + nright(1) * ncoarse(3) * ncoarse(2)

          WRITE(stdout,*) 'coord_cube_global(1)',coord_cube_global(1,ik)
          WRITE(stdout,*) 'coord_cube_global(2)',coord_cube_global(2,ik)
          WRITE(stdout,*) 'coord_cube_global(3)',coord_cube_global(3,ik)
          WRITE(stdout,*) 'coord_cube_global(4)',coord_cube_global(4,ik)
          WRITE(stdout,*) 'coord_cube_global(5)',coord_cube_global(5,ik)
          WRITE(stdout,*) 'coord_cube_global(6)',coord_cube_global(6,ik)
          WRITE(stdout,*) 'coord_cube_global(7)',coord_cube_global(7,ik)
          WRITE(stdout,*) 'coord_cube_global(8)',coord_cube_global(8,ik)


          !Find relative distance of the point between 2 points in coarse grid.
          do i=1,3
           k_interp%pos_cube(i,ik) = ( q(i) - dble(nleft(i))*dq(i) ) / dq(i)
          enddo

          do i=1,8
           diff_kpoints(coord_cube_global(i,ik)) = 1
          enddo
      enddo   

      !Next determine the number of kpoint from coarse grid needed to build Vnl in each local pro
      k_interp%nk_smooth_loc = sum(diff_kpoints)

      !WRITE(mpime+1000,*) 'nk_smooth_loc',k_interp%nk_smooth_loc
      !WRITE(mpime+1000,*) 'sh%nkb,sh%ntot_e,k%nk',sh%nkb,sh%ntot_e,k%nk
      !flush(mpime+1000)
      WRITE(stdout,*) 'k%nk',k%nk

      !global coarse k -> local for interp grid.
      counter = 0
      do i = 1,k%nk !Loop through coarse grid kpoint
        if (diff_kpoints(i) /= 0) then
          counter = counter + 1
          diff_kpoints(i) = counter !Now it be consecutive list of 1 -> nk_smooth_loc with 0 in between
                                     !if the points are not used for interpolation
        endif
      enddo

      WRITE(stdout,*) 'diff_kpoints',diff_kpoints
      !Now map global coordinate of the kpoint into locally store one.
      do ik=1,k_interp%nk_loc 
        do i=1,8
         k_interp%coord_cube(i,ik) = diff_kpoints(coord_cube_global(i,ik))
        enddo
      enddo

      !Next open the .hamiltonian_k file to read. We then read the Vnl part one by one.
      !If the index is in the list, then it's read into shirley data type.
      write(stdout,*)'obf_diag: opening file "hamiltonian_k" for coarse grid kpoint'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian_k', status='old',form='unformatted')
         write(stdout,*)'File opened'
      endif

      ! Hamiltonain file is list of kpoint. Read in one by one. If it matches the diff kpoints.
      ! Then it is read into type(shirley). Will read the whole k-dependent Hamiltonian then put into
      ! sh data type.
      if (sh%noncolin) then
        allocate(sh%beck_nc(sh%nkb,sh%npol,sh%ntot_e,k_interp%nk_smooth_loc)) ! (nkb,npol,ntot_e,nk)
        allocate(sum_beck_nc_tmp(sh%nkb,sh%npol,sh%ntot_e,k%nk))
      else
        allocate(sh%beckc(sh%nkb,sh%ntot_e,k_interp%nk_smooth_loc))
        allocate(sum_beckc_tmp(sh%nkb,sh%ntot_e,k%nk))
      endif
     ! if (sh%nonlocal_commutator) then
     !   allocate(sh%commut(sh%ntot_e,sh%ntot_e,3,k_interp%nk_smooth_loc))
     !   allocate(sum_commut_tmp(sh%ntot_e,sh%ntot_e,3,k%nk))
     ! endif

      !Now read required H into each proc
      do ik=1,k%nk
        WRITE(stdout,*) 'k-point index from coarse grid: ',ik
        ! WRITE(stdout,*) 'sh%noncolin',sh%noncolin
      !  WRITE(stdout,*) 'sh%nonlocal_commutator',sh%nonlocal_commutator
        !WRITE(stdout,*) 'Size of sh%beck',size(sh%beckc,1),size(sh%beckc,2),size(sh%beckc,3)
        !WRITE(stdout,*) 'Size of sum_beck_tmp',size(sum_beck_nc_tmp,1),size(sum_beck_nc_tmp,2),size(sum_beck_nc_tmp,3)
       ! WRITE(stdout,*) 'Size of sh%commut',size(sh%commut,1),size(sh%commut,2),size(sh%commut,3),size(sh%commut,4)
       ! WRITE(stdout,*) 'Size of sum_commut_tmp',size(sum_commut_tmp,1),size(sum_commut_tmp,2),&
       !         &size(sum_commut_tmp,3),size(sum_commut_tmp,4)
        if (sh%noncolin) then
          if (ionode) read(iun) sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e,ik) !Read from blcok nkpoints(3) to
        else                                                                       ! 1 at a time
          if (ionode) read(iun) sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e,ik)
        endif
        WRITE(stdout,*) 'read in sum_beckc_tmp for ik',ik
        !
        !if (sh%nonlocal_commutator) then
        !  if (ionode) read(iun) sum_commut_tmp(1:sh%ntot_e,1:sh%ntot_e,1:3,ik)  ! read commutator matrix
        !endif
        !WRITE(stdout,*) 'read in sum_commut_tmp'
        !WrITE(stdout,*) 'ionode_id',ionode_id
        !WRITE(stdout,*) 'world_comm',world_comm
        !
        if (sh%noncolin) then
          call mp_bcast(sum_beck_nc_tmp,ionode_id,world_comm)
        else
          call mp_bcast(sum_beckc_tmp,ionode_id,world_comm)
        endif
        !
        !WRITE(stdout,*) 'done with mp_bcast sum_bec_tmp'
        !if (sh%nonlocal_commutator) then
        !  call mp_bcast(sum_commut_tmp,ionode_id,world_comm)
        !endif

        !WRITE(stdout,*) 'done with mp_bcast sum_commut_tmp'

           if (diff_kpoints(ik) /= 0 ) then
                   if (sh%noncolin) then
                           sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e,diff_kpoints(ik)) =&
                                   &sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e,ik)
                   else
                           sh%beckc(1:sh%nkb,1:sh%ntot_e,diff_kpoints(ik)) = sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e,ik)
                   endif
                 !  if (sh%nonlocal_commutator) then
                  !         sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3,diff_kpoints(ik)) = sum_commut_tmp(1:sh%ntot_e,1:sh%ntot_e,1:3,ik)
                  ! endif
           endif


      enddo

      if(ionode) then
        close(iun)
        write(stdout,*)'File closed'
      endif
      !
      if (sh%noncolin) then
        deallocate(sum_beck_nc_tmp)
      else
        deallocate(sum_beckc_tmp)
      endif
      !
      !if (sh%nonlocal_commutator) then
      !  deallocate(sum_commut_tmp)
      !endif

      !Convert back to cart
      CALL cryst_to_cart( k_interp%nk, k_interp%qk, k_interp%bg, -1 )

    END SUBROUTINE read_shirley_k_interp

    SUBROUTINE read_shirley_k_flat_spin(obfdiag_in,sh,ene)
      USE input_obf_diag, ONLY : input_options_obf_diag
      USE mp,                   ONLY : mp_bcast, mp_barrier
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      !
      IMPLICIT NONE
      !
      TYPE(input_options_obf_diag) :: obfdiag_in
      TYPE(shirley) :: sh
      TYPE(energies) :: ene
      !
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun, ll, ii , jj , kk
      COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: sum_beckc_tmp
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: sum_beck_nc_tmp !, sum_commut_tmp
      !
      write(stdout,*)'obf_diag: opening file "hamiltonian_k" flat'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(obfdiag_in%prefix)//'.hamiltonian_k_down', status='old',form='unformatted')
         write(stdout,*)'File opened'
      endif
      !
      if (sh%noncolin) then
        allocate(sum_beck_nc_tmp(sh%nkb,sh%npol,sh%ntot_e,sh%nkp_tot)) !From sh%nkpoints(3) flatening to nkp_tot
        allocate(sh%beck_nc(sh%nkb,sh%npol,sh%ntot_e,ene%nk_loc))
      else
        allocate(sum_beckc_tmp(sh%nkb,sh%ntot_e,sh%nkp_tot))
        allocate(sh%beckc(sh%nkb,sh%ntot_e,ene%nk_loc))
      endif
      !
      !if (sh%nonlocal_commutator) then
      !  allocate(sum_commut_tmp(sh%ntot_e,sh%ntot_e,3,sh%nkp_tot))
      !  allocate(sh%commut(sh%ntot_e,sh%ntot_e,3,ene%nk_loc))
      !endif
      !
      write(stdout,*) 'Total number of k-point (for reading): ' , sh%nkp_tot
      do ii=1,sh%nkp_tot
        !
        write(stdout,*) 'k-point index: ' , ii
        if (sh%noncolin) then
          if (ionode) read(iun) sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e,ii) !Read from blcok nkpoints(3) to
        else                                                                       ! 1 at a time
          if (ionode) read(iun) sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e,ii)
        endif
        !
        !if (sh%nonlocal_commutator) then
        !  if (ionode) read(iun) sum_commut_tmp(1:sh%ntot_e,1:sh%ntot_e,1:3,ii)  ! read commutator matrix
        !endif
        !
        if (sh%noncolin) then
          call mp_bcast(sum_beck_nc_tmp,ionode_id,world_comm)
        else
          call mp_bcast(sum_beckc_tmp,ionode_id,world_comm)
        endif
        !
        !if (sh%nonlocal_commutator) then
        !  call mp_bcast(sum_commut_tmp,ionode_id,world_comm)
        !endif
        !
        !do kk=1,sh%nkpoints(3)  No longer this a loop since it's one by one kpoint now.
          !
          !ll = kk + (jj - 1)*sh%nkpoints(3) + (ii - 1)*sh%nkpoints(2)*sh%nkpoints(3)   ! global kindex
          !
          if (ene%ik_first <= ii .and. ii <= ene%ik_last) then !global index is now ll -> ii
            !
            if (sh%noncolin) then
              sh%beck_nc(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e,ii - ene%ik_first + 1) = &
              & sum_beck_nc_tmp(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e, ii) !SJ kk -> ii
            else
              sh%beckc(1:sh%nkb, 1:sh%ntot_e, ii - ene%ik_first + 1) = &
              & sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e, ii) !SJ kk -> ii
            endif
            !
            !if (sh%nonlocal_commutator) then
            !  sh%commut(1:sh%ntot_e, 1:sh%ntot_e, 1:3, ii - ene%ik_first + 1) = &
            !  & sum_commut_tmp(1:sh%ntot_e, 1:sh%ntot_e, 1:3, ii)
            !endif
          endif
          !
        !enddo
        !
        !enddo
      enddo
      !
      if(ionode) then
        close(iun)
        write(stdout,*)'File closed'
      endif
      !
      if (sh%noncolin) then
        deallocate(sum_beck_nc_tmp)
      else
        deallocate(sum_beckc_tmp)
      endif
      !
      !if (sh%nonlocal_commutator) then
      !  deallocate(sum_commut_tmp)
      !endif
      !
    END SUBROUTINE read_shirley_k_flat_spin
!SJ
    !


END MODULE obf_diag_objects
