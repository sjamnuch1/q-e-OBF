
subroutine dielectric(sh,input,kpnts,energy)
  !------------------------------------------------------------------------
  !
  ! This subroutine diagonalizes the k-dependent Hamiltonian for every k in interp_grid
  ! and calculates the IP dielectric function
  !
  USE io_files,  ONLY : tmp_dir !SJ
  USE obf_diag_objects
  USE input_obf_diag
  USE tetra_ip, ONLY : tetrahedra1, weights_delta1
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE io_files,  ONLY : create_directory
  USE constants, ONLY : rytoev, pi
  USE mp_world,  ONLY : world_comm,mpime
  USE mp,        ONLY : mp_sum
  !USE cell_base, ONLY : tpiba
  !
  IMPLICIT NONE
  !
  TYPE(shirley) :: sh
  TYPE(input_options_obf_diag) :: input
  TYPE(energies) :: energy
  TYPE(kpoints) :: kpnts
  TYPE(eigen) :: eig
  REAL(kind=DP) :: q(3)
  INTEGER :: ik, idir, ii , counting , iun , iunbnk, ibnk
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: tmp , tmp2
  !COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: commut_interp
  REAL(kind=DP), EXTERNAL :: efermig , efermit
  REAL(kind=DP), EXTERNAL ::  w0gauss
  REAL(kind=DP), EXTERNAL ::  wgauss
  INTEGER ::  is , iw , iband1, iband2 , num_energy_dos ,im
  REAL(kind=DP) :: nelec, ef, arg, tpiovera, alpha, diff_occ, delta_e, w
  REAL(kind=DP), DIMENSION(6) :: plasma_freq
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: wk
  INTEGER, DIMENSION(:), ALLOCATABLE :: isk
  REAL(kind=DP), DIMENSION(:,:), ALLOCATABLE :: energy_sum
  REAL(kind=DP), DIMENSION(:,:), ALLOCATABLE :: focc, focc_tmp, eps_im, eps_re, eps_tot_re, eps_tot_im, eels
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: wgrid, wgrid_dos, dos , jdos , eps_avg_re, eps_avg_im 
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: refractive_index , extinct_coeff , reflectivity
  REAL(kind=DP) :: delta_e_thr = 1.0E-6   ! in eV
  INTEGER :: num_tetra, iuneig
  INTEGER, allocatable :: tetra(:,:)
  REAL(kind=DP), allocatable :: weights(:) 
  CHARACTER(len=256), DIMENSION(6) :: headline
  CHARACTER(256) :: str, filename, save_dir,eigenfile
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=6), EXTERNAL :: int_to_char !SJ
  !COMPLEX(kind=DP) :: wfc(:,:), unk(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: wfc(:,:), unk(:,:)
  INTEGER :: igwx,npol,nbnd,nbasis,iglobal_k,nbnd_save
  !
  call start_clock('dielectric')
  !
  call initialize_eigen(eig)
  !
  !allocate(tmp(sh%ntot_e,sh%num_bands))
  !allocate(tmp2(sh%num_bands,sh%num_bands))

  !if (input%nonlocal_commutator .and. sh%nonlocal_commutator) then
  !  allocate(commut_interp(sh%ntot_e,sh%ntot_e,3))
  !endif

  allocate(focc(sh%num_bands,kpnts%nk),focc_tmp(sh%num_bands,energy%nk_loc))
  !allocate(wgrid(input%nw),jdos(input%nw),eps_im(input%nw,3),eps_re(input%nw,3),eels(input%nw,3))
  !
  !!!!!!! Interpolation of the bands
  counting = 0  !DEBUG
  tpiovera = 2.d0*pi/sh%alat
  write(stdout,*) ' '
  write(stdout,*) 'Computing the band structure...'
  write(stdout,*) ' '
  write(stdout,*) 'energy%nk_loc',energy%nk_loc
  !write(1000+mpime,*) mpime,energy%nk_loc
  !flush(1000+mpime)
  !write(2000+mpime,*) mpime,energy%ik_first
  !flush(2000+mpime)

  IF (ionode) THEN
  eigenfile='eigen.txt' 
  iuneig=find_free_unit()
  open(unit=iuneig,file=trim(eigenfile),status='unknown',form='formatted')
  write(iuneig,999) sh%num_bands,kpnts%nk,sh%nspin
  999 format(3I6)
  close(iuneig)
  ENDIF

  do ik=1,energy%nk_loc
     !
     !write(1000+mpime,*) mpime,ik
     !flush(1000+mpime) 
     q(1:3) = kpnts%qk(1:3,ik + energy%ik_first - 1) !ik + energy%ik_first - 1 
     !TODO modify this to take kpnts of interpolate grid
     call diagonalization(q,sh,input,eig,ik,kpnts)
     energy%energy(1:sh%num_bands,ik) = eig%energy(1:sh%num_bands) ! Save bands in energy%energy
     !
     counting = counting + 1  !DEBUG
     iglobal_k = ik + energy%ik_first - 1
     write(stdout,*)  '*****************************************************************'  !DEBUG
     write(stdout,*) 'k-point:', counting
     write(stdout,*) 'k-point coordinates:', q  !DEBUG
     write(stdout,*) 'Interpolated bands:' , energy%energy(1:sh%num_bands,ik)*rytoev !DEBUG
     write(stdout,*) ' '
     !write(iuneig,998) energy%energy(1:sh%num_bands,ik) !Should be moved to after diagonalization 
     !998 format(*(f25.16))                              !Now loop from ik=1 -> total kpoint
     write(stdout,*) 'Diagonalized matrix size ',size(eig%wave_func,1),sh%num_bands,size(eig%wave_func) !SJ
!     save_dir=trim(tmp_dir)//trim(input%prefix)//'.save/'
     save_dir=trim(input%unkdir)//'/'//trim(input%prefix)//'.save/'

     !Create the directory for saving unk
     CALL create_directory( trim(input%unkdir) ) 
     CALL create_directory( save_dir )
     if ( input%lcollect_bnk ) then
       write(stdout,*) 'Saving diagonalized matrix b_nk(i) to disk'
       write(stdout,*) 'Writing to ',trim(tmp_dir)//trim(input%prefix)//'.save/'
       !save_dir=trim(tmp_dir)//trim(input%prefix)//'.save/'
       iunbnk=find_free_unit()
       filename = '.bnk' // TRIM(int_to_char(iglobal_k))
       open( unit= iunbnk, file=trim(tmp_dir)//trim(input%prefix)//'.save/'//trim(filename), &
            &status='unknown',form='unformatted')
       write(iunbnk) sh%num_bands
       write(iunbnk) eig%wave_func(1:sh%ntot_e,1:sh%num_bands)
       close(iunbnk)
     endif
     !Determine how many bands to save
     nbnd=sh%num_bands

     IF (input%lshrink) THEN  !SJ
              IF ( input%nband_save < 1 ) THEN
                      write(stdout,*) 'Invalid or unspecified nband_save. Will save the OBF to the size of DFT bands'
                      nbnd_save = sh%num_val + sh%num_cond
              ELSE
                      nbnd_save = input%nband_save
              ENDIF
         if (nbnd_save > sh%ntot_e) THEN !Check if the requested nband_save is larger than the maximum possible
                  write(stdout,*) 'The requested final number of bands is larger than the number of band in OBF basis'
                  write(stdout,*) 'nbnd_save',nbnd_save,'nbasis OBF',sh%ntot_e
                  write(stdout,*) 'Will save to the size of nbasis OBF instead'
                  nbnd_save = sh%ntot_e     ! If yes then we reduce that to maximum possible.
         ENDIF
                  !nbnd_save = sh%num_val + sh%num_cond
     ELSE
                  nbnd_save = nbnd
     ENDIF
     iuneig=find_free_unit()
     filename= TRIM(int_to_char(iglobal_k))//'.txt'
     !write(1000+mpime,*) save_dir,tmp_dir,filename
     open(unit=iuneig, file=trim(save_dir)//'eig'//trim(filename), &
             &status='unknown',form='formatted')
     write(iuneig,998) energy%energy(1:nbnd_save,ik)
     998 format(*(f25.16))
     close(iuneig)
     !do  im = 1, sh%num_bands
     !    write(stdout,*) eig%wave_func(:,im)
     !enddo
  !SJ
  IF (input%build_unk) THEN
          WRITE(stdout,*) 'Building u_nk from diagonalized matrix and wfc from obf_basis.x'

          CALL read_igwx(iglobal_k,input,npol,nbasis,igwx) 
          WRITE(stdout,*) 'ik igwx npol nbasis sh%num_bands' !ik should be k-point index in each proc?
          WRITE(stdout,*) iglobal_k,igwx,npol,nbasis,sh%num_bands
          ALLOCATE(wfc(npol*igwx,nbasis))
          wfc=(0.d0,0.d0)
          CALL read_wfc(iglobal_k,input,npol,nbasis,igwx,wfc) !Read wfc(ik).dat from simple
          WRITE(stdout,*) 'Done reading wfc at ik',ik
          ALLOCATE(unk(npol*igwx,nbnd))
          unk=(0.d0,0.d0)
          CALL unk_from_wfc(eig%wave_func,sh,wfc,npol,nbasis,nbnd,igwx,unk)
          DEALLOCATE(wfc)
          CALL write_wfc(iglobal_k,input,npol,nbnd,igwx,tpiovera*q,unk,nbnd_save,1)
          !CALL write_wfc(iglobal_k,input,npol,nbnd,igwx,tpiovera*q,unk) 
          DEALLOCATE(unk)
  ENDIF
  enddo
  close(iuneig)
  !IF (input%kpoint_manual) THEN 
          WRITE(stdout,*) 'Doing only band diagonalization'
          return 
  !ENDIF !SJ
  !!!!!!! End interpolation of the bands
  !
  allocate(energy_sum(sh%num_bands,kpnts%nk))
  energy_sum = 0.d0
  ! Gather all the bands divided in the various processors by subgroups of k-points
  do ik =1, energy%nk_loc
    energy_sum(1:sh%num_bands,ik+energy%ik_first-1) = energy%energy(1:sh%num_bands, ik)
  enddo
  call mp_sum(energy_sum,world_comm)
  !
  !!!!!!! Fermi energy calculation
  if (input%fermi_energy == -1) then
    write(stdout,*) ' '
    write(stdout,*) 'Computing the Fermi energy...'
    nelec = sh%nelec
    allocate(wk(1:kpnts%nk),isk(1:kpnts%nk))  ! kpnts%nk = total number of k-points in interp_grid
    wk(1:kpnts%nk) = 2.d0/(sh%npol*dble(kpnts%nk))   ! uniform weights
    isk(1:kpnts%nk) = 1
    is = 0
    !
    if (input%tetrahedron_method) then
        ! Tetrahedra
        num_tetra = 6 * kpnts%nk
        allocate(tetra(4,num_tetra))
        allocate(weights(kpnts%nk))
        call tetrahedra1(kpnts%nkgrid(1), kpnts%nkgrid(2), kpnts%nkgrid(3), num_tetra, tetra)
        ef = efermit (energy_sum, sh%num_bands, kpnts%nk, nelec, sh%nspin, num_tetra, tetra, is, isk)
        write(stdout,'(a,f10.5,a)') ' Fermi energy (tetrahedra) = ', ef*rytoev , ' eV'
    else
        ! Broadening methods
        ef = efermig (energy_sum, sh%num_bands, kpnts%nk, nelec, wk, input%fermi_degauss, input%fermi_ngauss, is, isk)
        write(stdout,'(a,f10.5,a)') ' Fermi energy (broadening method) = ', ef*rytoev , ' eV'
    endif
    deallocate(wk,isk)
  !
  else
  !
    write(stdout,*) 'Reading the Fermi energy from input...'
    ef = input%fermi_energy
    write(stdout,*) 'Fermi energy = ', ef*rytoev
  endif
  !!!!!!! End Fermi energy calculation
  !
  ! Occupation of the states (with Fermi-Dirac distribution)
  focc = 0.0d0
  do ik=1,energy%nk_loc
    do ii=1,sh%num_bands
      arg = (ef - energy%energy(ii,ik))/input%elec_temp
      focc_tmp(ii,ik) = wgauss(arg,-99)
    enddo
  enddo
  do ik=1,energy%nk_loc
    focc(1:sh%num_bands,ik+energy%ik_first-1) = focc_tmp(1:sh%num_bands,ik)
  enddo
  call mp_sum(focc,world_comm) 
  !
  !!! Calculate DOS  (units: states/eV)
  headline(6) = "       Energy grid [eV]               DOS"
  num_energy_dos = int( ( maxval(energy_sum(sh%num_bands,:)) - minval(energy_sum(1,:)) )&
                    & /input%delta_energy_dos + 0.5) + 1
  allocate(wgrid_dos(1:num_energy_dos),dos(1:num_energy_dos))
  wgrid_dos = 0.0d0
  do iw = 1, num_energy_dos
      wgrid_dos(iw) = minval(energy_sum(1,:)) + (iw-1) * input%delta_energy_dos ! energy grid for DOS
  enddo
  ! 
  dos = 0.d0
  if (input%tetrahedron_method) then
    ! Calculate DOS (Tetrahedron method)
    write(stdout,*) 'Computing the DOS (Tetrahedron method)...'
    do iw = 1, num_energy_dos
      w = wgrid_dos(iw)
      do ii=1, sh%num_bands
        call weights_delta1(w, energy_sum(ii,:), num_tetra, tetra, kpnts%nk, weights)
        do ik=1,kpnts%nk
          dos(iw) = dos(iw) + weights(ik)
        enddo
      enddo
    enddo
    dos(1:num_energy_dos) = 2.d0 * dos(1:num_energy_dos) / sh%npol
    if (ionode) then
      write(stdout,"(/,1x, 'Writing output on file...' )")
      call writetofile(input%prefix,"dos_tetrahedra",headline(6),num_energy_dos,wgrid_dos,1,dos/rytoev)
    endif
  else
    ! Calculate DOS (Broadening methods)
    write(stdout,*) ' '
    write(stdout,*) 'Computing the DOS (Broadening method)...'
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ik,ii,iw,w,arg) &
    !$OMP FIRSTPRIVATE(eig) REDUCTION(+:dos)
    do ik=1,energy%nk_loc
        do ii=1, sh%num_bands
           do iw = 1, num_energy_dos
             w = wgrid_dos(iw)
             !! Approximate Dirac-delta with the derivative of the Fermi-Dirac function
             arg = (energy%energy(ii, ik) - w)/input%inter_broadening
             dos(iw) = dos(iw) + w0gauss(arg,input%drude_ngauss)/input%inter_broadening
           enddo
         enddo
    enddo
    !$OMP END PARALLEL DO
    call mp_sum(dos,world_comm)
    dos(1:num_energy_dos) = 2.d0*dos(1:num_energy_dos) / (sh%npol*dble(kpnts%nk))
    if (ionode) then
        write(stdout,"(/,1x, 'Writing output on file...' )")
        call writetofile(input%prefix,"dos_broadening",headline(6),num_energy_dos,wgrid_dos,1,dos/rytoev)
    endif
  endif
  deallocate(wgrid_dos,dos,energy_sum)
  !!! End calculation of DOS
  !
  !$OMP END PARALLEL DO
  !
  !
end subroutine dielectric






