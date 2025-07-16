program obf_diag
  !
  USE start_end
  USE obf_diag_objects
  USE input_obf_diag
  !
  implicit none
  !
  TYPE(input_options_obf_diag) :: din
  TYPE(shirley) :: sh
  TYPE(kpoints) :: k,k_interp
  TYPE(energies) :: e
  !
  call initialize_shirley(sh)
  call initialize_kpoints(k)
  call initialize_energies(e)
  !
  !setup MPI environment
  call startup
  !
  call start_clock('obf_diag')
  call start_clock('init (read)')
  !
  call read_input_obf_diag( din )  
  call read_shirley(din,sh)
  !SJ
  !if (din%kpoint_manual) then
  call read_ham_kgrid(din,k,sh)
  call create_energies_flat(sh,k,e)
  !else
  !        call kgrid_creation(din,k,sh)
  !        call create_energies(sh,k,e)
  !endif
  !

  if (din%nonlocal_interpolation) then
     write(stdout,*) 'Using interpolation for nonlocal part of H'
     call interpolated_kgrid_create(din,k,k_interp,sh)
     call create_energies_flat(sh,k_interp,e) 
     call read_shirley_k_interp(din,sh,e,k,k_interp)  ! trilinear interpolation
  else
   ! if (din%hamil_flat) then !SJ
     call create_energies_flat(sh,k,e)
     if (din%hamilk_byfile) then
             call read_shirley_k_byfile(din,sh,e)
     else
             call read_shirley_k_flat(din,sh,e) !SJ
     endif
   ! else       !SJ 
   !     call read_shirley_k(din,sh,e)
   ! endif !SJ
  endif
  !
  call stop_clock('init (read)')
  !
  ! calculate IP dielectric function (interband + intraband)
  if (din%nonlocal_interpolation) then
          call dielectric(sh,din,k_interp,e)
  else
          call dielectric(sh,din,k,e)
  endif
  !
  !SJ Another loop for spin case. Will now read spin down H and build u_nk down
  call start_clock('down spin')
  IF (sh%nspin == 2) THEN
          call deallocate_shirley(sh)
          call deallocate_energies(e)
          !if (din%kpoint_manual .and. din%hamil_flat) then
          !SJ Read in the down spin hamiltonian.
          call read_shirley_spin(din,sh)
          
          ! interpolation vs no interpolation
          if (din%nonlocal_interpolation) then
              call interpolated_kgrid_create(din,k,k_interp,sh)
              call create_energies_flat(sh,k_interp,e)
              call read_shirley_k_interp(din,sh,e,k,k_interp)  ! trilinear interpolation
          else
              call create_energies_flat(sh,k,e)
                   if (din%hamilk_byfile) then !Same routine since the .hamiltonian_k should be the same for
                                               !spin up/down cases. 
                       call read_shirley_k_byfile(din,sh,e)
                   else
                       call read_shirley_k_flat(din,sh,e) !SJ
                   endif
              !call read_shirley_k_flat(din,sh,e)

              !The usual routine should work since for NCPP k-dependent Hs are the same
              !between up and down spin. The obf_ham will no longer write 2 files for H(k)
              !!call read_shirley_k_flat_spin(din,sh,e)
          endif
          !endif
          call dielectric_spin(sh,din,k,e)
  ENDIF
  call stop_clock('down spin')

  !
  call stop_run
  call deallocate_shirley(sh)
  call deallocate_kpoints(k)
  call deallocate_energies(e)
  !
  call stop_clock('obf_diag')
  call print_clock('init (read)')
  call print_clock('diagonalization')
  call print_clock('diago_vnloc')
  call print_clock('diago_zheevx')
  call print_clock('optic_elements')
  call print_clock('dielectric')
  call print_clock('down_spin')
  call print_clock('obf_diag')
  stop
  !
end program obf_diag

