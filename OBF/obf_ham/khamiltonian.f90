subroutine khamiltonian
     !------------------------------------------------------------------------------------------
     !
     !  This subroutine calculates and saves to disk the elements needed for the construction of
     !  the k-dependent Hamiltonian in the optimal basis |e_i>:
     !  K1_ij, K0_ij, K1_ij, Vloc_ij, V_nonloc_ij (only the projectors <\beta_l|exp(ik.r)|e_j>)
     !  It calculates also the nonlocal commutator for the the velocity matrix elements
     !
     USE kinds, ONLY : DP
     USE input_ham
     USE input_parameters,         ONLY : qk, nksinterp_shirley, wk
     USE scf, ONLY : vrs
     USE wavefunctions, ONLY : psic
     USE cell_base, ONLY : tpiba2, tpiba, bg, alat, at
     USE gvect, ONLY : gg, g
     USE io_global, ONLY : stdout, ionode, ionode_id
     USE mp_world, ONLY : world_comm,mpime
     USE mp, ONLY : mp_sum,mp_barrier, mp_bcast
     USE noncollin_module, ONLY: npol, noncolin, domag
     USE io_files,  ONLY : prefix, tmp_dir, iunwfc, nwordwfc
     USE fft_base,         ONLY : dffts !smooth grid fft
     USE fft_interfaces,   ONLY : fwfft, invfft
     USE cell_base, ONLY : omega
     USE ions_base,     ONLY: nat, ntyp => nsp, ityp
     USE uspp_param,    ONLY: nh, nhm
     USE uspp,          ONLY: nkb, deeq, ofsbeta, deeq_nc
     USE lsda_mod,      ONLY : nspin
     USE becmod,        ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
     USE wvfct, ONLY : npwx
     USE klist, ONLY : nelec
     USE wannier_gw, ONLY : num_nbndv
     USE uspp_init,        ONLY : init_us_2
     use wavefunctions,  ONLY : evc
     !
     IMPLICIT NONE
     !
    COMPLEX(kind=DP), ALLOCATABLE :: g2e(:,:) , g2e_mat(:,:) , rwfc(:,:) , v_rwfc(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: deeqc(:,:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: d_wfc_e(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: sum_beck_nc(:,:,:,:) , sum_beckc(:,:,:) 
    COMPLEX(kind=DP) :: sup, sdwn
    INTEGER :: ipol, ig, itot_e, iun, iun1, idir, ii, ir, jj, kk, jtot_e, kindex, sum_kindex
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER, ALLOCATABLE :: igkk(:)
    REAL(kind=DP) :: qk_uniform(3)
    TYPE(bec_type) :: beck, beck2
    !
    call start_clock('khamiltonian')
    !
    write(stdout,*) '  '
    write(stdout,*) 'Compute and store (in files prefix.hamiltonian and prefix.hamiltonian_k) the'
    write(stdout,*) 'matrix elements needed for the k-dependent Hamiltonian in the optimal basis:'
    !
    call davcio (evc, 2*nwordwfc, iunwfc, 1, -1) !Get OBF. Only one file since OBF is reduced to gamma only.
                                                 !Now use evc instead of wfc_e
    !
    WRITE(stdout,*) 'nwordwfc',nwordwfc
    WRITE(stdout,*) 'npw_max',npw_max
    WRITE(stdout,*) 'npwx',npwx
    WRITE(stdout,*) 'ntot_e',ntot_e
    WRITE(stdout,*) 'npol',npol
    WRITE(stdout,*) 'nkb',nkb
    npw_max = npwx
   WRITE(stdout,*) 'npw_max',npw_max

    allocate(g2e(npw_max*npol,ntot_e),g2e_mat(ntot_e,ntot_e),d_wfc_e(npw_max*npol,ntot_e))
    !
    if (noncolin) then
      allocate(sum_beck_nc(nkb,npol,ntot_e,nksinterp_shirley))
    else
       allocate(sum_beckc(nkb,ntot_e,nksinterp_shirley))
    endif
    !
    !
    g2e(1:npw_max*npol,1:ntot_e) = evc(1:npw_max*npol,1:ntot_e)
    !
    ! K0_ij (kinetic part)
    call start_clock('K0')
    write(stdout,*) 'K0_ij'
    do itot_e=1,ntot_e
     do ipol=0,npol-1
        do ig=1,npw_max
           g2e(ig + ipol*npw_max,itot_e) = gg(ig)*g2e(ig + ipol*npw_max,itot_e)
        enddo
     enddo
    enddo
    ! Scalar product: at the end we have \sum_G e_i*(G) G**2 e_j(G)
    call ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max,(1.d0,0.d0),g2e,npol*npw_max,evc,npol*npw_max,(0.d0,0.d0),g2e_mat,ntot_e)
    call mp_sum(g2e_mat,world_comm)
    !
    g2e_mat = g2e_mat*tpiba2
    !
    if(ionode) then
      iun=find_free_unit()
      open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.hamiltonian', &
                &status='unknown',form='unformatted')
      write(iun) ntot_e
    endif
    if (ionode) write(iun) g2e_mat(1:ntot_e,1:ntot_e)
    call stop_clock('K0')
    ! End K0_ij
    !
    ! K1_ij (kinetic part)
    call start_clock('K1')
    write(stdout,*) 'K1_ij'
    do idir=1,3
      g2e(1:npw_max*npol,1:ntot_e) = evc(1:npw_max*npol,1:ntot_e)
      do itot_e=1,ntot_e
       do ipol=0,npol-1
          do ig=1,npw_max
             g2e(ig + ipol*npw_max,itot_e) = g(idir,ig)*g2e(ig + ipol*npw_max,itot_e)
          enddo
       enddo
      enddo
     ! Scalar product: at the end we have \sum_G e_i*(G) G e_j(G)
     call ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max,(1.d0,0.d0),g2e, &
             npol*npw_max,evc,npol*npw_max,(0.d0,0.d0),g2e_mat,ntot_e)
     call mp_sum(g2e_mat,world_comm)
     !
     g2e_mat = g2e_mat*tpiba*2.0 ! In the definition of K1_ij there is a factor 2 in front
     !
     if (ionode) write(iun) g2e_mat(1:ntot_e,1:ntot_e) ! write to disk
     !
    enddo
    call stop_clock('K1')
    ! End K1_ij
    !
    ! Vloc_ij
    call start_clock('Vloc')
    write(stdout,*) 'Vloc_ij'
    allocate(rwfc(dffts%nnr*npol,ntot_e))
    allocate(v_rwfc(dffts%nnr*npol,ntot_e))
    !
    !WRITE(15000+mpime,*) 'dffts%nl(1:npw_max)',dffts%nl(1:npw_max)
    do ii=1,ntot_e
         psic(1:dffts%nnr)=0.d0
         psic(dffts%nl(1:npw_max))=evc(1:npw_max,ii) !dffts%nl is the position of gvec on FFT grid
                                                       !psic is a collection of wfc_e across proc.
         CALL invfft ('Wave', psic, dffts)      ! Inverse fourier transform of |e_i> on the smooth grid
         rwfc(1:dffts%nnr,ii)=psic(1:dffts%nnr) ! Wavefunction in real space

         if(npol>1) then
            psic(1:dffts%nnr)=0.d0
            psic(dffts%nl(1:npw_max))=evc(npw_max+1:npw_max+npw_max,ii)
            CALL invfft ('Wave', psic, dffts)
            rwfc(dffts%nnr+1:2*dffts%nnr,ii)=psic(1:dffts%nnr)
         endif
    enddo
    !
    v_rwfc(1:dffts%nnr*npol,1:ntot_e) = rwfc(1:dffts%nnr*npol,1:ntot_e)

    !do ir=1,dffts%nnr
    !WRITE(10000+mpime,*) vrs(ir,1)
    !enddo

    WRITE(stdout,*) 'noncolin',noncolin
    WRITE(stdout,*) 'domag',domag
    if (noncolin) then
     do ii=1,ntot_e
       if (domag) then
               do ir=1, dffts%nnr
                  sup = v_rwfc(ir,ii) * (vrs(ir,1)+vrs(ir,4)) + &
                        v_rwfc(ir+ dffts%nnr,ii) * (vrs(ir,2)-(0.d0,1.d0)*vrs(ir,3))
                  sdwn =  v_rwfc(ir+ dffts%nnr,ii) * (vrs(ir,1)-vrs(ir,4)) + &
                         v_rwfc(ir,ii) * (vrs(ir,2)+(0.d0,1.d0)*vrs(ir,3))
                  v_rwfc(ir,ii)=sup
                  v_rwfc(ir+ dffts%nnr,ii)=sdwn
               enddo
       else
               do ir=1, dffts%nnr
                  v_rwfc(ir,ii) = v_rwfc(ir,ii) * vrs(ir,1)
                  v_rwfc(dffts%nnr+ir,ii) = v_rwfc(dffts%nnr+ir,ii) * vrs(ir,1)
               enddo
               !v_rwfc(dffts%nnr+1:dffts%nnr*npol,ii) = 0.d0
       endif
     enddo
    else
      do ii=1,ntot_e
        do ir=1,dffts%nnr
          v_rwfc(ir,ii) = v_rwfc(ir,ii)*vrs(ir,1)  ! V_loc(r) x e_i(r)
        enddo
      enddo
    end if


    !do ii=1,ntot_e
    !do ir=1,dffts%nnr
    !write(30000+mpime,*) REAL(rwfc(ir,ii)),AIMAG(rwfc(ir,ii))
    !write(50000+mpime,*) REAL(v_rwfc(ir,ii)),AIMAG(v_rwfc(ir,ii))
    !enddo
    !enddo

    !write(80000+mpime,*) 'nr1 nr2 nr3 nnr',dffts%nr1,dffts%nr2,dffts%nr3,dffts%nnr
    ! Calculate integral of e_j(r)* x V_loc(r) x e_i(r)
    call ZGEMM('C','N',ntot_e,ntot_e,npol*dffts%nnr,(1.d0,0.d0),rwfc, &
            npol*dffts%nnr,v_rwfc,npol*dffts%nnr,(0.d0,0.d0),g2e_mat,ntot_e)
    call mp_sum(g2e_mat,world_comm)
    !
    g2e_mat = g2e_mat / dble(dffts%nr1*dffts%nr2*dffts%nr3)
    !
    !write(stdout,*) 'g2e_mat(1,1)',g2e_mat(1,1)
    if (ionode) write(iun) g2e_mat(1:ntot_e,1:ntot_e)
    !
    deallocate(rwfc,v_rwfc)
    call stop_clock('Vloc')
    ! End Vloc_ij

    ! Vnonloc_ij
    write(stdout,*) 'Vnonloc_ij(k) (projectors)'
    !
    if(ionode) then
      write(iun) noncolin
      write(iun) nat
      write(iun) ntyp
      write(iun) nhm
      write(iun) nspin
      write(iun) nkb
      write(iun) npol
      write(iun) ityp(1:nat)
      write(iun) nh(1:ntyp)
      write(iun) ofsbeta(1:nat)
      write(iun) nksinterp_shirley   !This is the k-point grid. Array with size of 3. nkpoints(3)
      write(iun) nksinterp_shirley !nkp_total so originally nkp(1)*nkp(2)*nkp(3)

      !IF (kpoint_manual) THEN !SJ Also write total number of kpoints to file
      !        write(stdout,*) 'Manual kpoint grid read in will replace total kpoints' !SJ
      !        write(stdout,*) 'with nksinterp_shirley instead of nkpoints(1)*(2)*(3)' !SJ
      !        write(iun) nksinterp_shirley !SJ
      !ELSE !SJ
      !        write(iun) nkpoints(1)*nkpoints(2)*nkpoints(3) !SJ
      !ENDIF !SJ
    endif
    !
    allocate(deeqc(nhm,nhm,nat), vkb_max(npw_max,nkb), igkk(npw_max))
    if (noncolin) then
     if (ionode) write(iun) deeq_nc(1:nhm,1:nhm,1:nat,1:nspin)
    else
     deeqc(1:nhm,1:nhm,1:nat) = cmplx(deeq(1:nhm,1:nhm,1:nat,1), 0.d0, KIND=DP)
     !write(stdout,*) 'deeqc_up(1,1,1)',deeqc(1,1,1)!SJ
     !WRITE(stdout,*) 'deeq(1,1,1,1)',deeq(1:nhm,1:nhm,1:nat,1) !SJ up/dw are the same in norm-conserving PP. 
                                                               !Nonspin dependent
     !DO ii=1,nhm
     !DO jj=1,nhm
     !DO kk=1,nat
     !WRITE(12000+mpime,*) 'deeq(',ii,',',jj,',',kk,')',deeq(ii,jj,kk,1)
     !ENDDO
     !ENDDO
     !ENDDO
     if (ionode) write(iun) deeqc(1:nhm,1:nhm,1:nat)
    endif
    !
    do ig=1,npw_max
      igkk(ig) = ig
    enddo 
    !
    call allocate_bec_type(nkb,ntot_e,beck)
    call allocate_bec_type(nkb,ntot_e,beck2)
    !
    if(ionode) then
      iun1=find_free_unit()
      open( unit= iun1, file=trim(tmp_dir)//trim(prefix)//'.hamiltonian_k', &
                &status='unknown',form='unformatted')
    endif
    !
    !! Calculation of the non-local part of the pseudopotential on a uniform k-mesh
    !SJ
    write(stdout,*) 'Total number of k-points: ', nksinterp_shirley

    do ii=1,nksinterp_shirley
                write(stdout,*) qk(:,ii) !SJ kpoint in crystal coordinate
                !
                call start_clock('Vnloc')
                if (nkb>0) then
                   call init_us_2_max(npw_max,igkk,qk(:,ii),vkb_max)  ! get the projectors \beta_Ilm (k-dependent)
                endif
                !vkb_max(npwx+1:npw_max,1:nkb) = 0.d0 ! WARNING: HERE I PUT TO ZERO THE ELEMENTS OF BETA WiTH G > npwx
                !
                call calbec(npw_max,vkb_max,evc,beck,ntot_e) ! scalar product <\beta_Ilm|exp(ik.r)|e_j>
                !
                if (noncolin) then
                   sum_beck_nc(1:nkb,1:npol,1:ntot_e,ii) = beck%nc(1:nkb,1:npol,1:ntot_e)
                else
                   sum_beckc(1:nkb,1:ntot_e,ii) = beck%k(1:nkb,1:ntot_e)
                endif
                call stop_clock('Vnloc')
                !
                !
           !enddo
             !
             call start_clock('Vnloc_write')
             if (noncolin) then
               if (ionode) write(iun1) sum_beck_nc(1:nkb,1:npol,1:ntot_e,ii)
               if (ionode .and. ii == 1) write(stdout,*) 'Writing sum_beck_nc: Size ', size(sum_beck_nc,1),&
                       size(sum_beck_nc,2),size(sum_beck_nc,3),size(sum_beck_nc,4)
             else
               if (ionode) write(iun1) sum_beckc(1:nkb,1:ntot_e,ii)
               if (ionode .and. ii == 1) write(stdout,*) 'Writing sum_beckc: Size ', size(sum_beckc,1),&
                       size(sum_beckc,2),size(sum_beckc,3)
             endif
             call stop_clock('Vnloc_write')
             !
             !
    enddo

    !
    deallocate(deeqc,g2e_mat,d_wfc_e,g2e)
    if (noncolin) then
      deallocate(sum_beck_nc)
    else
      deallocate(sum_beckc)
    endif
    !
    !
    call deallocate_bec_type(beck)
    call deallocate_bec_type(beck2)
    deallocate(vkb_max,igkk)
    !End Vnonloc_ij
    !
    if(ionode) then
       write(iun) alat
       write(iun) bg(1:3,1:3)
       write(iun) at(1:3,1:3)
       write(iun) nelec
       write(iun) omega
       write(iun) num_val
       write(iun) num_cond
       write(iun) num_nbndv
       !write(iun) nonlocal_commutator
       write(iun) s_bands
    endif
    !
    if(ionode) then
     close(iun)
     close(iun1)
    endif
    !
    write(stdout,*) ' '
    write(stdout,*) '**********************************'
    write(stdout,*) 'Matrix elements computed and saved'
    write(stdout,*) '**********************************'
    !
    call stop_clock('khamiltonian')
    !
    ! Write in output the input parameters
    write(stdout,*) '                     '
    write(stdout,*) 'INPUT PARAMETERS:'
    write(stdout,'(a,f10.5)') ' s_bands [a.u.] = ', s_bands
    write(stdout,'(a,I6)')           ' number of k-point = ', nksinterp_shirley
    write(stdout,*) '                     '
    !
    call print_clock('K0')
    call print_clock('K1')
    call print_clock('Vloc')
    call print_clock('Vnloc')
    call print_clock('Vnloc_write')
    call print_clock('khamiltonian')
    !
end subroutine khamiltonian
