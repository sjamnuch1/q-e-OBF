!this subroutine creates a optimal basis for the periodic functions |u_{nk}>
MODULE f_basis

        USE kinds, ONLY : DP

        COMPLEX(DP), ALLOCATABLE, TARGET :: &
        wfc0(:,:)
contains
subroutine f_basis_obf
  USE kinds, ONLY : DP
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc,tmp_dir
  USE wavefunctions, ONLY : evc,psic
  USE wvfct, ONLY : nbnd,npw,npwx,et 
  USE mp, ONLY : mp_sum,mp_barrier,mp_max
  USE klist, ONLY : nks,ngk,xk, igk_k
  USE becmod,        ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE mp_world, ONLY : world_comm,mpime
  USE ions_base,  ONLY : nat, nsp, ityp
  USE io_global, ONLY : stdout, ionode
  USE input_basis
  USE wannier_gw, ONLY : num_nbndv
  USE fft_base,         ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, gstart,gg, g
  

  implicit none


  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:),omatj(:),omatij(:,:),dmat(:,:)
  INTEGER ::  ipol,ierr
  TYPE(BEC_TYPE), ALLOCATABLE :: bec0(:),bec1(:),bec2(:)
  INTEGER :: ijkb0,nt,na,jh,ih,i,j,ikb,jkb
  INTEGER :: ik,ig,nfound,nfound_gs
  INTEGER :: npw1,ntot,isca, npw_work
  COMPLEX(kind=DP), ALLOCATABLE :: wfc1(:,:),wfc3(:,:),wfc2(:,:) !wfc0(:,:),wfc1(:,:),wfc3(:,:)
  !COMPLEX(DP), ALLOCATABLE, TARGET :: &
  !wfc0(:,:)
       !! wavefunctions in the PW basis set.
       !! noncolinear case: first index is a combined PW + spin index
       !! same format as wfc0
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: ZDOTC

  COMPLEX(kind=DP), ALLOCATABLE :: wfc_t(:,:)
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: iun,iloop,nbasis
  COMPLEX(kind=DP), ALLOCATABLE :: rwfc(:,:), cwfc(:),rprod(:),prod_g(:,:)
  real(kind=DP), allocatable :: grids(:),gridd(:)
  INTEGER :: ii,jj
  COMPLEX(kind=DP), ALLOCATABLE :: valbec(:),zmat(:,:),valbecj(:,:,:)
  LOGICAL :: debug=.false.
  INTEGER :: ikk, jkk, kk

  COMPLEX(kind=DP), ALLOCATABLE :: wfc_e_tmp1(:),wfc_e_tmp2(:)
  TYPE(BEC_TYPE), ALLOCATABLE :: bec_e_tmp1(:),bec_e_tmp2(:)

  !Diagonalization variables
  REAL(kind=DP), DIMENSION(:), ALLOCATABLE :: rwork,energies_tmp
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwork , ifail
  INTEGER :: lwork,low_eig,high_eig,info,nwork
  COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: work
  REAL(kind=DP) :: abstol , vl , vu,trace,trace_ik,ctrace,rnorm
  COMPLEX(kind=DP) :: norm

  INTEGER :: iq(3)
  LOGICAL :: lnorm = .false.
  call start_clock('optimal_basis')
  
!determine npw_max
  if(nks>1) then
     !rewind (unit = iunigk)
     !READ( iunigk ) igk
     npw = ngk (1)
     npw_max=maxval(igk_k(1:npw,1))
     do ik=2,nks
        !READ( iunigk ) igk
        npw = ngk (ik)
        isca=maxval(igk_k(1:npw,ik))
        if(isca>npw_max) npw_max=isca
     end do
  else
     npw = ngk (1)
     npw_max=maxval(igk_k(1:npw,1))
     
  endif
  write(stdout,*) 'nks',nks
  write(stdout,*) 'NPWX NPW_MAX', npwx, npw_max

 ! Try forcing NPW_MAX to be the same across all proc !SJ
  CALL mp_max(npw_max,world_comm)  !SJ
  write(stdout,*) 'NPWX NPW_MAX after mp_max', npwx, npw_max !SJ
  

 ! if (nks>1) rewind (unit = iunigk)
   npw = ngk (1)
   !IF ( nks > 1 ) READ( iunigk ) igk
   !SJ size of evc
   WRITE(stdout,*) 'Size of evc',size(evc,1),size(evc,2),size(evc)
   !call davcio (evc, 2*nwordwfc, iunwfc, 1, -1) !read evc
      

!allocate basis 

   ntot=(num_val+num_cond)*nks !Load the whole thing into the array
   ntot_e=(num_val+num_cond)  
   !SJ Basis variables initialization
   IF ( size(evc,2) < ntot_e) THEN
           WRITE(stdout,*) 'Requested bands for OBF are greater than one read in from DFT stage'
           WRITE(stdout,*) 'You should rerun DFT stage with larger number of bands if this is the expected input bands'
           WRITE(stdout,*) 'nbnd(from DFT): ',size(evc,2), 'nobf(requested): ', ntot_e
           STOP
   ENDIF

   WRITE(stdout,*) 'npol npwx ntot_e ntot'
   WRITE(stdout,*) npol,npwx,ntot_e, ntot
   WRITE(stdout,*) 'Using s_bands threshold of ',s_bands
   !SJ
   !allocate(wfc_t(npol*npwx,ntot_e))
   

!put k=1 wfcs on basis
   !wfc_t=(0.d0,0.d0)
   !do ipol=0,npol-1
   !   wfc_t(1+ipol*npwx:npw+ipol*npwx,1:ntot_e)=evc(1+ipol*npwx:npw+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
   !enddo
   
!check orthonormality
!   allocate(dmat(ntot,ntot),zmat(ntot,ntot))
!   call ZGEMM('C','N',ntot,ntot,npol*npwx,(1.d0,0.d0),wfc_t,npol*npwx,wfc_t,npol*npwx,(0.d0,0.d0),dmat,ntot)
!   call mp_sum(dmat,world_comm)

   
   !Diagonalize the S matrix
 !        write(stdout,*) 'Diagonalization of Smatrix: at ik 1'
 !        allocate(rwork(7*ntot),iwork(5*ntot),ifail(ntot))
 !        allocate(work(1))
 !        allocate(energies_tmp(ntot))
 !        lwork = -1
 !        abstol = -1.d0
 !        vl = 0.d0
  !       vu = 0.d0
  !       low_eig = 1
   !      high_eig = ntot

    !     CALL ZHEEVX('V','I', 'U', ntot, -dmat, ntot, vl, vu, low_eig, high_eig, abstol, &
    !         & nwork, energies_tmp, zmat, ntot, work, lwork, rwork, iwork, ifail, info)

     !    lwork = int(work(1))
      !   deallocate(work)
       !  allocate(work(1:lwork))
  !
        ! CALL ZHEEVX('V','I', 'U', ntot, -dmat, ntot, vl, vu, low_eig, high_eig, abstol, &
!             & nwork, energies_tmp, zmat, ntot, work, lwork, rwork, iwork, ifail, info)

!       write(stdout,*) 'Diagonalized overlap matrix at ik 1'
!       do i=1,ntot
!       write(stdout,*) 'ik: i eigenvalue', "1",i,-energies_tmp(i)
!       enddo
!       deallocate(rwork,zmat,dmat,iwork,ifail,work,energies_tmp)
   
!   if(ionode) then
!      do i=1,ntot_e
!         !do j=1,ntot_e
!            write(stdout,*) 'K POINT 1, ORTHONORMALITY CHECK: ', i, i, omat(i,i)
         !enddo
!      end do
!   endif
!   deallocate(omat)

   allocate (wfc0(npw_max*npol,ntot))!,wfc1(npw_max*npol,ntot_e))
     wfc0=(0.d0,0.d0)
     !iloop=2
     npw_work=npw_max*npol
     
   do ik=1,nks
      call start_clock('wfc_loop')
      call davcio (evc, 2*nwordwfc, iunwfc, ik, -1) !read evc at ik
      npw1 = ngk (ik)
      write(stdout,*) 'ik npw1',ik,npw1
   do ipol=0,npol-1
   do ig=1,npw1
           wfc0(igk_k(ig,ik)+ipol*npw_max,1+(ik-1)*ntot_e:ntot_e+(ik-1)*ntot_e)=&
                  &  evc(ig+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
   enddo ! npw
   enddo ! ipol
   enddo ! ik read evc loop
     ik = nks
      !SJ
     ! WRITE(stdout,*) 'ik: Size of wfc1',ik, size(wfc1,1),size(wfc1,2),size(wfc1)
      WRITE(stdout,*) 'ik: Size of wfc0',ik, size(wfc0,1),size(wfc0,2),size(wfc0)
      !WRITE(stdout,*) 'npw_work',npw_work
      !WRITE(mpime+11000,*) 'Size of wfc1', size(wfc1,1),size(wfc1,2)
      !flush(mpime+11000)

     
!we need the bec factors only at k-point ik

    
     
       
!!allocate  overlap matrix
       allocate(omat(ntot,ntot))

!!at this point wfc_e are expressed by k-states up to ik-1

!!calculate overlap
 !      call start_clock('zgemm')
 !      call ZGEMM('C','N',ntot,ntot,npw_work,(1.d0,0.d0),wfc0,npw_work,wfc0,npw_work,(0.d0,0.d0),omat,ntot)
       !call ZGEMM('C','N',ntot_e,ntot,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
 !      call stop_clock('zgemm')
 !      call mp_sum(omat,world_comm)

 !      call mp_barrier(world_comm)
       !if(ionode.and.debug) then
 !         do i=1,ntot
 !            do j=1,ntot
 !               write(stdout,*) 'K POINT I J, ORTHONORMALITY CHECK: ', i,j, omat(i,j)
 !            enddo
 !         enddo
       !endif

       
!!project out
  !     call start_clock('zgemm')
  !     call ZGEMM('N','N',npw_work,ntot,ntot,(-1.d0,0.d0),wfc0,npw_work,omat,ntot,(1.d0,0.d0),wfc0,npw_work)
       !call ZGEMM('N','N',npw_max*npol,ntot,ntot_e,(-1.d0,0.d0),wfc0,npw_max*npol,omat,ntot_e,(1.d0,0.d0),wfc1,npw_max*npol)
  !     call stop_clock('zgemm')
!updates bec1
!loop on k-points up to ik

!DEBUG part check if they are really orthogonal
       !!caluclate overlap
  !     if(debug) then
  !        call start_clock('zgemm')
  !        call ZGEMM('C','N',ntot,ntot,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc0,npol*npw_max,(0.d0,0.d0),omat,ntot)
  !        call stop_clock('zgemm')
  !        call mp_sum(omat,world_comm)



  !        if(ionode.and.debug) then
  !           do i=1,ntot_e
  !              do j=1,ntot
  !                 write(stdout,*) 'K POINT I J, ORTHONORMALITY CHECK2: ', i,j, omat(i,j)
  !              enddo
  !           enddo
  !        endif !print loop
  !     endif!debug

!Check S matrix
       allocate(dmat(ntot,ntot),zmat(ntot,ntot))
       !write(8000+mpime,*) 'ntot',ntot
       !flush(8000+mpime)
       !call ZGEMM('C','N',ntot,ntot,npol*npw_max,(1.d0,0.d0),wfc1,npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),dmat,ntot)
       call ZGEMM('C','N',ntot,ntot,npw_work,(1.d0,0.d0),wfc0,npw_work,wfc0,npw_work,(0.d0,0.d0),dmat,ntot)
       !call ZGEMM('C','N',ntot,ntot,npw_work,(-1.d0,0.d0),wfc0,npw_work,wfc0,npw_work,(0.d0,0.d0),dmat,ntot)
       call mp_sum(dmat,world_comm)
       !write(stdout,*) 'Reached barrier after mp_sum(dmat)'
       call mp_barrier(world_comm)
       !write(stdout,*) 'Passed barrier after mp_sum(dmat)'

       !do i=1,ntot_e
       !write(1111+mpime,*) dmat(1,1)
       !flush(1111+mpime)
       !enddo

       !Diagonalize the S matrix
         write(stdout,*) 'Diagonalization of Smatrix:',ik
         allocate(rwork(7*ntot),iwork(5*ntot),ifail(ntot))
         allocate(work(1))
         allocate(energies_tmp(ntot))
         lwork = -1
         abstol = -1.d0
         vl = 0.d0
         vu = 0.d0
         low_eig = 1
         high_eig = ntot

         CALL ZHEEVX('V','I', 'U', ntot, -dmat, ntot, vl, vu, low_eig, high_eig, abstol, &
             & nwork, energies_tmp, zmat, ntot, work, lwork, rwork, iwork, ifail, info)

         lwork = int(work(1))
         deallocate(work)
         allocate(work(1:lwork))

         !write(9000+mpime,*) 'lwork',lwork
         !write(9000+mpime,*) 'ntot',ntot
         !write(9000+mpime,*) 'Size dmat',size(dmat,1),size(dmat,2)
         !flush(9000+mpime)
  !
         CALL ZHEEVX('V','I', 'U', ntot, -dmat, ntot, vl, vu, low_eig, high_eig, abstol, &
             & nwork, energies_tmp, zmat, ntot, work, lwork, rwork, iwork, ifail, info)

       do i=1,ntot
       write(stdout,*) 'ik: i eigenvalue', ik,i,-energies_tmp(i)
       enddo

       write(stdout,*) 'Done with printing eigenvalues at ', ik

       !Truncating the basis based on value of s_bands
       nfound = 0
       energies_tmp=energies_tmp**2
       trace = sum(energies_tmp) !Current ik trace
       !IF ( ik == 2 ) THEN !First one is from ntot - trace
       !    trace_ik = ntot-trace     !Remaining trace ~ (ideally this should sum up to ntot)
       !ELSE
       !    trace_ik = trace_ik -trace ! For anything after we subtract current trace out.
       !ENDIF
                                 
       write(stdout,*) 'Smatrix trace at ik',ik,trace
       write(stdout,*) 'Smatrix trace diff at ik', ik,ntot-trace
       !write(stdout,*) 'Smatrix trace difference at ik',ik,trace_ik

       !ctrace= 0.d0 ! ntot - trace
       !trace = ntot
       !do i=1,ntot
       !  ctrace = ctrace + energies_tmp(i)  ! 
       !  write(stdout,*) 'Percent coverage based on current trace',i,ctrace,ctrace/trace*100
       !enddo

       write(stdout,*) ' '

       !ctrace= 0.d0 ! ntot - trace
       !do i=1,ntot
       ! ctrace = ctrace + energies_tmp(i)  !
       !  write(stdout,*) 'Percent coverage based on remaining trace',i,ctrace,ctrace/(trace_ik)*100
       !enddo

       IF (trace_mode == 0) THEN

       ctrace= ntot-trace ! ntot - trace
       do i=1,ntot
        ctrace = ctrace + energies_tmp(i)  !
         write(stdout,*) 'Percent coverage based on difference band and trace',i,ctrace,ctrace/(ntot)*100
       enddo

       write(stdout,*) 'Basis based on ntot-trace'
       if ( s_bands > 0 ) then
               ctrace  = ntot-trace
               do i=1,ntot
                ctrace =ctrace + energies_tmp(i)
                 if ( abs(ctrace/(ntot)) > 1.d0-s_bands ) exit
               enddo
               nfound = min(i,ntot)
               write(stdout,*) 'Found ', nfound, 'basis functions at ik ',ik
               write(stdout,*) 'Trace at this point', ctrace
               write(stdout,*) 'Discarded trace', ntot-ctrace
       endif
       
       ELSEIF (trace_mode == 1 ) THEN

       ctrace= 0.d0 ! ntot - trace
       do i=1,ntot
         ctrace = ctrace + energies_tmp(i)  !
         write(stdout,*) 'Percent coverage based on current trace',i,ctrace,ctrace/trace*100
       enddo

       write(stdout,*) 'Basis based on current trace'
       if ( s_bands > 0 ) then
               ctrace  =0.d0
               do i=1,ntot
                ctrace =ctrace + energies_tmp(i)
                 if ( abs(ctrace/trace) > 1.d0 - s_bands ) exit
               enddo
               nfound = min(i,ntot)
               write(stdout,*) 'Found ', nfound, 'basis functions at ik ',ik
       endif

       ELSEIF (trace_mode == 2) THEN
          WRITE(stdout,*) 'Keep basis based on fraction of all input bands'
          nfound = s_bands * ntot
          write(stdout,*) 'Found ', nfound, 'basis functions at ik ',ik
       ELSE
          WRITE(stdout,*) 'trace_mode is', calc_mode
          WRITE(stdout,*) 'trace_mode 0 for ntot-trace and 1 for current trace'
          WRITE(stdout,*) 'Incorrect trace_mode. Exiting'
          STOP
       ENDIF


       !write(stdout,*) 'Basis based on trace_ik'
       !Do based on trace_ik
       !if ( s_bands > 0 ) then
       !        ctrace  =0.d0
       !        do i=1,ntot
       !         ctrace =ctrace + energies_tmp(i)
       !          if ( abs(ctrace/trace_ik) > 1.d0 - s_bands ) exit
       !        enddo
       !        nfound = min(i,ntot)
       !        write(stdout,*) 'Found ', nfound, 'basis functions at ik ',ik
       !endif


       !Get wfc1 in OBF basis from diagonalized matrix zmat

       !allocate(wfc3(npw_max*npol,nfound))

       !write(8000+mpime,*) 'npw_work',npw_work
       !write(8000+mpime,*) 'nfound ',nfound

       

       allocate(wfc3(npw_work,nfound))

       CALL MPI_BARRIER(world_comm,ierr)

       write(stdout,*) 'Finished with building basis at', ik ,'Performing GS.'
       !write(8000+mpime,*) 'Finished with building basis at', ik, 'Performing GS.'
       !flush(8000+mpime)

       !project zmat out in wfc1 to get obf basis
       !zmat is required up until 1->nfound. So we get wfc1[npw,nbnd] * zmat[nbnd,nfound] = obf[npw,nfound] 
       !obf is then added to the existing basis
       !call ZGEMM('N','N',npw_max*npol,nfound,ntot,(1.d0,0.d0),wfc1,npw_max*npol,zmat(:,1:nfound),&
       !        & ntot,(0.d0,0.d0),wfc3,npw_max*npol)
       call ZGEMM('N','N',npw_work,nfound,ntot,(1.d0,0.d0),wfc0,npw_work,zmat(:,1:nfound),&
               & ntot,(0.d0,0.d0),wfc3,npw_work)

       !Using GS to keep only orthonormal one.
       IF ( gsthres > 0 ) THEN
         !IF ( eigval_as_thres ) THEN ! Use threshold based on the last kept eigenvalue from diagonalized Smat
              ! call optimal_gram_schmidt_sm(nfound,wfc3,1,-energies_tmp(nfound),ik,nfound_gs)   
         !     call optimal_gram_schmidt_sm(nfound,wfc3,1,-energies_tmp(nfound),ik,npw_work,nfound_gs)
         !ELSEIF ( no_gs ) THEN
         !IF ( no_gs) THEN
         !        WRITE(stdout,*) 'Skipping Gram-schmidt orthonormal'
         !        nfound_gs = nfound
              call optimal_gram_schmidt_sm(nfound,wfc3,1,gsthres              ,ik,npw_work,nfound_gs)
       ELSE   !Use input specified threshold 
              !call optimal_gram_schmidt_sm(nfound,wfc3,1,gsthres,ik,nfound_gs)   
              nfound_gs = nfound
              !call optimal_gram_schmidt_sm(nfound,wfc3,1,gsthres              ,ik,npw_work,nfound_gs)
         !ENDIF
       ENDIF ! GS final loop if needed

       !nfound_gs = nfound
       write(stdout,*) 'Keep ',nfound_gs, 'from',nfound,'basis functions'
       deallocate(wfc0)
       allocate(wfc0(npw_work,nfound_gs))
       wfc0=(0.d0,0.d0)
       wfc0(1:npw_work,1:nfound_gs)=wfc3(1:npw_work,1:nfound_gs)
       !Now add the basis into the existing OBF.
       !Create tmp wfc2 for storing old obf + new obf
       !allocate(wfc2(npw_max*npol,ntot_e+nfound_gs))
       !allocate(wfc2(npw_work,nfound_gs))
       !wfc2(1:npw_max*npol,1:ntot_e)=wfc0(1:npw_max*npol,1:ntot_e)
       !wfc2(1:npw_max*npol,ntot_e+1:ntot_e+nfound_gs)=wfc3(1:npw_max*npol,1:nfound_gs)
       !wfc2(1:npw_work,1:ntot_e)=wfc0(1:npw_work,1:nfound_gs)
       !wfc2(1:npw_work,ntot_e+1:ntot_e+nfound_gs)=wfc3(1:npw_work,1:nfound_gs)
       !deallocate(wfc3)
       !ntot_e =ntot_e+nfound_gs
       !deallocate(wfc0)
       !allocate(wfc0(npw_max*npol,ntot_e))
       !allocate(wfc0(npw_work,ntot_e))
       !wfc0=(0.d0,0.d0)
       !wfc0(1:npw_max*npol,1:ntot_e)=wfc2(1:npw_max*npol,1:ntot_e)
       !wfc0(1:npw_work,1:ntot_e)=wfc2(1:npw_work,1:ntot_e)
       !deallocate(wfc2)
       ntot_e=nfound_gs
       !Normalize the basis 1/sqrt(norm)
       do i=1,ntot_e
        rnorm = zdotc(npw_work,wfc0(1,i),1,wfc0(1,i),1)
        !rnorm = zdotc(npw_max*npol,wfc0(1,i),1,wfc0(1,i),1)
        call mp_sum(rnorm,world_comm)
        rnorm=dsqrt(dble(rnorm))
        !wfc0(1:npw_max*npol,i)=wfc0(1:npw_max*npol,i)/rnorm
        wfc0(1:npw_work,i)=wfc0(1:npw_work,i)/rnorm
       enddo

       deallocate(dmat,zmat, rwork, work, iwork, ifail, energies_tmp)

!DEBUG part check if they are really orthogonal
       !!caluclate overlap
      ! if(debug) then
      !    call start_clock('zgemm')
      !    call ZGEMM('C','N',ntot_e,ntot,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
      !    call stop_clock('zgemm')
      !    call mp_sum(omat,world_comm)
         
          
          
      !    if(ionode.and.debug) then
      !       do i=1,ntot_e
      !          do j=1,ntot
      !             write(stdout,*) 'K POINT I J, ORTHONORMALITY CHECK2: ', i,j, omat(i,j)
      !          enddo
      !       enddo
     !     endif
     !  endif!debug
       
!!orhtonormalize them

!!add to basis and updates arrays and counters
       deallocate(omat)


       !call start_clock('wfc_optimal')
       !call optimal_gram_schmidt_z(ntot,wfc1,1,s_bands,ik,nfound)
       
       !call stop_clock('wfc_optimal')

!do the same copying for the bec factors
        write(stdout,*) 'DIMENSION OF BASIS at ik', ntot_e, ik !SJ nband new
!check orthonormality of basis

      if(debug) then
         write(stdout,*) 'ORTHONORMALITY check of basis: ', ik
          allocate(omat(ntot_e,ntot_e))
          call ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc0,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
          call mp_sum(omat,world_comm)
     
          do ii=1,ntot_e
             do jj=1,ntot_e
                write(stdout,*) 'ORTHONORMALITY of basis (ii,jj):', ii,jj, omat(ii,jj)
             enddo
          enddo
          deallocate(omat)
      endif


!recalculate bec's
      !deallocate(omat)
       write(stdout,*) 'DIMENSION OF BASIS', ntot_e !SJ nband new
       call stop_clock('wfc_loop')
       call print_clock('wfc_loop')
       call print_clock('wfc_optimal')
       call print_clock('zgemm')
  !end !ik loop
      !if(debug) then
      !Check the orthonormality of the basis
         write(stdout,*) 'ORTHONORMALITY check of basis. Attempting to print out large norm.'
         write(stdout,*) 'Consider increasing gsthres if basis is not orthonormal'
          allocate(omat(ntot_e,ntot_e))
          !call ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max,(1.d0,0.d0),wfc0,npol*npw_max,wfc0,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
          call ZGEMM('C','N',ntot_e,ntot_e,npw_work,(1.d0,0.d0),wfc0,npw_work,wfc0,npw_work,(0.d0,0.d0),omat,ntot_e)
          call mp_sum(omat,world_comm)

          do ii=1,ntot_e
             do jj=1,ntot_e
                if  ( ii /= jj ) THEN
                        if ( ABS(omat(ii,jj)) > 0.25 ) THEN
                        write(stdout,*) 'WARNING'
                        write(stdout,*) 'ORTHONORMALITY of basis (ii,jj):', ii,jj, omat(ii,jj)
                        lnorm = .true.
                        ENDIF
                ENDIF
             enddo
          enddo
          deallocate(omat)
          IF (lnorm) THEN
                  write(stdout,*) 'WARNING. Detected large norm. Expected repeated bands'
          ENDIF
      !endif

!copy results to common variable (sigh..)
     !IF ( allocated( wfc_e ) )      DEALLOCATE( wfc_e )
     !if( allocated( wfc_e ) ) deallocate( wfc_e )

    IF ( restart_obf) THEN
      deallocate(wfc_e)
    ENDIF

    !The wfc_e is stored in memory and is used when writing OBF to file.
    allocate(wfc_e(npw_work,ntot_e))
    !allocate(wfc_e(npw_max*npol,ntot_e))
    !wfc_e(1:npw_max*npol,1:ntot_e)=wfc0(1:npw_max*npol,1:ntot_e)
    wfc_e(1:npw_work,1:ntot_e)=wfc0(1:npw_work,1:ntot_e)
 

    do ik=1,nks
       write(stdout,*) ' IK', ik
       write(stdout,*) ' xk', xk(1:3,ik)
    enddo!ik
    write(stdout,*) ' '
    write(stdout,*) 'TOTAL NUMBER OF OPTIMAL BASIS VECTORS :', ntot_e
    write(stdout,*) ' '

    
   !SJ Final Optimal basis is stored in wfc0
   !SJ
   WRITE(stdout,*) 'Size of wfc0'
   WRITE(stdout,*) size(wfc0,1),size(wfc0,2),size(wfc0)
   !WRITE(mpime+1000,*) 'Size of wfc0 in s_band_module'
   !WRITE(mpime+1000,*) size(wfc0,1),size(wfc0,2),size(wfc0)
   !flush(mpime+1000)
      WRITE(stdout,*) 'Size of wfc_e'
   WRITE(stdout,*) size(wfc_e,1),size(wfc_e,2),size(wfc_e)
   WRITE(stdout,*) 'Size of evc'
   WRITE(stdout,*) size(evc,1),size(evc,2),size(evc)

   !WRITE wfc0 to file for testing
   !open( unit= 1000, file='wfc0', status='unknown',form='unformatted')
   !WRITE(1000+mpime,*) size(wfc0,1),size(wfc0,2)
   !do ii= 1,ntot_e
   !WRITE(1000+mpime,*) wfc0(1:npw_max*npol,ii)
   !enddo
   !close(1000+mpime)

   !SJ We need to put wfc0 it into evc file so that it can be saved to wfc#.dat
   !This is currently done in modified punch.f90
   !By calling punch at the end of run
   !The write_collected_wfc_simple() is called if we're running simple code (by setting lsimple =.true.)
   !The wfc files are written in wfc_shirley#.dat
   !deallocate(wfc1)
   deallocate(wfc0)!,wfc1)

  call stop_clock('optimal_basis')

  return


end subroutine f_basis_obf

!SUBROUTINE optimal_gram_schmidt_sm(num_in,wfcs,ithres,thres,ik,num_out)
SUBROUTINE optimal_gram_schmidt_sm(num_in,wfcs,ithres,thres,ik,npw_in,num_out)
!this subroutine performs a gram_schmidt orthonormalization and retains
!vectors which are above the given threshold
!For orthogonalizing the current OBF with new one we check if we are going to retain it or not.

  USE kinds,                ONLY : DP
  USE mp_world, ONLY : world_comm, mpime, nproc
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE noncollin_module, ONLY: npol, noncolin
  USE input_basis, ONLY : npw_max,vkb_max
  USE becmod,        ONLY : bec_type,calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE ions_base,  ONLY : nat, nsp, ityp



 implicit none

  INTEGER, INTENT(in) :: num_in!number of initial vectors
  !SJ
  INTEGER, INTENT(in) :: ik !index of ik
  !SJ
  !COMPLEX(kind=DP), INTENT(inout) :: wfcs(npw_max*npol,num_in)!in input non-orthonormal in output optimal basis
  INTEGER, INTENT(in) :: npw_in !npw_work
  COMPLEX(kind=DP), INTENT(inout) :: wfcs(npw_in,num_in)!in input non-orthonormal in output optimal basis
  INTEGER, INTENT(in) :: ithres!kind of threshold
  REAL(kind=DP), INTENT(in) :: thres!thrshold for the optimal basis
  !INTEGER, INTENT(in) :: npw_in !npw_work
  INTEGER, INTENT(out) :: num_out!final number of orthonormal basis functions



  INTEGER :: i,j
  COMPLEX(kind=DP), ALLOCATABLE :: prod(:)
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  REAL(kind=DP) :: sca

  TYPE(BEC_TYPE) :: bec0
  TYPE(BEC_TYPE), ALLOCATABLE :: bec1(:)
  INTEGER :: ijkb0,nt,na,jh,ih,ikb,jkb,ipol
  
  allocate(prod(num_in))
  num_out=0

  do i=1,num_in
     if(num_out >0) then 
        call zgemv('C',npw_max*npol,num_out,(1.d0,0.d0), wfcs,npw_max*npol,wfcs(1,i),1,(0.d0,0.d0),prod,1)
        call mp_sum(prod(1:num_out),world_comm)
        call start_clock('zgemm')
        !call zgemm('N','N',npw_max*npol,1,num_out,(-1.d0,0.d0),wfcs,npw_max*npol,prod,num_in,(1.d0,0.d0),wfcs(1,i),npw_max*npol)
        call zgemm('N','N',npw_in,1,num_out,(-1.d0,0.d0),wfcs,npw_in,prod,num_in,(1.d0,0.d0),wfcs(1,i),npw_in)
        call stop_clock('zgemm')
     endif
     !csca = zdotc(npw_max*npol,wfcs(1,i),1,wfcs(1,i),1)
     csca = zdotc(npw_in,wfcs(1,i),1,wfcs(1,i),1)
     call mp_sum(csca,world_comm)

     write(stdout,*) 'ik: ',ik,'num_in csca',i,dble(csca),dble(csca)-thres



     if(dble(csca) >= thres) then
        num_out=num_out+1
        sca=dsqrt(dble(csca))
        !wfcs(1:npw_max*npol,num_out)=wfcs(1:npw_max*npol,i)/sca
        wfcs(1:npw_in,num_out)=wfcs(1:npw_in,i)/sca
        !wfcs(1:npw_max*npol,num_out)=wfcs(1:npw_max*npol,i)/sca
        wfcs(1:npw_in,num_out)=wfcs(1:npw_in,i)/sca
     endif
  enddo


  deallocate(prod)
  return
END SUBROUTINE optimal_gram_schmidt_sm

end module
