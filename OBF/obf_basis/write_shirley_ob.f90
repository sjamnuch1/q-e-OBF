SUBROUTINE write_shirley_ob()
!This subroutine write wfc0 to wfc_shirley_ob.dat in QE format.
!First new ecut is determined to get to new npw_max
!Next we will deallocate all fft parameter and reallocate it for new ecut
!The new gspace and miller index mapping are created.
!The collection of miller index and OBF wfc is the same way how QE would normally do.
!Then all the data are written into QE format.

    !USE kinds, ONLY : DP
    USE kinds
    USE parallel_include
    USE constants, ONLY : pi , rytoev
    USE lsda_mod, ONLY : lsda
    USE control_flags, ONLY : smallmem
    USE input_basis
    USE io_global, ONLY : stdout, ionode
    USE ions_base,  ONLY : nsp
    USE io_files,  ONLY : tmp_dir,prefix,psfile
    USE io_files,  ONLY : create_directory
    USE clib_wrappers,        ONLY : f_copy
    !USE wfc_basis,            ONLY : wfc0 !Use wfc_e instead because
                                           !wfc0 is also used in s_band_wfc for Smatrix method
                                           !OBF is stored into wfc_e for both method 
    USE gvecw,                ONLY : gcutw, ecutwfc
    USE pwcom,                ONLY : nspin, nkstot, xk
    USE cell_base, ONLY : bg,at,tpiba2
    USE fft_base,           ONLY : dfftp, dffts
    USE gvect,              ONLY : g, gg, mill, gcutm,gstart,ngm_g,ngm,ig_l2g
    USE gvecs,              ONLY : gcutms, ngms
    USE recvec_subs,        ONLY : ggen, ggens
    use mp_world, ONLY: mpime, world_comm, nproc,root
    USE mp, ONLY : mp_sum,mp_barrier, mp_bcast,mp_max
    USE mp_wave, ONLY : mergekg
    USE fft_types,  ONLY : fft_type_deallocate
    !USE pw_restart_new,       ONLY : pw_write_schema
    USE obf_xml_write,        ONLY : obf_write_schema

    implicit none
    character (len=100)  :: filename,save_dir,writename,read_dir, source,cp_source,cp_target

    !integer, intent(IN) :: iread
    !integer, intent(IN) :: npol
    !integer, intent(IN) :: nbnd
    !integer, intent(IN) :: igwx
    !real(dp), intent(IN) :: xk(3)
    !complex(dp), INTENT(inout) :: unk(npol*igwx,nbnd)
    !complex(DP), INTENT(IN), allocatable :: unk(:,:)
    integer :: igwx_,igwx_gamma,igwx,ii,length


    integer  :: ibnd, iunwfc, iununk
    complex(DP), allocatable :: evc(:,:)
    !integer , allocatable :: mill(:,:)
    real(dp) :: xk_(3)
    integer  :: ik, nbnd_, ispin, npol_, ngw,nbnd ,cp_status
    integer  :: dummy_int, narg,iecutnew,ig,found_new_ecut
    logical  :: gamma_only
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    INTEGER, ALLOCATABLE  :: igk_gamma(:),igk_k_gamma(:),igk_l2g_kdip_gamma(:),mill_tmp(:,:)
    INTEGER, ALLOCATABLE  :: igk_l2g(:),igk_l2g_kdip(:),mill_k(:,:),ig_ip(:),mill_collected(:,:)
    REAL(dp) , ALLOCATABLE :: g2kin_gamma(:)
    COMPLEX(dp), ALLOCATABLE :: wfc0_collect(:,:),pw_ip(:)
    !
    integer :: iuni = 1111, i,ngk_g_gamma,npw_g,npw_tmp,npw_nexti !,npw_gamma
    integer :: itmp, ierr,ip,ngw_ip,ngw_lmax,npw_max_sent
    real(dp) :: scalef,ecutnew,bglen,xk_gamma(3,1)
    real(dp) :: b1(3), b2(3), b3(3), dummy_real

    INTEGER :: istatus(MPI_STATUS_SIZE)

    !Where outputs from pw.x are stored.
    read_dir=trim(tmp_dir)//trim(prefix)//'.save/'
    iunwfc=find_free_unit()
    iununk=find_free_unit()
    WRITE(stdout,*) 'nspin',nspin
    IF ( nspin == 2 ) THEN
            WRITE(stdout,*) 'Spin-polarized calculation found looking for wfcup1.dat for header'
            filename = TRIM(read_dir)//'wfcup1'//'.dat'
    ELSE
            WRITE(stdout,*) 'Non spin-polarized calculation. Looking for wfc1.dat for header'
            filename = TRIM(read_dir)//'wfc1'//'.dat'
    ENDIF
    !filename = TRIM(save_dir)//'wfc' // TRIM(int_to_char(iread)) // '.dat'
    !filename = TRIM(save_dir)//'wfc1'//'.dat' !TODO for ispin != 1 -> then it should read for wfc1up.dat ?
                                              !TODO maybe would be better to fetch this data from data-file-schema.xml? 
    open( unit= iunwfc, file=trim(filename), status='old',form='unformatted')

    !open(UNIT = iunwfc, FILE = trim(filename), FORM = 'unformatted', status = 'old')
    read( iunwfc) ik, xk_, ispin, gamma_only, scalef 
    read (iunwfc) ngw, igwx_, npol_, nbnd_
    read (iunwfc) b1, b2, b3


    !WRITE(stdout,*) 'ik',ik
    !WRITE(stdout,*) 'kpoint:',xk(:)
    !WRITE(stdout,*) 'Number of planewave ngw', ngw
    !WRITE(stdout,*) 'Max number of planewave igwx', igwx
    !WRITE(stdout,*) 'nbnd',nbnd
    !WRITE(stdout,*) 'b1',b1(:)
    !WRITE(stdout,*) 'b2',b2(:)
    !WRITE(stdout,*) 'b3',b3(:)
    !
    ! avoid reading miller indices of G vectors below E_cut for this kpoint
    ! if needed allocate  an integer array of dims (1:3,1:igwx)
    !
    close(iunwfc)

    WRITE(stdout,*) 'Successfully read in wfc.dat for header'

    !ALLOCATION AND INIT BLOCK
    ecutnew = gcutw
    bglen=dot_product(bg(:,1),bg(:,1))
    xk_gamma =0.d0
    npw_gamma = 0
    npw_tmp=npw_max
    !if( allocated(igk_l2g_kdip_gamma) ) deallocate( igk_l2g_kdip_gamma )
    !ALLOCATE ( igk_l2g_kdip_gamma( npw_max+1 ) )
    !if( allocated( igk_gamma ) ) deallocate( igk_gamma )
    !allocate( igk_gamma(npw_max+1) )
    !if( allocated( g2kin_gamma ) ) deallocate( g2kin_gamma )
    !allocate( g2kin_gamma(npw_max+1) )
    !igwx_gamma=0
    !igk_l2g_kdip_gamma = 0
    !g2kin_gamma=0
    !iecutnew=0
    CALL mp_max(npw_tmp,world_comm)
    !WRITE(stdout,*) 'MP_MAX(npw_max)',npw_tmp
    !WRITE(mpime+30,*) 'proc',mpime,'npw_max',npw_max,'ngm',ngm,'npw_tmp',npw_tmp
    !FLUSH(mpime+30)
    !found_new_ecut=0

    !if( allocated(igk_l2g_kdip_gamma) ) deallocate( igk_l2g_kdip_gamma )
    !ALLOCATE ( igk_l2g_kdip_gamma( npw_tmp+1 ) )
    if( allocated( igk_gamma ) ) deallocate( igk_gamma )
    allocate( igk_gamma(npw_tmp+1) )
    if( allocated( g2kin_gamma ) ) deallocate( g2kin_gamma )
    allocate( g2kin_gamma(npw_tmp+1) )
    igwx_gamma=0
    !igk_l2g_kdip_gamma = 0
    g2kin_gamma=0
    iecutnew=0


    !LOOP FOR NEW ECUT
    do while (npw_gamma < npw_max )!( igwx_gamma < igwx  .and. ecutnew <= gcutm )

      call gk_sort_simple(xk_gamma(1,1),ngm,g,ecutnew,npw_gamma, igk_gamma, g2kin_gamma)
      !WRITE(stdout,*) 'Size of igk_gamma',size(igk_gamma)

      !CALL gk_l2gmap( ngm, ig_l2g(1), npw_gamma, igk_gamma, &
       !                       igk_l2g_kdip_gamma )
      !igwx_gamma = maxval(igk_l2g_kdip_gamma)
      !call mp_max(igwx_gamma, world_comm)
      iecutnew=iecutnew+1
      write(stdout,*) 'i, igwx_gamma, old ecut, new ecut, npw_gamma, npw_max'
      write(stdout,*) iecutnew, igwx_gamma,  gcutw, ecutnew, npw_gamma,npw_max
      !write(mpime+40,*) iecutnew, igwx_gamma,  gcutw, ecutnew, npw_gamma,npw_max
      !flush(mpime+40)
      ecutnew = ecutnew + 1.d0*bglen
    enddo
    !WRITE(mpime+20,*) 'proc',mpime,'exited loop'
    !WRITE(mpime+20,*) 'proc',mpime,'ecutnew',ecutnew,'npw_gamma',npw_gamma
    !WRITE(mpime+20,*) 'proc',mpime,'Size of igk_gamma',size(igk_gamma)
    !flush(mpime+20)
    ecutnew = ecutnew - 1.d0*bglen

    CALL mp_max(ecutnew,world_comm)
    
    WRITE(stdout,*) 'New ecut in gk_sort', ecutnew


    !DEALLOCATION BLOCK
    WRITE(stdout,*) 'Deallocating fft gvectors part'
    CALL deallocate_fft()

    !Need fft_type_deallocate for dfttp?
    WRITE(stdout,*) 'Deallocating fft_type'
    CALL fft_type_deallocate( dfftp )
    CALL fft_type_deallocate( dffts )

    !REALLOCATE FFT
    WRITE(stdout,*) 'Calling data_structure'
    CALL data_structure_obf(ecutnew,gamma_only) !ecutnew is taken here and gcutm is updated
    WRITE(stdout,*) 'Calling allocate_fft()'
    CALL allocate_fft()

    !CALLING GGEN FOR NEW MILL AND GVECTORS
    WRITE(stdout,*) 'gcutm before calling ggen',gcutm
    WRITE(stdout,*) 'old ecutwfc', ecutwfc
    WRITE(stdout,*) 'new ecutwfc', gcutm/4*tpiba2

    !Update ecutwfc so that the xml routine will write the correct value for obf_ham later.
    ecutwfc = gcutm/4*tpiba2


    CALL MPI_BARRIER(world_comm,ierr)
    CALL MP_MAX(gcutm,world_comm)
    !WRITE(mpime+80,*) 'ngm_g ngm gcutm',ngm_g,ngm,gcutm
    !WRITE(mpime+80,*) 'smallmem',smallmem
    !flush(mpime+80)
    CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
            g, gg, mill, ig_l2g, gstart, no_global_sort = smallmem )! Try local only with flag = .true.? Didn't work
    CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )
    WRITE(stdout,*) 'ngm after ggen', ngm
    WRITE(stdout,*) 'npw_gamma',npw_gamma
    !Determining new igwx_
    ngk_g_gamma = npw_gamma
    CALL mp_sum( ngk_g_gamma, world_comm)
    !CALL mp_sum( ngk_g_gamma, intra_bgrp_comm)
    WRITE(stdout,*) 'After ggen barrier'
    CALL MPI_BARRIER(world_comm,ierr)

    WRITE(stdout,*) 'ngk_g_gamma after mp_sum of npw_gamma',ngk_g_gamma
    !Merge mill and wfc0 across different proc
     CALL MPI_ALLREDUCE( npw_gamma, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, world_comm, IERR )  !This should be npw_max so might not need.
    !ALLOCATE (mill_tmp(3,npw_max))

    !Further expand the ecut by one increment t o find npw_nexti.
    !This should prevent out of bound ig_l2g such as case where the 2000th miller index is not found.
    !Example case: In the processor the ig_igk is at 174 but the ngw_lmax is 172.
    if( allocated( igk_gamma ) ) deallocate( igk_gamma )
    allocate( igk_gamma(ngw_lmax+1) )
    if( allocated( g2kin_gamma ) ) deallocate( g2kin_gamma )
    allocate( g2kin_gamma(ngw_lmax+1) )
    call gk_sort_count(xk_gamma(1,1),ngm,g,ecutnew+1.d0*bglen,npw_nexti)
    !WRITE(mpime+5000,*) 'npw_nexti',npw_nexti,'ngw_lmax',ngw_lmax
    !flush(mpime+5000)

    if( allocated(igk_l2g_kdip_gamma) ) deallocate( igk_l2g_kdip_gamma )
    ALLOCATE ( igk_l2g_kdip_gamma( npw_nexti ) )
    if( allocated( igk_gamma ) ) deallocate( igk_gamma )
    allocate( igk_gamma(npw_nexti) )
    if( allocated( g2kin_gamma ) ) deallocate( g2kin_gamma )
    allocate( g2kin_gamma(npw_nexti) )
    igk_gamma=0

    call gk_sort_ngw_lmax(xk_gamma(1,1),ngm,g,ecutnew,npw_tmp, igk_gamma, g2kin_gamma,npw_nexti)
      !WRITE(stdout,*) 'Size of igk_gamma',size(igk_gamma)

    CALL gk_l2gmap( ngm, ig_l2g(1), npw_nexti, igk_gamma, & !ngw_lmax instead of npw_tmp
                              igk_l2g_kdip_gamma )

    !write(mpime+50,*) ngm,ngw_lmax,npw_tmp,size(igk_gamma),size(mill_k,2),size(mill,2)
    !flush(mpime+50)
    ALLOCATE ( igk_l2g( npw_nexti ) )
      !
      ! ... the igk_l2g_kdip local-to-global map yields the correspondence
      ! ... between the global order of k+G and the local index for k+G.
      !
    !ALLOCATE ( igk_l2g_kdip( npw_nexti ) )
      !
    ALLOCATE ( mill_k( 3, npw_nexti ) )
    mill_k = 0
    igk_l2g=0
    DO ig = 1, npw_nexti !ngw_lmax+1 !ngk (ik)
            IF ( igk_gamma(ig) == 0 ) CYCLE
            igk_l2g(ig) = ig_l2g(igk_gamma(ig)) !igk_gamma replaces igk_k. This is obtained from gk_sort
    END DO

    npw_g = MAXVAL( igk_l2g(1:npw_max)) !This should be ngw?
     CALL MPI_ALLREDUCE( npw_g, itmp, 1, MPI_INTEGER, MPI_MAX, world_comm, IERR ) !Find max among npw_g put into itmp
     igwx = itmp
    WRITE(stdout,*) 'itmp ngw_lmax',itmp,ngw_lmax
    DO ig = 1,npw_nexti
     mill_k(:,ig) = mill(:,igk_gamma(ig))
    ! mill_k(:,ig) = mill(:,ig_l2g(ig))
    ENDDO

    
    ALLOCATE ( mill_collected( 3, igwx ))!npw_g ) )
    mill_collected=0
    WRITE(stdout,*) 'Size of mill', size(mill,1),size(mill,2)
    !CALL mergekg( mill_k, mill_collected, npw_max, igk_l2g_kdip_gamma, mpime, &
    !       nproc, root, world_comm )
    !CALL mergekg( mill_k, mill_collected, npw_gamma, ig_l2g, mpime, &
    !       nproc, root, world_comm )
    
    !mill_k(:,1:npw_g) = mill(:,1:npw_g)

    !do i = 1,npw_nexti
    !write(500+mpime,*) mill_k(:,i)
    !enddo
    
    !do ig =1,npw_nexti
    !write(600+mpime,*) ig,igk_gamma(ig), ig_l2g(ig), mill_k(:,ig),igk_l2g(ig)
    !flush(600+mpime)
    !enddo 

    !do ig =1,size(mill,2)
    !write(800+mpime,*) mill(:,ig)
    !enddo

    !do ig=1,npw_nexti
    ! write(900+mpime,*) ig_l2g(ig)
    ! flush(900+mpime)
    !enddo

    !if( mpime == 0 ) then
    !do i = 1,npw_g
    !write(700+mpime,*) mill_collected(:,i)
    !enddo
    !endif
    
    !Call mergek and mergewf to gather across all proc?
     
    ik=1
    !IF (lshrink ) THEN
    !       nbnd = nbnd_
    !ELSE
           nbnd = ntot_e
    !ENDIF
    !do ii=1,nbnd
    !do i=1,npw_max*npol_
    !write(1000+mpime,*) wfc0(i,ii)
    !enddo
    !enddo
    !WRITE(mpime+1000,*) size(wfc_e,1),size(wfc_e,2)
    !FLUSH(mpime+1000)
    !Should now save to something like Out_obf/system_obf.save/
    !if the output files are in Out/system.save
    WRITE(stdout,*) 'tmp_dir ',trim(tmp_dir)
    WRITE(stdout,*) 'prefix ',trim(prefix)

    source=TRIM(tmp_dir)//TRIM(prefix)//'.save/'
    WRITE(stdout,*) 'source', trim(source)
    

    length = LEN_TRIM(tmp_dir)
    IF (tmp_dir(length:length) == '/' ) length=length-1
    tmp_dir=tmp_dir(1:length)

    save_dir=trim(tmp_dir)//'_obf/'//trim(prefix)//'.save/'
    tmp_dir=trim(tmp_dir)//'_obf/'
    !prefix=trim(prefix)//'_obf.save/'
    WRITE(stdout,*) 'save_dir ',trim(save_dir)
    WRITE(stdout,*) 'tmp_dir ',tmp_dir
    WRITE(stdout,*) 'prefix ',prefix

    CALL create_directory( tmp_dir )
    CALL create_directory( save_dir )
    !Copy pseudo and charge_density.dat to new directory.
    IF(ionode) THEN
        do ii=1,nsp
            cp_source = TRIM(source)//psfile(ii)
            cp_target = TRIM(save_dir)//psfile(ii)
            !WRITE(stdout,*) 'cp_source',TRIM(cp_source)
            !WRITE(stdout,*) 'cp_target',TRIM(cp_target)
            cp_status = f_copy(cp_source,cp_target)
            WRITE(stdout,*) 'cp_status',cp_status
        enddo
            cp_source = TRIM(source)//'charge-density.dat'
            cp_target = TRIM(save_dir)//'charge-density.dat'
            !WRITE(stdout,*) 'cp_source',TRIM(cp_source)
            !WRITE(stdout,*) 'cp_target',TRIM(cp_target)
            cp_status = f_copy(cp_source,cp_target)
            WRITE(stdout,*) 'cp_status',cp_status
    ENDIF

    !Set nkstot = 1 to gamma 
    nkstot = 1
    !xk=0.0




    !Write data-file-schema.xml file
    WRITE(stdout,*) 'Calling pw_write_schema'
    CALL obf_write_schema( .false. , .false.)
     
    WRITE(stdout,*) 'Done writing xmlfile'
    WRITE(stdout,*) 'Begin reading wfc1.dat for header info'
    !Write obf
    IF (mpime == 0) THEN 
     IF (lsda) THEN
             writename=TRIM(save_dir)//'wfcup1'//'.dat' !For spin case this would be up/dw for the sake of restart routine in qe            
     ELSE
             writename=TRIM(save_dir)//'wfc1'//'.dat' !Single point calculation within QE will generate wfc1.dat
     ENDIF
     open( unit= iununk, file=trim(writename), status='unknown',form='unformatted')
     write (iununk) ik, xk_gamma, ispin, gamma_only, scalef  !For OBF. ik =1 , xk = 0.d0, ispin = as read, 
                                                      !gamma_only = F, scalef = 1
     !write (iununk) npw_g, npw_g, npol_, nbnd            !ngw == igwx_ = size(wfc_e)/MAXVAL(igwx across proc) after merging
                                                      !npol as read, nbnd should be dependent on lshrink 
     !write (iununk) igwx,igwx, npol_, nbnd 
     write (iununk) ngk_g_gamma,ngk_g_gamma, npol_, nbnd
     write (iununk) b1, b2, b3                         !As read
     !write (iununk) mill_collected(1:3,1:npw_g)
     !write (iununk) mill_k(1:3,1:npw_g)                   !New from ggen. Should build some dummy like 
     flush (iununk) 
    ENDIF                                                  !mill_k(3,igk(1->npw_max) for write
    !enddo
    WRITE(stdout,*) 'Done writing the header part of wfc_shirley_ob.dat'
    WRITE(stdout,*) 'Reached first MPI_BARRIER'
    CALL MPI_BARRIER(world_comm,ierr)
    WRITE(stdout,*) 'Passed first MPI_BARRIER'
    !
    !WRITE(stdout,*) mill(:,1:10)

    !WRITE same thing from the head of wfc.dat file

    !Merge wfc_e across different proc
     !CALL MPI_ALLREDUCE( npw_gamma, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, world_comm, IERR )  !This should be npw_max so might not need. 
    !  Use ngw_lmax = npw_max? Internally this is now npw_gamma.
     !CALL MPI_ALLREDUCE( npw_g, itmp, 1, MPI_INTEGER, MPI_MAX, world_comm, IERR ) !Find max among npw_g put into itmp
     !igwx = itmp
     
     !WRITE(2000+mpime,*) 'itmp from npw_g',itmp
     !WRITE(2000+mpime,*) 'ngw_lmax',ngw_lmax
     !WRITE(2000+mpime,*) 'npw_max',npw_max
     !WRITE(2000+mpime,*) 'npw_g', npw_g
     !WRITE(2000+mpime,*) 'Size of ig_l2g',size(ig_l2g)
     !WRITE(2000+mpime,*) 'Size of mill_k',size(mill_k)
     !WRITE(2000+mpime,*) 'Size of g',size(g,1),size(g,2)
     !WRITE(2000+mpime,*) 'mpime', mpime
     !WRITE(2000+mpime,*) 'npw_gamma',npw_gamma
     !WRITE(2000+mpime,*) 'root',root
     !WRITE(2000+mpime,*) ((mpime -1) /= root )
     !WRITE(2000+mpime,*) 'npw_nexti',npw_nexti
     !flush(2000+mpime)
     !WRITE(3000+mpime,*) ig_l2g
     !flush(3000+mpime)
     
     ALLOCATE(wfc0_collect(igwx*npol_,nbnd))
     wfc0_collect=(0.d0,0.d0)

     !Try setting npw_nexti to global max. Should resolve error with larger nproc runs.
     CALL MP_MAX(npw_nexti,world_comm)

     !For mill
     WRITE(stdout,*) 'Begin mill index stuff'
     DO ip = 1,nproc
       IF (( ip -1 ) /= root) THEN
               IF ( mpime == (ip-1) ) THEN
                       CALL MPI_SEND(ig_l2g,npw_nexti,MPI_INTEGER,root,ip+1000,world_comm,IERR)
                       CALL MPI_SEND(mill_k,3*npw_nexti,MPI_INTEGER,root,ip+2000,world_comm,IERR)
                       CALL MPI_SEND(npw_nexti,1,MPI_INTEGER,root,ip+3000,world_comm,IERR)
               ENDIF
               IF ( mpime == root) THEN
                       ALLOCATE(ig_ip(npw_nexti))
                       ALLOCATE(mill_tmp(3,npw_nexti))
                       CALL MPI_RECV(ig_ip,npw_nexti,MPI_INTEGER,(ip-1),ip+1000,world_comm,istatus,ierr)
                       !WRITE(mpime+10000,*) ip,'sent ig_ip ',npw_nexti
                       !FLUSH(mpime+10000)
                       CALL MPI_GET_COUNT(istatus,MPI_INTEGER,ngw_ip,ierr)
                       CALL MPI_RECV(mill_tmp,3*npw_nexti,MPI_INTEGER,(ip-1),ip+2000,world_comm,istatus,ierr)
                       CALL MPI_RECV(npw_max_sent,1,MPI_INTEGER,(ip-1),ip+3000,world_comm,istatus,ierr)
                       DO i = 1,npw_max_sent
                            IF ( ig_ip(i) > igwx ) THEN
                                                      !This should skip any ig_l2g that is greater than saved 
                                                               !collected wfc0
                            !WRITE(mpime+1200,*) ip,i,ig_ip(i),mill_tmp(:,i),size(ig_ip),size(mill_tmp,2) 
                            CYCLE
                            ENDIF

                            !WRITE(ip+1300,*) i,ig_ip(i),mill_tmp(:,i)
                            !Need to zero out non-existence planewave in wfc0
                            !Example when npw_max < i. We don't have wfc0(npw_max*npol,bnd) for this.
                            !IF ( i > npw_max_sent ) THEN
                            !        WRITE(mpime+1400,*) 'proc',(ip-1),'index i:',i,'exceeds',npw_max_sent
                            !        WRITE(mpime+1400,*) 'ig_ip(i)',ig_ip(i)
                            !        mill_collected(:,ig_ip(i)) = (0.d0,0.d0)
                            !        CYCLE
                            !ELSE
                            mill_collected(:,ig_ip(i)) = mill_tmp(:,i)
                            !ENDIF
                       ENDDO
                       DEALLOCATE(ig_ip)
                       DEALLOCATE(mill_tmp)
               ENDIF
       ELSE
               IF (mpime == root) THEN
                       DO ii=1,npw_nexti !This should work?
                        IF (ig_l2g(ii) > igwx) THEN
                        !WRITE(mpime+1500,*) ip,i,ig_l2g(ii),mill_k(:,ii),size(ig_l2g),size(mill_k,2)
                        !FLUSH(mpime+1500)
                          CYCLE
                        ENDIF
                          mill_collected(:,ig_l2g(ii)) = mill_k(:,ii)
                       ENDDO
               ENDIF
       ENDIF
                CALL MPI_BARRIER(world_comm,ierr)
      ENDDO
    
     WRITE(stdout,*) 'Done getting miller index across procs'

     !if( mpime == 0 ) then
     !do i = 1,npw_g
     !write(700+mpime,*) mill_collected(:,i)
     !flush(700+mpime)
     !enddo
     !endif     

     WRITE(stdout,*) 'Collecting wfc_e'


     !Should be double loop over 1->nbnd then 1 -> nproc
     DO i=1,nbnd
      DO ip = 1,nproc
        IF  (( ip -1 ) /= root) THEN
                IF ( mpime == (ip-1) ) THEN
                     CALL MPI_SEND(ig_l2g,npw_gamma, MPI_INTEGER, root, ip,world_comm,IERR)
                     CALL MPI_SEND(wfc_e(1,i),npw_gamma,MPI_DOUBLE_COMPLEX,root,ip+nproc,world_comm,IERR)
                     CALL MPI_SEND(npw_max,1,MPI_INTEGER,root,ip+nproc*4,world_comm,IERR)
                ENDIF
                IF ( mpime == root) THEN
                     ALLOCATE(ig_ip(ngw_lmax))
                     ALLOCATE(pw_ip(ngw_lmax))
                     CALL MPI_RECV(ig_ip,ngw_lmax,MPI_INTEGER,(ip-1),ip,world_comm,istatus,ierr)
                     CALL MPI_RECV(pw_ip,ngw_lmax,MPI_DOUBLE_COMPLEX,(ip-1),ip+nproc,world_comm,istatus,ierr)
                     CALL MPI_GET_COUNT(istatus,MPI_DOUBLE_COMPLEX,ngw_ip,ierr)
                     CALL MPI_RECV(npw_max_sent,1,MPI_INTEGER,(ip-1),ip+nproc*4,world_comm,istatus,ierr)
                     !WRITE(mpime+1100,*) ip,ngw_ip,size(pw_ip),size(ig_ip),npw_max_sent
                     !flush(mpime+1100)
                     DO ii =1,ngw_ip
                      IF ( ii > npw_max_sent) CYCLE
                      wfc0_collect(ig_ip(ii),i) = pw_ip(ii)
                     ENDDO
                     DEALLOCATE(ig_ip)
                     DEALLOCATE(pw_ip)
                ENDIF
        ELSE
                IF (mpime == root) THEN
                        DO ii = 1,npw_gamma
                         IF ( ii > npw_max) CYCLE
                        !if ( i==1 .and. ii ==1) THEN
                        !        write(2000,*) nbnd,npw_gamma,npw_max
                        !        write(2000,*) size(wfc_e,1),size(wfc_e,2)
                        !        write(2000,*) size(wfc0_collect,1),size(wfc0_collect,2)
                        !        flush(2000)
                        !endif
                        wfc0_collect(ig_l2g(ii),i) = wfc_e(ii,i)
                        ENDDO
                ENDIF
        ENDIF
                CALL MPI_BARRIER(world_comm,ierr)
      ENDDO

     ENDDO                

     !allocate (wfc(npol*igwx,nbnd))
    !do i = 1, nbnd
    WRITE(stdout,*) 'Size of unk in write_wfc', size(wfc_e,1),size(wfc_e,2)
    WRITE(stdout,*) 'npol,igwx,nbnd,ngk_g_gamma'
    WRITE(stdout,*) npol_,igwx,nbnd,ngk_g_gamma
    IF (ionode) THEN
     WRITE(iununk) mill_collected(1:3,1:ngk_g_gamma) !igwx is too large... get only ngk_g_gamma 
    do i =1,nbnd
     WRITE(iununk) wfc0_collect (1:npol_*ngk_g_gamma,i)                          !Here as is wfc0
    end do
    ENDIF
    close(iunwfc)
    close(iununk)
       

    WRITE(stdout,*) 'Successfully wrote wfc(ik): ', writename
    IF (lsda) THEN
            WRITE(stdout,*) 'Now writing wfc for dw spin'
            cp_source = TRIM(save_dir)//'wfcup1'//'.dat'       
            cp_target = TRIM(save_dir)//'wfcdw1'//'.dat'
            cp_status = f_copy(cp_source,cp_target)
    ENDIF
   ! STOP

    !WRITE(stdout,*) 'Size of wfc',size(evc,1),size(evc,2),size(evc)
end SUBROUTINE write_shirley_ob
