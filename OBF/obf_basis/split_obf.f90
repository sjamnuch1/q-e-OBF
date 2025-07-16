SUBROUTINE split_obf
!This subroutine splits the obf wfc from the ionode into
!each proc. Same way as QE would do. ig_l2g will dictate
!Which proc has which miller index <-> PW.

      USE kinds
      USE parallel_include
      USE gvect,              ONLY : g, gg, mill, gcutm,gstart,ngm_g,ngm,ig_l2g
      USE input_basis
      USE io_global,          ONLY : stdout,ionode, ionode_id
      USE wvfct,              ONLY : npwx
      USE klist,              ONLY : ngk
      use mp_world, ONLY: mpime, world_comm, nproc,root
      USE mp, ONLY : mp_bcast

      IMPLICIT NONE
      INTEGER :: i,nbnd,igwx,ip,ngw_lmax,ierr,ngw_ip,j,ngwl
      COMPLEX(DP), allocatable :: PWT(:)
      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)
      REAL(DP) :: xk_gamma(3,1)
      INTEGER :: npw_tmp,npw_proc
      INTEGER, ALLOCATABLE :: igk_obf(:),igk_l2g_kdip_gamma(:)
      REAL(DP), ALLOCATABLE :: g2kin_gamma(:)


      integer istatus(MPI_STATUS_SIZE)

      nbnd = size(obf_e,2)
      write(stdout,*) 'Size of nbasis', nbnd
      CALL mp_bcast(nbnd, ionode_id , world_comm )
      !allocate(wfc_e(npwx,nbnd))  !Allocate the wfc_e which is the wfc at each proc
                                           !To the size of original npwx and nbasis
                                           !Will remove any expanded planewave
      !igwx = MAXVAL( ig_l2g(1:MAXVAL(ngk)))
      igwx = size(obf_e,1)
      CALL mp_bcast(igwx, ionode_id , world_comm )

      !write(stdout,*) 'Size of wfc_e', size(wfc_e,1),size(wfc_e,2)
      !write(mpime+3000,*) 'Size of wfc_e', size(wfc_e,1),size(wfc_e,2)
      write(stdout,*) 'Size of ig_l2g', size(ig_l2g)
      write(stdout,*) 'Size of ngk', size(ngk)
      write(stdout,*) 'Max of ngk (ngwl)', MAXVAL(ngk) 
      write(stdout,*) 'igwx from read in OBF',igwx
      write(stdout,*) 'Size of obf_e(:,i)', size(obf_e(:,1))
      write(stdout,*) 'ngm', ngm
      ngwl = MAXVAL(ngk)

      xk_gamma=0.d0
      npw_proc = igwx/nproc+1 !Not round down in case of integer

      write(stdout,*) 'Planewave per proc',npw_proc
      CALL mp_bcast(npw_proc, ionode_id , world_comm )
      allocate(wfc_e(npw_proc,nbnd))  !Allocate the wfc_e which is the wfc at each proc
                                           !To the size of original npwx and nbasis
                                           !Will remove any expanded planewave
      allocate(igk_obf(npw_proc))
      allocate(g2kin_gamma(npw_proc))
      allocate(igk_l2g_kdip_gamma(npw_proc))
      wfc_e=(0.d0,0.d0)

      write(stdout,*) 'Size of wfc_e', size(wfc_e,1),size(wfc_e,2)




      CALL gk_sort_ngw_lmax(xk_gamma(1,1),ngm,g,ecutobf,igwx, igk_obf, g2kin_gamma,npw_proc)
      !WRITE(stdout,*) 'Size of igk_gamma',size(igk_gamma)

      CALL gk_l2gmap( ngm, ig_l2g(1), npw_proc, igk_obf, & !ngw_lmax instead of npw_tmp
                              igk_l2g_kdip_gamma )
      !write(mpime+3000,*) ngwl
      !flush(mpime+3000)

      !write(mpime+5000,*) igk_obf
      !flush(mpime+5000)

      !write(mpime+7000,*) ig_l2g
      !flush(mpime+7000)

      !write(mpime+6000,*) igk_l2g_kdip_gamma
      !flush(mpime+6000)
      CALL MPI_ALLREDUCE(ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, world_comm, IERR )
      !CALL MPI_ALLREDUCE(ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, world_comm, IERR )

      write(stdout,*) 'ngw_lmax', ngw_lmax

      !Loop through each band
      do i=1,nbnd
        IF (mpime == root .AND. igwx > size(obf_e(:,i))) THEN
                WRITE(stdout,*) 'Number of planewave in obf_e is smaller than expected'
                WRITE(stdout,*) 'Size of obf_e', size(obf_e,1)
                WRITE(stdout,*) 'Expected igwx', igwx
                STOP
        ENDIF
       do ip = 1,nproc
        IF ( (ip-1) /= root) THEN 
          IF ( mpime == (ip-1) ) THEN
              !write(mpime+4000,*) ig_l2g
              !CALL MPI_SEND(ig_l2g,ngwl,MPI_INTEGER,root,ip,world_comm,ierr)
              !CALL MPI_RECV(wfc_e(:,i),ngwl,MPI_DOUBLE_COMPLEX,root,ip+nproc,world_comm,istatus,ierr)
              CALL MPI_SEND(igk_l2g_kdip_gamma,npw_proc,MPI_INTEGER,root,ip,world_comm,ierr)
              CALL MPI_RECV(wfc_e(:,i),npw_proc,MPI_DOUBLE_COMPLEX,root,ip+nproc,world_comm,istatus,ierr)
              !WRITE(mpime+2000,*) 'Received data from root to proc', ip
          ENDIF
          IF (mpime == root) THEN
              !write(mpime+4000,*) ig_l2g
              !allocate(ig_ip(ngw_lmax))
              !allocate(pw_ip(ngw_lmax))
              allocate(ig_ip(npw_proc))
              allocate(pw_ip(npw_proc))
              !CALL MPI_RECV(ig_ip,ngw_lmax,MPI_INTEGER,(ip-1),ip,world_comm,istatus,ierr)
              CALL MPI_RECV(ig_ip,npw_proc,MPI_INTEGER,(ip-1),ip,world_comm,istatus,ierr)
              CALL MPI_GET_COUNT(istatus, MPI_INTEGER,ngw_ip,ierr)
              !WRITE(ip+4000,*) 'npw_proc',npw_proc
              !WRITE(ip+4000,*) 'ngw_ip',ngw_ip
              !FLUSH(ip+4000)
              pw_ip=0.d0
              DO j = 1,ngw_ip
                 IF ( ig_ip(j) > igwx ) THEN
                         WRITE(mpime+4000,*) 'ig_ip(j)',j,ig_ip(j)
                         CYCLE
                 ENDIF
                 pw_ip(j) = obf_e(ig_ip(j),i) 
              ENDDO

              CALL MPI_SEND(pw_ip,ngw_ip,MPI_DOUBLE_COMPLEX,(ip-1),ip+nproc,world_comm,ierr)
              deallocate(ig_ip)
              deallocate(pw_ip)
              !WRITE(mpime+2000,*) 'Done sending data from root to proc', ip
          ENDIF
        ELSE
              IF ( mpime == root) THEN
              DO j = 1,npw_proc
                !wfc_e(j,i) = obf_e(ig_l2g(j),i)
                 IF ( igk_l2g_kdip_gamma(j) > igwx ) THEN
                    WRITE(mpime+4000,*) 'ig_ip(j)',j,igk_l2g_kdip_gamma(j)
                    CYCLE
                 ENDIF
                 wfc_e(j,i) = obf_e(igk_l2g_kdip_gamma(j),i)
              ENDDO
              ENDIF
        ENDIF
           CALL MPI_BARRIER(world_comm,ierr)
           !WRITE(stdout,*) 'Done with nbasis ip ',i,ip
        enddo !Loop nproc
      enddo !Loop nbnd

      !DO i=1,npw_proc
      !DO j=1,nbnd
      !WRITE(mpime+10000,*) wfc_e(i,j)
      !flush(mpime+10000)
      !ENDDO
      !ENDDO

      !Now that obf is read in and distributed. The full OBF can be deallocated from memory.
      !obf_e should only be in the root.
      IF (ionode) THEN
      DEALLOCATE(obf_e)
      ENDIF

      write(stdout,*) 'Done distributing OBF from previous run into each proc'




END SUBROUTINE split_obf
