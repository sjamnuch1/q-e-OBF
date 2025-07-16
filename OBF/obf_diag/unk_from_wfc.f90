SUBROUTINE unk_from_wfc(bnk,sh,wfc,npol,nbasis,nbnd,igwx,unk)
!This subroutine will build unk from diagonalized bnk and wfc from simple.x run.
!Loop through each pw 1 -> igwx from wfc[igwx,ntot_e]
!We will then do multiplication wfc[1,ntot_e] * bnk[ntot_e,ntot_e] 
!We have unk[1,ntot_e] from each pw.

    USE kinds, ONLY : DP
    USE constants, ONLY : pi , rytoev
    USE obf_diag_objects
    USE input_obf_diag
    USE io_global, ONLY : stdout
    USE io_files,  ONLY : tmp_dir

    TYPE(shirley) :: sh
    COMPLEX(kind=DP), INTENT(IN) :: bnk(nbasis,nbnd)
    INTEGER, INTENT(IN) :: npol
    INTEGER, INTENT(IN) :: nbasis
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(IN) :: igwx
    COMPLEX(kind=DP), INTENT(IN) :: wfc(npol*igwx,nbasis)
    complex(DP), INTENT(inout) :: unk(npol*igwx,nbnd)
    complex(DP), allocatable :: wfc_pw(:),wfc_tmp(:)

    INTEGER :: ntot_e
    INTEGER :: i

    !nbnd=sh%num_bands !total band can be either original band from nscf if lshrink
                      !or new optimal basis band ntot_e
    !ntot_e=sh%ntot_e

    WRITE(stdout,*) 'Building unk from shirley ob'
    !ALLOCATE(unk(igwx,ntot_e))
    !ALLOCATE(wfc_pw(nbasis),wfc_tmp(nbasis))
    WRITE(stdout,*) 'Size of bnk in unk_from_wfc.f90 ',size(bnk,1),size(bnk,2)
    WRITE(stdout,*) 'Size of wfc in unk_from_wfc.f90 ',size(wfc,1),size(wfc,2)
    WRITE(stdout,*) 'Size of unk in unk_from_wfc.f90 ',size(unk,1),size(unk,2)
    !unk=(0.d0,0.d0)
    !wfc_pw=(0.d0,0.d0)
    !wfc_tmp=(0.d0,0.d0)
    !TO DO check when do npol == 2. Need to check datastructure of bnk then wfc0.
    ! igwx * npol -> divided into each proc for paralellization
    CALL ZGEMM('N','N',igwx*npol,nbnd,nbasis,(1.d0,0.d0), &
        & wfc,igwx*npol,bnk,nbasis,(0.d0,0.d0),unk,igwx*npol)
       
    !DO i=1,igwx*npol
    !    wfc_pw=wfc(i,1:nbasis)
    !   call ZGEMM('N','N',1,nbnd,nbasis,(1.d0,0.d0), &
    !   & wfc_pw,1,bnk,nbasis,(0.d0,0.d0),wfc_tmp,1)
    !   unk(i,1:nbnd) = wfc_tmp
    !ENDDO

    !DEALLOCATE(wfc_pw,wfc_tmp)

END SUBROUTINE unk_from_wfc
