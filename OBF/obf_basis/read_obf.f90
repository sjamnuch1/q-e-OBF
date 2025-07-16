SUBROUTINE read_obf()
!This subroutine will read the wfc1.dat which is the obf from obf_basis.x
!the wfc is then stored and be used for building obf. This is for using multiple wfc.dat 
!with different bands size/mesh.
USE kinds, ONLY : DP
USE input_basis
USE io_files,  ONLY : prefix
USE io_global, ONLY : stdout
USE wvfct, ONLY : npw,npwx
USE gvecw,                ONLY : gcutw

    implicit none
    character (len=100)  :: filename,read_dir
    integer :: igwx,npol,nbnd
    complex(dp),allocatable :: wfc(:,:)
    integer , allocatable :: mill(:,:)
    real(dp) :: xk(3)
    real(dp) :: scalef
    real(dp) :: b1(3), b2(3), b3(3), dummy_real
    integer  :: ik, ispin, ngw,iunwfc
    integer  :: dummy_int, narg,i
    integer  :: dummy_int2
    integer  :: dummy_int3
    integer  :: dummy_int4
    integer  :: dummy_int5
    !integer  :: stdout =6
    logical  :: gamma_only
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER(LEN=6), EXTERNAL :: int_to_char

    !IF ( npw_obf <= 0) THEN
    !        WRITE(stdout,*) 'npw_obf is not specified in the input or less than zero. Exiting.'
    !        STOP
    !ENDIF

    IF ( ecutobf <= gcutw ) THEN
            WRITE(stdout,*) 'ecutobf is not specified in the input or less than dft gcutw. Exiting'
            STOP
    ENDIF

    read_dir=trim(obf_dir)//trim(prefix)//'.save/'
    filename=trim(read_dir)//'wfc1'//'.dat' 
    iunwfc=find_free_unit()
    open( unit= iunwfc, file=trim(filename), status='old',form='unformatted')
    read( iunwfc) ik, xk, ispin, gamma_only, scalef
    read ( iunwfc) ngw, igwx, npol, nbnd
    read (iunwfc) b1, b2, b3

    allocate(mill(3,igwx))
    read (iunwfc) mill(1:3,1:igwx)
    deallocate(mill)

    allocate(wfc(npol*igwx,nbnd))
    do i = 1, nbnd
            read(iunwfc) wfc (1:npol*igwx,i)
    end do
    close(iunwfc)
    allocate(obf_e(npol*igwx,nbnd))
    obf_e = wfc
    deallocate(wfc)

    WRITE(stdout,*) 'Successfully read in obf wfc'
    WRITE(stdout,*) 'Size of OBF wfc', size(obf_e,1),size(obf_e,2)
    WRITE(stdout,*) 'npw npwx', npw,npwx
    !WRITE(stdout,*) 'wfc(1,1) ', wfc_e(1,1)


END SUBROUTINE read_obf
