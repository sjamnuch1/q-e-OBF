SUBROUTINE read_wfc(iread,input,npol,nbnd,igwx,wfc)
!This subroutine will read wfc_(iread).dat from simple.x run.
!Then return wfc and igwx from that file.
!Rest of variables stored will not be collected.
!Will reread the wfc when doing write wfc.
        
    USE kinds, ONLY : DP
    USE constants, ONLY : pi , rytoev
    USE obf_diag_objects
    USE input_obf_diag
    USE io_global, ONLY : stdout
    USE io_files,  ONLY : tmp_dir 

    implicit none
    character (len=100)  :: filename,save_dir

    integer, intent(IN) :: iread
    TYPE(input_options_obf_diag) :: input
    integer, intent(IN) :: npol
    integer, intent(IN) :: nbnd
    integer, intent(IN) :: igwx
    complex(dp), INTENT(inout) :: wfc(npol*igwx,nbnd)
    !complex(DP), INTENT(inout), allocatable :: wfc(:,:)
    !integer, intent(OUT) :: igwx


    integer  :: ibnd, iunwfc
    complex(DP), allocatable :: evc(:,:)
    integer , allocatable :: mill(:,:)
    real(dp) :: xk(3)
    integer  :: ik, ispin, ngw
    integer  :: dummy_int, narg
    integer  :: dummy_int2
    integer  :: dummy_int3
    integer  :: dummy_int4
    integer  :: dummy_int5
    !integer  :: stdout =6
    logical  :: gamma_only
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    integer ::  i
    real(dp) :: scalef
    real(dp) :: b1(3), b2(3), b3(3), dummy_real

    save_dir=trim(tmp_dir)//trim(input%prefix)//'.save/'
    iunwfc=find_free_unit()
    filename = TRIM(save_dir)//'wfc1'//'.dat'
    !filename = TRIM(save_dir)//'wfc' // TRIM(int_to_char(iread)) // '.dat'
    open( unit= iunwfc, file=trim(filename), status='old',form='unformatted')

    !open(UNIT = iunwfc, FILE = trim(filename), FORM = 'unformatted', status = 'old')
    read( iunwfc) ik, xk, ispin, gamma_only, scalef
    read ( iunwfc) ngw, dummy_int, dummy_int2, dummy_int3
    read (iunwfc) b1, b2, b3

    !WRITE(stdout,*) 'ik',ik
    !WRITE(stdout,*) 'kpoint:',xk(:)
    !WRITE(stdout,*) 'Number of planewave ngw', ngw
    !WRITE(stdout,*) 'Max number of planewave igwx', igwx
    !WRITE(stdout,*) 'nbnd',nbnd
    WRITE(stdout,*) 'b1',b1(:)
    WRITE(stdout,*) 'b2',b2(:)
    WRITE(stdout,*) 'b3',b3(:)
    !
    ! avoid reading miller indices of G vectors below E_cut for this kpoint
    ! if needed allocate  an integer array of dims (1:3,1:igwx)
    !
    allocate (mill(3,igwx))
    !do i = 1,igwx
    read (iunwfc) mill(1:3,1:igwx)
    !enddo
    !
    !WRITE(stdout,*) mill(:,1:10)
    WRITE(stdout,*) 'Size of input wfc',size(wfc,1),size(wfc,2),size(wfc)

    !allocate (wfc(npol*igwx,nbnd))
    do i = 1, nbnd
          read(iunwfc) wfc (1:npol*igwx,i)
    end do
    close(iunwfc)
    !wfc = evc
    !deallocate(evc)
    

    WRITE(stdout,*) 'Successfully read in ',filename

    WRITE(stdout,*) 'Size of output wfc',size(wfc,1),size(wfc,2),size(wfc)
end SUBROUTINE read_wfc
