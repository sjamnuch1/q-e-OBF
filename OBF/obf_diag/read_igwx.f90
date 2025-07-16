SUBROUTINE read_igwx(iread,input,npol,nbnd,igwx)
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
    !complex(dp), INTENT(out) :: wfc(:,:)
    !complex(DP), INTENT(inout), allocatable :: wfc(:,:)
    integer, intent(OUT) :: npol
    integer, intent(OUT) :: nbnd
    integer, intent(OUT) :: igwx


    integer  :: ibnd, iunwfc
    complex(DP), allocatable :: evc(:,:)
    integer , allocatable :: mill(:,:)
    real(dp) :: xk(3)
    integer  :: ik, ispin, ngw
    integer  :: dummy_int, narg
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
    filename = TRIM(save_dir)//'wfc1' // '.dat'
    open( unit= iunwfc, file=trim(filename), status='old',form='unformatted')

    !open(UNIT = iunwfc, FILE = trim(filename), FORM = 'unformatted', status = 'old')
    read( iunwfc) ik, xk, ispin, gamma_only, scalef
    read ( iunwfc) ngw, igwx, npol, nbnd
    read (iunwfc) b1, b2, b3

    WRITE(stdout,*) 'ik',ik
    WRITE(stdout,*) 'kpoint:',xk(:)
    WRITE(stdout,*) 'Number of planewave ngw', ngw
    WRITE(stdout,*) 'Max number of planewave igwx', igwx
    WRITE(stdout,*) 'nbnd',nbnd
    !WRITE(stdout,*) 'b1',b1(:)
    !WRITE(stdout,*) 'b2',b2(:)
    !WRITE(stdout,*) 'b3',b3(:)
    !
    ! avoid reading miller indices of G vectors below E_cut for this kpoint
    ! if needed allocate  an integer array of dims (1:3,1:igwx)
    !
    close(iunwfc)

    !WRITE(stdout,*) 'Successfully read in ',filename

    !WRITE(stdout,*) 'Size of wfc',size(evc,1),size(evc,2),size(evc)
end SUBROUTINE read_igwx
