SUBROUTINE write_wfc(iread,input,npol,nbnd,igwx,xk,unk,bnd_unk,ispin)
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
    character (len=100)  :: filename,save_dir,writename,load_dir

    integer, intent(IN) :: iread
    TYPE(input_options_obf_diag) :: input
    integer, intent(IN) :: npol
    integer, intent(IN) :: nbnd
    integer, intent(IN) :: igwx,bnd_unk
    real(dp), intent(IN) :: xk(3)
    complex(dp), INTENT(inout) :: unk(npol*igwx,nbnd)
    integer, intent(IN) , optional :: ispin
    !complex(DP), INTENT(IN), allocatable :: unk(:,:)
    integer :: igwx_


    integer  :: ibnd, iunwfc, iununk
    complex(DP), allocatable :: evc(:,:)
    integer , allocatable :: mill(:,:)
    real(dp) :: xk_(3)
    integer  :: ik, nbnd_, npol_, ngw
    integer  :: dummy_int, narg, ispin_
    logical  :: gamma_only
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    integer :: iuni = 1111, i
    real(dp) :: scalef
    real(dp) :: b1(3), b2(3), b3(3), dummy_real


    save_dir=trim(input%unkdir)//'/'//trim(input%prefix)//'.save/'
    load_dir=trim(tmp_dir)//trim(input%prefix)//'.save/'
    iunwfc=find_free_unit()
    iununk=find_free_unit()
    !filename = TRIM(save_dir)//'wfc' // TRIM(int_to_char(iread)) // '.dat'
    filename = TRIM(load_dir)//'wfc1'//'.dat' !Don't worry about up/dw because OBF is always a single file wfc1.dat
    open( unit= iunwfc, file=trim(filename), status='old',form='unformatted')

    !open(UNIT = iunwfc, FILE = trim(filename), FORM = 'unformatted', status = 'old')
    read( iunwfc) ik, xk_, ispin_, gamma_only, scalef
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
    allocate (mill(3,igwx))
    !do i = 1,igwx
    read (iunwfc) mill(1:3,1:igwx)

    writename=TRIM(save_dir)//'wfc_ob' // TRIM(int_to_char(iread)) //'.dat'
    IF (ispin .eq. 2 ) THEN
            writename=TRIM(save_dir)//'wfc_obdw' // TRIM(int_to_char(iread)) //'.dat'
    ENDIF
    open( unit= iununk, file=trim(writename), status='unknown',form='unformatted')
    write (iununk) ik, xk, ispin, gamma_only, scalef
    write (iununk) ngw, igwx_, npol_, bnd_unk!nbnd
    write (iununk) b1, b2, b3
    write (iununk) mill(1:3,1:igwx)
    !enddo
    !
    !WRITE(stdout,*) mill(:,1:10)

    !WRITE same thing from the head of wfc.dat file




    !allocate (wfc(npol*igwx,nbnd))
    !do i = 1, nbnd
    WRITE(stdout,*) 'Size of unk in write_wfc', size(unk,1),size(unk,2)
    WRITE(stdout,*) 'npol,igwx,nbnd,bnd_unk'
    WRITE(stdout,*) npol,igwx,nbnd,bnd_unk
    do i =1,bnd_unk
     WRITE(iununk) unk (1:npol*igwx,i)
    end do
    close(iunwfc)
    close(iununk)
       

    WRITE(stdout,*) 'Successfully wrote wfc(ik): ', writename

    !WRITE(stdout,*) 'Size of wfc',size(evc,1),size(evc,2),size(evc)
end SUBROUTINE write_wfc
