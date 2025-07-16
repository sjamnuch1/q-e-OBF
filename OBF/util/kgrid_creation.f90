PROGRAM kgrid_creation
INTEGER :: nk(3),mode
REAL :: xk(3)
character(len=12), dimension(:), allocatable :: args
INTEGER :: narg,i,j,k,iglobal,nks
REAL , allocatable :: q(:,:)


narg = command_argument_count()
allocate(args(narg))

if (narg /= 3 ) then
WRITE(*,*) "Usage kgrid_creation.x nk1 nk2 nk3"
stop
endif

do i=1,3
call get_command_argument(i,args(i))
!WRITE(*,*) i,args(i)
read(args(i),'(I10)') nk(i)
!WRITE(*,*) i,nk(i)
enddo

call get_command_argument(4,args(4))
read(args(4),'(I10)') mode

!if (mode == 0) THEN

!DO i=0,nk(1)
!DO j=0,nk(2)
!DO k=0,nk(3)
!   iglobal = (k) + (j)*(1+nk(3)) + (i)*(1+nk(2))*(1+nk(3)) + 1
!  xk(1) = dble(i)/nk(1)
!   xk(2) = dble(j)/nk(2)
!   xk(3) = dble(k)/nk(3)
!  WRITE(*,*) iglobal,xk

!ENDDO
!ENDDO
!ENDDO

nks=nk(1)*nk(2)*nk(3) ! 7 other gamma points
ALLOCATE(q(nks+7,3))
q=0.d0

!ELSEIF (mode == 1 ) THEN
        DO i=1,nk(1)
        DO j=1,nk(2)
        DO k=1,nk(3)
           iglobal = (k-1) + (j-1)*nk(3) + (i-1)*nk(2)*nk(3) + 1
              xk(1) = dble(i-1)/nk(1)
              q(iglobal,1)=xk(1)
              xk(2) = dble(j-1)/nk(2)
              q(iglobal,2)=xk(2)
              xk(3) = dble(k-1)/nk(3)
              q(iglobal,3)=xk(3)
              !WRITE(*,*) iglobal,xk

        ENDDO
        ENDDO
        ENDDO
!ENDIF

DO i=nks+1,nks+3
j=i-nks
q(i,j)=1.d0
ENDDO
q(nks+4:nks+7,:)=1.d0
DO i=nks+4,nks+6
j=i-nks-3
q(i,j)=0.d0
ENDDO

DO i=1,nks+7
WRITE(*,*) i,q(i,:)
ENDDO

WRITE(*,*) 'Done'
DEALLOCATE(q,args)
STOP

END PROGRAM kgrid_creation
