subroutine kpoint_mesh

        USE kinds,        ONLY : DP
        USE input_parameters,        ONLY : qk, nksinterp_shirley, wk
        USE input_ham
        USE cell_base, ONLY : bg 
        use io_global, ONLY : stdout 
        !USE control_flags,        ONLY : lmanual_shirley_grid
        INTEGER :: iunkpt=999,ii,jj,kk,ll

        OPEN( UNIT = iunkpt, FILE = trim(kpoint_filename), FORM='FORMATTED', &
           STATUS = 'OLD', IOSTAT = ierr )
        IF ( ierr > 0 ) THEN
                WRITE(stdout, "('fatal error opening ',A)") TRIM(kpoint_filename)
                CALL exit(ierr)
        ENDIF
        !lmanual_shirley_grid=.true.
        read(iunkpt,*) nk(1),nk(2),nk(3) ,kshift(1) ,kshift(2), kshift(3)

        IF ( kpoint_wrap) THEN
        WRITE(stdout,*) 'kpoint with BZ from [0,1]. Will not consider any kshift/offset.'
        WRITE(stdout,*) 'Total number of kpoints:', (nk(1)+1)*(nk(2)+1)*(nk(3)+1)
        nksinterp_shirley = (nk(1)+1)*(nk(2)+1)*(nk(3)+1)
        ALLOCATE(qk(3,nksinterp_shirley),wk(nksinterp_shirley))
        wk=1.d0/nksinterp_shirley
        DO ii=0,nk(1)
        DO jj=0,nk(2)
        DO kk=0,nk(3)

        ll = (kk) + (jj)*(nk(3)+1) + (ii)*(nk(3)+1)*(nk(2)+1)+1
        qk(1,ll)=dble(ii)/nk(1) !+dble(kshift(1))/nk(1)/2
        qk(2,ll)=dble(jj)/nk(2) !+dble(kshift(2))/nk(2)/2
        qk(3,ll)=dble(kk)/nk(3) !+dble(kshift(3))/nk(3)/2
        WRITE(stdout,*) 'xk',ll, qk(:,ll),wk(ll)

        ENDDO
        ENDDO
        ENDDO

        ELSE

        WRITE(stdout,*) 'Default QE automatic kpoint'        
        WRITE(stdout,*) 'Total number of kpoints:', nk(1)*nk(2)*nk(3)
        nksinterp_shirley = nk(1)*nk(2)*nk(3)
        ALLOCATE(qk(3,nksinterp_shirley),wk(nksinterp_shirley))
        wk=1.d0/nksinterp_shirley        

        DO ii=1,nk(1)
        DO jj=1,nk(2)
        DO kk=1,nk(3)

        ll = (kk-1) + (jj-1)*nk(3) + (ii-1)*nk(3)*nk(2)+1
        qk(1,ll)=dble(ii-1)/nk(1)+dble(kshift(1))/nk(1)!Same way as OCEAN. QE do /2.
        qk(2,ll)=dble(jj-1)/nk(2)+dble(kshift(2))/nk(2)!Same way as OCEAN. QE do /2.
        qk(3,ll)=dble(kk-1)/nk(3)+dble(kshift(3))/nk(3)!Same way as OCEAN. QE do /2.
        
        WRITE(stdout,*) 'xk',ll, qk(:,ll),wk(ll)

        ENDDO
        ENDDO
        ENDDO


        ENDIF


        WRITE(stdout,*) 'bg',bg
        CALL cryst_to_cart( nksinterp_shirley, qk, bg, 1 )

        WRITE(stdout,*) 'kpoint in cart coord 2pi/alat'
        DO i=1,nksinterp_shirley
        WRITE(stdout,*) 'xk', qk(:,i),wk(i)
        ENDDO

        WRITE(stdout,*) 'Done generating k-point mesh'

end subroutine kpoint_mesh
