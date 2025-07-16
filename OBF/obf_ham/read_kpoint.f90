subroutine read_kpoint

        USE kinds,        ONLY : DP
        USE input_parameters,        ONLY : qk, nksinterp_shirley, wk
        USE input_ham
        USE cell_base, ONLY : bg 
        use io_global, ONLY : stdout 
        !USE control_flags,        ONLY : lmanual_shirley_grid
        INTEGER :: iunkpt=999

        OPEN( UNIT = iunkpt, FILE = trim(kpoint_filename), FORM='FORMATTED', &
           STATUS = 'OLD', IOSTAT = ierr )
        IF ( ierr > 0 ) THEN
                WRITE(stdout, "('fatal error opening ',A)") TRIM(kpoint_filename)
                CALL exit(ierr)
        ENDIF
        !lmanual_shirley_grid=.true.
        read(iunkpt,*) nksinterp_shirley
        WRITE(stdout,*) 'Total number of kpoints:', nksinterp_shirley
        ALLOCATE(qk(3,nksinterp_shirley),wk(nksinterp_shirley))
        DO i=1,nksinterp_shirley
        read(iunkpt,*) qk(1,i),qk(2,i),qk(3,i),wk(i)
        ENDDO

        WRITE(stdout,*) 'kpoint in crystal coord'
        DO i=1,nksinterp_shirley
        WRITE(stdout,*) 'xk', qk(:,i),wk(i)
        ENDDO

        WRITE(stdout,*) 'bg',bg
        CALL cryst_to_cart( nksinterp_shirley, qk, bg, 1 )

        WRITE(stdout,*) 'kpoint in cart coord 2pi/alat'
        DO i=1,nksinterp_shirley
        WRITE(stdout,*) 'xk', qk(:,i),wk(i)
        ENDDO

end subroutine read_kpoint
