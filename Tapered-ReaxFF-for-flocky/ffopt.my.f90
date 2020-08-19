!*********************************************************************
!                                                                    *
!     REAXFF Reactive force field program                            *
!                                                                    *
!     Developed and written by Adri van Duin, duin@wag.caltech.edu   *
!                                                                    *
!     Copyright (c) 2001-2010 California Institute of Technology     *
!                                                                    *
!     This is an open-source program. Feel free to modify its        *
!     contents. Please keep me informed of any useful modification   *
!     or addition that you made. Please do not distribute this       *
!     program to others; if people are interested in obtaining       *
!     a copy of this program let them contact me first.              *
!                                                                    *
!                          * * * *                                   *
!                                                                    *
!     Modified by David Furman, df398@cam.ac.uk                      *
!     University of Cambridge                                        *
!                                                                    *
!     New features:                                                  *
!     - Tapered ReaxFF model (2019)                                  *
!     - Low gradient (lg) correction of Liu et al. JCPA 2011 (2019)  *
!     - Numerically stable dihedrals formulation (2019)              *
!     - Added atomic forces to training set (2020)                   *
!     - Numerically stable lone pairs formulation (2020)             *
!     - Numerically stable sbo2 formulation in valence angles (2020) *
!                                                                    *
!*********************************************************************
!*********************************************************************

    program ffopt

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(33) :: qfile2
    character(33) :: qfilh
!*********************************************************************
!                                                                    *
!     Optimization routines for the reactive MD-force field          *
!                                                                    *
!*********************************************************************
    ustime=zero
    systime=zero
    big=1e+30
    starttime=secnds(0.0)
    vpmax=zero
    vpmin=zero
    read (35,*)dseed
    iagain=0
    do i1=1,nat
        do i2=1,mbond+3
            ia(i1,i2)=0
            iag(i1,i2)=0
        end do
    end do
    read (20,1000)iopt

    if (iopt == 0) then
        call version
        call reac
        call calcerr(1)
        if (nsurp >= 2) then
            rewind (90)
            rewind (98)
            do i1=1,nmollset
                call extractgeo(i1)
                call writebgf(90)
                call writegeo(98)
                 
                if (nsurp == 3) then
                    qfilh=qfile(i1)
                    qfile2=qfile(i1)
                    if (imodfile == 0) then
                        istart=1
                        qstrana1(1:33)=qfilh
                        call stranal(istart,iend,vout,iout,1)
                        qfile2=qfilh(istart:iend-1)//".bgf"
                    end if
                    open (88,file=qfile2,status='unknown')
                    if (imodfile == 0) call writebgf(88)
                    if (imodfile == 1 .AND. ngeofor == 1) call writebgf(88)
                    close (88)
                    if (imodfile == 0) qfile2=qfilh(istart:iend-1)//".geo"
                    open (88,file=qfile2,status='unknown')
                    if (imodfile == 0) call writegeo(88)
                    if (imodfile == 1 .AND. ngeofor == 0) call writegeo(88)
                    close (88)
                end if
            end do
        end if
         
        write (*,*) 'Normal end of MD-simulation'
        stop
    end if

    if (iopt == 2) then !force field extrapolation
        call ffext
        stop 'Normal end of program; extrapolated force field in unit 13'
    end if

    call version
    call readc
    call ffin2(4)
!     read (21,1100)ichn(1),ichn(2),ichn(3),vchange,vpmax,vpmin
!     read (21,*)ichn(1),ichn(2),ichn(3),vchange,vpmax,vpmin
    read (21,'(a200)')qstrana1
    istart=1
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    ichn(1)=iout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    ichn(2)=iout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    ichn(3)=iout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    vchange=vout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    vpmax=vout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    vpmin=vout
    if (vpmax < vpmin) then
        vph=vpmax
        vpmax=vpmin
        vpmin=vph
    end if
    if (vchange < 0.00001) vchange=0.00001
    if (vpmax < 0.001 .AND. vpmax > -0.001) vpmax=25000.0
    if (vpmin < 0.001 .AND. vpmin > -0.001) vpmin=-25000.0
    icount=0
    valpar=fpar(ichn(1),ichn(2),ichn(3))
    call optdat
    9 change=parc1*vchange
    open (25,file='fort.25',status='unknown')
    write (25,1200)ichn(1),ichn(2),ichn(3),valpar,vchange*parc1, &
    vpmax,vpmin
    close (25)
    call ffchng(1,4,big,-big)

    iflga=0
    do i1=1,3
        call ffchng(i1,4,big,-big)
        call reac
        call increm
        call rewind
    !     call bepsdy(i1)
        call calcerr(i1)
    end do
    errsav=sdy(3)
    call kmin(vpmax,vpmin)
    call ffchng(4,4,big,-big)
    iflga=1
    call reac
    call increm
    call rewind
    call rewind2
!     call bepsdy(3)
    call calcerr(3)
    if (sdy(3) > errsav+accincr .AND. icount < 4) then     !Do not accept parameter change
        call ffchng(3,4,big,-big)
        if (sdy(1) < errsav) call ffchng(1,4,vpmax,vpmin)
        if (sdy(2) < errsav) call ffchng(2,4,vpmax,vpmin)
        if (sdy(1) > errsav .AND. sdy(2) > errsav) then
            vchange=vchange*0.50
            icount=icount+1
            goto 9
        end if
    end if
    if (ingeo == 1 .AND. sdy(3) < errsav .AND. nsurp /= 2 &
    .and. ngeofor == 0) call cpfile(98,3)
    if (ingeo == 1 .AND. sdy(3) < errsav .AND. nsurp /= 2 &
    .and. ngeofor == 1) call cpfile(90,3)
    if (nsurp >= 2) then
        rewind (90)
        rewind (98)
        do i1=1,nmollset
            call extractgeo(i1)
            call writebgf(90)
            call writegeo(98)
            if (imodfile == 1) then
                open (88,file=qfile(i1),status='unknown')
                if (ngeofor == 0) call writegeo(88)
                if (ngeofor == 1) call writebgf(88)
                close(88)
            end if
        end do
        if (imodfile == 0 .AND. sdy(3) < errsav .AND. ingeo == 1 &
        .and. ngeofor == 0) call cpfile(98,3)
        if (imodfile == 0 .AND. sdy(3) < errsav .AND. ingeo == 1 &
        .and. ngeofor == 1) call cpfile(90,3)
    end if

    call ffsav
    10 continue

    read (21,'(a200)',end=20,err=20)qstrana1
    istart=1
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    ichn(1)=iout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    ichn(2)=iout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    ichn(3)=iout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    vchange=vout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    vpmax=vout
    call stranal(istart,iend,vout,iout,1)
    istart=iend
    vpmin=vout
    if (vpmax < vpmin) then
        vph=vpmax
        vpmax=vpmin
        vpmin=vph
    end if
    if (vpmax < 0.001 .AND. vpmax > -0.001) vpmax=25000.0
    if (vpmin < 0.001 .AND. vpmin > -0.001) vpmin=-25000.0
    if (vchange < 0.00001) vchange=0.00001
    iagain=0
    icount=0
    valpar=fpar(ichn(1),ichn(2),ichn(3))
    11 change=vchange*parc1
    open (25,file='fort.25',status='unknown')
    write (25,1200)ichn(1),ichn(2),ichn(3),valpar,vchange*parc1, &
    vpmax,vpmin
    close (25)
    iflga=0
    do i1=1,2+ingeo
    !     do i1=1,2
        if (i1 == 3) iflga=1
        call ffchng(i1,4,big,-big)
        call reac
        call increm
        call rewind
    !     call bepsdy(i1)
        call calcerr(i1)
    end do
    errsav=sdy(3)
    call kmin(vpmax,vpmin)
    call ffchng(4,4,big,-big)
    call rewind2
    call reac
    call rewind
    call increm
!     call bepsdy(3)
    call calcerr(3)
    if (sdy(3) > errsav+accincr .AND. icount < 4) then     !Do not accept parameter change
        call ffchng(3,4,big,-big)
        if (sdy(1) < errsav) call ffchng(1,4,vpmax,vpmin)
        if (sdy(2) < errsav) call ffchng(2,4,vpmax,vpmin)
        if (sdy(1) > errsav .AND. sdy(2) > errsav) then
            vchange=vchange*0.50
            icount=icount+1
            goto 11
        end if
    end if
    if (ingeo == 1 .AND. sdy(3) < errsav .AND. nsurp /= 2 .and. &
    ngeofor == 0) call cpfile(98,3)
    if (ingeo == 1 .AND. sdy(3) < errsav .AND. nsurp /= 2 .and. &
    ngeofor == 1) call cpfile(90,3)
    if (nsurp >= 2) then
        rewind (90)
        rewind (98)
        do i1=1,nmollset
            call extractgeo(i1)
            call writebgf(90)
            call writegeo(98)
            if (imodfile == 1) then
                open (88,file=qfile(i1),status='unknown')
                if (ngeofor == 0) call writegeo(88)
                if (ngeofor == 1) call writebgf(88)
                close(88)
            end if
        end do
        if (imodfile == 0 .AND. sdy(3) < errsav .AND. ingeo == 1 .and. &
        ngeofor == 0) call cpfile(98,3)
        if (imodfile == 0 .AND. sdy(3) < errsav .AND. ingeo == 1 .and. &
        ngeofor == 1) call cpfile(90,3)
    end if
    call ffsav
    goto 10
    20 continue
    call ffchng(4,4,big,-big)
    call finrep
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    1000 format (i3)
    1100 format (3i3,3f8.4)
    1200 format (3i3,4f12.4)
    END PROGRAM
!*********************************************************************
!*********************************************************************

    subroutine ffin2(nunit)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!                                                                    *
!     Read in force field for optimization                           *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In ffin2'
        call timer(65)
        close (65)
    end if
    read (nunit,1000)qff(1)
    read (nunit,1010)npar,qff(2)
    do i1=1,npar
        read (nunit,1015)fpar(1,i1,1),qff(2+i1)
    end do

    ! df398 read lg dispersion scaling factor. (when iopt=1 all ffield is read here, not in reac.f90)
    dispscale = fpar(1,36,1) 

    read (nunit,1010)nso,qff(2+npar+1)
    read (nunit,1000)qff(2+npar+2)
    read (nunit,1000)qff(2+npar+3)
    read (nunit,1000)qff(2+npar+4)
    do i1=1,nso
        read (nunit,1200)qas2(i1),(fpar(2,i1,i2),i2=1,8)
        read (nunit,1250)(fpar(2,i1,i2),i2=9,16)
        read (nunit,1250)(fpar(2,i1,i2),i2=17,24)
        read (nunit,1250)(fpar(2,i1,i2),i2=25,32)
        if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 lg included
           read (nunit,1260)(fpar(2,i1,i2),i2=33,34)
        endif
    end do
    ihulp=2+npar+4
    read (nunit,1010)nboty,qff(ihulp+1)
    read (nunit,1000)qff(ihulp+2)
    do i1=1,nboty
        read (nunit,1400)iboo(i1,1),iboo(i1,2),(fpar(3,i1,i2),i2=1,8)
        read (nunit,1450)(fpar(3,i1,i2),i2=9,16)
    end do
    read (nunit,1010)nodmty,qff(ihulp+3)
    do i1=1,nodmty
        if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 lg included
            read (nunit,1400)idmo(i1,1),idmo(i1,2),(fpar(4,i1,i2),i2=1,7)
        else
            read (nunit,1400)idmo(i1,1),idmo(i1,2),(fpar(4,i1,i2),i2=1,6)
        endif
    end do
    read (nunit,1010)nvaty,qff(ihulp+4)
    do i1=1,nvaty
        read (nunit,1500)ivao(i1,1),ivao(i1,2),ivao(i1,3), &
        (fpar(5,i1,i2),i2=1,7)
    end do
    read (nunit,1010)ntoty,qff(ihulp+5)
    do i1=1,ntoty
        read (nunit,1600)itoo(i1,1),itoo(i1,2),itoo(i1,3),itoo(i1,4), &
        (fpar(6,i1,i2),i2=1,7)
    end do
    read (nunit,1010)nhbty,qff(ihulp+6)
    do i1=1,nhbty
        read (nunit,1500)ihbo(i1,1),ihbo(i1,2),ihbo(i1,3), &
        (fpar(7,i1,i2),i2=1,4)
    end do
    return
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    1000 format (a80)
    1010 format (i3,a80)
    1015 format (f10.4,a80)
    1100 format (i3,2x,a2,3x,3d22.15)
    1200 format (1x,a2,10f9.4)
    1250 format (3x,10f9.4)
    1260 format (3x,2f9.4)
    1300 format (f10.4)
    1400 format (2i3,8f9.4)
    1450 format (6x,8f9.4)
    1500 format (3i3,7f9.4)
    1600 format (4i3,7f9.4)
    end subroutine ffin2
!*********************************************************************
!*********************************************************************

    subroutine ffchng(ichange,nunit,vhigh,vlow)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!                                                                    *
!     Change the force field                                         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In ffchng'
        call timer(65)
        close (65)
    end if
    if (ichange > 4) goto 30   !just write the force field
    if (valpar >= zero) then
        if (ichange == 1) vparn=valpar-change
        if (ichange == 2) vparn=valpar+change
        if (ichange == 3) vparn=valpar
        if (ichange == 4) vparn=valnew
    else
        if (ichange == 1) vparn=valpar+change
        if (ichange == 2) vparn=valpar-change
        if (ichange == 3) vparn=valpar
        if (ichange == 4) vparn=valnew
    end if

    if (vparn > vhigh) vparn=vhigh
    if (vparn < vlow) vparn=vlow
    fpar(ichn(1),ichn(2),ichn(3))=vparn
!*********************************************************************
!                                                                    *
!     Change linked parameters                                       *
!                                                                    *
!*********************************************************************
    10 read (23,'(4i3)',end=20)ihk1,ihk2,ihk3,ihk4
    if (ihk1 == ichn(1) .AND. ihk2 == ichn(2) .AND. ihk3 == ichn(3)) &
    then
        do i2=1,ihk4
            read (23,'(3i3)')ipk1,ipk2,ipk3
            fpar(ipk1,ipk2,ipk3)=vparn
        end do
    else
        do i2=1,ihk4
            read (23,'(3i3)')ipk1,ipk2,ipk3
        end do
    end if
    goto 10
    20 rewind (23)
!*********************************************************************
!                                                                    *
!     Write the force field                                          *
!                                                                    *
!*********************************************************************
    30 continue
    rewind (nunit)
    write (nunit,1000)qff(1)
    write (nunit,1010)npar,qff(2)
    do i1=1,npar
        write (nunit,1015)fpar(1,i1,1),qff(2+i1)
    end do
    write (nunit,1010)nso,qff(2+npar+1)
    write (nunit,1000)qff(2+npar+2)
    write (nunit,1000)qff(2+npar+3)
    write (nunit,1000)qff(2+npar+4)
    do i1=1,nso
        write (nunit,1200)qas2(i1),(fpar(2,i1,i2),i2=1,8)
        write (nunit,1250)(fpar(2,i1,i2),i2=9,16)
        write (nunit,1250)(fpar(2,i1,i2),i2=17,24)
        write (nunit,1250)(fpar(2,i1,i2),i2=25,32)
        if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 lg included
           write (nunit,1250)(fpar(2,i1,i2),i2=33,34)
        endif
    end do
    if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 lg included
        ihulp=2+npar+5
    else
        ihulp=2+npar+4
    endif
    write (nunit,1010)nboty,qff(ihulp+1)
    write (nunit,1000)qff(ihulp+2)
    do i1=1,nboty
        write (nunit,1400)iboo(i1,1),iboo(i1,2),(fpar(3,i1,i2),i2=1,8)
        write (nunit,1450)(fpar(3,i1,i2),i2=9,16)
    end do
    write (nunit,1010)nodmty,qff(ihulp+3)
    do i1=1,nodmty
        if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 lg included
            write (nunit,1400)idmo(i1,1),idmo(i1,2),(fpar(4,i1,i2),i2=1,7)
        else
           write (nunit,1400)idmo(i1,1),idmo(i1,2),(fpar(4,i1,i2),i2=1,6)
        endif
    end do
    write (nunit,1010)nvaty,qff(ihulp+4)
    do i1=1,nvaty
        write (nunit,1500)ivao(i1,1),ivao(i1,2),ivao(i1,3), &
        (fpar(5,i1,i2),i2=1,7)
    end do
    write (nunit,1010)ntoty,qff(ihulp+5)
    do i1=1,ntoty
        write (nunit,1600)itoo(i1,1),itoo(i1,2),itoo(i1,3),itoo(i1,4), &
        (fpar(6,i1,i2),i2=1,7)
    end do
    write (nunit,1010)nhbty,qff(ihulp+6)
    do i1=1,nhbty
        write (nunit,1500)ihbo(i1,1),ihbo(i1,2),ihbo(i1,3), &
        (fpar(7,i1,i2),i2=1,4)
    end do
    rewind (nunit)
    return
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    1000 format (a80)
    1010 format (i3,a65)
    1015 format (f10.4,a40)
    1100 format (i3,2x,a2,3x,3d22.15)
    1200 format (1x,a2,10f9.4)
    1250 format (3x,10f9.4)
    1300 format (f10.4)
    1400 format (2i3,8f9.4)
    1450 format (6x,8f9.4)
    1500 format (3i3,7f9.4)
    1600 format (4i3,7f9.4)
    end subroutine ffchng
!*********************************************************************
!*********************************************************************

    subroutine optdat

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!                                                                    *
!     Read force field optimization data                             *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In optdat'
        call timer(65)
        close (65)
    end if
    imam=1
    ndata=0
    iheatf=0
    10 read (9,100,end=20,err=20)ndatm(imam),qmdat(imam)
    do i1=1,ndatm(imam)
        read (9,200)idat(i1+ndata),datopt(i1+ndata),devi(i1+ndata)
        if (idat(i1+ndata) == 4) then
            iheatf=iheatf+1
            iheada(iheatf)=i1+ndata
        end if
    end do
    ndata=ndata+ndatm(imam)
    imam=imam+1
    goto 10
    20 continue
    rewind (9)
!*********************************************************************
!                                                                    *
!     Read in coupled data                                           *
!                                                                    *
!*********************************************************************
    read (22,300,err=30,end=30)nkop
    do i1=1,nkop
        read (22,300)ikop1(i1),ikop2(i1),mu1(i1),mu2(i1), &
        vkop(i1),devkop(i1)
    end do
    30 continue
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format (i3,a60)
    200 format (i3,12x,f10.2,f10.4)
    300 format (4i4,2f8.4)
    return
    end subroutine optdat
!*********************************************************************
!*********************************************************************

    subroutine bepsdy(inum)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!                                                                    *
!     Calculate deviation between training set and calculated data   *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In bepsdy'
        call timer(65)
        close (65)
    end if
    toterr=zero
    ihu=0
    do i1=1,imam
        write (99,100)qmdat(i1)
        do i2=1,ndatm(i1)
            ihu=ihu+1
            diff=datopt(ihu)-caldat(ihu)
            err=(diff/devi(ihu))**2
            toterr=toterr+err
            write (99,200)ihu,datopt(ihu),caldat(ihu),devi(ihu),err,toterr
        end do
    end do
    do i1=1,nkop
        diffh=caldat(ikop1(i1))/mu1(i1)-caldat(ikop2(i1))/mu2(i1)
        diff=diffh-vkop(i1)
        err=(diff/devkop(i1))**2
        toterr=toterr+err
        write (99,300)caldat(ikop1(i1))/mu1(i1), &
        caldat(ikop2(i1))/mu2(i1),diffh,vkop(i1),devkop(i1),err,toterr
    end do
    sdy(inum)=toterr
    open (13,file='fort.13',status='unknown',access='append')
    write (13,400)toterr
    close (13)
    rewind (99)
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format (a60)
    200 format (i4,5f20.6)
    300 format (7f20.6)
    400 format (f20.6)
    return
    end subroutine bepsdy
!*********************************************************************
!*********************************************************************

    subroutine calcerr(inum)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(80) :: qrom2
!*********************************************************************
!                                                                    *
!     Calculate error force field and generate output to unit 99     *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In calcerr'
        call timer(65)
        close (65)
    end if
    toterr=zero
    open (99,file='fort.99',status='unknown')
    write (99,100)
    if (nsurp >= 2) nprob=nmollset+1
    do i1=1,nprob-1
        ifound=0
        do i2=1,ndata2
            qrom2=qdatid(i2)
            if (qrom2(1:20) == qkeyw(i1)) then
                error=(caldat(i2)-compdat(i2))/weightdat(i2)
                toterr=toterr+error*error
                if (ifound == 0) then
                    write (99,'(a60,5f20.4)')qrom2(1:80),caldat(i2), &
                    compdat(i2),weightdat(i2),error*error,toterr
                else
                    write (99,'(a60,5f20.4)')qrom2(22:80),caldat(i2), &
                    compdat(i2),weightdat(i2),error*error,toterr
                end if

            !     if (ifound.eq.0) write (99,'(a10)')qkeyw(i1)
            !     write (99,'(a60,5f12.4)')qrom2(12:60),caldat(i2),
            !    $compdat(i2),weightdat(i2),error*error,toterr

                ifound=1
            end if
        end do
    end do
    write (99,*)
    do i2=1,ndata2
        qrom2=qdatid(i2)
        if (qrom2(1:6) == 'Energy') then
            error=(caldat(i2)-compdat(i2))/weightdat(i2)
            toterr=toterr+error*error
            write (99,'(a60,5f20.4)')qdatid(i2),caldat(i2), &
            compdat(i2),weightdat(i2),error*error,toterr
        end if
        if (qrom2(1:4) == 'Geo ') then
            error=(caldat(i2)-compdat(i2))/weightdat(i2)
            toterr=toterr+error*error
            write (99,'(a60,5f20.4)')qdatid(i2),caldat(i2), &
            compdat(i2),weightdat(i2),error*error,toterr
        end if
    end do
    sdy(inum)=toterr
    open (13,file='fort.13',status='unknown',access='append')
    write (13,400)toterr
    close (13)
    close(99)
    !rewind (99)

    100 format (71x,'FField value',7x,'QM/Lit value',11x,'Weight', &
    & 14x,'Error',13x,'Total error')
    400 format (f30.10)
    return
    end subroutine calcerr
!*********************************************************************
!*********************************************************************

    subroutine finrep

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension qlab(5)
    character(20) :: qlab
!*********************************************************************
!                                                                    *
!     Write out final report of minimization                         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In finrep'
        call timer(65)
        close (65)
    end if
    qlab(1)='Bond length      '
    qlab(2)='Valency angle    '
    qlab(3)='Torsion angle    '
    qlab(4)='Heat of formation'
    toterr=zero
    ihu=0
    do i1=1,imam-1
        write (67,100)qmdat(i1)
        write (67,150)
        write (67,*)
        do i2=1,ndatm(i1)
            ihu=ihu+1
            diff=datopt(ihu)-caldat(ihu)
            err=(diff/devi(ihu))**2
            toterr=toterr+err
            write (67,200)qlab(idat(ihu)),datopt(ihu), &
            caldat(ihu)
        end do
        write (67,*)
    end do
    do i1=1,nkop
        diffh=caldat(ikop1(i1))-caldat(ikop2(i1))
        diff=diffh-vkop(i1)
        err=(diff/devkop(i1))**2
        toterr=toterr+err
        write (67,300)caldat(ikop1(i1)),caldat(ikop2(i1)), &
        diffh,vkop(i1)
    end do
    ihu=0
    do i1=1,imam-1
        do i2=1,ndatm(i1)
            ihu=ihu+1
            if (idat(ihu) == 4) then
                write (66,310)qmdat(i1),datopt(ihu),caldat(ihu)
            end if
        end do
    end do

          
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format (a60)
    150 format ('Data type',15x,'Literature',2x,'Calculated')
    200 format (a20,2f12.4)
    300 format (20x,4f12.4)
    310 format (a40,4f12.4)
    return
    end subroutine finrep
!*********************************************************************
!*********************************************************************

    subroutine rewind

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!                                                                    *
!     Rewind files                                                   *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In rewind'
        call timer(65)
        close (65)
    end if
    rewind (2)
    rewind (3)
    rewind (9)
    rewind (10)
    rewind (18)
    rewind (28)
    rewind (38)
    rewind (55)
    rewind (56)
    rewind (58)
    rewind (67)
    rewind (76)
    rewind (86)
    rewind (90)
    rewind (98)
    rewind (81)
    rewind (91)
!     open (unit=68,file='xmolout',status='unknown')

    return
    end subroutine rewind
!*********************************************************************
!*********************************************************************

    subroutine rewind2

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!                                                                    *
!     Rewind other files                                             *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In rewind2'
        call timer(65)
        close (65)
    end if
    open (unit=71,file='fort.71',status='unknown')
    write (71,*)
    close (71)
    open (unit=73,file='fort.73',status='unknown')
    write (73,*)
    close (73)
    open (unit=74,file='fort.74',status='unknown')
    write (74,*)
    close (74)
    open (unit=62,file='fort.62',status='unknown')
    write (62,*)
    close (62)
    return
    end subroutine rewind2
!*********************************************************************
!*********************************************************************

    subroutine kmin(vpmax,vpmin)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension stem(3)
!*********************************************************************
!                                                                    *
!     Perform parabolic search                                       *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In kmin'
        call timer(65)
        close (65)
    end if
    if (valpar >= zero) then
        stem(1)=valpar-change
        stem(2)=valpar+change
        stem(3)=valpar
    else
        stem(1)=valpar+change
        stem(2)=valpar-change
        stem(3)=valpar
    end if

    if (abs(sdy(1)-sdy(2)) < 1e-5 .AND. abs(sdy(1)-sdy(3)) < 1e-5) &
    then
        temax=stem(3)
        goto 10
    end if

    ste32=stem(3)*stem(3)
    ste22=stem(2)*stem(2)
    ste12=stem(1)*stem(1)
    if ((ste22-ste12) < 1e-5) ste22=ste22+0.0010   !Avoid problems with parabolic search
    b1=(sdy(3)-sdy(1))*(ste22-ste12)-(sdy(2)-sdy(1))*(ste32-ste12)
    b2=(stem(3)-stem(1))*(ste22-ste12)-(stem(2)-stem(1)) &
    *(ste32-ste12)
    b=b1/b2
    a=(sdy(2)-sdy(1)-b*(stem(2)-stem(1)))/(ste22-ste12)
    ac=sdy(1)-a*ste12-b*stem(1)
    temax=-b/(2.0*a)
    expec=0.01*a*temax*temax+0.01*b*temax+0.01*ac
    if (abs(sdy(1)-sdy(2)) < 1e-5 .AND. abs(sdy(1)-sdy(3)) < 1e-5) &
    then
        temax=stem(3)
    end if
    10 open (79,file='fort.79',status='unknown',access='append')
    write (79,120)ichn(1),ichn(2),ichn(3),stem
    write (79,130)sdy
    write (79,140)a,b,ac
    write (79,150)temax
    write (79,155)100*expec
    hu1=1.00-parc2
    hu2=1.00+parc2
!     if (valpar.ge.zero) then
!     if (temax.lt.hu1*stem(1)) temax=hu1*stem(1)
!     if (temax.gt.hu2*stem(2)) temax=hu2*stem(2)
!     end if
!     if (valpar.lt.zero) then
!     if (temax.gt.hu1*stem(2)) temax=hu1*stem(2)
!     if (temax.lt.hu2*stem(1)) temax=hu2*stem(1)
!     end if
    if (valpar > zero) then
        if (temax < hu1*stem(1)) temax=hu1*stem(1)
        if (temax > hu2*stem(2)) temax=hu2*stem(2)
    end if
    if (valpar < zero) then
        if (temax > hu1*stem(1)) temax=hu1*stem(1)
        if (temax < hu2*stem(2)) temax=hu2*stem(2)
    end if

    if (a < zero) then
        iagain=iagain+1
        if (sdy(1) < sdy(2)) then
            temax=stem(1)
        else
            temax=stem(2)
        end if
    end if
    if (temax > vpmax) temax=vpmax
    if (temax < vpmin) temax=vpmin
    if (a > zero) iagain=0
    expec=0.01*a*temax*temax+0.01*b*temax+0.01*ac
    write (79,160)temax
    write (79,165)100*expec
    valnew=temax
    write (79,*)
    close (79)
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    120 format ('Values used for parameter',3i3,/,3d20.10)
    130 format ('Differences found        '/,3d20.10)
    140 format ('Parabol: a=',d20.10,1X,' b=',d20.10,1X,' c=',d20.10)
    150 format ('Minimum of the parabol',d20.10)
    155 format ('Difference belonging to minimum of parabol',d20.10)
    160 format ('New parameter value',d20.10)
    165 format ('Difference belonging to new parameter value',d20.10)
    return
    end subroutine kmin
!*********************************************************************
!*********************************************************************

    subroutine cpfile(nunit1,nunit2)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(80) :: qr1
!*********************************************************************
!                                                                    *
!     Copy file from nunit1 to nunit2                                *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In cpfile'
        call timer(65)
        close (65)
    end if
    rewind (nunit1)
    rewind (nunit2)
    10 read (nunit1,100,err=20,end=20)qr1
    write (nunit2,100)qr1
    goto 10
    20 continue
    rewind (nunit1)
    rewind (nunit2)
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format (a80)
    return
    end subroutine cpfile
!*********************************************************************
!*********************************************************************

    subroutine writegeo(nunit1)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!                                                                    *
!     Copy new geometries to unit nunit1                             *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In writegeo'
        call timer(65)
        close (65)
    end if
    if (axiss(1) < zero) then
        if (nrestra == 0 .AND. nrestrat == 0 .AND. &
        nrestrav == 0) &
        write (nunit1,300)qr,qmol
        if (nrestra > 0) write (nunit1,301)qr, &
        rrstra(1),qmol
        if (nrestrav > 0) write (nunit1,301)qr, &
        vrstra(1),qmol
        if (nrestrat > 0) write (nunit1,301)qr, &
        trstra(1),qmol
    else
        write (nunit1,310)qr,qmol
        write (nunit1,320)axiss(1),axiss(2),axiss(3)
        write (nunit1,320)angles(1),angles(2),angles(3)
    end if
    do i1=1,na
        if (nbiolab /= 2) write (nunit1,400)i1,qa(i1),(c(i1,i2),i2=1,3)
        if (nbiolab == 2) write (nunit1,401)i1,qa(i1),(c(i1,i2),i2=1,3) !Delphi-format
    end do
    if (nbiolab /= 2) write (nunit1,*)
         
    return

    300 format (2x,a1,1x,a60)
    301 format (2x,a1,1x,f6.2,a60)
    310 format (2x,a1,1x,a60)
    320 format (3f10.4)
    400 format (i4,1x,a2,3x,3(d21.14,1x),1x,a5,1x,i5)
    401 format (i3,2x,a2,3x,3(d21.14,1x),1x,a5,1x,i5)
    end subroutine writegeo
!*********************************************************************
!*********************************************************************

    subroutine extractgeo(nmol)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!*********************************************************************
!*********************************************************************
!                                                                    *
!     Extract geometry information from training set                 *
!                                                                    *
!*********************************************************************
    i1=nmol
    na=naset(i1)
    qmol=qmolset(i1)
    qr=qrset(i1)
    qremark(1)=qremset(i1)
    nfc=nfcset(i1)
    iredo=iredoset(i1)
    iexco=iexcoset(i1)
    ncellopt=ncellset(i1)
    iruid=iruidset(i1)
    icell=icellset(i1)
    icelo2=icelo2lset(i1)
    nmmax=nmmaxset(i1)
    ibity=ibityset(i1)
    vvol=vvolset(i1)
    endpo=endposet(i1)
    nrestra=nrestraset(i1)
    icgeo=icgeopt(i1)
    ifreq=ifreqset(i1)
    iremark=1
    do i2=1,nrestra
        rrstra(i2)=rrstraset(i1,i2)
        vkrstr(i2)=vkrstrset(i1,i2)
        vkrst2(i2)=vkrst2set(i1,i2)
        irstra(i2,1)=irstraset(i1,i2,1)
        irstra(i2,2)=irstraset(i1,i2,2)
        rrcha(i1)=rrchaset(i1,i2)
    end do
    nrestrav=nrestravset(i1)
    do i2=1,nrestrav
        irstrav(i2,1)=irstravset(i1,i2,1)
        irstrav(i2,2)=irstravset(i1,i2,2)
        irstrav(i2,3)=irstravset(i1,i2,3)
        vrstra(i2)=vrstraset(i1,i2)
        vkrv(i2)=vkrvset(i1,i2)
        vkr2v(i2)=vkr2vset(i1,i2)
    end do
    nrestrat=nrestratset(i1)
    do i2=1,nrestrat
        irstrat(i2,1)=irstratset(i1,i2,1)
        irstrat(i2,2)=irstratset(i1,i2,2)
        irstrat(i2,3)=irstratset(i1,i2,3)
        irstrat(i2,4)=irstratset(i1,i2,4)
        trstra(i2)=trstraset(i1,i2)
        vkrt(i2)=vkrtset(i1,i2)
        vkr2t(i2)=vkr2tset(i1,i2)
    end do
    do i2=1,3
        axiss(i2)=axisset(i1,i2)
        angles(i2)=anglesset(i1,i2)
    end do
    do i2=1,na
        qa(i2)=qaset(i1,i2)
        chgbgf(i2)=chaset(i1,i2)
        do i3=1,3
            c(i2,i3)=cset(i1,i2,i3)
        end do
    end do

    return
    end subroutine extractgeo
!*********************************************************************
!*********************************************************************

    subroutine writebgf(nunit1)

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension qdir(3)
    character(2) :: qt
    character(1) :: qdir
!*********************************************************************
!                                                                    *
!     Copy new Biograf-geometries to unit nunit1                     *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In newbgf'
        call timer(65)
        close (65)
    end if
    irom=1
    qdir(1)='x'
    qdir(2)='y'
    qdir(3)='z'
    ibgfversion=200
    if (ibity == 1) write (nunit1,1500)ibgfversion
    if (ibity == 2) write (nunit1,1600)ibgfversion
!     if (qr.ne.'F'.and.qr.ne.'5'.and.qr.ne.'Y')
!    $write (nunit1,1500)ibgfversion
!     if (qr.eq.'F'.or.qr.eq.'5'.or.qr.eq.'Y')
!    $write (nunit1,1600)ibgfversion
    write (nunit1,1700)qmol
!     write (nunit1,1700)qkeyw(nprob)
    do i1=1,iremark
        write (nunit1,1800)qremark(i1)
    end do
    qruid='NORMAL RUN'
    if (iruid == 0) then
        write (nunit1,2000)
    else
        if (abs(endpo-endpoold) > 1e-5) write (nunit1,2010)endpo
        if (nmmax /= nmmaxold) write (nunit1,2020)nmmax
        if (nfc /= nfcold) write (nunit1,2030)nfc
        if (ncha /= nchaold) write (nunit1,2036)ncha
        if (iredo > 1) write (nunit1,2035)iredo
        if (ipropt == 1) write (nunit1,2037)vpresopt,vprestax, &
        vprestay,vprestaz
        if (ipropt == 2) write (nunit1,2038)
        if (icell /= icellold) then
            if (icell == 0) write (nunit1,2033)
            if (icell > 0) write (nunit1,2034)ncellopt
        end if
    end if
    if (iexco /= 0 .AND. nsurp > 0) then
        write (nunit1,2040)vvol
        write (nunit1,3500)
        write (nunit1,*)
        return
    end if
    if (nmcharge > 0) then
        do i3=1,nmcharge
            write (nunit1,2050)iat1mc(i3),iat2mc(i3),vmcha(i3)
        end do
    end if

    ims=0
    do i1=1,na
        if (ims == 0 .AND. imove(i1) == 0) then
            if1=i1
            ims=1
        end if
        if (ims == 1 .AND. imove(i1) == 1) then
            write (nunit1,2060)if1,i1-1
            ims=0
        end if
    end do
    if (ims == 1) then
        write (nunit1,2060)if1,na
    end if

!     if (qr.eq.'F'.or.qr.eq.'5'.or.qr.eq.'Y')
    if (ibity == 2) &
    write (nunit1,2100)axiss(1),axiss(2),axiss(3),angles(1), &
    angles(2),angles(3)

    if (neqdis > 0) write (nunit1,2250)
    do i2=1,neqdis
        write (nunit1,2260) &
        ieqdis(i2,1),ieqdis(i2,2),ieqdis(i2,3),ieqdis(i2,4), &
        vkeqd1(i2),vkeqd2(i2)
    end do
    if (nrestra > 0) write (nunit1,2300)
    do i2=1,nrestra
        write (nunit1,2400) &
        irstra(i2,1),irstra(i2,2),rrstra(i2), &
        vkrstr(i2),vkrst2(i2),rrcha(i2),itstart(i2),itend(i2)
    end do

    if (nrestrav > 0) write (nunit1,2500)
    do i2=1,nrestrav
        write (nunit1,2600) &
        irstrav(i2,1),irstrav(i2,2),irstrav(i2,3), &
        vrstra(i2),vkrv(i2),vkr2v(i2),zero
    end do

    if (nrestrat > 0) write (nunit1,2700)
    do i2=1,nrestrat
        write (nunit1,2800) &
        irstrat(i2,1),irstrat(i2,2),irstrat(i2,3), &
        irstrat(i2,4),trstra(i2),vkrt(i2), &
        vkr2t(i2),zero
    end do

    if (nrestram > 0) write (nunit1,2810)
    do i2=1,nrestram
        write (nunit1,2820) &
        qdir(irstram(i2,1)),irstram(i2,2),irstram(i2,3), &
        rmstra1(i2),irstram(i2,4),irstram(i2,5),rmstra2(i2), &
        rmstra3(i2),rmcha(i2)
    end do

    if (icgeo == 0 .AND. ingeo == 0) write (nunit1,2830)
    if (icgeo == 1 .AND. ingeo == 1) write (nunit1,2840)
    if (ifreq == 1) write (nunit1,2850)
    write (nunit1,2900)
    do i2=1,na
        write (nunit1,3000)i2,qa(i2),c(i2,1),c(i2,2),c(i2,3), &
        qa(i2),ibgr1(i2),ibgr2(i2),chgbgf(i2)
    end do
    write (nunit1,3100)
    if (nsurp < 2) then
        do i1=1,na
            write (nunit1,3200)i1,(iag(i1,2+i2),i2=1,iag(i1,2))
        end do
        write (nunit1,3300)
        write (nunit1,3400)estrc
    end if

    write (nunit1,3500)
    write (nunit1,*)
         
    return
    1500 format ('BIOGRF',i4)
    1600 format ('XTLGRF',i4)
    1700 format ('DESCRP ',a60)
    1800 format ('REMARK ',a60)
    1900 format ('FFIELD ',a40)
    2000 format ('RUTYPE NORMAL RUN')
    2010 format ('RUTYPE ENDPO',f6.3)
    2020 format ('RUTYPE MAXIT',i6)
    2030 format ('RUTYPE MAXMOV',i9)
    2033 format ('RUTYPE NO CELL OPT')
    2034 format ('RUTYPE CELL OPT',i6)
    2035 format ('RUTYPE REDO',i6)
    2036 format ('RUTYPE CHARGEMET',i6)
    2037 format ('RUTYPE PRESS OPT',4f8.4)
    2038 format ('RUTYPE PRESS CALC')
    2040 format ('VCHANGE',f8.4)
    2050 format ('MOLCHARGE',2i4,f6.2)
    2060 format ('FIXATOMS',2i6)
    2100 format ('CRYSTX ',6f11.5)
    2200 format ('CELLS ',6i5)
    2250 format ('#              At1 At2 At3 At4    Force1  Force2  ')
    2260 format ('EQUIV DISTANCE ',4i4,f8.2,f8.4)
    2300 format ('#              At1 At2   R12    Force1  Force2  ', &
    'dR12/dIter(MD) Start (MD) End (MD)')
    2400 format ('BOND RESTRAINT ',2i4,f8.4,f8.2,f8.4,1x,f10.7,2i8)
    2500 format ('#               At1 At2 At3 Angle   Force1  Force2', &
    '  dAngle/dIteration (MD only)')
    2600 format ('ANGLE RESTRAINT ',3i4,2f8.2,f8.4,f9.6)
    2700 format ('#                 At1 At2 At3 At3 Angle   Force1  ', &
    'Force2  dAngle/dIteration (MD only)')
    2800 format ('TORSION RESTRAINT ',4i4,2f8.2,f8.4,f9.6)
    2810 format ('#              x/y/z At1 At2    R   At3 At4 Force1', &
    '  Force2  dR/dIteration (MD only)')
    2820 format ('MASCEN RESTRAINT ',a1,1x,2i4,f8.2,2i4,2f8.2,f9.6)
    2830 format ('GEOUPD')
    2840 format ('NO GEOUPD')
    2850 format ('FREQUENCY')
    2900 format ('FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,', &
    '3f15.10,1x,a5,i3,i2,1x,f8.5)')
    3000 format ('HETATM',1x,i5,1x,a2,3x,1x,3x,1x,1x,1x,5x,3f10.5,1x, &
    a5,i3,i2,1x,f8.5)
    3100 format ('FORMAT CONECT (a6,12i6)')
    3200 format ('CONECT',12i6)
    3300 format ('UNIT ENERGY   kcal')
    3400 format ('ENERGY',5x,f14.6)
    3500 format ('END')

    end subroutine writebgf
!*********************************************************************
!*********************************************************************

    subroutine ffsav

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(80) :: qr1
!*********************************************************************
!                                                                    *
!     Save changes in force field                                    *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In ffsav'
        call timer(65)
        close (65)
    end if
    rewind (4)
    write (83,100)sdy(3)
    10 read (4,200,err=20,end=20)qr1
    write (83,200)qr1
    goto 10
    20 continue
    rewind (4)
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format ('Error force field:',f12.4)
    200 format (a80)
    return
    end subroutine ffsav
!*********************************************************************
!*********************************************************************

    subroutine ffext

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension fpar1(7,80,80),fpar2(7,80,80)
    character(40) :: qff1,qff2
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In ffext'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Extrapolate between two force fields                           *
!                                                                    *
!*********************************************************************
    write (6,*)'First force field file (old force field)'
    read (5,'(a40)')qff1
    write (6,*)'Second force field file (new force field)'
    read (5,'(a40)')qff2
    write (6,*)'Extrapolation factor'
    read (5,*)vextr

    open (4,file=qff1,status='old')
    call ffin2(4)
    close (4)
    do i1=1,7
        do i2=1,80
            do i3=1,80
                fpar1(i1,i2,i3)=fpar(i1,i2,i3)
            end do
        end do
    end do

    open (4,file=qff2,status='old')
    call ffin2(4)
    close (4)
    do i1=1,7
        do i2=1,80
            do i3=1,80
                fpar2(i1,i2,i3)=fpar(i1,i2,i3)
            end do
        end do
    end do
          
    do i1=1,7
        do i2=1,80
            do i3=1,80
                change=fpar2(i1,i2,i3)-fpar1(i1,i2,i3)
                fpar(i1,i2,i3)=fpar2(i1,i2,i3)+vextr*change
            end do
        end do
    end do

    call ffchng(5,13,big,-big)
    return
    end subroutine ffext
!*********************************************************************
!*********************************************************************

    subroutine increm

!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension vn(nsort,maxmdat+1),vt(maxmdat,nsort),vtv(nsort,nsort)
    dimension delt(maxmdat+1),sdevhf(maxmdat),vinn(nsort),vb(nsort)
!*********************************************************************
!                                                                    *
!     Optimize heat of formation increments                          *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In increm'
        call timer(65)
        close (65)
    end if
    nso2=nso
!     nso2=2
    gcor=4.0*rgasc*xjouca*0.29815
    do i1=1,iheatf
        delt(i1)=compdat(iheada(i1))-caldat(iheada(i1))-gcor
        sdevhf(i1)=weightdat(iheada(i1))
    end do
    if (nincrop == 0) goto 10   !No increment optimisation
    do i1=1,iheatf
        i3=iheada2(i1)
        do i2=1,nso2
            vn(i2,i1)=float(molin(i3,i2))
        end do
    end do

    do i1=1,iheatf
        do i2=1,nso2
            vt(i1,i2)=vn(i2,i1)
        end do
    end do

    do i3=1,nso2
        do i2=1,nso2
            sum=zero
            do i1=1,iheatf
                sum=sum+vt(i1,i3)*vn(i2,i1)/(sdevhf(i1)**2)
            end do
            vtv(i2,i3)=sum
        end do
    end do

    do i1=1,nso2
        sum=zero
        do i2=1,iheatf
            sum=sum+vt(i2,i1)*delt(i2)/(sdevhf(i2)**2)
        end do
        vb(i1)=sum
    end do

    call matsym5(nso2,nso2,nso2,nso2,vtv,vinn,vb)
    open (24,file='fort.24',status='unknown',access='append')
    write (24,100)(vinn(i1),i1=1,nso2)
    close (24)
    do i1=1,nso
        vincr(i1)=vinn(i1)
        fpar(2,i1,19)=vinn(i1)
    end do
    10 if (nincrop == 0) then
        do i1=1,nso
            vinn(i1)=vincr(i1)
        end do
    end if
    do i1=1,iheatf
        i3=iheada2(i1)
        do i2=1,nso2
            caldat(iheada(i1))=caldat(iheada(i1))+molin(i3,i2)*vinn(i2)
        end do
        caldat(iheada(i1))=caldat(iheada(i1))+gcor
    end do
    return
!***********************************************************************
!                                                                      *
!     Format part                                                      *
!                                                                      *
!***********************************************************************
    100 format (20f10.4)
    end subroutine increm
!*********************************************************************
!*********************************************************************

    subroutine matsym5(n,ndimm,ndimvx,ndimvy,rmat,vecx,vecy)

!*********************************************************************
    include 'cbka.blk'
    DIMENSION RMAT(nsort,nsort),VECX(nsort),VECY(nsort)
!***********************************************************************
!                                                                      *
!     Construction of the lower triangle and backsubstitution          *
!     VERTICAL: ONLY for SYMMETRIC matrices                            *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In matsym5'
        call timer(65)
        close (65)
    end if
    DO  K=1,N
        DIAG=RMAT(K,K)
        DO  L=K+1,N
            VECX(L)=RMAT(L,K)
            RMAT(L,K)=RMAT(L,K)/DIAG
        END DO
        VECY(K)=VECY(K)/DIAG
        RMAT(K,K)=1.0D0
        DO  L=K+1,N
            FACTOR=VECX(L)
            DO  M=L,N
                RMAT(M,L)=RMAT(M,L)-FACTOR*RMAT(M,K)
            END DO
            VECY(L)=VECY(L)-FACTOR*VECY(K)
       END DO
    ENDDO
    DO  K=N,1,-1
        VECX(K)=VECY(K)
        DO  L=K+1,N
            VECX(K)=VECX(K)-RMAT(L,K)*VECX(L)
        END DO
    END DO

    RETURN

    end subroutine matsym5

!***********************************************************************
!*********************************************************************

    double precision function gauss(sigma,v0,dseed)

!*********************************************************************
    implicit double precision (a-h,o-z),integer (i-n)
    1 vr1=2.0*random(dseed)-1.0
    vr2=2.0*random(dseed)-1.0
    r=vr1*vr1+vr2*vr2
    if (r >= 1.0) goto 1
    fac=vr1*sqrt(-2.0*log(r)/r)
    gauss=v0+sigma*fac
    return
    end function gauss
!*********************************************************************
!*********************************************************************

    double precision function random(dseed)

!*********************************************************************
    implicit double precision (a-h,o-z),integer (i-n)
    data d2p31m/2147483647.d0/
    data d2p31/2147483711.d0/
    dseed=dmod(16807.d0*dseed,d2p31m)
    random=dseed/d2p31
    return
    end function random
!*********************************************************************
