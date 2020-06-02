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
     
    subroutine vibra(qfreqfile,vibreax,vibqc,imatch,errmatch, &
    itrain,klinear)
     
!*********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension vmode(navib*3,navib*3),vmodqc(navib*3,navib*3)
    dimension vibreax(navib*3),vibqc(navib*3),errmatch(3*navib)
    dimension errmatch2(3*navib)
    dimension csav(nat,3)
    dimension hulp1(navib*3),hulp2(navib*3)
    dimension htrm(3)
    dimension imatch(3*navib),imatch2(3*navib)
    character(40) :: qfreqfile
    character(200) :: qhulp
!*********************************************************************
!                                                                    *
!     Calculate vibrational frequencies                              *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In vibra'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Calculate vibrational frequencies                              *
!                                                                    *
!*********************************************************************
    CLIGHT=2.99792458D8        !             [m/s]
    conrcm=sqrt(caljou)*1d11/(2.0*pi*clight)
    if (na > navib) stop 'Nr. of atoms > navib; stop in vibra'
    do i1=1,na*3
        do i2=1,na*3
            d2(i1,i2)=zero
        end do
    end do
!*********************************************************************
!                                                                    *
!     Calculate numerical second derivatives                         *
!                                                                    *
!*********************************************************************
    delta=delvib
    ico=0
    do i1=1,na
        do i2=1,3
            csav(i1,i2)=c(i1,i2)
        end do
    end do

     do i1=1,na
        if (ndebug == 1) then
            open (51,file='thermo.out',status='unknown',access='append')
            write (51,*)i1
            close (51)
        end if
        do  k1=1,3
            ih1=3*(i1-1)+k1
            do  i2=i1,na
                do  k2=1,3
                    ihulp=0
                    if (i1 == i2 .AND. k1 == k2) ihulp=1
                    ih2=3*(i2-1)+k2

                    call distan
                    call vlist
                    call intcor
                    call boncor
                    call charvib
                    call covbon
                    call lonpar
                    call ovcor
                    call calval
                    call valang
                    call hbond
                    call torang
                    call nonbon
                    estrc=eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+ncha2*ech+efi
                    e0=estrc
                    ico=ico+1
                    if (ihulp == 0) then
                        c(i1,k1)=csav(i1,k1)+delta
                        c(i2,k2)=csav(i2,k2)+delta
                    else
                        c(i1,k1)=csav(i1,k1)+2.0d0*delta
                    end if

                    call distan
                    call vlist
                    call intcor
                    call boncor
                !     call charvib
                    call covbon
                    call lonpar
                    call ovcor
                    call calval
                    call valang
                    call hbond
                    call torang
                    call nonbon
                    estrc=eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+ncha2*ech+efi
                    e1=estrc
                    e2=zero
                    e3=zero

                    if (ihulp == 0) then
                        c(i1,k1)=csav(i1,k1)+delta
                        c(i2,k2)=csav(i2,k2)-delta
                    else
                        c(i1,k1)=csav(i1,k1)
                    end if

                    call distan
                    call vlist
                    call intcor
                    call boncor
                !     call charvib
                    call covbon
                    call lonpar
                    call ovcor
                    call calval
                    call valang
                    call hbond
                    call torang
                    call nonbon
                    estrc=eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+ncha2*ech+efi
                    e2=estrc

                    if (ihulp == 0) then
                        c(i1,k1)=csav(i1,k1)-delta
                        c(i2,k2)=csav(i2,k2)+delta
                    else
                        c(i1,k1)=csav(i1,k1)
                    end if

                    call distan
                    call vlist
                    call intcor
                    call boncor
                !     call charvib
                    call covbon
                    call lonpar
                    call ovcor
                    call calval
                    call valang
                    call hbond
                    call torang
                    call nonbon
                    estrc=eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+ncha2*ech+efi
                    e3=estrc
                    if (ihulp == 0) then
                        c(i1,k1)=csav(i1,k1)-delta
                        c(i2,k2)=csav(i2,k2)-delta
                    else
                        c(i1,k1)=csav(i1,k1)-2.0d0*delta
                    end if

                    call distan
                    call vlist
                    call intcor
                    call boncor
                !     call charvib
                    call covbon
                    call lonpar
                    call ovcor
                    call calval
                    call valang
                    call hbond
                    call torang
                    call nonbon
                    estrc=eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+ncha2*ech+efi
                    e4=estrc
                    if (ndebug >= 1) then
                        open (53,file='fort.53',status='unknown',access='append')
                        write (53,'(4i6,7f16.4)')i1,i2,k1,k2,e0,e1,e2,e3,e4
                        close (53)
                    end if
                    d2h=((e1-e2)-(e3-e4))/(4.0d0*delta*delta)
                    d2(ih1,ih2)=d2h/sqrt(xmasat(i1)*xmasat(i2))
                    d2(ih2,ih1)=d2h/sqrt(xmasat(i1)*xmasat(i2))

                    c(i1,k1)=csav(i1,k1)
                    c(i2,k2)=csav(i2,k2)
    enddo
    enddo
    enddo
    enddo
    if (ndebug >= 1) then
        open (54,file='fort.54',status='unknown',access='append')
        write (54,*)'Numerical 2nd derivatives:'
        do i1=1,na*3
            write(54,'(200f8.2)')(d2(i1,i2),i2=1,na*3)
        end do
        close (54)
    end if
!*********************************************************************
!                                                                    *
!     Diagonalize matrix                                             *
!                                                                    *
!*********************************************************************
!     call diag(d2,vmode,navib*3,vibreax)
    call EIGVAL(navib*3,na*3,d2,vibreax,hulp1,hulp2,ierr)
    do i1=1,na*3
    end do

    do i1=1,na*3
        do i2=1,na*3
            iathu=int(float(i1-1)/3.0)+1
            vmode(i1,i2)=d2(i1,i2)/sqrt(xmasat(iathu))
        end do
    end do

    if (iopt == 0 .AND. nsurp == 0) then
        open (51,file='thermo.out',status='unknown',access='append')
        write (51,*)qmol
    end if
    do i1=1,3*na
        nr=1
        if (vibreax(i1) < zero) nr=-1
        vibreax(i1)=float(nr)*sqrt(abs(vibreax(i1)))*conrcm
        if (iopt == 0 .AND. nsurp == 0) then
            write (51,*)i1,vibreax(i1)
            do i2=1,na
                write (51,'(a2,6(f10.5))')qa(i2), &
                (c(i2,i3),i3=1,3),(vmode((i2-1)*3+i3,i1),i3=1,3)
            end do
        end if
    end do

    CLIGHT=2.99792458D8        !             [m/s]
    HPLNCK=6.6260755D-34       ! [Js]      = [Nms]
    ABSNUL=-273.15D0           !             [K]
    TEMPC=25.0D0               ! [deg C]   = [K]
    PRESS=ONE                  !

    TEMPK=TEMPC-ABSNUL
    HAENCE=HPLNCK*avognr*CLIGHT
    HNC100=1.0D2*HAENCE
    call INERTI(htrm)
    XIXYZ=htrm(1)+htrm(2)+htrm(3)
!     XIXYZ=5.0
    K1=-1
    klinear=0
    if (htrm(1) < 1e-5 .OR. htrm(2) < 1e-5 .OR. htrm(3) < 1e-5) &
    then
        klinear=1  !Linear molecule
        if (iopt == 0 .AND. nsurp == 0) then
            write (51,*)'Linear molecule'
        end if
    end if

    if (iopt == 0 .AND. nsurp == 0) then
        WRITE(51,600)PRESS

        35 ERTE=TEMPK*rgasc
        ENULP=0.0
        HTCNT=0.0
        CEPE=0.0
        SVIBR1=0.0
        SVIBR2=0.0
        DO 7 K=1,klinear+3*na-6
            kmod=k+6-klinear
            ENULP=ENULP+vibreax(kmod)
            ARGUM=HNC100*vibreax(kmod)/ERTE
            IF (ARGUM > 85.0D0) GOTO 6
            EMACHT=EXP(ARGUM)
            EMM1=EMACHT-ONE
            FRGEM=vibreax(kmod)/EMM1
            HTCNT=HTCNT+FRGEM
            CEPE=CEPE+FRGEM*EMACHT*FRGEM
            SVIBR1=SVIBR1+FRGEM+vibreax(kmod)
            SVIBR2=SVIBR2+LOG(EMM1)
        7 END DO
         
        ENULP=HAENCE*ENULP*xjouca*0.5D-1
        HTCNT=xjouca*(4.0D-3*ERTE+1.0D-1*HAENCE*HTCNT)
        CEPE=xjouca*(4.0D0*rgasc+1.0D4*CEPE*(HAENCE/TEMPK)**2/rgasc)
        SVIBR=xjouca*(1.0D2*HAENCE*SVIBR1/TEMPK-rgasc*SVIBR2)
        SROT=rgasc*xjouca*(1.5D0+1.5D0*LOG(8.0D0*pi**(7.0D0/3.0D0)*ERTE/( &
        (avognr*HPLNCK)**2))+0.5D0*LOG(XIXYZ)-34.5D0*LOG(10.0D0))
        SIGMA=float(isymm)
        SROT=SROT-xjouca*rgasc*LOG(SIGMA)
        STRAN=xjouca*(2.5D0*rgasc+rgasc*LOG((2.0D0*pi*xmasmd*ERTE)**1.5D0 &
        *ERTE/((avognr*HPLNCK)**3*avognr*PRESS*1.01325D0*10.0D0**9.5D0)))
         
    !     IF (NCL.GT.0) THEN
    !     HTCNT=HTCNT-xjouca*4.0D-3*ERTE
    !     CEPE=CEPE-xjouca*4.0D0*rgasc
    !     SROT=ZERO
    !     STRAN=ZERO
    !     ENDIF
         
        STOT=STRAN+SROT+SVIBR
         
        TEDS=STOT*TEMPK*1.0D-3
        WRITE(51,700)TEMPK,ENULP,HTCNT,TEDS,CEPE,SVIBR,SROT,STRAN,STOT
        IF(K1 == -1) THEN
            WRITE(*,1100)
        END IF
        6 K1=K1+1
        TEMPK=23.15+K1*100.0D0
        IF(TEMPK < 5000.0)GOTO 35
        WRITE(51,800)

    end if

    if (iopt == 0 .AND. nsurp == 0) close (51)
    if (itrain == 0) return
!*********************************************************************
!                                                                    *
!     Read Jaguar vibrational data                                   *
!                                                                    *
!*********************************************************************
    ifh=0
    nufreq=klinear+3*na-6
    open (19,file=qfreqfile,status='old',err=60)
    !read (19,'(a200)')qhulp
    !read (19,'(a200)')qhulp
    !read (19,'(a200)')qhulp
    ijagvers=1
    !if (qhulp(21:23) == '4.2') ijagvers=2

    15 read (19,'(a200)',err=60,end=50)qhulp

    !if (ijagvers == 1) then
        ihulp=6
        if (ifh+6 > nufreq) ihulp=nufreq-ifh
        if (qhulp(3:13) == 'frequencies') then
            read (qhulp,'(17x,6(f8.2,1x))')(vibqc(i1+ifh),i1=1,ihulp)
            !read (19,*)
            !read (19,*)
            !do i1=1,na
            !    do i2=1,3
            !        read (19,'(17x,6(f8.2,1x))') &
            !        (vmodqc(i2+(i1-1)*3,i3+ifh),i3=1,ihulp)
            !    end do
            !end do
            ifh=ifh+6
        end if
        goto 15
    !end if

    !if (ijagvers == 2) then
    !    ihulp=7
    !    if (ifh+7 > nufreq) ihulp=nufreq-ifh
    !    if (qhulp(3:13) == 'frequencies') then
    !        read (qhulp,'(13x,7(f8.2,1x))')(vibqc(i1+ifh),i1=1,ihulp)
    !        read (19,*)
    !        do i1=1,na
    !            do i2=1,3
    !                read (19,'(13x,7(f8.2,1x))') &
    !                (vmodqc(i2+(i1-1)*3,i3+ifh),i3=1,ihulp)
    !            end do
    !        end do
    !        ifh=ifh+7
    !    end if
    !    goto 15
    !end if
     
    50 close (19)
!*********************************************************************
!                                                                    *
!     Compare Reax modes to Jaguar modes                             *
!                                                                    *
!*********************************************************************
    do i1=1,navib*3
        !errmatch2(i1)=5e+35
        !imatch(i1)=-1
        !imatch2(i1)=-1
        errmatch2(i1)=0.0d0
        imatch(i1)=i1
        imatch2(i1)=i1
    end do
    nagain=0

    !100 continue
    !nagain=nagain+1
    !iagain=0
    !do i1=1,klinear+3*na-6
    !    diffmin=5e+35
    !    if (imatch(i1) < 0) then
    !        do i2=1,klinear+3*na-6
    !            diffp=0.0
    !            diffm=0.0
    !            do i3=1,3*na
    !                diffph=vmode(i3,i1+6-klinear)-vmodqc(i3,i2)
    !                diffmh=vmode(i3,i1+6-klinear)+vmodqc(i3,i2)
    !                diffp=diffp+diffph*diffph*diffph*diffph
    !                diffm=diffm+diffmh*diffmh*diffmh*diffmh
    !            end do
    !            diffp=diffp/float(3*na)
    !            diffm=diffm/float(3*na)
    !            if (diffp < diffmin) then
    !                if (diffp < errmatch2(i2)) then
    !                    diffmin=diffp
    !                    imatch(i1)=i2
    !                    errmatch(i1)=diffp
    !                    if (imatch2(i2) < i1) then
    !                        errmatch(imatch2(i2))=5e+35
    !                        imatch(imatch2(i2))=-1
    !                        iagain=1
    !                    end if
    !                    imatch2(i2)=i1
    !                    errmatch2(i2)=diffp
    !                end if
    !            end if
    !            if (diffm < diffmin) then
    !                if (diffm < errmatch2(i2)) then
    !                    diffmin=diffm
    !                    imatch(i1)=i2
    !                    errmatch(i1)=diffm
    !                    if (imatch2(i2) < i1) then
    !                        errmatch(imatch2(i2))=5e+35
    !                        imatch(imatch2(i2))=-1
    !                        iagain=1
    !                    end if
    !                    imatch2(i2)=i1
    !                    errmatch2(i2)=diffm
    !                end if
    !            end if
    !        end do
    !    end if
    !end do
    ! 
    !do i1=1,navib*3
    !    errmatch2(i1)=5e+35
    !    imatch2(i1)=-1
    !end do

    !do i1=1,klinear+3*na-6
    !    imatch2(imatch(i1))=i1
    !    errmatch2(imatch(i1))=errmatch(i1)
    !end do

    !if (iopt == 0 .AND. nsurp < 2) then
    !    write (61,*)qmol,nagain,iagain
    !     
    !    do i1=1,klinear+3*na-6
    !        write (61,*)i1,imatch(i1),errmatch(i1)
    !        write (61,*)vibqc(imatch(i1)),vibreax(i1+6-klinear)
    !        do i2=1,na
    !            write (61,'(i4,a2,6(f8.2,1x))')i2,qa(i2), &
    !            (vmodqc((i2-1)*3+i3,imatch(i1)),i3=1,3), &
    !            (vmode((i2-1)*3+i3,i1+6-klinear),i3=1,3)
    !        end do
    !    end do
    !end if

    !if (iagain == 1 .AND. nagain < 15) goto 100

    !if (imfreq == -1) then   !Do not match modes
    !    do i1=1,klinear+3*na-6
    !        imatch(i1)=i1
    !        errmatch(i1)=0.0
    !    end do
    !end if
     
    return
    60 write (*,*)'Could not open Jaguar vibrational data'
    return
    600 FORMAT(1X,'THERMODYNAMIC QUANTITIES OF IDEAL GAS PHASE AT ',F6.2, &
    ' atm AND TEMPERATURE (K) AS INDICATED:'//T4,'TKELVIN',T18,'A',T33 &
    ,'B',T48,'C',T63,'D',T78,'E',T93,'F',T108,'G',T123,'H')
    700 FORMAT(1(1X,F9.3,3X,8(1X,F14.8)))
    800 FORMAT(/1X,'A = ZERO POINT ENERGY',T32,'B = HEAT CONTENT AT TKELVI &
    N',T62,'C = ENTROPY*TKELVIN',T92,'D = SPECIFIC HEAT CP'/1X,'E = EN &
    TROPY OF VIBRATION',T32,'F = ENTROPY OF ROTATION',T62,'G = ENTROPY &
    OF TRANSLATION',T92,'H = TOTAL ENTROPY E+F+G+CORRECTIONS'/1X,'A,B &
    AND C: kcal/mol, OTHER QUANTITIES: cal/(mol*K)'//1X,'THE ENTROPY &
    OF ROTATION IS CORRECTED FOR THE SYMMETRY NUMBER; NOTE THAT THE EN&
    &TROPY OF A 1:1 MIXTURE OF ENANTIOMERS IS RLN2 LARGER')
    1100 FORMAT(1X)
    end subroutine vibra
!***********************************************************************
!*********************************************************************

    subroutine charvib

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
    dimension elcvec(nat),char(nat)
!*********************************************************************
!                                                                    *
!     Determine charges on atoms: full system exact solution         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',access='append')
        write (65,*) 'In charvib'
        call timer(65)
        close (65)
    end if

    ech=zero
    cfield=332.0638

    do i2=1,na+1      !Zero EEM-matrix
        do i3=1,na+1
            xmortr(i3,i2)=zero
        end do
    end do

    do i2=1,na+1
        elcvec(i2)=zero
    end do
    chamol=syscha    !System charge

    do i1=1,na
        xmortr(i1,i1)=2.0*eta(ia(i1,1))  !Mortier method
        elcvec(i1)=-chi(ia(i1,1))        !Mortier method
              
        xmortr(na+1,i1)=1.0
    end do
    elcvec(na+1)=chamol     !Charge on molecule

    do 10 ivl=1,nvpair-nvlself      !Verlet-list
        i2=nvl1(ivl)
        i3=nvl2(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        dx=c(i2,1)-c(i3,1)+ix*tm11
        dy=c(i2,2)-c(i3,2)+ix*tm21+iy*tm22
        dz=c(i2,3)-c(i3,3)+ix*tm31+iy*tm32+iz*tm33
        dis2=sqrt(dx*dx+dy*dy+dz*dz)

        ity1=ia(i2,1)
        ity2=ia(i3,1)
        gamt=sqrt(gam(ity1)*gam(ity2))

        if (dis2 < swb .AND. dis2 > 0.001) then
            call taper(dis2,dis2*dis2)
            hulp1=(dis2**3+(1.0/(gamt**3)))
            hulp2=hulp1**(1.0/3.0)
            xmortr(i3,i2)=xmortr(i3,i2)+sw*14.40/hulp2

        end if

    10 END DO

    call matsym4(na+1,na+1,na+1,na+1,na+1,xmortr,char,elcvec)

    ech=zero
    do i2=1,na
        ch(i2)=char(i2)
        ech=ech+23.02*(chi(ia(i2,1))*ch(i2)+ &
        eta(ia(i2,1))*ch(i2)*ch(i2))
    end do
    chisys=char(na+1)

    return
    end subroutine charvib
!*********************************************************************
!***********************************************************************
     
    SUBROUTINE INERTI(htrm)
     
!***********************************************************************
!                                                                      *
!     CALCULATION OF MATRIX OF MOMENTS OF INERTIA                      *
!                                                                      *
!     TRANSFORMATION OF COORDINATES SUCH THAT PRINCIPLE AXES OF        *
!     INERTIA COINCIDE WITH CARTESIAN COORDINATE SYSTEM AXES           *
!                                                                      *
!     CALCULATION OF TOLERANCIES IN VALUES OF PRINCIPLE MOMENTS OF     *
!     INERTIA AND IN DIFFERENCES BETWEEN THESE VALUES                  *
!                                                                      *
!***********************************************************************
     
    INCLUDE 'cbka.blk'
    INCLUDE 'opt.blk'
     
    DIMENSION ctemp(nat,3)
    DIMENSION ELVEC(3,3),TRM(3,3),HULP1(3),HULP2(3),htrm(3)
     
    EQUIVALENCE (TRM(1,1),ELVEC(1,1))
     
    do i1=1,3
        do i2=1,3
            TRM(i1,i2)=0.0
        end do
    end do

!***********************************************************************
!                                                                      *
!     Place system at origin                                           *
!                                                                      *
!***********************************************************************
    ccx=0.0
    ccy=0.0
    ccz=0.0
    do i1=1,na
        ccx=ccx+c(i1,1)*xmasat(i1)/xmasmd
        ccy=ccy+c(i1,2)*xmasat(i1)/xmasmd
        ccz=ccz+c(i1,3)*xmasat(i1)/xmasmd
    end do
    xt2=-ccx
    yt2=-ccy
    zt2=-ccz
    do i1=1,na
        ctemp(i1,1)=c(i1,1)+xt2
        ctemp(i1,2)=c(i1,2)+yt2
        ctemp(i1,3)=c(i1,3)+zt2
    end do

    DO 1 I=1,na
        XMA=xmasat(I)
        CX=ctemp(I,1)
        CY=ctemp(I,2)
        CZ=ctemp(I,3)
        TRM(1,1)=TRM(1,1)+XMA*(CY*CY+CZ*CZ)
        TRM(2,2)=TRM(2,2)+XMA*(CX*CX+CZ*CZ)
        TRM(3,3)=TRM(3,3)+XMA*(CX*CX+CY*CY)
        TRM(2,1)=TRM(2,1)-XMA*CX*CY
        TRM(3,1)=TRM(3,1)-XMA*CX*CZ
        TRM(3,2)=TRM(3,2)-XMA*CY*CZ
    1 END DO
     
!***********************************************************************
!                                                                      *
!     Calculation of principle moments of inertia.                     *
!                                                                      *
!***********************************************************************
     
    CALL EIGVAL(3,3,TRM,HTRM,HULP1,HULP2,IERR)
!     IF (IERR.NE.0) WRITE(*,200)IERR
    if (iopt == 0 .AND. nsurp == 0) then
        WRITE(51,100)htrm(1),htrm(2),htrm(3)
    end if
     
!***********************************************************************
!                                                                      *
!     Transformation of coordinates such that principle axes of        *
!     inertia coincide with cartesian coordinate system axes.          *
!                                                                      *
!     ELVEC=TRM                                                        *
!                                                                      *
!***********************************************************************

!     DELI1=ZERO
!     DELI2=ZERO
!     DELI3=ZERO

!     DO 2 I=1,na
!     CX=c(I,1)
!     CY=c(I,2)
!     CZ=c(I,3)
!     c(I,1)=ELVEC(1,1)*CX+ELVEC(2,1)*CY+ELVEC(3,1)*CZ
!     c(I,2)=ELVEC(1,2)*CX+ELVEC(2,2)*CY+ELVEC(3,2)*CZ
!     c(I,3)=ELVEC(1,3)*CX+ELVEC(2,3)*CY+ELVEC(3,3)*CZ
!     XMA=XMASAT(I)
!     DELI1=DELI1+2*EPSC*XMA*(ABS(C(I,2))+ABS(C(I,3)))
!     DELI2=DELI2+2*EPSC*XMA*(ABS(C(I,1))+ABS(C(I,3)))
!   2 DELI3=DELI3+2*EPSC*XMA*(ABS(C(I,1))+ABS(C(I,2)))

!***********************************************************************
!                                                                      *
!     Calculation of tolerancies in values of principle moments of     *
!     intertia and in differences between these values.                *
!                                                                      *
!***********************************************************************

!     D32=DELI3+DELI2
!     D21=DELI2+DELI1
     
    RETURN
     
!***********************************************************************
!                                                                      *
!     FORMAT PART                                                      *
!                                                                      *
!***********************************************************************
     
    100 FORMAT(' PRINCIPAL MOMENTS OF INERTIA IN ATOMIC MASS UNITS AND ANG &
    STROMS:',3F15.8)
    200 FORMAT(1X,'IERR<>0 in INERTI !!!!!!!!!. Ask authors. IERR=',I4)
     
    END SUBROUTINE INERTI
     
!***********************************************************************
!***********************************************************************

    subroutine diag(C,A,NP,D)

!***********************************************************************
    implicit real*8 (a-h,o-z)
!		include "maxpar.inc"
!     Driver for routine TQLI
!***********************************************************************
    PARAMETER(TINY=1.0E-6)
    double precision :: A(NP,NP),C(NP,NP), &
    D(NP),E(NP),F(NP)
!      DATA C/
!     *            4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,
!     *            3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,
!     *            2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,
!     *            1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,
!     *            0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,
!     *            -1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,
!     *            -2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,
!     *            -3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,
!     *            -4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0/

!	write(6,*)'np', np
!		open(50,file='vel_corr.dat',status='old')
!	do i=1,np
!		write(6,*)(c(i,j),j=1,NP)
!	enddo
!		close(50)
    DO 12 I=1,NP
        DO 11 J=1,NP
            A(I,J)=C(I,J)
        11 END DO
    12 END DO
    CALL TRED2(A,NP,NP,D,E)
    CALL TQLI(D,E,NP,NP,A)
!     WRITE(*,'(/1X,A)') 'Eigenvectors for a real symmetric matrix'
    DO 16 I=1,NP
        DO 14 J=1,NP
            F(J)=0.0
            DO 13 K=1,NP
                F(J)=F(J)+C(J,K)*A(K,I)
            13 END DO
        14 END DO
    !           WRITE(*,'(/1X,A,I3,A,F10.6)') 'Eigenvalue',I,' =',D(I)
    !           WRITE(*,'(/1X,T7,A,T17,A,T31,A)') 'Vector','Mtrx*Vect.',
    !    &			'Ratio'
        DO 15 J=1,NP
        !                 IF (ABS(A(J,I)).LT.TINY) THEN
        !                       WRITE(*,'(1X,2F12.6,A12)') A(J,I),F(J),
        !    &						'div. by 0'
        !                 ELSE
        !                       WRITE(*,'(1X,2F12.6,E14.6)') A(J,I),F(J),
        !    *                              F(J)/A(J,I)
        !                 ENDIF
        15 END DO
    !            WRITE(*,'(/1X,A)') 'press ENTER to continue...'
    !            READ(*,*)
    16 END DO
    end subroutine diag

!***********************************************************************
!***********************************************************************

    SUBROUTINE TRED2(A,N,NP,D,E)

!***********************************************************************
    implicit real*8 (a-h,o-z)
    double precision :: A(NP,NP),D(NP),E(NP)
    IF(N > 1)THEN
        DO 18 I=N,2,-1
            L=I-1
            H=0.
            SCALE=0.
            IF(L > 1)THEN
                DO 11 K=1,L
                    SCALE=SCALE+ABS(A(I,K))
                11 END DO
                IF(SCALE == 0.)THEN
                    E(I)=A(I,L)
                ELSE
                    DO 12 K=1,L
                        A(I,K)=A(I,K)/SCALE
                        H=H+A(I,K)**2
                    12 END DO
                    F=A(I,L)
                    G=-SIGN(SQRT(H),F)
                    E(I)=SCALE*G
                    H=H-F*G
                    A(I,L)=F-G
                    F=0.
                    DO 15 J=1,L
                        A(J,I)=A(I,J)/H
                        G=0.
                        DO 13 K=1,J
                            G=G+A(J,K)*A(I,K)
                        13 END DO
                        IF(L > J)THEN
                            DO 14 K=J+1,L
                                G=G+A(K,J)*A(I,K)
                            14 END DO
                        ENDIF
                        E(J)=G/H
                        F=F+E(J)*A(I,J)
                    15 END DO
                    HH=F/(H+H)
                    DO 17 J=1,L
                        F=A(I,J)
                        G=E(J)-HH*F
                        E(J)=G
                        DO 16 K=1,J
                            A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                        16 END DO
                    17 END DO
                ENDIF
            ELSE
                E(I)=A(I,L)
            ENDIF
            D(I)=H
        18 END DO
    ENDIF
    D(1)=0.
    E(1)=0.
    DO 23 I=1,N
        L=I-1
        IF(D(I) /= 0.)THEN
            DO 21 J=1,L
                G=0.
                DO 19 K=1,L
                    G=G+A(I,K)*A(K,J)
                19 END DO
                DO 20 K=1,L
                    A(K,J)=A(K,J)-G*A(K,I)
                20 END DO
            21 END DO
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L >= 1)THEN
            DO 22 J=1,L
                A(I,J)=0.
                A(J,I)=0.
            22 END DO
        ENDIF
    23 END DO
    RETURN
    END SUBROUTINE TRED2

!***********************************************************************
!***********************************************************************

    SUBROUTINE TQLI(D,E,N,NP,Z)

!***********************************************************************
    implicit real*8 (a-h,o-z)
    double precision :: D(NP),E(NP),Z(NP,NP)
    IF (N > 1) THEN
        DO 11 I=2,N
            E(I-1)=E(I)
        11 END DO
        E(N)=0.
        DO 15 L=1,N
            ITER=0
            1 DO 12 M=L,N-1
                DD=ABS(D(M))+ABS(D(M+1))
                IF (ABS(E(M))+DD == DD) GO TO 2
            12 END DO
            M=N
            2 IF(M /= L)THEN
                IF(ITER == 30) stop 'too many iterations'
                ITER=ITER+1
                G=(D(L+1)-D(L))/(2.*E(L))
                R=SQRT(G**2+1.)
                G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
                S=1.
                C=1.
                P=0.
                DO 14 I=M-1,L,-1
                    F=S*E(I)
                    B=C*E(I)
                    IF(ABS(F) >= ABS(G))THEN
                        C=G/F
                        R=SQRT(C**2+1.)
                        E(I+1)=F*R
                        S=1./R
                        C=C*S
                    ELSE
                        S=F/G
                        R=SQRT(S**2+1.)
                        E(I+1)=G*R
                        C=1./R
                        S=S*C
                    ENDIF
                    G=D(I+1)-P
                    R=(D(I)-G)*S+2.*C*B
                    P=S*R
                    D(I+1)=G+P
                    G=C*R-B
                    DO 13 K=1,N
                        F=Z(K,I+1)
                        Z(K,I+1)=S*Z(K,I)+C*F
                        Z(K,I)=C*Z(K,I)-S*F
                    13 END DO
                14 END DO
                D(L)=D(L)-P
                E(L)=G
                E(M)=0.
                GO TO 1
            ENDIF
        15 END DO
    ENDIF
    RETURN
    END SUBROUTINE TQLI
!***********************************************************************
!***********************************************************************
     
    SUBROUTINE EIGVAL(NDIM,N,A,W,FV1,FV2,IERR)
     
!***********************************************************************
!                                                                      *
!     THIS ROUTINE FINDS EIGENVALUES AND EIGENVECTORS (IF DESIRED)     *
!     OF A REAL SYMMETRIC MATRIX                                       *
!                                                                      *
!     INPUT:  NDIM IS THE DIMENSION OF THE REAL SYMMETRIC MATRIX  A    *
!             N DENOTES THE PART OF THE MATRIX A WHICH IS USED         *
!                                                                      *
!             MATZ EQUAL TO ZERO ONLY EIGENVALUES ARE CALCULATED       *
!             MATZ NOT EQUAL TO ZERO BOTH EIGENVALUES AND EIGENVECTORS *
!             ARE CALCULATED                                           *
!                                                                      *
!     OUTPUT: W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER           *
!             A  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO         *
!             IERR  IS AN ERROR COMPLETION CODE, THE NORMAL COMPLETION *
!             CODE IS ZERO, FOR OTHER CODES SEE DOCUMENTATION          *
!             ROUTINES TOLRAT AND TQL2                                 *
!                                                                      *
!     FV1 AND FV2 ARE TEMPORARY STORAGE ARRAYS                         *
!                                                                      *
!***********************************************************************
     
    INCLUDE 'cbka.blk'
     
    DIMENSION A(NDIM,NDIM),W(NDIM),FV1(NDIM),FV2(NDIM)
     
!***********************************************************************
!                                                                      *
!     FIND BOTH EIGENVALUES AND EIGENVECTORS                           *
!                                                                      *
!***********************************************************************
     
    1 CALL TRED2B(NDIM,N,A,W,FV1)
    CALL TQL2(NDIM,N,W,FV1,A,IERR)
     
    RETURN
     
    END SUBROUTINE EIGVAL
     
!***********************************************************************
!***********************************************************************
     
    SUBROUTINE TRED2B(NDIM,N,Z,D,E)
     
!***********************************************************************
!                                                                      *
!     THIS ROUTINE IS BASED ON THE ALGOL PROCEDURE TRED2               *
!     NUM.MATH.11,181-195(1968) BY MARTIN, REINSCH, AND WILKINSON      *
!     HANDBOOK FOR AUTO.COMP.,VOL.II-LINEAR ALGEBRA,212-226(1971)      *
!                                                                      *
!     THIS ROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A SYMMETRIC      *
!     TRIDIAGONAL MATRIX USING AND ACCUMULATING ORTHOGONAL SIMILARITY  *
!     TRANSFORMATIONS                                                  *
!                                                                      *
!     INPUT:  N IS THE DIMENSION OF THE REAL SYMMETRIC MATRIX A        *
!             (ONLY THE LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED) *
!                                                                      *
!     OUTPUT: D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL      *
!             MATRIX                                                   *
!             E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL   *
!             MATRIX IN ITS LAST N-1 POSITIONS, E(1) IS SET TO ZERO    *
!             Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX          *
!             PRODUCED IN THE REDUCTION                                *
!             A IS UNALTERED                                           *
!             E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF *
!             E                                                        *
!                                                                      *
!***********************************************************************
     
    implicit double precision (a-h,o-z),integer(i-n)
    DIMENSION D(NDIM),E(NDIM),Z(NDIM,NDIM)
     
!***********************************************************************
!                                                                      *
!                                                                      *
!***********************************************************************
     
    ZERO=0.0
    ONE=1.0
    DO 1 I=1,N
        D(I)=Z(N,I)
    1 END DO
    IF(N == 1)GOTO 3
     
!***********************************************************************
!                                                                      *
!     FOR I=N STEP -1 UNTIL 2 DO -                                     *
!                                                                      *
!***********************************************************************
     
    DO 4 II= 2,N
        I=N+2-II
        L=I-1
        H=ZERO
        SCALE=ZERO
        IF(L < 2)GOTO 5
         
    !***********************************************************************
    !                                                                      *
    !     SCALE ROW (ALGOL TOL THEN NOT NEEDED)                            *
    !                                                                      *
    !***********************************************************************
        DO 6 K=1,L
            SCALE=SCALE+ABS(D(K))
        6 END DO
        IF(SCALE /= ZERO)GOTO 7
        5 E(I)=D(L)
        DO 8 J=1,L
            D(J)=Z(L,J)
            Z(I,J)=ZERO
            Z(J,I)=ZERO
        8 END DO
        GOTO 9
        7 DO 10 K=1,L
            D(K)=D(K)/SCALE
            H=H+D(K)*D(K)
        10 END DO
        F=D(L)
        G=-SIGN(SQRT(H),F)
        E(I)=SCALE*G
        H=H-F*G
        D(L)=F-G
         
    !***********************************************************************
    !                                                                      *
    !     FORM A*U                                                         *
    !                                                                      *
    !***********************************************************************
         
        DO 11 J=1,L
            E(J)=ZERO
        11 END DO
        DO 12 J=1,L
            F=D(J)
            Z(J,I)=F
            G=E(J)+Z(J,J)*F
            JP1=J+1
            IF(L < JP1)GOTO 13
            DO 14 K=JP1,L
                G=G+Z(K,J)*D(K)
                E(K)=E(K)+Z(K,J)*F
            14 END DO
            13 E(J)=G
        12 END DO
         
    !***********************************************************************
    !                                                                      *
    !     FORM P                                                           *
    !                                                                      *
    !***********************************************************************
         
        F=ZERO
        DO 15 J=1,L
            E(J)=E(J)/H
            F=F+E(J)*D(J)
        15 END DO
        HH=F/(H+H)
         
    !***********************************************************************
    !                                                                      *
    !     FORM Q                                                           *
    !                                                                      *
    !***********************************************************************
         
        DO 16 J=1,L
            E(J)=E(J)-HH*D(J)
        16 END DO
         
    !***********************************************************************
    !                                                                      *
    !     FORM REDUCED A                                                   *
    !                                                                      *
    !***********************************************************************
        DO 17 J=1,L
            F=D(J)
            G=E(J)
            DO 18 K=J,L
                Z(K,J)=Z(K,J)-F*E(K)-G*D(K)
            18 END DO
            D(J)=Z(L,J)
            Z(I,J)=ZERO
        17 END DO
        9 D(I)=H
    4 END DO
     
!***********************************************************************
!                                                                      *
!     ACCUMULATION OF TRANSFORMATION MATRICES                          *
!                                                                      *
!***********************************************************************
     
    DO  I=2,N
        L=I-1
        Z(N,L)=Z(L,L)
        Z(L,L)=ONE
        H=D(I)
        IF (H == ZERO) GOTO 20
          DO  K=1,L
              D(K)=Z(K,I)/H
          END DO
        DO J=1,L
            G=ZERO
              DO  K=1,L
                G=G+Z(K,I)*Z(K,J)
              END DO
              DO  K=1,L
               Z(K,J)=Z(K,J)-G*D(K)
              END DO
        END DO
      20 DO  K=1,L
            Z(K,I)=ZERO
         END DO
    END DO
  3 DO  I=1,N
        D(I)=Z(N,I)
        Z(N,I)=ZERO
    END DO
    Z(N,N)=ONE
    E(1)=ZERO
     
    RETURN
     
    END SUBROUTINE TRED2B
     
!***********************************************************************
!***********************************************************************
     
    SUBROUTINE TQL2(NDIM,N,DIAVAL,SUBDIA,EIGVEC,IERR)
     
!***********************************************************************
!                                                                      *
!     FINDS THE EIGENVALUES AND EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL*
!     MATRIX BY THE QL METHOD                                          *
!                                                                      *
!     INPUT:  NDIM IS THE DIMENSION OF THE MATRIX                      *
!             N IS THE NUMBER OF COLUMS AND ROWS FILLED                *
!             D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX     *
!             E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX  *
!             IN ITS LAST N-1 POSITIONS, E(1) IS ARBITRARY             *
!             Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY TRED2   *
!             (IF PERFORMED) IF THE EIGENVECTORS OF THE TRIDIAGONAL    *
!             MATRIX ARE DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX   *
!                                                                      *
!     OUTPUT: D CONTAINS THE EIGENVALUES IN ASCENDING ORDER (IF AN     *
!             ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT      *
!             UNORDERED FOR INDICES 1,2,...,IERR-1)                    *
!             E HAS BEEN DESTROYED                                     *
!             Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC     *
!             TRIDIAGONAL (OR FULL) MATRIX (IF AN ERROR EXIT IS MADE Z *
!             CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED     *
!             EIGENVALUES)                                             *
!             IERR IS SET TO ZERO FOR NORMAL RETURN                    *
!             IERR IS SET TO J IF THE J-TH EIGENVALUE HAS NOT BEEN     *
!             DETERMINED AFTER 30 ITERATIONS                           *
!                                                                      *
!***********************************************************************
     
    implicit double precision (a-h,o-z),integer(i-n)
    DIMENSION DIAVAL(NDIM),SUBDIA(NDIM),EIGVEC(NDIM,NDIM)
     
!***********************************************************************
     
    ZERO=0.0
    IERR=0
    IF(N == 1)GOTO 1
    DO 2 I=2,N
        SUBDIA(I-1)=SUBDIA(I)
    2 END DO
    F=0.0D0
    TST1=0.0D0
    SUBDIA(N)=0.0D0
    DO 3 L=1,N
        J=0
        H=ABS(DIAVAL(L))+ABS(SUBDIA(L))
        IF(TST1 < H)TST1=H
         
    !***********************************************************************
    !                                                                      *
    !     LOOK FOR SMALL SUB-DIAGONAL ELEMENT                              *
    !                                                                      *
    !***********************************************************************
         
        DO 4 M=L,N
            TST2=TST1+ABS(SUBDIA(M))
            IF(TST2 == TST1)GOTO 5
        4 END DO
         
    !***********************************************************************
    !                                                                      *
    !     SUBDIA(N) IS ALWAYS ZERO, SO THERE IS NO EXIT THROUGH THE BOTTOM *
    !     OF THE LOOP                                                      *
    !                                                                      *
    !***********************************************************************
         
        5 IF(M == L)GOTO 220
        130 IF(J == 30)GOTO 1000
        J=J+1
         
    !***********************************************************************
    !                                                                      *
    !     FORM SHIFT                                                       *
    !                                                                      *
    !***********************************************************************
         
        L1=L+1
        L2=L1+1
        G=DIAVAL(L)
        P=(DIAVAL(L1)-G)/(2.0D0*SUBDIA(L))
        R=PYTHAG(P,1.0D0)
        DIAVAL(L)=SUBDIA(L)/(P+SIGN(R,P))
        DIAVAL(L1)=SUBDIA(L)*(P+SIGN(R,P))
        DL1=DIAVAL(L1)
        H=G-DIAVAL(L)
        IF(L2 > N)GOTO 145
        DO 140 I=L2,N
            DIAVAL(I)=DIAVAL(I)-H
        140 END DO
        145 F=F+H
         
    !***********************************************************************
    !                                                                      *
    !     QL TRANSFORMATION                                                *
    !                                                                      *
    !***********************************************************************
        P=DIAVAL(M)
        C=1.0D0
        C2=C
        EL1=SUBDIA(L1)
        S=0.0D0
        MML=M-L
        DO 200 II=1,MML
            C3=C2
            C2=C
            S2=S
            I=M-II
            G=C*SUBDIA(I)
            H=C*P
            R=PYTHAG(P,SUBDIA(I))
            SUBDIA(I+1)=S*R
            S=SUBDIA(I)/R
            C=P/R
            P=C*DIAVAL(I)-S*G
            DIAVAL(I+1)=H+S*(C*G+S*DIAVAL(I))
             
        !***********************************************************************
        !                                                                      *
        !     FORM VECTOR                                                      *
        !                                                                      *
        !***********************************************************************
             
            DO 180 K=1,N
                H=EIGVEC(K,I+1)
                EIGVEC(K,I+1)=S*EIGVEC(K,I)+C*H
                EIGVEC(K,I)=C*EIGVEC(K,I)-S*H
            180 END DO
        200 END DO
        P=-S*S2*C3*EL1*SUBDIA(L)/DL1
        SUBDIA(L)=S*P
        DIAVAL(L)=C*P
        TST2=TST1+ABS(SUBDIA(L))
        IF(TST2 > TST1)GOTO 130
        220 DIAVAL(L)=DIAVAL(L)+F
    3 END DO
     
!***********************************************************************
!                                                                      *
!     ORDER EIGENVALUES AND EIGENVECTORS                               *
!                                                                      *
!***********************************************************************
     
    DO 300 II=2,N
        I=II-1
        K=I
        P=DIAVAL(I)
        DO 260 J=II,N
            IF(DIAVAL(J) >= P)GOTO 260
            K=J
            P=DIAVAL(J)
        260 END DO
        IF(K == I)GOTO 300
        DIAVAL(K)=DIAVAL(I)
        DIAVAL(I)=P
        DO 280 J=1,N
            P=EIGVEC(J,I)
            EIGVEC(J,I)=EIGVEC(J,K)
            EIGVEC(J,K)=P
        280 END DO
    300 END DO
    GOTO 1
!***********************************************************************
!                                                                      *
!     SET ERROR IF NO CONVERGENCE TO AN EIGENVALUE AFTER 30 ITERATIONS *
!                                                                      *
!***********************************************************************
     
    1000 IERR = L
     
    1 RETURN
     
    END SUBROUTINE TQL2
     
!***********************************************************************
!***********************************************************************
     
    DOUBLE PRECISION FUNCTION PYTHAG(A,B)
     
!***********************************************************************
!                                                                      *
!     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW  *
!                                                                      *
!***********************************************************************
     
    DOUBLE PRECISION :: A,B
    DOUBLE PRECISION :: P,R,S,T,U
     
    P=MAX(ABS(A),ABS(B))
    IF (P == 0.0D0) GOTO 20
    R=(MIN(ABS(A),ABS(B))/P)**2
    10 CONTINUE
    T=4.0D0+R
    IF (T == 4.0D0) GOTO 20
    S=R/T
    U=1.0D0+2.0D0*S
    P=U*P
    R=(S/U)**2*R
    GOTO 10
    20 PYTHAG=P
     
    RETURN
    END FUNCTION PYTHAG 
     
!***********************************************************************
!***********************************************************************
