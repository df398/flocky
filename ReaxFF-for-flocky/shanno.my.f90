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
!     - Evaluate energies only for (unique) training set             *
!       structures (2020)                                            *
!                                                                    *
!*********************************************************************
!***********************************************************************

    subroutine minim

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension csav(nat,3),yts(na1mx3),pts(na1mx3),gis(na1mx3), &
    ais(na1mx3),bis(na1mx3),yis(na1mx3),pns(na1mx3)
    character(6) :: minmet
    character(2) :: minopt,minop(2)
    data minop /'PI','MM'/
!***********************************************************************
!                                                                      *
!     Minimisation of the molecular energy                             *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In minim'
        call timer(65)
        close (65)
    end if
    1 nit=0
    nits=0
    irset=0
    gvc=0.0400
    fctor=0.5
    eold=zero
    eest=zero
    elowest=1.0e+08
    moi=0
    nfact=0
    nfcsav=nfc
    gdicmax=zero
    minmet='SHANNO'
    iimp=0
!***********************************************************************
!                                                                      *
!     Open output files (unit 57 and unit 58)                          *
!                                                                      *
!***********************************************************************
    if (nsurp == 0) then
        open (57,file='fort.57',status='unknown',position='append')
        write (57,50)qmol
        if (nfc >= 0) write (57,60)
        if (nfc < 0) write (57,65)
        close (57)
        open (58,file='fort.58',status='unknown',position='append')
        write (58,50)qmol
        write (58,70)
        close (58)
    end if
!***********************************************************************
!                                                                      *
!     Calculate initial energy                                         *
!                                                                      *
!***********************************************************************
    call distan
    call vlist
    if (imolde == 1) call readmol
    if (iconne == 2 .AND. nbon == 0) call intcon
    if (iconne == 0 .OR. nbon == 0) call srtbon1
    if (iconne > 0 .AND. nbon > 0) call srtbon2
    call inilp
    call intcor
    call encalc
    if (nsurp == 0) call mdsav (0,qfile(nprob))
    if (nfc < 0 .AND. nmmax > 1) then        !MD-energy minimization
        call minimmd(elowest)
    end if
!***********************************************************************
!                                                                      *
!     Calculate initial RMSG                                           *
!                                                                      *
!***********************************************************************
    rmsg=0.0
    nmovh=0
    do i1=1,na
        do i2=1,3
            rmsg=rmsg+imove(i1)*d(i2,i1)*d(i2,i1)
            nmovh=nmovh+imove(i1)
        end do
    end do
    rmsg=sqrt(rmsg/float(nmovh))
!***********************************************************************
!                                                                      *
!     Output intial energies to units 57 and 58                        *
!                                                                      *
!***********************************************************************
    if (nsurp == 0) then
        open (57,file='fort.57',status='unknown',position='append')
        write (57,100)nit,estrc,gdicmax,fctor,rmsg,nfc
        close (57)
        open (58,file='fort.58',status='unknown',position='append')
        write (58,200)nit,ea,elp,eb,emol,ev+epen,ecoa,ehb,et,eco, &
        ew,ep,ech,efi,eres
        close(58)
        !open(200,file='fort.200',STATUS='unknown',position='append')
        !write (200,*) estrc
        !write (200,*) eb
        !write (200,*) ea
        !write (200,*) elp
        !write (200,*) emol
        !write (200,*) ev+epen+ecoa
        !write (200,*) ehb
        !write (200,*) et+eco
        !write (200,*) eco
        !write (200,*) ew+ep
        !write (200,*) ncha2*ech
        !close(200)
    end if
    if (nfc < 0 .OR. nmmax <= 1) goto 3
    if (qr == '1' .OR. qr == '5') goto 3       !Single point
    rmsgo=rmsg
    estrco=estrc+eres
!***********************************************************************
!                                                                      *
!     Iterative loop                                                   *
!                                                                      *
!***********************************************************************
    2 continue
    nit=nit+1
!***********************************************************************
!                                                                      *
!     Update charge distribution (every nchaud iterations)             *
!                                                                      *
!***********************************************************************
    if (mod(nit,nchaud) == 0) then
        call chargess
    end if
!***********************************************************************
!                                                                      *
!     Save old coordinates and SHANNO-variables                        *
!                                                                      *
!***********************************************************************
    do i1=1,na
        do i2=1,3
            csav(i1,i2)=c(i1,i2)
        end do
    end do
    do i1=1,na1mx3
        yts(i1)=yt(i1)
        pts(i1)=pt(i1)
        gis(i1)=gi(i1)
        ais(i1)=ai(i1)
        bis(i1)=bi(i1)
        yis(i1)=yi(i1)
        pns(i1)=pn(i1)
        eolds=eold
        eests=eest
        gvcs=gvc
    end do
!***********************************************************************
!                                                                      *
!     Call minimisation routine                                        *
!                                                                      *
!***********************************************************************
    if (nfc > 0) CALL STEDES(EOLD,EEST,moi,NFACT,GVC)
    if (nfc == 0) call shanno(eold,eest,moi,nfact,gvc,nits)
!***********************************************************************
!                                                                      *
!     Calculate new energies                                           *
!                                                                      *
!***********************************************************************
    call intcor
    call encalc
!***********************************************************************
!                                                                      *
!     Keep lowest energy                                               *
!                                                                      *
!***********************************************************************
    if (estrc+erestra < elowest) elowest=estrc
!***********************************************************************
!                                                                      *
!     Update RMSG                                                      *
!                                                                      *
!***********************************************************************
    rmsg=0.0
    nmovh=0
    do i1=1,na
        do i2=1,3
            rmsg=rmsg+imove(i1)*d(i2,i1)*d(i2,i1)
            nmovh=nmovh+imove(i1)
        end do
    end do
    rmsg=sqrt(rmsg/float(nmovh))
!***********************************************************************
!                                                                      *
!     Update step size                                                 *
!                                                                      *
!***********************************************************************
    if (estrc+eres < estrco) then
        iimp=iimp+1
        if (iimp > 4) then
            iimp=0
        !     if (nfc.gt.0.and.nfc.lt.4) nfc=nfc+1
            if (nfc > 20) nfc=nfc+25
            if (nfc < 20 .AND. nfc > 1) nfc=nfc*2
        end if
        rmsgo=rmsg
        estrco=estrc+eres
    end if
    if (estrc+eres > estrco .AND. nfc > 25) then
    !***********************************************************************
    !                                                                      *
    !     Reject move; restore coordinates and SHANNO-variables;           *
    !     lower step size                                                  *
    !                                                                      *
    !***********************************************************************
        iimp=0
        nfc=nfc/2
        if (nfc < 1) nfc=1
        do i1=1,na
            do i2=1,3
                c(i1,i2)=csav(i1,i2)
            end do
        end do
        do i1=1,na1mx3
            yt(i1)=yts(i1)
            pt(i1)=pts(i1)
            gi(i1)=gis(i1)
            ai(i1)=ais(i1)
            bi(i1)=bis(i1)
            yi(i1)=yis(i1)
            pn(i1)=pns(i1)
            eold=eolds
            eest=eests
            gvc=gvcs
        end do
        if (nit < 3 .AND. iopt == 0) nits=nits-1
        nit=nit-1
    end if
!***********************************************************************
!                                                                      *
!     Output new energies to units 57 and 58                           *
!                                                                      *
!***********************************************************************
    if (nsurp == 0) then
        open (57,file='fort.57',status='unknown',position='append')
        write (57,100)nit,estrc,gdicmax,fctor,rmsg,nfc
        close (57)
        open (58,file='fort.58',status='unknown',position='append')
        write (58,200)nit,ea,elp,eb,emol,ev,ecoa,ehb,et,eco,ew,ep,ech,efi, &
        eres
        close (58)
        !open(200,file='fort.200',STATUS='unknown',position='append')
        !write (200,*) estrc
        !write (200,*) eb
        !write (200,*) ea
        !write (200,*) elp
        !write (200,*) emol
        !write (200,*) ev+epen+ecoa
        !write (200,*) ehb
        !write (200,*) et+eco
        !write (200,*) eco
        !write (200,*) ew+ep
        !close(200)
    end if
    IF (ABS(EEST) > 1.0D8) EEST=ZERO
    if (nsurp == 0 .AND. iopt == 0) WRITE(59,300)NIT,EOLD,EEST,MINMET, &
    NCONJ,MINOPT,NFACT,FCTOR,GVC,RMSG
!***********************************************************************
!                                                                      *
!     Cell parameter optimisation (icell=2)                            *
!                                                                      *
!***********************************************************************
    if ((qr == 'F' .OR. qr == '3' .OR. qr == 'P') &
     .AND. icell == 2) call optlpar(iiter)
    if (qr == '3' .AND. icell == 2) then
        icelo2=3
        call optlpar(iiter)
    end if
!***********************************************************************
!                                                                      *
!     Save intermediate geometry every nsav2 iterations                *
!                                                                      *
!***********************************************************************
    if (mod(nit,nsav2) == 0 .AND. nsurp == 0) call mdsav(0,qfile(nprob))

    endposav=endpo
    if (mod(nit,ncontrol) == 0) call readc
    if (iruid == 1) endpo=endposav
    if (nit >= nmmax) goto 3
    if (rmsg > endpo) goto 2
    if (nit < 1) goto 2
!***********************************************************************
!                                                                      *
!     End of minimisation loop                                         *
!                                                                      *
!***********************************************************************
    3 continue
!***********************************************************************
!                                                                      *
!     Lattice parameters optimisation (icell=1)                        *
!                                                                      *
!***********************************************************************
    if ((qr == 'F' .OR. qr == '3' .OR. qr == 'P') .AND. icell == 1 &
     .AND. iexco == 0) call optlpar(iiter)
    if (qr == '3' .AND. icell == 1) then
        icelo2=3
        call optlpar(iiter)
    end if
!***********************************************************************
!                                                                      *
!     Close units 57 and 58                                            *
!                                                                      *
!***********************************************************************
    if (nsurp == 0) then
        if (i5758 == 1 .OR. iopt == 1) then
            open (57,file='fort.57',status='unknown')
            write (57,*)
            close (57)
            open (58,file='fort.58',status='unknown')
            write (58,*)
            close (58)
        end if
    end if
!***********************************************************************
!                                                                      *
!     Output final partial energies to unit 73                         *
!                                                                      *
!***********************************************************************
! df398 no need for fort.73 during minimization (fort.58 is better)
    if (nsurp /= 2) then
        open (73,file='fort.73',status='unknown',position='append')
        write (73,'(i8,1x,19(f25.10,1x))')nit,estrc,eb,ea,elp,emol,ev+epen,ecoa, & !df398 ev --> ev+epen
        ehb,et,eco,ew,ep,ech,efi,rmsg,bo(1),d(1,1)
        close (73)
    end if
    nfc=nfcsav
!***********************************************************************
!                                                                      *
!     Check whether final energy equals lowest energy                  *
!                                                                      *
!***********************************************************************
    if (elowest < estrc .AND. ikeep == 1) then
        if (iopt == 0) write (*,*)estrc,elowest
        if (iopt == 0) write (*,*)'Replaced lowest energy for',qmol
        estrc=elowest
    end if
          
    return
!***********************************************************************
!                                                                      *
!     Format part                                                      *
!                                                                      *
!***********************************************************************
    50 format (a40)
    60 format(2x,'Iter.',10x,'Epot',12x,'Max.move',6x,'Factor', &
    & 8x,'RMSG',7x,'nfc')
    65 format(2x,'Iter.',10x,'Epot',12x,'Temp(K) ',6x,'Tset (K)', &
    & 6x,'RMSG',7x,'nfc')
    70 format(2x,'Iter.',4x,'Eatom',4x,'Elopa',8x,'Ebond',6x, &
    'Emol',5x,'Eval',5x,'Ecoa',5x,'Ehb',6x,'Etor',5x,'Econj',4x, &
    'Evdw',5x,'Ecoul')
    100 format(i6,2x,f20.10,3(1x,f12.6),2x,i6)
    200 format(i6,2x,2(f14.5,1x),f14.5,1x,6(f14.5,1x),6(f14.5,1x))
    300 FORMAT(1X,I4,2(1X,F17.10),1X,A6,I6,2X,A2,I6,2X,F10.6,F14.9,F14.8, &
    F10.2,F12.2,2X,I3,1X,I1)
    end subroutine minim
!***********************************************************************
!***********************************************************************

    SUBROUTINE STEDES(EOLD,EEST,moi,NFACT,GVC)

!***********************************************************************
!                                                                      *
!     CONTROL OF MINIMISATION BY STEEPEST DESCENT METHOD               *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
!***********************************************************************
!                                                                      *
!     ZERO-RESET AND CALCULATION OF ARRAY OF FIRST DERIVATIVES BY      *
!     ROUTINE  ENCALC                                                  *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In stedes'
        call timer(65)
        close (65)
    end if
    EOLD=ESTRC+eres
!***********************************************************************
!                                                                      *
!     DIC IS THE DESCENT DIRECTION                                     *
!     IN THE FIRST STEP DIC=-D (D REPRESENTS THE GRADIENT)             *
!                                                                      *
!     THE MAXIMUM INDIVIDUAL VALUE OF DIC IS LIMITED TO 0.01           *
!     FOR PROPER SCALING IN THE LINE SEARCH (ROUTINE  FACTOR)          *
!                                                                      *
!     WHEN RETURNED FROM THE LINE SEARCH DIC REPRESENTS THE CHANGES IN *
!     THE CARTESIAN COORDINATES                                        *
!                                                                      *
!***********************************************************************
    vhulp=0.01
    if (nfc == 0) vhulp=0.001
    DO  K=1,3
        DO  L=1,na
            DIC(K,L)=-D(K,L)
       END DO
    END DO
    CALL CAGDIC(GDIC)
    IF (GDIC <= vhulp) GOTO 3
    GD10M=vhulp/GDIC
    GDIC=vhulp
    DO  K=1,3
        DO  L=1,na
            DIC(K,L)=DIC(K,L)*GD10M
        END DO
    END DO
!***********************************************************************
!                                                                      *
!     DETERMINATION OF OPTIMAL SCALING FACTOR BY ROUTINE FACTOR        *
!                                                                      *
!***********************************************************************
    3 continue
    if (nfc > 0) CALL FACTOR1(EEST,NFACT,gdic,gvc)
    if (nfc == 0) CALL FACTOR(EEST,moi,NFACT,gdic,gvc)

    END SUBROUTINE STEDES
!***********************************************************************
!***********************************************************************

    SUBROUTINE CAGDIC(GDIC)

!***********************************************************************
!                                                                      *
!     Calculation of largest move in any of the cartesian coordinates  *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DICA(NA1MX3)
    EQUIVALENCE (DIC(1,1),DICA(1))
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In cagdic'
        call timer(65)
        close (65)
    end if
    GDIC=ZERO
    DO 1 K=1,3*na
        IF (ABS(DICA(K)) > GDIC) GDIC=ABS(DICA(K))
    1 END DO

    gdicmax=gdic
    RETURN

    END SUBROUTINE CAGDIC
!***********************************************************************
!***********************************************************************

    SUBROUTINE FACTOR1(EEST,NFACT,GDIC,gvc)

!***********************************************************************
!***********************************************************************
    INCLUDE 'cbka.blk'
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In factor'
        call timer(65)
        close (65)
    end if

    rmsg=0.0
    nmovh=0
    do i1=1,na
        do i2=1,3
            rmsg=rmsg+imove(i1)*d(i2,i1)*d(i2,i1)
            nmovh=nmovh+imove(i1)
        end do
    end do
    rmsg=sqrt(rmsg/float(nmovh))
!     if (rmsg.gt.1.0) rmsg=1.0

    fctor=rmsg*float(nfc)*gdicmax/1000000.0
    if (fctor < 0.00001) fctor=0.000001
    IF (FCTOR > 1.00) FCTOR=1.00
    if (fctor < zero) fctor=-fctor
!***********************************************************************
!                                                                      *
!     New coordinates are obtained with fctor                          *
!                                                                      *
!***********************************************************************
    GVC=GDIC
    DO  K1=1,na
        DO  K2=1,3
            DIC(K2,K1)=FCTOR*DIC(K2,K1)
            if (imove(k1) == 1) C(K1,K2)=C(K1,K2)+DIC(K2,K1)
       END DO
    ENDDO

    RETURN
    END SUBROUTINE FACTOR1

!***********************************************************************
!***********************************************************************

    subroutine optlpar(iiter)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension chax(3)
!***********************************************************************
!                                                                      *
!     Lattice parameter optimization                                   *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In optlpar'
        call timer(65)
        close (65)
    end if
    iiter=0
    call distan
    call vlist
    call intcor
    call chargess
    call encalc
    esavo=estrc
    esav1=estrc
    if (nsurp /= 2 .AND. iopt == 0) then
        open (58,file='fort.58',status='unknown',position='append')
        write (58,200)nit,estrc,ea,elp,eb,emol,ev,epen,ehb,et,eco,ew,ep, &
        axis(1),axis(2),axis(3)
        write (58,'(4f9.3)')pressu,presx,presy,presz
        close (58)
    end if
    1 prest=presx+presy+presz
    chax(1)=ccpar
    chax(2)=ccpar
    chax(3)=ccpar
    do i1=1,3
        if (icelo2 == 0) axis(i1)=axis(i1)*chax(i1)
    end do
    if (icelo2 > 0 .AND. icelo2 < 4) axis(icelo2)= &
    axis(icelo2)*chax(icelo2)
    if (icelo2 == 4) then
        axis(3)=axis(3)*chax(3)
        axis(1)=axis(1)/sqrt(chax(1))
        axis(2)=axis(2)/sqrt(chax(2))
    end if
    do i1=1,na
        if (icelo2 == 0) then
            do i2=1,3
                c(i1,i2)=c(i1,i2)*chax(i2)
            end do
        end if
        if (icelo2 > 0 .AND. icelo2 < 4) c(i1,icelo2)= &
        c(i1,icelo2)*chax(icelo2)
        if (icelo2 == 4) then
            c(i1,3)=c(i1,3)*chax(3)
            c(i1,1)=c(i1,1)/sqrt(chax(1))
            c(i1,2)=c(i1,2)/sqrt(chax(2))
        end if
    end do

    call distan
    call vlist
    call intcor
    call chargess
    call encalc
    if (nsurp /= 2 .AND. iopt == 0) then
        open (58,file='fort.58',status='unknown',position='append')
        write (58,200)nit,estrc,ea,elp,eb,emol,ev,epen,ehb,et,eco,ew,ep, &
        axis(1),axis(2),axis(3)
        write (58,'(4f9.3)')pressu,presx,presy,presz
        close (58)
    end if
    esav2=estrc
    if (esav2 < esav1) then
        esav1=esav2
        iiter=iiter+1
        if (icell == 2) goto 5
        goto 1
    else
        do i1=1,3
            if (icelo2 == 0) axis(i1)=axis(i1)/chax(i1)
        end do
        if (icelo2 > 0 .AND. icelo2 < 4) axis(icelo2)= &
        axis(icelo2)/chax(icelo2)
        if (icelo2 == 4) then
            axis(3)=axis(3)/chax(3)
            axis(1)=axis(1)*sqrt(chax(1))
            axis(2)=axis(2)*sqrt(chax(2))
        end if
        do i1=1,na
            if (icelo2 == 0) then
                do i2=1,3
                    c(i1,i2)=c(i1,i2)/chax(i2)
                end do
            end if
            if (icelo2 > 0 .AND. icelo2 < 4) c(i1,icelo2)= &
            c(i1,icelo2)/chax(icelo2)
            if (icelo2 == 4) then
                c(i1,3)=c(i1,3)/chax(3)
                c(i1,1)=c(i1,1)*sqrt(chax(1))
                c(i1,2)=c(i1,2)*sqrt(chax(2))
            end if
        end do
    end if
    2 prest=presx+presy+presz
    chax(1)=ccpar
    chax(2)=ccpar
    chax(3)=ccpar
    do i1=1,3
        if (icelo2 == 0) axis(i1)=axis(i1)/chax(i1)
    end do
    if (icelo2 > 0 .AND. icelo2 < 4) axis(icelo2)= &
    axis(icelo2)/chax(icelo2)
    if (icelo2 == 4) then
        axis(3)=axis(3)/chax(3)
        axis(1)=axis(1)*sqrt(chax(1))
        axis(2)=axis(2)*sqrt(chax(2))
    end if
    do i1=1,na
        if (icelo2 == 0) then
            do i2=1,3
                c(i1,i2)=c(i1,i2)/chax(i2)
            end do
        end if
        if (icelo2 > 0 .AND. icelo2 < 4) c(i1,icelo2)= &
        c(i1,icelo2)/chax(icelo2)
        if (icelo2 == 4) then
            c(i1,3)=c(i1,3)/chax(3)
            c(i1,1)=c(i1,1)*sqrt(chax(1))
            c(i1,2)=c(i1,2)*sqrt(chax(2))
        end if
    end do
    call distan
    call vlist
    call intcor
    call chargess
    call encalc
    if (nsurp /= 2 .AND. iopt == 0) then
        open (58,file='fort.58',status='unknown',position='append')
        write (58,200)nit,estrc,ea,elp,eb,emol,ev,epen,ehb,et,eco,ew,ep, &
        axis(1),axis(2),axis(3)
        write (58,'(4f9.3)')pressu,presx,presy,presz
        close (58)
    end if
    esav3=estrc
    if (esav3 < esav1) then
        esav1=esav3
        iiter=iiter+1
        if (icell == 2) goto 5
        goto 2
    else
        do i1=1,3
            if (icelo2 == 0) axis(i1)=axis(i1)*chax(i1)
        end do
        if (icelo2 > 0 .AND. icelo2 < 4) axis(icelo2)= &
        axis(icelo2)*chax(icelo2)
        if (icelo2 == 4) then
            axis(3)=axis(3)*chax(3)
            axis(1)=axis(1)/sqrt(chax(1))
            axis(2)=axis(2)/sqrt(chax(2))
        end if
        do i1=1,na
            if (icelo2 == 0) then
                do i2=1,3
                    c(i1,i2)=c(i1,i2)*chax(i2)
                end do
            end if
            if (icelo2 > 0 .AND. icelo2 < 4) c(i1,icelo2)= &
            c(i1,icelo2)*chax(icelo2)
            if (icelo2 == 4) then
                c(i1,3)=c(i1,3)*chax(3)
                c(i1,1)=c(i1,1)/sqrt(chax(1))
                c(i1,2)=c(i1,2)/sqrt(chax(2))
            end if
        end do
    end if

    5 do i1=1,3
        axiss(i1)=axis(i1)
    end do
    return
!***********************************************************************
!                                                                      *
!     Format part                                                      *
!                                                                      *
!***********************************************************************
    200 format(i6,2x,f14.5,1x,2(f14.5,1x),f14.5,1x,6(f14.5,1x),2(f14.5,1x), &
    & 3(f14.5,1x))
    end subroutine optlpar
!***********************************************************************
!***********************************************************************

    SUBROUTINE SHANNO(EOLD,EEST,MOI,NFACT,GVC,NITS)

!***********************************************************************
!                                                                      *
!     CONTROL OF MINIMISATION BY SHANNO'S CONJUGATE GRADIENT METHOD    *
!                                                                      *
!     SEE: S.J.WATOWICH, E.S.MEYER, R.HAGSTROM, R.JOSEPS               *
!          J.COMP.CHEM.9,650-661(1988)                                 *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In shanno'
        call timer(65)
        close (65)
    end if
!***********************************************************************
!                                                                      *
!     ZERO-RESET AND CALCULATION OF ARRAY OF FIRST DERIVATIVES BY      *
!     ROUTINE  ENCALC                                                  *
!                                                                      *
!***********************************************************************
    GOTO(2,3)NITS
    1 continue
    CALL SHANN1(EOLD,EEST,MOI,NFACT,GVC)
    NITS=NITS+1
    RETURN
    2 continue
    EOLD=ESTRC+eres
    CALL SHANN2(EEST,MOI,NFACT,GVC)
    NITS=NITS+1
    RETURN
    3 continue
    EOLD=ESTRC+eres
    CALL SHANN3(EEST,MOI,NFACT,GVC)

    END SUBROUTINE SHANNO
!***********************************************************************
!***********************************************************************

    SUBROUTINE SHANN1(EOLD,EEST,MOI,NFACT,GVC)

!***********************************************************************
!                                                                      *
!     GRADIENT AND CHANGES IN THE CARTESIAN COORDINATES ARE STORED FOR *
!     USE IN SUBSEQUENT ITERATIONS                                     *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DA(NA1MX3),DICA(NA1MX3)
    EQUIVALENCE (D(1,1),DA(1)),(DIC(1,1),DICA(1))
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In shann1'
        call timer(65)
        close (65)
    end if

    CALL STEDES(EOLD,EEST,MOI,NFACT,GVC)
    nelc3=3*na
    DO 1 K=1,NELC3
        PT(K)=DICA(K)
        GI(K)=DA(K)
    1 END DO

    END SUBROUTINE SHANN1
!***********************************************************************
!***********************************************************************

    SUBROUTINE SHANN2(EEST,MOI,NFACT,GVC)

!***********************************************************************
!                                                                      *
!     IN THE SECOND ITERATION SHANNO'S ALGORITHM IS APPLIED WITHOUT THE*
!     THE SO-CALLED RESTART VECTOR                                     *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DA(NA1MX3),DICA(NA1MX3)
    EQUIVALENCE (D(1,1),DA(1)),(DIC(1,1),DICA(1))
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In shann2'
        call timer(65)
        close (65)
    end if

    C1=ZERO
    C4=ZERO
    C9=ZERO
    C10=ZERO
    nelc3=3*na
    DO 1 K=1,NELC3
        HULP1=DA(K)
        YT(K)=HULP1-GI(K)
        GI(K)=HULP1
        C4=YT(K)*YT(K)+C4
        C1=PT(K)*YT(K)+C1
        C9=PT(K)*HULP1+C9
        C10=YT(K)*HULP1+C10
    1 END DO
    DO 2 K=1,NELC3
        DICA(K)=-DA(K)-((ONE+C4/C1)*C9/C1-C10/C1)*PT(K)+C9/C1*YT(K)
    2 END DO
    CALL CAGDIC(GDIC)
!***********************************************************************
!                                                                      *
!     THE SEARCH DIRECTION IS SCALED TO A MAXIMUM VALUE OF 0.2 BEFORE  *
!     THE LINE SEARCH IN ROUTINE  FACTOR IS ENTERED                    *
!                                                                      *
!***********************************************************************
    SCALE=0.20
    CALL PREFAC(SCALE,EEST,MOI,NFACT,GVC,GDIC)
!***********************************************************************
!                                                                      *
!     CHANGES IN THE CARTESIAN COORDINATES ARE STORED FOR SUBSEQUENT   *
!     ITERATIONS                                                       *
!                                                                      *
!***********************************************************************
    DO 3 K=1,NELC3
        PN(K)=DICA(K)
    3 END DO

    END SUBROUTINE SHANN2
!***********************************************************************
!***********************************************************************

    SUBROUTINE SHANN3(EEST,MOI,NFACT,GVC)

!***********************************************************************
!                                                                      *
!     STARTING WITH THE THIRD ITERATION SHANNO'S ALGORITHM IS FULLY    *
!     APPLIED                                                          *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DA(NA1MX3),DICA(NA1MX3)
    EQUIVALENCE (D(1,1),DA(1)),(DIC(1,1),DICA(1))
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In shann3'
        call timer(65)
        close (65)
    end if

    C2=ZERO
    C3=ZERO
    C2A=ZERO
    C3A=ZERO
    C5=ZERO
    C8=ZERO
    nelc3=3*na
    DO 1 K=1,NELC3
        HULP1=DA(K)
        YI(K)=HULP1-GI(K)
        C2=PT(K)*YI(K)+C2
        C3=YT(K)*YI(K)+C3
        C2A=PT(K)*HULP1+C2A
        C3A=YT(K)*HULP1+C3A
        C5=PN(K)*YI(K)+C5
        C8=PN(K)*HULP1+C8
    1 END DO
    C6=ZERO
    C7=ZERO
    H1=C1/C4
    H2=C2/C4
    H3=TWO*C2/C1-C3/C4
    H2A=C2A/C4
    H3A=TWO*C2A/C1-C3A/C4
    C4=ZERO
    DO 2 K=1,NELC3
        C4=YI(K)*YI(K)+C4
        AI(K)=H1*YI(K)-H2*YT(K)+H3*PT(K)
        C6=AI(K)*YI(K)+C6
        C7=AI(K)*DA(K)+C7
        BI(K)=H1*DA(K)-H2A*YT(K)+H3A*PT(K)
    2 END DO
    H1=C8/C5
    H2=(ONE+C6/C5)*C8/C5-C7/C5
    DO 3 K=1,NELC3
        DICA(K)=-BI(K)+H1*AI(K)-H2*PN(K)
    3 END DO
    CALL CAGDIC(GDIC)
!***********************************************************************
!                                                                      *
!     THE SEARCH DIRECTION IS SCALED TO A MAXIMUM VALUE OF 0.1 BEFORE  *
!     THE LINE SEARCH IN ROUTINE  FACTOR IS ENTERED                    *
!                                                                      *
!***********************************************************************
    SCALE=0.20
    CALL PREFAC(SCALE,EEST,MOI,NFACT,GVC,GDIC)
!***********************************************************************
!                                                                      *
!     GRADIENT AND CHANGES IN THE CARTESIAN COORDINATES ARE STORED FOR *
!     SUBSEQUENT ITERATIONS                                            *
!                                                                      *
!     SOME REORGANIZATION IS APPLIED WITH RESPECT TO THE STORAGE OF    *
!     SOME VECTORS THAT ARE USED IN SUBSEQUENT ITERATIONS              *
!                                                                      *
!***********************************************************************
    C1=C5
    DO 4 K=1,NELC3
        PT(K)=PN(K)
        PN(K)=DICA(K)
        YT(K)=YI(K)
        GI(K)=DA(K)
    4 END DO

    END SUBROUTINE SHANN3
!***********************************************************************
!***********************************************************************

    SUBROUTINE PREFAC(SCALE,EEST,MOI,NFACT,GVC,GDIC)

!***********************************************************************
!                                                                      *
!     THE SEARCH DIRECTION IS SCALED TO A MAXIMUM VALUE OF SCALE BEFORE*
!     THE LINE SEARCH IN ROUTINE  FACTOR IS ENTERED                    *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DICA(NA1MX3)
    EQUIVALENCE (DIC(1,1),DICA(1))
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In prefac'
        call timer(65)
        close (65)
    end if
    nelc3=na*3

    FCT=ONE
    IF (GDIC > SCALE) FCT=SCALE/GDIC

    DO 1 K=1,NELC3
        DICA(K)=DICA(K)*FCT
    1 END DO

    GDIC=FCT*GDIC
!     NFC=0
    CALL FACTOR(EEST,MOI,NFACT,GDIC,GVC)

    END SUBROUTINE PREFAC
!***********************************************************************
!***********************************************************************

    SUBROUTINE FACTOR(EEST,MOI,NFACT,GDIC,GVC)

!***********************************************************************
!                                                                      *
!     PERFORMS A LINE SEARCH TO OPTIMIZE THE  MINIMIZATION ALONG A     *
!     PARTICULAR SEARCH DIRECTION (STEEPEST DESCENT, FULL MATRIX       *
!     NEWTON-RAPHSON, TRUNCATED NEWTON-RAPHSON)                        *
!                                                                      *
!     THE LINE SEARCH EMPLOYS A QUADRATIC INTERPOLATION/EXTRAPOLATION  *
!     PROCEDURE                                                        *
!                                                                      *
!***********************************************************************

    INCLUDE 'cbka.blk'

    LOGICAL :: TRYSTD,OK

    DIMENSION DUMP(3,nat)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In factor'
        call timer(65)
        close (65)
    end if

    nelc=na
!***********************************************************************
!                                                                      *
!     Next is to check whether the C-array still contains real values. *
!     On our unix computer, during error situations the values can     *
!     become NaN (when formatfree printed), without stopping execution.*
!                                                                      *
!***********************************************************************

    OK= .FALSE. 
    IF (DIC(1,1) >= ZERO) OK= .TRUE. 
    IF (DIC(1,1) <= ZERO) OK= .TRUE. 
    IF ( .NOT. OK) THEN
        WRITE(6,*)'DIC IS NOT A NUMBER (NaN); CONTINUING WITH NEXT GEOMETRY'
        RETURN
    ENDIF
!***********************************************************************
!                                                                      *
!     Stepsize for parabolic interpolation procedure is initially 1.5. *
!     Using DFLAG(25) this delta is scaled to previous delta.          *
!     When no proper factor can be found, STEDES is tried. To do this  *
!     only once, TRYSTD is used.                                       *
!     NFACT is the number of times a parabolic interpolation is        *
!     performed. Initially, MOI=2, indicating no parabolic interpolat. *
!     Old cartesian coordinates C are kept in array DUMP.              *
!                                                                      *
!***********************************************************************

    DELTA=1.5
    NFACT=0
    TRYSTD= .FALSE. 
    MOI=2                                                        ! MM
    DO K1=1,na
        DO K2=1,3
            DUMP(K2,K1)=C(K1,K2)
        END DO
    END DO

!***********************************************************************
!                                                                      *
!     NFC=0:                                                           *
!     The largest element of the search direction is compared to the   *
!     maximum change in one of the coordinates in the previous step.   *
!     The elements in the search direction are adjusted if the largest *
!     element is larger than the maximum change in the previous step.  *
!     NFC>0:                                                           *
!     According to the value of NFC either the largest change in any   *
!     cartesian coordinate is limited to NFC/1000 Angstrom or a line   *
!     search (if NFC = zero) is carried out to optimize the solution   *
!     along the search direction.                                      *
!                                                                      *
!***********************************************************************

    rmsg=0.0
    nmovh=0
    do i1=1,na
        do i2=1,3
            rmsg=rmsg+imove(i1)*d(i2,i1)*d(i2,i1)
            nmovh=nmovh+imove(i1)
        end do
    end do
    rmsg=sqrt(rmsg/float(nmovh))
!     if (rmsg.gt.1.0) rmsg=1.0

    71 continue

    IF (NFC /= 0) GOTO 151
    IF (GVC > 0.4) THEN
        GVC=0.4
    ENDIF
    GDICH=GDIC
    IF (GDIC > 10.0*GVC) GDIC=10.0*GVC
    IF (GDIC > 0.6) GDIC=0.6
    IF (GDIC /= GDICH) THEN
        FACT=GDIC/GDICH
        DO K1=1,3
            DO K2=1,NELC
                DIC(K1,K2)=DIC(K1,K2)*FACT
            END DO
        END DO
    ENDIF
    GOTO 152
    151 FCTOR=rmsg*FLOAT(NFC)/(GDIC*1.0D6)
! 151 FCTOR=FLOAT(NFC)/(GDIC*1.0D6)
    75 IF (FCTOR > 1.0) FCTOR=1.0
    EEST=1.0D20
    GOTO 62
    152 continue

!***********************************************************************
!                                                                      *
!     An optimal scaling factor FCTOR is determined by a quadratic     *
!     interpolation/extrapolation procedure in routine  FACTO2.        *
!     22:
!     61: no factor obtained. Try STEDES.                              *
!     62:
!     63:
!     64:
!     65:
!                                                                      *
!***********************************************************************

    EOLD=ESTRC
    CALL FACTO2(DUMP,DELTA,EEST,F0,F1,F2,E0,E1,E2,MOI,NFACT,*61,*62, &
    *63,*64,*65)
    GOTO 22
    61 IF (TRYSTD) GOTO 55
    DO   K1=1,NELC
        DO   K2=1,3
            C(K1,K2)=DUMP(K2,K1)
       END DO
    END DO
!     IF (IFLAG(3)) THEN
!     FCTOR=FLOAT(100)/(GDIC*1.0D3)
!     GOTO 75
!     ENDIF
    if (nsurp == 0) WRITE(59,*)' APPLYING STEDES'
    ESTRC=EOLD
    GDIC=ZERO
    FCTOR=ONE
    TRYSTD= .TRUE. 

    DO  K1=1,3
        DO  K2=1,NELC
            DIC(K1,K2)=-D(K1,K2)
            IF(ABS(D(K1,K2)) > GDIC)GDIC=ABS(D(K1,K2))
       END DO
    END DO

    IF(GDIC > 0.1)THEN
        GD10M=0.1/GDIC
        GDIC=0.1
        DO   K1=1,3
            DO  K2=1,NELC
                DIC(K1,K2)=DIC(K1,K2)*GD10M
            END DO
        END DO
    ENDIF

    IF (DELTA < 0.75) DELTA=DELTA*TWO
    GOTO 71
    55 IF (GVC < 1.0D-5) GVC=0.01
    FCTOR=ONE
    IF (GDIC > 2*GVC) FCTOR=2*GVC/GDIC
    GOTO 62

!***********************************************************************
!                                                                      *
!     THE OPTIMUM SOLUTION IS APPLIED AND NEW CARTESIAN COORDINATES    *
!     ARE CALCULATED                                                   *
!                                                                      *
!***********************************************************************
!           ****************************************
!            RETURN FROM MM PROCEDURE IN FIRST STEP
!           ****************************************
    63 continue
!     WRITE(*,100)E1,DELTA,F1
    if (nsurp == 0) WRITE(59,100)E1,DELTA,F1
    FCTOR=ONE
    MOI=2
    EEST=1.0D20
    GOTO 62
!           **********************
!            RETURN WITH NFACT=11
!           **********************
    64 IF (FCTOR == ZERO) THEN
        FCTOR=ONE
        SCALE=0.05
        IF (GDIC > SCALE) FCTOR=SCALE/GDIC
    ENDIF
    GOTO 62
!           ***************************************
!            RETURN FROM MM PROCEDURE IN LATER STEP
!           ***************************************
    65 IF (F0 == ZERO) GOTO 63
!     WRITE(*,100)E0,E1,DELTA,F0,F1
    if (nsurp == 0) WRITE(59,100)E0,E1,DELTA,F0,F1
    MOI=2
    EEST=E1
    FCTOR=F1
    GOTO 62
!           ****************************
!            CHECK ON CALCULATED FACTOR
!           ****************************
    22 IF (FCTOR == ZERO) THEN
        FCTOR=F2
        IF (E1 < E2) FCTOR=F1
        IF (FCTOR == ZERO) FCTOR=0.05
    ENDIF
    IF (MOI == 2) GVC=0.4
!           *****************************************************
!            new coordinates are calculated with obtained factor
!           *****************************************************
    62 GVC=GDIC*ABS(FCTOR)
    DO  K1=1,NELC
        DO  K2=1,3
            DIC(K2,K1)=FCTOR*DIC(K2,K1)
            if (imove(k1) == 1) C(K1,K2)=DUMP(K2,K1)+DIC(K2,K1)
       END DO
    END DO

    RETURN

    100 FORMAT(F15.7,F8.4,F9.5)

    END SUBROUTINE FACTOR

!***********************************************************************
!***********************************************************************

    SUBROUTINE FACTO2(DUMP,DELTA,EEST,F0,F1,F2,E0,E1,E2,MOI,NFACT,*,*, &
    *,*,*)

!***********************************************************************
!                                                                      *
!                                                                      *
!     PERFORMS A LINE SEARCH TO OPTIMIZE THE  MINIMIZATION ALONG A     *
!     PARTICULAR SEARCH DIRECTION ( STEEPEST DESCENT, FULL MATRIX      *
!     NEWTON-RAPHSON, TRUNCATED NEWTON-RAPHSON )                       *
!                                                                      *
!     THE LINE SEARCH EMPLOYS A QUADRATIC INTERPOLATION/EXTRAPOLATION  *
!     PROCEDURE                                                        *
!                                                                      *
!***********************************************************************

    INCLUDE 'cbka.blk'
    LOGICAL :: LTEST
    DIMENSION DUMP(3,NAT)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In facto2'
        call timer(65)
        close (65)
    end if


    LTEST= .TRUE. 

!***********************************************************************
!                                                                      *
!     AN OPTIMAL SCALING FACTOR FCTOR IS DETERMINED BY A QUADRATIC     *
!     INTERPOLATION/EXTRAPOLATION PROCEDURE                            *
!                                                                      *
!     FIRST THE ENERGY IS CALCULATED FOR TWO POINTS ALONG              *
!     THE SEARCH DIRECTION                                             *
!                                                                      *
!***********************************************************************

    MOI=1
    E0=ESTRC
    F0=ZERO
    CALL FACTWO(DUMP,E0,E1,E2,F0,F1,F2,DELTA,*93)

!***********************************************************************
!                                                                      *
!     THE THREE AVAILABLE POINTS (ONE OF THESE MAY BE THE ORIGIN)      *
!     ALONG THE SEARCH DIRECTION ARE FITTED INTO A QUADRATIC           *
!     FUNCTION FROM WHICH THE OPTIMAL SOLUTION IS ESTIMATED            *
!                                                                      *
!     THE ALGORITHM IS PROTECTED FOR AF2 BEING ZERO (STRAIGHT LINE)    *
!                                                                      *
!***********************************************************************

    1 AF1=(E2-E0)/(F2-F0)
    AF2=(E0+E2-TWO*E1)/(DELTA*DELTA)
    IF (ABS(AF2) < 1.0D-10) THEN
        FCTOR=F2
        IF (E0 < E2) FCTOR=F0
    ELSE
        FCTOR=F1-AF1/AF2
        EMB=E1-HALF*AF1/AF2*AF1
    ENDIF
    NFACT=NFACT+1
    if (nsurp == 0) WRITE(59,100)E0,E1,E2,DELTA,EMB,F0,F1,F2,FCTOR

!***********************************************************************
!                                                                      *
!     THE ESTIMATE OF THE OPTIMUM SOLUTION IS ACCEPTED WHEN THIS       *
!     SOLUTION IS BRACKETED BETWEEN TWO OF THE THREE CALCULATED        *
!     POINTS AND WHEN FCTOR IS AT LEAST ONE TENTH OF DELTA             *
!                                                                      *
!***********************************************************************

    CALL FACTO3(AF2,DUMP,DELTA,E0,E1,E2,F0,F1,F2,EEST,EMB,NFACT, &
    LTEST,*1,*91,*92,*94,*95)

    RETURN

    91 continue
    if (nsurp == 0) WRITE(59,*)' FACTO2: RETURN 1'
    RETURN 1                            ! RETURN WITHOUT A NEW FACTOR
    92 continue
    if (nsurp == 0) WRITE(59,*)' FACTO2: RETURN 2'
    RETURN 2                           ! RETURN WITH BEST POINT FOUND
    93 continue
    if (nsurp == 0) WRITE(59,*)' FACTO2: RETURN 3'
    RETURN 3                               ! RETURN WITH MM PROCEDURE
    94 continue
    if (nsurp == 0) WRITE(59,*)' FACTO2: RETURN 4'
    RETURN 4                                   ! RETURN WITH NFACT=11
    95 continue
    if (nsurp == 0) WRITE(59,*)' FACTO2: RETURN 5'
    RETURN 5                 ! RETURN WITH MM PROCEDURE IN LATER STEP

!***********************************************************************
!                                                                      *
!     FORMAT PART                                                      *
!                                                                      *
!***********************************************************************


    100 FORMAT(3F15.7,F8.4,F15.7,3F8.4,F12.4)

    END SUBROUTINE FACTO2

!***********************************************************************
!***********************************************************************

    SUBROUTINE FACTWO(DUMP,E0,E1,E2,F0,F1,F2,DELTA,*)

!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DUMP(3,nat)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In factwo'
        call timer(65)
        close (65)
    end if

    if (nsurp == 0) WRITE(59,*)' FACTWO'

    F1=F0+DELTA
    CALL CFDIC(DUMP,F1,E1)

!***********************************************************************
!                                                                      *
!     THE ENERGY CALCULATED FOR THIS PARTICULAR POINT ON THE SEARCH    *
!     DIRECTION IS COMPARED WITH THE PREVIOUS ENERGY TO SEE IF THE     *
!     MINIMIZATION IS ALMOST DONE; IN THAT CASE NO LINE SEARCH IS      *
!     CARRIED OUT BUT THE SOLUTION IS APPLIED AS SUCH                  *
!                                                                      *
!***********************************************************************

    IF (ABS(E0-E1) < 1.0D-10) RETURN 1

!***********************************************************************
!                                                                      *
!     A SECOND POINT IS CALCULATED ALONG THE POSITIVE DIRECTION        *
!     FURTHER ALONG THIS DIRECTION OR HALF WAY                         *
!                                                                      *
!***********************************************************************

    IF (E1 > E0) GOTO 1
    F2=F0+TWO*DELTA
    CALL CFDIC(DUMP,F2,E2)

    RETURN

    1 F2=F1
    E2=E1
    DELTA=DELTA/TWO
    F1=F0+DELTA
    CALL CFDIC(DUMP,F1,E1)

    RETURN

    END SUBROUTINE FACTWO

!***********************************************************************
!***********************************************************************

    SUBROUTINE FACTO3(AF2,DUMP,DELTA,E0,E1,E2,F0,F1,F2,EEST,EMB, &
    NFACT,LTEST,*,*,*,*,*)

!***********************************************************************
    INCLUDE 'cbka.blk'
    LOGICAL :: LTEST
    DIMENSION DUMP(3,nat)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In facto3'
        call timer(65)
        close (65)
    end if

    CALL KEEPLO(E0,E1,E2,F0,F1,F2,EL,FL)
    IF (NFACT == 11) GOTO 5

    IF (FCTOR < ZERO) THEN
        if (nsurp == 0) WRITE(59,*)' NEGATIVE  FACTOR'
        IF (EL < E0) THEN      ! HILLPARABOL
            if (nsurp == 0) WRITE(59,*)' HILLPARABOL'
            IF (E2 < E1) GOTO 4
            CALL KEEPLO(E0,E1,E2,F0,F1,F2,EEST,FCTOR)
            if (nsurp == 0) WRITE(59,*)' LOWEST POINT=E1 ACCEPTED'
            RETURN
        ENDIF
    !     IF ((.NOT.IFLAG(4)).OR.DFLAG(33)) THEN
        if (nsurp == 0) WRITE(59,*)'IS NOT ACCEPTED'
        GOTO 1
    !     ELSE
    !     CONTINUE
    !     ENDIF
    ENDIF

    IF (AF2 < ZERO) THEN   ! hillparabool
        if (nsurp == 0) WRITE(59,*)' HILLPARABOL'
        IF (E2 < E0) GOTO 4
        GOTO 13
    ENDIF

    IF (E1 < E0 .AND. E1 <= E2) GOTO 2
    IF (E1 >= E0) GOTO 3
    IF (E2 < E1) GOTO 4
    if (nsurp == 0) WRITE(59,*)' MAG HIER NIET KOMEN, FACTO3'
!     WRITE(*,*)' MAG HIER NIET KOMEN, FACTO3'
    if (nsurp == 0) WRITE(59,'(6F18.7)')E0,E1,E2,F0,F1,F2
    WRITE(*,'(6F18.7)')E0,E1,E2,F0,F1,F2
    WRITE(6,*)'FACTO3 error, continuing with next geometry'
    RETURN

!***********************************************************************
!                                                                      *
!     NEGATIVE  FACTOR IS NOT ACCEPTED                                 *
!                                                                      *
!***********************************************************************

    1 IF (F0 /= ZERO) GOTO 5
    IF (DELTA < 0.001) GOTO 92
    GOTO 13

!***********************************************************************
!                                                                      *
!     E1 IS THE POINT WITH LOWEST ENERGY                               *
!     USING DFLAG(22) A CHECK IS DONE ON A POSSIBLE POTENTIAL WALL. IN *
!     THIS CASE THE ESTIMATED POINT IS HALFWAY BETWEEN POINT 0 AND 1   *
!     OR HALFWAY BETWEEN POINT 1 AND 2.                                *
!                                                                      *
!***********************************************************************

    2 IF(FCTOR < (DELTA/10.0) .AND. F0 == ZERO .AND. DELTA > 0.001)GOTO 13
    GOTO 14                          ! EMB ACCEPTEREN

!***********************************************************************
!                                                                      *
!     E1 HAS A HIGHER ENERGY THAN E0, --> INTERPOLATION OR ACCEPT      *
!                                                                      *
!***********************************************************************

    3 continue
!   3 IF (IFLAG(4).AND..NOT.DFLAG(33)) THEN
!     IF (FCTOR.LT.F0-DELTA) GOTO 11   ! PUNT VERDER ZOEKEN
!     GOTO 13                          ! EMB ACCEPTEREN
!     ENDIF
    IF (DELTA < 0.0001 .AND. F0 == ZERO) GOTO 25
    IF (DELTA < 0.01 .AND. F0 /= ZERO) GOTO 5
    GOTO 13

!***********************************************************************
!                                                                      *
!     E2 HAS A LOWER ENERGY THAN E1, --> EXTRAPOLATION OR ACCEPT       *
!                                                                      *
!***********************************************************************

    4 continue
!   4 IF (DFLAG(22)) THEN               ! IVM POTENTIAAL MUUR
!     T1=F1+(HALF+DREAL(2))*DELTA
!     T2=F1+(HALF-DREAL(2))*DELTA
!     IF (FCTOR.LT.T1.AND.FCTOR.GT.T2) GOTO 11      ! PUNT VERDER ZOEKEN
!     ENDIF
    IF (EMB > E2) GOTO 5            ! E2 ACCEPTEREN
    IF (FCTOR > F2+DELTA) GOTO 11   ! PUNT VERDER ZOEKEN
    GOTO 15                          ! EMB ACCEPTEREN

!***********************************************************************
!                                                                      *
!     ACCEPT E1                                                        *
!                                                                      *
!***********************************************************************

    25 eest=e1
    fctor=f1
    if (nsurp == 0) WRITE(59,*)' accept e1'
!     WRITE(*,*)' accept e1'
!     WRITE(*,100)E0,E1,E2,DELTA,EMB,F0,F1,F2,FCTOR

    RETURN

!***********************************************************************
!                                                                      *
!     ACCEPT LOWEST POINT                                              *
!                                                                      *
!***********************************************************************

    5 CALL KEEPLO(E0,E1,E2,F0,F1,F2,EEST,FCTOR)
    if (nsurp == 0) WRITE(59,*)' keeplo'
    IF (NFACT == 11) GOTO 94

    RETURN

!***********************************************************************
!                                                                      *
!     EXTRAPOLATION                                                    *
!                                                                      *
!***********************************************************************
    11 continue
    if (nsurp == 0) WRITE(59,*)' extrapolation'
    IF (LTEST) THEN
        DELTA=(FCTOR-F2)*0.8
        IF (FCTOR < ZERO) DELTA=(F0-FCTOR)*0.8
        DELTA=ABS(DELTA)
        IF (DELTA > 10.0) DELTA=SQRT(DELTA)*SQRT(10.0)
        IF (FCTOR < ZERO) DELTA=-DELTA
        E0=E1
        F0=F1
        GOTO 21
    ENDIF
    CALL FACEXT(E0,E1,E2,E3,F0,F1,F2,DELTA,DUMP,*22)
    GOTO 91
    22 EEST=E2
    FCTOR=F2

    RETURN

!***********************************************************************
!                                                                      *
!     INTERPOLATION (not called)                                       *
!                                                                      *
!***********************************************************************

!  12 continue
!     if (nsurp.eq.0) WRITE(59,*)' interpolation'
!     CALL FACINT(E0,E1,E2,E3,F0,F1,F2,DELTA,FCTOR,DUMP)
!     GOTO 91

!***********************************************************************
!                                                                      *
!     DIMINSH STEP AND CALCULATE AGAIN TWO POINTS                      *
!                                                                      *
!***********************************************************************

    13 continue
    if (nsurp == 0) WRITE(59,*)' delta=delta/three'
    DELTA=DELTA/3.0
    IF (AF2 > ZERO) LTEST= .FALSE. 

!***********************************************************************
!                                                                      *
!     TWO NEW POINTS ARE CALCULATED                                    *
!                                                                      *
!***********************************************************************
    21 continue
    if (nsurp == 0) WRITE(59,*)' two new points'
    CALL FACTWO(DUMP,E0,E1,E2,F0,F1,F2,DELTA,*95)
    GOTO 91

!***********************************************************************
!                                                                      *
!     ACCEPT ESTIMATED POINT                                           *
!     WITH DFLAG(29) THE BEST POINT IS SELECTED AFTER CHECKING ESTRC   *
!                                                                      *
!***********************************************************************

    14 continue
    if (nsurp == 0) WRITE(59,*)' accept estimated point'
!     IF (.NOT.DFLAG(29)) GOTO 15
!     CALL CFDIC(DUMP,FCTOR,EF)
!     if (nsurp.eq.0) WRITE(59,*)' DFlag(29) check:',EL,EF
!     IF (EL.LT.EF) THEN
!     if (nsurp.eq.0) WRITE(59,*)
!    $' ESTIMATED POINT REJECTED, BEST POINT TAKEN'
!     FCTOR=FL
!     EMB=EL
!     ENDIF
!     EEST=EMB
!     GOTO 93

    15 EEST=EMB

    RETURN

    91 continue
    if (nsurp == 0) WRITE(59,*)' FACTO3: RETURN 1'
    RETURN 1  ! RETURN AFTER SUCCESSFULL INTER- OR EXTRAPOLATION
    92 continue
    if (nsurp == 0) WRITE(59,*)' FACTO3: RETURN 2'
    RETURN 2  ! RETURN WITHOUT A NEW FACTOR
    94 continue
    if (nsurp == 0) WRITE(59,*)' FACTO3: RETURN 4'
    RETURN 4  ! RETURN WITH NFACT=11
    95 continue
    if (nsurp == 0) WRITE(59,*)' FACTO3: RETURN 5'
    RETURN 5  ! RETURN WITH MM

    100 FORMAT(3F15.7,F8.4,F15.7,3F8.4,F12.4)

    END SUBROUTINE FACTO3

!***********************************************************************
!***********************************************************************

    SUBROUTINE CFDIC(DUMP,FFACTOR,EF)

!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DUMP(3,nat),dsav(3,nat)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In cfdic'
        call timer(65)
        close (65)
    end if

    nelc=na
    DO K=1,NELC
        DO L=1,3
            dsav(l,k)=d(l,k)
            if (imove(k) == 1) C(K,L)=DUMP(L,K)+FFACTOR*DIC(L,K)
        ENDDO
    ENDDO
    call intcor
    CALL ENCALC

    do k=1,nelc
        do l=1,3
            d(l,k)=dsav(l,k)
        end do
    end do

    EF=ESTRC+eres

    RETURN

    END SUBROUTINE CFDIC
!***********************************************************************
!***********************************************************************

    SUBROUTINE KEEPLO(E0,E1,E2,F0,F1,F2,EL,FL)

!***********************************************************************
    INCLUDE 'cbka.blk'
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In keeplo'
        call timer(65)
        close (65)
    end if

    EL=E1    ! het beste punt wordt bewaard
    FL=F1
    IF (E1 < E0 .AND. E1 < E2) RETURN
    EL=E2
    FL=F2
    IF (E2 < E0) RETURN
    EL=E0
    FL=F0

    END SUBROUTINE KEEPLO

!***********************************************************************
!***********************************************************************

    SUBROUTINE FACEXT(E0,E1,E2,E3,F0,F1,F2,DELTA,DUMP,*)

!***********************************************************************
!                                                                      *
!     ANOTHER POINT IS CALCULATED FURTHER ALONG THE POSITIVE DIRECTION *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DUMP(3,nat)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In facext'
        call timer(65)
        close (65)
    end if

    if (nsurp == 0) WRITE(59,*)' FACEXT'

    IF (E2 < E1) GOTO 1
    GOTO 2

!***********************************************************************
!                                                                      *
!     POSITIVE  FACTOR                                                 *
!                                                                      *
!***********************************************************************

    1 F3=F2+DELTA
    CALL CFDIC(DUMP,F3,E3)
    IF (E3 >= E0) RETURN 1
    E0=E1
    E1=E2
    E2=E3
    F0=F1
    F1=F2
    F2=F3

    RETURN

!***********************************************************************
!                                                                      *
!     POSITIVE  FACTOR                                                 *
!                                                                      *
!***********************************************************************
    2 F3=F0-DELTA
    CALL CFDIC(DUMP,F3,E3)
    IF (E3 >= E2) RETURN 1
    E2=E1
    E1=E0
    E0=E3
    F2=F1
    F1=F0
    F0=F3

    RETURN

    END SUBROUTINE FACEXT

!***********************************************************************
!***********************************************************************

    SUBROUTINE FACINT(E0,E1,E2,E3,F0,F1,F2,DELTA,DUMP)

!***********************************************************************
!                                                                      *
!     ANOTHER POINT IS CALCULATED HALFWAY ALONG THE POSITIVE DIRECTION *
!                                                                      *
!***********************************************************************
    INCLUDE 'cbka.blk'
    DIMENSION DUMP(3,nat)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In facint'
        call timer(65)
        close (65)
    end if

    if (nsurp == 0) WRITE(59,*)' FACINT'

    DELTA=DELTA/TWO
    F3=F0+DELTA
    IF (FCTOR > F1) F3=F1+DELTA
    CALL CFDIC(DUMP,F3,E3)

    IF (F3 > F1) THEN
        E0=E1
        F0=F1
    ELSE
        E2=E1
        F2=F1
    ENDIF

    E1=E3
    F1=F3

    RETURN

    END SUBROUTINE FACINT

!***********************************************************************
!***********************************************************************

    subroutine minimmd(elowest)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!***********************************************************************
!                                                                      *
!     Minimization using MD-routines                                   *
!                                                                      *
!***********************************************************************
    tsetsav=tset
    5 rmsg=0.0
    nmovh=0
    do i1=1,na
        do i2=1,3
            rmsg=rmsg+imove(i1)*d(i2,i1)*d(i2,i1)
            nmovh=nmovh+imove(i1)
        end do
    end do
    rmsg=sqrt(rmsg/float(nmovh))
    if (nsurp == 0) then
        open (57,file='fort.57',status='unknown',position='append')
        write (57,100)nit,estrc,tempmd,tset,rmsg,nfc
        close (57)
        open (58,file='fort.58',status='unknown',position='append')
        write (58,200)nit,ea,elp,eb,emol,ev,ecoa,ehb,et,eco,ew,ep,ech,efi, &
        eres
        close (58)
    end if

    if (rmsg < endpo) then
        if (ipropt == 0) then
            tset=tsetsav
            return
        end if
    end if


    call verlet1
    call intcor
    if (ipropt > 0) icpres=1
    call encalc
    call verlet2
    call tdamp(scasum)
!     write (64,*)nmethod,ipropt,icpres,pressu,0.001*presx,
!    $0.001*presy,0.001*presz,pset,vprestax,vprestay,vprestaz,vpresopt
    if (ipropt == 1) call pdamp
    nit=nit+1
    if (estrc < elowest) elowest=estrc
    if (tset+tincr > 0.01) tset=tset+tincr
    endposav=endpo

    if (rmsg < endpo) then
        if (ipropt == 2) then
            tset=tsetav
            return
        end if
        testp=abs(0.001*presx-vprestax)+abs(0.001*presy-vprestay) &
        +abs(0.001*presz-vprestaz)
    !     write (64,*)pressu,pset,vpresopt,testp
        if (ipropt == 1 .AND. testp < vpresopt) then
            tset=tsetav
            return
        end if
    end if

    if (mod(nit,ncontrol) == 0) call readc
    if (iruid == 1) endpo=endposav
    if (nit < nmmax) goto 5
    10 continue
    tset=tsetsav
    return
    100 format(i6,2x,f20.10,3(1x,f12.6),2x,i6)
    200 format(i6,2x,2(f14.5,1x),f14.5,1x,6(f14.5,1x),4(f14.5,1x))
    end subroutine minimmd
!***********************************************************************
!***********************************************************************
