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
    subroutine reac
    include 'cbka.blk'
    include 'cbkm.blk'
    include 'opt.blk'
    dimension clast(nat,3)
    dimension totf1(3),totf2(3),totf3(3),avef1(3),avef2(3),avef3(3)
    dimension nmcb(nmolmax),nmce(nmolmax), &
    disres(mrestra)
    dimension vibreax(navib*3),vibqc(navib*3),errmatch(3*navib)
    dimension imatch (3*navib)
    dimension isumattreg(mtreg)
    real :: tarray(2)
    character(80) :: qromb
    character(40) :: qfreqfile
    character(25) :: qfileh
    character(2) :: qrom
    character(4) :: hdr
    character qmou
    integer :: icntrl(20)
    logical :: lperiod,lmolxtl,lcanon,ldefcel,lprtthrm,lnose, &
    lnpecan,ltmpdamp
    !df398
    junique = 0
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In reac'
        call timer(65)
        close (65)
    end if
    convmd=4.184d26
    pi=3.14159265358979323846D0 
    avognr=6.022140857D23 
    rdndgr=180.0d0/pi
    dgrrdn=1.0d0/rdndgr
    rgasc=8.3144598d0 !df398 molar gas constant to 7 decimal places
    caljou=4.184d0
    xjouca=1.0d0/caljou
    zero=0.0d0
    ech=zero
    one=1.0d0
    two=2.0d0
    three=3.0d0
    half=one/two
    nzero=0
    none=1
    ntwo=2
    nthree=3
    irun=1
    invt=0
    ndata2=0
    iheatf=0
    nradcount=0
    itemp=1
    xinh=zero
    ifieldx=0
    ifieldy=0
    ifieldz=0
    mdstep=0
    kx=0
    ky=0
    kz=0
    nit=0
    nbon=0
    angle(1)=90.0d0
    angle(2)=90.0d0
    angle(3)=90.0d0
    axiss(1)=zero
    axiss(2)=zero
    axiss(3)=zero
    do i1=1,nat
        id(i1,1)=0
        id(i1,2)=0
        id(i1,3)=0
    end do
    call readtrans    !Apply translations from previous simulations to id-array

    do i1=1,nmolset
        icgeopt(i1)=0
        ifreqset(i1)=0
    end do

    do i2=1,neem      !Zero EEM-matrix
        do i3=1,neem
            xmortr(i3,i2)=zero
        end do
    end do

    icgeo=0
    sumhe=zero
    call readc
    call ffinpt
    call tap7th
    n=0

    if (nmethod == 3 .AND. nmm == 0 .AND. nhop2 < maxstp) then
        write (*,*)'NVE-ensemble: removal of translational and '
        write (*,*)'rotational energy switched off'
    end if

!     call etime(tarray)
!     totime=tarray(1)
!     prvtme=totime
    call readtreg(isumattreg)
    call readvreg
    call readereg
    call readaddmol
    call readpiston
    nprob=0
    nuge=1
    nubgf=1
    imodfile=0
    5 continue
    ioldchg=0
    nprob=nprob+1
    if (nprob > nmolset) then
        write (*,*) nprob,nmolset
        write (*,*) 'Too many molecules in training set'
        stop 'Too many molecules in training set'
    end if
    naold=na
    na=0
    nmollset=nprob-1
    nrestra=0
    nrestrav=0
    nrestrat=0
    nrestram=0
    nrestras=0
    neqdis=0
    nrestraset(nprob)=0
    nrestratset(nprob)=0
    nrestravset(nprob)=0
    tset=tsetor
    tm11=axis(1)
    tm21=zero
    tm31=zero
    tm22=axis(2)
    tm32=zero
    tm33=axis(3)
    qruid='NORMAL RUN'
    Hug_chi = 0.d0
    Hug_chi_dot = 0.d0
    do i1=1,nat
        imove(i1)=1
    end do

    open (14,file='models.in',status='old',err=10)
    imodfile=1
    iline=0
    6 continue
    if (nsurp /= 2) read (14,'(a80)',end=900)qromb
    if (nsurp >= 2) read (14,'(a80)',end=50)qromb
    qstrana1(1:80)=qromb
    if (qromb(1:1) /= '#') iline=iline+1
    if (iline < nprob) goto 6
    7 continue
    istart=1
    call stranal(istart,iend,vout,iout,1)
    qfile(nprob)=qromb(istart:iend-1)
    istart=iend
    call stranal(istart,iend,vout,iout,1)
    qkeyw(nprob)=qromb(istart:iend-1)
    qfileh=qfile(nprob)
    if (qfileh(1:1) == '#') goto 7
    close (14)
    10 continue

    qfileh=qfile(nprob)
    if (qfileh(1:3) == 'END' .AND. nsurp /= 2) goto 900
    if (qfileh(1:3) == 'END' .AND. nsurp >= 2) goto 50

    if (ngeofor == 0) then
        call readdelphi (qfileh,iend,naold)
        namov=na
        if (imodfile == 1) then
            if (ingeo == 0) then ! do not overwrite old geometry files
                qromb=qfile(nprob)
                qstrana1(1:80)=qromb
                istart=1
                call stranal(istart,iend2,vout,iout,1)
                qfile(nprob)=qromb(istart:iend2-1)//".out"
            end if
        else
            istart=1
            qstrana1(1:60)=qmol
            call stranal(istart,iend2,vout,iout,1)
            if (iend2 > 30) iend2=29
            qfile(nprob)=qmol(istart:iend2-1)
            if (iend2 > 30) iend2=30
            qkeyw(nprob)=qmol(istart:iend2-1)
        end if
        if (iend == 1 .AND. nsurp /= 2) goto 900
        if (iend == 1 .AND. nsurp >= 2) goto 50
    end if

    if (ngeofor == 1) then
        !df398
        call gettrainstructs2
        call gettrainstructs
        call readbgf(qfile(nprob),iend,naold)

        if (imodfile == 1) then
            if (ingeo == 0) then ! do not overwrite old geometry files
                qromb=qfile(nprob)
                qstrana1(1:80)=qromb
                istart=1
                call stranal(istart,iend2,vout,iout,1)
                qfile(nprob)=qromb(istart:iend2-1)//".out"
            end if
        else
            istart=1
            qstrana1(1:60)=qmol
            call stranal(istart,iend2,vout,iout,1)
            if (iend2 > 30) iend2=30
            qfile(nprob)=qmol(istart:iend2-1)
            if (iend2 > 30) iend2=30
            qkeyw(nprob)=qmol(istart:iend2-1)
        end if
        if (iend == 1 .AND. nsurp /= 2) goto 900
        if (iend == 1 .AND. nsurp >= 2) goto 50
    end if


    if (ngeofor == 2) then
        if (nsurp >= 2) stop 'Xmol only with 1-geometry read-in'
    !*********************************************************************
    !                                                                    *
    !     Read in free format (xmol) geometry                            *
    !                                                                    *
    !*********************************************************************
        qr='1'
        read (3,'(i6)',end=900)na
        namov=na
        read (3,'(a60)',end=900)qmol
        do i1=1,na
            read (3,'(a80)')qromb
            ifirstchar=80
            do i2=1,80
                if (qromb(i2:i2) /= ' ' .AND. i2 < ifirstchar) ifirstchar=i2
            end do
            read (qromb(ifirstchar:80),'(a2)')qa(i1)
            read (qromb(ifirstchar+2:80),*)c(i1,1),c(i1,2),c(i1,3)
            qlabel(i1)=qa(i1)
            qresi1(i1)='   '
            qresi2(i1)=' '
            qresi3(i1)='     '
            qffty(i1)='     '
        end do
        ibity=1
        axiss(1)=-1.0d0
    end if

          
    if (ngeofor == 3) then
        if (nsurp >= 2) stop 'CC1 only with 1-geometry read-in'
    !*********************************************************************
    !                                                                    *
    !     Read in ChemDraw CC1-file                                      *
    !                                                                    *
    !*********************************************************************
        qr='1'
        read (3,*)na
        namov=na
        read (3,'(a60)')qmol
        do i1=1,na
            read (3,'(2x,a2,5x,3f12.6)')qa(i1),c(i1,1),c(i1,2),c(i1,3)
        end do
    end if

    if (ngeofor == 4) then
        if (nsurp >= 2) stop 'pdb only with 1-geometry read-in'
    !*********************************************************************
    !                                                                    *
    !     Read in .pdb-format                                            *
    !                                                                    *
    !*********************************************************************
        qr='C'
        call readpdb(iendf)
        namov=na
        ibity=1
        axiss(1)=-1.0d0
        qfile(nprob)=qmol
        if (iendf == 1) goto 900
    end if
!*********************************************************************
!                                                                    *
!     Check consistency temperature programs in tregime.in           *
!                                                                    *
!*********************************************************************
    if (ntrc > 0) then
        do i1=1,ntrc
            if (isumattreg(i1) /= na) then
                write (*,*)'Inconsistency in temperature regime nr.',i1
                write (*,*)'Number of atoms defined in tregime.in:', &
                isumattreg(i1)
                write (*,*)'Number of atoms in system:',na
                stop 'Inconsistency in tregime.in'
            end if
        end do
    end if

    if (nsurp == 0) write (*,2000)qmol
          
    naset(nprob)=na
    enmolset(nprob)=enmol
    formolset(nprob)=formol

    qrset(nprob)=qr
    qmolset(nprob)=qmol
    qremset(nprob)=qremark(1)

    kxset(nprob)=kx !! df398 previously was set to kz (probably a bug)
    kyset(nprob)=ky
    kzset(nprob)=kz

    nrestraset(nprob)=nrestra
    iredoset(nprob)=iredo
    iexcoset(nprob)=iexco
    ncellset(nprob)=ncellopt
    iruidset(nprob)=iruid
    icellset(nprob)=icell
    icelo2lset(nprob)=icelo2
    ibityset(nprob)=ibity
    nfcset(nprob)=nfc
    nmmaxset(nprob)=nmmax
    endposet(nprob)=endpo
    vvolset(nprob)=vvol
    do i1=1,nrestra
        irstraset(nprob,i1,1)=irstra(i1,1)
        irstraset(nprob,i1,2)=irstra(i1,2)
        rrstraset(nprob,i1)=rrstra(i1)
        vkrstrset(nprob,i1)=vkrstr(i1)
        vkrst2set(nprob,i1)=vkrst2(i1)
        rrchaset(nprob,i1)=rrcha(i1)
    end do

    nrestratset(nprob)=nrestrat
    do i1=1,nrestrat
        irstratset(nprob,i1,1)=irstrat(i1,1)
        irstratset(nprob,i1,2)=irstrat(i1,2)
        irstratset(nprob,i1,3)=irstrat(i1,3)
        irstratset(nprob,i1,4)=irstrat(i1,4)
        trstraset(nprob,i1)=trstra(i1)
        vkrtset(nprob,i1)=vkrt(i1)
        vkr2tset(nprob,i1)=vkr2t(i1)
    end do

    nrestravset(nprob)=nrestrav
    do i1=1,nrestrav
        irstravset(nprob,i1,1)=irstrav(i1,1)
        irstravset(nprob,i1,2)=irstrav(i1,2)
        irstravset(nprob,i1,3)=irstrav(i1,3)
        vrstraset(nprob,i1)=vrstra(i1)
        vkrvset(nprob,i1)=vkrv(i1)
        vkr2vset(nprob,i1)=vkr2v(i1)
    end do

    do i1=1,3
        axisset(nprob,i1)=axiss(i1)
        anglesset(nprob,i1)=angles(i1)
    end do

    do i2=1,na
        do i1=1,3
            cset(nprob,i2,i1)=c(i2,i1)
        end do
        qaset(nprob,i2)=qa(i2)
        chaset(nprob,i2)=chgbgf(i2)
    end do


    if (nsurp >= 2) goto 5
!*********************************************************************
!                                                                    *
!     Read in morphing geometry                                      *
!                                                                    *
!*********************************************************************
    imorph=0
    if (qr == 'W') then
        imorph=1
        read (39,*)
        do i1=1,na
            read (39,1100)ir,qrom,(cmo(i1,i2),i2=1,3)
            read (40,1110)ihulp,vmo1(i1),vmo2(i1)
        end do
    end if
!*********************************************************************
!                                                                    *
!     Read in fixed chargess from unit 26 (optional)                  *
!                                                                    *
!*********************************************************************
    if (ncha == 5) then   !fixed chargess
        read (26,1220)ntymol
        read (26,1220)(nmcb(i1),nmce(i1),i1=1,ntymol)

        nhulp=0
        do i1=1,ntymol
            read (26,*)
            read (26,1220)natmty

            do i2=1,natmty
                read (26,1240)irom,charge
                ch(nhulp+i2)=charge
                do i3=1,nmce(i1)-nmcb(i1)
                    ch(nhulp+i3*natmty+i2)=charge
                end do
            end do
            nhulp=nhulp+natmty*(nmce(i1)-nmcb(i1)+1)
        end do
        if (nhulp /= na) then
            write (*,*)nhulp,na
            stop 'INCONSISTENCY IN CHARGE DATA (UNIT 26)'
        end if
        totch=0.0d0
        do i1=1,na
            totch=totch+ch(i1)
        end do
        !write (*,*)'Number of atoms',na,' Total charge',totch
        !df398 added rewind(26) so that fixed charges option could be used
        !with more than one structure in geo file. Otherwise when new structures
        !are read from geo, the position of fort.26 is at END OF FILE and cannot
        !proceed to read the charges. This is not needed in REAXFFGMIN
        ! or NEWREAXFFOPTIM since we would not normally use more than one
        ! structure in fort.3 file. Added for consistency with the standalone
        ! reaxff code when running calculations on several structures in geo file.
        rewind(26)
    end if
    50 if (nsurp >= 2) then
        nmmax=nmmaxold
        ioldchg=0
        nfc=nfcold
        endpo=endpoold
        icell=icellold
        icelo2=icelo2old
        iruid=0
        nprob=nuge
        na=naset(nuge)
        iredo=iredoset(nuge)
        iexco=iexcoset(nuge)
        ncellopt=ncellset(nuge)
        iruid=iruidset(nuge)
        icell=icellset(nuge)
        icelo2=icelo2lset(nuge)
        ncellopt=ncellset(nuge)
        ibity=ibityset(nuge)
        nfc=nfcset(nuge)
        nmmax=nmmaxset(nuge)
        endpo=endposet(nuge)
        qr=qrset(nuge)
        qmol=qmolset(nuge)
        nrestra=nrestraset(nuge)
        nrestrat=nrestratset(nuge)
        nrestrav=nrestravset(nuge)
        kx=kxset(nuge)
        ky=kyset(nuge)
        kz=kzset(nuge)
        do i1=1,nrestra
            rrstra(i1)=rrstraset(nuge,i1)
            vkrstr(i1)=vkrstrset(nuge,i1)
            vkrst2(i1)=vkrst2set(nuge,i1)
            irstra(i1,1)=irstraset(nuge,i1,1)
            irstra(i1,2)=irstraset(nuge,i1,2)
            rrcha(i1)=rrchaset(nuge,i1)
        end do
        do i1=1,nrestrat
            irstrat(i1,1)=irstratset(nuge,i1,1)
            irstrat(i1,2)=irstratset(nuge,i1,2)
            irstrat(i1,3)=irstratset(nuge,i1,3)
            irstrat(i1,4)=irstratset(nuge,i1,4)
            trstra(i1)=trstraset(nuge,i1)
            vkrt(i1)=vkrtset(nuge,i1)
            vkr2t(i1)=vkr2tset(nuge,i1)
        end do
        do i1=1,nrestrav
            irstrav(i1,1)=irstravset(nuge,i1,1)
            irstrav(i1,2)=irstravset(nuge,i1,2)
            irstrav(i1,3)=irstravset(nuge,i1,3)
            vrstra(i1)=vrstraset(nuge,i1)
            vkrv(i1)=vkrvset(nuge,i1)
            vkr2v(i1)=vkr2vset(nuge,i1)
        end do
        do i1=1,3
            axiss(i1)=axisset(nuge,i1)
            angles(i1)=anglesset(nuge,i1)
        end do
        do i1=1,na
            qa(i1)=qaset(nuge,i1)
            do i2=1,3
                c(i1,i2)=cset(nuge,i1,i2)
            end do
        end do

        if (iexco == 1) then
            vred=(1.0+(vvolset(nprob)-1.0d0))**(0.333333333333d0)
            na=naset(nuge-1)
            do i1=1,3
                axiss(i1)=vred*axis(i1)
                angles(i1)=angle(i1)
            end do
            do i1=1,na
                qa(i1)=qaset(nuge-1,i1)
                do i2=1,3
                !     c(i1,i2)=vred*cset(nuge-1,i1,i2)
                    c(i1,i2)=vred*clast(i1,i2)
                end do
            end do
            halfa=angles(1)*dgrrdn
            hbeta=angles(2)*dgrrdn
            hgamma=angles(3)*dgrrdn
            sinalf=sin(halfa)
            cosalf=cos(halfa)
            sinbet=sin(hbeta)
            cosbet=cos(hbeta)
            cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
            if (cosphi > 1.0d0) cosphi=1.0d0
            sinphi=sqrt(one-cosphi*cosphi)
            tm11=axiss(1)*sinbet*sinphi
            tm21=axiss(1)*sinbet*cosphi
            tm31=axiss(1)*cosbet
            tm22=axiss(2)*sinalf
            tm32=axiss(2)*cosalf
            tm33=axiss(3)
            kx=int(2.0d0*swb/tm11)
            ky=int(2.0d0*swb/tm22)
            kz=int(2.0d0*swb/tm33)
            qr='Y'
        end if

        if (qr == 'M' .OR. qr == 'A') then
            nmmsav=nmm
            nmm=2
        end if
        if (qr == 'A') nmm=1
        if (qr == 'D') endpo=endpo/25
        if (qr == 'D') nmmax=nmmax*5
        if (qr == 'H') nmmax=nmmax/10
        if (qr == '1' .OR. qr == '5') then
            nmm=1
            nmmax=1
        end if

    end if
          
!     if (qr.eq.'F'.or.qr.eq.'Y'.or.qr.eq.'3'.or.qr.eq.'5'.
!    $or.qr.eq.'P'.or.qr.eq.'B'.or.qr.eq.'S') then
    axis(1)=axiss(1)
    axis(2)=axiss(2)
    axis(3)=axiss(3)
    angle(1)=angles(1)
    angle(2)=angles(2)
    angle(3)=angles(3)
    if (axiss(1) < zero) then
        axis(1)=axis1
        axis(2)=axis2
        axis(3)=axis3
        angle(1)=angle1
        angle(2)=angle2
        angle(3)=angle3
    end if
    halfa=angle(1)*dgrrdn
    hbeta=angle(2)*dgrrdn
    hgamma=angle(3)*dgrrdn
    sinalf=sin(halfa)
    cosalf=cos(halfa)
    sinbet=sin(hbeta)
    cosbet=cos(hbeta)
    cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
    if (cosphi > 1.0d0) cosphi=1.0d0
    sinphi=sqrt(one-cosphi*cosphi)
    tm11=axis(1)*sinbet*sinphi
    tm21=axis(1)*sinbet*cosphi
    tm31=axis(1)*cosbet
    tm22=axis(2)*sinalf
    tm32=axis(2)*cosalf
    tm33=axis(3)
!     end if
          
    ngnh=3*na
!*********************************************************************
!                                                                    *
!     Write file headers                                             *
!                                                                    *
!*********************************************************************
    if (nsurp /= 2) then
        open (unit=69,file='summary.txt',status='unknown')
        write (69,1025)
        close (69)
        open (71,file='fort.71',status='unknown',position='append')
        if (nmm /= 0) write (71,1150) qmol
        write (71,1006)
        close (71)
        open (73,file='fort.73',status='unknown',position='append')
        if (nmm /= 0) write (73,1150) qmol
        write (73,1010)
        close (73)
        if (ntrc > 0) then
            open (75,file='fort.75',status='unknown',position='append')
            write (75,1020)
            close (75)
        end if
        if (nmethod == 4) then
            open (59,file='fort.59',status='unknown',position='append')
            write (59,1030)
            close (59)
        end if
    end if
!*********************************************************************
!                                                                    *
!     Determine atom types in system                                 *
!                                                                    *
!*********************************************************************
    xmasmd=0.0d0
    do i1=1,nso
        molin(nprob,i1)=0
        nasort(i1)=0
    end do
    do i1=1,na
        ia(i1,1)=0
        iag(i1,1)=0
        do i2=1,nso
            if (qa(i1) == qas(i2)) then
                ia(i1,1)=i2
                iag(i1,1)=i2
                molin(nprob,i2)=molin(nprob,i2)+1
                xmasat(i1)=amas(i2)
                xmasmd=xmasmd+amas(i2)
                nasort(i2)=nasort(i2)+1
            end if
        end do
        if (ia(i1,1) == 0) then
            write (*,*)'Unknown atom type: ',qa(i1)
            stop
        end if
    end do
!*********************************************************************
!                                                                    *
!     icentr=1: Place centre of mass at centre periodic box          *
!     icentr=2: Place centre of mass at origin                       *
!                                                                    *
!*********************************************************************
    if (icentr > 0) then
        ccx=0.0d0
        ccy=0.0d0
        ccz=0.0d0
        do i1=1,na
            ccx=ccx+c(i1,1)*xmasat(i1)/xmasmd
            ccy=ccy+c(i1,2)*xmasat(i1)/xmasmd
            ccz=ccz+c(i1,3)*xmasat(i1)/xmasmd
        end do
        xt2=-ccx
        yt2=-ccy
        zt2=-ccz
        if (icentr == 1) xt2=0.50d0*axis(1)-ccx
        if (icentr == 1) yt2=0.50d0*axis(2)-ccy
        if (icentr == 1) zt2=0.50d0*axis(3)-ccz
        do i1=1,na
            c(i1,1)=c(i1,1)+xt2
            c(i1,2)=c(i1,2)+yt2
            c(i1,3)=c(i1,3)+zt2
        end do
    end if
!*********************************************************************

    if (nreac == 2) then
        call default
        stop 'Placed default atoms (unit 98)'
    end if
         
!     if (nbiolab.eq.1) then
!     write (*,*)'Warning: generating new Biograf-labels'
!     call intcor
!     do i1=1,nmolo
!     na1m=nmolat(i1,1)
!     do i2=1,na1m
!     ihu=nmolat(i1,i2+1)
!     imolsta(ihu)=i1
!     write (qromb,'(f6.3)')float(i2)/1d3
!     qlabel(ihu)=qa(ihu)(1:1)//qromb(4:6)
!     end do
!     end do
!     end if
!*********************************************************************
!                                                                    *
!     MD-simulation loop                                             *
!                                                                    *
!*********************************************************************
    mdstep=0
    do i1=1,nbomax
        do i2=1,3
            do i3=1,2
                dbodc(i1,i2,i3)=0.0d0
            end do
        end do
    end do
    do i1=1,na
        vlp(i1)=0.0d0
        do i2=1,3
            vel(i2,i1)=0.0d0
            accel(i2,i1)=0.0d0
            aold(i2,i1)=0.0d0
        end do
    end do
    tempmd=tset
    call inivel
    irun=irun+2
    if (nvel == 1) then
        do i1=1,na
            do i2=1,3
                vel(i2,i1)=0.0d0
                accel(i2,i1)=0.0d0
                aold(i2,i1)=0.0d0
            end do
        end do
    end if
    if (ntest > 0) then
        call testaf
        stop 'Testing derivatives'
    end if

    if (nmm == 1) then          !Energy minimization
        call minim
        if (iredo > 1) then
            sumin=0.0d0
            do i1=1,na
                sumin=sumin+vincr(ia(i1,1))
            end do
            heatfo=sumin+estrc+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
            open (74,file='fort.74',status='unknown',position='append')
            if (nmm == 2 .AND. sdeva > endpo2) then
                heatfo=sumin+eavn+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
                write (74,1404)qmol,eavn,mdstep,heatfo,sdeva, &
                sumtt/float(mdstep)
            else
                if (axiss(1) < zero) write (74,1405)qmol,estrc,nit,heatfo
                if (axiss(1) > zero) write (74,1406)qmol,estrc,nit,heatfo, &
                tm11*tm22*tm33,(xmasmd*1.0d+24)/(avognr*tm11*tm22*tm33)
                if (axiss(1) > zero) write (62,1407)qmol,estrc/na,nit, &
                tm11*tm22*tm33/na,(xmasmd*1.0d+24)/(avognr*tm11*tm22*tm33)
            end if
            close (74)
            do i45=1,iredo
                icelo2=0
                call minim
                if (iredo > 1 .AND. i45 < iredo) then
                    sumin=0.0d0
                    do i1=1,na
                        sumin=sumin+vincr(ia(i1,1))
                    end do
                    heatfo=sumin+estrc+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
                    open (74,file='fort.74',status='unknown',position='append')
                    if (nmm == 2 .AND. sdeva > endpo2) then
                        heatfo=sumin+eavn+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
                        write (74,1404)qmol,eavn,mdstep,heatfo,sdeva, &
                        sumtt/float(mdstep)
                    else
                        if (axiss(1) < zero) write (74,1405)qmol,estrc,nit,heatfo
                        if (axiss(1) > zero) write (74,1406)qmol,estrc,nit,heatfo, &
                        tm11*tm22*tm33,(xmasmd*1.0d+24)/(avognr*tm11*tm22*tm33)
                        if (axiss(1) > zero) write (62,1407)qmol,estrc/na,nit, &
                        tm11*tm22*tm33/na,(xmasmd*1.0d+24)/(avognr*tm11*tm22*tm33)
                    end if
                    close (74)
                end if
            end do
        end if

        goto 800
    end if

    call distan
    call vlist
    if (iconne == 2) call intcon
    if (imolde == 1) call readmol
    if (iconne == 0 .OR. nbon == 0) then
        call srtbon1
    else
        if (iconne > 0 .AND. nbon > 0) call srtbon2
    end if
    call inilp
    call intcor
    call chargess
    if (iopt > 1 .AND. nmolo > 1) stop '>1 mol and ff optimization'

    tottime=0.0d0
    sum1=0.0d0
    sdev=0.0d0
    sdeva=0.0d0
    sum12=0.0d0
    sumt=0.0d0
    sump=0.0d0
    sumtt=0.0d0
    tmax=0.0d0
    eaver=0.0d0
    eav2=0.0d0
    eav3=0.0d0
    etot2s=0.0d0
    ediff=0.0d0
    idone=0d0
    idone2=0d0

!*********************************************************************
!                                                                    *
!     Calculate and output initial energy                            *
!                                                                    *
!*********************************************************************
    call encalc

    nmovh=0
    do i2=1,na
        do i3=1,3
            rmsg=rmsg+imove(i2)*d(i3,i2)*d(i3,i2)
            nmovh=nmovh+imove(i2)
        end do
    end do
          
    rmsg=sqrt(rmsg/float(nmovh))

    do i1=1,na
        do i2=1,3
            ekin=ekin+imove(i1)*xmasat(i1)*vel(i2,i1)*vel(i2,i1)
        end do
    end do
    ekin=0.50d0*ekin/convmd

    call mdsav(0,qfile(nprob))
    open (71,file='fort.71',status='unknown',position='append')
    open (73,file='fort.73',status='unknown',position='append')
    write (71,'(i8,2i4,1x,19(f10.2,1x))')mdstep+nprevrun,nmolo, &
    nmolo5,estrc,ekin,estrc+ekin,tempmd,estrc,estrc, &
    tempmd,tempmd,sump/float(nrep1),sdev,sdeva,tset, &
    tstep*1d+15,rmsg,tottime
    write (73,'(i8,1x,14(f10.2,1x))')mdstep+nprevrun,eb,ea,elp, &
    emol,ev+epen,ecoa,ehb,et,eco,ew,ep,ech,efi
    close (71)
    close (73)

    if (nmethod == 7) then
    !        xImpVcm = 1.d0 !nm/ps
        xImpVcm = xImpVcm*1d+13 !A/s
    elseif (nmethod == 8) then
    !        shock_vel = 2.d0
        shock_vel = shock_vel*1d13
        if (ishock_type == 1) then
            vwall = -1.*shock_vel
        else
            vwall = 0.d0
        endif
        call trarot2
        do i1=1,na
            if (c(i1,3) > shock_z_sep) then
                vel(3,i1) = vel(3,i1) - shock_vel
            endif
        enddo
    endif
    100 continue
    mdstep=mdstep+1
    if (iflext == 1) call chdt
    goto (110,120,130,140,150,160,170,180,190) nmethod
    stop 'Unrecognised method (check control-file)'
    110 call verlet1 !NVT-dynamics with Berendsen thermostat
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    call encalc
    call verlet2
    if (ilavel == 1) call layvel
    if (ntrc == 0) call tdamp(scasum)
    if (ntrc > 0) call tregime
    if (nvrc > 0) call vregime
    if (nerc > 0) call eregime
    if (npist > 0) call pireg
    goto 200
    120 tseto=tset !NVT-dynamics with Nose-Hoover thermostat
    tstepo=tstep
    call vernh1
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    call encalc
    call vernh2
    tset=tseto
    tstep=tstepo
    goto 200
    130 call verlet1 !NVE-dynamics
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    call encalc
    call verlet2
    if (nvrc > 0) call vregime
    if (nerc > 0) call eregime
    if (npist > 0) call pireg
    goto 200
    140 call verlet1 !NPT-dynamics
    icpres=1
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    call encalc
    call verlet2
    call pdamp
    if (ntrc == 0) call tdamp(scasum)
    if (ntrc > 0) call tregime
    if (nvrc > 0) call vregime
    if (nerc > 0) call eregime
    if (npist > 0) call pireg
    goto 200
    150 call verlet1 !NVE-dynamics; switch to NVT at high-T
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    call encalc
    call verlet2
    if (nvrc > 0) call vregime
    if (nerc > 0) call eregime
    if (npist > 0) call pireg
    if (invt == 1) call tdamp(scasum)
    goto 200
    160 call verlet1 !Hugoniostat Maillet PRE 2000
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    call encalc
!      call Hugoniostat
    etot = ekin+estrc
!		Hugoniostat variables
!     Hug_E0 = -45323.2
!     Hug_P0 = 1318.49 !in MPa
!     Hug_V0 = 6020.62173388
    Hug_chi_dot = (etot-Hug_E0-0.5d0* &
    (presz+Hug_P0)/1000.d0* &
    (Hug_V0-tm11*tm22*tm33)*.1438388139d0)/abs(Hug_E0)
    Hug_chi_dot =  (1.d0/taut)**2*Hug_chi_dot
    write(6,*)'Chi dot', Hug_chi_dot, taut, tstep
    write(6,*)"ene, target", etot, Hug_E0+0.5d0* &
    (presz+Hug_P0)/1000.d0* &
    (Hug_V0-tm11*tm22*tm33)*.1438388139d0
    write(6,*)'stress: ', presx,presy,presz

    Hug_chi = Hug_chi + Hug_chi_dot*tstep
    write(6,*)'chi, chi_dot: ', Hug_chi, Hug_chi_dot
    write(6,*)'Force: ', d(1,1), vel(1,i1)*xmasat(i1)/convmd, &
    Hug_chi*vel(1,i1)*xmasat(i1)/convmd
    do i1=1,na
        d(1,i1) = d(1,i1) + Hug_chi*vel(1,i1)*xmasat(i1)/convmd
        d(2,i1) = d(2,i1) + Hug_chi*vel(2,i1)*xmasat(i1)/convmd
        d(3,i1) = d(3,i1) + Hug_chi*vel(3,i1)*xmasat(i1)/convmd
    enddo
    call verlet2
!      if (ntrc.eq.0) call tdamp(scasum)
    goto 200
    170 continue !shear
    write(6,*)'Shear'
    write(6,*)'Impulse vel:', xImpVcm, 1d+13
!		ImpVcm = 0.d0
    vcm1 = 0.d0
    xmass1 = 0.d0
    vcm2 = 0.d0
    xmass2 = 0.d0
    do i1=1,na
        if (c(i1,3) > tm33/2.d0) then
            vcm2 = vcm2+ xmasat(i1) * vel(1,i1)
            xmass2 = xmass2+ xmasat(i1)
        else
            vcm1 = vcm1+ xmasat(i1) * vel(1,i1)
            xmass1 = xmass1+ xmasat(i1)
        endif
    enddo
    vcm1 = vcm1/xmass1
    vcm2 = vcm2/xmass2
    dvcm1 = ximpVcm - vcm1
    dvcm2 = -ximpVcm - vcm2
    write(6,*)'Vcm: ', vcm1, vcm2
!      write(6,*)'mass: ', xmass1, xmass2
!      write(6,*)'calc: ', ximpVcm, vcm1, dvcm1
    write(6,*)'Delta Vcm: ', dvcm1, dvcm2
    do i1=1,na
        if (c(i1,3) > tm33/2.d0) then
            vel(1,i1) = vel(1,i1) + dvcm2
        else
            vel(1,i1) = vel(1,i1) + dvcm1
        endif
    enddo
    call verlet1 !NVE-dynamics
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    call encalc
    call verlet2
    goto 200
    180 call verlet1 !NEMD shock (compressing boundary conditions)
    if (nreac == 0) call intcor
    if (nreac == 1) call distan
    axis(3) = axis(3)+vwall*tstep
    axiss(3)=axis(3)
    write(6,*)'lattice parameter: ', axis(3)
    call encalc
    call verlet2
    if (nvrc > 0) call vregime
    if (nerc > 0) call eregime
    goto 200
    190 stop 'Unrecognised method (check control-file)'
    200 continue
!*********************************************************************
!                                                                    *
!     Change distance restraints                                     *
!                                                                    *
!*********************************************************************
    do i1=1,nrestra
        if (itend(i1) == 0 .OR. (mdstep > itstart(i1) .AND. mdstep < &
        itend(i1))) then
            rrstra(i1)=rrstra(i1)+rrcha(i1)
        end if
    end do
    do i1=1,nrestrav
        vrstra(i1)=vrstra(i1)+rvcha(i1)
    end do
    do i1=1,nrestrat
        trstra(i1)=trstra(i1)+rtcha(i1)
    end do
    do i1=1,nrestram
        rmstra1(i1)=rmstra1(i1)+rmcha(i1)
    end do
!*********************************************************************
!                                                                    *
!     Check temperature regime                                       *
!                                                                    *
!*********************************************************************
    tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
    if (ntrc == 0) then
        tset=tset+itemp*tincr
        if (tset < 0.1d0) tset=0.1d0
    end if
!*********************************************************************
!                                                                    *
!     Call runtime analysis routine                                  *
!                                                                    *
!*********************************************************************
    if (ianaly == 2) call runanal
!*********************************************************************
!                                                                    *
!     Add atom                                                       *
!                                                                    *
!*********************************************************************
    if (nrdd == 3) then
        if (mod(mdstep,nrddf) == 0 .OR. mdstep == 1) then
            ioldchg=0
            na=na+1
            c(na,1)=axiss(1)*random(dseed)
            c(na,2)=30.0d0
            c(na,3)=axiss(3)*random(dseed)
            ia(na,1)=ityrad
            iag(na,1)=ityrad
            qa(na)=qas(ityrad)
            xmasat(na)=amas(ia(na,1))
            xmasmd=xmasmd+amas(ia(na,1))
            molin(nprob,4)=molin(nprob,4)+1
            nasort(4)=nasort(4)+1
            vel(1,na)=random(dseed)-0.5d0
            vel(2,na)=-random(dseed)
            vel(3,na)=random(dseed)-0.5d0
            velsq=xmasat(na)*(vel(1,na)*vel(1,na)+vel(2,na)*vel(2,na)+ &
            vel(3,na)*vel(3,na))
            ekinat=half*velsq/convmd
            tempn=2.0d0*ekinat/(float(3)*rgasc*xjouca/1.0d3)
            factor=sqrt(tpnrad/tempn)
            vel(1,na)=vel(1,na)*factor
            vel(2,na)=vel(2,na)*factor
            vel(3,na)=vel(3,na)*factor
        end if
    end if
!*********************************************************************
!                                                                    *
!     Add molecule to system                                         *
!                                                                    *
!*********************************************************************
    if (iaddfreq > 0) then
        if (mod(mdstep,iaddfreq) == 0 .OR. mdstep == 5) then

            ioldchg=0  !get new chargess
            iattempt=0
            iaccplace=0
            32 iattempt=iattempt+1
                  
            if (icentp /= 1) then   !Put molecule at random location
                xaddh=axiss(1)*(random(dseed)-0.5d0)
                yaddh=axiss(2)*(random(dseed)-0.5d0)
                zaddh=axiss(3)*(random(dseed)-0.5d0)
            else    !Put molecule close to centre-of-mass
                ccx=0.0d0
                ccy=0.0d0
                ccz=0.0d0
                do i1=1,na
                    ccx=ccx+c(i1,1)*xmasat(i1)/xmasmd
                    ccy=ccy+c(i1,2)*xmasat(i1)/xmasmd
                    ccz=ccz+c(i1,3)*xmasat(i1)/xmasmd
                end do
            !     write (64,*)iattempt,ccx,ccy,ccz,xmasmd
                xaddh=ccx+dcentm*(random(dseed)-0.5d0)
                yaddh=ccy+dcentm*(random(dseed)-0.5d0)
                zaddh=ccz+dcentm*(random(dseed)-0.5d0)
            end if
                  
            if (xadd > -5000.0d0) xaddh=xadd
            if (yadd > -5000.0d0) yaddh=yadd
            if (zadd > -5000.0d0) zaddh=zadd

            do i1=1,naa
                c(na+i1,1)=cadd(i1,1)+xaddh
                c(na+i1,2)=cadd(i1,2)+yaddh
                c(na+i1,3)=cadd(i1,3)+zaddh
                qa(na+i1)=qadd(i1)
                ia(na+i1,1)=ityadd(i1)
                iag(na+i1,1)=ityadd(i1)
                xmasat(na+i1)=amas(ia(na+i1,1))
                xmasmd=xmasmd+xmasat(na+i1)
                vel(1,na+i1)=random(dseed)-0.5d0
                vel(2,na+i1)=random(dseed)-0.5d0
                vel(3,na+i1)=random(dseed)-0.5d0
            end do
        !*********************************************************************
        !                                                                    *
        !     Check distances between added molecule and rest of system      *
        !                                                                    *
        !*********************************************************************
            distmin=100.0d0
            distminc=100.0d0
            kxt=kx
            kyt=ky
            kzt=kz
            if (kxt == 0) kxt=1
            if (kyt == 0) kyt=1
            if (kzt == 0) kzt=1
            do i1=1,naa
                do i2=1,na
                    dx1=c(i2,1)-c(na+i1,1)
                    dy1=c(i2,2)-c(na+i1,2)
                    dz1=c(i2,3)-c(na+i1,3)
                    do k1=-kxt,kxt
                        do k2=-kyt,kyt
                            do k3=-kzt,kzt
                                a1=dx1+k1*tm11
                                a2=dy1+k1*tm21+k2*tm22
                                a3=dz1+k1*tm31+k2*tm32+k3*tm33
                                rr=sqrt(a1*a1+a2*a2+a3*a3)
                                if (rr < distmin) then
                                    distmin=rr
                                end if
                                if (ia(i2,1) == 1) then    !find distance to carbon atoms
                                    if (rr < distminc) then
                                        distminc=rr
                                    end if
                                end if

                            end do
                        end do
                    end do
                end do
            end do
                 
            if (distmin > addist .AND. distminc > addist*2.0d0) then
                iaccplace=1
            else
                if (iattempt < nadattempt) then
                    do i1=1,naa
                        xmasmd=xmasmd-xmasat(na+i1)
                    end do
                    goto 32
                else
                    write (*,*)'Could not place new molecule in iteration:',mdstep
                end if
            end if

            if (iaccplace == 1) then
                velsq=zero
                do i1=1,naa
                    velsq=velsq+xmasat(na+i1)*(vel(1,na+i1)*vel(1,na+i1)+ &
                    vel(2,na+i1)*vel(2,na+i1)+vel(3,na+i1)*vel(3,na+i1))
                end do

                ekinat=half*velsq/convmd
                tempadd=2.0d0*ekinat/(float(3*naa)*rgasc*xjouca/1.0d3)
                temph=taddmol
                if (temph < zero) temph=tempmd
                factor=sqrt(temph/tempadd)

                do i1=1,naa
                    vel(1,na+i1)=vel(1,na+i1)*factor
                    vel(2,na+i1)=vel(2,na+i1)*factor
                    vel(3,na+i1)=vel(3,na+i1)*factor
                end do

                if (iveladd == 2) then
                    do i1=1,naa
                        vel(1,na+i1)=veladd(1,i1)
                        vel(2,na+i1)=veladd(2,i1)
                        vel(3,na+i1)=veladd(3,i1)
                    end do
                end if

                do i1=1,naa
                    accel(1,na+i1)=zero
                    accel(2,na+i1)=zero
                    accel(3,na+i1)=zero
                end do

                do i1=1,naa
                    velsq=velsq+xmasat(na+i1)*(vel(1,na+i1)*vel(1,na+i1)+ &
                    vel(2,na+i1)*vel(2,na+i1)+vel(3,na+i1)*vel(3,na+i1))
                end do

                ekinat=half*velsq/convmd
                tempadd=2.0d0*ekinat/(float(3*naa)*rgasc*xjouca/1.0d3)

                write (*,*)'Placed new molecule in iteration:',mdstep
                write (*,*)'Location:'
                do i1=1,naa
                    write (*,'(i6,3f15.5)')i1,c(na+i1,1),c(na+i1,2),c(na+i1,3)
                end do
                write (*,*)'Velocity:'
                do i1=1,naa
                    write (*,'(i6,3d20.5)')i1,vel(1,na+i1),vel(2,na+i1),vel(3,na+i1)
                end do
                write (*,'(a25,f10.2)')'Temperature new molecule:',tempadd
                na=na+naa
                call vlist
            end if
                  
        end if
    end if
!*********************************************************************
!                                                                    *
!     Output MD-statistics and atom coordinates                      *
!                                                                    *
!*********************************************************************
    rmsg=0.0
    nmovh=0
          
    do k=1,3
        totf1(k)=0.0d0
        totf2(k)=0.0d0
        totf3(k)=0.0d0
    end do

    do i2=1,na
        do i3=1,3
            rmsg=rmsg+imove(i2)*d(i3,i2)*d(i3,i2)
            nmovh=nmovh+imove(i2)
            if(ibgr2(i2) <= 0) then
                totf1(i3)=totf1(i3)+d(i3,i2)
            else if (ibgr2(i2) > 0) then
                totf2(i3)=totf2(i3)+d(i3,i2)
            !     else if(iforc(i2).eq.2) then
            !     totf3(i3)=totf3(i3)+d(i3,i2)
            end if                                    !Qing
        end do
    end do

    rmsg=sqrt(rmsg/float(nmovh))

    do k=1,3
        avef1(k)=(avef1(k)*(mdstep-1)+totf1(k))/float(mdstep)
        avef2(k)=(avef2(k)*(mdstep+nprevrun-1)+ &
        totf2(k))/float(mdstep+nprevrun)
    !     avef3(k)=(avef3(k)*(mdstep-1)+totf3(k))/float(mdstep)
    end do

    tottime=tottime+tstep*1d+15
    sum1=sum1+estrc
    sum12=sum12+estrc*estrc
    sumt=sumt+tempmd
    sumtt=sumtt+tempmd
    sump=sump+pressu/1000.0d0
    eaver=eaver+estrc
    eav2=eav2+eaver/float(mdstep)
    eav3=eav3+(eaver/float(mdstep))*(eaver/float(mdstep))
    if (tempmd > tmax) tmax=tempmd
    if (mod(mdstep,nrep1) == 0) then
    !     scale=1.0+scasum/float(nrep1*na)
        if (nrep1 > 1) &
        sdev=sqrt((sum12-sum1*sum1/float(nrep1))/float(nrep1-1))
        eavn=eaver/float(mdstep)
        if (mdstep > 1) &
        sdeva=sqrt((eav3-eav2*eav2/float(mdstep))/float(mdstep-1))
        open (71,file='fort.71',status='unknown',position='append')
        open (73,file='fort.73',status='unknown',position='append')
        write (71,'(i8,2i4,1x,19(f10.2,1x))')mdstep+nprevrun,nmolo, &
        nmolo5,estrc,ekin,estrc+ekin,tempmd,sum1/float(nrep1),eavn, &
        sumt/float(nrep1),tmax,sump/float(nrep1),sdev,sdeva,tset, &
        tstep*1d+15,rmsg,tottime
        write (73,'(i8,1x,14(f10.2,1x))')mdstep+nprevrun,eb,ea,elp, &
        emol,ev+epen,ecoa,ehb,et,eco,ew,ep,ech,efi
        close (71)
        close (73)
        if (ilavel == 1) then
            open (68,file='force.out',status='unknown',position='append')
            write (68,*)mdstep,(totf1(k),k=1,3),(totf2(k),k=1,3), &
            (avef1(k),k=1,3),(avef2(k),k=1,3)
            close (68)
        end if

        if ((sumt/float(nrep1)) > tset .AND. nmethod == 5) then
            if (invt == 0) write (*,*)'Switched to NVT in iteration',mdstep
            invt=1
        end if
        if (nrestra > 0 .OR. nrestrat > 0) &
        open (76,file='fort.76',status='unknown',position='append')
        if (nrestra > 0) then
            do i2=1,nrestra
                call dista2(irstra(i2,1),irstra(i2,2),disres(i2),dx,dy,dz)
            end do
            open (76,file='fort.76',status='unknown',position='append')
            write (76,'(i8,1x,40f12.4)')mdstep,eres,estrc, &
            (rrstra(i2),disres(i2),i2=1,nrestra)
        end if

        if (nrestrat > 0) then
            open (76,file='fort.76',status='unknown',position='append')
            do i2=1,nrestrat
                do i3=1,ntor
                    ih1=irstrat(i2,1)
                    ih2=irstrat(i2,2)
                    ih3=irstrat(i2,3)
                    ih4=irstrat(i2,4)
                    if (ih1 == it(i3,2) .AND. ih2 == it(i3,3) .AND. ih3 == it(i3,4) &
                     .AND. ih4 == it(i3,5)) ittr=i3
                end do
                write (76,'(i8,1x,40f12.4)')mdstep,eres, &
                trstra(i2),thg(ittr)
            end do
        end if

        if (nrestra > 0 .OR. nrestrat > 0) close(76)

        if (nrestram > 0) then
            open (76,file='fort.76',status='unknown',position='append')
            do i2=1,nrestram
                write (76,'(2i8,1x,20f12.4)')mdstep,i2,eres,rmstra1(i2), &
                dismacen(i2)
            end do
            close (76)
        end if

        if (nmm == 2 .AND. iequi == 1 .AND. rmsg < endpo2) then
            sumin=0.0
            do i1=1,na
                sumin=sumin+vincr(ia(i1,1))
            end do
            heatfo=sumin+estrc+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
            open (74,file='fort.74',status='unknown',position='append')
            if (nrestra == 0) write (74,1400)qmol,estrc,sdeva,sumt,heatfo
            if (nrestra > 0) write (74,1410)(rrstra(i2),i2=1,nrestra),heatfo
            close (74)
            estrmin=estrc

            if (qr == 'M') call minim
            goto 800
        end if

        sum1=0.0d0
        sum12=0.0d0
        sumt=0.0d0
        sump=0.0d0
        tmax=0.0d0
    end if
    if (mod(mdstep,nhop2) == 0 .AND. na > 2 .AND. nmethod /= 3) &
    call trarot1
    if (mod(mdstep,ncontrol) == 0) call readc
    if (mod(mdstep,nsav) == 0) call mdsav(0,qfile(nprob))
    if (mod(mdstep,nsav3) == 0) call mdsav(2,qfile(nprob))
    if (mdstep > nequi .AND. iequi == 0) then
        mdstep=0
        iequi=1
        eaver=0.0d0
        eav2=0.0d0
        eav3=0.0d0
        sum1=0.0d0
        sum12=0.0d0
        sumtt=0.0d0
    end if
    if (mdstep < maxstp) goto 100

    800 if (nmm > 0) then
        if (qr == 'M') call minim
        sumin=0.0d0
        do i1=1,na
            sumin=sumin+vincr(ia(i1,1))
        end do

        if (iopt == 0 .AND. nsurp < 2) &
        call dipmom(naold,dpmm,xdip,ydip,zdip,xdir,ydir,zdir)

        heatfo=sumin+estrc+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
        open (74,file='fort.74',status='unknown',position='append')
        if (nmm == 2 .AND. sdeva > endpo2) then
            heatfo=sumin+eavn+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
            write (74,1404)qmol,eavn,mdstep,heatfo,sdeva, &
            sumtt/float(mdstep)
        else
            if (axiss(1) < zero) write (74,1405)qmol,estrc,nit,heatfo
            if (axiss(1) > zero) write (74,1406)qmol,estrc,nit,heatfo, &
            tm11*tm22*tm33,(xmasmd*1.0d+24)/(avognr*tm11*tm22*tm33)
            if (axiss(1) > zero) write (62,1407)qmol,estrc/na,nit, &
            tm11*tm22*tm33/na,(xmasmd*1.0d+24)/(avognr*tm11*tm22*tm33)
        end if
        close (74)
    !*********************************************************************
    !                                                                    *
    !     Output restraint information                                   *
    !                                                                    *
    !*********************************************************************
        if (nrestra+nrestrav+nrestrat+nrestram+nrestras+neqdis > 0 &
         .AND. iopt == 0) then
        !     open (76,file='fort.76',status='unknown')

            if (nrestra > 0 .AND. iopt == 0) then
                do i2=1,nrestra
                    call dista2(irstra(i2,1),irstra(i2,2),disres(i2),dx,dy,dz)
                end do
                write (76,'(i8,1x,40f12.4)')mdstep,eres,estrc, &
                (rrstra(i2),disres(i2),i2=1,nrestra)
            end if

            if (nrestrav > 0 .AND. iopt == 0) then
                do i2=1,nrestrav
                    ittr=0
                    do i3=1,nval
                        ih1=irstrav(i2,1)
                        ih2=irstrav(i2,2)
                        ih3=irstrav(i2,3)
                        if (ih1 == iv(i3,2) .AND. ih2 == iv(i3,3) .AND. ih3 == iv(i3,4) &
                        ) ittr=i3
                    end do
                    write (76,1411)rdndgr*h(ittr),heatfo
                end do
            end if

            if (nrestrat > 0 .AND. iopt == 0) then
                do i2=1,nrestrat
                    ittr=0
                    do i3=1,ntor
                        ih1=irstrat(i2,1)
                        ih2=irstrat(i2,2)
                        ih3=irstrat(i2,3)
                        ih4=irstrat(i2,4)
                        if (ih1 == it(i3,2) .AND. ih2 == it(i3,3) .AND. ih3 == it(i3,4) &
                         .AND. ih4 == it(i3,5)) ittr=i3
                    end do
                    write (76,1411)thg(ittr),heatfo
                end do
            end if

            if (nrestram > 0 .AND. iopt == 0) then
                do i2=1,nrestram
                    write (76,'(2i8,1x,20f12.4)')mdstep,i2,eres,rmstra1(i2), &
                    dismacen(i2)
                end do
            end if

            if (nrestras > 0 .AND. iopt == 0) then
                do i2=1,nrestras
                    write (76,'(2i8,1x,a6,2f12.4)')mdstep,i2,qrstras(i2),eres, &
                    sysdissum(i2)
                end do
            end if

            if (neqdis > 0 .AND. iopt == 0) then
                do i2=1,neqdis
                    call dista2(ieqdis(i2,1),ieqdis(i2,2),reqd1,dx,dy,dz)
                    call dista2(ieqdis(i2,3),ieqdis(i2,4),reqd2,dx,dy,dz)
                    write (76,'(i8,1x,2f12.4,4i8,2f12.4)')mdstep,eres,estrc, &
                    (ieqdis(i2,i3),i3=1,4),reqd1,reqd2
                end do
            end if

        !     close (76)
        end if

        estrmin=estrc
        enmolend(nprob)=estrc
        formolend(nprob)=rmsg

        itrain=0
        if (ifreq == 1 .AND. iopt == 0 .AND. nsurp < 2) &
        call vibra(qfreqfile,vibreax,vibqc,imatch,errmatch,itrain,klinear)
        if (nsurp /= 2 .AND. icgeo == 0) call mdsav(1,qfile(nprob))
        call outres2
        if (iopt == 0 .AND. nsurp /= 2 .AND. qr == 'C') call outint
        call readc
        if (qr == 'M' .OR. qr == 'A') then
            nmm=nmmsav
        end if
        do i1=1,na
            do i2=1,3
                clast(i1,i2)=c(i1,i2)
            end do
        end do
        if (nsurp /= 2) then
            if (icgeo == 0) then
                do i1=1,na
                    do i2=1,3
                        cset(nprob,i1,i2)=c(i1,i2)
                    end do
                end do
            else
                do i1=1,na
                    do i2=1,3
                        c(i1,i2)=cset(nprob,i1,i2)   !Restore start coordinates
                    end do
                end do
                do i1=1,3
                    axiss(i1)=axisset(nprob,i1)
                    angles(i1)=anglesset(nprob,i1)
                end do
                call mdsav(1,qfile(nprob))
            end if

            goto 5     !Next geometry
        end if

        if (nsurp >= 2) then
            if ((iopt == 0 .OR. iflga == 1) .AND. icgeopt(nuge) == 0) then
                do i1=1,na
                    do i2=1,3
                        cset(nuge,i1,i2)=c(i1,i2)
                    end do
                end do
                do i1=1,3
                    axisset(nuge,i1)=axiss(i1)
                    anglesset(nuge,i1)=angles(i1)
                end do
            end if
            nuge=nuge+1
            if (nuge > nmollset) goto 900
            goto 50
        end if

    end if
    if (nmm == 0) call mdsav(1,qfile(nprob))
    900 continue
    if (imodfile == 1 .AND. iline == 0) then
        write (*,*)'Could not find geometries in models.in'
        stop 'Could not find geometries in models.in'
    end if
    close (14)
    call outresend
    if (ianaly == 1) call analysis(naold)
    open (72,file='fort.72',status='unknown',position='append')
    call timer(72)
    close (72)
!     call etime(tarray)
!     totime=tarray(1)
!     eltme=totime-prvtme
!     write (72,1300)eltme
    return
! 999 stop 'Error writing binary file'
!     return
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    1000 format (2x,a1,1x,a60)
    1001 format (2x,a1,7x,a60)
    1005 format (3f10.4)
    1006 format ('   Iter. Nmol    Epot         Ekin      Etot ', &
    '       T(K)  Eaver(block) Eaver(total) Taver      Tmax   ', &
    ' Pres(GPa)   sdev(Epot)  sdev(Eaver)    Tset      Timestep', &
    '    RMSG     Totaltime')
    1010 format ('  Iter.      Etotal      Ebond       Eatom       Elp        Emol', &
    '       Eval       Ecoa       Ehbo       Etors      Econj', &
    '      Evdw      Ecoul    Echarge      Efield     RMSG     BO(1,1)    d(1,1)')
    1020 format ('   Iter.  Tsys    Tzone1  Tset1   Tzone2  Tset2')
    1025 format ('  Iteration Nmol    Time(fs)  Epot(kcal)', &
    '    Vol(A^3)     T(K)       Pres(MPa)    Dens(kg/dm3)')
    1030 format ('    Iter.     a             b          c        px', &
    '(GPa)    py(GPa)      pz(GPa)     pset(GPa)  Volume ')
    1050 format (f7.3)
    1055 format (f7.4)
    1060 format (i8)
    1080 format (i6,7x,f8.4)
    1100 format (i4,1x,a2,3x,3d22.15,1x,a5,1x,i5)
    1110 format (i4,2f10.2)
    1150 format (1x,a60)
    1200 format (3f12.4)
    1220 format (10i4)
    1240 format (i4,f12.6)
    1300 format ('Total elapsed time:',f12.2,' seconds')
    1400 format (1x,a20,'Emin: ',f12.3,' Sdev.: ',f10.3, &
    ' Temp.: ',f10.3,' Heatfo: ',f10.3)
    1404 format (1x,a30,'[E]: ',f12.2,' MD-Iterations:',i6 &
    ,' Heatfo: ',f8.2,' Error: ',f8.2,' [T]: ',f8.2)
    1405 format (1x,a20,'Emin: ',f12.3,' Iter.:',i6 &
    ,' Heatfo: ',f10.3)
    1406 format (1x,a20,'Emin: ',f12.3,' Iter.:',i6 &
    ,' Hf: ',f10.3,' Vol: ',f10.3, ' Dens): ',f10.3)
    1407 format (1x,a20,'Emin/atom: ',f12.3,' Iter.:',i6 &
    ,' Vol/atom: ',f10.3, ' Dens): ',f10.3)
    1410 format (f12.3,' Heatfo: ',f10.3)
    1411 format (f12.3,' Heatfo: ',f10.3)
    1415 format (f12.3,' Heatfo: ',f10.3,' Estrc: ',f10.3)
    2000 format ('Geometry identifier: ',a40)
    end subroutine reac
!*********************************************************************
!*********************************************************************

    subroutine verlet1

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Velocity Verlet time step part 1                               *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In verlet1'
        call timer(65)
        close (65)
    end if
    do i1=1,na
        do i2=1,3
            if (imove(i1) == 1) &
            c(i1,i2)=c(i1,i2)+vel(i2,i1)*tstep+accel(i2,i1)*ts22
            vel(i2,i1)=imove(i1)*vel(i2,i1)+ts2*accel(i2,i1)
        end do
    end do
    return
    end subroutine verlet1
!*********************************************************************
!*********************************************************************

    subroutine verlet2

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Velocity Verlet time step part 2                               *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In verlet2'
        call timer(65)
        close (65)
    end if
    do i1=1,na
        do i2=1,3
        !     aold(i2,i1)=accel(i2,i1)
            accel(i2,i1)=-convmd*d(i2,i1)/xmasat(i1)
            vel(i2,i1)=imove(i1)*vel(i2,i1)+accel(i2,i1)*ts2
            ekin=ekin+imove(i1)*xmasat(i1)*vel(i2,i1)*vel(i2,i1)
        end do
    end do
    ekin=0.50d0*ekin/convmd
    eksav=ekin
    return
    end subroutine verlet2
!*********************************************************************
!*********************************************************************

    subroutine chdt

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Change time-step as a function of maximum velocity             *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In chdt'
        call timer(65)
        close (65)
    end if
    vmax=zero
    do i1=1,na
        vhulp=zero
        do i2=1,3
            vhulp=vhulp+vel(i2,i1)*vel(i2,i1)
        end do
        if (vhulp > vmax) vmax=vhulp
    end do
        
    vmax2=sqrt(vmax)
    tstepmax=tstep0*1d-15
    tstepmin=1.0d-19
    tstep=dtmax/vmax2
    write (64,*)mdstep,vmax2,dtmax,tstepmax,tstepmin,tstep
    if (tstep > tstepmax) tstep=tstepmax
    if (tstep < tstepmin) tstep=tstepmin
    ts2=tstep/2.0d0
    ts22=tstep*ts2
    return
    end subroutine chdt
!*********************************************************************
!*********************************************************************

    subroutine tdamp(scasum)

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Temperature control using temperature damping                  *
!     (Berendsen et al., J. Chem. Phys. 1984,81,3684)                *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In tdamp'
        call timer(65)
        close (65)
    end if
    if (ncons == 2) then
        tempmds=tempmd
        5 tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
        sca1=1.0d0+(tstep/taut)*(tset/tempmd-1.0d0)
        if (sca1 < 0.0d0) sca1=0.0d0
        scale=sqrt(sca1)
        if (ntscale == 1) scale=sqrt(tset/(tempmd+0.1d0))
        if (scale > 2.0d0) scale=2.0d0
        do i1=1,na
            do i2=1,3
                vel(i2,i1)=imove(i1)*scale*vel(i2,i1)
                ekin=ekin+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        ekin=0.50d0*ekin/convmd
        tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
        if (tempmd > 2.0d0*tempmds .AND. mdstep > 100) then
            open (75,file='fort.75',status='unknown',position='append')
            write (75,*) mdstep,tempmd,tempmds
            close (75)
        end if
    end if
!*********************************************************************
!                                                                    *
!     Temperature control using temperature damping                  *
!     applied to individual atoms                                    *
!                                                                    *
!*********************************************************************
    if (ncons == 0) then
        scasum=1.0d0
        do i1=1,na
            if (imove(i1) == 1) then
                ekinat=0.0d0
                do i2=1,3
                    ekinat=ekinat+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
                end do
                ekinat=0.50d0*ekinat/convmd
                tempat=2.0d0*ekinat/(3.0*rgasc*xjouca/1.0d3)
                10 sca1=1.0d0+(tstep/taut)*(tset/tempat-1.0d0)
                if (sca1 < 0.0d0) sca1=0.0d0
                scale=sqrt(sca1)
                if (ntscale == 1) scale=sqrt(tset/(tempmd+0.1d0))
                if (scale > 2.0d0) scale=2.0d0
                scasum=scasum*scale
                ekinat=0.0d0
                do i2=1,3
                    vel(i2,i1)=vel(i2,i1)*scale
                    ekinat=ekinat+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
                end do
                ekinat=0.50d0*ekinat/convmd
                tempat=2.0d0*ekinat/(3.0d0*rgasc*xjouca/1.0d3)
                if (tempat > tset+5000.0d0) goto 10
                if (tempat < tset-5000.0d0) goto 10
            end if
        end do

        ekin=zero
        do i1=1,na
            do i2=1,3
                ekin=ekin+xmasat(i1)*imove(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        ekin=0.50d0*ekin/convmd
    !     tempmd=2.0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
        scasum=scasum*1.0d0/float(na)
    end if
!*********************************************************************
!                                                                    *
!     Temperature control using temperature damping                  *
!     applied to individual molecules                                *
!                                                                    *
!*********************************************************************
    if (ncons == 3) then
        if (mod(mdstep,nrep1) == 0) then
            open (75,file='fort.75',status='unknown')
        end if
        do i1=1,nmolo
            ekinmol=0.0
            na1m=nmolat(i1,1)
            do i2=1,na1m
                ihu=nmolat(i1,i2+1)
                do i3=1,3
                    ekinmol=ekinmol+xmasat(ihu)*vel(i3,ihu)*vel(i3,ihu)
                end do
            end do
            ekinmol=0.50d0*ekinmol/convmd
            tempmol=2.0d0*ekinmol/(float(3*na1m)*rgasc*xjouca/1.0d3)
            20 sca1=1.0d0+(tstep/taut)*(tset/tempmol-1.0d0)
            if (sca1 < 0.0d0) sca1=0.0d0
            scale=sqrt(sca1)
            if (mod(mdstep,nrep1) == 0) then
                write (75,'(2i8,4f12.6)')mdstep,i1,tempmol,tset, &
                taut*1e+15,scale
            end if
            ekinmol=0.0d0
            do i2=1,na1m
                ihu=nmolat(i1,i2+1)
                do i3=1,3
                    vel(i3,ihu)=vel(i3,ihu)*scale
                    ekinmol=ekinmol+xmasat(ihu)*vel(i3,ihu)*vel(i3,ihu)
                end do
            end do
            ekinmol=0.50d0*ekinmol/convmd
            tempmol=2.0d0*ekinmol/(float(3*na1m)*rgasc*xjouca/1.0d3)
        !     if (tempmol.gt.2.0*tset) goto 20
        end do

        ekin=zero
        do i1=1,na
            do i2=1,3
                ekin=ekin+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        ekin=0.50d0*ekin/convmd
    !     tempmd=2.0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
        if (mod(mdstep,nrep1) == 0) then
            close (75)
        end if
    end if
          
!*********************************************************************
!                                                                    *
!     Enforce energy conservation                                    *
!     ediff+estrc should remain constant during the simulation       *
!                                                                    *
!*********************************************************************
    if (ncons == 1) then
        stop 'Energy conservation does not work'
    !     if (abs(tempmd-tset).lt.2.5*sqrt(taut).and.idone.eq.0) then
    !     idone2=idone2+1
    !     if (idone2.gt.10) then
    !     etot2s=estrc+ekin
    !     idone=1
    !     end if
    !     end if
    !     ediff=ediff+eksav-ekin
    !     etot2=estrc+ekin+ediff
    !     if (idone.eq.1) then
    !     hulp1=etot2-etot2s
    !     hu12=(1.0-(hulp1/ekin))
    !     hulp2=1.0
    !     if (hu12.gt.zero) hulp2=sqrt(hu12)
    !     do i1=1,na
    !     do i2=1,3
    !     vel(i2,i1)=vel(i2,i1)*hulp2
    !     end do
    !     end do
    !     end if
    end if
!*********************************************************************
!                                                                    *
!     Temperature control using Anderson's thermostat                *
!     (Anderson, J. Chem. Phys. 1980,72,2384)                        *
!                                                                    *
!*********************************************************************
    if (ncons == 4) then
        hulp=(tstep/taut)
        do i1=1,na
            if (random(dseed) < hulp) then   !collision with bath
                v0=sqrt(3.0d0*convmd*rgasc*tset/(1.0d3*caljou*xmasat(i1)))
                sigma=sqrt(v0)
                velnew=gauss(sigma,v0,dseed)
                velold=sqrt(vel(1,i1)*vel(1,i1)+vel(2,i1)*vel(2,i1)+ &
                vel(3,i1)*vel(3,i1))
                scale=velnew/velold
                scale2=sqrt(scale)
                do i2=1,3
                    vel(i2,i1)=scale2*vel(i2,i1)
                end do
            end if

        end do

        ekin=zero
        do i1=1,na
            do i2=1,3
                ekin=ekin+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        ekin=0.50d0*ekin/convmd

    end if
!*********************************************************************
!                                                                    *
!     Temperature control using temperature damping                  *
!     applied to individual molecules with 2 different damping       *
!     constants                                                      *
!                                                                    *
!*********************************************************************
    if (ncons == 5) then
        if (mod(mdstep,nrep1) == 0) then
            open (75,file='fort.75',status='unknown')
        end if
        do i1=1,nmolo
            ekinmol=0.0
            na1m=nmolat(i1,1)
            do i2=1,na1m
                ihu=nmolat(i1,i2+1)
                do i3=1,3
                    ekinmol=ekinmol+xmasat(ihu)*vel(i3,ihu)*vel(i3,ihu)
                end do
            end do
            ekinmol=0.50d0*ekinmol/convmd
            tempmol=2.0d0*ekinmol/(float(3*na1m)*rgasc*xjouca/1.0d3)
            ekinmol=0.0d0
            do i2=1,na1m
                ihu=nmolat(i1,i2+1)
                tauhu=taut
                if (ihu > ndtau) tauhu=taut2
                sca1=1.0d0+(tstep/tauhu)*(tset/tempmol-1.0d0)
                if (sca1 < 0.0d0) sca1=0.0d0
                scale=sqrt(sca1)
                if (mod(mdstep,nrep1) == 0) then
                    write (75,'(2i8,4f18.6)')mdstep,i1,tempmol,tset, &
                    tauhu*1e+15,scale
                end if
                do i3=1,3
                    vel(i3,ihu)=vel(i3,ihu)*scale
                    ekinmol=ekinmol+xmasat(ihu)*vel(i3,ihu)*vel(i3,ihu)
                end do
            end do
            ekinmol=0.50d0*ekinmol/convmd
            tempmol=2.0d0*ekinmol/(float(3*na1m)*rgasc*xjouca/1.0d3)
        end do

        ekin=zero
        do i1=1,na
            do i2=1,3
                ekin=ekin+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        ekin=0.50d0*ekin/convmd
        if (mod(mdstep,nrep1) == 0) then
            close (75)
        end if

    end if

!*********************************************************************
!                                                                    *
!     Temperature control using temperature damping                  *
!     applied to individual atoms with two different temperature     *
!     regimes                                                        *
!                                                                    *
!*********************************************************************
    if (ncons == 6) then
        ekin1=zero
        ekin2=zero
        nhu=ndtau
        nhu2=na-ndtau
        if (nhu > na) then
            nhu=na
            nhu2=0
        end if
        do i1=1,nhu
            do i2=1,3
                ekin1=ekin1+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        do i1=nhu+1,na
            do i2=1,3
                ekin2=ekin2+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        ekin1=0.50d0*ekin1/convmd
        ekin2=0.50d0*ekin2/convmd
        tempmd1=2.0d0*ekin1/(float(3*nhu)*rgasc*xjouca/1.0d3)
        tempmd2=2.0d0*ekin2/(float(3*nhu2)*rgasc*xjouca/1.0d3)
        sca1=1.0d0+(tstep/taut)*(tset/tempmd1-1.0d0)
        sca2=1.0d0+(tstep/taut2)*(tset2/tempmd2-1.0d0)
        if (sca1 < 0.0d0) sca1=0.0d0
        if (sca2 < 0.0d0) sca2=0.0d0
        scale1=sqrt(sca1)
        scale2=sqrt(sca2)
        if (scale1 > 2.0d0) scale=2.0d0
        if (scale2 > 2.0d0) scale=2.0d0
        ekin1=zero
        ekin2=zero
        do i1=1,na
            do i2=1,3
                if (i1 <= nhu) vel(i2,i1)=scale1*vel(i2,i1)
                if (i1 > nhu) vel(i2,i1)=scale2*vel(i2,i1)
                ekin=ekin+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
                if (i1 <= nhu) ekin1=ekin1+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
                if (i1 > nhu) ekin2=ekin2+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
        end do
        ekin=0.50d0*ekin/convmd
        tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
        ekin1=0.50d0*ekin1/convmd
        ekin2=0.50d0*ekin2/convmd
        tempmd1=2.0d0*ekin1/(float(3*nhu)*rgasc*xjouca/1.0d3)
        tempmd2=2.0d0*ekin2/(float(3*nhu2)*rgasc*xjouca/1.0d3)
        write (75,'(i8,12f12.4)') mdstep,tempmd,tempmd1,tempmd2, &
        tset,tset2,scale1,scale2
    end if

    return
    end subroutine tdamp
!*********************************************************************
!*********************************************************************

    subroutine tregime

!*********************************************************************
    include 'cbka.blk'
    dimension tsetcur(mtzone),ekinz(mtzone),tempmdz(mtzone), &
    scaz(mtzone),scalez(mtzone)
!*********************************************************************
!                                                                    *
!     Temperature control based on tregime.in                        *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In tregime'
        call timer(65)
        close (65)
    end if

    goto (100,500,200) ncons+1
    stop 'Tregime: unknown temperature damping method'
!*********************************************************************
!                                                                    *
!     Temperature control based on individual atoms
!                                                                    *
!*********************************************************************
    100 continue
    ekin=zero
    mdsteph=mdstep+nprevrun
    do i1=1,ntrc
        if (mdsteph >= nittc(i1)) itcur=i1
    end do
    do i1=1,nntreg(itcur)
        tsetcur(i1)=tsettreg(itcur,i1)+(mdsteph-nittc(itcur))* &
        dttreg(itcur,i1)
        if (tsetcur(i1) < 0.01d0) tsetcur(i1)=0.01d0
    end do

    do i1=1,na
        ekinat=0
        do i2=1,3
            ekinat=ekinat+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
        end do
        ekinat=0.50d0*ekinat/convmd
        tempat=2.0d0*ekinat/(3.0d0*rgasc*xjouca/1.0d3)
        ifound=0
        do i2=1,nntreg(itcur)
            if (i1 >= ia1treg(itcur,i2) .AND. i1 <= ia2treg(itcur,i2)) then
                tseth=tsetcur(i2)
                tautz=1.0d-15*tdamptreg(itcur,i2)
                ifound=1
            end if
        end do
    !     if (ifound.eq.0) then
    !     write (*,*)'Error in temperature regime:',itcur,' for atom:',i1
    !     stop 'Wrong temperature regime; see log-file'
    !     end if
        sca1=1.0d0+(tstep/tautz)*(tseth/tempat-1.0d0)
        if (ifound == 0) sca1=1.0d0
        if (sca1 < 0.0d0) sca1=0.0d0
        scale=sqrt(sca1)
        do i2=1,3
            vel(i2,i1)=vel(i2,i1)*scale
        end do
    end do
          
    do i1=1,nntreg(itcur)
        ekinz(i1)=zero
        nahu=ia2treg(itcur,i1)-ia1treg(itcur,i1)+1
        do i2=ia1treg(itcur,i1),ia2treg(itcur,i1)
            do i3=1,3
                ekinz(i1)=ekinz(i1)+xmasat(i2)*vel(i3,i2)*vel(i3,i2)
            end do
        end do
        ekinz(i1)=0.50d0*ekinz(i1)/convmd
        tempmdz(i1)=2.0d0*ekinz(i1)/(float(3*nahu)*rgasc*xjouca/1.0d3)
    end do

    do i1=1,na
        do i2=1,3
            ekin=ekin+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
        end do
    end do

    ekin=0.50d0*ekin/convmd
    tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
    if (mod(mdstep,nrep1) == 0) then
        open (75,file='fort.75',status='unknown',position='append')
        write (75,'(i8,21f8.2)') mdsteph,tempmd,(tempmdz(i1), &
        tsetcur(i1),i1=1,nntreg(itcur))
        close (75)
    end if

    return
!*********************************************************************
!                                                                    *
!     Temperature control based on system temperatures               *
!                                                                    *
!*********************************************************************
    200 continue
    do i1=1,mtzone
        ekinz(i1)=zero
    end do
    ekin=zero
    do i1=1,ntrc
        if (mdstep >= nittc(i1)) itcur=i1
    end do
    do i1=1,nntreg(itcur)
        tsetcur(i1)=tsettreg(itcur,i1)+(mdstep-nittc(itcur))* &
        dttreg(itcur,i1)
        if (tsetcur(i1) < 0.01d0) tsetcur(i1)=0.01d0
    end do
          
    do i1=1,nntreg(itcur)
        do i2=ia1treg(itcur,i1),ia2treg(itcur,i1)
            do i3=1,3
                ekinz(i1)=ekinz(i1)+xmasat(i2)*vel(i3,i2)*vel(i3,i2)
            end do
        end do

        ekinz(i1)=0.50d0*ekinz(i1)/convmd
        nahu=ia2treg(itcur,i1)-ia1treg(itcur,i1)+1
        tautz=1.0d-15*tdamptreg(itcur,i1)
        tempmdz(i1)=2.0d0*ekinz(i1)/(float(3*nahu)*rgasc*xjouca/1.0d3)
        scaz(i1)=1.0d0+(tstep/tautz)*(tsetcur(i1)/tempmdz(i1)-1.0d0)
        if (scaz(i1) < zero) scaz(i1)=zero
        scalez(i1)=sqrt(scaz(i1))
        if (scalez(i1) > 2.0d0) scalez(i1)=2.0d0

        ekinz(i1)=zero
        do i2=ia1treg(itcur,i1),ia2treg(itcur,i1)
            do i3=1,3
                vel(i3,i2)=scalez(i1)*vel(i3,i2)
                ekinz(i1)=ekinz(i1)+xmasat(i2)*vel(i3,i2)*vel(i3,i2)
                ekin=ekin+xmasat(i2)*vel(i3,i2)*vel(i3,i2)
            end do
        end do
        ekinz(i1)=0.50d0*ekinz(i1)/convmd
        tempmdz(i1)=2.0d0*ekinz(i1)/(float(3*nahu)*rgasc*xjouca/1.0d3)
    end do

    ekin=0.50d0*ekin/convmd
    tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
    if (mod(mdstep,nrep1) == 0) then
        open (75,file='fort.75',status='unknown',position='append')
        write (75,'(i8,21f8.2)') mdstep,tempmd,(tempmdz(i1), &
        tsetcur(i1),i1=1,nntreg(itcur))
        close (75)
    end if
    return

    500 continue
    stop 'Tregime: unknown temperature damping method'

    return
    end subroutine tregime
!*********************************************************************
!*********************************************************************

    subroutine vregime

!*********************************************************************
    include 'cbka.blk'
    dimension fc(nat,3)
    character(5) :: qvtyp
!*********************************************************************
!                                                                    *
!     Volume control based on vregime.in                             *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In vregime'
        call timer(65)
        close (65)
    end if

    mdsteph=mdstep+nprevrun

    do i1=1,nvrc
        if (mdsteph >= nitvc(i1)) ivcur=i1
    end do

    do i1=1,nnvreg(ivcur)

        do i2=1,na
            fc(i2,1)=c(i2,1)/tm11
            fc(i2,2)=(c(i2,2)-tm21*fc(i2,1))/tm22
            fc(i2,3)=(c(i2,3)-tm31*fc(i2,1)-tm32*fc(i2,2))/tm33
        end do

        qvtyp=qvtype(ivcur,i1)

        if (qvtyp(1:1) == 'a' .AND. qvtyp(1:4) /= 'alfa') then
            axis(1)=axis(1)+dvvreg(ivcur,i1)
        end if
        if (qvtyp(1:1) == 'b' .AND. qvtyp(1:4) /= 'beta') then
            axis(2)=axis(2)+dvvreg(ivcur,i1)
        end if
        if (qvtyp(1:1) == 'c') then
            axis(3)=axis(3)+dvvreg(ivcur,i1)
        end if
        if (qvtyp(1:4) == 'alfa') then
            angle(1)=angle(1)+dvvreg(ivcur,i1)
        end if
        if (qvtyp(1:4) == 'beta') then
            angle(2)=angle(2)+dvvreg(ivcur,i1)
        end if
        if (qvtyp(1:5) == 'gamma') then
            angle(3)=angle(3)+dvvreg(ivcur,i1)
        end if

        axiss(1)=axis(1)
        axiss(2)=axis(2)
        axiss(3)=axis(3)
        halfa=angle(1)*dgrrdn
        hbeta=angle(2)*dgrrdn
        hgamma=angle(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axis(1)*sinbet*sinphi
        tm21=axis(1)*sinbet*cosphi
        tm31=axis(1)*cosbet
        tm22=axis(2)*sinalf
        tm32=axis(2)*cosalf
        tm33=axis(3)
        volume=tm11*tm22*tm33

        if (ivsca(ivcur,i1) == 1) then  !rescale coordinates
            do i2=1,na
                c(i2,1)=fc(i2,1)*tm11
                c(i2,2)=fc(i2,1)*tm21+fc(i2,2)*tm22
                c(i2,3)=fc(i2,1)*tm31+fc(i2,2)*tm32+fc(i2,3)*tm33
            end do
        end if

    end do

    if (mod(mdstep,nrep1) == 0) then
        open (77,file='fort.77',status='unknown',position='append')
        write (77,'(i8,21f8.2)') mdsteph,tm11*tm22*tm33,axis(1),axis(2), &
        axis(3),angle(1),angle(2),angle(3)
        close (77)
    end if
    return
    end subroutine vregime
!*********************************************************************
!*********************************************************************

    subroutine eregime

!*********************************************************************
    include 'cbka.blk'
    dimension fc(nat,3)
    character(5) :: qetyp
!*********************************************************************
!                                                                    *
!     Electric field control based on eregime.in                     *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In eregime'
        call timer(65)
        close (65)
    end if

    mdsteph=mdstep+nprevrun

    do i1=1,nerc
        if (mdsteph >= nitec(i1)) iecur=i1
    end do

    ifieldx=0
    ifieldy=0
    ifieldz=0
    vfieldx=zero
    vfieldy=zero
    vfieldz=zero

    do i1=1,nnereg(iecur)
        qetyp=qetype(iecur,i1)
        if (qetyp(1:1) == 'x') then
            vfieldx=ereg(iecur,i1)
            ifieldx=1
        end if
        if (qetyp(1:1) == 'y') then
            vfieldy=ereg(iecur,i1)
            ifieldy=1
        end if
        if (qetyp(1:1) == 'z') then
            vfieldz=ereg(iecur,i1)
            ifieldz=1
        end if
    end do

    if (mod(mdstep,nrep1) == 0) then
        open (78,file='fort.78',status='unknown',position='append')
        write (78,'(i8,3f10.6,4f12.2)') mdsteph,vfieldx,vfieldy,vfieldz, &
        efi,efix,efiy,efiz
        close (78)
    end if

    return
    end subroutine eregime
!*********************************************************************
!!*********************************************************************
! df398 commented this subroutine because it conflicts with
! another subroutine in GMIN with the same name and generates
! a compilation error when issuing 'make' when building gmin
    subroutine reflect
!
!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Add reflecting wall                                            *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In reflect'
        call timer(65)
        close (65)
    end if

    do i1=1,na

        if (ireflx == 1) then     !reflection in x-direction

            if (c(i1,1) > axiss(1)) then
                diplace=c(i1,1)-axiss(1)
                c(i1,1)=axiss(1)-displace
                vel(1,i1)=-vel(1,i1)
            end if

            if (c(i1,1) < 0.0) then
                diplace=c(i1,1)
                c(i1,1)=-displace
                vel(1,i1)=-vel(1,i1)
            end if

        end if

        if (irefly == 1) then    !reflection in y-direction

            if (c(i1,2) > axiss(2)) then
                diplace=c(i1,2)-axiss(2)
                c(i1,2)=axiss(2)-displace
                vel(2,i1)=-vel(2,i1)
            end if

            if (c(i1,2) < 0.0) then
                diplace=c(i1,2)
                c(i1,2)=-displace
                vel(2,i1)=-vel(2,i1)
            end if

        end if

        if (ireflz == 1) then    !reflection in z-direction

            if (c(i1,3) > axiss(3)) then
                diplace=c(i1,3)-axiss(3)
                c(i1,3)=axiss(3)-displace
                vel(3,i1)=-vel(3,i1)
            end if

            if (c(i1,3) < 0.0) then
            !     write (64,*)i1,c(i1,3)
                diplace=c(i1,3)
                c(i1,3)=-displace
                vel(3,i1)=-vel(3,i1)
            end if

        end if

    end do

    return
    end subroutine reflect
!*********************************************************************
!*********************************************************************

    subroutine pireg

!*********************************************************************
    include 'cbka.blk'
    dimension fc(nat,3)
    character(5) :: qptyp
!*********************************************************************
!                                                                    *
!     Piston control based on piston.in                              *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In pireg'
        call timer(65)
        close (65)
    end if

    mdsteph=mdstep+nprevrun

    do i1=1,npist
        if (mdsteph >= nitpist(i1)) ipcur=i1
    end do

    nsteph=mdsteph-nitpist(ipcur)
    edeep=edeep2(ipcur)
    rdeep=rdeep2(ipcur)
    pshft=pshft2(ipcur)+speedp2(ipcur)*nsteph
    rcut=rcut2(ipcur)
    qptyp=qptype(ipdir)
    if (qptyp(1:1) == 'x') ipdir=1
    if (qptyp(1:1) == 'y') ipdir=2
    if (qptyp(1:1) == 'z') ipdir=3

    if (mod(mdstep,nrep1) == 0) then
        open (79,file='fort.79',status='unknown',position='append')
        write (79,100) mdsteph,edeep,rdeep,pshft, &
        rcut,epist
        close (79)
    end if

    return
    100 format (i8,2x,5f12.2)
    end subroutine pireg
!*********************************************************************
!*********************************************************************

    subroutine calcpres

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Pressure calculation                                           *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In calpres'
        call timer(65)
        close (65)
    end if

    deda(1)=0.0d0
    deda(2)=0.0d0
    deda(3)=0.0d0
    pressu=0.0d0
!*********************************************************************
!                                                                    *
!     Construct total force matrix                                   *
!                                                                    *
!*********************************************************************
    do i1=1,na
        do i2=1,3
            do ix=1,3
                do iy=0,2
                    do iz=0,2
                        i3=ix+3*iy+9*iz
                        d(i2,i1)=d(i2,i1)+dcell(i2,i1,i3)
                    end do
                end do
            end do
        end do
    end do
!*********************************************************************
!                                                                    *
!     Calculate internal pressure                                    *
!                                                                    *
!*********************************************************************
    dedaintx=0.0d0
    dedainty=0.0d0
    dedaintz=0.0d0
    do i1=1,na
    !     fx=c(i1,1)/tm11
    !     fy=(c(i1,2)-tm21*fx)/tm22
    !     fz=(c(i1,3)-tm31*fx-tm32*fy)/tm33

    
    ! This seems to be the correct form for triclinic, APT 2/28/2008
    
        fx=c(i1,1)/tm11
        fy=c(i1,2)/tm22
        fz=c(i1,3)/tm33

        dedaintx=dedaintx+d(1,i1)*fx
        dedainty=dedainty+d(2,i1)*fy
        dedaintz=dedaintz+d(3,i1)*fz
    end do
!*********************************************************************
!                                                                    *
!     Calculate pressure contribution of interactions across         *
!     the periodic boundaries                                        *
!                                                                    *
!*********************************************************************
    dedaextx=0.0d0
    dedaexty=0.0d0
    dedaextz=0.0d0
    do i1=1,na
        do i2=1,3
            do ix=1,3
                do iy=0,2
                    do iz=0,2
                        i3=ix+3*iy+9*iz
                        dxda=0.0d0
                        dyda=0.0d0
                        dzda=0.0d0
                    !     if (i2.eq.1.and.ix.ne.2) then
                    !     dxda=-float(ix-2)
                    !     end if
                    !     if (i2.eq.2.and.iy.ne.1) then
                    !     dyda=-float(iy-1)
                    !     end if
                    !     if (i2.eq.3.and.iz.ne.1) then
                    !     dzda=-float(iz-1)
                    !     end if

                    
                    ! This seems to be the correct form for triclinic, APT 2/28/2008
                    
                        if (i2 == 1) then
                            dxda=-float(ix-2)
                        end if
                        if (i2 == 2) then
                            dyda=-float(iy-1)-float(ix-2)*tm21/tm22
                        end if
                        if (i2 == 3) then
                            dzda=-float(iz-1)-float(ix-2)*tm31/tm33-float(iy-1)*tm32/tm33
                        end if

                        dedaxh=dxda*dcell(i2,i1,i3)
                        dedayh=dyda*dcell(i2,i1,i3)
                        dedazh=dzda*dcell(i2,i1,i3)
                        dedaextx=dedaextx+dedaxh
                        dedaexty=dedaexty+dedayh
                        dedaextz=dedaextz+dedazh
                    end do
                end do
            end do
        end do
    end do

    if (ipresm == 0) then        !Both internal and external pressure; default
        deda(1)=dedaintx+dedaextx
        deda(2)=dedainty+dedaexty
        deda(3)=dedaintz+dedaextz
    end if
    if (ipresm == 1) then       !Only external pressure
        deda(1)=dedaextx
        deda(2)=dedaexty
        deda(3)=dedaextz
    end if
    if (ipresm == 2) then       !Only internal pressure
        deda(1)=dedaintx
        deda(2)=dedainty
        deda(3)=dedaintz
    end if

    cpres=1.0e-24*avognr*xjouca
    volu=tm11*tm22*tm33
    presx=-deda(1)/(tm22*tm33*cpres)+2.d0*ekin/(3.d0*volu*cpres)
    presy=-deda(2)/(tm11*tm33*cpres)+2.d0*ekin/(3.d0*volu*cpres)
    presz=-deda(3)/(tm11*tm22*cpres)+2.d0*ekin/(3.d0*volu*cpres) !GPa
    presx=presx*1d+3
    presy=presy*1d+3
    presz=presz*1d+3 !in MPa
    pressu=(presx+presy+presz)/3.d0

    volume=tm11*tm22*tm33
    if (mod(mdstep,nrep1) == 0) then
    !*********************************************************************
    !                                                                    *
    !     Pressure output to fort.59                                     *
    !                                                                    *
    !*********************************************************************
        open (59,file='fort.59',status='unknown',position='append')
        if (ipropt /= 1) write (59,'(i8,20(1x,f12.6))')mdstep+nprevrun, &
        axis(1),axis(2),axis(3),presx/1000.0,presy/1000.0,presz/1000.0, &
        pset,volume
        if (ipropt == 1) write (59,'(i8,22(1x,f12.6))')mdstep+nprevrun, &
        axis(1),axis(2),axis(3),presx/1000.0,presy/1000.0,presz/1000.0, &
        vprestax,vprestay,vprestaz,volume
        close (59)
    end if

    return
    end subroutine calcpres
!*********************************************************************
!*********************************************************************

    subroutine pdamp

!*********************************************************************
    include 'cbka.blk'
    dimension fc(nat,3)
!*********************************************************************
!                                                                    *
!     Pressure control                                               *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In pdamp'
        call timer(65)
        close (65)
    end if


    third=one/three
!**********************************************************************
!                                                                     *
!     Conversion to fractional coordinates                            *
!                                                                     *
!**********************************************************************
    if (ipresm == 0 .OR. ipresm == 2) then     !internal pressure included
        do i1=1,na
            fc(i1,1)=c(i1,1)/tm11
            fc(i1,2)=(c(i1,2)-tm21*fc(i1,1))/tm22
            fc(i1,3)=(c(i1,3)-tm31*fc(i1,1)-tm32*fc(i1,2))/tm33
        end do
    end if
     
    h1=tstep/taup
    difg=1.01d0
    difl=1.0d0/difg
    scalex=(one+h1*(0.001*presx-pset))**third
    if (ipropt == 1) scalex=(one+h1*(0.001d0*presx-vprestax))**third
    if (scalex > difg) scalex=difg
    if (scalex < difl) scalex=difl
    if (inpt == 1 .OR. inpt == 4) scalex=1.0d0
    scaley=(one+h1*(0.001d0*presy-pset))**third
    if (ipropt == 1) scaley=(one+h1*(0.001d0*presy-vprestay))**third
    if (scaley > difg) scaley=difg
    if (scaley < difl) scaley=difl
    if (inpt == 2 .OR. inpt == 4) scaley=1.0d0
    scalez=(one+h1*(0.001d0*presz-pset))**third
    if (ipropt == 1) scalez=(one+h1*(0.001d0*presz-vprestaz))**third
    if (inpt == 3) scalez=1.0d0
    if (scalez > difg) scalez=difg
    if (scalez < difl) scalez=difl
     
    if (inpt == 5) then     !scale all cell parameters by same amount
        scasum2=scalex+scaley+scalez
        scalex=third*scasum2
        scaley=third*scasum2
        scalez=third*scasum2
    end if

    axis(1)=scalex*axis(1)
    axis(2)=scaley*axis(2)
    axis(3)=scalez*axis(3)
     
    if (ipresm == 0 .OR. ipresm == 2) then     !internal pressure included
        do i1=1,na
            fc(i1,1)=scalex*fc(i1,1)
            fc(i1,2)=scaley*fc(i1,2)
            fc(i1,3)=scalez*fc(i1,3)
        end do
    end if

    axiss(1)=axis(1)
    axiss(2)=axis(2)
    axiss(3)=axis(3)
    halfa=angle(1)*dgrrdn
    hbeta=angle(2)*dgrrdn
    hgamma=angle(3)*dgrrdn
    sinalf=sin(halfa)
    cosalf=cos(halfa)
    sinbet=sin(hbeta)
    cosbet=cos(hbeta)
    cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
    if (cosphi > 1.0d0) cosphi=1.0d0
    sinphi=sqrt(one-cosphi*cosphi)
    tm11=axis(1)*sinbet*sinphi
    tm21=axis(1)*sinbet*cosphi
    tm31=axis(1)*cosbet
    tm22=axis(2)*sinalf
    tm32=axis(2)*cosalf
    tm33=axis(3)

!*********************************************************************
!                                                                    *
!     Conversion to cartesian coordinates                            *
!                                                                    *
!*********************************************************************
    if (ipresm == 0 .OR. ipresm == 2) then     !internal pressure included
        do i1=1,na
            c(i1,1)=fc(i1,1)*tm11
            c(i1,2)=fc(i1,1)*tm21+fc(i1,2)*tm22
            c(i1,3)=fc(i1,1)*tm31+fc(i1,2)*tm32+fc(i1,3)*tm33
        end do
    end if
         
    return
    end subroutine pdamp
!*********************************************************************
!*********************************************************************

    subroutine vernh1

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Velocity Verlet+Nose-Hoover thermostat part 1                  *
!                                                                    *
!     see Frenkel and Smit, Understanding Molecular Simulation       *
!     Academic Press Inc. San Diego, 1996 pp.388-391                 *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In vernh1'
        call timer(65)
        close (65)
    end if
    vsq=zero
    hulp=0.001d0*convmd*rgasc*xjouca
    do i1=1,na
        do i2=1,3
            fi=accel(i2,i1)-xinh*vel(i2,i1)
            c(i1,i2)=c(i1,i2)+tstep*(vel(i2,i1)+0.50d0*tstep*fi)
            vsq=vsq+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            vel(i2,i1)=vel(i2,i1)+0.50d0*tstep*fi
        end do
    end do
    fsnh=(vsq-hulp*ngnh*tset)/vqnh
    snh=snh+tstep*(xinh+0.50d0*tstep*fsnh)
    xinh=xinh+0.50d0*tstep*fsnh
    return
    end subroutine vernh1
!*********************************************************************
!*********************************************************************

    subroutine vernh2

!*********************************************************************
    include 'cbka.blk'
    dimension veln(3,nat),velo(3,nat),hnh(3,nat)
    logical :: ready
!*********************************************************************
!                                                                    *
!     Velocity Verlet+Nose-Hoover thermostat part 2                  *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In vernh2'
        call timer(65)
        close (65)
    end if
    hulp=0.001d0*convmd*rgasc*xjouca
    vsq=zero
    do i1=1,na
        do i2=1,3
            accel(i2,i1)=-convmd*d(i2,i1)/xmasat(i1)
            veln(i2,i1)=vel(i2,i1)
            vsq=vsq+xmasat(i1)*veln(i2,i1)*veln(i2,i1)
        end do
    end do
    xin=xinh
    ready= .FALSE. 
    iit=0
    10 if ( .NOT. ready) then
        iit=iit+1
        xio=xin
        delxi=zero
        do i1=1,na
            do i2=1,3
                velo(i2,i1)=veln(i2,i1)
                hnh(i2,i1)=velo(i2,i1)-vel(i2,i1)- &
                (accel(i2,i1)-xio*velo(i2,i1))*0.50d0*tstep
                binh=velo(i2,i1)*tstep/vqnh
                delxi=delxi-hnh(i2,i1)*binh
            end do
        end do
        dnh=-0.50d0*xio*tstep-1.0d0
        hnho=xio-xin+dnh*(-vsq+hulp*ngnh*tset)*tstep/(2.0d0*vqnh)
        cibi=-vsq*tstep*tstep/(2.0d0*vqnh)
        delxi=(delxi+hnho*dnh)/(dnh-cibi)
        xin=xio+delxi
        vsq=zero
        do i1=1,na
            do i2=1,3
                ci=velo(i2,i1)*0.50d0*tstep
                veln(i2,i1)=velo(i2,i1)+(hnh(i2,i1)-ci*delxi)/dnh
                vsq=vsq+xmasat(i1)*veln(i2,i1)*veln(i2,i1)
            end do
        end do
        i=0
        ready= .TRUE. 
        20 if (ready .AND. i <= na) then
            i=i+1
            if (i <= na) then
                do i2=1,3
                    if (abs((veln(i2,i)-velo(i2,i))/veln(i2,i)) > errnh) &
                    ready= .FALSE. 
                end do
            else
                if (abs((xin-xio)/xin) > errnh) ready= .FALSE. 
            end if
            goto 20
        end if
        goto 10
    end if

    write (67,*)mdstep,iit
    do i1=1,na
        do i2=1,3
            vel(i2,i1)=veln(i2,i1)
        end do
    end do
    xinh=xin
    xin2=xinh/convmd
    vsq2=vsq/convmd
    hulp2=0.001d0*rgasc*xjouca
    ham=estrc+0.50d0*vsq2+ &
    & 0.50d0*xin2*xin2*vqnh+hulp2*ngnh*tset*snh
    write (67,'(6d20.10)')ham,estrc,vsq2,xin2,snh,xinh
    ekin=0.50d0*vsq/convmd
    tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
    return
    end subroutine vernh2
!*********************************************************************
!*********************************************************************

    subroutine chargess

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
    dimension elcvec(nat),char(nat)
!*********************************************************************
!                                                                    *
!     Determine chargess on atoms: separate molecules                 *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In chargess'
        call timer(65)
        close (65)
    end if
    if (ncha == 5) return  !fixed chargess

    if (nmcharge > 0) then  !fixed molecular chargess
        if (ioldchg == 0) then
            do i1=1,na
                ch(i1)=zero
                chgold(i1)=zero
            end do
            call charmol
        else
            call charmol
        end if
        return
    end if

    if (kx > 0 .OR. ky > 0 .OR. kz > 0) then
        call charcrys
        return
    end if

    if (ncha == 4) then
        call chargesst    !chargess on full system
        return
    end if

    if (ncha == 0) then
        call chargesstnocg    !chargess on full system without CG solver
        return
    end if

    if (ncha == 6) then
        call charfrac    !Takes molecular fragmentation into account
        return
    end if

    ech=zero

!     swbo=swb
!     swb=5.0
!     call tap7th

    do i1=1,nmolo5
        na1m=nmolat2(i1,1)

        do i2=1,na1m+1
            do i3=1,na1m+1
                xmortr(i2,i3)=zero
            end do
        end do
        elnor=zero
        do i2=1,na1m
            ihu=nmolat2(i1,i2+1)
            elnor=elnor+stlp(ia(ihu,1))
        end do
        chamol=elnor-elmol2(i1)
        do i2=1,na1m
            ihu=nmolat2(i1,i2+1)
            xmortr(i2,i2)=2.0d0*eta(ia(ihu,1))  !Mortier method
            elcvec(i2)=-chi(ia(ihu,1))        !Mortier method
        !     xmortr(i2,i2)=vjqeq(ia(ihu,1))     !Qeq method
        !     elcvec(i2)=-chiqeq(ia(ihu,1))     !Qeq method
            xmortr(na1m+1,i2)=1.0d0
        end do
        elcvec(na1m+1)=chamol     !Charge on molecule
        do i2=1,na1m-1
            ihu2=nmolat2(i1,i2+1)
            ity1=ia(ihu2,1)
            do i3=i2+1,na1m
                ihu3=nmolat2(i1,i3+1)
                call dista2(ihu2,ihu3,dis,dx,dy,dz)
                if (dis*dis < swb*swb) then
                    call taper(dis,dis*dis)
                    ity2=ia(ihu3,1)
                !     gamt=sqrt(gam(ity1)*gam(ity2))
                !     sig=1.0/(rqeq(ity1)+rqeq(ity2))
                !     sig=0.75/(rqeq(ity1)+rqeq(ity2))
                !     xmortr(i3,i2)=sw*14.40*exp(-sig*dis)  !Qeq method
                    hulp1=(dis**3+gamcco(ity1,ity2))
                    hulp2=hulp1**(1.0d0/3.0d0)
                    if (ncha == 1) xmortr(i3,i2)=sw*14.40d0/dis     !Mortier method
                    if (ncha == 3) xmortr(i3,i2)=sw*14.40d0/hulp2   !Shielded Mortier method
                end if
            end do
        end do
        call matsym4(na1m+1,na1m+1,na1m+1,na1m+1,na1m+1,xmortr,char, &
        elcvec)

    !     write (65,'(a40,i6,f12.4)')qmol,na1m+1,char(1)
    !     do i2=1,na1m+1
    !     write (65,'(40f8.2)')(xmortr(i3,i2),i3=1,na1m+1)
    !     end do

        do i2=1,na1m
            ihu=nmolat2(i1,i2+1)
            ch(ihu)=char(i2)
            ech=ech+23.02d0*(chi(ia(ihu,1))*ch(ihu)+ &
            eta(ia(ihu,1))*ch(ihu)*ch(ihu))
        end do

      END DO
!     swb=swbo
!     call tap7th
    return
    end subroutine chargess
!*********************************************************************
!*********************************************************************

    subroutine charmol

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
    dimension elcvec(nat),char(nat)
!     Space for dsysv
!     integer ipiv(nat)
!     real*8 work(5*nat)
!*********************************************************************
!                                                                    *
!     Determine chargess on atoms: separate molecules with fixed      *
!     molecular chargess                                              *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In charmol'
        call timer(65)
        close (65)
    end if

    ech=zero
    cfield=332.0638d0
     
    nsum=0

    do i2=1,na+nmcharge
        elcvec(i2)=zero
        do i3=1,na+nmcharge
            xmortr(i2,i3)=zero
        end do
    end do

    do  i1=1,nmcharge
        nastart=iat1mc(i1)
        naend=iat2mc(i1)
        nhulp=naend-nastart+1
        nsum=nsum+nhulp
         
        do i2=1,nhulp                  !Fill in diagonal
            ihu=nastart+i2-1
            xmortr(ihu,ihu)=2.0*eta(ia(ihu,1))
            elcvec(ihu)=-chi(ia(ihu,1))
            xmortr(na+i1,ihu)=1.0d0
! df398 commented next block. cfield multiplication is errorneous. no conversion to kcal/mol/A is needed
! since the electric field should already contain the 1/(4*pi*epsilon_0) (=14.40 eV*A/e^2) factor when specified in eV/A. Second, the EEM system of
! linear equations is solved in the eV units, which means that no conversion from eV to kcal/mol should be done for the EEM right-hand side term due to
! the external electric field. see bug-fix r60538 in ADF website: https://www.scm.com/support/documentation-previous-versions/changelog-bug-fixes-adf-modeling-suite-2017/
!            if (ifieldx == 1 .AND. ipolar == 1) then
!                elcvec(ihu)=elcvec(ihu)-vfieldx*cfield*c(ihu,1)
!            end if
!            if (ifieldy == 1 .AND. ipolar == 1) then
!                elcvec(ihu)=elcvec(ihu)-vfieldy*cfield*c(ihu,2)
!            end if
!            if (ifieldz == 1 .AND. ipolar == 1) then
!                elcvec(ihu)=elcvec(ihu)-vfieldz*cfield*c(ihu,3)
!            end if
            if (ifieldx == 1 .AND. ipolar == 1) then
                elcvec(ihu)=elcvec(ihu)-vfieldx*c(ihu,1)
            end if
            if (ifieldy == 1 .AND. ipolar == 1) then
                elcvec(ihu)=elcvec(ihu)-vfieldy*c(ihu,2)
            end if
            if (ifieldz == 1 .AND. ipolar == 1) then
                elcvec(ihu)=elcvec(ihu)-vfieldz*c(ihu,3)
            end if
        end do
        elcvec(na+i1)=vmcha(i1)     !Charge on molecule

        do i2=1,nhulp-1                   !Fill in matrix for each molecule
            ihu2=nastart+i2-1
            ity1=ia(ihu2,1)
            do i3=i2+1,nhulp
                ihu3=nastart+i3-1
                call dista2(ihu2,ihu3,dis,dx,dy,dz)
                if (dis*dis < swb*swb) then
                    call taper(dis,dis*dis)
                    ity2=ia(ihu3,1)
                    hulp1=(dis**3+gamcco(ity1,ity2))
                    hulp2=hulp1**(1.0d0/3.0d0)
                    xmortr(ihu3,ihu2)=sw*14.40d0/hulp2
                end if
            end do
        end do

        do i2=1,nhulp                      !Add polarization by system
            ihu2=nastart+i2-1
            ity1=ia(ihu2,1)
            do i3=ihu2+1,na
                call dista2(ihu2,i3,dis,dx,dy,dz)
                if (dis*dis < swb*swb) then
                    call taper(dis,dis*dis)
                    ity2=ia(i3,1)
                    hulp1=(dis**3+gamcco(ity1,ity2))
                    hulp2=hulp1**(1.0d0/3.0d0)
                    xmortr(i3,ihu2)=sw*14.40d0/hulp2
                end if
            end do
        end do

    END DO
    call matsym4(na+nmcharge,neem,neem,neem,neem,xmortr,char,elcvec)
!     if (ioldchg .ne. 0) then
!        call fillmatrix(na+nmcharge,xmortr,neem) ! fill in top half of matrix
!        call dcopy(na+nmcharge,chgold,1,char,1)
!        call cgsolve(na+nmcharge,xmortr,neem,char,elcvec,mdstep,convg)
!     else
!     call matsym4(na+nmcharge,neem,neem,neem,neem,xmortr,char,elcvec)
!        call dcopy(na+nmcharge,elcvec,1,char,1)
!        call dsysv('L',na+nmcharge,1,xmortr,neem,ipiv,
!    $        char,nat,work,5*nat,info)
!        if (info .ne. 0) stop 'dsysv solved failed in charmol'
!        ioldchg=1              ! Only called once
!     endif
!     call dcopy(na+nmcharge,char,1,chgold,1) ! Save the chargess
     
    do i2=1,na
        ch(i2)=char(i2)
        ech=ech+23.02d0*(chi(ia(i2,1))*ch(i2)+ &
        eta(ia(i2,1))*ch(i2)*ch(i2))
    end do
     
    ioldchg=1
    if (nsum /= na) then
        write (*,*)nsum,na
        write (*,*)'Not all atoms defined in MOLCHARGE-keywords'
        stop 'Not all atoms defined in MOLCHARGE-keywords'
    end if

    return
    end subroutine charmol
!*********************************************************************
!*********************************************************************

    subroutine charcrys

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
    dimension elcvec(nat),char(nat)
    dimension ct(nat,3)

!     Space to hold old chargess
!     real*8 chgold(nat)        ! The chargess from the previous iteration
!     integer ioldchg           ! Flag that chgold has something interesting
!     common /cchgold/ ioldchg,chgold
!     data ioldchg /0/

!     Space for dsysv
    integer :: ipiv(nat)
    real*8 :: work(5*nat)

!*********************************************************************
!                                                                    *
!     Determine chargess on atoms: crystal                            *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In charcrys'
        call timer(65)
        close (65)
    end if
    ech=zero
    cfield=332.0638d0
    chmax=2.5d0
    chmin=-2.0d0
          
    kxt=kx
    kyt=ky
    kzt=kz
    if (kxt == 0) kxt=1
    if (kyt == 0) kyt=1
    if (kzt == 0) kzt=1

    do i1=1,na+1
        elcvec(i1)=zero
        do i2=1,na+1
            xmortr(i1,i2)=zero
        end do
    end do
    chamol=syscha       !System charge

    do i1=1,na
        ity=ia(i1,1)
        xmortr(i1,i1)=2.0d0*eta(ity)  !EEM method
        elcvec(i1)=-chi(ity)
        xmortr(na+1,i1)=1.0d0
! df398 commented next block. cfield multiplication is errorneous. no conversion to kcal/mol/A is needed
! since the electric field should already contain the 1/(4*pi*epsilon_0) (=14.40 eV*A/e^2) factor when specified in eV/A. Second, the EEM system of
! linear equations is solved in the eV units, which means that no conversion from eV to kcal/mol should be done for the EEM right-hand side term due to
! the external electric field. see bug-fix r60538 in ADF website: https://www.scm.com/support/documentation-previous-versions/changelog-bug-fixes-adf-modeling-suite-2017/
!        if (ifieldx == 1 .AND. ipolar == 1) then
!            elcvec(i1)=elcvec(i1)-vfieldx*cfield*c(i1,1)
!        end if
!        if (ifieldy == 1 .AND. ipolar == 1) then
!            elcvec(i1)=elcvec(i1)-vfieldy*cfield*c(i1,2)
!        end if
!        if (ifieldz == 1 .AND. ipolar == 1) then
!            elcvec(i1)=elcvec(i1)-vfieldz*cfield*c(i1,3)
!        end if
        if (ifieldx == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldx*c(i1,1)
        end if
        if (ifieldy == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldy*c(i1,2)
        end if
        if (ifieldz == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldz*c(i1,3)
        end if

        do k1=-kxt,kxt              !Sum over periodic self-images
            do k2=-kyt,kyt
                do k3=-kzt,kzt
                          
                    dx=k1*tm11
                    dy=k1*tm21+k2*tm22
                    dz=k1*tm31+k2*tm32+k3*tm33
                    dis=sqrt(dx*dx+dy*dy+dz*dz)

                    if (dis < swb .AND. dis > 0.001d0) then
                        call taper(dis,dis*dis)
                        hulp1=(dis**3+(1.0d0/(gam(ity)**3)))
                        hulp2=hulp1**(1.0d0/3.0d0)
                        xmortr(i1,i1)=xmortr(i1,i1)+1.0d0*sw*14.40d0/hulp2
                    end if
                ENDDO
           END DO
        END DO
    end do

    elcvec(na+1)=chamol     !Charge on system

    do  ivl=1,nvpair-nvlself      !Verlet-list
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

        if (dis2 < swb .AND. dis2 > 0.001d0) then
            call taper(dis2,dis2*dis2)
            hulp1=(dis2**3+(1.0d0/(gamt**3)))
            hulp2=hulp1**(1.0d0/3.0d0)
            xmortr(i3,i2)=xmortr(i3,i2)+sw*14.40d0/hulp2

        end if

    END DO

!     do i2=1,na-1                   !No Verlet-list
!     do i3=i2+1,na
!     ity1=ia(i2,1)
!     ity2=ia(i3,1)
!     gamt=sqrt(gam(ity1)*gam(ity2))
!     dx=c(i2,1)-c(i3,1)
!     dy=c(i2,2)-c(i3,2)
!     dz=c(i2,3)-c(i3,3)

!     do 10 k1=-kxt,kxt   !Sum over periodic images
!     dx2=dx+k1*tm11
!     do 10 k2=-kyt,kyt
!     dy2=dy+k1*tm21+k2*tm22
!     do 10 k3=-kzt,kzt
!     dz2=dz+k1*tm31+k2*tm32+k3*tm33
!     dis2=sqrt(dx2*dx2+dy2*dy2+dz2*dz2)

!     if (dis2.lt.swb.and.dis2.gt.0.001) then
!     call taper(dis2,dis2*dis2)
!     hulp1=(dis2**3+(1.0/(gamt**3)))
!     hulp2=hulp1**(1.0/3.0)
!     xmortr(i3,i2)=xmortr(i3,i2)+sw*14.40/hulp2
!     end if

!  10 continue

!     end do
!     end do

!     call matsym4(na+1,na+1,na+1,na+1,na+1,xmortr,char,
!    $elcvec)

!     if (ioldchg .ne. 0) then

! ---------- next block commented df39 --------------
!!    call fillmatrix(na+1,xmortr,neem) ! fill in top half of matrix
!!    call dcopy(na+1,chgold,1,char,1)
!!    if (ndebug == 1) then
!!        open (65,file='fort.65',status='unknown',position='append')
!!        write (65,*) 'Calling cgsolve'
!!        call timer(65)
!!        close (65)
!!    end if
!!    call cgsolve(na+1,xmortr,neem,char,elcvec,mdstep,convg)
!! -------------------------------------------------------

    call matsym4(na+1,na+1,na+1,na+1,na+1,xmortr,char,elcvec) !! df398 inserted (direct solver))

!     else
!        call matsym4(na+1,na+1,na+1,na+1,na+1,xmortr,char,
!    $   elcvec)
!        call dcopy(na+1,elcvec,1,char,1)
!        call dsysv('L',na+1,1,xmortr,neem,ipiv,
!    $        char,nat,work,5*nat,info)
!        if (info .ne. 0) stop 'dsysv solved failed in charcrys'
!        ioldchg=1              ! Only called once
!     endif


    ech=zero
    do i2=1,na
        ch(i2)=char(i2)
    !     if (ch(i2).gt.chmax) ch(i2)=chmax
    !     if (ch(i2).lt.chmin) ch(i2)=chmin
        ech=ech+23.02d0*(chi(ia(i2,1))*ch(i2)+ &
        eta(ia(i2,1))*ch(i2)*ch(i2))
    end do
    chisys=char(na+1)   !! df398 adapted from chargestnocg
 
    !call dcopy(nat,ch,1,chgold,1) ! Save the chargess !! df398 commented out

    return
    end subroutine charcrys
!*********************************************************************
!*********************************************************************

    subroutine chargessm (namo,nmohu,charge,chag)

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
    dimension elcvec(nat),chag(nat)
    dimension nmohu(nat)
!*********************************************************************
!                                                                    *
!     Determine chargess on a single molecule                         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In chargessm'
        call timer(65)
        close (65)
    end if
    do i1=1,namo+1
        do i2=1,namo+1
            xmortr(i1,i2)=zero
        end do
    end do

    do i2=1,namo
        ihu=nmohu(i2)
        ity=ia(ihu,1)
        xmortr(i2,i2)=2.0d0*eta(ity)  !Mortier method
        elcvec(i2)=-chi(ity)                              !Mortier method
        xmortr(namo+1,i2)=1.0d0
        xmortr(i2,namo+1)=1.0d0
    end do
    xmortr(namo+1,namo+1)=zero
    elcvec(namo+1)=charge     !Charge on molecule

    do i2=1,namo-1
        ihu2=nmohu(i2)
        ity1=ia(ihu2,1)
        do i3=i2+1,namo
            ihu3=nmohu(i3)
            ity2=ia(ihu3,1)
            call dista2 (ihu2,ihu3,dis,dx,dy,dz)

            if (dis*dis < swb*swb) then
                call taper(dis,dis*dis)
                gamt=sqrt(gam(ity1)*gam(ity2))
                hulp1=(dis**3+(1.0d0/(gamt**3)))
                hulp2=hulp1**(1.0d0/3.0d0)

            !     if (ncha.eq.1) xmortr(i3,i2)=sw*14.40/dis     !Mortier method

            !     if (ncha.eq.3) then   !Shielded EEM
                xmortr(i3,i2)=sw*14.40d0/hulp2
            !     end if

            end if

        end do
    end do

    call matsym4(namo+1,namo+1,namo+1,namo+1,namo+1,xmortr,chag, &
    elcvec)

    echh=zero
    do i2=1,namo
        ihu=nmohu(i2)
        echh=echh+23.02d0*(chi(ia(ihu,1))*chag(ihu)+ &
        eta(ia(ihu,1))*chag(ihu)*chag(ihu))
    end do

    chimol=-chag(namo+1)

    return
    end subroutine chargessm

!*********************************************************************
!*********************************************************************

    subroutine charfrac

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
!*********************************************************************
!                                                                    *
!     Determine chargess on atoms including molecular fragmentation   *
!     effects                                                        *
!                                                                    *
!*********************************************************************
    dimension chadevs(nat,20),velvec(nat),char(nat),charfull(nat), &
    chacs(nat),chi2(nat),eta2(nat),chag1(nat),c8(nat), &
    chag2(nat),bofrag(20),a(20),b(20),cf(20)

    dimension iam(nat,mbond+3),nmohu(nat),nmohu1(nat),nmohu2(nat), &
    icsort(nat)
    logical :: done(nat)
    chisyswat=-5.3539d0
!*********************************************************************
!                                                                    *
!     Find molecular fragments within molecules                      *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In charfrac'
        call timer(65)
        close (65)
    end if
    if (axiss(1) > zero) then
        write (*,*) 'Do not use fragment EEM with periodic systems'
        stop 'Do not use fragment EEM with periodic systems'
    end if
!*********************************************************************
!                                                                    *
!     Calculate full system charge distribution                      *
!                                                                    *
!*********************************************************************
    call chargesst
    do i1=1,na
        charfull(i1)=ch(i1)
    end do
    ech=zero

    do 10 i1=1,nmolo5

        chmole=zero
        elmole=zero
        nfrag=1
        na1m=nmolat2(i1,1)
        do i2=1,na1m
            nmohu(i2)=nmolat2(i1,1+i2)
            chmole=chmole+charfull(nmohu(i2))
            elmole=elmole+stlp(ia(nmohu(i2),1))
        end do

        if (2*int(elmole/2) == int(elmole)) then  !Check for odd-electron molecules
            vcharge=zero
        else
            if (chisys < chisyswat) vcharge=-1.0d0
            if (chisys > chisyswat) vcharge=+1.0d0
        end if
    !*********************************************************************
    !                                                                    *
    !     Calculate charge distribution unfragmented molecule            *
    !                                                                    *
    !*********************************************************************
        call chargessm (na1m,nmohu,vcharge,chag1)
        do i2=1,na1m
            chadevs(nmohu(i2),1)=chag1(i2)
        end do

        do i2=1,na1m-1
            do i3=i2+1,na1m
                j1=nmolat2(i1,i2+1)
                j2=nmolat2(i1,i3+1)
                nhu=nubon(j1,j2)
            !*********************************************************************
            !                                                                    *
            !     Identify possible fragmentations in molecule                   *
            !                                                                    *
            !*********************************************************************
                if (nhu /= 0 .AND. bo(nhu) < 0.75d0) then

                !*********************************************************************
                !                                                                    *
                !     Create temporary connection table (iam)                        *
                !                                                                    *
                !*********************************************************************
                    do i4=1,na1m
                        do i5=1,3+mbond
                            iam (i4,i5)=0
                        end do
                    end do
                    do i4=1,na1m-1
                        ih=nmolat2(i1,1+i4)
                        iam(i4,1)=ih
                        iam(i4,3+mbond)=0
                        do i5=i4,na1m
                            ih2=nmolat2(i1,1+i5)
                            nhu2=nubon(ih,ih2)
                        !     if (nhu2.ne.0.and.
                        !    $(ih.ne.j1.or.ih2.ne.j2).and.(ih.ne.j2.or.ih2.ne.j1)) then
                            if (nhu2 /= 0 .AND. bo(nhu2) > cutof3 .AND. &
                            (ih /= j1 .OR. ih2 /= j2) .AND. (ih /= j2 .OR. ih2 /= j1)) then
                                iam(i4,2)=iam(i4,2)+1
                                iam(i5,2)=iam(i5,2)+1
                                iam(i4,(iam(i4,2)+2))=i5
                                iam(i5,(iam(i5,2)+2))=i4
                            end if

                        end do
                    end do
                    iam(na1m,1)=nmolat2(i1,1+na1m)

                    nmoln=0
                !     write (65,*)j1,j2,nfrag+1,bo(nhu)
                !     do i4=1,na1m
                !     write (65,'(20i4)')(iam(i4,i5),i5=1,10)
                !     end do
                    nmoln=0
                !     write (65,*)na1m,nmoln
                    call sortmol(iam,na1m,nmoln)
                !     write (65,*)
                !     do i4=1,na1m
                !     write (65,'(20i4)')i4,(iam(i4,2+i5),i5=1,nsbmax),
                !    $iam(i4,3+mbond)
                !     end do
                !     write (65,*)na1m,nmoln
                !*********************************************************************
                !                                                                    *
                !     Calculate charge distributions in fragments                    *
                !                                                                    *
                !*********************************************************************
                    if (nmoln > 1) then
                        nfrag=nfrag+1
                        bofrag(nfrag)=bo(nhu)

                        ihu1=0
                        ihu2=0
                        chmole1=zero
                        elmole1=zero
                        chmole2=zero
                        elmole2=zero

                        do i4=1,na1m
                            if (iam(i4,3+mbond) == 1) then
                                ihu1=ihu1+1
                                nmohu1(ihu1)=iam(i4,1)
                                chmole1=chmole1+charfull(nmohu1(ihu1))
                                elmole1=elmole1+stlp(ia(nmohu1(ihu1),1))
                            end if
                            if (iam(i4,3+mbond) == 2) then
                                ihu2=ihu2+1
                                nmohu2(ihu2)=iam(i4,1)
                                chmole2=chmole2+charfull(nmohu2(ihu2))
                                elmole2=elmole2+stlp(ia(nmohu2(ihu2),1))
                            end if
                        end do

                        if (ihu1+ihu2 /= na1m) write (*,*) qmol,'Error in valec'
                        if (2*int(elmole1/2) == int(elmole1)) then  !Check for odd-electron molecules
                            chatot1=zero
                        else
                            if (chmole1 < zero) chatot1=-1.0d0
                            if (chmole1 >= zero) chatot1=1.0d0
                        end if

                        if (2*int(elmole2/2) == int(elmole2)) then
                            chatot2=zero
                        else
                            if (chmole2 < zero) chatot2=-1.0d0
                            if (chmole2 >= zero) chatot2=1.0d0
                        end if

                        call chargessm (ihu1,nmohu1,chatot1,chag1)
                        call chargessm (ihu2,nmohu2,chatot2,chag2)
                        do i4=1,ihu1
                            chadevs(nmohu1(i4),nfrag)=chag1(i4)
                        end do
                        do i4=1,ihu2
                            chadevs(nmohu2(i4),nfrag)=chag2(i4)
                        end do

                    end if

                end if

            end do
        end do
    !*********************************************************************
    !                                                                    *
    !     Determine final charge distribution                            *
    !                                                                    *
    !*********************************************************************

        do i2=1,nfrag-1
        !     a(i2)=1.0-bofrag(i2+1)
        !     b(i2)=bofrag(i2+1)
            a(i2)=1.0d0-sqrt(bofrag(i2+1))
            b(i2)=sqrt(bofrag(i2+1))
        end do
        do i2=1,nfrag-1
            xvelec(i2,1)=a(i2)
            xvelec(i2,1+i2)=-b(i2)
            velvec(i2)=1.0d0
            xvelec(nfrag,i2)=1.0d0
        end do
        xvelec(nfrag,nfrag)=1.0d0
        velvec(nfrag)=1.0d0
    !*********************************************************************
    !                                                                    *
    !     Calculate transposed matrix                                    *
    !                                                                    *
    !*********************************************************************
        do i2=1,nfrag
            do i3=1,nfrag
                xvelect(i2,i3)=xvelec(i3,i2)
            end do
        end do
    !*********************************************************************
    !                                                                    *
    !     Calculate transposed matrix*matrix                             *
    !                                                                    *
    !*********************************************************************
        do i2=1,nfrag
            do i3=1,nfrag
                sum=0.0d0
                do i4=1,nfrag
                    sum=sum+xvelect(i3,i4)*xvelec(i4,i2)
                end do
                xvt(i3,i2)=sum
            end do
        end do

        !df398 NAG xvt complaint call matsym4(nfrag,nfrag,nfrag,nfrag,nfrag,xvt,cf,velvec)
    !*********************************************************************
    !                                                                    *
    !     Multiply fragment charge distributions by coefficients         *
    !     to get final charge distribution                               *
    !                                                                    *
    !*********************************************************************
        echh=zero
        do i2=1,na1m
            ihu=nmohu(i2)
            ity=ia(ihu,1)
            charf=zero
            do i3=1,nfrag
                charf=charf+cf(i3)*chadevs(ihu,i3)
            end do
            ch(ihu)=charf
            echh=echh+23.02d0* &
            (chi(ity)*charf+eta(ity)*charf*charf)
        end do
          
        ech=ech+echh

    10 END DO


    return
    end subroutine charfrac
!*********************************************************************
!*********************************************************************

    subroutine chargesstnocg

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
    dimension elcvec(nat),char(nat)
!*********************************************************************
!                                                                    *
!     Determine chargess on atoms: full system exact solution         *
!     taken from vibra.f90                                           *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In chargesstnocg'
        call timer(65)
        close (65)
    end if

    ech=zero
    cfield=332.0638d0

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
        xmortr(i1,i1)=2.0d0*eta(ia(i1,1))  !Mortier method
        elcvec(i1)=-chi(ia(i1,1))        !Mortier method
        ! df398 added here a block to include external electric field effect
        ! into EEM if provided in eregime.in file (from the chargesst subroutine)
        ! -------
        if (ifieldx == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldx*c(i1,1)
        end if
        if (ifieldy == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldy*c(i1,2)
        end if
        if (ifieldz == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldz*c(i1,3)
        end if   
        ! -----         
        xmortr(na+1,i1)=1.0d0
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

        if (dis2 < swb .AND. dis2 > 0.001d0) then
            call taper(dis2,dis2*dis2)
            hulp1=(dis2**3+(1.0d0/(gamt**3)))
            hulp2=hulp1**(1.0d0/3.0d0)
            xmortr(i3,i2)=xmortr(i3,i2)+sw*14.40d0/hulp2

        end if

    10 END DO

    call matsym4(na+1,na+1,na+1,na+1,na+1,xmortr,char,elcvec)

    ech=zero
    do i2=1,na
        ch(i2)=char(i2)
        ech=ech+23.02d0*(chi(ia(i2,1))*ch(i2)+ &
        eta(ia(i2,1))*ch(i2)*ch(i2))
    end do
    chisys=char(na+1)

    return
    end subroutine chargesstnocg
!*********************************************************************
!*********************************************************************

    subroutine chargesst

!*********************************************************************
    include 'cbka.blk'
    include 'cbkm.blk'
    dimension elcvec(nat),char(nat)
!     Space for dsysv
    integer :: ipiv(nat)
    real*8 :: work(5*nat)
!*********************************************************************
!                                                                    *
!     Determine chargess on atoms: full system                        *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In chargesst'
        call timer(65)
        close (65)
    end if

    ech=zero
    cfield=332.0638d0

    if (nmethod /= 0) then
        do i2=1,na+1      !Zero EEM-matrix
            do i3=1,na+1
                xmortr(i3,i2)=zero
            end do
        end do
    end if

    do i2=1,na+1
        elcvec(i2)=zero
    end do
    chamol=syscha    !System charge

    do i1=1,na
        xmortr(i1,i1)=2.0d0*eta(ia(i1,1))  !Mortier method
        elcvec(i1)=-chi(ia(i1,1))        !Mortier method
    !     write (64,*)i1,xmortr(i1,i1),elcvec(i1)

! df398 commented next block. cfield multiplication is errorneous. no conversion to kcal/mol/A is needed
! since the electric field should already contain the 1/(4*pi*epsilon_0) (=14.40 eV*A/e^2) factor when specified in eV/A. Second, the EEM system of
! linear equations is solved in the eV units, which means that no conversion from eV to kcal/mol should be done for the EEM right-hand side term due to
! the external electric field. see bug-fix r60538 in ADF website: https://www.scm.com/support/documentation-previous-versions/changelog-bug-fixes-adf-modeling-suite-2017/
! 
!        if (ifieldx == 1 .AND. ipolar == 1) then
!            elcvec(i1)=elcvec(i1)-vfieldx*cfield*c(i1,1)
!        end if
!        if (ifieldy == 1 .AND. ipolar == 1) then
!            elcvec(i1)=elcvec(i1)-vfieldy*cfield*c(i1,2)
!        end if
!        if (ifieldz == 1 .AND. ipolar == 1) then
!            elcvec(i1)=elcvec(i1)-vfieldz*cfield*c(i1,3)
!        end if
        if (ifieldx == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldx*c(i1,1)
        end if
        if (ifieldy == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldy*c(i1,2)
        end if
        if (ifieldz == 1 .AND. ipolar == 1) then
            elcvec(i1)=elcvec(i1)-vfieldz*c(i1,3)
        end if
              
        xmortr(na+1,i1)=1.0d0
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

        if (dis2 < swb .AND. dis2 > 0.001d0) then
            call taper(dis2,dis2*dis2)
            hulp1=(dis2**3+(1.0d0/(gamt**3)))
            hulp2=hulp1**(1.0d0/3.0d0)
            xmortr(i3,i2)=xmortr(i3,i2)+sw*14.40d0/hulp2
        !     write (64,*)i2,i3,xmortr(i3,i2)

        end if

    10 END DO
    do i1=1,na+1
    !     write (64,*)i1,(xmortr(i1,i2),i2=1,na)
    end do

    call fillmatrix(na+1,xmortr,neem) ! fill in top half of matrix
    call dcopy(na+1,chgold,1,char,1)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In cgsolve'
        call timer(65)
        close (65)
    end if
    call cgsolve(na+1,xmortr,neem,char,elcvec,mdstep,convg)
    ioldchg=1              ! Only called once

    ech=zero
    do i2=1,na
        ch(i2)=char(i2)
        ech=ech+23.02d0*(chi(ia(i2,1))*ch(i2)+ &
        eta(ia(i2,1))*ch(i2)*ch(i2))
    end do
    chisys=char(na+1)

    call dcopy(nat,ch,1,chgold,1) ! Save the chargess

    do 20 ivl=1,nvpair-nvlself      !zero xmortr-matrix
        i2=nvl1(ivl)
        i3=nvl2(ivl)
        xmortr(i3,i2)=zero
    20 END DO

    return
    end subroutine chargesst
!*********************************************************************
!*********************************************************************

    subroutine ffinpt

!*********************************************************************
    include 'cbka.blk'
    dimension rcore2(nsort),ecore2(nsort),acore2(nsort)
! df398 08.10.18 added a new file unit (fort.1)) to store various ffield warnings
! instead writing them to std output
! df398 01.12.19 added support for low-gradient dispersion model (Liu et al. JPCA 2011)
    dimension dispc6i(nsort),disprei(nsort)
!*********************************************************************
!                                                                    *
!     Read in force field                                            *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In ffinpt'
        call timer(65)
        close (65)
    end if
    rewind (4)
    iline=0
    read (4,'(a40)',end=990,err=990)qffield
    if (nsurp == 0) then
        open (1,file='fort.1',status='unknown',position='append') ! df398 write reaxff warnings to fort.1
        write (1,100)qffield
        close(1)
    end if
    iline=iline+1
!*********************************************************************
!                                                                    *
!     Read in general force field parameters                         *
!                                                                    *
!*********************************************************************
    read (4,1100,end=990,err=990)npar
    iline=iline+1
    do i1=1,npar
        read (4,1300,end=990,err=990)vpar(i1)
        iline=iline+1
    end do
    cutoff=0.01d0*vpar(30)
    swa=vpar(12)
    if (abs(swa) > 0.01d0) then
        open (1,file='fort.1',status='unknown',position='append') ! df398 write reaxff warnings to fort.1
        write (1,*) 'Warning: non-zero value for lower Taper-radius cutoff'
        close (1)
    end if
    swb=vpar(13)
    if (swb < zero) stop &
    'Negative value for upper Taper-radius cutoff'
    if (swb < 5.0d0) then 
        open (1,file='fort.1',status='unknown',position='append')
        write (1,*) 'Warning: very low value for upper Taper-radius cutoff:',swb
        close(1)
    end if
    dispscale=vpar(36)
    if (dispscale.lt.0.5) then
        write (*,*) 'Warning: scale factor in dispersion too small:',dispscale
        if (dispscale == 0.0 .or. dispscale == 20.0) then
           write (*,*) 'dispscale is zero or ffield is not lg. Assuming no lg correction is needed'
        endif
    endif
!*********************************************************************
!                                                                    *
!     Read in atom type data                                         *
!                                                                    *
!*********************************************************************
    read (4,1100,end=990,err=990) nso
    iline=iline+1
    read (4,*,end=990,err=990)
    iline=iline+1
    read (4,*,end=990,err=990)
    iline=iline+1
    read (4,*,end=990,err=990)
    iline=iline+1
    if (nso > nsort) stop 'Maximum number of atom types exceeded'
    ivdwty=0    ! vdWaals type: 1:Shielded Morse, no inner-wall 2: inner wall, no shielding  3: inner wall+shielding 4: lg+inner wall, no shielding
    do i1=1,nso
        read (4,1200,end=990,err=990)qas(i1),rat(i1),aval(i1),amas(i1), &
        rvdw(i1),eps(i1),gam(i1),rapt(i1),stlp(i1)
        iline=iline+1
        read (4,1250,end=990,err=990)alf(i1),vop(i1),valf(i1), &
        valp1(i1),valp2(i1),chi(i1),eta(i1),vnphb
        iline=iline+1
        read (4,1250,end=990,err=990)vnq(i1),vlp1(i1),vincr(i1), &
        bo131(i1),bo132(i1),bo133(i1),sigqeq(i1),default
        iline=iline+1
        read (4,1250,end=990,err=990)vovun(i1),vval1(i1),vlp2(i1), &
        vval3(i1),vval4(i1),rcore2(i1),ecore2(i1),acore2(i1)
        if (dispscale /= 0.0 .and. dispscale /= 20.0) then
          iline=iline+1
          read (4,1260,end=990,err=990)dispc6i(i1),disprei(i1)
        endif
        iline=iline+1
        idef(i1)=int(default)
        nphb(i1)=int(vnphb)

        if (rcore2(i1) > 0.01 .AND. acore2(i1) > 0.01) then          !Inner wall-parameters present
            if (vop(i1) > 0.5) then                                     !Shielding vdWaals-parameters present
                if (ivdwty /= 0 .AND. ivdwty /= 3) then
                    open (1,file='fort.1',status='unknown',position='append')
                    write (1,*)'Warning: inconsistent vdWaals-parameters'
                    write (1,*)'Force field parameters for element:',qas(i1)
                    write (1,*)'indicate inner wall+shielding, but earlier'
                    write (1,*)'atoms indicate different vdWaals-method.'
                    write (1,*)'This may cause division-by-zero errors.'
                    write (1,*)'Keeping vdWaals-setting for earlier atoms.'
                    write (1,*)
                    close(1)
                else
                    if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 add lg dispersion
                        ivdwty=4
                        open (1,file='fort.1',status='unknown',position='append')
                        write (1,*) 'vdWaals type for element ',qas(i1), &
                        ': Shielding+inner-wall+lg'
                        write (1,*) 'Warning: shielding is ignored with lg correction'
                        close(1)
                    else
                        ivdwty=3
                        open (1,file='fort.1',status='unknown',position='append')
                        write (1,*)'vdWaals type for element ',qas(i1), &
                        ': Shielding+inner-wall'
                        close(1)
                    endif
                end if
            else                                                         !No shielding vdWaals parameters present
                if (ivdwty /= 0 .AND. ivdwty /= 2) then
                    open (1,file='fort.1',status='unknown',position='append')
                    write (1,*)'Warning: inconsistent vdWaals-parameters'
                    write (1,*)'Force field parameters for element:',qas(i1)
                    write (1,*)'indicate inner wall without shielding, but earlier'
                    write (1,*)'atoms indicate different vdWaals-method.'
                    write (1,*)'This may cause division-by-zero errors.'
                    write (1,*)'Keeping vdWaals-setting for earlier atoms.'
                    write (1,*)
                    close(1)
                else
                    ivdwty=2
                    if (dispscale /= 0.00 .and. dispscale /= 20.0) then !df398 add lg dispersion
                        ivdwty=4
                        open (1,file='fort.1',status='unknown',position='append')
                        write (1,*) 'vdWaals type for element ',qas(i1), &
                        ': No shielding+inner-wall+lg'
                        close(1)
                    else
                        open (1,file='fort.1',status='unknown',position='append')
                        write (1,*)'vdWaals type for element ',qas(i1), &
                        ': No shielding+inner-wall'
                        close(1)
                    end if
                end if
            end if
        else                                                         !No Inner wall parameters present
            if (dispscale /= 0.0 .and. dispscale /= 20.0) then
                open (1,file='fort.1',status='unknown',position='append')
                write (1,*) 'Error: lg correction requires inner-wall parameters to be present for element ', qas(i1)
                stop 'no inner-wall parameters found but lg correction was included for element'
                close(1)
            endif
            if (vop(i1) > 0.5) then                                     !Shielding vdWaals-parameters present
                if (ivdwty /= 0 .AND. ivdwty /= 1) then
                    open (1,file='fort.1',status='unknown',position='append')
                    write (1,*)'Warning: inconsistent vdWaals-parameters'
                    write (1,*)'Force field parameters for element:',qas(i1)
                    write (1,*)'indicate shielding without inner wall, but earlier'
                    write (1,*)'atoms indicate different vdWaals-method.'
                    write (1,*)'This may cause division-by-zero errors.'
                    write (1,*)'Keeping vdWaals-setting for earlier atoms.'
                    write (1,*)
                    close(1)
                else
                    ivdwty=1
                    open (1,file='fort.1',status='unknown',position='append')
                    write (*,*)'vdWaals type for element ',qas(i1), &
                    ': Shielding, no inner-wall'
                close(1)
                end if
            else
                open (1,file='fort.1',status='unknown',position='append')
                write (*,*)'Error: inconsistent vdWaals-parameters'
                write (*,*)'No shielding or inner-wall set for'
                write (*,*)'element:',qas(i1)
                close(1)
                stop 'Wrong vdWaals parameters, see run.log-output'
            end if
        end if

    end do
!*********************************************************************
!                                                                    *
!     Equate vval3 to valf for first-row elements (25/10/2004)       *
!                                                                    *
!*********************************************************************
    do i1=1,nso
        if (amas(i1) < 21.0) then
            if (vval3(i1) /= valf(i1)) then
                open (1,file='fort.1',status='unknown',position='append')
                write (1,*)'Warning:changed vval3 to valf for ',qas(i1)
                close(1)
                vval3(i1)=valf(i1)
            end if
        end if
    end do
!*********************************************************************
!                                                                    *
!     Calculate van der Waals and Coulomb pair-parameters            *
!                                                                    *
!*********************************************************************
    do i1=1,nso
        do i2=1,nso
            if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 include lg
                dispc6(i1,i2)=sqrt(dispc6i(i1)*dispc6i(i2))
                dispre(i1,i2)=sqrt(4.0*disprei(i1)*disprei(i2))
            endif
            p1co(i1,i2)=sqrt(4.0d0*rvdw(i1)*rvdw(i2))
            p2co(i1,i2)=sqrt(eps(i1)*eps(i2))
            p3co(i1,i2)=sqrt(alf(i1)*alf(i2))
            gamwh=sqrt(vop(i1)*vop(i2))
            gamwco(i1,i2)=1.0d0/gamwh**vpar(29)
            gamch=sqrt(gam(i1)*gam(i2))
            gamcco(i1,i2)=1.0d0/gamch**3
            rob1(i1,i2)=0.50d0*(rat(i1)+rat(i2))
            rob2(i1,i2)=0.50d0*(rapt(i1)+rapt(i2))
            rob3(i1,i2)=0.50d0*(vnq(i1)+vnq(i2))
            rcore(i1,i2)=sqrt(rcore2(i1)*rcore2(i2))
            ecore(i1,i2)=sqrt(ecore2(i1)*ecore2(i2))
            acore(i1,i2)=sqrt(acore2(i1)*acore2(i2))

        end do
    end do
!*********************************************************************
!                                                                    *
!     Read in bond type data                                         *
!                                                                    *
!*********************************************************************
    read (4,1100,end=990,err=990)nboty
    iline=iline+1
    read (4,*,end=990,err=990)
    iline=iline+1
    if (2*nboty > nbotym) stop 'Maximum nr. of bond types exceeded'
    ih=0
    do i1=1,nboty
        ih=ih+1
        read (4,1400,end=990,err=990)nbs(ih,1),nbs(ih,2),de1(ih), &
        de2(ih),de3(ih),psi(ih),pdo(ih),v13cor(ih),popi(ih),vover(ih)
        iline=iline+1
        read (4,1450,end=990,err=990)psp(ih),pdp(ih),ptp(ih), &
        bom(ih),bop1(ih),bop2(ih),ovc(ih),vuncor(ih)
        iline=iline+1
        if (nbs(ih,1) /= nbs(ih,2)) then
            ih=ih+1
            nbs(ih,1)=nbs(ih-1,2)
            nbs(ih,2)=nbs(ih-1,1)
            de1(ih)=de1(ih-1)
            de2(ih)=de2(ih-1)
            de3(ih)=de3(ih-1)
            psi(ih)=psi(ih-1)
            pdo(ih)=pdo(ih-1)
            v13cor(ih)=v13cor(ih-1)
            vover(ih)=vover(ih-1)
            psp(ih)=psp(ih-1)
            pdp(ih)=pdp(ih-1)
            ptp(ih)=ptp(ih-1)
            bop1(ih)=bop1(ih-1)
            bop2(ih)=bop2(ih-1)
        !     bop3(ih)=bop3(ih-1)
        !     bop4(ih)=bop4(ih-1)
            bom(ih)=bom(ih-1)
            popi(ih)=popi(ih-1)
            ovc(ih)=ovc(ih-1)
            vuncor(ih)=vuncor(ih-1)
        end if
    end do
    nboty2=ih

!     do i1=1,nboty2   !copy vpar(11) to bom if bom=zero (triple bond stabilization)
!     if (bom(i1).lt.0.001.and.bom(i1).gt.-0.001) then
!     bom(i1)=vpar(11)
!     end if
!     end do
!*********************************************************************
!                                                                    *
!     Read in off-diagonal parameters                                *
!                                                                    *
!*********************************************************************
    read (4,1100,end=990,err=990)nodmty
    iline=iline+1
    if (nodmty > nodmtym) &
    stop 'Maximum nr. of off-diagonal Morse types exceeded'
    ih=0
    do i1=1,nodmty
        ih=ih+1
        if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 lg included
            read (4,1400,end=990,err=990)nodm1,nodm2,deodmh,rodmh,godmh, &
            rsig,rpi,rpi2,dispc6ij
        else ! no lg
            read (4,1400,end=990,err=990)nodm1,nodm2,deodmh,rodmh,godmh, &
            rsig,rpi,rpi2
        endif
        iline=iline+1
        if (dispscale /= 0.0 .and. dispscale /= 20.0) then !df398 lg included
            if (dispc6ij > zero) dispc6(nodm1,nodm2)=dispc6ij
            if (dispc6ij > zero) dispc6(nodm2,nodm1)=dispc6ij
        endif
        if (rsig > zero) rob1(nodm1,nodm2)=rsig
        if (rsig > zero) rob1(nodm2,nodm1)=rsig
        if (rpi > zero) rob2(nodm1,nodm2)=rpi
        if (rpi > zero) rob2(nodm2,nodm1)=rpi
        if (rpi2 > zero) rob3(nodm1,nodm2)=rpi2
        if (rpi2 > zero) rob3(nodm2,nodm1)=rpi2
        if (rodmh > zero) p1co(nodm1,nodm2)=2.0d0*rodmh
        if (rodmh > zero) p1co(nodm2,nodm1)=2.0d0*rodmh
        if (deodmh > zero) p2co(nodm1,nodm2)=deodmh
        if (deodmh > zero) p2co(nodm2,nodm1)=deodmh
        if (godmh > zero) p3co(nodm1,nodm2)=godmh
        if (godmh > zero) p3co(nodm2,nodm1)=godmh
    end do
!*********************************************************************
!                                                                    *
!     Read in valency angle and conjugation type data                *
!                                                                    *
!*********************************************************************
    read (4,1100,end=990,err=990)nvaty
    iline=iline+1
    if (nvaty > nvatym) &
    stop 'Maximum nr. of valency angle types exceeded'
    do i1=1,nvaty
        read (4,1500,end=990,err=990)nvs(i1,1),nvs(i1,2), &
        nvs(i1,3),th0(i1),vka(i1),vka3(i1),vka8(i1),vkac(i1),vkap(i1), &
        vval2(i1)
        iline=iline+1
    end do
!*********************************************************************
!                                                                    *
!     Read in torsion angle type data                                *
!                                                                    *
!*********************************************************************
    read (4,1100,end=990,err=990)ntoty
    iline=iline+1
    if (ntoty > ntotym) &
    stop 'Maximum nr. of torsion angle types exceeded'
    do i1=1,ntoty
        read (4,1600,end=990,err=990)nts(i1,1),nts(i1,2),nts(i1,3), &
        nts(i1,4),v1(i1), &
        v2(i1),v3(i1),v4(i1),vconj(i1),v2bo(i1),v3bo(i1)
        iline=iline+1
    end do
!*********************************************************************
!                                                                    *
!     Read in hydrogen bond type data                                *
!                                                                    *
!*********************************************************************
    read (4,1100,end=990,err=990)nhbty
    iline=iline+1
    if (nhbty > nhbtym) &
    stop 'Maximum nr. of hydrogen bond types exceeded'
    do i1=1,nhbty
        read (4,1500,end=990,err=990)nhbs(i1,1),nhbs(i1,2), &
        nhbs(i1,3),rhb(i1),dehb(i1),vhb1(i1),vhb2(i1)
        iline=iline+1
    end do
!*********************************************************************
!                                                                    *
!     Calculate vdWaals interaction parameters                       *
!                                                                    *
!*********************************************************************
    do i1=1,nso
        do i2=1,nso
            rr=(rvdw(i1)+rvdw(i2))
            rr2=rr*rr
            eps2=sqrt(eps(i1)*eps(i2))
            rr6=rr2*rr2*rr2
            pvdw1(i1,i2)=eps2*rr6*rr6
            pvdw1(i2,i1)=eps2*rr6*rr6
            pvdw2(i1,i2)=2.0d0*eps2*rr6
            pvdw2(i2,i1)=2.0d0*eps2*rr6
        end do
    end do
!*********************************************************************
!                                                                    *
!     Error part                                                     *
!                                                                    *
!*********************************************************************
    goto 999
    990 write (*,*)'Error or end-of-file reading unit 4 on line:',iline
    stop
    999 continue
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format ('Force field identifier: ',a40)
    1100 format (i3,2x,a2,3x,3d22.15)
    1200 format (1x,a2,10d9.4)
    1250 format (3x,10f9.4)
    1260 format (3x,2f9.4)
    1300 format (d10.4)
    1400 format (2i3,8f9.4)
    1450 format (6x,8f9.4)
    1500 format (3i3,7f9.4)
    1600 format (4i3,7f9.4)
    return
    end subroutine ffinpt
!*********************************************************************
!*********************************************************************

    subroutine intcon

!*********************************************************************
    include 'cbka.blk'
    character(120) :: qhulp
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In intcon'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Input fixed connection table                                   *
!                                                                    *
!*********************************************************************
    do i1=1,na
        ia(i1,2)=0
        iag(i1,2)=0
        do i2=1,mbond
            nubon2(i1,i2)=0
        end do
    end do
    nsbmax=0
    nsbma2=0
    open (12,file='cnt.in',status='old',err=900)
    read (12,'(i4)',end=910,err=910)nah
    if (nah /= na) stop 'Inconsistent number of atoms in cnt-file'
    do i1=1,na
        read (12,'(a120)',end=910,err=910)qhulp
        read (qhulp,'(4x,i4)')ish
        if (ish /= ia(i1,1)) then
            write (*,*)'Wrong atom type in connection table: atom nr.',i1
            stop 'Wrong atom type in connection table'
        end if
        ihulp=9
        10 read (qhulp(ihulp:ihulp+4),'(i4)')icoh
        if (icoh > 0) then
            ia(i1,2)=ia(i1,2)+1
            ia(i1,ia(i1,2)+2)=icoh
            iag(i1,2)=iag(i1,2)+1
            iag(i1,ia(i1,2)+2)=icoh
            if (ia(i1,2) > nsbma2) nsbma2=ia(i1,2)
            if (iag(i1,2) > nsbmax) nsbmax=ia(i1,2)
            ihulp=ihulp+4
        end if
        if (icoh > 0) goto 10
    end do

    close (12)
    call sortmol(ia,na1m,nmoln)
    call sortmol(iag,na1m,nmoln)
!*********************************************************************
!                                                                    *
!     Find number of bonds and force field types                     *
!                                                                    *
!*********************************************************************
    nbon=0
    do i1=1,na
        ity1=ia(i1,1)
        do i2=1,ia(i1,2)
            i3=ia(i1,2+i2)
            if (i3 > i1) then
                nbon=nbon+1
                ib(nbon,2)=i1
                ib(nbon,3)=i3
                ity2=ia(i3,1)
                if (ity1 > ity2) then
                    ity2=ia(i1,1)
                    ity1=ia(i3,1)
                end if
                do i4=1,nboty2
                    if (ity1 == nbs(i4,1) .AND. ity2 == nbs(i4,2)) ibtyp=i4
                end do
                if (ibtyp == 0) then
                    write (*,*)qa(i1),'-',qa(i3),'Fatal: Unknown bond in molecule'
                    stop
                end if
                ib(nbon,1)=ibtyp
                nubon2(i1,ia(i1,2))=nbon
                nubon2(i3,ia(i3,2))=nbon
            end if

        end do
    end do


!     do i1=1,na
!     write (7,'(10i4)')(ia(i1,i2),i2=1,10)
!     end do
!     write (7,*)nbon
          
    return
    900 stop 'Error opening connection table (cnt.in)'
    910 stop 'Error or end of file in connection table (cnt.in)'
    end subroutine intcon
!*********************************************************************
!*********************************************************************

    subroutine readmol

!*********************************************************************
    include 'cbka.blk'
    character(120) :: qhulp
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readmol'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Read in fixed molecular definitions                            *
!                                                                    *
!*********************************************************************
    open (12,file='moldef.in',status='old',err=900)
    nmolo=0
    nmolo5=0
    do i1=1,na
        read (12,*,end=910,err=920)irom,imol
        ia (i1,3+mbond)=imol
        iag (i1,3+mbond)=imol
        if (imol > nmolo) nmolo=imol
        if (imol > nmolo5) nmolo5=imol
    end do
    return
    900 stop 'Error opening molecular definitions (moldef.in)'
    910 stop 'End of file reading molecular definitions (moldef.in)'
    920 stop 'Error reading molecular definitions (moldef.in)'
    end subroutine readmol
!*********************************************************************
!*********************************************************************

    subroutine intcor

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Determine internal coordinates in molecule                     *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In intcor'
        call timer(65)
        close (65)
    end if

    call distan

    if (ireflx == 1 .OR. irefly == 1 .OR. ireflz == 1) call reflect  !Reflective cell boundaries
    if (mod(mdstep,nvlist) == 0 .AND. mdstep > 1) then
        call vlist
    end if

    if (iconne == 0 .OR. nbon == 0) then
        call srtbon1
    else
        if (iconne > 0 .AND. nbon > 0) call srtbon2
    end if
!     call srtang
!     call srttor

!*********************************************************************
!                                                                    *
!     Determine changes in molecules                                 *
!                                                                    *
!*********************************************************************
    call molec
!*********************************************************************
!                                                                    *
!     Distribute electrons                                           *
!                                                                    *
!*********************************************************************
!     if (mod(mdstep,nrddf).eq.0.and.nrdd.eq.1) then
!     call distrad
!     end if

    return
    end subroutine intcor
!*********************************************************************
!*********************************************************************

    subroutine inilp

!*********************************************************************
    include 'cbka.blk'
    character(60) :: qm2
!*********************************************************************
!                                                                    *
!     Determine initial number of electrons in each molecule         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In inilp'
        call timer(65)
        close (65)
    end if
    do i1=1,nmolo
        elmol(i1)=0.0d0
        do i2=1,na
            if (ia(i2,3+mbond) == i1) then
                it1=ia(i2,1)
                elmol(i1)=elmol(i1)+stlp(it1)
            end if
        end do
    end do

    do i1=1,nmolo5
        elmol2(i1)=0.0d0
        do i2=1,na
            if (iag(i2,3+mbond) == i1) then
                it1=iag(i2,1)
                elmol2(i1)=elmol2(i1)+stlp(it1)
            end if
        end do
    end do
!*********************************************************************
!                                                                    *
!     Adjust normal amount of lone pairs (unit 2)                    *
!                                                                    *
!*********************************************************************
    iadj=0
    read (2,100,err=15,end=15)idev,qm2
!     write (95,*)qm2
!     write (95,*)qmol
    if (qm2 /= qmol) then
        write (*,*)qm2
        write (*,*)qmol
        stop 'Wrong molecule in inilp-file'
    end if
    do i1=1,idev
        read (2,200,err=98,end=98)imolnr,devlp
        if (devlp /= zero) then
            write (*,*)'Warning: deviation from normal electron amount'
            elmol(imolnr)=elmol(imolnr)+devlp
        end if
    end do
    15 continue
    sumelec=0.0d0
    do i1=1,nmolo
        sumelec=sumelec+elmol(i1)
    end do
    97 goto 99
    98 stop 'Error in inilp-file'
    99 continue
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format (i3,a60)
    200 format (i3,f8.4)
    return
    end subroutine inilp
!*********************************************************************
!*********************************************************************

    subroutine molec

!*********************************************************************
    include 'cbka.blk'
    dimension nmolo2(nat),iseen(nmolmax),isee2(nmolmax)
!*********************************************************************
!                                                                    *
!     Determine changes in molecules                                 *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In molec'
        call timer(65)
        close (65)
    end if
    npreac=0

    do i1=1,nmolo
        natmol=0
        do i2=1,na
            if (ia(i2,3+mbond) == i1) then
                natmol=natmol+1
                nmolat(i1,natmol+1)=i2
            end if
        end do
        nmolat(i1,1)=natmol
    end do

    if (nmolo5 < nmolo5o) nradcount=0     !reset reaction counter
    do i1=1,nmolo5
        natmol=0
        do i2=1,na
            if (iag(i2,3+mbond) == i1) then
                natmol=natmol+1
                nmolat2(i1,natmol+1)=i2
            end if
        end do
        nmolat2(i1,1)=natmol
    end do
    nmolo5o=nmolo5

    do i1=nmolo+1,nmoloold
        do i2=1,nmolat(i1,1)
            nmolat(i1,1+i2)=0
        end do
        nmolat(i1,1)=0
    end do

    do i1=1,nmolo
        elmol(i1)=0.0d0
        do i2=1,nmolat(i1,1)
            ihu=nmolat(i1,i2+1)
            ity=ia(ihu,1)
            elmol(i1)=elmol(i1)+stlp(ity)
        end do
    end do
     
    do i1=1,nmolo5
        elmol2(i1)=0.0d0
        do i2=1,nmolat2(i1,1)
            ihu=nmolat2(i1,i2+1)
            ity=iag(ihu,1)
            elmol2(i1)=elmol2(i1)+stlp(ity)
        end do
    end do
     
    return

!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format (' Iter.',i6,': mol.',i3,' and mol.',i3,' reacted to', &
    ' mol.',i3,' and mol.',i3)
    200 format (' Iter.',i6,': mol.',i3,' reacted to mol.',i3,' and', &
    ' mol.',i3)
    300 format (' Iter.',i6,': mol.',i3,' and mol.',i3,' reacted to', &
    ' mol.',i3)
    end subroutine molec
!*********************************************************************
!*********************************************************************

    subroutine encalc

!*********************************************************************
    include 'cbka.blk'
    real :: tarray(2)
    logical :: found
!*********************************************************************
!                                                                    *
!     Calculate energy and first derivatives                         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In encalc'
        call timer(65)
        close (65)
    end if
    estrc=0.0d0

    do i1=1,na
        estrain(i1)=0.0d0
        do i2=1,3
            d(i2,i1)=0.0d0
            if (icpres == 1) then
                do i3=1,27
                    dcell(i2,i1,i3)=0.0d0
                end do
            end if
        end do
    end do

    eb=zero
    ea=zero
    elp=zero
    emol=zero
    ev=zero
    ehb=zero
    ecoa=zero
    epen=zero
    et=zero
    eco=zero
    eres=zero
    epist=zero
    eradbo=zero
    efi=zero
    if (icpres == 1) then
    !*********************************************************************
    !                                                                    *
    !     Calculate minimum images for NPT                               *
    !                                                                    *
    !*********************************************************************
        do  i1=1,na-1
            do  i2=i1+1,na
                dx1=c(i1,1)-c(i2,1)
                dy1=c(i1,2)-c(i2,2)
                dz1=c(i1,3)-c(i2,3)
                dismin=1d+10
                do  k1=-1,1
                    do  k2=-1,1
                        do  k3=-1,1
                            a1=dx1+k1*tm11
                            a2=dy1+k1*tm21+k2*tm22
                            a3=dz1+k1*tm31+k2*tm32+k3*tm33
                            rr=sqrt(a1*a1+a2*a2+a3*a3)
                            if (rr < dismin) then
                                dismin=rr
                                ixsav=k1
                                iysav=k2
                                izsav=k3
                            end if
                        enddo
                    enddo
                 enddo
                nmpx(i1,i2)=ixsav
                nmpx(i2,i1)=-ixsav
                nmpy(i1,i2)=iysav
                nmpy(i2,i1)=-iysav
                nmpz(i1,i2)=izsav
                nmpz(i2,i1)=-izsav
            enddo
        enddo
    end if
     
    do i1=1,na
        nmpx(i1,i1)=0
        nmpy(i1,i1)=0
        nmpz(i1,i1)=0
    end do

    call boncor

    if (icobo == 1) then    !base molecule numbers on corrected bond orders

        do i1=1,na
            do i2=1,mbond+3
                iag(i1,i2)=0
            end do
        end do
    !***********************************************************************
    !                                                                      *
    !     Create connection table based on corrected bond orders           *
    !                                                                      *
    !***********************************************************************
        do i1=1,nbon
            if (bo(i1) > cutof3) then
                j1=ib(i1,2)
                j2=ib(i1,3)
                ity1=ia(j1,1)
                ity2=ia(j2,1)
                if (ity1 /= igno .AND. ity2 /= igno) then
                    iag(j1,2)=iag(j1,2)+1
                    iag(j1,2+iag(j1,2))=j2
                    iag(j2,2)=iag(j2,2)+1
                    iag(j2,2+iag(j2,2))=j1
                end if
            end if
        end do
    !*********************************************************************
    !                                                                    *
    !     Find molecules                                                 *
    !                                                                    *
    !*********************************************************************
        nmolo5=0
        found= .FALSE. 
        DO 61 k1=1,na
            IF (iag(K1,3+mbond) == 0) found= .TRUE. 
            IF (iag(K1,3+mbond) > nmolo5) nmolo5=iag(K1,3+mbond)
        61 END DO
        IF ( .NOT. FOUND) GOTO 62
    !***********************************************************************
    !                                                                      *
    !     Molecule numbers are assigned. No restrictions are made for the  *
    !     sequence of the numbers in the connection table.                 *
    !                                                                      *
    !***********************************************************************
        N3=1
        64 N2=N3
        nmolo5=nmolo5+1
        if (nmolo5 > nmolmax) stop 'Too many molecules in system'
        iag(N2,3+mbond)=nmolo5
        67 FOUND= .FALSE. 
        DO 66 N1=N2+1,na
            IF (iag(N1,3+mbond) /= 0) GOTO 66
            DO 65 L=1,mbond
                IF (iag(N1,l+2) == 0) GOTO 66
                IF (iag(iag(N1,l+2),3+mbond) == nmolo5) THEN
                    FOUND= .TRUE. 
                    iag(N1,3+mbond)=nmolo5
                    GOTO 66
                ENDIF
            65 END DO
        66 END DO
        IF (FOUND) GOTO 67
        DO 63 N3=N2+1,NA
            if (iag(N3,3+mbond) == 0) goto 64
        63 END DO
    !***********************************************************************
    !                                                                      *
    !     The assigned or input molecule numbers are checked for their     *
    !     consistency.                                                     *
    !                                                                      *
    !***********************************************************************
        62 FOUND= .FALSE. 
        DO 72 N1=1,NA
            DO 71 L=1,mbond
                IF (iag(N1,L+2) == 0) GOTO 72
                IF (iag(iag(N1,L+2),3+mbond) /= iag(N1,3+mbond)) THEN
                    FOUND= .TRUE. 
                ENDIF
            71 END DO
        72 END DO
        IF (FOUND) THEN
            write (7,'(i4,a40)')na,qmol
            do i1=1,na
                write (7,'(40i4)')i1,iag(i1,1),(iag(i1,2+i2),i2=1,nsbmax), &
                iag(i1,3+mbond)
            end do
            STOP 'Mol.nrs. not consistent; maybe wrong cell parameters'
        ENDIF

    !     write (64,*)mdstep,nmolo5
    end if

    call lonpar
    call covbon
    call ovcor

    if (nreac == 0) call molen

    call srtang   !Determine valency angles
    call srttor   !Determine torsion angles
!     call srtoop   !Determine out of plane angles
    call srthb    !Determine hydrogen bonds

    call calval
    call valang

!     call oopang

    call torang
    call hbond

    if (nmm /= 1 .AND. mod(mdstep,nchaud) == 0) then
        call chargess
    end if

    if (nmm == 1 .AND. mod(nit,nchaud) == 0) then
        call chargess
    end if

    call nonbon
    call efield

    if (nrdd == 1) call radbo

    call restraint
    if (edeep > zero) call piston

    if (icpres == 1) call calcpres

    estrc=eb+ea+elp+ev+ecoa+emol+epen+et+ehb+eco+ew+ep+epist+ &
    ncha2*ech+efi

    if (ndebug == 1) then    !write forces to unit 66
        open (66,file='fort.66',status='unknown',position='append')
        !write (66,*)'MD-step:',mdstep
        write (66,*)'ReaxFF> Total forces'
        write (66,*)'ReaxFF> Force matrix:'
        do i1=1,na
            write (66,'(i6,3f20.8)')i1,(d(i2,i1),i2=1,3)
        end do
        close (66)
    end if

    if (estrc > zero) return
    if (estrc <= zero) then
        goto 10
    else
        write (*,*)mdstep
        write (92,*)eb,ea,elp,ev,ecoa,emol,epen,eoop,et,eco,ew, &
        ep,ech,eres,eradbo
        write (92,*) 'bosia:',bosia, 'bopia:',bopia, 'bopi2a:', bopi2a
        write (92,*) 'estrc:', estrc, '-de1h:',-de1h, 'exphu1:', exphu1, '-de2h:', -de2h, '-deh3:',-de3h
        write (*,*)eb,ea,elp,ev,ecoa,emol,epen,eoop,et,eco,ew, &
        ep,ech,eres,eradbo
        write (*,*) bosia, bopia, bopi2a
        write (*,*) 'eb:', eb
        write (*,*) 'ea:', ea
        write (*,*) 'elp:', elp
        write (*,*) 'ev:', ev
        write (*,*) 'ecoa:', ecoa
        write (*,*) 'emol:', emol
        write (*,*) 'epen:', epen
        write (*,*) 'et:', et
        write (*,*) 'ehb:', ehb
        write (*,*) 'eco:', eco
        write (*,*) 'ew:', ew
        write (*,*) 'ep:', ep
        write (*,*) 'epist:', epist
        write (*,*) 'ncha2*ech:', ncha2*ech
        write (*,*) 'efi:', efi
        stop 'Energy not a number'
    end if

    10 continue
    return
    end subroutine encalc
!*********************************************************************
!*********************************************************************

    subroutine distan

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Set up interatomic distance calculation                        *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In distan'
        call timer(65)
        close (65)
    end if
    iortho=1
    if (qr == 'F' .OR. qr == 'Y' .OR. qr == '3' .OR. qr == '5' .or. &
    qr == 'P' .OR. qr == 'B' .OR. qr == 'S' .OR. qr == 'R') then
        if (angle(1) /= 90.0d0 .OR. angle(2) /= 90.0d0 .OR. angle(3) /= 90.0d0) &
        iortho=0           !Non-orthogonal
    !     if (iortho.eq.0) stop 'Only for orthogonal systems'
        halfa=angle(1)*dgrrdn
        hbeta=angle(2)*dgrrdn
        hgamma=angle(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axis(1)*sinbet*sinphi
        tm21=axis(1)*sinbet*cosphi
        tm31=axis(1)*cosbet
        tm22=axis(2)*sinalf
        tm32=axis(2)*cosalf
        tm33=axis(3)
    end if
    aaxh=sqrt(tm11*tm11+tm21*tm21+tm31*tm31)
    baxh=sqrt(tm22*tm22+tm32*tm32)
    caxh=tm33


    return
    end subroutine distan

!*********************************************************************
!*********************************************************************

    subroutine dista2 (n1,n2,dista,dx,dy,dz)

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Determine interatomic distances                                *
!                                                                    *
!*********************************************************************
!     if (ndebug.eq.1) then
!     open (65,file='fort.65',status='unknown',position='append')
!     write (65,*) 'In dista2'
!     call timer(65)
!     close (65)
!     end if
          
    if (iortho == 1) then
        dx=c(n1,1)-c(n2,1)
        dy=c(n1,2)-c(n2,2)
        dz=c(n1,3)-c(n2,3)
        dx=dx-anint(dx/aaxh)*aaxh
        dy=dy-anint(dy/baxh)*baxh
        dz=dz-anint(dz/caxh)*caxh
        dista=sqrt(dx*dx+dy*dy+dz*dz)

    else

        dismin=1d+10
        do  k1=-1,1
            dx=(c(n1,1)-c(n2,1)+k1*tm11)
            do  k2=-1,1
                dy=(c(n1,2)-c(n2,2)+k1*tm21+k2*tm22)
                do  k3=-1,1
                    dz=(c(n1,3)-c(n2,3)+k1*tm31+k2*tm32+k3*tm33)
                    dis=sqrt(dx*dx+dy*dy+dz*dz)
                    if ((dis < dismin) .AND. dis > 0.001d0) then
                        dxmin=dx
                        dymin=dy
                        dzmin=dz
                        dismin=dis
                    end if
                enddo
             enddo
         enddo
        dista=dismin
        dx=dxmin
        dy=dymin
        dz=dzmin
              
    end if

    return
    end subroutine dista2
!*********************************************************************
!*********************************************************************

    subroutine distrad

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Determine distance between radical centres and nearest double  *
!     bonds                                                          *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In distrad'
        call timer(65)
        close (65)
    end if
    do i1=1,20
        irap(i1)=0
    end do
    irac=0
    do i1=1,nmolo
        if (elmol(i1) > 2*int(elmol(i1)*0.50d0)) then
            irac=irac+1
            vlpmax=0.0d0
            do i2=1,nmolat(i1,1)
                ihu=nmolat(i1,i2+1)
                if (vlp(ihu) > vlpmax) then
                    vlpmax=vlp(ihu)
                    iradpos=ihu
                end if
            end do
            irap(irac)=iradpos
            dirmin=50.0d0  !Excessively large value

            do i2=1,na
                imol2=ia(i2,3+mbond)
                if (imol2 /= i1 .AND. &
                int(elmol(imol2)) == 2*int(elmol(imol2)*0.50d0)) then
                    do i3=1,ia(i2,2)
                        ihu=nubon2(i2,i3)
                        if (bo(ihu) > 1.3) then
                            call dista2(iradpos,i2,dirb1,dx,dy,dz)
                            call dista2(iradpos,ia(i2,2+i3),dirb2,dx,dy,dz)
                            if (dirb1 < dirmin) then
                                dirmin=dirb1
                                irdo(irac,1)=i2
                                irdo(irac,2)=ia(i2,2+i3)
                            end if
                            if (dirb2 < dirmin) then
                                dirmin=dirb2
                                irdo(irac,1)=i2
                                irdo(irac,2)=ia(i2,2+i3)
                            end if
                        end if
                    end do
                end if
            end do
                  
        end if
    end do
          
    return
    end subroutine distrad
!*********************************************************************
!*********************************************************************

    subroutine vlist

!*********************************************************************
    include 'cbka.blk'
    dimension fc(nat,3)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In vlist'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Update Verlet list                                             *
!                                                                    *
!*********************************************************************
    nvpair=0
    nvlself=0
    kxt=kx
    kyt=ky
    kzt=kz
    if (kxt == 0) kxt=1
    if (kyt == 0) kyt=1
    if (kzt == 0) kzt=1

    call trarot2  !Translate atoms back into unit cell

    do i1=1,na-1
        do i2=i1+1,na
        !*********************************************************************
        !                                                                    *
        !     Construct periodic images for each interaction                 *
        !                                                                    *
        !*********************************************************************
            dx1=c(i1,1)-c(i2,1)
            dy1=c(i1,2)-c(i2,2)
            dz1=c(i1,3)-c(i2,3)
            do k1=-kxt,kxt
                do k2=-kyt,kyt
                    do  k3=-kzt,kzt
                        a1=dx1+k1*tm11
                        a2=dy1+k1*tm21+k2*tm22
                        a3=dz1+k1*tm31+k2*tm32+k3*tm33
                        rr=sqrt(a1*a1+a2*a2+a3*a3)
                        if (rr < swb+vrange) then
                            nvpair=nvpair+1
                            nvl1(nvpair)=i1
                            nvl2(nvpair)=i2
                            nvlx(nvpair)=k1
                            nvly(nvpair)=k2
                            nvlz(nvpair)=k3
                            nvlbo(nvpair)=0
                            if (rr < vlbora) nvlbo(nvpair)=1
                        end if
                    enddo
                enddo
            enddo
         enddo
      enddo
!*********************************************************************
!                                                                    *
!     Add interaction of atoms with the corresponding atom           *
!     in the surrounding periodic cells                              *
!                                                                    *
!*********************************************************************
    do i1=1,na
        do k1=-kx,kx
            do k2=-ky,ky
                do k3=-kz,kz
                    if (k1 /= 0 .OR. k2 /= 0 .OR. k3 /= 0) then
                        a1=k1*tm11
                        a2=k1*tm21+k2*tm22
                        a3=k1*tm31+k2*tm32+k3*tm33
                        rr=sqrt(a1*a1+a2*a2+a3*a3)
                        if (rr < swb+vrange) then
                            nvpair=nvpair+1
                            nvlself=nvlself+1
                            nvl1(nvpair)=i1
                            nvl2(nvpair)=i1
                            nvlx(nvpair)=k1
                            nvly(nvpair)=k2
                            nvlz(nvpair)=k3
                            nvlbo(nvpair)=0
                            if (rr < vlbora) nvlbo(nvpair)=1
                        end if
                    end if
                enddo
            enddo
         enddo
    END DO

    return
    end subroutine vlist
!*********************************************************************
!*********************************************************************

    subroutine srtbon1

!*********************************************************************
    include 'cbka.blk'
    logical :: found
!*********************************************************************
!                                                                    *
!     Determine connections within the molecule                      *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In srtbon1'
        call timer(65)
        close (65)
    end if
    do i1=1,na
        abo(i1)=0.0d0
    end do
    nbon=0
    nbon2=0
    nsbmax=0
    nsbma2=0

    if (imolde == 0) then

        nmolo=0
        nmolo5=0
    end if
    if (imolde == 0) then
        do i1=1,na
            do i2=2,mbond+3
                ia(i1,i2)=0
                iag(i1,i2)=0
            end do
        end do
             
    else

        do i1=1,na
            do i2=2,mbond+2
                ia(i1,i2)=0
                iag(i1,i2)=0
            end do
        end do

    end if

    do i1=1,na
        do i2=1,mbond
            nubon1(i1,i2)=0
            nubon2(i1,i2)=0
        end do
    end do
     
!     do 10 ivl=1,nvpair
    do 10 ivl=1,nvpair
        if (nvlbo(ivl) == 0) goto 10  !not in bond order range
        i1=nvl1(ivl)
        i2=nvl2(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        dxm=c(i1,1)-c(i2,1)+ix*tm11
        dym=c(i1,2)-c(i2,2)+ix*tm21+iy*tm22
        dzm=c(i1,3)-c(i2,3)+ix*tm31+iy*tm32+iz*tm33
        dis=sqrt(dxm*dxm+dym*dym+dzm*dzm)
    !     do 10 i1=1,na-1
    !     do 10 i2=i1+1,na
    !     call dista2(i1,i2,dis,dxm,dym,dzm)
        nubon(i1,i2)=0
        nubon(i2,i1)=0
        ih1=ia(i1,1)
        ih2=ia(i2,1)
    !     if (dis.gt.5.0*rob) goto 10
        disdx=dxm/dis
        disdy=dym/dis
        disdz=dzm/dis
        itype=0
        if (ih1 > ih2) then
            ih1=ia(i2,1)
            ih2=ia(i1,1)
        end if
        do i3=1,nboty2
            if (ih1 == nbs(i3,1) .AND. ih2 == nbs(i3,2)) itype=i3
        end do
        if (itype == 0 .AND. rat(ih1) > zero .AND. rat(ih2) > zero) then
            call mdsav(1,qfile(nprob))
            write (*,*)qa(i1),'-',qa(i2),'Fatal: Unknown bond in molecule'
            stop
        end if

        rhulp=dis/rob1(ih1,ih2)
         
    !*********************************************************************
    !                                                                    *
    !     Determine bond orders                                          *
    !                                                                    *
    !*********************************************************************
        rh2=zero
        rh2p=zero
        rh2pp=zero
        ehulp=zero
        ehulpp=zero
        ehulppp=zero
        if (rapt(ih1) > zero .AND. rapt(ih2) > zero) then
            rhulp2=dis/rob2(ih1,ih2)
            rh2p=rhulp2**ptp(itype)
            ehulpp=exp(pdp(itype)*rh2p)
        end if
        if (vnq(ih1) > zero .AND. vnq(ih2) > zero) then
            rhulp3=dis/rob3(ih1,ih2)
            rh2pp=rhulp3**popi(itype)
            ehulppp=exp(pdo(itype)*rh2pp)
        end if

        if (rat(ih1) > zero .AND. rat(ih2) > zero) then
            rh2=rhulp**bop2(itype)
            ehulp=(1.0d0+cutoff)*exp(bop1(itype)*rh2)
        end if

        bor=ehulp+ehulpp+ehulppp
        borsi=ehulp
        borpi=ehulpp
        borpi2=ehulppp
        dbordrob=bop2(itype)*bop1(itype)*rh2*(1.0d0/dis)*ehulp+ &
        ptp(itype)*pdp(itype)*rh2p*(1.0d0/dis)*ehulpp+ &
        popi(itype)*pdo(itype)*rh2pp*(1.0d0/dis)*ehulppp
        
        dborsidrob=bop2(itype)*bop1(itype)*rh2*(1.0d0/dis)*ehulp
        dborpidrob=ptp(itype)*pdp(itype)*rh2p*(1.0d0/dis)*ehulpp
        dborpi2drob=popi(itype)*pdo(itype)*rh2pp*(1.0d0/dis)*ehulppp
             
        nbon2=nbon2+1
        j1=i1
        j2=i2

    !*********************************************************************
    !                                                                    *
    !     Determine bond orders                                          *
    !                                                                    *
    !*********************************************************************
        if (bor > cutoff) then
            nbon=nbon+1
            ib(nbon,1)=itype
            ib(nbon,2)=j1
            ib(nbon,3)=j2
            ibsym(nbon)=ivl
            drdc(nbon,1,1)=disdx
            drdc(nbon,2,1)=disdy
            drdc(nbon,3,1)=disdz
            drdc(nbon,1,2)=-disdx
            drdc(nbon,2,2)=-disdy
            drdc(nbon,3,2)=-disdz


            if (bor < 4.0D0*cutoff) then
                !call ataper(.false.,borsi,cutoff,6.0d0*cutoff,ft,dftdBO,d2ftdBO2,d3ftdBO3,.true.,.false.,.false.)
                call valtaper(.false.,borsi,cutoff,4.0D0*cutoff,ft,dftdBO)
                bosi(nbon)=ft*(ehulp-cutoff)
                dborsidrob=(ft + dftdBO*(ehulp - cutoff))*bop2(itype)*bop1(itype)*rh2*(1.0d0/dis)*ehulp
            else
                bosi(nbon) = ehulp - cutoff
            end if

            if (bor < 4.0D0*cutoff) then
                !call ataper(.false.,borpi,cutoff,6.0d0*cutoff,ft,dftdBO,d2ftdBO2,d3ftdBO3,.true.,.false.,.false.)
                call valtaper(.false.,borpi,cutoff,4.0D0*cutoff,ft,dftdBO)
                bopi(nbon)=ft*ehulpp
                dborpidrob=(ft + dftdBO*ehulpp)*ptp(itype)*pdp(itype)*rh2p*(1.0d0/dis)*ehulpp
            else
                bopi(nbon) = ehulpp
            end if

            if (bor < 4.0D0*cutoff) then
                !call ataper(.false.,borpi2,cutoff,6.0d0*cutoff,ft,dftdBO,d2ftdBO2,d3ftdBO3,.true.,.false.,.false.)
                call valtaper(.false.,borpi2,cutoff,4.0D0*cutoff,ft,dftdBO)
                bopi2(nbon)=ft*ehulppp
                dborpi2drob=(ft + dftdBO*ehulppp)*popi(itype)*pdo(itype)*rh2pp*(1.0d0/dis)*ehulppp
            else
                bopi2(nbon) = ehulppp
            end if

            if (bor < 4.0D0*cutoff) then
                !call ataper(.false.,bor,cutoff,6.0d0*cutoff,ft,dftdBO,d2ftdBO2,d3ftdBO3,.true.,.false.,.false.)
                call valtaper(.false.,bor,cutoff,4.0D0*cutoff,ft,dftdBO)
                abo(i1)=abo(i1)+ft*(bor-cutoff)
                if (i1 /= i2) abo(i2)=abo(i2)+ft*(bor-cutoff)
                bo(nbon)=ft*(bor-cutoff)
                bos(nbon)=ft*(bor-cutoff)
                dbordrob = dborsidrob + dborpidrob + dborpi2drob
            else
                abo(i1)=abo(i1)+bor-cutoff
                if (i1 /= i2) abo(i2)=abo(i2)+bor-cutoff
                bo(nbon)=bor-cutoff
                bos(nbon)=bor-cutoff
            end if

            rbo(nbon)=dis
            dbodr(nbon)=dbordrob
        !     dbosidr(nbon)=dborsidrob
            dbopidr(nbon)=dborpidrob
            dbopi2dr(nbon)=dborpi2drob
            dbodc(nbon,1,1)=dbodr(nbon)*drdc(nbon,1,1)
            dbodc(nbon,2,1)=dbodr(nbon)*drdc(nbon,2,1)
            dbodc(nbon,3,1)=dbodr(nbon)*drdc(nbon,3,1)
            dbodc(nbon,1,2)=dbodr(nbon)*drdc(nbon,1,2)
            dbodc(nbon,2,2)=dbodr(nbon)*drdc(nbon,2,2)
            dbodc(nbon,3,2)=dbodr(nbon)*drdc(nbon,3,2)
        !     dbosidc(nbon,1,1)=dbosidr(nbon)*drdc(nbon,1,1)
        !     dbosidc(nbon,2,1)=dbosidr(nbon)*drdc(nbon,2,1)
        !     dbosidc(nbon,3,1)=dbosidr(nbon)*drdc(nbon,3,1)
        !     dbosidc(nbon,1,2)=dbosidr(nbon)*drdc(nbon,1,2)
        !     dbosidc(nbon,2,2)=dbosidr(nbon)*drdc(nbon,2,2)
        !     dbosidc(nbon,3,2)=dbosidr(nbon)*drdc(nbon,3,2)
            dbopidc(nbon,1,1)=dbopidr(nbon)*drdc(nbon,1,1)
            dbopidc(nbon,2,1)=dbopidr(nbon)*drdc(nbon,2,1)
            dbopidc(nbon,3,1)=dbopidr(nbon)*drdc(nbon,3,1)
            dbopidc(nbon,1,2)=dbopidr(nbon)*drdc(nbon,1,2)
            dbopidc(nbon,2,2)=dbopidr(nbon)*drdc(nbon,2,2)
            dbopidc(nbon,3,2)=dbopidr(nbon)*drdc(nbon,3,2)
            dbopi2dc(nbon,1,1)=dbopi2dr(nbon)*drdc(nbon,1,1)
            dbopi2dc(nbon,2,1)=dbopi2dr(nbon)*drdc(nbon,2,1)
            dbopi2dc(nbon,3,1)=dbopi2dr(nbon)*drdc(nbon,3,1)
            dbopi2dc(nbon,1,2)=dbopi2dr(nbon)*drdc(nbon,1,2)
            dbopi2dc(nbon,2,2)=dbopi2dr(nbon)*drdc(nbon,2,2)
            dbopi2dc(nbon,3,2)=dbopi2dr(nbon)*drdc(nbon,3,2)
            ia(i1,2)=ia(i1,2)+1
            if (i1 /= i2) ia(i2,2)=ia(i2,2)+1
            ia(i1,ia(i1,2)+2)=i2
            ia(i2,ia(i2,2)+2)=i1
            if (ia(i1,2) > nsbma2) nsbma2=ia(i1,2)
            if (ia(i2,2) > nsbma2) nsbma2=ia(i2,2)
            if (bor > cutof3) then
                iag(i1,2)=iag(i1,2)+1
                iag(i2,2)=iag(i2,2)+1
                iag(i1,iag(i1,2)+2)=i2
                iag(i2,iag(i2,2)+2)=i1
                nubon1(i1,iag(i1,2))=nbon
                nubon1(i2,iag(i2,2))=nbon
                if (iag(i1,2) > nsbmax) nsbmax=iag(i1,2)
                if (iag(i2,2) > nsbmax) nsbmax=iag(i2,2)
            end if
            if (abs(de1(ib(nbon,1))) > -0.01d0) then
                nubon2(i1,ia(i1,2))=nbon
                nubon2(i2,ia(i2,2))=nbon
                nubon(i1,i2)=nbon     !Only valid for non-periodic systems
                nubon(i2,i1)=nbon
            else
                nbon=nbon-1      !Inorganics
            end if
        end if
    10 END DO

    if (nbon > nbomax) then
        write (95,*)nbon,nbomax
        stop 'Too many bonds; maybe wrong cell parameters'
    end if
!*********************************************************************
!                                                                    *
!     Sort molecules                                                 *
!                                                                    *
!*********************************************************************
    if (imolde == 1) return    !fixed molecular definitions

    FOUND= .FALSE. 
    DO 31 K1=1,NA
        IF (IA(K1,3+mbond) == 0) FOUND= .TRUE. 
        IF (IA(K1,3+mbond) > NMOLO) NMOLO=IA(K1,3+mbond)
    31 END DO
    IF ( .NOT. FOUND) GOTO 32
!***********************************************************************
!                                                                      *
!     Molecule numbers are assigned. No restrictions are made for the  *
!     sequence of the numbers in the connection table.                 *
!                                                                      *
!***********************************************************************
    N3=1
    34 N2=N3
    NMOLO=NMOLO+1
    if (nmolo > nmolmax) then
        write (*,*)nmolmax
        write (*,*)'Too many molecules in system; increase nmolmax'
        stop 'Too many molecules in system'
    end if
    IA(N2,3+mbond)=NMOLO
    37 FOUND= .FALSE. 
    DO 36 N1=N2+1,NA
        IF (IA(N1,3+mbond) /= 0.0d0) GOTO 36
        DO 35 L=1,mbond
            IF (IA(N1,l+2) == 0.0d0) GOTO 36
            IF (IA(IA(N1,l+2),3+mbond) == NMOLO) THEN
                FOUND= .TRUE. 
                IA(N1,3+mbond)=NMOLO
                GOTO 36
            ENDIF
        35 END DO
    36 END DO
    IF (FOUND) GOTO 37
    DO 33 N3=N2+1,NA
        IF (IA(N3,3+mbond) == 0) GOTO 34
    33 END DO
!***********************************************************************
!                                                                      *
!     The assigned or input molecule numbers are checked for their     *
!     consistency.                                                     *
!                                                                      *
!***********************************************************************
    32 FOUND= .FALSE. 
    DO 42 N1=1,NA
        DO 41 L=1,mbond
            IF (IA(N1,L+2) == 0) GOTO 42
            IF (IA(IA(N1,L+2),3+mbond) /= IA(N1,3+mbond)) THEN
                FOUND= .TRUE. 
            ENDIF
        41 END DO
    42 END DO
    IF (FOUND) THEN
        write (7,1000)na,qmol
        do i1=1,na
            write (7,1100)i1,ia(i1,1),(ia(i1,2+i2),i2=1,nsbmax), &
            ia(i1,3+mbond)
        end do
        write (7,*)tm11,tm22,tm33,angle(1),angle(2),angle(3)
        STOP 'Mol.nrs. not consistent; maybe wrong cell parameters'
    end if
!*********************************************************************
!                                                                    *
!     Sort molecules again                                           *
!                                                                    *
!*********************************************************************
    FOUND= .FALSE. 
    DO 61 K1=1,NA
        IF (IAG(K1,3+mbond) == 0) FOUND= .TRUE. 
        IF (IAG(K1,3+mbond) > NMOLO5) NMOLO5=IAG(K1,3+mbond)
    61 END DO
    IF ( .NOT. FOUND) GOTO 62
!***********************************************************************
!                                                                      *
!     Molecule numbers are assigned. No restrictions are made for the  *
!     sequence of the numbers in the connection table.                 *
!                                                                      *
!***********************************************************************
    N3=1
    64 N2=N3
    NMOLO5=NMOLO5+1
    if (nmolo5 > nmolmax) stop 'Too many molecules in system'
    IAG(N2,3+mbond)=NMOLO5
    67 FOUND= .FALSE. 
    DO 66 N1=N2+1,NA
        IF (IAG(N1,3+mbond) /= 0) GOTO 66
        DO 65 L=1,mbond
            IF (IAG(N1,l+2) == 0) GOTO 66
            IF (IAG(IAG(N1,l+2),3+mbond) == NMOLO5) THEN
                FOUND= .TRUE. 
                IAG(N1,3+mbond)=NMOLO5
                GOTO 66
            ENDIF
        65 END DO
    66 END DO
    IF (FOUND) GOTO 67
    DO 63 N3=N2+1,NA
        IF (IAG(N3,3+mbond) == 0) GOTO 64
    63 END DO
!***********************************************************************
!                                                                      *
!     The assigned or input molecule numbers are checked for their     *
!     consistency.                                                     *
!                                                                      *
!***********************************************************************
    62 FOUND= .FALSE. 
    DO 72 N1=1,NA
        DO 71 L=1,mbond
            IF (IAG(N1,L+2) == 0) GOTO 72
            IF (IAG(IAG(N1,L+2),3+mbond) /= IAG(N1,3+mbond)) THEN
                FOUND= .TRUE. 
            ENDIF
        71 END DO
    72 END DO
    IF (FOUND) THEN
        write (7,1000)na,qmol
        do i1=1,na
            write (7,1100)i1,iag(i1,1),(iag(i1,2+i2),i2=1,nsbmax), &
            iag(i1,3+mbond)
        end do
        write (7,*)tm11,tm22,tm33,angle(1),angle(2),angle(3)
        STOP 'Mol.nrs. not consistent; maybe wrong cell parameters'
    ENDIF

!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    1000 format (i3,2x,a60)
    1100 format (8i3)
    end subroutine srtbon1
!*********************************************************************
!*********************************************************************

    subroutine srtbon2

!*********************************************************************
    include 'cbka.blk'
    logical :: found
!*********************************************************************
!                                                                    *
!     Determine bond orders within the molecule from                 *
!     fixed connectivity                                             *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In srtbon2'
        call timer(65)
        close (65)
    end if
    do i1=1,na
        abo(i1)=0.0d0
    end do
    do 10 ib1=1,nbon
        i1=ib(ib1,2)
        i2=ib(ib1,3)
        itype=ib(ib1,1)
        call dista2(i1,i2,dis,dxm,dym,dzm)
        ih1=ia(i1,1)
        ih2=ia(i2,1)
        disdx=dxm/dis
        disdy=dym/dis
        disdz=dzm/dis

        rhulp=dis/rob1(ih1,ih2)
         
    !*********************************************************************
    !                                                                    *
    !     Determine bond orders                                          *
    !                                                                    *
    !*********************************************************************
        rh2=zero
        rh2p=zero
        rh2pp=zero
        ehulp=zero
        ehulpp=zero
        ehulppp=zero
        if (rapt(ih1) > zero .AND. rapt(ih2) > zero) then
            rhulp2=dis/rob2(ih1,ih2)
            rh2p=rhulp2**ptp(itype)
            ehulpp=exp(pdp(itype)*rh2p)
        end if
        if (vnq(ih1) > zero .AND. vnq(ih2) > zero) then
            rhulp3=dis/rob3(ih1,ih2)
            rh2pp=rhulp3**popi(itype)
            ehulppp=exp(pdo(itype)*rh2pp)
        end if

        if (rat(ih1) > zero .AND. rat(ih2) > zero) then
            rh2=rhulp**bop2(itype)
            ehulp=(1.0d0+cutoff)*exp(bop1(itype)*rh2)
        end if

        bor=ehulp+ehulpp+ehulppp
        borsi=ehulp
        borpi=ehulpp
        borpi2=ehulppp
        dbordrob=bop2(itype)*bop1(itype)*rh2*(1.0d0/dis)*ehulp+ &
        ptp(itype)*pdp(itype)*rh2p*(1.0d0/dis)*ehulpp+ &
        popi(itype)*pdo(itype)*rh2pp*(1.0d0/dis)*ehulppp
        dborsidrob=bop2(itype)*bop1(itype)*rh2*(1.0d0/dis)*ehulp
        dborpidrob=ptp(itype)*pdp(itype)*rh2p*(1.0d0/dis)*ehulpp
        dborpi2drob=popi(itype)*pdo(itype)*rh2pp*(1.0d0/dis)*ehulppp
             
    !*********************************************************************
    !                                                                    *
    !     Determine bond orders                                          *
    !                                                                    *
    !*********************************************************************
        drdc(ib1,1,1)=disdx
        drdc(ib1,2,1)=disdy
        drdc(ib1,3,1)=disdz
        drdc(ib1,1,2)=-disdx
        drdc(ib1,2,2)=-disdy
        drdc(ib1,3,2)=-disdz
        abo(i1)=abo(i1)+bor-cutoff
        abo(i2)=abo(i2)+bor-cutoff
        bo(ib1)=bor-cutoff
        bos(ib1)=bor-cutoff
        bosi(ib1)=borsi-cutoff
        bopi(ib1)=borpi
        bopi2(ib1)=borpi2
        rbo(ib1)=dis
        dbodr(ib1)=dbordrob
        dbopidr(ib1)=dborpidrob
        dbopi2dr(ib1)=dborpi2drob
        dbodc(ib1,1,1)=dbodr(ib1)*drdc(ib1,1,1)
        dbodc(ib1,2,1)=dbodr(ib1)*drdc(ib1,2,1)
        dbodc(ib1,3,1)=dbodr(ib1)*drdc(ib1,3,1)
        dbodc(ib1,1,2)=dbodr(ib1)*drdc(ib1,1,2)
        dbodc(ib1,2,2)=dbodr(ib1)*drdc(ib1,2,2)
        dbodc(ib1,3,2)=dbodr(ib1)*drdc(ib1,3,2)
        dbopidc(ib1,1,1)=dbopidr(ib1)*drdc(ib1,1,1)
        dbopidc(ib1,2,1)=dbopidr(ib1)*drdc(ib1,2,1)
        dbopidc(ib1,3,1)=dbopidr(ib1)*drdc(ib1,3,1)
        dbopidc(ib1,1,2)=dbopidr(ib1)*drdc(ib1,1,2)
        dbopidc(ib1,2,2)=dbopidr(ib1)*drdc(ib1,2,2)
        dbopidc(ib1,3,2)=dbopidr(ib1)*drdc(ib1,3,2)
        dbopi2dc(ib1,1,1)=dbopi2dr(ib1)*drdc(ib1,1,1)
        dbopi2dc(ib1,2,1)=dbopi2dr(ib1)*drdc(ib1,2,1)
        dbopi2dc(ib1,3,1)=dbopi2dr(ib1)*drdc(ib1,3,1)
        dbopi2dc(ib1,1,2)=dbopi2dr(ib1)*drdc(ib1,1,2)
        dbopi2dc(ib1,2,2)=dbopi2dr(ib1)*drdc(ib1,2,2)
        dbopi2dc(ib1,3,2)=dbopi2dr(ib1)*drdc(ib1,3,2)
    10 END DO

    return
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    1000 format (i3,2x,a60)
    1100 format (8i3)
    end subroutine srtbon2
!*********************************************************************
!*********************************************************************

    subroutine srtang

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),b(3),j(3)
    dimension ityva(100)
!*********************************************************************
!                                                                    *
!     Find valency angles in molecule                                *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In srtang'
        call timer(65)
        close (65)
    end if
    nval=0
    if (nvaty == 0) return
    do i1=1,nbon-1
        k4=ib(i1,2)
        k5=ib(i1,3)
        if (bo(i1) < cutof2) goto 51
        do i2=i1+1,nbon
            iju=0
            if (bo(i2) < cutof2) goto 50
            !! df398 hard coded cutoff might lead to discontinuity
            !! changed from original 0.00001d0 to cutof2*cutof2 to avoid discontinuity
            !! if (bo(i1)*bo(i2) < cutof2*cutof2) goto 50   !Scott Habershon

            k7=ib(i2,2)
            k8=ib(i2,3)

            if (k4 == k7 .AND. k5 == k8 .AND. k4 /= k8 .AND. k5 /= k7) then
                nval=nval+1
                iv(nval,2)=k5
                iv(nval,3)=k4
                iv(nval,4)=k8
                iv(nval,5)=i1
                iv(nval,6)=i2
                nval=nval+1
                iv(nval,2)=k4
                iv(nval,3)=k5
                iv(nval,4)=k7
                iv(nval,5)=i1
                iv(nval,6)=i2
                iju=2
            end if
            if (iju == 2) goto 50

            if (k4 == k8 .AND. k5 == k7 .AND. k4 /= k7 .AND. k5 /= k8) then
                nval=nval+1
                iv(nval,2)=k5
                iv(nval,3)=k4
                iv(nval,4)=k7
                iv(nval,5)=i1
                iv(nval,6)=i2
                nval=nval+1
                iv(nval,2)=k4
                iv(nval,3)=k5
                iv(nval,4)=k8
                iv(nval,5)=i1
                iv(nval,6)=i2
                iju=2
            end if
            if (iju == 2) goto 50

            if (k4 == k7) then
                nval=nval+1
                iv(nval,2)=k5
                iv(nval,3)=k4
                iv(nval,4)=k8
                iv(nval,5)=i1
                iv(nval,6)=i2
                iju=1
            end if
            if (iju == 1) goto 50

            if (k4 == k8) then
                nval=nval+1
                iv(nval,2)=k5
                iv(nval,3)=k4
                iv(nval,4)=k7
                iv(nval,5)=i1
                iv(nval,6)=i2
                iju=1
            end if
            if (iju == 1) goto 50

            if (k5 == k7) then
                nval=nval+1
                iv(nval,2)=k4
                iv(nval,3)=k5
                iv(nval,4)=k8
                iv(nval,5)=i1
                iv(nval,6)=i2
                iju=1
            end if
            if (iju == 1) goto 50

            if (k5 == k8) then
                nval=nval+1
                iv(nval,2)=k4
                iv(nval,3)=k5
                iv(nval,4)=k7
                iv(nval,5)=i1
                iv(nval,6)=i2
                iju=1
            end if
            50 continue
            if (iju > 0) then
            !*********************************************************************
            !                                                                    *
            !     Determine force field types of angles                          *
            !                                                                    *
            !*********************************************************************
                ityva(1)=0
                ityva(2)=0
                ityva(3)=0
                ityva(4)=0
                ityva(5)=0
                ih1=ia(iv(nval,2),1)
                ih2=ia(iv(nval,3),1)
                ih3=ia(iv(nval,4),1)
                if (ih3 < ih1) then
                    ih3=ia(iv(nval,2),1)
                    ih2=ia(iv(nval,3),1)
                    ih1=ia(iv(nval,4),1)
                end if

                nfound=0
                do i3=1,nvaty
                    if (ih1 == nvs(i3,1) .AND. ih2 == nvs(i3,2) .AND. &
                    ih3 == nvs(i3,3)) then
                        nfound=nfound+1
                        ityva(nfound)=i3
                    end if
                end do

                if (ityva(1) == 0 .OR. abs(vka(ityva(1))) < 0.001d0) then   !Valence angle does not exist in force field;ignore
                    nval=nval-1
                    ihul=0
                else
                    iv(nval,1)=ityva(1)
                    ihul=1

                    do i3=1,nfound-1           !Found multiple angles of the same type
                        nval=nval+1
                        iv(nval,1)=ityva(i3+1)
                        do i4=2,6
                            iv(nval,i4)=iv(nval-1,i4)
                        end do
                    end do

                end if

                if (iju == 2) then
                    ityva(1)=0
                    ih1=ia(iv(nval-ihul,2),1)
                    ih2=ia(iv(nval-ihul,3),1)
                    ih3=ia(iv(nval-ihul,4),1)
                    if (ih3 < ih1) then
                        ih3=ia(iv(nval-ihul,2),1)
                        ih2=ia(iv(nval-ihul,3),1)
                        ih1=ia(iv(nval-ihul,4),1)
                    end if

                    nfound=0
                    do i3=1,nvaty
                        if (ih1 == nvs(i3,1) .AND. ih2 == nvs(i3,2) .AND. &
                        ih3 == nvs(i3,3)) then
                            nfound=nfound+1
                            ityva(nfound)=i3
                        end if
                    end do

                    if (ityva(1) == 0 .OR. abs(vka(ityva(1))) < 0.001d0) then   !Valence angle does not exist in force field;ignore
                        if (ihul == 1) then
                            do i3=1,6
                                iv(nval-1,i3)=iv(nval,i3)
                            end do
                        end if
                        nval=nval-1
                    else
                        iv(nval-ihul,1)=ityva(1)

                        do i3=1,nfound-1             !Found multiple angles of the same type
                            nval=nval+1
                            iv(nval-ihul,1)=ityva(i3+1)
                            do i4=2,6
                                iv(nval-ihul,i4)=iv(nval-1-ihul,i4)
                            end do
                        end do

                    end if

                end if

            end if

        end do
        51 continue
    end do

    nbonop=0

    if (nval > nvamax) stop 'Too many valency angles'
         
    return
    end subroutine srtang
!*********************************************************************
!*********************************************************************

    subroutine srttor

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Find torsion angles in molecule                                *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In srttor'
        call timer(65)
        close (65)
    end if
    ntor=0
    if (ntoty == 0) return
    do i1=1,nbon
        k2=ib(i1,2)
        k3=ib(i1,3)
        iob1=ia(k2,2)
        iob2=ia(k3,2)
        do i2=1,iob1        !Atoms connected to k2
            k4=ia(k2,2+i2)
            ibo2=nubon2(k2,i2)
            do i3=1,iob2        !Atoms connected to k3
                k5=ia(k3,2+i3)
                ibo3=nubon2(k3,i3)
                bopr=bo(i1)*bo(ibo2)*bo(ibo3)
                !if (bopr > cutof2**3 .AND. k2 /= k5 .AND. k3 /= k4 .AND. k4 /= k5) then
                if (k2 /= k5 .AND. k3 /= k4 .AND. k4 /= k5) then
                    ntor=ntor+1
                    it(ntor,2)=k4
                    it(ntor,3)=k2
                    it(ntor,4)=k3
                    it(ntor,5)=k5
                    it(ntor,6)=ibo2
                    it(ntor,7)=i1
                    it(ntor,8)=ibo3
                !*********************************************************************
                !                                                                    *
                !     Determine force field types of torsion angles                  *
                !                                                                    *
                !*********************************************************************
                    ity=0
                    ih1=ia(it(ntor,2),1)
                    ih2=ia(it(ntor,3),1)
                    ih3=ia(it(ntor,4),1)
                    ih4=ia(it(ntor,5),1)

                    if (ih2 > ih3) then
                        ih1=ia(it(ntor,5),1)
                        ih2=ia(it(ntor,4),1)
                        ih3=ia(it(ntor,3),1)
                        ih4=ia(it(ntor,2),1)
                    end if

                    if (ih2 == ih3 .AND. ih4 < ih1) then
                        ih1=ia(it(ntor,5),1)
                        ih2=ia(it(ntor,4),1)
                        ih3=ia(it(ntor,3),1)
                        ih4=ia(it(ntor,2),1)
                    end if

                    do i4=1,ntoty
                        if (ih1 == nts(i4,1) .AND. ih2 == nts(i4,2) .AND. ih3 == nts(i4,3) &
                         .AND. ih4 == nts(i4,4)) ity=i4
                    end do

                    if (ity == 0) then
                        do i4=1,ntoty
                            if (nts(i4,1) == 0 .AND. ih2 == nts(i4,2) .AND. ih3 == nts(i4,3) &
                             .AND. nts(i4,4) == 0) ity=i4
                        end do
                    end if

                    if (ity == 0) then
                        ntor=ntor-1           !Torsion angle does not exist in force field: ignore
                    else
                        it(ntor,1)=ity
                    end if

                end if
            end do
        end do
    end do
    if (ntor > ntomax) stop 'Too many torsion angles'
!     do i1=1,ntor
!     write (41,'(20i4)')i1,it(i1,1),it(i1,2),it(i1,3),
!    $it(i1,4),it(i1,5),it(i1,6),it(i1,7),it(i1,8)
!     end do

    return
    end subroutine srttor
!*********************************************************************
!*********************************************************************

    subroutine srtoop

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In srtoop'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Find out of plane angles in molecule                           *
!                                                                    *
!*********************************************************************
    noop=0
    do i1=1,nval
        k2=iv(i1,2)
        k3=iv(i1,3)
        k4=iv(i1,4)
        k5=iv(i1,5)
        k6=iv(i1,6)
        do i2=1,nbon
            k7=ib(i2,2)
            k8=ib(i2,3)
            if (bo(i2) > cutof2) then
                if (k7 == k3 .AND. k8 /= k4 .AND. k8 /= k2) then
                    noop=noop+1
                    ioop(noop,2)=k8
                    ioop(noop,3)=k3
                    ioop(noop,4)=k2
                    ioop(noop,5)=k4
                    ioop(noop,6)=i2
                    ioop(noop,7)=iv(i1,5)
                    ioop(noop,8)=iv(i1,6)
                    ioop(noop,9)=i1
                end if
                if (k8 == k3 .AND. k7 /= k4 .AND. k7 /= k2) then
                    noop=noop+1
                    ioop(noop,2)=k7
                    ioop(noop,3)=k3
                    ioop(noop,4)=k2
                    ioop(noop,5)=k4
                    ioop(noop,6)=i2
                    ioop(noop,7)=iv(i1,5)
                    ioop(noop,8)=iv(i1,6)
                    ioop(noop,9)=i1
                end if
            end if
        end do
    end do
          
    do i1=1,noop
        call caltor(ioop(i1,2),ioop(i1,3),ioop(i1,4),ioop(i1,5),hoop)
    end do
          
!*********************************************************************
    return
    end subroutine srtoop
!*********************************************************************
!*********************************************************************

    subroutine srthb

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Find hydrogen bonds in molecule                                *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In srthb'
        call timer(65)
        close (65)
    end if
    nhb=0
!*********************************************************************
!                                                                    *
!     Locate donor/acceptor bonds                                    *
!                                                                    *
!*********************************************************************
    do 20 ivl=1,nvpair-nvlself      !Use Verlet-list to find donor-acceptor pairs
        j1=nvl1(ivl)
        j2=nvl2(ivl)
        k1=nvlx(ivl)
        k2=nvly(ivl)
        k3=nvlz(ivl)
        ity1=ia(j1,1)
        ity2=ia(j2,1)
        ihhb1=nphb(ia(j1,1))
        ihhb2=nphb(ia(j2,1))

        if (ihhb1 > ihhb2) then        !Make j1 donor atom and j2 acceptor atom
            j2=nvl1(ivl)
            j1=nvl2(ivl)
            ity1=ia(j1,1)
            ity2=ia(j2,1)
            ihhb1=nphb(ia(j1,1))
            ihhb2=nphb(ia(j2,1))
            k1=-nvlx(ivl)
            k2=-nvly(ivl)
            k3=-nvlz(ivl)
        end if

        if (ihhb1 == 1 .AND. ihhb2 == 2) then
            dx=c(j1,1)-c(j2,1)+k1*tm11
            dy=c(j1,2)-c(j2,2)+k1*tm21+k2*tm22
            dz=c(j1,3)-c(j2,3)+k1*tm31+k2*tm32+k3*tm33
            dishb=sqrt(dx*dx+dy*dy+dz*dz)
            if (dishb < 7.50d0) then        !df398 hydrogen bond cutoff (Angstrom)
                do 10 i23=1,ia(j1,2)        !Search for acceptor atoms bound to donor atom
                    j3=ia(j1,2+i23)
                    ity3=ia(j3,1)
                    nbohb=nubon(j1,j3)
                    if (nphb(ity3) == 2 .AND. j3 /= j2 .AND. bo(nbohb) > 0.01d0) then  !!! df398 bond order cutoff in HB-acceptor 
                    !*********************************************************************
                    !                                                                    *
                    !     Accept hydrogen bond and find hydrogen bond type               *
                    !                                                                    *
                    !*********************************************************************
                        nhb=nhb+1

                        if (nhb > nhbmax) then
                            write (*,*)nhb,nhbmax
                            write (*,*)'Maximum number of hydrogen bonds exceeded'
                            stop 'Maximum number of hydrogen bonds exceeded'
                        end if

                        ihb(nhb,1)=0

                        do i3=1,nhbty
                            if (ity3 == nhbs(i3,1) .AND. ity1 == nhbs(i3,2) .AND. ity2 == &
                            nhbs(i3,3)) ihb(nhb,1)=i3
                        end do

                        if (ihb(nhb,1) == 0) then    !Hydrogen bond not in force field
                            nhb=nhb-1
                        !     write (*,*)'Warning: added hydrogen bond ',ity3,ity1,ity2
                        !     nhbty=nhbty+1
                        !     nhbs(nhbty,1)=ity3
                        !     nhbs(nhbty,2)=ity1
                        !     nhbs(nhbty,3)=ity2
                        !     rhb(nhbty)=2.70
                        !     dehb(nhbty)=zero
                        !     vhb1(nhbty)=5.0
                        !     vhb2(nhbty)=20.0
                        !     ihb(nhb,1)=nhbty
                        end if

                        ihb(nhb,2)=j3
                        ihb(nhb,3)=j1
                        ihb(nhb,4)=j2
                        ihb(nhb,5)=nbohb
                        ihb(nhb,6)=k1
                        ihb(nhb,7)=k2
                        ihb(nhb,8)=k3
                    end if
                10 END DO

            end if
        end if
20 END DO

!     stop 'end in srthb'
    return
    end subroutine srthb
!*********************************************************************
!*********************************************************************

    subroutine testaf

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
    dimension esav(nat,3,2),ebas1(nat,3,2),ebas2(nat,3,2), &
    ebas3(nat,3,2),emole(nat,3,2), &
    evs(nat,3,2),ehbs(nat,3,2),etcos(nat,3,2),enbs(nat,3,2), &
    erads(nat,3,2), &
    dedc(3,nat),ds(3,nat),demodc(3,nat), &
    deb1adc(3,nat),deb2adc(3,nat),deb3adc(3,nat),debtadc(3,nat), &
    devdc(3,nat),dehbdc(3,nat),detcodc(3,nat),denbdc(3,nat), &
    deraddc(3,nat),dass(3,nat), &
    idum(mbond+3),bodum(mbond+3)
!*********************************************************************
!                                                                    *
!     Test first derivative calculation                              *
!     Lone pairs are not taken into account                          *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In testaf'
        call timer(65)
        close (65)
    end if
    do i1=1,mbond+3
        idum(i1)=nzero
        bodum(i1)=zero
    end do
!     do i1=1,nso
!     vlp1(i1)=0.0
!     vlp2(i1)=0.0
!     end do
    call distan
    call vlist
    call srtbon1
    call inilp
    call intcor
    call chargess

    write (7,100)na,qmol
    do i1=1,na
        bosum=0.0d0
        do i3=1,nsbmax
            if (iag(i1,2+i3) > 0) bosum=bosum+bo(nubon2(i1,i3))
        end do
        if (nsbmax < 5) then
            write (7,200)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,5-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,5-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 10) then
            write (7,210)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,10-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,10-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 15) then
            write (7,220)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,15-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,15-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 20) then
            write (7,230)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,20-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,20-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 25) then
            write (7,240)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,25-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,25-iag(i1,2)),bosum,vlp(i1),ch(i1)
        end if
    end do

    estrc=0.0d0
    do i1=1,na
        do i2=1,3
            d(i2,i1)=0.0d0
            ds(i2,i1)=0.0d0
        end do
    end do
    call boncor
    call srtang
    call srttor
    call srthb
    call calval
    call covbon
    write (41,*)'Bonds'
    do i1=1,na
        write (41,'(i3,3f12.6)')i1,(d(i2,i1),i2=1,3)
        do i2=1,3
            ds(i2,i1)=d(i2,i1)
        end do
    end do
    call lonpar
    write (49,*)'Lone pair'
    do i1=1,na
        write (49,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
        do i2=1,3
            ds(i2,i1)=d(i2,i1)
        end do
    end do
    call ovcor
    write (48,*)'Over/undercoordination'
    do i1=1,na
        write (48,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
        do i2=1,3
            ds(i2,i1)=d(i2,i1)
        end do
    end do
    call molen
    write (42,*)'Molecular'
    do i1=1,na
        write (42,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
        do i2=1,3
            ds(i2,i1)=d(i2,i1)
        end do
    end do
    call valang
    write (43,*)'Valency angles+valency conjugation'
    do i1=1,na
        write (43,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
        do i2=1,3
            ds(i2,i1)=d(i2,i1)
        end do
    end do
    call hbond
    write (44,*)'Hydrogen bonds'
    do i1=1,na
        write (44,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
        do i2=1,3
            ds(i2,i1)=d(i2,i1)
        end do
    end do
    call torang
    write (45,*)'Torsion angles+conjugation'
    do i1=1,na
        write (45,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
        do i2=1,3
            ds(i2,i1)=d(i2,i1)
        end do
    end do
    if (nrdd == 1) then
        call radbo
        write (47,*)'Radical/double bonds'
        do i1=1,na
            write (47,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
            do i2=1,3
                ds(i2,i1)=d(i2,i1)
            end do
        end do
    end if
    call nonbon
    write (46,*)'Nonbonded'
    do i1=1,na
        write (46,'(i3,3f12.6)')i1,(d(i2,i1)-ds(i2,i1),i2=1,3)
    end do
    write (40,*)'Total'
    do i1=1,na
        write (40,'(i3,3f12.6)')i1,(d(i2,i1),i2=1,3)
    end do
    do i1=1,na
        do i2=1,3
            dass(i2,i1)=d(i2,i1)
        end do
    end do

    write (7,100)na,qmol
    do i1=1,na
        bosum=0.0d0
        do i3=1,nsbmax
            if (nubon2(i1,i3) > 0) bosum=bosum+bo(nubon2(i1,i3))
        end do
        if (nsbmax < 5) then
            write (7,200)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,5-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,5-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 10) then
            write (7,210)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,10-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,10-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 15) then
            write (7,220)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,15-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,15-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 20) then
            write (7,230)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,20-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,20-iag(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbmax < 25) then
            write (7,240)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
            (idum(i2),i2=1,25-iag(i1,2)), &
            iag(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,iag(i1,2)), &
            (bodum(i2),i2=1,25-iag(i1,2)),bosum,vlp(i1),ch(i1)
        end if
    end do

!     stop
    open (8,file='fort.8',status='unknown')
    write (8,100)na,qmol
    do i1=1,na
        bosum=0.0d0
        do i3=1,nsbma2
            if (nubon2(i1,i3) > 0) bosum=bosum+bo(nubon2(i1,i3))
        end do
        if (nsbma2 < 5) then
            write (8,200)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
            (idum(i2),i2=1,5-ia(i1,2)), &
            ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
            (bodum(i2),i2=1,5-ia(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbma2 < 10) then
            write (8,210)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
            (idum(i2),i2=1,10-ia(i1,2)), &
            ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
            (bodum(i2),i2=1,10-ia(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbma2 < 15) then
            write (8,220)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
            (idum(i2),i2=1,15-ia(i1,2)), &
            ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
            (bodum(i2),i2=1,15-ia(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbma2 < 20) then
            write (8,230)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
            (idum(i2),i2=1,20-ia(i1,2)), &
            ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
            (bodum(i2),i2=1,20-ia(i1,2)),bosum,vlp(i1),ch(i1)
        else if (nsbma2 < 25) then
            write (8,240)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
            (idum(i2),i2=1,25-ia(i1,2)), &
            ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
            (bodum(i2),i2=1,25-ia(i1,2)),bosum,vlp(i1),ch(i1)
        end if
    end do
    close(8)

    write (7,*)'Nr. of bonds:',nbon
    write (7,*)'Nr. of valency angles:',nval
    write (7,*)'Nr. of torsion angles:',ntor
    write (7,*)'Nr. of hydrogen bonds:',nhb
    write (71,150)
    write (71,'(12f14.4)')eb,ea,elp,emol,ev,ecoa,ehb,et,ehb, &
    ew,ep,eradbo
    if (ntest == 2) then
        call mdsav(0,qfile(nprob))
        stop 'Single point'
    end if
!*********************************************************************
!                                                                    *
!     Numeric calculation of first derivatives for cell parameters   *
!                                                                    *
!*********************************************************************
    if (ntest == 3) then

        write (90,*)'Total'
        write (91,*)'Bonds'
        write (93,*)'Valency angles+valency conjugation'
        write (94,*)'Hydrogen bonds'
        write (95,*)'Torsion angles+conjugation'
        write (96,*)'Nonbonded'
        write (97,*)'Bonds+atoms+lone pairs'
        write (98,*)'Atoms'
        write (99,*)'Lone pairs'

        do i1=1,na
            do i2=1,3
                d(i2,i1)=0.0
                do i3=1,27
                    dcell(i2,i1,i3)=0.0d0
                end do
            end do
        end do
    !*********************************************************************
    !                                                                    *
    !     Calculate minimum images for NPT                               *
    !                                                                    *
    !*********************************************************************
        do i1=1,na-1
            do  i2=i1+1,na
                dx1=c(i1,1)-c(i2,1)
                dy1=c(i1,2)-c(i2,2)
                dz1=c(i1,3)-c(i2,3)
                dismin=1d+10
                do  k1=-1,1
                    do  k2=-1,1
                        do  k3=-1,1
                            a1=dx1+k1*tm11
                            a2=dy1+k1*tm21+k2*tm22
                            a3=dz1+k1*tm31+k2*tm32+k3*tm33
                            rr=sqrt(a1*a1+a2*a2+a3*a3)
                            if (rr < dismin) then
                                dismin=rr
                                ixsav=k1
                                iysav=k2
                                izsav=k3
                            end if
                        enddo
                END DO
                enddo
                nmpx(i1,i2)=ixsav
                nmpx(i2,i1)=-ixsav
                nmpy(i1,i2)=iysav
                nmpy(i2,i1)=-iysav
                nmpz(i1,i2)=izsav
                nmpz(i2,i1)=-izsav
          END DO
        enddo
         
        do i1=1,na
            nmpx(i1,i1)=0
            nmpy(i1,i1)=0
            nmpz(i1,i1)=0
        end do

        call distan
        call vlist
        call intcor
        call boncor
        call covbon
        call calcpres
        write (91,'(3f20.6)')(deda(i1),i1=1,3)
        deda1s=deda(1)
        deda2s=deda(2)
        deda3s=deda(3)
        call lonpar
        call calcpres
        write (99,'(3f20.6)')deda(1)-deda1s,deda(2)-deda2s, &
        deda(3)-deda3s
        deda1s=deda(1)
        deda2s=deda(2)
        deda3s=deda(3)
        call ovcor
        call calcpres
        write (98,'(3f20.6)')deda(1)-deda1s,deda(2)-deda2s, &
        deda(3)-deda3s
        write (97,'(3f20.6)')deda(1),deda(2),deda(3)
        deda1s=deda(1)
        deda2s=deda(2)
        deda3s=deda(3)
        call calval
        call valang
        call calcpres
        write (93,'(3f20.6)')deda(1)-deda1s,deda(2)-deda2s, &
        deda(3)-deda3s
        deda1s=deda(1)
        deda2s=deda(2)
        deda3s=deda(3)
        call hbond
        call calcpres
        write (94,'(3f20.6)')deda(1)-deda1s,deda(2)-deda2s, &
        deda(3)-deda3s
        deda1s=deda(1)
        deda2s=deda(2)
        deda3s=deda(3)
        call torang
        call calcpres
        write (95,'(3f20.6)')deda(1)-deda1s,deda(2)-deda2s, &
        deda(3)-deda3s
        deda1s=deda(1)
        deda2s=deda(2)
        deda3s=deda(3)
        call nonbon
        call calcpres
        write (90,'(3f20.6)')deda(1),deda(2),deda(3)
        write (96,'(3f20.6)')deda(1)-deda1s,deda(2)-deda2s, &
        deda(3)-deda3s
        deda1s=deda(1)
        deda2s=deda(2)
        deda3s=deda(3)

        delta=0.00001d0
        axis1=axis(1)
        axis2=axis(2)
        axis3=axis(3)
        angle1=angle(1)
        angle2=angle(2)
        angle3=angle(3)
        do i1=1,3
            axis(i1)=axis(i1)+delta
            call distan
            call vlist
            call intcor
            call boncor
            call covbon
            call lonpar
            call ovcor
            call calval
            call valang
            call hbond
            call torang
            call nonbon
            write (71,'(12f14.4)')eb,ea,elp,emol,ev,ecoa,ehb,et,ehb, &
            ew,ep,eradbo
            esav(i1,1,1)=eb+ea+elp+emol+ev+ecoa+epen+ehb+et+eco+ew+ep+ &
            eradbo
            ebas1(i1,1,1)=eb
            ebas2(i1,1,1)=ea
            ebas3(i1,1,1)=elp
            evs(i1,1,1)=ev+ecoa
            ehbs(i1,1,1)=ehb
            etcos(i1,1,1)=et+eco
            enbs(i1,1,1)=ew+ep

            axis(i1)=axis(i1)-2.0d0*delta
            call distan
            call vlist
            call intcor
            call boncor
            call covbon
            call lonpar
            call ovcor
            call calval
            call valang
            call hbond
            call torang
            call nonbon
            write (71,'(12f14.4)')eb,ea,elp,emol,ev,ecoa,ehb,et,ehb, &
            ew,ep,eradbo
            esav(i1,1,2)=eb+ea+elp+emol+ev+ecoa+epen+ehb+et+eco+ew+ep+ &
            eradbo
            ebas1(i1,1,2)=eb
            ebas2(i1,1,2)=ea
            ebas3(i1,1,2)=elp
            evs(i1,1,2)=ev+ecoa
            ehbs(i1,1,2)=ehb
            etcos(i1,1,2)=et+eco
            enbs(i1,1,2)=ew+ep
            axis(i1)=axis(i1)+delta

            dedc(i1,1)=(esav(i1,1,1)-esav(i1,1,2))/(2.0d0*delta)
            deb1adc(i1,1)=(ebas1(i1,1,1)-ebas1(i1,1,2))/(2.0d0*delta)
            deb2adc(i1,1)=(ebas2(i1,1,1)-ebas2(i1,1,2))/(2.0d0*delta)
            deb3adc(i1,1)=(ebas3(i1,1,1)-ebas3(i1,1,2))/(2.0d0*delta)
            debtadc(i1,1)=deb1adc(i1,1)+deb2adc(i1,1)+deb3adc(i1,1)
            devdc(i1,1)=(evs(i1,1,1)-evs(i1,1,2))/(2.0d0*delta)
            dehbdc(i1,1)=(ehbs(i1,1,1)-ehbs(i1,1,2))/(2.0d0*delta)
            detcodc(i1,1)=(etcos(i1,1,1)-etcos(i1,1,2))/(2.0d0*delta)
            denbdc(i1,1)=(enbs(i1,1,1)-enbs(i1,1,2))/(2.0d0*delta)
        end do

        write (90,'(3f20.6)')(dedc(i1,1),i1=1,3)
        write (91,'(3f20.6)')(deb1adc(i1,1),i1=1,3)
        write (93,'(3f20.6)')(devdc(i1,1),i1=1,3)
        write (94,'(3f20.6)')(dehbdc(i1,1),i1=1,3)
        write (95,'(3f20.6)')(detcodc(i1,1),i1=1,3)
        write (96,'(3f20.6)')(denbdc(i1,1),i1=1,3)
        write (97,'(3f20.6)')(debtadc(i1,1),i1=1,3)
        write (98,'(3f20.6)')(deb2adc(i1,1),i1=1,3)
        write (99,'(3f20.6)')(deb3adc(i1,1),i1=1,3)

        write (*,*) &
        'Numerical calculation of cell parameter derivatives'
        stop 'Numerical calculation of cell parameter derivatives'
    end if
!*********************************************************************
!                                                                    *
!     Numeric calculation of first derivatives                       *
!                                                                    *
!*********************************************************************
    delta=0.00001d0
        
    do i1=1,na
        do i2=1,3
                  
            c(i1,i2)=c(i1,i2)+delta
            call intcor
            call calval
            call boncor
            call lonpar
            call covbon
            call ovcor
            call srttor
            call molen
            call valang
            call hbond
            call torang
            call nonbon
            eradbo=zero
            if (nrdd == 1) call radbo
            esav(i1,i2,1)=eb+ea+elp+emol+ev+ecoa+epen+ehb+et+eco+ew+ep+ &
            eradbo
            ebas1(i1,i2,1)=eb
            ebas2(i1,i2,1)=ea
            ebas3(i1,i2,1)=elp
            emole(i1,i2,1)=emol
            evs(i1,i2,1)=ev+ecoa
            ehbs(i1,i2,1)=ehb
            etcos(i1,i2,1)=et+eco
            enbs(i1,i2,1)=ew+ep
            erads(i1,i2,1)=eradbo
            c(i1,i2)=c(i1,i2)-2.0d0*delta
            call intcor
            call calval
            call boncor
            call lonpar
            call covbon
            call ovcor
            call srttor
            call molen
            call valang
            call hbond
            call torang
            call nonbon
            eradbo=zero
            if (nrdd == 1) call radbo
            esav(i1,i2,2)=eb+ea+elp+emol+ev+ecoa+epen+ehb+et+eco+ew+ep+ &
            eradbo
            ebas1(i1,i2,2)=eb
            ebas2(i1,i2,2)=ea
            ebas3(i1,i2,2)=elp
            emole(i1,i2,2)=emol
            evs(i1,i2,2)=ev+ecoa
            ehbs(i1,i2,2)=ehb
            etcos(i1,i2,2)=et+eco
            enbs(i1,i2,2)=ew+ep
            erads(i1,i2,2)=eradbo
            dedc(i2,i1)=(esav(i1,i2,1)-esav(i1,i2,2))/(2.0d0*delta)
            deb1adc(i2,i1)=(ebas1(i1,i2,1)-ebas1(i1,i2,2))/(2.0d0*delta)
            deb2adc(i2,i1)=(ebas2(i1,i2,1)-ebas2(i1,i2,2))/(2.0d0*delta)
            deb3adc(i2,i1)=(ebas3(i1,i2,1)-ebas3(i1,i2,2))/(2.0d0*delta)
            demodc(i2,i1)=(emole(i1,i2,1)-emole(i1,i2,2))/(2.0d0*delta)
            debtadc(i2,i1)=deb1adc(i2,i1)+deb2adc(i2,i1)+deb3adc(i2,i1)
            devdc(i2,i1)=(evs(i1,i2,1)-evs(i1,i2,2))/(2.0d0*delta)
            dehbdc(i2,i1)=(ehbs(i1,i2,1)-ehbs(i1,i2,2))/(2.0d0*delta)
            detcodc(i2,i1)=(etcos(i1,i2,1)-etcos(i1,i2,2))/(2.0d0*delta)
            denbdc(i2,i1)=(enbs(i1,i2,1)-enbs(i1,i2,2))/(2.0d0*delta)
            deraddc(i2,i1)=(erads(i1,i2,1)-erads(i1,i2,2))/(2.0d0*delta)
            c(i1,i2)=c(i1,i2)+delta
        end do
    end do
! df398
!d = dedc ! equate analytical derivatives to numeric derivatives to sanity check via icheck=1 in control
    write (90,*)'Total'
    write (51,*)'Differences'
    write (91,*)'Bonds'
    write (92,*)'Molecular'
    write (93,*)'Valency angles+valency conjugation'
    write (94,*)'Hydrogen bonds'
    write (95,*)'Torsion angles+conjugation'
    write (96,*)'Nonbonded'
    write (97,*)'Bonds+atoms+lone pairs'
    write (98,*)'Atoms'
    write (99,*)'Lone pairs'
!     write (100,*)'Radical/double bonds'
    do i1=1,na
        write (90,'(i3,3f12.6)')i1,(dedc(i2,i1),i2=1,3)
        write (51,'(i3,3f12.6)')i1,(dedc(i2,i1)-dass(i2,i1),i2=1,3)
        write (91,'(i3,3f12.6)')i1,(deb1adc(i2,i1),i2=1,3)
        write (92,'(i3,3f12.6)')i1,(demodc(i2,i1),i2=1,3)
        write (93,'(i3,3f12.6)')i1,(devdc(i2,i1),i2=1,3)
        write (94,'(i3,3f12.6)')i1,(dehbdc(i2,i1),i2=1,3)
        write (95,'(i3,3f12.6)')i1,(detcodc(i2,i1),i2=1,3)
        write (96,'(i3,3f12.6)')i1,(denbdc(i2,i1),i2=1,3)
        write (97,'(i3,3f12.6)')i1,(debtadc(i2,i1),i2=1,3)
        write (98,'(i3,3f12.6)')i1,(deb2adc(i2,i1),i2=1,3)
        write (99,'(i3,3f12.6)')i1,(deb3adc(i2,i1),i2=1,3)
    !     write (100,'(i3,3f12.6)')i1,(deraddc(i2,i1),i2=1,3)
    end do
    return
    100 format (i3,2x,a60)
    150 format (6x,'Ebond',11x,'Eatom',9x,'Elopa',10x,'Emol',10x,'Eval' &
    ,10x,'Ecoa',10x,'Ehb',11x,'Etors',9x,'Econj',9x,'Evdw',10x, &
    'Epp')
    200 format (8i4,8f7.3)
    210 format (13i4,13f7.3)
    220 format (18i4,18f7.3)
    230 format (23i4,23f7.3)
    240 format (28i4,28f7.3)
    end subroutine testaf

!**********************************************************************
!**********************************************************************

    SUBROUTINE TRAROT1

!**********************************************************************
!                                                                     *
!    Remove translational or rotational kinetic energy.               *
!                                                                     *
!**********************************************************************

    INCLUDE 'cbka.blk'

    DIMENSION VMC(3),ANGMOM(3),TENSOR(3,3),ANGVEL(3),XMC(3)

    NSTART=1
    NEND=NA
!**********************************************************************
!                                                                     *
!     Compute the position, the velocity of and angular momentum about*
!     the centre of mass.                                             *
!                                                                     *
!**********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In trarot1'
        call timer(65)
        close (65)
    end if
    DO 1 K=1,3
        XMC(K)=0.0d0
        VMC(K)=0.0d0
        ANGMOM(K)=0.0d0
    1 END DO
    DO 3 K=NSTART,NEND
        XMASSI=XMASAT(K)
    !     cxh=c(k,1)-anint(c(k,1)/tm11)*tm11
    !     cyh=c(k,2)-anint(c(k,2)/tm22)*tm22
    !     czh=c(k,3)-anint(c(k,3)/tm33)*tm33
        cxh=c(k,1)
        cyh=c(k,2)
        czh=c(k,3)
        XMC(1)=XMC(1)+CXH*XMASSI
        XMC(2)=XMC(2)+CYH*XMASSI
        XMC(3)=XMC(3)+CZH*XMASSI
        DO 2 J=1,3
            VMC(J)=VMC(J)+VEL(J,K)*XMASSI
        2 END DO
        ANGMOM(1)=ANGMOM(1)+(CYH*VEL(3,K)-CZH*VEL(2,K))*XMASSI
        ANGMOM(2)=ANGMOM(2)+(CZH*VEL(1,K)-CXH*VEL(3,K))*XMASSI
        ANGMOM(3)=ANGMOM(3)+(CXH*VEL(2,K)-CYH*VEL(1,K))*XMASSI
    3 END DO
    EKTRAN=0.0
    DO 4 K=1,3
        XMC(K)=XMC(K)/XMASMD
        VMC(K)=VMC(K)/XMASMD
        EKTRAN=EKTRAN+VMC(K)*VMC(K)
    4 END DO
    EKTRAN=HALF*EKTRAN*XMASMD/CONVMD
    ANGMOM(1)=ANGMOM(1)-(XMC(2)*VMC(3)-XMC(3)*VMC(2))*XMASMD
    ANGMOM(2)=ANGMOM(2)-(XMC(3)*VMC(1)-XMC(1)*VMC(3))*XMASMD
    ANGMOM(3)=ANGMOM(3)-(XMC(1)*VMC(2)-XMC(2)*VMC(1))*XMASMD

!**********************************************************************
!                                                                     *
!     Calculate and then invert the inertial tensor.                  *
!                                                                     *
!**********************************************************************

    XX=0.0d0
    XY=0.0d0
    XZ=0.0d0
    YY=0.0d0
    YZ=0.0d0
    ZZ=0.0d0
    DO 5 I=NSTART,NEND
        XMASSI=XMASAT(I)
        XDIFF=C(I,1)-XMC(1)
        YDIFF=C(I,2)-XMC(2)
        ZDIFF=C(I,3)-XMC(3)
    !     xdiff=xdiff-anint(xdiff/tm11)*tm11
    !     ydiff=ydiff-anint(ydiff/tm22)*tm22
    !     zdiff=zdiff-anint(zdiff/tm33)*tm33
        XX=XX+XDIFF*XDIFF*XMASSI
        XY=XY+XDIFF*YDIFF*XMASSI
        XZ=XZ+XDIFF*ZDIFF*XMASSI
        YY=YY+YDIFF*YDIFF*XMASSI
        YZ=YZ+YDIFF*ZDIFF*XMASSI
        ZZ=ZZ+ZDIFF*ZDIFF*XMASSI
    5 END DO
    TENSOR(1,1)=YY+ZZ
    TENSOR(2,1)=-XY
    TENSOR(3,1)=-XZ
    TENSOR(1,2)=-XY
    TENSOR(2,2)=XX+ZZ
    TENSOR(3,2)=-YZ
    TENSOR(1,3)=-XZ
    TENSOR(2,3)=-YZ
    TENSOR(3,3)=XX+YY
    CALL INVMAT(3,3,TENSOR,3)

!**********************************************************************
!                                                                     *
!     Compute the angular velocity about the centre of mass.          *
!                                                                     *
!**********************************************************************

    EKROT=0.0d0
    DO  I=1,3
        ANGVEL(I)=0.0d0
        DO  J=1,3
            ANGVEL(I)=ANGVEL(I)+TENSOR(I,J)*ANGMOM(J)
        END DO
        EKROT=EKROT+ANGVEL(I)*ANGMOM(I)
    END DO
    EKROT=HALF*EKROT/CONVMD

!**********************************************************************
!                                                                     *
!     Eliminate rotation about the centre of mass.                    *
!                                                                     *
!**********************************************************************

    DO I=NSTART,NEND
        XDIFF=C(I,1)-XMC(1)
        YDIFF=C(I,2)-XMC(2)
        ZDIFF=C(I,3)-XMC(3)
    !     xdiff=xdiff-anint(xdiff/tm11)*tm11
    !     ydiff=ydiff-anint(ydiff/tm22)*tm22
    !     zdiff=zdiff-anint(zdiff/tm33)*tm33
        VEL(1,I)=VEL(1,I)-ANGVEL(2)*ZDIFF+ANGVEL(3)*YDIFF
        VEL(2,I)=VEL(2,I)-ANGVEL(3)*XDIFF+ANGVEL(1)*ZDIFF
        VEL(3,I)=VEL(3,I)-ANGVEL(1)*YDIFF+ANGVEL(2)*XDIFF
    END DO

!**********************************************************************
!                                                                     *
!     Eliminate translation of the centre of mass.                    *
!                                                                     *
!**********************************************************************

      DO I=1,NA
        DO  J=1,3
            VEL(J,I)=VEL(J,I)-VMC(J)
        END DO
      ENDDO

!**********************************************************************
!                                                                     *
!     Remove the translation occured during the last 100 steps by back*
!     translation of the entire sytem.                                *
!                                                                     *
!**********************************************************************

!     D1=CENTRE(1)-XMC(1)
!     D2=CENTRE(2)-XMC(2)
!     D3=CENTRE(3)-XMC(3)
!     DO 11 K=1,NA
!     C(K,1)=C(K,1)+D1
!     C(K,2)=C(K,2)+D2
!  11 C(K,3)=C(K,3)+D3

    return
    END SUBROUTINE TRAROT1

!**********************************************************************
!**********************************************************************

    SUBROUTINE TRAROT2

!**********************************************************************
    INCLUDE 'cbka.blk'
    dimension cmcent(nmolmax,3),natmol(nmolmax),fc(nat,3)
    dimension rang(3)
!**********************************************************************
!                                                                     *
!     Translate atoms back to unit cell                               *
!                                                                     *
!**********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In trarot2'
        call timer(65)
        close (65)
    end if
!**********************************************************************
!                                                                     *
!     Conversion to fractional coordinates                            *
!                                                                     *
!**********************************************************************
    do i1=1,na
        fc(i1,1)=c(i1,1)/tm11
        fc(i1,2)=(c(i1,2)-tm21*fc(i1,1))/tm22
        fc(i1,3)=(c(i1,3)-tm31*fc(i1,1)-tm32*fc(i1,2))/tm33
        rang(1)=range/tm11
        rang(2)=(range-tm21*rang(1))/tm22
        rang(3)=(range-tm31*rang(1)-tm32*rang(2))/tm22
        rangdis=(rang(1)+rang(2)+rang(3))/3.0d0
    end do

    if (itrans == 1 .OR. iortho == 0) then
    !**********************************************************************
    !                                                                     *
    !     Simple back-translation                                         *
    !                                                                     *
    !**********************************************************************
        do i1=1,na
            do i2=1,3
                if (fc(i1,i2) < -rang(i2)) then
                    itrans2=int(fc(i1,i2))
                    fc(i1,i2)=fc(i1,i2)+1.0d0-itrans2
                    id(i1,i2)=id(i1,i2)-1.0d0+itrans2
                end if
                if (fc(i1,i2) > 1.0d0+rang(i2)) then
                    itrans2=int(fc(i1,i2))
                    fc(i1,i2)=fc(i1,i2)-itrans2
                    id(i1,i2)=id(i1,i2)+itrans2
                end if
            end do
        end do

    end if

    if (itrans == 2 .OR. itrans == 3) then
    !**********************************************************************
    !                                                                     *
    !     Back translation based on molecular composition                 *
    !                                                                     *
    !**********************************************************************
        do i1=1,nmolo5
            natmol(i1)=0
            do i2=1,3
                cmcent(i1,i2)=0.0d0
            end do
        end do
        do i1=1,na
            imol=iag(i1,3+mbond)
            do i2=1,3
                cmcent(imol,i2)=cmcent(imol,i2)+fc(i1,i2)
            end do
            natmol(imol)=natmol(imol)+1
        end do
        do i1=1,nmolo5
            do i2=1,3
                cmcent(i1,i2)=cmcent(i1,i2)/float(natmol(i1))
            end do
        end do
        do i1=1,nmolo5
            do i2=1,3
                if (cmcent(i1,i2) > 1.0d0+rang(i2)) then
                !     open (58,file='fort.58',status='unknown',position='append')
                !     write (58,*)'Iteration',mdstep,' Moved molecule number',i1
                !     write (58,*)i2,-axis(i2)
                !     close (58)
                    do i3=1,na
                        if (iag(i3,3+mbond) == i1) then
                            fc(i3,i2)=fc(i3,i2)-1.0d0
                            id(i3,i2)=id(i3,i2)+1.0d0
                        end if
                    end do
                    cmcent(i1,i2)=cmcent(i1,i2)-1.0d0
                end if
                if (cmcent(i1,i2) < -rang(i2)) then
                !     open (58,file='fort.58',status='unknown',position='append')
                !     write (58,*)'Iteration',mdstep,' Moved molecule number',i1
                !     write (58,*)i2,axis(i2)
                !     close (58)
                    do i3=1,na
                        if (iag(i3,3+mbond) == i1) then
                            fc(i3,i2)=fc(i3,i2)+1.0d0
                            id(i3,i2)=id(i3,i2)-1.0d0
                        end if
                    end do
                    cmcent(i1,i2)=cmcent(i1,i2)+1.0d0
                end if
            end do
        end do
    !**********************************************************************
    !                                                                     *
    !     Translate atoms on same bond to the same side of the unit cell  *
    !                                                                     *
    !**********************************************************************
        if (itrans == 3) then
            do i1=1,nbon
                j1=ib(i1,2)
                j2=ib(i1,3)
                if (bo(i1) > cutof3) then
                    do i2=1,3
                        dis=fc(j1,i2)-fc(j2,i2)
                        if (int(2.0d0*(abs(dis)-rangdis)) /= 0) then
                            tra2=1.0d0
                            if (int(2.0d0*dis) < 0) tra2=-1.0d0
                            imol=iag(j1,3+mbond)
                            discen1=fc(j1,i2)-cmcent(imol,i2)
                            discen2=fc(j2,i2)-cmcent(imol,i2)
                            ihu1=int(2.0d0*discen1)
                            ihu2=int(2.0d0*discen2)
                             
                            if (ihu1 /= 0 .AND. ihu2 == 0) then
                            !     open (58,file='fort.58',status='unknown',position='append')
                            !     fc(j1,i2)=fc(j1,i2)-tra2
                            !     write (58,*)'Iteration',mdstep,' Moved atom number',j1
                            !     write (58,*)i2,j2,-tra2
                            !     close (58)
                            end if
                            if (ihu1 == 0 .AND. ihu2 /= 0) then
                            !     open (58,file='fort.58',status='unknown',position='append')
                            !     fc(j2,i2)=fc(j2,i2)+tra2
                            !     write (58,*)'Iteration',mdstep,' Moved atom number',j2
                            !     write (58,*)i2,j2,tra2
                            !     close (58)
                            end if
                            if (ihu1 /= 0 .AND. ihu2 /= 0) then
                            !     open (58,file='fort.58',status='unknown',position='append')
                            !     fc(j2,i2)=fc(j2,i2)+tra2
                            !     write (58,*)'Iteration',mdstep,' Moved atom number',j2
                            !     write (58,*)i2,j2,tra2
                            !     close (58)
                            end if
                             
                        end if
                    end do
                end if
            end do

        end if
    end if
!*********************************************************************
!                                                                    *
!     Conversion to cartesian coordinates                            *
!                                                                    *
!*********************************************************************
    do i1=1,na
        c(i1,1)=fc(i1,1)*tm11
        c(i1,2)=fc(i1,1)*tm21+fc(i1,2)*tm22
        c(i1,3)=fc(i1,1)*tm31+fc(i1,2)*tm32+fc(i1,3)*tm33
        cp(i1,1)=(fc(i1,1)+id(i1,1))*tm11
        cp(i1,2)=(fc(i1,1)+id(i1,1))*tm21+(fc(i1,2)+id(i1,2))*tm22
        cp(i1,3)=(fc(i1,1)+id(i1,1))*tm31+(fc(i1,2)+id(i1,2))*tm32+ &
        (fc(i1,3)+id(i1,3))*tm33
    end do

    RETURN

    END SUBROUTINE TRAROT2
!**********************************************************************
!**********************************************************************

    subroutine mdsav(ni,qfileh)

!**********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension idum(mbond+3),bodum(mbond+3),qat2(2)
    dimension atomtemp(nat)
    character(400) :: q11_1,q11_2,q11_3
    character(25) :: qfileh
    character(33) :: qfile2
    character(4) :: qext
    character(6) :: qmdfi
    character (7) :: var
    character (3) :: qat2,pepname
    character (1) :: qrtemp
!***********************************************************************
!                                                                      *
!     Save coordinates, velocities and accelerations of MD-system      *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In mdsav'
        call timer(65)
        close (65)
    end if
    qmdfi='moldyn'
    pepname='   '
    ipeptide=0
    if (ni == 2) qmdfi='molsav'

    if (iopt == 0) then

    !***********************************************************************
    !                                                                      *
    !     Write out connection table to units 7 and 8                      *
    !                                                                      *
    !***********************************************************************
        do i1=1,mbond+3
            idum(i1)=nzero
            bodum(i1)=zero
        end do
        if (napp == 1) &
        open (7,file='fort.7',status='unknown',position='append')
        if (napp /= 1) &
        open (7,file='fort.7',status='unknown')
        nsbmaxh=5*((nsbmax/5)+1)
        write (7,100)na,qmol,mdstep,nsbmaxh
        if (nbiolab == 2) write (67,101)na,qmol
        do i1=1,na
            bosum=0.0d0
            do i3=1,nsbmax
                if (iag(i1,2+i3) > 0) bosum=bosum+bo(nubon1(i1,i3))
            end do
            if (icolab == 0) then
                if (nsbmax < 5) then
                    write (7,200)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,5-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,5-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                    if (nbiolab == 2) then      !Delphi-connection table output
                        write (67,201)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2))
                    end if
                else if (nsbmax < 10) then
                    write (7,210)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,10-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,10-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax < 15) then
                    write (7,220)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,15-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,15-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax < 20) then
                    write (7,230)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,20-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,20-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax < 25) then
                    write (7,240)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,25-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,25-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax > 25) then
                    write (7,250)i1,iag(i1,1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,35-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,35-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                end if
            else
                if (nsbmax < 5) then
                    write (7,202)i1,qa(i1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,5-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,5-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                    if (nbiolab == 2) then      !Delphi-connection table output
                        write (67,201)i1,qa(i1),(iag(i1,2+i2),i2=1,iag(i1,2))
                    end if
                else if (nsbmax < 10) then
                    write (7,212)i1,qa(i1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,10-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,10-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax < 15) then
                    write (7,222)i1,qa(i1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,15-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,15-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax < 20) then
                    write (7,232)i1,qa(i1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,20-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,20-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax < 25) then
                    write (7,242)i1,qa(i1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,25-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,25-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbmax > 25) then
                    write (7,252)i1,qa(i1),(iag(i1,2+i2),i2=1,iag(i1,2)), &
                    (idum(i2),i2=1,35-iag(i1,2)), &
                    iag(i1,3+mbond),(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
                    (bodum(i2),i2=1,35-iag(i1,2)),abo(i1),vlp(i1),ch(i1)
                end if
            end if
        end do
        boss=zero
        vlps=0.0d0
        if (napp == 1) &
        open (8,file='fort.8',status='unknown',position='append')
        if (napp /= 1) &
        open (8,file='fort.8',status='unknown')
        nsbmaxh=5*((nsbma2/5)+1)
        write (8,100)na,qmol,mdstep,nsbmaxh
        chsum=0.0d0
        do i1=1,na
            bosum=0.0d0
            do i3=1,nsbma2
                if (ia(i1,2+i3) > 0) bosum=bosum+bo(nubon2(i1,i3))
            end do
            if (icolab == 0) then
                if (nsbma2 < 5) then
                    write (8,200)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,5-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,5-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 10) then
                    write (8,210)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,10-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,10-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 15) then
                    write (8,220)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,15-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,15-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 20) then
                    write (8,230)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,20-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,20-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 25) then
                    write (8,240)i1,ia(i1,1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,25-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,25-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                end if
            else
                if (nsbma2 < 5) then
                    write (8,202)i1,qa(i1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,5-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,5-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 10) then
                    write (8,212)i1,qa(i1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,10-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,10-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 15) then
                    write (8,222)i1,qa(i1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,15-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,15-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 20) then
                    write (8,232)i1,qa(i1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,20-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,20-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                else if (nsbma2 < 25) then
                    write (8,242)i1,qa(i1),(ia(i1,2+i2),i2=1,ia(i1,2)), &
                    (idum(i2),i2=1,25-ia(i1,2)), &
                    ia(i1,3+mbond),(bo(nubon2(i1,i2)),i2=1,ia(i1,2)), &
                    (bodum(i2),i2=1,25-ia(i1,2)),abo(i1),vlp(i1),ch(i1)
                end if
            end if

            boss=boss+bosum/2.0d0
            vlps=vlps+vlp(i1)
            chsum=chsum+ch(i1)
        end do
        write (7,*)2.0d0*boss,vlps,2.0d0*boss+2.0d0*vlps,chsum
        close(8)
        close(7)
    !***********************************************************************
    !                                                                      *
    !     Write out GRASP-type connection table to unit 11                 *
    !                                                                      *
    !***********************************************************************
    !    df398> supress output to fort.11 during flocky training
    !    if (napp == 1) &
    !    open (11,file='fort.11',status='unknown',position='append')
    !    if (napp /= 1) &
    !    open (11,file='fort.11',status='unknown')

    !    do i1=1,na
    !        q11_1=''
    !        q11_2=''
    !        q11_3=''
    !        write (q11_1,'(200i8)')i1,iag(i1,1),iag(i1,2), &
    !        (iag(i1,2+i2),i2=1,iag(i1,2)),iag(i1,3+mbond)
    !        iq11_1=(iag(i1,2)+4)*8
    !        write (q11_2,'(200f7.3)')(bo(nubon1(i1,i2)),i2=1,iag(i1,2)), &
    !        abo(i1),vlp(i1),ch(i1)
    !        iq11_2=(iag(i1,2)+3)*8
    !        q11_3(1:iq11_1)=q11_1(1:iq11_1)
    !        q11_3(iq11_1+1:iq11_1+iq11_2)=q11_2(1:iq11_2)
    !        ilhulp=iq11_1+iq11_2
    !        if (ilhulp < 50) write (11,'(a50)')q11_3(1:50)
    !        if (ilhulp >= 50 .AND. ilhulp < 100) &
    !        write (11,'(a100)')q11_3(1:100)
    !        if (ilhulp >= 100 .AND. ilhulp < 150) &
    !        write (11,'(a150)')q11_3(1:150)
    !        if (ilhulp >= 150 .AND. ilhulp < 200) &
    !        write (11,'(a200)')q11_3(1:200)
    !        if (ilhulp >= 200 .AND. ilhulp < 250) &
    !        write (11,'(a250)')q11_3(1:250)
    !        if (ilhulp >= 250 .AND. ilhulp < 300) &
    !        write (11,'(a300)')q11_3(1:300)
    !        if (ilhulp >= 300 .AND. ilhulp < 350) &
    !        write (11,'(a350)')q11_3(1:350)
    !        if (ilhulp >= 350 .AND. ilhulp < 400) &
    !        write (11,'(a400)')q11_3(1:400)
    !    end do
    !    close(11)
             
    end if

    if (noutpt == 0) then
        write (var,'(f7.4)')float(mdstep/nsav)/1d4
        if (ni == 0) open (unit=67,file=qmdfi//var(3:7), &
        status='unknown')
        write (67,300)qmol
        do i1=1,na
            write (67,400)i1,qa(i1),(c(i1,i2),i2=1,3)
        end do
        write (67,*)
        close(67)
    end if

    if (noutpt == 2) then
        open (88,file='moldyn.bgf',status='unknown',position='append')
        call writebgf(88)
        close (88)
    end if

    if ((ni == 1 .AND. iopt == 0) .OR. &
    (ni == 1 .AND. iopt == 1 .AND. iflga == 1)) then
        qrtemp=qr
        if (qr == 'I') qr='C'
        if (qfileh == ' ') then
            write (*,*)'Warning: no file name given; use Unknown'
            qfileh='Unknown'
        end if
        qfile2=qfileh
        if (imodfile == 0) then
            istart=1
            qstrana1(1:25)=qfileh
            call stranal(istart,iend,vout,iout,1)
            qfile2=qfileh(istart:iend-1)//".geo"
        end if
        !df398> supress output during flocky training
        !call writegeo(98)

        !if (imodfile == 1 .OR. iopt == 0 .AND. nbiolab == 0) then
        !    open (88,file=qfile2,status='unknown')
        !    call writegeo(88)
        !    close (88)
        !end if

        qr=qrtemp

        if (iopt == 0) then

            do i1=1,na
                write (56,410) i1,ch(i1)
                write (55,410) i1,chgbgf(i1)
            end do
        !*********************************************************************
        !                                                                    *
        !     Write .pdb output file                                         *
        !                                                                    *
        !*********************************************************************
        ! df398> supress output during flocky training
        !    open (unit=47,file='output.pdb',status='unknown')
        !    do i1=1,na
        !        write (47,412)'ATOM  ',i1,qa(i1),pepname,ipeptide,c(i1,1), &
        !        c(i1,2),c(i1,3),1.0,2.2,qa(i1)
        !    end do
        !    write (47,*) 'TER'
        !    write (47,*) 'END'
        !    close (47)

            if (nsurp == 0) then
                if (kx > 0 .OR. ky > 0 .OR. kz > 0) then
                    qrtemp=qr
                !*********************************************************************
                !                                                                    *
                !     Write crystal structure including periodic images              *
                !                                                                    *
                !*********************************************************************
                !     mux=(1+kx+kx)
                !     muy=(1+ky+ky)
                !     muz=(1+kz+kz)
                !     qr='F'
                !     write (86,'(2x,a1,1x,a60)')qr,qmol
                !     qr=qrtemp
                !     write (86,'(3f10.4)')mux*axiss(1),muy*axiss(2),muz*axiss(3)
                !     write (86,'(3f10.4)')angle(1),angle(2),angle(3)
                !     do i1=1,na
                !     write (86,'(i4,1x,a2,3x,3d22.15)')i1,qa(i1),(c(i1,i2),i2=1,3)
                !     end do
                !     nhulp=na+1
                !     do k1=-kx,kx
                !     do k2=-ky,ky
                !     do k3=-kz,kz
                !     if (k1.ne.0.or.k2.ne.0.or.k3.ne.0) then
                !     do i1=1,na
                !     cx=c(i1,1)+k1*tm11
                !     cy=c(i1,2)+k1*tm21+k2*tm22
                !     cz=c(i1,3)+k1*tm31+k2*tm32+k3*tm33
                !     write (86,'(i4,1x,a2,3x,3d22.15)')nhulp,qa(i1),cx,cy,cz
                !     nhulp=nhulp+1
                !     end do
                !     end if
                !     end do
                !     end do
                !     end do
                !     write (86,*)
                !*********************************************************************
                !                                                                    *
                !     Write crystal structure with extra unit cells                  *
                !                                                                    *
                !*********************************************************************
                    mux=1+iexx
                    muy=1+iexy
                    muz=1+iexz
                    qr='F'
                    write (85,'(2x,a1,1x,a60)')qr,qmol
                    qr=qrtemp
                    write (85,'(3f10.4)')mux*axiss(1),muy*axiss(2),muz*axiss(3)
                    write (85,'(3f10.4)')angle(1),angle(2),angle(3)
                    do i1=1,na
                        write (85,'(i4,1x,a2,3x,3d22.15)')i1,qa(i1),(c(i1,i2),i2=1,3)
                    end do
                    nhulp=na+1
                    do k1=0,iexx
                        do k2=0,iexy
                            do k3=0,iexz
                                if (k1 /= 0 .OR. k2 /= 0 .OR. k3 /= 0) then
                                    do i1=1,na
                                        cx=c(i1,1)+k1*tm11
                                        cy=c(i1,2)+k1*tm21+k2*tm22
                                        cz=c(i1,3)+k1*tm31+k2*tm32+k3*tm33
                                        write (85,'(i4,1x,a2,3x,3d22.15)')nhulp,qa(i1),cx,cy,cz
                                        nhulp=nhulp+1
                                    end do
                                end if
                            end do
                        end do
                    end do
                    write (85,*)

                end if
            end if
        end if

    end if

    if (ni == 0 .OR. ni == 2) then
    !*********************************************************************
    !                                                                    *
    !     Write ASCII trajectory file                                    *
    !                                                                    *
    !*********************************************************************
        close(66)
        if (ni == 0) open(unit=66,file=qmdfi//'.vel',status='unknown')
        ! if (ni == 0) open(unit=66,file='vels',status='unknown')
        if (ni == 2) then
            write (var,'(f7.4)')float(mdstep/nsav3)/1d4
            open (unit=66,file=qmdfi//var(3:7),status='unknown')
        end if
        write (66,500)axis(1),axis(2),axis(3)
        write (66,550)angle(1),angle(2),angle(3)
        write (66,600)na,((c(i,j),j=1,3),qlabel(i),i=1,na)
        write (66,700)((vel(j,i),j=1,3),i=1,na)
        write (66,800)((accel(j,i),j=1,3),i=1,na)
        write (66,900)((aold(j,i),j=1,3),i=1,na)
        write (66,1000)tempmd
        write (66,1050)
        close (66)
    end if
    if (ni /= 2 .AND. iopt == 0) then

        open (unit=68,file='xmolout',status='unknown',position='append')
        write (68,1200)na
        if (ixmolo /= 1) write (68,1300)qmol,mdstep+nit+nprevrun,estrc, &
        axis(1),axis(2),axis(3),angle(1),angle(2),angle(3)
        if (ixmolo == 1) write (68,1390)tm11,tm22,tm33,tm32,tm31, &
        tm21,tm11*tm22*tm33, &
        tstep*1d+15*(mdstep+nit)/1000.0d0,mdstep+nit,tstep*1d+15,tset, &
        tstep,tstep,tstep, &
        (estrc+ekin)/23.06d0,pressu/1000.0d0,tempmd

        if (ixmolo == 5) then   ! determine atom temperatures
            do i1=1,na
                ekinat=0
                do i2=1,3
                    ekinat=ekinat+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
                end do
                ekinat=0.50d0*ekinat/convmd
                atomtemp(i1)=2.0d0*ekinat/(3.0d0*rgasc*xjouca/1.0d3)
            end do
        end if

        do i1=1,na
            if (ixmolo == 0) write (68,1400)qa(i1),(c(i1,i2),i2=1,3)
            if (ixmolo == 1) write (68,1400)qa(i1),(c(i1,i2),i2=1,3), &
            (vel(i2,i1)/1e+10,i2=1,3),iag(i1,3+mbond)
            if (ixmolo == 2) write (68,1401)qa(i1),(c(i1,i2),i2=1,3), &
            iag(i1,3+mbond)
            if (ixmolo == 3) write (68,1405)qa(i1),(c(i1,i2),i2=1,3), &
            iag(i1,3+mbond),estrain(i1)+vincr(ia(i1,1))
            if (ixmolo == 4) write (68,1405)qa(i1),(c(i1,i2),i2=1,3), &
            iag(i1,3+mbond),abo(i1)
            if (ixmolo == 5) write (68,1406)qa(i1),(c(i1,i2),i2=1,3), &
            atomtemp(i1),iag(i1,3+mbond)
            if (ixmolo == 6) write (68,1400)qa(i1),(c(i1,i2),i2=1,3), &
            (d(i2,i1),i2=1,3),iag(i1,3+mbond)
        end do
        close (68)

        if (itrout /= 0) then
            iless=0

            if (ignotr > 0) then         !Find out number of ignotr-atoms
                do i1=1,na
                    if (ia(i1,1) == ignotr) iless=iless+1
                end do
            end if

            open (unit=69,file='diff_traj.xyz',status='unknown', &
            position='append')
            write (69,1200)na-iless
            write (69,1300)qmol,mdstep+nit+nprevrun,estrc, &
            axis(1),axis(2),axis(3),angle(1),angle(2),angle(3)
            do i1=1,na
                if (ia(i1,1) /= ignotr) then  ! ignore element
                    if (ixmolo == 0) write (69,1400)qa(i1),(cp(i1,i2),i2=1,3)
                    if (ixmolo == 1) write (69,1400)qa(i1),(cp(i1,i2),i2=1,3), &
                    (vel(i2,i1)/1d+10,i2=1,3),iag(i1,3+mbond)
                    if (ixmolo == 2) write (69,1401)qa(i1),(cp(i1,i2),i2=1,3), &
                    iag(i1,3+mbond)
                    if (ixmolo == 3) write (69,1405)qa(i1),(cp(i1,i2),i2=1,3), &
                    iag(i1,3+mbond),estrain(i1)+vincr(ia(i1,1))
                end if
            end do
            close (69)
            open (unit=69,file='translate.out',status='unknown')
            write (69,'(i6)')na
            do i1=1,na
                write (69,'(4i6)')i1,id(i1,1),id(i1,2),id(i1,3)
            end do
            close (69)
        end if

        idummy=0
        call molanal(idummy,cutof3,1,ntotmol)
        if (ignore > 0) call molanal(ignore,cutof3,2,ntotmol)
        if (cutmol1 > 0.01d0) call molanal(idummy,cutmol1,3,ntotmol)
        if (cutmol2 > 0.01d0) call molanal(idummy,cutmol2,4,ntotmol)
        if (cutmol3 > 0.01d0) call molanal(idummy,cutmol3,5,ntotmol)
        if (cutmol4 > 0.01d0) call molanal(idummy,cutmol4,6,ntotmol)
        if (cutmol5 > 0.01d0) call molanal(idummy,cutmol5,7,ntotmol)
    end if
!*********************************************************************
!                                                                    *
!     Write summary.txt output-file                                  *
!                                                                    *
!*********************************************************************
    if (ni /= 2 .AND. iopt == 0 .AND. mdstep > 0) then
        open (unit=69,file='summary.txt',status='unknown', &
        position='append')
        write (69,'(i9,i6,12f12.2)')mdstep+nit+nprevrun,ntotmol,tottime, &
        estrc,tm11*tm22*tm33,tempmd,pressu, &
        (xmasmd*1.0d+24)/(avognr*tm11*tm22*tm33)
        close (69)
    end if
!*********************************************************************
!                                                                    *
!     Generate BIOGRAF output-file                                   *
!                                                                    *
!*********************************************************************
    if ((ni == 1 .AND. iopt == 0) &
     .OR. (ni == 1 .AND. iopt == 1 .AND. iflga == 1)) then

        if (qfileh == ' ') then
            write (*,*)'Warning: no file name given; use Unknown'
            qfileh='Unknown'
        end if
        qfile2=qfileh
        if (imodfile == 0) then
            istart=1
            qstrana1(1:25)=qfileh
            call stranal(istart,iend,vout,iout,1)
            qfile2=qfileh(istart:iend-1)//".bgf"
        end if
        call writebgf(90)

        if (imodfile == 1 .OR. iopt == 0 .AND. nbiolab == 0) then
            open (88,file=qfile2,status='unknown')
            call writebgf(88)
            close (88)
        end if

    end if
         
    return
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format (i4,1x,a40,'Iteration:',i8,' #Bonds:',i4)
    101 format (i3,2x,a40)
    200 format (8i5,8f7.3)
    201 format (8i3)
    202 format (i5,2x,a2,1x,6i5,8f7.3)
    210 format (13i5,13f7.3)
    212 format (i5,2x,a2,1x,11i5,13f7.3)
    220 format (18i5,18f7.3)
    222 format (i5,2x,a2,1x,16i5,18f7.3)
    230 format (23i5,23f7.3)
    232 format (i5,2x,a2,1x,21i5,23f7.3)
    240 format (28i5,28f7.3)
    242 format (i5,2x,a2,1x,26i5,28f7.3)
    250 format (38i5,38f7.3)
    252 format (i5,2x,a2,1x,36i5,38f7.3)
    300 format (2x,a1,1x,a60)
    301 format (2x,a1,1x,f6.2,a60)
    302 format (2x,a1,1x,2f6.2,a60)
    310 format (2x,a1,1x,a60)
    320 format (3f10.4)
    400 format (i4,1x,a2,3x,3(d21.14,1x),1x,a5,1x,i5)
    410 format (i4,f12.6)
    412 format(A6,I5,1x,A2,3x,A3,2x,i4,4x,3f8.3,f6.2,f6.2,4x,2x,A6)
    500 format (1x,'Lattice parameters:',/(3f15.8))
    550 format (3f15.8)
    600 format (i4,1x,'Atom coordinates (Angstrom):',/ &
    (3d24.15,1x,a5))
    700 format (1x,'Atom velocities (Angstrom/s):',/(3d24.15))
    800 format (1x,'Atom accelerations (Angstrom/s**2):',/(3d24.15))
    900 format (1x,'Previous atom accelerations:',/(3d24.15))
    1000 format (1x,'MD-temperature (K):',/(1d24.15))
    1050 format (1x,'Connections, bond orders and lone pairs:')
    1100 format (8i3,8f8.4)
    1200 format (i4)
    1300 format (a40,i6,f12.2,6f7.2)
    1390 format (6(1x,f9.4),f11.4,f9.4,i8,f7.4,f10.2,3(f7.2),3f12.4)
    1400 format (a2,3f10.5,3f15.5,i6)
    1401 format (a2,3f10.5,i6)
    1405 format (a2,3f10.5,i6,2x,f15.5)
    1406 format (a2,3f10.5,f15.5,i6)
    1500 format ('BIOGRF',i4)
    1600 format ('XTLGRF',i4)
    1700 format ('DESCRP ',a60)
    1800 format ('REMARK ',a60)
    1900 format ('FFIELD ',a40)
    2000 format ('RUTYPE ',a40)
    2100 format ('CRYSTX ',6f11.5)
    2200 format ('CELLS ',6i5)
    2300 format ('#              At1 At2   R12    Force1  Force2  ', &
    'dR12/dIteration(MD only)')
    2400 format ('BOND RESTRAINT ',2i4,f8.4,f8.2,f8.5,f10.7)
    2500 format ('#               At1 At2 At3 Angle   Force1  Force2', &
    '  dAngle/dIteration (MD only)')
    2600 format ('ANGLE RESTRAINT ',3i4,2f8.2,f8.4,f9.6)
    2700 format ('#                 At1 At2 At3 At3 Angle   Force1  ', &
    'Force2  dAngle/dIteration (MD only)')
    2800 format ('TORSION RESTRAINT ',4i4,2f8.2,f8.4,f9.6)
    2900 format ('FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,', &
    '3f10.5,1x,a5,i3,i2,1x,f8.5)')
    3000 format ('HETATM',1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x, &
    a5,i3,i2,1x,f8.5)
    3100 format ('FORMAT CONECT (a6,12i6)')
    3200 format ('CONECT',12i6)
    3300 format ('UNIT ENERGY   kcal')
    3400 format ('ENERGY',5x,f14.6)
    3500 format ('END')
    end subroutine mdsav

!***********************************************************************
!***********************************************************************

    subroutine inivel

!***********************************************************************
    include 'cbka.blk'
!***********************************************************************
!                                                                      *
!     Read coordinates, velocities and accelerations of MD-system      *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In inivel'
        call timer(65)
        close (65)
    end if
    ivels=0
    if (irun > 1) goto 10
    if (nmm > 0) goto 10
    open(unit=66,file='moldyn.vel',status='old',err=10)
    ivels=1
    read (66,*)
    read (66,100)aaxis,baxis,caxis
    read (66,100)angles(1),angles(2),angles(3)
    if (qr == 'F' .OR. qr == 'P' .OR. ngeofor == 1) then
        axis(1)=aaxis
        axis(2)=baxis
        axis(3)=caxis
        axiss(1)=axis(1)
        axiss(2)=axis(2)
        axiss(3)=axis(3)
        angle(1)=angles(1)
        angle(2)=angles(2)
        angle(3)=angles(3)
        halfa=angle(1)*dgrrdn
        hbeta=angle(2)*dgrrdn
        hgamma=angle(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axis(1)*sinbet*sinphi
        tm21=axis(1)*sinbet*cosphi
        tm31=axis(1)*cosbet
        tm22=axis(2)*sinalf
        tm32=axis(2)*cosalf
        tm33=axis(3)
    end if
    if (aaxis /= axis(1) .OR. baxis /= axis(2) .OR. caxis /= axis(3)) &
    stop 'Wrong cell parameters in moldyn.vel'
    read (66,200)nan
    if (nan /= na) stop 'Wrong number of atoms in moldyn.vel-file'
    if (nbiolab == 2) write (*,*)'Warning: using labels in vels-file'
    read (66,250)((c(i,j),j=1,3),qlabel(i),i=1,na)
    read (66,*)
    read (66,300)((vel(j,i),j=1,3),i=1,na)
    read (66,*)
    read (66,300)((accel(j,i),j=1,3),i=1,na)
    read (66,*)
    read (66,300,end=10,err=10)((aold(j,i),j=1,3),i=1,na)
    read (66,*)
    read (66,300,end=10,err=10)tempmd
    read (66,*)
    read (66,350,end=10,err=10)nsbma2
    10 continue
    if (ivels == 0 .AND. nrand == 1 .AND. na > 1) then
    !***********************************************************************
    !                                                                      *
    !     Random initial velocities                                        *
    !                                                                      *
    !***********************************************************************
        sumx=0.0d0
        sumy=0.0d0
        sumz=0.0d0
        do i1=1,na
            xx=2.0d0*random(dseed)-1.0d0 !df398 replaced rand(0) to random(dseed). rand is not supported by nagfor.
            yy=2.0d0*random(dseed)-1.0d0
            zz=2.0d0*random(dseed)-1.0d0
            xyz=1/sqrt(xx*xx+yy*yy+zz*zz)
            vel(1,i1)=xx/xyz
            vel(2,i1)=yy/xyz
            vel(3,i1)=zz/xyz
            sumx=sumx+vel(1,i1)
            sumy=sumy+vel(2,i1)
            sumz=sumz+vel(3,i1)
        end do
        velsq=0.0d0
        do i1=1,na
            vel(1,i1)=vel(1,i1)-sumx/na
            vel(2,i1)=vel(2,i1)-sumy/na
            vel(3,i1)=vel(3,i1)-sumz/na
            velsq=velsq+xmasat(i1)*(vel(1,i1)*vel(1,i1)+ &
            vel(2,i1)*vel(2,i1)+vel(3,i1)*vel(3,i1))
        end do
    !     ekin=half*velsq/convmd
    !     tempn=2.0*ekin/(float(3*na)*rgasc*xjouca/1.0d3)
    !     factor=sqrt(tset/tempn)
        velsq=0.0d0
        do i1=1,na
            ekinat=0.0d0
            do i2=1,3
                ekinat=ekinat+xmasat(i1)*vel(i2,i1)*vel(i2,i1)
            end do
            ekinat=0.50d0*ekinat/convmd
            tempat=2.0d0*ekinat/(3.0d0*rgasc*xjouca/1.0d3)
            factor=sqrt(tset/tempat)
            vel(1,i1)=vel(1,i1)*factor
            vel(2,i1)=vel(2,i1)*factor
            vel(3,i1)=vel(3,i1)*factor
            velsq=velsq+xmasat(i1)*(vel(1,i1)*vel(1,i1)+ &
            vel(2,i1)*vel(2,i1)+vel(3,i1)*vel(3,i1))
        end do
        ekin=half*velsq/convmd
        tempmd=2.0d0*ekin/(float(3*namov)*rgasc*xjouca/1.0d3)
    end if
!*********************************************************************
!                                                                    *
!     Format part                                                    *
!                                                                    *
!*********************************************************************
    100 format(3d15.8)
    200 format(i4)
    250 format(3d24.15,1x,a5)
    300 format(3d24.15)
    350 format(i3)
    400 format (8i3,8f8.4)
    return
    end subroutine inivel
!***********************************************************************
!***********************************************************************

    subroutine readc

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(2) :: qig
    character(6) :: qident
    character(20) :: qhulp
!     dimension qident(100)
!***********************************************************************
!                                                                      *
!     Read control file                                                *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readc'
        call timer(65)
        close (65)
    end if
    if (mdstep > 0 .OR. nit > 0) nmmsav=nmm
!***********************************************************************
!                                                                      *
!     Set default values                                               *
!                                                                      *
!***********************************************************************
    nreac=0
    axis1=200.0d0
    axis2=200.0d0
    axis3=200.0d0
    angle1=90.0d0
    angle2=90.0d0
    angle3=90.0d0
    cutof2=0.001d0
    cutof3=0.300d0
    tsetor=298.0d0
    tset2=298.0d0
    pset=0.0d0
    tincr=0.0d0
    tstep0=0.25d0
!     swa=0.0   !Moved to force field
!     swb=12.5  !Moved to force field
    taut=2.5d0
    taut2=2.5d0
    ndtau=50000
    taup=500.0d0
    vqnd=100.0d0
    errnh=1.0d0
    range=2.5d0
    maxstp=1000
    nequi=0
    nmethod=3
    ncha=3
    ncha2=1
    nchaud=1
    nvlist=25
    nrep1=5
    nsav=50
    icheck=0
    ivels=0
    itfix=0
    ncontrol=25
    noutpt=0
    napp=0
    nsurp=0
    ncons=2
    nrand=0
    nmm=0
    endpo=1.0d0
    endpo2=1.0d0
    nfc=50
    nsav2=50
    nmmax=50
    i5758=0
    parc1=1.0d0
    parc2=0.001d0
    icell=0
    ingeo=1
    ccpar=1.0005d0
    icelo2=0
    nrdd=0
    nrddf=200000
    nbiolab=0
    ngeofor=0
    nincrop=0
    accerr=2.50d0
    vrange=2.50d0
    vlbora=5.00d0
    nsav3=1000
    nhop2=25
    nprevrun=0
    ndebug=0
    volcha=10.00d0
    ixmolo=0
    inpt=0
    iconne=0
    imolde=0
    ianaly=0
    icentr=0
    itrans=1
    itrout=0
    tpnrad=300.0d0
    ityrad=3
    iexx=1
    iexy=1
    iexz=1
    syscha=0.00d0
    inmov1=0
    inmov2=0
    vfield=0.00d0
    itstep=0
    ifreq=0
    isymm=1
    icpres=0
    ipresm=0
    ilavel=0
    exfx=0.0d0
    exfy=0.0d0
    exfz=0.0d0
    ignore=0
    ignotr=0
    ipolar=1
    ikeep=0
    cutmol1=0.0d0
    cutmol2=0.0d0
    cutmol3=0.0d0
    cutmol4=0.0d0
    cutmol5=0.0d0
    delvib=1.0d-4
    convg=1.0d-6  !Convergence criterion for EEM conjugate gradient solve
    edeep=0.0d0   !Piston magnitude
    rdeep=0.0d0   !0: infinite piston; 1: Edeep piston
    pshft=0.0d0   !Piston position
    rcut=1.0    !Piston width
    ipdir=1    !Piston direction (1:x; 2:y; 3:z)
    ireflx=0    !non-reflective x-boundary
    irefly=0    !non-reflective y-boundary
    ireflz=0    !non-reflective z-boundary
    imfreq=0    !-1: do not match Reax/Jaguar vibrational modes
    iflext=0    !0: fixed time step  1:velocity-dependent time step
    dtmax=1.0d-2  !Maximum atom displacement per step (with iflext=1)
    icobo=0     !0: use uncorrected bond orders in xmolout, fort.7 and fort.71  1: use corrected bond orders
    ibocha=0     !0: do not output bondchange.out  1: generate bondchange.out
    iremov=0     !0: do not remove molecules >0 frequency of molecule removal
    vmamin= 1.0d-1  !Removal mass criterion
    vmami2= 1.0d-2  !Removal mass criterion  2
    icolab=0     !0: use atom type numbers in connection table  1: use atom names in connection table
!     shock variables
    shock_vel = 2.0d0 ! impact velocity for shock simulations (nm/ps)
    shock_z_sep = 10.0d0 ! separation z value to apply initial velocities in shocks
    ishock_type = 0.0d0 ! shock type. 0: simple impact; 1: compressing c axis
!     Hugoniostat variables
    Hug_E0 = 0.d0 ! Reference energy
    Hug_P0 = 0.d0 ! Reference pressure
    Hug_V0 = 0.d0 ! Reference volume
!     Shear flow simulations for viscosity
    xImpVcm = 1.0d0 ! velocity applied in shear simulations (in nm/ps), left half mover at -xImpVcm and right at +xImpVcm
!***********************************************************************
!                                                                      *
!     Read control-file                                                *
!                                                                      *
!***********************************************************************
    open (unit=10, file='control', status='old')
    10 read (10,'(a20)',end=20,err=30)qhulp
    if (qhulp(1:1) == '#') goto 10
    read (qhulp,*,err=30)vhulp
    read (qhulp,'(8x,a6)',err=30)qident
    if (qident == 'Hug_V0') Hug_P0=vhulp
    if (qident == 'Hug_P0') Hug_V0=vhulp
    if (qident == 'Hug_E0') Hug_E0=vhulp
    if (qident == 'shea_v') xImpVcm=vhulp
    if (qident == 'shok_t') ishock_type=int(vhulp)
    if (qident == 'shok_z') shock_z_sep=vhulp
    if (qident == 'shok_v') shock_vel=vhulp
    if (qident == 'nreac') nreac=int(vhulp)
    if (qident == 'axis1') axis1=vhulp
    if (qident == 'axis2') axis2=vhulp
    if (qident == 'axis3') axis3=vhulp
    if (qident == 'angle1') angle1=vhulp
    if (qident == 'angle2') angle2=vhulp
    if (qident == 'angle3') angle3=vhulp
    if (qident == 'cutof2') cutof2=vhulp
    if (qident == 'cutof3') cutof3=vhulp
    if (qident == 'mdtemp') tsetor=vhulp
    if (qident == 'mdtem2') tset2=vhulp
    if (qident == 'mdpres') pset=vhulp*0.001d0
    if (qident == 'tincr') tincr=vhulp
    if (qident == 'tstep') tstep0=vhulp
!     if (qident.eq.'lowtap') swa=vhulp !Moved to force field
!     if (qident.eq.'uptap') swb=vhulp  !Moved to force field
    if (qident == 'tdamp1') taut=vhulp
    if (qident == 'tdamp2') taut2=vhulp
    if (qident == 'ntdamp') ndtau=int(vhulp)
    if (qident == 'pdamp1') taup=vhulp
    if (qident == 'tdhoov') vqnd=vhulp
    if (qident == 'achoov') errnh=vhulp/100.0d0
    if (qident == 'range') range=vhulp
    if (qident == 'nmdit') maxstp=int(vhulp)
    if (qident == 'nmdeqi') nequi=int(vhulp)
    if (qident == 'imdmet') nmethod=int(vhulp)
    if (qident == 'icharg') ncha=int(vhulp)
    nchaold=ncha
    if (qident == 'ichaen') ncha2=int(vhulp)
    if (qident == 'ichupd') nchaud=int(vhulp)
    if (qident == 'iout1') nrep1=int(vhulp)
    if (qident == 'iout2') nsav=int(vhulp)
    if (qident == 'icheck') ntest=int(vhulp)
    if (qident == 'ivels') nvel=int(vhulp)
    if (qident == 'itfix') ntscale=int(vhulp)
    if (qident == 'irecon') ncontrol=int(vhulp)
    if (qident == 'iout3') noutpt=int(vhulp)
    if (qident == 'iappen') napp=int(vhulp)
    if (qident == 'isurpr') nsurp=int(vhulp)
    if (qident == 'itdmet') ncons=int(vhulp)
    if (qident == 'iravel') nrand=int(vhulp)
    if (qident == 'imetho') nmm=int(vhulp)
    if (qident == 'endmm') endpo=vhulp
    endpoold=endpo
    if (qident == 'endmd') endpo2=vhulp
    if (qident == 'imaxmo') nfc=int(vhulp)
    nfcold=nfc
    if (qident == 'iout4') nsav2=int(vhulp)
    if (qident == 'imaxit') nmmax=int(vhulp)
    nmmaxold=nmmax
    if (qident == 'iout5') i5758=int(vhulp)
    if (qident == 'parsca') parc1=vhulp
    if (qident == 'parext') parc2=vhulp
    if (qident == 'icelop') icell=int(vhulp)
    icellold=icell
    if (qident == 'igeopt') ingeo=int(vhulp)
    if (qident == 'celopt') ccpar=vhulp
    if (qident == 'icelo2') icelo2=int(vhulp)
    icelo2old=icelo2
    if (qident == 'ideve1') nrdd=int(vhulp)
    if (qident == 'ideve2') nrddf=int(vhulp)
    if (qident == 'ibiola') nbiolab=int(vhulp)
    if (qident == 'igeofo') ngeofor=int(vhulp)
    if (qident == 'iincop') nincrop=int(vhulp)
    if (qident == 'accerr') accincr=vhulp
    if (qident == 'iout6') nsav3=int(vhulp)
    if (qident == 'irten') nhop2=int(vhulp)
    if (qident == 'npreit') nprevrun=int(vhulp)
    if (qident == 'idebug') ndebug=int(vhulp)
    if (qident == 'volcha') volcha=vhulp
    if (qident == 'ixmolo') ixmolo=int(vhulp)
    if (qident == 'inpt') inpt=int(vhulp)
    if (qident == 'iconne') iconne=int(vhulp)
    if (qident == 'imolde') imolde=int(vhulp)
    if (qident == 'ianaly') ianaly=int(vhulp)
    if (qident == 'icentr') icentr=int(vhulp)
    if (qident == 'itrans') itrans=int(vhulp)
    if (qident == 'itrout') itrout=int(vhulp)
    if (qident == 'nvlist') nvlist=int(vhulp)
    if (qident == 'vrange') vrange=vhulp
    if (qident == 'vlbora') vlbora=vhulp
    if (qident == 'tpnrad') tpnrad=vhulp
    if (qident == 'ityrad') ityrad=int(vhulp)
    if (qident == 'iexx') iexx=int(vhulp)
    if (qident == 'iexy') iexy=int(vhulp)
    if (qident == 'iexz') iexz=int(vhulp)
    if (qident == 'syscha') syscha=vhulp
    if (qident == 'inmov1') inmov1=int(vhulp)
    if (qident == 'inmov2') inmov2=int(vhulp)
    if (qident == 'itstep') itstep=int(vhulp)
    if (qident == 'ifreq') ifreq=int(vhulp)
    if (qident == 'isymm') isymm=int(vhulp)
    if (qident == 'icpres') icpres=int(vhulp)
    if (qident == 'ipresm') ipresm=int(vhulp)
    if (qident == 'ilavel') ilavel=int(vhulp)
    if (qident == 'delvib') delvib=vhulp
    if (qident == 'exfx') exfx=vhulp
    if (qident == 'exfy') exfy=vhulp
    if (qident == 'exfz') exfz=vhulp
    if (qident == 'edeep') edeep=vhulp
    if (qident == 'rdeep') rdeep=vhulp
    if (qident == 'pshft') pshft=vhulp
    if (qident == 'rcut') rcut=vhulp
    if (qident == 'ipdir') ipdir=int(vhulp)
    if (qident == 'convg') convg=vhulp
    if (qident == 'ignore') ignore=int(vhulp)
    if (qident == 'ignotr') ignotr=int(vhulp)
    if (qident == 'cutmo1') cutmol1=vhulp
    if (qident == 'cutmo2') cutmol2=vhulp
    if (qident == 'cutmo3') cutmol3=vhulp
    if (qident == 'cutmo4') cutmol4=vhulp
    if (qident == 'cutmo5') cutmol5=vhulp
    if (qident == 'ipolar') ipolar=int(vhulp)
    if (qident == 'ikeep') ikeep=int(vhulp)
    if (qident == 'ireflx') ireflx=int(vhulp)
    if (qident == 'irefly') irefly=int(vhulp)
    if (qident == 'ireflz') ireflz=int(vhulp)
    if (qident == 'imfreq') imfreq=int(vhulp)
    if (qident == 'iflext') iflext=int(vhulp)
    if (qident == 'dtmax') dtmax=vhulp
    if (qident == 'icobo') icobo=int(vhulp)
    if (qident == 'iremov') iremov=int(vhulp)
    if (qident == 'vmamin') vmamin=vhulp
    if (qident == 'vmami2') vmami2=vhulp
    if (qident == 'ibocha') ibocha=int(vhulp)
    if (qident == 'icolab') icolab=int(vhulp)
    goto 10
    20 continue
    close (10)
    axis(1)=axis1
    axis(2)=axis2
    axis(3)=axis3
    angle(1)=angle1
    angle(2)=angle2
    angle(3)=angle3
    if (axiss(1) > zero) then
        axis(1)=axiss(1)
        axis(2)=axiss(2)
        axis(3)=axiss(3)
        angle(1)=angles(1)
        angle(2)=angles(2)
        angle(3)=angles(3)
    end if
    if (tincr < 0.0001d0 .AND. tincr > -0.0001d0) tset=tsetor
    iequi=1
    if (nequi > 0) iequi=0
    if (iopt == 1 .AND. napp == 1) then
        stop 'No fort.7 and fort.8 append with iopt=1 !'
    end if
    if (mdstep > 0 .OR. nit > 0) nmm=nmmsav
    if (mdstep > 0 .AND. itstep == 1) then
        tstepmax=tstep
        tstep=tstep*(tsetor/tempmd)
        if (tstep > tstepmax) tstep=tstepmax
    end if
    if (iflext == 0) then
        tstep=1.0d-15*tstep0
        taus=taut
        taut=1.0d-15*taut
        taut2=1.0d-15*taut2
        taup=1.0d-15*taup
        ts2=tstep/2.0d0
        ts22=tstep*ts2
    end if
    return
    30 continue
    write (*,*)'Error reading control-file'
    stop 'Error reading control-file'
!***********************************************************************
!                                                                      *
!     Format part                                                      *
!                                                                      *
!***********************************************************************
    1050 format (f7.3)
    1055 format (f7.4)
    1056 format (f9.4)
    1060 format (i8)
    1070 format (f7.5)
    end subroutine readc
!***********************************************************************
!***********************************************************************

    subroutine valtaper(lnorm,R,rmin,rmax,f,dfdr)

    include 'cbka.blk'
!***********************************************************************
!
!  Passed variables
!
    logical,   intent(in)             :: lnorm
    real*8,    intent(in)             :: R
    real*8,    intent(in)             :: rmin
    real*8,    intent(in)             :: rmax
    real*8,    intent(out)            :: f
    real*8,    intent(out)            :: dfdr



! local variables
    real*8                   :: SWC7MY,SWC6MY,SWC5MY,SWC4MY, SWC3MY,SWC2MY
    real*8                   :: SWC1MY, SWC0MY,D7MY,D1MY
    real*8                   :: SWA2MY,SWA3MY,SWB2MY,SWB3MY
    real*8                   :: rrmax, rrmin, RR, R2, R3

!***********************************************************************
!                                                                      *
!     7th order taper function setup                                   *
!                                                                      *
!***********************************************************************
  if (ndebug == 1) then
      open (65,file='fort.65',status='unknown',position='append')
      write (65,*) 'In valtaper'
      call timer(65)
      close (65)
  end if


  if (lnorm) then
    if (R .le. rmin) then
!*************
!  r < rmin  *
!*************
      f = 1.0D0
      dfdr = 0.0D0
    elseif (R .ge. rmax) then
!*************
!  r > rmax  *
!*************
      f = 0.0D0
      dfdr = 0.0D0
    else
!********************
!  rmin < r < rmax  *
!********************
      D1MY=rmax-rmin
      D7MY=D1MY**7
      SWA2MY=rmin*rmin
      SWA3MY=SWA2MY*rmin
      SWB2MY=rmax*rmax
      SWB3MY=SWB2MY*rmax

      SWC7MY=  20.0D0!/D7MY
      SWC6MY= -70.0D0*(rmin+rmax)!/D7MY
      SWC5MY=  84.0D0*(SWA2MY+3.0D0*rmax*rmin+SWB2MY)!/D7MY
      SWC4MY= -35.0D0*(SWA3MY+9.0D0*SWA2MY*rmax+9.0D0*rmin*SWB2MY+SWB3MY)!/D7MY
      SWC3MY= 140.0D0*(SWA3MY*rmax+3.0D0*SWA2MY*SWB2MY+rmin*SWB3MY)!/D7MY
      SWC2MY=-210.0D0*(SWA3MY*SWB2MY+SWA2MY*SWB3MY)!/D7MY
      SWC1MY= 140.0D0*SWA3MY*SWB3MY!/D7MY
      SWC0MY=(-35.0D0*SWA3MY*SWB2MY*SWB2MY+21.0D0*SWA2MY*SWB3MY*SWB2MY- &
      & 7.0D0*rmin*SWB3MY*SWB3MY+SWB3MY*SWB3MY*rmax)!/D7MY

      R2 = R*R
      R3 = R2*R

      f = (SWC7MY*R3*R3*R+SWC6MY*R3*R3+SWC5MY*R3*R2+SWC4MY*R2*R2+SWC3MY*R3+SWC2MY*R2+ &
          SWC1MY*R+SWC0MY) / D7MY
      dfdr = (7.0D0*SWC7MY*R3*R3+6.0D0*SWC6MY*R3*R2+5.0D0*SWC5MY*R2*R2+ &
           & 4.0D0*SWC4MY*R3+3.0D0*SWC3MY*R2+2.0D0*SWC2MY*R+SWC1MY) / D7MY  
    endif
  else ! lnorm = false


    if (R .le. rmin) then
!*************
!  r < rmin  *
!*************
      f = 0.0D0
      dfdr = 0.0D0
    elseif (R .ge. rmax) then
!*************
!  r > rmax  *
!*************
      f = 1.0D0
      dfdr = 0.0D0
    else
!********************
!  rmin < r < rmax  *
!********************
      D1MY=rmin-rmax
      D7MY=D1MY**7
      SWA2MY=rmax*rmax
      SWA3MY=SWA2MY*rmax
      SWB2MY=rmin*rmin
      SWB3MY=SWB2MY*rmin

      SWC7MY= 20.0D0!/D7MY
      SWC6MY= -70.0D0*(rmin+rmax)!/D7MY
      SWC5MY=  84.0D0*(SWA2MY+3.0D0*rmin*rmax+SWB2MY)!/D7MY
      SWC4MY= -35.0D0*(SWA3MY+9.0D0*SWA2MY*rmin+9.0D0*rmax*SWB2MY+SWB3MY)!/D7MY
      SWC3MY= 140.0D0*(SWA3MY*rmin+3.0D0*SWA2MY*SWB2MY+rmax*SWB3MY)!/D7MY
      SWC2MY=-210.0D0*(SWA3MY*SWB2MY+SWA2MY*SWB3MY)!/D7MY
      SWC1MY= 140.0D0*SWA3MY*SWB3MY!/D7MY
      SWC0MY=(-35.0D0*SWA3MY*SWB2MY*SWB2MY+21.0D0*SWA2MY*SWB3MY*SWB2MY- &
            & 7.0D0*rmax*SWB3MY*SWB3MY+SWB3MY*SWB3MY*rmin)!/D7MY

      R2 = R*R
      R3 = R2*R

      f = (SWC7MY*R3*R3*R+SWC6MY*R3*R3+SWC5MY*R3*R2+SWC4MY*R2*R2+SWC3MY*R3+SWC2MY*R2+ &
          SWC1MY*R+SWC0MY) / D7MY
      dfdr = (7.0D0*SWC7MY*R3*R3+6.0D0*SWC6MY*R3*R2+5.0D0*SWC5MY*R2*R2+ &
           & 4.0D0*SWC4MY*R3+3.0D0*SWC3MY*R2+2.0D0*SWC2MY*R+SWC1MY) / D7MY
    endif
  endif
 return
end subroutine valtaper
!***********************************************************************
!***********************************************************************

    subroutine taper(r,r2)

!***********************************************************************
    include 'cbka.blk'
!***********************************************************************
!                                                                      *
!     Taper function for Coulomb interaction                           *
!                                                                      *
!***********************************************************************
!     if (ndebug.eq.1) write (65,*) 'In taper'

      r3=r2*r
      SW=SWC7*R3*R3*R+SWC6*R3*R3+SWC5*R3*R2+SWC4*R2*R2+SWC3*R3+SWC2*R2+ &
      SWC1*R+SWC0
      SW1=7.0D0*SWC7*R3*R3+6.0D0*SWC6*R3*R2+5.0D0*SWC5*R2*R2+ &
      & 4.0D0*SWC4*R3+THREE*SWC3*R2+TWO*SWC2*R+SWC1
    return
    end subroutine taper
!***********************************************************************
!***********************************************************************

    subroutine tap7th

!***********************************************************************
    include 'cbka.blk'
!***********************************************************************
!                                                                      *
!     7th order taper function setup                                   *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In tap7th'
        call timer(65)
        close (65)
    end if
    D1=SWB-SWA
    D7=D1**7.0D0
    SWA2=SWA*SWA
    SWA3=SWA2*SWA
    SWB2=SWB*SWB
    SWB3=SWB2*SWB
     
!***********************************************************************
!     7th order taper function                                         *
!***********************************************************************
     
    SWC7=  20.0D0/D7
    SWC6= -70.0D0*(SWA+SWB)/D7
    SWC5=  84.0D0*(SWA2+3.0D0*SWA*SWB+SWB2)/D7
    SWC4= -35.0D0*(SWA3+9.0D0*SWA2*SWB+9.0D0*SWA*SWB2+SWB3)/D7
    SWC3= 140.0D0*(SWA3*SWB+3.0D0*SWA2*SWB2+SWA*SWB3)/D7
    SWC2=-210.0D0*(SWA3*SWB2+SWA2*SWB3)/D7
    SWC1= 140.0D0*SWA3*SWB3/D7
    SWC0=(-35.0D0*SWA3*SWB2*SWB2+21.0D0*SWA2*SWB3*SWB2- &
    & 7.0D0*SWA*SWB3*SWB3+SWB3*SWB3*SWB)/D7
     
    return
    end subroutine tap7th
!***********************************************************************
!***********************************************************************
  subroutine ataper(lnorm,r,rmin,rmax,f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,lgrad3)
!
!  Subroutine to calculate the taper function and derivatives 
!  using interpolation.
!
!  On entry : 
!
!  lnorm           = if .true. then apply normal taper
!                    if .false. then make the taper go to zero at the other end
!  r               = current distance for which the taper function
!                    is to be calculated
!  rmin            = minimum distance of taper
!  rmax            = maximum distance of taper
!  lgrad1          = if .true. calculate the first derivative
!  lgrad2          = if .true. calculate the second derivative
!  lgrad3          = if .true. calculate the third derivative
!
!  On exit :
!
!  f               = the value of the taper function at r
!  dfdr            = first derivative of taper function
!  d2fwdr2         = second derivative of taper function
!  d3fwdr3         = third derivative of taper function

  implicit none
!
!  Passed variables
!
  real*8,    intent(in)             :: r
  real*8,    intent(in)             :: rmin
  real*8,    intent(in)             :: rmax
  real*8,    intent(out)            :: f
  real*8,    intent(out)            :: dfdr
  real*8,    intent(out)            :: d2fdr2
  real*8,    intent(out)            :: d3fdr3
  logical,     intent(in)           :: lgrad1
  logical,     intent(in)           :: lgrad2
  logical,     intent(in)           :: lgrad3
  logical,     intent(in)           :: lnorm
!
!  Local variables
!
  real*8                            :: rd
  real*8                            :: rterm
  real*8                            :: trm1
!

  if (lnorm) then
    if (r .le. rmin) then
!*************
!  r < rmin  *
!*************
      f = 1.0D0
      if (lgrad1) then
        dfdr = 0.0D0
        if (lgrad2) then
          d2fdr2 = 0.0D0
          if (lgrad3) then
            d3fdr3 = 0.0D0
          endif
        endif
      endif
    elseif (r .ge. rmax) then
!*************
!  r > rmax  *
!*************
      f = 0.0D0
      if (lgrad1) then
        dfdr = 0.0D0
        if (lgrad2) then
          d2fdr2 = 0.0D0
          if (lgrad3) then
            d3fdr3 = 0.0D0
          endif
        endif
      endif
    else
!********************
!  rmin < r < rmax  *
!********************
      rterm = 1.0D0/(rmax - rmin)**3
      rd = rmax - r
      trm1 = (rmax + 2.0D0*r - 3.0D0*rmin)
      f = rterm*rd*rd*trm1
      if (lgrad1) then
        dfdr = 2.0D0*rd*rterm*(rd - trm1)
        if (lgrad2) then
          d2fdr2 = - 2.0D0*rterm*(4.0D0*rd - trm1)
          if (lgrad3) then
            d3fdr3 = 12.0D0*rterm
          endif
        endif
      endif
    endif
  else
    if (r .le. rmin) then
!*************
!  r < rmin  *
!*************
      f = 0.0D0
      if (lgrad1) then
        dfdr = 0.0D0
        if (lgrad2) then
          d2fdr2 = 0.0D0
          if (lgrad3) then
            d3fdr3 = 0.0D0
          endif
        endif
      endif
    elseif (r .ge. rmax) then
!*************
!  r > rmax  *
!*************
      f = 1.0D0
      if (lgrad1) then
        dfdr = 0.0D0
        if (lgrad2) then
          d2fdr2 = 0.0D0
          if (lgrad3) then
            d3fdr3 = 0.0D0
          endif
        endif
      endif
    else
!********************
!  rmin < r < rmax  *
!********************
      rterm = 1.0D0/(rmin - rmax)**3
      rd = rmin - r
      trm1 = (rmin + 2.0D0*r - 3.0D0*rmax)
      f = rterm*rd*rd*trm1
      if (lgrad1) then
        dfdr = 2.0D0*rd*rterm*(rd - trm1)
        if (lgrad2) then
          d2fdr2 = - 2.0D0*rterm*(4.0D0*rd - trm1)
          if (lgrad3) then
            d3fdr3 = 12.0D0*rterm
          endif
        endif
      endif
    endif
  endif

  return
  end subroutine ataper
!***********************************************************************

    subroutine default

!***********************************************************************
    include 'cbka.blk'
!***********************************************************************
!                                                                      *
!     Add default atoms to a molecular structure                       *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In default'
        call timer(65)
        close (65)
    end if
    inplc=0
    call distan
    call srtbon1
    call intcor
    call mdsav(0,qfile(nprob))
    do i1=1,na
        ity1=ia(i1,1)
        icon=0
        do i2=1,nsbmax
            if (iag(i1,2+i2) > 0) icon=icon+1
        end do
        if (icon < int(aval(ity1)) .AND. idef(ity1) /= 0 .AND. &
        aval(ity1) < 5.0) then
            inplc=inplc+1
            ity2=idef(ity1)
            ia(na+inplc,1)=ity2
            qa(na+inplc)=qas(ity2)
            ity3=ia(iag(i1,3),1)
            rrnew=0.5d0*(rat(ity1)+rat(ity2))
            rrold=0.5d0*(rat(ity1)+rat(ity3))
            if (ity3 < ity2) then
                ity3=idef(ity1)
                ity2=ia(iag(i1,3),1)
            end if
            do i2=1,nvaty
                if (ity2 == nvs(i2,1) .AND. ity1 == nvs(i2,2) .AND. &
                ity3 == nvs(i2,3)) then
                    angnew=th0(i2)
                end if
            end do
            r3=sqrt(rrnew*rrnew+rrold*rrold- &
            & 2.0d0*rrnew*rrold/cos(dgrrdn*angnew))
            cxd=c(i1,1)-c(iag(i1,3),1)
            cyd=c(i1,2)-c(iag(i1,3),2)
            czd=c(i1,3)-c(iag(i1,3),3)
            call dista2(i1,iag(i1,3),dis,cxd,cyd,czd)
            dfg=rrnew/rrold
            c(na+inplc,1)=c(i1,1)+dfg*cxd
            c(na+inplc,2)=c(i1,2)+dfg*cyd
            c(na+inplc,3)=c(i1,3)+dfg*czd
                  
        end if

    end do
    na=na+inplc
    call distan
    call srtbon1
    call intcor
    call mdsav(1,qfile(nprob))

    return
    end subroutine default
!***********************************************************************
!***********************************************************************

    subroutine staint

!***********************************************************************
    include 'cbka.blk'
    dimension bvt(nat,4)
!***********************************************************************
!                                                                      *
!     Generate cartesian coordinates from internal coordinate input    *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In staint'
        call timer(65)
        close (65)
    end if
    k=0
    10 read (3,1200,end=20,err=20)(ijk(k+1,k1),k1=1,3),k2,qa(k+1), &
    bvt(k+1,3),bvt(k+1,2),bvt(k+1,1)
    qlabel(k+1)=qa(k+1)
    qresi1(k+1)='   '
    qresi2(k+1)=' '
    qresi3(k+1)='     '
    qffty(k+1)='     '
    if (k2 /= k+1) then
        write (*,*)'Wrong order in internal coordinates at atom:',k2
        goto 20
    !     stop 'Wrong order in internal coordinates'
    end if
    k=k+1
    if (k > nat) then
        write (*,*)na,nat
        stop 'Maximum number of atoms exceeded'
    end if
    goto 10
    20 continue
    na=k
          
!***********************************************************************
!                                                                      *
!     CALCULATION OF CARTESIAN COORDINATES FROM INTERNAL COORDINAATES  *
!                                                                      *
!***********************************************************************

    12 C(1,1)=ZERO
    C(1,2)=ZERO
    C(1,3)=ZERO
    C(2,1)=BVT(2,1)
    C(2,2)=ZERO
    C(2,3)=ZERO
    HR=(BVT(3,2)-90.0D0)*DGRRDN
    C(3,1)=C(2,1)+BVT(3,1)*SIN(HR)
    C(3,2)=BVT(3,1)*COS(HR)
    C(3,3)=ZERO
    DO  K1=4,NA
        J=IJK(K1,2)
        KB=K1-1
        XH=C(J,1)
        YH=C(J,2)
        ZH=C(J,3)
        DO  K2=1,KB
            C(K2,1)=C(K2,1)-XH
            C(K2,2)=C(K2,2)-YH
            C(K2,3)=C(K2,3)-ZH
            DO  K3=1,3
                IF (ABS(C(K2,K3)) < 1.0D-15) C(K2,K3)=1.0D-15
            END DO
        ENDDO
        K=IJK(K1,3)
        P2=C(K,2)*C(K,2)+C(K,3)*C(K,3)
        IF (P2 /= ZERO) THEN
            P=SQRT(P2)
            Q=SQRT(C(K,1)*C(K,1)+P2)
            SA=C(K,2)/P
            CA=C(K,3)/P
            SB=-C(K,1)/Q
            CB=P/Q
        ELSE
            SA=ZERO
            CA=ONE
            SB=ONE
            CB=ZERO
        ENDIF
        DO  K2=1,KB
            AZ=C(K2,1)
            BZ=C(K2,2)
            C(K2,1)=AZ*CB+BZ*SB*SA+C(K2,3)*SB*CA
            C(K2,2)=BZ*CA-C(K2,3)*SA
            C(K2,3)=-AZ*SB+BZ*CB*SA+C(K2,3)*CB*CA
        END DO
        IF (C(K,3) <= ZERO) THEN
            DO  K2=1,KB
                C(K2,3)=-C(K2,3)
            END DO
        ENDIF
        I=IJK(K1,1)
        IF (1.0D5*ABS(C(I,1)) <= ABS(C(I,2))) THEN
            T1=HALF*PI
        ELSE
            YX=ABS(C(I,2)/C(I,1))
            T1=ATAN(YX)
        ENDIF
        IF (C(I,1) >= ZERO .AND. C(I,2) < ZERO) T1=TWO*PI-T1
        IF (C(I,1) < ZERO .AND. C(I,2) >= ZERO) T1=PI-T1
        IF (C(I,1) < ZERO .AND. C(I,2) < ZERO) T1=T1+PI
        DO  K2=1,KB
            IF (C(K2,1) == ZERO .AND. C(K2,2) == ZERO) GOTO 31
            IF (1.0D5*ABS(C(K2,1)) < ABS(C(K2,2))) THEN
                T2=HALF*PI
            ELSE
                YX=ABS(C(K2,2)/C(K2,1))
                T2=ATAN(YX)
            ENDIF
            IF (C(K2,1) >= ZERO .AND. C(K2,2) < ZERO) T2=TWO*PI-T2
            IF (C(K2,1) < ZERO .AND. C(K2,2) >= ZERO) T2=PI-T2
            IF (C(K2,1) < ZERO .AND. C(K2,2) < ZERO) T2=T2+PI
            T3=T2-T1
            IF (T3 < ZERO)T3=T3+TWO*PI
            RZ=SQRT(C(K2,1)*C(K2,1)+C(K2,2)*C(K2,2))
            C(K2,1)=RZ*COS(T3)
            C(K2,2)=RZ*SIN(T3)
        31 END DO
        HR=(BVT(K1,2)-90.0D0)*DGRRDN
        HT=BVT(K1,3)*DGRRDN
        CHR=COS(HR)
        C(K1,1)=BVT(K1,1)*CHR*COS(HT)
        C(K1,2)=BVT(K1,1)*CHR*SIN(HT)
        C(K1,3)=C(IJK(K1,3),3)+BVT(K1,1)*SIN(HR)
    END DO
          
    return
    1200 FORMAT(4I3,1X,A2,3F10.5,4X,I1,F10.5)
    end subroutine staint
!***********************************************************************
!***********************************************************************

    subroutine outint

!***********************************************************************
    include 'cbka.blk'
!***********************************************************************
!                                                                      *
!     Output internal coordinates                                      *
!                                                                      *
!***********************************************************************
    dimension dvdc(3,3),dargdc(3,3)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In outint'
        call timer(65)
        close (65)
    end if
    write (91,50)qmol
    open (82,file='output.MOP',status='unknown')
    write (82,*)
    write (82,'(a40)')qmol
    write (82,*)
    close (82)

!     IF (NMOLO.GT.1) THEN
!     WRITE(6,*)' OUTPUT INTERNAL COORDINATES NOT POSSIBLE FOR CALCULATI
!    $ON ON MORE THAN ONE MOLECULE'
!     RETURN
!     END IF

!***********************************************************************
!                                                                      *
!     Output of internal coordinates.                                  *
!     First 3 atoms of other input file.                               *
!                                                                      *
!***********************************************************************
    N1=1
    N2=2
    N3=3
    open (82,file='output.MOP',status='unknown',position='append')
    write(91,100)N1,qa(n1)
    write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(n1), &
    zero,nzero,zero,nzero,zero,nzero,nzero,nzero,nzero
    call dista2(n1,n2,rr,dx,dy,dz)
    write(91,200)N1,N2,qa(n2),RR
    write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(n2), &
    rr,none,zero,nzero,zero,nzero,n1,nzero,nzero
    close (82)

    call dista2(n2,n3,rr,dx,dy,dz)
    hv=zero
    call calvalres(n1,n2,n3,arg,hv,dvdc,dargdc)
    WRITE(91,300)N1,N2,N3,qa(n3),rdndgr*HV,RR
    open (82,file='output.MOP',status='unknown',position='append')
    write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(n3), &
    rr,none,rdndgr*hv,none,zero,nzero,n2,n1,nzero
    close (82)

    naih=3
          
    do i1=naih+1,na
        bomax=zero
        j1=0
        do i2=1,ia(i1,2)
            iob=ia(i1,2+i2)
            ncubo=nubon2(i1,i2)
            if (bo(ncubo) > bomax .AND. iob < i1) then
                bomax=bo(ncubo)
                j1=iob
            end if
        end do
        if (j1 == 0) j1=i1-1
        call dista2(j1,i1,rr,dx,dy,dz)
              
        bomax=zero
        j2=0
        do i2=1,ia(j1,2)
            iob=ia(j2,2+i2)
            ncubo=nubon2(j1,i2)
            if (bo(ncubo) > bomax .AND. iob < i1 .AND. &
            abo(iob) > bo(ncubo)+0.2d0) then
                bomax=bo(ncubo)
                j2=iob
            end if
        end do
        if (j2 == 0) j2=i1-2
        if (j2 == j1) j2=j2+1

        call calvalres(j2,j1,i1,arg,hh,dvdc,dargdc)

        bomax=zero
        j3=0
        do i2=1,ia(j2,2)
            iob=ia(j2,2+i2)
            ncubo=nubon2(j2,i2)
            if (bo(ncubo) > bomax .AND. iob < i1 .AND. iob /= j1) then
                bomax=bo(ncubo)
                j3=iob
            end if
        end do
        if (j3 == 0) j3=i1-3
        if (j3 == j2 .AND. j3 /= j1-1) j3=j3+1
        if (j3 == j2 .AND. j3 /= j1-2) j3=j3+2
        if (j3 == j1 .AND. j3 /= j2-1) j3=j3+1
        if (j3 == j1 .AND. j3 /= j2-2) j3=j3+2

        call caltor(j3,j2,j1,i1,ht)

        write(91,400)j3,j2,j1,i1,qa(i1),ht,rdndgr*hh,rr
        open (82,file='output.MOP',status='unknown',position='append')
        write(82,'(2x,a2,f12.6,i3,f12.6,i3,f12.6,i3,1x,3i4)')qa(i1), &
        rr,none,rdndgr*hh,none,ht,none,j1,j2,j3
        close (82)
    end do

    close(82)
    return
    50 format ('  I',2x,a60)
    100 FORMAT(9X,I3,1x,a2)
    200 FORMAT(6X,2I3,1x,a2,20X,F10.5)
    300 FORMAT(3X,3I3,1x,a2,10X,2F10.5)
    400 FORMAT(4I3,1x,a2,3F10.5)

    end subroutine outint
!***********************************************************************
!***********************************************************************

    subroutine molanal(igno,cutofm,iname,ntotmol)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension iam(nat,mbond+3),nmolata(nmolmax,nat)
    dimension molfra(nmolmax,nsort),ndup(nmolmax)
    dimension vmolmas(nmolmax)
    character(20) :: qname
    character(40) :: qmolan1
    character(100) :: qmolan
    character(6) :: qremove
    character(7) :: var
    logical :: found
!***********************************************************************
!                                                                      *
!     Analyse and output molecular fragments                           *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In molanal'
        call timer(65)
        close (65)
    end if

    qremove='remove'
    if (iname == 1) qname="molfra.out"
    if (iname == 2) qname="molfra_ig.out"
    if (iname == 3) qname="molfra2.out"
    if (iname == 4) qname="molfra3.out"
    if (iname == 5) qname="molfra4.out"
    if (iname == 6) qname="molfra5.out"
    if (iname == 7) qname="molfra6.out"
    do i1=1,nmolmax
        do i2=1,nsort
            molfra(i1,i2)=0
        end do
        ndup(i1)=1
    end do

    do i1=1,na
        do i2=1,mbond+3
            iam(i1,i2)=0
        end do
    end do
!***********************************************************************
!                                                                      *
!     Create connection table based on corrected bond orders           *
!                                                                      *
!***********************************************************************
    do i1=1,nbon
        if (bo(i1) > cutofm) then
            j1=ib(i1,2)
            j2=ib(i1,3)
            ity1=ia(j1,1)
            ity2=ia(j2,1)
            if (ity1 /= igno .AND. ity2 /= igno) then
                iam(j1,2)=iam(j1,2)+1
                iam(j1,2+iam(j1,2))=j2
                iam(j2,2)=iam(j2,2)+1
                iam(j2,2+iam(j2,2))=j1
            end if
        end if
    end do
!*********************************************************************
!                                                                    *
!     Find molecules                                                 *
!                                                                    *
!*********************************************************************
    nmolo6=0
    found= .FALSE. 
    DO 61 k1=1,na
        IF (iam(K1,3+mbond) == 0) found= .TRUE. 
        IF (iam(K1,3+mbond) > nmolo6) nmolo6=iam(K1,3+mbond)
    61 END DO
    IF ( .NOT. FOUND) GOTO 62
!***********************************************************************
!                                                                      *
!     Molecule numbers are assigned. No restrictions are made for the  *
!     sequence of the numbers in the connection table.                 *
!                                                                      *
!***********************************************************************
    N3=1
    64 N2=N3
    nmolo6=nmolo6+1
    if (nmolo6 > nmolmax) stop 'Too many molecules in system'
    iam(N2,3+mbond)=nmolo6
    67 FOUND= .FALSE. 
    DO 66 N1=N2+1,na
        IF (iam(N1,3+mbond) /= 0) GOTO 66
        DO 65 L=1,mbond
            IF (iam(N1,l+2) == 0) GOTO 66
            IF (iam(iam(N1,l+2),3+mbond) == nmolo6) THEN
                FOUND= .TRUE. 
                iam(N1,3+mbond)=nmolo6
                GOTO 66
            ENDIF
        65 END DO
    66 END DO
    IF (FOUND) GOTO 67
    DO 63 N3=N2+1,NA
        if (iam(N3,3+mbond) == 0) goto 64
    63 END DO
!***********************************************************************
!                                                                      *
!     The assigned or input molecule numbers are checked for their     *
!     consistency.                                                     *
!                                                                      *
!***********************************************************************
    62 FOUND= .FALSE. 
    DO 72 N1=1,NA
        DO 71 L=1,mbond
            IF (iam(N1,L+2) == 0) GOTO 72
            IF (iam(iam(N1,L+2),3+mbond) /= iam(N1,3+mbond)) THEN
                FOUND= .TRUE. 
            ENDIF
        71 END DO
    72 END DO
    IF (FOUND) THEN
        write (7,'(i4,a40)')na,qmol
        do i1=1,na
            write (7,'(40i4)')i1,iam(i1,1),(iam(i1,2+i2),i2=1,nsbmax), &
            iam(i1,3+mbond)
        end do
        STOP 'Mol.nrs. not consistent; maybe wrong cell parameters'
    ENDIF

    do i1=1,nmolo6
        natmol=0
        do i2=1,na
            if (iam(i2,3+mbond) == i1) then
                natmol=natmol+1
                nmolata(i1,natmol+1)=i2
            end if
        end do
        nmolata(i1,1)=natmol
    end do
!***********************************************************************
!                                                                      *
!     Analyze molecules                                                *
!                                                                      *
!***********************************************************************
    do i1=1,nmolo6
        do i2=1,nmolata(i1,1)
            i3=nmolata(i1,1+i2)
            ityp=ia(i3,1)
            molfra(i1,ityp)=molfra(i1,ityp)+1
        end do
    end do

    do i1=1,nmolo6
        isee=0
        do i2=1,nmolo6
            isee2=1
            do i3=1,nso
                if (molfra(i1,i3) /= molfra(i2,i3)) isee2=0
            end do
            if (isee2 == 1 .AND. i1 > i2 .AND. isee == 0) then  !molecule type already exists
                ndup(i2)=ndup(i2)+1
                ndup(i1)=0
                isee=1
            end if

        end do
    end do

    open (45,file=qname,status='unknown',position='append')
    if (mdstep == 0 .AND. igno > 0) write (45,100)cutofm,qas(igno)
    if (mdstep == 0 .AND. igno == 0) write (45,101)cutofm
    write (45,110)
    ntotmol=0
    ntotat=0
    vtotmass=zero
    do i1=1,nmolo6
        if (ndup(i1) > 0) then
        !     write (45,110)i1,(molfra(i1,i2),i2=1,nso),ndup(i1)
            ntotmol=ntotmol+ndup(i1)
            qmolan=' '
            qmolan1=' '
            istart=-4
            ihulp=0
            vmass=zero
            do i2=1,nso
                vmass=vmass+molfra(i1,i2)*amas(i2)
                ntotat=ntotat+molfra(i1,i2)*ndup(i1)
                if (molfra(i1,i2) > 0) then
                    istart=istart+6
                    iend=istart+5
                    if (molfra(i1,i2) > 1) then
                        write (qmolan(istart:iend),'(a2,i3)')qas(i2),molfra(i1,i2)
                    else
                        write (qmolan(istart:iend-2),'(a2)')qas(i2)
                    end if
                end if
            end do
            ihulp=1
            do i2=1,iend
                if (qmolan(i2:i2) /= ' ') then
                    qmolan1(ihulp:ihulp)=qmolan(i2:i2)
                    ihulp=ihulp+1
                end if
            end do

        !     write (45,120)ndup(i1),qmolan(1:iend),vmass
            write (45,120)mdstep+nprevrun,ndup(i1),qmolan1,vmass
            vtotmass=vtotmass+ndup(i1)*vmass
        end if
    end do
    write (45,*)'Total number of molecules:',ntotmol
    write (45,*)'Total number of atoms:',ntotat
    write (45,*)'Total system mass:',vtotmass
    close (45)

    if (ibocha == 1) then
    !***********************************************************************
    !                                                                      *
    !     Analyze changes in connection table                              *
    !                                                                      *
    !***********************************************************************
        open (46,file='bondchange.out',status='unknown',position='append')
        do i1=1,na
            if (iag(i1,2) < iagold(i1,2) .AND. mdstep > 1) &
            write (46,130)mdstep,i1
            if (iag(i1,2) > iagold(i1,2) .AND. mdstep > 1) &
            write (46,140)mdstep,i1
            do i2=1,mbond+3
                iagold(i1,i2)=iag(i1,i2)
            end do
        end do
        close (46)
    end if

    if (mdstep > 0 .AND. iremov > 0 .AND. mod(mdstep,iremov) == 0) then
    !***********************************************************************
    !                                                                      *
    !     Remove molecular fragments from system                           *
    !                                                                      *
    !***********************************************************************
        write (var,'(f7.4)')float(mdstep/iremov)/1d4
    !***********************************************************************
    !                                                                      *
    !     Find molecular masses                                            *
    !                                                                      *
    !***********************************************************************
        do i1=1,nmolo6
            vmolmas(i1)=zero
        end do
        do i1=1,na
            imolh=iam(i1,3+mbond)
            vmolmas(imolh)=vmolmas(imolh)+xmasat(i1)
        end do
            
    !     do i1=1,nmolo6
    !     write (*,*)i1,vmolmas(i1)
    !     end do
            
        natrem=0
        vstrrem=zero
        do i1=1,nmolo6   !Check how many atoms will be removed for .xyz output header
            if (vmolmas(i1) <= vmamin .AND. vmolmas(i1) > vmami2) then
                do i2=1,na
                    if (iam(i2,3+mbond) == i1) natrem=natrem+1
                    if (iam(i2,3+mbond) == i1) vstrrem=vstrrem+estrain(i2)
                end do
            end if
        end do

        if (natrem > 0) then
            open (64,file=qremove//var(3:7),status='unknown',position='append')
            write (64,'(i8)')natrem
            write (64,*)'Atoms_removed_in_iteration         ',mdstep,vstrrem
        end if

        do i1=1,nmolo6   !Remove the atoms
            if (vmolmas(i1) <= vmamin .AND. vmolmas(i1) > vmami2) then
                write (*,*)'Iteration',mdstep,' :removing molecule', i1, &
                ' with mass',vmolmas(i1)
                35 continue
                do i2=1,na

                    if (iam(i2,3+mbond) == i1) then
                        write (64,*)qa(i2),c(i2,1),c(i2,2),c(i2,3),iam(i2,3+mbond), &
                        estrain(i2)+vincr(ia(i2,1))

                        do i3=i2,na-1
                            c(i3,1)=c(i3+1,1)
                            c(i3,2)=c(i3+1,2)
                            c(i3,3)=c(i3+1,3)
                            vel(1,i3)=vel(1,i3+1)
                            vel(2,i3)=vel(2,i3+1)
                            vel(3,i3)=vel(3,i3+1)
                            accel(1,i3)=accel(1,i3+1)
                            accel(2,i3)=accel(2,i3+1)
                            accel(3,i3)=accel(3,i3+1)
                            qa(i3)=qa(i3+1)
                            xmasat(i3)=xmasat(i3+1)
                            iam(i3,1)=iam(i3+1,1)
                            ia(i3,1)=ia(i3+1,1)
                            iam(i3,3+mbond)=iam(i3+1,3+mbond)
                            estrain(i3)=estrain(i3+1)
                        end do
                        na=na-1
                        namov=namov-1           !January 24 2009 modification
                        goto 35   !start with a new loop
                    end if

                end do
                call intcor  !Build new connection table
                call encalc

            !     stop 'Trying to remove fragments'
            end if
             
        end do
        close (64)

    end if

    return
    100 format('Bond order cutoff:',f6.4, &
    ' Ignoring bonding from atom type:',a2)
    101 format('Bond order cutoff:',f6.4)
    110 format('Iteration Freq. Molecular formula',15x,'Molecular mass')
    120 format(i8,i4,' x  ',a35,f10.4)
    130 format('MD-step:',i10,' Atom ',i6,' has lost bond(s)')
    140 format('MD-step:',i10,' Atom ',i6,' has gained bond(s)')
    end subroutine molanal
!***********************************************************************
!***********************************************************************

    subroutine runanal

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension iconn(10)
!***********************************************************************
!                                                                      *
!     Perform runtime analysis                                         *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In runanal'
        call timer(65)
        close (65)
    end if
!***********************************************************************
!                                                                      *
!     Output number of bonds around hydrogen atoms to unit 64          *
!                                                                      *
!***********************************************************************
    do i1=1,10
        iconn(i1)=0
    end do
    do i1=1,na
        if (qa(i1) == 'H ') then
            iconh=iag(i1,2)
            iconn(iconh)=iconn(iconh)+1
        end if
    end do

    return
    end subroutine runanal
!***********************************************************************
!***********************************************************************

    subroutine analysis(naold)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension sbotype(nat,nsort),nbotype(nat,nsort),sbot(nat)
    character(2) :: qaan,qva1,qva2,qva3
    character(1) :: qrom
    character(10) :: qfileout
!***********************************************************************
!                                                                      *
!     Perform post-run analysis                                        *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In analysis'
        call timer(65)
        close (65)
    end if
!***********************************************************************
!                                                                      *
!     Give options                                                     *
!                                                                      *
!***********************************************************************
    10 write (6,*)'Analysis options:'
    write (6,*)'1) Bond orders/total bond order around an atom type'
    write (6,*)'2) Valence angles'
    write (6,*)'3) End analysis'
    write (6,*)'Which analysis option do you want to use ?'
    read (5,*) ioption
    if (ioption < 1 .OR. ioption > 3) then
        write (6,*)'Invalid option'
        goto 10
    end if
    goto (100,200,300) ioption
!***********************************************************************
!                                                                      *
!     Bond order analysis                                              *
!                                                                      *
!***********************************************************************
    100 continue
    write (6,*)'Give atom type as in ffield-file'
    read (5,'(a2)') qaan
    write (6,*)'Bond order cutoff'
    read (5,*) bocu
    write (6,*)'Give output file name (< 10 characters)'
    read (5,'(a10)') qfileout
    iatyp=0
    do i1=1,nso
        if (qaan == qas(i1)) iatyp=i1
    end do
    if (iatyp == 0) then
        write (6,*)'Unknown atom type'
        goto 100
    end if

    if (na == 0) na=naold
    do i1=1,na
        do i2=1,nso
            nbotype(i1,i2)=0
            sbotype(i1,i2)=zero
        end do
        sbot(i1)=zero
    end do
    do i1=1,na
        do i2=1,ia(i1,2)
            iob=ia(i1,2+i2)
            ncubo=nubon2(i1,i2)
            ityp=ia(iob,1)
            if (bo(ncubo) > bocu) then
                nbotype(i1,ityp)=nbotype(i1,ityp)+1
                sbotype(i1,ityp)=sbotype(i1,ityp)+bo(ncubo)
            end if
            sbot(i1)=sbot(i1)+bo(ncubo)
        end do
    end do
    open (66,file=qfileout,status='new',err=500)
    write (66,1000)qas(iatyp)
    write (66,1005)bocu
    write (66,1010)(qas(i1),i1=1,nso),(qas(i1),i1=1,nso)
    do i2=1,na
        if (ia(i2,1) == iatyp) then
            if (nso == 1) write (66,'(i6,3f12.4,1i4,2f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 2) write (66,'(i6,3f12.4,2i4,3f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 3) write (66,'(i6,3f12.4,3i4,6f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 4) write (66,'(i6,3f12.4,4i4,5f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 5) write (66,'(i6,3f12.4,5i4,6f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 6) write (66,'(i6,3f12.4,6i4,7f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 7) write (66,'(i6,3f12.4,7i4,8f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 8) write (66,'(i6,3f12.4,8i4,9f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 9) write (66,'(i6,3f12.4,9i4,10f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 10) write (66,'(i6,3f12.4,10i4,11f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 11) write (66,'(i6,3f12.4,11i4,12f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 12) write (66,'(i6,3f12.4,12i4,13f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 13) write (66,'(i6,3f12.4,13i4,14f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 14) write (66,'(i6,3f12.4,14i4,15f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)
            if (nso == 15) write (66,'(i6,3f12.4,15i4,16f12.4)')i2,c(i2,1), &
            c(i2,2),c(i2,3),(nbotype(i2,i3),i3=1,nso), &
            (sbotype(i2,i3),i3=1,nso),sbot(i2)

            if (nso > 15) then
                write (*,*)'Atom types in force field > 15, causing problems.'
                write (*,*)'Modify subroutine analysis in reac.f'
                stop 'End in subroutine analysis'
            end if

        end if
    end do
    close (66)
    write (6,*)'Do you want to do another analysis (y/n)'
    read (5,'(a1)') qrom
    if (qrom == 'y') goto 10
    goto 300
!***********************************************************************
!                                                                      *
!     Valence angle analysis                                           *
!                                                                      *
!***********************************************************************
    200 continue
    write (6,*)'Give atom types in valence angle as in ffield-file'
    write (6,*)'(e.g. C C H for C-C-H angles)'
    read (5,'(a)')qstrana1
    write (6,*)'Bond order cutoff'
    read (5,*) bocu
    write (6,*)'Give output file name (< 10 characters)'
    read (5,'(a10)') qfileout
!***********************************************************************
!                                                                      *
!     Determine force field type angle                                 *
!                                                                      *
!***********************************************************************
    istart=1
    call stranal(istart,iend,vout,iout,1)
    qva1=qstrana2(1:2)
    istart=iend
    iatyp1=0
    do i1=1,nso
        if (qva1 == qas(i1)) iatyp1=i1
    end do
    if (iatyp1 == 0) then
        write (6,*)'Unknown atom type for first atom'
        goto 200
    end if
    call stranal(istart,iend,vout,iout,1)
    qva2=qstrana2(1:2)
    istart=iend
    iatyp2=0
    do i1=1,nso
        if (qva2 == qas(i1)) iatyp2=i1
    end do
    if (iatyp2 == 0) then
        write (6,*)'Unknown atom type for second atom'
        goto 200
    end if
    call stranal(istart,iend,vout,iout,1)
    qva3=qstrana2(1:2)
    iatyp3=0
    do i1=1,nso
        if (qva3 == qas(i1)) iatyp3=i1
    end do
    if (iatyp3 == 0) then
        write (6,*)'Unknown atom type for third atom'
        goto 200
    end if
    if (iatyp1 > iatyp3) then
        ihu=iatyp1
        iatyp1=iatyp3
        iatyp3=ihu
    end if
    ivty=0
    do i1=1,nvaty
        if (iatyp1 == nvs(i1,1) .AND. iatyp2 == nvs(i1,2) .AND. &
        iatyp3 == nvs(i1,3)) ivty=i1
    end do
    if (ivty == 0) then
        write (6,*)'Valence angle does not exist in force field'
        goto 200
    end if

    open (66,file=qfileout,status='new',err=500)
    write (66,1100)qva1,qva2,qva3
    write (66,1005)bocu
    write (66,1110)
    write (66,1120)
    do i1=1,nval
        if (iv(i1,1) == ivty) then
            la=iv(i1,5)
            lb=iv(i1,6)
            if (bo(la) > bocu .AND. bo(lb) > bocu) then
                write (66,'(3i5,10f12.4)')iv(i1,2),iv(i1,3),iv(i1,4), &
                c(iv(i1,3),1),c(iv(i1,3),2),c(iv(i1,3),3), &
                rdndgr*h(i1),bo(la),bo(lb)
            end if
        end if
    end do
    close (66)
    write (6,*)'Do you want to do another analysis (y/n)'
    read (5,'(a1)') qrom
    if (qrom == 'y') goto 10
    goto 300
    300 continue
    return
    500 stop 'Error opening file'
    1000 format ('Atom type:',a2,38x,'Nr. of bonds to',38x, &
    'Total bond order with',45x,'Total bond order')
    1005 format ('Bond order cutoff:',f10.4)
    1010 format ('   Atom#       x           y           z     ',8(a2,2x), &
    & 6x,8(a2,10x))
    1100 format ('Angles between atom types: ',3(a2,1x))
    1110 format ('                          position central atom')
    1120 format ('atom1 atom2 atom3       x           y           ', &
    'z        angle       BO(1-2)   BO(2-3)')
    end subroutine analysis
!***********************************************************************
!***********************************************************************

    subroutine outres

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
!***********************************************************************
!                                                                      *
!     Output molecular data                                            *
!                                                                      *
!***********************************************************************
    dimension isort(100),iad1(100),iad2(100),iad3(100),iad4(100)
    character(60) :: qm2
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In outres'
        call timer(65)
        close (65)
    end if
    read (9,100,end=50)idata,qm2
!     if (qm2.ne.qmol) then
!     write (*,*)'Wrong molecule in outres-file'
!     write (*,*)qmol
!     write (*,*)qm2
!     return
!     end if
    do 25 i1=1,idata
        read (9,200)isort(i1),iad1(i1),iad2(i1),iad3(i1),iad4(i1)
        ndata2=ndata2+1

        if (isort(i1) == 1) then
        !     do i2=1,nbon
        !     if (ib(i2,2).eq.iad1(i1).and.ib(i2,3).eq.iad2(i1)) then
        !     if (iopt.ne.1) write (81,*)iad1(i1),iad2(i1),rbo(i2)
        !     caldat(ndata2)=rbo(i2)
        !     end if
        !     end do
            call dista2(iad1(i1),iad2(i1),dish,dx,dy,dz)
            write (81,*)iad1(i1),iad2(i1),dish
            caldat(ndata2)=dish
        end if

        if (isort(i1) == 2) then
            do i2=1,nval
                if (iv(i2,2) == iad1(i1) .AND. iv(i2,3) == iad2(i1) .AND. &
                iv(i2,4) == iad3(i1)) then
                    if (iopt /= 1) write (81,*)iad1(i1),iad2(i1), &
                    iad3(i1),h(i2)*rdndgr
                    caldat(ndata2)=h(i2)*rdndgr
                end if
            end do
        end if

        if (isort(i1) == 3) then
            do i2=1,ntor
                if (it(i2,2) == iad1(i1) .AND. it(i2,3) == iad2(i1) .AND. &
                it(i2,4) == iad3(i1) .AND. it(i2,5) == iad4(i1)) then
                    if (iopt /= 1) write (81,*)iad1(i1),iad2(i1),iad3(i1),iad4(i1), &
                    abs(thg(i2))
                    caldat(ndata2)=abs(thg(i2))
                end if
            end do
        end if

        if (isort(i1) == 4) then
            if (iopt /= 1) write (81,*)estrmin
            caldat(ndata2)=estrmin
        end if

        if (isort(i1) == 5) then
            if (iopt /= 1) write (81,*)estrmin
            caldat(ndata2)=estrmin
        end if

        if (isort(i1) == 6) then
            if (iopt /= 1) write (81,*)iad1(i1),axiss(iad1(i1))
            caldat(ndata2)=axiss(iad1(i1))
        end if

        if (isort(i1) == 7) then
            if (iopt /= 1) write (81,*)eco
            caldat(ndata2)=eco
        end if

        if (isort(i1) == 8) then
            do i2=1,nbon
                if (ib(i2,2) == iad1(i1) .AND. ib(i2,3) == iad2(i1)) then
                    if (iopt /= 1) write (81,*)iad1(i1),iad2(i1),bo(i2)
                    caldat(ndata2)=bo(i2)
                end if
            end do
        end if

        if (isort(i1) == 9) then
            if (iopt /= 1) write (81,*)ch(iad1(i1))
            caldat(ndata2)=ch(iad1(i1))
        end if
              
        if (isort(i1) == 10) then
            rmsg=0.0
            nmovh=0
            do i2=1,na
                do i3=1,3
                    rmsg=rmsg+imove(i2)*d(i3,i2)*d(i3,i2)
                    nmovh=nmovh+imove(i2)
                end do
            end do
            rmsg=sqrt(rmsg/float(nmovh))

            if (iopt /= 1) write (81,*)rmsg
            caldat(ndata2)=rmsg
        end if

        if (isort(i1) == 11) then
            if (iopt /= 1) write (81,*)1000.0d0*pressu
            caldat(ndata2)=1000.0d0*pressu
        end if

    25 END DO
          
    50 return
!***********************************************************************
!                                                                      *
!     Format part                                                      *
!                                                                      *
!***********************************************************************
    100 format (i3,a60)
    200 format (5i3)
    end subroutine outres
!***********************************************************************
!***********************************************************************

   subroutine gettrainstructs2

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension ctemp(nat,3),dargdc(3,3),dhrdc(3,3),iat(5)
    dimension vibreax(navib*3),vibqc(navib*3),errmatch(3*navib)
    dimension imatch (3*navib)
    character(80) :: qhulp
    character(20) :: qkeyhulp
    character(40) :: qfreqfile
    character(20) :: qrom
    logical :: redundant

    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In gettrainstructs2'
        call timer(65)
        close (65)
    end if
!***********************************************************************
!                                                                      *
!     Output molecular data                                            *
!                                                                      *
!***********************************************************************
    open (63,file='trainset.in',status='old',err=9000)
    10 read (63,'(a80)',err=9010,end=500)qhulp
    qstrana1(1:80)=qhulp
    if (qhulp(1:1) == '#') goto 10
    !if (qhulp(1:6) == 'STRESS') then
    !***********************************************************************
    !                                                                      *
    !     Stress data                                                      *
    !                                                                      *
    !***********************************************************************
    !    11 read (63,'(a80)',err=9010,end=500)qhulp
    !    qstrana1(1:80)=qhulp
    !    if (qhulp(1:9) == 'ENDSTRESS') goto 10
    !    if (qhulp(1:1) == '#') goto 11
    !    istart=1
    !    call stranal(istart,iend,vout,iout,1)
    !    qkeyhulp=qstrana2(1:20)
    !    istart=iend
    !    call stranal(istart,iend,vout,iout,1)
    !    weight=vout
    !    istart=iend
    !        
    !    ihulp=-6
    !    vcomp=-350.0d0
    !    call stranal(istart,iend,vout,iout,1)
    !    if (istart >= iend) goto 12
    !    if (qstrana2 == 'x') ihulp=1
    !    if (qstrana2 == 'y') ihulp=2
    !    if (qstrana2 == 'z') ihulp=3
    !    istart=iend

    !    call stranal(istart,iend,vout,iout,1)
    !    if (istart >= iend) goto 12
    !    vcomp=vout
    !    12 continue

    !    if (qkeyhulp == qkeyw(nprob)) then
    !        ndata2=ndata2+1
    !        if (ihulp == 1) then
    !            caldat(ndata2)=presx/1000.0d0
    !            write (qdatid(ndata2),180)qkeyhulp,iatn
    !        end if
    !        if (ihulp == 2) then
    !            caldat(ndata2)=presy/1000.0d0
    !            write (qdatid(ndata2),181)qkeyhulp,iatn
    !        end if
    !        if (ihulp == 3) then
    !            caldat(ndata2)=presz/1000.0d0
    !            write (qdatid(ndata2),182)qkeyhulp,iatn
    !        end if
    !        compdat(ndata2)=vcomp
    !        weightdat(ndata2)=weight
    !    end if

    !    goto 11
    !end if
       
    if (qhulp(1:6) == 'CHARGE') then
    !***********************************************************************
    !                                                                      *
    !     Charge distribution data                                         *
    !                                                                      *
    !***********************************************************************
        !df398
        redundant = .false.
        15 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:9) == 'ENDCHARGE') goto 10
        if (qhulp(1:1) == '#') goto 15
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:30)
        !df398 add unique structures to trainset_structs array
        do j=1,nmolset
           !write(*,*) 'qkeyhulp = ', qkeyhulp
           if (trainset_structs(j) == qkeyhulp .or. (qkeyhulp(1:1) >= '0' .and. qkeyhulp(1:1) <= '9')) then
              redundant = .true.
              !write(*,*) 'redundant == true or not a structure description'
              exit
           endif
        enddo
        if (redundant .eqv. .false.) then
           junique = junique + 1
           trainset_structs(junique) = qkeyhulp
           !write(*,*) 'CHARGES >>> trainset_struct(',junique,') = ', trainset_structs(junique)
        endif
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        iatn=0
        vcomp=-250.0d0
        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 16
        iatn=iout
        istart=iend

        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 16
        vcomp=vout

        16 continue
              
        !if (qkeyhulp == qkeyw(nprob)) then
        !    if (iatn > 0) then
        !        ndata2=ndata2+1
        !        caldat(ndata2)=ch(iatn)
        !        compdat(ndata2)=chaset(nprob,iatn)
        !        if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),105)qkeyhulp,iatn
        !    else
        !        do i1=1,naset(nprob)
        !            ndata2=ndata2+1
        !            caldat(ndata2)=ch(i1)
        !            compdat(ndata2)=chaset(nprob,i1)
        !            if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !            weightdat(ndata2)=weight
        !            write (qdatid(ndata2),106)qkeyhulp,i1
        !        end do
        !    end if
        !end if

        goto 15
    end if

    if (qhulp(1:8) == 'GEOMETRY') then
    !***********************************************************************
    !                                                                      *
    !     Geometry data                                                    *
    !                                                                      *
    !***********************************************************************
        !df398
        redundant = .false.
        20 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:1) == '#') goto 20
        if (qhulp(1:11) == 'ENDGEOMETRY') goto 10
             
        vcomp=zero
        do i1=1,4
            iat(i1)=0
        end do

        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:20)
        !df398 add unique structures to trainset_structs array
        do j=1,nmolset
           !write(*,*) 'qkeyhulp = ', qkeyhulp
           if (trainset_structs(j) == qkeyhulp .or. (qkeyhulp(1:1) >= '0' .and. qkeyhulp(1:1) <= '9')) then
              redundant = .true.
              !write(*,*) 'redundant == true or not a structure description'
              exit
           endif
        enddo
        if (redundant .eqv. .false.) then
           junique = junique + 1
           trainset_structs(junique) = qkeyhulp
           !write(*,*) 'GEOMETRY >>> trainset_struct(',junique,') = ', trainset_structs(junique)
        endif

        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        do i1=1,5
            call stranal(istart,iend,vout,iout,1)
            if (istart >= iend) then
                iat(i1-1)=0
                goto 21
            end if

            iat(i1)=iout
            vcomp=vout
            istart=iend
        end do

        21 continue

        !if (qkeyhulp == qkeyw(nprob)) then

        !    if (iat(4) /= 0) then
        !        call caltor(iat(1),iat(2),iat(3),iat(4),htcalc)
        !        if (abs(vcomp) > zero) then
        !            htcomp=vcomp
        !        else
        !            do i1=1,na
        !                do i2=1,3
        !                    ctemp(i1,i2)=c(i1,i2)
        !                    c(i1,i2)=cset(nprob,i1,i2)
        !                end do
        !            end do
        !            call caltor(iat(1),iat(2),iat(3),iat(4),htcomp)
        !            do i1=1,na
        !                do i2=1,3
        !                    c(i1,i2)=ctemp(i1,i2)
        !                end do
        !            end do
        !        end if
        !        ndata2=ndata2+1
        !        caldat(ndata2)=htcalc
        !        compdat(ndata2)=htcomp
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),110)qkeyhulp,iat(1),iat(2),iat(3),iat(4)
        !    end if

        !    if (iat(4) == 0 .AND. iat(3) /= 0) then
        !        call calvalres(iat(1),iat(2),iat(3),arg,hr,dhrdc,dargdc)
        !        hvcalc=hr*rdndgr
        !        if (vcomp > zero) then
        !            hvcomp=vcomp
        !        else
        !            do i1=1,na
        !                do i2=1,3
        !                    ctemp(i1,i2)=c(i1,i2)
        !                    c(i1,i2)=cset(nprob,i1,i2)
        !                end do
        !            end do
        !            call calvalres(iat(1),iat(2),iat(3),arg,hr,dhrdc,dargdc)
        !            hvcomp=hr*rdndgr
        !            do i1=1,na
        !                do i2=1,3
        !                    c(i1,i2)=ctemp(i1,i2)
        !                end do
        !            end do
        !        end if
        !        ndata2=ndata2+1
        !        caldat(ndata2)=hvcalc
        !        compdat(ndata2)=hvcomp
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),120)qkeyhulp,iat(1),iat(2),iat(3)
        !    end if

        !    if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) /= 0) then
        !        call dista2(iat(1),iat(2),rcalc,dx1,dy1,dz1)
        !        if (vcomp > zero) then
        !            rcomp=vcomp
        !        else
        !            do i1=1,na
        !                do i2=1,3
        !                    ctemp(i1,i2)=c(i1,i2)
        !                    c(i1,i2)=cset(nprob,i1,i2)
        !                end do
        !            end do
        !            call dista2(iat(1),iat(2),rcomp,dx1,dy1,dz1)
        !            do i1=1,na
        !                do i2=1,3
        !                    c(i1,i2)=ctemp(i1,i2)
        !                end do
        !            end do
        !        end if
        !        ndata2=ndata2+1
        !        caldat(ndata2)=rcalc
        !        compdat(ndata2)=rcomp
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),130)qkeyhulp,iat(1),iat(2)
        !    end if

        !    if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) == 0 &
        !    .and. iat(1) > 0) then
        !        dx=c(iat(1),1)-cset(nprob,iat(1),1)
        !        dy=c(iat(1),2)-cset(nprob,iat(1),2)
        !        dz=c(iat(1),3)-cset(nprob,iat(1),3)
        !        diff=sqrt(dx*dx+dy*dy+dz*dz)
        !        ndata2=ndata2+1
        !        caldat(ndata2)=diff
        !        compdat(ndata2)=zero
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),140)qkeyhulp,iat(1)
        !    end if
        !     
        !    if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) == 0 .and. &
        !    iat(1) == -1) then
        !        crms=0.0
        !        do i1=1,na
        !            dx=c(i1,1)-cset(nprob,i1,1)
        !            dy=c(i1,2)-cset(nprob,i1,2)
        !            dz=c(i1,3)-cset(nprob,i1,3)
        !            diff=sqrt(dx*dx+dy*dy+dz*dz)
        !            crms=crms+diff
        !        end do
        !        crms=(crms/float(na))
        !        ndata2=ndata2+1
        !        caldat(ndata2)=crms
        !        compdat(ndata2)=zero
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),145)qkeyhulp
        !    end if
        !     
        !    if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) == 0 .and. &
        !    iat(1) == 0) then
        !        ndata2=ndata2+1
        !        if (vcomp > zero) then
        !            rmsgcomp=vcomp
        !        else
        !            rmsgcomp=zero
        !        end if

        !        rmsg=0.0d0
        !        nmovh=0
        !        do i2=1,na
        !            do i3=1,3
        !                rmsg=rmsg+imove(i2)*d(i3,i2)*d(i3,i2)
        !                nmovh=nmovh+imove(i2)
        !            end do
        !        end do
        !        rmsg=sqrt(rmsg/float(nmovh))

        !        caldat(ndata2)=rmsg
        !        compdat(ndata2)=rmsgcomp
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),150)qkeyhulp
        !    end if

        !end if

        goto 20
    end if

    if (qhulp(1:6) == 'HEATFO') then
    !***********************************************************************
    !                                                                      *
    !     Heat of formation data                                           *
    !                                                                      *
    !***********************************************************************
        redundant = .false.
        25 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:1) == '#') goto 25
        if (qhulp(1:9) == 'ENDHEATFO') goto 10
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:20)
        !df398 add unique structures to trainset_structs array
        do j=1,nmolset
           !write(*,*) 'qkeyhulp = ', qkeyhulp
           if (trainset_structs(j) == qkeyhulp .or. (qkeyhulp(1:1) >= '0' .and. qkeyhulp(1:1) <= '9')) then
              redundant = .true.
              !write(*,*) 'redundant == true or not a structure description'
              exit
           endif
        enddo
        if (redundant .eqv. .false.) then
           junique = junique + 1
           trainset_structs(junique) = qkeyhulp
           !write(*,*) 'HEATFO >>> trainset_struct(',junique,') = ', trainset_structs(junique)
        endif
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vlithf=vout
        !if (qkeyhulp == qkeyw(nprob)) then
        !    sumin=0.0
        !    do i1=1,na
        !        sumin=sumin+vincr(ia(i1,1))
        !    end do
        !    heatfo=sumin+estrc+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
        !    ndata2=ndata2+1
        !    caldat(ndata2)=heatfo
        !    if (iopt == 1) caldat(ndata2)=estrmin
        !    compdat(ndata2)=vlithf
        !    weightdat(ndata2)=weight
        !    write (qdatid(ndata2),160)qkeyhulp
        !    iheatf=iheatf+1
        !    iheada(iheatf)=ndata2
        !    iheada2(iheatf)=nprob
        !end if
        goto 25
    end if

    if (qhulp(1:15) == 'CELL PARAMETERS') then
    !***********************************************************************
    !                                                                      *
    !     Cell parameter data                                              *
    !                                                                      *
    !***********************************************************************
        !df398
        redundant = .false.
        30 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:1) == '#') goto 30
        if (qhulp(1:18) == 'ENDCELL PARAMETERS') goto 10
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:20)
        !df398 add unique structures to trainset_structs array
        do j=1,nmolset
           !write(*,*) 'qkeyhulp = ', qkeyhulp
           if (trainset_structs(j) == qkeyhulp .or. (qkeyhulp(1:1) >= '0' .and. qkeyhulp(1:1) <= '9')) then
              redundant = .true.
              !write(*,*) 'redundant == true or not a structure description'
              exit
           endif
        enddo
        if (redundant .eqv. .false.) then
           junique = junique + 1
           trainset_structs(junique) = qkeyhulp
           !write(*,*) 'CELL PARAMETERS >>> trainset_struct(',junique,') = ', trainset_structs(junique)
        endif
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        ihulp=-6
        vcomp=-350.0d0
        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 31
        if (qstrana2 == 'a') ihulp=1
        if (qstrana2 == 'b') ihulp=2
        if (qstrana2 == 'c') ihulp=3
        if (qstrana2 == 'alpha') ihulp=4
        if (qstrana2 == 'alfa') ihulp=4
        if (qstrana2 == 'beta') ihulp=5
        if (qstrana2 == 'gamma') ihulp=6
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 31
        vcomp=vout
        31 continue
        !if (qkeyhulp == qkeyw(nprob)) then
        !    if (ihulp == 1) then
        !        ndata2=ndata2+1
        !        compdat(ndata2)=axisset(nprob,1)
        !        if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !        caldat(ndata2)=axiss(1)
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),170)qkeyhulp
        !    end if
        !          
        !    if (ihulp == 2) then
        !        ndata2=ndata2+1
        !        compdat(ndata2)=axisset(nprob,2)
        !        if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !        caldat(ndata2)=axiss(2)
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),171)qkeyhulp
        !    end if
        !          
        !    if (ihulp == 3) then
        !        ndata2=ndata2+1
        !        compdat(ndata2)=axisset(nprob,3)
        !        if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !        caldat(ndata2)=axiss(3)
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),172)qkeyhulp
        !    end if
        !          
        !    if (ihulp == 4) then
        !        ndata2=ndata2+1
        !        compdat(ndata2)=anglesset(nprob,1)
        !        if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !        caldat(ndata2)=angles(1)
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),175)qkeyhulp
        !    end if
        !          
        !    if (ihulp == 5) then
        !        ndata2=ndata2+1
        !        compdat(ndata2)=anglesset(nprob,2)
        !        if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !        caldat(ndata2)=angles(2)
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),176)qkeyhulp
        !    end if
        !          
        !    if (ihulp == 6) then
        !        ndata2=ndata2+1
        !        compdat(ndata2)=anglesset(nprob,3)
        !        if (vcomp > -200.0d0) compdat(ndata2)=vcomp
        !        caldat(ndata2)=angles(3)
        !        weightdat(ndata2)=weight
        !        write (qdatid(ndata2),177)qkeyhulp
        !    end if

        !end if
              
        goto 30
    end if

    !if (qhulp(1:11) == 'FREQUENCIES') then
    !    40 read (63,'(a80)',err=9010,end=500)qhulp
    !    qstrana1(1:80)=qhulp
    !    if (qhulp(1:1) == '#') goto 40
    !    if (qhulp(1:14) == 'ENDFREQUENCIES') goto 10
    !    istart=1
    !    call stranal(istart,iend,vout,iout,1)
    !    qkeyhulp=qstrana2(1:20)
    !    istart=iend
    !    call stranal(istart,iend,vout,iout,1)
    !    weightf=vout
    !    istart=iend
    !    call stranal(istart,iend,vout,iout,1)
    !    weightm=vout
    !    istart=iend
    !    call stranal(istart,iend,vout,iout,1)
    !    qfreqfile=qstrana2
    !    istart=iend

    !    if (qkeyhulp == qkeyw(nprob)) then
    !        itrain=1
    !        call vibra(qfreqfile,vibreax,vibqc,imatch,errmatch,itrain,klinear)
    !        do i1=1,na*3-6+klinear
    !            ndata2=ndata2+1
    !            caldat(ndata2)=vibreax(i1+6-klinear)
    !            weightdat(ndata2)=weightf
    !            compdat(ndata2)=vibqc(imatch(i1))
    !            write (qdatid(ndata2),178)qkeyhulp,i1+6
    !            ndata2=ndata2+1
    !            caldat(ndata2)=1000.0d0*errmatch(i1)
    !            weightdat(ndata2)=weightm
    !            compdat(ndata2)=zero
    !            write (qdatid(ndata2),179)qkeyhulp,i1+6
    !        end do
    !    end if

    !    goto 40
    !end if
    goto 10
    500 continue
    close (63)
    return
    105 format (a20,' Charge atom:',i4)
    106 format (a20,' Charge atom:',i4)
    110 format (a20,' Torsion angle:',4i4)
    120 format (a20,' Valence angle:',3i4)
    130 format (a20,' Bond distance:',2i4)
    140 format (a20,' Position atom:',i4)
    145 format (a20,' CRMS:',i4)
    150 format (a20,' RMSG force:')
    160 format (a20,' Heat of formation:')
    170 format (a20,' a:')
    171 format (a20,' b:')
    172 format (a20,' c:')
    175 format (a20,' alfa:')
    176 format (a20,' beta:')
    177 format (a20,' gamma:')
    178 format (a20,' Frequency:',i4)
    179 format (a20,' Error in mode:',i4)
    180 format (a20,' pressure x:')
    181 format (a20,' pressure y:')
    182 format (a20,' pressure z:')
    9000 continue
    if (iopt == 1) stop 'Cannot open trainset.in'
!     write (*,*) 'Cannot open trainset.in'
    return
    9010 stop 'Error reading trainset.in'

    end subroutine gettrainstructs2
!***********************************************************************
!***********************************************************************

   subroutine outres2

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension ctemp(nat,3),dargdc(3,3),dhrdc(3,3),iat(5)
    dimension vibreax(navib*3),vibqc(navib*3),errmatch(3*navib)
    dimension imatch (3*navib)
    character(80) :: qhulp
    character(40) :: qkeyhulp
    character(40) :: qfreqfile
    character(40) :: qrom
    logical :: redundant

    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In outres2'
        call timer(65)
        close (65)
    end if
!***********************************************************************
!                                                                      *
!     Output molecular data                                            *
!                                                                      *
!***********************************************************************
    open (63,file='trainset.in',status='old',err=9000)
    10 read (63,'(a80)',err=9010,end=500)qhulp
    qstrana1(1:80)=qhulp
    if (qhulp(1:1) == '#') goto 10
    if (qhulp(1:6) == 'STRESS') then
    !***********************************************************************
    !                                                                      *
    !     Stress data                                                      *
    !                                                                      *
    !***********************************************************************
        11 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:9) == 'ENDSTRESS') goto 10
        if (qhulp(1:1) == '#') goto 11
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:30)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
            
        ihulp=-6
        vcomp=-350.0d0
        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 12
        if (qstrana2 == 'x') ihulp=1
        if (qstrana2 == 'y') ihulp=2
        if (qstrana2 == 'z') ihulp=3
        istart=iend

        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 12
        vcomp=vout
        12 continue

        if (qkeyhulp == qkeyw(nprob)) then
            ndata2=ndata2+1
            if (ihulp == 1) then
                caldat(ndata2)=presx/1000.0d0
                write (qdatid(ndata2),180)qkeyhulp,iatn
            end if
            if (ihulp == 2) then
                caldat(ndata2)=presy/1000.0d0
                write (qdatid(ndata2),181)qkeyhulp,iatn
            end if
            if (ihulp == 3) then
                caldat(ndata2)=presz/1000.0d0
                write (qdatid(ndata2),182)qkeyhulp,iatn
            end if
            compdat(ndata2)=vcomp
            weightdat(ndata2)=weight
        end if

        goto 11
    end if
       
    if (qhulp(1:6) == 'CHARGE') then
    !***********************************************************************
    !                                                                      *
    !     Charge distribution data                                         *
    !                                                                      *
    !***********************************************************************
        15 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:9) == 'ENDCHARGE') goto 10
        if (qhulp(1:1) == '#') goto 15
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:30)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        iatn=0
        vcomp=-250.0d0
        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 16
        iatn=iout
        istart=iend

        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 16
        vcomp=vout

        16 continue
              
        if (qkeyhulp == qkeyw(nprob)) then
            if (iatn > 0) then
                ndata2=ndata2+1
                caldat(ndata2)=ch(iatn)
                compdat(ndata2)=chaset(nprob,iatn)
                if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                weightdat(ndata2)=weight
                write (qdatid(ndata2),105)qkeyhulp,iatn
            else
                do i1=1,naset(nprob)
                    ndata2=ndata2+1
                    caldat(ndata2)=ch(i1)
                    compdat(ndata2)=chaset(nprob,i1)
                    if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                    weightdat(ndata2)=weight
                    write (qdatid(ndata2),106)qkeyhulp,i1
                end do
            end if
        end if

        goto 15
    end if

    if (qhulp(1:8) == 'GEOMETRY') then
    !***********************************************************************
    !                                                                      *
    !     Geometry data                                                    *
    !                                                                      *
    !***********************************************************************
        20 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:1) == '#') goto 20
        if (qhulp(1:11) == 'ENDGEOMETRY') goto 10
             
        vcomp=zero
        do i1=1,4
            iat(i1)=0
        end do

        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:30)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        do i1=1,5
            call stranal(istart,iend,vout,iout,1)
            if (istart >= iend) then
                iat(i1-1)=0
                goto 21
            end if

            iat(i1)=iout
            vcomp=vout
            istart=iend
        end do

        21 continue

        if (qkeyhulp == qkeyw(nprob)) then

            if (iat(4) /= 0) then
                call caltor(iat(1),iat(2),iat(3),iat(4),htcalc)
                if (abs(vcomp) > zero) then
                    htcomp=vcomp
                else
                    do i1=1,na
                        do i2=1,3
                            ctemp(i1,i2)=c(i1,i2)
                            c(i1,i2)=cset(nprob,i1,i2)
                        end do
                    end do
                    call caltor(iat(1),iat(2),iat(3),iat(4),htcomp)
                    do i1=1,na
                        do i2=1,3
                            c(i1,i2)=ctemp(i1,i2)
                        end do
                    end do
                end if
                ndata2=ndata2+1
                caldat(ndata2)=htcalc
                compdat(ndata2)=htcomp
                weightdat(ndata2)=weight
                write (qdatid(ndata2),110)qkeyhulp,iat(1),iat(2),iat(3),iat(4)
            end if

            if (iat(4) == 0 .AND. iat(3) /= 0) then
                call calvalres(iat(1),iat(2),iat(3),arg,hr,dhrdc,dargdc)
                hvcalc=hr*rdndgr
                if (vcomp > zero) then
                    hvcomp=vcomp
                else
                    do i1=1,na
                        do i2=1,3
                            ctemp(i1,i2)=c(i1,i2)
                            c(i1,i2)=cset(nprob,i1,i2)
                        end do
                    end do
                    call calvalres(iat(1),iat(2),iat(3),arg,hr,dhrdc,dargdc)
                    hvcomp=hr*rdndgr
                    do i1=1,na
                        do i2=1,3
                            c(i1,i2)=ctemp(i1,i2)
                        end do
                    end do
                end if
                ndata2=ndata2+1
                caldat(ndata2)=hvcalc
                compdat(ndata2)=hvcomp
                weightdat(ndata2)=weight
                write (qdatid(ndata2),120)qkeyhulp,iat(1),iat(2),iat(3)
            end if

            if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) /= 0) then
                call dista2(iat(1),iat(2),rcalc,dx1,dy1,dz1)
                if (vcomp > zero) then
                    rcomp=vcomp
                else
                    do i1=1,na
                        do i2=1,3
                            ctemp(i1,i2)=c(i1,i2)
                            c(i1,i2)=cset(nprob,i1,i2)
                        end do
                    end do
                    call dista2(iat(1),iat(2),rcomp,dx1,dy1,dz1)
                    do i1=1,na
                        do i2=1,3
                            c(i1,i2)=ctemp(i1,i2)
                        end do
                    end do
                end if
                ndata2=ndata2+1
                caldat(ndata2)=rcalc
                compdat(ndata2)=rcomp
                weightdat(ndata2)=weight
                write (qdatid(ndata2),130)qkeyhulp,iat(1),iat(2)
            end if

            if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) == 0 &
            .and. iat(1) > 0) then
                dx=c(iat(1),1)-cset(nprob,iat(1),1)
                dy=c(iat(1),2)-cset(nprob,iat(1),2)
                dz=c(iat(1),3)-cset(nprob,iat(1),3)
                diff=sqrt(dx*dx+dy*dy+dz*dz)
                ndata2=ndata2+1
                caldat(ndata2)=diff
                compdat(ndata2)=zero
                weightdat(ndata2)=weight
                write (qdatid(ndata2),140)qkeyhulp,iat(1)
            end if
             
            if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) == 0 .and. &
            iat(1) == -1) then
                crms=0.0
                do i1=1,na
                    dx=c(i1,1)-cset(nprob,i1,1)
                    dy=c(i1,2)-cset(nprob,i1,2)
                    dz=c(i1,3)-cset(nprob,i1,3)
                    diff=sqrt(dx*dx+dy*dy+dz*dz)
                    crms=crms+diff
                end do
                crms=(crms/float(na))
                ndata2=ndata2+1
                caldat(ndata2)=crms
                compdat(ndata2)=zero
                weightdat(ndata2)=weight
                write (qdatid(ndata2),145)qkeyhulp
            end if
             
            if (iat(4) == 0 .AND. iat(3) == 0 .AND. iat(2) == 0 .and. &
            iat(1) == 0) then
                ndata2=ndata2+1
                if (vcomp > zero) then
                    rmsgcomp=vcomp
                else
                    rmsgcomp=zero
                end if

                rmsg=0.0d0
                nmovh=0
                do i2=1,na
                    do i3=1,3
                        rmsg=rmsg+imove(i2)*d(i3,i2)*d(i3,i2)
                        nmovh=nmovh+imove(i2)
                    end do
                end do
                rmsg=sqrt(rmsg/float(nmovh))

                caldat(ndata2)=rmsg
                compdat(ndata2)=rmsgcomp
                weightdat(ndata2)=weight
                write (qdatid(ndata2),150)qkeyhulp
            end if

        end if

        goto 20
    end if

    if (qhulp(1:6) == 'HEATFO') then
    !***********************************************************************
    !                                                                      *
    !     Heat of formation data                                           *
    !                                                                      *
    !***********************************************************************
        25 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:1) == '#') goto 25
        if (qhulp(1:9) == 'ENDHEATFO') goto 10
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:30)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vlithf=vout
        if (qkeyhulp == qkeyw(nprob)) then
            sumin=0.0
            do i1=1,na
                sumin=sumin+vincr(ia(i1,1))
            end do
            heatfo=sumin+estrc+4.0d0*rgasc*xjouca*0.29815d0*nmolo5
            ndata2=ndata2+1
            caldat(ndata2)=heatfo
            if (iopt == 1) caldat(ndata2)=estrmin
            compdat(ndata2)=vlithf
            weightdat(ndata2)=weight
            write (qdatid(ndata2),160)qkeyhulp
            iheatf=iheatf+1
            iheada(iheatf)=ndata2
            iheada2(iheatf)=nprob
        end if
        goto 25
    end if

    if (qhulp(1:15) == 'CELL PARAMETERS') then
    !***********************************************************************
    !                                                                      *
    !     Cell parameter data                                              *
    !                                                                      *
    !***********************************************************************
        30 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:1) == '#') goto 30
        if (qhulp(1:18) == 'ENDCELL PARAMETERS') goto 10
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:30)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weight=vout
        istart=iend
        ihulp=-6
        vcomp=-350.0d0
        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 31
        if (qstrana2 == 'a') ihulp=1
        if (qstrana2 == 'b') ihulp=2
        if (qstrana2 == 'c') ihulp=3
        if (qstrana2 == 'alpha') ihulp=4
        if (qstrana2 == 'alfa') ihulp=4
        if (qstrana2 == 'beta') ihulp=5
        if (qstrana2 == 'gamma') ihulp=6
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        if (istart >= iend) goto 31
        vcomp=vout
        31 continue
        if (qkeyhulp == qkeyw(nprob)) then
            if (ihulp == 1) then
                ndata2=ndata2+1
                compdat(ndata2)=axisset(nprob,1)
                if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                caldat(ndata2)=axiss(1)
                weightdat(ndata2)=weight
                write (qdatid(ndata2),170)qkeyhulp
            end if
                  
            if (ihulp == 2) then
                ndata2=ndata2+1
                compdat(ndata2)=axisset(nprob,2)
                if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                caldat(ndata2)=axiss(2)
                weightdat(ndata2)=weight
                write (qdatid(ndata2),171)qkeyhulp
            end if
                  
            if (ihulp == 3) then
                ndata2=ndata2+1
                compdat(ndata2)=axisset(nprob,3)
                if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                caldat(ndata2)=axiss(3)
                weightdat(ndata2)=weight
                write (qdatid(ndata2),172)qkeyhulp
            end if
                  
            if (ihulp == 4) then
                ndata2=ndata2+1
                compdat(ndata2)=anglesset(nprob,1)
                if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                caldat(ndata2)=angles(1)
                weightdat(ndata2)=weight
                write (qdatid(ndata2),175)qkeyhulp
            end if
                  
            if (ihulp == 5) then
                ndata2=ndata2+1
                compdat(ndata2)=anglesset(nprob,2)
                if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                caldat(ndata2)=angles(2)
                weightdat(ndata2)=weight
                write (qdatid(ndata2),176)qkeyhulp
            end if
                  
            if (ihulp == 6) then
                ndata2=ndata2+1
                compdat(ndata2)=anglesset(nprob,3)
                if (vcomp > -200.0d0) compdat(ndata2)=vcomp
                caldat(ndata2)=angles(3)
                weightdat(ndata2)=weight
                write (qdatid(ndata2),177)qkeyhulp
            end if

        end if
              
        goto 30
    end if

    if (qhulp(1:11) == 'FREQUENCIES') then
        40 read (63,'(a80)',err=9010,end=500)qhulp
        qstrana1(1:80)=qhulp
        if (qhulp(1:1) == '#') goto 40
        if (qhulp(1:14) == 'ENDFREQUENCIES') goto 10
        istart=1
        call stranal(istart,iend,vout,iout,1)
        qkeyhulp=qstrana2(1:30)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weightf=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        weightm=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        qfreqfile=qstrana2
        istart=iend

        if (qkeyhulp == qkeyw(nprob)) then
            itrain=1
            call vibra(qfreqfile,vibreax,vibqc,imatch,errmatch,itrain,klinear)
            do i1=1,na*3-6+klinear
                ndata2=ndata2+1
                caldat(ndata2)=vibreax(i1+6-klinear)
                weightdat(ndata2)=weightf
                compdat(ndata2)=vibqc(imatch(i1))
                write (qdatid(ndata2),178)qkeyhulp,i1+6
                ndata2=ndata2+1
                caldat(ndata2)=1000.0d0*errmatch(i1)
                weightdat(ndata2)=weightm
                compdat(ndata2)=zero
                write (qdatid(ndata2),179)qkeyhulp,i1+6
            end do
        end if

        goto 40
    end if
    goto 10
    500 continue
    close (63)
    return
    105 format (a20,' Charge atom:',i4)
    106 format (a20,' Charge atom:',i4)
    110 format (a20,' Torsion angle:',4i4)
    120 format (a20,' Valence angle:',3i4)
    130 format (a20,' Bond distance:',2i4)
    140 format (a20,' Position atom:',i4)
    145 format (a20,' CRMS:',i4)
    150 format (a20,' RMSG force:')
    160 format (a20,' Heat of formation:')
    170 format (a20,' a:')
    171 format (a20,' b:')
    172 format (a20,' c:')
    175 format (a20,' alfa:')
    176 format (a20,' beta:')
    177 format (a20,' gamma:')
    178 format (a20,' Frequency:',i4)
    179 format (a20,' Error in mode:',i4)
    180 format (a20,' pressure x:')
    181 format (a20,' pressure y:')
    182 format (a20,' pressure z:')
    9000 continue
    if (iopt == 1) stop 'Cannot open trainset.in'
!     write (*,*) 'Cannot open trainset.in'
    return
    9010 stop 'Error reading trainset.in'
    end subroutine outres2
!***********************************************************************
!***********************************************************************

    subroutine gettrainstructs

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension ctemp(nat,3),dargdc(3,3),dhrdc(3,3),vdiv(5)
    dimension iat(2,5),ndiv(5),ist2(5),iend2(5),qkeyh(5), &
    qsign(5),qdiv(5)
    character(100) :: qhulp
    character(30) :: qhulp2,qhulp3
    character(20) :: qkeyh
    character(5) :: qdiv
    character(1) :: qsign
    logical :: redundant
!***********************************************************************
!                                                                      *
!     Output molecular energy data                                     *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In gettrainstructs'
        call timer(65)
        close (65)
    end if

    iline=0
    !df398
    open (63,file='trainset.in',status='old',err=9000)
    10 read (63,'(a100)',err=9010,end=500)qhulp
    iline=iline+1
    if (qhulp(1:1) == '#') goto 10
    if (qhulp(1:6) == 'ENERGY') then
        15 read (63,'(a100)',err=9010,end=500)qhulp
        qstrana1(1:100)=qhulp
        iline=iline+1
        if (qhulp(1:1) == '#') goto 15
        if (qhulp(1:9) == 'ENDENERGY') goto 10
        vcomp=zero
        istart=1
        sumcalc=zero
        sumcomp=zero

        call stranal(istart,iend,vout,iout,2)
        weight=vout
        istart=iend

        nhulp=0
        do i1=1,5
            qsign(i1)=' '
            qkeyh(i1)=' '
            qdiv(i1)=' '
            ist2(i1)=iend
            iend2(i1)=iend
        end do
        do i1=1,5
            call stranal(istart,iend,vout,iout,2)
            qsign(i1)=qstrana2(1:1)
            if (qsign(i1) == '+' .OR. qsign(i1) == '-') &
            istart=istart+1
            call stranal(istart,iend,vout,iout,2)
            ist2(i1)=istart
            iend2(i1)=iend-1
            qkeyh(i1)=qstrana2(1:30)
            !df398 add unique structures to trainset_structs
            redundant = .false.
            do j=1,nmolset
               if (trainset_structs(j) == qkeyh(i1) .or. (qkeyh(i1)(1:1) >= '0' .and. qkeyh(i1)(1:1) <= '9')) then
                  redundant = .true.
                  !write(*,*) 'redundant == true or not a structure description'
                  exit
               endif
            enddo
            if (redundant .eqv. .false.) then
               junique = junique + 1
               !write(*,*) 'SIZE of trainset_structs = ', size(trainset_structs)
               trainset_structs(junique) = qkeyh(i1)
               !write(*,*) 'ENERGY>>> trainset_struct(',junique,') = ', trainset_structs(junique)
            endif
            istart=iend
            call stranal(istart,iend,vout,iout,2)
            vdiv(i1)=vout
            if (abs(vdiv(i1)) > 0.0001d0) then
                nhulp=nhulp+1
                istart=iend
                istart2=iend
                !if (vdiv(i1) < 10.0d0) write (qdiv(nhulp), &
                !'(a1,f4.2,a1)')'/',vdiv(i1),' '
                if (vdiv(i1) < 10.0) write (qdiv(nhulp), '(a1,i3,a1)')'/',int(vdiv(i1)),' '

                if (qsign(i1) /= '+' .AND. qsign(i1) /= '-') then
                    qsign(i1)='+'
                    write (*,*)'Warning: forgot operator symbol on line:' &
                    ,iline,';assume +'
                end if

            else
                qsign(i1)=' '
                qkeyh(i1)=' '
                ist2(i1)=iend
                iend2(i1)=iend
                goto 20
            end if
        end do

        20 continue
        istart=istart2
        call stranal(istart,iend,vout,iout,2)
        vcomp=vout
        istart=iend

        !do i1=1,nhulp
        !    isign=1
        !    if (qsign(i1) == '-') isign=-1
        !    ifound=0
        !    do i2=1,nprob
        !        if (qkeyh(i1) == qkeyw(i2)) then
        !            sumcalc=sumcalc+isign*enmolend(i2)/vdiv(i1)
        !            sumcomp=sumcomp+isign*enmolset(i2)/vdiv(i1)
        !            ifound=1
        !        end if
        !    end do
        !    if (ifound == 0) then
        !        if (nsurp == 0) then
        !            write (*,*)'Could not find data linked to keyword:', &
        !            qkeyh(i1),';ignored'
        !        end if
        !        goto 15
        !    end if
        !end do
        !ndata2=ndata2+1
        !caldat(ndata2)=sumcalc
        !compdat(ndata2)=sumcomp
        !if (vcomp /= zero) compdat(ndata2)=vcomp
        !weightdat(ndata2)=weight
        !qdatid(ndata2)='Energy '//qsign(1)// &
        !qhulp(ist2(1):iend2(1))//qdiv(1)//qsign(2)// &
        !qhulp(ist2(2):iend2(2))//qdiv(2)//qsign(3)// &
        !qhulp(ist2(3):iend2(3))//qdiv(3)//qsign(4)// &
        !qhulp(ist2(4):iend2(4))//qdiv(4)//qsign(5)// &
        !qhulp(ist2(5):iend2(5))//qdiv(5)


        goto 15
    end if

    goto 10
    500 continue
    close (63)
    return
    9000 continue
    if (iopt == 1) stop 'Cannot open trainset.in'
    return
    9010 stop 'Error reading trainset.in'

    end subroutine gettrainstructs
!***********************************************************************
!***********************************************************************

    subroutine outresend

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    dimension ctemp(nat,3),dargdc(3,3),dhrdc(3,3),vdiv(5)
    dimension iat(2,5),ndiv(5),ist2(5),iend2(5),qkeyh(5), &
    qsign(5),qdiv(5)
    character(100) :: qhulp
    character(30) :: qhulp2,qhulp3
    character(30) :: qkeyh
    character(5) :: qdiv
    character(1) :: qsign
    logical :: redundant
!***********************************************************************
!                                                                      *
!     Output molecular energy data                                     *
!                                                                      *
!***********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In outresend'
        call timer(65)
        close (65)
    end if

    iline=0
    open (63,file='trainset.in',status='old',err=9000)
    10 read (63,'(a100)',err=9010,end=500)qhulp
    iline=iline+1
    if (qhulp(1:1) == '#') goto 10
    if (qhulp(1:6) == 'ENERGY') then
        15 read (63,'(a100)',err=9010,end=500)qhulp
        qstrana1(1:100)=qhulp
        iline=iline+1
        if (qhulp(1:1) == '#') goto 15
        if (qhulp(1:9) == 'ENDENERGY') goto 10
        vcomp=zero
        istart=1
        sumcalc=zero
        sumcomp=zero

        call stranal(istart,iend,vout,iout,2)
        weight=vout
        istart=iend

        nhulp=0
        do i1=1,5
            qsign(i1)=' '
            qkeyh(i1)=' '
            qdiv(i1)=' '
            ist2(i1)=iend
            iend2(i1)=iend
        end do
        do i1=1,5
            call stranal(istart,iend,vout,iout,2)
            qsign(i1)=qstrana2(1:1)
            if (qsign(i1) == '+' .OR. qsign(i1) == '-') &
            istart=istart+1
            call stranal(istart,iend,vout,iout,2)
            ist2(i1)=istart
            iend2(i1)=iend-1
            qkeyh(i1)=qstrana2(1:30)
            istart=iend
            call stranal(istart,iend,vout,iout,2)
            vdiv(i1)=vout
            if (abs(vdiv(i1)) > 0.0001d0) then
                nhulp=nhulp+1
                istart=iend
                istart2=iend
                !if (vdiv(i1) < 10.0d0) write (qdiv(nhulp), &
                !'(a1,f4.2,a1)')'/',vdiv(i1),' '
                if (vdiv(i1) < 10.0) write (qdiv(nhulp), '(a1,i3,a1)')'/',int(vdiv(i1)),' '

                if (qsign(i1) /= '+' .AND. qsign(i1) /= '-') then
                    qsign(i1)='+'
                    write (*,*)'Warning: forgot operator symbol on line:' &
                    ,iline,';assume +'
                end if

            else
                qsign(i1)=' '
                qkeyh(i1)=' '
                ist2(i1)=iend
                iend2(i1)=iend
                goto 20
            end if
        end do

        20 continue
        istart=istart2
        call stranal(istart,iend,vout,iout,2)
        vcomp=vout
        istart=iend

        do i1=1,nhulp
            isign=1
            if (qsign(i1) == '-') isign=-1
            ifound=0
            do i2=1,nprob
                if (qkeyh(i1) == qkeyw(i2)) then
                    sumcalc=sumcalc+isign*enmolend(i2)/vdiv(i1)
                    sumcomp=sumcomp+isign*enmolset(i2)/vdiv(i1)
                    ifound=1
                end if
            end do
            if (ifound == 0) then
                if (nsurp == 0) then
                    write (*,*)'Could not find data linked to keyword:', &
                    qkeyh(i1),';ignored'
                end if
                goto 15
            end if
        end do
        ndata2=ndata2+1
        caldat(ndata2)=sumcalc
        compdat(ndata2)=sumcomp
        if (vcomp /= zero) compdat(ndata2)=vcomp
        weightdat(ndata2)=weight
        qdatid(ndata2)='Energy '//qsign(1)// &
        qhulp(ist2(1):iend2(1))//qdiv(1)//qsign(2)// &
        qhulp(ist2(2):iend2(2))//qdiv(2)//qsign(3)// &
        qhulp(ist2(3):iend2(3))//qdiv(3)//qsign(4)// &
        qhulp(ist2(4):iend2(4))//qdiv(4)//qsign(5)// &
        qhulp(ist2(5):iend2(5))//qdiv(5)


        goto 15
    end if

    if (qhulp(1:7) == 'GEODIFF') then
        30 read (63,'(a100)',err=9010,end=500)qhulp
        qstrana1(1:100)=qhulp
        iline=iline+1
        if (qhulp(1:1) == '#') goto 30
        if (qhulp(1:10) == 'ENDGEODIFF') goto 10
        vcomp=zero
        istart=1
        nhulp=0
        do i1=1,2
            qsign(i1)=' '
            qkeyh(i1)=' '
            ist2(i1)=iend
            iend2(i1)=iend
        end do
        do i1=1,4
            iat(1,i1)=0
            iat(2,i1)=0
        end do

        sumcalc=zero
        sumcomp=zero
        call stranal(istart,iend,vout,iout,2)
        weight=vout
        istart=iend
        call stranal(istart,iend,vout,iout,2)
        qsign(1)=qstrana2(1:1)
        if (qsign(1) == '+' .OR. qsign(1) == '-') then
            istart=istart+1
        else
            write (*,*)'Warning: forgot first operator symbol on line:' &
            ,iline,';assume +'
            qsign(1)='+'
        end if
        call stranal(istart,iend,vout,iout,2)
        qkeyh(1)=qstrana2(1:30)
        ist2(1)=istart
        iend2(1)=iend-1
        istart=iend

        do i1=1,5
            call stranal(istart,iend,vout,iout,2)
            if (iout == 0) goto 31
            iat(1,i1)=iout
            istart=iend
        end do
        31 continue

        call stranal(istart,iend,vout,iout,2)
        qsign(2)=qstrana2(1:1)
        if (qsign(2) == '+' .OR. qsign(2) == '-') then
            istart=istart+1
        else
            write (*,*)'Warning: forgot second operator symbol on line:' &
            ,iline,';assume -'
            qsign(2)='-'
        end if
        call stranal(istart,iend,vout,iout,2)
        qkeyh(2)=qstrana2(1:30)
        ist2(2)=istart
        iend2(2)=iend-1
        istart=iend

        do i1=1,5
            call stranal(istart,iend,vout,iout,2)
            if (istart >= iend) then
                iat(2,i1-1)=0
                goto 32
            end if
            iat(2,i1)=iout
            vcomp=vout
            istart=iend
        end do

        32 continue

        do i1=1,2
            isign=1
            if (qsign(i1) == '-') isign=-1
            ifound=0
            do i2=1,nprob
                if (qkeyh(i1) == qkeyw(i2)) then
                    ifound=1

                    do i3=1,nat
                        do i4=1,3
                            ctemp(i3,i4)=c(i3,i4)
                            c(i3,i4)=cset(i2,i3,i4)
                        end do
                    end do

                    if (iat(i1,4) /= 0) then    !Torsion angle
                        call caltor(iat(i1,1),iat(i1,2),iat(i1,3),iat(i1,4),htcalc)
                        sumcalc=sumcalc+isign*htcalc
                        if (i1 == 1) write (qhulp2,100)qsign(i1),qkeyh(i1), &
                        iat(i1,1),iat(i1,2),iat(i1,3),iat(i1,4)
                        if (i1 == 2) write (qhulp3,100)qsign(i1),qkeyh(i1), &
                        iat(i1,1),iat(i1,2),iat(i1,3),iat(i1,4)
                        100 format(a1,a9,'Tor:',4i3)
                    end if

                    if (iat(i1,4) == 0 .AND. iat(i1,3) /= 0) then !Valence angle
                        call calvalres(iat(i1,1),iat(i1,2),iat(i1,3),arg,hr, &
                        dhrdc,dargdc)
                        hvcalc=hr*rdndgr
                        sumcalc=sumcalc+isign*hvcalc
                        if (i1 == 1) write (qhulp2,120)qsign(i1),qkeyh(i1), &
                        iat(i1,1),iat(i1,2),iat(i1,3)
                        if (i1 == 2) write (qhulp3,120)qsign(i1),qkeyh(i1), &
                        iat(i1,1),iat(i1,2),iat(i1,3)
                        120 format(a1,a9,'Angl:',3i3)
                    end if

                    if (iat(i1,4) == 0 .AND. iat(i1,3) == 0 .AND. iat(i1,2) /= 0) then !Bond
                        call dista2(iat(i1,1),iat(i1,2),rcalc,dx1,dy1,dz1)
                        sumcalc=sumcalc+isign*rcalc
                        if (i1 == 1) write (qhulp2,140)qsign(i1),qkeyh(i1), &
                        iat(i1,1),iat(i1,2)
                        if (i1 == 2) write (qhulp3,140)qsign(i1),qkeyh(i1), &
                        iat(i1,1),iat(i1,2)
                        140 format(a1,a9,'Bond:',2i3)
                    end if

                    do i3=1,na
                        do i4=1,3
                            c(i3,i4)=ctemp(i3,i4)
                        end do
                    end do

                end if
            end do
            if (ifound == 0) then
                write (*,*)'Could not find data linked to keyword:', &
                qkeyh(i1),';ignored'
                goto 30
            end if
        end do

        ndata2=ndata2+1
        caldat(ndata2)=sumcalc
        compdat(ndata2)=vcomp
        weightdat(ndata2)=weight
        qdatid(ndata2)='Geo '//qhulp2//qhulp3
        goto 30
    end if

    goto 10
    500 continue
    close (63)
    return
    9000 continue
    if (iopt == 1) stop 'Cannot open trainset.in'
    return
    9010 stop 'Error reading trainset.in'
    end subroutine outresend
!***********************************************************************
!***********************************************************************

    subroutine stranal(istart,iend,vout,iout,icheck)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(1) :: qchar
    dimension qchar(5)
!*********************************************************************
!                                                                    *
!     Analyze string for special characters; find words in string    *
!                                                                    *
!*********************************************************************
    qchar(1)=' '
    qchar(2)='/'
          
    ifound1=0
    do i1=istart,200
        ifound2=0
        do i2=1,icheck

            if (qstrana1(i1:i1) == qchar(i2)) then
                ifound2=1
                if (ifound1 == 1) then     !End of word
                    iend=i1
                    goto 10
                end if

            end if

        end do

        if (ifound2 == 0 .AND. ifound1 == 0) then     !Start of word
            istart2=i1
            ifound1=1
        end if

    end do

    10 continue
    qstrana2=' '
    vout=zero
    iout=0

    if (ifound1 == 1) then
        qstrana2=qstrana1(istart2:iend-1)
        istart=istart2
        vout=zero
        read (qstrana2,*,end=20,err=20) vout
        20 iout=int(vout)
    end if

    return
    end subroutine stranal
!***********************************************************************
!***********************************************************************

    subroutine sortmol(ias,na2,nmohulp)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    logical :: found
    dimension ias(nat,mbond+3)
!*********************************************************************
!                                                                    *
!     Sort molecules                                                 *
!                                                                    *
!*********************************************************************
    FOUND= .FALSE. 
    DO 31 K1=1,NA2
        IF (ias(K1,3+mbond) == 0) FOUND= .TRUE. 
        IF (ias(K1,3+mbond) > nmohulp) nmohulp=ias(k1,3+mbond)
    31 END DO
    IF ( .NOT. FOUND) GOTO 32
!***********************************************************************
!                                                                      *
!     Molecule numbers are assigned. No restrictions are made for the  *
!     sequence of the numbers in the connection table.                 *
!                                                                      *
!***********************************************************************
    N3=1
    34 N2=N3
    nmohulp=nmohulp+1
    if (nmohulp > nmolmax) stop 'Too many molecules in system'
    ias(N2,3+mbond)=nmohulp
    37 FOUND= .FALSE. 
    DO 36 N1=N2+1,NA2
        IF (ias(N1,3+mbond) /= 0) GOTO 36
        DO 35 L=1,mbond
            IF (ias(N1,l+2) == 0) GOTO 36
            IF (ias(ias(N1,l+2),3+mbond) == nmohulp) THEN
                FOUND= .TRUE. 
                ias(N1,3+mbond)=nmohulp
                GOTO 36
            ENDIF
        35 END DO
    36 END DO
    IF (FOUND) GOTO 37
    DO 33 N3=N2+1,NA2
        IF (ias(N3,3+mbond) == 0) GOTO 34
    33 END DO
!***********************************************************************
!                                                                      *
!     The assigned or input molecule numbers are checked for their     *
!     consistency.                                                     *
!                                                                      *
!***********************************************************************
    32 FOUND= .FALSE. 
    DO 42 N1=1,NA2
        DO 41 L=1,mbond
            IF (ias(N1,L+2) == 0) GOTO 42
            IF (ias(ias(N1,L+2),3+mbond) /= ias(N1,3+mbond)) THEN
                FOUND= .TRUE. 
            ENDIF
        41 END DO
    42 END DO
    IF (FOUND) THEN
        STOP 'Mol.nrs. not consistent; maybe wrong cell parameters'
    end if

    return
    end subroutine sortmol

!*********************************************************************
!***********************************************************************

    subroutine readdelphi (qfileh,iend,naold)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(25) :: qfileh
!*********************************************************************
!                                                                    *
!     Read in geometries in Delphi-format (xyz)                      *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readdelphi'
        call timer(65)
        close (65)
    end if

    if (imodfile == 1) then
        open (3,file=qfileh,status='old')
    end if
    nmmax=nmmaxold
    nfc=nfcold
    ibity=1
    iredo=1
    endpo=endpoold
    icell=icellold
    icelo2=icelo2old
    iend=0
    read (3,1000,end=900)qr,qmol
!*********************************************************************
!                                                                    *
!     Read in restraint information (optional)                       *
!                                                                    *
!*********************************************************************
    if (qr == 'R' .OR. qr == 'P' .OR. qr == 'X') then
        qmol=qmol(7:60)
        qmolset(nuge)=qmol
        read (18,1070,end=4,err=4) nrestra
        do i1=1,nrestra
            read (18,1090)irstra(i1,1),irstra(i1,2),rrstra(i1),vkrstr(i1), &
            vkrst2(i1),rrcha(i1)
        end do
        4 continue
    end if
!*********************************************************************
!                                                                    *
!     Read in torsion restraint information (optional)               *
!                                                                    *
!*********************************************************************
    if (qr == 'T' .OR. qr == 'X') then
        if (qr == 'T') then
            qmol=qmol(7:60)
            qmolset(nuge)=qmol
        end if
        read (28,1070,end=6,err=6) nrestrat
        do i1=1,nrestrat
            read (28,1091)irstrat(i1,1),irstrat(i1,2),irstrat(i1,3), &
            irstrat(i1,4),trstra(i1),vkrt(i1),vkr2t(i1),rtcha(i1)
        end do
        6 continue
    end if
!*********************************************************************
!                                                                    *
!     Read in valency angle restraint information (optional)         *
!                                                                    *
!*********************************************************************
    if (qr == 'V') then
        qmol=qmol(7:60)
        qmolset(nuge)=qmol
        read (38,1070,end=7,err=7) nrestrav
        do i1=1,nrestrav
            read (38,1092)irstrav(i1,1),irstrav(i1,2),irstrav(i1,3), &
            vrstra(i1),vkrv(i1),vkr2v(i1)
        end do
        7 continue
    end if
!*********************************************************************
!                                                                    *
!     Read in geometry                                               *
!                                                                    *
!*********************************************************************
    ibh2=0
    iequi=1
    iexco=0
    if (nequi > 0) iequi=0
    axiss(1)=-1.0d0

    if (qr == 'O' .OR. qr == 'L') stop 'Not xyz-format'

    if (qr == 'I') then      !Delphi internal coordinates
        if (nsurp >= 2) stop 'Int.coordinates only with 1 gemetry'
        call staint
        goto 20
    end if

    if (qr == 'B') then      !Previous geometry with volume reduction
        read (3,*)
        vred=(1.0d0-0.01d0*volcha)**(0.33333d0)
        iexco=1
        na=naold
        do i1=1,3
            qmol=qmol
            axiss(i1)=vred*axis(i1)
            angles(i1)=angle(i1)
            do i2=1,na
                c(i2,i1)=vred*c(i2,i1)
            end do
        end do

        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)
        ibity=2

        goto 20
    end if

    if (qr == 'S') then      !Previous geometry with volume expansion
        read (3,*)
        vexp=(1.0d0+0.01d0*volcha)**(0.333333333333d0)
        na=naold
        iexco=1
        do i1=1,3
            qmol=qmol
            axiss(i1)=vexp*axis(i1)
            angles(i1)=angle(i1)
            do i2=1,na
                c(i2,i1)=vexp*c(i2,i1)
            end do
        end do

        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)
        ibity=2

        goto 20
    end if

    if (qr == 'F' .OR. qr == 'Y' .OR. qr == '3' .OR. qr == '5' &
    .or. qr == 'P') then
        kx=0
        ky=0
        kz=0
        ibity=2
        read(3,1005)axiss(1),axiss(2),axiss(3)
        read(3,1005)angles(1),angles(2),angles(3)

        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)

    end if

    if (qr == 'M' .OR. qr == 'A') then
        nmmsav=nmm
        nmm=2
    end if

    if (qr == 'A') nmm=1

    if (qr == 'D') then
        endpo=endpo/25
        nmmax=nmmax*5
        qruid='HIGH PRECISION'
    end if

    if (qr == 'H') then
        nmmax=nmmax/10
        qruid='LOW PRECISION'
    end if

    if (qr == '1' .OR. qr == '5') then
        nmm=1
        nmmax=1
        qruid='SINGLE POINT'
    end if

    if (qr == 'Y') then
        icell=0
        qruid='NO CELL OPT'
    end if

    10 read (3,1100,end=20,err=20)ir,qa(na+1),(c(na+1,i2),i2=1,3)
    qlabel(na+1)=qa(na+1)
    qresi1(na+1)='   '
    qresi2(na+1)=' '
    qresi3(na+1)='     '
    qffty(na+1)='     '
    if (ir == 0) goto 20
    na=na+1

    if (na > nat) then
        write (*,*)'Maximum number of atom exceeded ',na,nat
        stop 'Maximum number of atoms exceeded'
    end if

    goto 10
    20 continue

    if (imodfile == 1) close (3)

    return
    900 iend=1
    return
    1000 format (2x,a1,1x,a60)
    1005 format (3f10.4)
    1070 format (i3)
    1090 format (2i4,2f8.4,f8.6,f10.8)
    1091 format (4i4,2f8.4,3f8.6)
    1092 format (3i4,2f8.4,2f8.6)
    1100 format (i4,1x,a2,3x,3d22.15,1x,a5,1x,i5)
    end subroutine readdelphi
!***********************************************************************
!***********************************************************************

    subroutine readbgf(qfileh,iendf,naold)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(80) :: qromb
    character(2) :: qrom
    character(5) :: quen
    character(5) :: qlabhulp
    character(25) :: qfileh
    character(200) :: qhulp
    !df398
    logical :: match
    character(200) :: temp
!*********************************************************************
!                                                                    *
!     Read in BIOGRAF-geometry                                       *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readbgf'
        call timer(65)
        close (65)
    end if

    iendf=0
    ienread=0
    iredo=0
    qremark(1)=' '
    enmol=zero
    formol=zero
    if (imodfile == 1) then
        open (3,file=qfileh,status='old')
    end if
    10 continue
    read (3,'(a40)',end=900)qromb
    ibity=0
    if (qromb(1:6) == 'BIOGRF') ibity=1
    if (qromb(1:6) == 'XTLGRF') ibity=2
    inobio=0
    if (ibity == 0) then
        !write (*,*)qromb(1:6)
    !     stop 'Unknown Biograf-file'
        inobio=1
    end if
    if (inobio == 1) goto 10

    read (qromb,'(6x,i4)')ibgfversion
    if (ibity == 1) qr='C'
    if (ibity == 2) qr='F'
    iremark=0
    iformat=0
    iline=0
    iexco=0
    iruid=1
    ipropt=0
    vvol=1.0d0
    nmcharge=0
    nmmax=nmmaxold
    nfc=nfcold
    ncha=nchaold
    endpo=endpoold
    icell=icellold
    icelo2=icelo2old
    axiss(1)=-1.0d0

    30 read (3,'(a200)',end=46,err=40)qhulp
    qstrana1(1:200)=qhulp
    iline=iline+1
    irecog=0

    if (qhulp(1:6) == 'DESCRP') then
        read (qhulp,'(7x,a40)',end=46,err=46)temp
        !irecog=1
        !df398> only read structures that are present in trainset.in
        match = .false.
        do j=1,junique
           if (temp == trainset_structs(j)) then
               !write(*,*) 'match found! ', temp,' is indeed defined in trainset.in'
               match = .true.
               irecog=1
               qmol = temp
               exit
           endif
        enddo
        if (match .eqv. .false.) then
           irecog=0
           goto 10
           !write(*,*) 'temp(1:3) = ', temp(1:3)
           !skip lines until we reach 'END' and then start reading
           !a new structure
           do while (temp(1:3) /= 'END')
              read(3,'(a200)',end=46,err=40)temp
              !write(*,*) 'temp(1:3) now is: ', temp(1:3)
              if (temp(1:3) == 'END') then
                 !write(*,*) 'temp = END. going to 30'
                 if (nsurp /= 2) goto 900
                 if (nsurp >= 2) then 
                    goto 45 
                 endif
              !else
              !   goto 45
              end if
           end do
        endif
    end if

    if (qhulp(1:6) == 'REMARK') then
        if (iremark < 20) iremark=iremark+1
        read (qhulp,'(7x,a40)',end=46,err=46)qremark(iremark)
        irecog=1
    end if

    if (qhulp(1:6) == 'FORMAT') then
        if (iformat < 20) iformat=iformat+1
        read(qhulp,'(7x,a40)',end=46,err=46)qformat(iformat)
        irecog=1
    end if

    if (qhulp(1:7) == 'VCHANGE') then
        read (qhulp(8:60),*)vvol
        vred=(1.0+(vvol-1.0))**(0.33333333333333d0)
        iexco=1
        na=naold
        qmol=qmol
        do i1=1,3
            axiss(i1)=vred*axis(i1)
            angles(i1)=angle(i1)
            do i2=1,na
                c(i2,i1)=vred*c(i2,i1)
            end do
        end do

        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0) cosphi=1.0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)
        ibity=2
        irecog=1
    end if

    if (qhulp(1:7) == 'VCHANGX') then
        read (qhulp(8:60),*)vvol
        vred=vvol
        iexco=1
        na=naold
        qmol=qmol
        do i1=1,3
            axiss(i1)=axis(i1)
            angles(i1)=angle(i1)
            do i2=1,na
                c(i2,i1)=c(i2,i1)
            end do
        end do

        axiss(1)=vred*axiss(1)
        do i2=1,na
            c(i2,1)=vred*c(i2,1)
        end do

        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)
        axis(1)=axiss(1)
        ibity=2
        irecog=1
    end if

    if (qhulp(1:7) == 'VCHANGY') then
        read (qhulp(8:60),*)vvol
        vred=vvol
        iexco=1
        na=naold
        qmol=qmol
        do i1=1,3
            axiss(i1)=axis(i1)
            angles(i1)=angle(i1)
            do i2=1,na
                c(i2,i1)=c(i2,i1)
            end do
        end do

        axiss(2)=vred*axiss(2)
        do i2=1,na
            c(i2,2)=vred*c(i2,2)
        end do

        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)
        axis(2)=axiss(2)
        ibity=2
        irecog=1
    end if

    if (qhulp(1:7) == 'VCHANGZ') then
        read (qhulp(8:60),*)vvol
        vred=vvol
        iexco=1
        na=naold
        qmol=qmol

        do i1=1,3
            axiss(i1)=axis(i1)
            angles(i1)=angle(i1)
            do i2=1,na
                c(i2,i1)=c(i2,i1)
            end do
        end do

        axiss(3)=vred*axiss(3)
        do i2=1,na
            c(i2,3)=vred*c(i2,3)
        end do

        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0) cosphi=1.0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)
        axis(3)=axiss(3)
        ibity=2
        irecog=1
    end if

    if (qhulp(1:6) == 'CRYSTX') then
        read (qhulp,'(8x,6f11.5)',end=46,err=46)axiss(1), &
        axiss(2),axiss(3),angles(1),angles(2),angles(3)
        kx=0
        ky=0
        kz=0
        halfa=angles(1)*dgrrdn
        hbeta=angles(2)*dgrrdn
        hgamma=angles(3)*dgrrdn
        sinalf=sin(halfa)
        cosalf=cos(halfa)
        sinbet=sin(hbeta)
        cosbet=cos(hbeta)
        cosphi=(cos(hgamma)-cosalf*cosbet)/(sinalf*sinbet)
        if (cosphi > 1.0d0) cosphi=1.0d0
        sinphi=sqrt(one-cosphi*cosphi)
        tm11=axiss(1)*sinbet*sinphi
        tm21=axiss(1)*sinbet*cosphi
        tm31=axiss(1)*cosbet
        tm22=axiss(2)*sinalf
        tm32=axiss(2)*cosalf
        tm33=axiss(3)
        kx=int(2.0d0*swb/tm11)
        ky=int(2.0d0*swb/tm22)
        kz=int(2.0d0*swb/tm33)
        qr='F'
        if (nmmax == 1 .AND. nmmaxold > 1) qr='5'
        if (icell == 0 .AND. icellold > 0) qr='Y'
        ibity=2
        irecog=1
    end if

    if (qhulp(1:6) == 'PERIOD') then
        read (qhulp,'(7x,i3)',end=46,err=46)iperiod
        irecog=1
    end if

    if (qhulp(1:4) == 'AXES') then
        read (qhulp,'(7x,a3)',end=46,err=46)qbgfaxes
        irecog=1
    end if

    if (qhulp(1:6) == 'SGNAME') then
        read (qhulp,'(7x,a3)',end=46,err=46)qbgfsgn
        irecog=1
    end if

!     if (qhulp(1:5).eq.'CELLS') then
!     read (qhulp,'(7x,*)',end=40,err=40)kx,ky,kz
!     irecog=1
!     end if

    if (qhulp(1:6) == 'HETATM') then
        if (ibgfversion < 400) then
            read (qhulp, &
            '(7x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f15.10,1x,a5,i3,i2,1x,f8.5)' &
            ,end=40,err=40) &
            ir,qlabel(na+1),qresi1(na+1),qresi2(na+1),qresi3(na+1), &
            c(na+1,1),c(na+1,2), &
            c(na+1,3),qffty(na+1),ibgr1(na+1),ibgr2(na+1),chgbgf(na+1)
        else
            stop 'Unsupported Biograf-version'
        end if
        qlabhulp=qlabel(na+1)
        if (qlabhulp(1:1) == ' ') qlabhulp=qlabhulp(2:5)
        if (qlabhulp(1:1) == ' ') qlabhulp=qlabhulp(2:4)
        if (qlabhulp(1:1) == ' ') qlabhulp=qlabhulp(2:3)
        qa(na+1)=qlabhulp(1:2)
        if (qlabhulp(1:1) == 'C') qa(na+1)='C '
        if (qlabhulp(1:2) == 'C ') qa(na+1)='C '
        if (qlabhulp(1:2) == 'Ca') qa(na+1)='Ca'
        if (qlabhulp(1:2) == 'Ce') qa(na+1)='Ce'
        if (qlabhulp(1:2) == 'Cl') qa(na+1)='Cl'
        if (qlabhulp(1:2) == 'Cu') qa(na+1)='Cu'
        if (qlabhulp(1:2) == 'Co') qa(na+1)='Co'
        if (qlabhulp(1:1) == 'H') qa(na+1)='H '
        if (qlabhulp(1:2) == 'H ') qa(na+1)='H '
        if (qlabhulp(1:2) == 'He') qa(na+1)='He'
        if (qlabhulp(1:2) == 'Hf') qa(na+1)='Hf'
        if (qlabhulp(1:1) == 'D') qa(na+1)='D '
        if (qlabhulp(1:2) == 'D ') qa(na+1)='D '
        if (qlabhulp(1:1) == 'N') qa(na+1)='N '
        if (qlabhulp(1:2) == 'N ') qa(na+1)='N '
        if (qlabhulp(1:2) == 'Na') qa(na+1)='Na'
        if (qlabhulp(1:2) == 'Ni') qa(na+1)='Ni'
        if (qlabhulp(1:2) == 'Nb') qa(na+1)='Nb'
        if (qlabhulp(1:2) == 'Ne') qa(na+1)='Ne'
        if (qlabhulp(1:1) == 'O') qa(na+1)='O '
        if (qlabhulp(1:2) == 'O ') qa(na+1)='O '
        if (qlabhulp(1:1) == 'B') qa(na+1)='B '
        if (qlabhulp(1:2) == 'B ') qa(na+1)='B '
        if (qlabhulp(1:2) == 'Ba') qa(na+1)='Ba'
        if (qlabhulp(1:2) == 'Bi') qa(na+1)='Bi'
        if (qlabhulp(1:1) == 'F') qa(na+1)='F '
        if (qlabhulp(1:2) == 'F ') qa(na+1)='F '
        if (qlabhulp(1:2) == 'Fe') qa(na+1)='Fe'
        if (qlabhulp(1:2) == 'P ') qa(na+1)='P '
        if (qlabhulp(1:1) == 'S') qa(na+1)='S '
        if (qlabhulp(1:2) == 'S ') qa(na+1)='S '
        if (qlabhulp(1:2) == 'Sb') qa(na+1)='Sb'
        if (qlabhulp(1:2) == 'Si') qa(na+1)='Si'
        if (qlabhulp(1:2) == 'Se') qa(na+1)='Se'
        if (qlabhulp(1:2) == 'Sn') qa(na+1)='Sn'
        if (qlabhulp(1:1) == 'Y') qa(na+1)='Y '
        if (qlabhulp(1:2) == 'Y ') qa(na+1)='Y '
        if (qlabhulp(1:1) == 'K') qa(na+1)='K '
        if (qlabhulp(1:1) == 'K ') qa(na+1)='K '
        if (qlabhulp(1:2) == 'Kr') qa(na+1)='Kr'
        if (qlabhulp(1:1) == 'V') qa(na+1)='V '
        if (qlabhulp(1:2) == 'V ') qa(na+1)='V '
        if (qlabhulp(1:2) == 'Al ') qa(na+1)='Al'
        if (qlabhulp(1:2) == 'Au ') qa(na+1)='Au'
        if (qlabhulp(1:2) == 'Pd') qa(na+1)='Pd'
        if (qlabhulp(1:2) == 'Pt') qa(na+1)='Pt'
        if (qlabhulp(1:2) == 'Mo') qa(na+1)='Mo'
        if (qlabhulp(1:2) == 'Mg') qa(na+1)='Mg'
        if (qlabhulp(1:2) == 'Ar') qa(na+1)='Ar'
        if (qlabhulp(1:2) == 'Ti') qa(na+1)='Ti'
        if (qlabhulp(1:2) == 'Ru') qa(na+1)='Ru'
        if (qlabhulp(1:2) == 'Rb') qa(na+1)='Rb'
        if (qlabhulp(1:2) == 'Te') qa(na+1)='Te'
        if (qlabhulp(1:2) == 'Li') qa(na+1)='Li'
        if (qlabhulp(1:2) == 'Zr') qa(na+1)='Zr'
        if (qlabhulp(1:1) == 'X') qa(na+1)='X '
        if (qlabhulp(1:2) == 'X ') qa(na+1)='X '
        if (qlabhulp(1:2) == 'Xe') qa(na+1)='Xe'
        if (qlabhulp(1:2) == 'In') qa(na+1)='In'
        if (qlabhulp(1:2) == 'As') qa(na+1)='As'
        if (qlabhulp(1:2) == 'Cs') qa(na+1)='Cs'
        if (qlabhulp(1:2) == 'Cr') qa(na+1)='Cr'
        if (qlabhulp(1:2) == 'Sr') qa(na+1)='Sr'
        na=na+1
        if (na > nat) then
            write (*,*)'Maximum number of atoms: ',nat
            stop 'Maximum number of atoms exceeded; increase nat in cbka.blk'
        end if
        irecog=1
    end if

    if (qhulp(1:6) == 'RUTYPE') then          !run-type identifiers
        irecrun=0
        read (qhulp,'(7x,a40)',end=46,err=46)qruid

        if (qruid(1:10) == 'NORMAL RUN') then
            iruid=0
            irecrun=1
        end if

        if (qruid(1:14) == 'HIGH PRECISION') then
            endpo=endpo/25
            nmmax=nmmax*5
            qr='D'
            iruid=1
            irecrun=1
        end if

        if (qruid(1:13) == 'LOW PRECISION') then
            nmmax=nmmax/10
            qr='H'
            iruid=1
            irecrun=1
        end if

        if (qruid(1:12) == 'SINGLE POINT') then
            iruid=1
            nmmax=1
            qr='1'
            if (ibity == 2) qr='5'
            irecrun=1
        end if

        if (qruid(1:11) == 'NO CELL OPT') then
            iruid=1
            icell=0
            if (ibity == 2) qr='Y'
            irecrun=1
        end if

        if (qruid(1:8) == 'CELL OPT') then
            iruid=1
            icell=1
            iexco=0   !Override from VCHANGE
            read (qruid,'(8x,i6)',end=46,err=46)ncellopt
            if (ncellopt == 2) icell=2 !cell optimisation during energy minimisation
            if (ncellopt == 3) icelo2=4  !c/a optimisation
            if (ncellopt == 4) icelo2=1  !only a optimisation
            if (ncellopt == 5) icelo2=2  !only b optimisation
            if (ncellopt == 6) icelo2=3  !only c optimisation
            if (ncellopt == 7) then
                icelo2=4  !c/a optimisation
                icell=2 !cell optimisation during energy minimisation
            end if
            if (ibity == 2) qr='F'
            irecrun=1
        end if

        if (qruid(1:9) == 'PRESS OPT') then
            iruid=1
            ipropt=1
            istart=10
            read (qruid,'(9x,4f8.4)')vpresopt,vprestax,vprestay,vprestaz
            if (ibity == 2) qr='F'
            irecrun=1
        end if

        if (qruid(1:10) == 'PRESS CALC') then
            iruid=1
            ipropt=2
            if (ibity == 2) qr='F'
            irecrun=1
        end if

        if (qruid(1:6) == 'MAXMOV') then
            iruid=1
            read (qruid,'(6x,i9)',end=46,err=46)nfc
            irecrun=1
        end if

        if (qruid(1:4) == 'REDO') then
            iruid=1
            read (qruid,'(4x,i6)',end=46,err=46)iredo
            irecrun=1
        end if

        if (qruid(1:5) == 'MAXIT') then
            iruid=1
            read (qruid,'(6x,i6)',end=46,err=46)nmmax
            if (qruid(14:18) == 'ENDPO') then
                read (qruid,'(18x,f6.3)',end=46,err=46)endpo
            end if
            irecrun=1
        end if
        if (qruid(1:5) == 'ENDPO') then
            iruid=1
            read (qruid,'(6x,f6.3)',end=46,err=46)endpo
            irecrun=1
        end if
              
        if (qruid(1:9) == 'CHARGEMET') then
            iruid=1
            read (qruid,'(9x,i6)',end=46,err=46)ncha
            irecrun=1
        end if

        if (irecrun == 0) then
            write (*,*)'Warning: ignored RUTYPE identifier ',qruid(1:12)
        end if

        irecog=1
    end if

    if (qhulp(1:14) == 'EQUIV DISTANCE') then
        neqdis=neqdis+1
        if (neqdis > mrestra) &
        stop 'Too many restraints; increase mrestra in cbka.blk'
        istart=15
        call stranal(istart,iend,vout,iout,1)
        istart=15
        call stranal(istart,iend,vout,iout,1)
        ieqdis(neqdis,1)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        ieqdis(neqdis,2)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        ieqdis(neqdis,3)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        ieqdis(neqdis,4)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkeqd1(neqdis)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkeqd2(neqdis)=iout
        istart=iend
        qr='R'
        irecog=1
    end if

    if (qhulp(1:14) == 'BOND RESTRAINT') then
        nrestra=nrestra+1
        if (nrestra > mrestra) &
        stop 'Too many restraints; increase mrestra in cbka.blk'
        istart=15
        call stranal(istart,iend,vout,iout,1)
        irstra(nrestra,1)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        irstra(nrestra,2)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        rrstra(nrestra)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkrstr(nrestra)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkrst2(nrestra)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        rrcha(nrestra)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        itstart(nrestra)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        itend(nrestra)=iout
        istart=iend
    !     read (qhulp,'(15x,2i4,f8.4,f8.2,f8.5,f10.7),end=46,err=46')
    !    $irstra(nrestra,1),irstra(nrestra,2),rrstra(nrestra),
    !    $vkrstr(nrestra),vkrst2(nrestra),rrcha(nrestra)
        qr='R'
        irecog=1
    end if

    if (qhulp(1:15) == 'ANGLE RESTRAINT') then
        nrestrav=nrestrav+1
        if (nrestrav > mrestra) &
        stop 'Too many restraints; increase mrestra in cbka.blk'
        istart=16
        call stranal(istart,iend,vout,iout,1)
        irstrav(nrestrav,1)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        irstrav(nrestrav,2)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        irstrav(nrestrav,3)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vrstra(nrestrav)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkrv(nrestrav)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkr2v(nrestrav)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        rvcha(nrestrav)=vout
        istart=iend
    !     read (qhulp,'(16x,3i4,2f8.2,f8.4,f9.6),end=46,err=46')
    !    $irstrav(nrestrav,1),irstrav(nrestrav,2),irstrav(nrestrav,3),
    !    $vrstra(nrestrav),vkrv(nrestrav),vkr2v(nrestrav),
    !    $rvcha(nrestrav)
        qr='V'
        irecog=1
    end if

    if (qhulp(1:17) == 'TORSION RESTRAINT') then
        nrestrat=nrestrat+1
        if (nrestrat > mrestra) &
        stop 'Too many restraints; increase mrestra in cbka.blk'
        istart=18
        call stranal(istart,iend,vout,iout,1)
        irstrat(nrestrat,1)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        irstrat(nrestrat,2)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        irstrat(nrestrat,3)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        irstrat(nrestrat,4)=iout
        istart=iend

    !     iatyt1=ia(irstrat(nrestrat,1),1)
    !     iatyt2=ia(irstrat(nrestrat,2),1)
    !     iatyt3=ia(irstrat(nrestrat,3),1)
    !     iatyt4=ia(irstrat(nrestrat,4),1)

    !     if (iatyt3.lt.iatyt2) then   !reorganize torsion restraint
    !     irstrat(nrestrat,1)=iatyt4
    !     irstrat(nrestrat,2)=iatyt3
    !     irstrat(nrestrat,3)=iatyt2
    !     irstrat(nrestrat,4)=iatyt1
    !     write (*,*)'Warning: reorganized torsion restraint data'
    !     end if
    
    !     if (iatyt3.eq.iatyt2.and.iatyt4.lt.iatyt1) then   !reorganize torsion restraint
    !     irstrat(nrestrat,1)=iatyt4
    !     irstrat(nrestrat,2)=iatyt3
    !     irstrat(nrestrat,3)=iatyt2
    !     irstrat(nrestrat,4)=iatyt1
    !     write (*,*)'Warning: reorganized torsion restraint data'
    !     end if
            
        call stranal(istart,iend,vout,iout,1)
        trstra(nrestrat)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkrt(nrestrat)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vkr2t(nrestrat)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        rtcha(nrestrat)=vout
        istart=iend
    !     read (qhulp,'(18x,4i4,2f8.2,f8.4,f9.6),end=46,err=46')
    !    $irstrat(nrestrat,1),irstrat(nrestrat,2),irstrat(nrestrat,3),
    !    $irstrat(nrestrat,4),trstra(nrestrat),vkrt(nrestrat),
    !    $vkr2t(nrestrat),rtcha(nrestrat)
        qr='T'
        irecog=1
    end if

    if (qhulp(1:16) == 'MASCEN RESTRAINT') then
        nrestram=nrestram+1
        istart=17
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        irstram(nrestram,1)=0
        if (qstrana2 == 'x') irstram(nrestram,1)=1
        if (qstrana2 == 'y') irstram(nrestram,1)=2
        if (qstrana2 == 'z') irstram(nrestram,1)=3
        if (qstrana2 == 'p') irstram(nrestram,1)=4 !fixed center of mass
        if (irstram(nrestram,1) == 0) &
        stop 'Error in mass centre restraint'
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        irstram(nrestram,2)=iout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        irstram(nrestram,3)=iout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        rmstra1(nrestram)=vout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        if (irstram(nrestram,1) <= 3) irstram(nrestram,4)=iout
        if (irstram(nrestram,1) == 4) rmstra2(nrestram)=vout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        if (irstram(nrestram,1) <= 3) irstram(nrestram,5)=iout
        if (irstram(nrestram,1) == 4) rmstra3(nrestram)=vout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        if (irstram(nrestram,1) <= 3) rmstra2(nrestram)=vout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        if (irstram(nrestram,1) <= 3) rmstra3(nrestram)=vout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        if (irstram(nrestram,1) <= 3) rmcha(nrestram)=vout
        irecog=1
    end if

    if (qhulp(1:14) == 'SYMM RESTRAINT') then
        nrestras=nrestras+1
        istart=15
        call stranal(istart,iend,vout,iout,1)
        irstras(nrestras,1)=iout
        istart=iend
        if (irstras(nrestras,1) > maxsrest) then
            write (*,*)'Too many atoms used in symmetry restraint'
            write (*,*)'Increase maxsrest in cbka.blk'
            stop 'Too many atoms used in symmetry restraint'
        end if
        do i1=1,irstras(nrestras,1)
            call stranal(istart,iend,vout,iout,1)
            irstras(nrestras,1+i1)=iout
            istart=iend
        end do
        call stranal(istart,iend,vout,iout,1)
        qrstras(nrestras)=qstrana1(iend-6:iend)
        write (64,*)qrstras(nrestras)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vksym1(nrestras)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        vksym2(nrestras)=vout
        istart=iend
        irecog=1
    end if

    if (qhulp(1:9) == 'MOLCHARGE') then
        nmcharge=nmcharge+1
        istart=10
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        iat1mc(nmcharge)=iout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        iat2mc(nmcharge)=iout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        vmcha(nmcharge)=vout
        irecog=1
    end if
          
    if (qhulp(1:8) == 'FIXATOMS') then
        istart=9
        call stranal(istart,iend,vout,iout,1)
        if1=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        if2=iout
        do i12=if1,if2
            imove(i12)=0
        end do
        irecog=1
    end if

    if (qhulp(1:11) == 'UNIT ENERGY') then
        eenconv=zero
        read (qhulp,'(14x,a5)',end=46,err=46)quen
        if (quen == 'eV') eenconv=23.0408d0
        if (quen == 'EV') eenconv=23.0408d0
        if (quen == 'ev') eenconv=23.0408d0
        if (quen == 'h') eenconv=627.5d0
        if (quen == 'H') eenconv=627.5d0
        if (quen == 'kcal') eenconv=1.0d0
        if (quen == 'kCal') eenconv=1.0d0
        if (quen == 'KCAL') eenconv=1.0d0
        if (eenconv == zero) then
            write (*,*)quen,': unknown energy unit; assuming kcal/mol'
            eenconv=1.0d0
        end if
        irecog=1
    end if

    if (qhulp(1:6) == 'ENERGY') then
        read (qhulp(7:80),*,end=46,err=46)enmol
        ienread=1
        irecog=1
    end if

    if (qhulp(1:6) == 'GEOUPD') then
        icgeopt(nprob)=0
        icgeo=0
        irecog=1
    end if
     
    if (qhulp(1:9) == 'NO GEOUPD') then
        icgeopt(nprob)=1
        icgeo=1
        irecog=1
    end if

    if (qhulp(1:9) == 'FREQUENCY') then
        ifreqset(nprob)=1
        ifreq=1
        irecog=1
    end if

!     if (qhulp(1:5).eq.'FORCE') then
!     read (qhulp(6:80),*,end=46,err=46)formol
!     ienread=1
!     irecog=1
!     end if

    if (qhulp(1:6) == 'FFIELD') goto 30
    if (qhulp(1:6) == 'CONECT') goto 30
    if (qhulp(1:5) == 'ORDER') goto 30
    if (qhulp(1:1) == '#') goto 30
    if (qhulp(1:3) == 'END') goto 45

    if (irecog == 0) then
        write (*,*)'Warning: ignored line starting with: ',qhulp(1:10)
    end if
     
    goto 30
          
    40 write (*,*)'Error on line ',iline+1,' of Biograf-input'
    stop
    45 read (3,*,err=46,end=46)
    46 continue
    if (ienread == 1) then
        if (eenconv == zero) then
            write (*,*)'No energy unit given; assuming kcal/mol'
            eenconv=1.0
        end if
        enmol=enmol*eenconv                !Convert energies to kcal/mol
    end if
          
    namov=0                            !calculate number of moving atoms
    do i1=1,na
        if (imove(i1) == 1) namov=namov+1
    end do

    if (imodfile == 1) close (3)
    return
    900 iendf=1
    return
    end subroutine readbgf
!***********************************************************************
!***********************************************************************

    subroutine readpdb (iendf)

!***********************************************************************
    include 'cbka.blk'
    include 'opt.blk'
    character(200) :: qhulp
!*********************************************************************
!                                                                    *
!     Read in .pdb-geometry                                          *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readpdb'
        call timer(65)
        close (65)
    end if
    iendf=1
    qmol='pdb_in'
    5 read (3,'(a200)',end=10,err=900) qhulp
    qstrana1(1:200)=qhulp
    istart=1
    call stranal(istart,iend,vout,iout,1)
    istart=iend

    if (qstrana2(1:6) == 'HEADER') then
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        qmol=qstrana2(1:30)
    end if

    if (qstrana2(1:6) == 'HETATM' .OR. qstrana2(1:4) == 'ATOM') then
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        qa(na+1)=qstrana2(1:2)
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        c(na+1,1)=vout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        c(na+1,2)=vout
        call stranal(istart,iend,vout,iout,1)
        istart=iend
        c(na+1,3)=vout
        na=na+1
    end if

    if (qstrana2(1:3) == 'END' .OR. qstrana2(2:4) == 'END') then
        iendf=0
        goto 10
    end if

    goto 5
    10 continue
    return
    900 write (*,*)'Error reading in .pdb-format'
    stop 'Error reading in .pdb-format'
    end subroutine readpdb
!***********************************************************************
!***********************************************************************

    subroutine readtreg(isumattreg)

!***********************************************************************
    include 'cbka.blk'
    dimension isumattreg(mtreg)
    character(200) :: qrom
!*********************************************************************
!                                                                    *
!     Read in temperature regime                                     *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readtreg'
        call timer(65)
        close (65)
    end if
    ntrc=0
    open (19,file='tregime.in',status='old',err=60)
    10 read (19,'(a200)',end=50,err=900)qrom
    qstrana1(1:200)=qrom
    if (qrom(1:1) == '#') goto 10
    istart=1
    ntrc=ntrc+1
    if (ntrc > mtreg) then
        write (*,*)'Too many temperature regimes in tregime.in;', &
        ' inrease mtreg in cbka.blk'
        stop 'Too many temperature regimes in tregime.in'
    end if
    call stranal(istart,iend,vout,iout,1)
    nittc(ntrc)=iout
    istart=iend

    if (ntrc > 1) then
        if (nittc(ntrc) < nittc(ntrc-1)) then
            ntrc=ntrc-1
            write (*,*)'Warning: wrong order or empty line in tregime.in'
            write (*,*)'Ignored lines below iteration:',nittc(ntrc)
            goto 50
        end if
    end if

    call stranal(istart,iend,vout,iout,1)
    nntreg(ntrc)=iout
    if (nntreg(ntrc) > mtzone) then
        write (*,*)'Too many temperature zones in tregime.in;', &
        ' inrease mtzone in cbka.blk'
        stop 'Too many temperature zones in tregime.in'
    end if
    istart=iend
    isumattreg(ntrc)=0
    do i1=1,nntreg(ntrc)
        call stranal(istart,iend,vout,iout,1)
        ia1treg(ntrc,i1)=iout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        ia2treg(ntrc,i1)=iout
        istart=iend
        isumattreg(ntrc)=isumattreg(ntrc)+1+ia2treg(ntrc,i1)- &
        ia1treg(ntrc,i1)
        call stranal(istart,iend,vout,iout,1)
        tsettreg(ntrc,i1)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        tdamptreg(ntrc,i1)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        dttreg(ntrc,i1)=vout
        istart=iend
    end do
    goto 10
    50 continue
    close (19)
    60 continue
    return
    900 stop 'Error reading tregime.in'
    end subroutine readtreg
!***********************************************************************
!***********************************************************************

    subroutine readvreg

!***********************************************************************
    include 'cbka.blk'
    character(200) :: qrom
!*********************************************************************
!                                                                    *
!     Read in volume regime                                          *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readvreg'
        call timer(65)
        close (65)
    end if
    nvrc=0
    open (19,file='vregime.in',status='old',err=60)
    10 read (19,'(a200)',end=50,err=900)qrom
    qstrana1(1:200)=qrom
    if (qrom(1:1) == '#') goto 10
    istart=1
    nvrc=nvrc+1
    if (nvrc > mvreg) then
        write (*,*)'Too many volume regimes in vregime.in;', &
        ' inrease mvreg in cbka.blk'
        stop 'Too many volume regimes in vregime.in'
    end if

    call stranal(istart,iend,vout,iout,1)
    nitvc(nvrc)=iout
    istart=iend

    if (nvrc > 1) then
        if (nitvc(nvrc) < nitvc(nvrc-1)) then
            nvrc=nvrc-1
            write (*,*)'Warning: wrong order or empty line in vregime.in'
            write (*,*)'Ignored lines below iteration:',nitvc(nvrc)
            goto 50
        end if
    end if

    call stranal(istart,iend,vout,iout,1)
    nnvreg(nvrc)=iout
    if (nnvreg(nvrc) > mvzone) then
        write (*,*)'Too many volume regimes in vregime.in;', &
        ' inrease mvzone in cbka.blk'
        stop 'Too many volume zones in vregime.in'
    end if
    istart=iend
    do i1=1,nnvreg(nvrc)
        call stranal(istart,iend,vout,iout,1)
        if (qstrana2(1:1) /= 'a' .AND. qstrana2(1:1) /= 'b' .AND. &
        qstrana2(1:1) /= 'c' .AND. qstrana2(1:4) /= 'alfa' .AND. &
        qstrana2(1:4) /= 'beta' .AND. qstrana2(1:5) /= 'gamma') then
            write (*,*)qstrana2
            write (*,*)'Invalid cell parameter type in vregime.in ;', &
            ' use a,b,c,alfa,beta or gamma'
            stop 'Invalid cell parameter type in vregime.in'
        end if
        qvtype(nvrc,i1)=qstrana2
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        dvvreg(nvrc,i1)=vout
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        ivsca(nvrc,i1)=1
        if (qstrana2(1:1) == 'n') ivsca(nvrc,i1)=0
        istart=iend
    end do
    goto 10
    50 continue
    close (19)
    60 continue
    return
    900 stop 'Error reading vregime.in'
    end subroutine readvreg
!***********************************************************************
!***********************************************************************

    subroutine readereg

!***********************************************************************
    include 'cbka.blk'
    character(200) :: qrom
!*********************************************************************
!                                                                    *
!     Read in electric field regime                                  *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readereg'
        call timer(65)
        close (65)
    end if
    nerc=0
    open (19,file='eregime.in',status='old',err=60)
    10 read (19,'(a200)',end=50,err=900)qrom
    qstrana1(1:200)=qrom
    if (qrom(1:1) == '#') goto 10
    istart=1
    nerc=nerc+1
    if (nerc > mereg) then
        write (*,*)'Too many electric field regimes in eregime.in;', &
        ' inrease mereg in cbka.blk'
        stop 'Too many electric field regimes in eregime.in'
    end if
    call stranal(istart,iend,vout,iout,1)
    nitec(nerc)=iout

    if (nerc > 1) then
        if (nitec(nerc) < nitec(nerc-1)) then
            nerc=nerc-1
            write (*,*)'Warning: wrong order or empty line in eregime.in'
            write (*,*)'Ignored lines below iteration:',nitec(nerc)
            goto 50
        end if
    end if
          
    istart=iend
    call stranal(istart,iend,vout,iout,1)
    nnereg(nerc)=iout
    if (nnereg(nerc) > mezone) then
        write (*,*)'Too many electric field zones in eregime.in;', &
        ' inrease mezone in cbka.blk'
        stop 'Too many electric field zones in eregime.in'
    end if
    istart=iend
    do i1=1,nnereg(nerc)
        call stranal(istart,iend,vout,iout,1)
        if (qstrana2(1:1) /= 'x' .AND. qstrana2(1:1) /= 'y' .AND. &
        qstrana2(1:1) /= 'z') then
            write (*,*)qstrana2
            write (*,*)'Invalid field direction in eregime.in ;', &
            ' use x,y or z'
            stop 'Invalid field direction in eregime.in'
        end if
        qetype(nerc,i1)=qstrana2
        istart=iend
        call stranal(istart,iend,vout,iout,1)
        ereg(nerc,i1)=vout
        istart=iend
    end do
    goto 10
    50 continue
    close (19)
    60 continue
    return
    900 stop 'Error reading eregime.in'
    end subroutine readereg
!***********************************************************************
!***********************************************************************

    subroutine readpiston

!***********************************************************************
    include 'cbka.blk'
    character(200) :: qrom
!*********************************************************************
!                                                                    *
!     Read in piston characteristics                                 *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readpiston'
        call timer(65)
        close (65)
    end if
    npist=0
    open (19,file='piston.in',status='old',err=60)
    10 read (19,'(a200)',end=50,err=900)qrom
    qstrana1(1:200)=qrom
    if (qrom(1:1) == '#') goto 10
    istart=1
    npist=npist+1
    if (npist > mpist) then
        write (*,*)'Too many piston regimes in piston.in;', &
        ' inrease mpist in cbka.blk'
        stop 'Too many piston regimes in piston.in'
    end if
    call stranal(istart,iend,vout,iout,1)
    nitpist(npist)=iout
    istart=iend

    if (npist > 1) then
        if (nitpist(npist) < nitpist(npist-1)) then
            npist=npist-1
            write (*,*)'Warning: wrong order or empty line in piston.in'
            write (*,*)'Ignored lines below iteration:',nitpist(npist)
            goto 50
        end if
    end if

    call stranal(istart,iend,vout,iout,1)

    if (qstrana2(1:1) /= 'x' .AND. qstrana2(1:1) /= 'y' .AND. &
    qstrana2(1:1) /= 'z') then
        write (*,*)qstrana2
        write (*,*)'Invalid piston direction in piston.in ;', &
        ' use x,y or z'
        stop 'Invalid piston direction in piston.in'
    end if

    qptype(npist)=qstrana2
    istart=iend

    call stranal(istart,iend,vout,iout,1)
    edeep2(npist)=vout
    istart=iend
    call stranal(istart,iend,vout,iout,1)
    rdeep2(npist)=vout
    istart=iend
    call stranal(istart,iend,vout,iout,1)
    pshft2(npist)=vout
    istart=iend
    call stranal(istart,iend,vout,iout,1)
    rcut2(npist)=vout
    istart=iend
    call stranal(istart,iend,vout,iout,1)
    speedp2(npist)=vout
    istart=iend

    goto 10
    50 continue
    close (19)
    60 continue
    return
    900 stop 'Error reading piston.in'
    end subroutine readpiston
!***********************************************************************
!***********************************************************************

    subroutine readaddmol

!***********************************************************************
    include 'cbka.blk'
    character(80) :: qromb
    character(200) :: qhulp
    character(5) :: qlabhulp
!*********************************************************************
!                                                                    *
!     Read in molecule coordinates. This molecule will be added to   *
!     the system at regular intervals                                *
!     Accepts only .bgf-format                                       *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readaddmol'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Set default values                                             *
!                                                                    *
!*********************************************************************
    iaddfreq=-1   !frequency of molecule addition; <0: no addition
    xadd=-9000.0d0  !x-coordinate for added molecule; <-5000.0: random
    yadd=-9000.0d0  !y-coordinate for added molecule; <-5000.0: random
    zadd=-9000.0d0  !z-coordinate for added molecule; <-5000.0: random
    iveladd=1     !1: random initial velocities; 2: read in velocities
! rom addmol.vel
    addist=-1.00d0  !Minimum distance between added molecule and rest
! f system. < 0.0: do not check
    nadattempt=10  !Number of attempts at adding the molecule
    taddmol=-1.0d0  !Temperature added molecule. <0.0: system temperature
    icentp=0      !1: place new molecule close to centre-of-mass
    dcentm=3.0d0    !distance from centre-of-mass
    open (19,file='addmol.bgf',status='old',err=60)
    read (19,'(a40)',end=900,err=900)qromb
    if (qromb(1:6) /= 'BIOGRF') then
        write (*,*)'addmol.bgf should start with BIOGRF'
        stop 'addmol.bgf should start with BIOGRF'
    end if
    naa=0
    iline=0
    30 read (19,'(a200)',end=900,err=900)qhulp
    irecog=0
    iline=iline+1

    if (qhulp(1:6) == 'DESCRP') then
        irecog=1
    end if
     
    if (qhulp(1:6) == 'FORMAT') then
        irecog=1
    end if

    if (qhulp(1:6) == 'REMARK') then
        irecog=1
    end if

    if (qhulp(1:6) == 'HETATM') then
        irecog=1
        read (qhulp,'(7x,i5,1x,a5,1x,3x,1x,1x,1x,5x,3f10.5)' &
        ,end=900,err=900) &
        ir,qlabhulp,cadd(naa+1,1),cadd(naa+1,2),cadd(naa+1,3)
        if (qlabhulp(1:1) == ' ') qlabhulp=qlabhulp(2:5)
        if (qlabhulp(1:1) == ' ') qlabhulp=qlabhulp(2:4)
        if (qlabhulp(1:1) == ' ') qlabhulp=qlabhulp(2:3)
        if (qlabhulp(1:1) == 'C ') qadd(naa+1)='C '
        if (qlabhulp(1:2) == 'Ca') qadd(naa+1)='Ca'
        if (qlabhulp(1:2) == 'Cl') qadd(naa+1)='Cl'
        if (qlabhulp(1:2) == 'Cu') qadd(naa+1)='Cu'
        if (qlabhulp(1:2) == 'Co') qadd(naa+1)='Co'
        if (qlabhulp(1:1) == 'H ') qadd(naa+1)='H '
        if (qlabhulp(1:2) == 'He') qadd(naa+1)='He'
        if (qlabhulp(1:1) == 'N ') qadd(naa+1)='N '
        if (qlabhulp(1:2) == 'Ni') qadd(naa+1)='Ni'
        if (qlabhulp(1:2) == 'Ne') qadd(naa+1)='Ne'
        if (qlabhulp(1:1) == 'O ') qadd(naa+1)='O '
        if (qlabhulp(1:1) == 'B ') qadd(naa+1)='B '
        if (qlabhulp(1:1) == 'F ') qadd(naa+1)='F '
        if (qlabhulp(1:2) == 'Fe') qadd(naa+1)='Fe'
        if (qlabhulp(1:1) == 'P ') qadd(naa+1)='P '
        if (qlabhulp(1:1) == 'S ') qadd(naa+1)='S '
        if (qlabhulp(1:1) == 'Y ') qadd(naa+1)='Y '
        if (qlabhulp(1:2) == 'Al') qadd(naa+1)='Al'
        if (qlabhulp(1:2) == 'Au') qadd(naa+1)='Au'
        if (qlabhulp(1:2) == 'Si') qadd(naa+1)='Si'
        if (qlabhulp(1:2) == 'Pt') qadd(naa+1)='Pt'
        if (qlabhulp(1:2) == 'Mo') qadd(naa+1)='Mo'
        if (qlabhulp(1:2) == 'Mg') qadd(naa+1)='Mg'
        if (qlabhulp(1:2) == 'Ar') qadd(naa+1)='Ar'
        if (qlabhulp(1:2) == 'Zr') qadd(naa+1)='Zr'
        if (qlabhulp(1:2) == 'Ba') qadd(naa+1)='Ba'
        if (qlabhulp(1:2) == 'X ') qadd(naa+1)='X '
        if (qlabhulp(1:2) == 'Xe') qadd(naa+1)='Xe'
        if (qlabhulp(1:2) == 'Cs') qadd(naa+1)='Cs'
        if (qlabhulp(1:2) == 'Sr') qadd(naa+1)='Sr'
        ityadd(naa+1)=0
        do i1=1,nso  !Find force field type
            if (qadd(naa+1) == qas(i1)) ityadd(naa+1)=i1
        end do
        if (ityadd(naa+1) == 0) then
            write (*,*) 'Unknown atom type:',qadd(naa+1)
            stop 'Unknown atom type'
        end if
        naa=naa+1
    end if

    if (qhulp(1:7) == 'FREQADD') then
        irecog=1
        read (qhulp,'(8x,i6)',end=900,err=900) iaddfreq
    end if

    if (qhulp(1:6) == 'VELADD') then
        irecog=1
        read (qhulp,'(8x,i6)',end=900,err=900) iveladd
    end if

    if (qhulp(1:6) == 'STARTX') then
        irecog=1
        read (qhulp,'(7x,f8.2)',end=900,err=900) xadd
    end if

    if (qhulp(1:6) == 'STARTY') then
        irecog=1
        read (qhulp,'(7x,f8.2)',end=900,err=900) yadd
    end if

    if (qhulp(1:6) == 'STARTZ') then
        irecog=1
        read (qhulp,'(7x,f8.2)',end=900,err=900) zadd
    end if

    if (qhulp(1:6) == 'ADDIST') then
        irecog=1
        read (qhulp,'(7x,f8.2)',end=900,err=900) addist
    end if

    if (qhulp(1:8) == 'NATTEMPT') then
        irecog=1
        read (qhulp,'(9x,i6)',end=900,err=900) nadattempt
    end if

    if (qhulp(1:7) == 'TADDMOL') then
        irecog=1
        read (qhulp,'(8x,f8.2)',end=900,err=900) taddmol
    end if

    if (qhulp(1:6) == 'ICENTP') then
        irecog=1
        read (qhulp,'(7x,i6)',end=900,err=900) icentp
    end if

    if (qhulp(1:6) == 'DCENTM') then
        irecog=1
        read (qhulp,'(7x,f8.2)',end=900,err=900) dcentm
    end if

    if (qhulp(1:6) == 'FFIELD') goto 30
    if (qhulp(1:6) == 'CONECT') goto 30
    if (qhulp(1:5) == 'ORDER') goto 30
    if (qhulp(1:1) == '#') goto 30
    if (qhulp(1:3) == 'END') goto 45

    if (irecog == 0) then
        write (*,*)'Warning: ignored line starting with: ',qhulp(1:10)
    end if

    goto 30

    45 continue
    close (19)
    if (iveladd == 2) then
        open (19,file='addmol.vel',status='old',err=800)
        read (19,*)
        read (19,'(3d24.15)',err=850,end=850) &
        ((veladd(j,i),j=1,3),i=1,naa)
        close (19)
    end if
!***********************************************************************
!                                                                      *
!     Place molecule at origin                                         *
!                                                                      *
!***********************************************************************
    ccx=0.0d0
    ccy=0.0d0
    ccz=0.0d0
    do i1=1,naa
        ccx=ccx+cadd(i1,1)/float(naa)
        ccy=ccy+cadd(i1,2)/float(naa)
        ccz=ccz+cadd(i1,3)/float(naa)
    end do
    do i1=1,naa
        cadd(i1,1)=cadd(i1,1)-ccx
        cadd(i1,2)=cadd(i1,2)-ccy
        cadd(i1,3)=cadd(i1,3)-ccz
    end do

    60 continue
    return
    800 stop 'Error opening addmol.vel'
    850 stop 'Error or end of file reading addmol.vel'
    900 write (*,*)'Error or end-of-file reading addmol.bgf on line:', &
    iline
    return
    end subroutine readaddmol
!***********************************************************************
!***********************************************************************

    subroutine readtrans

!***********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Read in translations from previous simulations                 *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In readtrans'
        call timer(65)
        close (65)
    end if
    open (19,file='translate.in',status='old',err=60)
    read (19,*)natrans

!     if (natrans.ne.na) then
!     write (*,*)natrans,na
!     write (*,*)'Wrong number of atoms in translate.in-file'
!     stop 'Wrong number of atoms in translate.in-file'
!     end if

    do i1=1,natrans
        read (19,*)ir,ix,iy,iz
        id(i1,1)=ix
        id(i1,2)=iy
        id(i1,3)=iz
    end do

    60 continue
    return
    end subroutine readtrans
!***********************************************************************
!***********************************************************************

    SUBROUTINE INVMAT(NDIM1,NDIM2,RMAT,N1)

!***********************************************************************
!     Matrix inverse                                                   *
!     N1=number of occupied elements                                   *
!                                                                      *
!***********************************************************************

    INCLUDE 'cbka.blk'

    dimension rmat(ndim1,ndim2)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In invmat'
        call timer(65)
        close (65)
    end if
    ONE=1.0D0
    ZERO=0.0D0

    DO  K=1,N1
        DIAG=RMAT(K,K)
        IF (DIAG == 0.0d0) STOP 'INVMAT'
        RMAT(K,K)=ONE/DIAG
    !---RIJ---
        DO  L=1,K-1
            RMAT(K,L)=RMAT(K,L)/DIAG
        END DO
        DO  L=K+1,N1
            RMAT(K,L)=RMAT(K,L)/DIAG
        END DO
    !---INVERSIE STUK----
        DO  M=1,K-1
            DO  L=1,K-1
                RMAT(M,L)=RMAT(M,L)-RMAT(M,K)*RMAT(K,L)
            END DO
        ENDDO
        DO  M=K+1,N1
            DO  L=1,K-1
                RMAT(M,L)=RMAT(M,L)-RMAT(M,K)*RMAT(K,L)
            END DO
        ENDDO
    !---GEWONE KOLOM---
        DO M=1,K-1
            RMAT(M,K)=-RMAT(M,K)/DIAG
        END DO
        DO  M=K+1,N1
            RMAT(M,K)=-RMAT(M,K)/DIAG
        END DO
    !---GEWONE STUK-----
        DO  M=1,K-1
            DO  L=K+1,N1
                RMAT(M,L)=RMAT(M,L)+RMAT(M,K)*RMAT(K,L)*DIAG
          END DO
        ENDDO
        DO  M=K+1,N1
            DO  L=K+1,N1
                RMAT(M,L)=RMAT(M,L)+RMAT(M,K)*RMAT(K,L)*DIAG
            END DO
        ENDDO
    !-----
    END DO

    RETURN

    END SUBROUTINE INVMAT

!***********************************************************************
!***********************************************************************

    SUBROUTINE MATSYM4(N,NDIMM1,NDIMM2,NDIMVX,NDIMVY,RMAT,VECX,VECY)

!***********************************************************************
!                                                                      *
!     Construction of the lower triangle and backsubstitution.         *
!     VERTICAL: ONLY for SYMMETRIC matrices not necessarily stored in  *
!     square arrays.                                                   *
!                                                                      *
!***********************************************************************

    INCLUDE 'cbka.blk'

    DIMENSION RMAT(neem,neem),VECX(neem),VECY(neem)

    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In matsym4'
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
    END DO

    DO  K=N,1,-1
        VECX(K)=VECY(K)
        DO  L=K+1,N
            VECX(K)=VECX(K)-RMAT(L,K)*VECX(L)
        END DO
    ENDDO

    RETURN

    END SUBROUTINE MATSYM4

!***********************************************************************
!***********************************************************************

    subroutine cgsolve(n,a,lda,x,b,mdstep,convg)

!***********************************************************************
    implicit real*8 (a-h,o-z),integer(i-n)

!     RPM. Conjugate gradient scheme to solve Ax=b

    INTEGER :: n                 ! INPUT: The order of the problem
    REAL*8 :: a(lda,n)           ! INPUT: The matrix to be solved
    INTEGER :: lda               ! INPUT: The leading dimension of A
    REAL*8 :: x(n)               ! OUTPUT: The solution x=A^-1*b
    REAL*8 :: b(n)               ! INPUT: The vector b in Ax=b
    INTEGER :: mdstep
    REAL*8 :: convg              ! Convergence criterion
    REAL*8 :: cptime,ustime,systime

!     Local temp storage
    dimension r(lda),p(lda),p_old(lda),q(lda)

    one = 1.0d0
    zero = 0.0d0
!WRITE (*,*) 'size of EEM matrix:', lda,'x', n
    maxiter = 150
!     Form r=b-Ax
    call dgemv('N',n,n,one,a,lda,x,1,zero,r,1) ! r = Ax
    call DSCAL(n,-one,r,1)
    call daxpy(n,one,b,1,r,1)
    rho_old = one

    do iter = 1,maxiter
        rho = ddot(n,r,1,r,1)
!WRITE(60,*) rho
        if (rho < convg) goto 666

        call dcopy(n,r,1,p,1) ! p=r
        if (iter > 1) then
            beta = rho/rho_old
            call daxpy(n,beta,p_old,1,p,1)
        endif
        call dgemv('N',n,n,one,a,lda,p,1,zero,q,1) ! q=Ap
        alpha = rho/ddot(n,p,1,q,1)
        call daxpy(n,alpha,p,1,x,1) ! x += alpha*p
        call daxpy(n,-alpha,q,1,r,1) ! r -= alpha*q
        call dcopy(n,p,1,p_old,1)
        rho_old = rho
    enddo
    !write(*,*)mdstep,'cgsolve:Warning:maxiter exceeded. Continuing.'
    666 continue
!     write (*,*)mdstep,iter

    return
    end subroutine cgsolve

!***********************************************************************
!***********************************************************************

    subroutine fillmatrix(n,a,lda)

!***********************************************************************
    implicit real*8 (a-h,o-z),integer(i-n)

!     fill in top half of matrix
    INTEGER :: n                 ! INPUT: The size of the matrix a
    REAL*8 :: a(lda,lda)         ! INPUT: The matrix to be filled in
! OUTPUT: Contains the full matrix
    INTEGER :: lda               ! INPUT: The leading dimension of a

    do i = 1, n
        do j = 1, i-1
            a(j,i) = a(i,j)
        enddo
    enddo

    return
    end subroutine fillmatrix

!***********************************************************************
!*********************************************************************

    subroutine dipmom(naold,dpmm,xdip,ydip,zdip,xdir,ydir,zdir)

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Calculate dipole moment                                        *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In dipmom'
        call timer(65)
        close (65)
    end if
!***********************************************************************
!                                                                      *
!     CONVERSION FACTOR TO DEBYE UNITS IS CALCULATED                   *
!     THE CALCULATION IS INITIALIZED                                   *
!                                                                      *
!***********************************************************************
     
    ELCHG=1.60217733D-19       ! [C]       = [As]
    CLIGHT=2.99792458D8        !             [m/s]
    DBCONV=ONE/(CLIGHT*ELCHG*1.0D11)
     
    CHCPX=ZERO
    CHCPY=ZERO
    CHCPZ=ZERO
    CHCMX=ZERO
    CHCMY=ZERO
    CHCMZ=ZERO
    XDIP=ZERO
    YDIP=ZERO
    ZDIP=ZERO
    XGRD=ZERO
    YGRD=ZERO
    ZGRD=ZERO
!***********************************************************************
!                                                                      *
!     CALCULATION OF MAGNITUDE AND CENTRES OF POSITIVE AND NEGATIVE    *
!     CHARGES                                                          *
!                                                                      *
!***********************************************************************
     
    if (na == 0) na=naold
    CHRG=ZERO
    DO 4 K1=1,NA
        CHK1=CH(K1)
        IF (CHK1 == ZERO) GOTO 4
        IF (CHK1 < ZERO) GOTO 3
        CHRG=CHRG+CHK1
        CHCPX=CHCPX+CHK1*C(K1,1)
        CHCPY=CHCPY+CHK1*C(K1,2)
        CHCPZ=CHCPZ+CHK1*C(K1,3)
        GOTO 4
        3 CHCMX=CHCMX-CHK1*C(K1,1)
        CHCMY=CHCMY-CHK1*C(K1,2)
        CHCMZ=CHCMZ-CHK1*C(K1,3)
    4 END DO
     
!***********************************************************************
!                                                                      *
!     CALCULATION OF DISTANCE BETWEEN CENTRES AND OF DIPOLE MOMENT     *
!     IN DEBIJE UNITS                                                  *
!                                                                      *
!***********************************************************************
     
    CHDSTX=CHCPX-CHCMX
    CHDSTY=CHCPY-CHCMY
    CHDSTZ=CHCPZ-CHCMZ
    DPMM=SQRT(CHDSTX*CHDSTX+CHDSTY*CHDSTY+CHDSTZ*CHDSTZ)/DBCONV

    xdip=zero
    ydip=zero
    zdip=zero
    xdir=zero
    ydir=zero
    zdir=zero
    IF(DPMM > 1.0D-4) then
        XDIP=HALF*(CHCPX+CHCMX)/CHRG
        YDIP=HALF*(CHCPY+CHCMY)/CHRG
        ZDIP=HALF*(CHCPZ+CHCMZ)/CHRG
        GRTST=MAX(CHDSTX,CHDSTY,CHDSTZ)
        XDIR=-CHDSTX/GRTST
        YDIR=-CHDSTY/GRTST
        ZDIR=-CHDSTZ/GRTST
    end if

    open (51,file='dipole.out',status='unknown',position='append')
    write (51,100)qmol
    write (51,110)dpmm,xdip,ydip,zdip,xdir,ydir,zdir
    close (51)
     
    100 format ('System identifier:',a40)
    110 format ('Dipole moment (Debye):',f12.4,' Location:',3f12.4, &
    ' Direction (-side):',3f12.4)
    return
    end subroutine dipmom
!***********************************************************************
!*********************************************************************

    subroutine layvel

!*********************************************************************
    include 'cbka.blk'
    dimension vent(3,nat),vmc1(3),vmc2(3),exf(3)
!*********************************************************************
!                                                                    *
!     Add velocity vector to certain atoms                           *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In layvel'
        call timer(65)
        close (65)
    end if

    exf(1)=exfx
    exf(2)=exfy
    exf(3)=exfz

    do i1=1,3
        vmc1(i1)=0.0d0
        vmc2(i1)=0.0d0
    end do
    xmasmd1=0.0d0
    xmasmd2=0.0d0

    do i1=1,na
        if (ibgr2(i1) <= 0) then
            xmasmd1=xmasmd1+amas(i1)
        else
            xmasmd2=xmasmd2+amas(i1)
        end if
    end do

    do i1=1,na
        do i2=1,3
            if(ibgr2(i1) <= 0) then
                vmc1(i2)=vmc1(i2)+vel(i2,i1)*xmasat(i1)
            else
                vmc2(i2)=vmc2(i2)+vel(i2,i1)*xmasat(i1)
            endif
        !      ekin=ekin+imove(i1)*xmasat(i1)*vel(i2,i1)*vel(i2,i1)
        end do
    end do

    do k=1,3
        vmc1(k)=vmc1(k)/xmasmd1
        vmc2(k)=vmc2(k)/xmasmd2
    end do

    do i1=1,na
        do i2=1,3
            if(ibgr2(i1) <= 0) then
                vent(i2,i1)=vel(i2,i1)-vmc1(i2)
            else
                vent(i2,i1)=vel(i2,i1)-vmc2(i2)
            end if
            ekin=ekin+imove(i1)*xmasat(i1)*vent(i2,i1)*vent(i2,i1)
        end do
    end do

    ekin=0.50d0*ekin/convmd
    eksav=ekin
    do i1=1,3
        do i2=1,na
            if(ibgr2(i2) == -1) then
                vel(i1,i2)=vel(i1,i2)-vmc1(i1)
            else if(ibgr2(i2) == 2) then
                vel(i1,i2)=vel(i1,i2)+(exf(i1)*1.0d+12)-vmc2(i1)
            end if
        end do
    end do

    return
    end subroutine layvel
!***********************************************************************
!***********************************************************************
     
    subroutine timer(nunit)
     
!***********************************************************************
    include 'cbka.blk'
!     real timear
!     real tarray(2)
!     real *4 function cptime()
!     call dtime (timear,tarray)
! df398 commented. secnds is not supported in nagfortran and this subroutine is not really needed.
!    cptime=secnds(0.0001)
!    ustime=systime
!    systime=cptime-starttime
!    write (nunit,100)systime-ustime,systime
!    return
!    100 format (' Time since last call:',f20.4, &
!    ' Total time:',f20.4)
    end subroutine timer
!***********************************************************************
!***********************************************************************
