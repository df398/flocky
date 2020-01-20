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
!*******************************************************************


    subroutine calval

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),b(3),j(3),dradc(3,3),drbdc(3,3),dtdc(3,3), &
    dargdc(3,3),dndc(3,3),dadc(3),dbdc(3)
!*********************************************************************
!                                                                    *
!     Calculate valency angles and their derivatives to cartesian    *
!     coordinates                                                    *
!     Valency angle energies are calculated in valang                *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In calval'
        call timer(65)
        close (65)
    end if

    third=1.0d0/3.0d0
    twothird=2.0d0/3.0d0
    dadc(1)=-1.0d0
    dadc(2)=1.0d0
    dadc(3)=0.0d0
    dbdc(1)=0.0d0
    dbdc(2)=1.0d0
    dbdc(3)=-1.0d0
    do k1=1,3
        do k2=1,3
            dradc(k1,k2)=0.0d0
            drbdc(k1,k2)=0.0d0
        end do
    end do
    if (nval == 0) return
     
    do  i1=1,nval
        ity=iv(i1,1)
        j(1)=iv(i1,2)
        j(2)=iv(i1,3)
        j(3)=iv(i1,4)
    !*********************************************************************
    !                                                                    *
    !     Determine valency angle                                        *
    !                                                                    *
    !*********************************************************************
        la=iv(i1,5)
        lb=iv(i1,6)
        ivl1=ibsym(la)
        ivl2=ibsym(lb)
        isign1=1
        isign2=1
        if (j(1) < j(2)) isign1=-1
        if (j(3) < j(2)) isign2=-1
        rla=rbo(la)
        rlb=rbo(lb)
        ix1=isign1*nvlx(ivl1)
        iy1=isign1*nvly(ivl1)
        iz1=isign1*nvlz(ivl1)
        ix2=isign2*nvlx(ivl2)
        iy2=isign2*nvly(ivl2)
        iz2=isign2*nvlz(ivl2)
         
    !     call dista2(j(2),j(1),dis,a(1),a(2),a(3))
    !     call dista2(j(2),j(3),dis,b(1),b(2),b(3))
        a(1)=c(j(2),1)-c(j(1),1)+ix1*tm11
        a(2)=c(j(2),2)-c(j(1),2)+ix1*tm21+iy1*tm22
        a(3)=c(j(2),3)-c(j(1),3)+ix1*tm31+iy1*tm32+iz1*tm33
        b(1)=c(j(2),1)-c(j(3),1)+ix2*tm11
        b(2)=c(j(2),2)-c(j(3),2)+ix2*tm21+iy2*tm22
        b(3)=c(j(2),3)-c(j(3),3)+ix2*tm31+iy2*tm32+iz2*tm33

        poem=rla*rlb
        tel=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
        arg=tel/poem
        arg2=arg*arg
        s1ma22=1.0d0-arg2
        if (s1ma22 < 1.0d-20) s1ma22=1.0d-20
        s1ma2=sqrt(s1ma22)
        if (arg > 1.0d0) arg=1.0d0
        if (arg < -1.0d0) arg=-1.0d0
        hl=acos(arg)
        h(i1)=hl

    !*********************************************************************
    !                                                                    *
    !     Calculate derivative valency angle to cartesian coordinates    *
    !                                                                    *
    !*********************************************************************
        if (j(1) == ib(la,2)) then
            do k1=1,3
                dradc(k1,1)=drdc(la,k1,1)
                dradc(k1,2)=drdc(la,k1,2)
            end do
        else
            do k1=1,3
                dradc(k1,1)=drdc(la,k1,2)
                dradc(k1,2)=drdc(la,k1,1)
            end do
        end if
        if (j(2) == ib(lb,2)) then
            do k1=1,3
                drbdc(k1,2)=drdc(lb,k1,1)
                drbdc(k1,3)=drdc(lb,k1,2)
            end do
        else
            do k1=1,3
                drbdc(k1,2)=drdc(lb,k1,2)
                drbdc(k1,3)=drdc(lb,k1,1)
            end do
        end if
        do k1=1,3
            do k2=1,3
                dndc(k1,k2)=rla*drbdc(k1,k2)+rlb*dradc(k1,k2)
                dtdc(k1,k2)=a(k1)*dbdc(k2)+b(k1)*dadc(k2)
                dargdc(k1,k2)=(dtdc(k1,k2)-arg*dndc(k1,k2))/poem
                dhdc(i1,k1,k2)=-dargdc(k1,k2)/s1ma2
            end do
        end do
              
    END DO

    return
    end subroutine calval
!*********************************************************************
!*********************************************************************

    subroutine boncor

!*********************************************************************
    include 'cbka.blk'
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In boncor'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Correction for overcoordination and 1-3 bond orders            *
!                                                                    *
!*********************************************************************
          
    do  i1=1,nbon
        ibt=ib(i1,1)
        j1=ib(i1,2)
        j2=ib(i1,3)
       ! df398: check for ovc and v13cor parameters (flags) in ffield bond section. 
       !        If both are smaller than 0.001 don't correct bond order for 
       !        over-coordination and 1-3 bond orders.
       if (ovc(ibt) < 0.001d0 .AND. v13cor(ibt) < 0.001d0) then
            idbo1(i1)=2
            idbo(i1,1)=j1
            idbo(i1,2)=j2
            do k1=1,3                                   
                dbondc(i1,k1,1)=dbodc(i1,k1,1)          ! df398 total b.o derivative atom 1
                dbondc(i1,k1,2)=dbodc(i1,k1,2)          ! total b.o derivative atom 2
                !dbosindc(i1,k1,1)=dbosidc(i1,k1,1)      ! sigma bond b.o derivative atom 1
                !dbosindc(i1,k1,2)=dbosidc(i1,k1,2)      !                                2
                dbopindc(i1,k1,1)=dbopidc(i1,k1,1)      ! pi                             1
                dbopindc(i1,k1,2)=dbopidc(i1,k1,2)      !                                2
                dbopi2ndc(i1,k1,1)=dbopi2dc(i1,k1,1)    ! pipi                           1
                dbopi2ndc(i1,k1,2)=dbopi2dc(i1,k1,2)    ! pipi                           2
            end do
            goto 10
       end if
        boo=bo(i1)
        bopio=bopi(i1)
        bopi2o=bopi2(i1)
        iti=ia(j1,1)
        itj=ia(j2,1)
        aboi=abo(j1)
        aboj=abo(j2)
        vp131=sqrt(bo131(iti)*bo131(itj))
        vp132=sqrt(bo132(iti)*bo132(itj))
        vp133=sqrt(bo133(iti)*bo133(itj))
        corrtot=1.0d0
        dbodsboi1=zero
        dbodsboj1=zero
        ! df398: check ovc flag to decide if overcoordination
        !        correction will be performed. minus => always correct.
	!        (original: 0.001d0)
        if (ovc(ibt) > 0.001d0) then
            ovi=aboi-aval(iti)            ! df398 calculation of b.o overcoordination for atom i
            ovj=aboj-aval(itj)
        !*********************************************************************
        !                                                                    *
        !     Correction for overcoordination                                *
        !                                                                    *
        !*********************************************************************
            exphu1=exp(-vpar(2)*ovi)                             ! df398 vpar(2) is p_boc2 ; vpar(1) is p_boc1
            exphu2=exp(-vpar(2)*ovj)
            exp11=exp(-vpar(1)*ovi)
            exp21=exp(-vpar(1)*ovj)
            exphu12=(exphu1+exphu2)                              
            ovcor=-(1.0d0/vpar(2))*log(0.50d0*exphu12)           ! equation f3 correction
        !     huli=((1.0/ovc(ibt))*aval(iti)+exp11+exp21)
        !     hulj=((1.0/ovc(ibt))*aval(itj)+exp11+exp21)
            huli=aval(iti)+exp11+exp21                           ! val_i + f2 correction
            hulj=aval(itj)+exp11+exp21                           
            corr1=huli/(huli+ovcor)                              ! (val_i + f2 correction) / (val_i + f2 + f3)
            corr2=hulj/(hulj+ovcor)
            corrtot=0.50d0*(corr1+corr2)                         ! equation f1 correction

            dbodsboi1=0.50d0*(-vpar(1)*exp11/(huli+ovcor)- &                   ! derviative of ... ??
            (corr1/(huli+ovcor))* &
            (-vpar(1)*exp11+exphu1/exphu12)-vpar(1)*exp11/(hulj+ovcor)- &
            (corr2/(hulj+ovcor))*(-vpar(1)*exp11+exphu1/exphu12))

            dbodsboj1=0.50d0*(-vpar(1)*exp21/(huli+ovcor)- &                   ! derivative of ... ??
            (corr1/(huli+ovcor))* &
            (-vpar(1)*exp21+exphu2/exphu12)-vpar(1)*exp21/(hulj+ovcor)- &
            (corr2/(hulj+ovcor))*(-vpar(1)*exp21+exphu2/exphu12))
        end if
    !*********************************************************************
    !                                                                    *
    !     Correction for 1-3 bond orders                                 *
    !                                                                    *
    !*********************************************************************
        dbodsboi2=zero
        dbodsboj2=zero
        bocor1=1.0d0
        bocor2=1.0d0
        ! df398: check v13cor flag if 1-3 bond-order correction is performed
	!        minus => always correct. (original: 0.001d0)
        if (v13cor(ibt) > 0.001d0) then
        !     ovi2=aboi-valf(iti)
        !     ovj2=aboj-valf(itj)
            ovi2=aboi-vval3(iti)                           !Modification for Pt
            ovj2=aboj-vval3(itj)
        !     ovi2=aboi-aval(iti)
        !     ovj2=aboj-aval(itj)
            cor1=vp131*boo*boo-ovi2
            cor2=vp131*boo*boo-ovj2
        !     exphu3=v13cor(ibt)*exp(-vp132*cor1+vp133)
        !     exphu4=v13cor(ibt)*exp(-vp132*cor2+vp133)
            exphu3=exp(-vp132*cor1+vp133)
            exphu4=exp(-vp132*cor2+vp133)
            bocor1=1.0d0/(1.0d0+exphu3)                     ! df398 equation f4 correction for atom i
            bocor2=1.0d0/(1.0d0+exphu4)                     ! equation f5 for atom j
            dbodsboi2=-bocor1*bocor1*bocor2*vp132*exphu3    ! derivative of f4 for atom i wrt to b.o
            dbodsboj2=-bocor1*bocor2*bocor2*vp132*exphu4    ! derivative of f5 for atom j wrt to b.o
        end if

        bo(i1)=boo*corrtot*bocor1*bocor2                    ! df398 corrected_sigma_b.o = uncorrected_sigma_b.o * f1 * f4 * f5
	!df398 comment following to prevent discontinuity        
	!------------------------------------------------
	!if (bo(i1) < 1.0d-10) bo(i1)=zero
	! -----------------------------------------------
        corrtot2=corrtot*corrtot
        bopi(i1)=bopio*corrtot2*bocor1*bocor2               ! corrected pi b.o
        bopi2(i1)=bopi2o*corrtot2*bocor1*bocor2             ! corrected double pi b.o
	!df398 comment following to prevent discontinuity
	!------------------------------------------
        !if (bopi(i1) < 1.0d-10) bopi(i1)=zero
        !if (bopi2(i1) < 1.0d-10) bopi2(i1)=zero
	!------------------------------------------

        dbodboo=corrtot*bocor1*bocor2+corrtot* &                  ! various b.o derivatives wrt to uncorrected b.o
        bocor1*bocor1*bocor2*boo*vp132*vp131*2.0d0*boo*exphu3+ &
        corrtot*bocor1*bocor2*bocor2*boo* &
        vp132*vp131*exphu4*2.0d0*boo

        dbopidbopio=corrtot2*bocor1*bocor2

        dbopidboo=corrtot2* &
        bocor1*bocor1*bocor2*boo*vp132*vp131*2.0d0*bopio*exphu3+ &
        corrtot2*bocor1*bocor2*bocor2*boo* &
        vp132*vp131*exphu4*2.0d0*bopio

        dbopi2dbopi2o=corrtot2*bocor1*bocor2

        dbopi2dboo=corrtot2* &
        bocor1*bocor1*bocor2*boo*vp132*vp131*2.0d0*bopi2o*exphu3+ &
        corrtot2*bocor1*bocor2*bocor2*boo* &
        vp132*vp131*exphu4*2.0d0*bopi2o

        dbodsboit=boo*dbodsboi1*bocor1*bocor2+boo*corrtot*dbodsboi2
        dbodsbojt=boo*dbodsboj1*bocor1*bocor2+boo*corrtot*dbodsboj2

        vhui=2.0d0*corrtot*dbodsboi1*bocor1*bocor2+corrtot2*dbodsboi2
        vhuj=2.0d0*corrtot*dbodsboj1*bocor1*bocor2+corrtot2*dbodsboj2
        dbopidsboit=bopio*vhui
        dbopidsbojt=bopio*vhuj

        dbopi2dsboit=bopi2o*vhui
        dbopi2dsbojt=bopi2o*vhuj

    !*********************************************************************
    !                                                                    *
    !     Calculate bond order derivatives                               *
    !                                                                    *
    !*********************************************************************
        idbo1(i1)=2+ia(j1,2)+ia(j2,2)
        idbo(i1,1)=j1
        idbo(i1,2)=j2
        nco=0
        do k1=1,3                                          ! df397 loop over cartesian coordinates for all bonds
            dbondc(i1,k1,1)=dbodc(i1,k1,1)*dbodboo
            dbondc(i1,k1,2)=dbodc(i1,k1,2)*dbodboo
        !     dbosindc(i1,k1,1)=dbosidc(i1,k1,1)*dbosidboo
        !     dbosindc(i1,k1,2)=dbosidc(i1,k1,2)*dbosidboo
            dbopindc(i1,k1,1) =dbopidc(i1,k1,1)*dbopidbopio+dbodc(i1,k1,1)*dbopidboo
            dbopindc(i1,k1,2) =dbopidc(i1,k1,2)*dbopidbopio+dbodc(i1,k1,2)*dbopidboo
            dbopi2ndc(i1,k1,1)=dbopi2dc(i1,k1,1)*dbopi2dbopi2o+dbodc(i1,k1,1)*dbopi2dboo
            dbopi2ndc(i1,k1,2)=dbopi2dc(i1,k1,2)*dbopi2dbopi2o+dbodc(i1,k1,2)*dbopi2dboo
        end do
        do i2=1,ia(j1,2)
            ihl=0
            iob=ia(j1,2+i2)
            if (iob < j1) ihl=1
            ncubo=nubon2(j1,i2)
            idbo(i1,2+nco+1)=iob
            do k1=1,3
                dbondc(i1,k1,1)=dbondc(i1,k1,1)+dbodc(ncubo,k1,1+ihl)*dbodsboit
                dbondc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbodsboit

            !     dbosindc(i1,k1,1)=dbosindc(i1,k1,1)+
            !    $dbodc(ncubo,k1,1+ihl)*dbosidsboit
            !     dbosindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbosidsboit

                dbopindc(i1,k1,1)=dbopindc(i1,k1,1)+dbodc(ncubo,k1,1+ihl)*dbopidsboit
                dbopindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopidsboit

                dbopi2ndc(i1,k1,1)=dbopi2ndc(i1,k1,1)+dbodc(ncubo,k1,1+ihl)*dbopi2dsboit
                dbopi2ndc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopi2dsboit

            end do
            nco=nco+1
        end do
        do i2=1,ia(j2,2)
            ihl=0
            iob=ia(j2,2+i2)
            if (iob < j2) ihl=1
            ncubo=nubon2(j2,i2)
            idbo(i1,2+nco+1)=iob
            do k1=1,3
                dbondc(i1,k1,2)=dbondc(i1,k1,2)+dbodc(ncubo,k1,1+ihl)*dbodsbojt
                dbondc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbodsbojt

            !     dbosindc(i1,k1,2)=dbosindc(i1,k1,2)+
            !    $dbodc(ncubo,k1,1+ihl)*dbosidsbojt
            !     dbosindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbosidsbojt

                dbopindc(i1,k1,2)=dbopindc(i1,k1,2)+dbodc(ncubo,k1,1+ihl)*dbopidsbojt
                dbopindc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopidsbojt

                dbopi2ndc(i1,k1,2)=dbopi2ndc(i1,k1,2)+dbodc(ncubo,k1,1+ihl)*dbopi2dsbojt
                dbopi2ndc(i1,k1,2+nco+1)=dbodc(ncubo,k1,2-ihl)*dbopi2dsbojt

            end do
            nco=nco+1
        end do

    10 END DO

    do i1=1,na
        abo(i1)=zero
    end do
!    do i1=1,na
!        do i2=1,ia(i1,2)
!           iob=ia(i1,2+i2)
!           ncubo=nubon2(i1,i2)
!           abo(i1)=abo(i1)+bo(ncubo)
!        end do
!    end do
    do i1=1,nbon
        j1=ib(i1,2)
        j2=ib(i1,3)
        abo(j1)=abo(j1)+bo(i1)
        if (j1 /= j2) abo(j2)=abo(j2)+bo(i1)
    end do

    15 continue
    return
    end subroutine boncor
!*********************************************************************
!*********************************************************************

    subroutine lonpar

!*********************************************************************
    include 'cbka.blk'
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In lonpar'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Calculate lone pair energy and first derivatives               *
!                                                                    *
!*********************************************************************
    elp=zero
    do i1=1,na
    !*********************************************************************
    !                                                                    *
    !     Determine number of lone pairs on atoms
    !                                                                    *
    !*********************************************************************
        ity=ia(i1,1)
    !     voptlp=0.50*(stlp(ity)-aval(ity))+vlp2(ity)       !For Si-lp
        voptlp=0.50d0*(stlp(ity)-aval(ity))
        vlp(i1)=zero
        vund=abo(i1)-stlp(ity)
        vlph=2.0d0*int(vund/2.0d0)
        vlpex=vund-vlph
        vp16h=vpar(16)-1.0d0

        expvlp=exp(-vpar(16)*(2.0d0+vlpex)*(2.0d0+vlpex))
        dvlpdsbo(i1)=-vpar(16)*2.0d0*(2.0d0+vlpex)*expvlp
        vlp(i1)=expvlp-int(vund/2.0d0)
    !*********************************************************************
    !                                                                    *
    !     Calculate lone pair energy                                     *
    !                                                                    *
    !*********************************************************************
        diffvlp=voptlp-vlp(i1)
        exphu1=exp(-75.0d0*diffvlp)
    !     exphu2=exp(75.0*diffvlp)
        hulp1=1.0d0/(1.0d0+exphu1)
    !     hulp2=1.0/(1.0+exphu2)
    !     elph=vlp1(ity)*diffvlp*hulp1+vlp2(ity)*diffvlp*hulp2   !Stabilize extra lone pair
        elph=vlp1(ity)*diffvlp*hulp1
        estrain(i1)=estrain(i1)+elph
    !     delpdvlp=-vlp1(ity)*hulp1-vlp1(ity)*diffvlp*hulp1*hulp1*   !Stabilize extra lone pair
    !    $75.0*exphu1-vlp2(ity)*hulp2+vlp2(ity)*diffvlp*hulp2*hulp2*
    !    $75.0*exphu2
        delpdvlp=-vlp1(ity)*hulp1-vlp1(ity)*diffvlp*hulp1*hulp1* &
        & 75.0d0*exphu1
        elp=elp+elph
        delpdsbo=delpdvlp*dvlpdsbo(i1)
    !*********************************************************************
    !                                                                    *
    !     Calculate first derivative of lone pair energy to              *
    !     cartesian coordinates                                          *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then
            do i3=1,ia(i1,2)
                iob=ia(i1,2+i3)
                ncubo=nubon2(i1,i3)
                do i4=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i4)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+delpdsbo*dbondc(ncubo,k1,i4)
                    end do
                end do
            end do
        else
            do i3=1,ia(i1,2)
                iob=ia(i1,2+i3)
                ncubo=nubon2(i1,i3)
                do i4=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i4)
                    ix=nmpx(i1,ihu)
                    iy=nmpy(i1,ihu)
                    iz=nmpz(i1,ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        delpdsbo*dbondc(ncubo,k1,i4)
                    end do
                end do
            end do
        end if

    end do

    return
    end subroutine lonpar
!*********************************************************************
!*********************************************************************

    subroutine covbon

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Calculate bond energy and first derivatives                    *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In covbon'
        call timer(65)
        close (65)
    end if
    eb=0.0d0
    if (nbon == 0) return
!*********************************************************************
!                                                                    *
!     Calculate bond energies                                        *
!                                                                    *
!*********************************************************************
    do i1=1,nbon

        boa=bo(i1)
    !     if (boa.lt.cutof2) goto 20
        j1=ib(i1,2)
        j2=ib(i1,3)
        vsymm=1.0d0
        if (j1 == j2) vsymm=0.5d0

        bopia=bopi(i1)
        bopi2a=bopi2(i1)
        bosia=boa-bopia-bopi2a

        if (bosia < zero) bosia=zero
        it1=ia(j1,1)
        it2=ia(j2,1)
        ibt=ib(i1,1)
        de1h=vsymm*de1(ibt)
        de2h=vsymm*de2(ibt)
        de3h=vsymm*de3(ibt)

        bopo1=bosia**psp(ibt)
        exphu1=exp(psi(ibt)*(1.0d0-bopo1))
        ebh=-de1h*bosia*exphu1-de2h*bopia-de3h*bopi2a
        estrain(j1)=estrain(j1)+0.50d0*ebh
        estrain(j2)=estrain(j2)+0.50d0*ebh

        debdbo=-de1h*exphu1+de1h*exphu1*psp(ibt)*psi(ibt)*bopo1
        debdbopi=-de2h
        debdbopi2=-de3h

        eb=eb+ebh

        if (icpres == 0) then
            do i2=1,idbo1(i1)
                ihu=idbo(i1,i2)
                do k1=1,3
                    d(k1,ihu)=d(k1,ihu)+debdbo*(dbondc(i1,k1,i2)-dbopindc(i1,k1,i2)- &
                    dbopi2ndc(i1,k1,i2))+ &
                    debdbopi*dbopindc(i1,k1,i2)+ &
                    debdbopi2*dbopi2ndc(i1,k1,i2)
                end do
            end do
        else
            do i2=1,idbo1(i1)
                ihu=idbo(i1,i2)
                ix=nmpx(j1,ihu)
                iy=nmpy(j1,ihu)
                iz=nmpz(j1,ihu)
                kcell=14+ix+3*iy+9*iz
                do k1=1,3
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+debdbo* &
                    (dbondc(i1,k1,i2)-dbopindc(i1,k1,i2)-dbopi2ndc(i1,k1,i2))+ &
                    debdbopi*dbopindc(i1,k1,i2)+ &
                    debdbopi2*dbopi2ndc(i1,k1,i2)
                end do
            end do
        end if
    !*********************************************************************
    !                                                                    *
    !     Stabilisation terminal triple bond                             *
    !                                                                    *
    !*********************************************************************
        if (boa < 1.00d0) goto 20
         if ((qa(j1) == 'C ' .AND. qa(j2) == 'O ') .OR. &
         (qa(j1) == 'O ' .AND. qa(j2) == 'C ')) then

            ba=(boa-2.50d0)*(boa-2.50d0)
            exphu=exp(-vpar(8)*ba)
            oboa=abo(j1)-boa
            obob=abo(j2)-boa
            exphua1=exp(-vpar(4)*oboa)
            exphub1=exp(-vpar(4)*obob)
            ovoab=abo(j1)-aval(it1)+abo(j2)-aval(it2)
            exphuov=exp(vpar(5)*ovoab)
            hulpov=1.0d0/(1.0d0+25.0d0*exphuov)
                  
            estriph=vpar(11)*exphu*hulpov*(exphua1+exphub1)
            estrain(j1)=estrain(j1)+0.50d0*estriph
            estrain(j2)=estrain(j2)+0.50d0*estriph
            eb=eb+estriph

            decobdbo=vpar(4)*vpar(11)*exphu*hulpov*(exphua1+exphub1) &
            -2.0d0*vpar(11)*vpar(8)*(boa-2.50d0)*hulpov*exphu* &
            (exphua1+exphub1)
            decobdboua=-25.0d0*vpar(5)*vpar(11)*exphu*exphuov*hulpov*hulpov* &
            (exphua1+exphub1)-vpar(11)*exphu*vpar(4)*hulpov*exphua1
            decobdboub=-25.0d0*vpar(5)*vpar(11)*exphu*exphuov*hulpov*hulpov* &
            (exphua1+exphub1)-vpar(11)*exphu*vpar(4)*hulpov*exphub1

            if (icpres == 0) then

                do i2=1,idbo1(i1)
                    ihu=idbo(i1,i2)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+decobdbo*dbondc(i1,k1,i2)
                    end do
                end do

                do i3=1,ia(j1,2)
                    iob=ia(j1,2+i3)
                    ncubo=nubon2(j1,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        do k1=1,3
                            d(k1,ihu)=d(k1,ihu)+decobdboua*dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do

                do i3=1,ia(j2,2)
                    iob=ia(j2,2+i3)
                    ncubo=nubon2(j2,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        do k1=1,3
                            d(k1,ihu)=d(k1,ihu)+decobdboub*dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do

            else

                do i2=1,idbo1(i1)
                    ihu=idbo(i1,i2)
                    ix=nmpx(j1,ihu)
                    iy=nmpy(j1,ihu)
                    iz=nmpz(j1,ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        decobdbo*dbondc(i1,k1,i2)
                    end do
                end do

                do i3=1,ia(j1,2)
                    iob=ia(j1,2+i3)
                    ncubo=nubon2(j1,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        ix=nmpx(j1,ihu)
                        iy=nmpy(j1,ihu)
                        iz=nmpz(j1,ihu)
                        kcell=14+ix+3*iy+9*iz
                        do k1=1,3
                            dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                            decobdboua*dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do

                do i3=1,ia(j2,2)
                    iob=ia(j2,2+i3)
                    ncubo=nubon2(j2,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        ix=nmpx(j1,ihu)
                        iy=nmpy(j1,ihu)
                        iz=nmpz(j1,ihu)
                        kcell=14+ix+3*iy+9*iz
                        do k1=1,3
                            dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                            decobdboub*dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do

            end if 
        end if !df398 close the commented if statement on C-O check
    !*********************************************************************
    !                                                                    *
    !     Calculate penalty for lone pairs sharing a bond                *
    !                                                                    *
    !*********************************************************************
    !     elpbh=de2(ibt)*(2.0-boa)*boa*vlp(j1)*vlp(j2)
    !     estrain(j1)=estrain(j1)+0.50*elpbh
    !     estrain(j2)=estrain(j2)+0.50*elpbh
    !     write (65,'(2i4,8f12.4)')j1,j2,boa,vlp(j1),vlp(j2),elpbh
    !     delpbdbo=-de2(ibt)*vlp(j1)*vlp(j2)*boa+de2(ibt)*vlp(j1)*
    !    $vlp(j2)*(2.0-boa)
    !     delpbdvlp1=de2(ibt)*(2.0-boa)*boa*vlp(j2)
    !     delpbdvlp2=de2(ibt)*(2.0-boa)*boa*vlp(j1)
    !     delpbdsbo1=delpbdvlp1*dvlpdsbo(j1)
    !     delpbdsbo2=delpbdvlp2*dvlpdsbo(j2)
    !     eb=eb+elpbh

    !     do i2=1,idbo1(i1)
    !     ihu=idbo(i1,i2)
    !     do k1=1,3
    !     d(k1,ihu)=d(k1,ihu)+delpbdbo*dbondc(i1,k1,i2)
    !     end do
    !     end do

    !     do i3=1,ia(j1,2)
    !     iob=ia(j1,2+i3)
    !     ncubo=nubon2(j1,i3)
    !     do i4=1,idbo1(ncubo)
    !     ihu=idbo(ncubo,i4)
    !     do k1=1,3
    !     d(k1,ihu)=d(k1,ihu)+delpbdsbo1*dbondc(ncubo,k1,i4)
    !     end do
    !     end do
    !     end do

    !     do i3=1,ia(j2,2)
    !     iob=ia(j2,2+i3)
    !     ncubo=nubon2(j2,i3)
    !     do i4=1,idbo1(ncubo)
    !     ihu=idbo(ncubo,i4)
    !     do k1=1,3
    !     d(k1,ihu)=d(k1,ihu)+delpbdsbo2*dbondc(ncubo,k1,i4)
    !     end do
    !     end do
    !     end do

    20 END DO
    return
    end subroutine covbon
!*********************************************************************
!*********************************************************************

    subroutine ovcor

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Calculate atom energy                                          *
!     Correction for over- and undercoordinated atoms                *
!                                                                    *
!*********************************************************************
    dimension vlptemp(nat)
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In ovcor'
        call timer(65)
        close (65)
    end if
    do i1=1,na
        ity1=ia(i1,1)
        vlptemp(i1)=vlp(i1)
        if (amas(ity1) > 21.0d0) vlptemp(i1)=0.50d0*(stlp(ity1)-aval(ity1))  !Only for 1st-row elements
    end do
    25 ea=zero
    eaot=zero
    eaut=zero

    do 30 i1=1,na
        ity1=ia(i1,1)
        dfvl=1.0d0
        if (amas(ity1) > 21.0d0) dfvl=0.0d0  !Only for 1st-row elements
    !*********************************************************************
    !                                                                    *
    !     Calculate overcoordination energy                              *
    !     Valency is corrected for lone pairs                            *
    !                                                                    *
    !*********************************************************************
        voptlp=0.50d0*(stlp(ity1)-aval(ity1))
        diffvlph=dfvl*(voptlp-vlptemp(i1))
    !*********************************************************************
    !                                                                    *
    !     Determine coordination neighboring atoms                       *
    !                                                                    *
    !*********************************************************************
        sumov=0.0d0
        sumov2=0.0d0
        do i3=1,ia(i1,2)
            iat2=ia(i1,2+i3)
            ity2=ia(iat2,1)
            ncubo=nubon2(i1,i3)
            ibt=ib(ncubo,1)
            voptlp2=0.50d0*(stlp(ity2)-aval(ity2))
            diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
            sumov=sumov+(bopi(ncubo)+bopi2(ncubo))* &
            (abo(iat2)-aval(ity2)-diffvlp2)
            sumov2=sumov2+vover(ibt)*de1(ibt)*bo(ncubo)
        end do

        exphu1=exp(vpar(32)*sumov)
        vho=1.0d0/(1.0d0+vpar(33)*exphu1)
        diffvlp=diffvlph*vho
              
        vov1=abo(i1)-aval(ity1)-diffvlp
        dvov1dsumov=diffvlph*vpar(32)*vpar(33)*vho*vho*exphu1
        exphuo=exp(vovun(ity1)*vov1)
        hulpo=1.0d0/(1.0d0+exphuo)

        hulpp=(1.0d0/(vov1+aval(ity1)+1d-8))

        eah=sumov2*hulpp*hulpo*vov1
        estrain(i1)=estrain(i1)+eah
        deadvov1=-sumov2*hulpp*hulpp*vov1*hulpo+ &
        sumov2*hulpp*hulpo-sumov2*hulpp*vov1*vovun(ity1)* &
        hulpo*hulpo*exphuo
        ea=ea+eah

    !*********************************************************************
    !                                                                    *
    !     Calculate first derivative of overcoordination energy to       *
    !     cartesian coordinates                                          *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then

            do i3=1,ia(i1,2)
                iob=ia(i1,2+i3)
                ncubo=nubon2(i1,i3)
                ibt=ib(ncubo,1)
                deadbo=vover(ibt)*de1(ibt)*hulpp*hulpo*vov1
                do i4=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i4)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+deadvov1*(1.0d0+dfvl*vho*dvlpdsbo(i1))* &
                        dbondc(ncubo,k1,i4)+deadbo*dbondc(ncubo,k1,i4)
                    end do
                end do
            end do
             
            do i2=1,ia(i1,2)
                 
                iat2=ia(i1,2+i2)
                ity2=ia(iat2,1)
                nbosa=nubon2(i1,i2)
                deadvov2=deadvov1*dvov1dsumov*(bopi(nbosa)+bopi2(nbosa))
                 
                voptlp2=0.50d0*(stlp(ity2)-aval(ity2))
                diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
                deadpibo=deadvov1*dvov1dsumov*(abo(iat2)-aval(ity2)-diffvlp2)
                 
                do i4=1,idbo1(nbosa)
                    ihu=idbo(nbosa,i4)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+deadpibo*(dbopindc(nbosa,k1,i4)+ &
                        dbopi2ndc(nbosa,k1,i4))
                    end do
                end do
                 
                do i3=1,ia(iat2,2)
                    iob=ia(iat2,2+i3)
                    ncubo=nubon2(iat2,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        do k1=1,3
                            d(k1,ihu)=d(k1,ihu)+deadvov2*(1.0d0+dfvl*dvlpdsbo(iat2))* &
                            dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do
                 
            end do

        else

            do i3=1,ia(i1,2)
                iob=ia(i1,2+i3)
                ncubo=nubon2(i1,i3)
                ibt=ib(ncubo,1)
                deadbo=vover(ibt)*de1(ibt)*hulpp*hulpo*vov1
                do i4=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i4)
                    ix=nmpx(i1,ihu)
                    iy=nmpy(i1,ihu)
                    iz=nmpz(i1,ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+deadvov1* &
                        (1.0d0+dfvl*vho*dvlpdsbo(i1))* &
                        dbondc(ncubo,k1,i4)+deadbo*dbondc(ncubo,k1,i4)
                    end do
                end do
            end do

            do i2=1,ia(i1,2)

                iat2=ia(i1,2+i2)
                ity2=ia(iat2,1)
                nbosa=nubon2(i1,i2)
                deadvov2=deadvov1*dvov1dsumov*(bopi(nbosa)+bopi2(nbosa))

                voptlp2=0.50d0*(stlp(ity2)-aval(ity2))
                diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
                deadpibo=deadvov1*dvov1dsumov*(abo(iat2)-aval(ity2)-diffvlp2)

                do i4=1,idbo1(nbosa)
                    ihu=idbo(nbosa,i4)
                    ix=nmpx(i1,ihu)
                    iy=nmpy(i1,ihu)
                    iz=nmpz(i1,ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+deadpibo* &
                        (dbopindc(nbosa,k1,i4)+dbopi2ndc(nbosa,k1,i4))
                    end do
                end do

                do i3=1,ia(iat2,2)
                    iob=ia(iat2,2+i3)
                    ncubo=nubon2(iat2,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        ix=nmpx(i1,ihu)
                        iy=nmpy(i1,ihu)
                        iz=nmpz(i1,ihu)
                        kcell=14+ix+3*iy+9*iz
                        do k1=1,3
                            dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+deadvov2* &
                            (1.0d0+dfvl*dvlpdsbo(iat2))*dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do

            end do
        end if

    !*********************************************************************
    !                                                                    *
    !     Calculate undercoordination energy                             *
    !                                                                    *
    !*********************************************************************
        if (valp1(ity1) < zero) goto 30  !skip undercoordination
        exphu2=exp(vpar(10)*sumov)
        vuhu1=1.0d0+vpar(9)*exphu2
        hulpu2=1.0d0/vuhu1

        exphu3=-exp(vpar(7)*vov1)
        hulpu3=-(1.0d0+exphu3)

        dise2=valp1(ity1)
        exphuu=exp(-vovun(ity1)*vov1)
        hulpu=1.0d0/(1.0d0+exphuu)
        eahu=dise2*hulpu*hulpu2*hulpu3
        estrain(i1)=estrain(i1)+eahu
        deaudvov1=dise2*hulpu2*vovun(ity1)*hulpu*hulpu*exphuu*hulpu3- &
        dise2*hulpu*hulpu2*vpar(7)*exphu3
        ea=ea+eahu
	deaudsumov=-dise2*hulpu*vpar(9)*vpar(10)*hulpu3*exphu2* &
        hulpu2*hulpu2

    !*********************************************************************
    !                                                                    *
    !     Calculate first derivative of atom energy to cartesian         *
    !     coordinates                                                    *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then

            do i3=1,ia(i1,2)
                iob=ia(i1,2+i3)
                ncubo=nubon2(i1,i3)
                do i4=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i4)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+deaudvov1*(1.0d0+dfvl*vho*dvlpdsbo(i1))* &
                        dbondc(ncubo,k1,i4)
                    end do
                end do
            end do
             
            do i2=1,ia(i1,2)
                 
                iat2=ia(i1,2+i2)
                ity2=ia(iat2,1)
                nbosa=nubon2(i1,i2)
                deadvov2=(deaudsumov+dvov1dsumov*deaudvov1)* &
                (bopi(nbosa)+bopi2(nbosa))
                 
                voptlp2=0.50d0*(stlp(ity2)-aval(ity2))
                diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
                deadpibo1=(dvov1dsumov*deaudvov1+deaudsumov)* &
                (abo(iat2)-aval(ity2)-diffvlp2)
                 
                do i4=1,idbo1(nbosa)
                    ihu=idbo(nbosa,i4)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+deadpibo1* &
                        (dbopindc(nbosa,k1,i4)+dbopi2ndc(nbosa,k1,i4))
                    end do
                end do
                 
                do i3=1,ia(iat2,2)
                    iob=ia(iat2,2+i3)
                    ncubo=nubon2(iat2,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        do k1=1,3
                            d(k1,ihu)=d(k1,ihu)+deadvov2*(1.0d0+dfvl*dvlpdsbo(iat2))* &
                            dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do
                 
            end do

        else

            do i3=1,ia(i1,2)
                iob=ia(i1,2+i3)
                ncubo=nubon2(i1,i3)
                do i4=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i4)
                    ix=nmpx(i1,ihu)
                    iy=nmpy(i1,ihu)
                    iz=nmpz(i1,ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+deaudvov1* &
                        (1.0d0+dfvl*vho*dvlpdsbo(i1))*dbondc(ncubo,k1,i4)
                    end do
                end do
            end do

            do i2=1,ia(i1,2)

                iat2=ia(i1,2+i2)
                ity2=ia(iat2,1)
                nbosa=nubon2(i1,i2)
                deadvov2=(deaudsumov+dvov1dsumov*deaudvov1)* &
                (bopi(nbosa)+bopi2(nbosa))

                voptlp2=0.50d0*(stlp(ity2)-aval(ity2))
                diffvlp2=dfvl*(voptlp2-vlptemp(iat2))
                deadpibo1=(dvov1dsumov*deaudvov1+deaudsumov)* &
                (abo(iat2)-aval(ity2)-diffvlp2)

                do i4=1,idbo1(nbosa)
                    ihu=idbo(nbosa,i4)
                    ix=nmpx(i1,ihu)
                    iy=nmpy(i1,ihu)
                    iz=nmpz(i1,ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+deadpibo1* &
                        (dbopindc(nbosa,k1,i4)+dbopi2ndc(nbosa,k1,i4))
                    end do
                end do

                do i3=1,ia(iat2,2)
                    iob=ia(iat2,2+i3)
                    ncubo=nubon2(iat2,i3)
                    do i4=1,idbo1(ncubo)
                        ihu=idbo(ncubo,i4)
                        ix=nmpx(i1,ihu)
                        iy=nmpy(i1,ihu)
                        iz=nmpz(i1,ihu)
                        kcell=14+ix+3*iy+9*iz
                        do k1=1,3
                            dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+deadvov2* &
                            (1.0d0+dfvl*dvlpdsbo(iat2))*dbondc(ncubo,k1,i4)
                        end do
                    end do
                end do

            end do

        end if

    30 END DO
     
!*********************************************************************
!                                                                    *
!     Calculate correction for MgH                                   *
!                                                                    *
!*********************************************************************
!     emgh=zero
!     do 39 i1=1,na
!     ity1=ia(i1,1)
!     vov4=abo(i1)-aval(ity1)

!     do i2=1,ia(i1,2)
!     iat2=ia(i1,2+i2)
!     nbohu=nubon2(i1,i2)

!     ibt=ib(nbohu,1)
!     vc2=vuncor(ibt)
!     elph=zero
!     deahu2dbo=zero
!     deahu2dsbo=zero
!     vov3=bo(nbohu)-vov4-vpar(14)*(vov4**4)
!     if (vov3.gt.1.0) then
!     elph=vc2*(vov3-1.0)*(vov3-1.0)
!     estrain(i1)=estrain(i1)+elph
!     deahu2dbo=2.0*vc2*(vov3-1.0)
!     deahu2dsbo=2.0*vc2*(vov3-1.0)*(-1.0-
!    $4.0*vpar(14)*(vov4**3))
!     end if

!     emgh=emgh+elph

!     do i3=1,idbo1(nbohu)
!     ihu=idbo(nbohu,i3)
!     do k1=1,3
!     d(k1,ihu)=d(k1,ihu)+deahu2dbo*dbondc(nbohu,k1,i3)
!     end do
!     end do

!     do i3=1,ia(i1,2)
!     iob=ia(i1,2+i3)
!     ncubo=nubon2(i1,i3)
!     do i4=1,idbo1(ncubo)
!     ihu=idbo(ncubo,i4)
!     do k1=1,3
!     d(k1,ihu)=d(k1,ihu)+deahu2dsbo*dbondc(ncubo,k1,i4)
!     end do
!     end do
!     end do

!     end do

!  39 continue
!     elp=elp+emgh
!*********************************************************************
!                                                                    *
!     Calculate correction for C2                                    *
!                                                                    *
!*********************************************************************

! df398 always perform correction
    if (abs(vpar(6)) > 0.0001d0) then
        do 40 i1=1,na
            ity1=ia(i1,1)
            vov4=abo(i1)-aval(ity1)
             
            do i2=1,ia(i1,2)


                iat2=ia(i1,2+i2)
                if (qa(i1) == 'C ' .AND. qa(iat2) == 'C ') then
                    nbohu=nubon2(i1,i2)
                     
                    ibt=ib(nbohu,1)
                    elph=zero
                    deahu2dbo=zero
                    deahu2dsbo=zero
                    vov3=bo(nbohu)-vov4-0.040d0*(vov4**4)
                    if (vov3 > 3.0d0) then
                        elph=vpar(6)*(vov3-3.0d0)*(vov3-3.0d0)
                        estrain(i1)=estrain(i1)+elph
                        deahu2dbo=2.0d0*vpar(6)*(vov3-3.0d0)
                        deahu2dsbo=2.0d0*vpar(6)*(vov3-3.0d0)*(-1.0d0- &
                        & 0.16d0*(vov4**3))
                    end if
                     
                    elp=elp+elph
                     
                    if (icpres == 0) then
                        do i3=1,idbo1(nbohu)
                            ihu=idbo(nbohu,i3)
                            do k1=1,3
                                d(k1,ihu)=d(k1,ihu)+deahu2dbo*dbondc(nbohu,k1,i3)
                            end do
                        end do
                         
                        do i3=1,ia(i1,2)
                            iob=ia(i1,2+i3)
                            ncubo=nubon2(i1,i3)
                            do i4=1,idbo1(ncubo)
                                ihu=idbo(ncubo,i4)
                                do k1=1,3
                                    d(k1,ihu)=d(k1,ihu)+deahu2dsbo*dbondc(ncubo,k1,i4)
                                end do
                            end do
                        end do

                    else

                        do i3=1,idbo1(nbohu)
                            ihu=idbo(nbohu,i3)
                            ix=nmpx(i1,ihu)
                            iy=nmpy(i1,ihu)
                            iz=nmpz(i1,ihu)
                            kcell=14+ix+3*iy+9*iz
                            do k1=1,3
                                dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                                deahu2dbo*dbondc(nbohu,k1,i3)
                            end do
                        end do
                         
                        do i3=1,ia(i1,2)
                            iob=ia(i1,2+i3)
                            ncubo=nubon2(i1,i3)
                            do i4=1,idbo1(ncubo)
                                ihu=idbo(ncubo,i4)
                                ix=nmpx(i1,ihu)
                                iy=nmpy(i1,ihu)
                                iz=nmpz(i1,ihu)
                                kcell=14+ix+3*iy+9*iz
                                do k1=1,3
                                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                                    deahu2dsbo*dbondc(ncubo,k1,i4)
                                end do
                            end do
                        end do

                    end if

                end if
            end do

             
        40 END DO
! df398 close IF for correction detection
    end if

    return
    end subroutine ovcor
!*********************************************************************
!*********************************************************************

    subroutine molen

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Calculate molecular energy and first derivatives               *
!     Only used to prevent creating virtual electrons                *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In molen'
        call timer(65)
        close (65)
    end if
    emol=zero
    return
    do i1=1,nmolo

        enelm=0.0d0
        do i2=1,na
            if (ia(i2,3+mbond) == i1) then
                it1=ia(i2,1)
                enelm=enelm+aval(it1)
            end if
        end do

        na1m=nmolat(i1,1)

        enelm=2*int(enelm*0.50d0)
    !     enelm=elmol(i1)
        bomsum=zero
        do i2=1,na1m
            ihu=nmolat(i1,i2+1)
            do i3=1,ia(ihu,2)
                ihu2=nubon2(ihu,i3)
                bomsum=bomsum+bo(ihu2)
            end do
        end do
        diff=(bomsum-enelm)
        exphu=exp(-vpar(37)*diff)
        exphu2=1.0d0/(1.0d0+15.0d0*exphu)
        emolh=zero
        demoldsbo=zero
        emolh=vpar(38)*exphu2
        emol=emol+emolh
        demoldsbo=vpar(38)*vpar(37)*15.0d0*exphu2*exphu2*exphu

        do i2=1,na1m
            ihu1=nmolat(i1,i2+1)
            do i3=1,ia(ihu1,2)
                iob=ia(ihu1,2+i3)
                ncubo=nubon2(ihu1,i3)
                do i4=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i4)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+demoldsbo*dbondc(ncubo,k1,i4)
                    end do
                end do
            end do
        end do


    end do


    return
    end subroutine molen
!*********************************************************************
!*********************************************************************

    subroutine valang

!*********************************************************************
    include 'cbka.blk'
    dimension j(3)
!*********************************************************************
!                                                                    *
!     Calculate valency angle energies and first derivatives         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In valang'
        call timer(65)
        close (65)
    end if
!     eco=0.0
    ev=0.0d0
    ecoa=0.0d0
    epen=0.0d0
    if (nval == 0) return

    do 10 i1=1,nval
        ity=iv(i1,1)
        j(1)=iv(i1,2)
        j(2)=iv(i1,3)
        j(3)=iv(i1,4)

        la=iv(i1,5)
        lb=iv(i1,6)
        boa=bo(la)-cutof2
        bob=bo(lb)-cutof2
        if (boa < zero .OR. bob < zero) goto 10

        hl=h(i1)     ! Calculated earlier in routine calval
    !*********************************************************************
    !                                                                    *
    !     Calculate valency angle energy                                 *
    !                                                                    *
    !*********************************************************************
        nbocen=ia(j(2),2)
        sbo2=0.0d0
        vmbo=1.0d0
             
        do i2=1,nbocen
            ibv=nubon2(j(2),i2)
            vmbo=vmbo*exp(-bo(ibv)**8)
            sbo2=sbo2+bopi(ibv)+bopi2(ibv)
        end do

        ity2=ia(j(2),1)
    !     exbo=abo(j(2))-stlp(ia(j(2),1))
        exbo=abo(j(2))-valf(ity2)
    !     if (exbo.gt.zero) exbo=zero
    !     expov=exp(vka8(ity)*exbo)
    !     expov2=exp(-vpar(13)*exbo)
    !     htov1=2.0+expov2
    !     htov2=1.0+expov+expov2
    !     evboadj=htov1/htov2
        evboadj=1.0d0
        expun=exp(-vkac(ity)*exbo)
        expun2=exp(vpar(15)*exbo)
        htun1=2.0d0+expun2
        htun2=1.0d0+expun+expun2
        evboadj2=vval4(ity2)-(vval4(ity2)-1.0d0)*htun1/htun2
    !*********************************************************************
    !                                                                    *
    !     Calculate number of lone pairs                                 *
    !                                                                    *
    !*********************************************************************
        dsbo2dvlp=(1.0d0-vmbo)
        vlpadj=zero
        exlp1=abo(j(2))-stlp(ia(j(2),1))
        exlp2=2.0d0*int(exlp1/2.0d0)
        exlp=exlp1-exlp2
        if (exlp < zero) then
        !     expvlp=exp(-vpar(16)*(2.0+exlp)*(2.0+exlp))
        !     vlpadj=expvlp-int(exlp1/2.0)
        !     dsbo2dvlp=(1.0-vmbo)*(1.0-vpar(34)*
        !    $2.0*(2.0+exlp)*vpar(16)*expvlp)
            vlpadj=vlp(j(2))
            dsbo2dvlp=(1.0d0-vmbo)*(1.0d0+vpar(34)*dvlpdsbo(j(2)))
        end if

        sbo2=sbo2+(1.0d0-vmbo)*(-exbo-vpar(34)*vlpadj)
        dsbo2dvmbo=exbo+vpar(34)*vlpadj

        sbo2h=sbo2
        powv=vpar(17)
        if (sbo2 <= 0.0d0) sbo2h=0.0d0
        if (sbo2 > 0.0d0 .AND. sbo2 <= 1.0d0) sbo2h=sbo2**powv
        if (sbo2 > 1.0d0 .AND. sbo2 < 2.0d0) sbo2h=2.0d0-(2.0d0-sbo2)**powv
        if (sbo2 > 2.0d0) sbo2h=2.0d0
        thba=th0(ity)
        expsbo=exp(-vpar(18)*(2.0d0-sbo2h))
        thetao=180.0d0-thba*(1.0d0-expsbo)

        thetao=thetao*dgrrdn
        thdif=(thetao-hl)
        thdi2=thdif*thdif
        dthsbo=dgrrdn*thba*vpar(18)*expsbo
        if (sbo2 < 0.0d0) dthsbo=zero
        if (sbo2 > 0.0d0 .AND. sbo2 <= 1.0d0) &
        dthsbo=powv*(sbo2**(powv-1.0d0))*dgrrdn*thba*vpar(18)*expsbo
        if (sbo2 > 1.0d0 .AND. sbo2 < 2.0d0) &
        dthsbo=powv*((2.0d0-sbo2)**(powv-1.0d0))*dgrrdn*thba*vpar(18)*expsbo
        if (sbo2 > 2.0d0) dthsbo=zero

        exphu=vka(ity)*exp(-vka3(ity)*thdi2)
        exphu2=vka(ity)-exphu
        if (vka(ity) < zero) exphu2=exphu2-vka(ity)             !To avoid linear Me-H-Me angles (6/6/06)
        boap=boa**vval2(ity)
        boap2=boa**(vval2(ity)-1.0d0)
        bobp=bob**vval2(ity)
        bobp2=bob**(vval2(ity)-1.0d0)
        exa=exp(-vval1(ity2)*boap)
        exb=exp(-vval1(ity2)*bobp)
        dexadboa=vval2(ity)*vval1(ity2)*exa*boap2
        dexbdbob=vval2(ity)*vval1(ity2)*exb*bobp2
        exa2=(1.0d0-exa)
        exb2=(1.0d0-exb)

	! df398 calculate taper functions        
	!call ataper(.false.,bo(la),cutof2,4.0d0*cutof2,fij,dfijdBOij,d2fijdBO2ij,d3fijdBO3ij,.true.,.false.,.false.)
	!call ataper(.false.,bo(lb),cutof2,4.0d0*cutof2,fik,dfikdBOik,d2fikdBO2ik,d3fikdBO3ik,.true.,.false.,.false.)
	call valtaper(.false.,bo(la),cutof2,4.0D0*cutof2,fij,dfijdBOij)
	call valtaper(.false.,bo(lb),cutof2,4.0D0*cutof2,fik,dfikdBOik)

        evh=evboadj2*evboadj*exa2*exb2*exphu2*fij*fik
        estrain(j(2))=estrain(j(2))+evh

	! df398 multiply f7 and f8 derivatives by taper + add new derivative term for fij and fik
        devdlb=fij*fik*evboadj2*evboadj*dexbdbob*exa2*exphu2 + dfikdBOik*fij*evboadj2*evboadj*exb2*exa2*exphu2
        devdla=fij*fik*evboadj2*evboadj*dexadboa*exb2*exphu2 + fik*dfijdBOij*evboadj2*evboadj*exa2*exb2*exphu2

        devdsbo=2.0d0*evboadj2*evboadj*dthsbo*exa2*exb2*vka3(ity)*thdif*exphu*fij*fik

        devdh=-2.0d0*evboadj2*evboadj*exa2*exb2*vka3(ity)*thdif*exphu*fij*fik

        devdsbo2= evboadj*exa2*exb2*exphu2*(vval4(ity2)-1.0d0)*(-vpar(15)*expun2/htun2+htun1*&
                  (vpar(15)*expun2-vkac(ity)*expun)/(htun2*htun2))*fij*fik

    !     devdsbo2=-evboadj2*exa2*exb2*exphu2*(vpar(13)*expov2/htov2+
    !    $htov1*(vka8(ity)*expov-vpar(13)*expov2)/(htov2*htov2))+
    !    $evboadj*exa2*exb2*exphu2*(vpar(14)-1.0)*(-vpar(15)*expun2/htun2
    !    $+htun1*(vpar(15)*expun2-vkac(ity)*expun)/(htun2*htun2))



        ev=ev+evh
    !     write (64,'(4i8,18f8.2)')mdstep,j(1),j(2),j(3),sbo2,sbo2h,
    !    $thetao*rdndgr,hl*rdndgr,bo(la),bo(lb),bopi(la),
    !    $vlp(j(2)),exbo,vlpadj,vmbo,evh,ev,vka(ity)
    !*********************************************************************
    !                                                                    *
    !     Calculate penalty for two double bonds in valency angle        *
    !                                                                    *
    !*********************************************************************
        exbo=abo(j(2))-aval(ia(j(2),1))
        expov=exp(vpar(22)*exbo)
        expov2=exp(-vpar(21)*exbo)
        htov1=2.0d0+expov2
        htov2=1.0d0+expov+expov2
        ecsboadj=htov1/htov2
        exphu1=exp(-vpar(20)*(boa-2.0d0)*(boa-2.0d0))
        exphu2=exp(-vpar(20)*(bob-2.0d0)*(bob-2.0d0))

        epenh=vkap(ity)*ecsboadj*exphu1*exphu2
	! df398 multiply by taper functions fij, fik
        estrain(j(2))=estrain(j(2))+fij*fik*epenh
        epen=epen+fij*fik*epenh
	! df398 multiply by taper functions fij, fik and add derivative for fij, fik
        decoadboa=-2.0d0*vpar(20)*epenh*(boa-2.0d0)*fij*fik + dfijdBOij*fik*epenh
        decoadbob=-2.0d0*vpar(20)*epenh*(bob-2.0d0)*fij*fik + dfikdBOik*fij*epenh
        decdsbo2=-vkap(ity)*exphu1*exphu2*(vpar(21)*expov2/htov2+htov1* &
        (vpar(22)*expov-vpar(21)*expov2)/(htov2*htov2))*fij*fik

    !     write (65,'(4i3,20f12.4)')mdstep,j(1),j(2),j(3),exbo,expov,
    !    $expov2,htov1,htov2,ecsboadj,exphu1,exphu2,epenh,decoadboa,
    !    $dcoadbob,decdsbo2
    !*********************************************************************
    !                                                                    *
    !     Calculate valency angle conjugation energy                     *
    !                                                                    *
    !*********************************************************************
        unda=abo(j(1))-boa
    !     ovb=abo(j(2))-valf(ia(j(2),1))
        ovb=abo(j(2))-vval3(ia(j(2),1))    !Modification for Ru  7/6/2004
        undc=abo(j(3))-bob
        ba=(boa-1.50d0)*(boa-1.50d0)
        bb=(bob-1.50d0)*(bob-1.50d0)
        exphua=exp(-vpar(31)*ba)
        exphub=exp(-vpar(31)*bb)
        exphuua=exp(-vpar(39)*unda*unda)
        exphuob=exp(vpar(3)*ovb)
        exphuuc=exp(-vpar(39)*undc*undc)
        hulpob=1.0d0/(1.0d0+exphuob)

        ecoah=vka8(ity)*exphua*exphub*exphuua*exphuuc*hulpob
	! df398 multiply by taper funcs fij, fik
        estrain(j(2))=estrain(j(2))+ecoah*fij*fik

        ! df398 multiply all derivatives by taper functions fij, fik and add new derivatives
        decodbola=-2.0d0*vka8(ity)*(boa-1.50d0)*vpar(31)*exphua*exphub*exphuua*exphuuc*hulpob*fij*fik &
          + vpar(39)*vka8(ity)*exphua*exphub*exphuua*exphuuc*hulpob*2.0d0*unda*fij*fik + dfijdBOij*fik*ecoah

        decodbolb=-2.0d0*vka8(ity)*(bob-1.50d0)*vpar(31)*exphua*exphub*exphuua*exphuuc*hulpob*fij*fik &
          + vpar(39)*vka8(ity)*exphua*exphub*exphuua*exphuuc*hulpob*2.0d0*undc*fij*fik + dfikdBOik*fij*ecoah
        ! df398 multiply by taper funcs fij, fik
        decodboua=-2.0d0*unda*vka8(ity)*vpar(39)*exphua*exphub*exphuua*exphuuc*hulpob*fij*fik

        decodbouc=-2.0d0*undc*vka8(ity)*vpar(39)*exphua*exphub*exphuua*exphuuc*hulpob*fij*fik
        ! df398 decodboob = devalddeltai in GULP reaxffmd
        decodboob=-vka8(ity)*exphua*exphub*exphuua*exphuuc*hulpob*hulpob*vpar(3)*exphuob*fij*fik

	! df398 multiply by taper funcs fij, fik
        ecoa=ecoa+ecoah*fij*fik
    !     write (64,'(4i3,20f12.4)')mdstep,j(1),j(2),j(3),vka8(ity),
    !    $ba,bb,ovb,exphua,exphub,exphuua,exphuob,exphuuc,hulpob,ecoah
    !*********************************************************************
    !                                                                    *
    !     Calculate derivative valency energy to cartesian coordinates   *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then

            do k1=1,3
                do k2=1,3
                    d(k1,j(k2))=d(k1,j(k2))+devdh*dhdc(i1,k1,k2)
                end do
            end do
             
            do i2=1,idbo1(la)
                ihu=idbo(la,i2)
                do k1=1,3
                    d(k1,ihu)=d(k1,ihu)+(devdla+decoadboa+decodbola)* &
                    dbondc(la,k1,i2)
                end do
            end do
             
            do i2=1,idbo1(lb)
                ihu=idbo(lb,i2)
                do k1=1,3
                    d(k1,ihu)=d(k1,ihu)+(devdlb+decoadbob+decodbolb)* &
                    dbondc(lb,k1,i2)
                end do
            end do
             
            do i2=1,nbocen
                j5=ia(j(2),2+i2)
                ibv=nubon2(j(2),i2)
                dvmbodbo=-vmbo*8.0d0*bo(ibv)**7
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+(-dsbo2dvlp*devdsbo+devdsbo2+decdsbo2 &
                        +dvmbodbo*dsbo2dvmbo*devdsbo)* &
                        dbondc(ibv,k1,i3)+devdsbo*(dbopindc(ibv,k1,i3)+ &
                        dbopi2ndc(ibv,k1,i3))
                    end do
                end do
            end do
             
            do i2=1,ia(j(1),2)
                j5=ia(j(1),2+i2)
                ibv=nubon2(j(1),i2)
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+decodboua*dbondc(ibv,k1,i3)
                    end do
                end do
            end do
             
            do i2=1,ia(j(2),2)
                j5=ia(j(2),2+i2)
                ibv=nubon2(j(2),i2)
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+decodboob*dbondc(ibv,k1,i3)
                    end do
                end do
            end do

            do i2=1,ia(j(3),2)
                j5=ia(j(3),2+i2)
                ibv=nubon2(j(3),i2)
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+decodbouc*dbondc(ibv,k1,i3)
                    end do
                end do
            end do

        else

            do k1=1,3
                do k2=1,3
                    ix=nmpx(j(2),j(k2))
                    iy=nmpy(j(2),j(k2))
                    iz=nmpz(j(2),j(k2))
                    kcell=14+ix+3*iy+9*iz
                    dcell(k1,j(k2),kcell)=dcell(k1,j(k2),kcell)+devdh*dhdc(i1,k1,k2)
                end do
            end do

            do i2=1,idbo1(la)
                ihu=idbo(la,i2)
                ix=nmpx(j(2),ihu)
                iy=nmpy(j(2),ihu)
                iz=nmpz(j(2),ihu)
                kcell=14+ix+3*iy+9*iz
                do k1=1,3
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                    (devdla+decoadboa+decodbola)*dbondc(la,k1,i2)
                end do
            end do

            do i2=1,idbo1(lb)
                ihu=idbo(lb,i2)
                ix=nmpx(j(2),ihu)
                iy=nmpy(j(2),ihu)
                iz=nmpz(j(2),ihu)
                kcell=14+ix+3*iy+9*iz
                do k1=1,3
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                    (devdlb+decoadbob+decodbolb)*dbondc(lb,k1,i2)
                end do
            end do

            do i2=1,nbocen
                j5=ia(j(2),2+i2)
                ibv=nubon2(j(2),i2)
                dvmbodbo=-vmbo*8.0d0*bo(ibv)**7
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    ix=nmpx(j(2),ihu)
                    iy=nmpy(j(2),ihu)
                    iz=nmpz(j(2),ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        (-dsbo2dvlp*devdsbo+devdsbo2+decdsbo2+ &
                        dvmbodbo*dsbo2dvmbo*devdsbo)* &
                        dbondc(ibv,k1,i3)+devdsbo*(dbopindc(ibv,k1,i3)+ &
                        dbopi2ndc(ibv,k1,i3))
                    end do
                end do
            end do

            do i2=1,ia(j(1),2)
                j5=ia(j(1),2+i2)
                ibv=nubon2(j(1),i2)
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    ix=nmpx(j(2),ihu)
                    iy=nmpy(j(2),ihu)
                    iz=nmpz(j(2),ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        decodboua*dbondc(ibv,k1,i3)
                    end do
                end do
            end do

            do i2=1,ia(j(2),2)
                j5=ia(j(2),2+i2)
                ibv=nubon2(j(2),i2)
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    ix=nmpx(j(2),ihu)
                    iy=nmpy(j(2),ihu)
                    iz=nmpz(j(2),ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        decodboob*dbondc(ibv,k1,i3)
                    end do
                end do
            end do

            do i2=1,ia(j(3),2)
                j5=ia(j(3),2+i2)
                ibv=nubon2(j(3),i2)
                do i3=1,idbo1(ibv)
                    ihu=idbo(ibv,i3)
                    ix=nmpx(j(2),ihu)
                    iy=nmpy(j(2),ihu)
                    iz=nmpz(j(2),ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        decodbouc*dbondc(ibv,k1,i3)
                    end do
                end do
            end do

        end if

    10 END DO

    return
    end subroutine valang
!*********************************************************************
!*********************************************************************

    subroutine hbond

!*********************************************************************
    include 'cbka.blk'
    dimension drda(3),j(3),dvdc(3,3),dargdc(3,3)
    virial=zero
    virx=zero
    viry=zero
    virz=zero
!*********************************************************************
!                                                                    *
!     Calculate hydrogen bond energies and first derivatives         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In hbond'
        call timer(65)
        close (65)
    end if
    ehb=zero
    do 10 i1=1,nhb
        ityhb=ihb(i1,1)
        j(1)=ihb(i1,2)
        j(2)=ihb(i1,3)
        j(3)=ihb(i1,4)
        la=ihb(i1,5)
        ix=ihb(i1,6)
        iy=ihb(i1,7)
        iz=ihb(i1,8)
        dxm=c(j(2),1)-c(j(3),1)+ix*tm11
        dym=c(j(2),2)-c(j(3),2)+ix*tm21+iy*tm22
        dzm=c(j(2),3)-c(j(3),3)+ix*tm31+iy*tm32+iz*tm33
        boa=bo(la)
        rda=sqrt(dxm*dxm+dym*dym+dzm*dzm)

	!call ataper(.false.,boa,0.01D0,4.0d0*0.01D0,fHB,dfHBdBO,d2fHBdBO2,d3fHBdBO3,.true.,.false.,.false.)
	!call ataper(.true.,rda,0.8D0*7.5D0,7.5D0,frHB,dfrHBdr,d2frHBdr2,d3frHBdr3,.true.,.false.,.false.)
	call valtaper(.false.,boa,0.01D0,4.0D0*0.01D0,fHB,dfHBdBO)
	call valtaper(.true.,rda,0.9D0*7.5D0,7.5D0,frHB,dfrHBdr)

        drda(1)=dxm/rda
        drda(2)=dym/rda
        drda(3)=dzm/rda
        call calvalhb(j(1),j(2),j(3),ix,iy,iz,arg,hhb(i1),dvdc,dargdc)
        rhu1=rhb(ityhb)/rda
        rhu2=rda/rhb(ityhb)
        sinhu=sin(hhb(i1)/2.0d0)
        sin2=sinhu*sinhu
        exphu1=exp(-vhb1(ityhb)*boa)
        exphu2=exp(-vhb2(ityhb)*(rhu1+rhu2-2.0d0))

        ! ehbh=(1.0d0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2*fHB*frHB   !!! df398 old sin^4 form
        ehbh=(1.0d0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2*sin2*sin2*fHB*frHB  !!! new sin^8 form

        estrain(j(2))=estrain(j(2))+0.50d0*ehbh
        estrain(j(3))=estrain(j(3))+0.50d0*ehbh

        ehb=ehb+ehbh
    !     write (63,'(9i4,10f12.4)')i1,ityhb,j(1),j(2),j(3),ix,iy,iz,la,
    !    $boa,rda,rdndgr*hhb(i1),1.0-exphu1,exphu2,sin2*sin2,
    !    $ehbh,ehb
    !*********************************************************************
    !                                                                    *
    !     Calculate first derivatives                                    *
    !                                                                    *
    !*********************************************************************

        dehbdbo  = fHB*frHB*vhb1(ityhb)*exphu1*dehb(ityhb)*exphu2*sin2*sin2*sin2*sin2 &    !!!  df398 derivs for new sin^8 form
                   + dfHBdBO*(1.0d0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2*sin2*sin2*frHB

        dehbdv   = fHB*frHB*(1.0d0-exphu1)*dehb(ityhb)*exphu2*4.0d0*sin2*sin2*sin2*sinhu*cos(hhb(i1)/2.0d0) 

        dehbdrda = fHB*frHB*(1.0d0-exphu1)*dehb(ityhb)*sin2*sin2*sin2*sin2* &
                   vhb2(ityhb)*(rhb(ityhb)/(rda*rda)-1.0d0/rhb(ityhb))*exphu2 &
                   + dfrHBdr*(1.0d0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2*sin2*sin2*fHB

        !dehbdbo  = fHB*frHB*vhb1(ityhb)*exphu1*dehb(ityhb)*exphu2*sin2*sin2 &          !!! df398 these derivatives belong to sin^4
        !           + dfHBdBO*(1.0d0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2*frHB

        !dehbdv   = fHB*frHB*(1.0d0-exphu1)*dehb(ityhb)*exphu2*2.0d0*sin2*sinhu*cos(hhb(i1)/2.0d0)                          

        !dehbdrda = fHB*frHB*(1.0d0-exphu1)*dehb(ityhb)*sin2*sin2*vhb2(ityhb)* &
        !           (rhb(ityhb)/(rda*rda)- 1.0d0/rhb(ityhb))*exphu2 &
        !           + dfrHBdr*(1.0d0-exphu1)*dehb(ityhb)*exphu2*sin2*sin2*fHB

        if (icpres == 0) then

            do k1=1,3
                d(k1,j(2))=d(k1,j(2))+dehbdrda*drda(k1)
                d(k1,j(3))=d(k1,j(3))-dehbdrda*drda(k1)
            end do
             
            do k1=1,3
                do k2=1,3
                    d(k1,j(k2))=d(k1,j(k2))+dehbdv*dvdc(k1,k2)
                end do
            end do
             
            do i2=1,idbo1(la)
                ihu=idbo(la,i2)
                do k1=1,3
                    d(k1,ihu)=d(k1,ihu)+dehbdbo*dbondc(la,k1,i2)
                end do
            end do

        else

            ix=nmpx(j(2),j(3))
            iy=nmpy(j(2),j(3))
            iz=nmpz(j(2),j(3))
            kcell=14+ix+3*iy+9*iz
            do k1=1,3
                dcell(k1,j(2),14)=dcell(k1,j(2),14)+dehbdrda*drda(k1)
                dcell(k1,j(3),kcell)=dcell(k1,j(3),kcell)-dehbdrda*drda(k1)
            end do
                 
            do k1=1,3
                do k2=1,3
                    ix=nmpx(j(2),j(k2))
                    iy=nmpy(j(2),j(k2))
                    iz=nmpz(j(2),j(k2))
                    kcell=14+ix+3*iy+9*iz
                    dcell(k1,j(k2),kcell)=dcell(k1,j(k2),kcell)+dehbdv*dvdc(k1,k2)
                end do
            end do
              
            do i2=1,idbo1(la)
                ihu=idbo(la,i2)
                ix=nmpx(j(2),ihu)
                iy=nmpy(j(2),ihu)
                iz=nmpz(j(2),ihu)
                kcell=14+ix+3*iy+9*iz
                do k1=1,3
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                    dehbdbo*dbondc(la,k1,i2)
                end do
            end do

        end if

           
    10 END DO
    return
    end subroutine hbond

!*********************************************************************
!*********************************************************************

    subroutine torang

!*********************************************************************
    include 'cbka.blk'
    DIMENSION  A(3),DRDA(3),DADC(4),DRADC(3,4),DRBDC(3,4), &
    DRCDC(3,4),DHDDC(3,4),DHEDC(3,4),DRVDC(3,4),DTDC(3,4), &
    DNDC(3,4)
    dimension j(4),dh1rdc(3,3),dh2rdc(3,3),dargdc(3,3)
!*********************************************************************
!                                                                    *
!     Calculate torsion angle energies and first derivatives         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In torang'
        call timer(65)
        close (65)
    end if
    do k1=1,3
        do k2=1,4
            dhddc(k1,k2)=0.0d0
            dhedc(k1,k2)=0.0d0
            dradc(k1,k2)=0.0d0
            drbdc(k1,k2)=0.0d0
            drcdc(k1,k2)=0.0d0
        end do
    end do
    et=0.0d0
    eth12=0.0d0
    eco=0.0d0
    dadc(1)=1.0d0
    dadc(2)=0.0d0
    dadc(3)=0.0d0
    dadc(4)=-1.0d0
    if (ntor == 0) return

    do 10 i1=1,ntor
        j(1)=it(i1,2)
        j(2)=it(i1,3)
        j(3)=it(i1,4)
        j(4)=it(i1,5)
        ity=it(i1,1)
        la=it(i1,6)
        lb=it(i1,7)
        lc=it(i1,8)
        call calvalres(j(1),j(2),j(3),arg1,ht1,dh1rdc,dargdc)
        call calvalres(j(2),j(3),j(4),arg2,ht2,dh2rdc,dargdc)
        boa=bo(la)-cutof2
        bob=bo(lb)-cutof2
        boc=bo(lc)-cutof2
        if (boa < zero .OR. bob < zero .OR. boc < zero) goto 10

	!call ataper(.false.,bo(la),cutof2,4.0d0*cutof2,fij,dfijdBOij,d2fijdBO2ij,d3fijdBO3ij,.true.,.false.,.false.)
	!call ataper(.false.,bo(lb),cutof2,4.0d0*cutof2,fik,dfikdBOik,d2fikdBO2ik,d3fikdBO3ik,.true.,.false.,.false.)
        !call ataper(.false.,bo(lc),cutof2,4.0d0*cutof2,fjl,dfjldBOjl,d2fjldBO2jl,d3fjldBO3jl,.true.,.false.,.false.)
        call valtaper(.false.,bo(la),cutof2,4.0D0*cutof2,fij,dfijdBOij)
	call valtaper(.false.,bo(lb),cutof2,4.0D0*cutof2,fik,dfikdBOik)
	call valtaper(.false.,bo(lc),cutof2,4.0D0*cutof2,fjl,dfjldBOjl)
        
       !df398 skip applying taper
        !fij = 1.0d0
        !fik = 1.0d0
        !fjl = 1.0d0
        !dfijdBOij = 0.0d0
        !dfikdBOik = 0.0d0
        !dfjldBOjl = 0.0d0

        r42=0.0d0
        ivl1=ibsym(la)
        ivl2=ibsym(lb)
        ivl3=ibsym(lc)
        isign1=1
        isign2=1
        isign3=1
        if (j(2) < j(1)) isign1=-1.0D0
        if (j(3) < j(2)) isign2=-1.0D0
        if (j(4) < j(3)) isign3=-1.0D0
        rla=rbo(la)
        rlb=rbo(lb)
        ix1=isign1*nvlx(ivl1)+isign2*nvlx(ivl2)+isign3*nvlx(ivl3)
        iy1=isign1*nvly(ivl1)+isign2*nvly(ivl2)+isign3*nvly(ivl3)
        iz1=isign1*nvlz(ivl1)+isign2*nvlz(ivl2)+isign3*nvlz(ivl3)
         
        a(1)=c(j(1),1)-c(j(4),1)+ix1*tm11
        a(2)=c(j(1),2)-c(j(4),2)+ix1*tm21+iy1*tm22
        a(3)=c(j(1),3)-c(j(4),3)+ix1*tm31+iy1*tm32+iz1*tm33
        r4=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
    !*********************************************************************
    !                                                                    *
    !     Determine torsion angle                                        *
    !                                                                    *
    !*********************************************************************
        d142=r4*r4
        rla=rbo(la)
        rlb=rbo(lb)
        rlc=rbo(lc)
        coshd=cos(ht1)
        coshe=cos(ht2)
        sinhd=sin(ht1)
        sinhe=sin(ht2)
        poem=2.0d0*rla*rlc*sinhd*sinhe
        poem2=poem*poem
        tel=rla*rla+rlb*rlb+rlc*rlc-d142-2.0d0*(rla*rlb*coshd-rla*rlc* &
            coshd*coshe+rlb*rlc*coshe)
        if (poem < 1.0d-20) poem=1.0d-20
        arg=tel/poem
        if (arg > 1.0d0) arg=1.0d0
        if (arg < -1.0d0) arg=-1.0d0
        arg2=arg*arg
        thg(i1)=acos(arg)*rdndgr
        k1=j(1)
        k2=j(2)
        k3=j(3)
        k4=j(4)
        call dista2(k3,k2,dis,x3,y3,z3)
        y32z32=y3*y3+z3*z3
        wort1=sqrt(y32z32)+1.0D-6
        wort2=sqrt(y32z32+x3*x3)+1.0D-6
       !  if (wort1.lt.1e-6) wort1=1e-6
       !  if (wort2.lt.1e-6) wort2=1e-6
        sinalf=y3/wort1
        cosalf=z3/wort1
        sinbet=x3/wort2
        cosbet=wort1/wort2
        call dista2(k1,k2,dis,x1,y1,z1)
        x1=x1*cosbet-y1*sinalf*sinbet-z1*cosalf*sinbet
        y1=y1*cosalf-z1*sinalf
        wort3=sqrt(x1*x1+y1*y1)+1.0D-6
       !  if (wort3.lt.1e-6) wort3=1e-6
        singam=y1/wort3
        cosgam=x1/wort3
        call dista2(k4,k2,dis,x4,y4,z4)
        x4=x4*cosbet-y4*sinalf*sinbet-z4*cosalf*sinbet
        y4=y4*cosalf-z4*sinalf
        y4=x4*singam-y4*cosgam
        if (y4 > 0.0d0) thg(i1)=-thg(i1)
        if (thg(i1) < -179.9999999999999d0) thg(i1)=-179.9999999999999d0
        if (thg(i1) > 179.9999999999999d0) thg(i1)=179.9999999999999d0
        th2=thg(i1)*dgrrdn

    !*********************************************************************
    !                                                                    *
    !     Calculate torsion angle energy                                 *
    !                                                                    *
    !*********************************************************************

        exbo1=abo(j(2))-valf(ia(j(2),1))
        exbo2=abo(j(3))-valf(ia(j(3),1))
        htovt=exbo1+exbo2
        expov=exp(vpar(26)*htovt)
        expov2=exp(-vpar(25)*(htovt))
        htov1=2.0d0+expov2
        htov2=1.0d0+expov+expov2
        etboadj=htov1/htov2

!etboadj = 0.0d0

     !    btb2=bopi(lb)-1.0+etboadj
     !    bo2t=1.0-btb2
        bo2t=2.0d0-bopi(lb)-etboadj
        bo2p=bo2t*bo2t
        bocor2=exp(v4(ity)*bo2p)

        !df398 This form (original ReaxFF) suffers from numerical instability
        !hsin=sinhd*sinhe
        !df398 Numericaly stable form. Third-power sine decays much faster and ensures torsion energy and force go to zero when valence angles go to zero or pi. 
        hsin=(sinhd**3)*(sinhe**3)


        ethhulp=0.50d0*v1(ity)*(1.0d0+arg)+v2(ity)*bocor2*(1.0d0-arg2)+ &
                v3(ity)*(0.50d0+2.0d0*arg2*arg-1.50d0*arg) ! df398 trigo identity for 0.5*cos(3x)


        exphua=exp(-vpar(24)*boa) !df398 original form
        exphub=exp(-vpar(24)*bob)
        exphuc=exp(-vpar(24)*boc)

        !df398 2013 correction new form for f10
!        exphua=exp(-2.0d0*vpar(24)*boa*boa) 
!        exphub=exp(-2.0d0*vpar(24)*bob*bob)
!        exphuc=exp(-2.0d0*vpar(24)*boc*boc)

        bocor4=(1.0d0-exphua)*(1.0d0-exphub)*(1.0d0-exphuc) ! f10 equation

        eth=hsin*ethhulp*bocor4

        estrain(j(2))=estrain(j(2))+0.50d0*eth*fij*fik*fjl
        estrain(j(3))=estrain(j(3))+0.50d0*eth*fij*fik*fjl

        !df398 multiply by taper funcs, fij, fik, fjl
        detdar=hsin*bocor4*(0.50d0*v1(ity)-2.0d0*v2(ity)*bocor2*arg+v3(ity)*(6.0d0*arg2-1.5d0))*fij*fik*fjl

        ! df398 This form (original ReaxFF) suffers from numerical instability
        !detdhd=coshd*sinhe*bocor4*ethhulp*fij*fik*fjl
        !detdhe=sinhd*coshe*bocor4*ethhulp*fij*fik*fjl
        ! This form is stable because sin^(x) decays much faster than sin(x)
        ! Otherwise force discontinuities appear
        detdhd=3.0d0*sinhd**2*coshd*sinhe**3*bocor4*ethhulp*fij*fik*fjl
        detdhe=sinhd**3*3.0d0*sinhe**2*coshe*bocor4*ethhulp*fij*fik*fjl

        !df398 multiply by taper funcs and add taper derivatives
        detdboa=vpar(24)*exphua*(1.0D0-exphub)*(1.0D0-exphuc)*ethhulp*hsin*fij*fik*fjl &  ! df398 original form
                + eth*fik*fjl*dfijdBOij

        ! df398 2013 correction new form
!        detdboa=4.0d0*boa*vpar(24)*exphua*(1.0d0-exphub)*(1.0d0-exphuc)*ethhulp*hsin  

        detdbopib=-bocor4*2.0d0*v4(ity)*v2(ity)*bo2t*bocor2*(1.0d0-arg2)*hsin*fij*fik*fjl

        detdbob=vpar(24)*exphub*(1.0D0-exphua)*(1.0D0-exphuc)*ethhulp*hsin*fij*fik*fjl &   ! df398 original forms
                + eth*fij*fjl*dfikdBOik

        detdboc=vpar(24)*exphuc*(1.0D0-exphua)*(1.0D0-exphub)*ethhulp*hsin*fij*fik*fjl & ! 
                + eth*fij*fik*dfjldBOjl

         !df398 2013 correction new forms for f10 derivative
!        detdbob=4.0d0*bob*vpar(24)*exphub*(1.0d0-exphua)*(1.0d0-exphuc)*ethhulp*hsin  
!        detdboc=4.0d0*boc*vpar(24)*exphuc*(1.0d0-exphua)*(1.0d0-exphub)*ethhulp*hsin 

        detdsbo1=-(detdbopib)*( vpar(25)*expov2/htov2 + htov1*(vpar(26)*expov-vpar(25)*expov2)/(htov2*htov2) )

!        detdsbo1 = 0.0d0


        et=et+eth*fij*fik*fjl

    !*********************************************************************
    !                                                                    *
    !     Calculate conjugation energy                                   *
    !                                                                    *
    !*********************************************************************
        ba=(boa-1.50D0)*(boa-1.50D0)  ! df398 original form
        bb=(bob-1.50D0)*(bob-1.50D0)  ! 
        bc=(boc-1.50D0)*(boc-1.50D0)  ! 
        exphua1=exp(-vpar(28)*ba)     ! 
        exphub1=exp(-vpar(28)*bb)     ! 
        exphuc1=exp(-vpar(28)*bc)     ! 

        !!!df398 2013 correction new form 
        !!!exphua1=sin(pi/3.0d0*boa)*sin(pi/3.0d0*boa)*sin(pi/3.0d0*boa)*sin(pi/3.0d0*boa)   
        !!!exphub1=sin(pi/3.0d0*bob)*sin(pi/3.0d0*bob)*sin(pi/3.0d0*bob)*sin(pi/3.0d0*bob) 
        !!!exphuc1=sin(pi/3.0d0*boc)*sin(pi/3.0d0*boc)*sin(pi/3.0d0*boc)*sin(pi/3.0d0*boc) 


        sbo=exphua1*exphub1*exphuc1 ! this is equation f12

	!df398 derivatives of sbo wrt BO
        dbohua=-2.0D0*(boa-1.50D0)*vpar(28)*exphua1*exphub1*exphuc1 ! df398 original form
        dbohub=-2.0D0*(bob-1.50D0)*vpar(28)*exphua1*exphub1*exphuc1 ! 
        dbohuc=-2.0D0*(boc-1.50D0)*vpar(28)*exphua1*exphub1*exphuc1 ! 

	!!!df398 2013 correction new form
        !!dbohua=4.0d0/3.0d0*pi*sin(pi/3.0d0*boa)*sin(pi/3.0d0*boa)*sin(pi/3.0d0*boa)*cos(pi/3.0d0*boa)*exphub1*exphuc1
        !!dbohub=4.0d0/3.0d0*pi*sin(pi/3.0d0*bob)*sin(pi/3.0d0*bob)*sin(pi/3.0d0*bob)*cos(pi/3.0d0*bob)*exphua1*exphuc1
        !!dbohuc=4.0d0/3.0d0*pi*sin(pi/3.0d0*boc)*sin(pi/3.0d0*boc)*sin(pi/3.0d0*boc)*cos(pi/3.0d0*boc)*exphua1*exphub1


        arghu0=(arg2-1.0d0)*sinhd*sinhe
        ehulp=vconj(ity)*(arghu0+1.0d0)

        ecoh=ehulp*sbo

        estrain(j(2))=estrain(j(2)) + 0.50d0*ecoh*fij*fik*fjl
        estrain(j(3))=estrain(j(3)) + 0.50d0*ecoh*fij*fik*fjl

	!df398 derivative wrt arg
        decodar=sbo*vconj(ity)*2.0d0*arg*sinhd*sinhe*fij*fik*fjl
        
        !df398 derivatives wrt BO
        decodbola=dbohua*ehulp*fij*fik*fjl + ecoh*fik*fjl*dfijdBOij
        decodbolb=dbohub*ehulp*fij*fik*fjl + ecoh*fij*fjl*dfikdBOik
        decodbolc=dbohuc*ehulp*fij*fik*fjl + ecoh*fij*fik*dfjldBOjl
        
        !df398 derivatives wrt ht1 and ht2
        decodhd=coshd*sinhe*vconj(ity)*sbo*(arg2-1.0d0)*fij*fik*fjl
        decodhe=coshe*sinhd*vconj(ity)*sbo*(arg2-1.0d0)*fij*fik*fjl

	
        eco=eco+ecoh*fij*fik*fjl

        1 continue
    !*********************************************************************
    !                                                                    *
    !     Calculate derivative torsion angle and conjugation energy      *
    !     to cartesian coordinates                                       *
    !                                                                    *
    !*********************************************************************
        SINTH=SIN(THG(i1)*DGRRDN)
        IF (SINTH >= 0.0D0 .AND. SINTH < 1.0D-20) SINTH=1.0D-20
        IF (SINTH < 0.0D0 .AND. SINTH > -1.0D-20) SINTH=-1.0D-20
        IF (j(1) == IB(LA,2)) THEN
            DO  K1=1,3
                DRADC(K1,1)=DRDC(LA,K1,1)
                DRADC(K1,2)=DRDC(LA,K1,2)
            end do
        ELSE
            DO  K1=1,3
                DRADC(K1,1)=DRDC(LA,K1,2)
                DRADC(K1,2)=DRDC(LA,K1,1)
            end do
        ENDIF
        IF (j(2) == IB(LB,2)) THEN
            DO  K1=1,3
                DRBDC(K1,2)=DRDC(LB,K1,1)
                DRBDC(K1,3)=DRDC(LB,K1,2)
            end do
        ELSE
            DO K1=1,3
                DRBDC(K1,2)=DRDC(LB,K1,2)
                DRBDC(K1,3)=DRDC(LB,K1,1)
            end do
        ENDIF
        IF (j(3) == IB(LC,2)) THEN
            DO K1=1,3
                DRCDC(K1,3)=DRDC(LC,K1,1)
                DRCDC(K1,4)=DRDC(LC,K1,2)
            end do
        ELSE
            DO K1=1,3
                DRCDC(K1,3)=DRDC(LC,K1,2)
                DRCDC(K1,4)=DRDC(LC,K1,1)
            end do
        ENDIF

        do k1=1,3
            dhddc(1,k1)=dh1rdc(1,k1)
            dhddc(2,k1)=dh1rdc(2,k1)
            dhddc(3,k1)=dh1rdc(3,k1)
            dhedc(1,k1+1)=dh2rdc(1,k1)
            dhedc(2,k1+1)=dh2rdc(2,k1)
            dhedc(3,k1+1)=dh2rdc(3,k1)
        end do

    !*********************************************************************
    !     write (64,*)j(1),j(2),j(3),j(4)
    !     do k1=1,3
    !     write (64,'(10f12.4)')(dh1rdc(k1,k2),k2=1,3),
    !    $(dhdc(ld,k1,k2),k2=1,3),(dhddc(k1,k2),k2=1,4)
    !     write (64,'(10f12.4)')(dh2rdc(k1,k2),k2=1,3),
    !    $(dhdc(le,k1,k2),k2=1,3),(dhedc(k1,k2),k2=1,4)
    !     end do
    !     write (64,*)
    !*********************************************************************
        HTRA=RLA+COSHD*(RLC*COSHE-RLB)
        HTRB=RLB-RLA*COSHD-RLC*COSHE
        HTRC=RLC+COSHE*(RLA*COSHD-RLB)
        HTHD=RLA*SINHD*(RLB-RLC*COSHE)
        HTHE=RLC*SINHE*(RLB-RLA*COSHD)
        HNRA=RLC*SINHD*SINHE
        HNRC=RLA*SINHD*SINHE
        HNHD=RLA*RLC*COSHD*SINHE
        HNHE=RLA*RLC*SINHD*COSHE

        if (icpres == 0) then

            DO  K1=1,3
                DRDA(K1)=A(K1)/R4
                DO  K2=1,4
                    DRVDC(K1,K2)=DRDA(K1)*DADC(K2)
                    DTDC(K1,K2)=2.0d0*(DRADC(K1,K2)*HTRA+DRBDC(K1,K2)*HTRB+DRCDC(K1,K2 &
                    )*HTRC-DRVDC(K1,K2)*R4+DHDDC(K1,K2)*HTHD+DHEDC(K1,K2)*HTHE)
                    DNDC(K1,K2)=2.0d0*(DRADC(K1,K2)*HNRA+DRCDC(K1,K2)*HNRC+DHDDC(K1,K2 &
                    )*HNHD+DHEDC(K1,K2)*HNHE)
                    DARGTDC(i1,K1,K2)=(DTDC(K1,K2)-ARG*DNDC(K1,K2))/POEM
                     
                    D(K1,J(K2))=D(K1,J(K2))+DARGTDC(i1,K1,K2)*detdar+ &
                    dargtdc(i1,k1,k2)*decodar+(detdhd+decodhd)*dhddc(k1,k2)+ &
                    (detdhe+decodhe)*dhedc(k1,k2)
                end do
            end do
             
            do i2=1,idbo1(la)
                ihu=idbo(la,i2)
                do k1=1,3
                    d(k1,ihu)=d(k1,ihu)+dbondc(la,k1,i2)*(detdboa+decodbola)
                end do
            end do
             
            do i2=1,idbo1(lb)
                ihu=idbo(lb,i2)
                do k1=1,3
                    d(k1,ihu)=d(k1,ihu)+dbondc(lb,k1,i2)*(detdbob+decodbolb)
                    d(k1,ihu)=d(k1,ihu)+dbopindc(lb,k1,i2)*detdbopib
                end do
            end do
             
            do i2=1,idbo1(lc)
                ihu=idbo(lc,i2)
                do k1=1,3
                    d(k1,ihu)=d(k1,ihu)+dbondc(lc,k1,i2)*(detdboc+decodbolc)
                end do
            end do
             
            do i2=1,ia(j(2),2)
                iob=ia(j(2),2+i2)
                ncubo=nubon2(j(2),i2)
                do i3=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i3)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+detdsbo1*dbondc(ncubo,k1,i3)
                    end do
                end do
            end do

            do i2=1,ia(j(3),2)
                iob=ia(j(3),2+i2)
                ncubo=nubon2(j(3),i2)
                do i3=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i3)
                    do k1=1,3
                        d(k1,ihu)=d(k1,ihu)+detdsbo1*dbondc(ncubo,k1,i3)
                    end do
                end do
            end do

        else

            DO  K1=1,3
                DRDA(K1)=A(K1)/R4
                DO  K2=1,4
                    DRVDC(K1,K2)=DRDA(K1)*DADC(K2)
                    DTDC(K1,K2)=2.0d0*(DRADC(K1,K2)*HTRA+DRBDC(K1,K2)*HTRB+DRCDC(K1,K2 &
                    )*HTRC-DRVDC(K1,K2)*R4+DHDDC(K1,K2)*HTHD+DHEDC(K1,K2)*HTHE)
                    DNDC(K1,K2)=2.0d0*(DRADC(K1,K2)*HNRA+DRCDC(K1,K2)*HNRC+DHDDC(K1,K2 &
                    )*HNHD+DHEDC(K1,K2)*HNHE)
                    DARGTDC(i1,K1,K2)=(DTDC(K1,K2)-ARG*DNDC(K1,K2))/POEM
                    ix=nmpx(j(2),j(k2))
                    iy=nmpy(j(2),j(k2))
                    iz=nmpz(j(2),j(k2))
                    kcell=14+ix+3*iy+9*iz
                    dcell(k1,j(k2),kcell)=dcell(k1,j(k2),kcell) &
                    +dargtdc(i1,K1,K2)*detdar+ &
                    dargtdc(i1,k1,k2)*decodar+(detdhd+decodhd)*dhddc(k1,k2)+ &
                    (detdhe+decodhe)*dhedc(k1,k2)
                end do
            end do

            do i2=1,idbo1(la)
                ihu=idbo(la,i2)
                ix=nmpx(j(2),ihu)
                iy=nmpy(j(2),ihu)
                iz=nmpz(j(2),ihu)
                kcell=14+ix+3*iy+9*iz
                do k1=1,3
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                    dbondc(la,k1,i2)*(detdboa+decodbola)
                end do
            end do

            do i2=1,idbo1(lb)
                ihu=idbo(lb,i2)
                ix=nmpx(j(2),ihu)
                iy=nmpy(j(2),ihu)
                iz=nmpz(j(2),ihu)
                kcell=14+ix+3*iy+9*iz
                do k1=1,3
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                    dbondc(lb,k1,i2)*(detdbob+decodbolb)
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                    dbopindc(lb,k1,i2)*detdbopib
                end do
            end do

            do i2=1,idbo1(lc)
                ihu=idbo(lc,i2)
                ix=nmpx(j(2),ihu)
                iy=nmpy(j(2),ihu)
                iz=nmpz(j(2),ihu)
                kcell=14+ix+3*iy+9*iz
                do k1=1,3
                    dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                    dbondc(lc,k1,i2)*(detdboc+decodbolc)
                end do
            end do

            do i2=1,ia(j(2),2)
                iob=ia(j(2),2+i2)
                ncubo=nubon2(j(2),i2)
                do i3=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i3)
                    ix=nmpx(j(2),ihu)
                    iy=nmpy(j(2),ihu)
                    iz=nmpz(j(2),ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        detdsbo1*dbondc(ncubo,k1,i3)
                    end do
                end do
            end do

            do i2=1,ia(j(3),2)
                iob=ia(j(3),2+i2)
                ncubo=nubon2(j(3),i2)
                do i3=1,idbo1(ncubo)
                    ihu=idbo(ncubo,i3)
                    ix=nmpx(j(2),ihu)
                    iy=nmpy(j(2),ihu)
                    iz=nmpz(j(2),ihu)
                    kcell=14+ix+3*iy+9*iz
                    do k1=1,3
                        dcell(k1,ihu,kcell)=dcell(k1,ihu,kcell)+ &
                        detdsbo1*dbondc(ncubo,k1,i3)
                    end do
                end do
            end do

        end if

    10 END DO

    return
    end subroutine torang
!*********************************************************************
!*********************************************************************

    subroutine nonbon

!*********************************************************************
    include 'cbka.blk'
!*********************************************************************
!                                                                    *
!     Determine which vdWaals method to call                         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In nonbon',ivdwty
        call timer(65)
        close (65)
    end if
    if (ivdwty == 1) call nonbonnoinS
    if (ivdwty == 2) call nonboninnoS
    if (ivdwty == 3) call nonboninS
    if (ivdwty == 4) call nonboninnoSlg
    if (ivdwty < 1 .OR. ivdwty > 4) stop 'Unknown vdWaals-type'
    return
    end subroutine nonbon
!*********************************************************************
!*********************************************************************

    subroutine nonbonnoinS

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),da(6)
!*********************************************************************
!                                                                    *
!     Calculate vdWaals and Coulomb energies and derivatives         *
!     vdWaals no inner-wall, with shielding                          *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In nonbonnoinS'
        call timer(65)
        close (65)
    end if

    ew=0.0d0
    ep=0.0d0

    c1c=332.0638d0
    third=1.0d0/3.0d0
    fothird=4.0d0*third
    twothird=2.0d0*third
    h15=(vpar(29)-1.0d0)/vpar(29)

    do 10 ivl=1,nvpair-nvlself
        i1=nvl1(ivl)
        i2=nvl2(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        a(1)=c(i1,1)-c(i2,1)+ix*tm11
        a(2)=c(i1,2)-c(i2,2)+ix*tm21+iy*tm22
        a(3)=c(i1,3)-c(i2,3)+ix*tm31+iy*tm32+iz*tm33
    !*********************************************************************
    !                                                                    *
    !     Construct periodic images for each interaction                 *
    !                                                                    *
    !*********************************************************************
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 10
        ity1=ia(i1,1)
        ity2=ia(i2,1)
        imol1=iag(i1,3+mbond)
        imol2=iag(i2,3+mbond)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !                                                                    *
    !*********************************************************************
        p1=p1co(ity1,ity2)
        p2=p2co(ity1,ity2)
        p3=p3co(ity1,ity2)
        hulpw=(rr**vpar(29)+gamwco(ity1,ity2)) ! df398 f13 = shielding
        rrw=hulpw**(1.0d0/vpar(29))
        h1=exp(p3*(1.0d0-rrw/p1))
        h2=exp(0.50d0*p3*(1.0d0-rrw/p1))

        ewh=p2*(h1-2.0d0*h2)
        rrhuw=rr**(vpar(29)-1.0d0)
        dewdr=(p2*p3/p1)*(h2-h1)*rrhuw*(hulpw**(-h15))

    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i2)
        hulp1=(rr2*rr+gamcco(ity1,ity2))
        eph=c1c*q1q2/(hulp1**third)
        depdr=-c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=ewh*sw
        dewdrtap=dewdr*sw+ewh*sw1

    !     write (64,*)i1,i2,p1,p2,p3,gamwco(ity1,ity2),vpar(29),rr,ewh,ew
        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+0.50d0*(ewhtap+ephtap)
        estrain(i2)=estrain(i2)+0.50d0*(ewhtap+ephtap)
    !*********************************************************************
    !                                                                    *
    !     Calculate derivatives vdWaals energy to cartesian              *
    !     coordinates                                                    *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then

            do k4=1,3
                d(k4,i1)=d(k4,i1)+(dewdrtap+depdrtap)*(a(k4)/rr)
                d(k4,i2)=d(k4,i2)-(dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        else

            kcell=14+ix+3*iy+9*iz
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i2,kcell)=dcell(k4,i2,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        end if

    10 END DO
!*********************************************************************
!                                                                    *
!     Add interaction of atoms with the corresponding atom           *
!     in the surrounding periodic cells                              *
!                                                                    *
!*********************************************************************
    do 20 ivl=nvpair-nvlself,nvpair
        i1=nvl1(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        kcell=14+ix+3*iy+9*iz
        a(1)=ix*tm11
        a(2)=ix*tm21+iy*tm22
        a(3)=ix*tm31+iy*tm32+iz*tm33
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 20
        ity1=ia(i1,1)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !                                                                    *
    !*********************************************************************
        p1=p1co(ity1,ity1)
        p2=p2co(ity1,ity1)
        p3=p3co(ity1,ity1)

        hulpw=(rr**vpar(29)+gamwco(ity1,ity1))  ! df398 f13 = shielding
        rrw=hulpw**(1.0d0/vpar(29))
        h1=exp(p3*(1.0d0-rrw/p1))
        h2=exp(0.50d0*p3*(1.0d0-rrw/p1))

        ewh=0.50d0*p2*(h1-2.0d0*h2)
        rrhuw=rr**(vpar(29)-1.0d0)
        dewdr=0.50d0*(p2*p3/p1)*(h2-h1)*rrhuw*(hulpw**(-h15))
    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i1)
        hulp1=(rr*rr2+gamcco(ity1,ity1))
        eph=0.50d0*c1c*q1q2/(hulp1**third)
        depdr=-0.50d0*c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=ewh*sw
        dewdrtap=dewdr*sw+ewh*sw1

        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+ewhtap+ephtap

        if (icpres == 1) then
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i1,kcell)=dcell(k4,i1,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do
        end if

    20 END DO

!*********************************************************************
!                                                                    *
!     Calculate derivatives Coulomb energy to cartesian              *
!     coordinates                                                    *
!                                                                    *
!*********************************************************************

!     do i1=1,na
!     do k1=1,3
!     sum=0.0
!     do 20 i2=1,na

!     do 20 i3=i2+1,na
!     if (dqdc(i3,i1,k1).le.zero.and.dqdc(i2,i1,k1).le.zero)
!    $goto 20
!     rr=rrs(i2,i3)
!     ity1=ia(i2,1)
!     ity2=ia(i3,1)
!     gamt=sqrt(gam(ity1)*gam(ity2))
!     hulp1=(rr**3+(1.0/(gamt**3)))
!     depdq1=-c1c*ch(i3)/(hulp1**third)
!     depdq2=-c1c*ch(i2)/(hulp1**third)
!     sum=sum+depdq1*dqdc(i2,i1,k1)+depdq2*dqdc(i3,i1,k1)
!  20 continue
!     d(i1,k1)=d(i1,k1)+sum

!     end do
!     end do

    return
    end subroutine nonbonnoinS

!*********************************************************************
!*********************************************************************

    subroutine nonboninnoS

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),da(6)
!*********************************************************************
!                                                                    *
!     Calculate vdWaals and Coulomb energies and derivatives         *
!     vdWaals include inner-wall, no shielding                       *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In nonboninnoS'
        call timer(65)
        close (65)
    end if

    ew=0.0d0
    ep=0.0d0

    c1c=332.0638d0
    third=1.0d0/3.0d0
    fothird=4.0d0*third
    twothird=2.0d0*third
    h15=(vpar(29)-1.0d0)/vpar(29)

    do 10 ivl=1,nvpair-nvlself
        i1=nvl1(ivl)
        i2=nvl2(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        a(1)=c(i1,1)-c(i2,1)+ix*tm11
        a(2)=c(i1,2)-c(i2,2)+ix*tm21+iy*tm22
        a(3)=c(i1,3)-c(i2,3)+ix*tm31+iy*tm32+iz*tm33
    !*********************************************************************
    !                                                                    *
    !     Construct periodic images for each interaction                 *
    !                                                                    *
    !*********************************************************************
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 10
        ity1=ia(i1,1)
        ity2=ia(i2,1)
        imol1=iag(i1,3+mbond)
        imol2=iag(i2,3+mbond)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !                                                                    *
    !*********************************************************************
        p1=p1co(ity1,ity2)
        p2=p2co(ity1,ity2)
        p3=p3co(ity1,ity2)

        h1=exp(p3*(1.0d0-rr/p1))
        h2=exp(0.50d0*p3*(1.0d0-rr/p1))
        ewh=p2*(h1-2.0d0*h2)
        dewdr=(p2*p3/p1)*(h2-h1)
    !*********************************************************************
    !                                                                    *
    !     Calculate inner core repulsion                                 *
    !                                                                    *
    !*********************************************************************
        pc1=rcore(ity1,ity2)
        pc2=ecore(ity1,ity2)
        pc3=acore(ity1,ity2)
        ecoreh=pc2*(exp(pc3*(1.0d0-(rr/pc1))))
        decoredr=-(pc3/pc1)*ecoreh
    !     rr10=rr2*rr2*rr2*rr2*rr2      !Power-10 function
    !     ecoreh=ecore(ity1,ity2)*(rcore(ity1,ity2)/rr10)
    !     decoredr=-10.0*ecoreh/rr
    !     write (65,*)i1,i2,rcore(ity1,ity2),ecore(ity1,ity2),
    !    $acore(ity1,ity2),rr,ecoreh
    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i2)
        hulp1=(rr2*rr+gamcco(ity1,ity2))
        eph=c1c*q1q2/(hulp1**third)
        depdr=-c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=(ewh+ecoreh)*sw
        dewdrtap=(dewdr+decoredr)*sw+(ewh+ecoreh)*sw1

    !     write (64,*)i1,i2,p1,p2,p3,gamwco(ity1,ity2),vpar(29),rr,ewh,ew
        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+0.50d0*(ewhtap+ephtap)
        estrain(i2)=estrain(i2)+0.50d0*(ewhtap+ephtap)
    !*********************************************************************
    !                                                                    *
    !     Calculate derivatives vdWaals energy to cartesian              *
    !     coordinates                                                    *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then

            do k4=1,3
                d(k4,i1)=d(k4,i1)+(dewdrtap+depdrtap)*(a(k4)/rr)
                d(k4,i2)=d(k4,i2)-(dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        else

            kcell=14+ix+3*iy+9*iz
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i2,kcell)=dcell(k4,i2,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        end if

    10 END DO
!*********************************************************************
!                                                                    *
!     Add interaction of atoms with the corresponding atom           *
!     in the surrounding periodic cells                              *
!                                                                    *
!*********************************************************************
    do 20 ivl=nvpair-nvlself,nvpair
        i1=nvl1(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        kcell=14+ix+3*iy+9*iz
        a(1)=ix*tm11
        a(2)=ix*tm21+iy*tm22
        a(3)=ix*tm31+iy*tm32+iz*tm33
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 20
        ity1=ia(i1,1)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !     !df398> 01.12.19 removed shielding. original subroutine had it *
    !             included (bug?)                                        *
    !*********************************************************************
        p1=p1co(ity1,ity1)
        p2=p2co(ity1,ity1)
        p3=p3co(ity1,ity1)
  
  ! Remove vdW shielding
  !      hulpw=(rr**vpar(29)+gamwco(ity1,ity1))
  !      rrw=hulpw**(1.0/vpar(29))
  !      h1=exp(p3*(1.0-rrw/p1))
  !      h2=exp(0.50*p3*(1.0-rrw/p1))
  
  !      ewh=0.50*p2*(h1-2.0*h2)
  !      rrhuw=rr**(vpar(29)-1.0)
  !      dewdr=0.50*(p2*p3/p1)*(h2-h1)*rrhuw*(hulpw**-h15)
  
        h1=exp(p3*(1.0-rr/p1))
        h2=exp(0.50*p3*(1.0-rr/p1))
        ewh=0.50*p2*(h1-2.0*h2)
        dewdr=0.50*(p2*p3/p1)*(h2-h1)
    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i1)
        hulp1=(rr*rr2+gamcco(ity1,ity1))
        eph=0.50d0*c1c*q1q2/(hulp1**third)
        depdr=-0.50d0*c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=ewh*sw
        dewdrtap=dewdr*sw+ewh*sw1

        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+ewhtap+ephtap

        if (icpres == 1) then
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i1,kcell)=dcell(k4,i1,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do
        end if

    20 END DO

    return
    end subroutine nonboninnoS

!*********************************************************************
!*********************************************************************

    subroutine nonboninnoSlg

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),da(6)
!*********************************************************************
!                                                                    *
!     Calculate vdWaals and Coulomb energies and derivatives         *
!     vdWaals include inner-wall + lg , no shielding                 *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In nonboninnoSlg'
        call timer(65)
        close (65)
    end if

    ew=0.0d0
    ep=0.0d0

    c1c=332.0638d0
    third=1.0d0/3.0d0
    fothird=4.0d0*third
    twothird=2.0d0*third
    h15=(vpar(29)-1.0d0)/vpar(29)

    do 10 ivl=1,nvpair-nvlself
        i1=nvl1(ivl)
        i2=nvl2(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        a(1)=c(i1,1)-c(i2,1)+ix*tm11
        a(2)=c(i1,2)-c(i2,2)+ix*tm21+iy*tm22
        a(3)=c(i1,3)-c(i2,3)+ix*tm31+iy*tm32+iz*tm33
    !*********************************************************************
    !                                                                    *
    !     Construct periodic images for each interaction                 *
    !                                                                    *
    !*********************************************************************
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 10
        ity1=ia(i1,1)
        ity2=ia(i2,1)
        imol1=iag(i1,3+mbond)
        imol2=iag(i2,3+mbond)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !                                                                    *
    !*********************************************************************
        p1=p1co(ity1,ity2)
        p2=p2co(ity1,ity2)
        p3=p3co(ity1,ity2)

        h1=exp(p3*(1.0d0-rr/p1))
        h2=exp(0.50d0*p3*(1.0d0-rr/p1))
        ewh=p2*(h1-2.0d0*h2)
        dewdr=(p2*p3/p1)*(h2-h1)
    !*********************************************************************
    !                                                                    *
    !     Calculate inner core repulsion                                 *
    !                                                                    *
    !*********************************************************************
        pc1=rcore(ity1,ity2)
        pc2=ecore(ity1,ity2)
        pc3=acore(ity1,ity2)
        ecoreh=pc2*(exp(pc3*(1.0d0-(rr/pc1))))
        decoredr=-(pc3/pc1)*ecoreh
    !     rr10=rr2*rr2*rr2*rr2*rr2      !Power-10 function
    !     ecoreh=ecore(ity1,ity2)*(rcore(ity1,ity2)/rr10)
    !     decoredr=-10.0*ecoreh/rr
    !     write (65,*)i1,i2,rcore(ity1,ity2),ecore(ity1,ity2),
    !    $acore(ity1,ity2),rr,ecoreh
    !*********************************************************************
    !                                                                    *
    !     Calculate long range dispersion (Edisp=-C/(r^6+b*re^6))        *
    !                                                                    *
    !*********************************************************************
        pdisp=dispc6(ity1,ity2)
        pdispre=dispre(ity1,ity2)
        pdispre6=pdispre**6
        rr6=rr**6
        rr5=rr6/rr
        edisp=-pdisp/(rr6+dispscale*pdispre6)
        dedispdr=6.0*edisp*edisp*rr5/pdisp
    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i2)
        hulp1=(rr2*rr+gamcco(ity1,ity2))
        eph=c1c*q1q2/(hulp1**third)
        depdr=-c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=(ewh+ecoreh+edisp)*sw
        dewdrtap=(dewdr+decoredr+dedispdr)*sw+(ewh+ecoreh+edisp)*sw1
  
  !     write (64,*)i1,i2,p1,p2,p3,gamwco(ity1,ity2),vpar(29),rr,ewh,ew
        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+0.50*(ewhtap+ephtap)
        estrain(i2)=estrain(i2)+0.50*(ewhtap+ephtap)

    !*********************************************************************
    !                                                                    *
    !     Calculate derivatives vdWaals energy to cartesian              *
    !     coordinates                                                    *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then

            do k4=1,3
                d(k4,i1)=d(k4,i1)+(dewdrtap+depdrtap)*(a(k4)/rr)
                d(k4,i2)=d(k4,i2)-(dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        else

            kcell=14+ix+3*iy+9*iz
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i2,kcell)=dcell(k4,i2,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        end if

    10 END DO
!*********************************************************************
!                                                                    *
!     Add interaction of atoms with the corresponding atom           *
!     in the surrounding periodic cells                              *
!                                                                    *
!*********************************************************************
    do 20 ivl=nvpair-nvlself,nvpair
        i1=nvl1(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        kcell=14+ix+3*iy+9*iz
        a(1)=ix*tm11
        a(2)=ix*tm21+iy*tm22
        a(3)=ix*tm31+iy*tm32+iz*tm33
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 20
        ity1=ia(i1,1)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !                                                                    *
    !*********************************************************************
        p1=p1co(ity1,ity1)
        p2=p2co(ity1,ity1)
        p3=p3co(ity1,ity1)
  
  ! Remove vdW shielding
  !      hulpw=(rr**vpar(29)+gamwco(ity1,ity1))
  !      rrw=hulpw**(1.0/vpar(29))
  !      h1=exp(p3*(1.0-rrw/p1))
  !      h2=exp(0.50*p3*(1.0-rrw/p1))
  
  !      ewh=0.50*p2*(h1-2.0*h2)
  !      rrhuw=rr**(vpar(29)-1.0)
  !      dewdr=0.50*(p2*p3/p1)*(h2-h1)*rrhuw*(hulpw**-h15)
  
        h1=exp(p3*(1.0-rr/p1))
        h2=exp(0.50*p3*(1.0-rr/p1))
        ewh=0.50*p2*(h1-2.0*h2)
        dewdr=0.50*(p2*p3/p1)*(h2-h1)
    !********************************************************************
    !                                                                    *
    !     Calculate inner core repulsion                                 *
    !                                                                    *
    !*********************************************************************
        pc1=rcore(ity1,ity1)
        pc2=ecore(ity1,ity1)
        pc3=acore(ity1,ity1)
        ecoreh=0.5*pc2*(exp(pc3*(1.0-(rr/pc1))))
        decoredr=-0.5*(pc3/pc1)*ecoreh
    !*********************************************************************
    !                                                                    *
    !     Calculate long range dispersion (Edisp=-C/(r^6+d*re^6))        *
    !                                                                    *
    !*********************************************************************
        pdisp=dispc6(ity1,ity1)
        pdispre=dispre(ity1,ity1)
        pdispre6=pdispre**6
        rr6=rr**6
        rr5=rr6/rr
        edisp=-0.5*pdisp/(rr6+dispscale*pdispre6)
        dedispdr=3.0*edisp*edisp*rr5/pdisp
    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i1)
        hulp1=(rr*rr2+gamcco(ity1,ity1))
        eph=0.50d0*c1c*q1q2/(hulp1**third)
        depdr=-0.50d0*c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=(ewh+ecoreh+edisp)*sw
        dewdrtap=(dewdr+decoredr+dedispdr)*sw+(ewh+ecoreh+edisp)*sw1

        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+ewhtap+ephtap

        if (icpres == 1) then
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i1,kcell)=dcell(k4,i1,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do
        end if

    20 END DO

    return
    end subroutine nonboninnoSlg

!*********************************************************************
!*********************************************************************

    subroutine nonboninS

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),da(6)
!*********************************************************************
!                                                                    *
!     Calculate vdWaals and Coulomb energies and derivatives         *
!     vdWaals include inner-wall with shielding                      *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In nonbonins'
        call timer(65)
        close (65)
    end if

    ew=0.0d0
    ep=0.0d0

    c1c=332.0638d0
    third=1.0d0/3.0d0
    fothird=4.0d0*third
    twothird=2.0d0*third
    h15=(vpar(29)-1.0d0)/vpar(29)

    do 10 ivl=1,nvpair-nvlself
        i1=nvl1(ivl)
        i2=nvl2(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        a(1)=c(i1,1)-c(i2,1)+ix*tm11
        a(2)=c(i1,2)-c(i2,2)+ix*tm21+iy*tm22
        a(3)=c(i1,3)-c(i2,3)+ix*tm31+iy*tm32+iz*tm33
    !*********************************************************************
    !                                                                    *
    !     Construct periodic images for each interaction                 *
    !                                                                    *
    !*********************************************************************
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 10
        ity1=ia(i1,1)
        ity2=ia(i2,1)
        imol1=iag(i1,3+mbond)
        imol2=iag(i2,3+mbond)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !                                                                    *
    !*********************************************************************
        p1=p1co(ity1,ity2)
        p2=p2co(ity1,ity2)
        p3=p3co(ity1,ity2)
        hulpw=(rr**vpar(29)+gamwco(ity1,ity2))  ! f13, vdWaals with shielding
        rrw=hulpw**(1.0d0/vpar(29))
        h1=exp(p3*(1.0d0-rrw/p1))
        h2=exp(0.50d0*p3*(1.0d0-rrw/p1))

        ewh=p2*(h1-2.0d0*h2)
        rrhuw=rr**(vpar(29)-1.0d0)
        dewdr=(p2*p3/p1)*(h2-h1)*rrhuw*(hulpw**(-h15))
    !*********************************************************************
    !                                                                    *
    !     Calculate inner core repulsion                                 *
    !                                                                    *
    !*********************************************************************
        pc1=rcore(ity1,ity2)
        pc2=ecore(ity1,ity2)
        pc3=acore(ity1,ity2)
        ecoreh=pc2*(exp(pc3*(1.0d0-(rr/pc1))))
        decoredr=-(pc3/pc1)*ecoreh
    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i2)
        hulp1=(rr2*rr+gamcco(ity1,ity2))
        eph=c1c*q1q2/(hulp1**third)
        depdr=-c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=(ewh+ecoreh)*sw
        dewdrtap=(dewdr+decoredr)*sw+(ewh+ecoreh)*sw1

    !     write (64,*)i1,i2,p1,p2,p3,gamwco(ity1,ity2),vpar(29),rr,ewh,ew
        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+0.50d0*(ewhtap+ephtap)
        estrain(i2)=estrain(i2)+0.50d0*(ewhtap+ephtap)
    !*********************************************************************
    !                                                                    *
    !     Calculate derivatives vdWaals energy to cartesian              *
    !     coordinates                                                    *
    !                                                                    *
    !*********************************************************************
        if (icpres == 0) then

            do k4=1,3
                d(k4,i1)=d(k4,i1)+(dewdrtap+depdrtap)*(a(k4)/rr)
                d(k4,i2)=d(k4,i2)-(dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        else

            kcell=14+ix+3*iy+9*iz
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i2,kcell)=dcell(k4,i2,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do

        end if

    10 END DO
!*********************************************************************
!                                                                    *
!     Add interaction of atoms with the corresponding atom           *
!     in the surrounding periodic cells                              *
!                                                                    *
!*********************************************************************
    do 20 ivl=nvpair-nvlself,nvpair
        i1=nvl1(ivl)
        ix=nvlx(ivl)
        iy=nvly(ivl)
        iz=nvlz(ivl)
        kcell=14+ix+3*iy+9*iz
        a(1)=ix*tm11
        a(2)=ix*tm21+iy*tm22
        a(3)=ix*tm31+iy*tm32+iz*tm33
        rr=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
        if (rr > swb .OR. rr < 0.001d0) goto 20
        ity1=ia(i1,1)
        rr2=rr*rr

        sw=1.0d0
        sw1=0.0d0
        call taper(rr,rr2)
    !*********************************************************************
    !                                                                    *
    !     Calculate vdWaals energy                                       *
    !                                                                    *
    !*********************************************************************
        p1=p1co(ity1,ity1)
        p2=p2co(ity1,ity1)
        p3=p3co(ity1,ity1)

        hulpw=(rr**vpar(29)+gamwco(ity1,ity1)) ! df398 vdWaals with shielding
        rrw=hulpw**(1.0d0/vpar(29))
        h1=exp(p3*(1.0d0-rrw/p1))
        h2=exp(0.50d0*p3*(1.0d0-rrw/p1))

        ewh=0.50d0*p2*(h1-2.0d0*h2)
        rrhuw=rr**(vpar(29)-1.0d0)
        dewdr=0.50d0*(p2*p3/p1)*(h2-h1)*rrhuw*(hulpw**(-h15))
    !*********************************************************************
    !                                                                    *
    !     Calculate Coulomb energy                                       *
    !                                                                    *
    !*********************************************************************
        q1q2=ch(i1)*ch(i1)
        hulp1=(rr*rr2+gamcco(ity1,ity1))
        eph=0.50d0*c1c*q1q2/(hulp1**third)
        depdr=-0.50d0*c1c*q1q2*rr2/(hulp1**fothird)
    !*********************************************************************
    !                                                                    *
    !     Taper correction                                               *
    !                                                                    *
    !*********************************************************************
        ephtap=eph*sw
        depdrtap=depdr*sw+eph*sw1
        ewhtap=ewh*sw
        dewdrtap=dewdr*sw+ewh*sw1

        ew=ew+ewhtap
        ep=ep+ephtap
        estrain(i1)=estrain(i1)+ewhtap+ephtap

        if (icpres == 1) then
            do k4=1,3
                dcell(k4,i1,14)=dcell(k4,i1,14)+ &
                (dewdrtap+depdrtap)*(a(k4)/rr)
                dcell(k4,i1,kcell)=dcell(k4,i1,kcell)- &
                (dewdrtap+depdrtap)*(a(k4)/rr)
            end do
        end if

    20 END DO

    return
    end subroutine nonboninS

!*********************************************************************
!*********************************************************************
     
    subroutine efield
     
!*********************************************************************
    include 'cbka.blk'
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In efield'
        call timer(65)
        close (65)
    end if
!*********************************************************************
!                                                                    *
!     Electric field                                                 *
!                                                                    *
!*********************************************************************
    efi=0.0d0
    efix=0.0d0
    efiy=0.0d0
    efiz=0.0d0
    c1c=332.0638d0       !Coulomb energy conversion
     
    if (ifieldx == 1) then
        do i1=1,na
            efih=vfieldx*23.02d0*c1c*ch(i1)*c(i1,1)
            estrain(i1)=estrain(i1)+efih
            efix=efix+efih
            defidc=23.02d0*c1c*vfieldx*ch(i1)
            d(1,i1)=d(1,i1)+defidc
        end do
    end if
     
    if (ifieldy == 1) then
        do i1=1,na
            efih=vfieldy*23.02d0*c1c*ch(i1)*c(i1,2)
            estrain(i1)=estrain(i1)+efih
            efiy=efiy+efih
            defidc=23.02d0*c1c*vfieldy*ch(i1)
            d(2,i1)=d(2,i1)+defidc
        end do
    end if
     
    if (ifieldz == 1) then
        do i1=1,na
            efih=vfieldz*23.02d0*c1c*ch(i1)*c(i1,3)
            estrain(i1)=estrain(i1)+efih
            efiz=efiz+efih
            defidc=23.02d0*c1c*vfieldz*ch(i1)
            d(3,i1)=d(3,i1)+defidc
        end do
    end if
     
    efi=efix+efiy+efiz
    return
    end subroutine efield
!*********************************************************************
!*********************************************************************

    subroutine radbo

!*********************************************************************
    include 'cbka.blk'
    dimension a1(3)
!*********************************************************************
!                                                                    *
!     Calculate radical/double bond energy (to increase reaction     *
!     rates)                                                         *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In radbo'
        call timer(65)
        close (65)
    end if
!     eradbo=zero
!     do i1=1,na
!     if (xmasat(i1).gt.2.0) qa(i1)='C '
!     end do
!     nradcount=nradcount+1
!     if (nradcount.lt.nrddf) return
!     if (nradcount.gt.2*nrddf) nradcount=0

!     do i1=1,nmolo5
!     if (elmol(i1).gt.2*int(elmol(i1)*0.50)) then
!     vlpmax=0.0
!     do i2=1,nmolat2(i1,1)
!     eradmin=50.0
!     ihu=nmolat2(i1,i2+1)
!     if (xmasat(ihu).gt.2.0.and.vlp(ihu).gt.zero) then !no H-atoms
!     vlps=vlp(ihu)

!     do 10 i3=1,na
!     if (xmasat(i3).lt.2.0) goto 10       !no H-atoms
!     imol2=iag(i3,3+mbond)
!     if (i1.eq.imol2) goto 10            !no intermolecular reactions
!     bopisum=zero
!     do i4=1,ia(i3,2)
!     ihu2=nubon2(i3,i4)
!     bopisum=bopisum+bopi(ihu2)
!     end do
!     if (bopisum.lt.0.25) goto 10        !only atoms in double bonds
!     dirb1=dista(ihu,i3)
!     hu1=vpar(22)-dirb1
!     expvl=exp(-5.0*vlps)
!     hexp=1.0-expvl
!     if (hu1.ge.zero) then
!     exphu1=exp(-vpar(25)*hu1*hu1)
!     eradboh1=hexp*vpar(24)*exphu1
!     else
!     exphu1=exp(-vpar(26)*hu1*hu1)
!     eradboh1=hexp*vpar(24)*exphu1
!     end if

!     if (eradboh1.lt.eradmin) then
!     i3s=i3
!     disrmin=dirb1
!     eradmin=eradboh1
!     end if

!  10 continue

!     qa(ihu)='S '
!     qa(i3s)='N '
!     hu1=vpar(22)-disrmin
!     eradboh1=zero
!     expvl=exp(-5.0*vlps)
!     hexp=1.0-expvl
!     if (hu1.ge.zero) then
!     exphu1=exp(-vpar(25)*hu1*hu1)
!     eradboh1=hexp*vpar(24)*exphu1
!     deradbodr1=hexp*vpar(24)*vpar(25)*2.0*hu1*exphu1
!     else
!     exphu1=exp(-vpar(26)*hu1*hu1)
!     eradboh1=hexp*vpar(24)*exphu1
!     deradbodr1=hexp*vpar(24)*vpar(26)*2.0*hu1*exphu1
!     end if
!     eradbo=eradbo+eradboh1

!     a1(1)=dxm(ihu,i3s)
!     a1(2)=dym(ihu,i3s)
!     a1(3)=dzm(ihu,i3s)

!     do k1=1,3
!     d(k1,ihu)=d(k1,ihu)+deradbodr1*(a1(k1)/disrmin)
!     d(k1,i3s)=d(k1,i3s)-deradbodr1*(a1(k1)/disrmin)
!     end do

!     end if
!     end do

!     end if
!     end do


    return
    end subroutine radbo
!*********************************************************************
!*********************************************************************

    subroutine restraint

!*********************************************************************
    include 'cbka.blk'
    parameter (msymgroups=20)
    parameter (msymopts=15)
    dimension drda(3),j(4),dhrdc(3,3),dargdc(3,3)
    dimension fc(nat,3),csym(msymopts*nat,3),ccsym(msymopts*nat,3)
    dimension cs(nat,3)
    dimension qatsym(msymopts*nat)
    dimension vsymadd(msymgroups,msymopts,3)
    dimension vsymmulx(msymgroups,msymopts,3)
    dimension vsymmuly(msymgroups,msymopts,3)
    dimension vsymmulz(msymgroups,msymopts,3)
    dimension nsymopt(msymgroups)
    dimension qsymm(msymgroups)
    character(2) :: qatsym
    character(6) :: qsymm
!*********************************************************************
!                                                                    *
!     Calculate restraint energies                                   *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In restraint'
        call timer(65)
        close (65)
    end if

    if (nrestras > 0) then
    !*********************************************************************
    !                                                                    *
    !     Define symmetry groups                                         *
    !                                                                    *
    !*********************************************************************
        nsymel=3
        qsymm(1)='    P1'
        nsymopt(1)=1
        vsymadd(1,1,1)=0.0d0          !x=1.0*x   y=1.0*y   z=1.0*z
        vsymadd(1,1,2)=0.0d0
        vsymadd(1,1,3)=0.0d0
        vsymmulx(1,1,1)=1.0d0
        vsymmulx(1,1,2)=0.0d0
        vsymmulx(1,1,3)=0.0d0
        vsymmuly(1,1,1)=0.0d0
        vsymmuly(1,1,2)=1.0d0
        vsymmuly(1,1,3)=0.0d0
        vsymmulz(1,1,1)=0.0d0
        vsymmulz(1,1,2)=1.0d0
        vsymmulz(1,1,3)=0.0d0

        qsymm(2)='   C2V'
        nsymopt(2)=2
        vsymadd(2,1,1)=0.0d0          !x=1.0d0*x   y=1.0*y  z=1.0*z
        vsymadd(2,1,2)=0.0d0
        vsymadd(2,1,3)=0.0d0
        vsymmulx(2,1,1)=1.0d0
        vsymmulx(2,1,2)=0.0d0
        vsymmulx(2,1,3)=0.0d0
        vsymmuly(2,1,1)=0.0d0
        vsymmuly(2,1,2)=1.0d0
        vsymmuly(2,1,3)=0.0d0
        vsymmulz(2,1,1)=0.0d0
        vsymmulz(2,1,2)=0.0d0
        vsymmulz(2,1,3)=1.0d0

        vsymadd(2,2,1)=0.5d0          !x=0.5-1.0d0*y   y=0.5-1.0*x  z=1.0*z
        vsymadd(2,2,2)=0.5d0
        vsymadd(2,2,3)=0.0d0
        vsymmulx(2,2,1)=0.0d0
        vsymmulx(2,2,2)=-1.0d0
        vsymmulx(2,2,3)=0.0d0
        vsymmuly(2,2,1)=-1.0d0
        vsymmuly(2,2,2)=0.0d0
        vsymmuly(2,2,3)=0.0d0
        vsymmulz(2,2,1)=0.0d0
        vsymmulz(2,2,2)=0.0d0
        vsymmulz(2,2,3)=1.0d0

        qsymm(3)='  Pmmm'
        nsymopt(3)=8
        vsymadd(3,1,1)=0.0d0          !x=1.0d0*x   y=1.0*y  z=1.0*z
        vsymadd(3,1,2)=0.0d0
        vsymadd(3,1,3)=0.0d0
        vsymmulx(3,1,1)=1.0d0
        vsymmulx(3,1,2)=0.0d0
        vsymmulx(3,1,3)=0.0d0
        vsymmuly(3,1,1)=0.0d0
        vsymmuly(3,1,2)=1.0d0
        vsymmuly(3,1,3)=0.0d0
        vsymmulz(3,1,1)=0.0d0
        vsymmulz(3,1,2)=0.0d0
        vsymmulz(3,1,3)=1.0d0

        vsymadd(3,2,1)=0.0d0          !x=1.0d0*x   y=-1.0*y  z=1.0*z
        vsymadd(3,2,2)=0.0d0
        vsymadd(3,2,3)=0.0d0
        vsymmulx(3,2,1)=1.0d0
        vsymmulx(3,2,2)=0.0d0
        vsymmulx(3,2,3)=0.0d0
        vsymmuly(3,2,1)=0.0d0
        vsymmuly(3,2,2)=-1.0d0
        vsymmuly(3,2,3)=0.0d0
        vsymmulz(3,2,1)=0.0d0
        vsymmulz(3,2,2)=0.0d0
        vsymmulz(3,2,3)=1.0d0

        vsymadd(3,3,1)=0.0d0          !x=-1.0d0*x   y=1.0*y  z=-1.0*z
        vsymadd(3,3,2)=0.0d0
        vsymadd(3,3,3)=0.0d0
        vsymmulx(3,3,1)=-1.0d0
        vsymmulx(3,3,2)=0.0d0
        vsymmulx(3,3,3)=0.0d0
        vsymmuly(3,3,1)=0.0d0
        vsymmuly(3,3,2)=1.0d0
        vsymmuly(3,3,3)=0.0d0
        vsymmulz(3,3,1)=0.0d0
        vsymmulz(3,3,2)=0.0d0
        vsymmulz(3,3,3)=-1.0d0

        vsymadd(3,4,1)=0.0d0          !x=-1.0d0*x   y=1.0*y  z=1.0*z
        vsymadd(3,4,2)=0.0d0
        vsymadd(3,4,3)=0.0d0
        vsymmulx(3,4,1)=-1.0d0
        vsymmulx(3,4,2)=0.0d0
        vsymmulx(3,4,3)=0.0d0
        vsymmuly(3,4,1)=0.0d0
        vsymmuly(3,4,2)=1.0d0
        vsymmuly(3,4,3)=0.0d0
        vsymmulz(3,4,1)=0.0d0
        vsymmulz(3,4,2)=0.0d0
        vsymmulz(3,4,3)=1.0d0

        vsymadd(3,5,1)=0.0d0          !x=1.0d0*x   y=-1.0*y  z=-1.0*z
        vsymadd(3,5,2)=0.0d0
        vsymadd(3,5,3)=0.0d0
        vsymmulx(3,5,1)=1.0d0
        vsymmulx(3,5,2)=0.0d0
        vsymmulx(3,5,3)=0.0d0
        vsymmuly(3,5,1)=0.0d0
        vsymmuly(3,5,2)=-1.0d0
        vsymmuly(3,5,3)=0.0d0
        vsymmulz(3,5,1)=0.0d0
        vsymmulz(3,5,2)=0.0d0
        vsymmulz(3,5,3)=-1.0d0

        vsymadd(3,6,1)=0.0d0          !x=1.0d0*x   y=1.0*y  z=-1.0*z
        vsymadd(3,6,2)=0.0d0
        vsymadd(3,6,3)=0.0d0
        vsymmulx(3,6,1)=1.0d0
        vsymmulx(3,6,2)=0.0d0
        vsymmulx(3,6,3)=0.0d0
        vsymmuly(3,6,1)=0.0d0
        vsymmuly(3,6,2)=1.0d0
        vsymmuly(3,6,3)=0.0d0
        vsymmulz(3,6,1)=0.0d0
        vsymmulz(3,6,2)=0.0d0
        vsymmulz(3,6,3)=-1.0d0

        vsymadd(3,7,1)=0.0d0          !x=-1.0d0*x   y=-1.0*y  z=1.0*z
        vsymadd(3,7,2)=0.0d0
        vsymadd(3,7,3)=0.0d0
        vsymmulx(3,7,1)=-1.0d0
        vsymmulx(3,7,2)=0.0d0
        vsymmulx(3,7,3)=0.0d0
        vsymmuly(3,7,1)=0.0d0
        vsymmuly(3,7,2)=-1.0d0
        vsymmuly(3,7,3)=0.0d0
        vsymmulz(3,7,1)=0.0d0
        vsymmulz(3,7,2)=0.0d0
        vsymmulz(3,7,3)=1.0d0

        vsymadd(3,8,1)=0.0d0          !x=-1.0d0*x   y=-1.0*y  z=-1.0*z
        vsymadd(3,8,2)=0.0d0
        vsymadd(3,8,3)=0.0d0
        vsymmulx(3,8,1)=-1.0d0
        vsymmulx(3,8,2)=0.0d0
        vsymmulx(3,8,3)=0.0d0
        vsymmuly(3,4,1)=0.0d0
        vsymmuly(3,8,2)=1.0d0
        vsymmuly(3,8,3)=0.0d0
        vsymmulz(3,8,1)=0.0d0
        vsymmulz(3,8,2)=0.0d0
        vsymmulz(3,8,3)=1.0d0

    end if

!*********************************************************************
!                                                                    *
!     Calculate distance restraint energy                            *
!                                                                    *
!*********************************************************************
    do i1=1,nrestra
        ih1=irstra(i1,1)
        ih2=irstra(i1,2)
        if (itend(i1) == 0 .OR. (mdstep > itstart(i1) .AND. mdstep < &
        itend(i1))) then
            call dista2(ih1,ih2,rr,dx,dy,dz)
            diffr=rr-rrstra(i1)
        !     diffr=rrstra(i1)
            exphu=exp(-vkrst2(i1)*(diffr*diffr))
            erh=vkrstr(i1)*(1.0d0-exphu)
            deresdr=2.0d0*vkrst2(i1)*diffr*vkrstr(i1)*exphu
        !     deresdr=-2.0*vkrst2(i1)*diffr*vkrstr(i1)*exphu
            eres=eres+erh
            drda(1)=dx/rr
            drda(2)=dy/rr
            drda(3)=dz/rr
            do k1=1,3
                d(k1,ih1)=d(k1,ih1)+deresdr*drda(k1)
                d(k1,ih2)=d(k1,ih2)-deresdr*drda(k1)
            end do
        end if
    end do
          
!*********************************************************************
!                                                                    *
!     Calculate equivalent distance restraint energy                 *
!                                                                    *
!*********************************************************************
    do i1=1,neqdis
        ih1=ieqdis(i1,1)
        ih2=ieqdis(i1,2)
        ih3=ieqdis(i1,3)
        ih4=ieqdis(i1,4)
        call dista2(ih1,ih2,rr1,dx1,dy1,dz1)
        call dista2(ih3,ih4,rr2,dx2,dy2,dz2)
        diffr=rr1-rr2
        exphu=exp(-vkeqd2(i1)*(diffr*diffr))
        erh=vkeqd1(i1)*(1.0d0-exphu)
        deresdr1=2.0d0*vkeqd2(i1)*diffr*vkeqd1(i1)*exphu
        deresdr2=-2.0d0*vkeqd2(i1)*diffr*vkeqd1(i1)*exphu
        eres=eres+erh
        drda(1)=dx1/rr1
        drda(2)=dy1/rr1
        drda(3)=dz1/rr1
        do k1=1,3
            d(k1,ih1)=d(k1,ih1)+deresdr1*drda(k1)
            d(k1,ih2)=d(k1,ih2)-deresdr1*drda(k1)
        end do
        drda(1)=dx2/rr2
        drda(2)=dy2/rr2
        drda(3)=dz2/rr2
        do k1=1,3
            d(k1,ih3)=d(k1,ih3)+deresdr2*drda(k1)
            d(k1,ih4)=d(k1,ih4)-deresdr2*drda(k1)
        end do
    end do
!*********************************************************************
!                                                                    *
!     Calculate angle restraint energy                               *
!                                                                    *
!*********************************************************************
    do i1=1,nrestrav
        j(1)=irstrav(i1,1)
        j(2)=irstrav(i1,2)
        j(3)=irstrav(i1,3)
        ittr=0
    !     do i2=1,nval
    !     if (j(1).eq.iv(i2,2).and.j(2).eq.iv(i2,3).and.j(3).eq.iv(i2,4))
    !    $ittr=i2
    !     end do
    !     if (ittr.eq.0) stop 'Wrong valence angle restraint'
        call calvalres(j(1),j(2),j(3),arg,hr,dhrdc,dargdc)
        vaval=hr*rdndgr
        diffv=-(vaval-vrstra(i1))*dgrrdn
        exphu=exp(-vkr2v(i1)*(diffv*diffv))
        erh=vkrv(i1)*(1.0d0-exphu)
        deresdv=-2.0d0*vkr2v(i1)*diffv*vkrv(i1)*exphu
        eres=eres+erh
        do k1=1,3
            do k2=1,3
                d(k1,j(k2))=d(k1,j(k2))+deresdv*dhrdc(k1,k2)
            end do
        end do

    end do

!*********************************************************************
!                                                                    *
!     Calculate torsion restraint energy                             *
!                                                                    *
!*********************************************************************
    do i1=1,nrestrat
        j(1)=irstrat(i1,1)
        j(2)=irstrat(i1,2)
        j(3)=irstrat(i1,3)
        j(4)=irstrat(i1,4)
        ittr=0
        do i2=1,ntor
            if (j(1) == it(i2,2) .AND. j(2) == it(i2,3) .AND. j(3) == it(i2,4) &
             .AND. j(4) == it(i2,5)) ittr=i2
            if (j(4) == it(i2,2) .AND. j(3) == it(i2,3) .AND. j(2) == it(i2,4) &
             .AND. j(1) == it(i2,5)) ittr=i2
        end do
        if (ittr == 0) then
            write (*,*)'Wrong torsion restraint'
            write (*,*)i1,j(1),j(2),j(3),j(4)
            stop 'Wrong torsion restraint'
        end if
        vtor=thg(ittr)
        difft=-(vtor-trstra(i1))*dgrrdn
        exphu=exp(-vkr2t(i1)*(difft*difft))
        erh=vkrt(i1)*(1.0d0-exphu)
        deresdt=2.0d0*vkr2t(i1)*difft*vkrt(i1)*exphu
        if (vtor < zero) deresdt=-deresdt
        eres=eres+erh
        do k1=1,3
            do k2=1,4
                d(k1,j(k2))=d(k1,j(k2))+deresdt*dargtdc(ittr,k1,k2)
            end do
        end do

    end do
!*********************************************************************
!                                                                    *
!     Calculate mass centre restraint energy                         *
!                                                                    *
!*********************************************************************
    do i1=1,nrestram
        j1=irstram(i1,2)
        j2=irstram(i1,3)
        j3=irstram(i1,4)
        j4=irstram(i1,5)
        kdir=irstram(i1,1)
        cmx1=0.0d0
        cmy1=0.0d0
        cmz1=0.0d0
        cmx2=0.0d0
        cmy2=0.0d0
        cmz2=0.0d0
        summas1=0.0d0
        summas2=0.0d0
        do i2=j1,j2
            cmx1=cmx1+c(i2,1)*xmasat(i2)
            cmy1=cmy1+c(i2,2)*xmasat(i2)
            cmz1=cmz1+c(i2,3)*xmasat(i2)
            summas1=summas1+xmasat(i2)
        end do
        cmx1=cmx1/summas1
        cmy1=cmy1/summas1
        cmz1=cmz1/summas1
        if (mdstep < 2) then
            rmstrax(i1)=cmx1
            rmstray(i1)=cmy1
            rmstraz(i1)=cmz1
        end if
        if (kdir <= 3) then
            do i2=j3,j4
                cmx2=cmx2+c(i2,1)*xmasat(i2)
                cmy2=cmy2+c(i2,2)*xmasat(i2)
                cmz2=cmz2+c(i2,3)*xmasat(i2)
                summas2=summas2+xmasat(i2)
            end do
            cmx2=cmx2/summas2
            cmy2=cmy2/summas2
            cmz2=cmz2/summas2
        end if
        if (kdir == 1) dist=cmx1-cmx2
        if (kdir == 2) dist=cmy1-cmy2
        if (kdir == 3) dist=cmz1-cmz2
        if (kdir == 4) then
            distx=cmx1-rmstrax(i1)
            disty=cmy1-rmstray(i1)
            distz=cmz1-rmstraz(i1)
            dist=sqrt(distx*distx+disty*disty+distz*distz)
        end if
        dismacen(i1)=dist
        dist=dist-rmstra1(i1)
        erh=rmstra2(i1)*dist*dist
        deresdr=2.0*dist*rmstra2(i1)
    !     exphu=exp(-rmstra3(i1)*(dist*dist))
    !     erh=rmstra2(i1)*(1.0-exphu)
    !     deresdr=2.0*rmstra3(i1)*dist*rmstra2(i1)*exphu
        eres=eres+erh
        if (kdir <= 3) then
            do i2=j1,j2
                d(kdir,i2)=d(kdir,i2)+deresdr*xmasat(i2)/summas1
            end do
            do i2=j3,j4
                d(kdir,i2)=d(kdir,i2)-deresdr*xmasat(i2)/summas2
            end do
        end if
        if (kdir == 4 .AND. mdstep > 5) then
            do i2=j1,j2
                d(1,i2)=d(1,i2)+deresdr*(distx/dist)*(xmasat(i2)/summas1)
                d(2,i2)=d(2,i2)+deresdr*(disty/dist)*(xmasat(i2)/summas1)
                d(3,i2)=d(3,i2)+deresdr*(distz/dist)*(xmasat(i2)/summas1)
            end do
        end if
    end do
!*********************************************************************
!                                                                    *
!     Calculate symmetry restraint energy                            *
!                                                                    *
!*********************************************************************
    if (nrestras > 0) then
    !*********************************************************************
    !                                                                    *
    !     Convert to fractional coordinates                              *
    !                                                                    *
    !*********************************************************************
    !     write (67,'(f12.8)')c(1,1)
        do i1=1,na
            fc(i1,1)=c(i1,1)/tm11
            fc(i1,2)=(c(i1,2)-tm21*fc(i1,1))/tm22
            fc(i1,3)=(c(i1,3)-tm31*fc(i1,1)-tm32*fc(i1,2))/tm33
        end do
    !*********************************************************************
    !                                                                    *
    !     Move system to origin                                          *
    !                                                                    *
    !*********************************************************************

        vxmin=1000.0d0
        vymin=1000.0d0
        vzmin=1000.0d0

        do i1=1,na
            if (fc(i1,1) < vxmin) vxmin=fc(i1,1)
            if (fc(i1,2) < vymin) vymin=fc(i1,2)
            if (fc(i1,3) < vzmin) vzmin=fc(i1,3)
        end do

        do i1=1,na
            fc(i1,1)=fc(i1,1)-vxmin
            fc(i1,2)=fc(i1,2)-vymin
            fc(i1,3)=fc(i1,3)-vzmin
        end do

        do i1=1,na
            fc(i1,1)=fc(i1,1)-int(fc(i1,1))
            fc(i1,2)=fc(i1,2)-int(fc(i1,2))
            fc(i1,3)=fc(i1,3)-int(fc(i1,3))
        end do

    !*********************************************************************
    !                                                                    *
    !     Convert back to cartesian coordinates                          *
    !                                                                    *
    !*********************************************************************
        do i2=1,na
            cs(i2,1)=fc(i2,1)*tm11
            cs(i2,2)=fc(i2,1)*tm21+fc(i2,2)*tm22
            cs(i2,3)=fc(i2,1)*tm31+fc(i2,2)*tm32+fc(i2,3)*tm33
        end do

    !     write (67,'(2f12.8)')c(1,1),cs(i1,1)
        do i1=1,nrestras
            isymel=0
            do i2=1,nsymel
                if (qrstras(i1) == qsymm(i2)) isymel=i2
            end do
            if (isymel == 0) then
                write (*,*)qrstras(i1),nrestras,irstras(1,1)
                write (*,*)'Unknown symmetry element'
                stop 'Unknown symmetry element'
            end if
        !*********************************************************************
        !                                                                    *
        !     Construct symmetry images                                      *
        !                                                                    *
        !*********************************************************************
            natsym=0
            do i2=1,nsymopt(isymel)
                do i3=1,irstras(i1,1)
                    iatsy=irstras(i1,1+i3)
                    csym(natsym+1,1)=vsymadd(isymel,i2,1)+ &
                    vsymmulx(isymel,i2,1)*fc(iatsy,1)+ &
                    vsymmulx(isymel,i2,2)*fc(iatsy,2)+ &
                    vsymmulx(isymel,i2,3)*fc(iatsy,3)
                    csym(natsym+1,2)=vsymadd(isymel,i2,2)+ &
                    vsymmuly(isymel,i2,1)*fc(iatsy,1)+ &
                    vsymmuly(isymel,i2,2)*fc(iatsy,2)+ &
                    vsymmuly(isymel,i2,3)*fc(iatsy,3)
                    csym(natsym+1,3)=vsymadd(isymel,i2,3)+ &
                    vsymmulz(isymel,i2,1)*fc(iatsy,1)+ &
                    vsymmulz(isymel,i2,2)*fc(iatsy,2)+ &
                    vsymmulz(isymel,i2,3)*fc(iatsy,3)
                    qatsym(natsym+1)=qa(iatsy)
                !*********************************************************************
                !                                                                    *
                !     Put image back in unit cell                                    *
                !                                                                    *
                !*********************************************************************
                    do i4=1,3
                        if (csym(natsym+1,i4) > 1.0d0) &
                        csym(natsym+1,i4)=csym(natsym+1,i4)-1.0d0
                        if (csym(natsym+1,i4) < zero) &
                        csym(natsym+1,i4)=csym(natsym+1,i4)+1.0d0
                        if (csym(natsym+1,i4) > 0.99d0) &
                        csym(natsym+1,i4)=csym(natsym+1,i4)-1.0d0
                    end do
                     
                !*********************************************************************
                !                                                                    *
                !     Check overlap with existing images                             *
                !                                                                    *
                !*********************************************************************
                    dismin=1000.00d0
                    do i4=1,natsym
                        dx=csym(natsym+1,1)-csym(i4,1)
                        dy=csym(natsym+1,2)-csym(i4,2)
                        dz=csym(natsym+1,3)-csym(i4,3)
                        dx=dx-int(dx)
                        dy=dy-int(dy)
                        dz=dz-int(dz)

                        dissym=sqrt(dx*dx+dy*dy+dz*dz)
                        if (dissym < dismin) dismin=dissym
                    end do
                    if (dismin > 0.10d0) then
                        natsym=natsym+1   !Accept new image
                    end if

                end do
            end do
        !*********************************************************************
        !                                                                    *
        !     Move system to origin                                          *
        !                                                                    *
        !*********************************************************************
        !     do i2=1,natsym
        !     write (67,'(a2,3f12.8)')qatsym(i2),
        !    $csym(i2,1)*tm11,csym(i2,2),csym(i2,3)
        !     end do
        !     write (67,*)'2'

            vxmin=1000.0d0
            vymin=1000.0d0
            vzmin=1000.0d0

            do i2=1,natsym
                if (csym(i2,1) < vxmin) vxmin=csym(i2,1)
                if (csym(i2,2) < vymin) vymin=csym(i2,2)
                if (csym(i2,3) < vzmin) vzmin=csym(i2,3)
            end do

            do i2=1,natsym
                csym(i2,1)=csym(i2,1)-vxmin
                csym(i2,2)=csym(i2,2)-vymin
                csym(i2,3)=csym(i2,3)-vzmin
            end do

            do i2=1,natsym
                csym(i2,1)=csym(i2,1)-int(csym(i2,1))
                csym(i2,2)=csym(i2,2)-int(csym(i2,2))
                csym(i2,3)=csym(i2,3)-int(csym(i2,3))
            end do

        !     do i2=1,natsym
        !     write (67,'(a2,3f12.8)')qatsym(i2),
        !    $csym(i2,1),csym(i2,2),csym(i2,3)
        !     end do
        !     write (67,*)'3'
        !*********************************************************************
        !                                                                    *
        !     Convert symmetric images to cartesian coordinates              *
        !                                                                    *
        !*********************************************************************
            do i2=1,natsym
                ccsym(i2,1)=csym(i2,1)*tm11
                ccsym(i2,2)=csym(i2,1)*tm21+csym(i2,2)*tm22
                ccsym(i2,3)=csym(i2,1)*tm31+csym(i2,2)*tm32+csym(i2,3)*tm33
            end do
        !*********************************************************************
        !                                                                    *
        !     Shift symmetric images back to original location               *
        !                                                                    *
        !*********************************************************************
        !     dxs=ccsym(irstras(i1,2),1)-cs(irstras(i1,2),1)
        !     dys=ccsym(irstras(i1,2),2)-cs(irstras(i1,2),2)
        !     dzs=ccsym(irstras(i1,2),3)-cs(irstras(i1,2),3)

        !     do i2=1,natsym
        !     ccsym(i2,1)=ccsym(i2,1)-dxs
        !     ccsym(i2,2)=ccsym(i2,2)-dys
        !     ccsym(i2,3)=ccsym(i2,3)-dzs
        !     end do

        !     write (67,'(3f12.8)')c(1,1),cs(i1,1),csym(1,1)
        !     do i2=1,natsym
        !     write (67,'(a2,3f12.8)')qatsym(i2),
        !    $ccsym(i2,1),ccsym(i2,2),ccsym(i2,3)
        !     end do
        !     write (67,*)
        !     do i2=1,na
        !     write (67,'(a2,3f12.8)')qa(i2),
        !    $cs(i2,1),cs(i2,2),cs(i2,3)
        !     end do
        !     write (67,*)
        !*********************************************************************
        !                                                                    *
        !     Find closest symmetric image to coordinates                    *
        !     Calculate restraint energy and force                           *
        !                                                                    *
        !*********************************************************************
            sysdissum(i1)=zero
            do i2=1,na
                distmin=1000.00d0
                nclose=0
                do i3=1,natsym
                    disx=cs(i2,1)-ccsym(i3,1)
                    disy=cs(i2,2)-ccsym(i3,2)
                    disz=cs(i2,3)-ccsym(i3,3)
                    dist=sqrt(disx*disx+disy*disy+disz*disz)
                    if (dist < distmin) then
                        nclose=i3
                        dxclose=disx
                        dyclose=disy
                        dzclose=disz
                        distmin=dist
                    end if
                end do

                sysdissum(i1)=sysdissum(i1)+distmin
                exphu=exp(-vksym2(i1)*(distmin*distmin))
                erh=vksym1(i1)*(1.0d0-exphu)
                deresdr=2.0d0*vksym2(i1)*distmin*vksym1(i1)*exphu
                eres=eres+erh
            !     write (67,*)i2,nclose,distmin,erh,eres,deresdr
                drda(1)=dxclose
                drda(2)=dyclose
                drda(3)=dzclose
                do k1=1,3
                    d(k1,i2)=d(k1,i2)+deresdr*drda(k1)
                end do

            end do


        end do

        do i1=1,na
            do i2=1,3
                c(i1,i2)=cs(i1,i2)
            end do
        end do

    end if
!*********************************************************************
!                                                                    *
!     Calculate morphing energy                                      *
!                                                                    *
!*********************************************************************
    if (imorph == 1) then
        distot=zero
        do i1=1,na
            dmx=c(i1,1)-cmo(i1,1)
            dmy=c(i1,2)-cmo(i1,2)
            dmz=c(i1,3)-cmo(i1,3)
            dism=sqrt(dmx*dmx+dmy*dmy+dmz*dmz)
            distot=distot+dism
        !     exphu=exp(-vmo2(i1)*(dism*dism))
        !     erh=vmo1(i1)*(1.0-exphu)
            erh=vmo1(i1)*dism
            eres=eres+erh
        !     deresddis=2.0*vmo2(i1)*dism*vmo1(i1)*exphu
            deresddis=vmo1(i1)
            drda1=dmx/dism
            drda2=dmy/dism
            drda3=dmz/dism
            d(1,i1)=d(1,i1)+deresddis*drda1
            d(2,i1)=d(2,i1)+deresddis*drda2
            d(3,i1)=d(3,i1)+deresddis*drda3
        end do
        write (65,'(i6,6f12.4)')mdstep,distot,eres
          
    end if
          
          
    return
    end subroutine restraint
!*********************************************************************
!*********************************************************************

    subroutine piston

!*********************************************************************
    include 'cbka.blk'
!     dimension drda(3),j(4),dhrdc(3,3),dargdc(3,3)
!*********************************************************************
!                                                                    *
!     Calculate interactions with flat piston                        *
!     edeep:  Piston magnitude U=8*Edeep/(Rdeep+R^6)                 *
!     rdeep:  1:barrier equals Edeep at distance R=0                 *
!             0:infinite barrier at distance R=0                     *
!     pshft:  piston shift from left side of supercell               *
!     azm:    Location left side of supercell                        *
!     rcut:   Piston width                                           *
!     ipdir:  1: piston in x;2: piston in y;3: piston in z-direction *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In piston'
        call timer(65)
        close (65)
    end if
         
    epist=zero
    azm=zero
    if (icentr == 1) then
        azm=zero
    else if (icentr == 2) then
        if (ipdir == 1) azm= -axiss(1)/2.0d0
        if (ipdir == 2) azm= -axiss(2)/2.0d0
        if (ipdir == 3) azm= -axiss(3)/2.0d0
    end if

    apist0=azm+pshft
    apist1=azm+pshft+rcut
    edeep6=edeep/6.0d0

!     write (64,*)mdstep,edeep6,ipdir,apist1,apist0
    do i1=1,na
        episth=zero
        depistdr=zero
        if (c(i1,ipdir) < apist1) then
            dist=c(i1,ipdir)-apist0
            if (dist > zero) then
                rcube=dist*dist*dist
                episth=2.0d0*edeep6/(rcube+rdeep)
                depistdr=-dist*dist*edeep/((rcube+rdeep)*(rcube+rdeep))
            end if
        end if
          
        depistdc=depistdr
        epist=epist+episth
        d(ipdir,i1)=d(ipdir,i1)+depistdc

    end do

    if (mod(mdstep,nrep1) == 0) then
        open (81,file='fort.81',status='unknown',position='append')
        write (81,100)mdstep,epist
        close (81)
    end if

    return
    100 format (i8,f20.5)
    end subroutine piston
!*******************************************************************
!*******************************************************************

    subroutine calvalres (ja1,ja2,ja3,arg,hr,dhrdc,dargdc)

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),b(3),j(3),dradc(3,3),drbdc(3,3),dtdc(3,3), &
    dargdc(3,3),dndc(3,3),dadc(3),dbdc(3),dhrdc(3,3)
!*********************************************************************
!                                                                    *
!     Calculate valency angles and their derivatives to cartesian    *
!     coordinates  for restraint calculations                        *
!                                                                    *
!*********************************************************************
!     if (ndebug.eq.1) then
!     open (65,file='fort.65',status='unknown',position='append')
!     write (65,*) 'In calvalres'
!     call timer(65)
!     close (65)
!     end if

    dadc(1)=-1.0D0
    dadc(2)=1.0D0
    dadc(3)=0.0D0
    dbdc(1)=0.0D0
    dbdc(2)=1.0D0
    dbdc(3)=-1.0D0
    do k1=1,3
        do k2=1,3
            dradc(k1,k2)=0.0D0
            drbdc(k1,k2)=0.0D0
        end do
    end do
!*********************************************************************
!                                                                    *
!     Determine valency angle                                        *
!                                                                    *
!*********************************************************************
    call dista2(ja1,ja2,rla,dx1,dy1,dz1)
    call dista2(ja2,ja3,rlb,dx2,dy2,dz2)
     
    a(1)=-dx1
    a(2)=-dy1
    a(3)=-dz1
    b(1)=dx2
    b(2)=dy2
    b(3)=dz2
    poem=rla*rlb
    tel=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    arg=tel/poem
    arg2=arg*arg
    s1ma22=1.0D0-arg2
    if (s1ma22 < 1.0d-20) s1ma22=1.0d-20
    s1ma2=sqrt(s1ma22)
    if (arg > 1.0D0) arg=1.0D0
    if (arg < -1.0D0) arg=-1.0D0
    hr=acos(arg)
!*********************************************************************
!                                                                    *
!     Calculate derivative valency angle to cartesian coordinates    *
!                                                                    *
!*********************************************************************
    do k1=1,3
        dradc(k1,1)=-a(k1)/rla
        dradc(k1,2)=a(k1)/rla
    end do

    do k1=1,3
        drbdc(k1,2)=b(k1)/rlb
        drbdc(k1,3)=-b(k1)/rlb
    end do

    do k1=1,3
        do k2=1,3
            dndc(k1,k2)=rla*drbdc(k1,k2)+rlb*dradc(k1,k2)
            dtdc(k1,k2)=a(k1)*dbdc(k2)+b(k1)*dadc(k2)
            dargdc(k1,k2)=(dtdc(k1,k2)-arg*dndc(k1,k2))/poem
            dhrdc(k1,k2)=-dargdc(k1,k2)/s1ma2
        end do
    end do
          
    10 continue

    return
    end subroutine calvalres
!*********************************************************************
!*******************************************************************

    subroutine calvalhb (ja1,ja2,ja3,ix,iy,iz,arg,hr,dhrdc,dargdc)

!*********************************************************************
    include 'cbka.blk'
    dimension a(3),b(3),j(3),dradc(3,3),drbdc(3,3),dtdc(3,3), &
    dargdc(3,3),dndc(3,3),dadc(3),dbdc(3),dhrdc(3,3)
!*********************************************************************
!                                                                    *
!     Calculate valency angles and their derivatives to cartesian    *
!     coordinates  for hydrogen bond calculations                    *
!                                                                    *
!*********************************************************************
!     if (ndebug.eq.1) then
!     open (65,file='fort.65',status='unknown',position='append')
!     write (65,*) 'In calvalhb'
!     call timer(65)
!     close (65)
!     end if

    dadc(1)=-1.0d0
    dadc(2)=1.0d0
    dadc(3)=0.0d0
    dbdc(1)=0.0d0
    dbdc(2)=1.0d0
    dbdc(3)=-1.0d0
    do k1=1,3
        do k2=1,3
            dradc(k1,k2)=0.0d0
            drbdc(k1,k2)=0.0d0
        end do
    end do
!*********************************************************************
!                                                                    *
!     Determine valency angle                                        *
!                                                                    *
!*********************************************************************
    call dista2(ja1,ja2,rla,dx1,dy1,dz1)
    dx2=c(ja2,1)-c(ja3,1)+ix*tm11
    dy2=c(ja2,2)-c(ja3,2)+ix*tm21+iy*tm22
    dz2=c(ja2,3)-c(ja3,3)+ix*tm31+iy*tm32+iz*tm33
    rlb=sqrt(dx2*dx2+dy2*dy2+dz2*dz2)

    a(1)=-dx1
    a(2)=-dy1
    a(3)=-dz1
    b(1)=dx2
    b(2)=dy2
    b(3)=dz2
    poem=rla*rlb
    tel=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    arg=tel/poem
    arg2=arg*arg
    s1ma22=1.0d0-arg2
    if (s1ma22 < 1.0d-20) s1ma22=1.0d-20
    s1ma2=sqrt(s1ma22)
    if (arg > 1.0d0) arg=1.0d0
    if (arg < -1.0d0) arg=-1.0d0
    hr=acos(arg)
!*********************************************************************
!                                                                    *
!     Calculate derivative valency angle to cartesian coordinates    *
!                                                                    *
!*********************************************************************
    do k1=1,3
        dradc(k1,1)=-a(k1)/rla
        dradc(k1,2)=a(k1)/rla
    end do

    do k1=1,3
        drbdc(k1,2)=b(k1)/rlb
        drbdc(k1,3)=-b(k1)/rlb
    end do

    do k1=1,3
        do k2=1,3
            dndc(k1,k2)=rla*drbdc(k1,k2)+rlb*dradc(k1,k2)
            dtdc(k1,k2)=a(k1)*dbdc(k2)+b(k1)*dadc(k2)
            dargdc(k1,k2)=(dtdc(k1,k2)-arg*dndc(k1,k2))/poem
            dhrdc(k1,k2)=-dargdc(k1,k2)/s1ma2
        end do
    end do
          
    10 continue

    return
    end subroutine calvalhb
!*********************************************************************
!*********************************************************************

    subroutine caltor(ja1,ja2,ja3,ja4,ht)

!*********************************************************************
    include 'cbka.blk'
    DIMENSION  A(3),DRDA(3),DADC(4),DRADC(3,4),DRBDC(3,4), &
    DRCDC(3,4),DHDDC(3,4),DHEDC(3,4),DRVDC(3,4),DTDC(3,4), &
    DNDC(3,4)
    dimension j(4),dvdc1(3,3),dargdc1(3,3),dvdc2(3,3),dargdc2(3,3)
!*********************************************************************
!                                                                    *
!     Calculate torsion angle (for internal coordinates output)      *
!                                                                    *
!*********************************************************************
    if (ndebug == 1) then
        open (65,file='fort.65',status='unknown',position='append')
        write (65,*) 'In caltor'
        call timer(65)
        close (65)
    end if
    do k1=1,3
        do k2=1,4
            dhddc(k1,k2)=0.0d0
            dhedc(k1,k2)=0.0d0
            dradc(k1,k2)=0.0d0
            drbdc(k1,k2)=0.0d0
            drcdc(k1,k2)=0.0d0
        end do
    end do
    et=0.0d0
    eco=0.0d0
    dadc(1)=1.0d0
    dadc(2)=0.0d0
    dadc(3)=0.0d0
    dadc(4)=-1.0d0
    call dista2(ja1,ja2,rla,dx1,dy1,dz1)
    call dista2(ja2,ja3,rlb,dx2,dy2,dz2)
    call dista2(ja3,ja4,rlc,dx2,dy2,dz2)
    call dista2(ja1,ja4,r4,dx2,dy2,dz2)
    call calvalres(ja1,ja2,ja3,arg1,h1,dvdc1,dargdc1)
    call calvalres(ja2,ja3,ja4,arg2,h2,dvdc2,dargdc2)
!*********************************************************************
!                                                                    *
!     Determine torsion angle                                        *
!                                                                    *
!*********************************************************************
    d142=r4*r4
    coshd=cos(h1)
    coshe=cos(h2)
    sinhd=sin(h1)
    sinhe=sin(h2)
    poem=2.0d0*rla*rlc*sinhd*sinhe
    poem2=poem*poem
    tel=rla*rla+rlb*rlb+rlc*rlc-d142-2.0d0*(rla*rlb*coshd-rla*rlc* &
    coshd*coshe+rlb*rlc*coshe)
    arg=tel/poem
    if (arg > 1.0d0) arg=1.0d0
    if (arg < -1.0d0) arg=-1.0d0
    arg2=arg*arg
    ht=acos(arg)*rdndgr
    k1=ja1
    k2=ja2
    k3=ja3
    k4=ja4
    call dista2(k3,k2,dis,x3,y3,z3)
    y32z32=y3*y3+z3*z3
    wort1=sqrt(y32z32)+1.0d-6
    wort2=sqrt(y32z32+x3*x3)+1.0d-6
    sinalf=y3/wort1
    cosalf=z3/wort1
    sinbet=x3/wort2
    cosbet=wort1/wort2
    call dista2(k1,k2,dis,x1,y1,z1)
    x1=x1*cosbet-y1*sinalf*sinbet-z1*cosalf*sinbet
    y1=y1*cosalf-z1*sinalf
    wort3=sqrt(x1*x1+y1*y1)+1.0d-6
    singam=y1/wort3
    cosgam=x1/wort3
    call dista2(k4,k2,dis,x4,y4,z4)
    x4=x4*cosbet-y4*sinalf*sinbet-z4*cosalf*sinbet
    y4=y4*cosalf-z4*sinalf
    y4=x4*singam-y4*cosgam
    if (y4 > 0.0d0) ht=-ht
    if (ht < -179.999999999d0) ht=-179.999999999d0
    if (ht > 179.999999999d0) ht=179.999999999d0

    return
    end subroutine caltor
!*********************************************************************
!*********************************************************************
