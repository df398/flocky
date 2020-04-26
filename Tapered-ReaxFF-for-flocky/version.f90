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
!     - Numerically stable sbo2 formulation in valence angles (2020) *                                                               *
!                                                                    *
!*********************************************************************
    subroutine version
    !include 'cbka.blk'
    !character(20) :: qhulp
    !qhulp='Compiled on:'
    !write (*,100)
    !write (*,110)qhulp
    return
    !100 format ('ReaxFF version 2.0')
    !110 format(a20,'Mon Sep 27 19:42:48 EDT 2010')
    end subroutine version
