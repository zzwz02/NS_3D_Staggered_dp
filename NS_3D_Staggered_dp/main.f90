    !  NS_3D_Staggered_dp.f90
    !
    !  FUNCTIONS:
    !  NS_3D_Staggered_dp - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: NS_3D_Staggered_dp
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************
    include 'mkl_dfti.f90'

    program main

    use big_domain_3D

    implicit none
    include 'omp_lib.h'

    !call big_periodic_3D

    call re_simulation

    write (*,*) "Press anykey to continue..."
    read (*,*)

    end program main

