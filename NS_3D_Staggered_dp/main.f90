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


    program main

    use NS_3D_Staggered_dp

    implicit none
    include 'omp_lib.h'

    integer :: i



    call big_periodic_3D

    !call re_simulation
    
    !call re_simulation_sub_tstep

    !call test_subroutine

    write (*,*) "Press anykey to continue..."
    !read (*,*)

    end program main

