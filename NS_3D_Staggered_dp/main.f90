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

    !write (*,*) "Choose:"
    !write (*,*) "(1): big_periodic_3D"
    !write (*,*) "(2): re_simulation"
    !write (*,*) "(3): both"
    !read (*,'(i1)') i

    !if (i==1) then
        !call big_periodic_3D
    !end if
    !
    !if (i==2) then
        !call re_simulation
    !end if
    call test_subroutine

    write (*,*) "Press anykey to continue..."
    read (*,*)

    end program main

