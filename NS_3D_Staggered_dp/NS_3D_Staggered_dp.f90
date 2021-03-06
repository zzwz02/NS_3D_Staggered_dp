    !include 'mkl_dfti.f90'
    !include 'mkl_trig_transforms.f90'
    include 'lapack.f90'
    module NS_3D_Staggered_dp

    use MKL_DFTI!, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R
    use mkl_trig_transforms
    USE lapack95
    !use f95_precision
    use FD_functions
    use NS_functions
    use coo_mod
    use csr_mod
    use other_utility
    !use ogpf
    use hdf5
    use h5lt
    use, intrinsic :: ieee_arithmetic

    implicit none
    
    include 'mkl_lapack.fi'
    include 'mkl_pardiso.fi'
    
    contains

    include 'big_periodic_3D.f90'
    include 're_simulation.f90'
    include 're_simulation_sub_tstep.f90'
    include 're_simulation_sub_tstep_interBC_combine.f90'
    include 'test_subroutine.f90'

    end module NS_3D_Staggered_dp

