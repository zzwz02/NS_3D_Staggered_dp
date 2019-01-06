    module csr_mod
    
    implicit none
    
    type csr
        real(8), dimension (:), allocatable :: value
        integer, dimension (:), allocatable :: ja, ia
        integer, dimension (2) :: shape
    contains
    procedure, nopass :: csr_init
    procedure, public :: csr_check_fail
    procedure, nopass :: add
    end type csr

    interface csr
    !module procedure dencoo
    end interface

    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(csr) function csr_init(csr_value, csr_ja, csr_ia, csr_shape)

    implicit none
    real(8), dimension (:), intent(in) :: csr_value
    integer, dimension (:), intent(in) :: csr_ja, csr_ia
    integer, dimension (2), intent(in) :: csr_shape

    integer :: nnz=0, idx=0, i, j

    nnz=count(csr_value/=0)
    idx=0

    csr_init%shape=csr_shape
    allocate(csr_init%value(size(csr_value)), csr_init%ja(size(csr_ja)), csr_init%ia(size(csr_ia)))
    csr_init%value=csr_value
    csr_init%ja=csr_ja
    csr_init%ia=csr_ia;

    if (csr_init%csr_check_fail()) then
        print *, "CSR sparse matrix size is incorrect. (SUBROUTINE csr_init)"
        !stop(1)
    end if

    end function csr_init
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical function csr_check_fail(self)
    implicit none
    
    class(csr), intent(in) :: self
    
    if (size(self%value)/=size(self%ja) .or. size(self%ia)/=self%shape(1)+1 .or. self%ia(self%shape(1)+1)/=size(self%value)+1) then
        csr_check_fail = .true.
    else
        csr_check_fail = .false.
    end if
    
    end function csr_check_fail
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(csr) function add(a,b)
    implicit none
    
    class(csr), intent(in) :: a,b
    
    real(8), dimension(:), allocatable :: value
    integer, dimension(:), allocatable :: ja
    integer, dimension(:), allocatable :: ia
    integer :: info=-111
    
    allocate( value(size(a%value)+size(b%value)), ja(size(a%value)+size(b%value)), ia(a%shape(1)+1) )
    value=0; ja=0; ia=0;
    
    call mkl_dcsradd("N", 0, 0, a%shape(1), a%shape(2), a%value, a%ja, a%ia, 1.0d0, b%value, b%ja, b%ia, value, ja, ia, int(2e9), info)
    
    print '("      nnz: " (I8))', maxval(ia)
    add%value=value(1:ia(a%shape(1)+1)-1)
    add%ja=ja(1:ia(a%shape(1)+1)-1)
    add%ia=ia
    add%shape=a%shape
    
    end function add
    
    end module csr_mod