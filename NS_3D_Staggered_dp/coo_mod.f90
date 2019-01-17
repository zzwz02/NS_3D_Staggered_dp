    module coo_mod
    
    use csr_mod
    implicit none
    
    type coo
        real(8), dimension (:), allocatable :: value
        integer, dimension (:), allocatable :: row, col
        integer, dimension (2) :: shape
    contains
    procedure, nopass :: coo_init_or_clean
    procedure :: coo_check_fail
    procedure :: to_den
    !procedure :: reshape
    procedure, nopass :: kron
    procedure :: to_csr
    procedure, nopass :: from_csr
    end type coo

    interface coo
    module procedure dencoo
    end interface

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(coo) function coo_init_or_clean(coo_value, coo_row, coo_col, coo_shape)

    implicit none
    real(8), dimension (:), intent(in) :: coo_value
    integer, dimension (:), intent(in) :: coo_row, coo_col
    integer, dimension (2), intent(in) :: coo_shape

    integer :: nnz=0, idx=0, i, j

    nnz=count(coo_value/=0)
    idx=0

    coo_init_or_clean%shape=coo_shape
    allocate(coo_init_or_clean%value(nnz), coo_init_or_clean%row(nnz), coo_init_or_clean%col(nnz))
    coo_init_or_clean%value=0; coo_init_or_clean%row=0; coo_init_or_clean%col=0;

    if (nnz/=size(coo_value)) then
        do i=1, size(coo_value)
            if (coo_value(i)/=0) then
                idx=idx+1
                coo_init_or_clean%value(idx)=coo_value(i)
                coo_init_or_clean%row(idx)=coo_row(i)
                coo_init_or_clean%col(idx)=coo_col(i)
            end if
        end do
    else
        coo_init_or_clean%value=coo_value; coo_init_or_clean%row=coo_row; coo_init_or_clean%col=coo_col;
    end if

    if (coo_init_or_clean%coo_check_fail()) then
        print *, "COO sparse matrix size is incorrect. (SUBROUTINE coo_init_or_clean)"
        !stop(1)
    end if

    end function coo_init_or_clean

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(coo) function dencoo(mat)

    implicit none
    real(8), dimension (:,:), intent(in) :: mat

    integer :: nnz=0, idx=0, i, j

    nnz=count(mat/=0)
    idx=0

    dencoo%shape=shape(mat)
    allocate(dencoo%value(nnz), dencoo%row(nnz), dencoo%col(nnz))
    dencoo%value=0; dencoo%row=0; dencoo%col=0;

    do j=1, dencoo%shape(2)
        do i=1, dencoo%shape(1)
            if (mat(i,j)/=0) then
                idx=idx+1
                dencoo%value(idx)=mat(i,j)
                dencoo%row(idx)=i
                dencoo%col(idx)=j
            end if
        end do
    end do

    if (dencoo%coo_check_fail()) then
        print *, "COO sparse matrix size is incorrect. (SUBROUTINE dencoo)"
        !stop(1)
    end if

    end function dencoo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Logical function coo_check_fail(self)
    implicit none

    class(coo), intent(in) :: self

    if (minval(self%row)<1 .or. minval(self%col)<1 .or. maxval(self%row)>self%shape(1) .or. maxval(self%col)>self%shape(2) .or. &
        size(self%value)/=size(self%row) .or. size(self%value)/=size(self%col)) then
        coo_check_fail = .true.
    else
        coo_check_fail = .false.
    end if

    end function coo_check_fail

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(csr) function to_csr(self)
    implicit none
    
    class(coo), intent(in) :: self
    
    class(coo), allocatable :: self0
    integer :: nnz=0, idx=0, i, j
    real(8), dimension(:), allocatable :: value
    integer, dimension(:), allocatable :: ja
    integer, dimension(:), allocatable :: ia
    integer :: info=-111, job(8)=0
    
    self0=coo_init_or_clean(self%value, self%row, self%col, self%shape)
    if (self0%coo_check_fail()) then
        print *, "COO sparse matrix size is incorrect. (SUBROUTINE to_csr)"
        !stop(1)
    end if
    
    nnz=size(self0%value)
    allocate( value(nnz), ja(nnz), ia(self0%shape(1)+1) )
    value=0; ja=0; ia=0;
    
    job= [2,1,1,0,1000000,0,0,0] ! COO to CSR
    call mkl_dcsrcoo(job, self0%shape(1), value, ja, ia, nnz, self0%value, self0%row, self0%col, info)
    
    !allocate( to_csr%value(size(value)), to_csr%ja(size(ja)), to_csr%ia(size(ia)) )
    to_csr=csr_init(value,ja,ia,self%shape)
    
    end function to_csr
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(coo) function from_csr(self)
    implicit none
    
    class(csr), intent(in) :: self
    
    class(coo), allocatable :: self0
    integer :: nnz=0, idx=0, i, j
    real(8), dimension(:), allocatable :: value
    integer, dimension(:), allocatable :: row, col
    integer :: info=-111, job(8)=0
    
    !self0=coo_init_or_clean(self%value, self%row, self%col, self%shape)
    if (self%csr_check_fail()) then
        print *, "CSR sparse matrix size is incorrect. (SUBROUTINE from_csr)"
        !stop(1)
    end if
    
    nnz=size(self%value)
    allocate( value(nnz), row(nnz), col(nnz) )
    value=0; row=0; col=0;
    
    job= [0,1,1,0,int(2e9),3,0,0] ! COO to CSR
    call mkl_dcsrcoo(job, self%shape(1), self%value, self%ja, self%ia, nnz, value, row, col, info)
    
    !allocate( to_csr%value(size(value)), to_csr%row(size(ja)), to_csr%ia(size(ia)) )
    from_csr=coo_init_or_clean(value,row,col,self%shape)
    
    end function from_csr
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function to_den(self) result(mat)
    implicit none

    class(coo), intent(in) :: self
    real(8), dimension (:,:), allocatable :: mat

    integer :: nnz=0, idx=0, i, j

    if (self%coo_check_fail()) then
        print *, "COO sparse matrix size is incorrect. (SUBROUTINE to_den)"
        !stop(1)
    end if
    allocate( mat(self%shape(1), self%shape(2)) )

    mat=0
    do i=1, size(self%value)
        mat(self%row(i),self%col(i))=self%value(i)
    end do

    end function to_den

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(coo) function kron(a, b)
    implicit none

    class(coo), intent(in) :: a, b
    !class(coo), allocatable, intent(out) :: c

    integer, dimension (2) :: c_shape
    real(8), dimension (:), allocatable :: c_coo
    integer, dimension (:), allocatable :: c_rowind, c_colind

    integer :: sizea, sizeb, i
    real(8)  ::  t1_coo(size(b%value), size(a%value))
    integer :: t1_row(size(b%value), size(a%value)), t1_col(size(b%value), size(a%value))


        if (a.coo_check_fail()) then
            print *, "COO sparse matrix size is incorrect. (SUBROUTINE kron, mat A)"
            !stop(1)
        end if

        if (a.coo_check_fail()) then
            print *, "COO sparse matrix size is incorrect. (SUBROUTINE kron, mat B)"
            !stop(1)
        end if

        sizea=size(a%value); sizeb=size(b%value);
        !A = coo_matrix(A)
        !output_shape = (A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])
        c_shape(1)=a%shape(1)*b%shape(1); c_shape(2)=a%shape(2)*b%shape(2);
        !
        !if A.nnz == 0 or B.nnz == 0:
        !# kronecker product is the zero matrix
        !return coo_matrix(output_shape)
        !
        !# expand entries of a into blocks
        !row = A.row.repeat(B.nnz)
        !col = A.col.repeat(B.nnz)
        !data = A.data.repeat(B.nnz)
        !
        !row *= B.shape[0]
        !col *= B.shape[1]

        allocate( c_coo(sizea*sizeb), c_rowind(sizea*sizeb), c_colind(sizea*sizeb) )

        do i=1,sizea
            c_rowind((i-1)*sizeb+1:(i-1)*sizeb+sizeb)=a%row(i)-1
            c_colind((i-1)*sizeb+1:(i-1)*sizeb+sizeb)=a%col(i)-1
            c_coo((i-1)*sizeb+1:(i-1)*sizeb+sizeb)=a%value(i)
        end do
        c_rowind=c_rowind*b%shape(1); c_colind=c_colind*b%shape(2);
        !
        !# increment block indices
        !row,col = row.reshape(-1,B.nnz),col.reshape(-1,B.nnz)
        !row += B.row
        !col += B.col
        !row,col = row.reshape(-1),col.reshape(-1)
        t1_row=reshape(c_rowind, [sizeb,sizea]); t1_col=reshape(c_colind, [sizeb,sizea]);
        do i=1,sizea
            t1_row(:,i)=t1_row(:,i)+b%row-1
            t1_col(:,i)=t1_col(:,i)+b%col-1
        end do
        c_rowind=reshape(t1_row+1, [sizea*sizeb]); c_colind=reshape(t1_col+1, [sizea*sizeb]);
        !
        !# compute block entries
        !data = data.reshape(-1,B.nnz) * B.data
        !data = data.reshape(-1)
        t1_coo=reshape(c_coo, [sizeb,sizea]);
        do i=1,sizea
            t1_coo(:,i)=t1_coo(:,i)*b%value
        end do
        c_coo=reshape(t1_coo, [sizea*sizeb]);

        kron=coo_init_or_clean(c_coo,c_rowind,c_colind,c_shape)
    
    end function kron

    end module coo_mod

