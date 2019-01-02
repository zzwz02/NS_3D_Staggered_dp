    module sparse_pack

    implicit none
        
    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine den_kron(a_mat, b_mat, c_mat)
    implicit none

    real(8), dimension (:,:), intent(in) :: a_mat, b_mat
    real(8), dimension (:,:), allocatable, intent(out) :: c_mat

    integer, dimension (2) :: a_shape, b_shape, c_shape
    real(8), dimension (:), allocatable :: a_coo, b_coo, c_coo
    integer, dimension (:), allocatable :: a_rowind, a_colind, b_rowind, b_colind, c_rowind, c_colind

    call den_coo(a_mat, a_coo, a_rowind, a_colind, a_shape)
    call den_coo(b_mat, b_coo, b_rowind, b_colind, b_shape)
    call coo_kron(a_coo, a_rowind, a_colind, a_shape, b_coo, b_rowind, b_colind, b_shape, c_coo, c_rowind, c_colind, c_shape)
    call coo_den(c_mat, c_coo, c_rowind, c_colind, c_shape)

    end subroutine den_kron

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine den_coo(mat, a_coo, a_rowind, a_colind, a_shape)
    implicit none

    real(8), dimension (:,:), intent(in) :: mat
    integer, dimension (2), intent(out) :: a_shape
    real(8), dimension (:), allocatable, intent(out) :: a_coo
    integer, dimension (:), allocatable, intent(out) :: a_rowind, a_colind

    integer :: nnz=0, idx=0, i, j

    idx=0
    nnz=count(mat/=0)
    a_shape=shape(mat)
    allocate(a_coo(nnz), a_rowind(nnz), a_colind(nnz))

    a_coo=0; a_rowind=0; a_colind=0;
    do j=1, a_shape(2)
        do i=1, a_shape(1)
            if (mat(i,j)/=0) then
                idx=idx+1
                a_coo(idx)=mat(i,j)
                a_rowind(idx)=i
                a_colind(idx)=j
            end if
        end do
    end do

    end subroutine den_coo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coo_clean(a_coo, a_rowind, a_colind, b_coo, b_rowind, b_colind)
    implicit none

    real(8), dimension (:), intent(in) :: a_coo
    integer, dimension (:), intent(in) :: a_rowind, a_colind
    !integer, dimension (2), intent(inout) :: a_shape
    real(8), dimension (:), allocatable, intent(out) :: b_coo
    integer, dimension (:), allocatable, intent(out) :: b_rowind, b_colind

    integer :: nnz=0, idx=0, i

    nnz=count(a_coo/=0)
    if (nnz/=size(a_coo)) then
        
        allocate(b_coo(nnz), b_rowind(nnz), b_colind(nnz))

        b_coo=0; b_rowind=0; b_colind=0;
        do i=1, size(a_coo)
            if (a_coo(i)/=0) then
                idx=idx+1
                b_coo(idx)=a_coo(i)
                b_rowind(idx)=a_rowind(i)
                b_colind(idx)=a_colind(i)
            end if
        end do
    else
        allocate(b_coo(nnz), b_rowind(nnz), b_colind(nnz))
        b_coo=a_coo; b_rowind=a_rowind; b_colind=a_colind;
    end if

    end subroutine coo_clean

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coo_den(mat, a_coo, a_rowind, a_colind, a_shape)
    implicit none

    real(8), dimension (:,:), allocatable, intent(out) :: mat
    integer, dimension (2), intent(in) :: a_shape
    real(8), dimension (:), intent(in) :: a_coo
    integer, dimension (:), intent(in) :: a_rowind, a_colind

    integer :: nnz=0, idx=0, i, j

    call coo_check(a_coo, a_rowind, a_colind, a_shape)
    allocate( mat(a_shape(1), a_shape(2)) )

    mat=0
    do i=1, size(a_coo)
        mat(a_rowind(i),a_colind(i))=a_coo(i)
    end do

    end subroutine coo_den

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coo_reshape(a_coo, a_rowind, a_colind, old_shape, new_shape)
    implicit none

    integer, dimension (2), intent(in) :: old_shape, new_shape
    real(8), dimension (:), intent(inout) :: a_coo
    integer, dimension (:), intent(inout) :: a_rowind, a_colind

    integer :: flat_ind(size(a_coo))
    !integer :: nnz=0

    call coo_check(a_coo, a_rowind, a_colind, old_shape)
    if (old_shape(1)*old_shape(2)/=new_shape(1)*new_shape(2)) then
        print *, "To RESHAPE the number of elements must not change."
        stop (1)
    end if

    !nnz=0
    flat_ind=a_rowind+(a_colind-1)*old_shape(1)

    a_colind=ceiling(dble(flat_ind)/dble(new_shape(1)))
    a_rowind=mod(flat_ind-1,new_shape(1))+1

    end subroutine coo_reshape

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coo_check(a_coo, a_rowind, a_colind, a_shape)
    implicit none

    integer, dimension (2), intent(in) :: a_shape
    real(8), dimension (:), intent(in) :: a_coo
    integer, dimension (:), intent(in) :: a_rowind, a_colind

    if (minval(a_rowind)<1 .or. minval(a_colind)<1 .or. maxval(a_rowind)>a_shape(1) .or. maxval(a_colind)>a_shape(2) .or. size(a_coo)/=size(a_rowind) .or. size(a_coo)/=size(a_colind)) then
        print *, "COO sparse matrix size is incorrect."
        stop (1)
    end if

    end subroutine coo_check

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coo_kron(a_coo, a_rowind, a_colind, a_shape, b_coo, b_rowind, b_colind, b_shape, c_coo, c_rowind, c_colind, c_shape)
    implicit none

    integer, dimension (2), intent(in) :: a_shape, b_shape
    real(8), dimension (:), intent(in) :: a_coo, b_coo
    integer, dimension (:), intent(in) :: a_rowind, a_colind, b_rowind, b_colind

    integer, dimension (2), intent(out) :: c_shape
    real(8), dimension (:), allocatable, intent(out) :: c_coo
    integer, dimension (:), allocatable, intent(out) :: c_rowind, c_colind

    !integer :: flat_ind(size(acoo))
    integer :: sizea, sizeb, i
    real(8) :: a1_coo(size(a_coo)), b1_coo(size(b_coo))
    integer :: a1_rowind(size(a_coo)), a1_colind(size(a_coo)), b1_rowind(size(b_coo)), b1_colind(size(b_coo))
    real(8)  ::  t1_coo(size(b_coo),size(a_coo))
    integer :: t1_row(size(b_coo),size(a_coo)), t1_col(size(b_coo),size(a_coo))

    a1_coo=a_coo; b1_coo=b_coo;
    a1_rowind=a_rowind; a1_colind=a_colind; b1_rowind=b_rowind; b1_colind=b_colind;
    call coo_check(a1_coo, a1_rowind, a1_colind, a_shape)
    call coo_check(b1_coo, b1_rowind, b1_colind, b_shape)

    sizea=size(a1_coo); sizeb=size(b1_coo);
    a1_rowind=a1_rowind-1; a1_colind=a1_colind-1; b1_rowind=b1_rowind-1; b1_colind=b1_colind-1;
    !A = coo_matrix(A)
    !output_shape = (A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])
    c_shape(1)=a_shape(1)*b_shape(1); c_shape(2)=a_shape(2)*b_shape(2);
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
        c_rowind((i-1)*sizeb+1:(i-1)*sizeb+sizeb)=a1_rowind(i)
        c_colind((i-1)*sizeb+1:(i-1)*sizeb+sizeb)=a1_colind(i)
        c_coo((i-1)*sizeb+1:(i-1)*sizeb+sizeb)=a1_coo(i)
    end do
    c_rowind=c_rowind*b_shape(1); c_colind=c_colind*b_shape(2);
    !
    !# increment block indices
    !row,col = row.reshape(-1,B.nnz),col.reshape(-1,B.nnz)
    !row += B.row
    !col += B.col
    !row,col = row.reshape(-1),col.reshape(-1)
    t1_row=reshape(c_rowind, (/sizeb,sizea/)); t1_col=reshape(c_colind, (/sizeb,sizea/));
    do i=1,sizea
        t1_row(:,i)=t1_row(:,i)+b1_rowind
        t1_col(:,i)=t1_col(:,i)+b1_colind
    end do
    c_rowind=reshape(t1_row+1, (/sizea*sizeb/)); c_colind=reshape(t1_col+1, (/sizea*sizeb/));
    !
    !# compute block entries
    !data = data.reshape(-1,B.nnz) * B.data
    !data = data.reshape(-1)
    t1_coo=reshape(c_coo, (/sizeb,sizea/));
    do i=1,sizea
        t1_coo(:,i)=t1_coo(:,i)*b_coo
    end do
    c_coo=reshape(t1_coo, (/sizea*sizeb/));

    end subroutine coo_kron
    
    end module sparse_pack

