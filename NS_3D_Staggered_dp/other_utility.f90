    module other_utility

    implicit none

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine print_mat(mat)
    implicit none

    real(8), dimension (:,:), intent(in) :: mat

    integer :: i, j, mat_shape(2)

    mat_shape=shape(mat)
    do i=1,mat_shape(1)
        write(*,'(100g15.5)') ( mat(i,j), j=1,mat_shape(2) )
    enddo
    print *, "**************************************"


    end subroutine print_mat

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine print_spmat(mat)
    implicit none

    real(8), dimension (:,:), intent(in) :: mat
    !char(*) :: mat0

    integer :: i, j, mat_shape(2)

    mat_shape=shape(mat)
    !allocate(mat0(mat_shape(2)))

    do i=1,mat_shape(1)
        do j=1,mat_shape(2)
            if (mat(i,j)==0) then
                !mat0(j)=''
                write(*,'(".")', ADVANCE="NO")
            else
                !mat0(i)='*'
                write(*,'("*")', ADVANCE="NO")
            end if
        end do
        write(*,'(/)')
    end do

    !do i=1,mat_shape(1)
    !    write(*,'(100I1)') ( mat0(i,j), j=1,mat_shape(2) )
    !enddo
    print *, "**************************************"


    end subroutine print_spmat

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function unique_sort(vec_in) result(vec_out)
    implicit none
    integer :: i = 0, min_val, max_val
    integer, dimension(:), intent(in) :: vec_in
    integer, dimension(:), allocatable :: unique, vec_out

    allocate(unique(size(vec_in)))
    min_val = minval(vec_in)-1
    max_val = maxval(vec_in)
    do while (min_val<max_val)
        i = i+1
        min_val = minval(vec_in, mask=vec_in>min_val)
        unique(i) = min_val
    enddo
    !allocate(vec_out(i), source=unique(1:i))   !<-- Or, just use unique(1:i)
    vec_out=unique(1:i)
    !print "(10I5:)", vec_out
    end function unique_sort

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) function timediff(ic0, ic1)
    CHARACTER (LEN = 12), intent(in) :: ic0, ic1
    real(8) :: hh0, mm0, ss0, hh1, mm1, ss1

    read(ic0(1:2),'(f)') hh0
    read(ic0(3:4),'(f)') mm0
    read(ic0(5:10),'(f)') ss0

    read(ic1(1:2),'(f)') hh1
    read(ic1(3:4),'(f)') mm1
    read(ic1(5:10),'(f)') ss1

    ss0=hh0*3600+mm0*60+ss0
    ss1=hh1*3600+mm1*60+ss1
    if ( ss0<=ss1 ) then
        timediff=ss1-ss0
    else
        timediff=ss1-ss0+86400.0d0
    end if

    end function timediff
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) function mean(a)
    real(8), intent(in) :: a(:,:,:)

    mean=sum(a)!/size(a)
    
    end function mean
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) function rms(a)
    real(8), intent(in) :: a(:,:,:)

    rms=sqrt(sum(a*a)/size(a))
    
    end function rms

    end module other_utility

