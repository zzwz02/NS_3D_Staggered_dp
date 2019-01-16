    module FD_functions

    implicit none

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function avg_old(A, k, bc) result(B)
    implicit none
    integer, intent(in) :: k, bc
    real(8), dimension (:,:,:), intent(in) :: A
    real(8), dimension (:,:,:), allocatable :: B

    integer :: nx, ny, nz, i

    !system_clock
    REAL(8) :: system_clock_rate,t1,t2
    INTEGER :: c1,c2,cr,cm
    logical :: abc

    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    system_clock_rate = REAL(cr)

    !if (size(A,1)==1) A=transpose(A)
    
    nx=size(A,1)
    ny=size(A,2)
    nz=size(A,3)

    if (bc==1) then
        if (k==1) then
            allocate( B(nx,ny,nz) )
            !CALL CPU_TIME(t1)
            !CALL SYSTEM_CLOCK(c1)
            !$omp parallel do
            do i=1, nx-1
                B(i+1,:,:)=(A(i+1,:,:)+A(i,:,:))/2.0d0
            end do
            !$omp end parallel do
            !CALL CPU_TIME(t2)
            !CALL SYSTEM_CLOCK(c2)
            !print '("    omp wall time: ", F8.4, " second")', (c2-c1)/system_clock_rate
            !print '("    omp cpu time: ", F8.4, " second")', t2-t1
            B(1,:,:)=B(nx,:,:)
        else if (k==2) then
            allocate( B(nx,ny,nz) )
            !$omp parallel do
            do i=1, ny-1
                B(:,i+1,:)=(A(:,i+1,:)+A(:,i,:))/2.0d0
            end do
            !$omp end parallel do
            B(:,1,:)=B(:,ny,:)
        else if (k==3) then
            allocate( B(nx,ny,nz) )
            !$omp parallel do
            do i=1, nz-1
                B(:,:,i+1)=(A(:,:,i+1)+A(:,:,i))/2.0d0
            end do
            !$omp end parallel do
            B(:,:,1)=B(:,:,nz)
        end if
    else
        if (k==1) then
            allocate( B(nx-1,ny,nz) )
            !$omp parallel do
            do i=1, nx-1
                B(i,:,:)=(A(i+1,:,:)+A(i,:,:))/2.0d0
            end do
            !$omp end parallel do
        else if (k==2) then
            allocate( B(nx,ny-1,nz) )
            !$omp parallel do
            do i=1, ny-1
                B(:,i,:)=(A(:,i+1,:)+A(:,i,:))/2.0d0
            end do
            !$omp end parallel do
        else if (k==3) then
            allocate( B(nx,ny,nz-1) )
            !$omp parallel do
            do i=1, nz-1
                B(:,:,i)=(A(:,:,i+1)+A(:,:,i))/2.0d0
            end do
            !$omp end parallel do
        end if
    end if

    end function avg_old

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function diff_old(A, n, k) result(B)
    implicit none
    integer, intent(in) :: k, n
    real(8), dimension (:,:,:), intent(in) :: A
    real(8), dimension (:,:,:), allocatable :: B

    integer :: nx, ny, nz, i

    !if (size(A,1)==1) A=transpose(A)

    nx=size(A,1)
    ny=size(A,2)
    nz=size(A,3)

    if (k==1) then
        allocate( B(nx-1,ny,nz) )
        !$omp parallel do
        do i=1, nx-1
            B(i,:,:)=A(i+1,:,:)-A(i,:,:)
        end do
        !$omp end parallel do
    else if (k==2) then
        allocate( B(nx,ny-1,nz) )
        !$omp parallel do
        do i=1, ny-1
            B(:,i,:)=A(:,i+1,:)-A(:,i,:)
        end do
        !$omp end parallel do
    else if (k==3) then
        allocate( B(nx,ny,nz-1) )
        !$omp parallel do
        do i=1, nz-1
            B(:,:,i)=A(:,:,i+1)-A(:,:,i)
        end do
        !$omp end parallel do
    end if

    end function diff_old

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function diff2_old(A, k, bc) result(B)
    implicit none
    integer, intent(in) :: k, bc
    real(8), dimension (:,:,:), intent(in) :: A
    real(8), dimension (:,:,:), allocatable :: B
    real(8), dimension (:,:,:), allocatable :: A1

    integer :: nx, ny, nz, i

    !if (size(A,1)==1) A=transpose(A)

    nx=size(A,1)
    ny=size(A,2)
    nz=size(A,3)

    if (bc==1) then
        if (k==1) then
            allocate( A1(nx+1,ny,nz) )
            A1(1,:,:)=A(nx-1,:,:)
            !$omp parallel do
            do i=2, nx+1
                A1(i,:,:)=A(i-1,:,:)
            end do
            !$omp end parallel do

            allocate( B(nx-1,ny,nz) )
            !$omp parallel do
            do i=1, nx-1
                B(i,:,:)=A1(i,:,:)-2.0d0*A1(i+1,:,:)+A1(i+2,:,:)
            end do
            !$omp end parallel do
        else if (k==2) then
            allocate( A1(nx,ny+1,nz) )
            A1(:,1,:)=A(:,ny-1,:)
            !$omp parallel do
            do i=2, ny+1
                A1(:,i,:)=A(:,i-1,:)
            end do
            !$omp end parallel do

            allocate( B(nx,ny-1,nz) )
            !$omp parallel do
            do i=1, ny-1
                B(:,i,:)=A1(:,i,:)-2.0d0*A1(:,i+1,:)+A1(:,i+2,:)
            end do
            !$omp end parallel do
        else if (k==3) then
            allocate( A1(nx,ny,nz+1) )
            A1(:,:,1)=A(:,:,nz-1)
            !$omp parallel do
            do i=2, nz+1
                A1(:,:,i)=A(:,:,i-1)
            end do
            !$omp end parallel do

            allocate( B(nx,ny,nz-1) )
            !$omp parallel do
            do i=1, nz-1
                B(:,:,i)=A1(:,:,i)-2.0d0*A1(:,:,i+1)+A1(:,:,i+2)
            end do
            !$omp end parallel do
        end if
    else
        if (k==1) then
            allocate( B(nx-2,ny,nz) )
            !$omp parallel do
            do i=1, nx-2
                B(i,:,:)=A(i,:,:)-2.0d0*A(i+1,:,:)+A(i+2,:,:)
            end do
            !$omp end parallel do
        else if (k==2) then
            allocate( B(nx,ny-2,nz) )
            !$omp parallel do
            do i=1, ny-2
                B(:,i,:)=A(:,i,:)-2.0d0*A(:,i+1,:)+A(:,i+2,:)
            end do
            !$omp end parallel do
        else if (k==3) then
            allocate( B(nx,ny,nz-2) )
            !$omp parallel do
            do i=1, nz-2
                B(:,:,i)=A(:,:,i)-2.0d0*A(:,:,i+1)+A(:,:,i+2)
            end do
            !$omp end parallel do
        end if
    end if

    end function diff2_old

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function avg(A, k, bc) result(C)
    implicit none
    integer, intent(in) :: k, bc
    real(8), dimension (:,:,:), intent(in) :: A
    real(8), target :: B(size(A,1),size(A,2),size(A,3))
    real(8), pointer :: C(:,:,:)

    integer :: nx, ny, nz, i, j, kk

    !system_clock
    !REAL(8) :: system_clock_rate,t1,t2
    !INTEGER :: c1,c2,cr,cm

    ! First initialize the system_clock
    !CALL system_clock(count_rate=cr)
    !CALL system_clock(count_max=cm)
    !system_clock_rate = REAL(cr)

    !if (size(A,1)==1) A=transpose(A)

    nx=size(A,1)
    ny=size(A,2)
    nz=size(A,3)

    if (bc==1) then
        if (k==1) then
            !CALL CPU_TIME(t1)
            !CALL SYSTEM_CLOCK(c1)
            !$omp parallel do
            do kk=1, nz
                do j=1, ny
                    do i=1, nx
                        if (i==1) then
                            B(1,j,kk)=(A(nx,j,kk)+A(nx-1,j,kk))/2.0d0
                        else
                            B(i,j,kk)=(A(i,j,kk)+A(i-1,j,kk))/2.0d0
                        end if
                    end do

                end do
            end do
            !$omp end parallel do
            !CALL CPU_TIME(t2)
            !CALL SYSTEM_CLOCK(c2)
            !print '("    omp wall time: ", F8.4, " second")', (c2-c1)/system_clock_rate
            !print '("    omp cpu time: ", F8.4, " second")', t2-t1
            C => B(:,:,:)
        else if (k==2) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny
                    do i=1, nx
                        if (j==1) then
                            B(i,1,kk)=(A(i,ny,kk)+A(i,ny-1,kk))/2.0d0
                        else
                            B(i,j,kk)=(A(i,j,kk)+A(i,j-1,kk))/2.0d0
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,:,:)
        else if (k==3) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny
                    do i=1, nx
                        if (kk==1) then
                            B(i,j,1)=(A(i,j,nz)+A(i,j,nz-1))/2.0d0
                        else
                            B(i,j,kk)=(A(i,j,kk)+A(i,j,kk-1))/2.0d0
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,:,:)
        end if
    else
        if (k==1) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny
                    do i=1, nx-1
                        B(i,j,kk)=(A(i+1,j,kk)+A(i,j,kk))/2.0d0
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(1:nx-1,:,:)
        else if (k==2) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny-1
                    do i=1, nx
                        B(i,j,kk)=(A(i,j+1,kk)+A(i,j,kk))/2.0d0
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,1:ny-1,:)
        else if (k==3) then
            !$omp parallel do
            do kk=1, nz-1
                do j=1, ny
                    do i=1, nx
                        B(i,j,kk)=(A(i,j,kk+1)+A(i,j,kk))/2.0d0
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,:,1:nz-1)
        end if
    end if

    end function avg

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function diff(A, n, k) result(C)
    implicit none
    integer, intent(in) :: k, n
    real(8), dimension (:,:,:), intent(in) :: A
    real(8), target :: B(size(A,1),size(A,2),size(A,3))
    real(8), pointer :: C(:,:,:)

    integer :: nx, ny, nz, i, j, kk

    !system_clock
    !REAL(8) :: system_clock_rate,t1,t2
    !INTEGER :: c1,c2,cr,cm

    ! First initialize the system_clock
    !CALL system_clock(count_rate=cr)
    !CALL system_clock(count_max=cm)
    !system_clock_rate = REAL(cr)

    !if (size(A,1)==1) A=transpose(A)

    nx=size(A,1)
    ny=size(A,2)
    nz=size(A,3)

    if (k==1) then
        !CALL SYSTEM_CLOCK(c1)
        !$omp parallel do
        do kk=1, nz
            do j=1, ny
                do i=1, nx-1
                    B(i,j,kk)=A(i+1,j,kk)-A(i,j,kk)
                end do
            end do
        end do
        !$omp end parallel do
        !CALL SYSTEM_CLOCK(c2)
        !print '("    omp wall time 2: ", F8.4, " second")', (c2-c1)/system_clock_rate

        C => B(1:nx-1,:,:)
    else if (k==2) then
        !$omp parallel do
        do kk=1, nz
            do j=1, ny-1
                do i=1, nx
                    B(i,j,kk)=A(i,j+1,kk)-A(i,j,kk)
                end do
            end do
        end do
        !$omp end parallel do

        C => B(:,1:ny-1,:)
    else if (k==3) then
        !$omp parallel do
        do kk=1, nz-1
            do j=1, ny
                do i=1, nx
                    B(i,j,kk)=A(i,j,kk+1)-A(i,j,kk)
                end do
            end do
        end do
        !$omp end parallel do

        C => B(:,:,1:nz-1)
    end if

    end function diff

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function diff2(A, k, bc) result(C)
    implicit none
    integer, intent(in) :: k, bc
    real(8), dimension (:,:,:), intent(in) :: A
    real(8), target :: B(size(A,1),size(A,2),size(A,3))
    real(8), pointer :: C(:,:,:)

    integer :: nx, ny, nz, i, j, kk

    !system_clock
    REAL(8) :: system_clock_rate,t1,t2
    INTEGER :: c1,c2,cr,cm

    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    system_clock_rate = REAL(cr)

    !if (size(A,1)==1) A=transpose(A)

    nx=size(A,1)
    ny=size(A,2)
    nz=size(A,3)

    if (bc==1) then
        if (k==1) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny
                    do i=1, nx-1
                        if (i==1) then
                            B(1,j,kk)=A(nx-1,j,kk)-2.0d0*A(1,j,kk)+A(2,j,kk)
                        else
                            B(i,j,kk)=A(i-1,j,kk)-2.0d0*A(i,j,kk)+A(i+1,j,kk)
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(1:nx-1,:,:)
        else if (k==2) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny-1
                    do i=1, nx
                        if (j==1) then
                            B(i,1,kk)=A(i,ny-1,kk)-2.0d0*A(i,1,kk)+A(i,2,kk)
                        else
                            B(i,j,kk)=A(i,j-1,kk)-2.0d0*A(i,j,kk)+A(i,j+1,kk)
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,1:ny-1,:)
        else if (k==3) then
            !$omp parallel do
            do kk=1, nz-1
                do j=1, ny
                    do i=1, nx
                        if (kk==1) then
                            B(i,j,1)=A(i,j,nz-1)-2.0d0*A(i,j,1)+A(i,j,2)
                        else
                            B(i,j,kk)=A(i,j,kk-1)-2.0d0*A(i,j,kk)+A(i,j,kk+1)
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,:,1:nz-1)
        end if
    else
        if (k==1) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny
                    do i=1, nx-2
                        B(i,j,kk)=A(i,j,kk)-2.0d0*A(i+1,j,kk)+A(i+2,j,kk)
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(1:nx-2,:,:)
        else if (k==2) then
            !$omp parallel do
            do kk=1, nz
                do j=1, ny-2
                    do i=1, nx
                        B(i,j,kk)=A(i,j,kk)-2.0d0*A(i,j+1,kk)+A(i,j+2,kk)
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,1:ny-2,:)
        else if (k==3) then
            !$omp parallel do
            do kk=1, nz-2
                do j=1, ny
                    do i=1, nx
                        B(i,j,kk)=A(i,j,kk)-2.0d0*A(i,j,kk+1)+A(i,j,kk+2)
                    end do
                end do
            end do
            !$omp end parallel do

            C => B(:,:,1:nz-2)
        end if
    end if

    end function diff2

    end module FD_functions

