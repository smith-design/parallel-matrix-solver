program lu_solver
    implicit none
    integer, parameter :: wp = 8
    integer, parameter :: N = 3005
    
    complex(kind=wp), allocatable :: M(:,:), r(:), z(:)
    real(kind=wp), allocatable :: buf(:)
    real(kind=wp) :: vbuf(2)
    integer :: i, j, k, funit
    complex(kind=wp) :: s, residual, pivot, factor
    integer, allocatable :: ipiv(:)
    integer :: max_row
    complex(kind=wp) :: max_val, temp
    
    allocate(M(N,N), r(N), z(N), buf(2*N), ipiv(N))
    
    write(*,*) 'Reading matrix...'
    funit = 10
    open(unit=funit, file='datasets/complex_matrix.txt', status='old', action='read')
    do i = 1, N
        read(funit, *) buf
        do j = 1, N
            M(i,j) = cmplx(buf(2*j-1), buf(2*j), kind=wp)
        end do
    end do
    close(funit)
    
    write(*,*) 'Reading RHS...'
    funit = 11
    open(unit=funit, file='datasets/rhs_vector.txt', status='old', action='read')
    do i = 1, N
        read(funit, *) vbuf
        r(i) = cmplx(vbuf(1), vbuf(2), kind=wp)
    end do
    close(funit)
    deallocate(buf)
    
    write(*,*) 'LU factorization with partial pivoting...'
    do i = 1, N
        ipiv(i) = i
    end do
    
    do k = 1, N-1
        if (mod(k, 500) == 0) write(*,*) 'Column', k, 'of', N
        
        ! Find pivot
        max_val = M(k,k)
        max_row = k
        do i = k+1, N
            if (abs(M(i,k)) > abs(max_val)) then
                max_val = M(i,k)
                max_row = i
            end if
        end do
        
        ! Swap rows
        if (max_row /= k) then
            do j = 1, N
                temp = M(k,j)
                M(k,j) = M(max_row,j)
                M(max_row,j) = temp
            end do
            i = ipiv(k)
            ipiv(k) = ipiv(max_row)
            ipiv(max_row) = i
            
            temp = r(k)
            r(k) = r(max_row)
            r(max_row) = temp
        end if
        
        ! Elimination
        if (abs(M(k,k)) > 1.0e-14_wp) then
            do i = k+1, N
                factor = M(i,k) / M(k,k)
                M(i,k) = factor  ! Store L
                do j = k+1, N
                    M(i,j) = M(i,j) - factor * M(k,j)
                end do
                r(i) = r(i) - factor * r(k)
            end do
        end if
    end do
    
    write(*,*) 'Back substitution...'
    z = (0.0_wp, 0.0_wp)
    do i = N, 1, -1
        s = r(i)
        do j = i+1, N
            s = s - M(i,j) * z(j)
        end do
        if (abs(M(i,i)) > 1.0e-14_wp) then
            z(i) = s / M(i,i)
        end if
    end do
    
    write(*,*) 'Solution computed.'
    write(*,*) 'First 5 elements of z:'
    do i = 1, 5
        write(*,*) i, z(i)
    end do
    
    ! Write solution
    funit = 20
    open(unit=funit, file='solution_lu.txt', status='replace', action='write')
    do i = 1, N
        write(funit, '(2E25.16)') real(z(i)), aimag(z(i))
    end do
    close(funit)
    write(*,*) 'Solution written to solution_lu.txt'
    
    deallocate(M, r, z, ipiv)
end program lu_solver
