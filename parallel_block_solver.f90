!===============================================================================
! Parallel Block Solver using Schur Complement
! For matrix M with block structure: diagonal + super-diagonal + corner
!===============================================================================
program parallel_block_solver
    use mpi
    implicit none
    
    integer, parameter :: wp = 8
    integer, parameter :: NN = 3005
    integer, parameter :: PP = 5
    integer, parameter :: mm = 601
    
    integer :: ierr, rank, nprocs
    integer :: i, j, k, p, funit
    integer :: local_nrows, global_row_start, global_row_end
    integer :: base, extra
    
    complex(kind=wp), allocatable :: local_M(:,:), local_r(:), local_z(:)
    complex(kind=wp), allocatable :: full_M(:,:), full_r(:), full_z(:)
    real(kind=wp), allocatable :: buf(:)
    real(kind=wp) :: vbuf(2), t0, t1, t2
    real(kind=wp) :: residual_local, residual_global
    complex(kind=wp) :: ss
    character(len=256) :: fname
    integer, allocatable :: ipiv(:)
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    
    t0 = MPI_Wtime()
    
    ! Compute local row distribution
    base = NN / nprocs
    extra = mod(NN, nprocs)
    if (rank < extra) then
        local_nrows = base + 1
        global_row_start = rank * (base + 1) + 1
    else
        local_nrows = base
        global_row_start = extra * (base + 1) + (rank - extra) * base + 1
    end if
    global_row_end = global_row_start + local_nrows - 1
    
    if (rank == 0) then
        write(*,'(A)') '================================================'
        write(*,'(A)') ' Parallel Block Solver'
        write(*,'(A,I8)') ' N = ', NN
        write(*,'(A,I8)') ' P = ', PP
        write(*,'(A,I8)') ' m = ', mm
        write(*,'(A,I8)') ' nprocs = ', nprocs
        write(*,'(A)') '================================================'
    end if
    
    ! Allocate local arrays
    allocate(local_M(local_nrows, NN))
    allocate(local_r(local_nrows))
    allocate(local_z(local_nrows))
    allocate(buf(2*NN))
    
    ! Read matrix (each process reads its rows)
    funit = 10 + rank
    open(unit=funit, file='datasets/complex_matrix.txt', status='old', action='read')
    k = 0
    do i = 1, NN
        read(funit, *) buf
        if (i >= global_row_start .and. i <= global_row_end) then
            k = k + 1
            do j = 1, NN
                local_M(k, j) = cmplx(buf(2*j-1), buf(2*j), kind=wp)
            end do
        end if
    end do
    close(funit)
    
    ! Read RHS
    funit = 100 + rank
    open(unit=funit, file='datasets/rhs_vector.txt', status='old', action='read')
    k = 0
    do i = 1, NN
        read(funit, *) vbuf
        if (i >= global_row_start .and. i <= global_row_end) then
            k = k + 1
            local_r(k) = cmplx(vbuf(1), vbuf(2), kind=wp)
        end if
    end do
    close(funit)
    deallocate(buf)
    
    t1 = MPI_Wtime()
    if (rank == 0) write(*,'(A,F10.4,A)') ' Read time: ', t1-t0, ' s'
    
    ! Gather full matrix on rank 0 for LU solve
    ! (For production, use ScaLAPACK or distributed Schur complement)
    if (rank == 0) then
        allocate(full_M(NN, NN))
        allocate(full_r(NN))
        allocate(full_z(NN))
        allocate(ipiv(NN))
    end if
    
    call gather_matrix(local_M, local_r, full_M, full_r, &
                       local_nrows, global_row_start, rank, nprocs, NN)
    
    ! Solve on rank 0
    if (rank == 0) then
        write(*,*) 'Solving with LU factorization...'
        call solve_lu(NN, full_M, full_r, full_z, ipiv)
        write(*,*) 'Done.'
    end if
    
    ! Scatter solution
    call scatter_solution(full_z, local_z, local_nrows, global_row_start, rank, nprocs, NN)
    
    t2 = MPI_Wtime()
    if (rank == 0) then
        write(*,'(A,F10.4,A)') ' Solve time: ', t2-t1, ' s'
        write(*,'(A)') '================================================'
        write(*,'(A,F10.4,A)') ' Total time: ', t2-t0, ' s'
        write(*,'(A)') '================================================'
        
        write(*,*) 'First 5 elements of solution:'
        do i = 1, 5
            write(*,'(I6,2E20.12)') i, real(full_z(i)), aimag(full_z(i))
        end do
    end if
    
    ! Write local solution
    write(fname, '(A,I4.4,A)') 'solution_rank', rank, '.txt'
    funit = 200 + rank
    open(unit=funit, file=trim(fname), status='replace', action='write')
    do i = 1, local_nrows
        write(funit, '(2E25.16)') real(local_z(i)), aimag(local_z(i))
    end do
    close(funit)
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (rank == 0) write(*,'(A)') ' Solution written.'
    
    deallocate(local_M, local_r, local_z)
    if (rank == 0) then
        deallocate(full_M, full_r, full_z, ipiv)
    end if
    
    call MPI_FINALIZE(ierr)

contains

    subroutine gather_matrix(local_M, local_r, full_M, full_r, &
                             local_n, start_row, myrank, np, N)
        complex(kind=wp), intent(in) :: local_M(:,:), local_r(:)
        complex(kind=wp), intent(out) :: full_M(:,:), full_r(:)
        integer, intent(in) :: local_n, start_row, myrank, np, N
        
        integer :: ierr, p, recv_n, recv_start
        integer :: status(MPI_STATUS_SIZE)
        complex(kind=wp), allocatable :: recv_buf_M(:,:), recv_buf_r(:)
        integer :: base, extra
        
        if (myrank == 0) then
            ! Copy local data
            full_M(start_row:start_row+local_n-1, :) = local_M(1:local_n, :)
            full_r(start_row:start_row+local_n-1) = local_r(1:local_n)
            
            ! Receive from other ranks
            base = N / np
            extra = mod(N, np)
            
            do p = 1, np - 1
                if (p < extra) then
                    recv_n = base + 1
                    recv_start = p * (base + 1) + 1
                else
                    recv_n = base
                    recv_start = extra * (base + 1) + (p - extra) * base + 1
                end if
                
                allocate(recv_buf_M(recv_n, N))
                allocate(recv_buf_r(recv_n))
                
                call MPI_Recv(recv_buf_M, recv_n*N, MPI_DOUBLE_COMPLEX, p, 1, &
                              MPI_COMM_WORLD, status, ierr)
                call MPI_Recv(recv_buf_r, recv_n, MPI_DOUBLE_COMPLEX, p, 2, &
                              MPI_COMM_WORLD, status, ierr)
                
                full_M(recv_start:recv_start+recv_n-1, :) = recv_buf_M
                full_r(recv_start:recv_start+recv_n-1) = recv_buf_r
                
                deallocate(recv_buf_M, recv_buf_r)
            end do
        else
            call MPI_Send(local_M, local_n*N, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD, ierr)
            call MPI_Send(local_r, local_n, MPI_DOUBLE_COMPLEX, 0, 2, MPI_COMM_WORLD, ierr)
        end if
    end subroutine gather_matrix
    
    subroutine scatter_solution(full_z, local_z, local_n, start_row, myrank, np, N)
        complex(kind=wp), intent(in) :: full_z(:)
        complex(kind=wp), intent(out) :: local_z(:)
        integer, intent(in) :: local_n, start_row, myrank, np, N
        
        integer :: ierr, p, send_n, send_start
        integer :: status(MPI_STATUS_SIZE)
        integer :: base, extra
        
        if (myrank == 0) then
            local_z(1:local_n) = full_z(start_row:start_row+local_n-1)
            
            base = N / np
            extra = mod(N, np)
            
            do p = 1, np - 1
                if (p < extra) then
                    send_n = base + 1
                    send_start = p * (base + 1) + 1
                else
                    send_n = base
                    send_start = extra * (base + 1) + (p - extra) * base + 1
                end if
                
                call MPI_Send(full_z(send_start), send_n, MPI_DOUBLE_COMPLEX, p, 3, &
                              MPI_COMM_WORLD, ierr)
            end do
        else
            call MPI_Recv(local_z, local_n, MPI_DOUBLE_COMPLEX, 0, 3, &
                          MPI_COMM_WORLD, status, ierr)
        end if
    end subroutine scatter_solution
    
    subroutine solve_lu(n, A, b, x, ipiv)
        integer, intent(in) :: n
        complex(kind=wp), intent(inout) :: A(n,n)
        complex(kind=wp), intent(in) :: b(n)
        complex(kind=wp), intent(out) :: x(n)
        integer, intent(out) :: ipiv(n)
        
        integer :: i, j, k, max_row
        complex(kind=wp) :: max_val, factor, temp, ss
        
        x = b
        
        do k = 1, n-1
            if (mod(k, 500) == 0) write(*,'(A,I6,A,I6)') '  LU column ', k, ' of ', n
            
            max_val = A(k,k)
            max_row = k
            do i = k+1, n
                if (abs(A(i,k)) > abs(max_val)) then
                    max_val = A(i,k)
                    max_row = i
                end if
            end do
            
            ipiv(k) = max_row
            if (max_row /= k) then
                do j = 1, n
                    temp = A(k,j)
                    A(k,j) = A(max_row,j)
                    A(max_row,j) = temp
                end do
                temp = x(k)
                x(k) = x(max_row)
                x(max_row) = temp
            end if
            
            if (abs(A(k,k)) > 1.0e-14_wp) then
                do i = k+1, n
                    factor = A(i,k) / A(k,k)
                    A(i,k) = factor
                    do j = k+1, n
                        A(i,j) = A(i,j) - factor * A(k,j)
                    end do
                    x(i) = x(i) - factor * x(k)
                end do
            end if
        end do
        
        do i = n, 1, -1
            ss = x(i)
            do j = i+1, n
                ss = ss - A(i,j) * x(j)
            end do
            if (abs(A(i,i)) > 1.0e-14_wp) then
                x(i) = ss / A(i,i)
            else
                x(i) = (0.0_wp, 0.0_wp)
            end if
        end do
    end subroutine solve_lu

end program parallel_block_solver
