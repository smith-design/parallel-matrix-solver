program read_matrix_mpi
    use mpi
    implicit none

    ! === 全局常量定义 ===
    integer, parameter :: wp = 8 ! Double Precision

    ! 变量声明
    integer :: ierr, rank, nprocs
    integer :: rows_global, cols_global
    integer :: rows_local
    integer :: i, j, k, file_unit
    integer :: start_row, end_row
    CHARACTER(LEN=60) :: filename_m, filename_a
    
    ! 用于计算分配的临时变量
    integer :: n_base, n_rem

    ! 假设矩阵维度 (需与 MATLAB 生成的一致)
    parameter (rows_global = 3005)
    parameter (cols_global = 3005)
    
    ! 动态分配数组
    complex(kind=wp), allocatable :: local_M(:,:)
    complex(kind=wp), allocatable :: local_A(:) ! 新增：本地向量 A
    
    ! 临时读取缓冲区
    real(kind=wp), allocatable :: temp_buffer(:) ! 用于矩阵的一行
    real(kind=wp) :: temp_vec_buf(2)             ! Used for vector (Real, Imag)

    ! MPI 初始化
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    ! === 1. 计算负载均衡 (Row Decomposition) ===
    n_base = rows_global / nprocs
    n_rem  = mod(rows_global, nprocs)

    if (rank < n_rem) then
        rows_local = n_base + 1
        start_row = rank * (n_base + 1) + 1
    else
        rows_local = n_base
        start_row = n_rem * (n_base + 1) + (rank - n_rem) * n_base + 1
    end if
    end_row = start_row + rows_local - 1
    
    ! 安全检查
    if (rows_local < 0) rows_local = 0
    
    ! === 2. 内存分配 ===
    if (rows_local > 0) then
        allocate(local_M(rows_local, cols_global))
        allocate(local_A(rows_local)) ! 分配向量内存
    else
        allocate(local_M(0, cols_global)) 
        allocate(local_A(0))
    end if
    
    allocate(temp_buffer(2 * cols_global))

    ! === 3. 并行读取矩阵文件 (complex_matrix.txt) ===
    file_unit = 10
    open(unit=file_unit, file='complex_matrix.txt', status='old', action='read')

    k = 0 
    do i = 1, rows_global
        ! 读取矩阵的一行
        read(file_unit, *) temp_buffer
        
        if (i >= start_row .and. i <= end_row) then
            k = k + 1
            do j = 1, cols_global
                local_M(k, j) = cmplx(temp_buffer(2*j-1), temp_buffer(2*j), kind=wp)
            end do
        end if
    end do
    close(file_unit)
    
    ! === 4. 并行读取向量文件 (rhs_vector.txt) ===
    ! 所有进程都需要扫描文件，但只存储属于自己的行
    file_unit = 11
    open(unit=file_unit, file='rhs_vector.txt', status='old', action='read')
    
    k = 0
    do i = 1, rows_global
        ! 读取向量的一行 (实部 虚部)
        read(file_unit, *) temp_vec_buf
        
        if (i >= start_row .and. i <= end_row) then
            k = k + 1
            local_A(k) = cmplx(temp_vec_buf(1), temp_vec_buf(2), kind=wp)
        end if
        ! 如果不是我的行，循环继续，文件指针自动下移
    end do
    close(file_unit)

    ! === 5. 写入本地结果进行验证 ===
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    ! 构造文件名
    WRITE(filename_m, '("M_Rank",I3.3,".txt")') rank
    WRITE(filename_a, '("A_Rank",I3.3,".txt")') rank
    
    if (rows_local > 0) then
        WRITE(*,*) 'Rank ', rank, ' writing M to ', TRIM(filename_m), ' and A to ', TRIM(filename_a)
        
        ! 写入矩阵 M
        CALL write_complex_matrix_to_file(filename_m, local_M)
        
        ! 写入向量 A
        CALL write_complex_vector_to_file(filename_a, local_A)
    else
        WRITE(*,*) 'Rank ', rank, ' is empty, skipping write.'
    end if

    ! 清理内存
    if (allocated(local_M)) deallocate(local_M)
    if (allocated(local_A)) deallocate(local_A)
    if (allocated(temp_buffer)) deallocate(temp_buffer)

    call MPI_FINALIZE(ierr)

CONTAINS

    ! 子程序：写入复数矩阵
    SUBROUTINE write_complex_matrix_to_file(fname, matrix_data)
        CHARACTER(len=*), INTENT(IN) :: fname
        COMPLEX(kind=wp), INTENT(IN) :: matrix_data(:, :)

        INTEGER :: f_unit, ii, jj, nz, mm
        CHARACTER(len=50) :: fmt_string 
        
        nz = SIZE(matrix_data, 1)
        mm = SIZE(matrix_data, 2)

        OPEN(NEWUNIT=f_unit, FILE=fname, STATUS='REPLACE', ACTION='WRITE')

        ! 动态构建格式字符串: ( n个 (2F30.16, 2X) )
        WRITE(fmt_string, '("(",I0,"(2F30.16, 2X))")') mm

        DO ii = 1, nz
            WRITE(f_unit, fmt_string) &
                 ( REAL(matrix_data(ii, jj)), AIMAG(matrix_data(ii, jj)), jj = 1, mm )
        END DO

        CLOSE(f_unit)
    END SUBROUTINE write_complex_matrix_to_file

    ! 子程序：写入复数向量
    SUBROUTINE write_complex_vector_to_file(fname, vector_data)
        CHARACTER(len=*), INTENT(IN) :: fname
        COMPLEX(kind=wp), INTENT(IN) :: vector_data(:)

        INTEGER :: f_unit, ii, nz

        nz = SIZE(vector_data)

        OPEN(NEWUNIT=f_unit, FILE=fname, STATUS='REPLACE', ACTION='WRITE')
        
        DO ii = 1, nz
            ! 每个复数写一行: 实部 虚部
            WRITE(f_unit, '(2F30.16)') REAL(vector_data(ii)), AIMAG(vector_data(ii))
        END DO
        
        CLOSE(f_unit)
    END SUBROUTINE write_complex_vector_to_file

end program read_matrix_mpi
