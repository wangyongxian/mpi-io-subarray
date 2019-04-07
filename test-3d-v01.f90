! 利用MPI IO对分布式3维数组进行并行IO
! Wang Yong-Xian
! 2019-04-07
!
! compile:  mpif90 test-3d-v01.f90 -o test-3d-v01.exe
! run:      yhrun -n 16 -N 2 -p TH_NET1 ./test-3d-v01.exe
!
program test
    use mpi
    implicit none

    integer, parameter :: m = 1024, n = 512, d = 16
    integer, parameter :: npi = 2, npj = 4, npk = 2, np = npi * npj * npk
    
    real(kind=4), allocatable :: dat(:, :, :), cpy(:, :, :)
    integer            :: local_array_size
    integer            :: nprocs, rank
    integer            :: i, j, k, count, ierr
    character(len=80)  :: filename
    integer            :: filetype, comm

    integer, parameter :: ndims = 3
    integer, dimension(1:ndims) ::  gsizes, lsizes, psizes, &
                                    coords, start_indices

    call initialize()

    call MPI_TYPE_CREATE_SUBARRAY(ndims, &
            gsizes, lsizes, start_indices, &
            MPI_ORDER_FORTRAN, MPI_FLOAT, filetype, ierr)
    call MPI_TYPE_COMMIT(filetype, ierr)

    call write_file_demo()

    call read_file_demo()

    do i = 0, lsizes(1)-1
    do j = 0, lsizes(2)-1
    do k = 0, lsizes(3)-1
        if ( dat(i,j,k) /= cpy(i,j,k) ) then
            write(*, "(4I8,2F12.5)") rank, i, j, k, dat(i,j,k), cpy(i,j,k)
        end if
    end do
    end do
    end do

    call MPI_TYPE_FREE(filetype, ierr)

    call finalize()
contains

subroutine write_file_demo()
    implicit none

    integer :: status(MPI_STATUS_SIZE)
    integer :: fh    

    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,        &
            MPI_MODE_WRONLY + MPI_MODE_CREATE,          &
            MPI_INFO_NULL,                              &
            fh, ierr )
    call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND,       &
            MPI_FLOAT, filetype, "native",              &
            MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE_ALL(fh, dat, local_array_size,  &
            MPI_FLOAT, status, ierr)
    call MPI_GET_COUNT( status, MPI_INTEGER, count, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

    write(*, *) 'process ', rank, 'write ', count, 'integers'

end subroutine

subroutine read_file_demo()
    implicit none

    integer :: status(MPI_STATUS_SIZE)
    integer :: fh

    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,        &
            MPI_MODE_RDONLY,                            &
            MPI_INFO_NULL,                              &
            fh, ierr )
    call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND,       &
            MPI_FLOAT, filetype, "native",              &
            MPI_INFO_NULL, ierr)
    call MPI_FILE_READ_ALL(fh, cpy, local_array_size,   &
            MPI_FLOAT, status, ierr)
    call MPI_GET_COUNT( status, MPI_INTEGER, count, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

    write(*, *) 'process ', rank, 'read ', count, 'integers'

end subroutine

subroutine initialize() 
    implicit none

    LOGICAL, PARAMETER :: reorder = .FALSE.
    LOGICAL, PARAMETER :: periodic(1:ndims) = .FALSE.

    filename = "dat_2d_v01.dat"

    gsizes(1) = m   ! # rows in global array
    gsizes(2) = n   ! # cols in global array
    gsizes(3) = d

    psizes(1) = npi ! # procs in i-dimension
    psizes(2) = npj ! # procs in j-dimension
    psizes(3) = npk

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank,   ierr)
    if (nprocs < np .or. mod(nprocs, npi) /= 0) stop " nprocs < np or nprocs % npi /= 0"

    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, psizes, &
            periodic, reorder, comm, ierr)
    call MPI_CART_COORDS(comm, rank, ndims, coords, ierr)

    lsizes(:) = gsizes(:) / psizes(:) ! dim of local array
    local_array_size = lsizes(1) * lsizes(2) * lsizes(3)

    ! global indices of first element of local array
    start_indices(:) = coords(:) * lsizes(:)

    allocate( dat(0:lsizes(1)-1, 0:lsizes(2)-1, 0:lsizes(3)-1) )
    allocate( cpy(0:lsizes(1)-1, 0:lsizes(2)-1, 0:lsizes(3)-1) )

    do i = 0, lsizes(1)-1
    do j = 0, lsizes(2)-1
    do k = 0, lsizes(3)-1
        dat(i, j, k) = i * lsizes(2) + j * lsizes(3) + 1 + rank * 1000
        cpy(i, j, k) = 0
    end do
    end do
    end do
end subroutine

subroutine finalize()
    implicit none

    deallocate( dat )
    deallocate( cpy )

    call MPI_FINALIZE(ierr)
end subroutine

end program