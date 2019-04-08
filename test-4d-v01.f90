! 利用MPI IO对分布式多维数组进行并行IO
! 值得注意的是，本例中，本质上是4维的数据，实际上只表示成3维数组，并按3维进行并行划分
!
! Wang Yong-Xian
! 2019-04-07
!
! compile:  mpif90 test-4d-v01.f90 -o test-4d-v01.exe
! run:      yhrun -n 12 -N 3 -p TH_NET1 ./test-4d-v01.exe
!
program test
    use mpi
    implicit none

    integer, parameter :: ndims = 4
!    integer, parameter, dimension(1:ndims) :: gsizes = (/1020, 6, 1020, 50/)
!    integer, parameter, dimension(1:ndims) :: psizes = (/   6, 6,   12,  1/)
!    integer, parameter, dimension(1:ndims) :: lsizes = (/ 170, 1,   85, 50/)
    integer, parameter, dimension(1:ndims) :: gsizes = (/ 40, 33, 30,  10/)
    integer, parameter, dimension(1:ndims) :: psizes = (/  2,  3,  2,   1/)
    integer, parameter, dimension(1:ndims) :: lsizes = (/ 20, 11, 15,  10/)
    integer, dimension(1:ndims) :: coords, start_indices
    
    real(kind=4), allocatable :: dat(:, :, :), cpy(:, :, :)
    integer            :: local_array_size
    integer            :: nprocs, rank
    integer            :: i, j, k, count, ierr
    character(len=80)  :: filename
    integer            :: filetype


    call initialize()

    ! write(*,"(2I8,3(',', 4I8))") rank, ndims, gsizes, lsizes, start_indices
    call MPI_TYPE_CREATE_SUBARRAY(ndims,                &
            gsizes, lsizes, start_indices,              &
            MPI_ORDER_FORTRAN, MPI_FLOAT, filetype, ierr)
    call MPI_TYPE_COMMIT(filetype, ierr)

    call write_file_demo()

    call read_file_demo()

    do i = 0, lsizes(1)*lsizes(2)-1
    do j = 0, lsizes(3)-1
    do k = 0, lsizes(4)-1
            if ( dat(i,j,k) /= cpy(i,j,k) ) then
            write(*, "(4I8,2F14.0)") rank, i, j, k, dat(i,j,k), cpy(i,j,k)
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
    integer :: comm

    filename = "dat_4d_v01.dat"

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank,   ierr)
    if (nprocs < product(psizes)) stop " nprocs < np !!"
    
    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, psizes, &
            periodic, reorder, comm, ierr)
    call MPI_CART_COORDS(comm, rank, ndims, coords, ierr)

    start_indices(:) = coords(:) * lsizes(:)
    local_array_size = product(lsizes)

    allocate( dat(0:lsizes(1)*lsizes(2)-1, 0:lsizes(3)-1, 0:lsizes(4)-1) )
    allocate( cpy(0:lsizes(1)*lsizes(2)-1, 0:lsizes(3)-1, 0:lsizes(4)-1) )

    do i = 0, lsizes(1)*lsizes(2)-1
    do j = 0, lsizes(3)-1
    do k = 0, lsizes(4)-1
        dat(i, j, k) = k * 1 + j * 100 + i * 10000 + rank * 1000000
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