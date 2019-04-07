! 利用MPI IO对分布式1维数组进行并行IO
! Wang Yong-Xian
! 2019-04-07
!
! compile:  mpif90 test-1d-v01.f90 -o test-1d-v01.exe
! run:      yhrun -n 8 -N 1 -p TH_NET1 ./test-1d-v01.exe
!
program test
    use mpi
    implicit none

    integer, parameter :: filesize = 1024 * 1024 * 16
    integer, parameter :: int_size = 4
    integer, parameter :: n = filesize / int_size
    integer(kind=4)    :: dat(0:n-1), cpy(0:n-1)
    integer :: nprocs, rank
    integer :: i, ierr, num_ints, count, offset_ints
    character(len=80)  :: filename
    integer(kind=MPI_OFFSET_KIND) :: offset_bytes

    filename = "dat_1d_v01.dat"
    do i = 0, n-1
        dat(i) = i
        cpy(i) = 0
    end do

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank,   ierr)

    num_ints     = n / nprocs
    offset_ints  = rank * num_ints
    offset_bytes = offset_ints * int_size

    call write_file_demo()

    call read_file_demo()

    do i = offset_ints, offset_ints + num_ints - 1
        if ( dat(i) /= cpy(i) ) then
            write(*, *) "rank / index / dat / cpy: ", rank, i, dat(i), cpy(i)
        end if
    end do

    call MPI_FINALIZE(ierr)

contains
subroutine write_file_demo()
    implicit none

    integer :: status(MPI_STATUS_SIZE)
    integer :: fh

    call MPI_FILE_OPEN( MPI_COMM_WORLD, filename,       &
            MPI_MODE_WRONLY + MPI_MODE_CREATE,          &
            MPI_INFO_NULL,                              &
            fh, ierr )
    call MPI_FILE_WRITE_AT( fh, offset_bytes,           &
            dat(offset_ints:offset_ints+num_ints-1),    &
            num_ints,                                   &
            MPI_INTEGER, status, ierr )
    call MPI_GET_COUNT( status, MPI_INTEGER, count, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

    write(*, *) 'process ', rank, 'write ', count, 'integers'

end subroutine

subroutine read_file_demo()
    implicit none

    integer :: status(MPI_STATUS_SIZE)
    integer :: fh

    call MPI_FILE_OPEN( MPI_COMM_WORLD, filename,       &
            MPI_MODE_RDONLY,                            &
            MPI_INFO_NULL,                              &
            fh, ierr )
    call MPI_FILE_READ_AT( fh, offset_bytes,            &
            cpy(offset_ints:offset_ints+num_ints-1),    &
            num_ints,                                   &
            MPI_INTEGER, status, ierr )
    call MPI_GET_COUNT( status, MPI_INTEGER, count, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

    write(*, *) 'process ', rank, 'read ', count, 'integers'

end subroutine
end program