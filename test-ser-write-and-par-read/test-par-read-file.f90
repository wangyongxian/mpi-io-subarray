!======================================================================
     program main
!======================================================================
     use mpi
     implicit none

     integer, parameter :: ndims = 4
     integer, parameter, dimension(1:ndims) :: gsizes = (/ 4, 3, 4, 2/)
     integer, parameter, dimension(1:ndims) :: psizes = (/ 2, 3, 2, 1/)
     integer, parameter, dimension(1:ndims) :: lsizes = (/ 2, 1, 2, 2/)

     integer, parameter :: il = lsizes(1) * lsizes(2)
     integer, parameter :: jl = lsizes(3)
     integer, parameter :: kl = lsizes(4)

     real*4, allocatable :: subarray(:, :, :)
     real*4, allocatable :: localdat(:, :, :)

     character*3 :: ithreadstr
     character*90 :: filename
!===================================================================================
    integer, dimension(1:ndims) :: coords, start_indices
    integer            :: local_array_size
    integer :: nprocs, rank
    integer :: ierr
    integer :: filetype, fh
    integer :: status(MPI_STATUS_SIZE)
    integer :: count, val
    LOGICAL, PARAMETER :: reorder = .FALSE.
    LOGICAL, PARAMETER :: periodic(1:ndims) = .FALSE.
    integer :: comm
    integer :: i, j, k, d, i0, j0, k0, d0, ith

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank,   ierr)
    if (nprocs < product(psizes)) stop " nprocs < np !!"


    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, psizes, &
    &    periodic, reorder, comm, ierr)
    call MPI_CART_COORDS(comm, rank, ndims, coords, ierr)

    start_indices(:) = coords(:) * lsizes(:)
    local_array_size = product(lsizes)

    call MPI_TYPE_CREATE_SUBARRAY(ndims,                &
    &       gsizes, lsizes, start_indices,              &
    &       MPI_ORDER_FORTRAN, MPI_FLOAT, filetype, ierr)
    call MPI_TYPE_COMMIT(filetype, ierr)


    allocate( subarray(il, jl, kl) )
    allocate( localdat(il, jl, kl) )

!------- 'ith' can be derived from coordinates of current MPI rank
        i = coords(1); j = coords(2)
        k = coords(3); d = coords(4)
        ith = 1 + i + psizes(1) * (k + psizes(3) * (j + psizes(2) * d))
        write(ithreadstr,'(i3.3)')  ith
        filename = "test."//ithreadstr//".dat"

        open(21, file=filename, form='binary',          &
                 action='read',   &
                 status='old')
            read(21) localdat
        close(21)

!        write(*, "(A, 5I8, A, I8, A, 4I8)") "--> rank = ", rank, coords(:), &
!            ", ith = ", ith, ', start_indices ', start_indices

!------- read file using MPI IO operation
    
        filename = "test.global.dat"
        call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,        &
                MPI_MODE_RDONLY,                            &
                MPI_INFO_NULL,                              &
                fh, ierr )
        call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND,       &
                MPI_FLOAT, filetype, "native",            &
                MPI_INFO_NULL, ierr)
        call MPI_FILE_READ_ALL(fh, subarray, local_array_size,   &
                MPI_FLOAT, status, ierr)
        call MPI_GET_COUNT( status, MPI_FLOAT, count, ierr)
        call MPI_FILE_CLOSE(fh, ierr)
        !write(*, *) 'process ', rank, 'read ', count, 'integers'

        if (rank == 0) then
            do k0 = 1, kl
            do j0 = 1, jl
            do i0 = 1, il
                if (subarray(i0, j0, k0) /= localdat(i0, j0, k0)) then
                    write(*, "(A, I8, A, 3I8, 2(x, Z8.8))")  &
                                "rank = ", ith,         &
                                "pos = ", i0, j0, k0,   &
                                subarray(i0, j0, k0),   &
                                localdat(i0, j0, k0)
                endif
            enddo
            enddo
            enddo
        endif

    deallocate(subarray)
    deallocate(localdat)
        
    call MPI_TYPE_FREE(filetype, ierr)
    call MPI_FINALIZE(ierr)

!======================================================================
    end program
!======================================================================
