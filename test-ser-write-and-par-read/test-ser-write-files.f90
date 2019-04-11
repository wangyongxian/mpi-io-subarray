!======================================================================
     program main
!======================================================================
     implicit none

     integer, parameter :: ndims = 4
     integer, parameter, dimension(1:ndims) :: gsizes = (/ 4, 3, 4, 2/)
     integer, parameter, dimension(1:ndims) :: psizes = (/ 2, 3, 2, 1/)
     integer, parameter, dimension(1:ndims) :: lsizes = (/ 2, 1, 2, 2/)

     integer, parameter :: il = lsizes(1) * lsizes(2)
     integer, parameter :: jl = lsizes(3)
     integer, parameter :: kl = lsizes(4)
     integer, parameter :: im = gsizes(1) * gsizes(2)
     integer, parameter :: jm = gsizes(3)
     integer, parameter :: km = gsizes(4)

     real*4, allocatable :: local(:, :, :)
     real*4, allocatable :: global(:, :, :)

     character*3 :: ithreadstr
     character*90 :: filename
!===================================================================================
     integer, dimension(1:ndims) :: coords, start_indices
     integer :: count, val, ith
     integer :: i, j, k, d
     integer :: i0, j0, k0, d0
     integer :: ig, jg, kg, dg

    allocate( local(il, jl, kl) )
    allocate( global(im, jm, km) )
    local = 0
    global = 0

    count = 0
    do i = 0, psizes(1)-1
    do j = 0, psizes(2)-1
    do k = 0, psizes(3)-1
    do d = 0, psizes(4)-1
        start_indices = lsizes(:) * [i, j, k, d]
        ith = 1 + i + psizes(1) * (k + psizes(3) * (j + psizes(2) * d))
        !write(*, *) "proc ", i, j, k, d, ith

        do d0 = 1, lsizes(4)
        do k0 = 1, lsizes(3)
        do j0 = 1, lsizes(2)
        do i0 = 1, lsizes(1)
            count = count +  1
            val = count

            ig = i0+(j0-1)*lsizes(1)
            local(ig, k0, d0) = val

            ig = (start_indices(1) + i0) + (start_indices(2) + j0 - 1) * gsizes(1)
            global( ig, start_indices(3) + k0, start_indices(4) + d0) = val

        enddo
        enddo
        enddo
        enddo

        write(ithreadstr,'(i3.3)')  ith
        filename = "test."//ithreadstr//".dat"

        open(21, file=filename, form='binary', action='write')
            write(21) local
        close(21)
    enddo
    enddo
    enddo
    enddo

    filename = "test.global.dat"
    open(22, file=filename, form='binary', action='write')
    write(22) global
    close(22)

    deallocate(local)
    deallocate(global)

!======================================================================
    end program
!======================================================================
