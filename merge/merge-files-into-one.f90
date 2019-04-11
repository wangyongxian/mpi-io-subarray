#ifdef REAL8

#define REAL_T      real*8
#define MPI_REAL_T  MPI_DOUBLE_PRECISION

#else

#define REAL_T      real*4
#define MPI_REAL_T  MPI_REAL

#endif

!======================================================================
program main
!======================================================================
    use mpi
    implicit none

    integer, parameter :: ndims = 4
    integer, parameter, dimension(1:ndims) :: gsizes = (/1020, 6, 1020, 50/)
    integer, parameter, dimension(1:ndims) :: lsizes = (/ 170, 1,   85, 50/)
    integer, parameter, dimension(1:ndims) :: psizes = (/   6, 6,   12,  1/)

    integer, parameter :: il = lsizes(1) * lsizes(2)
    integer, parameter :: jl = lsizes(3)
    integer, parameter :: kl = lsizes(4)

    REAL_T, allocatable :: subarray(:, :, :)

!===================================================================================
    character*3  :: ithreadstr
    character*10 :: seqstr
    character*120 :: local_file, global_file
    integer, dimension(1:ndims) :: coords, start_indices
    integer :: subarray_size
    integer :: nprocs, rank
    integer :: ierr, ipoint
    integer :: filetype
    LOGICAL, PARAMETER :: reorder = .FALSE.
    LOGICAL, PARAMETER :: periodic(1:ndims) = .FALSE.
    integer :: comm
    integer :: ith

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank,   ierr)
    if (nprocs < product(psizes)) stop " nprocs >= 432 is expected !!"


    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, psizes, &
    &       periodic, reorder, comm, ierr)
    call MPI_CART_COORDS(comm, rank, ndims, coords, ierr)

    start_indices(:) = coords(:) * lsizes(:)
    subarray_size = product(lsizes)

    call MPI_TYPE_CREATE_SUBARRAY(ndims,                &
    &       gsizes, lsizes, start_indices,              &
    &       MPI_ORDER_FORTRAN, MPI_REAL_T, filetype, ierr)
    call MPI_TYPE_COMMIT(filetype, ierr)

    !--- 'ith' can be derived from coordinates of current MPI rank
    ith = 1 + coords(1) +           &
          psizes(1) * (coords(3) +  &
          psizes(3) * (coords(2) +  &
          psizes(2) * (coords(4))))
    write(ithreadstr,'(i3.3)')  ith

    allocate( subarray(il, jl, kl) )

    do ipoint = 36, 36, 36

        call get_file_seq_no(ipoint, seqstr)

        local_file = '../../model_result/THETA.'//seqstr//'.'//ithreadstr//'.001.data'
        !global_file = '../../model_result/THETAfhl.'//seqstr//'.data'
        global_file = '../data/THETAfhl.'//seqstr//'.data'

        call merge_files(local_file, global_file, filetype)

    enddo

    deallocate(subarray)

    call MPI_TYPE_FREE(filetype, ierr)
    call MPI_FINALIZE(ierr)

!======================================================================
contains
!======================================================================

    subroutine merge_files(readfile, writefile, datatype)
        implicit none
        character(len=*), intent(in) :: readfile, writefile
        integer, intent(in) :: datatype

        integer :: fh, total
        integer :: status(MPI_STATUS_SIZE)

        !------- read a single file per MPI rank using normal IO

        open(21, file=readfile, &
        &        form='binary', action='read', status='old')
        read(21) subarray
        close(21)

        !------- write file using MPI IO

        call MPI_FILE_OPEN(MPI_COMM_WORLD, writefile,       &
        &       MPI_MODE_WRONLY + MPI_MODE_CREATE,          &
        &       MPI_INFO_NULL, fh, ierr )
        call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND,       &
        &       MPI_REAL_T, datatype, "native",             &
        &       MPI_INFO_NULL, ierr)
        call MPI_FILE_WRITE_ALL(fh, subarray,               &
        &       subarray_size, MPI_REAL_T, status, ierr)
        ! call MPI_GET_COUNT( status, MPI_REAL_T, total, ierr)
        ! write(*, *) 'process ', rank, 'write ', total, 'number'
        call MPI_FILE_CLOSE(fh, ierr)

    end subroutine

    subroutine get_file_seq_no(ipoint, seqstr)
        implicit none
        integer, intent(in) :: ipoint
        character*10, intent(out) :: seqstr

        integer*4,        external :: point
        double precision, external :: relday
        double precision           :: tim
        integer, parameter         :: dt = 100
        integer     :: idaypmon(0:12),numpoint
        integer     :: iyear,imonth,iday,ihour,imin,isec
        character*4 :: chyear
        character*2 :: chmonth,chday,chhour
        data idaypmon/0,31,28,31,30,31,30,31,31,30,31,30,31/

        tim=(relday(0,0,0,1,1,1980)+ipoint*dt)/(24.0d0*3600.0d0)
        call unjuln(isec,imin,ihour,iday,imonth,iyear,tim)
        ! write(*,*)isec,imin,ihour,iday,imonth,iyear,tim
        write(chyear,'(i4)')        iyear
        write(chmonth,'(i2.2)')     imonth
        write(chday,'(i2.2)')       iday
        write(chhour,'(i2.2)')      ihour
        numpoint=point(iyear,imonth,iday,ihour)
        write(seqstr,'(I10.10)')    numpoint

    end subroutine
end program
!======================================================================


integer*4 function point(iyear,imonth,iday,ihour)
    integer iyear,imonth,iday,ihour
    double precision,external :: relday
    double precision tim1,tim2,tim
    tim1=relday(0,0,0,1,1,1980)
    tim2=relday(0,0,ihour,iday,imonth,iyear)
    tim=tim2-tim1
    point=int(tim/100.D0)
!       write(*,*)tim1,tim2,tim,point
    return
end function
!======================================================================

double precision function relday(jsec,jmin,jhour,iday,imonth,iyear)
!     ==================================================
!
! purpose:
!  ------
!        produce the time relative to 0000z "julref" of the current date
!        and time. the units are second, so hours and minutes transform to
!        fractions of a day.
!
!  arguments:
!  --------
!        jmin     : minutes for the current time
!        jhour    : hour of the current time
!        iday     : day of the current date
!        imonth   : month of the current date
!        iyear    : year of the current date (e.g. 1961)
!          e.g. 1340z on 21 sep 1987 is called by (40,13,21,9,1987)
!
      parameter (igreg=15+31*(10+12*1582),julref=2415021)
      if (iyear.lt.0) iyear=iyear+1
      if (imonth.gt.2) then
        jy=iyear
        jm=imonth+1
      else
        jy=iyear-1
        jm=imonth+13
      endif
      jultmp=int(365.25d0*jy)+int(30.6001d0*jm)+iday+1720995
      if (iday+31*(imonth+12*iyear).ge.igreg) then
        ja=int(0.01d0*jy)
        jultmp=jultmp+2-ja+int(0.25d0*ja)
      endif
      relday =((jultmp-julref) + dble((60*jhour+jmin)*60+jsec)/86400.d0)*86400.d0
      return
end function
!======================================================================
subroutine unjuln( jsec,jminut,jhour,jday,jmonth,jyear, relday )
!     ================================================================
!
!  Purpose:
!     Take the relative time in days and re-express in terms of
!     seconds, minutes, hours, days, month, year.
!
!
      parameter (igreg=2299161,julref=2415021)
!
      double precision rday,relday
!
      rday = relday
      jsec   = nint( 86400.d0*mod(rday,1.d0) )
      if( jsec.lt.0. ) jsec = 86400.+jsec
      jhour  = jsec/3600
      jminut = (jsec-3600*jhour)/60
      jsec   = mod(jsec,60)
!
      julian = julref + int(rday)
      if( rday.lt.0. ) julian = julian-1
      if( julian.ge.igreg ) then
        j1 = int((dble(julian-1867216)-0.25d0)/36524.25d0)
        ja=julian+1+j1-int( (0.25d0*j1) )
      else
        ja = julian
      endif
!
      jb = ja+1524
      jc=int(6680.d0+(dble(jb-2439870)-122.1d0)/365.25d0)
      jd = 365*jc+int(0.25d0*jc)
      je=int(dble(jb-jd)/30.6001d0)
      jday=jb-jd-int(30.6001d0*je)
      jmonth=je-1
      if(jmonth.gt.12) jmonth=jmonth-12
      jyear=jc-4715
      if(jmonth.gt.2) jyear=jyear-1
      if(jyear.le.0) jyear=jyear-1
      return
end subroutine
!====================================================================
SUBROUTINE CALLEAPYR(YR,DM)
  IMPLICIT NONE
  INTEGER :: YR,DM
  IF(MOD(YR,100).EQ.0)THEN
    IF(MOD(YR,400).EQ.0)THEN
      DM=29
    ELSE
      DM=28
    ENDIF
  ELSE
    IF(MOD(YR,4).EQ.0)THEN
      DM=29
    ELSE
      DM=28
    ENDIF
  ENDIF
  RETURN
END SUBROUTINE CALLEAPYR
!====================================================================
