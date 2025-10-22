!> mpas_pme_forcing_streams.F
module ocn_precip_factor_data
   use mpas_timekeeping
   use mpas_kind_types
   use mpas_log
   use mpas_stream_manager
   use netcdf
   use ocn_diagnostics_variables
   use mpi
   use ocn_mesh
   use ocn_constants
   use mpas_global_sum_mod

   implicit none
   private
   public :: pme_init_from_streams, pme_get_at_time, pme_finalize

   type(MPAS_streamManager_type), pointer :: sMgr => null()
   integer                 :: ncid=-1, vid_time=-1, vid_pmeN=-1, vid_pmeE=-1, vid_pmeS=-1
   integer :: nTimes=0, myrank=-1, nprocs=-1
   integer, allocatable    :: timeYMD(:)        ! YYYYMMDD (int)
   character(len=StrKIND) :: pme_file
   real(RKIND),   allocatable :: tsec(:)          ! seconds since ref
   real(RKIND),   allocatable :: pN(:), pE(:), pS(:)

   character(len=*), parameter :: CAL = 'noleap'

contains
   subroutine pme_init_from_streams(stream_mgr_in, stream_name, comm)
      type(MPAS_streamManager_type), target, intent(inout) :: stream_mgr_in
      character(len=*),               intent(in)     :: stream_name
      integer,                        intent(in)     :: comm
      character(len=NF90_MAX_NAME) :: dim_name

      integer :: ierr, root, ncid, did_t, tlen
      integer :: vid_time, vid_pN, vid_pE, vid_pS

      sMgr  => stream_mgr_in
      root = 0
      call MPI_Comm_rank(comm, myrank, ierr)

      if (myrank == root) then
         call mpas_get_stream_filename(sMgr, streamID=trim(stream_name),filename=pme_file,ierr=ierr)
         if (len_trim(pme_file) == 0) then
            call mpas_log_write('PME: failed to resolve filename for stream '//trim(stream_name))
            call MPI_Abort(comm, 2, ierr)
         end if

         ierr = nf90_open(trim(pme_file), NF90_NOWRITE, ncid); call nc_chk(ierr,'open')
         ierr = nf90_inq_dimid(ncid,'time',did_t);            call nc_chk(ierr,'inq_dimid(time)')
         ierr = nf90_inquire_dimension(ncid,did_t,dim_name,tlen);             call nc_chk(ierr,'inq_dimlen(time)')
         nTimes = tlen

         ierr = nf90_inq_varid(ncid,'time', vid_time);        call nc_chk(ierr,'inq_varid(time)')
         ierr = nf90_inq_varid(ncid,'pmeN', vid_pN);          call nc_chk(ierr,'inq_varid(pmeN)')
         ierr = nf90_inq_varid(ncid,'pmeE', vid_pE);          call nc_chk(ierr,'inq_varid(pmeE)')
         ierr = nf90_inq_varid(ncid,'pmeS', vid_pS);          call nc_chk(ierr,'inq_varid(pmeS)')

         allocate(timeYMD(nTimes), pN(nTimes), pE(nTimes), pS(nTimes))
         ierr = nf90_get_var(ncid, vid_time, timeYMD); call nc_chk(ierr,'get_var(time)')
         ierr = nf90_get_var(ncid, vid_pN,   pN);      call nc_chk(ierr,'get_var(pmeN)')
         ierr = nf90_get_var(ncid, vid_pE,   pE);      call nc_chk(ierr,'get_var(pmeE)')
         ierr = nf90_get_var(ncid, vid_pS,   pS);      call nc_chk(ierr,'get_var(pmeS)')
         ierr =  nf90_close(ncid); call nc_chk(ierr,'close file')
      end if

      call MPI_Bcast(nTimes, 1, MPI_INTEGER, root, comm, ierr)
      if (myrank /= root) then
         allocate(timeYMD(nTimes), pN(nTimes), pE(nTimes), pS(nTimes))
      end if

      call MPI_Bcast(timeYMD, nTimes, MPI_INTEGER, root, comm, ierr)
      call MPI_Bcast(pN,     nTimes, mpi_real8, root, comm, ierr)
      call MPI_Bcast(pE,     nTimes, mpi_real8, root, comm, ierr)
      call MPI_Bcast(pS,     nTimes, mpi_real8, root, comm, ierr)

      allocate(tsec(nTimes))
      call build_seconds_axis(timeYMD, CAL, tsec)

      if (myrank == root) then
         call mpas_log_write('PME: broadcast '//trim(itoa(nTimes))//' records from '//trim(pme_file))
      end if
   end subroutine pme_init_from_streams

   subroutine pme_get_at_time(domain, forcingPool)
      type(domain_type), intent(in)     :: domain

      type(mpas_pool_type), intent(inout) :: forcingPool
      integer              :: i0, i1, nCells, iCell, ierr
      real(kind=RKIND)                 :: w, dt01, dt0
      real(kind=RKIND), dimension(:), pointer :: evaporationFlux, rainFlux
      real(kind=RKIND) :: pSumE, pSumN, pSumS, eSumE, eSumN, eSumS
      !real(kind=RKIND) :: pSumEg, pSumNg, pSumSg, eSumEg, eSumNg, eSumSg
      real(kind=RKIND) :: corrN, corrE, corrS
      integer, dimension(2) :: reproIndexN, reproIndexE, reproIndexS
      integer, parameter :: nSums = 2
      real (kind=RKIND), dimension(nSums) :: reductions

      real (kind=RKIND), dimension(:,:), allocatable :: &
        pEReproSum, pSReproSum, pNReproSum, eEReproSum, eSReproSum, eNReprosum

      type(MPAS_Time_type) :: curTime
      real(RKIND) :: pmeN_t, pmeE_t, pmeS_t
      real(RKIND) :: latMax, latMin, dtr

      integer :: Y,M,D, h,mi,s
      real(RKIND) :: tmodel, t0, t1
      integer :: Y1,M1,D1, Ya,Ma,Da

      dtr = pi/180.0_RKIND
      latMin = config_fosi_pmt_latmIn*dtr
      latMax = config_fosi_pmt_latmax*dtr

      curTime = MPAS_get_clock_time(domain % clock, MPAS_NOW, ierr=ierr)
      ! Extract components from MPAS time (this exists in your branch)
      call mpas_get_time(curTime, YYYY=Y)
      call mpas_get_time(curTime, MM=M)
      call mpas_get_time(curTime, DD=D)
      call mpas_get_time(curTime, H =h)
      call mpas_get_time(curTime, M =mi)
      call mpas_get_time(curTime, S =s)

      tmodel = ymdhms_to_seconds(Y,M,D,h,mi,s, CAL)

      if (tmodel < tsec(1)) then
         ! Before first: interpolate between (first − 1 month) and first
         call split_ymd(timeYMD(1), Y1,M1,D1)
         call add_months(Y1,M1,D1, -1, Ya,Ma,Da, CAL)
         t0 = ymdhms_to_seconds(Ya,Ma, min(D1,days_in_month(Ya,Ma,CAL)), 0,0,0, CAL)
         t1 = tsec(1)
         i0 = nTimes; i1 = 1
      else if (tmodel >= tsec(nTimes)) then
         ! After last: last to (first + 1 month)
         call split_ymd(timeYMD(1), Y1,M1,D1)
         call add_months(Y1,M1,D1, +1, Ya,Ma,Da, CAL)
         t0 = tsec(nTimes)
         t1 = ymdhms_to_seconds(Ya,Ma, min(D1,days_in_month(Ya,Ma,CAL)), 0,0,0, CAL)
         i0 = nTimes; i1 = 1
      else
         call lower_upper_bounds(tmodel, tsec, i0, i1)
         t0 = tsec(i0); t1 = tsec(i1)
      end if

      w = max(0.0_RKIND, min(1.0_RKIND, (tmodel - t0) / max(1.0_RKIND, t1 - t0)))
      pmeN_t = (1.0_RKIND - w)*pN(i0) + w*pN(i1)
      pmeE_t = (1.0_RKIND - w)*pE(i0) + w*pE(i1)
      pmeS_t = (1.0_RKIND - w)*pS(i0) + w*pS(i1)

      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)

      nCells = nCellsOwned

      pSumE = 0.0_RKIND
      pSumN = 0.0_RKIND
      pSumS = 0.0_RKIND
      eSumE = 0.0_RKIND
      eSumN = 0.0_RKIND
      eSumS = 0.0_RKIND

      allocate (pNReproSum(nCells,2),pSReproSum(nCells,2),pEReproSum(nCells,2))
      allocate (eNReproSum(nCells,2),eSReproSum(nCells,2),eEReproSum(nCells,2))

      pNReproSum(:,:) = 0.0_RKIND
      pEReproSum(:,:) = 0.0_RKIND
      pSReproSum(:,:) = 0.0_RKIND
      eNReproSum(:,:) = 0.0_RKIND
      eEReproSum(:,:) = 0.0_RKIND
      eSReproSum(:,:) = 0.0_RKIND

      reproIndexN(1) = 1
      reproIndexN(2) = 0
      reproIndexE(1) = 1
      reproIndexE(2) = 0
      reproIndexS(1) = 1
      reproIndexS(2) = 0

      do iCell=1,nCells
         if (latCell(iCell) >= latMin .and. latCell(iCell) <= latMax) then
             reproIndexE(2) = reproIndexE(2) + 1
             pEReproSum(reproIndexE(2),1) = rainFlux(iCell)*areaCell(iCell)
!             pEReproSum(reproIndexE(2),2) = areaCell(iCell)
             eEReproSum(reproIndexE(2),1) = -evaporationFlux(iCell)*areaCell(iCell)
!             eEReproSum(reproIndexE(2),2) = areaCell(iCell)
         else if (latCell(iCell) > latMax) then
             reproIndexN(2) = reproIndexN(2) + 1
             pNReproSum(reproIndexN(2),1) = rainFlux(iCell)*areaCell(iCell)
!             pNReproSum(reproIndexN(2),2) = areaCell(iCell)
             eNReproSum(reproIndexN(2),1) = -evaporationFlux(iCell)*areaCell(iCell)
!             eNReproSum(reproIndexN(2),2) = areaCell(iCell)
         else
             reproIndexS(2) = reproIndexS(2) + 1
             pSReproSum(reproIndexS(2),1) = rainFlux(iCell)*areaCell(iCell)
 !            pSReproSum(reproIndexS(2),2) = areaCell(iCell)
             eSReproSum(reproIndexS(2),1) = -evaporationFlux(iCell)*areaCell(iCell)
 !            eSReproSum(reproIndexS(2),2) = areaCell(iCell)
         end if
      end do

      !now need to get a global P and E for each region
      ! just need area once for each region
      reductions(1) = mpas_global_sum(pEReproSum(:,2),domain%dminfo%comm, reproIndexE)
      reductions(2) = mpas_global_sum(pEReproSum(:,1),domain%dminfo%comm,reproIndexE)
      pSumE = reductions(2)!/(reductions(1) + 1.0e-20_RKIND)

      reductions(2) = mpas_global_sum(eEReproSum(:,1),domain%dminfo%comm,reproIndexE)
      eSumE = reductions(2)!/(reductions(1) + 1.0e-20_RKIND)

      !now North region
      reductions(1) = mpas_global_sum(pNReproSum(:,2),domain%dminfo%comm,reproIndexN)
      reductions(2) = mpas_global_sum(pNReproSum(:,1),domain%dminfo%comm,reproIndexN)
      pSumN = reductions(2)!/(reductions(1) + 1.0e-20_RKIND)

      reductions(2) = mpas_global_sum(eNReproSum(:,1),domain%dminfo%comm,reproIndexN)
      eSumN = reductions(2)!/(reductions(1) + 1.0e-20_RKIND)

      ! now south region
      reductions(1) = mpas_global_sum(pSReproSum(:,2),domain%dminfo%comm,reproIndexS)
      reductions(2) = mpas_global_sum(pSReproSum(:,1),domain%dminfo%comm,reproIndexS)
      pSumS = reductions(2)!/(reductions(1) + 1.0e-20_RKIND)

      reductions(2) = mpas_global_sum(eSReproSum(:,1),domain%dminfo%comm,reproIndexS)
      eSumS = reductions(2)!/(reductions(1) + 1.0e-20_RKIND)

      !compute all the correction terms
      corrE = (eSumE - (pmeN_t + pmeS_t))/(1.0E-20_RKIND + pSumE)
      corrN = (eSumN - pmeN_t)/(pSumN + 1.0E-20_RKIND)
      corrS = (eSumS - pmeS_t)/(pSumS + 1.0E-20_RKIND)

      do iCell=1,nCells
         if(latCell(iCell) >= latMax) then
            rainFlux(iCell) = corrN*rainFlux(iCell)
         else if(latCell(iCell) <= latMin) then
            rainFlux(iCell) = corrS*rainFlux(iCell)
         else
            rainFlux(iCell) = corrE*rainFlux(iCell)
         end if
      end do

   end subroutine pme_get_at_time


   subroutine pme_finalize()
      integer :: ierr
      if (allocated(timeYMD)) deallocate(timeYMD)
      if (allocated(pN))      deallocate(pN)
      if (allocated(pE))      deallocate(pE)
      if (allocated(pS))      deallocate(pS)
      if (ncid /= -1) then
        ierr = nf90_close(ncid)
      end if

      ncid = -1
   end subroutine pme_finalize

     !----------------- helpers (portable calendar math) -----------------

   subroutine build_seconds_axis(ymd, cal, outsec)
      integer, intent(in) :: ymd(:)
      character(len=*), intent(in) :: cal
      real(RKIND), intent(out) :: outsec(:)
      integer :: k, y,m,d
      do k=1,size(ymd)
         call split_ymd(ymd(k), y,m,d)
         outsec(k) = ymdhms_to_seconds(y,m,d, 0,0,0, cal)
      end do
   end subroutine

   pure subroutine split_ymd(ival, y,m,d)
      integer, intent(in)  :: ival
      integer,        intent(out) :: y,m,d
      y = ival / 10000
      m = mod(ival/100, 100)
      d = mod(ival, 100)
   end subroutine

   pure subroutine lower_upper_bounds(t, axis, i0,i1)
      real(RKIND), intent(in)  :: t
      real(RKIND), intent(in)  :: axis(:)
      integer, intent(out) :: i0,i1
      integer :: k
      do k=2,size(axis)
         if (t < axis(k)) then
            i0 = k-1; i1 = k
            return
         end if
      end do
      i0 = size(axis); i1 = size(axis)   ! fallback
   end subroutine

   pure function ymdhms_to_seconds(y,m,d, hh,mm,ss, cal) result(secs)
      integer, intent(in) :: y,m,d, hh,mm,ss
      character(len=*), intent(in) :: cal
      integer :: yy, j
      real(RKIND) :: secs
      real(RKIND), parameter :: SEC_DAY = 86400.0_RKIND

      yy = y-1
      if (trim(cal) == 'noleap') then
         secs = SEC_DAY * (365.0_RKIND*real(yy,RKIND))
      else
         secs = SEC_DAY * ( real(365*yy,RKIND) + real(yy/4 - yy/100 + yy/400,RKIND) )
      end if
      do j=1,m-1
         secs = secs + SEC_DAY * real(days_in_month(y,j,cal), RKIND)
      end do
      secs = secs + SEC_DAY*real(d-1,RKIND) + real(hh*3600 + mm*60 + ss, RKIND)
   end function

   pure function days_in_month(y,m, cal) result(nd)
      integer, intent(in) :: y,m
      character(len=*), intent(in) :: cal
      integer :: nd
      select case (m)
      case (1,3,5,7,8,10,12); nd = 31
      case (4,6,9,11);        nd = 30
      case (2)
         if (trim(cal) == 'noleap') then
            nd = 28
         else
            if ( (mod(y,400)==0) .or. (mod(y,4)==0 .and. mod(y,100)/=0) ) then
               nd = 29
            else
               nd = 28
            end if
         end if
      end select
   end function

   pure subroutine add_months(y0,m0,d0, nmon, y1,m1,d1, cal)
      integer, intent(in)  :: y0,m0,d0, nmon
      integer, intent(out) :: y1,m1,d1
      character(len=*), intent(in) :: cal
      integer :: totm
      totm = (y0*12 + (m0-1)) + nmon
      y1 = totm / 12
      m1 = mod(totm,12) + 1
      d1 = min(d0, days_in_month(y1,m1,cal))
   end subroutine

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=32) :: s
      write(s,'(I0)') i
   end function

   subroutine nc_chk(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where
      if (ierr /= nf90_noerr) then
         call mpas_log_write('NetCDF error in '//trim(where)//': '//trim(nf90_strerror(ierr)))
         stop 2
      end if
   end subroutine

   end module ocn_precip_factor_data

