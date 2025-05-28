MODULE CHUNK
CONTAINS
  SUBROUTINE Make2DChunks(Nxchunks, Nychunks, &
    &   N, Niv, ivchunk, rangechunk, Nxptchunks, Nivxpt, &
    &   ivxptchunk, rangexptchunk)!, ixychunk, Nixy)
    Use Dim, ONLY: nx, ny, nxpt
    Use Indexes, ONLY: igyl
    Use Lsode, ONLY: neq
    Use Xpoint_indices, ONLY: iysptrx1, ixpt1, ixpt2
    IMPLICIT NONE
    INTEGER, INTENT(INOUT):: N, Nxptchunks(nxpt), Nxchunks, Nychunks
!    INTEGER, ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT):: ixychunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(OUT):: rangexptchunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT)::  ivxptchunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT)::  ivchunk, &
    &               rangechunk, Nivxpt
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)::  Niv !, Nixy
    INTEGER:: Nmax!, Nixymax
    integer:: ix, iy, nxi, nyi, ii, idx(2), idxl, ixpt, iix, iicut, iscut
    real:: dx, dy, dyxpt(2)
    integer, allocatable:: xlims(:,:), ylims(:,:), xlimsxpt(:,:,:,:), ylimsxpt(:,:,:,:)
    ! TODO: Resize & use local variables, allocate real variables at last step only 

    ! Ensure chunking setup is valid: between 1 and n(x/y) chunks
    Nxchunks = MAX(MIN(nx, Nxchunks),1)
    Nychunks = MAX(MIN(ny-1, Nychunks),1)
    ! Calculate the number of chunks needed maximum
    N = Nxchunks * Nychunks + 1
    do ixpt = 1, nxpt
        Nxptchunks(ixpt) = MIN(Nxptchunks(ixpt),iysptrx1(ixpt)-1) ! Number of X-point chunks
        dyxpt(ixpt) = real(iysptrx1(ixpt)-1)/(Nxptchunks(ixpt)+1)
    end do
    Nmax = neq
!    Nixymax = (nx+2)*(ny+2)
    ! Get the dx/dy per chunk
    dx = real(nx)/(Nxchunks)
    dy = real(ny-1)/(Nychunks+1)
    ! Allocate the necassary arrays
    allocate( xlims(Nxchunks,2), ylims(Nychunks,2), xlimsxpt(nxpt, 2, MAXVAL(Nxptchunks),2), &!ixychunk(N, Nixymax, 2), &
    &   ivchunk(N, Nmax), rangechunk(N, 4), Niv(N), ylimsxpt(nxpt, 2, MAXVAL(Nxptchunks),2), &
    &   ivxptchunk(nxpt,MAXVAL(Nxptchunks),Nmax), Nivxpt(nxpt,MAXVAL(Nxptchunks)), &
    &   rangexptchunk(nxpt,2,MAXVAL(Nxptchunks),4))!, Nixy(N))
    ! Calculate poloidal X-point intervals
    do ixpt = 1, nxpt ! Separate teams for each X-point
        ! Calculate Y-chunks
        ylimsxpt(ixpt,1,1,1) = 2; ylimsxpt(ixpt,1,1,2)=1+max(1,int(dyxpt(ixpt)))
        ylimsxpt(ixpt,2,1,1) = 2; ylimsxpt(ixpt,2,1,2)=1+max(1,int(dyxpt(ixpt)))
        do iix = 2, Nxptchunks(ixpt)-1 ! Number of threads per Xpt - 1 for now
            ! Left cut
            ylimsxpt(ixpt, 1, iix, 1) = ylimsxpt(ixpt, 1, iix-1, 2)+1
            ylimsxpt(ixpt, 1, iix, 2) = 2+int(dyxpt(ixpt)*iix)
            ! Right cut
            ylimsxpt(ixpt, 2, iix, 1) = ylimsxpt(ixpt, 2, iix-1, 2)+1
            ylimsxpt(ixpt, 2, iix, 2) = 2+int(dyxpt(ixpt)*iix)
        end do
        if (Nxptchunks(ixpt) .gt. 1) then
            ylimsxpt(ixpt, 1, Nxptchunks(ixpt),1) = ylimsxpt(ixpt,1,Nxptchunks(ixpt)-1,2)+1
            ylimsxpt(ixpt, 2, Nxptchunks(ixpt),1) = ylimsxpt(ixpt,2,Nxptchunks(ixpt)-1,2)+1
        endif
        ylimsxpt(ixpt, 1, Nxptchunks(ixpt),2) = iysptrx1(ixpt)
        ylimsxpt(ixpt, 2, Nxptchunks(ixpt),2) = iysptrx1(ixpt)
        ! Set X-chunks
        do iix = 1, Nxptchunks(ixpt)
            ! Left cut
            xlimsxpt(ixpt, 1, iix, 1) = ixpt1(ixpt)
            xlimsxpt(ixpt, 1, iix, 2) = ixpt1(ixpt)+1
            ! Right cut
            xlimsxpt(ixpt, 2, iix, 1) = ixpt2(ixpt)
            xlimsxpt(ixpt, 2, iix, 2) = ixpt2(ixpt)+1
        end do
        ! Create ranges 
        do iix = 1, Nxptchunks(ixpt)
                ! Left cut
                rangexptchunk(ixpt, 1, iix, 1) = xlimsxpt(ixpt, 1, iix, 1)
                rangexptchunk(ixpt, 1, iix, 2) = xlimsxpt(ixpt, 1, iix, 2)
                rangexptchunk(ixpt, 1, iix, 3) = ylimsxpt(ixpt, 1, iix, 1)
                rangexptchunk(ixpt, 1, iix, 4) = ylimsxpt(ixpt, 1, iix, 2)
                ! Right cut
                rangexptchunk(ixpt, 2, iix, 1) = xlimsxpt(ixpt, 2, iix, 1)
                rangexptchunk(ixpt, 2, iix, 2) = xlimsxpt(ixpt, 2, iix, 2)
                rangexptchunk(ixpt, 2, iix, 3) = ylimsxpt(ixpt, 2, iix, 1)
                rangexptchunk(ixpt, 2, iix, 4) = ylimsxpt(ixpt, 2, iix, 2)
        end do
    end do
    ! TODO: Ensure one parallel chunk spans the whole iy=0 boundary!
    ! Calculate poloidal chunking intervals
    ! Include protections for boundaries
    xlims(1,1) = 0; xlims(1,2)=max(1,int(dx))
    do ix = 2, Nxchunks-1
        xlims(ix,1) = xlims(ix-1,2)+1
        xlims(ix,2) = int(dx*ix)
    end do
    if (Nxchunks .gt. 1) xlims(Nxchunks,1) = xlims(Nxchunks-1,2)+1
    xlims(Nxchunks,2) = nx+1
    ! Calculate radial chunking intervals
    ! Include protections for boundaries
    ylims(1,1) = 2; ylims(1,2)=2+int(dy)
    do iy = 2, Nychunks-1
        ylims(iy,1) = ylims(iy-1,2)+1
        ylims(iy,2) = 2+int(dy*iy)
    end do
    if (Nychunks .gt. 1) ylims(Nychunks,1) = ylims(Nychunks-1,2)+1
    ylims(Nychunks,2) = ny+1
    ! Ravel the ranges into chunks 
    rangechunk(1,1) =   0
    rangechunk(1,2) =   nx+1
    rangechunk(1,3) =   0
    rangechunk(1,4) =   1
    do ix = 1, Nxchunks
        do iy = 1, Nychunks
            rangechunk(Nxchunks*(iy-1) + ix + 1, 1) = xlims(ix,1)
            rangechunk(Nxchunks*(iy-1) + ix + 1, 2) = xlims(ix,2)
            rangechunk(Nxchunks*(iy-1) + ix + 1, 3) = ylims(iy,1)
            rangechunk(Nxchunks*(iy-1) + ix + 1, 4) = ylims(iy,2)
        end do
    end do
    ! Get the yldot-indices for the chunks 
    Niv = 0
    Nivxpt = 0
    ivchunk = 0
    ivxptchunk = 0
    iscut = 0
    do ii = 1, neq
        idx = igyl(ii,:)
        ! Get the X-point cut jacobian indices
        ! These loops will typically be much smaller than the main chunks
        do ixpt = 1, nxpt
            do iix = 1, Nxptchunks(ixpt)
                do iicut = 1, 2 
                    if (    (idx(1).ge.rangexptchunk(ixpt,iicut,iix,1)) &
                    &  .and.(idx(1).le.rangexptchunk(ixpt,iicut,iix,2)) &
                    &  .and.(idx(2).ge.rangexptchunk(ixpt,iicut,iix,3)) &
                    &  .and.(idx(2).le.rangexptchunk(ixpt,iicut,iix,4)) &
                    & ) then
                        Nivxpt(ixpt,iix) = Nivxpt(ixpt,iix) + 1
                        ivxptchunk(ixpt,iix,Nivxpt(ixpt,iix)) = ii
                        iscut=1
                    endif
                end do
            end do
        end do
        do iix = 1, N
            if ((idx(1).ge.rangechunk(iix,1)).and.(idx(1).le.rangechunk(iix,2)) .and. &
            &   (idx(2).ge.rangechunk(iix,3)).and.(idx(2).le.rangechunk(iix,4)) &
            &   .and.(iscut.ne.1)) then
                Niv(iix) = Niv(iix) + 1
                ivchunk(iix, Niv(iix)) = ii
            endif
        end do
        iscut=0
    enddo
!    write(*,*) SUM(Nivxpt)+SUM(Niv), MAXVAL(Niv), MAXVAL(Nivxpt)
    ! Do spatial chunking - no longer needed
!    Nixy = 0
!    ixychunk = 0
!    do ii = 1, N
!        do ix = rangechunk(ii,1), rangechunk(ii,2)
!            do iy = rangechunk(ii,3), rangechunk(ii,4)
!                Nixy(ii) = Nixy(ii) + 1
!                ixychunk(ii, Nixy(ii),1) = ix
!                ixychunk(ii, Nixy(ii),2) = iy
!            end do
!        end do
!    end do
    ! Update max array size required
!    Nmax = MAXVAL(Niv)
!    Nixymax = MAXVAL(Nixy)
!    ! Trim overly large array
!    allocate( xlims(Nxchunks,2), ylims(Nychunks,2), ixychunk(N, Nixymax, 2), &
!    &   ivchunk(N, Nmax), rangechunk(N, 4), Niv(N), Nixy(N))

  END SUBROUTINE Make2DChunks
END MODULE CHUNK




SUBROUTINE InitOMPJac()
    USE OMPJacSettings
    USE OMPJac
    USE ParallelSettings, ONLY: OMPParallelPandf1,Nthreads
    USE Jacobian, ONLY: nnzmx
    USE Lsode, ONLY: neq
    IMPLICIT NONE

    integer:: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS
#ifndef _OPENMP
        call xerrab( &
            "UEDGE was not compiled with OpenMP. Cannot use OMP parallel features.")
#else
        if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp, &
                    ' Max number of threads available:',OMP_GET_MAX_THREADS()
            call OMP_SET_NUM_THREADS(OMP_GET_MAX_THREADS())
              !$omp parallel
                if (OMP_GET_THREAD_NUM().eq.0) then
                    if (Nthreads.gt.OMP_GET_MAX_THREADS()) then
                        if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp, &
                                ' Warning: # threads requested > # threads available'
                            Nthreads=OMP_GET_NUM_THREADS()
                    endif
                    if (Nthreads.le.0) then
                        call xerrab('Nthread must be >0')
                    endif
                    if (OMPJacVerbose.gt.0) write(*,'(a,a,i3)') &
                            OMPJacStamp,' Number of threads for omp calculations:', &
                            Nthreads
                endif
              !$omp END parallel
            if (OMPJacNchunks.eq.0) then
                NchunksJac=neq
            elseif (OMPJacNchunks.lt.0) then
                NchunksJac=Nthreads
            else
                NchunksJac=OMPJacNchunks
            endif
 
            if (Nthreads.gt.1) then
                nnzmxperchunk=ceiling(real(nnzmx)/real(NchunksJac))*omplenpfac !nnzmx=neq*lenfac
            else
                nnzmxperchunk=ceiling(real(nnzmx)/real(NchunksJac))*omplenpfac !nnzmx=neq*lenfac
            endif

            if (OMPJacVerbose.gt.0) &
                write(*,"(a,a,i4,a,i6,a,i5,a,i6)") TRIM(OMPJacStamp), &
                    ' Nthreads:', Nthreads, ' NchunksJac:', &
                    NchunksJac, ' nnzmxperchunk:',nnzmxperchunk,' neq:',neq
            call gchange('OMPJac',0)
#endif    
    RETURN
END SUBROUTINE InitOMPJAC


SUBROUTINE InitOMPPandf1
    USE OMPPandf1Settings, ONLY: OMPPandf1Nxchunks,OMPPandf1Nychunks,OMPPandf1Stamp,OMPPandf1Verbose
    USE OMPPandf1, ONLY: NchunksPandf1,Nxchunks,Nychunks
    USE Dim, ONLY: ny
    IMPLICIT NONE

    if (1.eq.0) then
    if (OMPPandf1Nychunks.lt.0) then
        call xerrab('Nychunks<0. Nxchunks must be >=0')
    endif
    if (OMPPandf1Nychunks.eq.0) then
        Nychunks=ny
    else
        Nychunks=OMPPandf1Nychunks
    endif
    if (OMPPandf1Nxchunks.ne.1) then
        call xerrab('OMPPandf1Nxchunks!=1. Only Nxchunks=1 is implemented for the moment...')
    else
        Nxchunks=1
    endif
    endif

    ! this is a placeholder for further parallelization but we need to implement handling of x-points discon. in ix indexing
    NchunksPandf1=Nychunks

    call gchange('OMPPandf1',0)

    if (OMPPandf1Verbose.gt.0) then
        write(*,*) OMPPandf1Stamp, ' NchunksPandf1 = ',NchunksPandf1
    endif
    RETURN
END SUBROUTINE InitOMPPandf1


SUBROUTINE InitZeroOMP
    IMPLICIT NONE
#ifdef _OPENMP
        call OmpInitZerobbb
        call OmpInitZeroapi
        call OmpInitZerocom
#endif
    RETURN
END SUBROUTINE InitZeroOMP


#ifdef _OPENMP

  SUBROUTINE OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)
    USE OMPJac, ONLY: iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum,OMPTimeJacRow,NchunksJac
    USE OMPJacSettings, ONLY: OMPJacVerbose,OMPJacStamp,OMPTimingJacRow
    USE ParallelDebug, ONLY: OMPJacDebug
    IMPLICIT NONE
    integer,intent(in):: neq
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian
    real,intent(out)   :: rcsc(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: icsc(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: jcsc(neq+1)   ! pointers to beginning of each row in jac,ja
    integer,intent(out):: nnzcumout
    integer ichunk
    integer:: iunit,iv
    nnzcum(1:NchunksJac)=-1
    nnzcum(1)=nnz(1)-1

    do ichunk=2,NchunksJac
        nnzcum(ichunk)=nnzcum(ichunk-1)+nnz(ichunk)-1
    enddo
    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' nnz:',nnz(1:NchunksJac)
        write(*,*) OMPJacStamp,' nnzcum:',nnzcum(1:NchunksJac)
    endif
    if (OMPJacVerbose.gt.0) write(*,'(a,i9)') '**** Number of non-zero Jacobian elems:',nnzcum(NchunksJac)

    if (nnzcum(NchunksJac).gt.nnzmx) then
        write(*,*) 'nnzcum=',nnzcum
        write(*,*) nnzmx
        call xerrab(' Problem: nnzcum > nnzmx...')
    endif

    jcsc(OMPivmin(1):OMPivmax(1))= iJacRow(OMPivmin(1):OMPivmax(1))
    do ichunk=2,NchunksJac
        jcsc(OMPivmin(ichunk):OMPivmax(ichunk))= iJacRow(OMPivmin(ichunk):OMPivmax(ichunk))+nnzcum(ichunk-1)
    enddo

    rcsc(1:nnz(1)-1)= rJacElem(1:nnz(1)-1,1)
    icsc(1:nnz(1)-1)= iJacCol(1:nnz(1)-1,1)
    do ichunk=2,NchunksJac
        rcsc(nnzcum(ichunk-1)+1:nnzcum(ichunk))=rJacElem(1:nnz(ichunk)-1,ichunk)
        icsc(nnzcum(ichunk-1)+1:nnzcum(ichunk))=iJacCol(1:nnz(ichunk)-1,ichunk)
    enddo
    nnzcumout=nnzcum(NchunksJac)
    if (OMPTimingJacRow.gt.0) then
        OPEN(newunit = iunit, file = 'omptiming.dat')
        do iv=1,neq
            write(iunit,*) iv,OMPTimeJacRow(iv)
        enddo
        CLOSE(iunit)
    endif
    RETURN
  END SUBROUTINE OMPCOllectJacobian


  SUBROUTINE jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    !   Calculate Jacobian matrix (derivatives with respect to each
    !   dependent variable of the right-hand side of each rate equation).
    !   Lower and upper bandwidths are used to select for computation
    !   only those Jacobian elements that may be nonzero.
    !   Estimates of Jacobian elements are computed by finite differences.
    !   The Jacobian is stored in compressed sparse row format.

    USE Timing, ONLY: istimingon,ttjstor,ttotjf,ttimpjf
    USE PandfTiming
    USE Grid, ONLY: ngrid,ig,ijac,ijactot
    USE Jacobian_csc, ONLY: rcsc,jcsc,icsc,yldot_pert,yldot_unpt
    USE OMPJac, ONLY: NchunksJac,nnzmxperchunk
    USE ParallelSettings, ONLY: Nthreads, OMPParallelPandf1
    USE OMPJacSettings, ONLY: OMPJacVerbose,OMPCheckNaN,&
            OMPLoadBalance,OMPAutoBalance,OMPJacStamp,OMPBalanceStrength
    USE ParallelDebug, ONLY: WriteJacobian,OMPJacDebug
    USE Flags, ONLY: iprint
    USE OMPJac, ONLY:iJacCol,rJacElem,iJacRow,OMPivmin,OMPivmax,nnz,nnzcum,OMPLoadWeight,OMPTimeLocalJac
    USE UEpar, ONLY: svrpkg
    USE Math_problem_size, ONLY:neqmx
    IMPLICIT NONE
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(in)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja

    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine
    integer,allocatable :: iJacConstructor(:,:)
    real,allocatable:: rJacConstructor(:,:)
    integer:: nnzcumout
    ! ... Functions
    logical tstguardc
    real tick,tock
    external tick,tock

    ! ... Local variables:
    real tsjstor, tsimpjf, dtimpjf,time0,time1
    integer:: i,thread,ichunk,iv,TID, OMP_GET_THREAD_NUM
    character(len = 80) ::  filename

    OMPTotTimeJacCalc = tick()
    ! Calculate load distribution for threads
    if (OMPLoadBalance.ne.1 .and. OMPAutoBalance.ne.1) then
        OMPLoadWeight(1:NchunksJac)=1.0
    endif
    if (OMPAutoBalance.eq.1) then
        !Check that time are not zero
        if (OMPBalanceStrength<=0) call xerrab('OMPBalanceStrength must be >0')
        if (minval(OMPTimeLocalJac).gt.0.0) then
            do i=1,NchunksJac
                OMPLoadWeight(i)=OMPLoadWeight(i)*1/(OMPTimeLocalJac(i) &
                    /sum(OMPTimeLocalJac)*real(NchunksJac))**OMPBalanceStrength
            enddo
        else
            OMPLoadWeight(1:NchunksJac)=1.0
        endif
    endif
    !   Get the range of the iv index for each thread
    call OMPSplitIndex(1,neq,NchunksJac,OMPivmin,OMPivmax,OMPLoadWeight)

    if (OMPJacVerbose.gt.0) then
        write(*,*) ' *OMPJac* neq=',neq,neqmx
        write(*,*) &
            ' *OMPJac* Ivmin(ichunk),Ivmax(ichunk), OMPLoadWeight(ichunk) : OMPTimeLocalJac(ichunk) ***'
        do ichunk=1,NchunksJac
            write(*,'(a,I3,a,I7,I7,f6.2,a,f6.2)') '  * ichunk ', ichunk,':', &
                OMPivmin(ichunk),OMPivmax(ichunk),OMPLoadWeight(ichunk), &
                ' : ',OMPTimeLocalJac(ichunk)
        enddo
    endif

    OMPTimeJacCalc= tick()
    !   Get initial value of system cpu timer.
    if (istimingon .eq. 1) tsjstor = tick()

    !   Count Jacobian evaluations, both for total and for this case
    ijactot = ijactot + 1
    ijac(ig) = ijac(ig) + 1
    if ((svrpkg.eq.'nksol') .and. (iprint .ne. 0)) then
        write(*,'(a,i4,a,i6,a,i9)') ' Updating OMP Jacobian [', &
                    Nthreads,'|',NchunksJac, ']: npe = ', ijac(ig)
    endif
    !   Set up diagnostic arrays for debugging
    do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  ! for diagnostic only
        yldot_pert(iv) = 0.
    enddo

    !   build jacobian ##############################################################
    OMPTimeBuild=tick()
    nnz(1:NchunksJac)=-1
    call OMPJacBuilder(neq, t, yl,yldot00, ml, mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    OMPTotTimebuild = OMPTotTimeBuild+tock(OMPTimeBuild)
    if (OMPJacVerbose.gt.0) write(*,*)OMPJacStamp,' Time to build jac:',OMPTimeBuild
    !   end build jacobian ##############################################################

    !   collect jacobian ##############################################################
    OMPTimeCollect=tick()
    call OMPCollectJacobian(neq,nnzmx,rcsc,icsc,jcsc,nnzcumout)
    OMPTotTimeCollect = OMPTotTimeCollect+tock(OMPTimeCollect)
    if (OMPJacVerbose.gt.0) write(*,*)OMPJacStamp,' Time to collect jac:',OMPTimeCollect
    !   end collect jacobian ##############################################################

    jcsc(neq+1) = nnzcumout+1 ! This is set here out of OMPJAcCollect for compatibility with hybrid jac_calc

    !   for Debug purpose
    if (WriteJacobian.eq.1) then
        write(filename,'(a,3i3,a)') "jac_omp_",ijac(ig),".txt"
        call jac_write(filename,neq, rcsc, icsc, jcsc)
    endif

    if (OMPCheckNaN.gt.0) then
        do i=1,nnzmx
            if (isnan(rcsc(i))) then
                write(*,*) 'rcsc is NaN at i=',i,rcsc(i)
            endif
            if (icsc(i).ne.icsc(i)) then
                write(*,*) 'icsc is NaN at i=',i,icsc(i)
            endif
        enddo
        do i=1,neq+1
            if (jcsc(i).ne.jcsc(i)) then
                write(*,*) 'jcsc is NaN at i=',i
            endif
        enddo
    endif

    !   Convert Jacobian from compressed sparse column to compressedsparse row format.
    time1=tick()
    call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

    if (istimingon .eq. 1) ttjstor = ttjstor + tock(tsjstor)

    if ((OMPJacVerbose.gt.0) .and. (iprint .ne. 0)) &
        write(*,'(a,1pe9.2,a)') '  OMP Jac timing:', tock(OMPTimeJacCalc), 's'
        
    RETURN
  END SUBROUTINE jac_calc_omp


  SUBROUTINE OMPJacBuilder(neq, t, yl,yldot00, ml,mu,wk,iJacCol,rJacElem,iJacRow,nnz)
    USE OMPJacSettings, ONLY: OMPJacStamp,OMPJacVerbose,OMPLoopJacNchunk
    USE ParallelDebug, ONLY: OMPJacDebug,OMPCopyArray,OMPCopyScalar
    USE OMPJac, ONLY: OMPivmin,OMPivmax,OMPTimeLocalJac,OMPTimeJacRow,nnzmxperchunk,NchunksJac
    USE Selec, ONLY: yinc,xlinc,xrinc ! these variables are threadprivate because modify in pandf1 parallel loop. We copy them in ech thread with the copyin clause
    USE OmpCopybbb
    USE OmpCopycom
    USE OmpCopyapi
    USE omp_lib
    USE PandfTiming
    IMPLICIT NONE

    integer,intent(inout)::nnz(NchunksJac)
    integer,intent(in):: neq      ! total number of equations (all grid points)
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(in):: t           ! physical time
    real,intent(in) ::yl(*)       ! dependent variables
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    real,intent(inout) :: wk(neq)
    integer,intent(out)::iJacCol(nnzmxperchunk,NchunksJac)
    integer,intent(out):: iJacRow(neq)
    real,intent(out):: rJacElem(nnzmxperchunk,NchunksJac)
    real ::wkcopy(neq)
    real::ylcopy(neq+2)
    real::tick,tock
    external tick,tock
    integer ::iJacColCopy(nnzmxperchunk),iJacRowCopy(neq)
    integer ::ivmincopy(NchunksJac),ivmaxcopy(NchunksJac)
    integer ::NchunksJaccopy,nnzmxperchunkcopy
    real :: rJacElemCopy(nnzmxperchunk),TimeJacRowcopy(neq)
    integer:: ichunk,tid,nnzlocal
    DOUBLE PRECISION :: TimeThread
    OMPTimeCopy=tick() 
    if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,' Copying data....'
    call pandf (-1, -1, neq, 0.0, yl, ylcopy)
    if (OMPCopyArray.gt.0) then
        if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,' Copying array....'
        call OmpCopyPointerbbb
        call OmpCopyPointercom
        call OmpCopyPointerapi
    endif

    if (OMPCopyScalar.gt.0) then
        if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,' Copying scalar....'
        call OmpCopyScalarbbb
        call OmpCopyScalarcom
        call OmpCopyScalarapi
    endif

    !   We cannot use variables in the parallel construct declarations below when these variables are not in the scope of the subroutine
    NchunksJaccopy=NchunksJac
    nnzmxperchunkcopy=nnzmxperchunk
    ivmincopy(1:NchunksJac)=OMPivmin(1:NchunksJac)
    ivmaxcopy(1:NchunksJac)=OMPivmax(1:NchunksJac)
    iJacColCopy(1:nnzmxperchunk)=0
    rJacElemCopy(1:nnzmxperchunk)=0.0
    TimeJacRowcopy(1:neq)=0
    iJacRowCopy(1:neq)=0
    ylcopy(1:neq+2)=yl(1:neq+2) ! a very barbarian use of yl(neq+1) is implemented as a switch in pandf... Error-prone!
    wkcopy(1:neq)=wk(1:neq) ! Could be set equal to zero as well. The worker wk is not an output...

    if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,' Starting parallel loop'
    tid=-1
    nnzlocal=-10000
    OMPTotTimeCopy=OMPTotTimeCopy+tock(OMPTimeCopy)
    OMPTimeLocal=tick()
      !$omp parallel do schedule(dynamic,OMPLoopJacNchunk) default(shared)&
      !$omp& firstprivate(ivmincopy,ivmaxcopy,tid,nnzlocal,ylcopy)&
      !$omp& firstprivate(NchunksJaccopy,iJacRowCopy,iJacColCopy,rJacElemCopy,TimeJacRowcopy)&
      !$omp& private(TimeThread)  copyin(yinc,xlinc,xrinc)

        loopthread: do ichunk=1,NchunksJac !ichunk from 1 to Nthread, tid from 0 to Nthread-1
            Timethread = omp_get_wtime()
            tid=omp_get_thread_num()
            if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,' Thread id:',tid,' <-> ichunk:',ichunk
            ! we keep all these parameters as it is easier to debug LocalJacBuilder 
            ! and deal wichunk private/shared attributes
            call LocalJacBuilder(ivmincopy(ichunk),ivmaxcopy(ichunk),neq, t, ylcopy,yldot00,ml,mu,&
                iJacColcopy,rJacElemcopy,iJacRowcopy,ichunk,nnzlocal,nnzmxperchunk,TimeJacRowcopy)
            if (OMPJacDebug.gt.0) write(*,*) OMPJacStamp,',',tid,' nzlocal:',nnzlocal

            !!!!$omp  critical
            iJacCol(1:nnzlocal,ichunk) = iJacColCopy(1:nnzlocal)
            rJacElem(1:nnzlocal,ichunk) = rJacElemCopy(1:nnzlocal)
            iJacRow(OMPivmin(ichunk):OMPivmax(ichunk)) = iJacRowCopy(OMPivmin(ichunk):OMPivmax(ichunk))
            OMPTimeJacRow(ivmincopy(ichunk):ivmaxcopy(ichunk)) = &
                TimeJacRowcopy(ivmincopy(ichunk):ivmaxcopy(ichunk))
            nnz(ichunk)=nnzlocal
            OMPTimeLocalJac(tid+1)=omp_get_wtime() - Timethread
            !!!$omp  end critical

            if (OMPJacVerbose.gt.1) write(*,*) OMPJacStamp, &
                    ' Time in thread #', tid,':',OMPTimeLocalJac(tid+1)
            if (OMPJacVerbose.gt.1) write(*,'(a,I3,a)') 'OMP thread ',tid,' exiting...'
        enddo loopthread
      !$omp  END PARALLEL DO
    OMPTotTimeLocal=OMPTotTimeLocal+tock(OMPTimeLocal)

    if (OMPJacDebug.gt.0) then
        write(*,*) OMPJacStamp,' End of parallel loop....'
    endif
    RETURN
  END SUBROUTINE OMPJacBuilder


#endif

SUBROUTINE  OMPSplitIndex(ieqmin,ieqmax,NchunksJac,ivmin,ivmax,weight)
    IMPLICIT NONE
    integer,intent(in) ::ieqmin,ieqmax,NchunksJac
    real::weight(NchunksJac)
    integer,intent(out)::ivmin(NchunksJac),ivmax(NchunksJac)
    integer:: Nsize(NchunksJac),Msize,R,i,imax
    if (ieqmax-ieqmin+1<2) call xerrab('Number of equations to solve <2')
    if (NchunksJac.eq.ieqmax-ieqmin+1) then
        do i=1,NchunksJac
            ivmin(i)=i
            ivmax(i)=i
        enddo
        return
    endif

    if (NchunksJac.gt.1) then
        do i=1,NchunksJac
            if (weight(i)<=0) call xerrab('OMPSplitIndex: weight <0')
        enddo

        ! Normalized weights
        weight(1:NchunksJac)=weight(1:NchunksJac)/sum(weight(1:NchunksJac))*real(NchunksJac)
        do i=1,NchunksJac
            Nsize(i)=int(real((ieqmax-ieqmin+1)/NchunksJac)*weight(i))
        enddo

        do i=1,NchunksJac
            if (Nsize(i)<0) call xerrab('Nsize<0')
            if (Nsize(i)<2) Nsize(i)=Nsize(i)+1
        enddo
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) then
            imax=1
            do i=2,NchunksJac
                if (Nsize(i)>Nsize(i-1)) then
                    imax=i
                endif
            enddo
            Nsize(imax) = Nsize(imax) + ((ieqmax-ieqmin+1)-sum(Nsize))
        endif
        !write(*,*) Nsize,neq,sum(Nsize)
        if (ieqmax-ieqmin+1.ne.sum(Nsize)) call xerrab('Nsize .ne. neq!!!')
        ivmin(1)=ieqmin
        ivmax(1)=ieqmin+Nsize(1)-1
        do i=2,NchunksJac
            ivmin(i)=ivmax(i-1)+1
            ivmax(i)=ivmin(i)+Nsize(i)-1
        enddo
        if (ivmax(NchunksJac)-ivmin(1)+1.ne.(ieqmax-ieqmin+1)) &
            call xerrab('ivmax(NchunksJac)!=neq')
    else
        ivmin(NchunksJac)=ieqmin
        ivmax(NchunksJac)=ieqmax
    endif
    RETURN
END SUBROUTINE OMPSplitIndex


#ifdef _OPENMP

  SUBROUTINE OMPPandf1Rhs_old(neq,time,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,yincchunk,xincchunk, &
            iychunk,ixchunk,NchunksPandf1
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY:nx,ny
    USE Imprad, ONLY: prad
    USE Selec, ONLY:yinc,xrinc,xlinc
    USE Grid, ONLY:ijactot
    USE Cdv, ONLY: comnfe
    USE Rhsides, ONLY: psorcxg, psorrg, psordis
    IMPLICIT NONE
 
    integer yinc_bkp,xrinc_bkp,xlinc_bkp,iv,tid
    integer,intent(in)::neq
    real,intent(in)::yl(*)
    real,intent(out)::yldot(*)
    real,intent(in)::time
    real::yldotcopy(1:neq)
    real yldotsave(1:neq),ylcopy(1:neq+2), yldottot(1:neq)
    character*80 ::FileName
    real time1,time2
    integer::ichunk
    real tmp_prad(0:nx+1, 0:ny+1)
    ylcopy(1:neq+1)=yl(1:neq+1)
    yldotcopy = 0
    yldottot = 0
    tmp_prad = 0

    if (ijactot.gt.0) then

        Time1=omp_get_wtime()
!        call MakeChunksPandf1
        call OmpCopyPointerup
          !$omp parallel do default(shared) schedule(dynamic,OMPPandf1LoopNchunk) &
          !$omp& private(iv,ichunk) firstprivate(ylcopy,yldotcopy) copyin(yinc,xlinc,xrinc) &
          !$omp& REDUCTION(+:yldottot, tmp_prad)
            loopthread: do ichunk=1,NchunksPandf1 !ichunk from 1 to Nthread, tid from 0 to Nthread-1
            ! we keep all these parameters as it is easier to debug LocalJacBuilder and deal wichunk private/shared attributes
                yinc_bkp=yinc
                xlinc_bkp=xlinc
                xrinc_bkp=xrinc
                yldotcopy = 0
                if (iychunk(ichunk).ne.-1) then
                    yinc=yincchunk(ichunk)
                endif
                if (ixchunk(ichunk).ne.-1) then
                    xrinc=xincchunk(ichunk)
                    xlinc=xincchunk(ichunk)
                endif
                ! Necessary initialization for icntnunk=1 evaluation
                psorcxg = 0
                psorrg = 0
                psordis = 0

                call pandf (ixchunk(ichunk),iychunk(ichunk), neq, time, ylcopy, yldotcopy)

                do iv=1,Nivchunk(ichunk)
                    yldottot(ivchunk(ichunk,iv)) = yldottot(ivchunk(ichunk,iv)) + yldotcopy(ivchunk(ichunk,iv))
                enddo

                tmp_prad(0:nx+1, iychunk(ichunk)) =  &
                    & tmp_prad(0:nx+1, iychunk(ichunk)) + prad(0:nx+1, iychunk(ichunk))

                yinc=yinc_bkp
                xlinc=xlinc_bkp
                xrinc=xrinc_bkp
            enddo loopthread
          !$omp  END PARALLEL DO
        Time1=omp_get_wtime()-Time1

        OMPTimeParallelPandf1=Time1+OMPTimeParallelPandf1

        yldot(:neq) = yldottot
        prad = tmp_prad

        if (CheckPandf1.gt.0) then
            Time2=omp_get_wtime()
            call pandf (-1, -1, neq, time, ylcopy, yldotsave)
            Time2=omp_get_wtime()-Time2
            OMPTimeSerialPandf1=Time2+OMPTimeSerialPandf1
            if (OMPPandf1Verbose.gt.0) then
                write(*,*) "Timing Pandf1 serial:",OMPTimeSerialPandf1, &
                    "(",Time2,")/parallel:",OMPTimeParallelPandf1,'(',Time1,')'
            endif
            call Compare(yldot,yldotsave,neq)
            write(*,'(a,i4)') "  Serial and parallel pandf are identical for nfe = ", comnfe
        endif
    else
       call pandf (-1,-1, neq, time, yl, yldot)
    endif
    RETURN
  END SUBROUTINE OMPPandf1Rhs_old


  SUBROUTINE OMPinitialize_ranges2D(limits)
    Use Selec
    Use Bcond, ONLY: xcnearrb, xcnearlb
    Use Dim, ONLY: nfsp, nhsp, nisp, nx, ny
    Use Imprad, ONLY: isimpon
    IMPLICIT NONE
    integer, intent(in):: limits(4)
    integer:: xs, xe, ys, ye
    xs = limits(1)
    xe = limits(2)
    ys = limits(3)
    ye = limits(4)
      nfsp = nisp
      if (isimpon .eq. 3 .or. isimpon .eq. 4) nfsp = nhsp

    i1 = max(0, xs)
    i2 = max(1, xs)
    i2p = max(1, xs)
    i3 = xs 
    i4 = max(0, xs)
    i5 = min(nx, xe)
    i5m = min(nx-1, xe)
    i6 = min(nx+1, xe)
    i7 = xe     
    i8 = min(nx+1, xe)
    j1 = max(0, ys)
    j2 = max(1, ys)
    j1p = max(0, ys)
    j2p = max(1, ys)
    j3 = ys
    j4 = max(0, ys)
    j5 = min(ny, ye)
    j5m = min(ny-1, ye)
    j6 = min(ny+1, ye)
    j5p = min(ny, ye)
    j6p = min(ny+1, ye)
    j7 = ye
    j8 = min(ny+1, ye)



    i1 = max(0, xs-3)
    i2 = max(1, xs-2)
    i2p = max(1, xs-3)
    i3 = xs-2
    i4 = max(0, xs-2)
    i5 = min(nx, xe+2)
    i5m = min(nx-1, xe+2)
    i6 = min(nx+1, xe+3)
    i7 = xe+2
    i8 = min(nx+1, xe+2)
    j1 = max(0, ys-3)
    j2 = max(1, ys-2)
    j1p = max(0, ys-4)
    j2p = max(1, ys-3)
    j3 = ys-2
    j4 = max(0, ys-2)
    j5 = min(ny, ye+2)
    j5m = min(ny-1, ye+2)
    j6 = min(ny+1, ye+2)
    j5p = min(ny, ye+3)
    j6p = min(ny+1, ye+3)
    j7 = ye+2
    j8 = min(ny+1, ye+2)

    ixs = i2
    ixf = i5
    iys = j2
    iyf = j5
    ixs1 = i1
    ixf6 = i6
    iys1 = j1
    iyf6 = j6

    xcnearrb = .FALSE.
    xcnearlb = .FALSE.
    if (xs.eq.0) xcnearrb = .TRUE.
    if (xe.eq.nx+1) xcnearlb = .TRUE.

  END SUBROUTINE OMPinitialize_ranges2D


  SUBROUTINE OMPconvsr_vo1(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: ne, phi, ti, tg, te, nit, nz2, ng, nm, lng, ni
    USE Compla, ONLY: ne,phi,ti,tg,te,nit,nz2,ng,nm,lng,ni

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: ne_tmp(0:nx+1,0:ny+1), phi_tmp(0:nx+1,0:ny+1), ti_tmp(0:nx+1,0:ny+1), &
    &      tg_tmp(0:nx+1,0:ny+1,1:ngsp), te_tmp(0:nx+1,0:ny+1), &
    &      nit_tmp(0:nx+1,0:ny+1), nz2_tmp(0:nx+1,0:ny+1), &
    &      ng_tmp(0:nx+1,0:ny+1,1:ngsp), nm_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      lng_tmp(0:nx+1,0:ny+1,1:ngsp), ni_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    ne_tmp=0.; phi_tmp=0.; ti_tmp=0.; tg_tmp=0.; te_tmp=0.; nit_tmp=0.
    nz2_tmp=0.; ng_tmp=0.; nm_tmp=0.; lng_tmp=0.; ni_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:ne_tmp, phi_tmp, ti_tmp, tg_tmp, te_tmp, nit_tmp, nz2_tmp, ng_tmp, nm_tmp, &
    !$OMP &         lng_tmp, ni_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;phi=phi_cp;ti=ti_cp;tg=tg_cp;te=te_cp;nit=nit_cp;nz2=nz2_cp;ng=ng_cp
        nm=nm_cp;lng=lng_cp;ni=ni_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call convsr_vo1(-1, -1, ylcopy)

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        ne_tmp(xc,yc)=ne_tmp(xc,yc)+ne(xc,yc)
        phi_tmp(xc,yc)=phi_tmp(xc,yc)+phi(xc,yc)
        ti_tmp(xc,yc)=ti_tmp(xc,yc)+ti(xc,yc)
        tg_tmp(xc,yc,:)=tg_tmp(xc,yc,:)+tg(xc,yc,:)
        te_tmp(xc,yc)=te_tmp(xc,yc)+te(xc,yc)
        nit_tmp(xc,yc)=nit_tmp(xc,yc)+nit(xc,yc)
        nz2_tmp(xc,yc)=nz2_tmp(xc,yc)+nz2(xc,yc)
        ng_tmp(xc,yc,:)=ng_tmp(xc,yc,:)+ng(xc,yc,:)
        nm_tmp(xc,yc,:)=nm_tmp(xc,yc,:)+nm(xc,yc,:)
        lng_tmp(xc,yc,:)=lng_tmp(xc,yc,:)+lng(xc,yc,:)
        ni_tmp(xc,yc,:)=ni_tmp(xc,yc,:)+ni(xc,yc,:)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    ne=ne_tmp; phi=phi_tmp; ti=ti_tmp; tg=tg_tmp; te=te_tmp; nit=nit_tmp
    nz2=nz2_tmp; ng=ng_tmp; nm=nm_tmp; lng=lng_tmp; ni=ni_tmp
    call OmpCopyPointerne; call OmpCopyPointerphi; call OmpCopyPointerti
    call OmpCopyPointertg; call OmpCopyPointerte; call OmpCopyPointernit
    call OmpCopyPointernz2; call OmpCopyPointerng; call OmpCopyPointernm
    call OmpCopyPointerlng; call OmpCopyPointerni
    ne_cp=ne;phi_cp=phi;ti_cp=ti;tg_cp=tg;te_cp=te;nit_cp=nit;nz2_cp=nz2;ng_cp=ng
    nm_cp=nm;lng_cp=lng;ni_cp=ni
  END SUBROUTINE OMPconvsr_vo1


  SUBROUTINE OMPconvsr_vo2(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: up
    USE Compla, ONLY: up,nm

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: up_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    up_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:up_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;nm=nm_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call convsr_vo2(-1, -1, ylcopy)

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        up_tmp(xc,yc,:)=up_tmp(xc,yc,:)+up(xc,yc,:)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    up=up_tmp
    call OmpCopyPointerup
    up_cp=up
  END SUBROUTINE OMPconvsr_vo2


  SUBROUTINE OMPconvsr_aux1(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: pg, pr, pre, pri, tg
    USE Compla, ONLY: ne,ti,tg,te,ng,ni,pg,pr,pre,pri

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: pg_tmp(0:nx+1,0:ny+1,1:ngsp), pr_tmp(0:nx+1,0:ny+1), &
    &      pre_tmp(0:nx+1,0:ny+1), pri_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      tg_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    pg_tmp=0.; pr_tmp=0.; pre_tmp=0.; pri_tmp=0.; tg_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:pg_tmp, pr_tmp, pre_tmp, pri_tmp, tg_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;ti=ti_cp;tg=tg_cp;te=te_cp;ng=ng_cp;ni=ni_cp;pg=pg_cp;pr=pr_cp
        pre=pre_cp;pri=pri_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call convsr_aux1(-1, -1)

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        pg_tmp(xc,yc,:)=pg_tmp(xc,yc,:)+pg(xc,yc,:)
        pr_tmp(xc,yc)=pr_tmp(xc,yc)+pr(xc,yc)
        pre_tmp(xc,yc)=pre_tmp(xc,yc)+pre(xc,yc)
        pri_tmp(xc,yc,:)=pri_tmp(xc,yc,:)+pri(xc,yc,:)
        tg_tmp(xc,yc,:)=tg_tmp(xc,yc,:)+tg(xc,yc,:)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    pg=pg_tmp; pr=pr_tmp; pre=pre_tmp; pri=pri_tmp; tg=tg_tmp
    call OmpCopyPointerpg; call OmpCopyPointerpr; call OmpCopyPointerpre
    call OmpCopyPointerpri; call OmpCopyPointertg
    pg_cp=pg;pr_cp=pr;pre_cp=pre;pri_cp=pri;tg_cp=tg
  END SUBROUTINE OMPconvsr_aux1


  SUBROUTINE OMPconvsr_aux2(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: pgy0, tgy0, niy1, nity1, tiy1, phiy0, tiy0s, ney0, zeff, tiy0, pgy1, phiy0s, phiy1, &
    &    ngy0, priy0, phiv, niy0s, tey1, tiy1s, priy1, tgy1, niy0, tiv, priv, ney1, tev, &
    &    phiy1s, prtv, prev, tey0, znot, niy1s, nity0, ngy1
    USE Gradients, ONLY: gpry, gpix, ex, ey, gtex, gtiy, gpiy, gtey, gtix, gpex, gpondpotx, gprx, gpey
    USE Gradients, ONLY: gpry,gpix,ex,ey,gtex,gtiy,gpiy,gtey,gtix,gpex,gpondpotx,gprx,gpey
    USE Compla, ONLY: ne,phi,ti,tg,te,ng,ni,pg,pre,pri,pgy0,tgy0,niy1,nity1,tiy1,phiy0,tiy0s,ney0,zeff, &
    &    tiy0,pgy1,phiy0s,phiy1,ngy0,priy0,phiv,niy0s,tey1,tiy1s,priy1,tgy1,niy0,tiv,priv, &
    &    ney1,tev,phiy1s,prtv,prev,tey0,znot,niy1s,nity0,ngy1

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: pgy0_tmp(0:nx+1,0:ny+1,1:ngsp), tgy0_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      gpry_tmp(0:nx+1,0:ny+1), niy1_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      nity1_tmp(0:nx+1,0:ny+1), gpix_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      ex_tmp(0:nx+1,0:ny+1), ey_tmp(0:nx+1,0:ny+1), tiy1_tmp(0:nx+1,0:ny+1), &
    &      phiy0_tmp(0:nx+1,0:ny+1), tiy0s_tmp(0:nx+1,0:ny+1), &
    &      gtex_tmp(0:nx+1,0:ny+1), ney0_tmp(0:nx+1,0:ny+1), &
    &      gtiy_tmp(0:nx+1,0:ny+1), zeff_tmp(0:nx+1,0:ny+1), &
    &      tiy0_tmp(0:nx+1,0:ny+1), pgy1_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      phiy0s_tmp(0:nx+1,0:ny+1), phiy1_tmp(0:nx+1,0:ny+1), &
    &      ngy0_tmp(0:nx+1,0:ny+1,1:ngsp), gpiy_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      priy0_tmp(0:nx+1,0:ny+1,1:nisp), gtey_tmp(0:nx+1,0:ny+1), &
    &      gtix_tmp(0:nx+1,0:ny+1), phiv_tmp(0:nx+1,0:ny+1), &
    &      niy0s_tmp(0:nx+1,0:ny+1,1:nisp), tey1_tmp(0:nx+1,0:ny+1), &
    &      tiy1s_tmp(0:nx+1,0:ny+1), priy1_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      tgy1_tmp(0:nx+1,0:ny+1,1:ngsp), gpex_tmp(0:nx+1,0:ny+1), &
    &      niy0_tmp(0:nx+1,0:ny+1,1:nisp), tiv_tmp(0:nx+1,0:ny+1), &
    &      priv_tmp(0:nx+1,0:ny+1,1:nisp), ney1_tmp(0:nx+1,0:ny+1), &
    &      gpondpotx_tmp(0:nx+1,0:ny+1), tev_tmp(0:nx+1,0:ny+1), &
    &      phiy1s_tmp(0:nx+1,0:ny+1), prtv_tmp(0:nx+1,0:ny+1), &
    &      prev_tmp(0:nx+1,0:ny+1), tey0_tmp(0:nx+1,0:ny+1), &
    &      gprx_tmp(0:nx+1,0:ny+1), gpey_tmp(0:nx+1,0:ny+1), &
    &      znot_tmp(0:nx+1,0:ny+1), niy1s_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      nity0_tmp(0:nx+1,0:ny+1), ngy1_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    pgy0_tmp=0.; tgy0_tmp=0.; gpry_tmp=0.; niy1_tmp=0.; nity1_tmp=0.
    gpix_tmp=0.; ex_tmp=0.; ey_tmp=0.; tiy1_tmp=0.; phiy0_tmp=0.; tiy0s_tmp=0.
    gtex_tmp=0.; ney0_tmp=0.; gtiy_tmp=0.; zeff_tmp=0.; tiy0_tmp=0.
    pgy1_tmp=0.; phiy0s_tmp=0.; phiy1_tmp=0.; ngy0_tmp=0.; gpiy_tmp=0.
    priy0_tmp=0.; gtey_tmp=0.; gtix_tmp=0.; phiv_tmp=0.; niy0s_tmp=0.
    tey1_tmp=0.; tiy1s_tmp=0.; priy1_tmp=0.; tgy1_tmp=0.; gpex_tmp=0.
    niy0_tmp=0.; tiv_tmp=0.; priv_tmp=0.; ney1_tmp=0.; gpondpotx_tmp=0.
    tev_tmp=0.; phiy1s_tmp=0.; prtv_tmp=0.; prev_tmp=0.; tey0_tmp=0.
    gprx_tmp=0.; gpey_tmp=0.; znot_tmp=0.; niy1s_tmp=0.; nity0_tmp=0.
    ngy1_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:pgy0_tmp, tgy0_tmp, gpry_tmp, niy1_tmp, nity1_tmp, gpix_tmp, ex_tmp, ey_tmp, &
    !$OMP &         tiy1_tmp, phiy0_tmp, tiy0s_tmp, gtex_tmp, ney0_tmp, gtiy_tmp, zeff_tmp, &
    !$OMP &         tiy0_tmp, pgy1_tmp, phiy0s_tmp, phiy1_tmp, ngy0_tmp, gpiy_tmp, priy0_tmp, &
    !$OMP &         gtey_tmp, gtix_tmp, phiv_tmp, niy0s_tmp, tey1_tmp, tiy1s_tmp, priy1_tmp, &
    !$OMP &         tgy1_tmp, gpex_tmp, niy0_tmp, tiv_tmp, priv_tmp, ney1_tmp, gpondpotx_tmp, &
    !$OMP &         tev_tmp, phiy1s_tmp, prtv_tmp, prev_tmp, tey0_tmp, gprx_tmp, gpey_tmp, &
    !$OMP &         znot_tmp, niy1s_tmp, nity0_tmp, ngy1_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;phi=phi_cp;ti=ti_cp;tg=tg_cp;te=te_cp;ng=ng_cp;ni=ni_cp;pg=pg_cp
        pre=pre_cp;pri=pri_cp;pgy0=pgy0_cp;tgy0=tgy0_cp;gpry=gpry_cp;niy1=niy1_cp
        nity1=nity1_cp;gpix=gpix_cp;ex=ex_cp;ey=ey_cp;tiy1=tiy1_cp;phiy0=phiy0_cp
        tiy0s=tiy0s_cp;gtex=gtex_cp;ney0=ney0_cp;gtiy=gtiy_cp;zeff=zeff_cp;tiy0=tiy0_cp
        pgy1=pgy1_cp;phiy0s=phiy0s_cp;phiy1=phiy1_cp;ngy0=ngy0_cp;gpiy=gpiy_cp
        priy0=priy0_cp;gtey=gtey_cp;gtix=gtix_cp;phiv=phiv_cp;niy0s=niy0s_cp
        tey1=tey1_cp;tiy1s=tiy1s_cp;priy1=priy1_cp;tgy1=tgy1_cp;gpex=gpex_cp
        niy0=niy0_cp;tiv=tiv_cp;priv=priv_cp;ney1=ney1_cp;gpondpotx=gpondpotx_cp
        tev=tev_cp;phiy1s=phiy1s_cp;prtv=prtv_cp;prev=prev_cp;tey0=tey0_cp;gprx=gprx_cp
        gpey=gpey_cp;znot=znot_cp;niy1s=niy1s_cp;nity0=nity0_cp;ngy1=ngy1_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call convsr_aux2(-1, -1)

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        pgy0_tmp(xc,yc,:)=pgy0_tmp(xc,yc,:)+pgy0(xc,yc,:)
        tgy0_tmp(xc,yc,:)=tgy0_tmp(xc,yc,:)+tgy0(xc,yc,:)
        gpry_tmp(xc,yc)=gpry_tmp(xc,yc)+gpry(xc,yc)
        niy1_tmp(xc,yc,:)=niy1_tmp(xc,yc,:)+niy1(xc,yc,:)
        nity1_tmp(xc,yc)=nity1_tmp(xc,yc)+nity1(xc,yc)
        gpix_tmp(xc,yc,:)=gpix_tmp(xc,yc,:)+gpix(xc,yc,:)
        ex_tmp(xc,yc)=ex_tmp(xc,yc)+ex(xc,yc)
        ey_tmp(xc,yc)=ey_tmp(xc,yc)+ey(xc,yc)
        tiy1_tmp(xc,yc)=tiy1_tmp(xc,yc)+tiy1(xc,yc)
        phiy0_tmp(xc,yc)=phiy0_tmp(xc,yc)+phiy0(xc,yc)
        tiy0s_tmp(xc,yc)=tiy0s_tmp(xc,yc)+tiy0s(xc,yc)
        gtex_tmp(xc,yc)=gtex_tmp(xc,yc)+gtex(xc,yc)
        ney0_tmp(xc,yc)=ney0_tmp(xc,yc)+ney0(xc,yc)
        gtiy_tmp(xc,yc)=gtiy_tmp(xc,yc)+gtiy(xc,yc)
        zeff_tmp(xc,yc)=zeff_tmp(xc,yc)+zeff(xc,yc)
        tiy0_tmp(xc,yc)=tiy0_tmp(xc,yc)+tiy0(xc,yc)
        pgy1_tmp(xc,yc,:)=pgy1_tmp(xc,yc,:)+pgy1(xc,yc,:)
        phiy0s_tmp(xc,yc)=phiy0s_tmp(xc,yc)+phiy0s(xc,yc)
        phiy1_tmp(xc,yc)=phiy1_tmp(xc,yc)+phiy1(xc,yc)
        ngy0_tmp(xc,yc,:)=ngy0_tmp(xc,yc,:)+ngy0(xc,yc,:)
        gpiy_tmp(xc,yc,:)=gpiy_tmp(xc,yc,:)+gpiy(xc,yc,:)
        priy0_tmp(xc,yc,:)=priy0_tmp(xc,yc,:)+priy0(xc,yc,:)
        gtey_tmp(xc,yc)=gtey_tmp(xc,yc)+gtey(xc,yc)
        gtix_tmp(xc,yc)=gtix_tmp(xc,yc)+gtix(xc,yc)
        phiv_tmp(xc,yc)=phiv_tmp(xc,yc)+phiv(xc,yc)
        niy0s_tmp(xc,yc,:)=niy0s_tmp(xc,yc,:)+niy0s(xc,yc,:)
        tey1_tmp(xc,yc)=tey1_tmp(xc,yc)+tey1(xc,yc)
        tiy1s_tmp(xc,yc)=tiy1s_tmp(xc,yc)+tiy1s(xc,yc)
        priy1_tmp(xc,yc,:)=priy1_tmp(xc,yc,:)+priy1(xc,yc,:)
        tgy1_tmp(xc,yc,:)=tgy1_tmp(xc,yc,:)+tgy1(xc,yc,:)
        gpex_tmp(xc,yc)=gpex_tmp(xc,yc)+gpex(xc,yc)
        niy0_tmp(xc,yc,:)=niy0_tmp(xc,yc,:)+niy0(xc,yc,:)
        tiv_tmp(xc,yc)=tiv_tmp(xc,yc)+tiv(xc,yc)
        priv_tmp(xc,yc,:)=priv_tmp(xc,yc,:)+priv(xc,yc,:)
        ney1_tmp(xc,yc)=ney1_tmp(xc,yc)+ney1(xc,yc)
        gpondpotx_tmp(xc,yc)=gpondpotx_tmp(xc,yc)+gpondpotx(xc,yc)
        tev_tmp(xc,yc)=tev_tmp(xc,yc)+tev(xc,yc)
        phiy1s_tmp(xc,yc)=phiy1s_tmp(xc,yc)+phiy1s(xc,yc)
        prtv_tmp(xc,yc)=prtv_tmp(xc,yc)+prtv(xc,yc)
        prev_tmp(xc,yc)=prev_tmp(xc,yc)+prev(xc,yc)
        tey0_tmp(xc,yc)=tey0_tmp(xc,yc)+tey0(xc,yc)
        gprx_tmp(xc,yc)=gprx_tmp(xc,yc)+gprx(xc,yc)
        gpey_tmp(xc,yc)=gpey_tmp(xc,yc)+gpey(xc,yc)
        znot_tmp(xc,yc)=znot_tmp(xc,yc)+znot(xc,yc)
        niy1s_tmp(xc,yc,:)=niy1s_tmp(xc,yc,:)+niy1s(xc,yc,:)
        nity0_tmp(xc,yc)=nity0_tmp(xc,yc)+nity0(xc,yc)
        ngy1_tmp(xc,yc,:)=ngy1_tmp(xc,yc,:)+ngy1(xc,yc,:)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    pgy0=pgy0_tmp; tgy0=tgy0_tmp; gpry=gpry_tmp; niy1=niy1_tmp
    nity1=nity1_tmp; gpix=gpix_tmp; ex=ex_tmp; ey=ey_tmp; tiy1=tiy1_tmp
    phiy0=phiy0_tmp; tiy0s=tiy0s_tmp; gtex=gtex_tmp; ney0=ney0_tmp
    gtiy=gtiy_tmp; zeff=zeff_tmp; tiy0=tiy0_tmp; pgy1=pgy1_tmp
    phiy0s=phiy0s_tmp; phiy1=phiy1_tmp; ngy0=ngy0_tmp; gpiy=gpiy_tmp
    priy0=priy0_tmp; gtey=gtey_tmp; gtix=gtix_tmp; phiv=phiv_tmp
    niy0s=niy0s_tmp; tey1=tey1_tmp; tiy1s=tiy1s_tmp; priy1=priy1_tmp
    tgy1=tgy1_tmp; gpex=gpex_tmp; niy0=niy0_tmp; tiv=tiv_tmp; priv=priv_tmp
    ney1=ney1_tmp; gpondpotx=gpondpotx_tmp; tev=tev_tmp; phiy1s=phiy1s_tmp
    prtv=prtv_tmp; prev=prev_tmp; tey0=tey0_tmp; gprx=gprx_tmp; gpey=gpey_tmp
    znot=znot_tmp; niy1s=niy1s_tmp; nity0=nity0_tmp; ngy1=ngy1_tmp
    call OmpCopyPointerpgy0; call OmpCopyPointertgy0
    call OmpCopyPointergpry; call OmpCopyPointerniy1
    call OmpCopyPointernity1; call OmpCopyPointergpix; call OmpCopyPointerex
    call OmpCopyPointerey; call OmpCopyPointertiy1; call OmpCopyPointerphiy0
    call OmpCopyPointertiy0s; call OmpCopyPointergtex
    call OmpCopyPointerney0; call OmpCopyPointergtiy
    call OmpCopyPointerzeff; call OmpCopyPointertiy0
    call OmpCopyPointerpgy1; call OmpCopyPointerphiy0s
    call OmpCopyPointerphiy1; call OmpCopyPointerngy0
    call OmpCopyPointergpiy; call OmpCopyPointerpriy0
    call OmpCopyPointergtey; call OmpCopyPointergtix
    call OmpCopyPointerphiv; call OmpCopyPointerniy0s
    call OmpCopyPointertey1; call OmpCopyPointertiy1s
    call OmpCopyPointerpriy1; call OmpCopyPointertgy1
    call OmpCopyPointergpex; call OmpCopyPointerniy0; call OmpCopyPointertiv
    call OmpCopyPointerpriv; call OmpCopyPointerney1
    call OmpCopyPointergpondpotx; call OmpCopyPointertev
    call OmpCopyPointerphiy1s; call OmpCopyPointerprtv
    call OmpCopyPointerprev; call OmpCopyPointertey0
    call OmpCopyPointergprx; call OmpCopyPointergpey
    call OmpCopyPointerznot; call OmpCopyPointerniy1s
    call OmpCopyPointernity0; call OmpCopyPointerngy1
    pgy0_cp=pgy0;tgy0_cp=tgy0;gpry_cp=gpry;niy1_cp=niy1;nity1_cp=nity1;gpix_cp=gpix
    ex_cp=ex;ey_cp=ey;tiy1_cp=tiy1;phiy0_cp=phiy0;tiy0s_cp=tiy0s;gtex_cp=gtex
    ney0_cp=ney0;gtiy_cp=gtiy;zeff_cp=zeff;tiy0_cp=tiy0;pgy1_cp=pgy1
    phiy0s_cp=phiy0s;phiy1_cp=phiy1;ngy0_cp=ngy0;gpiy_cp=gpiy;priy0_cp=priy0
    gtey_cp=gtey;gtix_cp=gtix;phiv_cp=phiv;niy0s_cp=niy0s;tey1_cp=tey1
    tiy1s_cp=tiy1s;priy1_cp=priy1;tgy1_cp=tgy1;gpex_cp=gpex;niy0_cp=niy0;tiv_cp=tiv
    priv_cp=priv;ney1_cp=ney1;gpondpotx_cp=gpondpotx;tev_cp=tev;phiy1s_cp=phiy1s
    prtv_cp=prtv;prev_cp=prev;tey0_cp=tey0;gprx_cp=gprx;gpey_cp=gpey;znot_cp=znot
    niy1s_cp=niy1s;nity0_cp=nity0;ngy1_cp=ngy1
  END SUBROUTINE OMPconvsr_aux2


  SUBROUTINE OMPcalc_plasma_diffusivities(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Conduc, ONLY: dutm_use, difp_use, dif_use, kyi_use, trax_use, kxbohm, kxe_use, kxi_use, &
    &    kye_use, vy_use, kybohm, tray_use, dif2_use
    USE Compla, ONLY: betap
    USE OMPTiming
    USE Conduc, ONLY: dutm_use,difp_use,dif_use,kyi_use,trax_use,kxbohm,kxe_use,kxi_use,kye_use,vy_use, &
    &    kybohm,tray_use,dif2_use
    USE Compla, ONLY: te,ni,pr,betap,v2

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: dutm_use_tmp(0:nx+1,0:ny+1,1:nisp), difp_use_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      dif_use_tmp(0:nx+1,0:ny+1,1:nisp), kyi_use_tmp(0:nx+1,0:ny+1), &
    &      trax_use_tmp(0:nx+1,0:ny+1,1:nisp), kxbohm_tmp(0:nx+1,0:ny+1), &
    &      kxe_use_tmp(0:nx+1,0:ny+1), kxi_use_tmp(0:nx+1,0:ny+1), &
    &      kye_use_tmp(0:nx+1,0:ny+1), vy_use_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      betap_tmp(0:nx+1,0:ny+1), kybohm_tmp(0:nx+1,0:ny+1), &
    &      tray_use_tmp(0:nx+1,0:ny+1,1:nisp), dif2_use_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    dutm_use_tmp=0.; difp_use_tmp=0.; dif_use_tmp=0.; kyi_use_tmp=0.
    trax_use_tmp=0.; kxbohm_tmp=0.; kxe_use_tmp=0.; kxi_use_tmp=0.
    kye_use_tmp=0.; vy_use_tmp=0.; betap_tmp=0.; kybohm_tmp=0.
    tray_use_tmp=0.; dif2_use_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:dutm_use_tmp, difp_use_tmp, dif_use_tmp, kyi_use_tmp, trax_use_tmp, &
    !$OMP &         kxbohm_tmp, kxe_use_tmp, kxi_use_tmp, kye_use_tmp, vy_use_tmp, betap_tmp, &
    !$OMP &         kybohm_tmp, tray_use_tmp, dif2_use_tmp)
    DO ichunk = 1, NchunksPandf1
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_diffusivities

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        dutm_use_tmp(xc,yc,:)=dutm_use_tmp(xc,yc,:)+dutm_use(xc,yc,:)
        difp_use_tmp(xc,yc,:)=difp_use_tmp(xc,yc,:)+difp_use(xc,yc,:)
        dif_use_tmp(xc,yc,:)=dif_use_tmp(xc,yc,:)+dif_use(xc,yc,:)
        kyi_use_tmp(xc,yc)=kyi_use_tmp(xc,yc)+kyi_use(xc,yc)
        trax_use_tmp(xc,yc,:)=trax_use_tmp(xc,yc,:)+trax_use(xc,yc,:)
        kxbohm_tmp(xc,yc)=kxbohm_tmp(xc,yc)+kxbohm(xc,yc)
        kxe_use_tmp(xc,yc)=kxe_use_tmp(xc,yc)+kxe_use(xc,yc)
        kxi_use_tmp(xc,yc)=kxi_use_tmp(xc,yc)+kxi_use(xc,yc)
        kye_use_tmp(xc,yc)=kye_use_tmp(xc,yc)+kye_use(xc,yc)
        vy_use_tmp(xc,yc,:)=vy_use_tmp(xc,yc,:)+vy_use(xc,yc,:)
        betap_tmp(xc,yc)=betap_tmp(xc,yc)+betap(xc,yc)
        kybohm_tmp(xc,yc)=kybohm_tmp(xc,yc)+kybohm(xc,yc)
        tray_use_tmp(xc,yc,:)=tray_use_tmp(xc,yc,:)+tray_use(xc,yc,:)
        dif2_use_tmp(xc,yc,:)=dif2_use_tmp(xc,yc,:)+dif2_use(xc,yc,:)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    dutm_use=dutm_use_tmp; difp_use=difp_use_tmp; dif_use=dif_use_tmp
    kyi_use=kyi_use_tmp; trax_use=trax_use_tmp; kxbohm=kxbohm_tmp
    kxe_use=kxe_use_tmp; kxi_use=kxi_use_tmp; kye_use=kye_use_tmp
    vy_use=vy_use_tmp; betap=betap_tmp; kybohm=kybohm_tmp
    tray_use=tray_use_tmp; dif2_use=dif2_use_tmp
    call OmpCopyPointerdutm_use; call OmpCopyPointerdifp_use
    call OmpCopyPointerdif_use; call OmpCopyPointerkyi_use
    call OmpCopyPointertrax_use; call OmpCopyPointerkxbohm
    call OmpCopyPointerkxe_use; call OmpCopyPointerkxi_use
    call OmpCopyPointerkye_use; call OmpCopyPointervy_use
    call OmpCopyPointerbetap; call OmpCopyPointerkybohm
    call OmpCopyPointertray_use; call OmpCopyPointerdif2_use
    dutm_use_cp=dutm_use;difp_use_cp=difp_use;dif_use_cp=dif_use;kyi_use_cp=kyi_use
    trax_use_cp=trax_use;kxbohm_cp=kxbohm;kxe_use_cp=kxe_use;kxi_use_cp=kxi_use
    kye_use_cp=kye_use;vy_use_cp=vy_use;betap_cp=betap;kybohm_cp=kybohm
    tray_use_cp=tray_use;dif2_use_cp=dif2_use
  END SUBROUTINE OMPcalc_plasma_diffusivities


  SUBROUTINE OMPinitialize_driftterms(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Conduc, ONLY: eta1, dclass_e, rtaue, dclass_i
    USE UEpar, ONLY: ctaui, ctaue
    USE Compla, ONLY: loglambda
    USE UEpar, ONLY: ctaui,ctaue
    USE Conduc, ONLY: eta1,dclass_e,rtaue,dclass_i
    USE Compla, ONLY: ne,ti,te,nm,loglambda

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: eta1_tmp(0:nx+1,0:ny+1), ctaui_tmp(0:nx+1,0:ny+1,nisp), &
    &      dclass_e_tmp(0:nx+1,0:ny+1), loglambda_tmp(0:nx+1,0:ny+1), &
    &      rtaue_tmp(0:nx+1,0:ny+1), ctaue_tmp(0:nx+1,0:ny+1,nisp), &
    &      dclass_i_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    eta1_tmp=0.; ctaui_tmp=0.; dclass_e_tmp=0.; loglambda_tmp=0.; rtaue_tmp=0.
    ctaue_tmp=0.; dclass_i_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:eta1_tmp, ctaui_tmp, dclass_e_tmp, loglambda_tmp, rtaue_tmp, ctaue_tmp, &
    !$OMP &         dclass_i_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;ti=ti_cp;te=te_cp;nm=nm_cp;eta1=eta1_cp;ctaui=ctaui_cp
        dclass_e=dclass_e_cp;loglambda=loglambda_cp;rtaue=rtaue_cp;ctaue=ctaue_cp
        dclass_i=dclass_i_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call initialize_driftterms

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        eta1_tmp(xc,yc)=eta1_tmp(xc,yc)+eta1(xc,yc)
        ctaui_tmp(xc,yc,:)=ctaui_tmp(xc,yc,:)+ctaui(xc,yc,:)
        dclass_e_tmp(xc,yc)=dclass_e_tmp(xc,yc)+dclass_e(xc,yc)
        loglambda_tmp(xc,yc)=loglambda_tmp(xc,yc)+loglambda(xc,yc)
        rtaue_tmp(xc,yc)=rtaue_tmp(xc,yc)+rtaue(xc,yc)
        ctaue_tmp(xc,yc,:)=ctaue_tmp(xc,yc,:)+ctaue(xc,yc,:)
        dclass_i_tmp(xc,yc)=dclass_i_tmp(xc,yc)+dclass_i(xc,yc)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    eta1=eta1_tmp; ctaui=ctaui_tmp; dclass_e=dclass_e_tmp
    loglambda=loglambda_tmp; rtaue=rtaue_tmp; ctaue=ctaue_tmp
    dclass_i=dclass_i_tmp
    call OmpCopyPointereta1; call OmpCopyPointerctaui
    call OmpCopyPointerdclass_e; call OmpCopyPointerloglambda
    call OmpCopyPointerrtaue; call OmpCopyPointerctaue
    call OmpCopyPointerdclass_i
    eta1_cp=eta1;ctaui_cp=ctaui;dclass_e_cp=dclass_e;loglambda_cp=loglambda
    rtaue_cp=rtaue;ctaue_cp=ctaue;dclass_i_cp=dclass_i
  END SUBROUTINE OMPinitialize_driftterms


  SUBROUTINE OMPcalc_driftterms1(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: vydd, veycb, vyrd, vycf, vygp, vy, vycr, vyce, veycp, vycb, vycp
    USE Comtra, ONLY: coll_fe, diffusivwrk, coll_fi
    USE Conduc, ONLY: vyte_cft, vyti_cft, vy_cft
    USE Locflux, ONLY: cony
    USE Comtra, ONLY: coll_fe,diffusivwrk,coll_fi
    USE Conduc, ONLY: difp_use,dif_use,vy_use,eta1,rtaue,vyte_cft,vyti_cft,vy_cft
    USE Gradients, ONLY: gpry,ex,ey,gpiy,gtey,gpey
    USE Compla, ONLY: up,ne,ti,te,nit,ni,pr,niy1,tiy1,ney0,tiy0,phiv,tey1,niy0,priv,ney1,prev,tey0, &
    &    loglambda,vydd,veycb,vyrd,vycf,vygp,vy,vycr,vyce,veycp,vycb,vycp

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: vydd_tmp(0:nx+1,0:ny+1,1:nisp), veycb_tmp(0:nx+1,0:ny+1), &
    &      coll_fe_tmp(0:nx+1,0:ny+1), vyte_cft_tmp(0:nx+1,0:ny+1), &
    &      vyrd_tmp(0:nx+1,0:ny+1,1:nisp), vycf_tmp(0:nx+1,0:ny+1), &
    &      vygp_tmp(0:nx+1,0:ny+1,1:nisp), diffusivwrk_tmp(0:nx+1,0:ny+1), &
    &      vy_tmp(0:nx+1,0:ny+1,1:nisp), vycr_tmp(0:nx+1,0:ny+1), &
    &      vyce_tmp(0:nx+1,0:ny+1,1:nisp), coll_fi_tmp(0:nx+1,0:ny+1), &
    &      veycp_tmp(0:nx+1,0:ny+1), vyti_cft_tmp(0:nx+1,0:ny+1), &
    &      vycb_tmp(0:nx+1,0:ny+1,1:nisp), vycp_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      vy_cft_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    vydd_tmp=0.; veycb_tmp=0.; coll_fe_tmp=0.; vyte_cft_tmp=0.; vyrd_tmp=0.
    vycf_tmp=0.; vygp_tmp=0.; diffusivwrk_tmp=0.; vy_tmp=0.; vycr_tmp=0.
    vyce_tmp=0.; coll_fi_tmp=0.; veycp_tmp=0.; vyti_cft_tmp=0.; vycb_tmp=0.
    vycp_tmp=0.; vy_cft_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:vydd_tmp, veycb_tmp, coll_fe_tmp, vyte_cft_tmp, vyrd_tmp, vycf_tmp, vygp_tmp, &
    !$OMP &         diffusivwrk_tmp, vy_tmp, vycr_tmp, vyce_tmp, coll_fi_tmp, veycp_tmp, &
    !$OMP &         vyti_cft_tmp, vycb_tmp, vycp_tmp, vy_cft_tmp)
    DO ichunk = 1, NchunksPandf1
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_driftterms1

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        vydd_tmp(xc,yc,:)=vydd_tmp(xc,yc,:)+vydd(xc,yc,:)
        veycb_tmp(xc,yc)=veycb_tmp(xc,yc)+veycb(xc,yc)
        coll_fe_tmp(xc,yc)=coll_fe_tmp(xc,yc)+coll_fe(xc,yc)
        vyte_cft_tmp(xc,yc)=vyte_cft_tmp(xc,yc)+vyte_cft(xc,yc)
        vyrd_tmp(xc,yc,:)=vyrd_tmp(xc,yc,:)+vyrd(xc,yc,:)
        vycf_tmp(xc,yc)=vycf_tmp(xc,yc)+vycf(xc,yc)
        vygp_tmp(xc,yc,:)=vygp_tmp(xc,yc,:)+vygp(xc,yc,:)
        diffusivwrk_tmp(xc,yc)=diffusivwrk_tmp(xc,yc)+diffusivwrk(xc,yc)
        vy_tmp(xc,yc,:)=vy_tmp(xc,yc,:)+vy(xc,yc,:)
        vycr_tmp(xc,yc)=vycr_tmp(xc,yc)+vycr(xc,yc)
        vyce_tmp(xc,yc,:)=vyce_tmp(xc,yc,:)+vyce(xc,yc,:)
        coll_fi_tmp(xc,yc)=coll_fi_tmp(xc,yc)+coll_fi(xc,yc)
        veycp_tmp(xc,yc)=veycp_tmp(xc,yc)+veycp(xc,yc)
        vyti_cft_tmp(xc,yc)=vyti_cft_tmp(xc,yc)+vyti_cft(xc,yc)
        vycb_tmp(xc,yc,:)=vycb_tmp(xc,yc,:)+vycb(xc,yc,:)
        vycp_tmp(xc,yc,:)=vycp_tmp(xc,yc,:)+vycp(xc,yc,:)
        vy_cft_tmp(xc,yc,:)=vy_cft_tmp(xc,yc,:)+vy_cft(xc,yc,:)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    vydd=vydd_tmp; veycb=veycb_tmp; coll_fe=coll_fe_tmp
    vyte_cft=vyte_cft_tmp; vyrd=vyrd_tmp; vycf=vycf_tmp; vygp=vygp_tmp
    diffusivwrk=diffusivwrk_tmp; vy=vy_tmp; vycr=vycr_tmp; vyce=vyce_tmp
    coll_fi=coll_fi_tmp; veycp=veycp_tmp; vyti_cft=vyti_cft_tmp
    vycb=vycb_tmp; vycp=vycp_tmp; vy_cft=vy_cft_tmp
    call OmpCopyPointervydd; call OmpCopyPointerveycb
    call OmpCopyPointercoll_fe; call OmpCopyPointervyte_cft
    call OmpCopyPointervyrd; call OmpCopyPointervycf
    call OmpCopyPointervygp; call OmpCopyPointerdiffusivwrk
    call OmpCopyPointervy; call OmpCopyPointervycr; call OmpCopyPointervyce
    call OmpCopyPointercoll_fi; call OmpCopyPointerveycp
    call OmpCopyPointervyti_cft; call OmpCopyPointervycb
    call OmpCopyPointervycp; call OmpCopyPointervy_cft
    vydd_cp=vydd;veycb_cp=veycb;coll_fe_cp=coll_fe;vyte_cft_cp=vyte_cft;vyrd_cp=vyrd
    vycf_cp=vycf;vygp_cp=vygp;diffusivwrk_cp=diffusivwrk;vy_cp=vy;vycr_cp=vycr
    vyce_cp=vyce;coll_fi_cp=coll_fi;veycp_cp=veycp;vyti_cft_cp=vyti_cft;vycb_cp=vycb
    vycp_cp=vycp;vy_cft_cp=vy_cft
  END SUBROUTINE OMPcalc_driftterms1


  SUBROUTINE OMPcalc_driftterms2(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: v2dd, ve2cb, vytan, v2rd, ve2cd, vy, v2xgp, v2cd, vyavis, v2ce, q2cd, v2cb, v2
    USE Comflo, ONLY: fdiaxlb, fdiaxrb
    USE Xpoint_indices, ONLY: ixlb, ixrb
    USE Comflo, ONLY: fdiaxlb,fdiaxrb,fqya
    USE Conduc, ONLY: dif_use,dif2_use
    USE Gradients, ONLY: ey,gprx,gpey
    USE Compla, ONLY: up,ne,te,ni,pr,niy1,phiv,niy0,tiv,priv,tev,prev,loglambda,vy,v2dd,ve2cb,vytan, &
    &    v2rd,ve2cd,v2xgp,v2cd,vyavis,v2ce,q2cd,v2cb,v2

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii, jx
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: v2dd_tmp(0:nx+1,0:ny+1,1:nisp), ve2cb_tmp(0:nx+1,0:ny+1), &
    &      vytan_tmp(0:nx+1,0:ny+1,1:nisp), v2rd_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      fdiaxlb_tmp(0:ny+1,1:nxpt), ve2cd_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      vy_tmp(0:nx+1,0:ny+1,1:nisp), v2xgp_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      v2cd_tmp(0:nx+1,0:ny+1,1:nisp), vyavis_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      v2ce_tmp(0:nx+1,0:ny+1,1:nisp), q2cd_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      v2cb_tmp(0:nx+1,0:ny+1,1:nisp), v2_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      fdiaxrb_tmp(0:ny+1,1:nxpt)

    ! Initialize arrays to zero
    v2dd_tmp=0.; ve2cb_tmp=0.; vytan_tmp=0.; v2rd_tmp=0.; fdiaxlb_tmp=0.
    ve2cd_tmp=0.; vy_tmp=0.; v2xgp_tmp=0.; v2cd_tmp=0.; vyavis_tmp=0.
    v2ce_tmp=0.; q2cd_tmp=0.; v2cb_tmp=0.; v2_tmp=0.; fdiaxrb_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:v2dd_tmp, ve2cb_tmp, vytan_tmp, v2rd_tmp, fdiaxlb_tmp, ve2cd_tmp, vy_tmp, &
    !$OMP &         v2xgp_tmp, v2cd_tmp, vyavis_tmp, v2ce_tmp, q2cd_tmp, v2cb_tmp, v2_tmp, &
    !$OMP &         fdiaxrb_tmp)
    DO ichunk = 1, NchunksPandf1
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_driftterms2


        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)

        do jx = 1, nxpt
            if (xc .eq. ixlb(jx)) &
            &   fdiaxlb_tmp(yc,:)=fdiaxlb_tmp(yc,:)+fdiaxlb(yc,:)
            if (xc .eq. ixrb(jx)) &
            &   fdiaxrb_tmp(yc,:)=fdiaxrb_tmp(yc,:)+fdiaxrb(yc,:)
        end do

         ! Update locally calculated variables
        v2dd_tmp(xc,yc,:)=v2dd_tmp(xc,yc,:)+v2dd(xc,yc,:)
        ve2cb_tmp(xc,yc)=ve2cb_tmp(xc,yc)+ve2cb(xc,yc)
        vytan_tmp(xc,yc,:)=vytan_tmp(xc,yc,:)+vytan(xc,yc,:)
        v2rd_tmp(xc,yc,:)=v2rd_tmp(xc,yc,:)+v2rd(xc,yc,:)
        ve2cd_tmp(xc,yc,:)=ve2cd_tmp(xc,yc,:)+ve2cd(xc,yc,:)
        vy_tmp(xc,yc,:)=vy_tmp(xc,yc,:)+vy(xc,yc,:)
        v2xgp_tmp(xc,yc,:)=v2xgp_tmp(xc,yc,:)+v2xgp(xc,yc,:)
        v2cd_tmp(xc,yc,:)=v2cd_tmp(xc,yc,:)+v2cd(xc,yc,:)
        vyavis_tmp(xc,yc,:)=vyavis_tmp(xc,yc,:)+vyavis(xc,yc,:)
        v2ce_tmp(xc,yc,:)=v2ce_tmp(xc,yc,:)+v2ce(xc,yc,:)
        q2cd_tmp(xc,yc,:)=q2cd_tmp(xc,yc,:)+q2cd(xc,yc,:)
        v2cb_tmp(xc,yc,:)=v2cb_tmp(xc,yc,:)+v2cb(xc,yc,:)
        v2_tmp(xc,yc,:)=v2_tmp(xc,yc,:)+v2(xc,yc,:)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    v2dd=v2dd_tmp; ve2cb=ve2cb_tmp; vytan=vytan_tmp; v2rd=v2rd_tmp
    fdiaxlb=fdiaxlb_tmp; ve2cd=ve2cd_tmp; vy=vy_tmp; v2xgp=v2xgp_tmp
    v2cd=v2cd_tmp; vyavis=vyavis_tmp; v2ce=v2ce_tmp; q2cd=q2cd_tmp
    v2cb=v2cb_tmp; v2=v2_tmp; fdiaxrb=fdiaxrb_tmp
    call OmpCopyPointerv2dd; call OmpCopyPointerve2cb
    call OmpCopyPointervytan; call OmpCopyPointerv2rd
    call OmpCopyPointerfdiaxlb; call OmpCopyPointerve2cd
    call OmpCopyPointervy; call OmpCopyPointerv2xgp; call OmpCopyPointerv2cd
    call OmpCopyPointervyavis; call OmpCopyPointerv2ce
    call OmpCopyPointerq2cd; call OmpCopyPointerv2cb; call OmpCopyPointerv2
    call OmpCopyPointerfdiaxrb
    v2dd_cp=v2dd;ve2cb_cp=ve2cb;vytan_cp=vytan;v2rd_cp=v2rd;fdiaxlb_cp=fdiaxlb
    ve2cd_cp=ve2cd;vy_cp=vy;v2xgp_cp=v2xgp;v2cd_cp=v2cd;vyavis_cp=vyavis
    v2ce_cp=v2ce;q2cd_cp=q2cd;v2cb_cp=v2cb;v2_cp=v2;fdiaxrb_cp=fdiaxrb
  END SUBROUTINE OMPcalc_driftterms2


  SUBROUTINE OMPcalc_currents(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Comflo, ONLY: fq2d, fmity, fqym, fqyai, fqyb, fqydti, fqymi, fqydt, fqyao, fqya, fqyae, fqygp, &
    &    fq2, fqyd, fqy
    USE Comflo, ONLY: fq2d,fmity,fqym,fqyai,fqyb,fqydti,fqymi,fqydt,fqyao,fqya,fqyae,fqygp,fq2,fqyd, &
    &    fqy
    USE Conduc, ONLY: dutm_use
    USE Gradients, ONLY: ey,gpiy
    USE Compla, ONLY: ne,ti,te,niy1,tiy1,phiy0,tiy0s,ney0,zeff,tiy0,phiy0s,phiy1,niy0s,tey1,tiy1s,niy0, &
    &    ney1,phiy1s,prtv,tey0,niy1s,veycb,vycf,vy,vyce,vycb,vycp,vyavis

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: fq2d_tmp(0:nx+1,0:ny+1), fmity_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      fqym_tmp(0:nx+1,0:ny+1), fqyai_tmp(0:nx+1,0:ny+1), &
    &      fqyb_tmp(0:nx+1,0:ny+1), fqydti_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      fqymi_tmp(0:nx+1,0:ny+1,1:nisp), fqydt_tmp(0:nx+1,0:ny+1), &
    &      fqyao_tmp(0:nx+1,0:ny+1), fqya_tmp(0:nx+1,0:ny+1), &
    &      fqyae_tmp(0:nx+1,0:ny+1), fqygp_tmp(0:nx+1,0:ny+1), &
    &      fq2_tmp(0:nx+1,0:ny+1), fqyd_tmp(0:nx+1,0:ny+1), fqy_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    fq2d_tmp=0.; fmity_tmp=0.; fqym_tmp=0.; fqyai_tmp=0.; fqyb_tmp=0.
    fqydti_tmp=0.; fqymi_tmp=0.; fqydt_tmp=0.; fqyao_tmp=0.; fqya_tmp=0.
    fqyae_tmp=0.; fqygp_tmp=0.; fq2_tmp=0.; fqyd_tmp=0.; fqy_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fq2d_tmp, fmity_tmp, fqym_tmp, fqyai_tmp, fqyb_tmp, fqydti_tmp, fqymi_tmp, &
    !$OMP &         fqydt_tmp, fqyao_tmp, fqya_tmp, fqyae_tmp, fqygp_tmp, fq2_tmp, fqyd_tmp, &
    !$OMP &         fqy_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;ti=ti_cp;te=te_cp;niy1=niy1_cp;ey=ey_cp;tiy1=tiy1_cp;phiy0=phiy0_cp
        tiy0s=tiy0s_cp;ney0=ney0_cp;zeff=zeff_cp;tiy0=tiy0_cp;phiy0s=phiy0s_cp
        phiy1=phiy1_cp;gpiy=gpiy_cp;niy0s=niy0s_cp;tey1=tey1_cp;tiy1s=tiy1s_cp
        niy0=niy0_cp;ney1=ney1_cp;phiy1s=phiy1s_cp;prtv=prtv_cp;tey0=tey0_cp
        niy1s=niy1s_cp;dutm_use=dutm_use_cp;veycb=veycb_cp;vycf=vycf_cp;vy=vy_cp
        vyce=vyce_cp;vycb=vycb_cp;vycp=vycp_cp;vyavis=vyavis_cp;fq2d=fq2d_cp
        fmity=fmity_cp;fqym=fqym_cp;fqyai=fqyai_cp;fqyb=fqyb_cp;fqydti=fqydti_cp
        fqymi=fqymi_cp;fqydt=fqydt_cp;fqyao=fqyao_cp;fqya=fqya_cp;fqyae=fqyae_cp
        fqygp=fqygp_cp;fq2=fq2_cp;fqyd=fqyd_cp;fqy=fqy_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_currents

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        fq2d_tmp(xc,yc)=fq2d_tmp(xc,yc)+fq2d(xc,yc)
        fmity_tmp(xc,yc,:)=fmity_tmp(xc,yc,:)+fmity(xc,yc,:)
        fqym_tmp(xc,yc)=fqym_tmp(xc,yc)+fqym(xc,yc)
        fqyai_tmp(xc,yc)=fqyai_tmp(xc,yc)+fqyai(xc,yc)
        fqyb_tmp(xc,yc)=fqyb_tmp(xc,yc)+fqyb(xc,yc)
        fqydti_tmp(xc,yc,:)=fqydti_tmp(xc,yc,:)+fqydti(xc,yc,:)
        fqymi_tmp(xc,yc,:)=fqymi_tmp(xc,yc,:)+fqymi(xc,yc,:)
        fqydt_tmp(xc,yc)=fqydt_tmp(xc,yc)+fqydt(xc,yc)
        fqyao_tmp(xc,yc)=fqyao_tmp(xc,yc)+fqyao(xc,yc)
        fqya_tmp(xc,yc)=fqya_tmp(xc,yc)+fqya(xc,yc)
        fqyae_tmp(xc,yc)=fqyae_tmp(xc,yc)+fqyae(xc,yc)
        fqygp_tmp(xc,yc)=fqygp_tmp(xc,yc)+fqygp(xc,yc)
        fq2_tmp(xc,yc)=fq2_tmp(xc,yc)+fq2(xc,yc)
        fqyd_tmp(xc,yc)=fqyd_tmp(xc,yc)+fqyd(xc,yc)
        fqy_tmp(xc,yc)=fqy_tmp(xc,yc)+fqy(xc,yc)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    fq2d=fq2d_tmp; fmity=fmity_tmp; fqym=fqym_tmp; fqyai=fqyai_tmp
    fqyb=fqyb_tmp; fqydti=fqydti_tmp; fqymi=fqymi_tmp; fqydt=fqydt_tmp
    fqyao=fqyao_tmp; fqya=fqya_tmp; fqyae=fqyae_tmp; fqygp=fqygp_tmp
    fq2=fq2_tmp; fqyd=fqyd_tmp; fqy=fqy_tmp
    call OmpCopyPointerfq2d; call OmpCopyPointerfmity
    call OmpCopyPointerfqym; call OmpCopyPointerfqyai
    call OmpCopyPointerfqyb; call OmpCopyPointerfqydti
    call OmpCopyPointerfqymi; call OmpCopyPointerfqydt
    call OmpCopyPointerfqyao; call OmpCopyPointerfqya
    call OmpCopyPointerfqyae; call OmpCopyPointerfqygp
    call OmpCopyPointerfq2; call OmpCopyPointerfqyd; call OmpCopyPointerfqy
    fq2d_cp=fq2d;fmity_cp=fmity;fqym_cp=fqym;fqyai_cp=fqyai;fqyb_cp=fqyb
    fqydti_cp=fqydti;fqymi_cp=fqymi;fqydt_cp=fqydt;fqyao_cp=fqyao;fqya_cp=fqya
    fqyae_cp=fqyae;fqygp_cp=fqygp;fq2_cp=fq2;fqyd_cp=fqyd;fqy_cp=fqy
  END SUBROUTINE OMPcalc_currents


  SUBROUTINE OMPcalc_fqp1(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Poten, ONLY: dphi_iy1
    USE Compla, ONLY: vy, netap, vyavis
    USE Bcond, ONLY: fqpsatrb, fqpsatlb
    USE Comflo, ONLY: fqp, fqx, fqxb
    USE Xpoint_indices, ONLY: ixlb, ixrb
    USE Comflo, ONLY: fqp
    USE Compla, ONLY: ne,phi,te,pre,zeff,netap

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii, jx
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: netap_tmp(0:nx+1,0:ny+1), fqp_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    netap_tmp=0.; fqp_tmp=0.; 

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fqp_tmp,netap_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;phi=phi_cp;te=te_cp;pre=pre_cp;zeff=zeff_cp;netap=netap_cp;fqp=fqp_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_fqp1

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        fqp_tmp(xc,yc)=fqp_tmp(xc,yc)+fqp(xc,yc)
        netap_tmp(xc,yc)=netap_tmp(xc,yc)+netap(xc,yc)
            end do
    END DO
    !$OMP END PARALLEL DO
    fqp=fqp_tmp;netap=netap_tmp
    call OmpCopyPointernetap; call OmpCopyPointerfqp
    netap_cp=netap;fqp_cp=fqp
  END SUBROUTINE OMPcalc_fqp1


  SUBROUTINE OMPcalc_fqp2(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Poten, ONLY: dphi_iy1
    USE Compla, ONLY: vy, netap, vyavis
    USE Bcond, ONLY: fqpsatrb, fqpsatlb
    USE Comflo, ONLY: fqp, fqx, fqxb
    USE Xpoint_indices, ONLY: ixlb, ixrb
    USE Bcond, ONLY: fqpsatrb,fqpsatlb
    USE Poten, ONLY: dphi_iy1
    USE Comflo, ONLY: fdiaxlb,fdiaxrb,fq2,fqp,fqx,fqxb
    USE Compla, ONLY: up,ne,phi,te,ni,pre,zeff,ve2cb,v2ce,v2cb

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii, jx
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: dphi_iy1_tmp(0:nx+1), vy_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      fqpsatrb_tmp(0:ny+1,nxpt), &
    &      vyavis_tmp(0:nx+1,0:ny+1,1:nisp), fqp_tmp(0:nx+1,0:ny+1), &
    &      fqx_tmp(0:nx+1,0:ny+1), fqpsatlb_tmp(0:ny+1,nxpt), &
    &      fqxb_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    dphi_iy1_tmp=0.; vy_tmp=0.; fqpsatrb_tmp=0.; vyavis_tmp=0.
    fqp_tmp=0.; fqx_tmp=0.; fqpsatlb_tmp=0.; fqxb_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:dphi_iy1_tmp, vy_tmp, fqpsatrb_tmp, vyavis_tmp, fqp_tmp, fqx_tmp, &
    !$OMP &         fqpsatlb_tmp, fqxb_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ne=ne_cp;phi=phi_cp;te=te_cp;ni=ni_cp;pre=pre_cp;zeff=zeff_cp
        ve2cb=ve2cb_cp;fdiaxlb=fdiaxlb_cp;v2ce=v2ce_cp;v2cb=v2cb_cp;fdiaxrb=fdiaxrb_cp
        fq2=fq2_cp;fqp=fqp_cp;dphi_iy1=dphi_iy1_cp;fqpsatrb=fqpsatrb_cp;fqx=fqx_cp
        fqpsatlb=fqpsatlb_cp;fqxb=fqxb_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_fqp2

        do ii = 1, Nixychunk(ichunk)
   
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
 
        if (yc .eq. 1) &
        & dphi_iy1_tmp(xc)=dphi_iy1_tmp(xc)+dphi_iy1(xc)

        do jx = 1, nxpt
            if (xc .eq. ixlb(jx)) &
            &   fqpsatlb_tmp(yc, jx)=fqpsatlb_tmp(yc,jx)+fqpsatlb(yc,jx)
            if (xc .eq. ixrb(jx)) &
            &   fqpsatrb_tmp(yc, jx)=fqpsatrb_tmp(yc,jx)+fqpsatrb(yc,jx)
        end do


         ! Update locally calculated variables
        vy_tmp(xc,yc,:)=vy_tmp(xc,yc,:)+vy(xc,yc,:)
        vyavis_tmp(xc,yc,:)=vyavis_tmp(xc,yc,:)+vyavis(xc,yc,:)
        fqp_tmp(xc,yc)=fqp_tmp(xc,yc)+fqp(xc,yc)
        fqx_tmp(xc,yc)=fqx_tmp(xc,yc)+fqx(xc,yc)
        fqxb_tmp(xc,yc)=fqxb_tmp(xc,yc)+fqxb(xc,yc)
            end do
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    dphi_iy1=dphi_iy1_tmp; vy=vy_tmp; fqpsatrb=fqpsatrb_tmp
    vyavis=vyavis_tmp; fqp=fqp_tmp; fqx=fqx_tmp; fqpsatlb=fqpsatlb_tmp
    fqxb=fqxb_tmp
    call OmpCopyPointerdphi_iy1; call OmpCopyPointervy
    call OmpCopyPointerfqpsatrb
    call OmpCopyPointervyavis; call OmpCopyPointerfqp
    call OmpCopyPointerfqx; call OmpCopyPointerfqpsatlb
    call OmpCopyPointerfqxb
    dphi_iy1_cp=dphi_iy1;vy_cp=vy;fqpsatrb_cp=fqpsatrb;vyavis_cp=vyavis;fqp_cp=fqp
    fqx_cp=fqx;fqpsatlb_cp=fqpsatlb;fqxb_cp=fqxb
  END SUBROUTINE OMPcalc_fqp2


  SUBROUTINE OMPcalc_friction(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: upi, uz, uu, uup
    USE Cfric, ONLY: frici, frice
    USE Gradients, ONLY: ex
    USE UEpar, ONLY: cs
    USE Cfric, ONLY: frici,frice
    USE Comflo, ONLY: fqp
    USE Gradients, ONLY: ex,gtex,gpex,gpondpotx
    USE Compla, ONLY: up,ne,phi,ti,te,ni,vytan,v2cd,v2ce,v2,netap,upi,uz,uu,uup

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: upi_tmp(0:nx+1,0:ny+1,1:nisp), uz_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      uu_tmp(0:nx+1,0:ny+1,1:nisp), frici_tmp(0:nx+1,0:ny+1,nisp), &
    &      ex_tmp(0:nx+1,0:ny+1), uup_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      frice_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    upi_tmp=0.; uz_tmp=0.; uu_tmp=0.; frici_tmp=0.; ex_tmp=0.; uup_tmp=0.
    frice_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:upi_tmp, uz_tmp, uu_tmp, frici_tmp, ex_tmp, uup_tmp, frice_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ne=ne_cp;phi=phi_cp;ti=ti_cp;te=te_cp;ni=ni_cp;ex=ex_cp;gtex=gtex_cp
        gpex=gpex_cp;gpondpotx=gpondpotx_cp;vytan=vytan_cp;v2cd=v2cd_cp;v2ce=v2ce_cp
        v2=v2_cp;netap=netap_cp;fqp=fqp_cp;upi=upi_cp;uz=uz_cp;uu=uu_cp;frici=frici_cp
        uup=uup_cp;frice=frice_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_friction(-1)

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        upi_tmp(xc,yc,:)=upi_tmp(xc,yc,:)+upi(xc,yc,:)
        uz_tmp(xc,yc,:)=uz_tmp(xc,yc,:)+uz(xc,yc,:)
        uu_tmp(xc,yc,:)=uu_tmp(xc,yc,:)+uu(xc,yc,:)
        frici_tmp(xc,yc,:)=frici_tmp(xc,yc,:)+frici(xc,yc,:)
        ex_tmp(xc,yc)=ex_tmp(xc,yc)+ex(xc,yc)
        uup_tmp(xc,yc,:)=uup_tmp(xc,yc,:)+uup(xc,yc,:)
        frice_tmp(xc,yc)=frice_tmp(xc,yc)+frice(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    upi=upi_tmp; uz=uz_tmp; uu=uu_tmp; frici=frici_tmp; ex=ex_tmp; uup=uup_tmp
    frice=frice_tmp
    call OmpCopyPointerupi; call OmpCopyPointeruz; call OmpCopyPointeruu
    call OmpCopyPointerfrici; call OmpCopyPointerex; call OmpCopyPointeruup
    call OmpCopyPointerfrice
    upi_cp=upi;uz_cp=uz;uu_cp=uu;frici_cp=frici;ex_cp=ex;uup_cp=uup;frice_cp=frice
  END SUBROUTINE OMPcalc_friction


  SUBROUTINE OMPcalc_elec_velocities(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: vey, vex, upe
    USE Comflo, ONLY: fqy,fqp
    USE Gradients, ONLY: ex,ey
    USE Compla, ONLY: ne,ni,niy1,ney0,niy0,ney1,vydd,veycb,vy,vyce,ve2cb,vytan,ve2cd,v2ce,upi,vey,vex, &
    &    upe

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: vey_tmp(0:nx+1,0:ny+1), vex_tmp(0:nx+1,0:ny+1), upe_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    vey_tmp=0.; vex_tmp=0.; upe_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:vey_tmp, vex_tmp, upe_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;ni=ni_cp;niy1=niy1_cp;ex=ex_cp;ey=ey_cp;ney0=ney0_cp;niy0=niy0_cp
        ney1=ney1_cp;vydd=vydd_cp;veycb=veycb_cp;vy=vy_cp;vyce=vyce_cp;ve2cb=ve2cb_cp
        vytan=vytan_cp;ve2cd=ve2cd_cp;v2ce=v2ce_cp;fqy=fqy_cp;fqp=fqp_cp;upi=upi_cp
        vey=vey_cp;vex=vex_cp;upe=upe_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_elec_velocities

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        vey_tmp(xc,yc)=vey_tmp(xc,yc)+vey(xc,yc)
        vex_tmp(xc,yc)=vex_tmp(xc,yc)+vex(xc,yc)
        upe_tmp(xc,yc)=upe_tmp(xc,yc)+upe(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    vey=vey_tmp; vex=vex_tmp; upe=upe_tmp
    call OmpCopyPointervey; call OmpCopyPointervex; call OmpCopyPointerupe
    vey_cp=vey;vex_cp=vex;upe_cp=upe
  END SUBROUTINE OMPcalc_elec_velocities


  SUBROUTINE OMPcalc_volumetric_sources(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: psorrgc, psorgc, psordis, psorbgg, psordisg, smoc, psorc, snic, seic, psor, &
    &    psorbgz, seec, msor, psorg, psorcxg, psorrg, msorxr, psorxr, psorxrc
    USE Conduc, ONLY: nuelg, nurc, nuvl, nuiz, nucxi, nucx, nueli, nuix
    USE Compla, ONLY: rtauy, rtaux, rtau
    USE OMPTiming, ONLY: ParaTime, SerialTime
    USE Conduc, ONLY: nuelg,nurc,nuvl,nuiz,nucxi,nucx,nueli,nuix
    USE Rhsides, ONLY: psorrgc,psorgc,psordis,psorbgg,psordisg,smoc,psorc,snic,seic,psor,psorbgz,seec, &
    &    msor,psorg,psorcxg,psorrg,msorxr,psorxr,psorxrc,seik,seid,seidh
    USE Cfric, ONLY: frici,frice
    USE Comflo, ONLY: fqp
    USE Gradients, ONLY: gpix,ex,ey,gpiy,gpex,gpondpotx,gpey
    USE Compla, ONLY: up,ne,ti,te,nz2,ng,ni,vygp,vy,vycb,v2xgp,v2cb,upi,uu,vey,vex,upe,rtauy,rtaux, &
    &    rtau,vyg,uuxg

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
    real tick,tock, tsfe, tsjf, ttotfe, ttotjf, tserial, tpara
! Define local variables
    real:: psorrgc_tmp(0:nx+1,0:ny+1,1:ngsp), psorgc_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      nuelg_tmp(0:nx+1,0:ny+1,ngsp), nurc_tmp(0:nx+1,0:ny+1,ngsp), &
    &      psordis_tmp(0:nx+1,0:ny+1,1:nisp), psorbgg_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      psordisg_tmp(0:nx+1,0:ny+1,1:ngsp), smoc_tmp(0:nx+1,0:ny+1,1:nusp), &
    &      nuvl_tmp(0:nx+1,0:ny+1,nisp), psorc_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      nuiz_tmp(0:nx+1,0:ny+1,ngsp), nucxi_tmp(0:nx+1,0:ny+1,nisp), &
    &      snic_tmp(0:nx+1,0:ny+1,1:nisp), seic_tmp(0:nx+1,0:ny+1), &
    &      psor_tmp(0:nx+1,0:ny+1,1:nisp), rtauy_tmp(0:nx+1,0:ny+1), &
    &      psorbgz_tmp(0:nx+1,0:ny+1), seec_tmp(0:nx+1,0:ny+1), &
    &      msor_tmp(0:nx+1,0:ny+1,1:nisp), psorg_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      rtaux_tmp(0:nx+1,0:ny+1), psorcxg_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      psorrg_tmp(0:nx+1,0:ny+1,1:ngsp), msorxr_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      rtau_tmp(0:nx+1,0:ny+1), nucx_tmp(0:nx+1,0:ny+1,ngsp), &
    &      nueli_tmp(0:nx+1,0:ny+1,nisp), psorxr_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      psorxrc_tmp(0:nx+1,0:ny+1,1:nisp), nuix_tmp(0:nx+1,0:ny+1,ngsp)

    ! Initialize arrays to zero
    psorrgc_tmp=0.; psorgc_tmp=0.; nuelg_tmp=0.; nurc_tmp=0.; psordis_tmp=0.
    psorbgg_tmp=0.; psordisg_tmp=0.; smoc_tmp=0.; nuvl_tmp=0.; psorc_tmp=0.
    nuiz_tmp=0.; nucxi_tmp=0.; snic_tmp=0.; seic_tmp=0.; psor_tmp=0.
    rtauy_tmp=0.; psorbgz_tmp=0.; seec_tmp=0.; msor_tmp=0.; psorg_tmp=0.
    rtaux_tmp=0.; psorcxg_tmp=0.; psorrg_tmp=0.; msorxr_tmp=0.; rtau_tmp=0.
    nucx_tmp=0.; nueli_tmp=0.; psorxr_tmp=0.; psorxrc_tmp=0.; nuix_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


        tpara = tick()
    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:psorrgc_tmp, psorgc_tmp, nuelg_tmp, nurc_tmp, psordis_tmp, psorbgg_tmp, &
    !$OMP &         psordisg_tmp, smoc_tmp, nuvl_tmp, psorc_tmp, nuiz_tmp, nucxi_tmp, snic_tmp, &
    !$OMP &         seic_tmp, psor_tmp, rtauy_tmp, psorbgz_tmp, seec_tmp, msor_tmp, psorg_tmp, &
    !$OMP &         rtaux_tmp, psorcxg_tmp, psorrg_tmp, msorxr_tmp, rtau_tmp, nucx_tmp, &
    !$OMP &         nueli_tmp, psorxr_tmp, psorxrc_tmp, nuix_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ne=ne_cp;ti=ti_cp;te=te_cp;nz2=nz2_cp;ng=ng_cp;ni=ni_cp;gpix=gpix_cp
        ex=ex_cp;ey=ey_cp;gpiy=gpiy_cp;gpex=gpex_cp;gpondpotx=gpondpotx_cp;gpey=gpey_cp
        vygp=vygp_cp;vy=vy_cp;vycb=vycb_cp;v2xgp=v2xgp_cp;v2cb=v2cb_cp;fqp=fqp_cp
        upi=upi_cp;uu=uu_cp;frici=frici_cp;frice=frice_cp;vey=vey_cp;vex=vex_cp
        upe=upe_cp;psorrgc=psorrgc_cp;psorgc=psorgc_cp;nuelg=nuelg_cp;nurc=nurc_cp
        psordis=psordis_cp;psorbgg=psorbgg_cp;psordisg=psordisg_cp;smoc=smoc_cp
        nuvl=nuvl_cp;psorc=psorc_cp;nuiz=nuiz_cp;nucxi=nucxi_cp;snic=snic_cp
        seic=seic_cp;psor=psor_cp;rtauy=rtauy_cp;psorbgz=psorbgz_cp;seec=seec_cp
        msor=msor_cp;psorg=psorg_cp;rtaux=rtaux_cp;psorcxg=psorcxg_cp;psorrg=psorrg_cp
        msorxr=msorxr_cp;rtau=rtau_cp;nucx=nucx_cp;nueli=nueli_cp;psorxr=psorxr_cp
        psorxrc=psorxrc_cp;nuix=nuix_cp;vyg=vyg_cp;uuxg=uuxg_cp;seik=seik_cp
        seid=seid_cp;seidh=seidh_cp
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_volumetric_sources(1,1)

        do ii = 1, Nixychunk(ichunk)
        xc = ixychunk(ichunk,ii,1)
        yc = ixychunk(ichunk,ii,2)
        ! Update locally calculated variables
        psorrgc_tmp(xc,yc,:)=psorrgc_tmp(xc,yc,:)+psorrgc(xc,yc,:)
        psorgc_tmp(xc,yc,:)=psorgc_tmp(xc,yc,:)+psorgc(xc,yc,:)
        nuelg_tmp(xc,yc,:)=nuelg_tmp(xc,yc,:)+nuelg(xc,yc,:)
        nurc_tmp(xc,yc,:)=nurc_tmp(xc,yc,:)+nurc(xc,yc,:)
        psordis_tmp(xc,yc,:)=psordis_tmp(xc,yc,:)+psordis(xc,yc,:)
        psorbgg_tmp(xc,yc,:)=psorbgg_tmp(xc,yc,:)+psorbgg(xc,yc,:)
        psordisg_tmp(xc,yc,:)=psordisg_tmp(xc,yc,:)+psordisg(xc,yc,:)
        smoc_tmp(xc,yc,:)=smoc_tmp(xc,yc,:)+smoc(xc,yc,:)
        nuvl_tmp(xc,yc,:)=nuvl_tmp(xc,yc,:)+nuvl(xc,yc,:)
        psorc_tmp(xc,yc,:)=psorc_tmp(xc,yc,:)+psorc(xc,yc,:)
        nuiz_tmp(xc,yc,:)=nuiz_tmp(xc,yc,:)+nuiz(xc,yc,:)
        nucxi_tmp(xc,yc,:)=nucxi_tmp(xc,yc,:)+nucxi(xc,yc,:)
        snic_tmp(xc,yc,:)=snic_tmp(xc,yc,:)+snic(xc,yc,:)
        seic_tmp(xc,yc)=seic_tmp(xc,yc)+seic(xc,yc)
        psor_tmp(xc,yc,:)=psor_tmp(xc,yc,:)+psor(xc,yc,:)
        rtauy_tmp(xc,yc)=rtauy_tmp(xc,yc)+rtauy(xc,yc)
        psorbgz_tmp(xc,yc)=psorbgz_tmp(xc,yc)+psorbgz(xc,yc)
        seec_tmp(xc,yc)=seec_tmp(xc,yc)+seec(xc,yc)
        msor_tmp(xc,yc,:)=msor_tmp(xc,yc,:)+msor(xc,yc,:)
        psorg_tmp(xc,yc,:)=psorg_tmp(xc,yc,:)+psorg(xc,yc,:)
        rtaux_tmp(xc,yc)=rtaux_tmp(xc,yc)+rtaux(xc,yc)
        psorcxg_tmp(xc,yc,:)=psorcxg_tmp(xc,yc,:)+psorcxg(xc,yc,:)
        psorrg_tmp(xc,yc,:)=psorrg_tmp(xc,yc,:)+psorrg(xc,yc,:)
        msorxr_tmp(xc,yc,:)=msorxr_tmp(xc,yc,:)+msorxr(xc,yc,:)
        rtau_tmp(xc,yc)=rtau_tmp(xc,yc)+rtau(xc,yc)
        nucx_tmp(xc,yc,:)=nucx_tmp(xc,yc,:)+nucx(xc,yc,:)
        nueli_tmp(xc,yc,:)=nueli_tmp(xc,yc,:)+nueli(xc,yc,:)
        psorxr_tmp(xc,yc,:)=psorxr_tmp(xc,yc,:)+psorxr(xc,yc,:)
        psorxrc_tmp(xc,yc,:)=psorxrc_tmp(xc,yc,:)+psorxrc(xc,yc,:)
        nuix_tmp(xc,yc,:)=nuix_tmp(xc,yc,:)+nuix(xc,yc,:)
        end do
    END DO
    !$OMP  END PARALLEL DO
!        ParaTime = ParaTime + tock(tpara)

    ! Update global variables
    psorrgc=psorrgc_tmp; psorgc=psorgc_tmp; nuelg=nuelg_tmp; nurc=nurc_tmp
    psordis=psordis_tmp; psorbgg=psorbgg_tmp; psordisg=psordisg_tmp
    smoc=smoc_tmp; nuvl=nuvl_tmp; psorc=psorc_tmp; nuiz=nuiz_tmp
    nucxi=nucxi_tmp; snic=snic_tmp; seic=seic_tmp; psor=psor_tmp
    rtauy=rtauy_tmp; psorbgz=psorbgz_tmp; seec=seec_tmp; msor=msor_tmp
    psorg=psorg_tmp; rtaux=rtaux_tmp; psorcxg=psorcxg_tmp; psorrg=psorrg_tmp
    msorxr=msorxr_tmp; rtau=rtau_tmp; nucx=nucx_tmp; nueli=nueli_tmp
    psorxr=psorxr_tmp; psorxrc=psorxrc_tmp; nuix=nuix_tmp
    call OmpCopyPointerpsorrgc; call OmpCopyPointerpsorgc
    call OmpCopyPointernuelg; call OmpCopyPointernurc
    call OmpCopyPointerpsordis; call OmpCopyPointerpsorbgg
    call OmpCopyPointerpsordisg; call OmpCopyPointersmoc
    call OmpCopyPointernuvl; call OmpCopyPointerpsorc
    call OmpCopyPointernuiz; call OmpCopyPointernucxi
    call OmpCopyPointersnic; call OmpCopyPointerseic
    call OmpCopyPointerpsor; call OmpCopyPointerrtauy
    call OmpCopyPointerpsorbgz; call OmpCopyPointerseec
    call OmpCopyPointermsor; call OmpCopyPointerpsorg
    call OmpCopyPointerrtaux; call OmpCopyPointerpsorcxg
    call OmpCopyPointerpsorrg; call OmpCopyPointermsorxr
    call OmpCopyPointerrtau; call OmpCopyPointernucx
    call OmpCopyPointernueli; call OmpCopyPointerpsorxr
    call OmpCopyPointerpsorxrc; call OmpCopyPointernuix
    psorrgc_cp=psorrgc;psorgc_cp=psorgc;nuelg_cp=nuelg;nurc_cp=nurc
    psordis_cp=psordis;psorbgg_cp=psorbgg;psordisg_cp=psordisg;smoc_cp=smoc
    nuvl_cp=nuvl;psorc_cp=psorc;nuiz_cp=nuiz;nucxi_cp=nucxi;snic_cp=snic
    seic_cp=seic;psor_cp=psor;rtauy_cp=rtauy;psorbgz_cp=psorbgz;seec_cp=seec
    msor_cp=msor;psorg_cp=psorg;rtaux_cp=rtaux;psorcxg_cp=psorcxg;psorrg_cp=psorrg
    msorxr_cp=msorxr;rtau_cp=rtau;nucx_cp=nucx;nueli_cp=nueli;psorxr_cp=psorxr
    psorxrc_cp=psorxrc;nuix_cp=nuix
  END SUBROUTINE OMPcalc_volumetric_sources


  SUBROUTINE OMPneudifpg(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Locflux, ONLY: floxg, conyg, conxg, floyg
    USE Compla, ONLY: vy, vyg, uu, v2, vygtan, uug, uuxg
    USE Comflo, ONLY: fngy4ord, fngy, fngxy, fngx, fngx4ord
    USE Comflo, ONLY: fngy4ord,fngy,fngxy,fngx,fngx4ord
    USE Locflux, ONLY: floxg,conyg,conxg,floyg,floy,flox
    USE Conduc, ONLY: nuiz,nuix
    USE Compla, ONLY: up,tg,ng,pg,pgy0,pgy1,ngy0,ngy1,vy,v2,uu,vyg,vygtan,uug,uuxg

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: floxg_tmp(0:nx+1,0:ny+1), vy_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      conyg_tmp(0:nx+1,0:ny+1), vyg_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      fngy4ord_tmp(0:nx+1,0:ny+1,1:ngsp), fngy_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      fngxy_tmp(0:nx+1,0:ny+1,1:ngsp), uu_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      v2_tmp(0:nx+1,0:ny+1,1:nisp), conxg_tmp(0:nx+1,0:ny+1), &
    &      fngx_tmp(0:nx+1,0:ny+1,1:ngsp), vygtan_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      floyg_tmp(0:nx+1,0:ny+1), uug_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      fngx4ord_tmp(0:nx+1,0:ny+1,1:ngsp), uuxg_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    floxg_tmp=0.; vy_tmp=0.; conyg_tmp=0.; vyg_tmp=0.
    fngy4ord_tmp=0.; fngy_tmp=0.; fngxy_tmp=0.; uu_tmp=0.; v2_tmp=0.
    conxg_tmp=0.; fngx_tmp=0.; vygtan_tmp=0.; floyg_tmp=0.
    uug_tmp=0.; fngx4ord_tmp=0.; uuxg_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:floxg_tmp, vy_tmp, conyg_tmp, vyg_tmp, fngy4ord_tmp, &
    !$OMP &         fngy_tmp, fngxy_tmp, uu_tmp, v2_tmp, conxg_tmp, fngx_tmp, vygtan_tmp, &
    !$OMP &         floyg_tmp, uug_tmp, fngx4ord_tmp, uuxg_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;tg=tg_cp;ng=ng_cp;pg=pg_cp;pgy0=pgy0_cp;pgy1=pgy1_cp;ngy0=ngy0_cp
        ngy1=ngy1_cp;vy=vy_cp;v2=v2_cp;uu=uu_cp;nuiz=nuiz_cp;nuix=nuix_cp;floxg=floxg_cp
        conyg=conyg_cp;vyg=vyg_cp;fngy4ord=fngy4ord_cp;fngy=fngy_cp;fngxy=fngxy_cp
        conxg=conxg_cp;fngx=fngx_cp;vygtan=vygtan_cp;floyg=floyg_cp;uug=uug_cp
        fngx4ord=fngx4ord_cp;uuxg=uuxg_cp;floy=floy_cp;flox=flox_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call neudifpg

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        floxg_tmp(xc,yc)=floxg_tmp(xc,yc)+floxg(xc,yc)
        vy_tmp(xc,yc,:)=vy_tmp(xc,yc,:)+vy(xc,yc,:)
        conyg_tmp(xc,yc)=conyg_tmp(xc,yc)+conyg(xc,yc)
        vyg_tmp(xc,yc,:)=vyg_tmp(xc,yc,:)+vyg(xc,yc,:)
        fngy4ord_tmp(xc,yc,:)=fngy4ord_tmp(xc,yc,:)+fngy4ord(xc,yc,:)
        fngy_tmp(xc,yc,:)=fngy_tmp(xc,yc,:)+fngy(xc,yc,:)
        fngxy_tmp(xc,yc,:)=fngxy_tmp(xc,yc,:)+fngxy(xc,yc,:)
        uu_tmp(xc,yc,:)=uu_tmp(xc,yc,:)+uu(xc,yc,:)
        v2_tmp(xc,yc,:)=v2_tmp(xc,yc,:)+v2(xc,yc,:)
        conxg_tmp(xc,yc)=conxg_tmp(xc,yc)+conxg(xc,yc)
        fngx_tmp(xc,yc,:)=fngx_tmp(xc,yc,:)+fngx(xc,yc,:)
        vygtan_tmp(xc,yc,:)=vygtan_tmp(xc,yc,:)+vygtan(xc,yc,:)
        floyg_tmp(xc,yc)=floyg_tmp(xc,yc)+floyg(xc,yc)
        uug_tmp(xc,yc,:)=uug_tmp(xc,yc,:)+uug(xc,yc,:)
        fngx4ord_tmp(xc,yc,:)=fngx4ord_tmp(xc,yc,:)+fngx4ord(xc,yc,:)
        uuxg_tmp(xc,yc,:)=uuxg_tmp(xc,yc,:)+uuxg(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    floxg=floxg_tmp; vy=vy_tmp; conyg=conyg_tmp
    vyg=vyg_tmp; fngy4ord=fngy4ord_tmp; fngy=fngy_tmp; fngxy=fngxy_tmp
    uu=uu_tmp; v2=v2_tmp; conxg=conxg_tmp; fngx=fngx_tmp; vygtan=vygtan_tmp
    floyg=floyg_tmp; uug=uug_tmp; fngx4ord=fngx4ord_tmp
    uuxg=uuxg_tmp
    call OmpCopyPointerfloxg; call OmpCopyPointervy
    call OmpCopyPointerconyg; call OmpCopyPointervyg
    call OmpCopyPointerfngy4ord; call OmpCopyPointerfngy
    call OmpCopyPointerfngxy; call OmpCopyPointeruu; call OmpCopyPointerv2
    call OmpCopyPointerconxg; call OmpCopyPointerfngx
    call OmpCopyPointervygtan; 
    call OmpCopyPointerfloyg; call OmpCopyPointeruug
    call OmpCopyPointerfngx4ord; call OmpCopyPointeruuxg
    floxg_cp=floxg;vy_cp=vy;conyg_cp=conyg;vyg_cp=vyg;fngy4ord_cp=fngy4ord
    fngy_cp=fngy;fngxy_cp=fngxy;uu_cp=uu;v2_cp=v2;conxg_cp=conxg;fngx_cp=fngx
    vygtan_cp=vygtan;floyg_cp=floyg;uug_cp=uug;fngx4ord_cp=fngx4ord;uuxg_cp=uuxg
  END SUBROUTINE OMPneudifpg


  SUBROUTINE OMPcalc_srcmod(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: smoc, seic, seec, snic
    USE Rhsides, ONLY: smoc,snic,seic,seec
    USE Cfric, ONLY: frici,frice
    USE Comflo, ONLY: fqp
    USE Gradients, ONLY: gpix,ex,ey,gpiy,gpex,gpondpotx,gpey
    USE Compla, ONLY: up,ne,nz2,ni,vygp,vy,vycb,v2xgp,v2cb,upi,vey,vex,upe

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: smoc_tmp(0:nx+1,0:ny+1,1:nusp), seic_tmp(0:nx+1,0:ny+1), &
    &      seec_tmp(0:nx+1,0:ny+1), snic_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    smoc_tmp=0.; seic_tmp=0.; seec_tmp=0.; snic_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:smoc_tmp, seic_tmp, seec_tmp, snic_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ne=ne_cp;nz2=nz2_cp;ni=ni_cp;gpix=gpix_cp;ex=ex_cp;ey=ey_cp
        gpiy=gpiy_cp;gpex=gpex_cp;gpondpotx=gpondpotx_cp;gpey=gpey_cp;vygp=vygp_cp
        vy=vy_cp;vycb=vycb_cp;v2xgp=v2xgp_cp;v2cb=v2cb_cp;fqp=fqp_cp;upi=upi_cp
        frici=frici_cp;frice=frice_cp;vey=vey_cp;vex=vex_cp;upe=upe_cp;smoc=smoc_cp
        snic=snic_cp;seic=seic_cp;seec=seec_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_srcmod

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        smoc_tmp(xc,yc,:)=smoc_tmp(xc,yc,:)+smoc(xc,yc,:)
        seic_tmp(xc,yc)=seic_tmp(xc,yc)+seic(xc,yc)
        seec_tmp(xc,yc)=seec_tmp(xc,yc)+seec(xc,yc)
        snic_tmp(xc,yc,:)=snic_tmp(xc,yc,:)+snic(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    smoc=smoc_tmp; seic=seic_tmp; seec=seec_tmp; snic=snic_tmp
    call OmpCopyPointersmoc; call OmpCopyPointerseic
    call OmpCopyPointerseec; call OmpCopyPointersnic
    smoc_cp=smoc;seic_cp=seic;seec_cp=seec;snic_cp=snic
  END SUBROUTINE OMPcalc_srcmod


  SUBROUTINE OMPcalc_plasma_viscosities(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE UEpar, ONLY: ctaui
    USE Conduc, ONLY: alfneo, visy, ktneo, visxneo, nuiistar, nuii, visx, k2neo
    USE Wkspace, ONLY: w
    USE Wkspace, ONLY: w
    USE UEpar, ONLY: ctaui
    USE Conduc, ONLY: trax_use,tray_use,eta1,alfneo,visy,ktneo,visxneo,nuiistar,nuii,visx,k2neo
    USE Compla, ONLY: up,ti,tg,ng,nm,ni,loglambda,upi

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: alfneo_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      visy_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      ctaui_tmp(0:nx+1,0:ny+1,nisp), ktneo_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      visxneo_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      nuiistar_tmp(0:nx+1,0:ny+1,1:nisp), nuii_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      w_tmp(0:nx+1,0:ny+1), visx_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      k2neo_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    alfneo_tmp=0.; visy_tmp=0.
    ctaui_tmp=0.; ktneo_tmp=0.; visxneo_tmp=0.; nuiistar_tmp=0.
    nuii_tmp=0.; w_tmp=0.; visx_tmp=0.; k2neo_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:alfneo_tmp, visy_tmp, ctaui_tmp, &
    !$OMP &         ktneo_tmp, visxneo_tmp, nuiistar_tmp, nuii_tmp, w_tmp, visx_tmp, &
    !$OMP &         k2neo_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ti=ti_cp;tg=tg_cp;ng=ng_cp;nm=nm_cp;ni=ni_cp;trax_use=trax_use_cp
        tray_use=tray_use_cp;eta1=eta1_cp;ctaui=ctaui_cp;loglambda=loglambda_cp
        upi=upi_cp;alfneo=alfneo_cp;visy=visy_cp;ktneo=ktneo_cp;visxneo=visxneo_cp
        nuiistar=nuiistar_cp;nuii=nuii_cp;w=w_cp;visx=visx_cp;k2neo=k2neo_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_viscosities

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        alfneo_tmp(xc,yc,:)=alfneo_tmp(xc,yc,:)+alfneo(xc,yc,:)
        visy_tmp(xc,yc,:)=visy_tmp(xc,yc,:)+visy(xc,yc,:)
        ctaui_tmp(xc,yc,:)=ctaui_tmp(xc,yc,:)+ctaui(xc,yc,:)
        ktneo_tmp(xc,yc,:)=ktneo_tmp(xc,yc,:)+ktneo(xc,yc,:)
        visxneo_tmp(xc,yc,:)=visxneo_tmp(xc,yc,:)+visxneo(xc,yc,:)
        nuiistar_tmp(xc,yc,:)=nuiistar_tmp(xc,yc,:)+nuiistar(xc,yc,:)
        nuii_tmp(xc,yc,:)=nuii_tmp(xc,yc,:)+nuii(xc,yc,:)
        w_tmp(xc,yc)=w_tmp(xc,yc)+w(xc,yc)
        visx_tmp(xc,yc,:)=visx_tmp(xc,yc,:)+visx(xc,yc,:)
        k2neo_tmp(xc,yc,:)=k2neo_tmp(xc,yc,:)+k2neo(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    alfneo=alfneo_tmp; 
    visy=visy_tmp; ctaui=ctaui_tmp; ktneo=ktneo_tmp; visxneo=visxneo_tmp
    nuiistar=nuiistar_tmp; nuii=nuii_tmp; w=w_tmp; visx=visx_tmp
    k2neo=k2neo_tmp
    call OmpCopyPointeralfneo; 
    call OmpCopyPointervisy
    call OmpCopyPointerctaui; call OmpCopyPointerktneo
    call OmpCopyPointervisxneo; 
    call OmpCopyPointernuiistar; call OmpCopyPointernuii
    call OmpCopyPointerw; call OmpCopyPointervisx; call OmpCopyPointerk2neo
    alfneo_cp=alfneo;visy_cp=visy;ctaui_cp=ctaui;ktneo_cp=ktneo;visxneo_cp=visxneo
    nuiistar_cp=nuiistar;nuii_cp=nuii;w_cp=w;visx_cp=visx;k2neo_cp=k2neo
  END SUBROUTINE OMPcalc_plasma_viscosities


  SUBROUTINE OMPcalc_plasma_heatconductivities(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE UEpar, ONLY: ctaui, ctaue
    USE Conduc, ONLY: hcxineo, hcyij, hcxij, hcyi, hcyn, hcye, hcxi, hcxe, hcxn
    USE Wkspace, ONLY: w2, w1
    USE Comflo, ONLY: qipar
    USE Comflo, ONLY: qipar
    USE Wkspace, ONLY: w2,w1
    USE Comtra, ONLY: diffusivwrk
    USE UEpar, ONLY: ctaui,ctaue
    USE Conduc, ONLY: kyi_use,kxe_use,kxi_use,kye_use,dclass_e,dclass_i,nucx,k2neo,hcxineo,hcyij,hcxij, &
    &    hcyi,hcyn,hcye,hcxi,hcxe,hcxn
    USE Compla, ONLY: ne,ti,tg,te,ng,ni,tgy0,niy1,tiy1,zeff,tiy0,ngy0,tgy1,niy0,ngy1,loglambda

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: hcxineo_tmp(0:nx+1,0:ny+1), &
    &      ctaui_tmp(0:nx+1,0:ny+1,nisp), w2_tmp(0:nx+1,0:ny+1), &
    &      hcyij_tmp(0:nx+1,0:ny+1,1:nisp), hcxij_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      hcyi_tmp(0:nx+1,0:ny+1), ctaue_tmp(0:nx+1,0:ny+1,nisp), &
    &      hcyn_tmp(0:nx+1,0:ny+1), hcye_tmp(0:nx+1,0:ny+1), &
    &      qipar_tmp(0:nx+1,0:ny+1,nisp), hcxi_tmp(0:nx+1,0:ny+1), &
    &      w1_tmp(0:nx+1,0:ny+1), hcxe_tmp(0:nx+1,0:ny+1), hcxn_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    hcxineo_tmp=0.; ctaui_tmp=0.; w2_tmp=0.
    hcyij_tmp=0.; hcxij_tmp=0.; hcyi_tmp=0.; ctaue_tmp=0.; hcyn_tmp=0.
    hcye_tmp=0.; qipar_tmp=0.; hcxi_tmp=0.; w1_tmp=0.
    hcxe_tmp=0.; hcxn_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:hcxineo_tmp, ctaui_tmp, w2_tmp, hcyij_tmp, hcxij_tmp, &
    !$OMP &         hcyi_tmp, ctaue_tmp, hcyn_tmp, hcye_tmp, qipar_tmp, hcxi_tmp, &
    !$OMP &         w1_tmp, hcxe_tmp, hcxn_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;ti=ti_cp;tg=tg_cp;te=te_cp;ng=ng_cp;ni=ni_cp;tgy0=tgy0_cp;niy1=niy1_cp
        tiy1=tiy1_cp;zeff=zeff_cp;tiy0=tiy0_cp;ngy0=ngy0_cp;tgy1=tgy1_cp;niy0=niy0_cp
        ngy1=ngy1_cp;kyi_use=kyi_use_cp;kxe_use=kxe_use_cp;kxi_use=kxi_use_cp
        kye_use=kye_use_cp;ctaui=ctaui_cp;dclass_e=dclass_e_cp;loglambda=loglambda_cp
        ctaue=ctaue_cp;dclass_i=dclass_i_cp;diffusivwrk=diffusivwrk_cp;nucx=nucx_cp
        k2neo=k2neo_cp;hcxineo=hcxineo_cp;w2=w2_cp;hcyij=hcyij_cp;hcxij=hcxij_cp
        hcyi=hcyi_cp;hcyn=hcyn_cp;hcye=hcye_cp;qipar=qipar_cp;hcxi=hcxi_cp;w1=w1_cp
        hcxe=hcxe_cp;hcxn=hcxn_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_heatconductivities

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        hcxineo_tmp(xc,yc)=hcxineo_tmp(xc,yc)+hcxineo(xc,yc)
        ctaui_tmp(xc,yc,:)=ctaui_tmp(xc,yc,:)+ctaui(xc,yc,:)
        w2_tmp(xc,yc)=w2_tmp(xc,yc)+w2(xc,yc)
        hcyij_tmp(xc,yc,:)=hcyij_tmp(xc,yc,:)+hcyij(xc,yc,:)
        hcxij_tmp(xc,yc,:)=hcxij_tmp(xc,yc,:)+hcxij(xc,yc,:)
        hcyi_tmp(xc,yc)=hcyi_tmp(xc,yc)+hcyi(xc,yc)
        ctaue_tmp(xc,yc,:)=ctaue_tmp(xc,yc,:)+ctaue(xc,yc,:)
        hcyn_tmp(xc,yc)=hcyn_tmp(xc,yc)+hcyn(xc,yc)
        hcye_tmp(xc,yc)=hcye_tmp(xc,yc)+hcye(xc,yc)
        qipar_tmp(xc,yc,:)=qipar_tmp(xc,yc,:)+qipar(xc,yc,:)
        hcxi_tmp(xc,yc)=hcxi_tmp(xc,yc)+hcxi(xc,yc)
        w1_tmp(xc,yc)=w1_tmp(xc,yc)+w1(xc,yc)
        hcxe_tmp(xc,yc)=hcxe_tmp(xc,yc)+hcxe(xc,yc)
        hcxn_tmp(xc,yc)=hcxn_tmp(xc,yc)+hcxn(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    hcxineo=hcxineo_tmp; ctaui=ctaui_tmp; w2=w2_tmp
    hcyij=hcyij_tmp; hcxij=hcxij_tmp; hcyi=hcyi_tmp; ctaue=ctaue_tmp
    hcyn=hcyn_tmp; hcye=hcye_tmp; qipar=qipar_tmp; 
    hcxi=hcxi_tmp; w1=w1_tmp; hcxe=hcxe_tmp; hcxn=hcxn_tmp
    call OmpCopyPointerhcxineo
    call OmpCopyPointerctaui; call OmpCopyPointerw2
    call OmpCopyPointerhcyij; call OmpCopyPointerhcxij
    call OmpCopyPointerhcyi; call OmpCopyPointerctaue
    call OmpCopyPointerhcyn; call OmpCopyPointerhcye
    call OmpCopyPointerqipar; 
    call OmpCopyPointerhcxi; call OmpCopyPointerw1; call OmpCopyPointerhcxe
    call OmpCopyPointerhcxn
    hcxineo_cp=hcxineo;ctaui_cp=ctaui;w2_cp=w2;hcyij_cp=hcyij;hcxij_cp=hcxij
    hcyi_cp=hcyi;ctaue_cp=ctaue;hcyn_cp=hcyn;hcye_cp=hcye;qipar_cp=qipar
    hcxi_cp=hcxi;w1_cp=w1;hcxe_cp=hcxe;hcxn_cp=hcxn
  END SUBROUTINE OMPcalc_plasma_heatconductivities


  SUBROUTINE OMPcalc_plasma_equipartition(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Conduc, ONLY: eqp
    USE Conduc, ONLY: eqp
    USE Compla, ONLY: ne,ti,te,ni,loglambda

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: eqp_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    eqp_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:eqp_tmp)
    DO ichunk = 1, NchunksPandf1
        ne=ne_cp;ti=ti_cp;te=te_cp;ni=ni_cp;loglambda=loglambda_cp;eqp=eqp_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_equipartition

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        eqp_tmp(xc,yc)=eqp_tmp(xc,yc)+eqp(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO
    eqp=eqp_tmp

    ! Update global variables
    call OmpCopyPointereqp
    eqp_cp=eqp
  END SUBROUTINE OMPcalc_plasma_equipartition


  SUBROUTINE OMPcalc_gas_heatconductivities(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Conduc, ONLY: hcxg, hcyg
    USE Conduc, ONLY: hcyn,hcxn,hcxg,hcyg
    USE Compla, ONLY: tg,ng,ni,tgy0,niy1,ngy0,tgy1,niy0,ngy1

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: hcxg_tmp(0:nx+1,0:ny+1,1:ngsp), hcyg_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    hcxg_tmp=0.; hcyg_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:hcxg_tmp, hcyg_tmp)
    DO ichunk = 1, NchunksPandf1
        tg=tg_cp;ng=ng_cp;ni=ni_cp;tgy0=tgy0_cp;niy1=niy1_cp;ngy0=ngy0_cp;tgy1=tgy1_cp
        niy0=niy0_cp;ngy1=ngy1_cp;hcyn=hcyn_cp;hcxn=hcxn_cp;hcxg=hcxg_cp;hcyg=hcyg_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_gas_heatconductivities

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        hcxg_tmp(xc,yc,:)=hcxg_tmp(xc,yc,:)+hcxg(xc,yc,:)
        hcyg_tmp(xc,yc,:)=hcyg_tmp(xc,yc,:)+hcyg(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    hcxg=hcxg_tmp; hcyg=hcyg_tmp
    call OmpCopyPointerhcxg; call OmpCopyPointerhcyg
    hcxg_cp=hcxg;hcyg_cp=hcyg
  END SUBROUTINE OMPcalc_gas_heatconductivities


  SUBROUTINE OMPengbalg(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Locflux, ONLY: floxge, floyge, conxge, conyge
    USE Comflo, ONLY: fegx, fegy, fegxy
    USE Rhsides, ONLY: segc
    USE Locflux, ONLY: floxge,floyge,conxge,conyge
    USE Comflo, ONLY: fngy,fngx,fegx,fegxy,fegy
    USE Rhsides, ONLY: seic,segc
    USE Conduc, ONLY: nuiz,nucxi,nueli,nuix,hcxg,hcyg
    USE Compla, ONLY: ti,tg,ng,ni,pg,pgy0,pgy1,vyg,uuxg

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: floxge_tmp(0:nx+1,0:ny+1,1:ngsp), segc_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      floyge_tmp(0:nx+1,0:ny+1,1:ngsp), conxge_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      conyge_tmp(0:nx+1,0:ny+1,1:ngsp), fegx_tmp(0:nx+1,0:ny+1,1:ngsp), &
    &      fegxy_tmp(0:nx+1,0:ny+1,1:ngsp), fegy_tmp(0:nx+1,0:ny+1,1:ngsp)



    ! Initialize arrays to zero
    floxge_tmp=0.; segc_tmp=0.; floyge_tmp=0.; conxge_tmp=0.; conyge_tmp=0.
    fegx_tmp=0;fegxy_tmp=0;fegy_tmp=0

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:floxge_tmp, segc_tmp, floyge_tmp, conxge_tmp, & 
    !$OMP &      conyge_tmp, fegx_tmp, fegxy_tmp, fegy_tmp)
    DO ichunk = 1, NchunksPandf1
        ti=ti_cp;tg=tg_cp;ng=ng_cp;ni=ni_cp;pg=pg_cp;pgy0=pgy0_cp;pgy1=pgy1_cp
        nuiz=nuiz_cp;nucxi=nucxi_cp;seic=seic_cp;nueli=nueli_cp;nuix=nuix_cp;vyg=vyg_cp
        fngy=fngy_cp;fngx=fngx_cp;uuxg=uuxg_cp;hcxg=hcxg_cp;hcyg=hcyg_cp
        floxge=floxge_cp;segc=segc_cp;floyge=floyge_cp;conxge=conxge_cp;conyge=conyge_cp
        fegx=fegx_cp;fegxy=fegxy_cp;fegy=fegy_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call engbalg

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        floxge_tmp(xc,yc,:)=floxge_tmp(xc,yc,:)+floxge(xc,yc,:)
        segc_tmp(xc,yc,:)=segc_tmp(xc,yc,:)+segc(xc,yc,:)
        fegx_tmp(xc,yc,:)=fegx_tmp(xc,yc,:)+fegx(xc,yc,:)
        fegxy_tmp(xc,yc,:)=fegxy_tmp(xc,yc,:)+fegxy(xc,yc,:)
        fegy_tmp(xc,yc,:)=fegy_tmp(xc,yc,:)+fegy(xc,yc,:)
        floyge_tmp(xc,yc,:)=floyge_tmp(xc,yc,:)+floyge(xc,yc,:)
        conxge_tmp(xc,yc,:)=conxge_tmp(xc,yc,:)+conxge(xc,yc,:)
        conyge_tmp(xc,yc,:)=conyge_tmp(xc,yc,:)+conyge(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    floxge=floxge_tmp; segc=segc_tmp; floyge=floyge_tmp; conxge=conxge_tmp
    conyge=conyge_tmp; fegx=fegx_tmp; fegy=fegy_tmp;fegxy=fegxy_tmp
    call OmpCopyPointerfloxge; call OmpCopyPointersegc
    call OmpCopyPointerfloyge; call OmpCopyPointerconxge
    call OmpCopyPointerconyge; call OmpCopyPointerfegx
    call OmpCopyPointerfegxy; call OmpCopyPointerfegy
    floxge_cp=floxge;segc_cp=segc;floyge_cp=floyge;conxge_cp=conxge;conyge_cp=conyge
    fegx_cp=fegx;fegxy_cp=fegxy;fegy_cp=fegy
  END SUBROUTINE OMPengbalg


  SUBROUTINE OMPcalc_plasma_transport(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Comflo, ONLY: fniy, fniy4ord, fnixcb, fniycbo, fniycb, fnix
    USE Comflo, ONLY: fngy,fngx,fniy,fniy4ord,fnixcb,fniycb,fnix
    USE Compla, ONLY: ni,niy1,niy0,vy,vycb,vytan,v2cb,uu

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: fniy_tmp(0:nx+1,0:ny+1,1:nisp), fniy4ord_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      fnixcb_tmp(0:nx+1,0:ny+1,1:nisp), fniycbo_tmp(0:nx+1,1:nisp), &
    &      fniycb_tmp(0:nx+1,0:ny+1,1:nisp), fnix_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    fniy_tmp=0.; fniy4ord_tmp=0.; fnixcb_tmp=0.; fniycbo_tmp=0.; fniycb_tmp=0.
    fnix_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fniy_tmp, fniy4ord_tmp, fnixcb_tmp, fniycbo_tmp, fniycb_tmp, fnix_tmp)
    DO ichunk = 1, NchunksPandf1
        ni=ni_cp;niy1=niy1_cp;niy0=niy0_cp;vy=vy_cp;vycb=vycb_cp;vytan=vytan_cp
        v2cb=v2cb_cp;uu=uu_cp;fngy=fngy_cp;fngx=fngx_cp;fniy=fniy_cp
        fniy4ord=fniy4ord_cp;fnixcb=fnixcb_cp;fniycb=fniycb_cp;fnix=fnix_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_transport

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        fniy_tmp(xc,yc,:)=fniy_tmp(xc,yc,:)+fniy(xc,yc,:)
        fniy4ord_tmp(xc,yc,:)=fniy4ord_tmp(xc,yc,:)+fniy4ord(xc,yc,:)
        fnixcb_tmp(xc,yc,:)=fnixcb_tmp(xc,yc,:)+fnixcb(xc,yc,:)
        fniycbo_tmp(xc,yc)=fniycbo_tmp(xc,yc)+fniycbo(xc,yc)
        fniycb_tmp(xc,yc,:)=fniycb_tmp(xc,yc,:)+fniycb(xc,yc,:)
        fnix_tmp(xc,yc,:)=fnix_tmp(xc,yc,:)+fnix(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    fniy=fniy_tmp; fniy4ord=fniy4ord_tmp; fnixcb=fnixcb_tmp
    fniycbo=fniycbo_tmp; fniycb=fniycb_tmp; fnix=fnix_tmp
    call OmpCopyPointerfniy; call OmpCopyPointerfniy4ord
    call OmpCopyPointerfnixcb; call OmpCopyPointerfniycbo
    call OmpCopyPointerfniycb; call OmpCopyPointerfnix
    fniy_cp=fniy;fniy4ord_cp=fniy4ord;fnixcb_cp=fnixcb;fniycbo_cp=fniycbo
    fniycb_cp=fniycb;fnix_cp=fnix
  END SUBROUTINE OMPcalc_plasma_transport


  SUBROUTINE OMPcalc_plasma_momentum_coeffs(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Locflux, ONLY: conx, floy, flox, cony
    USE Locflux, ONLY: conx,floy,flox,cony
    USE Conduc, ONLY: visy,visx
    USE Compla, ONLY: up,nm,vy,vytan,v2,uu

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii, iusp, ixpt
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: conx_tmp(0:nx+1,0:ny+1,1:nusp), floy_tmp(0:nx+1,0:ny+1,1:nusp), &
    &      flox_tmp(0:nx+1,0:ny+1,1:nusp), cony_tmp(0:nx+1,0:ny+1,1:nusp)

    ! Initialize arrays to zero
    conx_tmp=0.; floy_tmp=0.; flox_tmp=0.; cony_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:conx_tmp, floy_tmp, flox_tmp, cony_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;nm=nm_cp;vy=vy_cp;vytan=vytan_cp;v2=v2_cp;uu=uu_cp;visy=visy_cp
        visx=visx_cp;conx=conx_cp;floy=floy_cp;flox=flox_cp;cony=cony_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_momentum_coeffs
        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        conx_tmp(xc,yc,:)=conx_tmp(xc,yc,:)+conx(xc,yc,:)
        floy_tmp(xc,yc,:)=floy_tmp(xc,yc,:)+floy(xc,yc,:)
        flox_tmp(xc,yc,:)=flox_tmp(xc,yc,:)+flox(xc,yc,:)
        cony_tmp(xc,yc,:)=cony_tmp(xc,yc,:)+cony(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    conx=conx_tmp; floy=floy_tmp;flox=flox_tmp;cony=cony_tmp
    call OmpCopyPointerconx; call OmpCopyPointerfloy
    call OmpCopyPointerflox; call OmpCopyPointercony
    conx_cp=conx;floy_cp=floy;flox_cp=flox;cony_cp=cony
  END SUBROUTINE OMPcalc_plasma_momentum_coeffs


  SUBROUTINE OMPcalc_plasma_momentum(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: fmivxpt, fmihxpt, vyvxpt, nixpt, vyhxpt, visyxpt, upxpt, up
    USE Comflo, ONLY: fmixy, fmix, fmiy
    USE Rhsides, ONLY: smoc
    USE UEpar, ONLY: methu, isupon
    USE Locflux, ONLY: flox, floy, conx, cony
    USE Comflo, ONLY: fmixy,fmix,fmiy
    USE Conduc, ONLY: visy,visx
    USE Rhsides, ONLY: smoc
    USE Compla, ONLY: up,ti,tg,nm,ni,vy,fmivxpt,fmihxpt,vyvxpt,nixpt,vyhxpt,visyxpt,upxpt

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii, iusp, ixpt, ifld
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: fmivxpt_tmp(1:nusp,1:nxpt), fmihxpt_tmp(1:nusp,1:nxpt), &
    &      fmixy_tmp(0:nx+1,0:ny+1,1:nusp), &
    &      vyvxpt_tmp(1:nusp,1:nxpt), nixpt_tmp(1:nusp,1:nxpt), &
    &      vyhxpt_tmp(1:nusp,1:nxpt), visyxpt_tmp(1:nusp,1:nxpt), &
    &      upxpt_tmp(1:nusp,1:nxpt), smoc_tmp(0:nx+1,0:ny+1,1:nusp), &
    &      fmix_tmp(0:nx+1,0:ny+1,1:nusp), &
    &      fmiy_tmp(0:nx+1,0:ny+1,1:nusp)

    ! Initialize arrays to zero
    fmivxpt_tmp=0.; fmihxpt_tmp=0.
    fmixy_tmp=0.; vyvxpt_tmp=0.; nixpt_tmp=0.; vyhxpt_tmp=0.
    visyxpt_tmp=0.; upxpt_tmp=0.; smoc_tmp=0.
    fmix_tmp=0.;fmiy_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fmivxpt_tmp, fmihxpt_tmp, fmixy_tmp, &
    !$OMP &         vyvxpt_tmp, nixpt_tmp, vyhxpt_tmp, visyxpt_tmp, upxpt_tmp, smoc_tmp, &
    !$OMP &         fmix_tmp, fmiy_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ti=ti_cp;tg=tg_cp;nm=nm_cp;ni=ni_cp;vy=vy_cp;smoc=smoc_cp;visy=visy_cp
        visx=visx_cp;fmivxpt=fmivxpt_cp;fmihxpt=fmihxpt_cp;fmixy=fmixy_cp
        vyvxpt=vyvxpt_cp;nixpt=nixpt_cp;vyhxpt=vyhxpt_cp;visyxpt=visyxpt_cp
        upxpt=upxpt_cp;fmix=fmix_cp;fmiy=fmiy_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_momentum(xc, yc)

        do iusp = 1, nusp
            do ixpt = 1, nxpt
                fmivxpt_tmp(iusp,ixpt)=fmivxpt_tmp(iusp,ixpt)+fmivxpt(iusp,ixpt)
                fmihxpt_tmp(iusp,ixpt)=fmihxpt_tmp(iusp,ixpt)+fmihxpt(iusp,ixpt)
                vyvxpt_tmp(iusp,ixpt)=vyvxpt_tmp(iusp,ixpt)+vyvxpt(iusp,ixpt)
                nixpt_tmp(iusp,ixpt)=nixpt_tmp(iusp,ixpt)+nixpt(iusp,ixpt)
                vyhxpt_tmp(iusp,ixpt)=vyhxpt_tmp(iusp,ixpt)+vyhxpt(iusp,ixpt)
                visyxpt_tmp(iusp,ixpt)=visyxpt_tmp(iusp,ixpt)+visyxpt(iusp,ixpt)
                upxpt_tmp(iusp,ixpt)=upxpt_tmp(iusp,ixpt)+upxpt(iusp,ixpt)
            end do
        end do 

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        fmixy_tmp(xc,yc,:)=fmixy_tmp(xc,yc,:)+fmixy(xc,yc,:)
        fmix_tmp(xc,yc,:)=fmix_tmp(xc,yc,:)+fmix(xc,yc,:)
        fmiy_tmp(xc,yc,:)=fmiy_tmp(xc,yc,:)+fmiy(xc,yc,:)
        smoc_tmp(xc,yc,:)=smoc_tmp(xc,yc,:)+smoc(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    fmivxpt=fmivxpt_tmp
    fmihxpt=fmihxpt_tmp; fmixy=fmixy_tmp; vyvxpt=vyvxpt_tmp
    nixpt=nixpt_tmp; vyhxpt=vyhxpt_tmp; visyxpt=visyxpt_tmp; upxpt=upxpt_tmp
    smoc=smoc_tmp; fmix=fmix_tmp; fmiy=fmiy_tmp
    ! TODO: Come up with a more sustainable solution for the fd2tra call
    call initialize_ranges(-1, -1, 0, 0, 0)
      do ifld = 1, nusp
      if(isupon(ifld) .ne. 0) then
         call fd2tra (nx,ny,flox(:,:,ifld),floy(:,:,ifld),conx(:,:,ifld),cony(:,:,ifld), &
         &            up(0:nx+1,0:ny+1,ifld),fmix(0:nx+1,0:ny+1,ifld), &
         &            fmiy(0:nx+1,0:ny+1,ifld),1, methu)
      end if
    end do

    call OmpCopyPointerfmivxpt; call OmpCopyPointerfmihxpt
    call OmpCopyPointerfmixy
    call OmpCopyPointervyvxpt; call OmpCopyPointernixpt
    call OmpCopyPointervyhxpt; call OmpCopyPointervisyxpt
    call OmpCopyPointerupxpt; call OmpCopyPointersmoc
    call OmpCopyPointerfmix; call OmpCopyPointerfmiy
    fmivxpt_cp=fmivxpt;fmihxpt_cp=fmihxpt;fmixy_cp=fmixy;vyvxpt_cp=vyvxpt
    nixpt_cp=nixpt;vyhxpt_cp=vyhxpt;visyxpt_cp=visyxpt;upxpt_cp=upxpt;smoc_cp=smoc
    fmix_cp=fmix;fmiy_cp=fmiy
  END SUBROUTINE OMPcalc_plasma_momentum


  SUBROUTINE OMPcalc_plasma_energy(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nzspmx, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: wvh, vsoree, edisse, seak, seik, emolia, pwribkg, psicx, seid, sead, seit, &
    &    pwrebkg, vsoreec, seidh, seadh, erliz
    USE Imprad, ONLY: prad, pwrzec, pwrze, ntau, pradcff, na, pradc, pradzc, pradz, nratio, nzloc
    USE Conduc, ONLY: eeli, pradhyd
    USE Comflo, ONLY: floxibgt, feexy, feeycbo, feey4ord, feixy, qipar, feey, feiy4ord, feiycbo, &
    &    floxebgt, feex, feiy, feix
    USE Wkspace, ONLY: w1, w0
    USE Locflux, ONLY: floye, conxi, floxe, conxe, floyi, conye, conyi, floxi
    USE MCN_sources, ONLY: seg_ue
    USE MCN_dim, ONLY: nfl
    USE MCN_sources, ONLY: seg_ue
    USE Locflux, ONLY: floye,conxi,floxe,conxe,floyi,conye,conyi,floxi
    USE Imprad, ONLY: prad,pwrzec,pwrze,ntau,pradcff,na,pradc,pradzc,pradz,nratio,nzloc
    USE Wkspace, ONLY: w1,w0
    USE Rhsides, ONLY: psorrgc,psordis,psordisg,psorc,psor,psorrg,vsoree,edisse,seak,seik,emolia, &
    &    pwribkg,psicx,seid,sead,seit,pwrebkg,vsoreec,seidh,seadh,erliz
    USE Comflo, ONLY: fqp,fngy,fngx,qipar,fniy,fnix,floxibgt,feexy,feeycbo,feey4ord,feixy,feey, &
    &    feiy4ord,feiycbo,floxebgt,feex,feiy,feix
    USE Conduc, ONLY: kyi_use,kye_use,vyte_cft,vyti_cft,nucx,hcyi,hcyn,hcye,hcxi,hcxe,eeli,pradhyd
    USE Gradients, ONLY: ex,ey,gtex
    USE Compla, ONLY: up,ne,ti,tg,te,nit,ng,ni,niy1,ney0,niy0,tiv,ney1,tev,vy,v2,vey,vex,rtau

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: wvh_tmp(0:nx+1,0:ny+1,1:nusp), prad_tmp(0:nx+1,0:ny+1), &
    &      pwrzec_tmp(0:nx+1,0:ny+1), eeli_tmp(0:nx+1,0:ny+1), &
    &      floxibgt_tmp(0:nx+1,0:ny+1,1:nisp), vsoree_tmp(0:nx+1,0:ny+1), &
    &      w1_tmp(0:nx+1,0:ny+1), pwrze_tmp(0:nx+1,0:ny+1), &
    &      edisse_tmp(0:nx+1,0:ny+1), feexy_tmp(0:nx+1,0:ny+1), feeycbo_tmp(0:nx+1), &
    &      feey4ord_tmp(0:nx+1,0:ny+1), seak_tmp(0:nx+1,0:ny+1), &
    &      seik_tmp(0:nx+1,0:ny+1), emolia_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      ntau_tmp(0:nx+1,0:ny+1), pradcff_tmp(0:nx+1,0:ny+1), &
    &      pwribkg_tmp(0:nx+1,0:ny+1), psicx_tmp(0:nx+1,0:ny+1), &
    &      floye_tmp(0:nx+1,0:ny+1), na_tmp(0:nx+1,0:ny+1), &
    &      feixy_tmp(0:nx+1,0:ny+1), pradc_tmp(0:nx+1,0:ny+1), &
    &      pradzc_tmp(0:nx+1,0:ny+1,0:nzspmx,1:nzspmx+1), &
    &      qipar_tmp(0:nx+1,0:ny+1,nisp), pradhyd_tmp(0:nx+1,0:ny+1), &
    &      seid_tmp(0:nx+1,0:ny+1), conxi_tmp(0:nx+1,0:ny+1), &
    &      floxe_tmp(0:nx+1,0:ny+1), sead_tmp(0:nx+1,0:ny+1), &
    &      conxe_tmp(0:nx+1,0:ny+1), seit_tmp(0:nx+1,0:ny+1), &
    &      pwrebkg_tmp(0:nx+1,0:ny+1), &
    &      pradz_tmp(0:nx+1,0:ny+1,0:nzspmx,1:nzspmx+1), feey_tmp(0:nx+1,0:ny+1), &
    &      seg_ue_tmp(0:nx+1,0:ny+1,1:nfl), vsoreec_tmp(0:nx+1,0:ny+1), &
    &      floyi_tmp(0:nx+1,0:ny+1), seidh_tmp(0:nx+1,0:ny+1), &
    &      seadh_tmp(0:nx+1,0:ny+1), conye_tmp(0:nx+1,0:ny+1), &
    &      conyi_tmp(0:nx+1,0:ny+1), feiy4ord_tmp(0:nx+1,0:ny+1), &
    &      feiycbo_tmp(0:nx+1), floxebgt_tmp(0:nx+1,0:ny+1), &
    &      feex_tmp(0:nx+1,0:ny+1), nratio_tmp(0:nx+1,0:ny+1), &
    &      feiy_tmp(0:nx+1,0:ny+1), feix_tmp(0:nx+1,0:ny+1), &
    &      floxi_tmp(0:nx+1,0:ny+1), erliz_tmp(0:nx+1,0:ny+1), &
    &      w0_tmp(0:nx+1,0:ny+1), nzloc_tmp(0:nzspmx)

    ! Initialize arrays to zero
    wvh_tmp=0.; prad_tmp=0.; pwrzec_tmp=0.; eeli_tmp=0.; floxibgt_tmp=0.
    vsoree_tmp=0.; w1_tmp=0.; pwrze_tmp=0.; edisse_tmp=0.; feexy_tmp=0.
    feeycbo_tmp=0.; feey4ord_tmp=0.; seak_tmp=0.; seik_tmp=0.; emolia_tmp=0.
    ntau_tmp=0.; pradcff_tmp=0.; pwribkg_tmp=0.; psicx_tmp=0.; floye_tmp=0.
    na_tmp=0.; feixy_tmp=0.; pradc_tmp=0.; pradzc_tmp=0.; qipar_tmp=0.
    pradhyd_tmp=0.; seid_tmp=0.; conxi_tmp=0.; floxe_tmp=0.; sead_tmp=0.
    conxe_tmp=0.; seit_tmp=0.; pwrebkg_tmp=0.; pradz_tmp=0.; feey_tmp=0.
    seg_ue_tmp=0.; vsoreec_tmp=0.; floyi_tmp=0.; seidh_tmp=0.; seadh_tmp=0.
    conye_tmp=0.; conyi_tmp=0.; feiy4ord_tmp=0.; feiycbo_tmp=0.
    floxebgt_tmp=0.; feex_tmp=0.; nratio_tmp=0.; feiy_tmp=0.; feix_tmp=0.
    floxi_tmp=0.; erliz_tmp=0.; w0_tmp=0.; nzloc_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:wvh_tmp, prad_tmp, pwrzec_tmp, eeli_tmp, floxibgt_tmp, vsoree_tmp, w1_tmp, &
    !$OMP &         pwrze_tmp, edisse_tmp, feexy_tmp, feeycbo_tmp, feey4ord_tmp, seak_tmp, &
    !$OMP &         seik_tmp, emolia_tmp, ntau_tmp, pradcff_tmp, pwribkg_tmp, psicx_tmp, &
    !$OMP &         floye_tmp, na_tmp, feixy_tmp, pradc_tmp, pradzc_tmp, qipar_tmp, pradhyd_tmp, &
    !$OMP &         seid_tmp, conxi_tmp, floxe_tmp, sead_tmp, conxe_tmp, seit_tmp, pwrebkg_tmp, &
    !$OMP &         pradz_tmp, feey_tmp, seg_ue_tmp, vsoreec_tmp, floyi_tmp, seidh_tmp, &
    !$OMP &         seadh_tmp, conye_tmp, conyi_tmp, feiy4ord_tmp, feiycbo_tmp, floxebgt_tmp, &
    !$OMP &         feex_tmp, nratio_tmp, feiy_tmp, feix_tmp, floxi_tmp, erliz_tmp, w0_tmp, &
    !$OMP &         nzloc_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ne=ne_cp;ti=ti_cp;tg=tg_cp;te=te_cp;nit=nit_cp;ng=ng_cp;ni=ni_cp
        niy1=niy1_cp;ex=ex_cp;ey=ey_cp;gtex=gtex_cp;ney0=ney0_cp;niy0=niy0_cp;tiv=tiv_cp
        ney1=ney1_cp;tev=tev_cp;kyi_use=kyi_use_cp;kye_use=kye_use_cp
        vyte_cft=vyte_cft_cp;vy=vy_cp;vyti_cft=vyti_cft_cp;v2=v2_cp;fqp=fqp_cp
        vey=vey_cp;vex=vex_cp;psorrgc=psorrgc_cp;psordis=psordis_cp;psordisg=psordisg_cp
        psorc=psorc_cp;psor=psor_cp;psorrg=psorrg_cp;rtau=rtau_cp;nucx=nucx_cp
        fngy=fngy_cp;fngx=fngx_cp;hcyi=hcyi_cp;hcyn=hcyn_cp;hcye=hcye_cp;qipar=qipar_cp
        hcxi=hcxi_cp;w1=w1_cp;hcxe=hcxe_cp;fniy=fniy_cp;fnix=fnix_cp;prad=prad_cp
        pwrzec=pwrzec_cp;eeli=eeli_cp;floxibgt=floxibgt_cp;vsoree=vsoree_cp
        pwrze=pwrze_cp;edisse=edisse_cp;feexy=feexy_cp;feeycbo=feeycbo_cp
        feey4ord=feey4ord_cp;seak=seak_cp;seik=seik_cp;emolia=emolia_cp;ntau=ntau_cp
        pradcff=pradcff_cp;pwribkg=pwribkg_cp;psicx=psicx_cp;floye=floye_cp;na=na_cp
        feixy=feixy_cp;pradc=pradc_cp;pradzc=pradzc_cp;pradhyd=pradhyd_cp;seid=seid_cp
        conxi=conxi_cp;floxe=floxe_cp;sead=sead_cp;conxe=conxe_cp;seit=seit_cp
        pwrebkg=pwrebkg_cp;pradz=pradz_cp;feey=feey_cp;seg_ue=seg_ue_cp
        vsoreec=vsoreec_cp;floyi=floyi_cp;seidh=seidh_cp;seadh=seadh_cp;conye=conye_cp
        conyi=conyi_cp;feiy4ord=feiy4ord_cp;feiycbo=feiycbo_cp;floxebgt=floxebgt_cp
        feex=feex_cp;nratio=nratio_cp;feiy=feiy_cp;feix=feix_cp;floxi=floxi_cp
        erliz=erliz_cp;w0=w0_cp;nzloc=nzloc_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_energy(-1,-1)


        nzloc_tmp=nzloc_tmp+nzloc
        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
        if (yc .eq. 0) then
            feeycbo_tmp(xc)=feeycbo_tmp(xc)+feeycbo(xc)
            feiycbo_tmp(xc)=feiycbo_tmp(xc)+feiycbo(xc)
        end if
         ! Update locally calculated variables
        wvh_tmp(xc,yc,:)=wvh_tmp(xc,yc,:)+wvh(xc,yc,:)
        prad_tmp(xc,yc)=prad_tmp(xc,yc)+prad(xc,yc)
        pwrzec_tmp(xc,yc)=pwrzec_tmp(xc,yc)+pwrzec(xc,yc)
        eeli_tmp(xc,yc)=eeli_tmp(xc,yc)+eeli(xc,yc)
        floxibgt_tmp(xc,yc,:)=floxibgt_tmp(xc,yc,:)+floxibgt(xc,yc,:)
        vsoree_tmp(xc,yc)=vsoree_tmp(xc,yc)+vsoree(xc,yc)
        w1_tmp(xc,yc)=w1_tmp(xc,yc)+w1(xc,yc)
        pwrze_tmp(xc,yc)=pwrze_tmp(xc,yc)+pwrze(xc,yc)
        edisse_tmp(xc,yc)=edisse_tmp(xc,yc)+edisse(xc,yc)
        feexy_tmp(xc,yc)=feexy_tmp(xc,yc)+feexy(xc,yc)
        feey4ord_tmp(xc,yc)=feey4ord_tmp(xc,yc)+feey4ord(xc,yc)
        seak_tmp(xc,yc)=seak_tmp(xc,yc)+seak(xc,yc)
        seik_tmp(xc,yc)=seik_tmp(xc,yc)+seik(xc,yc)
        emolia_tmp(xc,yc,:)=emolia_tmp(xc,yc,:)+emolia(xc,yc,:)
        ntau_tmp(xc,yc)=ntau_tmp(xc,yc)+ntau(xc,yc)
        pradcff_tmp(xc,yc)=pradcff_tmp(xc,yc)+pradcff(xc,yc)
        pwribkg_tmp(xc,yc)=pwribkg_tmp(xc,yc)+pwribkg(xc,yc)
        psicx_tmp(xc,yc)=psicx_tmp(xc,yc)+psicx(xc,yc)
        floye_tmp(xc,yc)=floye_tmp(xc,yc)+floye(xc,yc)
        na_tmp(xc,yc)=na_tmp(xc,yc)+na(xc,yc)
        feixy_tmp(xc,yc)=feixy_tmp(xc,yc)+feixy(xc,yc)
        pradc_tmp(xc,yc)=pradc_tmp(xc,yc)+pradc(xc,yc)
        pradzc_tmp(xc,yc,:,:)=pradzc_tmp(xc,yc,:,:)+pradzc(xc,yc,:,:)
        qipar_tmp(xc,yc,:)=qipar_tmp(xc,yc,:)+qipar(xc,yc,:)
        pradhyd_tmp(xc,yc)=pradhyd_tmp(xc,yc)+pradhyd(xc,yc)
        seid_tmp(xc,yc)=seid_tmp(xc,yc)+seid(xc,yc)
        conxi_tmp(xc,yc)=conxi_tmp(xc,yc)+conxi(xc,yc)
        floxe_tmp(xc,yc)=floxe_tmp(xc,yc)+floxe(xc,yc)
        sead_tmp(xc,yc)=sead_tmp(xc,yc)+sead(xc,yc)
        conxe_tmp(xc,yc)=conxe_tmp(xc,yc)+conxe(xc,yc)
        seit_tmp(xc,yc)=seit_tmp(xc,yc)+seit(xc,yc)
        pwrebkg_tmp(xc,yc)=pwrebkg_tmp(xc,yc)+pwrebkg(xc,yc)
        pradz_tmp(xc,yc,:,:)=pradz_tmp(xc,yc,:,:)+pradz(xc,yc,:,:)
        feey_tmp(xc,yc)=feey_tmp(xc,yc)+feey(xc,yc)
        seg_ue_tmp(xc,yc,:)=seg_ue_tmp(xc,yc,:)+seg_ue(xc,yc,:)
        vsoreec_tmp(xc,yc)=vsoreec_tmp(xc,yc)+vsoreec(xc,yc)
        floyi_tmp(xc,yc)=floyi_tmp(xc,yc)+floyi(xc,yc)
        seidh_tmp(xc,yc)=seidh_tmp(xc,yc)+seidh(xc,yc)
        seadh_tmp(xc,yc)=seadh_tmp(xc,yc)+seadh(xc,yc)
        conye_tmp(xc,yc)=conye_tmp(xc,yc)+conye(xc,yc)
        conyi_tmp(xc,yc)=conyi_tmp(xc,yc)+conyi(xc,yc)
        feiy4ord_tmp(xc,yc)=feiy4ord_tmp(xc,yc)+feiy4ord(xc,yc)
        floxebgt_tmp(xc,yc)=floxebgt_tmp(xc,yc)+floxebgt(xc,yc)
        feex_tmp(xc,yc)=feex_tmp(xc,yc)+feex(xc,yc)
        nratio_tmp(xc,yc)=nratio_tmp(xc,yc)+nratio(xc,yc)
        feiy_tmp(xc,yc)=feiy_tmp(xc,yc)+feiy(xc,yc)
        feix_tmp(xc,yc)=feix_tmp(xc,yc)+feix(xc,yc)
        floxi_tmp(xc,yc)=floxi_tmp(xc,yc)+floxi(xc,yc)
        erliz_tmp(xc,yc)=erliz_tmp(xc,yc)+erliz(xc,yc)
        w0_tmp(xc,yc)=w0_tmp(xc,yc)+w0(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    wvh=wvh_tmp; prad=prad_tmp; pwrzec=pwrzec_tmp; eeli=eeli_tmp
    floxibgt=floxibgt_tmp; vsoree=vsoree_tmp; w1=w1_tmp; pwrze=pwrze_tmp
    edisse=edisse_tmp; feexy=feexy_tmp; feeycbo=feeycbo_tmp
    feey4ord=feey4ord_tmp; seak=seak_tmp; seik=seik_tmp; emolia=emolia_tmp
    ntau=ntau_tmp; pradcff=pradcff_tmp; pwribkg=pwribkg_tmp; psicx=psicx_tmp
    floye=floye_tmp; na=na_tmp; feixy=feixy_tmp; pradc=pradc_tmp
    pradzc=pradzc_tmp; qipar=qipar_tmp; pradhyd=pradhyd_tmp; seid=seid_tmp
    conxi=conxi_tmp; floxe=floxe_tmp; sead=sead_tmp; conxe=conxe_tmp
    seit=seit_tmp; pwrebkg=pwrebkg_tmp; pradz=pradz_tmp; feey=feey_tmp
    seg_ue=seg_ue_tmp; vsoreec=vsoreec_tmp; floyi=floyi_tmp; seidh=seidh_tmp
    seadh=seadh_tmp; conye=conye_tmp; conyi=conyi_tmp; feiy4ord=feiy4ord_tmp
    feiycbo=feiycbo_tmp; floxebgt=floxebgt_tmp; feex=feex_tmp
    nratio=nratio_tmp; feiy=feiy_tmp; feix=feix_tmp; floxi=floxi_tmp
    erliz=erliz_tmp; w0=w0_tmp; nzloc=nzloc_tmp
    call OmpCopyPointerwvh; call OmpCopyPointerprad
    call OmpCopyPointerpwrzec; call OmpCopyPointereeli
    call OmpCopyPointerfloxibgt; call OmpCopyPointervsoree
    call OmpCopyPointerw1; call OmpCopyPointerpwrze
    call OmpCopyPointeredisse; call OmpCopyPointerfeexy
    call OmpCopyPointerfeeycbo; call OmpCopyPointerfeey4ord
    call OmpCopyPointerseak; call OmpCopyPointerseik
    call OmpCopyPointeremolia; call OmpCopyPointerntau
    call OmpCopyPointerpradcff; call OmpCopyPointerpwribkg
    call OmpCopyPointerpsicx; call OmpCopyPointerfloye
    call OmpCopyPointerna; call OmpCopyPointerfeixy
    call OmpCopyPointerpradc; call OmpCopyPointerpradzc
    call OmpCopyPointerqipar; call OmpCopyPointerpradhyd
    call OmpCopyPointerseid; call OmpCopyPointerconxi
    call OmpCopyPointerfloxe; call OmpCopyPointersead
    call OmpCopyPointerconxe; call OmpCopyPointerseit
    call OmpCopyPointerpwrebkg; call OmpCopyPointerpradz
    call OmpCopyPointerfeey; call OmpCopyPointerseg_ue
    call OmpCopyPointervsoreec; call OmpCopyPointerfloyi
    call OmpCopyPointerseidh; call OmpCopyPointerseadh
    call OmpCopyPointerconye; call OmpCopyPointerconyi
    call OmpCopyPointerfeiy4ord; call OmpCopyPointerfeiycbo
    call OmpCopyPointerfloxebgt; call OmpCopyPointerfeex
    call OmpCopyPointernratio; call OmpCopyPointerfeiy
    call OmpCopyPointerfeix; call OmpCopyPointerfloxi
    call OmpCopyPointererliz; call OmpCopyPointerw0
    call OmpCopyPointernzloc
    wvh_cp=wvh;prad_cp=prad;pwrzec_cp=pwrzec;eeli_cp=eeli;floxibgt_cp=floxibgt
    vsoree_cp=vsoree;w1_cp=w1;pwrze_cp=pwrze;edisse_cp=edisse;feexy_cp=feexy
    feeycbo_cp=feeycbo;feey4ord_cp=feey4ord;seak_cp=seak;seik_cp=seik
    emolia_cp=emolia;ntau_cp=ntau;pradcff_cp=pradcff;pwribkg_cp=pwribkg
    psicx_cp=psicx;floye_cp=floye;na_cp=na;feixy_cp=feixy;pradc_cp=pradc
    pradzc_cp=pradzc;qipar_cp=qipar;pradhyd_cp=pradhyd;seid_cp=seid;conxi_cp=conxi
    floxe_cp=floxe;sead_cp=sead;conxe_cp=conxe;seit_cp=seit;pwrebkg_cp=pwrebkg
    pradz_cp=pradz;feey_cp=feey;seg_ue_cp=seg_ue;vsoreec_cp=vsoreec;floyi_cp=floyi
    seidh_cp=seidh;seadh_cp=seadh;conye_cp=conye;conyi_cp=conyi;feiy4ord_cp=feiy4ord
    feiycbo_cp=feiycbo;floxebgt_cp=floxebgt;feex_cp=feex;nratio_cp=nratio
    feiy_cp=feiy;feix_cp=feix;floxi_cp=floxi;erliz_cp=erliz;w0_cp=w0;nzloc_cp=nzloc
  END SUBROUTINE OMPcalc_plasma_energy


  SUBROUTINE OMPcalc_gas_energy(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: seic, eiamoldiss
    USE Conduc, ONLY: eqpg
    USE Conduc, ONLY: eqpg
    USE Rhsides, ONLY: psordis,psordisg,seic,eiamoldiss
    USE Compla, ONLY: up,ti,tg,ng,ni,pg,vy,v2,vyg,uuxg

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: seic_tmp(0:nx+1,0:ny+1), eqpg_tmp(0:nx+1,0:ny+1,ngsp), &
    &      eiamoldiss_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    seic_tmp=0.; eqpg_tmp=0.; eiamoldiss_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:seic_tmp, eqpg_tmp, eiamoldiss_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ti=ti_cp;tg=tg_cp;ng=ng_cp;ni=ni_cp;pg=pg_cp;vy=vy_cp;v2=v2_cp
        psordis=psordis_cp;psordisg=psordisg_cp;seic=seic_cp;vyg=vyg_cp;uuxg=uuxg_cp
        eqpg=eqpg_cp;eiamoldiss=eiamoldiss_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_gas_energy

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        seic_tmp(xc,yc)=seic_tmp(xc,yc)+seic(xc,yc)
        eqpg_tmp(xc,yc,:)=eqpg_tmp(xc,yc,:)+eqpg(xc,yc,:)
        eiamoldiss_tmp(xc,yc,:)=eiamoldiss_tmp(xc,yc,:)+eiamoldiss(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    seic=seic_tmp; eqpg=eqpg_tmp; eiamoldiss=eiamoldiss_tmp
    call OmpCopyPointerseic; call OmpCopyPointereqpg
    call OmpCopyPointereiamoldiss
    seic_cp=seic;eqpg_cp=eqpg;eiamoldiss_cp=eiamoldiss
  END SUBROUTINE OMPcalc_gas_energy


  SUBROUTINE OMPcalc_plasma_particle_residuals(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE MCN_sources, ONLY: sng_ue
    USE MCN_dim, ONLY: nfl
    USE Rhsides, ONLY: resco
    USE MCN_sources, ONLY: sng_ue
    USE Comflo, ONLY: fniy,fnix
    USE Rhsides, ONLY: snic,psor,psorxr,resco
    USE Conduc, ONLY: nuvl
    USE Compla, ONLY: ti,ng,ni

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: sng_ue_tmp(0:nx+1,0:ny+1,1:nfl), resco_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    sng_ue_tmp=0.; resco_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:sng_ue_tmp, resco_tmp)
    DO ichunk = 1, NchunksPandf1
        ti=ti_cp;ng=ng_cp;ni=ni_cp;nuvl=nuvl_cp;snic=snic_cp;psor=psor_cp
        psorxr=psorxr_cp;fniy=fniy_cp;fnix=fnix_cp;sng_ue=sng_ue_cp;resco=resco_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_particle_residuals

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        sng_ue_tmp(xc,yc,:)=sng_ue_tmp(xc,yc,:)+sng_ue(xc,yc,:)
        resco_tmp(xc,yc,:)=resco_tmp(xc,yc,:)+resco(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    sng_ue=sng_ue_tmp; resco=resco_tmp
    call OmpCopyPointersng_ue; call OmpCopyPointerresco
    sng_ue_cp=sng_ue;resco_cp=resco
  END SUBROUTINE OMPcalc_plasma_particle_residuals


  SUBROUTINE OMPcalc_gas_continuity_residuals(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: resng
    USE MCN_sources, ONLY: sng_ue
    USE MCN_dim, ONLY: nfl
    USE MCN_sources, ONLY: sng_ue
    USE Comflo, ONLY: fngy,fngx,fnix
    USE Conduc, ONLY: nuiz,nucx
    USE Rhsides, ONLY: psordis,psorg,psorcxg,psorrg,resng
    USE Compla, ONLY: ti,tg,ng

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: resng_tmp(0:nx+1,0:ny+1,1:ngsp), sng_ue_tmp(0:nx+1,0:ny+1,1:nfl)

    ! Initialize arrays to zero
    resng_tmp=0.; sng_ue_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:resng_tmp, sng_ue_tmp)
    DO ichunk = 1, NchunksPandf1
        ti=ti_cp;tg=tg_cp;ng=ng_cp;psordis=psordis_cp;nuiz=nuiz_cp;psorg=psorg_cp
        psorcxg=psorcxg_cp;psorrg=psorrg_cp;nucx=nucx_cp;fngy=fngy_cp;fngx=fngx_cp
        fnix=fnix_cp;sng_ue=sng_ue_cp;resng=resng_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_gas_continuity_residuals

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        resng_tmp(xc,yc,:)=resng_tmp(xc,yc,:)+resng(xc,yc,:)
        sng_ue_tmp(xc,yc,:)=sng_ue_tmp(xc,yc,:)+sng_ue(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    resng=resng_tmp; sng_ue=sng_ue_tmp
    call OmpCopyPointerresng; call OmpCopyPointersng_ue
    resng_cp=resng;sng_ue_cp=sng_ue
  END SUBROUTINE OMPcalc_gas_continuity_residuals


  SUBROUTINE OMPcalc_plasma_momentum_residuals(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Wkspace, ONLY: w0, w2
    USE Cfric, ONLY: fricnrl
    USE Rhsides, ONLY: resmo
    USE Cfric, ONLY: fricnrl
    USE Comflo, ONLY: fmixy,fmix,fmiy
    USE Wkspace, ONLY: w2,w0
    USE Rhsides, ONLY: smoc,msor,msorxr,resmo
    USE Conduc, ONLY: nurc,nuiz,nucxi,nucx,nueli
    USE Compla, ONLY: up,ti,tg,ng,nm,ni,loglambda

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: w0_tmp(0:nx+1,0:ny+1), fricnrl_tmp(0:nx+1,0:ny+1,nusp), &
    &      resmo_tmp(0:nx+1,0:ny+1,1:nusp), w2_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    w0_tmp=0.; fricnrl_tmp=0.; resmo_tmp=0.; w2_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:w0_tmp, fricnrl_tmp, resmo_tmp, w2_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ti=ti_cp;tg=tg_cp;ng=ng_cp;nm=nm_cp;ni=ni_cp;loglambda=loglambda_cp
        nurc=nurc_cp;smoc=smoc_cp;nuiz=nuiz_cp;nucxi=nucxi_cp;msor=msor_cp
        msorxr=msorxr_cp;nucx=nucx_cp;nueli=nueli_cp;w2=w2_cp;fmixy=fmixy_cp
        fmix=fmix_cp;fmiy=fmiy_cp;w0=w0_cp;fricnrl=fricnrl_cp;resmo=resmo_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_momentum_residuals

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        w0_tmp(xc,yc)=w0_tmp(xc,yc)+w0(xc,yc)
        fricnrl_tmp(xc,yc,:)=fricnrl_tmp(xc,yc,:)+fricnrl(xc,yc,:)
        resmo_tmp(xc,yc,:)=resmo_tmp(xc,yc,:)+resmo(xc,yc,:)
        w2_tmp(xc,yc)=w2_tmp(xc,yc)+w2(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    w0=w0_tmp; fricnrl=fricnrl_tmp; resmo=resmo_tmp; w2=w2_tmp
    call OmpCopyPointerw0; call OmpCopyPointerfricnrl
    call OmpCopyPointerresmo; call OmpCopyPointerw2
    w0_cp=w0;fricnrl_cp=fricnrl;resmo_cp=resmo;w2_cp=w2
  END SUBROUTINE OMPcalc_plasma_momentum_residuals


  SUBROUTINE OMPcalc_gas_energy_residuals(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: reseg
    USE MCN_sources, ONLY: seg_ue
    USE Comflo, ONLY: fegx,fegy
    USE Conduc, ONLY: visy,visx,eqpg
    USE Rhsides, ONLY: psordis,psorg,segc,wvh,seak,sead,seit,seadh,eiamoldiss,reseg
    USE Compla, ONLY: up,ti,tg,ng,pg,vy,v2,upi,vyg,uuxg

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: reseg_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    reseg_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:reseg_tmp)
    DO ichunk = 1, NchunksPandf1
        up=up_cp;ti=ti_cp;tg=tg_cp;ng=ng_cp;pg=pg_cp;vy=vy_cp;v2=v2_cp;upi=upi_cp
        psordis=psordis_cp;psorg=psorg_cp;vyg=vyg_cp;uuxg=uuxg_cp;visy=visy_cp
        visx=visx_cp;segc=segc_cp;fegx=fegx_cp;fegy=fegy_cp;wvh=wvh_cp;seak=seak_cp
        sead=sead_cp;seit=seit_cp;seg_ue=seg_ue_cp;seadh=seadh_cp;eqpg=eqpg_cp
        eiamoldiss=eiamoldiss_cp;reseg=reseg_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_gas_energy_residuals

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        reseg_tmp(xc,yc,:)=reseg_tmp(xc,yc,:)+reseg(xc,yc,:)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    reseg=reseg_tmp
    call OmpCopyPointerreseg
    reseg_cp=reseg
  END SUBROUTINE OMPcalc_gas_energy_residuals


  SUBROUTINE OMPcalc_plasma_energy_residuals(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Wkspace, ONLY: w0
    USE Rhsides, ONLY: resei, wjdote, reseg, resee, wvh
    USE Wkspace, ONLY: w0
    USE MCN_sources, ONLY: seg_ue
    USE Imprad, ONLY: pwrze
    USE Conduc, ONLY: nuvl,nucxi,nucx,nueli,visy,visx,eqp
    USE Rhsides, ONLY: psordis,seic,psor,seec,wvh,vsoree,seik,emolia,pwribkg,seid,seit,pwrebkg,seidh, &
    &    resei,resee
    USE Comflo, ONLY: fqygp,fq2,fqy,fqp,fqx,feey,feex,feiy,feix
    USE Gradients, ONLY: ex,ey
    USE Compla, ONLY: up,ne,phi,ti,tg,te,ng,ni,vy,v2,upi

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: resei_tmp(0:nx+1,0:ny+1), &
    &      reseg_tmp(0:nx+1,0:ny+1), &
    &      resee_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    resei_tmp=0.; reseg_tmp=0.; resee_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:resei_tmp, reseg_tmp, resee_tmp)
    DO ichunk = 1, NchunksPandf1
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_plasma_energy_residuals(xc,yc)

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        resei_tmp(xc,yc)=resei_tmp(xc,yc)+resei(xc,yc)
        reseg_tmp(xc,yc)=reseg_tmp(xc,yc)+reseg(xc,yc,1)
        resee_tmp(xc,yc)=resee_tmp(xc,yc)+resee(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    resei=resei_tmp; reseg(:,:,1)=reseg_tmp
    resee=resee_tmp; 
    call OmpCopyPointerresei
    call OmpCopyPointerreseg
    call OmpCopyPointerresee
    resei_cp=resei;reseg_cp=reseg;resee_cp=resee
  END SUBROUTINE OMPcalc_plasma_energy_residuals


  SUBROUTINE OMPcalc_potential_residuals(neq, yl, yldot)
    USE PandfCopies
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: resphi
    USE Rhsides, ONLY: resphi
    USE Comflo, ONLY: fqy,fqx
    USE Compla, ONLY: phi

    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: resphi_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    resphi_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:resphi_tmp)
    DO ichunk = 1, NchunksPandf1
        phi=phi_cp;fqy=fqy_cp;fqx=fqx_cp;resphi=resphi_cp
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_potential_residuals

        do ii = 1, Nixychunk(ichunk)
         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
         ! Update locally calculated variables
        resphi_tmp(xc,yc)=resphi_tmp(xc,yc)+resphi(xc,yc)
            end do
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    resphi=resphi_tmp
    call OmpCopyPointerresphi
    resphi_cp=resphi
  END SUBROUTINE OMPcalc_potential_residuals


  SUBROUTINE OMPbouncon(neq,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,NchunksPandf1, rangechunk, &
    &       Nxchunks, Nychunks
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE PandfCopies
    USE Dim, ONLY:nx,ny,nxpt, nusp
    USE Math_problem_size, ONLY: numvar
    USE Xpoint_indices, ONLY: ixrb, ixlb, ixpt2
    USE UEpar, ONLY: igas
    USE Indices_domain_dcl, ONLY: iymnbcl,iymxbcl, ixmnbcl, ixmxbcl
    USE Share, ONLY: isudsym, geometry, islimon, ix_lim, nxc
    USE Bcond, ONLY: isfixlb
    USE Indexes, ONLY: idxu
    USE Selec, ONLY: j3
    
    USE Wkspace, ONLY: w
    USE Conduc, ONLY: nuiz,nuix,visx
    USE Bcond, ONLY: fqpsatrb,fqpsatlb
    USE Poten, ONLY: dphi_iy1
    USE Comflo, ONLY: fdiaxlb,fdiaxrb,fqym,fqyb,fqyao,fqya,fqyd,fqy,fqp,fqx,fngy,fngx,fegx,fegy,fniy, &
    &    fniycbo,fnix,fmix,fmiy,floxibgt,feeycbo,feey,feiycbo,floxebgt,feex,feiy,feix
    USE Gradients, ONLY: ex,ey,gpiy
    USE Compla, ONLY: up,ne,phi,ti,tg,te,ng,nm,ni,niy1,niy0,vy,vytan,v2cd,v2ce,upi,uz,uu,vey,vex

    IMPLICIT NONE
 
    integer,intent(in)::neq
    real,intent(in)::yl(*)
    real,intent(out)::yldot(*)
    real::yldotcopy(1:neq)
    real ylcopy(1:neq+2), yldottot(1:neq), Nxchunks_old, Nychunks_old
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, ii, idx, ifld, jx

        yldotcopy = yldot(1:neq)
        yldottot = 0

        ! TODO: Figure out all these chunking issues w/ bouncon

        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, NchunksPandf1
            ! TODO: Figure out more generalized chunking routines for BCs
            call OMPinitialize_ranges2D(rangechunk(ichunk,:))
            if (rangechunk(ichunk,3).eq.0) j3=0
            call bouncon(neq, yldotcopy)

            do ii=1,Nivchunk(ichunk)
                yldottot(ivchunk(ichunk,ii)) = yldottot(ivchunk(ichunk,ii)) + yldotcopy(ivchunk(ichunk,ii))
            end do

        END DO
        !$OMP END PARALLEL DO
        yldot(1:neq) = yldottot(1:neq)

        ! TODO: Resolve these issues with the core boundary...
        call initialize_ranges(-1,-1,2,2,2)
        call iwall_boundary(neq, yldot)

        RETURN
  END SUBROUTINE OMPbouncon


  SUBROUTINE OMPcalc_rhs(neq,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,NchunksPandf1, rangechunk
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY:nx,ny,nisp
    USE Math_problem_size, ONLY: numvar
    IMPLICIT NONE
 
    integer,intent(in)::neq
    real,intent(in)::yl(*)
    real,intent(out)::yldot(*)
    real::yldotcopy(1:neq)
    real ylcopy(1:neq+2), yldottot(1:neq)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, ii

        yldotcopy = yldot(1:neq)
        yldottot = 0


        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, NchunksPandf1
            call OMPinitialize_ranges2D(rangechunk(ichunk,:))
            call calc_rhs(yldotcopy)
            
            do ii=1,Nivchunk(ichunk)
                yldottot(ivchunk(ichunk,ii)) = yldottot(ivchunk(ichunk,ii)) + yldotcopy(ivchunk(ichunk,ii))
            end do

        END DO
        !$OMP END PARALLEL DO
        yldot(1:neq) = yldottot(1:neq)


    RETURN
  END SUBROUTINE OMPcalc_rhs


  SUBROUTINE OMPrscalf(neq,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,NchunksPandf1, rangechunk
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY:nx,ny,nisp
    USE Math_problem_size, ONLY: numvar
    IMPLICIT NONE
 
    integer,intent(in)::neq
    real,intent(in)::yl(*)
    real,intent(out)::yldot(*)
    real::yldotcopy(1:neq)
    real ylcopy(1:neq+2), yldottot(1:neq)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, ii

        yldotcopy(1:neq) = yldot(1:neq)
        ylcopy(1:neq) = yl(1:neq)
        yldottot = 0


        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, NchunksPandf1
            call OMPinitialize_ranges2D(rangechunk(ichunk,:))
            call rscalf(ylcopy,yldotcopy)

            do ii=1,Nivchunk(ichunk)
                yldottot(ivchunk(ichunk,ii)) = yldottot(ivchunk(ichunk,ii)) + yldotcopy(ivchunk(ichunk,ii))
            end do

        END DO
        !$OMP END PARALLEL DO
        yldot(1:neq) = yldottot(1:neq)


    RETURN
  END SUBROUTINE OMPrscalf


  SUBROUTINE OMPadd_timestep(neq,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,NchunksPandf1, rangechunk
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY:nx,ny,nisp
    USE Math_problem_size, ONLY: numvar
    USE UEpar, ONLY: isbcwdt
    IMPLICIT NONE
 
    integer,intent(in)::neq
    real,intent(in)::yl(*)
    real,intent(out)::yldot(*)
    real::yldotcopy(1:neq)
    real ylcopy(1:neq+2), yldottot(1:neq)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, ii

        yldotcopy = yldot(1:neq)
        ylcopy = yl(1:neq)
        yldottot = 0


        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, NchunksPandf1
            call OMPinitialize_ranges2D(rangechunk(ichunk,:))
            call add_timestep(neq, ylcopy, yldotcopy)
            
            do ii=1,Nivchunk(ichunk)
                yldottot(ivchunk(ichunk,ii)) = yldottot(ivchunk(ichunk,ii)) + yldotcopy(ivchunk(ichunk,ii))
            end do

        END DO
        !$OMP END PARALLEL DO
        yldot(1:neq) = yldottot(1:neq)


    RETURN
  END SUBROUTINE OMPadd_timestep


  SUBROUTINE OMPPandf1Rhs(neq,time,yl,yldot)
! Recreates Pandf using parallel structure
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE Grid, ONLY:ijactot
    USE Cdv, ONLY: comnfe
    USE Rhsides, ONLY: psorcxg, psorrg, psordis
    USE Time_dep_nwt, ONLY: dtreal, nufak
    USE Ynorm, ONLY: isflxvar, isrscalf
    USE MCN_sources, ONLY: ismcnon
    USE PandfTiming, ONLY: TimePandf, TotTimePandf, TimingPandfOn, TimeNeudif, &
    &   TotTimeNeudif
    USE OMPTiming, ONLY: ParaTime, SerialTime
    USE Coefeq, ONLY: cfvisxneov, cfvisxneoq    
    USE Math_problem_size, ONLY: numvar
    USE ParallelEval, ONLY: ParallelPandfCall

    USE UEpar, ONLY: igas
    USE Xpoint_indices, ONLY: ixrb, ixlb
    USE Indices_domain_dcl, ONLY: iymnbcl,iymxbcl, ixmnbcl, ixmxbcl
    USE Share, ONLY: isudsym, geometry, islimon, ix_lim, nxc
    USE Bcond, ONLY: isfixlb
    USE OMPPandf1, ONLY: Nxchunks, Nychunks
    USE CHUNK

    USE Indexes, ONLY: igyl
    IMPLICIT NONE
    INTEGER:: Nchunks, Nchunks_convert
!    INTEGER, ALLOCATABLE, DIMENSION(:,:,:):: ixychunk, ixychunk_convert
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: rangexptchunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) ::  ivxptchunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  ivchunk, rangechunk, ivchunk_convert, &
    &   rangechunk_convert, Nivxptchunk
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  Nivchunk, Nivchunk_convert
!    &   Nixychunk_convert
!    INTEGER, ALLOCATABLE, DIMENSION(:)::  Nivchunk, Nixychunk, Nivchunk_convert, &
    INTEGER tid, Nxptchunks(nxpt)

 
    integer yinc_bkp,xrinc_bkp,xlinc_bkp,iv
    integer,intent(in)::neq
    real,intent(in)::yl(neq+1)
    real,intent(out)::yldot(neq)
    real,intent(in)::time
    real:: yldotcopy(neq),yldotxpt1(neq),yldotxpt2(neq) 
    real:: yldottot(neq), yldotxpt1tot(neq), yldotxpt2tot(neq) 
    real yldotsave(neq),ylcopy(neq+2)
    character*80 ::FileName
    real time1,time2
    INTEGER:: ichunk, xc, yc, ii, icut, ixpt
    real tmp_prad(0:nx+1, 0:ny+1)
    real tick,tock, tsfe, tsjf, ttotfe, ttotjf, tserial, tpara
    external tick, tock
    integer dummy(4)
    dummy(1) = 0
    dummy(2) = nx+1
    dummy(3) = 0
    dummy(4) = ny+1

    ParallelPandfCall = 1
    ylcopy = yl
    yldotcopy = yldot
    yldottot = 0
    yldotxpt1 = 0
    yldotxpt2 = 0

!        tpara = tick()
!        ParaTime = ParaTime + tock(tpara)

!        tserial = tick()
!        call initialize_ranges(xc, yc, xlinc, xrinc, yinc)
!        SerialTime = SerialTime + tock(tserial)

    Nxptchunks(1) = 1
    if (ijactot.gt.0) then
    ! TODO: move to initializer
        Time1=omp_get_wtime()
        call Make2DChunks(Nxchunks, Nychunks, Nchunks, Nivchunk, &
        &   ivchunk, rangechunk, Nxptchunks, Nivxptchunk, ivxptchunk, &
        &   rangexptchunk)


!        write(*,*) '==========='
!        do ii = 1, Nchunks
!            write(*,*) rangechunk(ii,:)
!        end do
!        call Make2DChunks(Nxchunks, Nychunks, Nchunks, Nivchunk, &
!        &   ivchunk, rangechunk, ixychunk, Nixychunk)

!        call Make2DChunks(nx, ny, Nchunks_convert, Nivchunk_convert, &
!        &   ivchunk_convert, rangechunk_convert, ixychunk_convert, & 
!        &   Nixychunk_convert)

        tpara = tick()
        !$OMP TEAMS NUM_TEAMS(nxpt+1) PRIVATE(tid)
        tid = omp_get_team_num()
        if ( omp_get_num_teams() .ne. nxpt+1 ) stop "Too few teams allocated"
        ! Do the first X-point - assumed to always be present
        if(tid .eq. 0) then
        !$OMP PARALLEL DO PRIVATE(ichunk) SCHEDULE(dynamic) &
        !$OMP &      FIRSTPRIVATE(ylcopy, yldotcopy) REDUCTION(+:yldotxpt1)
            DO ichunk = 1, Nxptchunks(1)
                psorcxg = 0
                psorrg = 0
                psordis = 0
                ! Only one chunk for the time being
                call OMPPandf_XPT(neq, ylcopy, yldotcopy, rangexptchunk(1, :, ichunk,:))

                do iv=1,Nivxptchunk(1, ichunk)
                    yldotxpt1(ivxptchunk(1,ichunk,iv)) = yldotxpt1(ivxptchunk(1,ichunk,iv)) &
                    &       + yldotcopy(ivxptchunk(1,ichunk,iv))
                enddo
            END DO
        !$OMP END PARALLEL DO
        ! Do the second X-point if present
        elseif ((nxpt.eq.2).and.(tid.eq.1)) then
        !$OMP PARALLEL DO PRIVATE(ichunk) SCHEDULE(dynamic) &
        !$OMP &      FIRSTPRIVATE(ylcopy, yldotcopy) REDUCTION(+:yldotxpt1)
            DO ichunk = 1, Nxptchunks(2)
                psorcxg = 0
                psorrg = 0
                psordis = 0
                ! Only one chunk for the time being
                call OMPPandf_XPT(neq, ylcopy, yldotcopy, rangexptchunk(2, :, ichunk,:))

                do iv=1,Nivxptchunk(2, ichunk)
                    yldotxpt1(ivxptchunk(2,ichunk,iv)) = yldotxpt1(ivxptchunk(2,ichunk,iv)) &
                    &       + yldotcopy(ivxptchunk(2,ichunk,iv))
                enddo
            END DO
        !$OMP END PARALLEL DO


        else 
            !$OMP PARALLEL
            !$OMP DO    PRIVATE(ichunk, iv) SCHEDULE(dynamic) &
            !$OMP &      FIRSTPRIVATE(ylcopy, yldotcopy) REDUCTION(+:yldottot)
            DO ichunk= 1, Nchunks

                psorcxg = 0
                psorrg = 0
                psordis = 0
                call OMPPandf(neq, ylcopy, yldotcopy,rangechunk(ichunk,:))

                do iv=1,Nivchunk(ichunk)
                    yldottot(ivchunk(ichunk,iv)) = yldottot(ivchunk(ichunk,iv)) &
                    &       + yldotcopy(ivchunk(ichunk,iv))
                enddo
            END DO
            !$OMP END DO
            !$OMP END PARALLEL
            END IF
            !$OMP END TEAMS

        yldot = yldottot !+ yldotcut
!        write(*,*) "===========", SUM(Nivchunk),SUM(Nivxptchunk)
        do ii=1,Nivxptchunk(1, 1)
            iv = ivxptchunk(1,1,ii)
            yldot(iv) = yldotxpt1(iv)
        end do
        if (nxpt.eq.2) then
            do ii=1,Nivxptchunk(2, 1)
                iv = ivxptchunk(2,1,ii)
                yldot(iv) = yldotxpt1(iv)
            end do
        endif
        do ixpt = 1, nxpt
            do ii=1,Nivxptchunk(ixpt, 1)
                    iv = ivxptchunk(ixpt,1,ii)
                    yldot(iv) = yldotxpt1(iv)
!                    write(*,*) iv, ii
!                    if (ABS((yldottot(iv)-yldotcut(iv))/yldottot(iv)).gt.1e-6) &
!                    &   write(*,*) igyl(iv,:), iv!yldotcut(iv), yldottot(iv)
            end do
        enddo
            
        
        ParaTime = ParaTime + tock(tpara)

!        call OMPPandf(neq, yl, yldot, rangechunk(1,:))
        Time1=omp_get_wtime()-Time1

        OMPTimeParallelPandf1=Time1+OMPTimeParallelPandf1
        if (CheckPandf1.gt.0) then
            Time2=omp_get_wtime()
            call pandf (-1, -1, neq, time, ylcopy, yldotsave)
            Time2=omp_get_wtime()-Time2
            OMPTimeSerialPandf1=Time2+OMPTimeSerialPandf1
            if (OMPPandf1Verbose.gt.0) then
                write(*,*) "Timing Pandf1 serial:",OMPTimeSerialPandf1, &
                    "(",Time2,")/parallel:",OMPTimeParallelPandf1,'(',Time1,')'
            endif
            call Compare(yldot,yldotsave,neq)
            write(*,'(a,i4)') "  Serial and parallel pandf are identical for nfe = ", comnfe
        endif
    else
       call pandf (-1,-1, neq, time, yl, yldot)
    endif
    ParallelPandfCall = 0
    RETURN

    ! TODO: add serial pandf call outside of loops to update all arrays
    END SUBROUTINE OMPPandf1Rhs


    SUBROUTINE OMPPandf(neq, yl, yldot, range)
      USE UEpar, ONLY: isphion, svrpkg, isphiofft
      USE PandfTiming, ONLY: TimePandf, TotTimePandf, TimingPandfOn, &
      &     TimeNeudif, TotTimeNeudif
      USE Ynorm, ONLY: isflxvar, isrscalf
      USE Time_dep_nwt, ONLY: dtreal
      USE Dim, ONLY: nxpt, nx
      USE Xpoint_indices, ONLY: ixpt1, ixpt2, iysptrx1, ixlb, ixrb
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: neq
      INTEGER, INTENT(IN), DIMENSION(4) :: range
      REAL, INTENT(IN), DIMENSION(neq+2) :: yl
      REAL, INTENT(OUT), DIMENSION(neq) :: yldot
      INTEGER :: ichunk, xc, yc, ii, locrange(4), ixpt, corechunk = 1
      REAL :: tick,tock!, tsfe, tsjf, ttotfe, ttotjf, tserial, tpara

        ! Patch cells that locally intercept the X-point
        corechunk = 0
        locrange = range
        do ixpt = 1, nxpt
            ! Chunk spanning left cut
            if (locrange(3).le.iysptrx1(ixpt)) then
                ! This is likely required by the potential equation:
                ! Investigate whether the equation can be collapsed
                ! to avoid unnexxesary repeats of the whole core region
                locrange(4) = max(locrange(4), iysptrx1(ixpt))
                if (     (locrange(1).lt.ixpt1(ixpt)) &
                &   .and.(locrange(2).gt.ixpt1(ixpt)) &
                &   .and.(locrange(2).le.ixpt2(ixpt)) ) then
                    locrange(2) = ixpt2(ixpt)+1
                    corechunk = 1
                ! Chunk spanning right cut
                elseif ( (locrange(1).lt.ixpt2(ixpt)) &
                &   .and.(locrange(2).gt.ixpt2(ixpt)) &
                &   .and.(locrange(1).gt.ixpt1(ixpt)) ) then
                    locrange(1) = ixpt1(ixpt)
                    corechunk = 2
                elseif ( (locrange(1).gt.ixpt1(ixpt)) &
                &   .and.(locrange(2).le.ixpt2(ixpt)) ) then
                    locrange(1) = ixpt1(ixpt)
                    locrange(2) = ixpt2(ixpt)+1
                    corechunk = 3
                ! Not sure why the core region is needed for legs: potential again?
                elseif (    (locrange(2).le.ixpt1(ixpt)) &
                &       .or.(locrange(1).gt.ixpt2(ixpt)) ) then
                    corechunk = -2
                endif
            endif
            ! Core-intersecting chunks
            if (corechunk.gt.0) then
                call OMPinitialize_ranges2d(locrange)
                call convsr_vo1 (xc, yc, yl)
                call convsr_vo2 (xc, yc, yl) 
                call convsr_aux1 (xc, yc)
            ! PF-intersecting chunks
            elseif (corechunk.eq.-2) then
                ! Inner leg
                locrange(1) = ixlb(ixpt)
                locrange(2) = ixpt1(ixpt)+1 
                call OMPinitialize_ranges2d(locrange)
                call convsr_vo1 (xc, yc, yl)
                call convsr_vo2 (xc, yc, yl) 
                call convsr_aux1 (xc, yc)

                ! Outer leg
                locrange(1) = ixpt2(ixpt)
                locrange(2) = ixrb(ixpt)+1
                call OMPinitialize_ranges2d(locrange)
                call convsr_vo1 (xc, yc, yl)
                call convsr_vo2 (xc, yc, yl) 
                call convsr_aux1 (xc, yc)
            end if
            ! TODO: not sure why whole core region is
            ! necessary - likely related to the potentials
            locrange(1) = ixlb(ixpt)
            locrange(2) = ixrb(ixpt)+1
            call OMPinitialize_ranges2d(locrange)
            call convsr_aux2 (xc, yc)
        end do


        ! Initialize local thread ranges
        call OMPinitialize_ranges2d(locrange)
        xc=-1; yc=-1

        ! Calculate plasma variables from yl
        call convsr_vo1 (xc, yc, yl)
        call convsr_vo2 (xc, yc, yl) 

        ! Calculate derived quantities frequently used
        call convsr_aux1 (xc, yc)
        call convsr_aux2 (xc, yc)
        ! Calculate the plasma diffusivities and drift velocities
        call calc_plasma_diffusivities
        call initialize_driftterms  
        call calc_driftterms1
        call calc_driftterms2
        ! Calculate currents and potential
        if(isphion+isphiofft .eq. 1) then
            call calc_currents
            call calc_fqp1
            call calc_fqp2
        endif
        ! Get friction and electron velocities
        call calc_friction(xc)
        call calc_elec_velocities
        ! Volumetric plasma and gas sinks & sources
        call calc_volumetric_sources(xc, yc)
        if (TimingPandfOn.gt.0) TimeNeudif=tick()
        call neudifpg
        if (TimingPandfOn.gt.0) TotTimeNeudif=TotTimeNeudif+tock(TimeNeudif)
        call calc_srcmod
        ! Calculate plasma & gas conductivities etc.
        call calc_plasma_viscosities
        call calc_plasma_heatconductivities
        call calc_plasma_equipartition
        call calc_gas_heatconductivities
        call engbalg
        call calc_plasma_transport

        call calc_fniycbo 
        call calc_plasma_momentum_coeffs
        call calc_plasma_momentum(xc, yc)
        call calc_plasma_energy(xc, yc)
        call calc_gas_energy
        call calc_feeiycbo 
        call calc_atom_seic 

        call calc_plasma_particle_residuals
        call calc_gas_continuity_residuals
        call calc_plasma_momentum_residuals
        call calc_plasma_energy_residuals(xc, yc)
        call calc_gas_energy_residuals
        if (isphion.eq.1) call calc_potential_residuals


        ! Calculate yldot vector
        call calc_rhs(yldot)
        ! Set boundary conditions directly in yldot
        call bouncon(neq, yldot)
        if (TimingPandfOn.gt.0) & 
        &      TotTimePandf=TotTimePandf+tock(TimePandf)

        ! ================ BEGIN OLD PANDF1 ===================

        ! If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and
        ! the ODEs need to be modified as original equations 
        ! are for d(nv)/dt, etc If isflxvar=2, variables are 
        ! ni,v,nTe,nTi,ng. Boundary equations and potential 
        ! equations are not reordered.
        if(isflxvar.ne.1 .and. isrscalf.eq.1) call rscalf(yl,yldot)

        if(dtreal < 1.e15) then
            if ( &
            &   (svrpkg=='nksol' .and. yl(neq+1)<0) &
            &   .or. svrpkg == 'petsc' &
            & ) then
                call add_timestep(neq, yl, yldot)
            endif   !if-test on svrpkg and ylcopy(neq+1)
        endif    !if-test on dtreal
    END SUBROUTINE OMPPandf


    SUBROUTINE OMPPandf_XPT(neq, yl, yldot, ranges)
      USE UEpar, ONLY: isphion, svrpkg, isphiofft
      USE PandfTiming, ONLY: TimePandf, TotTimePandf, TimingPandfOn, &
      &     TimeNeudif, TotTimeNeudif
      USE Ynorm, ONLY: isflxvar, isrscalf
      USE Time_dep_nwt, ONLY: dtreal
      USE Selec, ONLY: yinc, xrinc, xlinc, j3, i4, i8
      USE Indices_domain_dcl, ONLY: iymnbcl
      USE Dim, ONLY: nx

        USE Rhsides, ONLY: resco
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: neq
      INTEGER, INTENT(IN), DIMENSION(2,4) :: ranges
      REAL, INTENT(IN), DIMENSION(neq+2) :: yl
      REAL, INTENT(OUT), DIMENSION(neq) :: yldot
      INTEGER, DIMENSION(4) :: corerange
      INTEGER :: ichunk, xc, yc, ii
      REAL :: tick,tock!, tsfe, tsjf, ttotfe, ttotjf, tserial, tpara
        ! Initialize local thread ranges
        xc=-1; yc=-1
        corerange(1) = ranges(1,1)
        corerange(2) = ranges(2,2)
        corerange(3) = ranges(1,3)
        corerange(4) = ranges(2,4)
!        call initialize_ranges(xc, yc, xlinc, xrinc, yinc)
        ! Calculate plasma variables from yl
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call convsr_vo1 (xc, yc, yl)
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call convsr_vo2 (xc, yc, yl)
        end do
        ! Calculate derived quantities frequently used
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call convsr_aux1 (xc, yc)
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call convsr_aux2 (xc, yc)
        end do
        ! Calculate the plasma diffusivities and drift velocities
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_diffusivities
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call initialize_driftterms  
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_driftterms1
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_driftterms2
        end do
        if(isphion+isphiofft .eq. 1) then
            ! Calculate currents and potential
            ! The potential is required over the whole core area
            call OMPinitialize_ranges2d(corerange)
            call calc_currents
            do ii = 1, 2
                call OMPinitialize_ranges2d(ranges(ii,:))
                call calc_fqp1
            end do
            do ii = 1, 2
                call OMPinitialize_ranges2d(ranges(ii,:))
                call calc_fqp2
            end do
        endif
        ! Get friction and electron velocities
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_friction(xc)
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_elec_velocities
        end do
        ! Volumetric plasma and gas sinks & sources
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_volumetric_sources(xc, yc)
        end do
        if (TimingPandfOn.gt.0) TimeNeudif=tick()
            do ii = 1, 2
                call OMPinitialize_ranges2d(ranges(ii,:))
                call neudifpg
            end do
        if (TimingPandfOn.gt.0) TotTimeNeudif=TotTimeNeudif+tock(TimeNeudif)
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_srcmod
        end do
        ! Calculate plasma & gas conductivities etc.
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_viscosities
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_heatconductivities
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_equipartition
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_gas_heatconductivities
        end do
        ! TODO: Fix ENGBALG
        call OMPinitialize_ranges2d(corerange)
        call engbalg
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_transport
        end do
        call calc_fniycbo 
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_momentum_coeffs
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_momentum(xc, yc)
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_energy(xc, yc)
        end do
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_gas_energy
        end do
        call calc_feeiycbo 
        ! TODO: FIX ATOM_SEIC 
        call OMPinitialize_ranges2d(corerange)
        call calc_atom_seic 
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_plasma_particle_residuals
            call calc_gas_continuity_residuals
            call calc_plasma_momentum_residuals
            call calc_plasma_energy_residuals(xc, yc)
            call calc_gas_energy_residuals
            if (isphion.eq.1) call calc_potential_residuals
        end do

        ! Calculate yldot vector
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            call calc_rhs(yldot)
        end do
        ! TODO: Use separate parallel chunk to set core boundary
!        ! Set boundary conditions directly in yldot
!        do ii = 1, 2
!            call OMPinitialize_ranges2d(ranges(ii,:))
!            if (iymnbcl .ne. 0) call iwall_boundary(neq, yldot)
!        end do
!        call initialize_ranges(xc, yc, xlinc, xrinc, yinc)
!        if (iymnbcl .ne. 0) call iwall_boundary(neq, yldot)

        if (TimingPandfOn.gt.0) & 
        &      TotTimePandf=TotTimePandf+tock(TimePandf)

            ! ================ BEGIN OLD PANDF1 ===================

            ! If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and
            ! the ODEs need to be modified as original equations 
            ! are for d(nv)/dt, etc If isflxvar=2, variables are 
            ! ni,v,nTe,nTi,ng. Boundary equations and potential 
            ! equations are not reordered.
        do ii = 1, 2
            call OMPinitialize_ranges2d(ranges(ii,:))
            if(isflxvar.ne.1 .and. isrscalf.eq.1) call rscalf(yl,yldot)
        end do
        if(dtreal < 1.e15) then
            if ( &
            &   (svrpkg=='nksol' .and. yl(neq+1)<0) &
            &   .or. svrpkg == 'petsc' &
            & ) then
                do ii = 1, 2
                    call OMPinitialize_ranges2d(ranges(ii,:))
                    call add_timestep(neq, yl, yldot)
                end do
            endif   !if-test on svrpkg and ylcopy(neq+1)
        endif    !if-test on dtreal



    END SUBROUTINE OMPPandf_XPT




  SUBROUTINE CreateBin(ieqmin,ieqmax,ichunkmin,ichunkmax,ichunktot,Padding,iCenterBin,iLeftBin,iRightBin,inc)
    IMPLICIT NONE
    integer,intent(in):: ieqmin, ieqmax,Padding,ichunkmin,ichunkmax, ichunktot
    integer,intent(out):: iCenterbin(ichunktot),inc(ichunktot),iLeftBin(ichunktot),iRightBin(ichunktot)
    integer ::N,SizeBin,Nchunk, i
    N=ieqmax-ieqmin+1
    Nchunk=ichunkmax-ichunkmin+1

    if (N>Nchunk) then
        SizeBin=int((N/Nchunk))
    else
        SizeBin=1
    endif

    iLeftBin(ichunkmin) = ieqmin
    iRightBin(ichunkmin) = iLeftBin(ichunkmin)+SizeBin-1
    iCenterBin(ichunkmin) = int((iLeftBin(ichunkmin)+iRightBin(ichunkmin))/2)
    inc(ichunkmin) = max(iCenterBin(ichunkmin)-iLeftBin(ichunkmin), &
                iRightBin(ichunkmin)-iCenterBin(ichunkmin)) + padding
    if (ichunkmax.gt.ichunkmin) then
        do i=ichunkmin+1,ichunkmax-1
            iLeftBin(i) = iRightBin(i-1)+1
            iRightBin(i) = iLeftBin(i)+SizeBin-1
            iCenterBin(i) = int((iLeftBin(i)+iRightBin(i))/2)
            inc(i) = max(iCenterBin(i)-iLeftBin(i),iRightBin(i)-iCenterBin(i))+padding
        enddo
        iLeftBin(ichunkmax) = iRightBin(ichunkmax-1)+1
        iRightBin(ichunkmax) = ieqmax
        iCenterBin(ichunkmax) = int((iLeftBin(ichunkmax)+iRightBin(ichunkmax))/2)
        inc(ichunkmax) = max(iCenterBin(ichunkmax)-iLeftBin(ichunkmax), &
                iRightBin(ichunkmax)-iCenterBin(ichunkmax))+padding
    endif
    RETURN
  END SUBROUTINE CreateBin


  SUBROUTINE MakeChunksPandf1()
    USE Indexes, ONLY: igyl
    USE OMPPandf1Settings, ONLY: xpadding,ypadding,OMPPandf1Verbose
    USE OMPPandf1, ONLY: NchunksPandf1,yincchunk,xincchunk,iychunk, &
            ixchunk,Nychunks_old,Nxchunks_old,neq_old,ivchunk,&
            Nivchunk,Nxchunks,Nychunks,iymaxchunk,ixmaxchunk,iyminchunk,ixminchunk
    USE Lsode, ONLY: neq
    USE Dim, ONLY: nx,ny
    IMPLICIT NONE

    integer:: remakechunk,i,ii,ichunk,iv,ix,iy
    integer:: iyCenterBin(Nychunks),iyRightBin(Nychunks),iyLeftBin(Nychunks),incy(Nychunks)
    integer::ixCenterBin(Nxchunks),ixRightBin(Nxchunks),ixLeftBin(Nxchunks),incx(Nxchunks)
    remakechunk=0
    if ((Nxchunks.ne.Nxchunks_old).or.(Nychunks.ne.Nychunks_old)) then
        if (Nychunks.gt.1) then
            if (Nychunks.eq.ny) then
                iyLeftBin(1)=0
                iyRightBin(1)=1
                iyCenterBin(1)=1
                incy(1)=ypadding
                call CreateBin(   2,ny-1,2,Nychunks-1, Nychunks, ypadding, &
                                iyCenterBin, iyLeftBin, iyRightBin,&
                                incy &
                )
                iyLeftBin(Nychunks)=ny
                iyRightBin(Nychunks)=ny+1
                iyCenterBin(Nychunks)=ny
                incy(Nychunks)=ypadding
                if (OMPPandf1Verbose.gt.1) then
                    write(*,*) '----- Bins in y direction: ', Nychunks, ny+2
                    do iy=1,Nychunks
                        write(*,*) iyCenterBin(iy),iyLeftBin(iy),iyRightBin(iy),incy(iy)
                    enddo
                endif
            else
                call CreateBin(0,ny+1,1,Nychunks,Nychunks,ypadding,iyCenterBin, &
                    iyLeftBin,iyRightBin,incy)
                if (OMPPandf1Verbose.gt.1) then
                    write(*,*) '----- Bins in y direction: ', Nychunks, ny+2
                    do iy=1,Nychunks
                        write(*,*) iyCenterBin(iy),iyLeftBin(iy),iyRightBin(iy),incy(iy)
                    enddo
                endif
                ! now we check the first and last bins to check that ypadding is 3 if iyCenterBin=2
                if (iyCenterBin(1)==0) then
                    incy(1)=incy(1)+1
                endif
                if (iyCenterBin(Nychunks)==ny+1) then
                    incy(Nychunks)=incy(Nychunks)+1
                endif
            endif
        else
            iyCenterBin(1)=-1
            iyLeftBin(1)=0
            iyLeftBin(1)=ny+1
            incy(1)=0 !not used
        endif
        if (Nxchunks.gt.1) then
            call CreateBin(0,nx+1,1,Nxchunks,Nxchunks,xpadding,ixCenterBin,ixLeftBin,ixRightBin,incx)
        else
            ixCenterBin(1)=-1
            ixLeftBin(1)=0
            ixRightBin(1)=nx+1
            incx(1)=0 !not used
        endif
        ichunk=1
        ! Build NchunksPandf1
        do iy=1,Nychunks
            if (iy==1) then
                iychunk(ichunk)=iyCenterBin(iy)
                iyminchunk(ichunk)=iyLeftBin(iy)
                iymaxchunk(ichunk)=iyRightBin(iy)
                yincchunk(ichunk)=incy(iy)
                ixchunk(ichunk)=-1
                ixminchunk(ichunk)=0
                ixmaxchunk(ichunk)=nx+1
                xincchunk(ichunk)=0
                ichunk=ichunk+1
                CYCLE
            endif
            if (iy==ny+1) then
                iychunk(ichunk)=iyCenterBin(iy)
                iyminchunk(ichunk)=iyLeftBin(iy)
                iymaxchunk(ichunk)=iyRightBin(iy)
                yincchunk(ichunk)=incy(iy)
                ixchunk(ichunk)=-1
                ixminchunk(ichunk)=0
                ixmaxchunk(ichunk)=nx+1
                xincchunk(ichunk)=0
                ichunk=ichunk+1
                CYCLE
            endif
            do ix=1,Nxchunks
                iychunk(ichunk)=iyCenterBin(iy)
                iyminchunk(ichunk)=iyLeftBin(iy)
                iymaxchunk(ichunk)=iyRightBin(iy)
                yincchunk(ichunk)=incy(iy)
                ixchunk(ichunk)=ixCenterBin(ix)
                ixminchunk(ichunk)=ixLeftBin(ix)
                ixmaxchunk(ichunk)=ixRightBin(ix)
                xincchunk(ichunk)=incx(ix)
                ichunk=ichunk+1
            enddo
        enddo
        Nychunks_old=Nychunks
        Nxchunks_old=Nxchunks
        remakechunk=1
    endif
    if ((neq.ne.neq_old) .or. remakechunk.gt.0) then
        do i=1,NchunksPandf1
            ii=1
            do iv=1,neq
                if ((igyl(iv,2).le.iymaxchunk(i) .and. igyl(iv,2).ge.iyminchunk(i)) &
                .and. (igyl(iv,1).le.ixmaxchunk(i) .and. igyl(iv,1).ge.ixminchunk(i))) then
                    ivchunk(i,ii)=iv
                    Nivchunk(i)=ii
                    ii=ii+1
                endif
            enddo
        enddo
        neq_old=neq
        if (OMPPandf1Verbose.gt.1) then
            write(*,"('Nychunks = ',I3,'; Nxchunks = ',I3)") Nychunks,Nxchunks
            write(*,"('Nchunks:',I3)") NchunksPandf1
            do i=1,NchunksPandf1
                write(*,'("ichunk: ",I3," | iyc = [",I3,";",I3,";",I3,"] ; ixc = ",I3, &
                " || xinc = ", I3,"; yinc = ",I3," | iv = [",I5,";",I5,"]")') &
                i,iyminchunk(i),iychunk(i),iymaxchunk(i), ixchunk(i),xincchunk(i),  &
                yincchunk(i), ivchunk(i,1),ivchunk(i,Nivchunk(i))
            enddo
        endif
    endif
    RETURN
  END SUBROUTINE MakeChunksPandf1

#endif

