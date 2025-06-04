MODULE CHUNK
CONTAINS
  SUBROUTINE Make2DChunks(Nxchunks, Nychunks, &
    &   N, Niv, ivchunk, rangechunk, Nxptchunks, Nivxpt, &
    &   ivxptchunk, rangexptchunk, Nmax, Nxptmax, Nivxptmax)!, ixychunk, Nixy)
    Use Dim, ONLY: nx, ny, nxpt
    Use Indexes, ONLY: igyl
    Use Lsode, ONLY: neq
    Use Xpoint_indices, ONLY: iysptrx1, ixpt1, ixpt2
    IMPLICIT NONE
    INTEGER, INTENT(OUT):: N, Nxptchunks(nxpt), Nxchunks, Nychunks
!    INTEGER, ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT):: ixychunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(OUT):: rangexptchunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT)::  ivxptchunk
    INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT)::  ivchunk, &
    &               rangechunk, Nivxpt
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)::  Niv !, Nixy
    INTEGER, INTENT(OUT) :: Nmax, Nxptmax, Nivxptmax
    integer:: ix, iy, nxi, nyi, ii, idx(2), idxl, ixpt, iix, iicut, iscut
    real:: dx, dy, dyxpt(2)
    integer, allocatable:: xlims(:,:), ylims(:,:), xlimsxpt(:,:,:,:), ylimsxpt(:,:,:,:)

    ! Ensure chunking setup is valid: between 1 and n(x/y) chunks
    Nxchunks = MAX(MIN(nx, Nxchunks),1)
    Nychunks = MAX(MIN(ny-1, Nychunks),1)
    ! Calculate the number of chunks needed maximum
    N = Nxchunks * Nychunks + 3*nxpt
    do ixpt = 1, nxpt
        Nxptchunks(ixpt) = MIN(Nxptchunks(ixpt),iysptrx1(ixpt)-1) ! Number of X-point chunks
        dyxpt(ixpt) = real(iysptrx1(ixpt)-1)/(Nxptchunks(ixpt))
    end do
    Nxptmax = MAXVAL(Nxptchunks)
    Nmax = neq
!    Nixymax = (nx+2)*(ny+2)
    ! Get the dx/dy per chunk
    dx = real(nx)/(Nxchunks)
    dy = real(ny-1)/(Nychunks)
    ! Allocate the necassary arrays
    allocate( &
    &       xlims(Nxchunks,2), &
    &       ylims(Nychunks,2), &
    &       xlimsxpt(nxpt, 2, MAXVAL(Nxptchunks),2), &
    &       ivchunk(N, Nmax), &
    &       rangechunk(N, 4), &
    &       Niv(N), &
    &       ylimsxpt(nxpt, 2, MAXVAL(Nxptchunks),2), &
    &       ivxptchunk(nxpt,MAXVAL(Nxptchunks),Nmax), & 
    &       Nivxpt(nxpt,MAXVAL(Nxptchunks)), &
    &       rangexptchunk(nxpt,2,MAXVAL(Nxptchunks),4))
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
        ylimsxpt(ixpt, 1, Nxptchunks(ixpt),2) = iysptrx1(ixpt)+1
        ylimsxpt(ixpt, 2, Nxptchunks(ixpt),2) = iysptrx1(ixpt)+1
        ! Set X-chunks
        do iix = 1, Nxptchunks(ixpt)
            ! Left cut
            xlimsxpt(ixpt, 1, iix, 1) = ixpt1(ixpt)-1
            xlimsxpt(ixpt, 1, iix, 2) = ixpt1(ixpt)+2
            ! Right cut
            xlimsxpt(ixpt, 2, iix, 1) = ixpt2(ixpt)-1
            xlimsxpt(ixpt, 2, iix, 2) = ixpt2(ixpt)+2
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
    ! Make special chunks for iy=0 boundary
    do ixpt = 1, nxpt
        rangechunk(1+3*(ixpt-1),1) =   0
        rangechunk(1+3*(ixpt-1),2) =   ixpt1(ixpt)-2
        rangechunk(1+3*(ixpt-1),3) =   0
        rangechunk(1+3*(ixpt-1),4) =   1
        rangechunk(2+3*(ixpt-1),1) =   ixpt1(ixpt)-1
        rangechunk(2+3*(ixpt-1),2) =   ixpt2(ixpt)+2
        rangechunk(2+3*(ixpt-1),3) =   0
        rangechunk(2+3*(ixpt-1),4) =   1
        rangechunk(3+3*(ixpt-1),1) =   ixpt2(ixpt)+3
        rangechunk(3+3*(ixpt-1),2) =   nx+1
        rangechunk(3+3*(ixpt-1),3) =   0
        rangechunk(3+3*(ixpt-1),4) =   1
    end do
    ! Ravel the ranges into chunks 
    do ix = 1, Nxchunks
        do iy = 1, Nychunks
            rangechunk(Nxchunks*(iy-1) + ix + 3*nxpt, 1) = xlims(ix,1)
            rangechunk(Nxchunks*(iy-1) + ix + 3*nxpt, 2) = xlims(ix,2)
            rangechunk(Nxchunks*(iy-1) + ix + 3*nxpt, 3) = ylims(iy,1)
            rangechunk(Nxchunks*(iy-1) + ix + 3*nxpt, 4) = ylims(iy,2)
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
    Nmax = MAXVAL(Niv)
    Nivxptmax = MAXVAL(Nivxpt)

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
    USE OMPPandf1, ONLY: Nxchunks,Nychunks, Nxptchunks, rangechunk, &
    &   rangexptchunk, Nchunks, Nxptchunks, ivchunk, ivxptchunk, &
    &   Nchunksmax, Nxptchunksmax, Nivxptchunk, Nivchunk, Nivxptchunksmax, &
    &   isnionxy_old, isngonxy_old, isuponxy_old, istionxy_old, &
    &   isteonxy_old, istgonxy_old, isphionxy_old, nisp_old, ngsp_old, &
    &   nx_old, ny_old, neq_old, Nxchunks_old, Nychunks_old, Nxptchunks_old
    USE Dim, ONLY: nx, ny, nisp, ngsp, nxpt
    USE Lsode, ONLY: neq
    USE Xpoint_indices, ONLY: iysptrx1
    USE UEpar, ONLY: isnionxy, isngonxy, isuponxy, istionxy, isteonxy, &
    &   istgonxy, isphionxy
    USE CHUNK 
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: rangexptchunk_tmp
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) ::  ivxptchunk_tmp
    INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  ivchunk_tmp, rangechunk_tmp, &
    &       Nivxptchunk_tmp
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  Nivchunk_tmp
    INTEGER :: rechunk, ixpt
    rechunk = 0

    if (nx_old.ne.nx) then
        rechunk = 1
    elseif (ny_old.ne.ny) then
        rechunk = 1
    elseif (nisp_old.ne.nisp) then
        rechunk = 1
    elseif (ngsp_old.ne.ngsp) then
        rechunk = 1
    elseif (MAXVAL(ABS(Nxptchunks_old-Nxptchunks)).ne.0) then
        rechunk = 1
    elseif (Nxchunks_old.ne.Nxchunks) then
        rechunk = 1
    elseif (Nychunks_old.ne.Nychunks) then
        rechunk = 1
    elseif (neq_old.ne.neq) then
        rechunk = 1
    elseif (MAXVAL(ABS(isnionxy_old-isnionxy)).ne.0) then
        rechunk = 1
    elseif (MAXVAL(ABS(isngonxy_old-isngonxy)).ne.0) then
        rechunk = 1
    elseif (MAXVAL(ABS(isuponxy_old-isuponxy)).ne.0) then
        rechunk = 1
    elseif (MAXVAL(ABS(istionxy_old-istionxy)).ne.0) then
        rechunk = 1
    elseif (MAXVAL(ABS(isteonxy_old-isteonxy)).ne.0) then
        rechunk = 1
    elseif (MAXVAL(ABS(istgonxy_old-istgonxy)).ne.0) then
        rechunk = 1
    elseif (MAXVAL(ABS(isphionxy_old-isphionxy)).ne.0) then
        rechunk = 1
    end if

    ! TODO: Add checks wheter anything has changed

    if (rechunk.ne.0) then
        call gchange('OMPPandf1',0)

        if (Nychunks.lt.0) then
            call xerrab('Nychunks<0. Nxchunks must be >=0')
        endif
        if (Nxchunks.lt.0) then
            call xerrab('Nxchunks<0. Nxchunks must be >=0')
        endif
        ! Set default to be approx 5x5 chunks
        if (Nychunks.eq.0) then
            Nychunks=int(ny/5)
        endif
        if (Nxchunks.eq.0) then
            Nxchunks=int(ny/5)
        endif
        do ixpt = 1, nxpt
            if (Nxptchunks(ixpt).eq.0) then
                Nxptchunks(ixpt) = int(iysptrx1(ixpt)/5)
            endif
        enddo


        call Make2DChunks(Nxchunks, Nychunks, Nchunks, Nivchunk_tmp, &
        &   ivchunk_tmp, rangechunk_tmp, Nxptchunks, Nivxptchunk_tmp, &
        &   ivxptchunk_tmp, rangexptchunk_tmp, Nchunksmax, Nxptchunksmax, &
        &   Nivxptchunksmax)

        call gchange('OMPPandf1',0)

        rangechunk=rangechunk_tmp; 
        ivchunk=ivchunk_tmp(:,:Nchunksmax)
        Nivchunk=Nivchunk_tmp
        rangexptchunk=rangexptchunk_tmp
        ivxptchunk=ivxptchunk_tmp(:,:,:Nivxptchunksmax)
        Nivxptchunk=Nivxptchunk_tmp;
    endif

    isnionxy_old = isnionxy; isngonxy_old = isngonxy; isuponxy_old = isuponxy
    istionxy_old = istionxy; isteonxy_old = isteonxy; istgonxy_old = istgonxy
    isphionxy_old = isphionxy;  nisp_old = nisp; ngsp_old = ngsp
    nx_old = nx; ny_old = ny; neq_old = neq; Nxptchunks_old = Nxptchunks
    Nxchunks_old = Nxchunks; Nychunks_old=Nychunks

    RETURN
END SUBROUTINE InitOMPPandf1


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

  SUBROUTINE OMPinitialize_ranges2D(limits)
    Use Selec
    Use Bcond, ONLY: xcnearrb, xcnearlb
    Use Dim, ONLY: nfsp, nhsp, nisp, nx, ny, nxpt
    Use Imprad, ONLY: isimpon
    Use Xpoint_indices, ONLY: ixlb, ixrb
    IMPLICIT NONE
    integer, intent(in):: limits(4)
    integer:: xs, xe, ys, ye, ixpt
    xs = limits(1)
    xe = limits(2)
    ys = limits(3)
    ye = limits(4)
      nfsp = nisp
      if (isimpon .eq. 3 .or. isimpon .eq. 4) nfsp = nhsp


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

    i1omp = xs
    i2omp = max(1, xs-0)
    i3omp = xs
    i4omp = max(0, xs-0)
    i5omp = min(nx, xe+0)
    i6omp = xe
    i7omp = xe
    i8omp = min(nx+1, xe+0)

    j1omp = max(ys, 0)
    j1pomp = max(ys-1,0)
    j1omp1 = max(ys-1, 0)
    j2omp = max(ys,1)
    j3omp = ys
    j4omp = max(ys,0)
    j5omp = min(ye+1, ny)
    j5pomp = min(ye+1, ny)
    j6omp = min(ye, ny+1)
    j6pomp = min(ye+1, ny+1)
    j7omp = ye
    j8omp = min(ye+1,ny+1)

    ixs = i2
    ixf = i5
    iys = j2omp
    iyf = j5omp
    ixs1 = i1
    ixf6 = i6
    iys1 = j1omp1
    iyf6 = j8omp


    xcnearrb = .FALSE.
    xcnearlb = .FALSE.
    do ixpt = 1, nxpt
        if (xs.eq.ixlb(ixpt)) xcnearrb = .TRUE.
        if (xe.eq.ixrb(ixpt)+1) xcnearlb = .TRUE.
    end do

  END SUBROUTINE OMPinitialize_ranges2D


  SUBROUTINE OMPPandf1Rhs(neq,time,yl,yldot)
! Recreates Pandf using parallel structure
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads, CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY: nxpt
    USE Grid, ONLY:ijactot
    USE Cdv, ONLY: comnfe
    USE Rhsides, ONLY: psorcxg, psorrg, psordis
    USE Time_dep_nwt, ONLY: dtreal, nufak
    USE Ynorm, ONLY: isflxvar, isrscalf
    USE PandfTiming, ONLY: TimePandf, TotTimePandf, TimingPandfOn, TimeNeudif, &
    &   TotTimeNeudif
    USE OMPTiming, ONLY: ParaTime, SerialTime
    USE ParallelEval, ONLY: ParallelPandfCall

    USE UEpar, ONLY: igas
    USE Xpoint_indices, ONLY: ixrb, ixlb
    USE Indices_domain_dcl, ONLY: iymnbcl,iymxbcl, ixmnbcl, ixmxbcl
    USE Share, ONLY: isudsym, geometry, islimon, ix_lim, nxc
    USE Bcond, ONLY: isfixlb
    USE OMPPandf1, ONLY: Nxptchunks, rangechunk, &
    &   rangexptchunk, Nchunks, ivchunk, ivxptchunk, &
    &   Nivxptchunk, Nivchunk
    USE OMPPandf1, ONLY: rangechunk, rangexptchunk, Nchunks, Nxptchunks, &
    &   ivchunk, ivxptchunk

    USE Indexes, ONLY: igyl
    IMPLICIT NONE
    INTEGER tid
    integer,intent(in)::neq
    real,intent(in)::yl(neq+1)
    real,intent(out)::yldot(neq)
    real,intent(in)::time
    real:: yldotcopy(neq), yldottot(neq)
    real yldotsave(neq),ylcopy(neq+2)
    INTEGER:: ichunk, xc, yc, ii, ixpt, iv
    real tick, tock
    external tick, tock

    ParallelPandfCall = 1
    ylcopy = yl
    yldotcopy = yldot
    yldottot = 0


    if (ijactot.gt.0) then
        !$OMP PARALLEL DO PRIVATE(ichunk) SCHEDULE(dynamic) &
        !$OMP &      FIRSTPRIVATE(ylcopy, yldotcopy) REDUCTION(+:yldottot)
            DO ichunk = 1, SUM(Nxptchunks) + Nchunks
                psorcxg = 0
                psorrg = 0
                psordis = 0
                if (ichunk .le. Nxptchunks(1)) then
                    call OMPPandf_XPT(neq, ylcopy, yldotcopy, rangexptchunk(1, :, ichunk,:))

                    do iv=1,Nivxptchunk(1, ichunk)
                        yldottot(ivxptchunk(1,ichunk,iv)) = yldottot(ivxptchunk(1,ichunk,iv)) &
                        &       + yldotcopy(ivxptchunk(1,ichunk,iv))
                    enddo
                elseif ((nxpt.eq.2).and.(ichunk.le.SUM(Nxptchunks))) then
                    call OMPPandf_XPT(neq, ylcopy, yldotcopy, rangexptchunk(2, :, ichunk-Nxptchunks(1),:))

                    do iv=1,Nivxptchunk(2, ichunk-Nxptchunks(1))
                        yldottot(ivxptchunk(2,ichunk-Nxptchunks(1),iv)) = &
                        &       yldottot(ivxptchunk(2,ichunk-Nxptchunks(1),iv)) &
                        &       + yldotcopy(ivxptchunk(2,ichunk-Nxptchunks(1),iv))
                    enddo
                else
                    call OMPPandf(neq, ylcopy, yldotcopy,rangechunk(ichunk-SUM(Nxptchunks),:))

                    do iv=1,Nivchunk(ichunk-SUM(Nxptchunks))
                        yldottot(ivchunk(ichunk-SUM(Nxptchunks),iv)) = &
                        &       yldottot(ivchunk(ichunk-SUM(Nxptchunks),iv)) &
                        &       + yldotcopy(ivchunk(ichunk-SUM(Nxptchunks),iv))
                    enddo
                endif
            END DO
            !$OMP END PARALLEL DO
        yldot = yldottot 

        if (CheckPandf1.gt.0) then
            call pandf (-1, -1, neq, time, ylcopy, yldotsave)
            call Compare(yldot,yldotsave,neq)
            write(*,'(a,i4)') "  Serial and parallel pandf are identical for nfe = ", comnfe
        endif
    else
       call pandf (-1,-1, neq, time, yl, yldot)
    endif
    ParallelPandfCall = 0
    RETURN

    END SUBROUTINE OMPPandf1Rhs


    SUBROUTINE OMPPandf(neq, yl, yldot, range)
      USE UEpar, ONLY: isphion, svrpkg, isphiofft
      USE PandfTiming, ONLY: TimePandf, TotTimePandf, TimingPandfOn, &
      &     TimeNeudif, TotTimeNeudif
      USE Ynorm, ONLY: isflxvar, isrscalf
      USE Time_dep_nwt, ONLY: dtreal
      USE Dim, ONLY: nxpt, nx, ny
      USE Xpoint_indices, ONLY: ixpt1, ixpt2, iysptrx1, ixlb, ixrb
      USE Bcond, ONLY: openbox
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: neq
      INTEGER, INTENT(IN), DIMENSION(4) :: range
      REAL, INTENT(IN), DIMENSION(neq+2) :: yl
      REAL, INTENT(OUT), DIMENSION(neq) :: yldot
      INTEGER :: ichunk, xc, yc, ii, ixpt, auxrange(4)
      REAL :: tick,tock!, tsfe, tsjf, ttotfe, ttotjf, tserial, tpara

        auxrange=range
        ! For some reason the right BC in bouncon only
        ! works robustly with the whole flux-tube. Since
        ! the boundary only calculates in the vicinity of the
        ! boundary, there is no additional cost to this.
        ! Just do it!
        do ixpt = 1, nxpt
            ! TODO: Figure out why this is needed??
            if ((auxrange(1).gt.ixpt1(ixpt)).and.(auxrange(2).le.ixpt2(ixpt))) then
                auxrange(1)=ixpt1(ixpt)+1
                auxrange(2)=ixpt2(ixpt)
            elseif (auxrange(2).le.ixpt1(ixpt)) then
                auxrange(1)=ixlb(ixpt)
                auxrange(2)=ixpt1(ixpt)
            endif
        end do



        ! Initialize local thread ranges
        xc=-1; yc=-1
        openbox = .true.
        
        call OMPinitialize_ranges2d(range)
        ! Calculate plasma variables from yl
        call convsr_vo1 (xc, yc, yl)
        call convsr_vo2 (xc, yc, yl) 
        ! Calculate derived quantities frequently used
        call convsr_aux1 (xc, yc)


        call OMPinitialize_ranges2d(auxrange)
        call convsr_aux2 (xc, yc)
        call OMPinitialize_ranges2d(range)


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

        call calc_plasma_momentum_coeffs
        call calc_plasma_momentum(xc, yc)
        call calc_plasma_energy(xc, yc)
        call calc_gas_energy
        call calc_atom_seic 

        if (isphion.eq.1) call calc_potential_residuals
        call calc_plasma_particle_residuals
        call calc_gas_continuity_residuals
        call calc_plasma_momentum_residuals
        call calc_plasma_energy_residuals(xc, yc)
        call calc_gas_energy_residuals

        ! Calculate yldot vector
        call calc_rhs(yldot)

        ! TODO: add more robust checks on ixrb: primary X-point
        ! RB may be in the middle of a chunk, rather than on the 
        ! boundary. Alternatively, the chunks should be defined to
        ! ensure a chunk on each LB/RB boundary

        ! Set boundary conditions directly in yldot
        do ixpt = 1, nxpt
            if (    (range(1).eq.ixlb(ixpt)) &
            &   .or.(range(2).eq.ixrb(ixpt)+1) &
            &   .or.(range(3).eq.0) &
            &   .or.(range(4).eq.ny+1) &
            & ) then
                call calc_fniycbo 
                call calc_feeiycbo 
                call bouncon(neq, yldot)
            endif
        end do
        !NOTE: This logic loop might be obsolete, due to internal checks
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
            do ii = 1, 2
                call OMPinitialize_ranges2d(ranges(ii,:))
                call calc_currents
            end do
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


    SUBROUTINE OMPilut (n,a,ja,ia,lfil,tol,alu,jlu,ju,iwk, &
    &               wu,wl,jr,jwl,jwu,ierr)
       IMPLICIT NONE
       integer n, ju0, j, ii, j1, j2, k, lenu, lenl, jj, nl, jrow
       integer jpos, len
       real tnorm, t, s, fact
       real a(*), alu(*), wu(n+1), wl(n), tol
       integer ja(*),ia(n+1),jlu(*),ju(n),jr(n), jwu(n), jwl(n), lfil, iwk, ierr
    !----------------------------------------------------------------------*
    !                      *** ILUT preconditioner ***                     *
    !                      ---------------------------                     *
    !      incomplete LU factorization with dual truncation mechanism      *
    !      VERSION 2 : sorting  done for both L and U.                     *
    !                                                                      *
    ! Bug Fix:  Version of 2-25-93.                                        *
    ! Modernizes loops: A. Holm, 6-2-25                                    *
    !----------------------------------------------------------------------*
    !---- coded by Youcef Saad May, 5, 1990. ------------------------------*
    !---- Dual drop-off strategy works as follows.                         *
    !                                                                      *
    !     1) Theresholding in L and U as set by tol. Any element whose size*
    !        is less than some tolerance (relative to the norm of current  *
    !        row in u) is dropped.                                         *
    !                                                                      *
    !     2) Keeping only the largest lfil elements in L and the largest   *
    !        lfil elements in U.                                           *
    !                                                                      *
    ! Flexibility: one can use tol=0 to get a strategy based on keeping the*
    ! largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n*
    ! will give the usual threshold strategy (however, fill-in is then     *
    ! impredictible).                                                      *
    !                                                                      *
    !----------------------------------------------------------------------*
    ! PARAMETERS
    !-----------
    !
    ! on entry:
    !==========
    ! n       = integer. The dimension of the matrix A.
    !
    ! a,ja,ia = matrix stored in Compressed Sparse Row format.
    !
    ! lfil    = integer. The fill-in parameter. Each row of L and
    !           each row of U will have a maximum of lfil elements
    !           in addition to the original number of nonzero elements.
    !           Thus storage can be determined beforehand.
    !           lfil must be .ge. 0.
    !
    ! iwk     = integer. The minimum (MAX??) length of arrays alu and jlu
    !
    ! On return:
    !===========
    !
    ! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
    !           the L and U factors together. The diagonal (stored in
    !           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
    !           contains the i-th row of L (excluding the diagonal entry=1)
    !           followed by the i-th row of U.
    !
    ! ju      = integer array of length n containing the pointers to
    !           the beginning of each row of U in the matrix alu,jlu.
    !
    ! ierr    = integer. Error message with the following meaning.
    !           ierr  = 0    --> successful return.
    !           ierr .gt. 0  --> zero pivot encountered at step number ierr.
    !           ierr  = -1   --> Error. input matrix may be wrong.
    !                            (The elimination process has generated a
    !                            row in L or U whose length is .gt.  n.)
    !           ierr  = -2   --> The matrix L overflows the array al.
    !           ierr  = -3   --> The matrix U overflows the array alu.
    !           ierr  = -4   --> Illegal value for lfil.
    !           ierr  = -5   --> zero row encountered.
    !
    ! work arrays:
    !=============
    ! jr,jwu,jwl 	  = integer work arrays of length n.
    ! wu, wl          = real work arrays of length n+1, and n resp.
    !
    ! Notes:
    ! ------
    ! A must have all nonzero diagonal elements.
    !-----------------------------------------------------------------------
    if (lfil .lt. 0) then
        ierr = -4
        RETURN
    endif 
    !-------------------------------
    ! initialize ju0 (points to next element to be added to alu,jlu)
    ! and pointer.
    !-----------------------------------------------------------------------
    !---- Initialize alu and jlu for diagnostic clarity --- TDR 1/24/00 ---
    alu(:iwk) = 0.
    jlu(:iwk) = 0.
    !------------------------------------------------------------------------
    ju0 = n+2
    jlu(1) = ju0
    !  integer double pointer array.
    jr(:j) = 0

    !-----------------------------------------------------------------------
    ! beginning of main loop.
    !-----------------------------------------------------------------------
    DO ii = 1, n
        j1 = ia(ii)
        j2 = ia(ii+1) - 1
        tnorm = 0.
        do k=j1,j2
            tnorm = tnorm+abs(a(k))
        end do
        if (tnorm .eq. 0.) then
            ierr = -5
            RETURN
        endif
        tnorm = tnorm/(j2-j1+1)
        ! unpack L-part and U-part of row of A in arrays wl, wu
        lenu = 1
        lenl = 0
        jwu(1) = ii
        wu(1) = 0.0
        jr(ii) = 1
        ! unpack lower and upper parts of row ii, in jwl-wl and
        ! jwu-wu compressed rows respectively. Ignore element if small
        do j = j1, j2
            k = ja(j)
            t = a(j)
            if ((abs(t) .lt. tol*tnorm).and.(k .ne. ii)) CYCLE
            if (k .lt. ii) then
                lenl = lenl+1
                jwl(lenl) = k
                wl(lenl) = t
                jr(k) = lenl
            else if (k .eq. ii) then
                wu(1) = t
            else
                lenu = lenu+1
                jwu(lenu) = k
                wu(lenu) = t
                jr(k) = lenu
            endif
        end do
        tnorm = tnorm/(j2-j1+1)
        !---------------------------------------------------------------
        jj = 1
        nl = 0
        ! eliminate previous rows
        DO WHILE( jj.le.lenl)
            !-------------------------------------------------------------------
            ! in order to do the elimination in the correct order we need to
            ! exchange the current row number with the one that has
            ! smallest column number, among jj,jj+1,...,lenl.
            !-------------------------------------------------------------------
            jrow = jwl(jj)
            k = jj
            ! determine smallest column index
            do j=jj+1,lenl
                if (jwl(j) .lt. jrow) then
                    jrow = jwl(j)
                    k = j
                endif
            end do
            ! exchange in jwl
            if (k .ne. jj) then
                j = jwl(jj)
                jwl(jj) = jwl(k)
                jwl(k) = j
                ! exchange in jr
                jr(jrow) = jj
                jr(j) = k
                ! exchange in wl
                s = wl(jj)
                wl(jj) = wl(k)
                wl(k) = s
            endif
            if (jrow .ge. ii) CYCLE
            ! get the multiplier for row to be eliminated: jrow
            fact = wl(jj)*alu(jrow)
            ! zero out element in row by setting jr(jrow) = 0
            jr(jrow) = 0
            if (abs(fact)*wu(n+2-jrow) .gt. tol*tnorm) then
                ! combine current row and row jrow
                do k = ju(jrow), jlu(jrow+1)-1
                    s = fact*alu(k)
                    j = jlu(k)
                    jpos = jr(j)
                    ! if fill-in element is small then disregard:
                    if (abs(s) .lt. tol*tnorm .and. jpos .eq. 0) CYCLE
                    if (j .ge. ii) then
                        ! dealing with upper part.
                        if (jpos .eq. 0) then
                            ! this is a fill-in element
                            lenu = lenu+1
                            if (lenu .gt. n) then
                                ierr = -1
                                RETURN
                            end if
                            jwu(lenu) = j
                            jr(j) = lenu
                            wu(lenu) = - s
                        else
                            ! no fill-in element --
                            wu(jpos) = wu(jpos) - s
                        endif
                    else
                        ! dealing with lower part.
                        if (jpos .eq. 0) then
                            ! this is a fill-in element
                            lenl = lenl+1
                            if (lenl .gt. n) then
                                ierr = -1
                                RETURN
                            end if
                            jwl(lenl) = j
                            jr(j) = lenl
                            wl(lenl) = - s
                        else
                            ! no fill-in element --
                            wl(jpos) = wl(jpos) - s
                        endif
                    endif
                end do 
                nl = nl+1
                wl(nl) = fact
                jwl(nl)  = jrow
            end if
            jj = jj+1
        END DO 
        ! update l-matrix
        len = min0(nl,lfil)
        call qsplit (wl,jwl,nl,len)
        do k=1, len
            if (ju0 .gt. iwk) then
                ierr = -2
                RETURN
            end if 
            alu(ju0) =  wl(k)
            jlu(ju0) =  jwl(k)
            ju0 = ju0+1
        end do
!
!     save pointer to beginning of row ii of U
!
        ju(ii) = ju0
!
!     reset double-pointer jr to zero (L-part - except first
!     jj-1 elements which have already been reset)
!
        do k= jj, lenl
            jr(jwl(k)) = 0
        end do

        len = min0(lenu,lfil)
        call qsplit (wu(2), jwu(2), lenu-1,len)
        ! update u-matrix
        t = abs(wu(1))
        if (len + ju0 .gt. iwk) then
            ierr = -3
            RETURN
        end if
        do k=2, len
            jlu(ju0) = jwu(k)
            alu(ju0) = wu(k)
            t = t + abs(wu(k) )
            ju0 = ju0+1
        end do
        ! save norm (in fact the average abs value) in wu (backwards)
        wu(n+2-ii) = t / (len+1)
        ! store inverse of diagonal element of u
        if (wu(1) .eq. 0.0) wu(1) = (0.0001 + tol)*tnorm
        alu(ii) = 1.0 / wu(1)
        ! update pointer to beginning of next row of U.
        jlu(ii+1) = ju0
        ! reset double-pointer jr to zero (U-part)
        do k=1, lenu
            jr(jwu(k)) = 0
        end do
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
    END DO 
    ierr = 0
    RETURN
    END SUBROUTINE OMPilut

#endif

