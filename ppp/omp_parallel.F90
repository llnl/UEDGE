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

    ! this is a placeholder for further parallelization but we need to implement handling of x-points discon. in ix indexing
    NchunksPandf1=Nychunks

    call gchange('OMPPandf1',0)

    if (OMPPandf1Verbose.gt.0) then
        write(*,*) OMPPandf1Stamp, ' NchunksPandf1 = ',NchunksPandf1
    endif
    RETURN
END SUBROUTINE InitOMPPandf1
!-------------------------------------------------------------------------------------------------

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
!-------------------------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------------------------
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
  SUBROUTINE OMPPandf1Rhs(neq,time,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,yincchunk,xincchunk, &
            iychunk,ixchunk,NchunksPandf1
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY:ny
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
    ylcopy(1:neq+1)=yl(1:neq+1)
    yldotcopy = 0
    yldottot = 0

    if (ijactot.gt.0) then

        Time1=omp_get_wtime()
        call MakeChunksPandf1
        call OmpCopyPointerup
          !$omp parallel do default(shared) schedule(dynamic,OMPPandf1LoopNchunk) &
          !$omp& private(iv,ichunk) firstprivate(ylcopy,yldotcopy) copyin(yinc,xlinc,xrinc) &
          !$omp& REDUCTION(+:yldottot)
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

                yinc=yinc_bkp
                xlinc=xlinc_bkp
                xrinc=xrinc_bkp
            enddo loopthread
          !$omp  END PARALLEL DO
        Time1=omp_get_wtime()-Time1

        OMPTimeParallelPandf1=Time1+OMPTimeParallelPandf1

        yldot(:neq) = yldottot

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
  END SUBROUTINE OMPPandf1Rhs

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


  SUBROUTINE chunk2d(nxl, nxu, nyl, nyu, chunks, nchunks)
  Use Lsode, only: neq
  IMPLICIT NONE
  integer, intent(in) :: nxl, nxu, nyl, nyu
  integer, intent(out) :: nchunks, chunks(neq,3)
    
    call chunk3d(nxl, nxu, nyl, nyu, 0, 0, chunks, nchunks)

  RETURN
  END SUBROUTINE chunk2d
  

  SUBROUTINE chunk3d(nxl, nxu, nyl, nyu, nzl, nzu, chunks, nchunks)
  Use Lsode, only: neq
  IMPLICIT NONE
  integer, intent(in) :: nxl, nxu, nyl, nyu, nzl, nzu
  integer, intent(out) :: nchunks, chunks(neq,3)
  integer ii, dx, dy, dz
    
    dx = nxu -nxl +1
    dy = nyu -nyl +1
    dz = nzu -nzl +1
    nchunks = dx*dy*dz

!   TODO OMP parallelize this loop?
    DO ii = 0, nchunks-1
        chunks(ii+1,3) = INT(ii/(dx*dy))
        chunks(ii+1,2) = nyl + INT( (ii -(chunks(ii+1,3)*dx*dy))/dx)
        chunks(ii+1,1) = nxl + MOD(ii, dx)
        chunks(ii+1,3) = chunks(ii+1,3) + nzl
    END DO
  RETURN
  END SUBROUTINE chunk3d
  


