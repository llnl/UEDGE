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

    ixs = i2
    ixf = i5
    iys = j2
    iyf = j5
    ixs1 = i1
    ixf6 = i6
    iys1 = j1
    iyf6 = j6


  END SUBROUTINE OMPinitialize_ranges2D

  SUBROUTINE OMPinitialize_ranges(xc, yc)
    Use Selec 
    Use Bcond, ONLY: xcnearrb, xcnearlb
    IMPLICIT NONE
    integer, intent(in):: xc, yc

      i1 = xc
      i2 = xc
      i2p = xc
      i3 = xc
      i4 = xc
      i5 = xc
      i5m = xc
      i6 = xc
      i7 = xc
      i8 = xc
      j1 = yc
      j1p = yc
      j2 = yc
      j2p = yc
      j3 = yc
      j4 = yc
      j5 = yc
      j5m = yc
      j5p = yc
      j6 = yc
      j6p = yc
      j7 = yc
      j8 = yc

            ixs = xc
            ixf = xc
            iys = yc
            iyf = yc
            ixs1 = xc
            ixf6 = xc
            iys1 = yc
            iyf6 = yc

      xcnearrb = .TRUE.
      xcnearlb = .TRUE.
      END SUBROUTINE OMPinitialize_ranges


  SUBROUTINE OMPconvsr_vo1(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: ne, phi, ti, tg, te, nit, nz2, ng, nm, lng, ni
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

  END SUBROUTINE OMPconvsr_vo1



  SUBROUTINE OMPconvsr_vo2(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: up
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

  END SUBROUTINE OMPconvsr_vo2

  SUBROUTINE OMPconvsr_aux1(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: pg, pr, pre, pri, tg
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

  END SUBROUTINE OMPconvsr_aux1


  SUBROUTINE OMPconvsr_aux2(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: pgy0, tgy0, niy1, nity1, tiy1, phiy0, tiy0s, ney0, zeff, tiy0, pgy1, phiy0s, phiy1, &
    &    ngy0, priy0, phiv, niy0s, tey1, tiy1s, priy1, tgy1, niy0, tiv, priv, ney1, tev, &
    &    phiy1s, prtv, prev, tey0, znot, niy1s, nity0, ngy1
    USE Gradients, ONLY: gpry, gpix, ex, ey, gtex, gtiy, gpiy, gtey, gtix, gpex, gpondpotx, gprx, gpey
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

  END SUBROUTINE OMPconvsr_aux2


  SUBROUTINE OMPcalc_plasma_diffusivities(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Conduc, ONLY: dutm_use, difp_use, dif_use, kyi_use, trax_use, kxbohm, kxe_use, kxi_use, &
    &    kye_use, vy_use, kybohm, tray_use, dif2_use
    USE Compla, ONLY: betap
    USE OMPTiming
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

  END SUBROUTINE OMPcalc_plasma_diffusivities

  SUBROUTINE OMPinitialize_driftterms(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Conduc, ONLY: eta1, dclass_e, rtaue, dclass_i
    USE UEpar, ONLY: ctaui, ctaue
    USE Compla, ONLY: loglambda
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

  END SUBROUTINE OMPinitialize_driftterms

  SUBROUTINE OMPcalc_driftterms1(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: vydd, veycb, vyrd, vycf, vygp, vy, vycr, vyce, veycp, vycb, vycp
    USE Comtra, ONLY: coll_fe, diffusivwrk, coll_fi
    USE Conduc, ONLY: vyte_cft, vyti_cft, vy_cft
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

  END SUBROUTINE OMPcalc_driftterms1


  SUBROUTINE OMPcalc_driftterms2(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: v2dd, ve2cb, vytan, v2rd, ve2cd, vy, v2xgp, v2cd, vyavis, v2ce, q2cd, v2cb, v2
    USE Comflo, ONLY: fdiaxlb, fdiaxrb
    USE Xpoint_indices, ONLY: ixlb, ixrb
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

  END SUBROUTINE OMPcalc_driftterms2


  SUBROUTINE OMPcalc_currents(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Comflo, ONLY: fq2d, fmity, fqym, fqyai, fqyb, fqydti, fqymi, fqydt, fqyao, fqya, fqyae, fqygp, &
    &    fq2, fqyd, fqy
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

  END SUBROUTINE OMPcalc_currents

  SUBROUTINE OMPcalc_fqp(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Poten, ONLY: dphi_iy1
    USE Compla, ONLY: vy, netap, vyavis
    USE Bcond, ONLY: fqpsatrb, fqpsatlb
    USE Comflo, ONLY: fqp, fqx, fqxb
    USE Xpoint_indices, ONLY: ixlb, ixrb
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: ichunk, xc, yc, ii, jx
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: dphi_iy1_tmp(0:nx+1), vy_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      netap_tmp(0:nx+1,0:ny+1), fqpsatrb_tmp(0:ny+1,nxpt), &
    &      vyavis_tmp(0:nx+1,0:ny+1,1:nisp), fqp_tmp(0:nx+1,0:ny+1), &
    &      fqx_tmp(0:nx+1,0:ny+1), fqpsatlb_tmp(0:ny+1,nxpt), &
    &      fqxb_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    dphi_iy1_tmp=0.; vy_tmp=0.; netap_tmp=0.; fqpsatrb_tmp=0.; vyavis_tmp=0.
    fqp_tmp=0.; fqx_tmp=0.; fqpsatlb_tmp=0.; fqxb_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0


    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fqp_tmp,netap_tmp)
    DO ichunk = 1, NchunksPandf1
        
        
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
    fqp_tmp = 0;netap_tmp=0
    call OmpCopyPointernetap; call OmpCopyPointerfqp
    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:dphi_iy1_tmp, vy_tmp, netap_tmp, fqpsatrb_tmp, vyavis_tmp, fqp_tmp, fqx_tmp, &
    !$OMP &         fqpsatlb_tmp, fqxb_tmp)
    DO ichunk = 1, NchunksPandf1
        
        
        call OMPinitialize_ranges2d(rangechunk(ichunk,:))
        call calc_fqp2

        do ii = 1, Nixychunk(ichunk)
    
        if (yc .eq. 1) &
        & dphi_iy1_tmp(xc)=dphi_iy1_tmp(xc)+dphi_iy1(xc)

        do jx = 1, nxpt
            if (xc .eq. ixlb(jx)) &
            &   fqpsatlb_tmp(yc, jx)=fqpsatlb_tmp(yc,jx)+fqpsatlb(yc,jx)
            if (xc .eq. ixrb(jx)) &
            &   fqpsatrb_tmp(yc, jx)=fqpsatrb_tmp(yc,jx)+fqpsatrb(yc,jx)
        end do


         xc = ixychunk(ichunk,ii,1)
         yc = ixychunk(ichunk,ii,2)
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

  END SUBROUTINE OMPcalc_fqp


  SUBROUTINE OMPcalc_friction(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: upi, uz, uu, uup
    USE Cfric, ONLY: frici, frice
    USE Gradients, ONLY: ex
    USE UEpar, ONLY: cs
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

  END SUBROUTINE OMPcalc_friction

  SUBROUTINE OMPcalc_elec_velocities(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Compla, ONLY: vey, vex, upe
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

  END SUBROUTINE OMPcalc_elec_velocities


  SUBROUTINE OMPcalc_volumetric_sources(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OMPPandf1, ONLY: NchunksPandf1, Nixychunk, ixychunk, rangechunk
    USE OmpCopybbb
    USE Rhsides, ONLY: psorrgc, psorgc, psordis, psorbgg, psordisg, smoc, psorc, snic, seic, psor, &
    &    psorbgz, seec, msor, psorg, psorcxg, psorrg, msorxr, psorxr, psorxrc
    USE Conduc, ONLY: nuelg, nurc, nuvl, nuiz, nucxi, nucx, nueli, nuix
    USE Compla, ONLY: rtauy, rtaux, rtau
    USE OMPTiming, ONLY: ParaTime, SerialTime
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

  END SUBROUTINE OMPcalc_volumetric_sources

  SUBROUTINE OMPneudifpg(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Locflux, ONLY: floxg, conyg, conxg, floyg
    USE Compla, ONLY: vy, vyg, uu, v2, vygtan, uug, uuxg
    USE Comflo, ONLY: fngy4ord, fngy, fngxy, fngx, fngx4ord
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:floxg_tmp, vy_tmp, conyg_tmp, vyg_tmp, fngy4ord_tmp, &
    !$OMP &         fngy_tmp, fngxy_tmp, uu_tmp, v2_tmp, conxg_tmp, fngx_tmp, vygtan_tmp, &
    !$OMP &         floyg_tmp, uug_tmp, fngx4ord_tmp, uuxg_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call neudifpg

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

  END SUBROUTINE OMPneudifpg

  SUBROUTINE OMPcalc_srcmod(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Rhsides, ONLY: smoc, seic, seec, snic
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: smoc_tmp(0:nx+1,0:ny+1,1:nusp), seic_tmp(0:nx+1,0:ny+1), &
    &      seec_tmp(0:nx+1,0:ny+1), snic_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    smoc_tmp=0.; seic_tmp=0.; seec_tmp=0.; snic_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:smoc_tmp, seic_tmp, seec_tmp, snic_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_srcmod

        ! Update locally calculated variables
        smoc_tmp(xc,yc,:)=smoc_tmp(xc,yc,:)+smoc(xc,yc,:)
        seic_tmp(xc,yc)=seic_tmp(xc,yc)+seic(xc,yc)
        seec_tmp(xc,yc)=seec_tmp(xc,yc)+seec(xc,yc)
        snic_tmp(xc,yc,:)=snic_tmp(xc,yc,:)+snic(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    smoc=smoc_tmp; seic=seic_tmp; seec=seec_tmp; snic=snic_tmp
    call OmpCopyPointersmoc; call OmpCopyPointerseic
    call OmpCopyPointerseec; call OmpCopyPointersnic

  END SUBROUTINE OMPcalc_srcmod


  SUBROUTINE OMPcalc_plasma_viscosities(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE UEpar, ONLY: ctaui
    USE Conduc, ONLY: alfneo, visy, ktneo, visxneo, nuiistar, nuii, visx, k2neo
    USE Wkspace, ONLY: w
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:alfneo_tmp, visy_tmp, ctaui_tmp, &
    !$OMP &         ktneo_tmp, visxneo_tmp, nuiistar_tmp, nuii_tmp, w_tmp, visx_tmp, &
    !$OMP &         k2neo_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_plasma_viscosities

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

  END SUBROUTINE OMPcalc_plasma_viscosities

  SUBROUTINE OMPcalc_plasma_heatconductivities(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE UEpar, ONLY: ctaui, ctaue
    USE Conduc, ONLY: hcxineo, hcyij, hcxij, hcyi, hcyn, hcye, hcxi, hcxe, hcxn
    USE Wkspace, ONLY: w2, w1
    USE Comflo, ONLY: qipar
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:hcxineo_tmp, ctaui_tmp, w2_tmp, hcyij_tmp, hcxij_tmp, &
    !$OMP &         hcyi_tmp, ctaue_tmp, hcyn_tmp, hcye_tmp, qipar_tmp, hcxi_tmp, &
    !$OMP &         w1_tmp, hcxe_tmp, hcxn_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_plasma_heatconductivities

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

  END SUBROUTINE OMPcalc_plasma_heatconductivities

  SUBROUTINE OMPcalc_plasma_equipartition(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Conduc, ONLY: eqp
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: eqp_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    eqp_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:eqp_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_plasma_equipartition

        ! Update locally calculated variables
        eqp_tmp(xc,yc)=eqp_tmp(xc,yc)+eqp(xc,yc)
    END DO
    !$OMP  END PARALLEL DO
    eqp=eqp_tmp

    ! Update global variables
    call OmpCopyPointereqp

  END SUBROUTINE OMPcalc_plasma_equipartition

  SUBROUTINE OMPcalc_gas_heatconductivities(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Conduc, ONLY: hcxg, hcyg
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: hcxg_tmp(0:nx+1,0:ny+1,1:ngsp), hcyg_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    hcxg_tmp=0.; hcyg_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:hcxg_tmp, hcyg_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_gas_heatconductivities

        ! Update locally calculated variables
        hcxg_tmp(xc,yc,:)=hcxg_tmp(xc,yc,:)+hcxg(xc,yc,:)
        hcyg_tmp(xc,yc,:)=hcyg_tmp(xc,yc,:)+hcyg(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    hcxg=hcxg_tmp; hcyg=hcyg_tmp
    call OmpCopyPointerhcxg; call OmpCopyPointerhcyg

  END SUBROUTINE OMPcalc_gas_heatconductivities

  SUBROUTINE OMPengbalg(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Locflux, ONLY: floxge, floyge, conxge, conyge
    USE Comflo, ONLY: fegx, fegy, fegxy
    USE Rhsides, ONLY: segc
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:floxge_tmp, segc_tmp, floyge_tmp, conxge_tmp, & 
    !$OMP &      conyge_tmp, fegx_tmp, fegxy_tmp, fegy_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call engbalg

        ! Update locally calculated variables
        floxge_tmp(xc,yc,:)=floxge_tmp(xc,yc,:)+floxge(xc,yc,:)
        segc_tmp(xc,yc,:)=segc_tmp(xc,yc,:)+segc(xc,yc,:)
        fegx_tmp(xc,yc,:)=fegx_tmp(xc,yc,:)+fegx(xc,yc,:)
        fegxy_tmp(xc,yc,:)=fegxy_tmp(xc,yc,:)+fegxy(xc,yc,:)
        fegy_tmp(xc,yc,:)=fegy_tmp(xc,yc,:)+fegy(xc,yc,:)
        floyge_tmp(xc,yc,:)=floyge_tmp(xc,yc,:)+floyge(xc,yc,:)
        conxge_tmp(xc,yc,:)=conxge_tmp(xc,yc,:)+conxge(xc,yc,:)
        conyge_tmp(xc,yc,:)=conyge_tmp(xc,yc,:)+conyge(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    floxge=floxge_tmp; segc=segc_tmp; floyge=floyge_tmp; conxge=conxge_tmp
    conyge=conyge_tmp; fegx=fegx_tmp; fegy=fegy_tmp;fegxy=fegxy_tmp
    call OmpCopyPointerfloxge; call OmpCopyPointersegc
    call OmpCopyPointerfloyge; call OmpCopyPointerconxge
    call OmpCopyPointerconyge; call OmpCopyPointerfegx
    call OmpCopyPointerfegxy; call OmpCopyPointerfegy
  END SUBROUTINE OMPengbalg

  SUBROUTINE OMPcalc_plasma_transport(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Comflo, ONLY: fniy, fniy4ord, fnixcb, fniycbo, fniycb, fnix
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: fniy_tmp(0:nx+1,0:ny+1,1:nisp), fniy4ord_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      fnixcb_tmp(0:nx+1,0:ny+1,1:nisp), fniycbo_tmp(0:nx+1,1:nisp), &
    &      fniycb_tmp(0:nx+1,0:ny+1,1:nisp), fnix_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    fniy_tmp=0.; fniy4ord_tmp=0.; fnixcb_tmp=0.; fniycbo_tmp=0.; fniycb_tmp=0.
    fnix_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fniy_tmp, fniy4ord_tmp, fnixcb_tmp, fniycbo_tmp, fniycb_tmp, fnix_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_plasma_transport

        ! Update locally calculated variables
        fniy_tmp(xc,yc,:)=fniy_tmp(xc,yc,:)+fniy(xc,yc,:)
        fniy4ord_tmp(xc,yc,:)=fniy4ord_tmp(xc,yc,:)+fniy4ord(xc,yc,:)
        fnixcb_tmp(xc,yc,:)=fnixcb_tmp(xc,yc,:)+fnixcb(xc,yc,:)
        fniycbo_tmp(xc,yc)=fniycbo_tmp(xc,yc)+fniycbo(xc,yc)
        fniycb_tmp(xc,yc,:)=fniycb_tmp(xc,yc,:)+fniycb(xc,yc,:)
        fnix_tmp(xc,yc,:)=fnix_tmp(xc,yc,:)+fnix(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    fniy=fniy_tmp; fniy4ord=fniy4ord_tmp; fnixcb=fnixcb_tmp
    fniycbo=fniycbo_tmp; fniycb=fniycb_tmp; fnix=fnix_tmp
    call OmpCopyPointerfniy; call OmpCopyPointerfniy4ord
    call OmpCopyPointerfnixcb; call OmpCopyPointerfniycbo
    call OmpCopyPointerfniycb; call OmpCopyPointerfnix

  END SUBROUTINE OMPcalc_plasma_transport

  SUBROUTINE OMPcalc_plasma_momentum_coeffs(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Locflux, ONLY: conx, floy, flox, cony
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, iusp, ixpt
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: conx_tmp(0:nx+1,0:ny+1,1:nusp), floy_tmp(0:nx+1,0:ny+1,1:nusp), &
    &      flox_tmp(0:nx+1,0:ny+1,1:nusp), cony_tmp(0:nx+1,0:ny+1,1:nusp)

    ! Initialize arrays to zero
    conx_tmp=0.; floy_tmp=0.; flox_tmp=0.; cony_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:conx_tmp, floy_tmp, flox_tmp, cony_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_plasma_momentum_coeffs
        ! Update locally calculated variables
        conx_tmp(xc,yc,:)=conx_tmp(xc,yc,:)+conx(xc,yc,:)
        floy_tmp(xc,yc,:)=floy_tmp(xc,yc,:)+floy(xc,yc,:)
        flox_tmp(xc,yc,:)=flox_tmp(xc,yc,:)+flox(xc,yc,:)
        cony_tmp(xc,yc,:)=cony_tmp(xc,yc,:)+cony(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    conx=conx_tmp; floy=floy_tmp;flox=flox_tmp;cony=cony_tmp
    call OmpCopyPointerconx; call OmpCopyPointerfloy
    call OmpCopyPointerflox; call OmpCopyPointercony

  END SUBROUTINE OMPcalc_plasma_momentum_coeffs


  SUBROUTINE OMPcalc_plasma_momentum(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Compla, ONLY: fmivxpt, fmihxpt, vyvxpt, nixpt, vyhxpt, visyxpt, upxpt, up
    USE Comflo, ONLY: fmixy, fmix, fmiy
    USE Rhsides, ONLY: smoc
    USE UEpar, ONLY: methu, isupon
    USE Locflux, ONLY: flox, floy, conx, cony
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, iusp, ixpt, ifld
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fmivxpt_tmp, fmihxpt_tmp, fmixy_tmp, &
    !$OMP &         vyvxpt_tmp, nixpt_tmp, vyhxpt_tmp, visyxpt_tmp, upxpt_tmp, smoc_tmp, &
    !$OMP &         fmix_tmp, fmiy_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
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

        ! Update locally calculated variables
        fmixy_tmp(xc,yc,:)=fmixy_tmp(xc,yc,:)+fmixy(xc,yc,:)
        fmix_tmp(xc,yc,:)=fmix_tmp(xc,yc,:)+fmix(xc,yc,:)
        fmiy_tmp(xc,yc,:)=fmiy_tmp(xc,yc,:)+fmiy(xc,yc,:)
        smoc_tmp(xc,yc,:)=smoc_tmp(xc,yc,:)+smoc(xc,yc,:)
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

  END SUBROUTINE OMPcalc_plasma_momentum


  SUBROUTINE OMPcalc_plasma_energy(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nzspmx, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
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
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

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
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_plasma_energy(xc,yc)

        if (yc .eq. 0) then
            feeycbo_tmp(xc)=feeycbo_tmp(xc)+feeycbo(xc)
            feiycbo_tmp(xc)=feiycbo_tmp(xc)+feiycbo(xc)
        end if

        nzloc_tmp=nzloc_tmp+nzloc
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

  END SUBROUTINE OMPcalc_plasma_energy


  SUBROUTINE OMPcalc_gas_energy(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Rhsides, ONLY: seic, eiamoldiss
    USE Conduc, ONLY: eqpg
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: seic_tmp(0:nx+1,0:ny+1), eqpg_tmp(0:nx+1,0:ny+1,ngsp), &
    &      eiamoldiss_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    seic_tmp=0.; eqpg_tmp=0.; eiamoldiss_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:seic_tmp, eqpg_tmp, eiamoldiss_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_gas_energy

        ! Update locally calculated variables
        seic_tmp(xc,yc)=seic_tmp(xc,yc)+seic(xc,yc)
        eqpg_tmp(xc,yc,:)=eqpg_tmp(xc,yc,:)+eqpg(xc,yc,:)
        eiamoldiss_tmp(xc,yc,:)=eiamoldiss_tmp(xc,yc,:)+eiamoldiss(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    seic=seic_tmp; eqpg=eqpg_tmp; eiamoldiss=eiamoldiss_tmp
    call OmpCopyPointerseic; call OmpCopyPointereqpg
    call OmpCopyPointereiamoldiss

  END SUBROUTINE OMPcalc_gas_energy

  SUBROUTINE OMPcalc_plasma_particle_residuals(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE MCN_sources, ONLY: sng_ue
    USE MCN_dim, ONLY: nfl
    USE Rhsides, ONLY: resco
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: sng_ue_tmp(0:nx+1,0:ny+1,1:nfl), resco_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    sng_ue_tmp=0.; resco_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:sng_ue_tmp, resco_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_plasma_particle_residuals

        ! Update locally calculated variables
        sng_ue_tmp(xc,yc,:)=sng_ue_tmp(xc,yc,:)+sng_ue(xc,yc,:)
        resco_tmp(xc,yc,:)=resco_tmp(xc,yc,:)+resco(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    sng_ue=sng_ue_tmp; resco=resco_tmp
    call OmpCopyPointersng_ue; call OmpCopyPointerresco

  END SUBROUTINE OMPcalc_plasma_particle_residuals

  SUBROUTINE OMPcalc_gas_continuity_residuals(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Rhsides, ONLY: resng
    USE MCN_sources, ONLY: sng_ue
    USE MCN_dim, ONLY: nfl
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: resng_tmp(0:nx+1,0:ny+1,1:ngsp), sng_ue_tmp(0:nx+1,0:ny+1,1:nfl)

    ! Initialize arrays to zero
    resng_tmp=0.; sng_ue_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:resng_tmp, sng_ue_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_gas_continuity_residuals

        ! Update locally calculated variables
        resng_tmp(xc,yc,:)=resng_tmp(xc,yc,:)+resng(xc,yc,:)
        sng_ue_tmp(xc,yc,:)=sng_ue_tmp(xc,yc,:)+sng_ue(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    resng=resng_tmp; sng_ue=sng_ue_tmp
    call OmpCopyPointerresng; call OmpCopyPointersng_ue

  END SUBROUTINE OMPcalc_gas_continuity_residuals

  SUBROUTINE OMPcalc_plasma_momentum_residuals(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Wkspace, ONLY: w0, w2
    USE Cfric, ONLY: fricnrl
    USE Rhsides, ONLY: resmo
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: w0_tmp(0:nx+1,0:ny+1), fricnrl_tmp(0:nx+1,0:ny+1,nusp), &
    &      resmo_tmp(0:nx+1,0:ny+1,1:nusp), w2_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    w0_tmp=0.; fricnrl_tmp=0.; resmo_tmp=0.; w2_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:w0_tmp, fricnrl_tmp, resmo_tmp, w2_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_plasma_momentum_residuals

        ! Update locally calculated variables
        w0_tmp(xc,yc)=w0_tmp(xc,yc)+w0(xc,yc)
        fricnrl_tmp(xc,yc,:)=fricnrl_tmp(xc,yc,:)+fricnrl(xc,yc,:)
        resmo_tmp(xc,yc,:)=resmo_tmp(xc,yc,:)+resmo(xc,yc,:)
        w2_tmp(xc,yc)=w2_tmp(xc,yc)+w2(xc,yc)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    w0=w0_tmp; fricnrl=fricnrl_tmp; resmo=resmo_tmp; w2=w2_tmp
    call OmpCopyPointerw0; call OmpCopyPointerfricnrl
    call OmpCopyPointerresmo; call OmpCopyPointerw2

  END SUBROUTINE OMPcalc_plasma_momentum_residuals

  SUBROUTINE OMPcalc_gas_energy_residuals(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Rhsides, ONLY: reseg
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: reseg_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    reseg_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:reseg_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_gas_energy_residuals

        ! Update locally calculated variables
        reseg_tmp(xc,yc,:)=reseg_tmp(xc,yc,:)+reseg(xc,yc,:)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    reseg=reseg_tmp
    call OmpCopyPointerreseg

  END SUBROUTINE OMPcalc_gas_energy_residuals

  SUBROUTINE OMPcalc_plasma_energy_residuals(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt, nusp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Wkspace, ONLY: w0
    USE Rhsides, ONLY: resei, wjdote, reseg, resee, wvh
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: resei_tmp(0:nx+1,0:ny+1), &
    &      reseg_tmp(0:nx+1,0:ny+1), &
    &      resee_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    resei_tmp=0.; reseg_tmp=0.; resee_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:resei_tmp, reseg_tmp, resee_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_plasma_energy_residuals(xc,yc)

        ! Update locally calculated variables
        resei_tmp(xc,yc)=resei_tmp(xc,yc)+resei(xc,yc)
        reseg_tmp(xc,yc)=reseg_tmp(xc,yc)+reseg(xc,yc,1)
        resee_tmp(xc,yc)=resee_tmp(xc,yc)+resee(xc,yc)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    resei=resei_tmp; reseg(:,:,1)=reseg_tmp
    resee=resee_tmp; 
    call OmpCopyPointerresei
    call OmpCopyPointerreseg
    call OmpCopyPointerresee

  END SUBROUTINE OMPcalc_plasma_energy_residuals


  SUBROUTINE OMPcalc_potential_residuals(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Rhsides, ONLY: resphi
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: resphi_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    resphi_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:resphi_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call OMPinitialize_ranges(xc, yc)
        call calc_potential_residuals

        ! Update locally calculated variables
        resphi_tmp(xc,yc)=resphi_tmp(xc,yc)+resphi(xc,yc)
    END DO
    !$OMP  END PARALLEL DO

    ! Update global variables
    resphi=resphi_tmp
    call OmpCopyPointerresphi

  END SUBROUTINE OMPcalc_potential_residuals




  SUBROUTINE OMPbouncon(neq,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,yincchunk,xincchunk, &
            iychunk,ixchunk,NchunksPandf1, rangechunk, Nxchunks, Nychunks
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY:nx,ny,nxpt, nusp
    USE Math_problem_size, ONLY: numvar
    USE Xpoint_indices, ONLY: ixrb, ixlb, ixpt2
    USE UEpar, ONLY: igas
    USE Indices_domain_dcl, ONLY: iymnbcl,iymxbcl, ixmnbcl, ixmxbcl
    USE Share, ONLY: isudsym, geometry, islimon, ix_lim, nxc
    USE Bcond, ONLY: isfixlb
    USE Indexes, ONLY: idxu
    
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
        Nxchunks_old = Nxchunks
        Nychunks_old = Nychunks
        Nxchunks = 1
        Nychunks = 1
        call Make2DChunks
        

        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, NchunksPandf1
            ! TODO: Figure out more generalized chunking routines for BCs
            call OMPinitialize_ranges2D(rangechunk(ichunk,:))
            call bouncon(neq, yldotcopy)

            do ii=1,Nivchunk(ichunk)
                yldottot(ivchunk(ichunk,ii)) = yldottot(ivchunk(ichunk,ii)) + yldotcopy(ivchunk(ichunk,ii))
            end do

        END DO
        !$OMP END PARALLEL DO
        yldot(1:neq) = yldottot(1:neq)

        Nxchunks = Nxchunks_old
        Nychunks = Nychunks_old
        call Make2DChunks

    RETURN
  END SUBROUTINE OMPbouncon




  SUBROUTINE OMPcalc_rhs(neq,yl,yldot)
    USE omp_lib
    USE OmpCopybbb
    USE ParallelSettings, ONLY: Nthreads,CheckPandf1
    USE OMPPandf1Settings, ONLY: OMPTimeParallelPandf1,OMPTimeSerialPandf1, &
            OMPPandf1Stamp,OMPPandf1Verbose,OMPPandf1Debug
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,yincchunk,xincchunk, &
            iychunk,ixchunk,NchunksPandf1
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

        call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, Nchunks
            xc = chunks(ichunk,1)
            yc = chunks(ichunk,2)

            call OMPinitialize_ranges(xc, yc)
            call calc_rhs(yldotcopy)
            
            do ii = 1, numvar
                yldottot((ichunk-1)*numvar + ii) = yldottot((ichunk-1)*numvar + ii) &
                &       + yldotcopy((ichunk-1)*numvar + ii)
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
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,yincchunk,xincchunk, &
            iychunk,ixchunk,NchunksPandf1
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

        call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, Nchunks
            xc = chunks(ichunk,1)
            yc = chunks(ichunk,2)

            call OMPinitialize_ranges(xc, yc)
            if ((xc.gt.0).and.(xc.lt.nx+1).and.(yc.gt.0).and.(yc.lt.ny+1)) then
                call rscalf(ylcopy,yldotcopy)
            
            end if

            do ii = 1, numvar
                yldottot((ichunk-1)*numvar + ii) = yldottot((ichunk-1)*numvar + ii) &
                &       + (yldotcopy((ichunk-1)*numvar + ii))
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
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,yincchunk,xincchunk, &
            iychunk,ixchunk,NchunksPandf1
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

        call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

        !$OMP    PARALLEL DO &
        !$OMP &      default(shared) &
        !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
        !$OMP &      private(ichunk,xc,yc) &
        !$OMP &      firstprivate(ylcopy, yldotcopy) &
        !$OMP &      REDUCTION(+:yldottot)
        DO ichunk = 1, Nchunks
            xc = chunks(ichunk,1)
            yc = chunks(ichunk,2)

            call OMPinitialize_ranges(xc, yc)
            if (isbcwdt .eq. 1) then
                call add_timestep(neq, ylcopy, yldotcopy)
            else
                if ((xc.gt.0).and.(xc.lt.nx+1).and.(yc.gt.0).and.(yc.lt.ny+1)) then
                    call add_timestep(neq, ylcopy, yldotcopy)
                end if
            end if 
            
            do ii = 1, numvar
                yldottot((ichunk-1)*numvar + ii) = yldottot((ichunk-1)*numvar + ii) &
                &       + yldotcopy((ichunk-1)*numvar + ii)
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
    USE OMPPandf1, ONLY: Nivchunk,ivchunk,yincchunk,xincchunk, &
            iychunk,ixchunk,NchunksPandf1
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE Dim, ONLY: nx, ny
    USE Selec, ONLY:yinc,xrinc,xlinc, i1, i6, j1, j6
    USE Grid, ONLY:ijactot
    USE Cdv, ONLY: comnfe
    USE Rhsides, ONLY: psorcxg, psorrg, psordis
    USE Time_dep_nwt, ONLY: dtreal, nufak
    USE Ynorm, ONLY: isflxvar, isrscalf
    USE MCN_sources, ONLY: ismcnon
    USE UEpar, ONLY: isphion, svrpkg, isphiofft
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
    USE OMPPandf1, ONLY: Nchunks, chunks
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
    INTEGER:: ichunk, xc, yc, ii
    real tmp_prad(0:nx+1, 0:ny+1)
    real tick,tock, tsfe, tsjf, ttotfe, ttotjf, tserial, tpara
    external tick, tock

    real yldot1(1:neq), yldot2(1:neq)

    ParallelPandfCall = 1
    xc=-1; yc=-1

!        tpara = tick()
!        ParaTime = ParaTime + tock(tpara)

!        tserial = tick()
!        call initialize_ranges(xc, yc, xlinc, xrinc, yinc)
!        SerialTime = SerialTime + tock(tserial)


    if (ijactot.gt.0) then
        Time1=omp_get_wtime()
!        call MakeChunksPandf1
        call Make2DChunks
        call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

        call OMPconvsr_vo1 (neq, yl, yldot) 
        call OMPconvsr_vo2 (neq, yl, yldot) 
        call OMPconvsr_aux1 (neq, yl, yldot) 
        call OMPconvsr_aux2 (neq, yl, yldot) 

        call OMPcalc_plasma_diffusivities (neq, yl, yldot) 
        call OMPinitialize_driftterms (neq, yl, yldot) 
        call OMPcalc_driftterms1(neq, yl, yldot)
        call OMPcalc_driftterms2(neq, yl, yldot)
        if(isphion+isphiofft .eq. 1) then
            call OMPcalc_currents(neq, yl, yldot)
            call OMPcalc_fqp(neq, yl, yldot)
        endif
        call OMPcalc_friction(neq, yl, yldot)
        call OMPcalc_elec_velocities(neq, yl, yldot)


!        tpara = tick()
        call OMPcalc_volumetric_sources(neq, yl, yldot)
!        ParaTime = ParaTime + tock(tpara)

!        call initialize_ranges(xc, yc, xlinc, xrinc, yinc)
!        tserial = tick()
!        call calc_volumetric_sources(xc, yc)
!        SerialTime = SerialTime + tock(tserial)

        if (TimingPandfOn.gt.0) TimeNeudif=tick()
        call OMPneudifpg(neq, yl, yldot)
        if (TimingPandfOn.gt.0) TotTimeNeudif=TotTimeNeudif+tock(TimeNeudif)
        call OMPcalc_srcmod(neq, yl, yldot)
        call OMPcalc_plasma_viscosities(neq, yl, yldot)
        call OMPcalc_plasma_heatconductivities(neq, yl, yldot)
        call OMPcalc_plasma_equipartition(neq, yl, yldot)
        call OMPcalc_gas_heatconductivities(neq, yl, yldot)
        call OMPengbalg(neq, yl, yldot)
        call OMPcalc_plasma_transport(neq, yl, yldot)
        call calc_fniycbo ! Nothing much to parallelize here, just do serial
        call OmpCopyPointerfniycbo
        call OMPcalc_plasma_momentum_coeffs(neq, yl, yldot)
        call OMPcalc_plasma_momentum(neq, yl, yldot) 
        call OMPcalc_plasma_energy(neq, yl, yldot)
        call OMPcalc_gas_energy(neq, yl, yldot)
        call calc_feeiycbo ! Nothing much to parallelize here, just do serial
        call OmpCopyPointerfeeycbo
        call OmpCopyPointerfeiycbo
        call OMPcalc_plasma_particle_residuals(neq, yl, yldot)
        call OMPcalc_gas_continuity_residuals(neq, yl, yldot)
        call OMPcalc_plasma_momentum_residuals(neq, yl, yldot)
        call OMPcalc_gas_energy_residuals(neq, yl, yldot)
        call calc_atom_seic ! Nothing much to parallelize here, just do serial
        call OmpCopyPointerseic
        call OMPcalc_plasma_energy_residuals(neq, yl, yldot)
        if (isphion.eq.1) call OMPcalc_potential_residuals(neq, yl, yldot)
        call OMPcalc_rhs(neq, yl, yldot)

        call OMPbouncon(neq, yl, yldot)

!        call initialize_ranges(xc, yc, xlinc, xrinc, yinc)
!        call bouncon(neq, yldot)


!        call OMPbouncon(neq, yl, yldot)

        if (TimingPandfOn.gt.0) & 
        &      TotTimePandf=TotTimePandf+tock(TimePandf)

        ! ================ BEGIN OLD PANDF1 ===================

        ! If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and
        ! the ODEs need to be modified as original equations 
        ! are for d(nv)/dt, etc If isflxvar=2, variables are 
        ! ni,v,nTe,nTi,ng. Boundary equations and potential 
        ! equations are not reordered.
        if(isflxvar.ne.1 .and. isrscalf.eq.1) call OMPrscalf(neq, yl,yldot)

        if(dtreal < 1.e15) then
            if ( &
            &   (svrpkg=='nksol' .and. yl(neq+1)<0) &
            &   .or. svrpkg == 'petsc' &
            & ) then
                yldot1 = yldot(1:neq)
                call OMPadd_timestep(neq, yl, yldot)
            endif   !if-test on svrpkg and ylcopy(neq+1)
        endif    !if-test on dtreal



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
  

  SUBROUTINE Make2DChunks
    Use Dim, ONLY: nx, ny
    Use Indexes, ONLY: igyl
    Use Lsode, ONLY: neq
    Use OMPPandf1, ONLY: Nxchunks, Nychunks, NchunksPandf1, Nivchunk, &
    &   ivchunk, Nchunksmax, rangechunk, ixychunk, Nixychunk, Nixychunksmax

    IMPLICIT NONE
    integer:: ix, iy, nxi, nyi, ii, idx(2), idxl
    real:: dx, dy
    integer, allocatable:: xlims(:,:), ylims(:,:), Nivchunks(:,:), Niv(:)
    Nxchunks = MAX(MIN(nx, Nxchunks),1)
    Nychunks = MAX(MIN(ny, Nychunks),1)
!    if (Nxchunks .gt. nx) Nxchunks = nx ! Limit nx to ensure 
    ! boundary chunks include guard cells
!    if (Nychunks .gt. ny) Nychunks = ny ! is the same necessary for Y-chunks?
    NchunksPandf1 = Nxchunks * Nychunks
    Nchunksmax = neq
    Nixychunksmax = (nx+2)*(ny+2)
    dx = real(nx)/Nxchunks
    dy = real(ny)/Nychunks
    allocate(xlims(Nxchunks,2), ylims(Nychunks,2))!, Nivchunks(NchunksPandf1,neq), Niv(NchunksPandf1))
    call gchange("OMPPandf1", 0)


    xlims(1,1) = 0; xlims(1,2)=max(1,int(dx))
    do ix = 2, Nxchunks-1
        xlims(ix,1) = xlims(ix-1,2)+1
        xlims(ix,2) = int(dx*ix)
    end do
    if (Nxchunks .gt. 1) xlims(Nxchunks,1) = xlims(Nxchunks-1,2)+1
    xlims(Nxchunks,2) = nx+1
    ylims(1,1) = 0; ylims(1,2)=max(1,int(dy))
    do iy = 2, Nychunks-1
        ylims(iy,1) = ylims(iy-1,2)+1
        ylims(iy,2) = int(dy*iy)
    end do
    if (Nychunks .gt. 1) ylims(Nychunks,1) = ylims(Nychunks-1,2)+1
    ylims(Nychunks,2) = ny+1

    do ix = 1, Nxchunks
        do iy = 1, Nychunks
            rangechunk(Nxchunks*(iy-1) + ix, 1) = xlims(ix,1)
            rangechunk(Nxchunks*(iy-1) + ix, 2) = xlims(ix,2)
            rangechunk(Nxchunks*(iy-1) + ix, 3) = ylims(iy,1)
            rangechunk(Nxchunks*(iy-1) + ix, 4) = ylims(iy,2)
        end do
    end do
    Nivchunk = 0
    ivchunk = 0
    Nixychunk = 0
    ixychunk = 0
    do ix = 1, Nxchunks
        do iy = 1, Nychunks
            idxl = Nxchunks*(iy-1) + ix
            do ii = 1, neq
                idx = igyl(ii,:)
                if ((idx(1).ge.rangechunk(idxl,1)).and.(idx(1).le.rangechunk(idxl,2)) .and. &
                &   (idx(2).ge.rangechunk(idxl,3)).and.(idx(2).le.rangechunk(idxl,4)) &
                &   ) then
                    Nivchunk(idxl) = Nivchunk(idxl) + 1
                    ivchunk(idxl, Nivchunk(idxl)) = ii
                endif
            end do    
        end do
    enddo
    do ii = 1, NchunksPandf1
        do ix = rangechunk(ii,1), rangechunk(ii,2)
            do iy = rangechunk(ii,3), rangechunk(ii,4)
                Nixychunk(ii) = Nixychunk(ii) + 1
                ixychunk(ii, Nixychunk(ii),1) = ix
                ixychunk(ii, Nixychunk(ii),2) = iy
            end do
        end do
    end do


    Nchunksmax = MAXVAL(Nivchunk)
    Nixychunksmax = MAXVAL(Nixychunk)
    ! Trim overly large array
    call gchange("OMPPandf1",0)

  END SUBROUTINE Make2DChunks



