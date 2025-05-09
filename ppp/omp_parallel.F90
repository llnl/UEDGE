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
        call MakeChunksPandf1
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


  SUBROUTINE OMPconvsr_vo1(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Compla, ONLY: ne, phi, ti, tg, te, nit, nz2, ng, nm, lng, ni
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:ne_tmp, phi_tmp, ti_tmp, tg_tmp, te_tmp, nit_tmp, nz2_tmp, ng_tmp, nm_tmp, &
    !$OMP &         lng_tmp, ni_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call convsr_vo1(xc, yc, ylcopy)

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
    USE OmpCopybbb
    USE Compla, ONLY: up
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: up_tmp(0:nx+1,0:ny+1,1:nisp)

    ! Initialize arrays to zero
    up_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:up_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call convsr_vo2(xc, yc, ylcopy)

        ! Update locally calculated variables
        up_tmp(xc,yc,:)=up_tmp(xc,yc,:)+up(xc,yc,:)
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    up=up_tmp
    call OmpCopyPointerup

  END SUBROUTINE OMPconvsr_vo2

  SUBROUTINE OMPconvsr_aux1(neq, yl, yldot)
    USE Dim, ONLY: nx, ny, ngsp, nisp
    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk
    USE OmpCopybbb
    USE Compla, ONLY: pg, pr, pre, pri, tg
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: pg_tmp(0:nx+1,0:ny+1,1:ngsp), pr_tmp(0:nx+1,0:ny+1), &
    &      pre_tmp(0:nx+1,0:ny+1), pri_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      tg_tmp(0:nx+1,0:ny+1,1:ngsp)

    ! Initialize arrays to zero
    pg_tmp=0.; pr_tmp=0.; pre_tmp=0.; pri_tmp=0.; tg_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:pg_tmp, pr_tmp, pre_tmp, pri_tmp, tg_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call convsr_aux1(xc, yc)

        ! Update locally calculated variables
        pg_tmp(xc,yc,:)=pg_tmp(xc,yc,:)+pg(xc,yc,:)
        pr_tmp(xc,yc)=pr_tmp(xc,yc)+pr(xc,yc)
        pre_tmp(xc,yc)=pre_tmp(xc,yc)+pre(xc,yc)
        pri_tmp(xc,yc,:)=pri_tmp(xc,yc,:)+pri(xc,yc,:)
        tg_tmp(xc,yc,:)=tg_tmp(xc,yc,:)+tg(xc,yc,:)
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
    USE OmpCopybbb
    USE Compla, ONLY: pgy0, tgy0, niy1, nity1, tiy1, phiy0, tiy0s, ney0, zeff, tiy0, pgy1, phiy0s, phiy1, &
    &    ngy0, priy0, phiv, niy0s, tey1, tiy1s, priy1, tgy1, niy0, tiv, priv, ney1, tev, &
    &    phiy1s, prtv, prev, tey0, znot, niy1s, nity0, ngy1
    USE Gradients, ONLY: gpry, gpix, ex, ey, gtex, gtiy, gpiy, gtey, gtix, gpex, gpondpotx, gprx, gpey
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

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
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call convsr_aux2(xc, yc)

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
    USE OmpCopybbb
    USE Conduc, ONLY: dutm_use, difp_use, dif_use, kyi_use, trax_use, kxbohm, kxe_use, kxi_use, &
    &    kye_use, vy_use, kybohm, tray_use, dif2_use
    USE Compla, ONLY: betap
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:dutm_use_tmp, difp_use_tmp, dif_use_tmp, kyi_use_tmp, trax_use_tmp, &
    !$OMP &         kxbohm_tmp, kxe_use_tmp, kxi_use_tmp, kye_use_tmp, vy_use_tmp, betap_tmp, &
    !$OMP &         kybohm_tmp, tray_use_tmp, dif2_use_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_plasma_diffusivities

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
    USE OmpCopybbb
    USE Conduc, ONLY: eta1, dclass_e, rtaue, dclass_i
    USE UEpar, ONLY: ctaui, ctaue
    USE Compla, ONLY: loglambda
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:eta1_tmp, ctaui_tmp, dclass_e_tmp, loglambda_tmp, rtaue_tmp, ctaue_tmp, &
    !$OMP &         dclass_i_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call initialize_driftterms

        ! Update locally calculated variables
        eta1_tmp(xc,yc)=eta1_tmp(xc,yc)+eta1(xc,yc)
        ctaui_tmp(xc,yc,:)=ctaui_tmp(xc,yc,:)+ctaui(xc,yc,:)
        dclass_e_tmp(xc,yc)=dclass_e_tmp(xc,yc)+dclass_e(xc,yc)
        loglambda_tmp(xc,yc)=loglambda_tmp(xc,yc)+loglambda(xc,yc)
        rtaue_tmp(xc,yc)=rtaue_tmp(xc,yc)+rtaue(xc,yc)
        ctaue_tmp(xc,yc,:)=ctaue_tmp(xc,yc,:)+ctaue(xc,yc,:)
        dclass_i_tmp(xc,yc)=dclass_i_tmp(xc,yc)+dclass_i(xc,yc)
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
    USE OmpCopybbb
    USE Compla, ONLY: vydd, veycb, vyrd, vycf, vygp, vy, vycr, vyce, veycp, vycb, vycp
    USE Comtra, ONLY: coll_fe, diffusivwrk, coll_fi
    USE Conduc, ONLY: vyte_cft, vyti_cft, vy_cft
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:vydd_tmp, veycb_tmp, coll_fe_tmp, vyte_cft_tmp, vyrd_tmp, vycf_tmp, vygp_tmp, &
    !$OMP &         diffusivwrk_tmp, vy_tmp, vycr_tmp, vyce_tmp, coll_fi_tmp, veycp_tmp, &
    !$OMP &         vyti_cft_tmp, vycb_tmp, vycp_tmp, vy_cft_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_driftterms1

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
    USE OmpCopybbb
    USE Compla, ONLY: v2dd, ve2cb, vytan, v2rd, ve2cd, vy, v2xgp, v2cd, vyavis, v2ce, q2cd, v2cb, v2
    USE Comflo, ONLY: fdiaxlb, fdiaxrb
    USE Xpoint_indices, ONLY: ixlb, ixrb
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, jx
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:v2dd_tmp, ve2cb_tmp, vytan_tmp, v2rd_tmp, fdiaxlb_tmp, ve2cd_tmp, vy_tmp, &
    !$OMP &         v2xgp_tmp, v2cd_tmp, vyavis_tmp, v2ce_tmp, q2cd_tmp, v2cb_tmp, v2_tmp, &
    !$OMP &         fdiaxrb_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_driftterms2


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
    USE OmpCopybbb
    USE Comflo, ONLY: fq2d, fmity, fqym, fqyai, fqyb, fqydti, fqymi, fqydt, fqyao, fqya, fqyae, fqygp, &
    &    fq2, fqyd, fqy
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fq2d_tmp, fmity_tmp, fqym_tmp, fqyai_tmp, fqyb_tmp, fqydti_tmp, fqymi_tmp, &
    !$OMP &         fqydt_tmp, fqyao_tmp, fqya_tmp, fqyae_tmp, fqygp_tmp, fq2_tmp, fqyd_tmp, &
    !$OMP &         fqy_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_currents

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
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc, jx
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

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:fqp_tmp,netap_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_fqp1

        ! Update locally calculated variables
        fqp_tmp(xc,yc)=fqp_tmp(xc,yc)+fqp(xc,yc)
        netap_tmp(xc,yc)=netap_tmp(xc,yc)+netap(xc,yc)
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
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_fqp2
    
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
    USE OmpCopybbb
    USE UEpar, ONLY: cs
    USE Compla, ONLY: upi, uu, uz, uup
    USE Cfric, ONLY: frici, frice
    USE Gradients, ONLY: ex
    IMPLICIT NONE
    INTEGER, INTENT(IN):: neq
    REAL, INTENT(IN):: yl(*)
    REAL, INTENT(OUT):: yldot(*)
    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc
    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)
! Define local variables
    real:: upi_tmp(0:nx+1,0:ny+1,1:nisp), &
    &      uu_tmp(0:nx+1,0:ny+1,1:nisp), frici_tmp(0:nx+1,0:ny+1,nisp), &
    &      uz_tmp(0:nx+1,0:ny+1,1:nisp), frice_tmp(0:nx+1,0:ny+1), &
    &      uup_tmp(0:nx+1,0:ny+1,1:nisp), ex_tmp(0:nx+1,0:ny+1)

    ! Initialize arrays to zero
    upi_tmp=0.; uu_tmp=0.; frici_tmp=0.; uz_tmp=0.; frice_tmp=0.
    uup_tmp=0.; ex_tmp=0.

    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0

    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)

    !$OMP    PARALLEL DO &
    !$OMP &      default(shared) &
    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &
    !$OMP &      private(ichunk,xc,yc) &
    !$OMP &      firstprivate(ylcopy, yldotcopy) &
    !$OMP &      REDUCTION(+:upi_tmp, uu_tmp, frici_tmp, uz_tmp, frice_tmp, uup_tmp, ex_tmp)
    DO ichunk = 1, Nchunks
        xc = chunks(ichunk,1)
        yc = chunks(ichunk,2)
        call initialize_ranges(xc, yc, 0, 0, 0)
        call calc_friction(xc)

        ! Update locally calculated variables
        upi_tmp(xc,yc,:)=upi_tmp(xc,yc,:)+upi(xc,yc,:)
        uu_tmp(xc,yc,:)=uu_tmp(xc,yc,:)+uu(xc,yc,:)
        frici_tmp(xc,yc,:)=frici_tmp(xc,yc,:)+frici(xc,yc,:)
        uz_tmp(xc,yc,:)=uz_tmp(xc,yc,:)+uz(xc,yc,:)
        frice_tmp(xc,yc)=frice_tmp(xc,yc)+frice(xc,yc)
        uup_tmp(xc,yc,:)=uup_tmp(xc,yc,:)+uup(xc,yc,:)
        ex_tmp(xc,yc)=ex_tmp(xc,yc)+ex(xc,yc)
    END DO
    !$OMP END PARALLEL DO

    ! Update global variables
    upi=upi_tmp; uu=uu_tmp; frici=frici_tmp; uz=uz_tmp
    frice=frice_tmp; uup=uup_tmp; ex=ex_tmp
    call OmpCopyPointerupi; call OmpCopyPointeruu
    call OmpCopyPointerfrici; call OmpCopyPointeruz
    call OmpCopyPointerfrice; call OmpCopyPointeruup; call OmpCopyPointerex

  END SUBROUTINE OMPcalc_friction

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
    USE Imprad, ONLY: prad
    USE Selec, ONLY:yinc,xrinc,xlinc
    USE Grid, ONLY:ijactot
    USE Cdv, ONLY: comnfe
    USE Rhsides, ONLY: psorcxg, psorrg, psordis
    USE Time_dep_nwt, ONLY: dtreal, nufak
    USE Ynorm, ONLY: isflxvar, isrscalf
    USE MCN_sources, ONLY: ismcnon
    USE UEpar, ONLY: isphion, svrpkg, isphiofft
    USE PandfTiming, ONLY: TimePandf, TotTimePandf, TimingPandfOn
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
    integer::ichunk, xc, yc
    real tmp_prad(0:nx+1, 0:ny+1)
    real tick,tock, tsfe, tsjf, ttotfe, ttotjf
    external tick, tock
    ylcopy(1:neq+1)=yl(1:neq+1)
    yldotcopy = 0
    yldottot = 0
    tmp_prad = 0
    xc=-1; yc=-1

    if (ijactot.gt.0) then
        Time1=omp_get_wtime()
        call MakeChunksPandf1

                ! ... Get initial value of system cpu timer.
                if(xc .lt. 0) then
                    tsfe = tick()
                else
                     tsjf = tick()
                endif
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

                call initialize_ranges(xc, yc, xlinc, xrinc, yinc)
!                if(isphion+isphiofft .eq. 1) call calc_fqp

                ! TODO: gather variables calculated in calc driftterms
                !       v2 needed by calc_friction
                ! TODO: Break out conditionals, move to top
                call calc_friction(xc)
                !******************************************************
                !*     Calculate the currents fqx, fqy, fq2 and fqp, if
                !*     isphion = 1 or if isphiofft = 1.
                !******************************************************

                call calc_elec_velocities
                ! Add checks on ishosor and ispsorave: parallel only works for == 0
                call calc_volumetric_sources(xc, yc)




                call calc_plasma_viscosities

                call calc_plasma_heatconductivities

                call calc_plasma_equipartition

                call calc_gas_heatconductivities
                ! ... Call routine to evaluate gas energy fluxes
                !******************************************************
                call engbalg

                call calc_plasma_transport

                !------------------------------------------------------
                !   SCALE SOURCE TERMS FROM MONTE-CARLO-NEUTRALS MODEL
                !
                ! These sources are used in the residuals (resco,resmo,
                ! resee,resei) so the call to scale_mcn must occur 
                ! BEFORE these residuals are evaluated.  Since they 
                ! scale with fnix at the divertor plates, the call to 
                ! scale_mcn must occur AFTER fnix has been calculated.

                if (ismcnon .ne. 0) call scale_mcnsor

                !******************************************************
                !  Here we do the neutral gas diffusion model
                !  The diffusion is flux limited using the thermal flux
                !******************************************************

                call calc_plasma_momentum(xc, yc)


                call calc_plasma_energy

                call calc_plasma_particle_residuals
                call calc_gas_continuity_residuals
                call calc_plasma_momentum_residuals()
                call calc_gas_energy_residuals
                !  Requires gas energy residuals
                call calc_plasma_energy_residuals(xc, yc)
                call calc_rhs(yldot)

                !  POTEN calculates the electrostatic potential, and 
                !  BOUNCON calculates the equations for the boundaries.
                !  For the vodpk solver, the B.C. are ODEs in time 
                !  (rate equations).  Both bouncon and poten must be 
                !  called before the perturbed variables are reset 
                !  below to get Jacobian correct

                if (isphion.eq.1) call calc_potential_residuals (neq, yl, yldot)

                call bouncon (neq, yldot)

                ! Accumulate cpu time spent here.
                if(xc .lt. 0) then
                    ttotfe = ttotfe + tock(tsfe)
                else
                    ttotjf = ttotjf + tock(tsjf)
                endif
                if (TimingPandfOn.gt.0) & 
                &      TotTimePandf=TotTimePandf+tock(TimePandf)


                ! ================ BEGIN OLD PANDF1 ===================

                ! If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and
                ! the ODEs need to be modified as original equations 
                ! are for d(nv)/dt, etc If isflxvar=2, variables are 
                ! ni,v,nTe,nTi,ng. Boundary equations and potential 
                ! equations are not reordered.

                if(isflxvar.ne.1 .and. isrscalf.eq.1) call rscalf(yl,yldot)

                ! Now add psuedo or real timestep for nksol method, but not both
                if (nufak.gt.1.e5 .and. dtreal.lt.1.e-5) then
                    call xerrab('***Both 1/nufak and dtreal < 1.e5 - illegal***')
                endif


                ! Add a real timestep, dtreal, to the nksol equations 
                ! NOTE!! condition yl(neq+1).lt.0 means a call from nksol, not jac_calc
                if(dtreal < 1.e15) then
                    if ( &
                    &   (svrpkg=='nksol' .and. yl(neq+1)<0) &
                    &   .or. svrpkg == 'petsc' &
                    & ) then
                        call add_timestep(neq, yl, yldot)
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
  


