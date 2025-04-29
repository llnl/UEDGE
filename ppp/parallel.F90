#include "../bbb/bbb.h"
#include "../com/com.h"
!-------------------------------------------------------------------------------------------------
SUBROUTINE InitParallel

    USE ParallelSettings,ONLY: OMPParallelJac,OMPParallelPandf1,ParallelWARNING
    USE ParallelEval, ONLY: ParallelJac,ParallelPandf1
    IMPLICIT NONE

    if (OMPParallelJac==1) then
        if (ParallelWARNING.gt.0) then
            write(*,*) "<<< WARNING >>> : You are using openmp routines to evaluate the Jacobian or/and pandf1."
            write(*,*) "<<< WARNING >>> : Read first the documentation about how to use these options."
            write(*,*) "<<< WARNING >>> : If this is your first run with this configuration, you should verify whether&
            the parallel evaluation of pandf1 and jacobian are correct by setting ppp.CheckJac=1 and ppp.OMPCheckPandf1=1.&
             If not error is thrown, turn these settings off (=0)."
            write(*,*) "<<< WARNING >>> : You can turn off this warning by setting ppp.ParallelWarning=0."
        endif

        if (OMPParallelJac.gt.0) then
            call InitOMPJac
            ParallelJac=1
        else
            ParallelJac=0
        endif
    else
        ParallelJac=0
    endif

    if (OMPParallelPandf1.gt.0) then
        if (OMPParallelJac==0) then
            write(*,*) "Parallel Pandf1  requires Parallel Jac: activating..."
            OMPParallelJac=1
            call InitOMPJac
            ParallelJac=1
        endif
        call InitOMPPandf1()
        ParallelPandf1=1
    else
        ParallelPandf1=0
    endif
    RETURN
END SUBROUTINE InitParallel
!-------------------------------------------------------------------------------------------------
SUBROUTINE jac_calc_parallel(neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    USE ParallelSettings, ONLY: OMPParallelJac
    USE ParallelSettings, ONLY: CheckJac
    USE ParallelDebug, ONLY: iidebugprint,ivdebugprint,DebugJac,ForceSerialCheck,DumpFullJac,DumpJac
    USE Cdv, ONLY:exmain_aborted
    IMPLICIT NONE
    ! ... Input arguments:
    integer,intent(in):: neq      !      total number of equations (all grid points)
    real,intent(in)   :: t              ! physical time
    real,intent(inout)   :: yl(*)          ! dependent variables
    real,intent(in)   :: yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu         ! lower and upper bandwidths
    integer,intent(in):: nnzmx          ! maximum number of nonzeros in Jacobian

    ! ... Output arguments:
    real,intent(out)   :: jac(nnzmx)     ! nonzero Jacobian elements
    integer,intent(out):: ja(nnzmx)   ! col indices of nonzero Jacobian elements
    integer,intent(out):: ia(neq+1)   ! pointers to beginning of each row in jac,ja
    real   :: jaccopy(nnzmx)     ! nonzero Jacobian elements
    integer:: jacopy(nnzmx)   ! col indices of nonzero Jacobian elements
    integer:: iacopy(neq+1)   ! pointers to beginning of each row in jac,ja
    integer i,iv
    real t_start
    real,external::tick,tock
    ! ... Work-array argument:
    real wk(neq)     ! work space available to this subroutine
    if (exmain_aborted) call xerrab('exmain aborted...')
    if (OMPParallelJac==1) then
        call jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
    else
        call xerrab('Cannot call calc_jac_parallel OMP jac calc is not enabled')
    endif

    if (DumpFullJac.gt.0) then
        write(*,*) 'Dumping full serial jacobian for analysis of bandwidth'
        call jac_calc (neq, t, yl, yldot00, neq, neq, wk,neq*neq, jaccopy, jacopy, iacopy)
        call jac_write('serialfulljac.dat',neq, jaccopy, jacopy, iacopy)
    endif

    if (ForceSerialCheck.gt.0) then
        write(*,*) 'Force check of serial evaluation of jacobian'
        write(*,*) '----- Performing first serial Evaluation of jacobian...'
        call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
        write(*,*) '----- Performing second serial Evaluation of jacobian...'
        call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jaccopy, jacopy, iacopy)
        write(*,*) '----- Comparing jacobians...: nzmax=',ia(neq)-1
        do i=1,ia(neq)-1
            if (abs(jaccopy(i)-jac(i)).gt.1e-14) then
                write(*,*) ' ****** diff between jacs :',i,jac(i),jaccopy(i),'|',ja(i),jacopy(i)
                do iv=1,neq
                    if (i<iacopy(iv)) then
                        write(*,*) 'idx row serial:',iv-1
                        exit
                    endif
                enddo
                do iv=1,neq
                    if (i<ia(iv)) then
                        write(*,*) 'idx row parallel:',iv-1
                        exit
                    endif
                enddo
                call xerrab('Parallel evaluation of Jacobian differs from serial evaluation... ')
                exit
            endif
        enddo
    endif

    if (CheckJac.gt.0) then
        write(*,*) '--- Checking whether parallel and serial evaluations of jacobian are the same...'
        write(*,*) '----- Performing serial Evaluation of jacobian...'
        t_start=tick()
        call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jaccopy, jacopy, iacopy)
        write(*,*) '----- Serial Evaluation of jacobian done in: ',tock(t_start) ,' s'
        write(*,*) '----- Comparing jacobians (nzmax=',ia(neq)-1,')'
        do i=1,ia(neq)-1
            if (abs(jaccopy(i)-jac(i)).gt.1e-14) then
                write(*,*) ' ****** diff between jacs :',i,jac(i),jaccopy(i),'|',ja(i),jacopy(i)
                do iv=1,neq
                    if (i<iacopy(iv)) then
                        write(*,*) 'idx row serial:',iv-1
                        exit
                    endif
                enddo
                do iv=1,neq
                    if (i<ia(iv)) then
                        write(*,*) 'idx row parallel:',iv-1
                        exit
                    endif
                enddo
                write(*,*) '----- Writing jacobian into serialjac.dat and paralleljac.dat'
                call jac_write('paralleljac.dat',neq, jac, ja, ia)
                call jac_write('serialjac.dat',neq, jaccopy, jacopy, iacopy)

                if (DumpJac.gt.0) then
                    iidebugprint=ja(i)
                    ivdebugprint=iv-1
                    write(*,*) 'recalculating jacobian to dump debug data for iv=',ivdebugprint,' ii=',iidebugprint
                    if (OMPParallelJac==1) then
                        call jac_calc_omp (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jac, ja, ia)
                    else
                        call xerrab('Cannot call calc_jac_parallel OMP jac calc is not enabled')
                    endif
                    call jac_calc (neq, t, yl, yldot00, ml, mu, wk,nnzmx, jaccopy, jacopy, iacopy)
                endif
                call xerrab('Parallel evaluation of Jacobian differs from serial evaluation... ')
                exit
            endif
        enddo
        write(*,*) '--- Parallel and serial evaluations of jacobian are identical...'
    endif

END SUBROUTINE jac_calc_parallel

!-------------------------------------------------------------------------------------------------
SUBROUTINE LocalJacBuilder(ivmin,ivmax,neq, t, yl, yldot00,ml, mu,iJacCol,rJacElem,iJacRow,ichunk,nnz,&
                            nnzmxperchunk,TimeJacRow)
    ! ... Calculate Jacobian matrix (derivatives with respect to each
    ! ... Calculate Jacobian matrix (derivatives with respect to each
    !     dependent variable of the right-hand side of each rate equation).
    !     Lower and upper bandwidths are used to select for computation
    !     only those Jacobian elements that may be nonzero.
    !     Estimates of Jacobian elements are computed by finite differences.
    !     The Jacobian is stored in compressed sparse row format.
    USE Dim, ONLY: nx
    USE Timing, ONLY: istimingon,ttjstor,ttotjf,ttimpjf
    USE Math_problem_size, ONLY: neqmx,numvar
    USE Indexes, ONLY: igyl,iseqalg,idxphi
    USE Variable_perturbation, ONLY: delperturb,dylconst,isjacreset
    USE Jacobian_clipping, ONLY: jaccliplim,istopjac,irstop,icstop
    USE Ynorm, ONLY: suscal,sfscal
    USE UEpar, ONLY: isphion,isnewpot,svrpkg,isbcwdt
    USE Model_choice, ONLY: iondenseqn
    USE Bcond, ONLY: isextrnpf,isextrtpf,isextrngc,isextrnw,isextrtw
    USE Time_dep_nwt, ONLY: nufak,dtreal,ylodt,dtuse,dtphi
    USE OMPJacSettings,ONLY: OMPTimingJacRow
    USE ParallelDebug, ONLY: iidebugprint,ivdebugprint
    USE Jacaux, ONLY: ExtendedJacPhi
    USE omp_lib
    IMPLICIT NONE

    integer,intent(in):: ichunk,nnzmxperchunk,ivmin,ivmax,neq
    integer,intent(inout)::nnz
    real,intent(in):: t           ! physical time
    real,intent(inout) ::yl(*)       ! dependent variables
    integer,intent(inout)::iJacCol(nnzmxperchunk),iJacRow(neq)
    real,intent(inout):: rJacElem(nnzmxperchunk)
    real,intent(in) ::yldot00(neq+2) ! right-hand sides evaluated at yl
    integer,intent(in):: ml, mu   ! lower and upper bandwidths
    real,intent(out)::TimeJacRow(neq)
    ! ... Work-array argument:
    real :: wk(neq)     ! work space available to this subroutine
    real :: yl_check(neq)     ! work space available to this subroutine
    ! ... Functions:
    logical ::tstguardc

    ! ... Local variables:
    real ::yold, dyl, jacelem,Time
    real(kind=4) sec4, tsjstor, tsimpjf, dtimpjf
    integer:: ii, iv, jv,ii1, ii2, xc, yc, ix, iy,tid


    ! ... Set beginning and ending indices of right-hand sides that might be
    !     perturbed.
    external:: xerrab, kaboom

    nnz=1
    loopiv: do iv=ivmin,ivmax

        if (OMPTimingJacRow.gt.0) Time=omp_get_wtime()
        ii1 = max(iv-mu, 1)
        ii2 = min(iv+ml, neq)
        ! ... Reset range if this is a potential perturbation with isnewpot=1
        if (ExtendedJacPhi.eq.1) then
            if  (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
                ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
                ii2 = min(iv+4*numvar*nx, neq)    ! 3*nx may be excessive
            endif
        else if (ExtendedJacPhi.eq.2) then
            if(isphion*isnewpot.eq.1) then
                ii1 = max(iv-4*numvar*nx, 1)      ! 3*nx may be excessive
                ii2 = min(iv+4*numvar*nx, neq)   ! 3*nx may be excessive
            endif
        endif

        ! ... Reset range if extrapolation boundary conditions are used
        if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
            ii1 = max(iv-2*numvar*(nx+3), 1)      ! guess to include extrap. bc
            ii2 = min(iv+2*numvar*(nx+3), neq)    ! guess to include extrap. bc
        endif

        ! ... Initialize all of those right-hand sides to their unperturbed
        !     values.
        do ii = ii1, ii2   ! a) below wk is reset, but only over limited range
            wk(ii) = yldot00(ii)
        enddo

        ! ... Set spatial-location indices for this dependent variable.
        xc = igyl(iv,1)
        yc = igyl(iv,2)

        ! ... Save value of dependent variable, then perturb it.
        !     The perturbation to the variable is proportional to parameter
        !     del and to a measure of the size of the variable.  That measure
        !     increases with the absolute value of the variable if it exceeds
        !     the typical size given by dylconst/suscal but can never be less
        !     than that typical size.

        yold = yl(iv)
        dyl = delperturb * (abs(yold) + dylconst / suscal(iv))
        yl(iv) = yold + dyl

        !Calculate right-hand sides near location of perturbed variable.
        call pandf (xc, yc, neq, t, yl, wk)

        !Calculate possibly nonzero Jacobian elements for this variable,
        !and store nonzero elements in compressed sparse column format.
        iJacRow(iv) = nnz      ! sets index for first Jac. elem. of var. iv
                                       ! the index is set to start at 0

        loopii: do ii = ii1, ii2
            jacelem = (wk(ii) - yldot00(ii)) / dyl

            !Add diagonal 1/dt for nksol
            ifdiagonal: if (iv.eq.ii) then
                if (iseqalg(iv)*(1-isbcwdt).eq.0) then
                    jacelem = jacelem - 1/dtuse(iv)
                endif
                ix = igyl(iv,1)
                iy = igyl(iv,2)
                if (idxphi(ix,iy)==iv .and. dtphi<1e10) then ! selects phi eqn
                    jacelem = jacelem - 1/dtphi
                endif
            endif ifdiagonal

            ! ...  Add a pseudo timestep to the diagonal #! if eqn is not algebraic
            if (nufak .gt. 0) then
                if (iv.eq.ii .and. yl(neq+1).eq.1) jacelem = jacelem - nufak 
            endif

            ! Debug
            debug: if (ii==iidebugprint.and.iv==ivdebugprint) then
                tid=omp_get_thread_num()
                write(*,'(a,i3,a,i7,i7,i7, E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)') '#', tid,&
                    ' : ',iv,nnz,ii,jacelem,wk(ii),yldot00(ii),sfscal(iv),jaccliplim,nufak
                call DebugHelper('dumpomp.txt') ! See at the end of this file
            endif debug

            storage: if (abs(jacelem*sfscal(iv)) .gt. jaccliplim) then
                if (nnz .gt. nnzmxperchunk) then
                    write(*,*) 'nnz=',nnz,'nnzmxperchunk=',nnzmxperchunk
                    write(*,*) 'chunk',ichunk,'*** jac_calc -- More storage needed for Jacobian. Increase omplenpfac.'
                    call xerrab("")
                endif
                iJacCol(nnz)=ii
                rJacElem(nnz)=jacelem
                nnz = nnz + 1
            endif storage

            check: if (istopjac.gt.0 .and. ii.eq.irstop .and. iv.eq.icstop) then
                if (istopjac == 2) then
                    yl(iv) = yold
                    call pandf (xc, yc, neq, t, yl, wk)
                endif
                call remark("***** non-zero jac_elem at irstop,icstop")
                write(*,*) 'irstop = ', irstop, ', icstop = ', icstop
                call xerrab("")
            endif check

        enddo loopii
        ! ... Restore dependent variable and plasma variables near its location.
        yl(iv) = yold

        call pandf (xc, yc, neq, t, yl, wk)

        !...  If this is the last variable before jumping to new cell, reset pandf
        reset: if( isjacreset.ge.1) then
            yl_check(1:neq)=wk(1:neq)
            !JG this call to pandf1 can be safely ignored with ijacreset=0 (and save some time...)
            if (mod(iv,numvar).eq.0) then
                call pandf (xc, yc, neq, t, yl, wk)
            endif

            do ii=ii1,ii2
                if (isnan(yl_check(ii))) then
                    write(*,*) 'NaN at ii=',ii,'equation #', mod(ii,numvar)+1
                    call xerrab('*** Nan in wk array in jac_calc ***')
                endif
                if (yl_check(ii).ne.wk(ii)) then
                    write(*,'(a,i5,e20.12,e20.12)') ' *** wk modified on second call to pandf1 at ii=', ii,yl_check(ii),wk(ii)
                    call xerrab('*** Stop ***')
                endif

            enddo
        endif reset
        if (OMPTimingJacRow.gt.0) then
            Time=omp_get_wtime()-Time
            TimeJacRow(iv)=Time
        endif
    enddo loopiv
! ... End loop over dependent variables and finish Jacobian storage.
END SUBROUTINE LocalJacBuilder
!-------------------------------------------------------------------------------------------------

SUBROUTINE Compare(yldot,yldotsave,neq)
    integer:: iv,neq
    integer::istop
    real yldot(neq+2),yldotsave(neq+2)
    istop=0
    do iv=1,neq
        if (abs(yldotsave(iv)-yldot(iv)).gt.1e-14) then
            if (max(abs(yldot(iv)),abs(yldotsave(iv)))>0) then
                if (abs(yldotsave(iv)-yldot(iv))/max(abs(yldot(iv)),abs(yldotsave(iv)))>1e-14) then
                    write(*,*) '>>>>',iv,yldotsave(iv),yldot(iv),abs(yldotsave(iv)-yldot(iv))/max(abs(yldot(iv)),abs(yldotsave(iv)))
                    istop=1
                endif
            endif
        endif
    enddo
    if (istop.gt.0) call xerrab('diff in rhsnk')

END SUBROUTINE Compare



