c
c $Id: nksol.m,v 7.4 2021/03/23 18:59:33 meyer8 Exp $
c

c!include "../mppl.h"
c!include "../sptodp.h"

      subroutine nksol (n, u, savf, f, jac, su, sf, ftol, stptol,
     *                  rwork, lrw, iwork, liw, iopt, iterm, pset,
     *                  psol, mf, mdif, ipflag, icflag, icnstr, rlx,
     *                  epscon1, epscon2, icntnu, adjf1)
c-----------------------------------------------------------------------
c this is the august 7, 1995 version of
c nksol.. a nonlinear krylov solver for nonlinear systems of
c         equations of the form f(u) = 0.
c-- this version special for use with uedge.
c
c this version is in single precision.
c
c authors..
c
c         peter n. brown
c         computing and mathematics research division
c         lawrence livermore national laboratory
c         livermore, ca  94550
c and
c         youcef saad
c         center for supercomputer research and development
c         university of illinois
c         urbana, il  61801
c
c references..
c 1. peter n. brown and youcef saad
c    hybrid krylov methods for nonlinear systems of equations
c    llnl report ucrl-97645, november 1987.
c
c 2. john e. dennis and robert b. schnabel
c    numerical methods for unconstrained optimization and nonlinear
c    equations, prentice-hall, englewood cliffs, nj, 1983.
c    (this reference contains the basic algorithms for the
c    stopping tests, dogleg and linesearch strategies used
c    in nksol.)
c
c revision history...
c  2-12-93  constraint checking for the approximate solution u was added.
c           this is only active when mf=1.  see descriptions for icflag
c           and array icnstr.
c  3-01-93  added wmp and iwmp to calling sequence for subroutine jac.
c           changed (1) to (*) in dimension statements to allow for
c           array out of bounds checking.
c  3-03-93  fixed possible infinite loop when using constraints.
c           if a constraint violation occurred, then trgupd strategy
c           would still try to increase the trust region size for the
c           current step.
c  3-05-93  generated an informational output routine.  informational
c           messages are printed when iermsg=0.
c  3-10-93  fixed a bug in the lnsrch module that left unew possibly
c           not satisfying the alpha condition when the beta condition
c           could not be satisfied.
c  3-15-93  fixed a bug that caused double scaling of j*v by sf in the
c           atv routine when using the user-supplied jac routine.
c           changed comments below to note that the preconditioner p
c           must well approximate the scaled jacobian matrix sf*j.
c  4-23-93  added logic to update the preconditioner less frequently.
c           new optional input, incpset, can be set to control
c           the frequency of preconditioner evaluations.
c  5-18-93  added rlx to restrict update to del u / u < rlx in cnstrt
c 10-12-93  corrected bug fix of 3-10-93 with proper fix.  the savf
c           and f1nrmp values were not always correct when exiting
c           the lnsrch routine.
c 10-18-93  changed intrinsic names amax1 and amin1 to generic
c           intrinsic names max and min (Gary R. Smith).
c           eliminated implicit variable typing in all routines by
c           adding appropriate integer and/or real statements as well
c           as implicit none (Gary R. Smith).
c 10-19-93  changed function name vnorm to vnormnk to avoid conflict
c           with vnorm included in some versions of lsode source
c           (Gary R. Smith).
c 11-16-93  Added epscon1 and epscon2 as input variables to define the
c           tolerance level for the linear iteration: 
c           epsfac = epscon1 * min(epscon2, frnm) (Tom Rognlien)
c  8-07-95  added arrays su and sf to the calling sequence for psol.
c           the calling sequence for dogstp needed to add sf.
c  2-01-95  added flag icntnu to determine if this is a continuation
c           call to nksol that makes use of old values, e.g., the
c           routine pset is not called initially to reform the
c           preconditioner, and various counters in common blocks
c           are saved.
c  6-15-98  added variable adjf1 which is passed to subroutine lnsrch
c           for adjusting the acceptance test for the mfnksol=3 option
c           which uses global strategy; test is fnrm_new/adjf1 >= frnm_old
c-----------------------------------------------------------------------
c introduction.
c
c nksol solves nonlinear systems f(u)=0 rewritten in the form
c
c              sf*f((su-inverse)*ubar) = 0, ubar = su*u,
c
c where f and u are n-vectors, and sf and su are diagonal scaling
c matrices with positive diagonal entries.  nksol uses an inexact
c newton method as the basic nonlinear iteration, where the newton
c equations are solved only approximately by a linear krylov iteration,
c coupled with either a linesearch or dogleg global strategy.  the
c user may optionally choose either arnoldi-s method (with linesearch
c backtracking) or the generalized minimum residual method (gmres) 
c (with either the linesearch or dogleg strategy) as the krylov
c iteration technique, with or without preconditioning.
c
c the scaling matrices su and sf should be chosen so that the vector
c su*u has all its components with roughly the same magnitude when u
c is close to a root of f, and similarly for sf*f(u) when u is not too
c near a root of f.
c
c
c nksol generates a sequence of approximations u(k) to a root of f.
c the stopping criteria for the nonlinear iteration are the
c following..
c
c 1.  maxnorm( sf*f(u(k)) ) .le. ftol, 
c
c where maxnorm() is the maximum norm function and ftol is a
c user-supplied stopping tolerance.
c
c 2.  max(abs(u(k,i)-u(k-1,i))/max(abs(u(k,i)),1/su(i))) .le. stptol,
c    i=1,...,n
c
c where u(j,i) is the i-th component of u(j) (j=k,k-1) , su(i) is
c the i-th diagonal element of the matrix su, and stptol is a user-
c supplied stopping tolerance.  the above test measures the length
c of the scaled distance between the last two iterates.
c
c
c the scaled/preconditioned newton equations are of the form
c
c     (sf*j*(p-inverse)*(su-inverse))*(su*p*x) = sf*b,
c
c where j = j(u) = df/du is the system jacobian matrix evaluated at the
c current iterate u, b = -f(u), and p is a preconditioner matrix
c (applied on the right).  the solution x is the approximate newton
c step.  note that when using a nontrivial scaling array sf, the
c preconditioner p must well approximate the scaled jacobian matrix
c sf*j.  The idea is that sf*j*(p-inverse) should be close to the
c identity matrix.
c
c the linear krylov iteration generates a sequence x(m) of
c approximate solutions to the newton equations within an m-dimensional
c krylov subspace (m = 1,2,...).  the stopping test for the linear
c krylov iteration is
c
c     norm(sf*(b - j*x(m))) .le. etak*norm(sf*b),
c
c where norm() is the euclidean norm, b = -f(u(k)), j = j(u(k)),
c and etak = 0.5**k.  the counter ncfl counts the number of times
c this condition is not met within mmax linear iterations.  (see the
c optional inputs and outputs sections below for the definitions of
c mmax and ncfl.)
c
c-----------------------------------------------------------------------
c full description of user interface to nksol.
c
c the user interface to nksol consists of the following parts.
c
c i.   the call sequence to subroutine nksol, which is a driver
c      routine for the solver.  this includes descriptions of both
c      the call sequence arguments and of user-supplied routines.
c      following these descriptions is a description of
c      optional inputs available through the call sequence, and then
c      a description of optional outputs (in the work arrays).
c
c ii.  description of two routines in the nksol package which
c      the user may replace with his own versions, if desired.
c      these relate to the above stopping criteria.
c
c-----------------------------------------------------------------------
c part i.  call sequence.
c
c the call sequence parameters used for input only are
c  n, f, jac, su, sf, ftol, stptol, lrw, liw, iopt, pset, psol, mf,
c  mdif, ipflag,
c those used for both input and output are
c  u,
c and those used only for output are
c  savf, iterm.
c the work arrays rwork and iwork are also used for conditional and
c optional inputs and optional outputs.  (the term output here refers
c to the return from subroutine nksol to the user-s calling program.)
c
c the legality of input parameters will be thoroughly checked on each
c call to nksol.
c
c the descriptions of the call arguments are as follows.
c
c n      = integer scalar containing the size of the nonlinear system,
c          and the length of arrays u, savf, su and sf.
c
c u      = real array of length n containing the initial guess to a
c          root of f(u) = 0 on input, and the approximate solution
c          on return.
c
c          this array is passed as the u argument in all calls to
c          f, jac, pset and psol.  hence its length may exceed n, and
c          locations u(n+1),... may be used to store other real data
c          and pass it to f (or jac), pset, and/or psol.
c
c savf   = real array of length n containing f(u) on return from
c          nksol, and undefined on input.
c
c f      = the name of the user-supplied subroutine for defining the
c          nonlinear system f(u) = 0.  f is a vector-valued function
c          of the vector u.  subroutine f is to compute the function
c          f.  it is to have the form 
c           subroutine f(n, u, savf)
c           dimension u(*), savf(n)
c          where n and u are input and the array savf = f(u) is output.
c          subroutine f should not alter u.  f must be declared
c          external in the calling program.
c
c          subroutine f may access user-defined quantities in
c          u(n+1),... if u has length exceeding n in the calling
c          program.  see description of u above.
c
c jac    = the name of the optional user-supplied subroutine for
c          calculating j(u)*v, where j(u) is the jacobian matrix of
c          f evaluated at u.  subroutine jac returns j(u)*v for a
c          given vector v.  it is to have the form
c           subroutine jac (n, u, savf, v, z, wmp, iwmp)
c           dimension u(*), savf(n), v(n), z(n), wmp(*), iwmp(*)
c          where n, u, savf and v are input, and z = j(u)*v as output.
c          the array u contains the last newton iterate and
c          savf = f(u).  subroutine jac should not alter u or savf.
c          jac must be declared external in the calling program, and
c          is only called when mdif=1.  when mdif=0, a dummy routine
c          may be supplied for jac.  any matrix data used by this
c          routine must be stored in the real and integer work
c          arrays wmp and iwmp.  these arrays are passed unaltered
c          to subroutine jac.
c
c          note: the lengths of arrays wmp and iwmp must be supplied
c                as optional inputs in the iwork array.  see the
c                optional inputs section below.
c
c          subroutine jac may also access user-defined quantities in
c          u(n+1),... if u has length exceeding n in the calling
c          program.  see description of u above.
c
c su     = real array of length n containing scale factors for the
c          solution vector u.
c
c sf     = real array of length n containing scale factors for the
c          function value f(u).
c
c          arrays su and sf can be the same array in the user-s
c          calling program.  if the user does not want to use scaling,
c          then su(i) and sf(i) should be loaded with 1 for all i.
c
c ftol   = real scalar containing the stopping tolerance on
c          maxnorm(sf*f(u)).  if ftol=zero on input, then
c          a default value of epsmch**(1/3) will be used, where
c          epsmch is the machine epsilon (or unit roundoff).
c
c stptol = real scalar containing the stopping tolerance on the
c          minimum scaled step u(k) - u(k-1).  if stptol=zero on
c          input, then a default value of epsmch**(2/3) will be
c          used, where epsmch is the machine epsilon (or unit
c          roundoff).
c
c mf     = integer scalar containing the nonlinear iteration
c          method flag.
c          mf = 1 means the dogleg global strategy with gmres is used.
c          mf = 2 means the linesearch backtracking with descent
c                 direction chosen by arnoldi is used.
c          mf = 3 means the linesearch backtracking with descent
c                 direction chosen by gmres is used.
c          mf negative implies no global strategy is used:
c          mf =-1 means inexact newton with gmres is used.
c          mf =-2 means inexact newton with arnoldi is used.
c          mf =-3 means inexact newton with gmres is used.
c
c rwork  = a real working array (single precision).
c          the length of rwork must be at least
c                4 + 4*n + lenk + lenwmp,
c          where
c           lenk = 4 + (4+mmax)*n + mmax**2,               if |mf|=2,
c           lenk = 6 + 4*n + (n+1)*mmax + 2*mmax*(mmax+1), if |mf|=1,3.
c          lenwmp is the length of the user-defined workspace for the
c          preconditioner routines pset and psol.  for the default
c          value of mmax(=10),
c           lenk = 104 + 14*n,  if |mf| = 2, and
c           lenk = 236 + 14*n,  if |mf| = 1 or 3.
c
c lrw    = integer scalar containing the length of the array rwork,
c          as declared by the user.  (this will be checked by the
c          solver).
c        
c iwork  = an integer work array.  the length of iwork must be at least
c                20 + mmax + lenimp,
c          where lenimp is the length of the user-defined integer work
c          space for the preconditioner routines pset and psol.
c          for the default value of mmax(=10), the length of iwork
c          must be at least 30 + lenimp.
c
c liw    = integer scalar containing the length of the array iwork,
c          as declared by the user.  (this will be checked by the
c          solver.)
c
c iopt   = an integer flag to specify whether or not any optional
c          inputs are being used.
c          iopt=0 means there are no optional inputs.
c                 default values will be used in all cases.
c          iopt .ne. 0 means one or more optional inputs are being
c                      used.  see the optional inputs section below
c                      for more information.
c
c iterm  = output flag.
c          iterm= 1 means maxnorm(sf*f(u)) .le. ftol, where
c                   maxnorm() is the maximum norm function.  u
c                   is probably an approximate root of f.
c          iterm= 2 means the scaled distance between the last two
c                   steps is less than stptol.  u may be an
c                   approximate root of f, but it is also possible
c                   that the algorithm is making very slow progress
c                   and is not near a root, or that stptol is too
c                   large.
c          iterm= 3 means iret=1 was returned from either the dogdrv
c                   or lnsrch module, and that the last global step
c                   failed to reduce norm(f) sufficiently.  either u
c                   is close to a root of f and no more accuracy is
c                   possible, or the finite-difference approximation
c                   to j*v is inaccurate, or stptol is too large.
c                   if the ncfl optional output value (see below) is
c                   close to the nni value, it may be the case that
c                   the krylov iteration is converging very slowly.
c                   in this case, the user may want to use
c                   preconditioning and/or increase the mmax
c                   value (i.e., increase the maximum dimension of the
c                   krylov subspace.)
c          iterm= 4 means that the maximum allowable number of
c                   nonlinear iterations has been reached.  this
c                   is currently set at 200, but may be changed as
c                   an optional input.  see below.
c          iterm= 5 means 5 consecutive steps of length stepmx (the
c                   maximum stepsize limit) have been taken.  either
c                   norm(f) asymptotes from above to a finite value
c                   in some direction, or stepmx is too small.  stepmx
c                   is computed internally as 
c                      stepmx = 1000*max(norm(su*u0),norm(su)),
c                   where u0 is the initial guess, and norm() is the
c                   euclidean norm.  norm(su) here means the euclidean
c                   norm of the n-dimensional array su.  stepmx may
c                   also be set as an optional input.  see below.
c          iterm= 6 means that more than 10 failures occurred when
c                   trying to satisfy the beta-condition in the
c                   linesearch algorithm.  it is likely that the
c                   iteration is making poor progress.
c          iterm= 7 means that there was a breakdown in the krylov
c                   iteration.  this will likely only occur when
c                   the jacobian matrix j or j*(p-inverse) is ill-
c                   conditioned.  if this error return occurs with
c                   mf=2, try either mf=1 or mf=3 instead.
c          iterm= 8 means there was a nonrecoverable error in pset
c                   causing the iteration to halt.
c          iterm= 9 means there was a nonrecoverable error in psol
c                   causing the iteration to halt.
c          iterm=-1 means the mf value was illegal.  the allowable
c                   mf values are 1, 2 or 3.
c          iterm=-2 means the mdif value was illegal.  the allowable
c                   mdif values are 0 and 1.
c          iterm=-3 means the ipflag value was illegal.  the allowable
c                   ipflag values are 0 and 1.
c          iterm=-4 means that an optional input value in iwork was
c                   negative.  only nonnegative values are allowed.
c          iterm=-5 means that the stepmx optional input value was
c                   negative.  stepmx must be .ge. 0.
c          iterm=-6 means that the eta optional input value was
c                   negative.  eta must be .ge. 0.
c          iterm=-7 means that the tau optional input value was
c                   negative.  tau must be .ge. 0.
c          iterm=-8 means that there was insufficient length declared
c                   for the rwork array.  see the lenrw optional output
c                   value for the minimum needed length of rwork.
c          iterm=-9 means that there was insufficient length declared
c                   for the iwork array.  see the leniw optional output
c                   value for the minimum needed length of iwork.
c          iterm=-10 means that the initial guess u did not
c                    satisfy the constraints.
c                  
c pset   = the name of the optional user-supplied subroutine for
c          calculating any matrix data associated with the
c          preconditioner p.  it is to have the form
c           subroutine pset (n, u, savf, su, sf, wk, f, wmp, iwmp, ier)
c           dimension u(*), savf(n), su(n), sf(n), wmp(*), iwmp(*)
c          where n, u, and savf are input.  f is the name of the
c          user-supplied subroutine defining f(u).  any matrix data
c          loaded by this routine must be stored in the real and integer
c          arrays wmp and iwmp.  these arrays are passed unaltered
c          to subroutine psol.  the array u contains the last newton
c          iterate and savf = f(u).  wk is a work array of length n
c          for use in this routine.  subroutine pset should not alter
c          u or savf, and is called infrequently in an effort to use
c          preconditioner data evaluated at an earlier newton iterate.
c          pset must be declared external in the calling program,
c          and is only called when ipflag=1.  when ipflag=0, a
c          dummy routine may be supplied for pset.  on return pset
c          should set the error flag as follows..
c             ier = 0 if pset was successful.
c             ier.ne.0 if an unrecoverable error occurred.  for
c                      example, the preconditioner may be singular.
c                      in this case the nonlinear iteration is halted.
c          note: the lengths of arrays wmp and iwmp must be supplied
c                as optional inputs in the iwork array.  see the
c                optional inputs section below.
c
c          subroutine pset may access user-defined quantities in
c          u(n+1),... if u has length exceeding n in the calling
c          program.  see description of u above.
c
c psol   = the name of the optional user-supplied subroutine for solving
c          the linear system p*x = c.  it uses matrix data loaded
c          in subroutine pset and passed via the arrays wmp and iwmp.
c          it is to have the form
c           subroutine psol (n, u, savf, su, sf, f, jac, wk, wmp, iwmp,
c          *                 x, ier)
c           dimension u(*),savf(n),su(n),sf(n),wk(n),wmp(*),iwmp(*)
c          where n, u, savf, su, sf, wmp and iwmp are input.  x contains
c          the right hand side on input and the solution on exit.  f is the
c          name of the user-supplied routine defining f(u), and jac is
c          the name of the user-supplied routine for calculating j(u)*v.
c          wk is a real work array of length n available for use in
c          psol.  psol should not alter u or savf, and is only called
c          when ipflag=1.  when ipflag=0, a dummy routine may be
c          supplied for psol.  on return, psol should set the error flag
c          ier as follows..
c             ier = 0 if psol was successful.
c             ier.lt.0 if an unrecoverable error occurred.  in this
c                      case the nonlinear iteration is halted.
c             ier.gt.0 if a recoverable error occurred, possibly caused
c                      by the preconditioner being out of date.  a call
c                      to pset is done, and the solution of the current
c                      linear system is retried.
c          note: the lengths of arrays wmp and iwmp must be supplied
c                as optional inputs in the iwork array.  see the
c                optional inputs section below.
c
c          subroutine psol may access user-defined quantities in
c          u(n+1),... if u has length exceeding n in the calling
c          program.  see description of u above.
c
c mdif   = integer scalar containing the method flag for the
c          j*v calculation.
c          mdif=1 means a user-supplied routine jac is used for
c                 calculating the product j*v in the krylov iteration.
c          mdif=0 means an internally generated j*v is used in the
c                 krylov iteration.  this is done by using an f
c                 evaluation and a difference-quotient.
c
c ipflag = integer scalar indicating if preconditioning is being used.
c          ipflag=0 means no preconditioning is used.
c          ipflag=1 means preconditioning is done (on the right only).
c
c
c icflag = integer scalar indicating if constraints violations will
c          be checked for in the approximate solution u.
c          icflag = 0 means no constraint checking.
c          icflag = 1 means constraint checking is done.
c
c icnstr = integer array of length n containing flags indicating
c          which entries in u are to be constrained.
c          if icnstr(i) .gt. 0, then u(i) must be .gt. 0,
c          if icnstr(i) .lt. 0, then u(i) must be .lt. 0, while
c          if icnstr(i) .eq. 0, then u(i) is not constrained.
c          icnstr is not referenced when icflag=0.
c
c  rlx   = real scalar which restricts update to del(u)/u < rlx in cnstrt
c
c  epscon1 real scalar for defining tolerance of linear solve, epsfac
c
c  epscon2 real scalar for defining tolerance of linear solve, epsfac
c          epsfac = epscon1 * min(epscon2, frnm)
c
c  adjf1   real scalar for global strategy test: frnm_new/adjf1>=frnm_old
c
c-----------------------------------------------------------------------
c optional inputs..
c
c the following is a list of the optional inputs provided for in the
c call sequence.  for each such input variable, this table lists its
c name as used in this documentation, its location in the call sequence,
c its meaning, and the default value.  the use of any of these inputs
c requires iopt .ne. 0, and in that case all of these inputs are
c examined.  a value of zero for any of these optional inputs will
c cause the default value to be used.  thus to use a subset of the
c optional inputs, simply preload locations 1 to 8 in iwork to 0, and
c locations 1 to 3 in rwork to 0.0.  then set those of interest
c to nonzero values.
c
c name    location      meaning and default value
c
c mmax    iwork(1)  maximum dimension of krylov subspace in the linear
c                   iteration.  the default value of mmax is 10.
c
c ------  iwork(2)  currently not used.
c
c lenwmp  iwork(3)  length of the real work array wmp for preconditioner
c                   matrix data storage used in subroutines pset and
c                   psol.  the default value of lenwmp is 0.
c
c lenimp  iwork(4)  length of the integer work array iwmp for
c                   preconditioner matrix data storage used in
c                   subroutines pset and psol.  the default value
c                   of lenimp is 0.
c
c iprint  iwork(5)  flag indicating whether optional statistics are
c                   desired.  
c                   iprint=0 means no statistics are printed.
c                            this is the default.
c                   iprint=1 means the nonlinear iteration count,
c                            norm(sf*f(u)) and number of calls to f are
c                            printed for each nonlinear iterate.
c                   iprint=2 includes the statistics for iprint=1, plus
c                            additional statistics from the krylov
c                            inexct, dogleg, or lnsrch strategies.
c                            for the dogleg statistics, tau is the
c                            trust region size, cpl is the length of
c                            of the step to the cauchy point, and gml
c                            is the length of the gmres step.
c                            for the lnsrch statistics, f1 refers to
c                            0.5*norm(sf*f(u))**2 and f1new refers to
c                            0.5*norm(sf*f(unew))**2.
c                   iprint=3 includes the statistics for iprint=2, plus
c                            additional statistics from the krylov
c                            iteration.
c
c iunit   iwork(6)  i/o unit number for optional statistics and error
c                   messages.  the default value of iunit is 6.
c
c iermsg  iwork(7)  flag indicating whether or not error messages are
c                   desired.  iermsg .ne. 0 turns off the printing
c                   of error messages.  the default value of iermsg
c                   is 0.
c
c itermx  iwork(8)  maximum allowable number of nonlinear iterations.
c                   the default value of itermx is 200.
c
c incpset iwork(9)  maximum number of nonlinear iterations before
c                   the preconditioner is reevaluated.
c                   the default value of incpset is 10.
c
c stepmx  rwork(1)  maximum allowable length of a newton step.  the
c                   default value of stepmx is
c                     stepmx = 1000*max(norm(su*u0),norm(su)),
c                   where u0 is the initial guess, and norm() is the
c                   euclidean norm.  norm(su) here means the euclidean
c                   norm of the n-dimensional array su.
c
c eta     rwork(2)  the relative error in computing f(u).  this
c                   value is used in the difference-quotient option
c                   for calculating j*v.  if the user-s f routine has
c                   large errors associated with the calculation of
c                   f(u), then supplying an accurate eta value will
c                   enhance the accuracy of the difference-quotient.
c                   the default value of eta is the machine epsilon.
c
c tau     rwork(3)  real scalar containing the initial trust region
c                   size for the dogleg strategy.  tau is the radius
c                   of the m-dimensional ball (m = dimension of the
c                   krylov subspace) in which the nonlinear function
c                       f1(x) = 0.5*norm( sf*f(u+x) )**2
c                   is well approximated by the local quadratic model
c                       g(x) = 0.5*norm( sf*(f(u) + j(u)*x) )**2
c                   for x in the krylov subspace, norm() the
c                   euclidean norm, and u the current newton iterate.
c                   if tau=zero on input, then an initial value for
c                   it will be chosen internally.
c                   (tau is not used with the linesearch strategy).
c
c-----------------------------------------------------------------------
c optional outputs.
c
c as optional additional output from nksol, the variables listed
c below are quantities related to the performance of nksol
c which are available to the user.  these are communicated by way of
c the work arrays, but also have internal mnemonic names as shown.
c except where stated otherwise, all of these outputs are defined
c on return from nksol.  on an illegal input return (iterm .lt. 0)
c they will be unchanged from their existing values (if any),
c except possibly for lenrw and leniw, and the integer counter values
c will have been reset to 0.
c
c name    location      meaning
c
c nni     iwork(10) the total number of nonlinear iterations.
c
c nli     iwork(11) the total number of linear krylov iterations.
c
c nfe     iwork(12) the total number of calls to the f routine.
c
c nje     iwork(13) the total number of calls to the jac routine.
c
c npe     iwork(14) the total number of calls to pset.
c
c nps     iwork(15) the total number of calls to psol.
c
c ncfl    iwork(16) the total number of times the error tolerance
c                   for the linear krylov iteration was not met
c                   in mmax iterations.
c
c nbcf    iwork(17) the total number of times the beta condition
c                   could not be met in the linesearch algorithm.
c                   the nonlinear iteration is halted if this value
c                   ever exceeds 10.
c
c nb      iwork(18) the total number of backtracks in the linesearch
c                   algorithm, and the number of extra f evaluations
c                   used by the dogleg strategy. nb is calculated as
c                   nb = nfe - nni - (nli - nje) - 1.
c
c lenrw   iwork(19) required minimum length of the rwork array.
c
c leniw   iwork(20) required minimum length of the iwork array.
c
c stepmx  rwork(1)  maximum allowable length of a newton step.  the
c                   default value of stepmx is
c                     stepmx = 1000*max(norm(su*u0),norm(su)),
c                   where u0 is the initial guess, and norm() is the
c                   euclidean norm.  if a value of stepmx .gt. 0 was
c                   input, it is not changed on output.
c
c fnrm    rwork(2)  the scaled norm norm(sf*f(u(k))) for the final
c                   computed iterate u(k).
c
c tau     rwork(3)  real scalar containing the last used trust region
c                   size for the dogleg strategy.
c
c-----------------------------------------------------------------------
c part ii.  optionally replaceable solver routines.
c
c below are descriptions of two routines in the nksol package which
c relate to the stopping criteria for the nonlinear iteration.  either
c routine may be replaced by a user-supplied version.  however, since
c such a replacement may have a major impact on performance, it should
c be done only when absolutely necessary, and only with great caution.
c (note.. the means by which the package version of a routine is
c superseded by the user-s version may be system-dependent.)
c
c (a) snrmf.
c the following subroutine is called to compute the scaled norm of
c f(u) at the current iterate u..
c     subroutine snrmf(n, savf, sf, fnorm)
c     dimension savf(n), sf(n)
c where n and sf are as in the nksol call sequence, savf contains
c f(u) with u the current iterate, and fnorm = norm(sf*f(u)) where
c norm() is the euclidean norm.  the package version of this routines
c uses the scaled norm of sf*f(u) defined in stopping criterion 1 above.
c
c if the user supplies this subroutine, it must return in fnorm
c the user-defined scaled norm of f(u).  if the user-s scaled norm
c of sf*f(u) is simply norm(sf*f(u)), then there is nothing for the
c new routine to compute since fnorm already has the desired value
c on input.
c
c (b) slngth.
c the following subroutine is called to compute the relative scaled
c length of the current step from u to unew = u + s..
c     subroutine slngth(n, u, s, su, rlngth)
c     dimension u(*), s(n), su(n)
c where n, u and su are as in the nksol call sequence, u is the
c current newton iterate, and s contains the current step from u
c to unew.  rlngth is undefined on input.  the package version of
c this routine computes the relative scaled length defined in stopping
c criterion 2 above.
c
c if the user supplies this subroutine, it must return in rlngth
c the user-defined relative scaled length of the step s.
c
c-----------------------------------------------------------------------
c other routines in the package..
c
c nkstop  stopping routine for the nonlinear iteration.
c nkstp0  stopping routine for the initial guess.
c slngth  computes the relative scaled length of a step in maxnorm.
c snrmf   computes the maxnorm of sf*f(u).
c model   driver for approximate solution of newton equations.
c solpk   interfaces to either spiom or spigmr.
c spiom   contains the main iteration loop for arnoldi.
c spigmr  contains the main iteration loop for gmres.
c atv     calculates j*v for a given vector v.
c svrorthog  orthogonalizes a new vector against older basis vectors.
c shefa   computes an lu decomposition of a hessenberg matrix.
c shesl   solves a hessenberg linear system, using lu factors.
c sheqr   computes a qr decomposition of a hessenberg matrix.
c shels   solves a hessenberg least-squares system, using qr factors.
c epsilon computes machine epsilon in double precision.
c vnormnk computes a scaled euclidean norm of a vector.
c dogdrv  driver for the dogleg strategy.
c dogstp  computes the dogleg step for the current trust region.
c trgupd  decides if step is successfull, and updates trust region size.
c lnsrch  driver for linesearch backtracking strategy.
c cnstrt0 checks for constraint violations in initial u.
c cnstrt  checks for constraint violations in dogleg strategy.
c errgen  outputs error messages.
c snrm2   computes the l2 norm of a vector.
c scopy   copies a vector to another vector.
c sswap   swaps two vectors.
c sdot    computes the dot product of two vectors.
c saxpy   adds a scalar times a vector to a vector.
c sscal   scales a vector by a scalar.
c isamax  finds the element of largest magnitude in a vector.
c note.. snrm2, scopy, sswap, sdot, saxpy, sscal, and isamax are from
c the blas collection (basic linear algebra modules).
c intrinsic fortran routines.. abs, sqrt, min, max, float, sign.
c-----------------------------------------------------------------------

c local variables
      implicit none
      integer n, lrw, iwork, liw, iopt, iterm, mf, mdif, ipflag
      integer locwmp, locimp, iersl, kmp, mmax, methn, methk, ipflg
      integer mfdif, nfe, nje, nni, nli, npe, nps, ncfl, nbcf
      integer iprint, iunit, iermsg, ipcur, nnipset, incpset, ierr
      integer lenwmp, lenimp, itermx, i, lup, lx, luprv, lfprv, lwm
      integer lenk, lenwm, lenrw, leniw, liwm, leniwm, nbcfmx, iter
c+ icntnu change
c      integer iret, ncscmx
      integer iret, ncscmx, icntnu
c- icntnu change
      external f, jac, pset, psol
      real savf,u,rwork,su,sf
      real stptol,epsmch,fnrm,f1nrm,f1nrmp,unrm
      real ftol,stepmx,epsilon,vnormnk,epsfac,sunrm,snrm2,tau
      dimension savf(n),u(*),rwork(lrw),su(n),sf(n)
      dimension iwork(liw)
      real rlx, epscon1, epscon2, adjf1
      integer icflag, icnstr, ivar
      dimension icnstr(n)
c
      real zero,one,two,three
      logical mxtkn
       Use(Cdv)  # comnfe, pandfitertag
c+pnb
C
C Type declaration for the additional diagnostics set up to inquire ----
C during the running of nksol. -----------------------------------------
C
      character*80 msgjm
      integer nrcv, ierrjm, ijmgetmr
c-pnb
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real eps, sqteta, rhom
      common /nks001/ eps, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
c-----------------------------------------------------------------------
c     nks003 common block.
c-----------------------------------------------------------------------
      real pthrsh
      common /nks003/ pthrsh, ipcur, nnipset, incpset
c
      save
      zero=0.
      one=1.0
      two=2.0
      three=3.0
c-----------------------------------------------------------------------
c zero counters in common block nks001.
c-----------------------------------------------------------------------
      nni = 0
      nli = 0
      npe = 0
      nps = 0
      ncfl = 0
      nbcf = 0
      nfe = 0
      nje = 0
c-----------------------------------------------------------------------
c load values in nks002 common block.
c-----------------------------------------------------------------------
      iprint = 0
      iunit = 6
      iermsg = 0
c+ icntnu change
c-----------------------------------------------------------------------
c load pthrsh value in nks003 common block if this is not a continuation
c call (i.e., icntnu = 0).
c-----------------------------------------------------------------------
c      pthrsh = two
      if (icntnu .eq. 0) then
	pthrsh = two
      else	
c set pthrsh = 0
	pthrsh = zero
c set ipcur = 0 to indicate that the preconditioner is from an earlier
c call.
	ipcur = 0
c set nnipset = 0 for a continuation call.
	nnipset = 0
      endif
c- icntnu change
c=======================================================================
c block a: initialization section.
c=======================================================================
c-----------------------------------------------------------------------
c     compute machine epsilon.
c-----------------------------------------------------------------------
      epsmch = epsilon(zero)
c-----------------------------------------------------------------------
c     initialize parameters.  check for illegal input.
c-----------------------------------------------------------------------

c:MVU 15-jan-2020      if ( (iabs(mf) .lt. 1) .or. (iabs(mf) .gt. 3) ) then
      if ( (iabs(mf) .lt. 1) .or. (iabs(mf) .gt. 4) ) then
c illegal value of method flag mf.
        iterm = -1
        ierr = 10
        call errgen(ierr,zero,zero,0,0)
        go to 500
        endif
c-----------------------------------------------------------------------
c load common variables methn and methk.
c methn =  0, no global strategy (inexact newton method used).
c          1, dogleg global strategy used.
c          2, linesearch backtracking global strategy used.
c methk =  1, arnoldi algorithm used.
c          2, gmres algorithm used.
c-----------------------------------------------------------------------
      if (mf .eq. 1) then
          methn = 1
          methk = 2
        elseif (mf .eq. -1) then
          methn = 0
          methk = 2
ccc MVU: 15-jan-2020					  
c        elseif (mf .ge. 2) then
c          methn = 2
c          methk = mf - 1
        elseif (mf .eq. 2) then
          methn = 2
          methk = 1
        elseif (mf .eq. 3) then
          methn = 2
          methk = 2
        elseif (mf .eq. 4) then
c........ using direct solver
          methn = 2
          methk = 3
ccc MVU: 15-jan-2020
        elseif (mf .le. -2) then
          methn = 0
          methk = -mf - 1
       endif
      if ( (mdif .lt. 0) .or. (mdif .gt. 1) ) then
c illegal value of method flag mdif.
        iterm = -2
        ierr = 20
        call errgen(ierr,zero,zero,0,0)
        go to 500
        endif
c load common variable mfdif.
      mfdif = 2 - mdif
      if ( (ipflag .lt. 0) .or. (ipflag .gt. 1) ) then
c illegal value of ipflag.
        iterm = -3
        ierr = 30
        call errgen(ierr,zero,zero,0,0)
        go to 500
        endif
c load common variable ipflg.
      ipflg = ipflag
c
      lenwmp = 0
      lenimp = 0
      mmax = 10
      kmp = mmax
      stepmx = zero
      itermx = 0
      sqteta = zero
      tau = zero
      incpset = 10
c check for optional inputs.
      if (iopt .ne. 0) then
c first check for illegal values.
        do 10 i = 1,9
          if (iwork(i) .lt. 0) then
            iterm = -4
            ierr = 40
            call errgen(ierr,zero,zero,i,iwork(i))
            go to 500
            endif
 10       continue
        if (rwork(1) .lt. zero) then
          iterm = -5
          ierr = 50
          call errgen(ierr,rwork(i),zero,i,0)
          go to 500
          endif
        if (rwork(2) .lt. zero) then
          iterm = -6
          ierr = 50
          call errgen(ierr,rwork(i),zero,i,0)
          go to 500
          endif
        if (rwork(3) .lt. zero) then
          iterm = -7
          ierr = 50
          call errgen(ierr,rwork(i),zero,i,0)
          go to 500
          endif
        if (iwork(1) .gt. 0) mmax = iwork(1)
        if (iwork(2) .gt. 0) kmp  = iwork(2)
        if (iwork(3) .gt. 0) lenwmp = iwork(3)
        if (iwork(4) .gt. 0) lenimp = iwork(4)
        if (iwork(5) .gt. 0) iprint = iwork(5)
        if (iwork(6) .gt. 0) iunit = iwork(6)
        if (iwork(7) .gt. 0) iermsg = iwork(7)
        if (iwork(8) .gt. 0) itermx = iwork(8)
        if (iwork(9) .gt. 0) incpset = iwork(9)
        if (rwork(1) .gt. zero) stepmx = rwork(1)
        if (rwork(2) .gt. zero) sqteta = sqrt(rwork(2))
        if (rwork(3) .gt. zero) tau = rwork(3)
        endif
      mmax = min0(mmax,n)
      kmp = mmax
      lup = 4
      lx = lup + n
      luprv = lx + n
      lfprv = luprv + n
      lwm = lfprv + n
      if (methk .eq. 1) then
          lenk = 1 + 3*n + 3 + n*mmax + n + mmax**2
        else if (methk .eq. 2) then
          lenk = 1 + 3*n + 3 + n*mmax + n + 2 + 2*mmax*(mmax+1) + mmax
ccc MVU: 15-jan-2020
        else if (methk .eq. 3) then
          lenk = 1 + 3*n + 3 + n*mmax + n + 2 + 2*mmax*(mmax+1) + mmax
ccc MVU: 15-jan-2020
        endif
      lenwm = lenk + lenwmp
      if (lenwmp .ne. 0) then
          locwmp = lenk + 1
        else
          locwmp = lenk
        endif
      lenrw = lwm + lenwm
      if (lrw .lt. lenrw) then
c insufficient storage in rwork.
        iterm = -8
        ierr = 100
        call errgen(ierr,zero,zero,lrw,lenrw)
        go to 500
        endif
c
      leniw = 20 + mmax + lenimp
      if (liw .lt. leniw) then
c insufficient storage in iwork
        iterm = -9
        ierr = 110
        call errgen(ierr,zero,zero,liw,leniw)
        go to 500
        endif
      liwm = 21
      leniwm = mmax + lenimp
      if (lenimp .ne. 0) then
          locimp = mmax + 1
        else
          locimp = mmax
        endif
      if (itermx .eq. 0) itermx = 200
      nbcfmx = 10
      if (iprint.gt.1) write(*,*)'0) sptol,epsmch', stptol,epsmch  
      if (stptol .eq. 0.0) stptol = epsmch**(2.0/3.0)
      if (iprint.gt.1) write(*,*)'1) sptol', stptol
      if (stepmx .eq. zero) then
        unrm = vnormnk(n,u,su)
        sunrm = snrm2(n,su,1)
        stepmx = 1.e3*max(unrm,sunrm)
        endif
      if (ftol .eq. zero) ftol = epsmch**(one/three)
      if (sqteta .eq. zero) sqteta = sqrt(epsmch)
c load initial trust region size tau.
      if (tau .le. zero) then
        tau = -one
        endif
      iter = 0
c if icflag .ne. 0, then call cnstrt0 to check that u satisfies
c the constraints.
      if (icflag .ne. 0) then
         call cnstrt0(n, u, icnstr, rlx, iret, ivar)
         if (iret .gt. 0) then
c           constraint violation occurred.
            iterm = -10
            ierr = 130
            call errgen(ierr,zero,zero,iret,ivar)
            go to 500
         endif
      endif
c-----------------------------------------------------------------------
c     call subroutine f to evaluate f at the initial guess.
c-----------------------------------------------------------------------
      call f(n, u, savf)
      nfe = nfe + 1
      comnfe = comnfe + 1
      fnrm = vnormnk(n,savf,sf)
      f1nrm = fnrm*fnrm/two
      if (iprint .ge. 1) write(iunit,400) TRIM(pandfitertag), iter,fnrm,nfe
c-----------------------------------------------------------------------
c     test to see if initial conditions satisfy stopping criterion.
c-----------------------------------------------------------------------
      call nkstp0(n, savf, sf, ftol, fnrm, iterm)
      if (iterm .ne. 0) go to 500
      ncscmx = 0
c=======================================================================
c block b: iteration section.
c=======================================================================
 100  continue
c+pnb
C-----------------------------------------------------------------------
C This section is here to obtain information while the solver ----------
C is running,  The keyword is:
C              kaboom  returns to the user prompt, this option is useful
C                      when the solution failes to converge.
C This section was added by PNB 4/2/92
C-----------------------------------------------------------------------
ccc      ierrjm = ijmgetmr(msgjm,80,1,nrcv)
ccc      if (ierrjm .eq. 0) then
ccc         if (msgjm(1:nrcv).eq.'kaboom' .or. msgjm(1:nrcv).eq.'k')then
ccc            call xerrab("")
ccc         endif
ccc      endif
C-----------------------------------------------------------------------
c-pnb
      iter = iter + 1
c-----------------------------------------------------------------------
c     call model to calculate the approximate newton step using the
c     appropriate krylov algorithm.
c-----------------------------------------------------------------------
c      epsfac = 0.5e0**(iter)
c     epsfac = min(1.e-5,0.01e0*fnrm)
      epsfac = epscon1 * min(epscon2, fnrm)
      eps = (epsmch + epsfac)*fnrm
      call model(n,rwork(lwm),lenwm,iwork(liwm),leniwm,u,savf,rwork(lx),
     *            f,jac,su,sf,pset,psol)
      if (iersl .ne. 0) then
        iterm = iersl
        go to 500
        endif
c-----------------------------------------------------------------------
c     call dogdrv or lnsrch to get an acceptable step x.
c-----------------------------------------------------------------------
      if (methn .eq. 0) then
          call inexct(n,u,savf,f1nrm,rwork(lx),su,sf,stepmx,stptol,
     *                iret,rwork(lup),f1nrmp,mxtkn,f,icflag,
     *                  icnstr, rlx)
        elseif (methn .eq. 1) then
          call dogdrv (n,rwork(lwm),lenwm,iwork(liwm),leniwm,u,savf,
     *                 f1nrm,rwork(lx),su,sf,stepmx,stptol,tau,iret,
     *                 rwork(luprv),rwork(lfprv),rwork(lup),f1nrmp,
     *                 mxtkn,f,jac,psol,icflag,icnstr,rlx)
          if (iersl .ne. 0) then
            iterm = iersl
            go to 500
            endif
        else
          call lnsrch(n,u,savf,f1nrm,rwork(lx),su,sf,stepmx,stptol,
     *                iret,rwork(lup),f1nrmp,mxtkn,f,jac,icflag,icnstr,
     *                                                       rlx,adjf1)
          if (nbcf .gt. nbcfmx) then
            iterm = 6
            ierr = 120
            call errgen(ierr,zero,zero,nbcf,nbcfmx)
            go to 500
            endif
        endif
c-----------------------------------------------------------------------
c     call nkstop to check if tolerances are met.
c-----------------------------------------------------------------------
      fnrm = sqrt(two*f1nrmp)
      # we add here a trap for aborted exmain when 'stop' is called. Added by J.Guterl
      if (exmain_aborted) call xerrab('exmain aborted...')
      call nkstop(n,u,rwork(lup),savf,fnrm,su,sf,stptol,rwork(lx),
     *            ftol,iret,iter,itermx,mxtkn,ncscmx,iterm)
c
      do 300 i = 1,n
        u(i) = rwork(i+lup-1)
 300    continue
      f1nrm = f1nrmp
      if (iprint .ge. 1) write(iunit,400) TRIM(pandfitertag), iter,fnrm,nfe
 400  format(' ', a,i4,' fnrm= ',g26.16,' nfe= ',i6)
      if (iterm .eq. 0) go to 100
c-----------------------------------------------------------------------
c     load optional outputs into iwork array and return.
c-----------------------------------------------------------------------
 500  continue
      call infgen (iterm,zero,zero,0,0)
      iwork(10) = nni
      iwork(11) = nli
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = npe
      iwork(15) = nps
      iwork(16) = ncfl
      iwork(17) = nbcf
      iwork(18) = nfe - nni - (nli - nje) - 1
      iwork(19) = lenrw
      iwork(20) = leniw
      rwork(1) = stepmx
      rwork(2) = fnrm
      rwork(3) = tau
      return
c----------------------- end of subroutine nksol -----------------------
      end
      subroutine nkstop(n, u, unew, savf, fnorm, su, sf, stptol, wk,
     *    ftol, iret, iter, itermx, mxtkn, ncscmx, iterm)
c-----------------------------------------------------------------------
c this routine decides whether to terminate nksol based on one of
c several stopping criteria described below.  (see description of
c iterm.)
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   u      = real array containing the previous approximate root.
c
c   unew   = real array containing the current approximate root.
c
c   savf   = real array containing the n-vector f(unew).
c
c   fnorm  = real scalar containing norm(sf*f(unew)), where norm()
c            denotes the euclidean norm.  (currently, not used, but
c            available if the user supplies his own stopping routine.)
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c   stptol = user-supplied tolerance on the minimum step length
c            s = unew-u.
c
c   wk     = real array of length n used for work space by nkstop.
c
c   ftol   = user-supplied tolerance on norm(sf*f(unew)).
c
c   iret   = output flag from either the inexct, dogdrv or lnsrch
c            module.
c
c   iter   = current nonlinear iteration count.
c
c   itermx = maximum allowable number of nonlinear iterations.
c
c   mxtkn  = logical flag returned from either lnsrch or dogdrv
c            module indicating if the last step taken was of maximum
c            length stepmx.
c
c   ncscmx = number of consecutive steps on which the stepsize
c            was of maximum length stepmx.
c
c on return
c
c   ncscmx = number of consecutive steps on which the maximum
c            allowable stepsize length stepmx was taken.  if
c            mxtkn=.true., then ncscmx has been incremented by
c            one, and ncscmx is reset to zero if mxtkn=.false.
c
c   iterm  = output flag.
c            iterm=0 means no termination criteria satisfied,
c                    and the iteration continues.
c            iterm=1 means maxnorm(sf*f(unew)) .le. ftol, where
c                    maxnorm() is the maximum norm function.  unew
c                    is probably an approximate root of f.
c            iterm=2 means the scaled distance between the last two
c                    steps is less than stptol.  unew may be an
c                    approximate root of f, but it is also possible
c                    that the algorithm is making very slow progress
c                    and is not near a root, or that stptol is too
c                    large.
c            iterm=3 means iret=1 was returned from either the dogdrv
c                    or lnsrch module, and that the last global step
c                    failed to reduce norm(f) sufficiently.  either u
c                    is close to a root of f and no more accuracy is
c                    possible, or the finite-difference approximation
c                    to j*v is inaccurate, or stptol is too large.
c            iterm=4 means that the maximum allowable number of
c                    nonlinear iterations has been exceeded.
c            iterm=5 means 5 consecutive steps of length stepmx (the
c                    maximum stepsize limit) have been taken.  either
c                    norm(f) asymptotes from above to a finite value
c                    in some direction, or stepmx is too small.
c                    
c-----------------------------------------------------------------------
      implicit none
       Use(Cdv)
      integer n, iret, iter, itermx, ncscmx, iterm, locwmp, locimp
      integer iersl, kmp, mmax, methn, methk, ipflg, mfdif, nfe, nje
      integer nni, nli, npe, nps, ncfl, nbcf, ipcur, nnipset
      integer incpset, i
      real u,unew,savf,fnorm,su,sf,stptol,ftol
      real fmax,rlngth,wk,thsnd,two
      dimension u(*),unew(n),savf(n),su(n),sf(n),wk(n)
      logical mxtkn
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real eps, sqteta, rhom
      common /nks001/ eps, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c-----------------------------------------------------------------------
c     nks003 common block.
c-----------------------------------------------------------------------
      real pthrsh
      common /nks003/ pthrsh, ipcur, nnipset, incpset
      save
c      data thsnd/1000.0e0/
      data thsnd/1.0e10/
      data two/2.0e0/
c
      iterm = 0
c check for error return from inexct, dogdrv or lnsrch.
      if (iret .eq. 1) then
c if preconditioner is out of date, then reevaluate, and continue
c the iteration.
         if (ipcur .eq. 0) then
            pthrsh = two
            go to 999
            endif
        iterm = 3
        go to 999
        endif
c check tolerance on norm(sf*f(unew)).
      fmax = fnorm
      call snrmf(n, savf, sf, fmax)
      if (fmax .le. ftol) then
        iterm = 1
        go to 999
        endif
c check for scaled distance between last two steps too small.
      do 200 i = 1,n
        wk(i) = unew(i) - u(i)
 200    continue
      call slngth(n, u, wk, su, rlngth)
      if (rlngth .le. stptol) then
c if preconditioner is out of date, then reevaluate, and continue
c the iteration.
         if (ipcur .eq. 0) then
            pthrsh = two
            go to 999
            endif
         endif
      if (rlngth .le. stptol) then
        iterm = 2
        go to 999
        endif
c check for maximum number of iterates exceeded.
      if (iter .ge. itermx) then
        iterm = 4
        go to 999
        endif
c check for consecutive number of steps taken of size stepmx.
c if mstkn=.false., set ncscmx = 0.
      if (mxtkn) then
          ncscmx = ncscmx + 1
        else
          ncscmx = 0
        endif
      if (ncscmx .eq. 5) iterm = 5
c load threshold for re-evaluating preconditioner.
c      pthrsh = fmax/(thsnd*ftol) + rlngth
ccc Modified to prevent new Jac eval for velocities; change;step TDR 5/26/99
ccc      pthrsh = rlngth
      pthrsh = 1e-20*rlngth
      if ((nni - nnipset) .ge. incpset) pthrsh = two
c
 999  continue
      return
c----------------------- end of subroutine nkstop ----------------------
      end
      subroutine slngth(n, u, s, su, rlngth)
c-----------------------------------------------------------------------
c this routine calculates the scaled steplength of the current step s.
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   u      = real array containing the previous approximate root.
c
c   s      = real array containing the current step from u to unew.
c
c   su     = real array containing scaling factors for u.
c
c on return
c
c   rlngth = real scalar containing the scaled steplength.
c-----------------------------------------------------------------------
      real u, s, su, one, rlngth, temp, zero
      integer i, n
      dimension u(n), s(n), su(n)
c
      zero=0.0
      one=1.0
      rlngth = zero
      do 10 i = 1,n
        temp = one/su(i)
        temp = max(abs(u(i)),temp)
        rlngth = max(abs(s(i))/temp,rlngth)
 10     continue
      return
c----------------------- end of subroutine slngth ----------------------
      end
      subroutine snrmf(n, savf, sf, fnorm)
c-----------------------------------------------------------------------
c this routine calculates the scaled norm of f(u) for the current
c iterate u.
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            savf and sf.
c
c   savf   = real array containing the n-vector f(u).
c
c   sf     = real array containing scaling factors for f(u).
c
c   fnorm  = real scalar containing norm(sf*f(unew)), where norm()
c            denotes the euclidean norm.  (currently, not used, but
c            available if the user supplies his own version of this
c            routine.)
c
c on return
c
c   fnorm = real scalar containing the scaled norm of f(u).
c-----------------------------------------------------------------------
      implicit none
      integer n, i
      real savf,sf,zero,temp,fnorm,fmax
      dimension savf(n),sf(n)
      save
      data zero/0.0e0/
c
      fmax = zero
      do 10 i = 1,n
        temp = sf(i)*abs(savf(i))
        fmax = max(temp,fmax)
 10     continue
      fnorm = fmax
      return
c----------------------- end of subroutine snrmf -----------------------
      end
      subroutine nkstp0(n, savf, sf, ftol, fnorm, iterm)
c-----------------------------------------------------------------------
c this routine decides whether to terminate nksol at iteration zero
c because u0 is an approximate root of f(u).
c
c on entry
c
c   n     = size of the nonlinear system, and the length of arrays
c           savf and sf.
c
c   savf  = real array containing the n-vector f(u0).
c
c   ftol  = real scalar containing the user-s stopping tolerance.
c
c   fnorm = real scalar containing norm(sf*f(u0)), where norm()
c           denotes the euclidean norm.  (currently, not used, but
c           available if the user supplies his own stopping routine.)
c
c on return
c
c   iterm = output flag.
c           iterm=0 means u0 is not an approximate root of f.
c           iterm=1 means u0 is an approximate root of f.
c
c-----------------------------------------------------------------------
      implicit none
      integer n, iterm
      real savf,sf,ftol,fmax,p01,fnorm
      dimension savf(n),sf(n)
      save
      data p01/1.e-2/
c
      fmax = fnorm
      call snrmf(n, savf, sf, fmax)
      iterm = 0
      if (fmax .le. p01*ftol) iterm = 1
c
      return
c----------------------- end of subroutine nkstp0 ----------------------
      end
      real function vnormnk (n, v, s)
c-----------------------------------------------------------------------
c this function routine computes the scaled euclidean norm
c of the vector of length n contained in the array v, with scale factors
c contained in the array s of length n..
c   vnormnk = sqrt( sum( v(i)*s(i) )**2 )
c-----------------------------------------------------------------------
      integer n, i
      real v, s, sum
      dimension v(n), s(n)
      sum = 0.0e0
      do 10 i = 1,n
 10     sum = sum + (v(i)*s(i))**2
      vnormnk = sqrt(sum)
      return
c----------------------- end of function vnormnk -------------------------
      end
      subroutine model(n, wm, lenwm, iwm, leniwm, u, savf, x, f, jac,
     *                 su, sf, pset, psol)
c-----------------------------------------------------------------------
c this routine interfaces to subroutine solpk for the approximate 
c solution of the newton equations in the newton iteration.
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   u      = real array containing the previous approximate root.
c
c   savf   = real array containing the n-vector f(u).
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c   wm     = real work space containing data for the krylov algorithm
c            (krylov basis vectors, hessenberg matrix, user-supplied
c            work space for preconditioning, etc.).
c
c   lenwm  = length of work array wm.
c
c   iwm    = integer work space for the krylov algorithm (pivot
c            information for the lu-factorization of the hessenberg
c            matrix, user-supplied work space for preconditioning,
c            etc.).
c
c   leniwm = length of work array iwm.
c
c on return
c
c   x      = the approximate solution vector of the newton equations.
c
c   iersl  = output flag (in common)..
c            iersl = 0 means no trouble occured.
c            iersl = 7 means the krylov solver suffered a breakdown, and
c                      so the solution x is undefined.
c            iersl = 8 means there was a nonrecoverable error in pset
c                      and the nonlinear iteration is halted.
c            iersl = 9 means there was a nonrecoverable error in psol
c                      and the nonlinear iteration is halted.
c
c-----------------------------------------------------------------------
      implicit none
       Use(Cdv)
      integer n, lenwm, leniwm, locwmp, locimp, iersl, kmp, mmax
      integer methn, methk, ipflg, mfdif, nfe, nje, nni, nli, npe
      integer nps, ncfl, nbcf, ipcur, nnipset, incpset, ier
      external f, jac, pset, psol
      integer iwm
      integer i
      real wm, x, su, sf, savf, u, onept5
      dimension wm(lenwm), iwm(leniwm), x(n)
      dimension u(*), savf(n), su(n), sf(n)
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real eps, sqteta, rhom
      common /nks001/ eps, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c-----------------------------------------------------------------------
c     nks003 common block.
c-----------------------------------------------------------------------
      real pthrsh
      common /nks003/ pthrsh, ipcur, nnipset, incpset
c
      save
      data onept5/1.5e0/
c-----------------------------------------------------------------------
c     call user-supplied routine pset to load preconditioner matrix
c     data if threshold achieved.  set ipcur=1 if p is re-evaluated.
c-----------------------------------------------------------------------
 10   continue
      ipcur = 0
      if ( (pthrsh .gt. onept5) .and. (ipflg .ne. 0) ) then
        ier = 0
        call pset (n, u, savf, su, sf, x, f, wm(locwmp), iwm(locimp),
     *             ier)
        npe = npe + 1
        ipcur = 1
        nnipset = nni
        if (ier .ne. 0) then
          iersl = 8
          return
          endif
        endif
c-----------------------------------------------------------------------
c     load x with -f(u).
c-----------------------------------------------------------------------
      do 100 i = 1,n
 100    x(i) = -savf(i)
c-----------------------------------------------------------------------
c     call solpk to solve j*x = -f using the appropriate krylov
c     algorithm.
c-----------------------------------------------------------------------
      call solpk (n,wm,lenwm,iwm,leniwm,u,savf,x,su,sf,f,jac,psol)
      if (iersl .lt. 0) then
c nonrecoverable error from psol.  set iersl and return.
         iersl = 9
         return
         endif
      if ( (iersl .gt. 0) .and. (ipflg .ne. 0) ) then
        if (ipcur .eq. 0) go to 10
        endif
c iersl=1 on return from solpk means there was a breakdown in the
c krylov iteration.  set iersl=7 to halt iteration and return.
      if (iersl .eq. 1) iersl = 7
c
      return
c----------------------- end of subroutine model -----------------------
      end
      subroutine solpk (n, wm, lenwm, iwm, leniwm, u, savf, x, su, sf,
     *                  f, jac, psol)
      implicit none
       Use(Cdv)
      integer lenwm, leniwm, locwmp, locimp, iersl, kmp, mmax
      integer methn, methk, ipflg, mfdif, nfe, nje, nni, nli, npe
      integer nps, ncfl, nbcf, iwk, npsl, mmaxp1, ihsv, iq, mgmr
      external f, jac, psol
      integer n, iwm
      real wm, x, su, sf, savf, u
      dimension u(*), savf(n), su(n), sf(n), wm(lenwm), iwm(leniwm),
     *          x(n)
c-----------------------------------------------------------------------
c this routine interfaces to either subroutine spiom or spigmr for
c the solution of the linear system j*x = -f arising from a newton
c iteration.
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   u      = real array containing the previous approximate root.
c
c   savf   = real array containing the n-vector f(u).
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c   wm     = real work space containing data for the krylov algorithm
c            (krylov basis vectors, hessenberg matrix, user-supplied
c            work space for preconditioning, etc.).
c
c   lenwm  = length of work array wm.
c
c   iwm    = integer work space for the krylov algorithm (pivot
c            information for the lu-factorization of the hessenberg
c            matrix, user-supplied work space for preconditioning,
c            etc.).
c
c   leniwm = length of work array iwm.
c
c on return
c
c   x      = the approximate solution vector of the newton equations.
c
c   iersl  = output flag (in common)..
c            iersl = 0 means no trouble occured.
c            iersl = 1 means the krylov solver suffered a breakdown, and
c                      so the solution x is undefined.
c            iersl =-1 means there was a nonrecoverable error in psol.
c                      this forces the iteration to halt.
c
c this routine also uses the common variables
c eps, locwmp, locimp, kmp, mmax, methk, nni, nli, npe, nps, ncfn, ncfl
c-----------------------------------------------------------------------
c methods used so far:
c 1. methk = 1      ------> spiom
c 2. methk = 2      ------> gmres
ccc MVU: 15-jan-2020
c 3. methk = 3      ------> direct
ccc MVU: 15-jan-2020
c-----------------------------------------------------------------------
      integer i, ib, iflag, ihes, iv, miom
      real bnrm
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real eps, sqteta, rhom
      common /nks001/ eps, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
      save
c
      iersl = 0
      go to (100,200,300), methk
c-----------------------------------------------------------------------
c use the spiom algorithm to solve the linear system j*x = -f.
c-----------------------------------------------------------------------
 100  continue
      iv = 3
      ib = iv + n*mmax
      ihes = ib + n
      iwk = ihes + mmax*mmax
      do 110 i = 1,n
 110    wm(i+ib-1) = x(i)
      call spiom (n, u, savf, wm(ib), su, sf, mmax, kmp, eps,
     *   f, jac, psol, npsl, x, wm(iv), wm(ihes), iwm, miom,
     *   wm(locwmp), iwm(locimp), wm(iwk), ipflg, iflag, rhom)
      nni = nni + 1
      nli = nli + miom
      nps = nps + npsl
      if (iflag .ne. 0) ncfl = ncfl + 1
      if (iflag .ge. 2) iersl = 1
      if (iflag .lt. 0) iersl = -1
      return
c-----------------------------------------------------------------------
c use the spigmr algorithm to solve the linear system j*x = -f.
c-----------------------------------------------------------------------
 200  continue
      mmaxp1 = mmax + 1
      iv = 3
      ib = iv + n*mmax
      ihes = ib + n + 1
      ihsv = ihes + mmax*(mmaxp1+1) + 1
      iwk = ihsv + mmax*mmaxp1
      iq = iwk + n
      do 210 i = 1,n
 210     wm(i+ib-1) = x(i)
      call spigmr (n, u, savf, wm(ib), su, sf, mmax, mmaxp1, kmp,
     *   eps, f, jac, psol, npsl, x, wm(iv), wm(ihes), wm(iq),
     *   wm(ihsv), mgmr, wm(locwmp), iwm(locimp), wm(iwk), methn,
     *   bnrm, ipflg, iflag, rhom)
      nni = nni + 1
      nli = nli + mgmr
      nps = nps + npsl
      if (iflag .ne. 0) ncfl = ncfl + 1
      if (iflag .ge. 2) iersl = 1
      if (iflag .lt. 0) iersl = -1
      if (iersl .eq. 0) then
        call scopy (mmax, wm(ib), 1, wm(ihes), 1)
        wm(1) = bnrm
        iwm(1) = mgmr
        endif
      return

ccc MVU: 15-jan-2020
c-----------------------------------------------------------------------
c use the direct algorithm to solve the linear system j*x = -f.
c-----------------------------------------------------------------------
 300  continue
      mmaxp1 = mmax + 1
      iv = 3
      ib = iv + n*mmax
      ihes = ib + n + 1
      ihsv = ihes + mmax*(mmaxp1+1) + 1
      iwk = ihsv + mmax*mmaxp1
      iq = iwk + n
      do 310 i = 1,n
 310     wm(i+ib-1) = x(i)

ccc - this is direct solver
      do i=1,n
       x(i)=x(i)*sf(i)/su(i)
      enddo
      call psol(n, u, savf, su, sf, f, jac, wm(iwk), wm(locwmp),
     *          iwm(locimp), x, iflag)
      npsl = 1
      bnrm = 0.
      mgmr = 1
      rhom = 0.

ccc - this way mfnksol=4 would identically match mfnksol=3
c$$$      call spigmr (n, u, savf, wm(ib), su, sf, mmax, mmaxp1, kmp,
c$$$     *   eps, f, jac, psol, npsl, x, wm(iv), wm(ihes), wm(iq),
c$$$     *   wm(ihsv), mgmr, wm(locwmp), iwm(locimp), wm(iwk), methn,
c$$$     *   bnrm, ipflg, iflag, rhom)

      nni = nni + 1
      nli = nli + mgmr
      nps = nps + npsl
      if (iflag .ne. 0) ncfl = ncfl + 1
      if (iflag .ge. 2) iersl = 1
      if (iflag .lt. 0) iersl = -1
      if (iersl .eq. 0) then
        call scopy (mmax, wm(ib), 1, wm(ihes), 1)
        wm(1) = bnrm
        iwm(1) = mgmr
        endif
      return
ccc MVU: 15-jan-2020




c----------------------- end of subroutine solpk -----------------------
      end
      subroutine spiom (n, u, savf, b, su, sf, mmax, iomp, eps,
     *            f, jac, psol, npsl, x, v, hes, ipvt, miom, wmp,
     *            iwmp, wk, ipflg, iflag, rho)
      implicit none
      integer npsl, iwmp, ipflg, iprint, iunit, iermsg, ier
      external f, jac, psol
      integer n, mmax, iomp, ipvt, miom, iflag
      real u, savf, b, su, sf, x, eps, v, hes, wmp, wk
      dimension  u(*), savf(n), b(n), su(n), sf(n), x(n), v(n,mmax+1),
     *           hes(mmax,mmax), wmp(*), wk(n), iwmp(*), ipvt(mmax)
c-----------------------------------------------------------------------
c this routine solves the linear system a * x = b using a scaled
c preconditioned version of the incomplete orthogonalization method.
c an initial guess of x = 0 is assumed.
c-----------------------------------------------------------------------
c
c on entry
c
c   n      = problem size, passed to f, jac and psol.  also, the order
c            of the matrix a, and the lengths of the vectors u, savf,
c            b, su, sf, and x.
c
c   u      = array containing current approximate solution to f(u) = 0.
c
c   savf   = array containing current value of f(u).
c
c   b      = the right hand side of the system a*x = b.
c            b is also used as work space when computing the
c            final approximation.
c
c   su     = array of length n containing scale factors for u.
c
c   sf     = array of length n containing scale factors for f.
c
c   mmax   = the maximum allowable order of the matrix hes.
c
c   iomp   = the number of previous vectors the new vector vnew
c            must be made orthogonal to.  iomp .le. mmax.
c
c   eps    = tolerance on residuals b-a*x in scaled euclidean norm.
c
c   wmp    = real work array used by preconditioner psol.
c
c   iwmp   = integer work array used by preconditioner psol.
c
c   ipflg  = integer flag
c            ipflg=0 means no preconditioning is done.
c            ipflg.ne.0 means preconditioning is being performed.
c
c on return
c
c   x      = the final computed approximation to the solution
c            of the system a*x = b.
c
c   rho    = the scaled euclidean norm of the residual associated
c            with x.
c
c   v      = the n x (miom+1) array containing the miom
c            orthogonal vectors v(*,1) to v(*,miom).
c
c   hes    = the lu factorization of the miom x miom upper
c            hessenberg matrix whose entries are the
c            scaled inner products of a*v(*,k) and v(*,i).
c
c   ipvt   = an integer array containg pivoting information.  it
c            is loaded in shefa and used in shesl.
c
c   miom   = the number of iterations performed, and current
c            order of the upper hessenberg matrix hes.
c
c   npsl   = the number of calls to psol.
c
c   iflag  = integer output flag..
c            iflag=0 means convergence in miom iterations, miom.le.mmax.
c            iflag=1 means the convergence test did not pass in mmax
c                    iterations.
c            iflag=2 means there was a recoverable error in psol
c                    caused by the preconditioner being out of date.
c            iflag=-1 means there was a nonrecoverable error in psol.
c
c note: the x array is used as work space in routines psol and atv.
c-----------------------------------------------------------------------
      integer i, info, j, k, m, mm1
      real bnrm, prod, rho, snormw, snrm2, tem
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
      save
c
      iflag = 0
      miom = 0
      npsl = 0
c zero out the hes array. ----------------------------------------------
      do 15 j = 1,mmax
        do 10 i = 1,mmax
 10       hes(i,j) = 0.0e0
 15     continue
c-----------------------------------------------------------------------
c the initial residual is the vector b.  apply scaling to b, and
c then normalize as v(*,1).
c-----------------------------------------------------------------------
      do 20 i = 1,n
 20     v(i,1) = b(i)*sf(i)
      bnrm = snrm2(n, v, 1)
      tem = 1.0e0/bnrm
      call sscal (n, tem, v(1,1), 1)
c-----------------------------------------------------------------------
c main loop to compute the vectors v(*,2) to v(*,mmax).
c the running product prod is needed for the convergence test.
c-----------------------------------------------------------------------
      prod = 1.0e0
      if (iprint .gt. 2) write(iunit,30)
 30   format(' ------ in routine spiom ------ ')
      do 90 m = 1,mmax
        miom = m
c-----------------------------------------------------------------------
c call routine atv to compute vnew = a*v(m).
c call routine svrorthog to orthogonalize the new vector vnew = v(*,m+1).
c call routine shefa to update the factors of hes.
c-----------------------------------------------------------------------
        call atv (n, u, savf, v(1,m), su, sf, x, f, jac, psol,
     *            v(1,m+1), wk, wmp, iwmp, ier, npsl)
        if (ier .ne. 0) go to 300
        call svrorthog (v(1,m+1),v,hes,n,m,mmax,iomp,snormw)
        call shefa (hes, mmax, m, ipvt, info, m)
        mm1 = m - 1
        if (m .gt. 1 .and. ipvt(mm1) .eq. mm1) prod = prod*hes(m,mm1)
        if (info .ne. m) go to 50
c-----------------------------------------------------------------------
c the last pivot in hes was found to be zero.
c if vnew = 0 or l = mmax, take an error return with iflag = 2.
c otherwise, continue the iteration without a convergence test.
c-----------------------------------------------------------------------
        if (snormw .eq. 0.0e0) go to 120
        if (m .eq. mmax) go to 120
        go to 60
c-----------------------------------------------------------------------
c update rho, the estimate of the norm of the residual b - a*x(m).
c test for convergence.  if passed, compute approximation x(m).
c if failed and m .lt. mmax, then continue iterating.
c-----------------------------------------------------------------------
 50     continue
        rho = bnrm*snormw*abs(prod/hes(m,m))
        if (iprint .gt. 2) write(iunit,55) m, rho, eps
 55     format(' m , res, eps ',i4,e25.16,1x,e25.16)
        if (rho .le. eps) go to 200
        if (m .eq. mmax) go to 150
c if l .lt. mmax, store hes(m+1,m) and normalize the vector v(*,m+1).
 60     continue
        hes(m+1,m) = snormw
        tem = 1.0e0/snormw
        call sscal (n, tem, v(1,m+1), 1)
 90     continue
c-----------------------------------------------------------------------
c l has reached mmax without passing the convergence test..
c compute a solution anyway and return with iflag = 1.
c otherwise, return with iflag = 2.
c-----------------------------------------------------------------------
 120  continue
      iflag = 2
      return
 150  iflag = 1
c-----------------------------------------------------------------------
c compute the approximation x(m) to the solution.
c since the vector x was used as work space in routine orthog,
c and since initial guess of the newton correction is zero, x must be
c set back to zero.
c-----------------------------------------------------------------------
 200  continue
      m = miom
      do 210 k = 1,m
 210    b(k) = 0.0e0
      b(1) = bnrm
      call shesl(hes,mmax,m,ipvt,b)
      do 220 k = 1,n
 220    x(k) = 0.0e0
      do 230 i = 1,m
        call saxpy (n, b(i), v(1,i), 1, x, 1)
 230    continue
      do 240 i = 1,n
 240    x(i) = x(i)/su(i)
      if (ipflg .eq. 1) then
        ier = 0
        call psol (n, u, savf, su, sf, f, jac, wk, wmp, iwmp, x, ier)
        npsl = npsl + 1
        if (ier .ne. 0) go to 300
        endif
      return
c-----------------------------------------------------------------------
c this block handles error returns forced by routine psol.
c-----------------------------------------------------------------------
 300  continue
      if (ier .lt. 0) iflag = -1
      if (ier .gt. 0) iflag = 3
      return
c----------------------- end of subroutine spiom -----------------------
      end
      subroutine atv (n, u, savf, v, su, sf, ftem, f, jac, psol, z,
     *                vtem, wmp, iwmp, ier, npsl)
      implicit none
       Use(Cdv)
      integer iwmp, ier, npsl, locwmp, locimp, iersl, kmp, mmax
      integer methn, methk, ipflg, mfdif, nfe, nje, nni, nli, npe
      integer nps, ncfl, nbcf, i
      external f, jac, psol
      integer n
      real u, savf, v, su, sf, ftem, vtem, z, wmp
      dimension  u(*), savf(n), v(n), su(n), sf(n), ftem(n), vtem(n),
     *           z(n)
      dimension iwmp(*), wmp(*)
c-----------------------------------------------------------------------
c this routine computes the product
c
c        sf*j*(p-inverse)*(su-inverse*v),
c
c where su and sf are diagonal scaling matrices, j is the jacobian
c matrix, and p is a preconditioning matrix.  the product is stored
c in z.  this is computed either by a call to jac and psol, or by a
c difference quotient, which involves calls to f and psol.
c-----------------------------------------------------------------------
c
c on entry
c
c   n      = problem size, passed to f, jac and psol.
c
c   u      = array containing current approximate to a solution of
c            f(u) = 0.
c
c   savf   = array containing current value of f(u).
c
c   v      = real array of length n.
c
c   ftem   = work array of length n.
c
c   vtem   = work array of length n used to store the unscaled
c            version of v.
c
c   su     = array of length n containing scale factors for u.
c
c   sf     = array of length n containing scale factors for f.
c
c   wmp    = real work array used by preconditioner psol.
c
c   iwmp   = integer work array used by preconditioner psol.
c
c on return
c
c   z      = array of length n containing desired scaled
c            matrix-vector product.
c
c in addition, this routine uses the common variables
c mfdif, nfe, and nje.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real delt, sqteta, rhom
      common /nks001/ delt, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c
      real sigma,tmp,utv,vtv,suitv,sdot,zero
      save
      data zero/0.0e0/
c
      if (mfdif .eq. 2) go to 200
c-----------------------------------------------------------------------
c mfdif = 1.  call user-supplied routine for computing j*v.
c-----------------------------------------------------------------------
      do 110 i = 1,n
 110    vtem(i) = v(i)/su(i)
      if (ipflg .eq. 1) then
        ier = 0
        call psol (n, u, savf, su, sf, f, jac, z, wmp, iwmp, vtem, ier)
        npsl = npsl + 1
        if (ier .ne. 0) return
        endif
      call jac (n, u, savf, vtem, z, wmp, iwmp)
      nje = nje + 1
c ... Add a 1/(pseudo time-step) carried in u(n+2): TDR 8/21/95
      do 118 i = 1,n
         z(i) = z(i) + vtem(i)*sf(i)*u(n+2)  # sf in numerator or dem??
 118  continue
c 3-15-93.  this scaling should already be performed in the jac routine.
c      do 120 i = 1,n
c 120    z(i) = z(i)*sf(i)
      return
c-----------------------------------------------------------------------
c mfdif = 2.  internally generated vector j*v.
c-----------------------------------------------------------------------
 200  continue
c set vtem = (su-inverse) * v.
      do 210 i = 1,n
        vtem(i) = v(i)/su(i)
 210    continue
      if (ipflg .eq. 0) then
c ipflg = 0.  save u in z and increment u by sigma*vtem.
          do 215 i = 1,n
 215        z(i) = u(i)*su(i)
          utv = sdot(n, z, 1, v, 1)
          suitv = zero
          do 220 i = 1,n
 220        suitv = suitv + abs(v(i))
c vtv = 1 in this case.
          sigma = sqteta*max(abs(utv),suitv)
          sigma = sign(sigma,utv)
          do 230 i = 1,n
 230        z(i) = u(i)
          do 240 i = 1,n
 240        u(i) = z(i) + sigma*vtem(i)
        else
c ipflg = 1. apply inverse of right preconditioner
          ier = 0
          call psol (n, u, savf, su, sf, f, jac, ftem, wmp, iwmp,
     *               vtem, ier)
          npsl = npsl + 1
          if (ier .ne. 0) return
          vtv = zero
          suitv = zero
          do 260 i = 1,n
            tmp = vtem(i)*su(i)
            z(i) = tmp*su(i)
            vtv = vtv + tmp*tmp
 260        suitv = suitv + abs(tmp)
          utv = sdot(n, u, 1, z, 1)
          sigma = sqteta*max(abs(utv),suitv)/vtv
          sigma = sign(sigma,utv)
          do 270 i = 1,n
 270        z(i) = u(i)
          do 280 i = 1,n
 280        u(i) = z(i) + sigma*vtem(i)
        endif
      call f(n, u, ftem)
      nfe = nfe + 1
      comnfe = comnfe + 1
ccc      do 281 i = 1, n                  # Begin 2nd order Jac
ccc         u(i) = z(i) - sigma*vtem(i)   # change sign of perturbation
ccc 281  continue
ccc      call f(n, u, vtem)		  # here vtem is second ftem
ccc      nfe = nfe + 1                    # End 2nd order Jac
      do 290 i = 1,n
 290    u(i) = z(i)
c ... u(n+2) contains the 1/(pseudo time-step) nufak: TDR 8/20/95
      do 300 i = 1,n
 300    z(i) = (ftem(i) - savf(i))/sigma - vtem(i)*u(n+2)
ccc 300    z(i) = (ftem(i) - vtem(i))/(2*sigma) - vtem(i)*u(n+2) # 2nd order Jac
cccc                                       vtem(i)*su(i)*u(n+2)/sf(i)
      do 310 i = 1,n
 310    z(i) = z(i)*sf(i)
      return
c----------------------- end of subroutine atv -------------------------
      end
      subroutine shefa (a, lda, n, ipvt, info, job)
      integer lda, n, ipvt(n), info, job
      real a(lda,1)
c-----------------------------------------------------------------------
c     this routine is a modification of the linpack routine sgefa and
c     performs an lu decomposition of an upper hessenberg matrix a.
c     there are two options available..
c
c          (1)  performing a fresh factorization
c          (2)  updating the lu factors by adding a row and a
c               column to the matrix a.
c-----------------------------------------------------------------------
c     shefa factors an upper hessenberg matrix by elimination.
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        job    integer
c                job = 1    means that a fresh factorization of the
c                           matrix a is desired.
c                job .ge. 2 means that the current factorization of a
c                           will be updated by the addition of a row
c                           and a column.
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that shesl will divide by zero if called.
c
c     modification of linpack. this version dated 07/20/83 .
c     peter brown, university of houston, lawrence livermore  lab.
c
c     subroutines and functions
c
c     blas saxpy,isamax
c
c-----------------------------------------------------------------------
c     internal variables
c
      real t
      integer isamax, j, k, km1, kp1, l, nm1
      save
c
      if (job .gt. 1) go to 80
c-----------------------------------------------------------------------
c a new facorization is desired.  this is essentially the linpack
c code with the exception that we know there is only one nonzero
c element below the main diagonal.
c-----------------------------------------------------------------------
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(2,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            a(k+1,k) = a(k+1,k)*t
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
c-----------------------------------------------------------------------
c the old factorization of a will be updated.  a row and a column
c has been added to the matrix a.
c n-1 is now the old order of the matrix.
c-----------------------------------------------------------------------
  80  continue
      nm1 = n - 1
c-----------------------------------------------------------------------
c perform row interchanges on the elements of the new column, and
c perform elimination operations on the elements using the multipliers.
c-----------------------------------------------------------------------
      if (nm1 .le. 1) go to 105
      do 100 k = 2,nm1
        km1 = k - 1
        l = ipvt(km1)
        t = a(l,n)
        if (l .eq. km1) go to 90
          a(l,n) = a(km1,n)
          a(km1,n) = t
  90    continue
        a(k,n) = a(k,n) + a(k,km1)*t
 100    continue
 105  continue
c-----------------------------------------------------------------------
c complete update of factorization by decomposing last 2x2 block.
c-----------------------------------------------------------------------
      info = 0
c
c        find l = pivot index
c
         l = isamax(2,a(nm1,nm1),1) + nm1 - 1
         ipvt(nm1) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,nm1) .eq. 0.0e0) go to 140
c
c           interchange if necessary
c
            if (l .eq. nm1) go to 110
               t = a(l,nm1)
               a(l,nm1) = a(nm1,nm1)
               a(nm1,nm1) = t
  110       continue
c
c           compute multipliers
c
            t = -1.0e0/a(nm1,nm1)
            a(n,nm1) = a(n,nm1)*t
c
c           row elimination with column indexing
c
               t = a(l,n)
               if (l .eq. nm1) go to 120
                  a(l,n) = a(nm1,n)
                  a(nm1,n) = t
  120          continue
               a(n,n) = a(n,n) + t*a(n,nm1)
         go to 150
  140    continue
            info = nm1
  150    continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
c----------------------- end of subroutine shefa -----------------------
      end
      subroutine shesl(a,lda,n,ipvt,b)
      integer lda, n, ipvt(n)
      real a(lda,1), b(n)
c-----------------------------------------------------------------------
c this is essentially the linpack routine sgesl except for changes
c due to the fact that a is an upper hessenberg matrix.
c-----------------------------------------------------------------------
c     shesl solves the real system a * x = b
c     using the factors computed by shefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from shefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from shefa.
c
c        b       real(n)
c                the right hand side vector.
c
c
c     on return
c
c        b       the solution vector  x .
c
c
c     modification of linpack. this version dated 07/20/83 .
c     peter brown, university of houston, lawrence livermore  lab.
c
c     subroutines and functions
c
c     blas saxpy
c-----------------------------------------------------------------------
c     internal variables
c
      real t
      integer k,kb,l,nm1
      save
c
      nm1 = n - 1
c
c         solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            b(k+1) = b(k+1) + t*a(k+1,k)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      return
c----------------------- end of subroutine shesl -----------------------
      end
      subroutine spigmr (n, u, savf, b, su, sf, mmax, mmaxp1, igmp,
     *                   eps, f, jac, psol, npsl, x, v, hes, q, hessv,
     *                   mgmr, wmp, iwmp, wk, methn, bnrm, ipflg, iflag,
     *                   rho)
      implicit none
      integer npsl, methn, ipflg, iprint, iunit, iermsg, j, ier, mp1
      external f, jac, psol
c-----------------------------------------------------------------------
c this routine solves the linear system a * x = b using a scaled
c preconditioned version of the generalized minimum residual method.
c-----------------------------------------------------------------------
c
c on entry
c
c   n      = problem size, passed to f, jac and psol.  also, the order
c            of the matrix a, and the lengths of the vectors u, savf,
c            b, su, sf, and x.
c
c   u      = array containing current approximate solution to
c            f(u) = 0.
c
c   savf   = array containing current value of f(u).
c
c   b      = the right hand side of the system a*x = b.
c            b is also used as work space when computing the
c            final approximation.
c
c   su     = array of length n containing scale factors for u.
c
c   sf     = array of length n containing scale factors for f.
c
c   x      = the initial guess of the solution of the system
c            a*x = b.
c
c   eps    = tolerance on residuals b-a*x in scaled euclidean norm.
c
c   mmax   = the maximum allowable order of the matrix h.
c
c   igmp   = the number of previous vectors the new vector vnew
c            must be made orthogonal to.  igmp .le. mmax.
c
c   wmp    = real work array used by preconditioner psol.
c
c   iwmp   = integer work array used by preconditioner psol.
c
c   methn  = integer flag indicating which global strategy is
c            being used (=0 inexact newton only, =1 dogleg,
c            =2 linesearch backtrack).
c            when methn=1, the x array is computed in dogstp, and
c            only the ygm value is returned in b.  otherwise, x
c            contains the gmres solution.
c
c   ipflg  = integer flag
c            ipflg=0 means no preconditioning is done.
c            ipflg.ne.0 means preconditioning is being performed.
c
c on return
c
c   x      = the final computed approximation to the solution
c            of the system a*x = b (when methn = 0 or 2).
c
c   rho    = the scaled euclidean norm of the residual associated
c            with x.
c
c   b      = the ygm value (when methn = 1).
c
c   bnrm   = the scaled euclidean norm of the right hand
c            side in a*x = b.
c
c   mgmr   = the number of iterations performed and
c            the current order of the upper hessenberg
c            matrix hes.
c
c   npsl   = the number of calls to psol.
c
c   v      = the n x (mgmr+1) array containing the mgmr
c            orthogonal vectors v(*,1) to v(*,mgmr).
c
c   hes    = the upper triangular factor of the qr decomposition
c            of the (mgmr+1) x mgmr upper hessenberg matrix whose
c            entries are the scaled inner-products of a*v(*,i)
c            and v(*,k).
c
c   hessv  = the (mgmr+1) x mgmr upper hessenberg matrix whose
c            entries are the scaled inner-products of a*v(*,i)
c            and v(*,k).
c
c   q      = a real array containing the components of the givens
c            rotations used in the qr decomposition of hes.  it
c            is loaded in sheqr and used in shels.
c
c   iflag  = integer output flag..
c            iflag=0 means convergence in mgmr iterations, mgmr.le.mmax.
c            iflag=1 means the convergence test did not pass in mmax
c                    iterations.
c            iflag=2 means there was a recoverable error in psol
c                    caused by the preconditioner being out of date.
c            iflag=-1 means there was a nonrecoverable error in psol.
c
c note: the x array is used as work space in routines psol and atv.
c-----------------------------------------------------------------------
      real v, hes, hessv, su, sf, q, b, savf, x, u, wk
      real wmp, snrm2
      integer n, iwmp, mmax, mmaxp1
      dimension v(n,mmaxp1), hes(mmaxp1,mmax), hessv(mmaxp1,mmax)
      dimension su(n), sf(n), q(*), b(n)
      dimension savf(n), x(n), u(*), wmp(*), wk(n), iwmp(*)
      real bnrm,eps,prod,rho,snormw,tem
      integer i,iflag,info,igmp,k,mgmr,m
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
      save
c
      iflag = 0
      mgmr = 0
      npsl = 0
c zero out the hes and hessv arrays ------------------------------------
      do 15 j = 1,mmax
        do 10 i = 1,mmaxp1
          hes(i,j) = 0.0e0
          hessv(i,j) = 0.0e0
 10       continue
 15     continue
c-----------------------------------------------------------------------
c the initial residual is the vector b.  apply scaling to b, and
c then normalize as v(*,1).
c-----------------------------------------------------------------------
      do 20 i = 1,n
 20     v(i,1) = b(i)*sf(i)
      bnrm = snrm2(n, v, 1)
      tem = 1.0e0/bnrm
      call sscal (n, tem, v(1,1), 1)
c-----------------------------------------------------------------------
c main loop to compute the vectors v(*,2) to v(*,mmax).
c the running product prod is needed for the convergence test.
c-----------------------------------------------------------------------
      prod = 1.0e0
      if (iprint .gt. 2) write(iunit,30)
 30   format(' ------ in routine spigmr ------')
      do 100 m = 1,mmax
        mgmr = m
c-----------------------------------------------------------------------
c call routine atv to compute vnew = a*v(m).
c call routine svrorthog to orthogonalize the new vector vnew = v(*,m+1).
c call routine sheqr to update the factors of hes.
c-----------------------------------------------------------------------
        call atv (n, u, savf, v(1,m), su, sf, x, f, jac, psol,
     *            v(1,m+1), wk, wmp, iwmp, ier, npsl)
        if (ier .ne. 0) go to 300
        call svrorthog (v(1,m+1),v,hes,n,m,mmaxp1,igmp,snormw)
        hes(m+1,m) = snormw
        hessv(m+1,m) = snormw
        do 60 i = 1,mgmr
          hessv(i,m) = hes(i,m)
 60       continue
        call sheqr(hes,mmaxp1,m,q,info,m)
        if (info .eq. m) go to 105
c-----------------------------------------------------------------------
c update rho, the estimate of the norm of the residual b - a*xl.
c-----------------------------------------------------------------------
        prod = prod*q(2*m)
        rho = abs(prod*bnrm)
c-----------------------------------------------------------------------
c test for convergence.  if passed, compute approximation xl.
c if failed and m .lt. mmax continue iterating.
c-----------------------------------------------------------------------
        if (iprint .gt. 2) write(iunit,65) m, rho, eps
 65     format(' m , res, eps ',i4,e25.16,1x,e25.16)
        if (rho .le. eps) go to 200
        if (m .eq. mmax) go to 110
c-----------------------------------------------------------------------
c rescale so that the norm of v(1,m+1) is one.
c-----------------------------------------------------------------------
        tem = 1.0e0/snormw
        call sscal (n, tem, v(1,m+1), 1)
 100    continue
 105  continue
      iflag = 2
      return
 110  iflag = 1
c-----------------------------------------------------------------------
c compute the approximation xl to the solution.
c since the vector x was used as work space in routine svrorthog and since
c the initial guess of the newton correction is zero, x must be set
c back to zero.
c-----------------------------------------------------------------------
 200  continue
      m = mgmr
      mp1 = m + 1
      do 210 k = 1,mp1
 210    b(k) = 0.0e0
      b(1) = bnrm
      call shels(hes,mmaxp1,m,q,b)
      if ((methn .eq. 0) .or. (methn .eq. 2)) then
        do 220 k = 1,n
 220      x(k) = 0.0e0
        do 230 i = 1,m
          call saxpy(n,b(i),v(1,i),1,x,1)
 230      continue
        do 240 i = 1,n
 240      x(i) = x(i)/su(i)
        if (ipflg .eq. 1) then
          ier = 0
          call psol (n, u, savf, su, sf, f, jac, wk, wmp, iwmp, x, ier)
          npsl = npsl + 1
          if (ier .ne. 0) go to 300
          endif
        endif
      return
c-----------------------------------------------------------------------
c this block handles error returns forced by routine psol.
c-----------------------------------------------------------------------
 300  continue
      if (ier .lt. 0) iflag = -1
      if (ier .gt. 0) iflag = 3
c
      return
c----------------------- end of subroutine spigmr ----------------------
      end
      subroutine dogdrv (n, wm, lenwm, iwm, leniwm, u, savf, f1nrm, x,
     *                   su, sf, stepmx, stptol, tau, iret, uprev,
     *                   fprev, unew, f1new, mxtkn, f, jac, psol,
     *                   icflag, icnstr, rlx)
c-----------------------------------------------------------------------
c this is the real version of subroutine dogdrv, which is the driver for
c the dogleg step.  its purpose is to find a unew on the dogleg curve
c such that f(unew) .le. f(u) + alpha*gt(unew-u) (alpha=1.e-4 used), 
c and scaled steplength = tau, starting with the input tau but
c increasing or decreasing tau if necessary.  also, it produces the
c starting trust region size tau for the next iteration.
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   wm     = real array of length lenwm used as work space in the
c            spigmr routine.  it contains the basis vectors v(i)
c            as well as the hessenberg matrix and its qr factors.
c
c   lenwm  = length of wm array.
c
c   iwm    = integer array of length leniwm used as work space.
c            it contains the current size m of the krylov subspace
c            (in iwm(1)), and the integer work space used by the
c            psol routine.
c
c   leniwm = length of iwm array.
c
c   u      = real array containing the previous approximate root.
c
c   f1nrm  = real scalar containing 0.5*norm(sf*f(u))**2, where
c            norm() denotes the euclidean norm.
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c   stepmx = real scalar containing the maximum allowable stepsize
c            length.
c
c   stptol = user-supplied tolerance on the minimum step length
c            s = unew-u.
c
c   tau    = the current size of the trust region.  it is the trust
c            region size tried at the beginning of the dogleg step.
c
c   uprev  = real array of length n used to hold a previous value
c            of u in the dogleg step.
c
c   fprev  = real array of length n used to hold f(uprev).
c
c   icflag = integer scalar.  if nonzero, then constraint violations
c            in the proposed new approximate solution will be checked
c            for, and the trust region size, tau, will be adjusted
c            accordingly.
c
c   icnstr = integer array of length n containing flags for checking
c            constraints.  see full description in subroutine cnstrt.
c
c   rlx    = real scalar restricting update to del(u)/u<rlx in cnstrt
c
c on return
c
c   unew   = real array containing the current approximate root.
c
c   savf   = real array containing the n-vector f(unew).
c
c   f1new  = real scalar containing 0.5*norm(sf*f(unew))**2, where
c            norm() denotes the euclidean norm.
c
c   x      = real array of length n containing the step computed
c            by the dogleg strategy.
c
c   mxtkn  = logical flag indicating if the step just taken was of
c            maximum length stepmx.
c
c   iret   = output flag.
c            iret=0 means that a satisfactory unew was found.
c            iret=1 means that the routine failed to find a unew
c                   that was sufficiently distinct from u.  this
c                   failure causes the nonlinear iteration to halt.
c
c-----------------------------------------------------------------------
      implicit none
       Use(Cdv)
      integer n, lenwm, iwm, leniwm, iret, locwmp, locimp, iersl
      integer kmp, mmax, methn, methk, ipflg, mfdif, nfe, nje, nni
      integer nli, npe, nps, ncfl, nbcf, iprint, iunit, iermsg, np1
      integer mmaxp1, m, mp1, iv, iwk, iygm, iycp, iynew, ihsv
      external f, jac, psol
      real u, uprev, fprev, unew, f1nrm, f1prv, f1new,
     *                 x, su, sf
      real beta, tau, stepmx, stptol, cpl, gml, wm, xl, savf
      dimension wm(lenwm), iwm(leniwm), x(n)
      dimension u(*), uprev(n), unew(n), savf(n), fprev(n), su(n), sf(n)
      real rlx
      integer icflag, icnstr, ivio, ivar
      dimension icnstr(n)
      logical mxtkn, dog1, nwttkn
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real delt, sqteta, rhom
      common /nks001/ delt, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
      save
c
      iret = 4
      dog1 = .true.
      if (iprint .gt. 1) write(iunit,5)
 5    format(' ------ in routine dogdrv ------ ')
c-----------------------------------------------------------------------
c initialize pointers into the wm array.
c-----------------------------------------------------------------------
      np1 = n + 1
      mmaxp1 = mmax + 1
      beta = wm(1)
      m = iwm(1)
      mp1 = m + 1
      iv = 3
      iwk = iv + n*mmax
      iygm = iwk + n + 1
      iycp = iygm + mmax
      iynew = iycp + mmax
      ihsv = iwk + n + 1 + mmax*(mmaxp1+1) + 1
      iersl = 0
      ivio = 0
 10   continue
c-----------------------------------------------------------------------
c find new step by the dogleg algorithm.
c-----------------------------------------------------------------------
      call dogstp (m, mp1, mmaxp1, wm(iygm), wm(iycp), beta, wm(ihsv),
     *             tau, wm(iynew), stepmx, dog1, nwttkn, cpl, gml, n,
     *             wm(iv), x, xl, wm(iwk), wm(locwmp), iwm(locimp),
     *             u, su, sf, savf, f, jac, psol)
      if (iersl .lt. 0) then
c nonrecoverable error from psol.  set iersl and return.
         iersl = 9
         return
         endif
c-----------------------------------------------------------------------
c check new point for constraint violations.
c ivio is set to 1 if there was a constraint violation during this call
c to dogdrv.  ivio is passed to trgupd to disable increasing the
c trust region size tau.
c-----------------------------------------------------------------------
      if (icflag .ne. 0) then
         call cnstrt (n, u, x, icnstr, tau, rlx, iret, ivar)
         if (iret .eq. 1) then
            ivio = 1
            if (iprint .gt. 1) write(iunit,9010) tau,ivio
            go to 10
            endif
         endif
c-----------------------------------------------------------------------
c check new point and update trust radius.
c-----------------------------------------------------------------------
      call trgupd (m, mp1, mmaxp1, n, np1, u, savf, f1nrm, x, xl,
     *             wm(iynew), su, sf, nwttkn, stepmx, beta, wm(ihsv),
     *             stptol,  mxtkn, tau, uprev, fprev, f1prv, unew,
     *             f1new, wm(iwk), ivio, iret, f)
c-----------------------------------------------------------------------
c if iret .eq. 2 or 3, go to 10 to repeat the step.
c if iret .eq. 0 or 1, return to calling program.
c-----------------------------------------------------------------------
      if (iprint .gt. 1) write(iunit,9000) tau,cpl,gml,iret
      if (iret .ge. 2) go to 10
      return
 9000 format(' tau= ',g12.4,' cpl= ',g12.4,' gml= ',g12.4,' iret= ',i2)
 9010 format(' tau= ',g12.4,' ivio= ',i2)
c----------------------- end of subroutine dogdrv ----------------------
      end
      subroutine dogstp (m, mp1, mmaxp1, ygm, ycp, beta, hes, tau, ynew,
     *                  stepmx, dog1, nwttkn, cpl, gml, n, v, xnew,
     *                  xnewl, wk, wmp, iwmp, u, su, sf, savf, f, jac,
     *                  psol)
      implicit none
       Use(Cdv)
      integer m, mp1, mmaxp1, n, iwmp, locwmp, locimp, iersl, kmp
      integer mmax, methn, methk, ipflg, mfdif, nfe, nje, nni, nli
      integer npe, nps, ncfl, nbcf, i, ier
      external f, jac, psol
      real ygm, ycp, beta, hes, tau, ynew, xnew, xnewl, v
      real wk, wmp, u, su, sf, savf
      dimension ygm(m), ycp(m), hes(mmaxp1,m), ynew(mp1), xnew(n),
     *          v(n,m), wk(n), wmp(*), iwmp(*), u(*), su(n), sf(n), 
     *          savf(n)
      logical dog1, nwttkn
c-----------------------------------------------------------------------
c this is subroutine dogstp, which computes the dogleg step for a
c given trust region size tau.
c
c on entry
c
c   m      = the current dimension of the krylov subspace returned
c            by routine spigmr.
c
c   mp1    = m + 1.
c
c   mmaxp1 = mmax + 1, where mmax is the maximum allowable krylov
c            subspace dimension.
c
c   ygm    = real array of length m containing the solution of the
c            m-dimensional least-squares problem computed in spigmr.
c
c   hes    = real mmaxp1-by-m array containing the hessenberg
c            matrix used in spigmr.
c
c   tau    = real scalar containing the current trust region size.
c
c   stepmx = real scalar containing the maximum allowable stepsize
c            length.
c
c   dog1   = logical variable.
c            dog1=.true. means that this is the first call to
c                        this routine by dogdrv for the current step.
c                        the length of ygm is calculated and ycp
c                        is calculated.
c            dog1=.false. means this is not the first call to this
c                         routine during the current step.
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, xnew, savf, su and sf.
c
c   v      = real n-by-m array containing the basis vectors for
c            the krylov space.
c   wmp    = real work array passed directly to psol routine.
c
c   iwmp   = integer work array passed directly to psol routine.
c
c   wk     = real array of length n used as work space.
c
c   u      = real array containing the previous approximate root.
c
c   beta   = real scalar containing norm(sf*f(u)), where norm() is
c            the euclidean norm.
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c on return
c
c   ycp    = real array of length m containing the cauchy point on
c            the dogleg curve. (only calculated when nwttkn=.false.)
c
c   cpl    = real scalar containing the length of the vector ycp.
c            (cpl = 0.e0  when nwttkn=.true.)
c
c   gml    = real scalar containing the length as the vector ygm.
c
c   ynew   = real array of length m+1 used in calculating xnew and ycp.
c
c   xnew   = real array containing the current dogleg step.
c
c   xnewl  = real scalar containing the scaled length of the
c            vector p*xnew, where p is the preconditioner matrix, and
c            is equal to norm(ynew).
c
c   nwttkn = logical variable.
c            nwttkn=.true. means that tau .ge. the length of the
c                          gmres step, and so the full gmres step
c                          is returned in xnew.
c            nwttkn=.false. means that tau .lt. the length of the
c                           gmres step, and so the cauchy point
c                           ycp must be calculated and the point
c                           on the dogleg curve found.
c
c-----------------------------------------------------------------------
      real a, b, c, cpl, gml, s, t, tem1, tem2, stepmx
      real sdot, snrm2
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real delt, sqteta, rhom
      common /nks001/ delt, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c
      save
c
      if (dog1) then
        gml = snrm2(m, ygm, 1)
        cpl = 0.0e0
        endif
      if (gml .le. tau) then
          nwttkn = .true.
          do 10 i = 1,m
            ynew(i) = ygm(i)
 10         continue
          tau = gml
c         return
        else
          if (dog1) then
c calculate cauchy point.
            dog1 = .false.
            do 20 i = 1,m
              ycp(i) = hes(1,i)
 20           continue
            tem1 = sdot(m, ycp, 1, ycp, 1)
            do 25 i = 1,mp1
              ynew(i) = 0.0e0
 25           continue
            do 30 i = 1,m
              call saxpy(mp1, ycp(i), hes(1,i), 1, ynew, 1)
 30           continue
            tem2 = sdot(mp1, ynew, 1, ynew, 1)
            t = beta*(tem1/tem2)
            call sscal(m, t, ycp, 1)
            cpl = snrm2(m, ycp, 1)
c if this is the first iteration and no initial trust region was
c provided by the user, then set tau = min(cpl,stepmx).
            if (tau .eq. -1.0e0) tau = min(cpl, stepmx)
            endif
          if (tau .le. cpl) then
              t = tau/cpl
              do 40 i = 1,m
                ynew(i) = t*ycp(i)
 40             continue
c             return
            else
              do 50 i = 1,m
                ynew(i) = ycp(i) - ygm(i)
 50             continue
              a = sdot(m, ynew, 1, ynew, 1)
              b = sdot(m, ygm, 1, ynew, 1)
              c = gml*gml - tau*tau
              s = ( -b - sqrt(max(b*b-a*c,0.e0)) )/a
              do 60 i = 1,m
                ynew(i) = ygm(i) + s*ynew(i)
 60             continue
c             return
            endif
        endif
c-----------------------------------------------------------------------
c calulate the dogleg step xnew = [v1,v2,...,vm]*ynew.
c-----------------------------------------------------------------------
      xnewl = snrm2(m, ynew, 1)
      do 90 i = 1,n
        xnew(i) = 0.e0
 90     continue
      do 100 i = 1,m
        call saxpy (n, ynew(i), v(1,i), 1, xnew, 1)
 100    continue
      do 110 i = 1,n
 110    xnew(i) = xnew(i)/su(i)
      if (ipflg .ge. 1) then
        ier = 0
        call psol (n, u, savf, su, sf, f, jac, wk, wmp, iwmp, xnew, ier)
        nps = nps + 1
        if (ier .ne. 0) iersl = -1
        endif
      return
c----------------------- end of subroutine dogstp ----------------------
      end
      subroutine trgupd (m, mp1, mmaxp1, n, np1, u, savf, f1nrm, x, xl,
     *                   ynew, su, sf, nwttkn, stepmx, beta, hes,
     *                   stptol, mxtkn, tau, uprev, fprev, f1prv, upls, 
     *                   f1pls, wk, ivio, iret, f)
      implicit none
       Use(Cdv)
      integer m, mp1, mmaxp1, n, np1, ivio, iret, locwmp, locimp
      integer iersl, kmp, mmax, methn, methk, ipflg, mfdif, nfe, nje
      integer nni, nli, npe, nps, ncfl, nbcf, i
      external f
      real u, f1nrm, x, xl, ynew, su, sf, stepmx, beta, hes,
     *                 stptol, tau, uprev, fprev, f1prv, upls, f1pls,
     *                 wk, savf
      dimension u(*), x(n), ynew(m), savf(n), su(n), sf(n),
     *          hes(mmaxp1,m), uprev(n), fprev(n), upls(n), wk(np1)
c-----------------------------------------------------------------------
c this is the real version of subroutine trgupd, which determines
c if the x(tau) returned by dogstp satisfies
c
c   f(unew) .le. f(u) + alpha*gt(unew-u) (alpha=1.e-4 used),
c
c where unew = u + x(tau).  if not, then a new tau is computed which
c hopefully will give a satisfactory x(tau)  (iret=0).  if this cannot
c be accomplished, then iret=1 is returned, and the nonlinear iteration
c is halted.
c
c on entry
c
c   m      = the current dimension of the krylov subspace returned
c            by routine spigmr.
c
c   mp1    = m + 1.
c
c   mmaxp1 = mmax + 1, where mmax is the maximum allowable krylov
c            subspace dimension.
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   np1    = n + 1.
c
c   u      = real array containing the previous approximate root.
c
c   savf   = real array containing f(u).
c
c   f1nrm  = real scalar containing 0.5*norm(sf*f(u))**2, where
c            norm() denotes the euclidean norm.
c
c   x      = real array of length n containing the current dogleg
c            step x(tau).
c
c   xl     = the scaled length of p*x(tau), where p is the
c            preconditioner matrix.
c
c   ynew   = real array of length m used in calculating x(tau).
c            x(tau) = [v1,...,vm]*ynew, where v1,...,vm are the
c            krylov basis vectors.
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c   nwttkn = logical variable.
c            nwttkn=.true. means that tau .ge. the length of the
c                          gmres step, and so the full gmres step
c                          is in x.
c            nwttkn=.false. means that tau .lt. the length of the
c                           gmres step, and so the cauchy point
c                           ycp has been calculated in dogstp.
c
c   stepmx = real scalar containing the maximum allowable stepsize
c            length.
c
c   beta   = real scalar containing norm(sf*f(u)), where norm() is
c            the euclidean norm.
c
c   hes    = real mmaxp1-by-m array containing the hessenberg
c            matrix used in spigmr.
c
c   stptol = user-supplied tolerance on the minimum step length
c            s = unew-u.
c
c   tau    = the current size of the trust region.  it is the trust
c            region size tried at the beginning of the dogleg step.
c
c   uprev  = real array of length n used to hold a previous value
c            of u in the dogleg step.
c
c   fprev  = real array of length n used to hold f(uprev).
c
c   f1prv  = real scalar used to hold 0.5*norm(sf*f(uprev))**2, where
c            norm() is the euclidean norm.
c
c   wk     = real array of length n used as work space.
c
c   ivio   = integer scalar indicating if there was a constraint
c            violation during the current call to dogdrv.
c            ivio=1 means there was a constraint violation,
c                =0 means there was not.
c            when ivio=1, the trust region size tau is not
c            allowed to increase.
c
c on return
c
c   upls   = real array containing the current approximate root.
c
c   savf   = real array containing the n-vector f(upls).
c
c   f1pls  = real scalar containing 0.5*norm(sf*f(upls))**2, where
c            norm() denotes the euclidean norm.
c
c   mxtkn  = logical flag indicating if the step just taken was of
c            maximum length stepmx.
c
c   iret   = output flag.
c            iret=0 means that a satisfactory upls was found.
c            iret=1 means that the routine failed to find a upls
c                   that was sufficiently distinct from u.
c            iret=2 means that the current tau was too big.  tau
c                   is reduced and the global step is repeated.
c            iret=3 means that the current tau may be too small.
c                   tau is increased and the global step is retried.
c
c   tau    = the new trust region size to be used at the beginning
c            of the next call to dogdrv, or for repeating the
c            current global step.
c
c-----------------------------------------------------------------------
      real alpha, pt1, pt5, pt75, pt99, zero, two
      real sdot, delf, dfpred, tautmp, fnrmp,
     *                 rlngth, slpi
      real vnormnk
      logical mxtkn, nwttkn
      save
      data alpha/1.e-4/, pt1/0.1e0/, pt5/0.5e0/, pt75/0.75e0/,
     *     pt99/0.99e0/, zero/0.0e0/, two/2.0e0/
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real delt, sqteta, rhom
      common /nks001/ delt, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c
      mxtkn = .false.
c-----------------------------------------------------------------------
c calculate upls = u + x.
c-----------------------------------------------------------------------
      call scopy (n, u, 1, upls, 1)
      do 20 i = 1,n
        u(i) = upls(i) + x(i)
 20     continue
      call f (n, u, savf)
      nfe = nfe + 1
      comnfe = comnfe + 1
      call sswap(n, u, 1, upls, 1)
      fnrmp = vnormnk(n, savf, sf)
      f1pls = pt5*fnrmp*fnrmp
      delf = f1pls - f1nrm
      slpi = -beta*sdot(m, hes, mmaxp1, ynew, 1)
      if (iret .ne. 3) f1prv = zero
      if ( (iret.eq.3) .and.
     *     ((f1pls.ge.f1prv) .or. (delf.gt.alpha*slpi)) ) then
c reset upls to uprev and terminate global step.
           iret = 0
           call scopy (n, uprev, 1, upls, 1)
           call scopy (n, fprev, 1, savf, 1)
           f1pls = f1prv
           return
        elseif (delf .ge. alpha*slpi) then
c f1(upls) too large.
           call slngth(n, u, x, su, rlngth)
           if (rlngth .lt. stptol) then
c (upls - u) too small, terminate global step.
               iret = 1
               call scopy (n, u, 1, upls, 1)
               return
             else
c reduce tau, continue global step.
               iret = 2
               tautmp = -slpi*xl/(two*(delf-slpi))
               if (tautmp .lt. pt1*tau) then
                    tau = pt1*tau
                 elseif (tautmp .gt. pt5*tau) then
                    tau = pt5*tau
                 else
                    tau = tautmp
                 endif
             return
             endif
        else
c f1(upls) sufficiently small.
           do 45 i = 1,mp1
             wk(i) = zero
 45          continue
           do 50 i = 1,m
             call saxpy (mp1, ynew(i), hes(1,i), 1, wk, 1)
 50          continue
           dfpred = slpi + pt5*sdot (mp1, wk, 1, wk, 1)
           if ( (iret.ne.2) .and. 
     *          ( (abs(dfpred-delf).le.pt1*abs(delf)) .or.
     *                (delf.le.slpi) ) .and. (.not.nwttkn) .and.
     *                (tau.le.pt99*stepmx) .and. (ivio.eq.0) ) then
c double tau and continue global step.
               iret = 3
               call scopy (n, upls, 1, uprev, 1)
               call scopy (n, savf, 1, fprev, 1)
               f1prv = f1pls
               tau = min( two*tau, stepmx )
               return
             else
c accept upls as new iterate.  choose new tau.
               iret = 0
               if (xl .gt. pt99*stepmx) mxtkn = .true.
               if (delf .ge. pt1*dfpred) then
c decrease tau for next iteration.
                   tau = tau/two
                 elseif (delf .le. pt75*dfpred) then
c increase tau for next iteration.
                   tau = min( two*tau, stepmx )
                 else
c leave tau unchanged.
                 endif
             endif
        endif
c     return
c----------------------- end of subroutine trgupd ----------------------
      end
      subroutine lnsrch(n, u, savf, f1nrm, p, su, sf, stepmx, stptol,
     *                  iret, unew, f1nrmp, mxtkn, f, jac, icflag,
     *                  icnstr, rlx,  adjf1)
c-----------------------------------------------------------------------
c this routine is the main driver for the linesearch algorithm.
c its purpose is to find a unew = u + rl*p in the direction p from u
c such that
c                               t
c    f(unew) .le. f(u) + alpha*g (unew-u) (alpha=1.e-4 used),
c
c and
c                              t
c    f(unew) .ge. f(u) + beta*g (unew-u)  (beta=0.9e0 used),
c
c where 0 .lt. rl .le. 1.
c-----------------------------------------------------------------------
c
c ON ENTRY ******************************
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   u      = real array containing the previous approximate root.
c
c   f1nrm  = real scalar containing 0.5*norm(sf*f(u))**2, where
c            norm() denotes the euclidean norm.
c
c   p      = step direction returned by routine model.
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c   stepmx = real scalar containing the maximum allowable stepsize
c            length.
c
c   stptol = user-supplied tolerance on the minimum step length
c            s = unew-u.
c
c   icflag = integer scalar.  if nonzero, then constraint violations
c            in the proposed new approximate solution will be checked
c            for, and the trust region size, tau, will be adjusted
c            accordingly.
c
c   icnstr = integer array of length n containing flags for checking
c            constraints.  see full description in subroutine cnstrt.
c
c   rlx    = real scalar restricting update to del(u)/u<rlx in cnstrt
c
c   adjf1  = real scalar for global strategy test: frnm_new/adjf1>=frnm_old
c
c ON RETURN  ******************************
c
c   unew   = real array containing the current approximate root.
c
c   savf   = real array containing the n-vector f(unew).
c
c   f1nrmp = real scalar containing 0.5*norm(sf*f(unew))**2, where
c            norm() denotes the euclidean norm.
c
c   iret   = output flag.
c            iret=0 means that a satisfactory unew was found.
c            iret=1 means that the routine failed to find a unew
c                   that was sufficiently distinct from u.
c
c   mxtkn  = logical flag indicating if the step just taken was of
c            maximum length stepmx.
c
c   nbcf   = number of steps for which beta-condition was not met
c            (in common).  this value is incremented by one if
c            the beta-condition could not be met on this call.
c
c-----------------------------------------------------------------------
      implicit none
       Use(Cdv)
      integer n, iret, icflag, icnstr, locwmp, locimp, iersl, kmp
      integer mmax, methn, methk, ipflg, mfdif, nfe, nje, nni, nli
      integer npe, nps, ncfl, nbcf, iprint, iunit, iermsg, i, ivio
      external f, jac
      real savf,u,unew,f1nrm,f1nrmp,p,su,sf
      real stptol,stepmx,pnrm,ratio,ratio1,slpi,rlngth,rl
      real fnrmp,rltmp,vnormnk
      real beta,f1lo,f1nprv,rldiff,rlincr,rlmin,rllo,
     1                 rlmax,rlprev
      real rlx, adjf1
      integer ivar
      dimension savf(n),u(*),unew(n),p(n),su(n),sf(n),icnstr(n)
      real pt1,pt1trl,pt99,one,two,alpha,acond,mcond,bcond
      logical mxtkn
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real delt, sqteta, rhom, tau
      common /nks001/ delt, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
c
      save
      data pt1/0.1e0/,alpha/1.e-4/,one/1.0e0/,two/2.0e0/
      data pt99/0.99e0/,beta/0.9e0/
c
      mxtkn = .false.
      pnrm = vnormnk(n,p,su)
      ratio = one
      if (pnrm .gt. stepmx) then
        ratio = stepmx/pnrm
        do 10 i = 1,n
          p(i) = p(i)*ratio
 10       continue
        pnrm = stepmx
        endif
      tau = pnrm
      ivio = 0
      if (iprint .gt. 1) write(iunit,9) pnrm
 9    format(' ------ in routine lnsrch (pnrm=',e12.4,') ---------- ')
      if (icflag .ne. 0) then
 11      continue
         call cnstrt (n, u, p, icnstr, tau, rlx, iret, ivar)
         if (iret .eq. 1) then
            ivio = 1
            ratio1 = tau/pnrm
            ratio = ratio*ratio1
            do 15 i = 1,n
              p(i) = p(i)*ratio1
 15           continue
            pnrm = tau
            if (iprint .gt. 1) write(iunit,16) pnrm, ivar
 16         format(
     .    ' ------ in routine lnsrch (ivio=1,pnrm=',e12.4,
     .    ', var index = ',i6,') -- ')
            go to 11
            endif
         endif
      slpi = -two*f1nrm*ratio
CCC should an adjf1 be used for f1nrm above??
      if (methk .eq. 2) slpi = slpi + rhom*rhom*ratio
ccc MVU: 15-jan-2020
      if (methk .eq. 3) slpi = slpi + rhom*rhom*ratio
ccc MVU: 15-jan-2020
      call slngth(n, u, p, su, rlngth)
      rlmin = stptol/rlngth
      rl = one
       if (iprint .gt. 1) write(iunit,20) rlmin,stptol,rlngth
 20   format(' -- in routine lnsrch (min lambda=',3e12.4,') -- ')
c-----------------------------------------------------------------------
c begin iteration to find rl value satisfying alpha- and beta-
c conditions.  if rl becomes .lt. rlmin, then terminate with iret=1.
c-----------------------------------------------------------------------
 115  continue
      call scopy(n, u, 1, unew, 1)
      do 120 i = 1,n
        u(i) = unew(i) + rl*p(i)
 120    continue
      call f(n, u, savf)
      nfe = nfe + 1
      comnfe = comnfe + 1
      call sswap(n, u, 1, unew, 1)
      fnrmp = vnormnk(n,savf,sf)
      f1nrmp = fnrmp*fnrmp/two
      acond=f1nrmp/adjf1 - f1nrm + alpha*slpi*rl
 
      if (iprint .gt. 1) then
       write(iunit,125) rl,f1nrm,f1nrmp,acond,nfe
      endif
 125  format(' lambda,f1,f1new,acon,nfe',4d28.16,I4)
 126  format(' lambda,f1,f1new,bcon,nfe',4d28.16,I4)
 127  format(' lambda,f1,f1new,acon,bcon,mcon,nfe',3d28.16,3d12.2,I4)
      if (f1nrmp/adjf1 .gt. f1nrm + alpha*slpi*rl) go to 200
c alpha-condition satisfied.  now check for beta-condition.
        if (f1nrmp/adjf1 .lt. f1nrm + beta*slpi*rl) then
          if ( (rl.eq.one) .and. (pnrm.lt.stepmx) ) then
            rlmax = stepmx/pnrm
            if (ivio .eq. 1) rlmax = one
 130        continue
            rlprev = rl
            f1nprv = f1nrmp
            rl = min(two*rl,rlmax)
            call scopy(n, u, 1, unew, 1)
            do 140 i = 1,n
 140          u(i) = unew(i) + rl*p(i)
            call f(n, u, savf)
            nfe = nfe + 1
            comnfe = comnfe + 1
            call sswap(n, u, 1, unew, 1)
            fnrmp = vnormnk(n,savf,sf)
            f1nrmp = fnrmp*fnrmp/two
            if (iprint .gt. 1) then
	     bcond=f1nrmp/adjf1 - f1nrm + beta*slpi*rl
              write(iunit,126) rl,f1nrm,f1nrmp,bcond,nfe
            endif
	    if ( (f1nrmp/adjf1 .le. f1nrm + alpha*slpi*rl) .and.
     *           (f1nrmp/adjf1 .lt. f1nrm + beta*slpi*rl)  .and.
     *           (rl .lt. rlmax) ) go to 130
            endif
          if ( (rl.lt.one) .or.
     *       ((rl.gt.one).and.(f1nrmp/adjf1.gt.f1nrm+alpha*slpi*rl)) ) 
     *                                                    then
            rllo = min(rl,rlprev)
            rldiff = abs(rlprev-rl)
            if (rl .lt. rlprev) then
                f1lo = f1nrmp
              else
                f1lo = f1nprv
              endif
 150        continue
            rlincr = rldiff/two
            rl = rllo + rlincr
            call scopy(n, u, 1, unew, 1)
            do 160 i = 1,n
 160          u(i) = unew(i) + rl*p(i)
            call f(n, u, savf)
            nfe = nfe + 1
            comnfe = comnfe + 1
            call sswap(n, u, 1, unew, 1)
            fnrmp = vnormnk(n,savf,sf)
            f1nrmp = fnrmp*fnrmp/two

            if (f1nrmp/adjf1 .gt.f1nrm + alpha*slpi*rl) then
                rldiff = rlincr
              elseif (f1nrmp/adjf1 .lt. f1nrm + beta*slpi*rl) then
                rllo = rl
                rldiff = rldiff - rlincr
                f1lo = f1nrmp
              endif
	    
	    if (iprint .gt. 1) then
	     mcond=rldiff-rlmin
	     acond=f1nrmp/adjf1 - f1nrm + alpha*slpi*rl
             bcond=f1nrmp/adjf1 - f1nrm + beta*slpi*rl
            write(iunit,127) rl,f1nrm,f1nrmp,acond,bcond,mcond,nfe
            endif
            if ( ( (f1nrmp/adjf1 .gt. f1nrm + alpha*slpi*rl) .or.
     *             (f1nrmp/adjf1 .lt. f1nrm + beta*slpi*rl) )  .and.
     *           (rldiff .gt. rlmin) ) go to 150
            if (f1nrmp/adjf1 .lt. f1nrm + beta*slpi*rl) then
c beta condition could not be satisfied.  set unew to last
c u value that satisfied alpha-condition and continue.
c increment counter on number of steps beta-condition not satisfied.
              call scopy(n, u, 1, unew, 1)
              do 170 i = 1,n
 170            u(i) = unew(i) + rllo*p(i)
              call f(n, u, savf)
              nfe = nfe + 1
              comnfe = comnfe + 1
              call sswap(n, u, 1, unew, 1)
              fnrmp = vnormnk(n,savf,sf)
              f1nrmp = fnrmp*fnrmp/two
              nbcf = nbcf + 1
              endif
            endif
          endif
        iret = 0
        if (rl*pnrm .gt. pt99*stepmx) mxtkn = .true.
        return
c-----------------------------------------------------------------------
c alpha-condition not satisfied.  perform quadratic backtrack to
c compute new rl value.
c-----------------------------------------------------------------------
 200  continue
      if (rl .lt. rlmin) then
c no satisfactory unew can be found sufficiently distinct from u.
c copy u into unew, set return code and return.
        iret = 1
        call scopy(n, u, 1, unew, 1)
c load savf and f1nrmp with current values.
        call f(n, u, savf)
        nfe = nfe + 1
        comnfe = comnfe + 1
        fnrmp = vnormnk(n,savf,sf)
        f1nrmp = fnrmp*fnrmp/two
        return
        endif
c compute new rl using a quadratic backtrack.
      rltmp = -slpi/(two*(f1nrmp - f1nrm - slpi))
      if (rltmp .gt. rl/two) rltmp = rl/two
      rlprev = rl
      f1nprv = f1nrmp
      pt1trl = pt1*rl
      if (rltmp .le. pt1trl) then
          rl = pt1trl
        else
          rl = rltmp
        endif
      go to 115
c----------------------- end of subroutine lnsrch ----------------------
      end
      subroutine errgen (ierr, v1, v2, i1, i2)
c-----------------------------------------------------------------------
c this routine prints error messages from the driver nksol.  the output
c from this routine can be turned off by setting a flag in iwork.
c-----------------------------------------------------------------------
      implicit none
      integer ierr, i1, i2, iprint, iunit, iermsg
      real v1, v2
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
      save
c
      if (iermsg .gt. 0) return
      if (ierr .eq. 10) then
        write(iunit,9000)
        return
        endif
      if (ierr .eq. 20) then
        write(iunit,9010)
        return
        endif
      if (ierr .eq. 30) then
        write(iunit,9020)
        return
        endif
      if (ierr .eq. 40) then
        write(iunit,9030) i1,i2
        return
        endif
      if (ierr .eq. 50) then
        write(iunit,9040)  i1,v1
        return
        endif
      if (ierr .eq. 100) then
        write(iunit,9100) i1,i2
        return
        endif
      if (ierr .eq. 110) then
        write(iunit,9110) i1,i2
        return
        endif
      if (ierr .eq. 120) then
        write(iunit,9120) i1,i2
        return
        endif
      if (ierr .eq. 130) then
        write(iunit,9130) i2
        return
        endif
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 9000   format(//
     * ' nksol ---  illegal value for mf.  mf must be between '
     */'            1 and 4, or between -3 and -1.'
     *)
 9010   format(//
     * ' nksol ---  illegal value for mdif.  mdif must either 0 or 1. '
     *)
 9020   format(//
     * ' nksol ---  illegal value for ipflag.  ipflag must be either '
     */'            0 or 1.'
     *)
 9030   format(//
     * ' nksol ---  illegal value for optional input in iwork. '
     */'            iwork(',i2,') = ',i3,' must be nonnegative.'
     *)
 9040   format(//
     * ' nksol ---  illegal value for optional input in rwork. '
     */'            rwork(',i2,') = ',e12.4,' must be nonnegative.'
     *)
 9100   format(//
     * ' nksol ---  insufficient length for rwork.'
     */'            rwork length given  - ',i8,
     */'            rwork length needed - ',i8
     *)
 9110   format(//
     * ' nksol ---  insufficient length for iwork.'
     */'            iwork length given  - ',i8,
     */'            iwork length needed - ',i8
     *)
 9120   format(//
     * ' nksol ---  maximum number of beta-condition test failures',
     */'            exceeded.'
     */'            number of beta-condition failures  - ',i8,
     */'            maximum number of failures allowed - ',i8
     *)
 9130   format(//
     * ' nksol ---  initial u failed to satisfy all constraints',
     */'            u(',i8,') violated its constraint.'
     *)
      return
c-------------  end of subroutine errgen  ------------------------------
      end
      subroutine cnstrt (n, u, x, icnstr, tau, rlx, iret, ivar)
c-----------------------------------------------------------------------
c this is the real version of subroutine cnstrt, which checks for
c constraint violations in the proposed new approximate solution u+x.
c if a constraint violation occurs, then a new trust region size, tau,
c is calculated, and this value is to be given to routine dogstp to
c calculate a new approximate solution u+x.
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, x and icnstr.
c
c   u      = real array containing the previous approximate root.
c
c   x      = real array containing the current dogleg step.
c
c   icnstr = integer array of length n containing flags indicating
c            which entries in unew are to be constrained.
c            if icnstr(i) .gt. 0, then u(i)+x(i) must be .gt. 0,
c            if icnstr(i) .lt. 0, then u(i)+x(i) must be .lt. 0, while
c            if icnstr(i) .eq. 0, then u(i)+x(i) is not constrained.
c
c   rlx    = real scalar restricting update to abs(x/u) < fac2*rlx 
c
c   tau    = the current size of the trust region.  it is the trust
c            region size tried at the beginning of the dogleg step.
c            equivalently, the current steplength for lnsrch.
c
c on return
c
c   tau    = the adjusted size of the trust region if a constraint
c            violation occurred (otherwise, it is unchanged).  it is
c            the trust region size to give to routine dogstp to
c            calculate a new dogleg step or linesearch steplength.
c
c   iret   = output flag.
c            iret=0 means that unew satisfied all constraints.
c            iret=1 means that unew failed to satisfy all the
c                   constraints, and a new dogleg or linesearch step
c                   must be computed.
c
c   ivar   = index of variable causing the constraint to be violated
c
c-----------------------------------------------------------------------
      implicit none
      integer n
      real u, x
      real tau, fac, zero
      real rlx, mxchng, fac2, cutlo
      integer icnstr, i, iret, ivar
      dimension u(*), x(n), icnstr(n)
      data fac /0.6e0/, zero/0.e0/, fac2 /0.9e0/, cutlo /1.0e-200/
      save
c-----------------------------------------------------------------------
c check constraints for proposed new step u+x.  if a constraint has
c been violated, then calculate a new trust region size, tau, to be
c used in dogstp, or equivalently a new step length to be used in
c lnsrch.
c-----------------------------------------------------------------------
      iret = 0
      mxchng = 0.
      ivar = 0
      do 100 i = 1,n
         if (icnstr(i) .gt. 0) then
            if ( abs(x(i)/(u(i)+cutlo)) .gt. mxchng) then
               mxchng = abs(x(i)/(u(i)+cutlo))
               ivar = i
            endif
            if ( (u(i) + x(i)) .le. zero) then
               tau = fac*tau
               ivar = i
               iret = 1
               return
            endif
         elseif (icnstr(i) .lt. 0) then
            if ( abs(x(i)/(u(i)+cutlo)) .gt. mxchng) then
	       mxchng = abs(x(i)/(u(i)+cutlo))
               ivar = i
            endif
            if ( (u(i) + x(i)) .ge. zero) then
               tau = fac*tau
               ivar = i
               iret = 1
               return
            endif
         endif
 100  continue

      if(mxchng .ge. rlx) then
         tau = fac2*tau*rlx/mxchng
         iret = 1
      endif
c
      return
c----------------------- end of subroutine cnstrt ----------------------
      end
      subroutine cnstrt0 (n, u, icnstr, rlx, iret, ivar)
c-----------------------------------------------------------------------
c this is the real version of subroutine cnstrt, which checks for
c constraint violations in the initial approximate solution u.
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u and icnstr.
c
c   u      = real array containing the initial approximate root.
c
c   icnstr = integer array of length n containing flags indicating
c            which entries in unew are to be constrained.
c            if icnstr(i) .gt. 0, then u(i) must be .gt. 0,
c            if icnstr(i) .lt. 0, then u(i) must be .lt. 0, while
c            if icnstr(i) .eq. 0, then u(i) is not constrained.
c
c   rlx    = real scalar restricting update to del(u)/u<rlx in cnstrt0
c
c on return
c
c   iret   = output flag.
c            iret=0    means that u satisfied all constraints.
c            iret.ne.0 means that u(iret) failed to satisfy its
c                      constraint.
c
c-----------------------------------------------------------------------
      implicit none
      integer n
      real u, zero
      real rlx
      integer icnstr, i, iret, ivar
      dimension u(*), icnstr(n)
      save
      data zero/0.e0/
c-----------------------------------------------------------------------
c check constraints for initial u.  if a constraint has been violated,
c set iret = 1 to signal an error return to call routine.
c-----------------------------------------------------------------------
      iret = 0
      do 100 i = 1,n
         if (icnstr(i) .gt. 0) then
            if ( u(i) .le. zero) then
               iret = 1
               ivar = i
               return
            endif
         elseif (icnstr(i) .lt. 0) then
            if ( u(i) .ge. zero) then
               iret = 1
               ivar = i
               return
            endif 
        endif
 100  continue
      return
c----------------------- end of subroutine cnstrt0 ---------------------
      end
      subroutine infgen (iterm, v1, v2, i1, i2)
c-----------------------------------------------------------------------
c this routine prints informational messages from the driver nksol.  
c the output from this routine can be turned off by setting a flag in
c iwork.
c-----------------------------------------------------------------------
      implicit none
      integer iterm, i1, i2, iprint, iunit, iermsg
      real v1, v2
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
      save
c
      if (iermsg .gt. 0) return
      if (iterm .eq. 1) then
        if (iprint .ge. 1) write(iunit,9000)
        return
        endif
      if (iterm .eq. 2) then
        if (iprint .ge. 1) write(iunit,9010)
        return
        endif
      if (iterm .eq. 3) then
        if (iprint .ge. 1) write(iunit,9020)
        return
        endif
      if (iterm .eq. 4) then
        if (iprint .ge. 1) write(iunit,9030) 
        return
        endif
      if (iterm .eq. 5) then
        if (iprint .ge. 1) write(iunit,9040) 
        return
        endif
      if (iterm .eq. 6) then
        if (iprint .ge. 1) write(iunit,9050) 
        return
        endif
      if (iterm .eq. 7) then
        if (iprint .ge. 1) write(iunit,9060) 
        return
        endif
      if (iterm .eq. 8) then
        if (iprint .ge. 1) write(iunit,9070) 
        return
        endif
      if (iterm .eq. 9) then
        if (iprint .ge. 1) write(iunit,9080) 
        return
        endif
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 9000   format(//
     * ' nksol ---  iterm = 1.'
     */'            maxnorm(sf*f(u)) .le. ftol, where maxnorm() is'
     */'            the maximum norm function.  u is probably an'
     */'            approximate root of f.'
     *)
 9010   format(//
     * ' nksol ---  iterm = 2.'
     */'            the scaled distance between the last two'
     */'            steps is less than stptol.  u may be an'
     */'            approximate root of f, but it is also possible'
     */'            that the algorithm is making very slow progress'
     */'            and is not near a root, or that stptol is too'
     */'            large.'
     *)
 9020   format(//
     * ' nksol ---  iterm = 3.'
     */'            the last global step failed to reduce norm(f)'
     */'            sufficiently.  either u is close to a root of f'
     */'            and no more accuracy is possible, or the'
     */'            finite-difference approximation to j*v is'
     */'            inaccurate, or stptol is too large.  if the ncfl'
     */'            optional output value is close to the nni value,'
     */'            it may be the case that the krylov iteration is'
     */'            converging very slowly.  in this case, the user'
     */'            may want to use preconditioning and/or increase'
     */'            the mmax value (i.e., increase the maximum'
     */'            dimension of the krylov subspace.)'
     *)
 9030   format(//
     * ' nksol ---  iterm = 4.'
     */'            the maximum allowable number of nonlinear'
     */'            iterations has been reached.'
     *)
 9040   format(//
     * ' nksol ---  iterm = 5.'
     */'            5 consecutive steps of length stepmx (the'
     */'            maximum stepsize limit) have been taken.  either'
     */'            norm(f) asymptotes from above to a finite value'
     */'            in some direction, or stepmx is too small.'
     *)
 9050   format(//
     * ' nksol ---  iterm = 6.'
     */'            more than 10 failures occurred when trying to'
     */'            satisfy the beta-condition in the linesearch'
     */'            algorithm.  it is likely that the iteration is'
     */'            making poor progress.'
     *)
 9060   format(//
     * ' nksol ---  iterm = 7.'
     */'            there was a breakdown in the krylov'
     */'            iteration.  this will likely only occur when'
     */'            the jacobian matrix j or j*(p-inverse) is ill-'
     */'            conditioned.  if this error return occurs with'
     */'            mf=2, try either mf=1 or mf=3 instead.'
     *)
 9070   format(//
     * ' nksol ---  iterm = 8.'
     */'            there was a nonrecoverable error in pset'
     */'            causing the iteration to halt.'
     *)
 9080   format(//
     * ' nksol ---  iterm = 9.'
     */'            there was a nonrecoverable error in psol'
     */'            causing the iteration to halt.'
     *)
      return
c-------------  end of subroutine infgen  ------------------------------
      end
      subroutine inexct(n, u, savf, f1nrm, p, su, sf, stepmx, stptol,
     *                  iret, unew, f1nrmp, mxtkn, f, icflag,
     *                  icnstr, rlx)
c-----------------------------------------------------------------------
c this routine is the main driver for the inexact newton algorithm.
c its purpose is to compute a unew = u + p in the direction p from u,
c taking the full inexact newton step p.  The step may be contrained if
c the constraint conditions are violated, or if the norm of p is
c greater than stepmx.
c-----------------------------------------------------------------------
c
c on entry
c
c   n      = size of the nonlinear system, and the length of arrays
c            u, unew, savf, su and sf.
c
c   u      = real array containing the previous approximate root.
c
c   f1nrm  = real scalar containing 0.5*norm(sf*f(u))**2, where
c            norm() denotes the euclidean norm.
c
c   p      = step direction returned by routine model.
c
c   su     = real array containing scaling factors for u.
c
c   sf     = real array containing scaling factors for f(u).
c
c   stepmx = real scalar containing the maximum allowable stepsize
c            length.
c
c   stptol = user-supplied tolerance on the minimum step length
c            s = unew-u.
c
c   icflag = integer scalar.  if nonzero, then constraint violations
c            in the proposed new approximate solution will be checked
c            for, and the trust region size, tau, will be adjusted
c            accordingly.
c
c   icnstr = integer array of length n containing flags for checking
c            constraints.  see full description in subroutine cnstrt.
c
c   rlx    = real scalar restricting update to del(u)/u<rlx in cnstrt
c
c on return
c
c   unew   = real array containing the current approximate root.
c
c   savf   = real array containing the n-vector f(unew).
c
c   f1nrmp = real scalar containing 0.5*norm(sf*f(unew))**2, where
c            norm() denotes the euclidean norm.
c
c   iret   = output flag.
c            iret=0 means that a satisfactory unew was found.
c            iret=1 means that the routine failed to find a unew
c                   that was sufficiently distinct from u.
c
c   mxtkn  = logical flag indicating if the step just taken was of
c            maximum length stepmx.
c
c-----------------------------------------------------------------------
      implicit none
       Use(Cdv)
      integer n, iret, icflag, icnstr, locwmp, locimp, iersl, kmp
      integer mmax, methn, methk, ipflg, mfdif, nfe, nje, nni, nli
      integer npe, nps, ncfl, nbcf, iprint, iunit, iermsg, i, ivio, ivar
      real beta
      external f
      real savf,u,unew,f1nrm,f1nrmp,p,su,sf
      real stepmx,stptol,pnrm,ratio,ratio1
      real fnrmp,vnormnk
      real rlx
      dimension savf(n),u(*),unew(n),p(n),su(n),sf(n),icnstr(n)
      real pt99,one,two
      logical mxtkn
c-----------------------------------------------------------------------
c     nks001 common block.
c-----------------------------------------------------------------------
      real delt, sqteta, rhom, tau
      common /nks001/ delt, rhom, sqteta, locwmp, locimp, iersl, kmp,
     *                mmax, methn, methk, ipflg, mfdif, nfe, nje, nni,
     *                nli,  npe, nps, ncfl, nbcf
c-----------------------------------------------------------------------
c     nks002 common block.
c-----------------------------------------------------------------------
      common /nks002/ iprint, iunit, iermsg
c
      save
      data one/1.0e0/,two/2.0e0/
      data pt99/0.99e0/,beta/0.9e0/
c
      mxtkn = .false.
      pnrm = vnormnk(n,p,su)
      ratio = one
      if (pnrm .gt. stepmx) then
        ratio = stepmx/pnrm
        do 10 i = 1,n
          p(i) = p(i)*ratio
 10       continue
        pnrm = stepmx
        endif
      tau = pnrm
      ivio = 0
      if (iprint .gt. 1) write(iunit,9) pnrm
 9    format(' ------ in routine inexct (pnrm=',e12.4,') ---------- ')
      if (icflag .ne. 0) then
 11      continue
         iret = 0
         call cnstrt (n, u, p, icnstr, tau, rlx, iret, ivar)
         if (iret .eq. 1) then
            ivio = 1
            ratio1 = tau/pnrm
            ratio = ratio*ratio1
            do 15 i = 1,n
              p(i) = p(i)*ratio1
 15           continue
            pnrm = tau
            if (iprint .gt. 1) write(iunit,16) pnrm, ivar
 16         format(
     .    ' ------ in routine inexct (ivio=1,pnrm=',e12.4,
     .    ', var index = ',i6,') -- ')
            if (pnrm .le. stptol) then
               iret = 1
               return
               endif
            go to 11
            endif
         endif
      call scopy(n, u, 1, unew, 1)
      do 20 i = 1,n
        u(i) = unew(i) + p(i)
 20     continue
      call f(n, u, savf)
      nfe = nfe + 1
      comnfe = comnfe + 1
      call sswap(n, u, 1, unew, 1)
      fnrmp = vnormnk(n,savf,sf)
      f1nrmp = fnrmp*fnrmp/two
      if (pnrm .gt. pt99*stepmx) mxtkn = .true.
      return
c----------------------- end of subroutine inexct ----------------------
      end
