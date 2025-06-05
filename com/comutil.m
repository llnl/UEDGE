c!include "../mppl.h"
c!include "../sptodp.h"
      subroutine com_set_dims
      implicit none
c ... Compute dimension variables nx, ny, nxm, and nym

c ... Common blocks:
      Use(Dim)     # nx,ny,nxm,nym,nxpt
      Use(Share)   # nycore,nysol,nxleg,nxcore,nxomit,geometry

c ... Compute variables used for grid and xpoint index dimensions.

      nxpt = 1
      nxm = nxleg(1) + nxcore(1) + nxcore(2)
     .                     + nxleg(2)  + 4*nxxpt
      nym = nycore + nysol
      if (geometry=="dnull" .or. geometry(1:9)=="snowflake" .or.
     .    geometry=="isoleg") then
	 nxpt = 2   # for geometry=isoleg, set after flxrun/grdrun
         if(geometry=="dnull") then  #Symmetric double-null config
           nxm = 2*nxm - 2
           nym = nycore + nysol + nyout
         endif
      endif
      nx = nxm - abs(nxomit)
      ny = nym - nyomitmx

      return
      end

c----------------------------------------------------------------------c

      subroutine wspltim

c...  Writes timing data for splines

      implicit none

      Use(Timespl)

      write(*,900) 'Total in B2VAhL spline routine = ', totb2val
      write(*,900) 'Total in INTRhV spline routine = ', totintrv
 900  format(a36,f10.4,20x,' sec')
      totb2val = 0.0e0
      totintrv = 0.0e0

      return
      end

c----------------------------------------------------------------------c

*DECK BSPDOC
      SUBROUTINE BSPDOC
C***BEGIN PROLOGUE  BSPDOC
C***PURPOSE  Documentation for BSPLINE, a package of subprograms for
C            working with piecewise polynomial functions
C            in B-representation.
C***LIBRARY   SLATEC
C***CATEGORY  E, E1A, K, Z
C***TYPE      ALL (BSPDOC-A)
C***KEYWORDS  B-SPLINE, DOCUMENTATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C           Boisvert, R. F., (NIST)
C***DESCRIPTION
C
C                          ****  Abstract  ****
C
C  BSPDOC is a non-executable documentation routine. Its purpose is  to
C  describe a suite of subprograms for performing basic operations with
C  piecewise polynomial (PP) functions. The operations included are (a)
C  evaluation  of PP functions and their derivatives, (b) interpolation
C  using PP functions, (c) numerical integration of PP  functions,  and
C  (d)  utility  operations such as evaluating B-spline basis functions
C  and converting between various  representations  for  PP  functions.
C  Also included are routines for interpolating gridded data in two and
C  three  dimensions  using  tensor  products  of  one-dimensional   PP
C  functions.
C
C                         ****  Background  ****
C
C  The basic package described here is a modification of  that  of  (de
C  Boor, 1977).  It differs from that package in the following ways:
C
C    * Names have been altered to prevent duplication  and  conflicts
C      with routines from (de Boor, 1978).
C    * The call lists used here are different.
C    * Work vectors were  added  to  ensure  portability  and  proper
C      execution  in  an overlay environment These work arrays can be
C      used for other purposes except as noted in BSPVhN.
C    * While most of the original routines in (de  Boor,  1977)  were
C      restricted  to orders 20 or less, this restriction was removed
C      from all routines except the quadrature  routine  BSQAD.  (See
C      the  section  below  on  differentiation  and  integration for
C      details.)
C
C  In 1990 the package was revised and  expanded.   Among  the  changes
C  introduced were:
C
C    * The existing B-spline evaluation routine, BVALU, was  replaced
C      by B1VAhL.
C    * Routines for  parametric  interpolation  of  arbitrary  planar
C      curves were added.
C    * Routines for interpolation of gridded data in  two  and  three
C      dimensions were added.
C
C  All the routines in the package are available  in  both  single  and
C  double  precision  versions. The double precision routines names are
C  obtained from the single precision names by prefixing with a D.  For
C  example,  B1VAhL  and  DB1VAhL  are  the  single  and double precision
C  routines for evaluating a B-spline or any of its derivatives in  the
C  B-representation.
C
C               ****  Piecewise Polynomial Functions  ****
C
C  Let x(1), x(2), ..., x(m) be a  strictly  increasing  sequence  with
C  x(1)=a and x(m)=b. This defines a subdivision of the interval (a,b).
C  A set of m-1 polynomials of fixed degree  k-1,  each  defined  on  a
C  unique  subinterval  (x(j),x(j+1)),  j=1,2,  ...,  m-1,  is called a
C  piecewise polynomial (PP) function of order k. If  the  PP  function
C  has  k-2 continuous derivatives on (a,b), then it is called a spline
C  of order k.
C
C  A PP function is uniquely determined by
C
C    (a) the number of subintervals m
C    (b) the order (degree+1) of the polynomial pieces k
C    (c) an increasing breakpoint sequence x(1),...,x(m)
C    (d) m-1 sets of k polynomial coefficients  (c(i,j),  i=1,...,k),
C        one for each subinterval, j=1,...,m-1
C
C  On the jth subinterval (x(j),x(j+1)) the function may be written
C
C                   k                  i-1
C         p(j,x) = SUM ( c(i,j)(x-x(j))    / (i-1)! )
C                  i=1
C
C  This is the so-called PP-representation for PP functions.  Note that
C  the  "polynomial  coefficients"  given here are actually the (right)
C  derivatives at the breakpoints.
C
C  Each of the m-1 polynomial pieces has k coefficients, making a total
C  of  k(m-1) parameters. In the most general case, the PP function and
C  its derivatives have m-2 jumps, one at each of the breakpoints x(j),
C  j=2,...,m-1.   Continuity   requirements   at  the  breakpoints  add
C  constraints and reduce the  number  of  free  parameters.  If  a  PP
C  function  is  continuous  at  each  of the m-2 interior breakpoints,
C  there are k(m-1)-(m-2) free parameters; if in addition  the  PP  has
C  continuous   first   derivatives,   there   are  k(m-1)-2(m-2)  free
C  parameters, etc., until we get to a spline where we have  k(m-1)-(k-
C  1)(m-2)  =  m+k-2  free  parameters.  Thus,  the  principle  is that
C  increasing the continuity of derivatives  decreases  the  number  of
C  free parameters and conversely.
C
C  The PP-representation is easy to understand and  fast  to  evaluate.
C  However,  when  trying to find a PP function which satisfies special
C  properties  (e.g.,   continuity,   interpolation)   an   alternative
C  representation is preferred. This is the so-called B-representation.
C  Here we write the PP function in the form
C
C                                 n
C                       p(x)  =  SUM a(i) b(i,x)
C                                i=1
C
C  That is, the PP function  is  expanded  in  a  fixed  basis  b(i,x),
C  i=1,...,n, and the a(i) are the free parameters. The basis functions
C  are splines which are zero except on (at most) k adjoining intervals
C  where  each  b(i,x)  is  positive  and,  in most cases, hat or bell-
C  shaped.
C
C  Since each basis function is a spline,  any  linear  combination  of
C  them  must  also  be  a spline, and hence the sum for p(x) above can
C  only produce splines.  A convenient way to  represent  PP  functions
C  with  less  continuity  is  to allow more than one breakpoint at the
C  same location.  Since we  defined  the  breakpoint  sequence  to  be
C  strictly  increasing  we  use  the  term knot sequence to denote the
C  nondecreasing sequence of breakpoints with each multiple  breakpoint
C  counted  multiple  times.  If two knots are allowed to come together
C  at some x(j), then we say that we have  a  knot  of  multiplicity  2
C  there,  and  the knot values are the x(j) values.  If we reverse the
C  procedure of the first paragraph, we find  that  adding  a  knot  to
C  increase  multiplicity  increases the number of free parameters and,
C  according to the principle above, thereby introduces a discontinuity
C  in  what  was  the highest continuous derivative at that knot. Thus,
C  the number of free parameters is n = l+k-2, where l is  the  sum  of
C  multiplicities at the x(i) values with x(1) and x(m) of multiplicity
C  1.  If l = m, all knots are simple, and we have a spline.  Each knot
C  can have a multiplicity of at most k.
C
C  In order for the nonzero part of b(1,x)  to  be  a  spline  covering
C  (x(1),x(2)), it is necessary to put k-1 additional knots to the left
C  of a.  Similarly, it is necessary to place k-1 additional  knots  to
C  the  right of b.  We see that the number of knots is l+2(k-1) = n+k.
C  We denote the knot  sequence  by  t(j),  j=1,  ...,  n+k.   In  many
C  problems where extrapolation beyond a or b is not anticipated, it is
C  common practice to set t(1)=t(2)=...=t(k)=a  and  t(n+1)=t(n+2)=...=
C  t(n+k)=b.
C
C  In summary, the B-representation of a PP function is given by,
C
C    (a) the order (degree+1) of the polynomial pieces k
C    (b) an nondecreasing knot sequence t(1),...,t(n+k)
C    (c) the coefficients a(i), i = 1,...,n
C
C  A PP function written using the B-representation is known  as  a  B-
C  spline.
C
C  The nonzero part of the basis function b(i,x) lies in  the  interval
C  (t(i),t(i+k)).  The fact that each basis function is nonzero over at
C  most k intervals means that for a given x value, there are at most k
C  nonzero  terms  of  the  sum  defining  p(x).  This  leads to banded
C  matrices in linear algebra problems, and references (de Boor,  1978)
C  and  (Hanson,  1979)  take  advantage of this in constructing higher
C  level routines to achieve speed and avoid ill-conditioning.
C
C                       ****  Basic Routines  ****
C
C  (Note: in the following we use upper case names  to  denote  Fortran
C  variables used in calling the routines in the B-spline package.)
C
C  The basic routines which most  casual  users  will  need  are  those
C  concerned   with   direct   evaluation   of  PP  functions.  The  B-
C  representation is preferred because of  numerical  stability.   Here
C  one  must  supply   the number of coefficients N, the order K of the
C  polynomial pieces (of degree K-1), the  knots  (T(j),  j=1,...,N+K),
C  and  the  B-spline  coefficients (A(i), i=1,...,N).  To evaluate the
C  B-spline or any of its derivatives at a point XVAL, one can use
C
C               Y = B1VAhL(XVAL,ID,T,N,K,A,INBV,WORK,IFLAG)
C
C  where ID  is  a  non-negative  integer  indicating  that  the  ID-th
C  derivative  is  desired (ID=0 gives function value itself).  INBV is
C  an initialization parameter which is set to 1 by  users  before  the
C  first call; distinct splines require distinct INBV parameters.  WORK
C  is a scratch array of length at least 3K,  and  IFLAG  is  a  return
C  code.
C
C  When more conventional  communication  is  needed  for  publication,
C  physical interpretation, etc., the B-representation can be converted
C  to to the PP-representation by
C
C                  CALL BSPPP(T,A,N,K,LDC,C,X,LX,WORK)
C
C  This  call  returns  the  number  of  polynomial  pieces   LX,   the
C  breakpoints  (distinct knots) (X(j),  j=1,...,LX+1), and the (right)
C  derivatives C(*,j) at each breakpoint X(j).  Here  LDC.GE.K  is  the
C  leading  dimension  of  the  matrix C and WORK is a scratch array of
C  length at least (N+3)K. To evaluate  the  B-spline  or  any  of  its
C  derivatives in this representation, one uses
C
C                 Y = PPVAL(LDC,C,X,LX,K,ID,XVAL,INPPV)
C
C  INPPV is an initialization parameter which is  set  to  1  by  users
C  before the first call.
C
C  A major reason for considering PPVAL  is  that  evaluation  is  much
C  faster than that from B1VAhL.  However, one must view this conversion
C  from the B- to the PP-representation with  some  skepticism  because
C  the conversion may result in the loss of significant digits when the
C  B-spline varies in an almost discontinuous fashion.  To determine to
C  what  extent  the  conversion  process  loses  digits,  compute  the
C  relative difference ABS((Y1-Y2)/Y2) for a set of x's on  (a,b)  with
C  Y1 from PPVAL and Y2 from B1VAhL.
C
C  No  routines  are  supplied  to  explicitly  convert  from  the  PP-
C  representation to the B-representation.  However, one can obtain the
C  B-representation by evaluating the PP function (using PPVAL)  at  an
C  appropriate  number  of  abscissas  to  create  data  for use by the
C  interpolation routine  BINThK  (see  below),  which  returns  the  B-
C  representation.   In order to correctly generate a knot sequence for
C  the B-representation one must know the smoothness of the PP function
C  at  each  of  the  breakpoints.  If the PP function is of continuity
C  class C(j) in the neighborhood of a  given  breakpoint,  then  k-1-j
C  knots  are required there.  If the function itself is discontinuous,
C  k-2 knots are required.  The total number of knots  must  be  n+k  =
C  l+2(k-1),  where  l  is  the  sum  of the multiplicities required to
C  obtain the required level of smoothness at the breakpoints.  If  the
C  PP  function is a spline, then l=m; otherwise l must be greater (see
C  the description of l above).
C
C  Recall that jump discontinuities in a PP function or its derivatives
C  can  occur  at  breakpoints.   B1VAhL and PPVAL return right-limiting
C  values in these cases.  An  exception  occurs  for  XVAL=X(LX+1)  in
C  PPVAL  or  for  XVAL=T(N+K) for B1VAhL where left-limiting values are
C  returned. Note also that  a  computed  XVAL  which,  mathematically,
C  would  be  a  breakpoint  value  may differ from the breakpoint by a
C  round-off error.  When this happens in evaluating a discontinuous PP
C  function   or  some  discontinuous  derivative,  the  value  at  the
C  breakpoint and the value at XVAL can be radically different. In this
C  case, setting XVAL to a T or X value makes the computation precise.
C
C                       ****  Interpolation  ****
C
C  Given a set of data points ((Y(i),F(i)), i=1,...,N),  the  order  K,
C  and  the knot sequence (T(j), j=1,...,N+K), BINThK generates B-spline
C  coefficients (A(i),i=1,...,N) which define a PP function which  will
C  interpolate the data by calls to B1VAhL.  A similar interpolation can
C  also be done for cubic splines using BINTh4.
C
C              ****  Differentiation and Integration  ****
C
C  Derivatives of PP  functions  are  obtained  from  B1VAhL  or  PPVAL.
C  Integrals  are  obtained  from  BSQAD using the B-representation and
C  PPQAD  using  the  PP-representation.  More  complicated   integrals
C  involving  the  product  of a function F and some derivative of a PP
C  function can be evaluated with BFQAD or PFQAD using the  B-  or  PP-
C  representations  respectively.  All  quadrature routines, except for
C  PPQAD, are limited in accuracy to 18 digits  or  working  precision,
C  whichever is smaller. PPQAD is limited to working precision only. In
C  addition, the order K for BSQAD is limited to 20 or less. If  orders
C  greater than 20 are required, use BFQAD with F(X) = 1.
C
C                       ****  Extrapolation  ****
C
C  Extrapolation outside the interval (a,b) can be accomplished  easily
C  by  the  PP-representation  using PPVAL.  However, caution should be
C  exercised, especially when multiple knots are located at a or  b  or
C  when  the  extrapolation  is carried significantly beyond a or b. On
C  the other hand, B1VAhL returns zero for XVAL  outside  (T(1),T(N+K)),
C  and   hence  some  manipulation  of  the  knots  may  be  needed  to
C  extrapolate with B1VAhL.
C
C                ****  Curve Fitting and Smoothing  ****
C
C  Unless one has many accurate data points,  direct  interpolation  is
C  not  recommended  for summarizing data. The results are often not in
C  accordance with intuition since the fitted curve tends to  oscillate
C  through  the  set  of  points.  Monotone splines (see SLATEC routine
C  PCHDOC (Fritsch and Carlson, 1980)) can help  curb  this  undulating
C  tendency  but  constrained  least  squares is more likely to give an
C  acceptable fit  with  fewer  parameters.   The  SLATEC  routine  FC,
C  described  in  (Hanson,  1978), is recommended for this purpose. The
C  output from this fitting process is the B-representation.
C
C                   **** Parametric Interpolation ****
C
C  Given a set of data points ((X(i),Y(i)), i=1,...,N),  the  order  K,
C  and  the  parameter  breakpoint  sequence  (P(i),  i=1,...,N), BINKP
C  generates a knot sequence (T(j),  j=1,...,N+K),  and  parametric  B-
C  spline    coefficients   ((AX(i),AY(i)),   i=1,...,N)   which   will
C  interpolate the data by calls to B1VAhL.  That  is,  the  x-component
C  function  has B-representation T,N,K,AX and the y-component function
C  has B-representation T,N,K,AY determined so that
C
C              x(P(i)) = X(i),  y(P(i)) = Y(i), i=1,...,N.
C
C  Parametric interpolation is suitable when  the  data  are  arbitrary
C  points in the plane, and do not necessarily represent a functional
C  dependence of Y on X.
C
C  A similar interpolation  can  also  be  done  for  parametric  cubic
C  splines using BIN4P.  An auxiliary routine PARGN is available to aid
C  in the determination of a  suitable  parameter  breakpoint  sequence
C  (P(i),i=1,...,N).
C
C                  ****  Tensor Product B-splines  ****
C
C  One can easily build up multidimensional B-splines  (on  rectangular
C  regions)  using  tensor  products  of one-dimensional B-splines.  We
C  illustrate this process in two dimensions.  Here such a function may
C  be written as
C
C                            nx   ny
C                p(x,y)  =  SUM  SUM  a(i,j) b(i,x) b(j,y)
C                           i=1  j=1
C
C  where the functions b(i,x) and b(j,y) are  one-dimensional  B-spline
C  basis  functions.  Distinct orders and knot sequences may be used in
C  each direction.  Note that for each fixed value of y p(x,y) is a  PP
C  function of x alone, and for each fixed x p(x,y) is a PP function of
C  y alone.  The function B2VAhL is designed to evaluate PP functions of
C  this  type.   B3VAL  performs  the  corresponding operation in three
C  dimensions.
C
C               ****  Multidimensional Interpolation  ****
C
C  The routines B2INhT and B3INT produce tensor product B-splines  which
C  interpolate  gridded data in two and three dimensions, respectively.
C  By gridded data  in  two  dimensions  we  mean  ((X(i),Y(j),F(i,j)),
C  i=1,...,NX, j=1,...,NY).
C
C                          ****  Summary  ****
C
C  Below we list the routines in this package according  to  what  jobs
C  they  perform.  The notation (D)SUB indicates that the routine comes
C  in both single and double precision. The single  precision  name  is
C  SUB and the double precision name is DSUB.
C
C  One dimension
C
C    Interpolation
C
C      (D)BINTh4 - interpolates with B-splines of order 4 (cubic)
C      (D)BINThK - interpolates with B-splines of order k
C
C    Parametric interpolation
C
C      (D)BIN4P - interpolates with parametric B-splines of order 4
C                 (cubic)
C      (D)BINKP - interpolates with parametric B-splines of order k
C      (D)PARGN - computes a parametrization from the data points.
C
C    Evaluation
C
C      (D)B1VAhL - evaluates the B-representation or a derivative
C      (D)PPVAL - evaluates the PP-representation or a derivative
C
C    Quadrature
C
C      (D)BFQAD - integrates the product of a function F and any spline
C                 derivative in the B-representation
C      (D)BSQAD - integrates the B-representation on subintervals
C      (D)PFQAD - integrates the product of a function F and any spline
C                 derivative in the PP-representation
C      (D)PPQAD - integrates the PP-representation
C
C    Utility
C
C      (D)INTRhV - gets the largest index of the knot to the left of x
C      (D)BSPDR - sets up difference array for BSPEV
C      (D)BSPEV - evaluates the B-representation and derivatives at x
C      (D)BSPPP - converts from B- to PP-representation
C      (D)BSPVhD - computes all nonzero B-spline basis functions and
C                 their derivatives at x
C      (D)BSPVhN - called by BSPEV, BSPVhD, BSPPP and BINThK for basis
C                 function and derivative evaluations
C
C  More than one dimension
C
C    Interpolation
C
C      (D)B2INhT - interpolates gridded data in 2D with tensor products
C                 of one-dimensional B-splines
C      (D)B3INT - interpolates gridded data in 3D with tensor products
C                 of one-dimensional B-splines
C
C    Evaluation
C
C      (D)B2VAhL - evaluates a 2D tensor product B-spline or derivative
C      (D)B3VAL - evaluates a 3D tensor product B-spline or derivative
C
C  Subsidiary Routines
C
C    (D)BKCHK, (D)BKNOT, (D)BNSLhV, (D)BNFAhC, (D)BSGQ8,
C    (D)BTPCF, (D)GNDIS, (D)GNMET, (D)PARTF, (D)PPGQ8
C
C  Error Handler
C
C    XERMShG
C
C  Machine-Dependent Routines
C
C    R1MACH9 - Replaced by EPSILON (JRC 27 Apr 08)
C
C***REFERENCES  D. E. Amos, Computation with splines and B-splines,
C                 Report SAND78-1968, Sandia Laboratories, March 1979.
C               D. E. Amos, Quadrature subroutines for splines and
C                 B-splines, Report SAND79-1825, Sandia Laboratories,
C                 December 1979.
C               Carl de Boor, On calculating with B-splines, Journal of
C                  Approximation Theory 6, (1972), pp. 50-62.
C               Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C               F. N. Fritsch and R. E. Carlson, Monotone piecewise
C                 cubic interpolation, SIAM Journal on Numerical Ana-
C                 lysis 17, 2 (April 1980), pp. 238-246.
C               R. J. Hanson, Constrained least squares curve fitting
C                 to discrete data using B-splines, a user's guide,
C                 Report SAND78-1291, Sandia Laboratories, December
C                 1978.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900723  PURPOSE section revised.  (WRB)
C   901120  DESCRIPTION rewritten.  (RFB)
C   930525  Changed name PARAM to the less common PARGN.  (FNF)
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  BSPDOC
C***FIRST EXECUTABLE STATEMENT  BSPDOC
      RETURN
      END

*DECK B3INT
      SUBROUTINE B3INT (X, NX, Y, NY, Z, NZ, KX, KY, KZ, TX, TY, TZ,
     +   FCN, LDF1, LDF2, WORK, IFLAG)
C***BEGIN PROLOGUE  B3INT
C***PURPOSE  B3INT determines a piecewise polynomial function that
C            interpolates three-dimensional gridded data. Users  specify
C            the polynomial order (degree+1) of the interpolant and
C            (optionally) the knot sequence.
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (B3INT-S, DB3INT-D)
C***KEYWORDS  INTERPOLATION, THREE DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   B3INT determines the parameters of a function that interpolates the
C   three-dimensional  gridded data (X(i),Y(j),Z(k),FCN(i,j,k)) for
C   i=1,..,NX, j=1,..,NY, and k=1,..,NZ. The interpolating function and
C   its derivatives may subsequently be evaluated by function B3VAL.
C
C   The interpolating function is a piecewise polynomial (pp) function
C   represented as a tensor product of one-dimensional B-splines.  The
C   form of this function is
C
C                        NX   NY   NZ
C          S(x,y,z)  =  SUM  SUM  SUM  a    U (x) V (y) W (z)
C                       i=1  j=1  k=1   ijk  i     j     k
C
C   where the functions U(i), V(j), and W(k) are one-dimensional B-
C   spline basis functions. The coefficients a(i,j,k) are chosen so
C
C   S(X(i),Y(j),Z(k))=FCN(i,j,k) for i=1,..,NX, j=1,..,NY, k=1,..,NZ.
C
C   Note that for fixed y and z S is a pp function of x alone, for
C   fixed x and z S is a pp function of y alone, and for fixed x and y
C   S is a pp function of z alone.  In one dimension a pp function may
C   be created by partitioning a given interval into subintervals and
C   defining a distinct polynomial on each one.  The points where
C   adjacent subintervals meet are called knots.  Each of the functions
C   U(i), V(j), and W(k) above is a piecewise polynomial.
C
C   Users of B3INT choose the order (degree+1) of the polynomial pieces
C   used to define the interpolant in each of the x, y and z directions
C   (KX, KY, and KZ). Users also may define their own knot sequences in
C   x, y and z separately (TX, TY, and TZ).  If IFLAG=0, however, B3INT
C   will choose knots that result in a pp interpolant with KX-2, KY-2
C   and KZ-2 continuous partial derivatives in x, y and z respectively.
C   The interpolating function is identically zero outside the rectan-
C   gular region defined by the knots.  See below for more information
C   on knot selection.
C
C   After a call to B3INT, all information necessary to define the
C   interpolating function is contained in the parameters NX, NY, NZ,
C   KX, KY, KZ, TX, TY, TZ, and FCN. These quantities should not be
C   altered until after the last call of the evaluation routine B3VAL.
C
C
C   I N P U T
C   ---------
C
C   X       Real 1D array (size NX)
C           Array of x abscissae. Must be strictly increasing.
C
C   NX      Integer scalar (2 .LE. NX .LE. LDF1)
C           Number of x abscissae.
C
C   Y       Real 1D array (size NY)
C           Array of y abscissae. Must be strictly increasing.
C
C   NY      Integer scalar (2 .LE. NY .LE. LDF2)
C           Number of y abscissae.
C
C   Z       Real 1D array (size NZ)
C           Array of z abscissae. Must be strictly increasing.
C
C   NZ      Integer scalar (NZ .GE. 2)
C           Number of z abscissae.
C
C   KX      Integer scalar (2 .LE. KX .LE. NX)
C           The order (degree + 1) of polynomial pieces in x.
C
C   KY      Integer scalar (2 .LE. KY .LE. NY)
C           The order (degree + 1) of polynomial pieces in y.
C
C   KZ      Integer scalar (2 .LE. KZ .LE. NZ)
C           The order (degree + 1) of polynomial pieces in z.
C
C
C   I N P U T   O R   O U T P U T
C   -----------------------------
C
C   TX      Real 1D array (size NX+KX)
C           The knots in the x direction. If IFLAG=0 these are chosen
C           by B3INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TX(i).LT.X(i).LT.TX(i+KX),  i=1,..,NX.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NX. (See information below on
C           knot placement.)
C
C   TY      Real 1D array (size NY+KY)
C           The knots in the y direction. If IFLAG=0 these are chosen
C           by B3INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TY(i).LT.Y(i).LT.TY(i+KY),  i=1,..,NY.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NY. (See information below on
C           knot placement.)
C
C   TZ      Real 1D array (size NZ+KZ)
C           The knots in the y direction. If IFLAG=0 these are chosen
C           by B3INT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TZ(i).LT.Z(i).LT.TZ(i+KZ),  i=1,..,NZ.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NZ. (See information below on
C           knot placement.)
C
C
C   I N P U T   A N D   O U T P U T
C   -------------------------------
C
C   FCN     Real 3D array (size LDF1 by LDF2 by NZ)
C           Input : Array of function values to interpolate. FCN(I,J,K)
C                   should contain the function value at the point
C                   (X(I),Y(J),Z(K)).
C           Output: Array of coefficients of the B-spline interpolant.
C
C   LDF1    Integer scalar (LDF1 .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C   LDF2    Integer scalar (LDF2 .GE. NY)
C           The actual second dimension of FCN used in the calling
C           program.
C
C
C   O T H E R
C   ---------
C
C   WORK    Real 1D array
C           Array of working storage. Must be dimensioned of length at
C           least NX*NY*NZ + 2*max( KX*(NX+1), KY*(NY+1), KZ*(NZ+1))
C
C   IFLAG   Integer scalar.
C           Must be set by the user before B3INT is called.
C           On return IFLAG indicates the status of the output.
C
C           Input :  0 : knot sequence chosen by user
C                    1 : knot sequence chosen by B3INT
C
C           Output:  0 : successful execution
C                    2 : IFLAG out of range
C                    3 : NX or LDF1 out of range
C                    4 : KX out of range
C                    5 : X not strictly increasing
C                    6 : TX is an illegal knot sequence
C                    7 : NY or LDF2 out of range
C                    8 : KY out of range
C                    9 : Y not strictly increasing
C                   10 : TY is an illegal knot sequence
C                   11 : NZ out of range
C                   12 : KZ out of range
C                   13 : Z not strictly increasing
C                   14 : TZ is an illegal knot sequence
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C           After a successful return, B3INT may be recalled without
C           resetting IFLAG, provided that only the array FCN has
C           changed.
C
C
C   K N O T   S E L E C T I O N
C   ---------------------------
C
C   In this section we describe the relationship between data points
C   and knots.  Users who choose to let B3INT select the knot sequences
C   TX, TY and TZ (by setting IFLAG=1) may skip this discussion.
C
C   To describe the relationship between knots and data points we first
C   consider the simpler case of one-dimensional interpolation; in
C   particular, we consider interpolating the data X(i), i=1,..,N with
C   a pp function of order K.
C
C   Knots are the points where the individual polynomial pieces join up,
C   and hence where the pp function may suffer a loss of smoothness.  To
C   define a pp function, one needs N+K knots. If the knots are distinct
C   the interpolant will be as smooth as possible (continuous, with K-2
C   continuous derivatives). If two adjacent knots come together, the
C   smoothness of the function is reduced at that point. In general, if
C   M knots are assigned to the same point, the pp function will have
C   K-M-1 continuous derivatives there. If K knots are taken at the same
C   point, then the pp function itself will be discontinuous there.
C
C   Typically, K knots are taken at or to the left of the leftmost data
C   point, K knots at or to the right of the rightmost data point, with
C   the remaining N-K knots in between.  In order for there to be a
C   solution to the interpolation problem the knots T(i), must satisfy
C   certain additional constraints. That is,
C
C                T(i) .LT. X(i) .LT. T(i+K)   i=1,..,N
C
C   Equality is permitted on the left for i=1 and on the right for i=N.
C
C   The three-dimensional interpolant computed by this routine is a
C   tensor product of one-dimensional pp interpolants.  The knots form
C   a grid (TX(i),TY(j),TZ(k)) in the same way that the data points do.
C   Along lines parallel to the coordinate axes, the interpolant reduces
C   to a one-dimensional pp function with order and knots (KX,TX),
C   (KY,TY), or (KZ,TZ). In this case the appropriate constraints on the
C   knots become
C
C            TX(i) .LT. X(i) .LT. TX(i+KX)   i=1,..,NX
C            TY(i) .LT. Y(i) .LT. TY(i+KY)   i=1,..,NY
C            TZ(i) .LT. Z(i) .LT. TZ(i+KZ)   i=1,..,NZ
C
C   with equality on the left permitted when i=1 and equality on the
C   right permitted when i=NX, i=NY and i=NZ respectively.
C
C   If these conditions are violated, then B3INT returns with IFLAG
C   equal to 6, 10 or 14.  The default knot sequence selected by B3INT
C   always satisfies these conditions.
C
C   When the user sets IFLAG=1 B3INT selects knots as follows. KX knots
C   are taken at each endpoint in the x direction, not-a-knot end
C   conditions (see references) are used, and the remaining knots are
C   placed at data points if KX is even or at midpoints between data
C   points if KX is odd.  The y and z directions are treated similarly.
C   This yields a three-dimensional pp function with KX-2, KY-2 and
C   KZ-2 continuous partial derivatives in x, y and z respectively. The
C   interpolant is zero outside the rectangular region defined by the
C   data points, and discontinuous along the boundary of this region.
C
C***SEE ALSO  B3VAL
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  BKCHK, BKNOT, BTPCF, BUPCK, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  B3INT
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          NX, NY, NZ, KX, KY, KZ, LDF1, LDF2, IFLAG
      real             X(NX), Y(NY), Z(NZ), TX(NX+KX), TY(NY+KY),
     +                 TZ(NZ+KZ), FCN(LDF1,LDF2,NZ), WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'SLATEC', SUBROU = 'B3INT')
      CHARACTER*50     MESSAG
      INTEGER          I, IW, J, K, LOC
      LOGICAL          BKCHK
C
      EXTERNAL         BKCHK, BKNOT, BTPCF, BUPCK, XERMShG
C
C***FIRST EXECUTABLE STATEMENT  B3INT
C
C  -----------------------
C  CHECK VALIDITY OF INPUT
C  -----------------------
C
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 9020
C
      IF ((NX .LT. 2) .OR. (NX .GT. LDF1))  GO TO 9030
      IF ((KX .LT. 2) .OR. (KX .GT. NX))  GO TO 9040
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 9050
   10 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. BKCHK(X,NX,KX,TX)) GO TO 9060
      ENDIF
C
      IF ((NY .LT. 2) .OR. (NY .GT. LDF2))  GO TO 9070
      IF ((KY .LT. 2) .OR. (KY .GT. NY))  GO TO 9080
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 9090
   20 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. BKCHK(Y,NY,KY,TY)) GO TO 9100
      ENDIF
C
      IF (NZ .LT. 2)  GO TO 9110
      IF ((KZ .LT. 2) .OR. (KZ .GT. NZ))  GO TO 9120
      DO 30 I=2,NZ
         IF (Z(I) .LE. Z(I-1))  GO TO 9130
   30 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. BKCHK(Z,NZ,KZ,TZ)) GO TO 9140
      ENDIF
C
C  ------------
C  CHOOSE KNOTS
C  ------------
C
      IF (IFLAG .EQ. 1) THEN
         CALL BKNOT(X,NX,KX,TX)
         CALL BKNOT(Y,NY,KY,TY)
         CALL BKNOT(Z,NZ,KZ,TZ)
      ENDIF
C
C  -------------------------------
C  CONSTRUCT B-SPLINE COEFFICIENTS
C  -------------------------------
C
      IW = NX*NY*NZ + 1
C
C     ... COPY FCN TO WORK IN PACKED FORM FOR BTPCF
C
      LOC = 0
      DO 300 K=1,NZ
         DO 200 J=1,NY
            DO 100 I=1,NX
               WORK(I+LOC) = FCN(I,J,K)
  100       CONTINUE
            LOC = LOC + NX
  200    CONTINUE
  300 CONTINUE
C
C     ... TENSOR-PRODUCT INTERPOLATION
C
      CALL BTPCF(X,NX,WORK,NX,NY*NZ,TX,KX,FCN, NY*NZ,WORK(IW))
      CALL BTPCF(Y,NY,FCN, NY,NX*NZ,TY,KY,WORK,NX*NZ,WORK(IW))
      CALL BTPCF(Z,NZ,WORK,NZ,NX*NY,TZ,KZ,FCN, NX*NY,WORK(IW))
C
C     ... UNPACK FCN
C
      CALL BUPCK(FCN,NX,NY,NZ,LDF1,LDF2)
C
      IFLAG = 0
      GO TO 9999
C
C  -----------
C  ERROR EXITS
C  -----------
C
 9020 CONTINUE
      IFLAG = 2
      MESSAG = 'IFLAG IS OUT OF RANGE'
      GO TO 9900
C
 9030 CONTINUE
      IFLAG = 3
      MESSAG = 'NX OR LDF1 IS OUT OF RANGE'
      GO TO 9900
C
 9040 CONTINUE
      IFLAG = 4
      MESSAG = 'KX IS OUT OF RANGE'
      GO TO 9900
C
 9050 CONTINUE
      IFLAG = 5
      MESSAG = 'X ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9060 CONTINUE
      IFLAG = 6
      MESSAG = 'TX IS AN ILLEGAL KNOT SEQUENCE'
      GO TO 9900
C
 9070 CONTINUE
      IFLAG = 7
      MESSAG = 'NY OR LDF2 IS OUT OF RANGE'
      GO TO 9900
C
 9080 CONTINUE
      IFLAG = 8
      MESSAG = 'KY IS OUT OF RANGE'
      GO TO 9900
C
 9090 CONTINUE
      IFLAG = 9
      MESSAG = 'Y ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9100 CONTINUE
      IFLAG = 10
      MESSAG = 'TY IS AN ILLEGAL KNOT SEQUENCE'
      GO TO 9900
C
 9110 CONTINUE
      IFLAG = 11
      MESSAG = 'NZ IS OUT OF RANGE'
      GO TO 9900
C
 9120 CONTINUE
      IFLAG = 12
      MESSAG = 'KZ IS OUT OF RANGE'
      GO TO 9900
C
 9130 CONTINUE
      IFLAG = 13
      MESSAG = 'Z ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9140 CONTINUE
      IFLAG = 14
      MESSAG = 'TZ IS AN ILLEGAL KNOT SEQUENCE'
      GO TO 9900
C
 9900 CONTINUE
      CALL XERMShG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
C
C  ----
C  EXIT
C  ----
C
 9999 CONTINUE
      RETURN
      END

*DECK BUPCK
      SUBROUTINE BUPCK (A, NX, NY, NZ, LD1, LD2)
C***BEGIN PROLOGUE  BUPCK
C***SUBSIDIARY
C***PURPOSE  Converts an array packed as though dimensioned A(NX,NY,NZ)
C            to one dimensioned as A(LD1,LD2,NZ).
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BUPCK-S, DBUPCK-D)
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   BUPCK is a subsidiary routine called by B3INT.
C
C   BUPCK converts an array packed as though dimensioned as A(NX,NY,NZ)
C   to one dimensioned as A(LD1,LD2,NZ), where LD1.GE.NX and LD2.GE.NY.
C
C   The routine assumes that the input array is packed, that is, the
C   first NX*NY*NZ locations contain data.  It then unpacks the data
C   so that the calling program can dimension it as A(LD1,LD2,NZ).
C   The operation is done in place.
C
C
C   I N P U T
C   ---------
C
C   A       Real 1D array (length at least LD1*LD2*NZ)
C           On input, contains the array to be converted.
C           On output, contains the converted array.
C
C   NX      Integer (NX.GT.0)
C           First dimension of input array.
C
C   NY      Integer (NY.GT.0)
C           Second dimension of input array.
C
C   NZ      Integer (NZ.GT.0)
C           Third dimension of input array.
C
C   LD1     Integer (LD1.GE.NX)
C           First dimension of output array.
C
C   LD2     Integer (LD2.GE.NY)
C           Second dimension of output array.
C
C
C   CAUTION: The constraints on the input variables are not checked by
C            this routine.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   900910  DATE WRITTEN
C   930510  Clarified purpose.  (FNF)
C   930618  Reformatted the AUTHOR section.  (WRB)
C***END PROLOGUE  BUPCK
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          LD1, LD2, NX, NY, NZ
      real             A(*)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I, IF2, IF3, IT2, IT3, J, K
C
C***FIRST EXECUTABLE STATEMENT  BUPCK
      IF (( NZ .EQ.  1) .AND. (LD1 .EQ. NX))  GO TO 900
      IF ((LD1 .EQ. NX) .AND. (LD2 .EQ. NY))  GO TO 900
C
      DO 300 K=NZ,1,-1
         IF3 = (K-1)*NY - 1
         IT3 = (K-1)*LD2 - 1
         DO 200 J=NY,1,-1
            IF2 = (J + IF3)*NX
            IT2 = (J + IT3)*LD1
            DO 100 I=NX,1,-1
               A(I+IT2) = A(I+IF2)
  100       CONTINUE
  200    CONTINUE
  300 CONTINUE
C
  900 CONTINUE
      RETURN
      END

*DECK B3VAL
      real FUNCTION B3VAL (XVAL, YVAL, ZVAL, IDX, IDY, IDZ,
     +   TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, FCN, LDF1, LDF2,
     +   ICONT, IWORK, WORK, IFLAG)
C***BEGIN PROLOGUE  B3VAL
C***PURPOSE  B3VAL evaluates the three-dimensional piecewise polynomial
C            interpolating function constructed by the routine B3INT.
C            Either function values or partial derivative values may be
C            be requested.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (B3VAL-S, DB3VAL-D)
C***KEYWORDS  THREE DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS, EVALUATION, DIFFERENTIATION
C***AUTHOR  Boisvert, R. F., (NIST)
C             Computing and Applied Mathematics Laboratory
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   B3VAL evaluates the tensor product piecewise polynomial interpolant
C   constructed by the routine B3INT or one of its derivatives at the
C   point (XVAL,YVAL,ZVAL).  The variables IDX, IDY and IDZ indicate
C   which function is to be evaluated.  If B(x,y,z) is the interpolant,
C
C                            IDX+IDY+IDZ
C                           d
C             B1VAhL  =   ----------------- B (XVAL,YVAL,ZVAL)
C                          IDX   IDY   IDZ
C                        dx    dy    dz
C
C   Thus, to evaluate the interpolant itself, set IDX=IDY=IDZ=0. To get
C   the first partial derivative with respect to x, set IDX=1 and
C   IDY=IDZ=0, and so on.

C   Since B is a piecewise polynomial of degree KX-1 in x, KY-1 in y and
C   KZ-1 in z, B3VAL returns zero whenever IDX.GE.KX or IDY.GE.KY or
C   IDZ.GE.KZ.  B3VAL also returns zero if (XVAL,YVAL,ZVAL) is out of
C   range, that is, if
C
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR.
C            ZVAL.LT.TZ(1) .OR. ZVAL.GT.ZY(NZ+KZ).
C
C   If the knots TX, TY, and TZ were chosen by B3INT, then this is
C   equivalent to
C
C              XVAL.LT.X(1) .OR. XVAL.GT.X(NX) .OR.
C              YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY) .OR.
C              ZVAL.LT.Z(1) .OR. ZVAL.GT.Z(NZ).
C
C   The input quantities TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, and FCN
C   should be unchanged since the call to B3INT which produced FCN.
C
C   Note that the derivative functions computed by B3VAL may be
C   discontinuous when x, y or z correspond to knots.  (When B3INT
C   selects knots this occurs only when IDX=KX-1, IDY=KY-1 or IDZ=KZ-1.)
C   In these cases B3VAL returns right limiting values (right
C   derivatives), except at the rightmost knot where left limiting
C   values are returned.
C
C
C   I N P U T
C   ---------
C
C   XVAL    Real scalar
C           X coordinate of evaluation point.
C
C   YVAL    Real scalar
C           Y coordinate of evaluation point.
C
C   ZVAL    Real scalar
C           Z coordinate of evaluation point.
C
C   IDX     Integer scalar  (IDX .GE. 0)
C           Indicates the x derivative of piecewise polynomial to
C           evaluate: IDX=J for the Jth partial derivative with respect
C           to x.  B3VAL will return 0 if IDX.GE.KX.
C
C   IDY     Integer scalar  (IDY .GE. 0)
C           Indicates the y derivative of piecewise polynomial to
C           evaluate: IDY=J for the Jth partial derivative with respect
C           to y.  B3VAL will return 0 if IDY.GE.KY.
C
C   IDZ     Integer scalar  (IDZ .GE. 0)
C           Indicates the z derivative of piecewise polynomial to
C           evaluate: IDZ=J for the Jth partial derivative with respect
C           to z.  B3VAL will return 0 if IDZ.GE.KZ.
C
C   TX      Real 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Unchanged since the call to B3INT which
C           produced FCN.)
C
C   TY      Real 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Unchanged since the call to B3INT which
C           produced FCN.)
C
C   TZ      Real 1D array (size NZ+KZ)
C           Sequence of knots defining the piecewise polynomial in
C           the z direction.  (Unchanged since the call to B3INT which
C           produced FCN.)
C
C   NX      Integer scalar (NX .GE. KX)
C           The number of interpolation points in x. (Unchanged since
C           the call to B3INT which produced FCN.)
C
C   NY      Integer scalar (NY .GE. KY)
C           The number of interpolation points in y. (Unchanged since
C           the call to B3INT which produced FCN.)
C
C   NZ      Integer scalar (NZ .GE. KZ)
C           The number of interpolation points in z. (Unchanged since
C           the call to B3INT which produced FCN.)
C
C   KX      Integer scalar (KX .GE. 2)
C           Order of polynomial pieces in x. (Unchanged since the call
C           to B3INT which produced FCN.)
C
C   KY      Integer scalar (KY .GE. 2)
C           Order of polynomial pieces in y. (Unchanged since the call
C           to B3INT which produced FCN.)
C
C   KZ      Integer scalar (KZ .GE. 2)
C           Order of polynomial pieces in z. (Unchanged since the call
C           to B3INT which produced FCN.)
C
C   FCN     Real 3D array (size LDF1 by LDF2 by NZ)
C           The B-spline coefficients computed by B3INT.
C
C   LDF1    Integer scalar (LDF1 .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C   LDF2    Integer scalar (LDF2 .GE. NY)
C           The actual second dimension of FCN used in the calling
C           program.
C
C
C   I N P U T / O U T P U T
C   -----------------------
C
C   ICONT   Integer scalar
C           An input flag.  Set ICONT=0 on the first call. Subsequent
C           calls evaluating the same piecewise polynomial function
C           (i.e., the same TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, and
C           FCN) should have ICONT=1.  (See notes on efficiency below.)
C           As a convenience, B3VAL sets ICONT=1 on output.
C
C   IWORK   Integer 1D array (size 10)
C           A working storage array.  NOTE: This array is used for
C           communication between successive calls of B3VAL.  It must
C           not be modified between successive calls when ICONT=1.
C           Different interpolants require different work arrays.
C
C   WORK    Real 1D array (size KY*KZ + 3*max(KX,KY,KZ) + KZ + 2)
C           A working storage array.  NOTE: This array is used for
C           communication between successive calls of B3VAL.  It must
C           not be modified between successive calls when ICONT=1.
C           Different interpolants require different work arrays.
C
C
C   O U T P U T
C   -----------
C
C   IFLAG   Integer scalar.
C           An output flag.  Possible values are
C           0 : successful execution
C           1 : KX out of range
C           2 : NX or LDF1 out of range
C           3 : KY out of range
C           4 : NY or LDF2 out of range
C           5 : KZ out of range
C           6 : NZ out of range
C           7 : IDX, IDY or IDZ out of range
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C
C   N O T E S   O N   E F F I C I E N C Y
C   -------------------------------------
C
C   B3VAL is designed so that it will be reasonably efficient when it is
C   being used to evaluate the fitted function on a grid.  The most
C   favorable situation occurs when B3VAL is called repeatedly with the
C   same XVAL, the same IDX and IDY, and slowly increasing YVAL and ZVAL.
C   When calling B3VAL in a triple loop over x, y and z :
C
C     (a) vary x in the OUTER loop, and
C     (b) vary z in the INNER loop, and
C     (c) use separate loops or distinct ICONT,IWORK,WORK for
C         each (IDX,IDY) pair
C
C   A typical example of this usage is as follows.
C
C      ICONT = 0
C      DO 30 I=1,NX
C         DO 20 J=1,NY
C            DO 10 K=1,NZ
C               G(I,J) = B3VAL(X(I),Y(J),Z(K),IDX,IDY,IDZ,TX,TY,TZ,
C     +                        NX,NY,NZ,KX,KY,KZ,FCN,LDF1,LDF2,
C     +                        ICONT,IWORK,WORK,IFLAG)
C               IF (IFLAG .NE. 0) GO TO 9999
C   10       CONTINUE
C   20    CONTINUE
C   30 CONTINUE
C
C   Note that ICONT=0 initially.  This signals B3VAL that it is starting
C   a new function.  B3VAL sets ICONT=1 on output, so ICONT=1 on all
C   subsequent calls, which tells B3VAL to attempt to reuse information
C   from the previous call to speed the evaluation process.
C
C
C   B A C K G R O U N D
C   -------------------
C
C   B3VAL evaluates the following function or its partial derivatives
C
C                      NX   NY   NZ
C        B(x,y,z)  =  SUM  SUM  SUM  FCN(i,j,k) U (x) V (y) W (z)
C                     i=1  j=1  k=1              i     j     k
C   where
C
C     U (x)  is the ith (one-dimensional) B-spline basis function
C      i     defined by NX, KX, and TX,
C
C     V (y)  is the jth (one-dimensional) B-spline basis function
C      j     defined by NY, KY, and TY, and
C
C     W (z)  is the kth (one-dimensional) B-spline basis function
C      k     defined by NZ, KZ, and TZ.
C
C   See (de Boor, 1978) for a description of the B-spline basis.
C
C  The summation above can be rewritten as
C
C                                  NZ
C  (1)               B(x,y,z)  =  SUM  b  V (z),
C                                 k=1   k  k
C  where
C                         NY
C  (2)            b  =   SUM  c   V (y),   k = 1,..,NZ,
C                  k     j=1   jk  j
C  and
C                     NX
C  (3)        c  =   SUM  FCN(i,j,k) U (x),   j = 1,..,NY,
C              jk    i=1              i       k = 1,..,NZ.
C
C  Note that each summation is the evaluation of a one-dimensional
C  B-spline, which can be done by calls to B1VAhL.  At most KZ basis
C  functions in (1) are nonzero.  The indices of these functions are
C  determined, and only the coefficients b(k) which multiply the
C  nonzero basis functions are be computed using (2).  We also
C  observe that at most KY basis functions in (2) are nonzero.  The
C  indices of these functions also are determined and only those
C  coefficients c(j,k) which multiply nonzero basis functions are
C  computed using (3).
C
C
C   A C K N O W L E D G E M E N T
C   -----------------------------
C
C   Thanks to Fred Fritsch of Lawrence Livermore National Laboratory
C   for his critical review of this code which lead to many
C   improvements.
C
C***SEE ALSO  B3INT
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  B1VAhL, INTRhV, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930510  Eliminated multiple calls to XERMShG.  (FNF)
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   940202  Substantial revision.  (RFB)
C           The revision allows for increased efficiency when evaluating
C           on a fine grid.  This is accomplished by reusing information
C           computed on previous calls.  Changes visible to users are :
C           (a) new arguments (ICONT and IWORK), and
C           (b) increase in size of WORK.
C***END PROLOGUE  B3VAL
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          IDX, IDY, IDZ, NX, NY, NZ, KX, KY, KZ, LDF1,
     +                 LDF2, ICONT, IWORK(10), IFLAG
      real             XVAL, YVAL, ZVAL, TX(NX+KX), TY(NY+KY),
     +                 TZ(NZ+KZ), FCN(LDF1,LDF2,NZ), WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'SLATEC', SUBROU = 'B3VAL')
      CHARACTER*50     MESSAG
      LOGICAL          NEEDCO
      INTEGER          I, IDLAST, IERR, IHAVEX, IHAVEY, ILOY, ILOZ,
     +                 INX, INY, INZ, IPOS, IW, IXY, IZ, J, JDLAST,
     +                 JMAX, JMIN, KMAX, KMIN, K, LEFTY, LEFTZ,
     +                 LYLAST, LZLAST, MFLAG
      real             B1VAhL, XLAST, YLAST
C
      EXTERNAL         INTRhV, B1VAhL, XERMShG
C
C***FIRST EXECUTABLE STATEMENT  B3VAL
      B3VAL = 0
C
C  ------------------------------
C  RESET STATE FROM PREVIOUS CALL
C  ------------------------------
C
C  These state variables indicate state of B3VAL at end of last call.
C  They tell us the range of validity of the coefficients b(k) in (2)
C  and c(j,k) in (3) computed in the previous call and saved in WORK.
C  We want to avoid recomputing these if possible.
C
C  IWORK( 1) == IHAVEX == 1 if coefficients c(j,k) in (3) are available
C  IWORK( 2) == IHAVEY == 1 if coefficients b(k) in (2) are available
C  IWORK( 3) == INX    == knot interval of XVAL (for B1VAhL)
C  IWORK( 4) == ILOY   == knot interval of YVAL (for INTRhV)
C  IWORK( 5) == ILOZ   == knot interval of ZVAL (for INTRhV)
C  IWORK( 6) == LYLAST == value of LEFTY from INTRhV for YVAL
C  IWORK( 7) == LZLAST == value of LEFTZ from INTRhV for ZVAL
C  IWORK( 8) == IDLAST == value of IDX
C  IWORK( 9) == JDLAST == value of IDY
C  IWORK(10) == KMIN   == smallest  index k of nonzero basis b(k) at ZVAL
C   WORK( 1) == XLAST  == value of x variable
C   WORK( 2) == YLAST  == value of y variable
C
      IF ( ICONT .NE. 1 ) THEN
         IHAVEX = 0
         IHAVEY = 0
         INX = 1
         ILOY = 1
         ILOZ = 1
         JMIN = 0
         KMIN = 0
      ELSE
         IHAVEX = IWORK(1)
         IHAVEY = IWORK(2)
         INX    = IWORK(3)
         ILOY   = IWORK(4)
         ILOZ   = IWORK(5)
         LYLAST = IWORK(6)
         LZLAST = IWORK(7)
         IDLAST = IWORK(8)
         JDLAST = IWORK(9)
         KMIN   = IWORK(10)
         XLAST  = WORK(1)
         YLAST  = WORK(2)
      ENDIF
C
C     ... set LEFTY, LEFTZ so that they have values when state is
C         saved after premature termination
C
      LEFTY = 0
      LEFTZ = 0
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IFLAG = 0
      IF (KX .LT. 1) THEN
         IFLAG = 1
         MESSAG = 'KX IS OUT OF RANGE'
      ELSE IF ((NX .LT. KX) .OR. (NX .GT. LDF1)) THEN
         IFLAG = 2
         MESSAG = 'NX OR LDF1 IS OUT OF RANGE'
      ELSE IF (KY .LT. 1) THEN
         IFLAG = 3
         MESSAG = 'KY IS OUT OF RANGE'
      ELSE IF ((NY .LT. KY) .OR. (NY .GT. LDF2)) THEN
         IFLAG = 4
         MESSAG = 'NY OR LDF2 IS OUT OF RANGE'
      ELSE IF (KZ .LT. 1) THEN
         IFLAG = 5
         MESSAG = 'KZ IS OUT OF RANGE'
      ELSE IF (NZ .LT. KZ) THEN
         IFLAG = 6
         MESSAG = 'NZ IS OUT OF RANGE'
      ELSE IF ((IDX .LT. 0) .OR. (IDY .LT. 0) .OR. (IDZ .LT. 0)) THEN
         IFLAG = 7
         MESSAG = 'IDX, IDY OR IDZ IS OUT OF RANGE'
      ELSE
         IFLAG = 0
      ENDIF
C
      IF ( IFLAG .NE. 0 ) THEN
         CALL XERMShG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
         GO TO 9999
      ENDIF
C
C     ... QUICK RETURN WHEN SPLINE EVALUATES TO ZERO
C
      IF ( (IDX  .GE. KX       ) .OR.
     +     (IDY  .GE. KY       ) .OR.
     +     (IDZ  .GT. KZ       ) .OR.
     +     (XVAL .LT. TX(1)    ) .OR.
     +     (XVAL .GT. TX(NX+KX)) .OR.
     +     (YVAL .LT. TY(1)    ) .OR.
     +     (YVAL .GT. TY(NY+KY)) .OR.
     +     (ZVAL .LT. TZ(1)    ) .OR.
     +     (ZVAL .GT. TZ(NZ+KZ))      ) GO TO 9999
C
C  ---------------------
C  EVALUATE THE B-SPLINE
C  ---------------------
C
C     ... PARTITION WORK ARRAY
C
      IXY = 3
      IZ = IXY + KY*KZ
      IW = IZ  + KZ
C
C     ... FIND KNOT INTERVAL CONTAINING ZVAL
C
      CALL INTRhV(TZ,NZ+KZ,ZVAL,ILOZ,LEFTZ,MFLAG)
      IF (MFLAG .NE. 0) THEN
C
C        ... ZVAL .EQ. T(NZ+KZ),  ADJUST TO GET LEFT LIMITING VALUE
C
   10    CONTINUE
         LEFTZ = LEFTZ - 1
         IF (ZVAL .EQ. TZ(LEFTZ)) GO TO 10
      ENDIF
C
C     ... FIND KNOT INTERVAL CONTAINING YVAL
C
      CALL INTRhV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0) THEN
C
C        ... YVAL .EQ. T(NY+KY),  ADJUST TO GET LEFT LIMITING VALUE
C
   20    CONTINUE
         LEFTY = LEFTY - 1
         IF (YVAL .EQ. TY(LEFTY)) GO TO 20
      ENDIF
C
C     ... CHECK IF WE ALREADY HAVE COEFFICIENTS c(j,k) IN (3)
C
      NEEDCO = (IHAVEX .EQ. 0     ) .OR.
     +         (LEFTZ  .NE. LZLAST) .OR.
     +         (LEFTY  .NE. LYLAST) .OR.
     +         (IDX    .NE. IDLAST) .OR.
     +         (XVAL   .NE. XLAST )
C
      IF ( NEEDCO ) THEN
C
C        ... FIND RANGE OF INDICES OF NONZERO BASIS FUNCTIONS IN (1)
C            (KMIN,KMAX) = (SMALLEST,LARGEST)
C
         IF (LEFTZ .LT. KZ) THEN
            KMIN = 1
            KMAX = KZ
         ELSEIF (LEFTZ .GT. NZ) THEN
            KMIN = NZ - KZ + 1
            KMAX = NZ
         ELSE
            KMIN = LEFTZ-KZ+1
            KMAX = LEFTZ
         ENDIF
C
C        ... FIND RANGE OF INDICES OF NONZERO BASIS FUNCTIONS IN (2)
C            (JMIN,JMAX) = (SMALLEST,LARGEST)
C
         IF (LEFTY .LT. KY) THEN
            JMIN = 1
            JMAX = KY
         ELSEIF (LEFTY .GT. NY) THEN
            JMIN = NY - KY + 1
            JMAX = NY
         ELSE
            JMIN = LEFTY-KY+1
            JMAX = LEFTY
         ENDIF
C
C        ... GET COEFFICIENTS OF NONZERO BASIS FUNCTIONS IN (3)
C
         IPOS = IXY - 1
         DO 50 K=KMIN,KMAX
            DO 50 J=JMIN,JMAX
               IPOS = IPOS + 1
               WORK(IPOS) = B1VAhL(XVAL,IDX,TX,NX,KX,FCN(1,J,K),INX,
     +                            WORK(IW),IERR)
   50    CONTINUE
         IHAVEX = 1
C
      ENDIF
C
C     ... CHECK IF WE ALREADY HAVE COEFFICIENTS b(j) IN (2)
C
      NEEDCO = NEEDCO               .OR.
     +         (IHAVEY .EQ. 0     ) .OR.
     +         (IDY    .NE. JDLAST) .OR.
     +         (YVAL   .NE. YLAST )
C
      IF ( NEEDCO ) THEN
C
C        ... FIND MIN INDEX OF NONZERO BASIS FUNCTIONS IN (2)
C
         IF (LEFTY .LT. KY) THEN
            JMIN = 1
         ELSEIF (LEFTY .GT. NY) THEN
            JMIN = NY - KY + 1
         ELSE
            JMIN = LEFTY-KY+1
         ENDIF
C
C        ... GET COEFFICIENTS OF NONZERO BASIS FUNCTIONS IN (2)
C
         INY = 1
         I = IXY - KY
         J = IZ - 1
         DO 60 K=1,KZ
            I = I + KY
            J = J + 1
            WORK(J) = B1VAhL(YVAL,IDY,TY(JMIN),KY,KY,WORK(I),INY,
     +                      WORK(IW),IERR)
  60     CONTINUE
         IHAVEY = 1
C
      ENDIF
C
C     ... EVALUATE INTERPOLANT USING (1)
C
      INZ = 1
      B3VAL = B1VAhL(ZVAL,IDZ,TZ(KMIN),KZ,KZ,WORK(IZ),INZ,WORK(IW),
     +              IERR)
C
C  -------------------
C  SAVE STATE AND EXIT
C  -------------------
C
 9999 CONTINUE
      IWORK(1)  = IHAVEX
      IWORK(2)  = IHAVEY
      IWORK(3)  = INX
      IWORK(4)  = ILOY
      IWORK(5)  = ILOZ
      IWORK(6)  = LEFTY
      IWORK(7)  = LEFTZ
      IWORK(8)  = IDX
      IWORK(9)  = IDY
      IWORK(10) = KMIN
      WORK(1) = XVAL
      WORK(2) = YVAL
      ICONT = 1
C
      RETURN
      END

*DECK BINTh4
      SUBROUTINE BINTh4 (X, Y, NDATA, IBCL, IBCR, FBCL, FBCR, KNTOPT, T,
     +   BCOEF, N, K, W)
C***BEGIN PROLOGUE  BINTh4
C***PURPOSE  Compute the B-representation of a cubic spline
C            which interpolates given data.
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (BINTh4-S, DBINTh4-D)
C***KEYWORDS  B-SPLINE, CUBIC SPLINES, DATA FITTING, INTERPOLATION
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         BINTh4 computes the B representation (T,BCOEF,N,K) of a
C         cubic spline (K=4) which interpolates data (X(I)),Y(I))),
C         I=1,NDATA.  Parameters IBCL, IBCR, FBCL, FBCR allow the
C         specification of the spline first or second derivative at
C         both X(1) and X(NDATA).  When this data is not specified
C         by the problem, it is common practice to use a natural
C         spline by setting second derivatives at X(1) and X(NDATA)
C         to zero (IBCL=IBCR=2,FBCL=FBCR=0.0).  The spline is defined on
C         T(4) .LE. X .LE. T(N+1) with (ordered) interior knots at X(I))
C         values where N=NDATA+2.  The knots T(1), T(2), T(3) lie to
C         the left of T(4)=X(1) and the knots T(N+2), T(N+3), T(N+4)
C         lie to the right of T(N+1)=X(NDATA) in increasing order.  If
C         no extrapolation outside (X(1),X(NDATA)) is anticipated, the
C         knots T(1)=T(2)=T(3)=T(4)=X(1) and T(N+2)=T(N+3)=T(N+4)=
C         T(N+1)=X(NDATA) can be specified by KNTOPT=1.  KNTOPT=2
C         selects a knot placement for T(1), T(2), T(3) to make the
C         first 7 knots symmetric about T(4)=X(1) and similarly for
C         T(N+2), T(N+3), T(N+4) about T(N+1)=X(NDATA).  KNTOPT=3
C         allows the user to make his own selection, in increasing
C         order, for T(1), T(2), T(3) to the left of X(1) and T(N+2),
C         T(N+3), T(N+4) to the right of X(NDATA) in the work array
C         W(1) through W(6).  In any case, the interpolation on
C         T(4) .LE. X .LE. T(N+1) by using function BVALU is unique
C         for given boundary conditions.
C
C     Description of Arguments
C         Input
C           X      - X vector of abscissae of length NDATA, distinct
C                    and in increasing order
C           Y      - Y vector of ordinates of length NDATA
C           NDATA  - number of data points, NDATA .GE. 2
C           IBCL   - selection parameter for left boundary condition
C                    IBCL = 1 constrain the first derivative at
C                             X(1) to FBCL
C                         = 2 constrain the second derivative at
C                             X(1) to FBCL
C           IBCR   - selection parameter for right boundary condition
C                    IBCR = 1 constrain first derivative at
C                             X(NDATA) to FBCR
C                    IBCR = 2 constrain second derivative at
C                             X(NDATA) to FBCR
C           FBCL   - left boundary values governed by IBCL
C           FBCR   - right boundary values governed by IBCR
C           KNTOPT - knot selection parameter
C                    KNTOPT = 1 sets knot multiplicity at T(4) and
C                               T(N+1) to 4
C                           = 2 sets a symmetric placement of knots
C                               about T(4) and T(N+1)
C                           = 3 sets TNP)=WNP) and T(N+1+I)=w(3+I),I=1,3
C                               where WNP),I=1,6 is supplied by the user
C           W      - work array of dimension at least 5*(NDATA+2)
C                    if KNTOPT=3, then W(1),W(2),W(3) are knot values to
C                    the left of X(1) and W(4),W(5),W(6) are knot
C                    values to the right of X(NDATA) in increasing
C                    order to be supplied by the user
C
C         Output
C           T      - knot array of length N+4
C           BCOEF  - B-spline coefficient array of length N
C           N      - number of coefficients, N=NDATA+2
C           K      - order of spline, K=4
C
C     Error Conditions
C         Improper  input is a fatal error
C         Singular system of equations is a fatal error
C
C***REFERENCES  D. E. Amos, Computation with splines and B-splines,
C                 Report SAND78-1968, Sandia Laboratories, March 1979.
C               Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  BNFAhC, BNSLhV, BSPVhD, EPSILON, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BINTh4
C
      INTEGER I, IBCL, IBCR, IFLAG, ILB, ILEFT, IT, IUB, IW, IWP, J,
     1 JW, K, KNTOPT, N, NDATA, NDM, NP, NWROW
      real BCOEF,FBCL,FBCR,T, TOL,TXN,TX1,VNIKX,W,WDTOL,WORK,X, XL,
     1 Y
C      real EPSILON
      DIMENSION X(*), Y(*), T(*), BCOEF(*), W(5,*), VNIKX(4,4), WORK(15)
C***FIRST EXECUTABLE STATEMENT  BINTh4
C      WDTOL = R1MACH9(4)
      WDTOL = EPSILON(BCOEF)
      TOL = SQRT(WDTOL)
      IF (NDATA.LT.2) GO TO 200
      NDM = NDATA - 1
      DO 10 I=1,NDM
        IF (X(I).GE.X(I+1)) GO TO 210
   10 CONTINUE
      IF (IBCL.LT.1 .OR. IBCL.GT.2) GO TO 220
      IF (IBCR.LT.1 .OR. IBCR.GT.2) GO TO 230
      IF (KNTOPT.LT.1 .OR. KNTOPT.GT.3) GO TO 240
      K = 4
      N = NDATA + 2
      NP = N + 1
      DO 20 I=1,NDATA
        T(I+3) = X(I)
   20 CONTINUE
      GO TO (30, 50, 90), KNTOPT
C     SET UP KNOT ARRAY WITH MULTIPLICITY 4 AT X(1) AND X(NDATA)
   30 CONTINUE
      DO 40 I=1,3
        T(4-I) = X(1)
        T(NP+I) = X(NDATA)
   40 CONTINUE
      GO TO 110
C     SET UP KNOT ARRAY WITH SYMMETRIC PLACEMENT ABOUT END POINTS
   50 CONTINUE
      IF (NDATA.GT.3) GO TO 70
      XL = (X(NDATA)-X(1))/3.0E0
      DO 60 I=1,3
        T(4-I) = T(5-I) - XL
        T(NP+I) = T(NP+I-1) + XL
   60 CONTINUE
      GO TO 110
   70 CONTINUE
      TX1 = X(1) + X(1)
      TXN = X(NDATA) + X(NDATA)
      DO 80 I=1,3
        T(4-I) = TX1 - X(I+1)
        T(NP+I) = TXN - X(NDATA-I)
   80 CONTINUE
      GO TO 110
C     SET UP KNOT ARRAY LESS THAN X(1) AND GREATER THAN X(NDATA) TO BE
C     SUPPLIED BY USER IN WORK LOCATIONS W(1) THROUGH W(6) WHEN KNTOPT=3
   90 CONTINUE
      DO 100 I=1,3
        T(4-I) = W(4-I,1)
        JW = MAX(1,I-1)
        IW = MOD(I+2,5)+1
        T(NP+I) = W(IW,JW)
        IF (T(4-I).GT.T(5-I)) GO TO 250
        IF (T(NP+I).LT.T(NP+I-1)) GO TO 250
  100 CONTINUE
  110 CONTINUE
C
      DO 130 I=1,5
        DO 120 J=1,N
          W(I,J) = 0.0E0
  120   CONTINUE
  130 CONTINUE
C     SET UP LEFT INTERPOLATION POINT AND LEFT BOUNDARY CONDITION FOR
C     RIGHT LIMITS
      IT = IBCL + 1
      CALL BSPVhD(T, K, IT, X(1), K, 4, VNIKX, WORK)
      IW = 0
      IF (ABS(VNIKX(3,1)).LT.TOL) IW = 1
      DO 140 J=1,3
        W(J+1,4-J) = VNIKX(4-J,IT)
        W(J,4-J) = VNIKX(4-J,1)
  140 CONTINUE
      BCOEF(1) = Y(1)
      BCOEF(2) = FBCL
C     SET UP INTERPOLATION EQUATIONS FOR POINTS I=2 TO I=NDATA-1
      ILEFT = 4
      IF (NDM.LT.2) GO TO 170
      DO 160 I=2,NDM
        ILEFT = ILEFT + 1
        CALL BSPVhD(T, K, 1, X(I), ILEFT, 4, VNIKX, WORK)
        DO 150 J=1,3
          W(J+1,3+I-J) = VNIKX(4-J,1)
  150   CONTINUE
        BCOEF(I+1) = Y(I)
  160 CONTINUE
C     SET UP RIGHT INTERPOLATION POINT AND RIGHT BOUNDARY CONDITION FOR
C     LEFT LIMITS(ILEFT IS ASSOCIATED WITH T(N)=X(NDATA-1))
  170 CONTINUE
      IT = IBCR + 1
      CALL BSPVhD(T, K, IT, X(NDATA), ILEFT, 4, VNIKX, WORK)
      JW = 0
      IF (ABS(VNIKX(2,1)).LT.TOL) JW = 1
      DO 180 J=1,3
        W(J+1,3+NDATA-J) = VNIKX(5-J,IT)
        W(J+2,3+NDATA-J) = VNIKX(5-J,1)
  180 CONTINUE
      BCOEF(N-1) = FBCR
      BCOEF(N) = Y(NDATA)
C     SOLVE SYSTEM OF EQUATIONS
      ILB = 2 - JW
      IUB = 2 - IW
      NWROW = 5
      IWP = IW + 1
      CALL BNFAhC(W(IWP,1), NWROW, N, ILB, IUB, IFLAG)
      IF (IFLAG.EQ.2) GO TO 190
      CALL BNSLhV(W(IWP,1), NWROW, N, ILB, IUB, BCOEF)
      RETURN
C
C
  190 CONTINUE
      CALL XERMShG ('SLATEC', 'BINTh4',
     +   'THE SYSTEM OF EQUATIONS IS SINGULAR', 2, 1)
      RETURN
  200 CONTINUE
      CALL XERMShG ('SLATEC', 'BINTh4', 'NDATA IS LESS THAN 2', 2, 1)
      RETURN
  210 CONTINUE
      CALL XERMShG ('SLATEC', 'BINTh4',
     +   'X VALUES ARE NOT DISTINCT OR NOT ORDERED', 2, 1)
      RETURN
  220 CONTINUE
      CALL XERMShG ('SLATEC', 'BINTh4', 'IBCL IS NOT 1 OR 2', 2, 1)
      RETURN
  230 CONTINUE
      CALL XERMShG ('SLATEC', 'BINTh4', 'IBCR IS NOT 1 OR 2', 2, 1)
      RETURN
  240 CONTINUE
      CALL XERMShG ('SLATEC', 'BINTh4', 'KNTOPT IS NOT 1, 2, OR 3', 2,1)
      RETURN
  250 CONTINUE
      CALL XERMShG ('SLATEC', 'BINTh4',
     +   'KNOT INPUT THROUGH W ARRAY IS NOT ORDERED PROPERLY', 2, 1)
      RETURN
      END

*DECK BSPVhD
      SUBROUTINE BSPVhD (T, K, NDERIV, X, ILEFT, LDVNIK, VNIKX, WORK)
C***BEGIN PROLOGUE  BSPVhD
C***PURPOSE  Calculate the value and all derivatives of order less than
C            NDERIV of all basis functions which do not vanish at X.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (BSPVhD-S, DBSPVhD-D)
C***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C         BSPVhD is the BSPLVD routine of the reference.
C
C         BSPVhD calculates the value and all derivatives of order
C         less than NDERIV of all basis functions which do not
C         (possibly) vanish at X.  ILEFT is input such that
C         T(ILEFT) .LE. X .LT. T(ILEFT+1).  A call to INTRhV(T,N+1,X,
C         ILO,ILEFT,MFLAG) will produce the proper ILEFT.  The output of
C         BSPVhD is a matrix VNIKX(I,J) of dimension at least (K,NDERIV)
C         whose columns contain the K nonzero basis functions and
C         their NDERIV-1 right derivatives at X, I=1,K, J=1,NDERIV.
C         These basis functions have indices ILEFT-K+I, I=1,K,
C         K .LE. ILEFT .LE. N. The nonzero part of the I-th basis
C         function lies in (T(I),T(I+K)), I=1,N.
C
C         If X=T(ILEFT+1) then VNIKX contains left limiting values
C         (left derivatives) at T(ILEFT+1).  In particular, ILEFT = N
C         produces left limiting values at the right end point
C         X=T(N+1). To obtain left limiting values at T(I), I=K+1,N+1,
C         set X= next lower distinct knot, call INTRhV to get ILEFT,
C         set X=T(I), and then call BSPVhD.
C
C     Description of Arguments
C         Input
C          T       - knot vector of length N+K, where
C                    N = number of B-spline basis functions
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          NDERIV  - number of derivatives = NDERIV-1,
C                    1 .LE. NDERIV .LE. K
C          X       - argument of basis functions,
C                    T(K) .LE. X .LE. T(N+1)
C          ILEFT   - largest integer such that
C                    T(ILEFT) .LE. X .LT. T(ILEFT+1)
C          LDVNIK  - leading dimension of matrix VNIKX
C
C         Output
C          VNIKX   - matrix of dimension at least (K,NDERIV) contain-
C                    ing the nonzero basis functions at X and their
C                    derivatives columnwise.
C          WORK    - a work vector of length (K+1)*(K+2)/2
C
C     Error Conditions
C         Improper input is a fatal error
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  BSPVhN, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BSPVhD
C
      implicit none
      INTEGER I,IDERIV,ILEFT,IPKMD,J,JJ,JLOW,JM,JP1MID,K,KMD, KP1, L,
     1 LDUMMY, M, MHIGH, NDERIV, LDVNIK, IWORK
      real FACTOR, FKMD, T, V, VNIKX, WORK, X
C     DIMENSION T(ILEFT+K), WORK((K+1)*(K+2)/2)
C     A(I,J) = WORK(I+J*(J+1)/2),  I=1,J+1  J=1,K-1
C     A(I,K) = W0RK(I+K*(K-1)/2)  I=1.K
C     WORK(1) AND WORK((K+1)*(K+2)/2) ARE NOT USED.
      DIMENSION T(*), VNIKX(LDVNIK,*), WORK(*)
C***FIRST EXECUTABLE STATEMENT  BSPVhD
      IF(K.LT.1) GO TO 200
      IF(NDERIV.LT.1 .OR. NDERIV.GT.K) GO TO 205
      IF(LDVNIK.LT.K) GO TO 210
      IDERIV = NDERIV
      KP1 = K + 1
      JJ = KP1 - IDERIV
      CALL BSPVhN(T, JJ, K, 1, X, ILEFT, VNIKX, WORK, IWORK)
      IF (IDERIV.EQ.1) GO TO 100
      MHIGH = IDERIV
      DO 20 M=2,MHIGH
        JP1MID = 1
        DO 10 J=IDERIV,K
          VNIKX(J,IDERIV) = VNIKX(JP1MID,1)
          JP1MID = JP1MID + 1
   10   CONTINUE
        IDERIV = IDERIV - 1
        JJ = KP1 - IDERIV
        CALL BSPVhN(T, JJ, K, 2, X, ILEFT, VNIKX, WORK, IWORK)
   20 CONTINUE
C
      JM = KP1*(KP1+1)/2
      DO 30 L = 1,JM
        WORK(L) = 0.0E0
   30 CONTINUE
C     A(I,I) = WORK(I*(I+3)/2) = 1.0       I = 1,K
      L = 2
      J = 0
      DO 40 I = 1,K
        J = J + L
        WORK(J) = 1.0E0
        L = L + 1
   40 CONTINUE
      KMD = K
      DO 90 M=2,MHIGH
        KMD = KMD - 1
        FKMD = KMD
        I = ILEFT
        J = K
        JJ = J*(J+1)/2
        JM = JJ - J
        DO 60 LDUMMY=1,KMD
          IPKMD = I + KMD
          FACTOR = FKMD/(T(IPKMD)-T(I))
          DO 50 L=1,J
            WORK(L+JJ) = (WORK(L+JJ)-WORK(L+JM))*FACTOR
   50     CONTINUE
          I = I - 1
          J = J - 1
          JJ = JM
          JM = JM - J
   60   CONTINUE
C
        DO 80 I=1,K
          V = 0.0E0
          JLOW = MAX(I,M)
          JJ = JLOW*(JLOW+1)/2
          DO 70 J=JLOW,K
            V = WORK(I+JJ)*VNIKX(J,M) + V
            JJ = JJ + J + 1
   70     CONTINUE
          VNIKX(I,M) = V
   80   CONTINUE
   90 CONTINUE
  100 RETURN
C
C
  200 CONTINUE
      CALL XERMShG ('SLATEC', 'BSPVhD', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  205 CONTINUE
      CALL XERMShG ('SLATEC', 'BSPVhD',
     +   'NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K', 2, 1)
      RETURN
  210 CONTINUE
      CALL XERMShG ('SLATEC', 'BSPVhD',
     +   'LDVNIK DOES NOT SATISFY LDVNIK.GE.K', 2, 1)
      RETURN
      END

c----------------------------------------------------------------------c

*DECK B2INhT
      SUBROUTINE B2INhT (X, NX, Y, NY, KX, KY, TX, TY, FCN, LDF, WORK,
     +   IFLAG)
C***BEGIN PROLOGUE  B2INhT
C***PURPOSE  B2INhT determines a piecewise polynomial function that
C            interpolates two-dimensional gridded data. Users specify
C            the polynomial order (degree+1) of the interpolant and
C            (optionally) the knot sequence.
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (B2INhT-S, DB2INhT-D)
C***KEYWORDS  INTERPOLATION, TWO DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   B2INhT determines the parameters of a function that interpolates the
C   two-dimensional gridded data (X(i),Y(j),FCN(i,j)) for i=1,..,NX and
C   j=1,..,NY.  The interpolating function and its derivatives may be
C   subsequently evaluated by the function B2VAhL.
C
C   The interpolating function is a piecewise polynomial (pp) function
C   represented as a tensor product of one-dimensional B-splines.  The
C   form of this function is
C
C                          NX   NY
C              S(x,y)  =  SUM  SUM  a   U (x) V (y)
C                         i=1  j=1   ij  i     j
C
C   where the functions U(i) and V(j) are one-dimensional B-spline
C   basis functions. The coefficients a(i,j) are chosen so that
C
C         S(X(i),Y(j)) = FCN(i,j)   for i=1,..,NX and j=1,..,NY.
C
C   Note that for each fixed value of y S(x,y) is a pp function of x
C   alone, and for each fixed x S(x,y) is a pp function of y alone.  In
C   one dimension a piecewise polynomial may be created by partitioning
C   a given interval into subintervals and defining a distinct poly-
C   nomial on each one. The points where adjacent subintervals meet are
C   called  knots.  Each of the functions U(i) and V(j) above is a
C   piecewise polynomial.
C
C   Users of B2INhT choose the order (degree+1) of the polynomial pieces
C   used to define the interpolant in each of the x and y directions
C   (KX and KY). Users also may define their own knot sequence in x and
C   y separately (TX and TY).  If IFLAG=0, however, B2INhT will choose
C   knots that result in a pp interpolant with KX-2 continuous partial
C   derivatives in x and KY-2 continuous partial derivatives in y.  The
C   interpolating function is identically zero outside the rectangular
C   region defined by the knots. See below for more information on knot
C   selection.
C
C   After a call to B2INhT, all information necessary to define the
C   interpolating function is contained in the parameters NX, NY, KX,
C   KY, TX, TY, and FCN.  These quantities should not be altered until
C   after the last call of the evaluation routine B2VAhL.
C
C
C   I N P U T
C   ---------
C
C   X       Real 1D array (size NX)
C           Array of x abscissae. Must be strictly increasing.
C
C   NX      Integer scalar (2 .LE. NX .LE. LDF)
C           Number of x abscissae.
C
C   Y       Real 1D array (size NY)
C           Array of y abscissae. Must be strictly increasing.
C
C   NY      Integer scalar (NY .GE. 2)
C           Number of y abscissae.
C
C   KX      Integer scalar (2 .LE. KX .LE. NX)
C           The order (degree + 1) of polynomial pieces in x.
C
C   KY      Integer scalar (2 .LE. KY .LE. NY)
C           The order (degree + 1) of polynomial pieces in y.
C
C
C   I N P U T   O R   O U T P U T
C   -----------------------------
C
C   TX      Real 1D array (size NX+KX)
C           The knots in the x direction. If IFLAG=0 these are chosen
C           by B2INhT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TX(i).LT.X(i).LT.TX(i+KX),  i=1,..,NX.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NX. (See information below on
C           knot placement.)
C
C   TY      Real 1D array (size NY+KY)
C           The knots in the y direction. If IFLAG=0 these are chosen
C           by B2INhT. If IFLAG=1 they are specified by the user as
C           input. The knots must be non-decreasing and also satisfy
C                  TY(i).LT.Y(i).LT.TY(i+KY),  i=1,..,NY.
C           Equality on the left is permitted for i=1 and equality on
C           the right is permitted for i=NY. (See information below on
C           knot placement.)
C
C
C   I N P U T   A N D   O U T P U T
C   -------------------------------
C
C   FCN     Real 2D array (size LDF by NY)
C           Input : Array of function values to interpolate. FCN(I,J)
C                   should contain the function value at the point
C                   (X(I),Y(J)).
C           Output: Array of coefficients of the B-spline interpolant.
C
C   LDF     Integer scalar (LDF .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C
C   O T H E R
C   ---------
C
C   WORK    Real 1D array
C           Array of working storage. Must be dimensioned of length
C           at least NX*NY + 2*max(KX*(NX+1),KY*(NY+1)).
C
C   IFLAG   Integer scalar.
C           Must be set by the user before B2INhT is called.
C           On return IFLAG indicates the status of the output.
C
C           Input :  0 : knot sequence chosen by user
C                    1 : knot sequence chosen by B2INhT
C
C           Output:  0 : successful execution
C                    2 : IFLAG out of range
C                    3 : NX or LDF out of range
C                    4 : KX out of range
C                    5 : X not strictly increasing
C                    6 : TX is an illegal knot sequence
C                    7 : NY out of range
C                    8 : KY out of range
C                    9 : Y not strictly increasing
C                   10 : TY is an illegal knot sequence
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C           After a successful return, B2INhT may be recalled without
C           resetting IFLAG, provided that only the array FCN has
C           changed.
C
C
C   K N O T   S E L E C T I O N
C   ---------------------------
C
C   In this section we describe the relationship between data points
C   and knots.  Users who choose to let B2INhT select the knot sequences
C   TX and TY (by setting IFLAG=1) may skip this discussion.
C
C   To describe the relationship between knots and data points we  first
C   consider the simpler case of one-dimensional interpolation; in
C   particular, we consider interpolating the data X(i), i=1,..,N with
C   a pp function of order K.
C
C   Knots are the points where the individual polynomial pieces join up,
C   and hence where the pp function may suffer a loss of smoothness.  To
C   define a pp function, one needs N+K knots. If the knots are distinct
C   the interpolant will be as smooth as possible (continuous, with K-2
C   continuous derivatives).  If two adjacent knots come together, the
C   smoothness of the function is reduced at that point.  In general, if
C   M knots are assigned to the same point, the pp function will have
C   K-M-1 continuous derivatives there. If K knots are taken at the same
C   point, then the pp function itself will be discontinuous there.
C
C   Typically, K knots are taken at or to the left of the leftmost data
C   point, K knots at or to the right of the rightmost data point, with
C   the remaining N-K knots in between.  In order for there to be a
C   solution to the interpolation problem the knots T(i), must satisfy
C   certain additional constraints. That is,
C
C                T(i) .LT. X(i) .LT. T(i+K)   i=1,..,N
C
C   Equality is permitted on the left for I=1 and on the right for I=N.
C
C   The two-dimensional interpolant computed by this routine is a tensor
C   product of one-dimensional pp interpolants.  The knots form a grid
C   (TX(i),TY(j)) in the same way that the data points do.  Along lines
C   parallel to the coordinate axes, the interpolant reduces to a one-
C   dimensional pp function with order and knots (KX,TX) or (KY,TY). In
C   this case the appropriate constraints on the knots become
C
C            TX(i) .LT. X(i) .LT. TX(i+KX)   i=1,..,NX
C            TY(i) .LT. Y(i) .LT. TY(i+KY)   i=1,..,NY
C
C   with equality on the left permitted when i=1 and equality on the
C   right permitted when i=NX and i=NY, respectively.
C
C   If these conditions are violated, then B2INhT returns with IFLAG
C   equal to 6 or 10.  The default knot sequence selected by B2INhT
C   always satisfies these conditions.
C
C   When the user sets IFLAG=1 B2INhT selects knots as follows. KX knots
C   are taken at each endpoint in the x direction, not-a-knot end
C   conditions (see references) are used, and the remaining knots are
C   placed at data points if KX is even or at midpoints between data
C   points if KX is odd.  The y direction is treated similarly.  This
C   yields a two-dimensional pp function with KX-2 continuous partial
C   derivatives in x and KY-2 continuous partial derivatives in y.  The
C   interpolant is zero outside the rectangular region defined by the
C   data points, and discontinuous along the boundary of this region.
C
C***SEE ALSO  B2VAhL
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  BKCHK, BKNOT, BTPCF, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  B2INhT
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          NX, NY, KX, KY, LDF, IFLAG
      real             X(NX), Y(NY), TX(NX+KX), TY(NY+KY), FCN(LDF,NY),
     +                 WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'SLATEC', SUBROU = 'B2INhT')
      CHARACTER*50     MESSAG
      INTEGER          I, IW
      LOGICAL          BKCHK
C
      EXTERNAL         BKCHK, BKNOT, BTPCF, XERMShG
C
C***FIRST EXECUTABLE STATEMENT  B2INhT
C
C  -----------------------
C  CHECK VALIDITY OF INPUT
C  -----------------------
C
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 9020
C
      IF ((NX .LT. 2) .OR. (NX .GT. LDF))  GO TO 9030
      IF ((KX .LT. 2) .OR. (KX .GT. NX))  GO TO 9040
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 9050
   10 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. BKCHK(X,NX,KX,TX)) GO TO 9060
      ENDIF
C
      IF (NY .LT. 2)  GO TO 9070
      IF ((KY .LT. 2) .OR. (KY .GT. NY))  GO TO 9080
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 9090
   20 CONTINUE
      IF (IFLAG .EQ. 0) THEN
         IF(.NOT. BKCHK(Y,NY,KY,TY)) GO TO 9100
      ENDIF
C
C  ------------
C  CHOOSE KNOTS
C  ------------
C
      IF (IFLAG .EQ. 1) THEN
         CALL BKNOT(X,NX,KX,TX)
         CALL BKNOT(Y,NY,KY,TY)
      ENDIF
C
C  -------------------------------
C  CONSTRUCT B-SPLINE COEFFICIENTS
C  -------------------------------
C
      IW = NX*NY + 1
      CALL BTPCF(X,NX,FCN,LDF,NY,TX,KX,WORK,NY,WORK(IW))
      CALL BTPCF(Y,NY,WORK,NY,NX,TY,KY,FCN,LDF,WORK(IW))
      IFLAG = 0
      GO TO 9999
C
C  -----------
C  ERROR EXITS
C  -----------
C
 9020 CONTINUE
      IFLAG = 2
      MESSAG = 'IFLAG IS OUT OF RANGE'
      GO TO 9900
C
 9030 CONTINUE
      IFLAG = 3
      MESSAG = 'NX OR LDF IS OUT OF RANGE'
      GO TO 9900
C
 9040 CONTINUE
      IFLAG = 4
      MESSAG = 'KX IS OUT OF RANGE'
      GO TO 9900
C
 9050 CONTINUE
      IFLAG = 5
      MESSAG = 'X ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9060 CONTINUE
      IFLAG = 6
      MESSAG = 'TX IS AN ILLEGAL KNOT SEQUENCE'
      GO TO 9900
C
 9070 CONTINUE
      IFLAG = 7
      MESSAG = 'NY IS OUT OF RANGE'
      GO TO 9900
C
 9080 CONTINUE
      IFLAG = 8
      MESSAG = 'KY IS OUT OF RANGE'
      GO TO 9900
C
 9090 CONTINUE
      IFLAG = 9
      MESSAG = 'Y ARRAY MUST BE STRICTLY INCREASING'
      GO TO 9900
C
 9100 CONTINUE
      IFLAG = 10
      MESSAG = 'TY IS AN ILLEGAL KNOT SEQUENCE'
      GO TO 9900
C
 9900 CONTINUE
      CALL XERMShG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
C
C  ----
C  EXIT
C  ----
C
 9999 CONTINUE
      RETURN
      END

*DECK BKCHK
      LOGICAL FUNCTION BKCHK (X, N, K, T)
C***BEGIN PROLOGUE  BKCHK
C***SUBSIDIARY
C***PURPOSE  BKCHK checks whether piecewise polynomial interpolation
C            of order K at the interpolation points X(i), i=1,..,N
C            is possible using the knot sequence T(i), i=1,..,N+K.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BKCHK-S, DBKCHK-D)
C***KEYWORDS  KNOT SELECTION, INTERPOLATION, PIECEWISE POLYNOMIALS,
C             SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   BKCHK is a subsidiary routine called by B2INhT and B3INT.
C
C   BKCHK checks whether piecewise polynomial interpolation of order K
C   at the interpolation points X(i), i=1,..,N is possible using the
C   knot sequence T(i), i=1,..,N+K.
C
C   Specifically, the following criteria are checked:
C
C     (1) T is non-decreasing
C     (2) T(1) .LE. X(1) .LT. T(1+K)
C     (3) T(i) .LT. X(i) .LT. T(i+K)   i=2,..,N-1
C     (4) T(N) .LT. X(N) .LE. T(N+K)
C
C   BKCHK returns .true. if the knot sequence is legal and .false. if
C   it is not legal.
C
C
C   I N P U T
C   ---------
C
C   X       Real 1D array (size N)
C           Array of data points. Must be strictly increasing.
C
C   N       Integer scalar (N .GE. 2)
C           Number of data points.
C
C   K       Integer scalar (2 .LE. K .LE. N)
C           The order (degree + 1) of polynomial pieces.
C
C   T       Real 1D array (size N+K)
C           The selected knot sequence. Non-decreasing.
C
C
C   CAUTION: The constraints on the variables X, N, and K are not
C            checked by this routine.
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   900910  DATE WRITTEN
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  BKCHK
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, K
      real             X(N), T(N+K)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I
C
C***FIRST EXECUTABLE STATEMENT  BKCHK
      BKCHK = .true.
C
C     ... CHECK WHETHER T IS NON-DECREASING
C
      DO 10 I=2,N+K
         IF (T(I-1) .GT. T(I))  GO TO 100
   10 CONTINUE
C
C     ... CHECK WHETHER THE X ARE PROPERLY DISTRIBUTED
C
      IF ((T(1) .GT. X(1)) .OR. (X(1) .GE. T(1+K)))  GO TO 100
      DO 20 I=2,N-1
         IF ((T(I) .GE. X(I)) .OR. (X(I) .GE. T(I+K)))  GO TO 100
   20 CONTINUE
      IF ((T(N) .GE. X(N)) .OR. (X(N) .GT. T(N+K))) GO TO 100
      GO TO 9999
C
C     ... EXIT THROUGH HERE IF ILLEGAL
C
  100 CONTINUE
      BKCHK = .false.
C
 9999 CONTINUE
      RETURN
      END

*DECK BKNOT
      SUBROUTINE BKNOT (X, N, K, T)
C***BEGIN PROLOGUE  BKNOT
C***SUBSIDIARY
C***PURPOSE  BKNOT selects a sequence of N+K knots for use in spline
C            interpolation of order K at a given set of N data points.
C            The selection yields an interpolant with K-2 continuous
C            derivatives, except at the two endpoints.  Not-a-knot end
C            conditions are used.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BKNOT-S, DBKNOT-D)
C***KEYWORDS  KNOT SELECTION, INTERPOLATION, PIECEWISE POLYNOMIALS,
C             SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   BKNOT is a subsidiary routine called by B2INhT and B3INT.
C
C   BKNOT selects a sequence of N+K knots for use in piecewise
C   polynomial interpolation of order K at a given set of N data points.
C
C   The knots T(i), i=1,..,N+K, are selected as follows:
C
C     number                    location
C     ----------------------------------------------------------------
C     K     at X(1)
C     N-K   between X(1+K/2) and X(N-K+K/2); these are placed at data
C           points if N is even and midway between data points if N is
C           odd
C     K     at X(N)
C
C   where X(i), i=1,..,N are the data points.
C
C   This selection yields an interpolant with K-2 continuous deriva-
C   tives, except at X(1) and X(N); the interpolant will be zero outside
C   (X(1),X(N)).  Knot-a-knot end conditions are used, resulting in K-1
C   continuous derivatives at T(K+1) and T(N).
C
C
C   I N P U T
C   ---------
C
C   X       Real 1D array (size N)
C           Array of data points. Must be strictly increasing.
C
C   N       Integer scalar (N .GE. 2)
C           Number of data points.
C
C   K       Integer scalar (2 .LE. K .LE. N)
C           The order (degree + 1) of polynomial pieces.
C
C
C   O U T P U T
C   -----------
C
C   T       Real 1D array (size N+K)
C           The selected knot sequence. Non-decreasing.
C
C
C   CAUTION: The constraints on the input variables are not checked by
C            this routine.
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  BKNOT
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, K
      real             X(N), T(N+K)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I, IPJ, J
      real             HALF
C
      PARAMETER ( HALF = 0.50E0 )
C
C  ----------------------------
C  PUT K KNOTS AT EACH ENDPOINT
C  ----------------------------
C
C***FIRST EXECUTABLE STATEMENT  BKNOT
      DO 100 J=1,K
         T(J)   = X(1)
         T(N+J) = X(N)
  100 CONTINUE
C
C  --------------------------
C  DISTRIBUTE REMAINING KNOTS
C  --------------------------
C
      IF (MOD(K,2) .EQ. 0) THEN
C
C        ... EVEN K --  KNOTS AT DATA POINTS
C
         I = (K/2) - K
         DO 120 J=K+1,N
            T(J) = X(I+J)
  120    CONTINUE
C
      ELSE
C
C        ... ODD K --  KNOTS BETWEEN DATA POINTS
C
         I = (K-1)/2 - K
         DO 160 J=K+1,N
            IPJ = I + J
            T(J) = HALF*( X(IPJ) + X(IPJ+1) )
  160    CONTINUE
C
      ENDIF
C
      RETURN
      END

*DECK BTPCF
      SUBROUTINE BTPCF (X, N, FCN, LDF, NF, T, K, BCOEF, LDB, WORK)
C***BEGIN PROLOGUE  BTPCF
C***SUBSIDIARY
C***PURPOSE  BTPCF computes NF sets of B-spline coefficients which
C            determine NF distinct piecewise polynomial interpolating
C            functions.  The functions share a common order, set of
C            interpolation points, and knot sequence.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BTPCF-S, DBTPCF-D)
C***KEYWORDS  INTERPOLATION, PIECEWISE POLYNOMIALS, SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   BTPCF is a subsidiary routine called by B2INhT and B3INT.
C
C   BTPCF computes NF sets of B-spline coefficients which determine
C   NF distinct piecewise polynomial interpolating functions of
C   order K.  The functions share a common set of interpolation
C   points X(i), i=1,..,N and knot sequence T(i), i=1,..,N+K.
C
C   The NF sets of N function values to be interpolated are stored in
C   the first NF columns of the array FCN.  On output, the
C   corresponding B-spline coefficients are stored in the ROWS of
C   the array BCOEF.  This is done to facilitate multidimensional
C   interpolation by tensor-products of one-dimensional B-splines.
C   See (de Boor, 1979) for details.
C
C
C   I N P U T
C   ---------
C
C   X       Real 1D array (size N)
C           Array of abscissae. Must be strictly increasing.
C
C   N       Integer scalar (N .GE. 2)
C           Number of abscissae.
C
C   FCN     Real 2D array (size LDF by NF)
C           Array of function values to interpolate. Each column of
C           FCN contains an independent set of N function values
C           corresponding to the array X of abscissae.  That is, the
C           jth data set is (X(i),FCN(i,j)), i=1,..,N.
C
C   LDF     Integer scalar (LDF .GE. N)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C   NF      Integer scalar (NF .GE. 1)
C           Number of data sets (columns of FCN).
C
C   T       Real 1D array (size N+K)
C           The knot sequence.
C
C   K       Integer scalar (2 .LE. K .LE. N)
C           The order (degree + 1) of the piecewise polynomial
C           interpolant.
C
C
C   O U T P U T
C   -----------
C
C   BCOEF   Real 2D array (size NF by N)
C           Array of coefficients of the B-spline interpolants.
C           The coefficients of the B-spline interpolating the
C           data in the jth column of FCN are found in the jth ROW
C           of BCOEF.  That is, the B-spline coefficients defining
C           an interpolant to (X(i),FCN(i,j)), i=1,..,N are
C           stored in BCOEF(j,i), i=1,..,N.
C
C   LDB     Integer scalar (LDB .GE. NF)
C           The actual leading dimension of BCOEF used in the calling
C           program.
C
C
C   M I S C E L L A N E O U S
C   -------------------------
C
C   WORK    Real 1D array (size 2*K*(N+1))
C           Array of working storage.
C
C
C   CAUTION: The constraints on the input variables are not checked by
C            this routine.
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  BINThK, BNSLhV
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  BTPCF
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, LDF, NF, K, LDB
      real             X(N), FCN(LDF,NF), T(N+K), BCOEF(LDB,N), WORK(*)
C
C  .. LOCAL VARIABLES
C
      INTEGER          I, IQ, IW, J, K1, K2
C
      EXTERNAL         BINThK, BNSLhV
C
C  ---------------------------------------------
C  CHECK FOR NULL INPUT AND PARTITION WORK ARRAY
C  ---------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  BTPCF
      IF (NF .LE. 0)  GO TO 500
      K1 = K - 1
      K2 = K1 + K
      IQ = 1 + N
      IW = IQ + K2*N + 1
C
C  -----------------------------
C  COMPUTE B-SPLINE COEFFICIENTS
C  -----------------------------
C
C     ... FIRST DATA SET
C
      CALL BINThK(X,FCN,T,N,K,WORK,WORK(IQ),WORK(IW))
      DO 20 I=1,N
         BCOEF(1,I) = WORK(I)
   20 CONTINUE
C
C     ... ALL REMAINING DATA SETS BY BACK-SUBSTITUTION
C
      IF (NF .GT. 1)  THEN
         DO 100 J=2,NF
            DO 50 I=1,N
               WORK(I) = FCN(I,J)
   50       CONTINUE
            CALL BNSLhV(WORK(IQ),K2,N,K1,K1,WORK)
            DO 60 I=1,N
               BCOEF(J,I) = WORK(I)
   60       CONTINUE
  100    CONTINUE
      ENDIF
C
C  ----
C  EXIT
C  ----
C
  500 CONTINUE
      RETURN
      END

*DECK BINThK
      SUBROUTINE BINThK (X, Y, T, N, K, BCOEF, Q, WORK)
C***BEGIN PROLOGUE  BINThK
C***PURPOSE  Compute the B-representation of a spline which interpolates
C            given data.
C***LIBRARY   SLATEC (BSPLINE)
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (BINThK-S, DBINThK-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C         BINThK is a modification of the SPLINT routine of reference
C         [de Boor 1979].
C
C         BINThK produces the B-spline coefficients, BCOEF, of the
C         B-spline of order K with knots T(I), I=1,...,N+K, which
C         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
C         any of its derivatives can be evaluated by calls to B1VAhL.
C
C         The I-th equation of the linear system A*BCOEF = B for the
C         coefficients of the interpolant enforces interpolation at
C         X(I), I=1,...,N.  Hence, B(I) = Y(I), for all I, and A is
C         a band matrix with 2K-1 bands if A is invertible.  The matrix
C         A is generated row by row and stored, diagonal by diagonal,
C         in the rows of Q, with the main diagonal going into row K.
C         The banded system is then solved by a call to BNFAhC (which
C         constructs the triangular factorization for A and stores it
C         again in Q), followed by a call to BNSLhV (which then
C         obtains the solution BCOEF by substitution).  BNFAhC does no
C         pivoting, since the total positivity of the matrix A makes
C         this unnecessary.  The linear system to be solved is
C         (theoretically) invertible if and only if
C                 T(I) .LT. X(I) .LT. T(I+K),        for all I.
C         Equality is permitted on the left for I=1 and on the right
C         for I=N when K knots are used at X(1) or X(N).  Otherwise,
C         violation of this condition is certain to lead to an error.
C
C     Description of Arguments
C
C         Input
C           X      - vector of length N containing data point abscissae
C                    in strictly increasing order.
C           Y      - corresponding vector of length N containing data
C                    point ordinates.
C           T      - knot vector of length N+K
C                    Since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
C                    .GE. X(N), this leaves only N-K knots (not nec-
C                    essarily X(I) values) interior to (X(1),X(N))
C           N      - number of data points, N .GE. K
C           K      - order of the spline, K .GE. 1
C
C         Output
C           BCOEF  - a vector of length N containing the B-spline
C                    coefficients
C           Q      - a work vector of length (2*K-1)*N, containing
C                    the triangular factorization of the coefficient
C                    matrix of the linear system being solved.  The
C                    coefficients for the interpolant of another data
C                    set (X(I),YY(I)), I=1,...,N, with the same
C                    abscissae can be obtained by loading YY into
C                    BCOEF and then executing
C                        CALL BNSLhV(Q,2K-1,N,K-1,K-1,BCOEF)
C           WORK   - work vector of length 2*K
C
C     Error Conditions
C         Improper input is a fatal error, as is a singular system of
C         equations.  Conditions checked and their XERMShG error numbers:
C               1   K is less than 1
C               2   N is less than K
C               3   X values are not distinct or not ordered
C               4   Some abscissa was not in the support of the corres-
C                   ponding basis function, and the system is singular.
C               8   The system solver has detected a singular system,
C                   although the theoretical conditions for a solution
C                   were satisfied.
C         Note: These numbers are printed by XERMShG, but not returned
C         to the calling program.
C
C***REFERENCES  D. E. Amos, Computation with splines and B-splines,
C                 Report SAND78-1968, Sandia Laboratories, March 1979.
C               Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  BNFAhC, BNSLhV, BSPVhN, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900910  Changed BVALU references to B1VAhL in DESCRIPTION.  (RFB)
C   901114  Miscellaneous housekeeping changes.  (FNF)
C           1. Improved readability of prologue and reduced single/
C              double differences.
C           2. Revised declarations.
C           3. Changed XERMShG calls so that different errors have
C              different error numbers, and recorded the numbers in
C              the prologue.
C   920501  Reformatted the REFERENCES section.  (WRB)
C   930527  Modified error return coding to look more like BINKP.  (FNF)
C***END PROLOGUE  BINThK
C
C  Implementation note:  The error numbers occur in the order of
C     detection.  The number 8 is historical.  Prior to 901114, the
C     rest all printed number 2.  The errors at statements 80 and 90
C     were detected by de Boor's original SPLINT and the rest were
C     added later by Amos.  (FNF)
C
C  Additional note: This routine is now ready to have IFLAG added to
C     its argument list, as in BINKP.
C**End
C
      INTEGER  N, K
      real  X(*), Y(*), T(*), BCOEF(*), Q(*), WORK(*)
C     DIMENSION Q(2*K-1,N), T(N+K)
C
      EXTERNAL  BNFAhC, BNSLhV, BSPVhN, XERMShG
C
      INTEGER  I, IFLAG, ILP1MX, IWORK, J, JJ, KM1, KPKM2, LEFT, LENQ,
     +   NP1
      real  XI
C
C***FIRST EXECUTABLE STATEMENT  BINThK
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      JJ = N - 1
      IF(JJ.EQ.0) GO TO 6
      DO 5 I=1,JJ
      IF(X(I).GE.X(I+1)) GO TO 110
    5 CONTINUE
    6 CONTINUE
      NP1 = N + 1
      KM1 = K - 1
      KPKM2 = 2*KM1
      LEFT = K
C                ZERO OUT ALL ENTRIES OF Q
      LENQ = N*(K+KM1)
      DO 10 I=1,LENQ
        Q(I) = 0
   10 CONTINUE
C
C  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO 50 I=1,N
        XI = X(I)
        ILP1MX = MIN(I+K,NP1)
C        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
C                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
C        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        LEFT = MAX(LEFT,I)
        IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
        LEFT = LEFT + 1
        IF (LEFT.LT.ILP1MX) GO TO 20
        LEFT = LEFT - 1
        IF (XI.GT.T(LEFT+1)) GO TO 80
C        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
C        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
C        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
C        ARE RETURNED, IN  BCOEF (USED FOR TEMP. STORAGE HERE), BY THE
C        FOLLOWING
   30   CALL BSPVhN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
C        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
C        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
C        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
C        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
C        BNFAhC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
C        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
C        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
C        ENTRY
C            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
C                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
C        OF  Q .
        JJ = I - LEFT + 1 + (LEFT-K)*(K+KM1)
        DO 40 J=1,K
          JJ = JJ + KPKM2
          Q(JJ) = BCOEF(J)
   40   CONTINUE
   50 CONTINUE
C
C     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
      CALL BNFAhC(Q, K+KM1, N, KM1, KM1, IFLAG)
      GO TO (60, 90), IFLAG
C     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60 DO 70 I=1,N
        BCOEF(I) = Y(I)
   70 CONTINUE
      CALL BNSLhV(Q, K+KM1, N, KM1, KM1, BCOEF)
      IFLAG = 0
      RETURN
C
C     ERROR RETURNS
C
   80 CONTINUE
      IFLAG = 4
      CALL XERMShG ('SLATEC', 'BINThK',
     +   'SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE CORRESPONDING ' //
     +   'BASIS FUNCTION AND THE SYSTEM IS SINGULAR.', IFLAG, 1)
      RETURN
   90 CONTINUE
      IFLAG = 8
      CALL XERMShG ('SLATEC', 'BINThK',
     +   'THE SYSTEM OF SOLVER DETECTS A SINGULAR SYSTEM, ALTHOUGH ' //
     +   'THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATISFIED.',
     +   IFLAG, 1)
      RETURN
  100 CONTINUE
      IFLAG = 1
      CALL XERMShG ('SLATEC', 'BINThK',
     +   'K DOES NOT SATISFY K.GE.1', IFLAG, 1)
      RETURN
  105 CONTINUE
      IFLAG = 2
      CALL XERMShG ('SLATEC', 'BINThK',
     +   'N DOES NOT SATISFY N.GE.K', IFLAG, 1)
      RETURN
  110 CONTINUE
      IFLAG = 3
      CALL XERMShG ('SLATEC', 'BINThK',
     +   'X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SOME I', IFLAG, 1)
      RETURN
      END

*DECK BNFAhC
      SUBROUTINE BNFAhC (W, NROWW, NROW, NBANDL, NBANDU, IFLAG)
C***BEGIN PROLOGUE  BNFAhC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BINTh4 and BINThK
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BNFAhC-S, DBNFAhC-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  BNFAhC is the BANFAC routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  Returns in  W  the lu-factorization (without pivoting) of the banded
C  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
C  onals in the work array  W .
C
C *****  I N P U T  ******
C  W.....Work array of size  (NROWW,NROW)  containing the interesting
C        part of a banded matrix  A , with the diagonals or bands of  A
C        stored in the rows of  W , while columns of  A  correspond to
C        columns of  W . This is the storage mode used in  LINPACK  and
C        results in efficient innermost loops.
C           Explicitly,  A  has  NBANDL  bands below the diagonal
C                            +     1     (main) diagonal
C                            +   NBANDU  bands above the diagonal
C        and thus, with    MIDDLE = NBANDU + 1,
C          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
C                                              J=1,...,NROW .
C        For example, the interesting entries of A (1,2)-banded matrix
C        of order  9  would appear in the first  1+1+2 = 4  rows of  W
C        as follows.
C                          13 24 35 46 57 68 79
C                       12 23 34 45 56 67 78 89
C                    11 22 33 44 55 66 77 88 99
C                    21 32 43 54 65 76 87 98
C
C        All other entries of  W  not identified in this way with an en-
C        try of  A  are never referenced .
C  NROWW.....Row dimension of the work array  W .
C        must be  .GE.  NBANDL + 1 + NBANDU  .
C  NBANDL.....Number of bands of  A  below the main diagonal
C  NBANDU.....Number of bands of  A  above the main diagonal .
C
C *****  O U T P U T  ******
C  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
C     If  IFLAG = 1, then
C  W.....contains the LU-factorization of  A  into a unit lower triangu-
C        lar matrix  L  and an upper triangular matrix  U (both banded)
C        and stored in customary fashion over the corresponding entries
C        of  A . This makes it possible to solve any particular linear
C        system  A*X = B  for  X  by A
C              CALL BNSLhV ( W, NROWW, NROW, NBANDL, NBANDU, B )
C        with the solution X  contained in  B  on return .
C     If  IFLAG = 2, then
C        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
C        one of the potential pivots was found to be zero indicating
C        that  A  does not have an LU-factorization. This implies that
C        A  is singular in case it is totally positive .
C
C *****  M E T H O D  ******
C     Gauss elimination  W I T H O U T  pivoting is used. The routine is
C  intended for use with matrices  A  which do not require row inter-
C  changes during factorization, especially for the  T O T A L L Y
C  P O S I T I V E  matrices which occur in spline calculations.
C     The routine should not be used for an arbitrary banded matrix.
C
C***SEE ALSO  BINTh4, BINThK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  BNFAhC
C
      INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K,
     1 KMAX, MIDDLE, MIDMK, NROWM1
      real W(NROWW,*), FACTOR, PIVOT
C
C***FIRST EXECUTABLE STATEMENT  BNFAhC
      IFLAG = 1
      MIDDLE = NBANDU + 1
C                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
      NROWM1 = NROW - 1
      IF (NROWM1) 120, 110, 10
   10 IF (NBANDL.GT.0) GO TO 30
C                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
      DO 20 I=1,NROWM1
        IF (W(MIDDLE,I).EQ.0.0E0) GO TO 120
   20 CONTINUE
      GO TO 110
   30 IF (NBANDU.GT.0) GO TO 60
C              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
C                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
      DO 50 I=1,NROWM1
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0E0) GO TO 120
        JMAX = MIN(NBANDL,NROW-I)
        DO 40 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   40   CONTINUE
   50 CONTINUE
      RETURN
C
C        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
   60 DO 100 I=1,NROWM1
C                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0E0) GO TO 120
C                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
C                     BELOW THE DIAGONAL .
        JMAX = MIN(NBANDL,NROW-I)
C              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
        DO 70 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   70   CONTINUE
C                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
C                     THE RIGHT OF THE DIAGONAL .
        KMAX = MIN(NBANDU,NROW-I)
C                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
C                  (BELOW ROW  I ) .
        DO 90 K=1,KMAX
          IPK = I + K
          MIDMK = MIDDLE - K
          FACTOR = W(MIDMK,IPK)
          DO 80 J=1,JMAX
            W(MIDMK+J,IPK) = W(MIDMK+J,IPK) - W(MIDDLE+J,I)*FACTOR
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
C                                       CHECK THE LAST DIAGONAL ENTRY .
  110 IF (W(MIDDLE,NROW).NE.0.0E0) RETURN
  120 IFLAG = 2
      RETURN
      END

*DECK BNSLhV
      SUBROUTINE BNSLhV (W, NROWW, NROW, NBANDL, NBANDU, B)
C***BEGIN PROLOGUE  BNSLhV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BINTh4 and BINThK
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BNSLhV-S, DBNSLhV-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  BNSLhV is the BANSLV routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  Companion routine to  BNFAhC . It returns the solution  X  of the
C  linear system  A*X = B  in place of  B , given the LU-factorization
C  for  A  in the work array  W from BNFAhC.
C
C *****  I N P U T  ******
C  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a
C        banded matrix  A  of order  NROW  as constructed in  BNFAhC .
C        For details, see  BNFAhC .
C  B.....Right side of the system to be solved .
C
C *****  O U T P U T  ******
C  B.....Contains the solution  X , of order  NROW .
C
C *****  M E T H O D  ******
C     (With  A = L*U, as stored in  W,) the unit lower triangular system
C  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
C  upper triangular system  U*X = Y  is solved for  X  . The calcul-
C  ations are so arranged that the innermost loops stay within columns.
C
C***SEE ALSO  BINTh4, BINThK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  BNSLhV
C
      INTEGER NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
      real W(NROWW,*), B(*)
C***FIRST EXECUTABLE STATEMENT  BNSLhV
      MIDDLE = NBANDU + 1
      IF (NROW.EQ.1) GO TO 80
      NROWM1 = NROW - 1
      IF (NBANDL.EQ.0) GO TO 30
C                                 FORWARD PASS
C            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
      DO 20 I=1,NROWM1
        JMAX = MIN(NBANDL,NROW-I)
        DO 10 J=1,JMAX
          B(I+J) = B(I+J) - B(I)*W(MIDDLE+J,I)
   10   CONTINUE
   20 CONTINUE
C                                 BACKWARD PASS
C            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
C            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW).
   30 IF (NBANDU.GT.0) GO TO 50
C                                A  IS LOWER TRIANGULAR .
      DO 40 I=1,NROW
        B(I) = B(I)/W(1,I)
   40 CONTINUE
      RETURN
   50 I = NROW
   60 B(I) = B(I)/W(MIDDLE,I)
      JMAX = MIN(NBANDU,I-1)
      DO 70 J=1,JMAX
        B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
   70 CONTINUE
      I = I - 1
      IF (I.GT.1) GO TO 60
   80 B(1) = B(1)/W(MIDDLE,1)
      RETURN
      END

*DECK BSPVhN
      SUBROUTINE BSPVhN (T, JHIGH, K, INDEX, X, ILEFT, VNIKX, WORK,
     +   IWORK)
C***BEGIN PROLOGUE  BSPVhN
C***PURPOSE  Calculate the value of all (possibly) nonzero basis
C            functions at X.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (BSPVhN-S, DBSPVhN-D)
C***KEYWORDS  EVALUATION OF B-SPLINE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C         BSPVhN is the BSPLVN routine of the reference.
C
C         BSPVhN calculates the value of all (possibly) nonzero basis
C         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where
C         T(K) .LE. X .LE. T(N+1) and J=IWORK is set inside the routine
C         on the first call when INDEX=1.  ILEFT is such that T(ILEFT)
C         .LE. X .LT. T(ILEFT+1).  A call to INTRhV(T,N+1,X,ILO,ILEFT,
C         MFLAG) produces the proper ILEFT.  BSPVhN calculates using the
C         basic algorithm needed in BSPVhD.  If only basis functions are
C         desired, setting JHIGH=K and INDEX=1 can be faster than
C         calling BSPVhD, but extra coding is required for derivatives
C         (INDEX=2) and BSPVhD is set up for this purpose.
C
C         Left limiting values are set up as described in BSPVhD.
C
C     Description of Arguments
C         Input
C          T       - knot vector of length N+K, where
C                    N = number of B-spline basis functions
C                    N = sum of knot multiplicities-K
C          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
C          K       - highest possible order
C          INDEX   - INDEX = 1 gives basis functions of order JHIGH
C                          = 2 denotes previous entry with WORK, IWORK
C                              values saved for subsequent calls to
C                              BSPVhN.
C          X       - argument of basis functions,
C                    T(K) .LE. X .LE. T(N+1)
C          ILEFT   - largest integer such that
C                    T(ILEFT) .LE. X .LT. T(ILEFT+1)
C
C         Output
C          VNIKX   - vector of length K for spline values.
C          WORK    - a work vector of length 2*K
C          IWORK   - a work parameter.  Both WORK and IWORK contain
C                    information necessary to continue for INDEX = 2.
C                    When INDEX = 1 exclusively, these are scratch
C                    variables and can be used for other purposes.
C
C     Error Conditions
C         Improper input is a fatal error.
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  XERMShG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BSPVhN
C
      INTEGER ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
      real T, VM, VMPREV, VNIKX, WORK, X
C     DIMENSION T(ILEFT+JHIGH)
      DIMENSION T(*), VNIKX(*), WORK(*)
C     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
C     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
C***FIRST EXECUTABLE STATEMENT  BSPVhN
      IF(K.LT.1) GO TO 90
      IF(JHIGH.GT.K .OR. JHIGH.LT.1) GO TO 100
      IF(INDEX.LT.1 .OR. INDEX.GT.2) GO TO 105
      IF(X.LT.T(ILEFT) .OR. X.GT.T(ILEFT+1)) GO TO 110
      GO TO (10, 20), INDEX
   10 IWORK = 1
      VNIKX(1) = 1.0E0
      IF (IWORK.GE.JHIGH) GO TO 40
C
   20 IPJ = ILEFT + IWORK
      WORK(IWORK) = T(IPJ) - X
      IMJP1 = ILEFT - IWORK + 1
      WORK(K+IWORK) = X - T(IMJP1)
      VMPREV = 0.0E0
      JP1 = IWORK + 1
      DO 30 L=1,IWORK
        JP1ML = JP1 - L
        VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
        VNIKX(L) = VM*WORK(L) + VMPREV
        VMPREV = VM*WORK(K+JP1ML)
   30 CONTINUE
      VNIKX(JP1) = VMPREV
      IWORK = JP1
      IF (IWORK.LT.JHIGH) GO TO 20
C
   40 RETURN
C
C
   90 CONTINUE
      CALL XERMShG ('SLATEC', 'BSPVhN', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  100 CONTINUE
      CALL XERMShG ('SLATEC', 'BSPVhN',
     +   'JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K', 2, 1)
      RETURN
  105 CONTINUE
      CALL XERMShG ('SLATEC', 'BSPVhN', 'INDEX IS NOT 1 OR 2', 2, 1)
      RETURN
  110 CONTINUE
      CALL XERMShG ('SLATEC', 'BSPVhN',
     +   'X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT+1)', 2, 1)
      RETURN
      END

*DECK INTRhV
      SUBROUTINE INTRhV (XT, LXT, X, ILO, ILEFT, MFLAG)
C***BEGIN PROLOGUE  INTRhV
C***PURPOSE  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
C            such that XT(ILEFT) .LE. X where XT(*) is a subdivision
C            of the X interval.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (INTRhV-S, DINTRhV-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C         INTRhV is the INTERV routine of the reference.
C
C         INTRhV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
C         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
C         the X interval.  Precisely,
C
C                      X .LT. XT(1)                1         -1
C         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
C           XT(LXT) .LE. X                         LXT        1,
C
C         That is, when multiplicities are present in the break point
C         to the left of X, the largest index is taken for ILEFT.
C
C     Description of Arguments
C         Input
C          XT      - XT is a knot or break point vector of length LXT
C          LXT     - length of the XT vector
C          X       - argument
C          ILO     - an initialization parameter which must be set
C                    to 1 the first time the spline array XT is
C                    processed by INTRhV.
C
C         Output
C          ILO     - ILO contains information for efficient process-
C                    ing after the initial call, and ILO must not be
C                    changed by the user.  Distinct splines require
C                    distinct ILO parameters.
C          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
C          MFLAG   - signals when X lies out of bounds
C
C     Error Conditions
C         None
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  INTRhV
C
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      real X, XT
      DIMENSION XT(*)
C***FIRST EXECUTABLE STATEMENT  INTRhV
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
C
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
C
C *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
C *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
C
C *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
C *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END

*DECK B2VAhL
      real FUNCTION B2VAhL (XVAL, YVAL, IDX, IDY, TX, TY, NX, NY,KX,KY,
     *   FCN, LDF, WORK, IFLAG)
C***BEGIN PROLOGUE  B2VAhL
C***PURPOSE  B2VAhL evaluates the two-dimensional piecewise polynomial
C            interpolating function constructed by the routine B2INhT.
C            Either function values or partial derivative values may be
C            be requested.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (B2VAhL-S, DB2VAhL-D)
C***KEYWORDS  TWO DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS, EVALUATION, DIFFERENTIATION
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   B2VAhL evaluates the tensor product piecewise polynomial interpolant
C   constructed by the routine B2INhT or one of its derivatives at the
C   point (XVAL,YVAL).  The variables IDX and IDY indicate which
C   function is to be evaluated.  If B(x,y) is the interpolant,
C
C                            IDX+IDY
C                           d
C              B2VAhL  =   ----------- B (XVAL,YVAL)
C                           IDX   IDY
C                         dx    dy
C
C   Thus, to evaluate the interpolant itself, set IDX=IDY=0. To get the
C   first partial derivative with respect to x, set IDX=1 and IDY=0, and
C   so on.
C
C   Since B is a piecewise polynomial of degree KX-1 in x and KY-1 in y,
C   B2VAhL returns zero whenever IDX.GE.KX or IDY.GE.KY.  B2VAhL also
C   returns zero if (XVAL,YVAL) is out of range, that is, if
C
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY).
C
C   If the knots TX and TY were chosen by B2INhT, then this is
C   equivalent to
C
C              XVAL.LT.X(1) .OR. XVAL.GT.X(NX) .OR.
C              YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY).
C
C   The input quantities TX, TY, NX, NY, KX, KY, and FCN should be
C   unchanged since the call to B2INhT which produced FCN.
C
C   Note that the derivative functions computed by B2VAhL may be
C   discontinuous when x or y correspond to knots.  (When B2INhT selects
C   knots this occurs only when IDX=KX-1 or IDY=KY-1.)  In these cases
C   B2VAhL returns right limiting values (right derivatives), except at
C   the rightmost knot where left limiting values are returned.
C
C   I N P U T
C   ---------
C
C   XVAL    Real scalar
C           X coordinate of evaluation point.
C
C   YVAL    Real scalar
C           Y coordinate of evaluation point.
C
C   IDX     Integer scalar  (IDX .GE. 0)
C           Indicates the x derivative of piecewise polynomial to
C           evaluate: IDX=J for the Jth partial derivative with respect
C           to x. B2VAhL will return 0 if IDX.GE.KX.
C
C   IDY     Integer scalar  (IDY .GE. 0)
C           Indicates the y derivative of piecewise polynomial to
C           evaluate: IDY=J for the Jth partial derivative with respect
C           to y. B2VAhL will return 0 if IDY.GE.KY.
C
C   TX      Real 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Unchanged since the call to B2INhT which
C           produced FCN.)
C
C   TY      Real 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Unchanged since the call to B2INhT which
C           produced FCN.)
C
C   NX      Integer scalar (NX .GE. KX)
C           The number of interpolation points in x. (Unchanged since
C           the call to B2INhT which produced FCN.)
C
C   NY      Integer scalar (NY .GE. KY)
C           The number of interpolation points in y. (Unchanged since
C           the call to B2INhT which produced FCN.)
C
C   KX      Integer scalar (KX .GE. 2)
C           Order of polynomial pieces in x. (Unchanged since the call
C           to B2INhT which produced FCN.)
C
C   KY      Integer scalar (KY .GE. 2)
C           Order of polynomial pieces in y. (Unchanged since the call
C           to B2INhT which produced FCN.)
C
C   FCN     Real 2D array (size LDF by NY)
C           The B-spline coefficients computed by B2INhT.
C
C   LDF     Integer scalar (LDF .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C
C   O T H E R
C   ---------
C
C   WORK    Real 1D array (size 3*max(KX,KY) + KY)
C           A working storage array.
C
C   IFLAG   Integer scalar.
C           On return IFLAG indicates the status of the output.
C
C           Output:  0 : successful execution
C                    1 : KX out of range
C                    2 : NX or LDF out of range
C                    3 : KY out of range
C                    4 : NY out of range
C                    5 : IDX or IDY out of range
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C
C   B A C K G R O U N D
C   -------------------
C
C   B2VAhL evaluates the following function or its partial derivatives
C
C                        NX   NY
C            B(x,y)  =  SUM  SUM  FCN(i,j) U (x) V (y)
C                       i=1  j=1            i     j
C
C   where
C
C     U (x)  is the ith (one-dimensional) B-spline basis function
C      i     defined by NX, KX, and TX, and
C
C     V (y)  is the jth (one-dimensional) B-spline basis function
C      j     defined by NY, KY, and TY.
C
C   See (de Boor, 1978) for a description of the B-spline basis.
C
C
C***SEE ALSO  B2INhT
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  B1VAhL, INTRhV, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930510  Eliminated multiple calls to XERMShG.  (FNF)
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  B2VAhL
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          IDX, IDY, NX, NY, KX, KY, LDF, IFLAG
      real             XVAL, YVAL, TX(NX+KX), TY(NY+KY), FCN(LDF,NY),
     +                 WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'SLATEC', SUBROU = 'B2VAhL')
      CHARACTER*50     MESSAG
      INTEGER          I, IERR, JMIN, JMAX, INX, INY, IPOS, IW, LEFTY,
     +                 MFLAG
      real             B1VAhL
C
      EXTERNAL         INTRhV, B1VAhL, XERMShG

      Use(Timespl)
      real tick, tock, tsint
      external tick, tock

C
C***FIRST EXECUTABLE STATEMENT  B2VAhL
      B2VAhL = 0
      IFLAG = 0
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IF (KX .LT. 1)  GO TO 9001
      IF ((NX .LT. KX) .OR. (NX .GT. LDF)) GO TO 9002
      IF (KY .LT. 1)  GO TO 9003
      IF (NY .LT. KY)  GO TO 9004
      IF ((IDX .LT. 0) .OR. (IDY .LT. 0)) GO TO 9005
C
C     ... QUICK RETURN WHEN SPLINE EVALUATES TO ZERO
C
      IF ((IDX .GE. KX) .OR. (IDY .GE. KY)) GO TO 9999
      IF ((XVAL .LT. TX(1)) .OR. (XVAL .GT. TX(NX+KX))) GO TO 9999
      IF ((YVAL .LT. TY(1)) .OR. (YVAL .GT. TY(NY+KY))) GO TO 9999
C
C  ---------------------
C  EVALUATE THE B-SPLINE
C  ---------------------
C
C  The summation in the Background section above can be rewritten as
C
C                                NY
C  (1)               B(x,y)  =  SUM  b  V (y),
C                               j=1   j  j
C  where
C                         NX
C  (2)             b  =  SUM  FCN(i,j) U (x),   j = 1,..,NY.
C                   j    i=1            i
C
C  Note that each summation is the evaluation of a one-dimensional
C  B-spline, which can be done by calls to B1VAhL.  At most KY basis
C  functions in (1) are nonzero.  The indices of these functions are
C  determined, and only the coefficients b-sub-j which multiply nonzero
C  basis functions are computed using (2).
C
C
C     ... FIND KNOT INTERVAL CONTAINING YVAL
C
      INY = 1
      tsint = tick()
      CALL INTRhV(TY,NY+KY,YVAL,INY,LEFTY,MFLAG)
      totintrv = totintrv + tock(tsint)
      IF (MFLAG .NE. 0) THEN
C
C        ... YVAL .EQ. T(NY+KY),  ADJUST TO GET LEFT LIMITING VALUE
C
   10    CONTINUE
         LEFTY = LEFTY - 1
         IF (YVAL .EQ. TY(LEFTY)) GO TO 10
      ENDIF
C
C     ... FIND RANGE OF INDICIES OF NONZERO BASIS FUNCTIONS IN (1)
C         (JMIN,JMAX) = (SMALLEST,LARGEST)
C
      IF (LEFTY .LT. KY) THEN
         JMIN = 1
         JMAX = KY
      ELSEIF (LEFTY .GT. NY) THEN
         JMIN = NY - KY + 1
         JMAX = NY
      ELSE
         JMIN = LEFTY-KY+1
         JMAX = LEFTY
      ENDIF
C
C     ... GET COEFFICIENTS OF NONZERO BASIS FUNCTIONS IN (1) USING (2)
C
      IPOS = 0
      INX = 1
      IW = KY + 1
      DO 30 I=JMIN,JMAX
         IPOS = IPOS + 1
         WORK(IPOS) = B1VAhL(XVAL,IDX,TX,NX,KX,FCN(1,I),INX,WORK(IW),
     +                      IERR)
   30 CONTINUE
C
C     ... EVALUATE INTERPOLANT USING (1)
C
      INY = KY - 1
      B2VAhL = B1VAhL(YVAL,IDY,TY(JMIN),KY,KY,WORK,INY,WORK(IW),IERR)
      GO TO 9999
C
C  -----
C  EXITS
C  -----
C
 9001 CONTINUE
      IFLAG = 1
      MESSAG = 'KX IS OUT OF RANGE'
      GO TO 9900
C
 9002 CONTINUE
      IFLAG = 2
      MESSAG = 'NX OR LDF IS OUT OF RANGE'
      GO TO 9900
C
 9003 CONTINUE
      IFLAG = 3
      MESSAG = 'KY IS OUT OF RANGE'
      GO TO 9900
C
 9004 CONTINUE
      IFLAG = 4
      MESSAG = 'NY IS OUT OF RANGE'
      GO TO 9900
C
 9005 CONTINUE
      IFLAG = 5
      MESSAG = 'IDX OR IDY IS OUT OF RANGE'
      GO TO 9900
C
 9900 CONTINUE
      CALL XERMShG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
C
C     ... EXIT
C
 9999 CONTINUE
      RETURN
      END

*DECK B1VAhL
      real FUNCTION B1VAhL (X, IDERIV, T, N, K, A, INBV, WORK, IFLAG)
C***BEGIN PROLOGUE  B1VAhL
C***PURPOSE  Evaluates the B-representation of a spline at X for the
C            function value or any of its derivatives.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (B1VAhL-S, DB1VAhL-D)
C***KEYWORDS  B-SPLINE, DIFFERENTIATION, EVALUATION, SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   B1VAhL evaluates the B-representation (T,N,K,A) of a spline or one
C   of its derivatives at X.  Variable IDERIV indicates which function
C   is to be evaluated.  If B(x) is the spline,
C
C                              IDERIV
C                             d
C                  B1VAhL  =   -------- B(X)
C                               IDERIV
C                             dx
C
C   Thus, to evaluate the spline itself, set IDERIV=0, to get the first
C   derivative with respect to x, set IDERIV=1, and so on.
C
C   Since B is a piecewise polynomial of degree K-1 B1VAhL returns zero
C   whenever IDERIV.GE.K.  B1VAhL also returns zero if X is out of
C   range, that is, if
C
C                    X.LT.T(1) .OR. X.GT.T(N+K).
C
C   Note that the derivative functions computed by B1VAhL may be
C   discontinuous when X corresponds to a knot.   In these cases B1VAhL
C   returns right limiting values (right derivatives), except at T(N+K)
C   where left limiting values are returned.
C
C
C   I N P U T
C   ---------
C
C   X       Real scalar
C           The point at which the B-spline is to be evaluated.
C
C   IDERIV  Integer scalar  (IDERIV .GE. 0)
C           Indicates the derivative of the B-spline to evaluate:
C           IDERIV=0 for the function itself, IDERIV=J for the Jth
C           derivative. B1VAhL will return 0 if IDERIV.GE.K.
C
C   T       Real 1D array (size N+K)
C           Sequence of knots defining the B-spline.
C
C   N       Integer scalar
C           The number of B-spline coefficients. (N = sum of knot
C           multiplicities - K.)
C
C   K       Integer scalar
C           Order of the B-spline.
C
C   A       Real 1D array (size N)
C           The B-spline coefficients.
C
C
C   I N P U T   A N D   O U T P U T
C   -------------------------------
C
C   INBV    Integer scalar
C           An initialization parameter which must be set to 1 the first
C           time B1VAhL is called. On output, contains information for
C           efficient processing after the initial call and must not be
C           changed by the user.  Distinct splines require distinct INBV
C           parameters.
C
C
C   O T H E R
C   ---------
C
C   WORK    Real 1D array (size 3*K)
C           Workspace.
C
C   IFLAG   Integer scalar.
C           On return IFLAG indicates the status of the output.
C
C           Output:  0 : successful execution
C                    1 : K out of range
C                    2 : N out of range
C                    3 : IDERIV out of range
C                    4 : Unable to find left-limiting value at T(N+K)
C
C
C  B1VAhL is a version of the routine BVALUE written by Carl de Boor and
C  given in the reference.
C
C  This routine replaces the former SLATEC routine BVALU.  It differs
C  from BVALU in that
C
C    (1) the order of arguments in the calling sequence has been
C        changed and a new parameter, IFLAG, has been added,
C    (2) X may be any real number,
C    (2) IDERIV may be any non-negative integer,
C    (3) right-limiting values are returned when X=T(N+1) and
C        T(N+K).GT.T(N+1), and
C    (4) the prescription given in BVALU for obtaining left-limiting
C        values at interior knots does not work.
C
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  INTRhV, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   900726  DATE WRITTEN
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C***END PROLOGUE  B1VAhL
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          N, K, IDERIV, INBV, IFLAG
      real             T(N+K), A(N), X, WORK(3*K)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'SLATEC', SUBROU = 'B1VAhL')
      CHARACTER*50     MESSAG
      INTEGER          I, IDL, IDR, ILO, IMK, IP1, J, JJ, JMIN, JMAX,
     +                 KM1, KMJ, MFLAG, NPK
      real             FKMJ
C
      EXTERNAL         INTRhV, XERMShG

      Use(Timespl)
      real tick, tock, tsint
      external tick, tock

C
C***FIRST EXECUTABLE STATEMENT  B1VAhL
      B1VAhL = 0
      NPK = N + K
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IF (K .LT. 1)  GO TO 9001
      IF (N .LT. K)  GO TO 9002
      IF (IDERIV .LT. 0)  GO TO 9003
C
C     ... QUICK RETURN WHEN SPLINE EVALUATES TO ZERO
C
      IF (IDERIV .GE. K)  GO TO 9999
      IF ((X .LT. T(1)) .OR. (X .GT. T(NPK)))  GO TO 9999
C
C  -------------------------------
C  FIND KNOT INTERVAL CONTAINING X
C  -------------------------------
C
C     ... FIND LARGEST I IN (1,N+K) SUCH THAT T(I) .LE. X .LT. T(I+1)
C
      tsint = tick()
      CALL INTRhV(T, NPK, X, INBV, I, MFLAG)
      totintrv = totintrv + tock(tsint)
      IF (MFLAG .NE. 0) THEN
C
C        ... X .EQ. T(N+K),  USE LEFT LIMITING VALUE
C
   10    CONTINUE
         IF (I .EQ. 1)  GO TO 9004
         I = I - 1
         IF (X .EQ. T(I)) GO TO 10
      ENDIF
C
C  ------------------------------------------------------
C  STORE COEFFICIENTS OF NON-ZERO BASIS FUNCTIONS IN WORK
C  ------------------------------------------------------
C
C     ... IF THERE ARE LESS THAN K BASIS FUNCTIONS, WE ASSUME FICTITIOUS
C         ONES WITH ZERO COEFFICIENTS AND KNOTS AT T(1) OR T(N+K)
C
      IMK = I - K
      JMIN = MAX(1-IMK,1)
      JMAX = MIN(NPK-I,K)
      DO 20 J=1,JMIN-1
         WORK(J) = 0
   20 CONTINUE
      DO 30 J=JMIN,JMAX
         WORK(J) = A(IMK+J)
   30 CONTINUE
      DO 40 J=JMAX+1,K
         WORK(J) = 0
   40 CONTINUE
C
C  ----------------------------
C  COMPUTE AUXILIARY QUANTITIES
C  ----------------------------
C
C     ... DL(J) = X - T(I+1-J) = WORK(IDL+J), J=1,..,K-1
C     ... DR(J) = T(I+J) - X   = WORK(IDR+J), J=1,..,K-1
C
      IDL = K
      IDR = IDL + K
      KM1 = K - 1
      IP1 = I + 1
      JMAX = MIN(I,KM1)
      DO 50 J=1,JMAX
         WORK(IDL+J) = X - T(IP1-J)
   50 CONTINUE
      DO 60 J=JMAX+1,KM1
         WORK(IDL+J) = WORK(IDL+JMAX)
   60 CONTINUE
      JMAX = MIN(NPK-I,KM1)
      DO 70 J=1,JMAX
         WORK(IDR+J) = T(I+J) - X
   70 CONTINUE
      DO 80 J=JMAX+1,KM1
         WORK(IDR+J) = WORK(IDR+JMAX)
   80 CONTINUE
C
C  ----------------------------------------
C  DIFFERENCE THE COEFFICIENTS IDERIV TIMES
C  ----------------------------------------
C
      DO 100 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        ILO = KMJ
        DO 90 JJ=1,KMJ
           WORK(JJ) = FKMJ*(WORK(JJ+1)-WORK(JJ))
     +                /(WORK(IDL+ILO)+WORK(IDR+JJ))
           ILO = ILO - 1
   90   CONTINUE
  100 CONTINUE
C
C  ---------------------------------
C  COMPUTE IDERIV-TH DERIVATIVE AT X
C  ---------------------------------
C
C     ... B-SPLINE COEFFICIENTS ARE IN WORK(1),..,WORK(K-IDERIV)
C
      DO 120 J=IDERIV+1,KM1
         KMJ = K - J
         ILO = KMJ
         DO 110 JJ=1,KMJ
            WORK(JJ) = (WORK(JJ+1)*WORK(IDL+ILO)+WORK(JJ)*WORK(IDR+JJ))
     +                 /(WORK(IDL+ILO)+WORK(IDR+JJ))
            ILO = ILO - 1
  110    CONTINUE
  120 CONTINUE
      B1VAhL = WORK(1)
      IFLAG = 0
      GO TO 9999
C
C  -----------
C  ERROR EXITS
C  -----------
C
 9001 CONTINUE
      IFLAG = 1
      MESSAG = 'K DOES NOT SATISFY K.GE.1'
      GO TO 9900
C
 9002 CONTINUE
      IFLAG = 2
      MESSAG = 'N DOES NOT SATISFY N.GE.K'
      GO TO 9900
C
 9003 CONTINUE
      IFLAG = 3
      MESSAG = 'IDERIV IS LESS THAN ZERO'
      GO TO 9900
C
 9004 CONTINUE
      IFLAG = 4
      MESSAG = 'NO LEFT LIMITING VALUE AT X=T(N+K)'
      GO TO 9900
C
 9900 CONTINUE
      CALL XERMShG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
C
C  ----
C  EXIT
C  ----
C
 9999 CONTINUE
      RETURN
      END

*DECK B2VA1
      real FUNCTION B2VA1 (XVAL, YVAL, IX, IY, IDX, IDY, TX, TY, NX, NY,
     *     KX, KY, FCN, LDF, WORK, IFLAG)
C***BEGIN PROLOGUE  B2VA1
C***PURPOSE  B2VA1 evaluates the two-dimensional piecewise polynomial
C            interpolating function constructed by the routine B2INhT.
C            Either function values or partial derivative values may be
C            be requested.  Input index values avoid table searches.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (B2VA1-S, DB2VA1-D)
C***KEYWORDS  TWO DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS, EVALUATION, DIFFERENTIATION
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   B2VA1 evaluates the tensor product piecewise polynomial interpolant
C   constructed by the routine B2INhT or one of its derivatives at the
C   point (XVAL,YVAL).
C   This routine is a modification of B2VAhL, and differs in the use
C   of the input indices IX and IY.  It assumes that B2INhT was called
C   with IFLAG = 1, and KX and KY even, in making use of IX and IY.
C
C   The variables IDX and IDY indicate which
C   function is to be evaluated.  If B(x,y) is the interpolant,
C
C                            IDX+IDY
C                           d
C              B2VA1  =   ----------- B (XVAL,YVAL)
C                           IDX   IDY
C                         dx    dy
C
C   Thus, to evaluate the interpolant itself, set IDX=IDY=0. To get the
C   first partial derivative with respect to x, set IDX=1 and IDY=0, and
C   so on.
C
C   Since B is a piecewise polynomial of degree KX-1 in x and KY-1 in y,
C   B2VA1 returns zero whenever IDX.GE.KX or IDY.GE.KY.  B2VA1 also
C   returns zero if (XVAL,YVAL) is out of range, that is, if
C
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY).
C
C   If the knots TX and TY were chosen by B2INhT, then this is
C   equivalent to
C
C              XVAL.LT.X(1) .OR. XVAL.GT.X(NX) .OR.
C              YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY).
C
C   The input quantities TX, TY, NX, NY, KX, KY, and FCN should be
C   unchanged since the call to B2INhT which produced FCN.
C
C   Note that the derivative functions computed by B2VA1 may be
C   discontinuous when x or y correspond to knots.  (When B2INhT selects
C   knots this occurs only when IDX=KX-1 or IDY=KY-1.)  In these cases
C   B2VA1 returns right limiting values (right derivatives), except at
C   the rightmost knot where left limiting values are returned.
C
C   I N P U T
C   ---------
C
C   XVAL    Real scalar
C           X coordinate of evaluation point.
C
C   YVAL    Real scalar
C           Y coordinate of evaluation point.
C
C   IX      Integer scalar (1 .LE. IX .LE. NX)
C           Index of interval containing XVAL in the x data array X
C           input to B2INhT, when X(1) .LE. XVAL .LE. X(NX):
C           X(IX) .LE. XVAL .LT. X(IX+1), with IX = NX if XVAL = X(NX).
C           If XVAL .LT. X(1) or XVAL .GT. X(NX), IX is not used.
C
C   IY      Integer scalar (1 .LE. IY .LE. NY)
C           Index of interval containing YVAL in the y data array Y
C           input to B2INhT, when Y(1) .LE. YVAL .LE. Y(NY):
C           Y(IY) .LE. YVAL .LT. Y(IY+1), with IY = NY if YVAL = Y(NY).
C           If YVAL .LT. Y(1) or YVAL .GT. Y(NY), IY is not used.
C
C   IDX     Integer scalar  (IDX .GE. 0)
C           Indicates the x derivative of piecewise polynomial to
C           evaluate: IDX=J for the Jth partial derivative with respect
C           to x. B2VA1 will return 0 if IDX.GE.KX.
C
C   IDY     Integer scalar  (IDY .GE. 0)
C           Indicates the y derivative of piecewise polynomial to
C           evaluate: IDY=J for the Jth partial derivative with respect
C           to y. B2VA1 will return 0 if IDY.GE.KY.
C
C   TX      Real 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Unchanged since the call to B2INhT which
C           produced FCN.)
C
C   TY      Real 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Unchanged since the call to B2INhT which
C           produced FCN.)
C
C   NX      Integer scalar (NX .GE. KX)
C           The number of interpolation points in x. (Unchanged since
C           the call to B2INhT which produced FCN.)
C
C   NY      Integer scalar (NY .GE. KY)
C           The number of interpolation points in y. (Unchanged since
C           the call to B2INhT which produced FCN.)
C
C   KX      Integer scalar (KX .GE. 2)
C           Order of polynomial pieces in x. (Unchanged since the call
C           to B2INhT which produced FCN.)
C
C   KY      Integer scalar (KY .GE. 2)
C           Order of polynomial pieces in y. (Unchanged since the call
C           to B2INhT which produced FCN.)
C
C   FCN     Real 2D array (size LDF by NY)
C           The B-spline coefficients computed by B2INhT.
C
C   LDF     Integer scalar (LDF .GE. NX)
C           The actual leading dimension of FCN used in the calling
C           program.
C
C
C   O T H E R
C   ---------
C
C   WORK    Real 1D array (size 3*max(KX,KY) + KY)
C           A working storage array.
C
C   IFLAG   Integer scalar.
C           On return IFLAG indicates the status of the output.
C
C           Output:  0 : successful execution
C                    1 : KX out of range
C                    2 : NX or LDF out of range
C                    3 : KY out of range
C                    4 : NY out of range
C                    5 : IDX or IDY out of range
C
C           Note: Error conditions are checked in the order shown here.
C                 Only the first error detected is reported.
C
C
C   B A C K G R O U N D
C   -------------------
C
C   B2VA1 evaluates the following function or its partial derivatives
C
C                        NX   NY
C            B(x,y)  =  SUM  SUM  FCN(i,j) U (x) V (y)
C                       i=1  j=1            i     j
C
C   where
C
C     U (x)  is the ith (one-dimensional) B-spline basis function
C      i     defined by NX, KX, and TX, and
C
C     V (y)  is the jth (one-dimensional) B-spline basis function
C      j     defined by NY, KY, and TY.
C
C   See (de Boor, 1978) for a description of the B-spline basis.
C
C
C***SEE ALSO  B2INhT
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C               Carl de Boor, Efficient computer manipulation of tensor
C                 products, ACM Transactions on Mathematical Software 5,
C                 (1979), pp. 173-182.
C***ROUTINES CALLED  B1VA1, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   820525  DATE WRITTEN
C   900910  Modified to meet SLATEC standards.
C   930510  Eliminated multiple calls to XERMShG.  (FNF)
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   931217  Renamed and added IX, IY arguments. (ACH)
C***END PROLOGUE  B2VA1
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          IDX, IDY, IX, IY, NX, NY, KX, KY, LDF, IFLAG
      real             XVAL, YVAL, TX(NX+KX), TY(NY+KY), FCN(LDF,NY),
     +                 WORK(*)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'SLATEC', SUBROU = 'B2VA1')
      CHARACTER*50     MESSAG
      INTEGER          I, IERR, JMIN, JMAX, IPOS, IW, LEFTX, LEFTY
      real             B1VA1
C
      EXTERNAL         INTRhV, B1VA1, XERMShG
C
C***FIRST EXECUTABLE STATEMENT  B2VA1
      B2VA1 = 0
      IFLAG = 0
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IF (KX .LT. 1)  GO TO 9001
      IF ((NX .LT. KX) .OR. (NX .GT. LDF)) GO TO 9002
      IF (KY .LT. 1)  GO TO 9003
      IF (NY .LT. KY)  GO TO 9004
      IF ((IDX .LT. 0) .OR. (IDY .LT. 0)) GO TO 9005
C
C     ... QUICK RETURN WHEN SPLINE EVALUATES TO ZERO
C
      IF ((IDX .GE. KX) .OR. (IDY .GE. KY)) GO TO 9999
      IF ((XVAL .LT. TX(1)) .OR. (XVAL .GT. TX(NX+KX))) GO TO 9999
      IF ((YVAL .LT. TY(1)) .OR. (YVAL .GT. TY(NY+KY))) GO TO 9999
C
C  ---------------------
C  EVALUATE THE B-SPLINE
C  ---------------------
C
C  The summation in the Background section above can be rewritten as
C
C                                NY
C  (1)               B(x,y)  =  SUM  b  V (y),
C                               j=1   j  j
C  where
C                         NX
C  (2)             b  =  SUM  FCN(i,j) U (x),   j = 1,..,NY.
C                   j    i=1            i
C
C  Note that each summation is the evaluation of a one-dimensional
C  B-spline, which can be done by calls to B1VA1.  At most KY basis
C  functions in (1) are nonzero.  The indices of these functions are
C  determined, and only the coefficients b-sub-j which multiply nonzero
C  basis functions are computed using (2).
C
C
C     ... SET KNOT INTERVAL CONTAINING YVAL
C
      LEFTY = IY + 2
      IF (LEFTY .LT. KY) LEFTY = KY
      IF (LEFTY .GT. NY) LEFTY = NY
C
C     ... SET RANGE OF INDICES OF NONZERO BASIS FUNCTIONS IN (1)
C         (JMIN,JMAX) = (SMALLEST,LARGEST)
C
         JMIN = LEFTY - KY + 1
         JMAX = LEFTY
C
C     ... GET COEFFICIENTS OF NONZERO BASIS FUNCTIONS IN (1) USING (2)
C
      IPOS = 0
      IW = KY + 1
      LEFTX = IX + KX/2
      IF (LEFTX .LT. KX) LEFTX = KX
      IF (LEFTX .GT. NX) LEFTX = NX
      DO 30 I=JMIN,JMAX
         IPOS = IPOS + 1
         WORK(IPOS) = B1VA1(XVAL,LEFTX,IDX,TX,NX,KX,FCN(1,I),
     +                      WORK(IW),IERR)
   30 CONTINUE
C
C     ... EVALUATE INTERPOLANT USING (1)
C
      B2VA1 = B1VA1(YVAL,KY,IDY,TY(JMIN),KY,KY,WORK,WORK(IW),IERR)
      GO TO 9999
C
C  -----
C  EXITS
C  -----
C
 9001 CONTINUE
      IFLAG = 1
      MESSAG = 'KX IS OUT OF RANGE'
      GO TO 9900
C
 9002 CONTINUE
      IFLAG = 2
      MESSAG = 'NX OR LDF IS OUT OF RANGE'
      GO TO 9900
C
 9003 CONTINUE
      IFLAG = 3
      MESSAG = 'KY IS OUT OF RANGE'
      GO TO 9900
C
 9004 CONTINUE
      IFLAG = 4
      MESSAG = 'NY IS OUT OF RANGE'
      GO TO 9900
C
 9005 CONTINUE
      IFLAG = 5
      MESSAG = 'IDX OR IDY IS OUT OF RANGE'
      GO TO 9900
C
 9900 CONTINUE
      CALL XERMShG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
C
C     ... EXIT
C
 9999 CONTINUE
      RETURN
      END

*DECK B1VA1
      real FUNCTION B1VA1 (X, I, IDERIV, T, N, K, A, WORK, IFLAG)
C***BEGIN PROLOGUE  B1VA1
C***PURPOSE  Evaluates the B-representation of a spline at X for the
C            function value or any of its derivatives.
C            Uses input index to avoid table search.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (B1VA1-S, DB1VA1-D)
C***KEYWORDS  B-SPLINE, DIFFERENTIATION, EVALUATION, SPLINES
C***AUTHOR  Boisvert, R. F., (NIST)
C             Center for Computing and Applied Mathematics
C             National Institute of Standards and Technology
C             Gaithersburg, MD 20899
C***DESCRIPTION
C
C   B1VA1 evaluates the B-representation (T,N,K,A) of a spline or one
C   of its derivatives at X.  Variable IDERIV indicates which function
C   is to be evaluated.  If B(x) is the spline,
C
C                              IDERIV
C                             d
C                  B1VA1  =   -------- B(X)
C                               IDERIV
C                             dx
C
C   Thus, to evaluate the spline itself, set IDERIV=0, to get the first
C   derivative with respect to x, set IDERIV=1, and so on.
C
C   Since B is a piecewise polynomial of degree K-1 B1VA1 returns zero
C   whenever IDERIV.GE.K.  B1VA1 also returns zero if X is out of
C   range, that is, if
C
C                    X.LT.T(1) .OR. X.GT.T(N+K).
C
C   Note that the derivative functions computed by B1VA1 may be
C   discontinuous when X corresponds to a knot.   In these cases B1VA1
C   returns right limiting values (right derivatives), except at T(N+K)
C   where left limiting values are returned.
C
C
C   I N P U T
C   ---------
C
C   X       Real scalar
C           The point at which the B-spline is to be evaluated.
C
C   I       Integer scalar (1 .LE. I .LE. N+K)
C           Index of interval in data containing X.  This is the
C           largest I in (1,N+K) such that T(I) .LE. X .LT. T(I+1).
C
C   IDERIV  Integer scalar  (IDERIV .GE. 0)
C           Indicates the derivative of the B-spline to evaluate:
C           IDERIV=0 for the function itself, IDERIV=J for the Jth
C           derivative. B1VA1 will return 0 if IDERIV.GE.K.
C
C   T       Real 1D array (size N+K)
C           Sequence of knots defining the B-spline.
C
C   N       Integer scalar
C           The number of B-spline coefficients. (N = sum of knot
C           multiplicities - K.)
C
C   K       Integer scalar
C           Order of the B-spline.
C
C   A       Real 1D array (size N)
C           The B-spline coefficients.
C
C
C   O T H E R
C   ---------
C
C   WORK    Real 1D array (size 3*K)
C           Workspace.
C
C   IFLAG   Integer scalar.
C           On return IFLAG indicates the status of the output.
C
C           Output:  0 : successful execution
C                    1 : K out of range
C                    2 : N out of range
C                    3 : IDERIV out of range
C                    4 : Unable to find left-limiting value at T(N+K)
C
C
C  B1VA1 is a version of the routine BVALUE written by Carl de Boor and
C  given in the reference.
C
C  This routine replaces the former SLATEC routine BVALU.  It differs
C  from BVALU in that
C
C    (1) the order of arguments in the calling sequence has been
C        changed and a new parameter, IFLAG, has been added,
C    (2) X may be any real number,
C    (2) IDERIV may be any non-negative integer,
C    (3) right-limiting values are returned when X=T(N+1) and
C        T(N+K).GT.T(N+1), and
C    (4) the prescription given in BVALU for obtaining left-limiting
C        values at interior knots does not work.
C
C
C***REFERENCES  Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  INTRhV, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   900726  DATE WRITTEN
C   930618  Reformatted the AUTHOR and REFERENCES sections.  (WRB)
C   931217  Renamed, added input index I, dropped INBV. (ACH)
C***END PROLOGUE  B1VA1
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  .. ARGUMENTS
C
      INTEGER          I, N, K, IDERIV, IFLAG
      real             T(N+K), A(N), X, WORK(3*K)
C
C  .. LOCAL VARIABLES
C
      CHARACTER*6      LIBRAR, SUBROU
      PARAMETER        (LIBRAR = 'SLATEC', SUBROU = 'B1VA1')
      CHARACTER*50     MESSAG
      INTEGER          IDL, IDR, ILO, IMK, IP1, J, JJ, JMAX,
     +                 KM1, KMJ, NPK
      real             FKMJ
C
      EXTERNAL         INTRhV, XERMShG
C
C***FIRST EXECUTABLE STATEMENT  B1VA1
      B1VA1 = 0
      NPK = N + K
C
C  -------------
C  SPECIAL CASES
C  -------------
C
C     ... CHECK INPUT FOR ERRORS
C
      IF (K .LT. 1)  GO TO 9001
      IF (N .LT. K)  GO TO 9002
      IF (IDERIV .LT. 0)  GO TO 9003
C
C     ... QUICK RETURN WHEN SPLINE EVALUATES TO ZERO
C
      IF (IDERIV .GE. K)  GO TO 9999
      IF ((X .LT. T(1)) .OR. (X .GT. T(NPK)))  GO TO 9999
C
C  ------------------------------------------------------
C  STORE COEFFICIENTS OF NON-ZERO BASIS FUNCTIONS IN WORK
C  ------------------------------------------------------
C
      IMK = I - K
      DO 30 J=1,K
         WORK(J) = A(IMK+J)
   30 CONTINUE
C
C  ----------------------------
C  COMPUTE AUXILIARY QUANTITIES
C  ----------------------------
C
C     ... DL(J) = X - T(I+1-J) = WORK(IDL+J), J=1,..,K-1
C     ... DR(J) = T(I+J) - X   = WORK(IDR+J), J=1,..,K-1
C
      IDL = K
      IDR = IDL + K
      KM1 = K - 1
      IP1 = I + 1
      JMAX = MIN(I,KM1)
      DO 50 J=1,JMAX
         WORK(IDL+J) = X - T(IP1-J)
   50 CONTINUE
      DO 60 J=JMAX+1,KM1
         WORK(IDL+J) = WORK(IDL+JMAX)
   60 CONTINUE
      JMAX = MIN(NPK-I,KM1)
      DO 70 J=1,JMAX
         WORK(IDR+J) = T(I+J) - X
   70 CONTINUE
      DO 80 J=JMAX+1,KM1
         WORK(IDR+J) = WORK(IDR+JMAX)
   80 CONTINUE
C
C  ----------------------------------------
C  DIFFERENCE THE COEFFICIENTS IDERIV TIMES
C  ----------------------------------------
C
      DO 100 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        ILO = KMJ
        DO 90 JJ=1,KMJ
           WORK(JJ) = FKMJ*(WORK(JJ+1)-WORK(JJ))
     +                /(WORK(IDL+ILO)+WORK(IDR+JJ))
           ILO = ILO - 1
   90   CONTINUE
  100 CONTINUE
C
C  ---------------------------------
C  COMPUTE IDERIV-TH DERIVATIVE AT X
C  ---------------------------------
C
C     ... B-SPLINE COEFFICIENTS ARE IN WORK(1),..,WORK(K-IDERIV)
C
      DO 120 J=IDERIV+1,KM1
         KMJ = K - J
         ILO = KMJ
         DO 110 JJ=1,KMJ
            WORK(JJ) = (WORK(JJ+1)*WORK(IDL+ILO)+WORK(JJ)*WORK(IDR+JJ))
     +                 /(WORK(IDL+ILO)+WORK(IDR+JJ))
            ILO = ILO - 1
  110    CONTINUE
  120 CONTINUE
      B1VA1 = WORK(1)
      IFLAG = 0
      GO TO 9999
C
C  -----------
C  ERROR EXITS
C  -----------
C
 9001 CONTINUE
      IFLAG = 1
      MESSAG = 'K DOES NOT SATISFY K.GE.1'
      GO TO 9900
C
 9002 CONTINUE
      IFLAG = 2
      MESSAG = 'N DOES NOT SATISFY N.GE.K'
      GO TO 9900
C
 9003 CONTINUE
      IFLAG = 3
      MESSAG = 'IDERIV IS LESS THAN ZERO'
      GO TO 9900
C
 9004 CONTINUE
      IFLAG = 4
      MESSAG = 'NO LEFT LIMITING VALUE AT X=T(N+K)'
      GO TO 9900
C
 9900 CONTINUE
      CALL XERMShG(LIBRAR,SUBROU,MESSAG,IFLAG,1)
C
C  ----
C  EXIT
C  ----
C
 9999 CONTINUE
      RETURN
      END

c----------------------------------------------------------------------c

*DECK FCh
      SUBROUTINE FCh (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, BKPT,
     +   NCONST, XCONST, YCONST, NDERIV, MODE, COEFF, W, IW)
c ... This is a special version of slatec FCh called FChh where changes
c ... have been made to allow conversion from real to double precision
c ..  using mmm - T. Rognlien 2/3/96 (REAL --> real, etc.)
C***BEGIN PROLOGUE  FCh
C***PURPOSE  Fit a piecewise polynomial curve to discrete data.
C            The piecewise polynomials are represented as B-splines.
C            The fitting is done in a weighted least squares sense.
C            Equality and inequality constraints can be imposed on the
C            fitted curve.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A1, K1A2A, L8A3
C***TYPE      SINGLE PRECISION (FCh-S, DFC-D)
C***KEYWORDS  B-SPLINE, CONSTRAINED LEAST SQUARES, CURVE FITTING,
C             WEIGHTED LEAST SQUARES
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C      This subprogram fits a piecewise polynomial curve
C      to discrete data.  The piecewise polynomials are
C      represented as B-splines.
C      The fitting is done in a weighted least squares sense.
C      Equality and inequality constraints can be imposed on the
C      fitted curve.
C
C      For a description of the B-splines and usage instructions to
C      evaluate them, see
C
C      C. W. de Boor, Package for Calculating with B-Splines.
C                     SIAM J. Numer. Anal., p. 441, (June, 1977).
C
C      For further documentation and discussion of constrained
C      curve fitting using B-splines, see
C
C      R. J. Hanson, Constrained Least Squares Curve Fitting
C                   to Discrete Data Using B-Splines, a User's
C                   Guide. Sandia Labs. Tech. Rept. SAND-78-1291,
C                   December, (1978).
C
C  Input..
C      NDATA,XDATA(*),
C      YDATA(*),
C      SDDATA(*)
C                         The NDATA discrete (X,Y) pairs and the Y value
C                         standard deviation or uncertainty, SD, are in
C                         the respective arrays XDATA(*), YDATA(*), and
C                         SDDATA(*).  No sorting of XDATA(*) is
C                         required.  Any non-negative value of NDATA is
C                         allowed.  A negative value of NDATA is an
C                         error.  A zero value for any entry of
C                         SDDATA(*) will weight that data point as 1.
C                         Otherwise the weight of that data point is
C                         the reciprocal of this entry.
C
C      NORD,NBKPT,
C      BKPT(*)
C                         The NBKPT knots of the B-spline of order NORD
C                         are in the array BKPT(*).  Normally the
C                         problem data interval will be included between
C                         the limits BKPT(NORD) and BKPT(NBKPT-NORD+1).
C                         The additional end knots BKPT(I),I=1,...,
C                         NORD-1 and I=NBKPT-NORD+2,...,NBKPT, are
C                         required to compute the functions used to fit
C                         the data.  No sorting of BKPT(*) is required.
C                         Internal to  FCh( ) the extreme end knots may
C                         be reduced and increased respectively to
C                         accommodate any data values that are exterior
C                         to the given knot values.  The contents of
C                         BKPT(*) is not changed.
C
C                         NORD must be in the range 1 .LE. NORD .LE. 20.
C                         The value of NBKPT must satisfy the condition
C                         NBKPT .GE. 2*NORD.
C                         Other values are considered errors.
C
C                         (The order of the spline is one more than the
C                         degree of the piecewise polynomial defined on
C                         each interval.  This is consistent with the
C                         B-spline package convention.  For example,
C                         NORD=4 when we are using piecewise cubics.)
C
C      NCONST,XCONST(*),
C      YCONST(*),NDERIV(*)
C                         The number of conditions that constrain the
C                         B-spline is NCONST.  A constraint is specified
C                         by an (X,Y) pair in the arrays XCONST(*) and
C                         YCONST(*), and by the type of constraint and
C                         derivative value encoded in the array
C                         NDERIV(*).  No sorting of XCONST(*) is
C                         required.  The value of NDERIV(*) is
C                         determined as follows.  Suppose the I-th
C                         constraint applies to the J-th derivative
C                         of the B-spline.  (Any non-negative value of
C                         J < NORD is permitted.  In particular the
C                         value J=0 refers to the B-spline itself.)
C                         For this I-th constraint, set
C                          XCONST(I)=X,
C                          YCONST(I)=Y, and
C                          NDERIV(I)=ITYPE+4*J, where
C
C                          ITYPE = 0,      if (J-th deriv. at X) .LE. Y.
C                                = 1,      if (J-th deriv. at X) .GE. Y.
C                                = 2,      if (J-th deriv. at X) .EQ. Y.
C                                = 3,      if (J-th deriv. at X) .EQ.
C                                             (J-th deriv. at Y).
C                          (A value of NDERIV(I)=-1 will cause this
C                          constraint to be ignored.  This subprogram
C                          feature is often useful when temporarily
C                          suppressing a constraint while still
C                          retaining the source code of the calling
C                          program.)
C
C        MODE
C                         An input flgag that directs the least squares
C                         solution method used by FCh( ).
C
C                         The variance function, referred to below,
C                         defines the square of the probable error of
C                         the fitted curve at any point, XVAL.
C                         This feature of  FCh( ) allows one to use the
C                         square root of this variance function to
C                         determine a probable error band around the
C                         fitted curve.
C
C                         =1  a new problem.  No variance function.
C
C                         =2  a new problem.  Want variance function.
C
C                         =3  an old problem.  No variance function.
C
C                         =4  an old problem.  Want variance function.
C
C                         Any value of MODE other than 1-4 is an error.
C
C                         The user with a new problem can skip directly
C                         to the description of the input parameters
C                         IW(1), IW(2).
C
C                         If the user correctly specifies the new or old
C                         problem status, the subprogram FCh( ) will
C                         perform more efficiently.
C                         By an old problem it is meant that subprogram
C                         FCh( ) was last called with this same set of
C                         knots, data points and weights.
C
C                         Another often useful deployment of this old
C                         problem designation can occur when one has
C                         previously obtained a Q-R orthogonal
C                         decomposition of the matrix resulting from
C                         B-spline fitting of data (without constraints)
C                         at the breakpoints BKPT(I), I=1,...,NBKPT.
C                         For example, this matrix could be the result
C                         of sequential accumulation of the least
C                         squares equations for a very large data set.
C                         The user writes this code in a manner
C                         convenient for the application.  For the
C                         discussion here let
C
C                                      N=NBKPT-NORD, and K=N+3
C
C                         Let us assume that an equivalent least squares
C                         system
C
C                                      RC=D
C
C                         has been obtained.  Here R is an N+1 by N
C                         matrix and D is a vector with N+1 components.
C                         The last row of R is zero.  The matrix R is
C                         upper triangular and banded.  At most NORD of
C                         the diagonals are nonzero.
C                         The contents of R and D can be copied to the
C                         working array W(*) as follows.
C
C                         The I-th diagonal of R, which has N-I+1
C                         elements, is copied to W(*) starting at
C
C                                      W((I-1)*K+1),
C
C                         for I=1,...,NORD.
C                         The vector D is copied to W(*) starting at
C
C                                      W(NORD*K+1)
C
C                         The input value used for NDATA is arbitrary
C                         when an old problem is designated.  Because
C                         of the feature of FCh( ) that checks the
C                         working storage array lengths, a value not
C                         exceeding NBKPT should be used.  For example,
C                         use NDATA=0.
C
C                         (The constraints or variance function request
C                         can change in each call to FCh( ).)  A new
C                         problem is anything other than an old problem.
C
C      IW(1),IW(2)
C                         The amounts of working storage actually
C                         allocated for the working arrays W(*) and
C                         IW(*).  These quantities are compared with the
C                         actual amounts of storage needed in FCh( ).
C                         Insufficient storage allocated for either
C                         W(*) or IW(*) is an error.  This feature was
C                         included in FCh( ) because misreading the
C                         storage formulas for W(*) and IW(*) might very
C                         well lead to subtle and hard-to-find
C                         programming bugs.
C
C                         The length of W(*) must be at least
C
C                           NB=(NBKPT-NORD+3)*(NORD+1)+
C                               2*MAX(NDATA,NBKPT)+NBKPT+NORD**2
C
C                         Whenever possible the code uses banded matrix
C                         processors BNDAChC( ) and BNDSOhL( ).  These
C                         are utilized if there are no constraints,
C                         no variance function is required, and there
C                         is sufficient data to uniquely determine the
C                         B-spline coefficients.  If the band processors
C                         cannot be used to determine the solution,
C                         then the constrained least squares code LSEIh
C                         is used.  In this case the subprogram requires
C                         an additional block of storage in W(*).  For
C                         the discussion here define the integers NEQCON
C                         and NINCON respectively as the number of
C                         equality (ITYPE=2,3) and inequality
C                         (ITYPE=0,1) constraints imposed on the fitted
C                         curve.  Define
C
C                           L=NBKPT-NORD+1
C
C                         and note that
C
C                           NCONST=NEQCON+NINCON.
C
C                         When the subprogram FCh( ) uses LSEIh( ) the
C                         length of the working array W(*) must be at
C                         least
C
C                           LW=NB+(L+NCONST)*L+
C                              2*(NEQCON+L)+(NINCON+L)+(NINCON+2)*(L+6)
C
C                         The length of the array IW(*) must be at least
C
C                           IW1=NINCON+2*L
C
C                         in any case.
C
C  Output..
C      MODE
C                         An output flag that indicates the status
C                         of the constrained curve fit.
C
C                         =-1  a usage error of FCh( ) occurred.  The
C                              offending condition is noted with the
C                              SLATEC library error processor, XWEMShG.
C                              In case the working arrays W(*) or IW(*)
C                              are not long enough, the minimal
C                              acceptable length is printed.
C
C                         = 0  successful constrained curve fit.
C
C                         = 1  the requested equality constraints
C                              are contradictory.
C
C                         = 2  the requested inequality constraints
C                              are contradictory.
C
C                         = 3  both equality and inequality constraints
C                              are contradictory.
C
C      COEFF(*)
C                         If the output value of MODE=0 or 1, this array
C                         contains the unknowns obtained from the least
C                         squares fitting process.  These N=NBKPT-NORD
C                         parameters are the B-spline coefficients.
C                         For MODE=1, the equality constraints are
C                         contradictory.  To make the fitting process
C                         more robust, the equality constraints are
C                         satisfied in a least squares sense.  In this
C                         case the array COEFF(*) contains B-spline
C                         coefficients for this extended concept of a
C                         solution.  If MODE=-1,2 or 3 on output, the
C                         array COEFF(*) is undefined.
C
C  Working Arrays..
C      W(*),IW(*)
C                         These arrays are respectively typed real and
C                         INTEGER.
C                         Their required lengths are specified as input
C                         parameters in IW(1), IW(2) noted above.  The
C                         contents of W(*) must not be modified by the
C                         user if the variance function is desired.
C
C  Evaluating the
C  Variance Function..
C                         To evaluate the variance function (assuming
C                         that the uncertainties of the Y values were
C                         provided to  FCh( ) and an input value of
C                         MODE=2 or 4 was used), use the function
C                         subprogram  CV( )
C
C                           VAR=CV(XVAL,NDATA,NCONST,NORD,NBKPT,
C                                  BKPT,W)
C
C                         Here XVAL is the point where the variance is
C                         desired.  The other arguments have the same
C                         meaning as in the usage of FCh( ).
C
C                         For those users employing the old problem
C                         designation, let MDATA be the number of data
C                         points in the problem.  (This may be different
C                         from NDATA if the old problem designation
C                         feature was used.)  The value, VAR, should be
C                         multiplied by the quantity
C
C                         real(MAX(NDATA-N,1))/MAX(MDATA-N,1)
C
C                         The output of this subprogram is not defined
C                         if an input value of MODE=1 or 3 was used in
C                         FCh( ) or if an output value of MODE=-1, 2, or
C                         3 was obtained.  The variance function, except
C                         for the scaling factor noted above, is given
C                         by
C
C                           VAR=(transpose of B(XVAL))*C*B(XVAL)
C
C                         The vector B(XVAL) is the B-spline basis
C                         function values at X=XVAL.
C                         The covariance matrix, C, of the solution
C                         coefficients accounts only for the least
C                         squares equations and the explicitly stated
C                         equality constraints.  This fact must be
C                         considered when interpreting the variance
C                         function from a data fitting problem that has
C                         inequality constraints on the fitted curve.
C
C  Evaluating the
C  Fitted Curve..
C                         To evaluate derivative number IDER at XVAL,
C                         use the function subprogram BVALU( ).
C
C                           F = BVALU(BKPT,COEFF,NBKPT-NORD,NORD,IDER,
C                                      XVAL,INBV,WORKB)
C
C                         The output of this subprogram will not be
C                         defined unless an output value of MODE=0 or 1
C                         was obtained from  FCh( ), XVAL is in the data
C                         interval, and IDER is nonnegative and .LT.
C                         NORD.
C
C                         The first time BVALU( ) is called, INBV=1
C                         must be specified.  This value of INBV is the
C                         overwritten by BVALU( ).  The array WORKB(*)
C                         must be of length at least 3*NORD, and must
C                         not be the same as the W(*) array used in
C                         the call to FCh( ).
C
C                         BVALU( ) expects the breakpoint array BKPT(*)
C                         to be sorted.
C
C***REFERENCES  R. J. Hanson, Constrained least squares curve fitting
C                 to discrete data using B-splines, a users guide,
C                 Report SAND78-1291, Sandia Laboratories, December
C                 1978.
C***ROUTINES CALLED  FChMN
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Convert references to XERRWV to references to XERMShG.  (RWC)
C   900607  Editorial changes to Prologue to make Prologues for EFCh,
C           DEFC, FCh, and DFCh look as much the same as possible.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960203  Changed name to FCh, REAL --> real, etc. (TDR)
C***END PROLOGUE  FCh
      implicit none
      real             BKPT(*), COEFF(*), SDDATA(*), W(*), XCONST(*),
     *   XDATA(*), YCONST(*), YDATA(*)
      INTEGER IW(*), MODE, NBKPT, NCONST, NDATA, NDERIV(*), NORD
C
      EXTERNAL FCMNh
C
      INTEGER I1, I2, I3, I4, I5, I6, I7, MDG, MDW
C
C***FIRST EXECUTABLE STATEMENT  FCh
      MDG = NBKPT - NORD + 3
      MDW = NBKPT - NORD + 1 + NCONST
C                         USAGE IN FCMNh( ) OF W(*)..
C     I1,...,I2-1      G(*,*)
C
C     I2,...,I3-1      XTEMP(*)
C
C     I3,...,I4-1      PTEMP(*)
C
C     I4,...,I5-1      BKPT(*) (LOCAL TO FCMNh( ))
C
C     I5,...,I6-1      BF(*,*)
C
C     I6,...,I7-1      W(*,*)
C
C     I7,...           WORK(*) FOR LSEIh( )
C
      I1 = 1
      I2 = I1 + MDG*(NORD+1)
      I3 = I2 + MAX(NDATA,NBKPT)
      I4 = I3 + MAX(NDATA,NBKPT)
      I5 = I4 + NBKPT
      I6 = I5 + NORD*NORD
      I7 = I6 + MDW*(NBKPT-NORD+1)
      CALL FCMNh(NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT, BKPT, NCONST,
     1   XCONST, YCONST, NDERIV, MODE, COEFF, W(I5), W(I2), W(I3),
     2   W(I4), W(I1), MDG, W(I6), MDW, W(I7), IW)
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK FCMNh
      SUBROUTINE FCMNh (NDATA, XDATA, YDATA, SDDATA, NORD, NBKPT,BKPTIN,
     +   NCONST, XCONST, YCONST, NDERIV, MODE, COEFF, BF, XTEMP, PTEMP,
     +   BKPT, G, MDG, W, MDW, WORK, IWORK)
C***BEGIN PROLOGUE  FCMNh
C***SUBSIDIARY
C***PURPOSE  Subsidiary to FCh
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (FCMNh-S, DFCMN-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This is a companion subprogram to FCh( ).
C     The documentation for FCh( ) has complete usage instructions.
C
C***SEE ALSO  FCh
C***ROUTINES CALLED  BNDAChC, BNDSOhL, BSPLVhD, BSPLVhN, LSEIh, SAXPY,
C                    SCOPY, SSCAL, SSORhT, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and extensively revised (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Convert XERRWV calls to XERMShG calls.  (RWC)
c   960203  Changed name to FCMNh with REAL --> real, etc. for mmm (TDR)
c   960214  Eliminated "broadcast" usage of SCOPY by using do loops (GRS)
c   960214  Made DUMMY an array because it is passed to a subroutine
c           expecting an array (GRS)
C***END PROLOGUE  FCMNh
      implicit none
      INTEGER IWORK(*), MDG, MDW, MODE, NBKPT, NCONST, NDATA, NDERIV(*),
     *   NORD
      real             BF(NORD,*), BKPT(*), BKPTIN(*), COEFF(*),
     *   G(MDG,*), PTEMP(*), SDDATA(*), W(MDW,*), WORK(*),
     *   XCONST(*), XDATA(*), XTEMP(*), YCONST(*), YDATA(*)
C
      EXTERNAL BNDAChC, BNDSOhL, BSPLVhD, BSPLVhN, LSEIh, SAXPY, SCOPY,
     *    SSCAL, SSORhT, XERMShG
C
      real            DUMMY(1), PRGOPT(10), RNORM, RNORME, RNORML, XMAX,
     *   XMIN, XVAL, YVAL
      INTEGER I, IDATA, IDERIV, ILEFT, INTRVL, INTW1, IP, IR, IROW,
     *   ITYPE, IW1, IW2, L, LW, MT, N, NB, NEQCON, NINCON, NORDM1,
     *   NORDP1, NP1
      integer i0
      LOGICAL BAND, NEW, VAR
      CHARACTER*8 XERN1
C
C***FIRST EXECUTABLE STATEMENT  FCMNh
C
C     Analyze input.
C
      IF (NORD.LT.1 .OR. NORD.GT.20) THEN
         CALL XERMShG ('SLATEC', 'FCMNh',
     +      'IN FCh, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',
     +      2, 1)
         MODE = -1
         RETURN
C
      ELSEIF (NBKPT.LT.2*NORD) THEN
         CALL XERMShG ('SLATEC', 'FCMNh',
     +      'IN FCh, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE ' //
     +      'THE B-SPLINE ORDER.', 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
      IF (NDATA.LT.0) THEN
         CALL XERMShG ('SLATEC', 'FCMNh',
     +      'IN FCh, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.',
     +      2, 1)
         MODE = -1
         RETURN
      ENDIF
C
C     Amount of storage allocated for W(*), IW(*).
C
      IW1 = IWORK(1)
      IW2 = IWORK(2)
      NB = (NBKPT-NORD+3)*(NORD+1) + 2*MAX(NDATA,NBKPT) + NBKPT +
     +     NORD**2
C
C     See if sufficient storage has been allocated.
C
      IF (IW1.LT.NB) THEN
         WRITE (XERN1, '(I8)') NB
         CALL XERMShG ('SLATEC', 'FCMNh',
     *      'IN FCh, INSUFFICIENT STORAGE FOR W(*).  CHECK NB = ' //
     *      XERN1, 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
      IF (MODE.EQ.1) THEN
         BAND = .true.
         VAR = .false.
         NEW = .true.
      ELSEIF (MODE.EQ.2) THEN
         BAND = .false.
         VAR = .true.
         NEW = .true.
      ELSEIF (MODE.EQ.3) THEN
         BAND = .true.
         VAR = .false.
         NEW = .false.
      ELSEIF (MODE.EQ.4) THEN
         BAND = .false.
         VAR = .true.
         NEW = .false.
      ELSE
         CALL XERMShG ('SLATEC', 'FCMNh',
     +      'IN FCh, INPUT VALUE OF MODE MUST BE 1-4.', 2, 1)
         MODE = -1
         RETURN
      ENDIF
      MODE = 0
C
C     Sort the breakpoints.
C
      CALL SCOPY (NBKPT, BKPTIN, 1, BKPT, 1)
      CALL SSORhT (BKPT, DUMMY, NBKPT, 1)
C
C     Initialize variables.
C
      NEQCON = 0
      NINCON = 0
      DO 100 I = 1,NCONST
         L = NDERIV(I)
         ITYPE = MOD(L,4)
         IF (ITYPE.LT.2) THEN
            NINCON = NINCON + 1
         ELSE
            NEQCON = NEQCON + 1
         ENDIF
  100 CONTINUE
C
C     Compute the number of variables.
C
      N = NBKPT - NORD
      NP1 = N + 1
      LW = NB + (NP1+NCONST)*NP1 + 2*(NEQCON+NP1) + (NINCON+NP1) +
     +     (NINCON+2)*(NP1+6)
      INTW1 = NINCON + 2*NP1
C
C     Save interval containing knots.
C
      XMIN = BKPT(NORD)
      XMAX = BKPT(NP1)
C
C     Find the smallest referenced independent variable value in any
C     constraint.
C
      DO 110 I = 1,NCONST
         XMIN = MIN(XMIN,XCONST(I))
         XMAX = MAX(XMAX,XCONST(I))
  110 CONTINUE
      NORDM1 = NORD - 1
      NORDP1 = NORD + 1
C
C     Define the option vector PRGOPT(1-10) for use in LSEIh( ).
C
      PRGOPT(1) = 4
C
C     Set the covariance matrix computation flag.
C
      PRGOPT(2) = 1
      IF (VAR) THEN
         PRGOPT(3) = 1
      ELSE
         PRGOPT(3) = 0
      ENDIF
C
C     Increase the rank determination tolerances for both equality
C     constraint equations and least squares equations.
C
      PRGOPT(4) = 7
      PRGOPT(5) = 4
      PRGOPT(6) = 1.E-4
C
      PRGOPT(7) = 10
      PRGOPT(8) = 5
      PRGOPT(9) = 1.E-4
C
      PRGOPT(10) = 1
C
C     Turn off work array length checking in LSEIh( ).
C
      IWORK(1) = 0
      IWORK(2) = 0
C
C     Initialize variables and analyze input.
C
      IF (NEW) THEN
C
C        To process least squares equations sort data and an array of
C        pointers.
C
         CALL SCOPY (NDATA, XDATA, 1, XTEMP, 1)
         DO 120 I = 1,NDATA
            PTEMP(I) = I
  120    CONTINUE
C
         IF (NDATA.GT.0) THEN
            CALL SSORhT (XTEMP, PTEMP, NDATA, 2)
            XMIN = MIN(XMIN,XTEMP(1))
            XMAX = MAX(XMAX,XTEMP(NDATA))
         ENDIF
C
C        Fix breakpoint array if needed.
C
         DO 130 I = 1,NORD
            BKPT(I) = MIN(BKPT(I),XMIN)
  130    CONTINUE
C
         DO 140 I = NP1,NBKPT
            BKPT(I) = MAX(BKPT(I),XMAX)
  140    CONTINUE
C
C        Initialize parameters of banded matrix processor, BNDAChC( ).
C
         MT = 0
         IP = 1
         IR = 1
         ILEFT = NORD
         DO 160 IDATA = 1,NDATA
C
C           Sorted indices are in PTEMP(*).
C
            L = PTEMP(IDATA)
            XVAL = XDATA(L)
C
C           When interval changes, process equations in the last block.
C
            IF (XVAL.GE.BKPT(ILEFT+1)) THEN
               CALL BNDAChC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
               MT = 0
C
C              Move pointer up to have BKPT(ILEFT).LE.XVAL,
C                 ILEFT.LT.NP1.
C
  150          IF (XVAL.GE.BKPT(ILEFT+1) .AND. ILEFT.LT.N) THEN
                  ILEFT = ILEFT + 1
                  GO TO 150
               ENDIF
            ENDIF
C
C           Obtain B-spline function value.
C
            CALL BSPLVhN (BKPT, NORD, 1, XVAL, ILEFT, BF)
C
C           Move row into place.
C
            IROW = IR + MT
            MT = MT + 1
            CALL SCOPY (NORD, BF, 1, G(IROW,1), MDG)
            G(IROW,NORDP1) = YDATA(L)
C
C           Scale data if uncertainty is nonzero.
C
            IF (SDDATA(L).NE.0.E0) CALL SSCAL (NORDP1, 1.E0/SDDATA(L),
     +                                  G(IROW,1), MDG)
C
C           When staging work area is exhausted, process rows.
C
            IF (IROW.EQ.MDG-1) THEN
               CALL BNDAChC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
               MT = 0
            ENDIF
  160    CONTINUE
C
C        Process last block of equations.
C
         CALL BNDAChC (G, MDG, NORD, IP, IR, MT, ILEFT-NORDM1)
C
C        Last call to adjust block positioning.
C
cmer (12 Sep 2006) REINSTATE SCOPY deleted by Gary Smith
cmer (11 Nov 2021) Change to SFILL as SCOPY can't be called with a scalar
         CALL SFILL (NORDP1, 0.E0, G(IR,1), MDG)
cmer         do i0 = 1, nordp1*mdg, mdg
cmer            g(ir-1+i0,1) = 0.
ctdr can also fix by using:           if (ir-1+i0 .le. mdg) g(ir-1+i0,1) = 0.
cmer         enddo
         CALL BNDAChC (G, MDG, NORD, IP, IR, 1, NP1)
      ENDIF
C
      BAND = BAND .AND. NCONST.EQ.0
      DO 170 I = 1,N
         BAND = BAND .AND. G(I,1).NE.0.E0
  170 CONTINUE
C
C     Process banded least squares equations.
C
      IF (BAND) THEN
         CALL BNDSOhL (1, G, MDG, NORD, IP, IR, COEFF, N, RNORM)
         RETURN
      ENDIF
C
C     Check further for sufficient storage in working arrays.
C
      IF (IW1.LT.LW) THEN
         WRITE (XERN1, '(I8)') LW
         CALL XERMShG ('SLATEC', 'FCMNh',
     *      'IN FCh, INSUFFICIENT STORAGE FOR W(*).  CHECK LW = ' //
     *      XERN1, 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
      IF (IW2.LT.INTW1) THEN
         WRITE (XERN1, '(I8)') INTW1
         CALL XERMShG ('SLATEC', 'FCMNh',
     *      'IN FCh, INSUFFICIENT STORAGE FOR IW(*).  CHECK IW1 = ' //
     *      XERN1, 2, 1)
         MODE = -1
         RETURN
      ENDIF
C
C     Write equality constraints.
C     Analyze constraint indicators for an equality constraint.
C
      NEQCON = 0
      DO 220 IDATA = 1,NCONST
         L = NDERIV(IDATA)
         ITYPE = MOD(L,4)
         IF (ITYPE.GT.1) THEN
            IDERIV = L/4
            NEQCON = NEQCON + 1
            ILEFT = NORD
            XVAL = XCONST(IDATA)
C
  180       IF (XVAL.LT.BKPT(ILEFT+1) .OR. ILEFT.GE.N) GO TO 190
            ILEFT = ILEFT + 1
            GO TO 180
C
  190       CALL BSPLVhD (BKPT, NORD, XVAL, ILEFT, BF, IDERIV+1)
cgrs        CALL SCOPY (NP1, 0.E0, 0, W(NEQCON,1), MDW)
            do i0 = 1, np1*mdw, mdw
               w(neqcon-1+i0,1) = 0.
            enddo
            CALL SCOPY (NORD, BF(1,IDERIV+1), 1, W(NEQCON,ILEFT-NORDM1),
     +                  MDW)
C
            IF (ITYPE.EQ.2) THEN
               W(NEQCON,NP1) = YCONST(IDATA)
            ELSE
               ILEFT = NORD
               YVAL = YCONST(IDATA)
C
  200          IF (YVAL.LT.BKPT(ILEFT+1) .OR. ILEFT.GE.N) GO TO 210
               ILEFT = ILEFT + 1
               GO TO 200
C
  210          CALL BSPLVhD (BKPT, NORD, YVAL, ILEFT, BF, IDERIV+1)
               CALL SAXPY (NORD, -1.E0, BF(1, IDERIV+1), 1,
     +                     W(NEQCON, ILEFT-NORDM1), MDW)
            ENDIF
         ENDIF
  220 CONTINUE
C
C     Transfer least squares data.
C
      DO 230 I = 1,NP1
         IROW = I + NEQCON
cgrs     CALL SCOPY (N, 0.E0, 0, W(IROW,1), MDW)
         do i0 = 1, n*mdw, mdw
            w(irow-1+i0,1) = 0.
         enddo
         CALL SCOPY (MIN(NP1-I, NORD), G(I,1), MDG, W(IROW,I), MDW)
         W(IROW,NP1) = G(I,NORDP1)
  230 CONTINUE
C
C     Write inequality constraints.
C     Analyze constraint indicators for inequality constraints.
C
      NINCON = 0
      DO 260 IDATA = 1,NCONST
         L = NDERIV(IDATA)
         ITYPE = MOD(L,4)
         IF (ITYPE.LT.2) THEN
            IDERIV = L/4
            NINCON = NINCON + 1
            ILEFT = NORD
            XVAL = XCONST(IDATA)
C
  240       IF (XVAL.LT.BKPT(ILEFT+1) .OR. ILEFT.GE.N) GO TO 250
            ILEFT = ILEFT + 1
            GO TO 240
C
  250       CALL BSPLVhD (BKPT, NORD, XVAL, ILEFT, BF, IDERIV+1)
            IROW = NEQCON + NP1 + NINCON
cgrs        CALL SCOPY (N, 0.E0, 0, W(IROW,1), MDW)
            do i0 = 1, n*mdw, mdw
               w(irow-1+i0,1) = 0.
            enddo
            INTRVL = ILEFT - NORDM1
            CALL SCOPY (NORD, BF(1, IDERIV+1), 1, W(IROW, INTRVL), MDW)
C
            IF (ITYPE.EQ.1) THEN
               W(IROW,NP1) = YCONST(IDATA)
            ELSE
               W(IROW,NP1) = -YCONST(IDATA)
               CALL SSCAL (NORD, -1.E0, W(IROW, INTRVL), MDW)
            ENDIF
         ENDIF
  260 CONTINUE
C
C     Solve constrained least squares equations.
C
      CALL LSEIh(W, MDW, NEQCON, NP1, NINCON, N, PRGOPT, COEFF, RNORME,
     +          RNORML, MODE, WORK, IWORK)
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK BNDAChC
      SUBROUTINE BNDAChC (G, MDG, NB, IP, IR, MT, JT)
C***BEGIN PROLOGUE  BNDAChC
C***PURPOSE  Compute the LU factorization of a banded matrices using
C            sequential accumulation of rows of the data matrix.
C            Exactly one right-hand side vector is permitted.
C***LIBRARY   SLATEC
C***CATEGORY  D9
C***TYPE      SINGLE PRECISION (BNDAChC-S, DBNDAC-D)
C***KEYWORDS  BANDED MATRIX, CURVE FITTING, LEAST SQUARES
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     These subroutines solve the least squares problem Ax = b for
C     banded matrices A using sequential accumulation of rows of the
C     data matrix.  Exactly one right-hand side vector is permitted.
C
C     These subroutines are intended for the type of least squares
C     systems that arise in applications such as curve or surface
C     fitting of data.  The least squares equations are accumulated and
C     processed using only part of the data.  This requires a certain
C     user interaction during the solution of Ax = b.
C
C     Specifically, suppose the data matrix (A B) is row partitioned
C     into Q submatrices.  Let (E F) be the T-th one of these
C     submatrices where E = (0 C 0).  Here the dimension of E is MT by N
C     and the dimension of C is MT by NB.  The value of NB is the
C     bandwidth of A.  The dimensions of the leading block of zeros in E
C     are MT by JT-1.
C
C     The user of the subroutine BNDAChC provides MT,JT,C and F for
C     T=1,...,Q.  Not all of this data must be supplied at once.
C
C     Following the processing of the various blocks (E F), the matrix
C     (A B) has been transformed to the form (R D) where R is upper
C     triangular and banded with bandwidth NB.  The least squares
C     system Rx = d is then easily solved using back substitution by
C     executing the statement CALL BNDSOhL(1,...). The sequence of
C     values for JT must be nondecreasing.  This may require some
C     preliminary interchanges of rows and columns of the matrix A.
C
C     The primary reason for these subroutines is that the total
C     processing can take place in a working array of dimension MU by
C     NB+1.  An acceptable value for MU is
C
C                       MU = MAX(MT + N + 1),
C
C     where N is the number of unknowns.
C
C     Here the maximum is taken over all values of MT for T=1,...,Q.
C     Notice that MT can be taken to be a small as one, showing that
C     MU can be as small as N+2.  The subprogram BNDAChC processes the
C     rows more efficiently if MU is large enough so that each new
C     block (C F) has a distinct value of JT.
C
C     The four principle parts of these algorithms are obtained by the
C     following call statements
C
C     CALL BNDAChC(...)  Introduce new blocks of data.
C
C     CALL BNDSOhL(1,...)Compute solution vector and length of
C                       residual vector.
C
C     CALL BNDSOhL(2,...)Given any row vector H solve YR = H for the
C                       row vector Y.
C
C     CALL BNDSOhL(3,...)Given any column vector W solve RZ = W for
C                       the column vector Z.
C
C     The dots in the above call statements indicate additional
C     arguments that will be specified in the following paragraphs.
C
C     The user must dimension the array appearing in the call list..
C     G(MDG,NB+1)
C
C     Description of calling sequence for BNDAChC..
C
C     The entire set of parameters for BNDAChC are
C
C     Input..
C
C     G(*,*)            The working array into which the user will
C                       place the MT by NB+1 block (C F) in rows IR
C                       through IR+MT-1, columns 1 through NB+1.
C                       See descriptions of IR and MT below.
C
C     MDG               The number of rows in the working array
C                       G(*,*).  The value of MDG should be .GE. MU.
C                       The value of MU is defined in the abstract
C                       of these subprograms.
C
C     NB                The bandwidth of the data matrix A.
C
C     IP                Set by the user to the value 1 before the
C                       first call to BNDAChC.  Its subsequent value
C                       is controlled by BNDAChC to set up for the
C                       next call to BNDAChC.
C
C     IR                Index of the row of G(*,*) where the user is
C                       to place the new block of data (C F).  Set by
C                       the user to the value 1 before the first call
C                       to BNDAChC.  Its subsequent value is controlled
C                       by BNDAChC. A value of IR .GT. MDG is considered
C                       an error.
C
C     MT,JT             Set by the user to indicate respectively the
C                       number of new rows of data in the block and
C                       the index of the first nonzero column in that
C                       set of rows (E F) = (0 C 0 F) being processed.
C
C     Output..
C
C     G(*,*)            The working array which will contain the
C                       processed rows of that part of the data
C                       matrix which has been passed to BNDAChC.
C
C     IP,IR             The values of these arguments are advanced by
C                       BNDAChC to be ready for storing and processing
C                       a new block of data in G(*,*).
C
C     Description of calling sequence for BNDSOhL..
C
C     The user must dimension the arrays appearing in the call list..
C
C     G(MDG,NB+1), X(N)
C
C     The entire set of parameters for BNDSOhL are
C
C     Input..
C
C     MODE              Set by the user to one of the values 1, 2, or
C                       3.  These values respectively indicate that
C                       the solution of AX = B, YR = H or RZ = W is
C                       required.
C
C     G(*,*),MDG,       These arguments all have the same meaning and
C      NB,IP,IR         contents as following the last call to BNDAChC.
C
C     X(*)              With mode=2 or 3 this array contains,
C                       respectively, the right-side vectors H or W of
C                       the systems YR = H or RZ = W.
C
C     N                 The number of variables in the solution
C                       vector.  If any of the N diagonal terms are
C                       zero the subroutine BNDSOhL prints an
C                       appropriate message.  This condition is
C                       considered an error.
C
C     Output..
C
C     X(*)              This array contains the solution vectors X,
C                       Y or Z of the systems AX = B, YR = H or
C                       RZ = W depending on the value of MODE=1,
C                       2 or 3.
C
C     RNORM             If MODE=1 RNORM is the Euclidean length of the
C                       residual vector AX-B.  When MODE=2 or 3 RNORM
C                       is set to zero.
C
C     Remarks..
C
C     To obtain the upper triangular matrix and transformed right-hand
C     side vector D so that the super diagonals of R form the columns
C     of G(*,*), execute the following Fortran statements.
C
C     NBP1=NB+1
C
C     DO 10 J=1, NBP1
C
C  10 G(IR,J) = 0.E0
C
C     MT=1
C
C     JT=N+1
C
C     CALL BNDAChC(G,MDG,NB,IP,IR,MT,JT)
C
C***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974, Chapter 27.
C***ROUTINES CALLED  H12h, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960203  Changed name to BNDAChC with REAL --> real, etc. for mmm (TDR)
C***END PROLOGUE  BNDAChC
      implicit none
      integer mdg, nb, ip, ir, mt, jt, nerr, iopt, i, nbp1, ig1, ig2,
     .        ie, j, ig, mu, l, k, lp1, jg, kh, mh
      real G(MDG,*), zero, rho
C***FIRST EXECUTABLE STATEMENT  BNDAChC
      ZERO=0.
C
C              ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE.
C
      NBP1=NB+1
      IF (MT.LE.0.OR.NB.LE.0) RETURN
C
      IF(.NOT.MDG.LT.IR) GO TO 5
      NERR=1
      IOPT=2
      CALL XERMShG ('SLATEC', 'BNDAChC', 'MDG.LT.IR, PROBABLE ERROR.',
     +   NERR, IOPT)
      RETURN
    5 CONTINUE
C
C                                             ALG. STEP 5
      IF (JT.EQ.IP) GO TO 70
C                                             ALG. STEPS 6-7
      IF (JT.LE.IR) GO TO 30
C                                             ALG. STEPS 8-9
      DO 10 I=1,MT
        IG1=JT+MT-I
        IG2=IR+MT-I
        DO 10 J=1,NBP1
        G(IG1,J)=G(IG2,J)
   10 CONTINUE
C                                             ALG. STEP 10
      IE=JT-IR
      DO 20 I=1,IE
        IG=IR+I-1
        DO 20 J=1,NBP1
        G(IG,J)=ZERO
   20 CONTINUE
C                                             ALG. STEP 11
      IR=JT
C                                             ALG. STEP 12
   30 MU=MIN(NB-1,IR-IP-1)
      IF (MU.EQ.0) GO TO 60
C                                             ALG. STEP 13
      DO 50 L=1,MU
C                                             ALG. STEP 14
        K=MIN(L,JT-IP)
C                                             ALG. STEP 15
        LP1=L+1
        IG=IP+L
        DO 40 I=LP1,NB
          JG=I-K
          G(IG,JG)=G(IG,I)
   40 CONTINUE
C                                             ALG. STEP 16
        DO 50 I=1,K
        JG=NBP1-I
        G(IG,JG)=ZERO
   50 CONTINUE
C                                             ALG. STEP 17
   60 IP=JT
C                                             ALG. STEPS 18-19
   70 MH=IR+MT-IP
      KH=MIN(NBP1,MH)
C                                             ALG. STEP 20
      DO 80 I=1,KH
        CALL H12h (1,I,MAX(I+1,IR-IP+1),MH,G(IP,I),1,RHO,
     1            G(IP,I+1),1,MDG,NBP1-I)
   80 CONTINUE
C                                             ALG. STEP 21
      IR=IP+KH
C                                             ALG. STEP 22
      IF (KH.LT.NBP1) GO TO 100
C                                             ALG. STEP 23
      DO 90 I=1,NB
        G(IR-1,I)=ZERO
   90 CONTINUE
C                                             ALG. STEP 24
  100 CONTINUE
C                                             ALG. STEP 25
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK BNDSOhL
      SUBROUTINE BNDSOhL (MODE, G, MDG, NB, IP, IR, X, N, RNORM)
C***BEGIN PROLOGUE  BNDSOhL
C***PURPOSE  Solve the least squares problem for a banded matrix using
C            sequential accumulation of rows of the data matrix.
C            Exactly one right-hand side vector is permitted.
C***LIBRARY   SLATEC
C***CATEGORY  D9
C***TYPE      SINGLE PRECISION (BNDSOhL-S, DBNDSL-D)
C***KEYWORDS  BANDED MATRIX, CURVE FITTING, LEAST SQUARES
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     These subroutines solve the least squares problem Ax = b for
C     banded matrices A using sequential accumulation of rows of the
C     data matrix.  Exactly one right-hand side vector is permitted.
C
C     These subroutines are intended for the type of least squares
C     systems that arise in applications such as curve or surface
C     fitting of data.  The least squares equations are accumulated and
C     processed using only part of the data.  This requires a certain
C     user interaction during the solution of Ax = b.
C
C     Specifically, suppose the data matrix (A B) is row partitioned
C     into Q submatrices.  Let (E F) be the T-th one of these
C     submatrices where E = (0 C 0).  Here the dimension of E is MT by N
C     and the dimension of C is MT by NB.  The value of NB is the
C     bandwidth of A.  The dimensions of the leading block of zeros in E
C     are MT by JT-1.
C
C     The user of the subroutine BNDAChC provides MT,JT,C and F for
C     T=1,...,Q.  Not all of this data must be supplied at once.
C
C     Following the processing of the various blocks (E F), the matrix
C     (A B) has been transformed to the form (R D) where R is upper
C     triangular and banded with bandwidth NB.  The least squares
C     system Rx = d is then easily solved using back substitution by
C     executing the statement CALL BNDSOhL(1,...). The sequence of
C     values for JT must be nondecreasing.  This may require some
C     preliminary interchanges of rows and columns of the matrix A.
C
C     The primary reason for these subroutines is that the total
C     processing can take place in a working array of dimension MU by
C     NB+1.  An acceptable value for MU is
C
C                       MU = MAX(MT + N + 1),
C
C     where N is the number of unknowns.
C
C     Here the maximum is taken over all values of MT for T=1,...,Q.
C     Notice that MT can be taken to be a small as one, showing that
C     MU can be as small as N+2.  The subprogram BNDAChC processes the
C     rows more efficiently if MU is large enough so that each new
C     block (C F) has a distinct value of JT.
C
C     The four principle parts of these algorithms are obtained by the
C     following call statements
C
C     CALL BNDAChC(...)  Introduce new blocks of data.
C
C     CALL BNDSOhL(1,...)Compute solution vector and length of
C                       residual vector.
C
C     CALL BNDSOhL(2,...)Given any row vector H solve YR = H for the
C                       row vector Y.
C
C     CALL BNDSOhL(3,...)Given any column vector W solve RZ = W for
C                       the column vector Z.
C
C     The dots in the above call statements indicate additional
C     arguments that will be specified in the following paragraphs.
C
C     The user must dimension the array appearing in the call list..
C     G(MDG,NB+1)
C
C     Description of calling sequence for BNDAChC..
C
C     The entire set of parameters for BNDAChC are
C
C     Input..
C
C     G(*,*)            The working array into which the user will
C                       place the MT by NB+1 block (C F) in rows IR
C                       through IR+MT-1, columns 1 through NB+1.
C                       See descriptions of IR and MT below.
C
C     MDG               The number of rows in the working array
C                       G(*,*).  The value of MDG should be .GE. MU.
C                       The value of MU is defined in the abstract
C                       of these subprograms.
C
C     NB                The bandwidth of the data matrix A.
C
C     IP                Set by the user to the value 1 before the
C                       first call to BNDAChC.  Its subsequent value
C                       is controlled by BNDAChC to set up for the
C                       next call to BNDAChC.
C
C     IR                Index of the row of G(*,*) where the user is
C                       the user to the value 1 before the first call
C                       to BNDAChC.  Its subsequent value is controlled
C                       by BNDAChC. A value of IR .GT. MDG is considered
C                       an error.
C
C     MT,JT             Set by the user to indicate respectively the
C                       number of new rows of data in the block and
C                       the index of the first nonzero column in that
C                       set of rows (E F) = (0 C 0 F) being processed.
C     Output..
C
C     G(*,*)            The working array which will contain the
C                       processed rows of that part of the data
C                       matrix which has been passed to BNDAChC.
C
C     IP,IR             The values of these arguments are advanced by
C                       BNDAChC to be ready for storing and processing
C                       a new block of data in G(*,*).
C
C     Description of calling sequence for BNDSOhL..
C
C     The user must dimension the arrays appearing in the call list..
C
C     G(MDG,NB+1), X(N)
C
C     The entire set of parameters for BNDSOhL are
C
C     Input..
C
C     MODE              Set by the user to one of the values 1, 2, or
C                       3.  These values respectively indicate that
C                       the solution of AX = B, YR = H or RZ = W is
C                       required.
C
C     G(*,*),MDG,       These arguments all have the same meaning and
C      NB,IP,IR         contents as following the last call to BNDAChC.
C
C     X(*)              With mode=2 or 3 this array contains,
C                       respectively, the right-side vectors H or W of
C                       the systems YR = H or RZ = W.
C
C     N                 The number of variables in the solution
C                       vector.  If any of the N diagonal terms are
C                       zero the subroutine BNDSOhL prints an
C                       appropriate message.  This condition is
C                       considered an error.
C
C     Output..
C
C     X(*)              This array contains the solution vectors X,
C                       Y or Z of the systems AX = B, YR = H or
C                       RZ = W depending on the value of MODE=1,
C                       2 or 3.
C
C     RNORM             If MODE=1 RNORM is the Euclidean length of the
C                       residual vector AX-B.  When MODE=2 or 3 RNORM
C                       is set to zero.
C
C     Remarks..
C
C     To obtain the upper triangular matrix and transformed right-hand
C     side vector D so that the super diagonals of R form the columns
C     of G(*,*), execute the following Fortran statements.
C
C     NBP1=NB+1
C
C     DO 10 J=1, NBP1
C
C  10 G(IR,J) = 0.E0
C
C     MT=1
C
C     JT=N+1
C
C     CALL BNDAChC(G,MDG,NB,IP,IR,MT,JT)
C
C***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974, Chapter 27.
C***ROUTINES CALLED  XERMShG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960203  Changed name to BNDSOhL with REAL --> real, etc. for mmm (TDR)
C***END PROLOGUE  BNDSOhL
      implicit none
      integer j,mode,mdg,nb,ip,ir,n,np1,irm1,ii,l,ie,ix,i,jg,i1,i2,nerr,
     .        iopt
      real G(MDG,*),X(*),zero,rnorm,rsq,s
C***FIRST EXECUTABLE STATEMENT  BNDSOhL
      ZERO=0.
C
      RNORM=ZERO
      GO TO (10,90,50), MODE
C                                   ********************* MODE = 1
C                                   ALG. STEP 26
   10      DO 20 J=1,N
           X(J)=G(J,NB+1)
   20 CONTINUE
      RSQ=ZERO
      NP1=N+1
      IRM1=IR-1
      IF (NP1.GT.IRM1) GO TO 40
           DO 30 J=NP1,IRM1
           RSQ=RSQ+G(J,NB+1)**2
   30 CONTINUE
      RNORM=SQRT(RSQ)
   40 CONTINUE
C                                   ********************* MODE = 3
C                                   ALG. STEP 27
   50      DO 80 II=1,N
           I=N+1-II
C                                   ALG. STEP 28
           S=ZERO
           L=MAX(0,I-IP)
C                                   ALG. STEP 29
           IF (I.EQ.N) GO TO 70
C                                   ALG. STEP 30
           IE=MIN(N+1-I,NB)
                DO 60 J=2,IE
                JG=J+L
                IX=I-1+J
                S=S+G(I,JG)*X(IX)
   60 CONTINUE
C                                   ALG. STEP 31
   70      IF (G(I,L+1)) 80,130,80
   80      X(I)=(X(I)-S)/G(I,L+1)
C                                   ALG. STEP 32
      RETURN
C                                   ********************* MODE = 2
   90      DO 120 J=1,N
           S=ZERO
           IF (J.EQ.1) GO TO 110
           I1=MAX(1,J-NB+1)
           I2=J-1
                DO 100 I=I1,I2
                L=J-I+1+MAX(0,I-IP)
                S=S+X(I)*G(I,L)
  100 CONTINUE
  110      L=MAX(0,J-IP)
           IF (G(J,L+1)) 120,130,120
  120      X(J)=(X(J)-S)/G(J,L+1)
      RETURN
C
  130 CONTINUE
      NERR=1
      IOPT=2
      CALL XERMShG ('SLATEC', 'BNDSOhL',
     +   'A ZERO DIAGONAL TERM IS IN THE N BY N UPPER TRIANGULAR ' //
     +   'MATRIX.', NERR, IOPT)
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK BSPLVhD
      SUBROUTINE BSPLVhD (T, K, X, ILEFT, VNIKX, NDERIV)
C***BEGIN PROLOGUE  BSPLVhD
C***SUBSIDIARY
C***PURPOSE  Subsidiary to FCh
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BSPLVhD-S, DFSPVD-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C Calculates value and deriv.s of all B-splines which do not vanish at X
C
C  Fill VNIKX(J,IDERIV), J=IDERIV, ... ,K  with nonzero values of
C  B-splines of order K+1-IDERIV , IDERIV=NDERIV, ... ,1, by repeated
C  calls to BSPLVhN
C
C***SEE ALSO  FCh
C***ROUTINES CALLED  BSPLVhN
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
c   960203  Changed name to BSPLVhD with REAL --> real, etc. for mmm (TDR)
C***END PROLOGUE  BSPLVhD
      implicit none
      integer i,ideriv,nderiv,k,ileft,idervm,j,kmd,m,ipkmd,jm1,
     .        l,jlow
      real T(*),VNIKX(K,*),A(20,20),x,fkmd,diff,v
C***FIRST EXECUTABLE STATEMENT  BSPLVhD
      CALL BSPLVhN(T,K+1-NDERIV,1,X,ILEFT,VNIKX(NDERIV,NDERIV))
      IF (NDERIV .LE. 1)               GO TO 99
      IDERIV = NDERIV
      DO 15 I=2,NDERIV
         IDERVM = IDERIV-1
         DO 11 J=IDERIV,K
   11       VNIKX(J-1,IDERVM) = VNIKX(J,IDERIV)
         IDERIV = IDERVM
         CALL BSPLVhN(T,0,2,X,ILEFT,VNIKX(IDERIV,IDERIV))
   15    CONTINUE
C
      DO 20 I=1,K
         DO 19 J=1,K
   19       A(I,J) = 0.
   20    A(I,I) = 1.
      KMD = K
      DO 40 M=2,NDERIV
         KMD = KMD-1
         FKMD = KMD
         I = ILEFT
         J = K
   21       JM1 = J-1
            IPKMD = I + KMD
            DIFF = T(IPKMD) - T(I)
            IF (JM1 .EQ. 0)            GO TO 26
            IF (DIFF .EQ. 0.)          GO TO 25
            DO 24 L=1,J
   24          A(L,J) = (A(L,J) - A(L,J-1))/DIFF*FKMD
   25       J = JM1
            I = I - 1
                                       GO TO 21
   26    IF (DIFF .EQ. 0.)             GO TO 30
         A(1,1) = A(1,1)/DIFF*FKMD
C
   30    DO 40 I=1,K
            V = 0.
            JLOW = MAX(I,M)
            DO 35 J=JLOW,K
   35          V = A(I,J)*VNIKX(J,M) + V
   40       VNIKX(I,M) = V
   99                                  RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK BSPLVhN
      SUBROUTINE BSPLVhN (T, JHIGH, INDEX, X, ILEFT, VNIKX)
C***BEGIN PROLOGUE  BSPLVhN
C***SUBSIDIARY
C***PURPOSE  Subsidiary to FCh
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BSPLVhN-S, DFSPVN-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C Calculates the value of all possibly nonzero B-splines at *X* of
C  order MAX(JHIGH,(J+1)(INDEX-1)) on *T*.
C
C***SEE ALSO  FCh
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
c   960203  Changed name to BSPLVhN with REAL --> real, etc. for mmm (TDR)
C***END PROLOGUE  BSPLVhN
      implicit none
      integer j,i,jhigh,ipj,ileft,index,imjp1,jp1,l,jp1ml
      real T(*),VNIKX(*),DELTAM(20),DELTAP(20),x,vmprev,vm
      SAVE J, DELTAM, DELTAP
      DATA J/1/,(DELTAM(I),I=1,20),(DELTAP(I),I=1,20)/40*0./
C***FIRST EXECUTABLE STATEMENT  BSPLVhN
                                       GO TO (10,20),INDEX
   10 J = 1
      VNIKX(1) = 1.
      IF (J .GE. JHIGH)                GO TO 99
C
   20    IPJ = ILEFT+J
         DELTAP(J) = T(IPJ) - X
         IMJP1 = ILEFT-J+1
         DELTAM(J) = X - T(IMJP1)
         VMPREV = 0.
         JP1 = J+1
         DO 26 L=1,J
            JP1ML = JP1-L
            VM = VNIKX(L)/(DELTAP(L) + DELTAM(JP1ML))
            VNIKX(L) = VM*DELTAP(L) + VMPREV
   26       VMPREV = VM*DELTAM(JP1ML)
         VNIKX(JP1) = VMPREV
         J = JP1
         IF (J .LT. JHIGH)             GO TO 20
C
   99                                  RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK H12h
      SUBROUTINE H12h (MODE, LPIVOT, L1, M, U, IUE, up, C, ICE, ICV,NCV)
C***BEGIN PROLOGUE  H12h
C***SUBSIDIARY
C***PURPOSE  Subsidiary to HFTIh, LSEIh and WNNLhS
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (H12h-S, DH12-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
C     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
C
C     Construction and/or application of a single
C     Householder transformation..     Q = I + U*(U**T)/B
C
C     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
C     LPIVOT is the index of the pivot element.
C     L1,M   If L1 .LE. M   the transformation will be constructed to
C            zero elements indexed from L1 through M.   If L1 GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,up    On entry to H1 U() contains the pivot vector.
C                   IUE is the storage increment between elements.
C                                       On exit from H1 U() and up
C                   contain quantities defining the vector U of the
C                   Householder transformation.   On entry to H2 U()
C                   and up should contain quantities previously computed
C                   by H1.  These will not be modified by H2.
C     C()    On entry to H1 or H2 C() contains a matrix which will be
C            regarded as a set of vectors to which the Householder
C            transformation is to be applied.  On exit C() contains the
C            set of transformed vectors.
C     ICE    Storage increment between elements of vectors in C().
C     ICV    Storage increment between vectors in C().
C     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0
C            no operations will be done on C().
C
C***SEE ALSO  HFTIh, LSEIh, WNNhLS
C***ROUTINES CALLED  SAXPY, SDOT, SSWAP
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
c   960203  Changed name to H12h.f with REAL --> real, etc. for mmm (TDR)
C***END PROLOGUE  H12h
      implicit none
      integer lpivot,l1,m,mode,j,ncv,mml1p2,i2,icv,ice,i3,i4,l1m1,
     .        kl2,kl1,klp,iue,incr,i
      real U(IUE,*), C(*),cl,sm,one,clinv,up,b,ul1m1,sdot
C***FIRST EXECUTABLE STATEMENT  H12h
      ONE=1.
C
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL=ABS(U(1,LPIVOT))
      IF (MODE.EQ.2) GO TO 60
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M
   10     CL=MAX(ABS(U(1,J)),CL)
      IF (CL) 130,130,20
   20 CLINV=ONE/CL
      SM=(U(1,LPIVOT)*CLINV)**2
          DO 30 J=L1,M
   30     SM=SM+(U(1,J)*CLINV)**2
      CL=CL*SQRT(SM)
      IF (U(1,LPIVOT)) 50,50,40
   40 CL=-CL
   50 up=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GO TO 70
C            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN
      B=up*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C
      IF (B) 80,130,130
   80 B=ONE/B
      MML1P2=M-L1+2
      IF (MML1P2.GT.20) GO TO 140
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
          DO 120 J=1,NCV
          I2=I2+ICV
          I3=I2+INCR
          I4=I3
          SM=C(I2)*up
              DO 90 I=L1,M
              SM=SM+C(I3)*U(1,I)
   90         I3=I3+ICE
          IF (SM) 100,120,100
  100     SM=SM*B
          C(I2)=C(I2)+SM*up
              DO 110 I=L1,M
              C(I4)=C(I4)+SM*U(1,I)
  110         I4=I4+ICE
  120     CONTINUE
  130 RETURN
  140 CONTINUE
      L1M1=L1-1
      KL1=1+(L1M1-1)*ICE
      KL2=KL1
      KLP=1+(LPIVOT-1)*ICE
      UL1M1=U(1,L1M1)
      U(1,L1M1)=up
      IF (LPIVOT.EQ.L1M1) GO TO 150
      CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
  150 CONTINUE
          DO 160 J=1,NCV
          SM=SDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
          SM=SM*B
          CALL SAXPY (MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
          KL1=KL1+ICV
  160 CONTINUE
      U(1,L1M1)=UL1M1
      IF (LPIVOT.EQ.L1M1) RETURN
      KL1=KL2
      CALL SSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK HFTIh
      SUBROUTINE HFTIh (A, MDA, M, N, B, MDB, NB, TAU, KRANK, RNORM, H,
     +   G, IP)
C***BEGIN PROLOGUE  HFTIh
C***PURPOSE  Solve a linear least squares problems by performing a QR
C            factorization of the matrix using Householder
C            transformations.
C***LIBRARY   SLATEC
C***CATEGORY  D9
C***TYPE      SINGLE PRECISION (HFTIh-S, DHFTI-D)
C***KEYWORDS  CURVE FITTING, LINEAR LEAST SQUARES, QR FACTORIZATION
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N)
C
C     This subroutine solves a linear least squares problem or a set of
C     linear least squares problems having the same matrix but different
C     right-side vectors.  The problem data consists of an M by N matrix
C     A, an M by NB matrix B, and an absolute tolerance parameter TAU
C     whose usage is described below.  The NB column vectors of B
C     represent right-side vectors for NB distinct linear least squares
C     problems.
C
C     This set of problems can also be written as the matrix least
C     squares problem
C
C                       AX = B,
C
C     where X is the N by NB solution matrix.
C
C     Note that if B is the M by M identity matrix, then X will be the
C     pseudo-inverse of A.
C
C     This subroutine first transforms the augmented matrix (A B) to a
C     matrix (R C) using premultiplying Householder transformations with
C     column interchanges.  All subdiagonal elements in the matrix R are
C     zero and its diagonal elements satisfy
C
C                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)),
C
C                       I = 1,...,L-1, where
C
C                       L = MIN(M,N).
C
C     The subroutine will compute an integer, KRANK, equal to the number
C     of diagonal terms of R that exceed TAU in magnitude. Then a
C     solution of minimum Euclidean length is computed using the first
C     KRANK rows of (R C).
C
C     To be specific we suggest that the user consider an easily
C     computable matrix norm, such as, the maximum of all column sums of
C     magnitudes.
C
C     Now if the relative uncertainty of B is EPS, (norm of uncertainty/
C     norm of B), it is suggested that TAU be set approximately equal to
C     EPS*(norm of A).
C
C     The user must dimension all arrays appearing in the call list..
C     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This
C     permits the solution of a range of problems in the same array
C     space.
C
C     The entire set of parameters for HFTIh are
C
C     INPUT..
C
C     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N
C                       matrix A of the least squares problem AX = B.
C                       The first dimensioning parameter of the array
C                       A(*,*) is MDA, which must satisfy MDA.GE.M
C                       Either M.GE.N or M.LT.N is permitted.  There
C                       is no restriction on the rank of A.  The
C                       condition MDA.LT.M is considered an error.
C
C     B(*),MDB,NB       If NB = 0 the subroutine will perform the
C                       orthogonal decomposition but will make no
C                       references to the array B(*).  If NB.GT.0
C                       the array B(*) must initially contain the M by
C                       NB matrix B of the least squares problem AX =
C                       B.  If NB.GE.2 the array B(*) must be doubly
C                       subscripted with first dimensioning parameter
C                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may
C                       be either doubly or singly subscripted.  In
C                       the latter case the value of MDB is arbitrary
C                       but it should be set to some valid integer
C                       value such as MDB = M.
C
C                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N)
C                       is considered an error.
C
C     TAU               Absolute tolerance parameter provided by user
C                       for pseudorank determination.
C
C     H(*),G(*),IP(*)   Arrays of working space used by HFTIh.
C
C     OUTPUT..
C
C     A(*,*)            The contents of the array A(*,*) will be
C                       modified by the subroutine. These contents
C                       are not generally required by the user.
C
C     B(*)              On return the array B(*) will contain the N by
C                       NB solution matrix X.
C
C     KRANK             Set by the subroutine to indicate the
C                       pseudorank of A.
C
C     RNORM(*)          On return, RNORM(J) will contain the Euclidean
C                       norm of the residual vector for the problem
C                       defined by the J-th column vector of the array
C                       B(*,*) for J = 1,...,NB.
C
C     H(*),G(*)         On return these arrays respectively contain
C                       elements of the pre- and post-multiplying
C                       Householder transformations used to compute
C                       the minimum Euclidean length solution.
C
C     IP(*)             Array in which the subroutine records indices
C                       describing the permutation of column vectors.
C                       The contents of arrays H(*),G(*) and IP(*)
C                       are not generally required by the user.
C
C***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974, Chapter 14.
C***ROUTINES CALLED  H12h, EPSILON, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   901005  Replace usage of DIFF with usage of R1MACH9.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960203  Changed name to HFTIh with REAL --> real, etc. for mmm (TDR)
C   080427  Replace usage of R1MACH9 with usage of EPSILON.  (JRC)
C***END PROLOGUE  HFTIh
      implicit none
      INTEGER IP(*),k,ldiag,mda,m,nerr,iopt,n,mdb,nb,krank,j,lmax,l,
     .        kp1,jb,i,ii,ip1,jj
      real A(MDA,*),B(MDB,*),H(*),G(*),RNORM(*),releps,szero,factor,
     .     hmax,tmp,sm1,tau
C     .     hmax,tmp,sm1,tau,r1mach9
      DOUBLE PRECISION SM,DZERO
      SAVE RELEPS
      DATA RELEPS /0.E0/
C***FIRST EXECUTABLE STATEMENT  HFTIh
C      IF (RELEPS.EQ.0) RELEPS = EPSILON(4)
      IF (RELEPS.EQ.0) RELEPS = EPSILON(releps)
      SZERO=0.
      DZERO=0.D0
      FACTOR=0.001
C
      K=0
      LDIAG=MIN(M,N)
      IF (LDIAG.LE.0) GO TO 270
      IF (.NOT.MDA.LT.M) GO TO 5
      NERR=1
      IOPT=2
      CALL XERMShG ('SLATEC', 'HFTIh', 'MDA.LT.M, PROBABLE ERROR.',
     +   NERR, IOPT)
      RETURN
    5 CONTINUE
C
      IF (.NOT.(NB.GT.1.AND.MAX(M,N).GT.MDB)) GO TO 6
      NERR=2
      IOPT=2
      CALL XERMShG ('SLATEC', 'HFTIh',
     +   'MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.', NERR, IOPT)
      RETURN
    6 CONTINUE
C
          DO 80 J=1,LDIAG
          IF (J.EQ.1) GO TO 20
C
C     UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C    ..
          LMAX=J
              DO 10 L=J,N
              H(L)=H(L)-A(J-1,L)**2
              IF (H(L).GT.H(LMAX)) LMAX=L
   10         CONTINUE
          IF (FACTOR*H(LMAX) .GT. HMAX*RELEPS) GO TO 50
C
C     COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C    ..
   20     LMAX=J
              DO 40 L=J,N
              H(L)=0.
                  DO 30 I=J,M
   30             H(L)=H(L)+A(I,L)**2
              IF (H(L).GT.H(LMAX)) LMAX=L
   40         CONTINUE
          HMAX=H(LMAX)
C    ..
C     LMAX HAS BEEN DETERMINED
C
C     DO COLUMN INTERCHANGES IF NEEDED.
C    ..
   50     CONTINUE
          IP(J)=LMAX
          IF (IP(J).EQ.J) GO TO 70
              DO 60 I=1,M
              TMP=A(I,J)
              A(I,J)=A(I,LMAX)
   60         A(I,LMAX)=TMP
          H(LMAX)=H(J)
C
C     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A AND B.
C    ..
   70     CALL H12h (1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,N-J)
   80     CALL H12h (2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
C
C     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
C    ..
          DO 90 J=1,LDIAG
          IF (ABS(A(J,J)).LE.TAU) GO TO 100
   90     CONTINUE
      K=LDIAG
      GO TO 110
  100 K=J-1
  110 KP1=K+1
C
C     COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
C
      IF (NB.LE.0) GO TO 140
          DO 130 JB=1,NB
          TMP=SZERO
          IF (KP1.GT.M) GO TO 130
              DO 120 I=KP1,M
  120         TMP=TMP+B(I,JB)**2
  130     RNORM(JB)=SQRT(TMP)
  140 CONTINUE
C                                           SPECIAL FOR PSEUDORANK = 0
      IF (K.GT.0) GO TO 160
      IF (NB.LE.0) GO TO 270
          DO 150 JB=1,NB
              DO 150 I=1,N
  150         B(I,JB)=SZERO
      GO TO 270
C
C     IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
C     DECOMPOSITION OF FIRST K ROWS.
C    ..
  160 IF (K.EQ.N) GO TO 180
          DO 170 II=1,K
          I=KP1-II
  170     CALL H12h (1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  180 CONTINUE
C
C
      IF (NB.LE.0) GO TO 270
          DO 260 JB=1,NB
C
C     SOLVE THE K BY K TRIANGULAR SYSTEM.
C    ..
              DO 210 L=1,K
              SM=DZERO
              I=KP1-L
              IF (I.EQ.K) GO TO 200
              IP1=I+1
                  DO 190 J=IP1,K
  190             SM=SM+A(I,J)*DBLE(B(J,JB))
  200         SM1=SM
  210         B(I,JB)=(B(I,JB)-SM1)/A(I,I)
C
C     COMPLETE COMPUTATION OF SOLUTION VECTOR.
C    ..
          IF (K.EQ.N) GO TO 240
              DO 220 J=KP1,N
  220         B(J,JB)=SZERO
              DO 230 I=1,K
  230         CALL H12h (2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,MDB,1)
C
C      RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
C      COLUMN INTERCHANGES.
C    ..
  240         DO 250 JJ=1,LDIAG
              J=LDIAG+1-JJ
              IF (IP(J).EQ.J) GO TO 250
              L=IP(J)
              TMP=B(L,JB)
              B(L,JB)=B(J,JB)
              B(J,JB)=TMP
  250         CONTINUE
  260     CONTINUE
C    ..
C     THE SOLUTION VECTORS, X, ARE NOW
C     IN THE FIRST  N  ROWS OF THE ARRAY B(,).
C
  270 KRANK=K
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK LPDPh
      SUBROUTINE LPDPh (A, MDA, M, N1, N2, PRGOPT, X, WNORM, MODE, WS,
     +   IS)
C***BEGIN PROLOGUE  LPDPh
C***SUBSIDIARY
C***PURPOSE  Subsidiary to LSEIh
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LPDPh-S, DLPDP-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1),
C     where N=N1+N2.  This is a slight overestimate for WS(*).
C
C     Determine an N1-vector W, and
C               an N2-vector Z
C     which minimizes the Euclidean length of W
C     subject to G*W+H*Z .GE. Y.
C     This is the least projected distance problem, LPDPh.
C     The matrices G and H are of respective
C     dimensions M by N1 and M by N2.
C
C     Called by subprogram LSIh( ).
C
C     The matrix
C                (G H Y)
C
C     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
C
C     The solution (W) is returned in X(*).
C                  (Z)
C
C     The value of MODE indicates the status of
C     the computation after returning to the user.
C
C          MODE=1  The solution was successfully obtained.
C
C          MODE=2  The inequalities are inconsistent.
C
C***SEE ALSO  LSEIh
C***ROUTINES CALLED  SCOPY, SDOT, SNRM2, SSCAL, WNNLhs
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910408  Updated the AUTHOR section.  (WRB)
c   960203  Changed name to LPDPh with REAL --> real, etc. for mmm (TDR)
c   960214  Eliminated "broadcast" usage of SCOPY by using do loops (GRS)
C***END PROLOGUE  LPDPh
C
C     SUBROUTINES CALLED
C
C     WNNLhs         SOLVES A NONNEGATIVELY CONSTRAINED LINEAR LEAST
C                   SQUARES PROBLEM WITH LINEAR EQUALITY CONSTRAINTS.
C                   PART OF THIS PACKAGE.
C
C++
C     SDOT,         SUBROUTINES FROM THE BLAS PACKAGE.
C     SSCAL,SNRM2,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308.
C     SCOPY
C
      implicit none
      INTEGER IS(*),mda,m,n1,n2,mode,n,i,np1,j,iw,ix,modew,l
      integer i0
      real             A(MDA,*), PRGOPT(*), WS(*), WNORM, X(*)
      real             FAC, ONE, RNORM, SC, YNORM, ZERO
      real             SDOT, SNRM2
      SAVE ZERO, ONE, FAC
      DATA ZERO, ONE /0.E0,1.E0/, FAC /0.1E0/
C***FIRST EXECUTABLE STATEMENT  LPDPh
      N = N1 + N2
      MODE = 1
      IF (.NOT.(M.LE.0)) GO TO 20
      IF (.NOT.(N.GT.0)) GO TO 10
cgrs  X(1) = ZERO
cgrs  CALL SCOPY(N, X, 0, X, 1)
      do i0 = 1, n
         x(i0) = zero
      enddo
   10 WNORM = ZERO
      RETURN
   20 NP1 = N + 1
C
C     SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
      DO 40 I=1,M
        SC = SNRM2(N,A(I,1),MDA)
        IF (.NOT.(SC.NE.ZERO)) GO TO 30
        SC = ONE/SC
        CALL SSCAL(NP1, SC, A(I,1), MDA)
   30   CONTINUE
   40 CONTINUE
C
C     SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
      YNORM = SNRM2(M,A(1,NP1),1)
      IF (.NOT.(YNORM.NE.ZERO)) GO TO 50
      SC = ONE/YNORM
      CALL SSCAL(M, SC, A(1,NP1), 1)
C
C     SCALE COLS OF MATRIX H.
   50 J = N1 + 1
   60 IF (.NOT.(J.LE.N)) GO TO 70
      SC = SNRM2(M,A(1,J),1)
      IF (SC.NE.ZERO) SC = ONE/SC
      CALL SSCAL(M, SC, A(1,J), 1)
      X(J) = SC
      J = J + 1
      GO TO 60
   70 IF (.NOT.(N1.GT.0)) GO TO 130
C
C     COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
      IW = 0
      DO 80 I=1,M
C
C     MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
        CALL SCOPY(N2, A(I,N1+1), MDA, WS(IW+1), 1)
        IW = IW + N2
C
C     MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
        CALL SCOPY(N1, A(I,1), MDA, WS(IW+1), 1)
        IW = IW + N1
C
C     MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
        WS(IW+1) = A(I,NP1)
        IW = IW + 1
   80 CONTINUE
cgrs  WS(IW+1) = ZERO
cgrs  CALL SCOPY(N, WS(IW+1), 0, WS(IW+1), 1)
      do i0 = 1, n
         ws(iw+i0) = zero
      enddo
      IW = IW + N
      WS(IW+1) = ONE
      IW = IW + 1
C
C     SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE
C     MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
C     F = TRANSPOSE OF (0,...,0,1).
      IX = IW + 1
      IW = IW + M
C
C     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLhs( ).
      IS(1) = 0
      IS(2) = 0
      CALL WNNLhs(WS, NP1, N2, NP1-N2, M, 0, PRGOPT, WS(IX), RNORM,
     1 MODEW, IS, WS(IW+1))
C
C     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
      SC = ONE - SDOT(M,A(1,NP1),1,WS(IX),1)
      IF (.NOT.(ONE+FAC*ABS(SC).NE.ONE .AND. RNORM.GT.ZERO)) GO TO 110
      SC = ONE/SC
      DO 90 J=1,N1
        X(J) = SC*SDOT(M,A(1,J),1,WS(IX),1)
   90 CONTINUE
C
C     COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS VECTOR.
      DO 100 I=1,M
        A(I,NP1) = A(I,NP1) - SDOT(N1,A(I,1),MDA,X,1)
  100 CONTINUE
      GO TO 120
  110 MODE = 2
      RETURN
  120 CONTINUE
  130 IF (.NOT.(N2.GT.0)) GO TO 180
C
C     COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
      IW = 0
      DO 140 I=1,M
        CALL SCOPY(N2, A(I,N1+1), MDA, WS(IW+1), 1)
        IW = IW + N2
        WS(IW+1) = A(I,NP1)
        IW = IW + 1
  140 CONTINUE
cgrs  WS(IW+1) = ZERO
cgrs  CALL SCOPY(N2, WS(IW+1), 0, WS(IW+1), 1)
      do i0 = 1, n2
         ws(iw+i0) = zero
      enddo
      IW = IW + N2
      WS(IW+1) = ONE
      IW = IW + 1
      IX = IW + 1
      IW = IW + M
C
C     SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE
C     OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
C     OF (0,...,0,1)).
C
C     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLhs( ).
      IS(1) = 0
      IS(2) = 0
      CALL WNNLhs(WS, N2+1, 0, N2+1, M, 0, PRGOPT, WS(IX), RNORM, MODEW,
     1 IS, WS(IW+1))
C
C     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
      SC = ONE - SDOT(M,A(1,NP1),1,WS(IX),1)
      IF (.NOT.(ONE+FAC*ABS(SC).NE.ONE .AND. RNORM.GT.ZERO)) GO TO 160
      SC = ONE/SC
      DO 150 J=1,N2
        L = N1 + J
        X(L) = SC*SDOT(M,A(1,L),1,WS(IX),1)*X(L)
  150 CONTINUE
      GO TO 170
  160 MODE = 2
      RETURN
  170 CONTINUE
C
C     ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
  180 CALL SSCAL(N, YNORM, X, 1)
      WNORM = SNRM2(N1,X,1)
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK LSEIh
      SUBROUTINE LSEIh (W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME,RNORML,
     +   MODE, WS, IP)
C***BEGIN PROLOGUE  LSEIh
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality and inequality constraints, and optionally compute
C            a covariance matrix.
C***LIBRARY   SLATEC
C***CATEGORY  K1A2A, D9
C***TYPE      SINGLE PRECISION (LSEIh-S, DLSEI-D)
C***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
C             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
C             QUADRATIC PROGRAMMING
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem with both equality and inequality constraints, and, if the
C     user requests, obtains a covariance matrix of the solution
C     parameters.
C
C     Suppose there are given matrices E, A and G of respective
C     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
C     respective lengths ME, MA and MG.  This subroutine solves the
C     linearly constrained least squares problem
C
C                   EX = F, (E ME by N) (equations to be exactly
C                                       satisfied)
C                   AX = B, (A MA by N) (equations to be
C                                       approximately satisfied,
C                                       least squares sense)
C                   GX .GE. H,(G MG by N) (inequality constraints)
C
C     The inequalities GX .GE. H mean that every component of the
C     product GX must be .GE. the corresponding component of H.
C
C     In case the equality constraints cannot be satisfied, a
C     generalized inverse solution residual vector length is obtained
C     for F-EX.  This is the minimal length possible for F-EX.
C
C     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
C     rank of the matrix E is estimated during the computation.  We call
C     this value KRANKE.  It is an output parameter in IP(1) defined
C     below.  Using a generalized inverse solution of EX=F, a reduced
C     least squares problem with inequality constraints is obtained.
C     The tolerances used in these tests for determining the rank
C     of E and the rank of the reduced least squares problem are
C     given in Sandia Tech. Rept. SAND-78-1290.  They can be
C     modified by the user if new values are provided in
C     the option list of the array PRGOPT(*).
C
C     The user must dimension all arrays appearing in the call list..
C     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
C     where K=MAX(MA+MG,N).  This allows for a solution of a range of
C     problems in the given working space.  The dimension of WS(*)
C     given is a necessary overestimate.  Once a particular problem
C     has been run, the output parameter IP(3) gives the actual
C     dimension required for that problem.
C
C     The parameters for LSEIh( ) are
C
C     Input..
C
C     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
C     ME,MA,MG,N    first dimensioning parameter equal to MDW.
C                   For this discussion let us call M = ME+MA+MG.  Then
C                   MDW must satisfy MDW .GE. M.  The condition
C                   MDW .LT. M is an error.
C
C                   The array W(*,*) contains the matrices and vectors
C
C                                  (E  F)
C                                  (A  B)
C                                  (G  H)
C
C                   in rows and columns 1,...,M and 1,...,N+1
C                   respectively.
C
C                   The integers ME, MA, and MG are the
C                   respective matrix row dimensions
C                   of E, A and G.  Each matrix has N columns.
C
C     PRGOPT(*)    This real-valued array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case, LINK=1 and the values KEY and DATA SET
C                  are not referenced.  The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1) = LINK1 (link to first entry of next group)
C               .  PRGOPT(2) = KEY1 (key to the option change)
C               .  PRGOPT(3) = data value (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
C               .  PRGOPT(LINK1+2) = data value
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK) = 1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK .GT. NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array, a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000, an error
C                  message is printed and the subprogram returns.
C
C                  Options..
C
C                  KEY=1
C                         Compute in W(*,*) the N by N
C                  covariance matrix of the solution variables
C                  as an output parameter.  Nominally the
C                  covariance matrix will not be computed.
C                  (This requires no user input.)
C                  The data set for this option is a single value.
C                  It must be nonzero when the covariance matrix
C                  is desired.  If it is zero, the covariance
C                  matrix is not computed.  When the covariance matrix
C                  is computed, the first dimensioning parameter
C                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N).
C
C                  KEY=10
C                         Suppress scaling of the inverse of the
C                  normal matrix by the scale factor RNORM**2/
C                  MAX(1, no. of degrees of freedom).  This option
C                  only applies when the option for computing the
C                  covariance matrix (KEY=1) is used.  With KEY=1 and
C                  KEY=10 used as options the unscaled inverse of the
C                  normal matrix is returned in W(*,*).
C                  The data set for this option is a single value.
C                  When it is nonzero no scaling is done.  When it is
C                  zero scaling is done.  The nominal case is to do
C                  scaling so if option (KEY=1) is used alone, the
C                  matrix will be scaled on output.
C
C                  KEY=2
C                         Scale the nonzero columns of the
C                         entire data matrix.
C                  (E)
C                  (A)
C                  (G)
C
C                  to have length one.  The data set for this
C                  option is a single value.  It must be
C                  nonzero if unit length column scaling
C                  is desired.
C
C                  KEY=3
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  (G)
C
C                  with a user-provided diagonal matrix.
C                  The data set for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=4
C                         Change the rank determination tolerance for
C                  the equality constraint equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity SRELPR is the
C                  largest positive number such that T=1.+SRELPR
C                  satisfies T .EQ. 1.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  KEY=5
C                         Change the rank determination tolerance for
C                  the reduced least squares equations from
C                  the nominal value of SQRT(SRELPR).  This quantity can
C                  be no smaller than SRELPR, the arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The data set for this option
C                  is the new tolerance.
C
C                  For example, suppose we want to change
C                  the tolerance for the reduced least squares
C                  problem, compute the covariance matrix of
C                  the solution parameters, and provide
C                  column scaling for the data matrix.  For
C                  these options the dimension of PRGOPT(*)
C                  must be at least N+9.  The Fortran statements
C                  defining these options would be as follows:
C
C                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
C                  PRGOPT(2)=1 (covariance matrix key)
C                  PRGOPT(3)=1 (covariance matrix wanted)
C
C                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
C                  PRGOPT(5)=5 (least squares equas.  tolerance key)
C                  PRGOPT(6)=... (new value of the tolerance)
C
C                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
C                  PRGOPT(8)=3 (user-provided column scaling key)
C
C                  CALL SCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
C                    scaling factors from the user array D(*)
C                    to PRGOPT(9)-PRGOPT(N+8))
C
C                  PRGOPT(N+9)=1 (no more options to change)
C
C                  The contents of PRGOPT(*) are not modified
C                  by the subprogram.
C                  The options for WNNLhs( ) can also be included
C                  in this array.  The values of KEY recognized
C                  by WNNLhs( ) are 6, 7 and 8.  Their functions
C                  are documented in the usage instructions for
C                  subroutine WNNLhs( ).  Normally these options
C                  do not need to be modified when using LSEIh( ).
C
C     IP(1),       The amounts of working storage actually
C     IP(2)        allocated for the working arrays WS(*) and
C                  IP(*), respectively.  These quantities are
C                  compared with the actual amounts of storage
C                  needed by LSEIh( ).  Insufficient storage
C                  allocated for either WS(*) or IP(*) is an
C                  error.  This feature was included in LSEIh( )
C                  because miscalculating the storage formulas
C                  for WS(*) and IP(*) might very well lead to
C                  subtle and hard-to-find execution errors.
C
C                  The length of WS(*) must be at least
C
C                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
C
C                  where K = max(MA+MG,N)
C                  This test will not be made if IP(1).LE.0.
C
C                  The length of IP(*) must be at least
C
C                  LIP = MG+2*N+2
C                  This test will not be made if IP(2).LE.0.
C
C     Output..
C
C     X(*),RNORME,  The array X(*) contains the solution parameters
C     RNORML        if the integer output flag MODE = 0 or 1.
C                   The definition of MODE is given directly below.
C                   When MODE = 0 or 1, RNORME and RNORML
C                   respectively contain the residual vector
C                   Euclidean lengths of F - EX and B - AX.  When
C                   MODE=1 the equality constraint equations EX=F
C                   are contradictory, so RNORME .NE. 0.  The residual
C                   vector F-EX has minimal Euclidean length.  For
C                   MODE .GE. 2, none of these parameters is defined.
C
C     MODE          Integer flag that indicates the subprogram
C                   status after completion.  If MODE .GE. 2, no
C                   solution has been computed.
C
C                   MODE =
C
C                   0  Both equality and inequality constraints
C                      are compatible and have been satisfied.
C
C                   1  Equality constraints are contradictory.
C                      A generalized inverse solution of EX=F was used
C                      to minimize the residual vector length F-EX.
C                      In this sense, the solution is still meaningful.
C
C                   2  Inequality constraints are contradictory.
C
C                   3  Both equality and inequality constraints
C                      are contradictory.
C
C                   The following interpretation of
C                   MODE=1,2 or 3 must be made.  The
C                   sets consisting of all solutions
C                   of the equality constraints EX=F
C                   and all vectors satisfying GX .GE. H
C                   have no points in common.  (In
C                   particular this does not say that
C                   each individual set has no points
C                   at all, although this could be the
C                   case.)
C
C                   4  Usage error occurred.  The value
C                      of MDW is .LT. ME+MA+MG, MDW is
C                      .LT. N and a covariance matrix is
C                      requested, or the option vector
C                      PRGOPT(*) is not properly defined,
C                      or the lengths of the working arrays
C                      WS(*) and IP(*), when specified in
C                      IP(1) and IP(2) respectively, are not
C                      long enough.
C
C     W(*,*)        The array W(*,*) contains the N by N symmetric
C                   covariance matrix of the solution parameters,
C                   provided this was requested on input with
C                   the option vector PRGOPT(*) and the output
C                   flag is returned with MODE = 0 or 1.
C
C     IP(*)         The integer working array has three entries
C                   that provide rank and working array length
C                   information after completion.
C
C                      IP(1) = rank of equality constraint
C                              matrix.  Define this quantity
C                              as KRANKE.
C
C                      IP(2) = rank of reduced least squares
C                              problem.
C
C                      IP(3) = the amount of storage in the
C                              working array WS(*) that was
C                              actually used by the subprogram.
C                              The formula given above for the length
C                              of WS(*) is a necessary overestimate.
C                              If exactly the same problem matrices
C                              are used in subsequent executions,
C                              the declared dimension of WS(*) can
C                              be reduced to this output value.
C     User Designated
C     Working Arrays..
C
C     WS(*),IP(*)              These are respectively type real
C                              and type integer working arrays.
C                              Their required minimal lengths are
C                              given above.
C
C***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Report SAND77-0552, Sandia
C                 Laboratories, June 1978.
C               K. H. Haskell and R. J. Hanson, Selected algorithms for
C                 the linearly constrained least squares problem - a
C                 users guide, Report SAND78-1290, Sandia Laboratories,
C                 August 1979.
C               K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Mathematical Programming
C                 21 (1981), pp. 98-118.
C               R. J. Hanson and K. H. Haskell, Two algorithms for the
C                 linearly constrained least squares problem, ACM
C                 Transactions on Mathematical Software, September 1982.
C***ROUTINES CALLED  H12h, LSIh, EPSILON, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
C                    SSCAL, SSWAP, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and extensively revised (WRB & RWC)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900510  Convert XERRWV calls to XERMShG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960214  Eliminated "broadcast" usage of SCOPY by using do loops (GRS)
C***END PROLOGUE  LSEIh
      implicit none
      INTEGER IP(3), MA, MDW, ME, MG, MODE, N
      real             PRGOPT(*), RNORME, RNORML, W(MDW,*), WS(*), X(*)
C
C      EXTERNAL H12h, LSIh, EPSILON, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
      EXTERNAL H12h, LSIh, SASUM, SAXPY, SCOPY, SDOT, SNRM2,
     *   SSCAL, SSWAP, XERMShG
C      real             EPSILON, SASUM, SDOT, SNRM2
      real             SASUM, SDOT, SNRM2
C
      real             ENORM, FNORM, GAM, RB, RN, RNMAX, SIZE, SN,
     *   SNMAX, SRELPR, T, TAU, UJ, up, VJ, XNORM, XNRME
      INTEGER I, IMAX, J, JP1, K, KEY, KRANKE, LAST, LCHK, LINK, M,
     *   MAPKE1, MDEQC, MEND, MEP1, N1, N2, NEXT, NLINK, NOPT, NP1,
     *   NTIMES
      integer i0
      LOGICAL COV, FIRST
      CHARACTER*8 XERN1, XERN2, XERN3, XERN4
      SAVE FIRST, SRELPR
C
      DATA FIRST /.true./
C***FIRST EXECUTABLE STATEMENT  LSEIh
C
C     Set the nominal tolerance used in the code for the equality
C     constraint equations.
C
C      IF (FIRST) SRELPR = R1MACH9(4)
      IF (FIRST) SRELPR = EPSILON(ENORM)
      FIRST = .false.
      TAU = SQRT(SRELPR)
C
C     Check that enough storage was allocated in WS(*) and IP(*).
C
      MODE = 4
      IF (MIN(N,ME,MA,MG) .LT. 0) THEN
         WRITE (XERN1, '(I8)') N
         WRITE (XERN2, '(I8)') ME
         WRITE (XERN3, '(I8)') MA
         WRITE (XERN4, '(I8)') MG
         CALL XERMShG ('SLATEC','LSEIh','ALL OF THE VARIABLES N, ME,' //
     *      ' MA, MG MUST BE .GE. 0$$ENTERED ROUTINE WITH' //
     *      '$$N  = ' // XERN1 //
     *      '$$ME = ' // XERN2 //
     *      '$$MA = ' // XERN3 //
     *      '$$MG = ' // XERN4, 2, 1)
         RETURN
      ENDIF
C
      IF (IP(1).GT.0) THEN
         LCHK = 2*(ME+N) + MAX(MA+MG,N) + (MG+2)*(N+7)
         IF (IP(1).LT.LCHK) THEN
            WRITE (XERN1, '(I8)') LCHK
            CALL XERMShG ('SLATEC', 'LSEIh', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR WS(*), NEED LW = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
C
      IF (IP(2).GT.0) THEN
         LCHK = MG + 2*N + 2
         IF (IP(2).LT.LCHK) THEN
            WRITE (XERN1, '(I8)') LCHK
            CALL XERMShG ('SLATEC', 'LSEIh', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR IP(*), NEED LIP = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
C
C     Compute number of possible right multiplying Householder
C     transformations.
C
      M = ME + MA + MG
      IF (N.LE.0 .OR. M.LE.0) THEN
         MODE = 0
         RNORME = 0
         RNORML = 0
         RETURN
      ENDIF
C
      IF (MDW.LT.M) THEN
         CALL XERMShG ('SLATEC', 'LSEIh', 'MDW.LT.ME+MA+MG IS AN ERROR',
     +      2, 1)
         RETURN
      ENDIF
C
      NP1 = N + 1
      KRANKE = MIN(ME,N)
      N1 = 2*KRANKE + 1
      N2 = N1 + N
C
C     Set nominal values.
C
C     The nominal column scaling used in the code is
C     the identity scaling.
C
cgrs  CALL SCOPY (N, 1.E0, 0, WS(N1), 1)
      do i0 = 1, n
         ws(n1-1+i0) = 1.
      enddo
C
C     No covariance matrix is nominally computed.
C
      COV = .false.
C
C     Process option vector.
C     Define bound for number of options to change.
C
      NOPT = 1000
      NTIMES = 0
C
C     Define bound for positive values of LINK.
C
      NLINK = 100000
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.EQ.0 .OR. LINK.GT.NLINK) THEN
         CALL XERMShG ('SLATEC', 'LSEIh',
     +      'THE OPTION VECTOR IS UNDEFINED', 2, 1)
         RETURN
      ENDIF
C
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
            CALL XERMShG ('SLATEC', 'LSEIh',
     +         'THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 2, 1)
            RETURN
         ENDIF
C
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) THEN
            COV = PRGOPT(LAST+2) .NE. 0.E0
         ELSEIF (KEY.EQ.2 .AND. PRGOPT(LAST+2).NE.0.E0) THEN
            DO 110 J = 1,N
               T = SNRM2(M,W(1,J),1)
               IF (T.NE.0.E0) T = 1.E0/T
               WS(J+N1-1) = T
  110       CONTINUE
         ELSEIF (KEY.EQ.3) THEN
            CALL SCOPY (N, PRGOPT(LAST+2), 1, WS(N1), 1)
         ELSEIF (KEY.EQ.4) THEN
            TAU = MAX(SRELPR,PRGOPT(LAST+2))
         ENDIF
C
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
         CALL XERMShG ('SLATEC', 'LSEIh',
     +      'THE OPTION VECTOR IS UNDEFINED', 2, 1)
            RETURN
         ENDIF
C
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
      DO 120 J = 1,N
         CALL SSCAL (M, WS(N1+J-1), W(1,J), 1)
  120 CONTINUE
C
      IF (COV .AND. MDW.LT.N) THEN
         CALL XERMShG ('SLATEC', 'LSEIh',
     +      'MDW .LT. N WHEN COV MATRIX NEEDED, IS AN ERROR', 2, 1)
         RETURN
      ENDIF
C
C     Problem definition and option vector OK.
C
      MODE = 0
C
C     Compute norm of equality constraint matrix and right side.
C
      ENORM = 0.E0
      DO 130 J = 1,N
         ENORM = MAX(ENORM,SASUM(ME,W(1,J),1))
  130 CONTINUE
C
      FNORM = SASUM(ME,W(1,NP1),1)
      SNMAX = 0.E0
      RNMAX = 0.E0
      DO 150 I = 1,KRANKE
C
C        Compute maximum ratio of vector lengths. Partition is at
C        column I.
C
         DO 140 K = I,ME
            SN = SDOT(N-I+1,W(K,I),MDW,W(K,I),MDW)
            RN = SDOT(I-1,W(K,1),MDW,W(K,1),MDW)
            IF (RN.EQ.0.E0 .AND. SN.GT.SNMAX) THEN
               SNMAX = SN
               IMAX = K
            ELSEIF (K.EQ.I .OR. SN*RNMAX.GT.RN*SNMAX) THEN
               SNMAX = SN
               RNMAX = RN
               IMAX = K
            ENDIF
  140    CONTINUE
C
C        Interchange rows if necessary.
C
         IF (I.NE.IMAX) CALL SSWAP (NP1, W(I,1), MDW, W(IMAX,1), MDW)
         IF (SNMAX.GT.RNMAX*TAU**2) THEN
C
C        Eliminate elements I+1,...,N in row I.
C
            CALL H12h (1, I, I+1, N, W(I,1), MDW, WS(I), W(I+1,1), MDW,
     +                1, M-I)
         ELSE
            KRANKE = I - 1
            GO TO 160
         ENDIF
  150 CONTINUE
C
C     Save diagonal terms of lower trapezoidal matrix.
C
  160 CALL SCOPY (KRANKE, W, MDW+1, WS(KRANKE+1), 1)
C
C     Use Householder transformation from left to achieve
C     KRANKE by KRANKE upper triangular form.
C
      IF (KRANKE.LT.ME) THEN
         DO 170 K = KRANKE,1,-1
C
C           Apply transformation to matrix cols. 1,...,K-1.
C
            CALL H12h (1, K, KRANKE+1, ME, W(1,K), 1, up, W, 1, MDW,K-1)
C
C           Apply to rt side vector.
C
            CALL H12h (2, K, KRANKE+1, ME, W(1,K), 1, up, W(1,NP1), 1,1,
     +                1)
  170    CONTINUE
      ENDIF
C
C     Solve for variables 1,...,KRANKE in new coordinates.
C
      CALL SCOPY (KRANKE, W(1, NP1), 1, X, 1)
      DO 180 I = 1,KRANKE
         X(I) = (X(I)-SDOT(I-1,W(I,1),MDW,X,1))/W(I,I)
  180 CONTINUE
C
C     Compute residuals for reduced problem.
C
      MEP1 = ME + 1
      RNORML = 0.E0
      DO 190 I = MEP1,M
         W(I,NP1) = W(I,NP1) - SDOT(KRANKE,W(I,1),MDW,X,1)
         SN = SDOT(KRANKE,W(I,1),MDW,W(I,1),MDW)
         RN = SDOT(N-KRANKE,W(I,KRANKE+1),MDW,W(I,KRANKE+1),MDW)
cgrs     IF (RN.LE.SN*TAU**2 .AND. KRANKE.LT.N)
cgrs *      CALL SCOPY (N-KRANKE, 0.E0, 0, W(I,KRANKE+1), MDW)
         if (RN.LE.SN*TAU**2 .AND. KRANKE.LT.N) then
            do i0 = 1, (n-kranke)*mdw, mdw
               w(i-1+i0,kranke+1) = 0.
            enddo
         endif
  190 CONTINUE
C
C     Compute equality constraint equations residual length.
C
      RNORME = SNRM2(ME-KRANKE,W(KRANKE+1,NP1),1)
C
C     Move reduced problem data upward if KRANKE.LT.ME.
C
      IF (KRANKE.LT.ME) THEN
         DO 200 J = 1,NP1
            CALL SCOPY (M-ME, W(ME+1,J), 1, W(KRANKE+1,J), 1)
  200    CONTINUE
      ENDIF
C
C     Compute solution of reduced problem.
C
      CALL LSIh(W(KRANKE+1, KRANKE+1), MDW, MA, MG, N-KRANKE, PRGOPT,
     +         X(KRANKE+1), RNORML, MODE, WS(N2), IP(2))
C
C     Test for consistency of equality constraints.
C
      IF (ME.GT.0) THEN
         MDEQC = 0
         XNRME = SASUM(KRANKE,W(1,NP1),1)
         IF (RNORME.GT.TAU*(ENORM*XNRME+FNORM)) MDEQC = 1
         MODE = MODE + MDEQC
C
C        Check if solution to equality constraints satisfies inequality
C        constraints when there are no degrees of freedom left.
C
         IF (KRANKE.EQ.N .AND. MG.GT.0) THEN
            XNORM = SASUM(N,X,1)
            MAPKE1 = MA + KRANKE + 1
            MEND = MA + KRANKE + MG
            DO 210 I = MAPKE1,MEND
               SIZE = SASUM(N,W(I,1),MDW)*XNORM + ABS(W(I,NP1))
               IF (W(I,NP1).GT.TAU*SIZE) THEN
                  MODE = MODE + 2
                  GO TO 290
               ENDIF
  210       CONTINUE
         ENDIF
      ENDIF
C
C     Replace diagonal terms of lower trapezoidal matrix.
C
      IF (KRANKE.GT.0) THEN
         CALL SCOPY (KRANKE, WS(KRANKE+1), 1, W, MDW+1)
C
C        Reapply transformation to put solution in original coordinates.
C
         DO 220 I = KRANKE,1,-1
            CALL H12h (2, I, I+1, N, W(I,1), MDW, WS(I), X, 1, 1, 1)
  220    CONTINUE
C
C        Compute covariance matrix of equality constrained problem.
C
         IF (COV) THEN
            DO 270 J = MIN(KRANKE,N-1),1,-1
               RB = WS(J)*W(J,J)
               IF (RB.NE.0.E0) RB = 1.E0/RB
               JP1 = J + 1
               DO 230 I = JP1,N
                  W(I,J) = RB*SDOT(N-J,W(I,JP1),MDW,W(J,JP1),MDW)
  230          CONTINUE
C
               GAM = 0.5E0*RB*SDOT(N-J,W(JP1,J),1,W(J,JP1),MDW)
               CALL SAXPY (N-J, GAM, W(J,JP1), MDW, W(JP1,J), 1)
               DO 250 I = JP1,N
                  DO 240 K = I,N
                     W(I,K) = W(I,K) + W(J,I)*W(K,J) + W(I,J)*W(J,K)
                     W(K,I) = W(I,K)
  240             CONTINUE
  250          CONTINUE
               UJ = WS(J)
               VJ = GAM*UJ
               W(J,J) = UJ*VJ + UJ*VJ
               DO 260 I = JP1,N
                  W(J,I) = UJ*W(I,J) + VJ*W(J,I)
  260          CONTINUE
               CALL SCOPY (N-J, W(J, JP1), MDW, W(JP1,J), 1)
  270       CONTINUE
         ENDIF
      ENDIF
C
C     Apply the scaling to the covariance matrix.
C
      IF (COV) THEN
         DO 280 I = 1,N
            CALL SSCAL (N, WS(I+N1-1), W(I,1), MDW)
            CALL SSCAL (N, WS(I+N1-1), W(1,I), 1)
  280    CONTINUE
      ENDIF
C
C     Rescale solution vector.
C
  290 IF (MODE.LE.1) THEN
         DO 300 J = 1,N
            X(J) = X(J)*WS(N1+J-1)
  300    CONTINUE
      ENDIF
C
      IP(1) = KRANKE
      IP(3) = IP(3) + 2*KRANKE + N
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK LSIh
      SUBROUTINE LSIh (W, MDW, MA, MG, N, PRGOPT, X, RNORM, MODE, WS,IP)
C***BEGIN PROLOGUE  LSIh
C***SUBSIDIARY
C***PURPOSE  Subsidiary to LSEIh
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LSIh-S, DLSI-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to LSEIh.  The documentation for
C     LSEIh has complete usage instructions.
C
C     Solve..
C              AX = B,  A  MA by N  (least squares equations)
C     subject to..
C
C              GX.GE.H, G  MG by N  (inequality constraints)
C
C     Input..
C
C      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
C                       (G H)
C
C     MDW,MA,MG,N
C              contain (resp) var. dimension of W(*,*),
C              and matrix dimensions.
C
C     PRGOPT(*),
C              Program option vector.
C
C     OUTPUT..
C
C      X(*),RNORM
C
C              Solution vector(unless MODE=2), length of AX-B.
C
C      MODE
C              =0   Inequality constraints are compatible.
C              =2   Inequality constraints contradictory.
C
C      WS(*),
C              Working storage of dimension K+N+(MG+2)*(N+7),
C              where K=MAX(MA+MG,N).
C      IP(MG+2*N+1)
C              Integer working storage
C
C***ROUTINES CALLED  H12h, HFTIh, LPDPh, EPSILON, SASUM, SAXPY, SCOPY,
C                    SDOT, SSCAL, SSWAP
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and extensively revised (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   920422  Changed CALL to HFTIh to include variable MA.  (WRB)
c   960203  Changed name to HFTIh with REAL -> real, etc. for mmm (TDR)
c   960214  Eliminated "broadcast" usage of SCOPY by using do loops (GRS)
c   960214  Added array rnorma(1) to pass to a subroutine expecting
c           an array argument (GRS)
C***END PROLOGUE  LSIh
      implicit none
      INTEGER IP(*), MA, MDW, MG, MODE, N
      real             PRGOPT(*), RNORM, W(MDW,*), WS(*), X(*)
C
C      EXTERNAL H12h, HFTIh, LPDPh, EPSILON, SASUM, SAXPY, SCOPY, SDOT,
      EXTERNAL H12h, HFTIh, LPDPh, SASUM, SAXPY, SCOPY, SDOT,
     *   SSCAL, SSWAP
C      real             EPSILON, SASUM, SDOT
      real             SASUM, SDOT
C
      real             ANORM, FAC, GAM, RB, SRELPR, TAU, TOL, XNORM
      INTEGER I, J, K, KEY, KRANK, KRM1, KRP1, L, LAST, LINK, M, MAP1,
     *   MDLPDP, MINMAN, N1, N2, N3, NEXT, NP1
      integer i0
      real rnorma(1)
      LOGICAL COV, FIRST, SCLCOV
C
      SAVE SRELPR, FIRST
      DATA FIRST /.true./
C
C***FIRST EXECUTABLE STATEMENT  LSIh
C
C     Set the nominal tolerance used in the code.
C
C      IF (FIRST) SRELPR = R1MACH9(4)
      IF (FIRST) SRELPR = EPSILON(ANORM)
      FIRST = .false.
      TOL = SQRT(SRELPR)
C
      MODE = 0
      RNORM = 0.E0
      M = MA + MG
      NP1 = N + 1
      KRANK = 0
      IF (N.LE.0 .OR. M.LE.0) GO TO 370
C
C     To process option vector.
C
      COV = .false.
      SCLCOV = .true.
      LAST = 1
      LINK = PRGOPT(1)
C
  100 IF (LINK.GT.1) THEN
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) COV = PRGOPT(LAST+2) .NE. 0.E0
         IF (KEY.EQ.10) SCLCOV = PRGOPT(LAST+2) .EQ. 0.E0
         IF (KEY.EQ.5) TOL = MAX(SRELPR,PRGOPT(LAST+2))
         NEXT = PRGOPT(LINK)
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
C     Compute matrix norm of least squares equations.
C
      ANORM = 0.E0
      DO 110 J = 1,N
         ANORM = MAX(ANORM,SASUM(MA,W(1,J),1))
  110 CONTINUE
C
C     Set tolerance for HFTIh( ) rank test.
C
      TAU = TOL*ANORM
C
C     Compute Householder orthogonal decomposition of matrix.
C
cgrs  CALL SCOPY (N, 0.E0, 0, WS, 1)
      do i0 = 1, n
         ws(i0) = 0.
      enddo
      CALL SCOPY (MA, W(1, NP1), 1, WS, 1)
      K = MAX(M,N)
      MINMAN = MIN(MA,N)
      N1 = K + 1
      N2 = N1 + N
      CALL HFTIh (W, MDW, MA, N, WS, MA, 1, TAU, KRANK, rnorma, WS(N2),
     +           WS(N1), IP)
      RNORM = rnorma(1)
      FAC = 1.E0
      GAM = MA - KRANK
      IF (KRANK.LT.MA .AND. SCLCOV) FAC = RNORM**2/GAM
C
C     Reduce to LPDPh and solve.
C
      MAP1 = MA + 1
C
C     Compute inequality rt-hand side for LPDPh.
C
      IF (MA.LT.M) THEN
         IF (MINMAN.GT.0) THEN
            DO 120 I = MAP1,M
               W(I,NP1) = W(I,NP1) - SDOT(N,W(I,1),MDW,WS,1)
  120       CONTINUE
C
C           Apply permutations to col. of inequality constraint matrix.
C
            DO 130 I = 1,MINMAN
               CALL SSWAP (MG, W(MAP1,I), 1, W(MAP1,IP(I)), 1)
  130       CONTINUE
C
C           Apply Householder transformations to constraint matrix.
C
            IF (KRANK.GT.0 .AND. KRANK.LT.N) THEN
               DO 140 I = KRANK,1,-1
                  CALL H12h (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),
     +                      W(MAP1,1), MDW, 1, MG)
  140          CONTINUE
            ENDIF
C
C           Compute permuted inequality constraint matrix times r-inv.
C
            DO 160 I = MAP1,M
               DO 150 J = 1,KRANK
                  W(I,J) = (W(I,J)-SDOT(J-1,W(1,J),1,W(I,1),MDW))/W(J,J)
  150          CONTINUE
  160       CONTINUE
         ENDIF
C
C        Solve the reduced problem with LPDPh algorithm,
C        the least projected distance problem.
C
         CALL LPDPh(W(MAP1,1), MDW, MG, KRANK, N-KRANK, PRGOPT, X,
     +             XNORM, MDLPDP, WS(N2), IP(N+1))
C
C        Compute solution in original coordinates.
C
         IF (MDLPDP.EQ.1) THEN
            DO 170 I = KRANK,1,-1
               X(I) = (X(I)-SDOT(KRANK-I,W(I,I+1),MDW,X(I+1),1))/W(I,I)
  170       CONTINUE
C
C           Apply Householder transformation to solution vector.
C
            IF (KRANK.LT.N) THEN
               DO 180 I = 1,KRANK
                  CALL H12h (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),
     +                      X, 1, 1, 1)
  180          CONTINUE
            ENDIF
C
C           Repermute variables to their input order.
C
            IF (MINMAN.GT.0) THEN
               DO 190 I = MINMAN,1,-1
                  CALL SSWAP (1, X(I), 1, X(IP(I)), 1)
  190          CONTINUE
C
C              Variables are now in original coordinates.
C              Add solution of unconstrained problem.
C
               DO 200 I = 1,N
                  X(I) = X(I) + WS(I)
  200          CONTINUE
C
C              Compute the residual vector norm.
C
               RNORM = SQRT(RNORM**2+XNORM**2)
            ENDIF
         ELSE
            MODE = 2
         ENDIF
      ELSE
         CALL SCOPY (N, WS, 1, X, 1)
      ENDIF
C
C     Compute covariance matrix based on the orthogonal decomposition
C     from HFTIh( ).
C
      IF (.NOT.COV .OR. KRANK.LE.0) GO TO 370
      KRM1 = KRANK - 1
      KRP1 = KRANK + 1
C
C     Copy diagonal terms to working array.
C
      CALL SCOPY (KRANK, W, MDW+1, WS(N2), 1)
C
C     Reciprocate diagonal terms.
C
      DO 210 J = 1,KRANK
         W(J,J) = 1.E0/W(J,J)
  210 CONTINUE
C
C     Invert the upper triangular QR factor on itself.
C
      IF (KRANK.GT.1) THEN
         DO 230 I = 1,KRM1
            DO 220 J = I+1,KRANK
               W(I,J) = -SDOT(J-I,W(I,I),MDW,W(I,J),1)*W(J,J)
  220       CONTINUE
  230    CONTINUE
      ENDIF
C
C     Compute the inverted factor times its transpose.
C
      DO 250 I = 1,KRANK
         DO 240 J = I,KRANK
            W(I,J) = SDOT(KRANK+1-J,W(I,J),MDW,W(J,J),MDW)
  240    CONTINUE
  250 CONTINUE
C
C     Zero out lower trapezoidal part.
C     Copy upper triangular to lower triangular part.
C
      IF (KRANK.LT.N) THEN
         DO 260 J = 1,KRANK
            CALL SCOPY (J, W(1,J), 1, W(J,1), MDW)
  260    CONTINUE
C
         DO 270 I = KRP1,N
cgrs        CALL SCOPY (I, 0.E0, 0, W(I,1), MDW)
            do i0 = 1, i*mdw, mdw
               w(i-1+i0,1) = 0.
            enddo
  270    CONTINUE
C
C        Apply right side transformations to lower triangle.
C
         N3 = N2 + KRP1
         DO 330 I = 1,KRANK
            L = N1 + I
            K = N2 + I
            RB = WS(L-1)*WS(K-1)
C
C           If RB.GE.0.E0, transformation can be regarded as zero.
C
            IF (RB.LT.0.E0) THEN
               RB = 1.E0/RB
C
C              Store unscaled rank one Householder update in work array.
C
cgrs           CALL SCOPY (N, 0.E0, 0, WS(N3), 1)
               do i0 = 1, n
                  ws(n3-1+i0) = 0.
               enddo
               L = N1 + I
               K = N3 + I
               WS(K-1) = WS(L-1)
C
               DO 280 J = KRP1,N
                  WS(N3+J-1) = W(I,J)
  280          CONTINUE
C
               DO 290 J = 1,N
                  WS(J) = RB*(SDOT(J-I,W(J,I),MDW,WS(N3+I-1),1)+
     +                    SDOT(N-J+1,W(J,J),1,WS(N3+J-1),1))
  290          CONTINUE
C
               L = N3 + I
               GAM = 0.5E0*RB*SDOT(N-I+1,WS(L-1),1,WS(I),1)
               CALL SAXPY (N-I+1, GAM, WS(L-1), 1, WS(I), 1)
               DO 320 J = I,N
                  DO 300 L = 1,I-1
                     W(J,L) = W(J,L) + WS(N3+J-1)*WS(L)
  300             CONTINUE
C
                  DO 310 L = I,J
                     W(J,L) = W(J,L) + WS(J)*WS(N3+L-1)+WS(L)*WS(N3+J-1)
  310             CONTINUE
  320          CONTINUE
            ENDIF
  330    CONTINUE
C
C        Copy lower triangle to upper triangle to symmetrize the
C        covariance matrix.
C
         DO 340 I = 1,N
            CALL SCOPY (I, W(I,1), MDW, W(1,I), 1)
  340    CONTINUE
      ENDIF
C
C     Repermute rows and columns.
C
      DO 350 I = MINMAN,1,-1
         K = IP(I)
         IF (I.NE.K) THEN
            CALL SSWAP (1, W(I,I), 1, W(K,K), 1)
            CALL SSWAP (I-1, W(1,I), 1, W(1,K), 1)
            CALL SSWAP (K-I-1, W(I,I+1), MDW, W(I+1,K), 1)
            CALL SSWAP (N-K, W(I, K+1), MDW, W(K, K+1), MDW)
         ENDIF
  350 CONTINUE
C
C     Put in normalized residual sum of squares scale factor
C     and symmetrize the resulting covariance matrix.
C
      DO 360 J = 1,N
         CALL SSCAL (J, FAC, W(1,J), 1)
         CALL SCOPY (J, W(1,J), 1, W(J,1), MDW)
  360 CONTINUE
C
  370 IP(1) = KRANK
      IP(2) = N + MAX(M,N) + (MG+2)*(N+7)
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK WNLIhT
      SUBROUTINE WNLIhT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE,RNORM,
     +   IDOPE, DOPE, done)
C***BEGIN PROLOGUE  WNLIhT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNNLhs
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (WNLIhT-S, DWNLIhT-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to WNNLhs( ).
C     The documentation for WNNLhs( ) has complete usage instructions.
C
C     Note  The M by (N+1) matrix W( , ) contains the rt. hand side
C           B as the (N+1)st col.
C
C     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
C     col interchanges.
C
C***SEE ALSO  WNNLhS
C***ROUTINES CALLED  H12h, ISAMAX, SCOPY, SROThM, SROTMhG, SSCAL, SSWAP,
C                    WNLTh1, WNLTh2, WNLTh3
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   890620  Revised to make WNLTh1, WNLTh2, and WNLhT3 subroutines.(RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
c   960203  Change name to WNLIhT with REAL -> real, etc. for mmm (TDR)
c   960214  Eliminated "broadcast" usage of SCOPY by using do loops (GRS)
C***END PROLOGUE  WNLIhT
      implicit none
      INTEGER IDOPE(*), IPIVOT(*), ITYPE(*), L, M, MDW, N
      real             DOPE(*), H(*), RNORM, SCALE(*), W(MDW,*)
      LOGICAL done
C
      EXTERNAL H12h, ISAMAX, SCOPY, SROThM, SROTMhG, SSCAL, SSWAP,WNLTh1,
     *   WNLTh2, WNLTh3
      INTEGER ISAMAX
      LOGICAL WNLTh2
C
      real             ALSQ, AMAX, EANORM, FACTOR, HBAR, RN, SPARAM(5),
     *   T, TAU
      INTEGER I, I1, IMAX, IR, J, J1, JJ, JP, KRANK, L1, LB, LEND, ME,
     *   MEND, NIV, NSOLN
      integer i0
      LOGICAL INDEP, RECALC
C
C***FIRST EXECUTABLE STATEMENT  WNLIhT
      ME    = IDOPE(1)
      NSOLN = IDOPE(2)
      L1    = IDOPE(3)
C
      ALSQ   = DOPE(1)
      EANORM = DOPE(2)
      TAU    = DOPE(3)
C
      LB     = MIN(M-1,L)
      RECALC = .true.
      RNORM  = 0.E0
      KRANK  = 0
C
C     We set FACTOR=1.0 so that the heavy weight ALAMDA will be
C     included in the test for column independence.
C
      FACTOR = 1.E0
      LEND = L
      DO 180 I=1,LB
C
C        Set IR to point to the I-th row.
C
         IR = I
         MEND = M
         CALL WNLTh1 (I, LEND, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,
     +                W)
C
C        Update column SS and find pivot column.
C
         CALL WNLTh3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C        Perform column interchange.
C        Test independence of incoming column.
C
  130    IF (WNLTh2(ME, MEND, IR, FACTOR, TAU, SCALE, W(1,I))) THEN
C
C           Eliminate I-th column below diagonal using modified Givens
C           transformations applied to (A B).
C
C           When operating near the ME line, use the largest element
C           above it as the pivot.
C
            DO 160 J=M,I+1,-1
               JP = J-1
               IF (J.EQ.ME+1) THEN
                  IMAX = ME
                  AMAX = SCALE(ME)*W(ME,I)**2
                  DO 150 JP=J-1,I,-1
                     T = SCALE(JP)*W(JP,I)**2
                     IF (T.GT.AMAX) THEN
                        IMAX = JP
                        AMAX = T
                     ENDIF
  150             CONTINUE
                  JP = IMAX
               ENDIF
C
               IF (W(J,I).NE.0.E0) THEN
                  CALL SROTMhG (SCALE(JP), SCALE(J), W(JP,I), W(J,I),
     +                         SPARAM)
                  W(J,I) = 0.E0
                  CALL SROThM (N+1-I, W(JP,I+1), MDW, W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  160       CONTINUE
         ELSE IF (LEND.GT.I) THEN
C
C           Column I is dependent.  Swap with column LEND.
C           Perform column interchange,
C           and find column in remaining set with largest SS.
C
            CALL WNLTh3 (I, LEND, M, MDW, IPIVOT, H, W)
            LEND = LEND - 1
            IMAX = ISAMAX(LEND-I+1, H(I), 1) + I - 1
            HBAR = H(IMAX)
            GO TO 130
         ELSE
            KRANK = I - 1
            GO TO 190
         ENDIF
  180 CONTINUE
      KRANK = L1
C
  190 IF (KRANK.LT.ME) THEN
         FACTOR = ALSQ
         DO 200 I=KRANK+1,ME
cgrs        CALL SCOPY (L, 0.E0, 0, W(I,1), MDW)
            do i0 = 1, l*mdw, mdw
               w(i-1+i0,1) = 0.
            enddo
  200    CONTINUE
C
C        Determine the rank of the remaining equality constraint
C        equations by eliminating within the block of constrained
C        variables.  Remove any redundant constraints.
C
         RECALC = .true.
         LB = MIN(L+ME-KRANK, N)
         DO 270 I=L+1,LB
            IR = KRANK + I - L
            LEND = N
            MEND = ME
            CALL WNLTh1 (I, LEND, ME, IR, MDW, RECALC, IMAX, HBAR, H,
     +                   SCALE, W)
C
C           Update col ss and find pivot col
C
            CALL WNLTh3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C           Perform column interchange
C           Eliminate elements in the I-th col.
C
            DO 240 J=ME,IR+1,-1
               IF (W(J,I).NE.0.E0) THEN
                 CALL SROTMhG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),
     +                        SPARAM)
                  W(J,I) = 0.E0
                  CALL SROThM (N+1-I, W(J-1,I+1), MDW,W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  240       CONTINUE
C
C           I=column being eliminated.
C           Test independence of incoming column.
C           Remove any redundant or dependent equality constraints.
C
            IF (.NOT.WNLTh2(ME, MEND, IR, FACTOR,TAU,SCALE,W(1,I))) THEN
               JJ = IR
               DO 260 IR=JJ,ME
cgrs              CALL SCOPY (N, 0.E0, 0, W(IR,1), MDW)
                  do i0 = 1, n*mdw, mdw
                     w(ir-1+i0,1) = 0.
                  enddo
                  RNORM = RNORM + (SCALE(IR)*W(IR,N+1)/ALSQ)*W(IR,N+1)
                  W(IR,N+1) = 0.E0
                  SCALE(IR) = 1.E0
C
C                 Reclassify the zeroed row as a least squares equation.
C
                  ITYPE(IR) = 1
  260          CONTINUE
C
C              Reduce ME to reflect any discovered dependent equality
C              constraints.
C
               ME = JJ - 1
               GO TO 280
            ENDIF
  270    CONTINUE
      ENDIF
C
C     Try to determine the variables KRANK+1 through L1 from the
C     least squares equations.  Continue the triangularization with
C     pivot element W(ME+1,I).
C
  280 IF (KRANK.LT.L1) THEN
         RECALC = .true.
C
C        Set FACTOR=ALSQ to remove effect of heavy weight from
C        test for column independence.
C
         FACTOR = ALSQ
         DO 350 I=KRANK+1,L1
C
C           Set IR to point to the ME+1-st row.
C
            IR = ME+1
            LEND = L
            MEND = M
            CALL WNLTh1 (I, L, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,
     +                   W)
C
C           Update column SS and find pivot column.
C
            CALL WNLTh3 (I, IMAX, M, MDW, IPIVOT, H, W)
C
C           Perform column interchange.
C           Eliminate I-th column below the IR-th element.
C
            DO 320 J=M,IR+1,-1
               IF (W(J,I).NE.0.E0) THEN
                 CALL SROTMhG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),
     +                        SPARAM)
                  W(J,I) = 0.E0
                  CALL SROThM (N+1-I, W(J-1,I+1),  MDW, W(J,I+1), MDW,
     +                        SPARAM)
               ENDIF
  320       CONTINUE
C
C           Test if new pivot element is near zero.
C           If so, the column is dependent.
C           Then check row norm test to be classified as independent.
C
            T = SCALE(IR)*W(IR,I)**2
            INDEP = T .GT. (TAU*EANORM)**2
            IF (INDEP) THEN
               RN = 0.E0
               DO 340 I1=IR,M
                  DO 330 J1=I+1,N
                     RN = MAX(RN, SCALE(I1)*W(I1,J1)**2)
  330             CONTINUE
  340          CONTINUE
               INDEP = T .GT. RN*TAU**2
            ENDIF
C
C           If independent, swap the IR-th and KRANK+1-th rows to
C           maintain the triangular form.  Update the rank indicator
C           KRANK and the equality constraint pointer ME.
C
            IF (.NOT.INDEP) GO TO 360
            CALL SSWAP(N+1, W(KRANK+1,1), MDW, W(IR,1), MDW)
            CALL SSWAP(1, SCALE(KRANK+1), 1, SCALE(IR), 1)
C
C           Reclassify the least square equation as an equality
C           constraint and rescale it.
C
            ITYPE(IR) = 0
            T = SQRT(SCALE(KRANK+1))
            CALL SSCAL(N+1, T, W(KRANK+1,1), MDW)
            SCALE(KRANK+1) = ALSQ
            ME = ME+1
            KRANK = KRANK+1
  350    CONTINUE
      ENDIF
C
C     If pseudorank is less than L, apply Householder transformation.
C     from right.
C
  360 IF (KRANK.LT.L) THEN
         DO 370 J=KRANK,1,-1
            CALL H12h (1, J, KRANK+1, L, W(J,1), MDW, H(J), W, MDW, 1,
     +                J-1)
  370    CONTINUE
      ENDIF
C
      NIV = KRANK + NSOLN - L
      IF (L.EQ.N) done = .true.
C
C     End of initial triangularization.
C
      IDOPE(1) = ME
      IDOPE(2) = KRANK
      IDOPE(3) = NIV
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK WNLShM
      SUBROUTINE WNLShM (W, MDW, MME, MA, N, L, PRGOPT, X, RNORM, MODE,
     +   IPIVOT, ITYPE, WD, H, SCALE, Z, TEMP, D)
C***BEGIN PROLOGUE  WNLShM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNNLhs
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (WNLShM-S, DWNLSM-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to WNNLhs.
C     The documentation for WNNLhs has complete usage instructions.
C
C     In addition to the parameters discussed in the prologue to
C     subroutine WNNLhs, the following work arrays are used in
C     subroutine WNLShM  (they are passed through the calling
C     sequence from WNNLhs for purposes of variable dimensioning).
C     Their contents will in general be of no interest to the user.
C
C         IPIVOT(*)
C            An array of length N.  Upon completion it contains the
C         pivoting information for the cols of W(*,*).
C
C         ITYPE(*)
C            An array of length M which is used to keep track
C         of the classification of the equations.  ITYPE(I)=0
C         denotes equation I as an equality constraint.
C         ITYPE(I)=1 denotes equation I as a least squares
C         equation.
C
C         WD(*)
C            An array of length N.  Upon completion it contains the
C         dual solution vector.
C
C         H(*)
C            An array of length N.  Upon completion it contains the
C         pivot scalars of the Householder transformations performed
C         in the case KRANK.LT.L.
C
C         SCALE(*)
C            An array of length M which is used by the subroutine
C         to store the diagonal matrix of weights.
C         These are used to apply the modified Givens
C         transformations.
C
C         Z(*),TEMP(*)
C            Working arrays of length N.
C
C         D(*)
C            An array of length N that contains the
C         column scaling for the matrix (E).
C                                       (A)
C
C***SEE ALSO  WNNLhs
C***ROUTINES CALLED  H12h, ISAMAX, EPSILON, SASUM, SAXPY, SCOPY, SNRM2,
C                    SROThM, SROTMhG, SSCAL, SSWAP, WNLIhT, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and revised.  (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900328  Added TYPE section.  (WRB)
C   900510  Fixed an error message.  (RWC)
c   960214  Eliminated "broadcast" usage of SCOPY by using do loops (GRS)
C***END PROLOGUE  WNLShM
      implicit none
      INTEGER IPIVOT(*), ITYPE(*), L, MA, MDW, MME, MODE, N
      real             D(*), H(*), PRGOPT(*), RNORM, SCALE(*), TEMP(*),
     *   W(MDW,*), WD(*), X(*), Z(*)
C
C      EXTERNAL H12h, ISAMAX, EPSILON, SASUM, SAXPY, SCOPY, SNRM2,SROThM,
      EXTERNAL H12h, ISAMAX, SASUM, SAXPY, SCOPY, SNRM2,SROThM,
     *   SROTMhG, SSCAL, SSWAP, WNLIhT, XERMShG
C      real             EPSILON, SASUM, SNRM2
      real             SASUM, SNRM2
      INTEGER ISAMAX
C
      real             ALAMDA, ALPHA, ALSQ, AMAX, BLOWUP, BNORM,
     *   DOPE(3), EANORM, FAC, SM, SPARAM(5), SRELPR, T, TAU, WMAX, Z2,
     *   ZZ
      INTEGER I, IDOPE(3), IMAX, ISOL, ITEMP, ITER, ITMAX, IWMAX, J,
     *   JCON, JP, KEY, KRANK, L1, LAST, LINK, M, ME, NEXT, NIV, NLINK,
     *   NOPT, NSOLN, NTIMES
      integer i0
      LOGICAL done, FEASBL, FIRST, HITCON, POS
C
      SAVE SRELPR, FIRST
      DATA FIRST /.true./
C***FIRST EXECUTABLE STATEMENT  WNLShM
C
C     Initialize variables.
C     SRELPR is the precision for the particular machine
C     being used.  This logic avoids resetting it every entry.
C
C      IF (FIRST) SRELPR = R1MACH9(4)
      IF (FIRST) SRELPR = EPSILON(ALAMDA)
      FIRST = .false.
C
C     Set the nominal tolerance used in the code.
C
      TAU = SQRT(SRELPR)
C
      M = MA + MME
      ME = MME
      MODE = 2
C
C     To process option vector
C
      FAC = 1.E-4
C
C     Set the nominal blow up factor used in the code.
C
      BLOWUP = TAU
C
C     The nominal column scaling used in the code is
C     the identity scaling.
C
cgrs  CALL SCOPY (N, 1.E0, 0, D, 1)
      do i0 = 1, n
         d(i0) = 1.
      enddo
C
C     Define bound for number of options to change.
C
      NOPT = 1000
C
C     Define bound for positive value of LINK.
C
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.LE.0 .OR. LINK.GT.NLINK) THEN
         CALL XERMShG ('SLATEC', 'WNLShM',
     +      'WNNLhs, THE OPTION VECTOR IS UNDEFINED', 3, 1)
         RETURN
      ENDIF
C
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
         CALL XERMShG ('SLATEC', 'WNLShM',
     +      'WNNLhs, THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 3, 1)
            RETURN
         ENDIF
C
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.0.E0) THEN
            DO 110 J = 1,N
               T = SNRM2(M,W(1,J),1)
               IF (T.NE.0.E0) T = 1.E0/T
               D(J) = T
  110       CONTINUE
         ENDIF
C
         IF (KEY.EQ.7) CALL SCOPY (N, PRGOPT(LAST+2), 1, D, 1)
         IF (KEY.EQ.8) TAU = MAX(SRELPR,PRGOPT(LAST+2))
         IF (KEY.EQ.9) BLOWUP = MAX(SRELPR,PRGOPT(LAST+2))
C
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
            CALL XERMShG ('SLATEC', 'WNLShM',
     +         'WNNLhS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
            RETURN
         ENDIF
C
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
      DO 120 J = 1,N
         CALL SSCAL (M, D(J), W(1,J), 1)
  120 CONTINUE
C
C     Process option vector
C
      done = .false.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      NSOLN = L
      L1 = MIN(M,L)
C
C     Compute scale factor to apply to equality constraint equations.
C
      DO 130 J = 1,N
         WD(J) = SASUM(M,W(1,J),1)
  130 CONTINUE
C
      IMAX = ISAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = SASUM(M,W(1,N+1),1)
      ALAMDA = EANORM/(SRELPR*FAC)
C
C     Define scaling diagonal matrix for modified Givens usage and
C     classify equation types.
C
      ALSQ = ALAMDA**2
      DO 140 I = 1,M
C
C        When equation I is heavily weighted ITYPE(I)=0,
C        else ITYPE(I)=1.
C
         IF (I.LE.ME) THEN
            T = ALSQ
            ITEMP = 0
         ELSE
            T = 1.E0
            ITEMP = 1
         ENDIF
         SCALE(I) = T
         ITYPE(I) = ITEMP
  140 CONTINUE
C
C     Set the solution vector X(*) to zero and the column interchange
C     matrix to the identity.
C
cgrs  CALL SCOPY (N, 0.E0, 0, X, 1)
      do i0 = 1, n
         x(i0) = 0.
      enddo
      DO 150 I = 1,N
         IPIVOT(I) = I
  150 CONTINUE
C
C     Perform initial triangularization in the submatrix
C     corresponding to the unconstrained variables.
C     Set first L components of dual vector to zero because
C     these correspond to the unconstrained variables.
C
cgrs  CALL SCOPY (L, 0.E0, 0, WD, 1)
      do i0 = 1, l
         wd(i0) = 0.
      enddo
C
C     The arrays IDOPE(*) and DOPE(*) are used to pass
C     information to WNLIhT().  This was done to avoid
C     a long calling sequence or the use of COMMON.
C
      IDOPE(1) = ME
      IDOPE(2) = NSOLN
      IDOPE(3) = L1
C
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = TAU
      CALL WNLIhT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,
     +            IDOPE, DOPE, done)
      ME    = IDOPE(1)
      KRANK = IDOPE(2)
      NIV   = IDOPE(3)
C
C     Perform WNNLhs algorithm using the following steps.
C
C     Until(done)
C        compute search direction and feasible point
C        when (HITCON) add constraints
C        else perform multiplier test and drop a constraint
C        fin
C     Compute-Final-Solution
C
C     To compute search direction and feasible point,
C     solve the triangular system of currently non-active
C     variables and store the solution in Z(*).
C
C     To solve system
C     Copy right hand side into TEMP vector to use overwriting method.
C
  160 IF (done) GO TO 330
      ISOL = L + 1
      IF (NSOLN.GE.ISOL) THEN
         CALL SCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 170 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.E0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL SAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  170    CONTINUE
      ENDIF
C
C     Increment iteration counter and check against maximum number
C     of iterations.
C
      ITER = ITER + 1
      IF (ITER.GT.ITMAX) THEN
         MODE = 1
         done = .true.
      ENDIF
C
C     Check to see if any constraints have become active.
C     If so, calculate an interpolation factor so that all
C     active constraints are removed from the basis.
C
      ALPHA = 2.E0
      HITCON = .false.
      DO 180 J = L+1,NSOLN
         ZZ = Z(J)
         IF (ZZ.LE.0.E0) THEN
            T = X(J)/(X(J)-ZZ)
            IF (T.LT.ALPHA) THEN
               ALPHA = T
               JCON = J
            ENDIF
            HITCON = .true.
         ENDIF
  180 CONTINUE
C
C     Compute search direction and feasible point
C
      IF (HITCON) THEN
C
C        To add constraints, use computed ALPHA to interpolate between
C        last feasible solution X(*) and current unconstrained (and
C        infeasible) solution Z(*).
C
         DO 190 J = L+1,NSOLN
            X(J) = X(J) + ALPHA*(Z(J)-X(J))
  190    CONTINUE
         FEASBL = .false.
C
C        Remove column JCON and shift columns JCON+1 through N to the
C        left.  Swap column JCON into the N th position.  This achieves
C        upper Hessenberg form for the nonactive constraints and
C        leaves an upper Hessenberg matrix to retriangularize.
C
  200    DO 210 I = 1,M
            T = W(I,JCON)
            CALL SCOPY (N-JCON, W(I, JCON+1), MDW, W(I, JCON), MDW)
            W(I,N) = T
  210    CONTINUE
C
C        Update permuted index vector to reflect this shift and swap.
C
         ITEMP = IPIVOT(JCON)
         DO 220 I = JCON,N - 1
            IPIVOT(I) = IPIVOT(I+1)
  220    CONTINUE
         IPIVOT(N) = ITEMP
C
C        Similarly permute X(*) vector.
C
         CALL SCOPY (N-JCON, X(JCON+1), 1, X(JCON), 1)
         X(N) = 0.E0
         NSOLN = NSOLN - 1
         NIV = NIV - 1
C
C        Retriangularize upper Hessenberg matrix after adding
C        constraints.
C
         I = KRANK + JCON - L
         DO 230 J = JCON,NSOLN
            IF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.E0) THEN
                  CALL SROTMhG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.E0
                  CALL SROThM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) THEN
C
C              Zero IP1 to I in column J
C
               IF (W(I+1,J).NE.0.E0) THEN
                  CALL SROTMhG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.E0
                  CALL SROThM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0) THEN
               CALL SSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
               CALL SSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
               ITEMP = ITYPE(I+1)
               ITYPE(I+1) = ITYPE(I)
               ITYPE(I) = ITEMP
C
C              Swapped row was formerly a pivot element, so it will
C              be large enough to perform elimination.
C              Zero IP1 to I in column J.
C
               IF (W(I+1,J).NE.0.E0) THEN
                  CALL SROTMhG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),
     +                         SPARAM)
                  W(I+1,J) = 0.E0
                  CALL SROThM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1) THEN
               IF (SCALE(I)*W(I,J)**2/ALSQ.GT.(TAU*EANORM)**2) THEN
C
C                 Zero IP1 to I in column J
C
                  IF (W(I+1,J).NE.0.E0) THEN
                     CALL SROTMhG (SCALE(I), SCALE(I+1), W(I,J),
     +                            W(I+1,J), SPARAM)
                     W(I+1,J) = 0.E0
                     CALL SROThM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,
     +                           SPARAM)
                  ENDIF
               ELSE
                  CALL SSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
                  CALL SSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
                  ITEMP = ITYPE(I+1)
                  ITYPE(I+1) = ITYPE(I)
                  ITYPE(I) = ITEMP
                  W(I+1,J) = 0.E0
               ENDIF
            ENDIF
            I = I + 1
  230    CONTINUE
C
C        See if the remaining coefficients in the solution set are
C        feasible.  They should be because of the way ALPHA was
C        determined.  If any are infeasible, it is due to roundoff
C        error.  Any that are non-positive will be set to zero and
C        removed from the solution set.
C
         DO 240 JCON = L+1,NSOLN
            IF (X(JCON).LE.0.E0) GO TO 250
  240    CONTINUE
         FEASBL = .true.
  250    IF (.NOT.FEASBL) GO TO 200
      ELSE
C
C        To perform multiplier test and drop a constraint.
C
         CALL SCOPY (NSOLN, Z, 1, X, 1)
cgrs     IF (NSOLN.LT.N) CALL SCOPY (N-NSOLN, 0.E0, 0, X(NSOLN+1), 1)
         if (NSOLN.LT.N) then
            do i0 = 1, n-nsoln
               x(nsoln+i0) = 0.
            enddo
         endif
C
C        Reclassify least squares equations as equalities as necessary.
C
         I = NIV + 1
  260    IF (I.LE.ME) THEN
            IF (ITYPE(I).EQ.0) THEN
               I = I + 1
            ELSE
               CALL SSWAP (N+1, W(I,1), MDW, W(ME,1), MDW)
               CALL SSWAP (1, SCALE(I), 1, SCALE(ME), 1)
               ITEMP = ITYPE(I)
               ITYPE(I) = ITYPE(ME)
               ITYPE(ME) = ITEMP
               ME = ME - 1
            ENDIF
            GO TO 260
         ENDIF
C
C        Form inner product vector WD(*) of dual coefficients.
C
         DO 280 J = NSOLN+1,N
            SM = 0.E0
            DO 270 I = NSOLN+1,M
               SM = SM + SCALE(I)*W(I,J)*W(I,N+1)
  270       CONTINUE
            WD(J) = SM
  280    CONTINUE
C
C        Find J such that WD(J)=WMAX is maximum.  This determines
C        that the incoming column J will reduce the residual vector
C        and be positive.
C
  290    WMAX = 0.E0
         IWMAX = NSOLN + 1
         DO 300 J = NSOLN+1,N
            IF (WD(J).GT.WMAX) THEN
               WMAX = WD(J)
               IWMAX = J
            ENDIF
  300    CONTINUE
         IF (WMAX.LE.0.E0) GO TO 330
C
C        Set dual coefficients to zero for incoming column.
C
         WD(IWMAX) = 0.E0
C
C        WMAX .GT. 0.E0, so okay to move column IWMAX to solution set.
C        Perform transformation to retriangularize, and test for near
C        linear dependence.
C
C        Swap column IWMAX into NSOLN-th position to maintain upper
C        Hessenberg form of adjacent columns, and add new column to
C        triangular decomposition.
C
         NSOLN = NSOLN + 1
         NIV = NIV + 1
         IF (NSOLN.NE.IWMAX) THEN
            CALL SSWAP (M, W(1,NSOLN), 1, W(1,IWMAX), 1)
            WD(IWMAX) = WD(NSOLN)
            WD(NSOLN) = 0.E0
            ITEMP = IPIVOT(NSOLN)
            IPIVOT(NSOLN) = IPIVOT(IWMAX)
            IPIVOT(IWMAX) = ITEMP
         ENDIF
C
C        Reduce column NSOLN so that the matrix of nonactive constraints
C        variables is triangular.
C
         DO 320 J = M,NIV+1,-1
            JP = J - 1
C
C           When operating near the ME line, test to see if the pivot
C           element is near zero.  If so, use the largest element above
C           it as the pivot.  This is to maintain the sharp interface
C           between weighted and non-weighted rows in all cases.
C
            IF (J.EQ.ME+1) THEN
               IMAX = ME
               AMAX = SCALE(ME)*W(ME,NSOLN)**2
               DO 310 JP = J - 1,NIV,-1
                  T = SCALE(JP)*W(JP,NSOLN)**2
                  IF (T.GT.AMAX) THEN
                     IMAX = JP
                     AMAX = T
                  ENDIF
  310          CONTINUE
               JP = IMAX
            ENDIF
C
            IF (W(J,NSOLN).NE.0.E0) THEN
               CALL SROTMhG (SCALE(JP), SCALE(J), W(JP,NSOLN),
     +                      W(J,NSOLN), SPARAM)
               W(J,NSOLN) = 0.E0
               CALL SROThM (N+1-NSOLN, W(JP,NSOLN+1), MDW, W(J,NSOLN+1),
     +                     MDW, SPARAM)
            ENDIF
  320    CONTINUE
C
C        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
C        this is nonpositive or too large.  If this was true or if the
C        pivot term was zero, reject the column as dependent.
C
         IF (W(NIV,NSOLN).NE.0.E0) THEN
            ISOL = NIV
            Z2 = W(ISOL,N+1)/W(ISOL,NSOLN)
            Z(NSOLN) = Z2
            POS = Z2 .GT. 0.E0
            IF (Z2*EANORM.GE.BNORM .AND. POS) THEN
               POS = .NOT. (BLOWUP*Z2*EANORM.GE.BNORM)
            ENDIF
C
C           Try to add row ME+1 as an additional equality constraint.
C           Check size of proposed new solution component.
C           Reject it if it is too large.
C
         ELSEIF (NIV.LE.ME .AND. W(ME+1,NSOLN).NE.0.E0) THEN
            ISOL = ME + 1
            IF (POS) THEN
C
C              Swap rows ME+1 and NIV, and scale factors for these rows.
C
               CALL SSWAP (N+1, W(ME+1,1), MDW, W(NIV,1), MDW)
               CALL SSWAP (1, SCALE(ME+1), 1, SCALE(NIV), 1)
               ITEMP = ITYPE(ME+1)
               ITYPE(ME+1) = ITYPE(NIV)
               ITYPE(NIV) = ITEMP
               ME = ME + 1
            ENDIF
         ELSE
            POS = .false.
         ENDIF
C
         IF (.NOT.POS) THEN
            NSOLN = NSOLN - 1
            NIV = NIV - 1
         ENDIF
         IF (.NOT.(POS.OR.done)) GO TO 290
      ENDIF
      GO TO 160
C
C     Else perform multiplier test and drop a constraint.  To compute
C     final solution.  Solve system, store results in X(*).
C
C     Copy right hand side into TEMP vector to use overwriting method.
C
  330 ISOL = 1
      IF (NSOLN.GE.ISOL) THEN
         CALL SCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 340 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
C
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.E0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL SAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  340    CONTINUE
      ENDIF
C
C     Solve system.
C
      CALL SCOPY (NSOLN, Z, 1, X, 1)
C
C     Apply Householder transformations to X(*) if KRANK.LT.L
C
      IF (KRANK.LT.L) THEN
         DO 350 I = 1,KRANK
            CALL H12h (2, I, KRANK+1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
  350    CONTINUE
      ENDIF
C
C     Fill in trailing zeroes for constrained variables not in solution.
C
cgrs  IF (NSOLN.LT.N) CALL SCOPY (N-NSOLN, 0.E0, 0, X(NSOLN+1), 1)
      if (NSOLN.LT.N) then
         do i0 = 1, n-nsoln
            x(nsoln+i0) = 0.
         enddo
      endif
C
C     Permute solution vector to natural order.
C
      DO 380 I = 1,N
         J = I
  360    IF (IPIVOT(J).EQ.I) GO TO 370
         J = J + 1
         GO TO 360
C
  370    IPIVOT(J) = IPIVOT(I)
         IPIVOT(I) = J
         CALL SSWAP (1, X(J), 1, X(I), 1)
  380 CONTINUE
C
C     Rescale the solution using the column scaling.
C
      DO 390 J = 1,N
         X(J) = X(J)*D(J)
  390 CONTINUE
C
      DO 400 I = NSOLN+1,M
         T = W(I,N+1)
         IF (I.LE.ME) T = T/ALAMDA
         T = (SCALE(I)*T)*T
         RNORM = RNORM + T
  400 CONTINUE
C
      RNORM = SQRT(RNORM)
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK WNLTh1
      SUBROUTINE WNLTh1 (I, LEND, MEND, IR, MDW, RECALC, IMAX, HBAR, H,
     +   SCALE, W)
C***BEGIN PROLOGUE  WNLTh1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIhT
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (WNLTh1-S, DWNLT1-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     To update the column Sum Of Squares and find the pivot column.
C     The column Sum of Squares Vector will be updated at each step.
C     When numerically necessary, these values will be recomputed.
C
C***SEE ALSO  WNLIhT
C***ROUTINES CALLED  ISAMAX
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIhT and made a subroutine.  (RWC))
c   960203  Changed name to WNLTh1 with REAL -> real, etc. for mmm (TDR)
C***END PROLOGUE  WNLTh1
      implicit none
      INTEGER I, IMAX, IR, LEND, MDW, MEND
      real             H(*), HBAR, SCALE(*), W(MDW,*)
      LOGICAL RECALC
C
      EXTERNAL ISAMAX
      INTEGER ISAMAX
C
      INTEGER J, K
C
C***FIRST EXECUTABLE STATEMENT  WNLTh1
      IF (IR.NE.1 .AND. (.NOT.RECALC)) THEN
C
C        Update column SS=sum of squares.
C
         DO 10 J=I,LEND
            H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
   10    CONTINUE
C
C        Test for numerical accuracy.
C
         IMAX = ISAMAX(LEND-I+1, H(I), 1) + I - 1
         RECALC = (HBAR+1.E-3*H(IMAX)) .EQ. HBAR
      ENDIF
C
C     If required, recalculate column SS, using rows IR through MEND.
C
      IF (RECALC) THEN
         DO 30 J=I,LEND
            H(J) = 0.E0
            DO 20 K=IR,MEND
               H(J) = H(J) + SCALE(K)*W(K,J)**2
   20       CONTINUE
   30    CONTINUE
C
C        Find column with largest SS.
C
         IMAX = ISAMAX(LEND-I+1, H(I), 1) + I - 1
         HBAR = H(IMAX)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK WNLTh2
      LOGICAL FUNCTION WNLTh2 (ME, MEND, IR, FACTOR, TAU, SCALE, WIC)
C***BEGIN PROLOGUE  WNLTh2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIhT
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (WNLTh2-S, DWNLT2-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     To test independence of incoming column.
C
C     Test the column IC to determine if it is linearly independent
C     of the columns already in the basis.  In the initial tri. step,
C     we usually want the heavy weight ALAMDA to be included in the
C     test for independence.  In this case, the value of FACTOR will
C     have been set to 1.E0 before this procedure is invoked.
C     In the potentially rank deficient problem, the value of FACTOR
C     will have been set to ALSQ=ALAMDA**2 to remove the effect of the
C     heavy weight from the test for independence.
C
C     Write new column as partitioned vector
C           (A1)  number of components in solution so far = NIV
C           (A2)  M-NIV components
C     And compute  SN = inverse weighted length of A1
C                  RN = inverse weighted length of A2
C     Call the column independent when RN .GT. TAU*SN
C
C***SEE ALSO  WNILT
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLIhT and made a subroutine.  (RWC))
c   960203  Changed name to WNLTh2 with REAL -> real, etc. for mmm (TDR)
C***END PROLOGUE  WNLTh2
      implicit none
      real             FACTOR, SCALE(*), TAU, WIC(*)
      INTEGER IR, ME, MEND
C
      real            RN, SN, T
      INTEGER J
C
C***FIRST EXECUTABLE STATEMENT  WNLTh2
      SN = 0.E0
      RN = 0.E0
      DO 10 J=1,MEND
         T = SCALE(J)
         IF (J.LE.ME) T = T/FACTOR
         T = T*WIC(J)**2
C
         IF (J.LT.IR) THEN
            SN = SN + T
         ELSE
            RN = RN + T
         ENDIF
   10 CONTINUE
      WNLTh2 = RN .GT. SN*TAU**2
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK WNLTh3
      SUBROUTINE WNLTh3 (I, IMAX, M, MDW, IPIVOT, H, W)
C***BEGIN PROLOGUE  WNLTh3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to WNLIhT
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (WNLTh3-S, DWNLT3-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Perform column interchange.
C     Exchange elements of permuted index vector and perform column
C     interchanges.
C
C***SEE ALSO  WNLIhT
C***ROUTINES CALLED  SSWAP
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890620  Code extracted from WNLT and made a subroutine.  (RWC))
c   960203  Changed name to WNLTh3 with REAL -> real, etc. for mmm (TDR)
C***END PROLOGUE  WNLTh3
      implicit none
      INTEGER I, IMAX, IPIVOT(*), M, MDW
      real             H(*), W(MDW,*)
C
      EXTERNAL SSWAP
C
      real             T
      INTEGER ITEMP
C
C***FIRST EXECUTABLE STATEMENT  WNLTh3
      IF (IMAX.NE.I) THEN
         ITEMP        = IPIVOT(I)
         IPIVOT(I)    = IPIVOT(IMAX)
         IPIVOT(IMAX) = ITEMP
C
         CALL SSWAP(M, W(1,IMAX), 1, W(1,I), 1)
C
         T       = H(IMAX)
         H(IMAX) = H(I)
         H(I)    = T
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
*DECK WNNLhS
      SUBROUTINE WNNLhS (W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE,
     +   IWORK, WORK)
C***BEGIN PROLOGUE  WNNLhS
C***PURPOSE  Solve a linearly constrained least squares problem with
C            equality constraints and nonnegativity constraints on
C            selected variables.
C***LIBRARY   SLATEC
C***CATEGORY  K1A2A
C***TYPE      SINGLE PRECISION (WNNLhS-S, DWNNLS-D)
C***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
C             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
C             NONNEGATIVITY CONSTRAINTS, QUADRATIC PROGRAMMING
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, K. H., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     This subprogram solves a linearly constrained least squares
C     problem.  Suppose there are given matrices E and A of
C     respective dimensions ME by N and MA by N, and vectors F
C     and B of respective lengths ME and MA.  This subroutine
C     solves the problem
C
C               EX = F, (equations to be exactly satisfied)
C
C               AX = B, (equations to be approximately satisfied,
C                        in the least squares sense)
C
C               subject to components L+1,...,N nonnegative
C
C     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
C
C     The problem is reposed as problem WNNLhS
C
C               (WT*E)X = (WT*F)
C               (   A)    (   B), (least squares)
C               subject to components L+1,...,N nonnegative.
C
C     The subprogram chooses the heavy weight (or penalty parameter) WT.
C
C     The parameters for WNNLhS are
C
C     INPUT..
C
C     W(*,*),MDW,  The array W(*,*) is double subscripted with first
C     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
C                  discussion let us call M = ME + MA.  Then MDW
C                  must satisfy MDW.GE.M.  The condition MDW.LT.M
C                  is an error.
C
C                  The array W(*,*) contains the matrices and vectors
C
C                       (E  F)
C                       (A  B)
C
C                  in rows and columns 1,...,M and 1,...,N+1
C                  respectively.  Columns 1,...,L correspond to
C                  unconstrained variables X(1),...,X(L).  The
C                  remaining variables are constrained to be
C                  nonnegative. The condition L.LT.0 or L.GT.N is
C                  an error.
C
C     PRGOPT(*)    This real-valued array is the option vector.
C                  If the user is satisfied with the nominal
C                  subprogram features set
C
C                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
C
C                  Otherwise PRGOPT(*) is a linked list consisting of
C                  groups of data of the following form
C
C                  LINK
C                  KEY
C                  DATA SET
C
C                  The parameters LINK and KEY are each one word.
C                  The DATA SET can be comprised of several words.
C                  The number of items depends on the value of KEY.
C                  The value of LINK points to the first
C                  entry of the next group of data within
C                  PRGOPT(*).  The exception is when there are
C                  no more options to change.  In that
C                  case LINK=1 and the values KEY and DATA SET
C                  are not referenced. The general layout of
C                  PRGOPT(*) is as follows.
C
C               ...PRGOPT(1)=LINK1 (link to first entry of next group)
C               .  PRGOPT(2)=KEY1 (key to the option change)
C               .  PRGOPT(3)=DATA VALUE (data value for this change)
C               .       .
C               .       .
C               .       .
C               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
C               .                       next group)
C               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
C               .  PRGOPT(LINK1+2)=DATA VALUE
C               ...     .
C               .       .
C               .       .
C               ...PRGOPT(LINK)=1 (no more options to change)
C
C                  Values of LINK that are nonpositive are errors.
C                  A value of LINK.GT.NLINK=100000 is also an error.
C                  This helps prevent using invalid but positive
C                  values of LINK that will probably extend
C                  beyond the program limits of PRGOPT(*).
C                  Unrecognized values of KEY are ignored.  The
C                  order of the options is arbitrary and any number
C                  of options can be changed with the following
C                  restriction.  To prevent cycling in the
C                  processing of the option array a count of the
C                  number of options changed is maintained.
C                  Whenever this count exceeds NOPT=1000 an error
C                  message is printed and the subprogram returns.
C
C                  OPTIONS..
C
C                  KEY=6
C                         Scale the nonzero columns of the
C                  entire data matrix
C                  (E)
C                  (A)
C                  to have length one. The DATA SET for
C                  this option is a single value.  It must
C                  be nonzero if unit length column scaling is
C                  desired.
C
C                  KEY=7
C                         Scale columns of the entire data matrix
C                  (E)
C                  (A)
C                  with a user-provided diagonal matrix.
C                  The DATA SET for this option consists
C                  of the N diagonal scaling factors, one for
C                  each matrix column.
C
C                  KEY=8
C                         Change the rank determination tolerance from
C                  the nominal value of SQRT(SRELPR).  This quantity
C                  can be no smaller than SRELPR, The arithmetic-
C                  storage precision.  The quantity used
C                  here is internally restricted to be at
C                  least SRELPR.  The DATA SET for this option
C                  is the new tolerance.
C
C                  KEY=9
C                         Change the blow-up parameter from the
C                  nominal value of SQRT(SRELPR).  The reciprocal of
C                  this parameter is used in rejecting solution
C                  components as too large when a variable is
C                  first brought into the active set.  Too large
C                  means that the proposed component times the
C                  reciprocal of the parameter is not less than
C                  the ratio of the norms of the right-side
C                  vector and the data matrix.
C                  This parameter can be no smaller than SRELPR,
C                  the arithmetic-storage precision.
C
C                  For example, suppose we want to provide
C                  a diagonal matrix to scale the problem
C                  matrix and change the tolerance used for
C                  determining linear dependence of dropped col
C                  vectors.  For these options the dimensions of
C                  PRGOPT(*) must be at least N+6.  The FORTRAN
C                  statements defining these options would
C                  be as follows.
C
C                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
C                  PRGOPT(2)=7 (user-provided scaling key)
C
C                  CALL SCOPY(N,D,1,PRGOPT(3),1) (copy the N
C                  scaling factors from a user array called D(*)
C                  into PRGOPT(3)-PRGOPT(N+2))
C
C                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
C                  PRGOPT(N+4)=8 (linear dependence tolerance key)
C                  PRGOPT(N+5)=... (new value of the tolerance)
C
C                  PRGOPT(N+6)=1 (no more options to change)
C
C
C     IWORK(1),    The amounts of working storage actually allocated
C     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
C                  respectively.  These quantities are compared with
C                  the actual amounts of storage needed for WNNLhS( ).
C                  Insufficient storage allocated for either WORK(*)
C                  or IWORK(*) is considered an error.  This feature
C                  was included in WNNLhS( ) because miscalculating
C                  the storage formulas for WORK(*) and IWORK(*)
C                  might very well lead to subtle and hard-to-find
C                  execution errors.
C
C                  The length of WORK(*) must be at least
C
C                  LW = ME+MA+5*N
C                  This test will not be made if IWORK(1).LE.0.
C
C                  The length of IWORK(*) must be at least
C
C                  LIW = ME+MA+N
C                  This test will not be made if IWORK(2).LE.0.
C
C     OUTPUT..
C
C     X(*)         An array dimensioned at least N, which will
C                  contain the N components of the solution vector
C                  on output.
C
C     RNORM        The residual norm of the solution.  The value of
C                  RNORM contains the residual vector length of the
C                  equality constraints and least squares equations.
C
C     MODE         The value of MODE indicates the success or failure
C                  of the subprogram.
C
C                  MODE = 0  Subprogram completed successfully.
C
C                       = 1  Max. number of iterations (equal to
C                            3*(N-L)) exceeded. Nearly all problems
C                            should complete in fewer than this
C                            number of iterations. An approximate
C                            solution and its corresponding residual
C                            vector length are in X(*) and RNORM.
C
C                       = 2  Usage error occurred.  The offending
C                            condition is noted with the error
C                            processing subprogram, XERMShG( ).
C
C     User-designated
C     Working arrays..
C
C     WORK(*)      A real-valued working array of length at least
C                  M + 5*N.
C
C     IWORK(*)     An integer-valued working array of length at least
C                  M+N.
C
C***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Report SAND77-0552, Sandia
C                 Laboratories, June 1978.
C               K. H. Haskell and R. J. Hanson, Selected algorithms for
C                 the linearly constrained least squares problem - a
C                 users guide, Report SAND78-1290, Sandia Laboratories,
C                 August 1979.
C               K. H. Haskell and R. J. Hanson, An algorithm for
C                 linear least squares problems with equality and
C                 nonnegativity constraints, Mathematical Programming
C                 21 (1981), pp. 98-118.
C               R. J. Hanson and K. H. Haskell, Two algorithms for the
C                 linearly constrained least squares problem, ACM
C                 Transactions on Mathematical Software, September 1982.
C               C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974.
C***ROUTINES CALLED  WNLSM, XERMShG
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   890618  Completely restructured and revised.  (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   900510  Convert XERRWV calls to XERMShG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960203  Changed name to WNNLhS with REAL -> real,etc for mmm (TDR)
c   960216  Moved declaration of mdw before its first use (GRS)
C***END PROLOGUE  WNNLhS
      implicit none
      integer mdw
      real              PRGOPT(*), RNORM, W(MDW,*), WORK(*), X(*)
      INTEGER IWORK(*),me,ma,n,l,mode,lw,liw,l1,l2,l3,l4,l5
      CHARACTER*8 XERN1
C
C
C***FIRST EXECUTABLE STATEMENT  WNNLhS
      MODE = 0
      IF (MA+ME.LE.0 .OR. N.LE.0) RETURN
      IF (IWORK(1).GT.0) THEN
         LW = ME + MA + 5*N
         IF (IWORK(1).LT.LW) THEN
            WRITE (XERN1, '(I8)') LW
            CALL XERMShG ('SLATEC', 'WNNLhS', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR WORK(*), NEED LW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
C
      IF (IWORK(2).GT.0) THEN
         LIW = ME + MA + N
         IF (IWORK(2).LT.LIW) THEN
            WRITE (XERN1, '(I8)') LIW
            CALL XERMShG ('SLATEC', 'WNNLhS', 'INSUFFICIENT STORAGE ' //
     *         'ALLOCATED FOR IWORK(*), NEED LIW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
C
      IF (MDW.LT.ME+MA) THEN
         CALL XERMShG ('SLATEC', 'WNNLhS',
     *      'THE VALUE MDW.LT.ME+MA IS AN ERROR', 1, 1)
         MODE = 2
         RETURN
      ENDIF
C
      IF (L.LT.0 .OR. L.GT.N) THEN
         CALL XERMShG ('SLATEC', 'WNNLhS',
     *      'L.GE.0 .AND. L.LE.N IS REQUIRED', 2, 1)
         MODE = 2
         RETURN
      ENDIF
C
C     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
C     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
C     REQUIRED BY THE MAIN SUBROUTINE WNLSM( ).
C
      L1 = N + 1
      L2 = L1 + N
      L3 = L2 + ME + MA
      L4 = L3 + N
      L5 = L4 + N
C
      CALL WNLShM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK,
     *           IWORK(L1), WORK(1), WORK(L1), WORK(L2), WORK(L3),
     *           WORK(L4), WORK(L5))
      RETURN
      END
c----------------------------------------------------------------------------
c----------------------------------------------------------------------------
*DECK SSORhT
      SUBROUTINE SSORhT (X, Y, N, KFLAG)
C***BEGIN PROLOGUE  SSORhT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      SINGLE PRECISION (SSORhT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   SSORhT sorts array X and optionally makes the same interchanges in
C   array Y.  The array X may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      X - array of values to be sorted   (usually abscissas)
C      Y - array to be (optionally) carried along
C      N - number of values in array X to be sorted
C      KFLAG - control parameter
C            =  2  means sort X in increasing order and carry Y along.
C            =  1  means sort X in increasing order (ignoring Y)
C            = -1  means sort X in decreasing order (ignoring Y)
C            = -2  means sort X in decreasing order and carry Y along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMShG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891024  Changed category.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMShG.  (THJ)
C   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
c   960207  Changed name to SSORhT with REAL -> real, etc. for mmm (TDR)
c   960216  Removed INTRINSIC statement (GRS)
C***END PROLOGUE  SSORhT
      implicit none
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      real X(*), Y(*)
C     .. Local Scalars ..
      real R, T, TT, TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
      EXTERNAL XERMShG
C***FIRST EXECUTABLE STATEMENT  SSORhT
      NN = N
      IF (NN .LT. 1) THEN
         CALL XERMShG ('SLATEC', 'SSORhT',
     +      'The number of values to be sorted is not positive.', 1, 1)
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         CALL XERMShG ('SLATEC', 'SSORhT',
     +      'The sort control parameter, K, is not 2, 1, -1, or -2.', 2,
     +      1)
         RETURN
      ENDIF
C
C     Alter array X to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort X only
C
      M = 1
      I = 1
      J = NN
      R = 0.375E0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = X(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I
C
   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80
C
C     Sort X and carry Y along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I
C
  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
*DECK SROThM
      SUBROUTINE SROThM (N, SX, INCX, SY, INCY, SPARAM)
C***BEGIN PROLOGUE  SROThM
C***PURPOSE  Apply a modified Givens transformation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A8
C***TYPE      SINGLE PRECISION (SROThM-S, DROTM-D)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C   SPARAM  5-element vector. SPARAM(1) is SFLAG described below.
C           Locations 2-5 of SPARAM contain elements of the
C           transformation matrix H described below.
C
C     --Output--
C       SX  rotated vector (unchanged if N .LE. 0)
C       SY  rotated vector (unchanged if N .LE. 0)
C
C     Apply the modified Givens transformation, H, to the 2 by N matrix
C     (SX**T)
C     (SY**T) , where **T indicates transpose.  The elements of SX are
C     in SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY.
C
C     With SPARAM(1)=SFLAG, H has one of the following forms:
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C
C     See SROThMG for a description of data storage in SPARAM.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960207  Change name to SROThM with REAL -> real for mppl (TDR)
C***END PROLOGUE  SROThM
      implicit none
      integer n,incx,incy,i,nsteps,kx,ky
      real SX(*), SY(*), SPARAM(5),zero,two,sflag,sh12,sh21,sh11,sh22,
     .     w,z
      SAVE ZERO, TWO
      DATA ZERO, TWO /0.0E0, 2.0E0/
C***FIRST EXECUTABLE STATEMENT  SROThM
      SFLAG=SPARAM(1)
      IF (N.LE.0 .OR. (SFLAG+TWO.EQ.ZERO)) GO TO 140
          IF (.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
C
               NSTEPS=N*INCX
               IF (SFLAG) 50,10,30
   10          CONTINUE
               SH12=SPARAM(4)
               SH21=SPARAM(3)
                    DO 20 I = 1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W+Z*SH12
                    SY(I)=W*SH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               SH11=SPARAM(2)
               SH22=SPARAM(5)
                    DO 40 I = 1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z
                    SY(I)=-W+SH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               SH11=SPARAM(2)
               SH12=SPARAM(4)
               SH21=SPARAM(3)
               SH22=SPARAM(5)
                    DO 60 I = 1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z*SH12
                    SY(I)=W*SH21+Z*SH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF (INCX .LT. 0) KX = 1+(1-N)*INCX
          IF (INCY .LT. 0) KY = 1+(1-N)*INCY
C
          IF (SFLAG) 120,80,100
   80     CONTINUE
          SH12=SPARAM(4)
          SH21=SPARAM(3)
               DO 90 I = 1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W+Z*SH12
               SY(KY)=W*SH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          SH11=SPARAM(2)
          SH22=SPARAM(5)
               DO 110 I = 1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z
               SY(KY)=-W+SH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          SH11=SPARAM(2)
          SH12=SPARAM(4)
          SH21=SPARAM(3)
          SH22=SPARAM(5)
               DO 130 I = 1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z*SH12
               SY(KY)=W*SH21+Z*SH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
      END
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
*DECK SROTMhG
      SUBROUTINE SROTMhG (SD1, SD2, SX1, SY1, SPARAM)
C***BEGIN PROLOGUE  SROTMhG
C***PURPOSE  Construct a modified Givens transformation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      SINGLE PRECISION (SROTMhG-S, DROTMG-D)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C      SD1  single precision scalar
C      SD2  single precision scalar
C      SX1  single precision scalar
C      SY2  single precision scalar
C   SPARAM  S.P. 5-vector. SPARAM(1)=SFLAG defined below.
C           Locations 2-5 contain the rotation matrix.
C
C     --Output--
C      SD1  changed to represent the effect of the transformation
C      SD2  changed to represent the effect of the transformation
C      SX1  changed to represent the effect of the transformation
C      SY2  unchanged
C
C     Construct the modified Givens transformation matrix H which zeros
C     the second component of the 2-vector  (SQRT(SD1)*SX1,SQRT(SD2)*
C     SY2)**T.
C     With SPARAM(1)=SFLAG, H has one of the following forms:
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C
C     Locations 2-5 of SPARAM contain SH11, SH21, SH12, and SH22,
C     respectively.  (Values of 1.E0, -1.E0, or 0.E0 implied by the
C     value of SPARAM(1) are not stored in SPARAM.)
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780301  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920316  Prologue corrected.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
c   960207  Changed name to SROTMhG with REAL -> real, etc for mppl (TDR)
C***END PROLOGUE  SROTMhG
      implicit none
      integer igo
      real SPARAM(5),sd1,sd2,sx1,sy1,zero,one,two,gam,gamsq,rgamsq,sq1,
     .     sq2,sp1,sp2,sflag,stemp,sh11,sh21,sh12,sh22,su
      SAVE ZERO, ONE, TWO, GAM, GAMSQ, RGAMSQ
      DATA ZERO, ONE, TWO /0.0E0, 1.0E0, 2.0E0/
      DATA GAM, GAMSQ, RGAMSQ /4096.0E0, 1.67772E7, 5.96046E-8/
C***FIRST EXECUTABLE STATEMENT  SROTMhG
      IF (.NOT. SD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-SX1..
          GO TO 60
   10 CONTINUE
C     CASE-SD1-NONNEGATIVE
      SP2=SD2*SY1
      IF (.NOT. SP2 .EQ. ZERO) GO TO 20
          SFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      SP1=SD1*SX1
      SQ2=SP2*SY1
      SQ1=SP1*SX1
C
      IF (.NOT. ABS(SQ1) .GT. ABS(SQ2)) GO TO 40
          SH21=-SY1/SX1
          SH12=SP2/SP1
C
          SU=ONE-SH12*SH21
C
          IF (.NOT. SU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   30     CONTINUE
               SFLAG=ZERO
               SD1=SD1/SU
               SD2=SD2/SU
               SX1=SX1*SU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF (.NOT. SQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   50     CONTINUE
               SFLAG=ONE
               SH11=SP1/SP2
               SH22=SX1/SY1
               SU=ONE+SH11*SH22
               STEMP=SD2/SU
               SD2=SD1/SU
               SD1=STEMP
               SX1=SY1*SU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-SX1..
   60 CONTINUE
          SFLAG=-ONE
          SH11=ZERO
          SH12=ZERO
          SH21=ZERO
          SH22=ZERO
C
          SD1=ZERO
          SD2=ZERO
          SX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF (.NOT. SFLAG .GE. ZERO) GO TO 90
C
          IF (.NOT. SFLAG .EQ. ZERO) GO TO 80
          SH11=ONE
          SH22=ONE
          SFLAG=-ONE
          GO TO 90
   80     CONTINUE
          SH21=-ONE
          SH12=ONE
          SFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF (.NOT. SD1 .LE. RGAMSQ) GO TO 130
               IF (SD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               SD1=SD1*GAM**2
               SX1=SX1/GAM
               SH11=SH11/GAM
               SH12=SH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF (.NOT. SD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               SD1=SD1/GAM**2
               SX1=SX1*GAM
               SH11=SH11*GAM
               SH12=SH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF (.NOT. ABS(SD2) .LE. RGAMSQ) GO TO 190
               IF (SD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               SD2=SD2*GAM**2
               SH21=SH21/GAM
               SH22=SH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF (.NOT. ABS(SD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               SD2=SD2/GAM**2
               SH21=SH21*GAM
               SH22=SH22*GAM
          GO TO 200
  220 CONTINUE
          IF (SFLAG) 250,230,240
  230     CONTINUE
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               GO TO 260
  240     CONTINUE
               SPARAM(2)=SH11
               SPARAM(5)=SH22
               GO TO 260
  250     CONTINUE
               SPARAM(2)=SH11
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               SPARAM(5)=SH22
  260 CONTINUE
          SPARAM(1)=SFLAG
          RETURN
      END

c*************************************************************************
c  Subroutine tanh_multi is from Osborne/Groebner for fitting DIII-D data
c*************************************************************************
      subroutine tanh_multi (ncoef, coeff, nrad, x, cfit_param, result)

      implicit none

      integer ncoef, nrad, ii, ir
      real coeff, x, result, c, rcfit_param
      dimension coeff(*),x(*), result(nrad), c(0:ncoef-1)
      character*8 cfit_param
cc      dimension cfit_param(*)

CC   If cfit_param is anything but 'none', we must convert it into a real
CC   called rfit_param and use it.
      real rfit_param, z, pz1, pz2, z0, cder
      dimension z(nrad), pz1(nrad), pz2(nrad)
      integer num

CC   If cfit_param is anything but 'none', we must convert it into a real
CC   and use it.

      do ii=1, ncoef
        c(ii-1) = coeff(ii)
      enddo

      do ir = 1, nrad
        z(ir) = 2.*(c(0)-x(ir))/c(1)
      enddo

      if (cfit_param == "none") then
        do ir = 1, nrad
          pz1(ir) = 1.+ c(4)*z(ir) + c(5)*z(ir)**2 + c(6)*z(ir)**3
        enddo
      else     ##  Convert the cfit_param string into a real number
        read (cfit_param, fmt=*) rcfit_param
        z0 = 2.*( c(0) - rfit_param )/c(1)
        cder = -( 2.0*c(4)*z0 + 3.0*c(5)*z0*z0 + 4.0*c(6)*z0*z0*z0 )
        do ir = 1, nrad
          pz1(ir) = 1 + cder*z(ir) + c(4)*z(ir)**2 +
     .                  c(5)*z(ir)**3 + c(6)*z(ir)**4
        enddo
      endif

CC  Here use pz2 = 1 + ( c(7)*z + (c(8)*z*z) ) depending if 7,8,or 9 coeff.
      num = ncoef    # size(c)

      do ir = 1, nrad
        pz2(ir) = 1.0
        if (num == 8) pz2(ir) = 1.0 + c(7)*z(ir)
        if (num == 9) pz2(ir) = 1.0 + c(7)*z(ir) + c(8)*z(ir)**2

        result(ir) =  0.5*(c(2)-c(3))* ( pz1(ir)*exp(z(ir)) -
     .                                       pz2(ir)*exp(-z(ir) ) )/
     .                   ( exp(z(ir)) + exp(-z(ir)) ) + 0.5*(c(2)+c(3))
      enddo

      return
      end
c*****************************
c  End of subroutine tanh_multi
c*****************************
c----------------------------------------------------------------------------

c******************************************************************************
c  Subroutine to read ne_tanh data and evaluate for input locations
c******************************************************************************

      subroutine readne_dat(fname)
       
      implicit none
      character*(*) fname
      Use(Dim)
      Use(Fitdata)
     
c     local variables --
      integer ios, n, nget, ii
      character*30 str1, str2, ptnam, fit_nam, dumm
      character*20 w1,w2,w3,w4,w5,w6,w7
      character*7 str3
      character*2 str4
      character*10 yunits
      character*80 line

c     procedures --
      external freeus, gallot, readrt1, xerrab

c----------------------------------------------------------------------c
c     Read tanh fit data from file (T. Osborne)
c----------------------------------------------------------------------c

ccc      call freeus(nget)
      data nget /55/
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call xerrab('**** netanh.dat file not found')
      endif

      read (nget,*) str1, str2, ptnam
      read (nget,'(a7,a2,a10)') str1, str2, yunits
      read (nget,*) str1, str2, fit_nam

      read (nget,*) str1, str2, fit_paramne_tanh
      read(nget,*) str1, str2, ncoefne_tanh
      read(nget,*) dumm

      call gchange("Fitdata",0)         # allocate fcoeff(numc)
      if (isdatfmtnew == 0) then  # old format with no text
         do ii = 1, ncoefne_tanh
           read(nget,*) fcoefne_tanh(ii)
         enddo
      elseif(isdatfmtnew == 1)
         read(nget,*) w1,w2,fcoefne_tanh(1)
         read(nget,*) w1,w2,fcoefne_tanh(2)
         read(nget,*) w1,fcoefne_tanh(3)
         read(nget,*) w1,fcoefne_tanh(4)
         read(nget,*) w1,w2,w3,w4,w5,w6,w7,fcoefne_tanh(5)
         read(nget,*) w1,w2,w3,w4,w5,w6,w7,fcoefne_tanh(6)
         read(nget,*) w1,w2,w3,w4,w5,w6,w7,fcoefne_tanh(7)
         if(ncoefne_tanh >= 8) read(nget,*) w1,w2,fcoefne_tanh(8)
         if(ncoefne_tanh == 9) read(nget,*) w1,w2,fcoefne_tanh(9)
      elseif(isdatfmtnew == 2)
         read(nget,*) fcoefne_tanh(1),w1,w2
         read(nget,*) fcoefne_tanh(2),w1,w2
         read(nget,*) fcoefne_tanh(3),w1
         read(nget,*) fcoefne_tanh(4),w1
         read(nget,*) fcoefne_tanh(5),w1,w2,w3,w4,w5,w6,w7
         read(nget,*) fcoefne_tanh(6),w1,w2,w3,w4,w5,w6,w7
         read(nget,*) fcoefne_tanh(7),w1,w2,w3,w4,w5,w6,w7
         if(ncoefne_tanh >= 8) read(nget,*) fcoefne_tanh(8),w1,w2
         if(ncoefne_tanh == 9) read(nget,*) fcoefne_tanh(9),w1,w2
      endif
      close (nget)

      return
      end
c-----------------------------------------------------------------------
c******************************************************************************
c  Subroutine to read te_tanh data and evaluate for input locations
c******************************************************************************

      subroutine readte_dat(fname)
       
      implicit none
      character*(*) fname
      Use(Dim)
      Use(Fitdata)

c     local variables --
      integer ios, n, nget, ii
      character*30 str1, str2, ptnam, fit_nam, dumm
      character*20 w1,w2,w3,w4,w5,w6,w7
      character*7 str3
      character*2 str4
      character*10 yunits
      character*80 line

c     procedures --
      external freeus, gallot, readrt1, xerrab

c----------------------------------------------------------------------c
c     Read tanh fit data from file (T. Osborne)
c----------------------------------------------------------------------c

ccc      call freeus(nget)
      data nget /55/
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call xerrab('**** tetanh.dat file not found')
      endif

      read (nget,*) str1, str2, ptnam
      read (nget,'(a7,a2,a10)') str1, str2, yunits
      read (nget,*) str1, str2, fit_nam

      read (nget,*) str1, str2, fit_paramte_tanh
      read(nget,*) str1, str2, ncoefte_tanh
      read(nget,*) dumm

      call gchange("Fitdata",0)         # allocate fcoeff(numc)
      if (isdatfmtnew == 0) then  # old format with no text
         do ii = 1, ncoefte_tanh
           read(nget,*) fcoefte_tanh(ii)
         enddo
      elseif(isdatfmtnew == 1)
         read(nget,*) w1,w2,fcoefte_tanh(1)
         read(nget,*) w1,w2,fcoefte_tanh(2)
         read(nget,*) w1,fcoefte_tanh(3)
         read(nget,*) w1,fcoefte_tanh(4)
         read(nget,*) w1,w2,w3,w4,w5,w6,w7,fcoefte_tanh(5)
         read(nget,*) w1,w2,w3,w4,w5,w6,w7,fcoefte_tanh(6)
         read(nget,*) w1,w2,w3,w4,w5,w6,w7,fcoefte_tanh(7)
         if(ncoefte_tanh >= 9) read(nget,*) w1,w2,fcoefte_tanh(8)
         if(ncoefte_tanh == 9) read(nget,*) w1,w2,fcoefte_tanh(9)
      elseif(isdatfmtnew == 2)
         read(nget,*) fcoefte_tanh(1),w1,w2
         read(nget,*) fcoefte_tanh(2),w1,w2
         read(nget,*) fcoefte_tanh(3),w1
         read(nget,*) fcoefte_tanh(4),w1
         read(nget,*) fcoefte_tanh(5),w1,w2,w3,w4,w5,w6,w7
         read(nget,*) fcoefte_tanh(6),w1,w2,w3,w4,w5,w6,w7
         read(nget,*) fcoefte_tanh(7),w1,w2,w3,w4,w5,w6,w7
         if(ncoefte_tanh >= 8) read(nget,*) fcoefte_tanh(8),w1,w2
         if(ncoefte_tanh == 9) read(nget,*) fcoefte_tanh(9),w1,w2
      endif

      close (nget)

      return
      end
c-----------------------------------------------------------------------
c******************************************************************************
c  Subroutine to read ti_spline data and evaluate for input locations
c******************************************************************************

      subroutine readti_dat(fname)
       
      implicit none
      character*(*) fname
      Use(Dim)
      Use(Fitdata)

c     local variables --
      integer ios, n, nget, ii
      character*30 str1, str2, ptnam, fit_nam, dumm
      character*7 str3
      character*2 str4
      character*10 yunits
      character*80 line

c     procedures --
      external freeus, gallot, readrt1, xerrab

c----------------------------------------------------------------------c
c     Read spline fit data from file (T. Osborne)
c----------------------------------------------------------------------c

ccc      call freeus(nget)
      data nget /55/
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call xerrab('**** tispline.dat file not found')
      endif

      read (nget,*) str1, str2, ptnam
      read (nget,'(a7,a2,a10)') str1, str2, yunits
      read (nget,*) str1, str2, fit_nam

      read(nget,*) str1, str2, numt_bs
      read(nget,*) dumm

      call gchange("Fitdata",0)         # allocate fit_t(numc)
      do ii = 1, numt_bs
        read(nget,*) fit_t_bs(ii)
      enddo
      read(nget,*) str1, str2, numc_bs
      read(nget,*) dumm

      call gchange("Fitdata",0)         # allocate fit_coef(numc)
      do ii = 1, numc_bs
        read(nget,*) fcoef_bs(ii)
      enddo
      read(nget,*) str1, str2, numk_bs
 
      close (nget)

      return
      end
c----------------------------------------------------------------------c

c*********************************************************
c Subroutine to call tanh fit for ne,te & spline for ti
c*********************************************************

      subroutine fit_neteti
       
      implicit none
      Use(Dim)              # ny
      Use(Xpoint_indices)   # iysptrx
      Use(Comgeo)           # psinormc
      Use(Fitdata)          # psishift,dumfit,nefit,tefit,tifit

c...  Local variables
      integer ir, inbv, iflag1
      real wrk1(3*(numk_bs+1))
      real psi_u(0:ny+1)

c...  Procedures
      real B1VAhL
      external B1VAhL

c...  Initialize inbv on entry for B1VAhL, then dont change
      inbv = 1

c...  Shift psi used by psishift to allow control of separatrix
      do ir = 0, ny+1
         psi_u(ir) = psishift + psinormc(ir)
      enddo

c...  Use tanh function to evaluate ne and Te
      call tanh_multi(ncoefne_tanh, fcoefne_tanh, ny+2, psi_u,
     &                fit_paramne_tanh, dumfit) 
c...  Copy & scale nefit by 1e20 to get correct units (1/m**3)
      do ir = 0, ny+1
        nefit(ir,ifitset) = 1e20*dumfit(ir)
      enddo
      call tanh_multi(ncoefte_tanh, fcoefte_tanh, ny+2, psi_u,
     &                fit_paramte_tanh, dumfit)
      do ir = 0, ny+1
        tefit(ir,ifitset) = dumfit(ir)
      enddo

c...  Use B_spline to evaluate Ti
      do ir = 0, iysptrx
         tifit(ir,ifitset) = B1VAhL(psi_u(ir),0,fit_t_bs,numc_bs,
     &                       numk_bs+1,fcoef_bs,inbv,wrk1,iflag1)
      enddo 

      return
      end   
c ***  End of subroutine fit_neteti ------------------------------c
c-----------------------------------------------------------------------c

c******************************************************************************
c  Subroutine to read ne, te, ti data and evaluate for input locations
c******************************************************************************

      subroutine read_exp_fit(fname)
       
      implicit none
      character*(*) fname
      Use(Dim)
      Use(Fitdata)     #epsi_fit,eprofile_fit,isprofvspsi,yyc_fit
     
c     local variables --
      integer ios, n, nget, ii
      character*30 prof_name

c     procedures --
      external gallot, xerrab
c----------------------------------------------------------------------c
c     Read expt profiles from file
c----------------------------------------------------------------------c

      data nget /55/
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call xerrab('**** D3D_fit file not found')
      endif

      read (nget,*) prof_name
      read (nget,*) num_elem

      call gchange("Fitdata",0)      # allocate epsi_fit,eprofile_fit
      do ii = 1, num_elem
        if (isprofvspsi == 1) then
          read(nget,*) epsi_fit(ii), eprofile_fit(ii)
        else
          read(nget,*) yyc_fit(ii), eprofile_fit(ii)
        endif
      enddo
     
      close (nget)

      return
      end
c-----------------------------------------------------------------------



