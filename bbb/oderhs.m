c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"

c-----------------------------------------------------------------------
      subroutine fd2tra (nx, ny, flox, floy, difx, dify, phi,
     .                   trax, tray, pos, meth)

*//documentation//
*
*
*  1. purpose 
*
*     FD2TRA computes the two-dimensional field of flow of some
*     quantity that is transported by convection and conduction.
*
*
*  2. specification
*
*     subroutine fd2tra (nx, ny, flox, floy, difx, dify, phi,
*    .                   trax, tray, pos, meth)  
*
*     integer nx, ny, pos, meth
*     (0:*) 'real' flox, floy, difx, dify, phi, trax, tray
*
*
*  3. description
*
*     This routine is part of the COCONUT package, ref. /1/.
*     This is a modification of the original B2 routine by Bas Braams
*
*
*  4. references
*
*
*  5. arguments
*
*
*  6. error indicators
*
*     If an error in the inputs is detected the routine will abort
*     through a call to subroutine xerrab.
*
*
*  7. auxiliary routines
*
*     none
*
*
*=======================================================================
*//declarations//

      implicit none

*  -- input arguments --
      integer nx, ny, pos, meth
      real flox(0:nx+1,0:ny+1), floy(0:nx+1,0:ny+1)
      real difx(0:nx+1,0:ny+1), dify(0:nx+1,0:ny+1)
      real phi(0:nx+1,0:ny+1)

*  -- common blocks --
      Use(Selec)   # i1,i4,i5,i8,j1,j4,j5,j8,ixp1
      Use(Noggeo)  # fxm,fx0,fxp,fxmy,fxpy
      Use(Share)   # isnonog,cutlo
      Use(PandfTiming)

      real tick,tock
      external tick tock


*  -- output arguments --
      real trax(0:nx+1,0:ny+1), tray(0:nx+1,0:ny+1)

*  -- local scalars --
      integer ix, iy, ix1, ix2, posx, posy, methx, methy
      real f, p1, p2, tpv, py0, py1

*  -- procedures --
      real upwind
      upwind(f, p1, p2) = max(f, 0.0e0) * p1 + min(f, 0.0e0) * p2

c..   note: dim(a,b) = max((a-b),0)
*=======================================================================
*//computation//
*=======================================================================

      if (TimingPandfOn.gt.0) Timefd2tra=tick()
*  ---------------------------------------------------------------------
*  -- auxiliaries --
*  ---------------------------------------------------------------------
      posx   = mod (pos, 10)
      posy   = pos / 10
      methx  = mod (meth, 10)
      methy  = meth / 10

*  ---------------------------------------------------------------------
*  compute the parallel transport.
*  ---------------------------------------------------------------------

*  -- simulate a CASE statement --
      goto (1, 10, 20, 30, 40, 50, 60, 70), abs(methx) + 1
*  -- if code drops through this goto, improper value of methx
        call xerrab('*** methx has improper value in fd2tra ***')

*  --------------------------------------------------------------------
*  -- /meth/ = 0 --
*  Assume flo = 0 and use central differencing.
    1 continue
      do 3 iy = j4, j8
         do 2 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            trax(ix2,iy) = -difx(ix2,iy)*
     .                          (phi(ix1,iy)-phi(ix,iy))
    2    continue
    3 continue
      goto 100

*  ---------------------------------------------------------------------
*  -- /meth/ = 1 --
*  Assume dif = 0 and use upwind differencing.
   10 continue
      do 14 iy = j4, j8
         do 12 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            trax(ix2,iy) = upwind(flox(ix2,iy),
     .                                phi(ix,iy),phi(ix1,iy))
   12    continue
   14 continue
      goto 100

*  ---------------------------------------------------------------------
*  -- /meth/ = 2 --
*  Regular central differencing.
   20 continue
      do 23 iy = j4, j8
         do 22 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            trax(ix2,iy) = flox(ix2,iy) *
     .                        (phi(ix1,iy)+phi(ix,iy))/2. -
     .                       difx(ix2,iy)*(phi(ix1,iy)-phi(ix,iy))
   22    continue
   23 continue
      goto 100

*  ------------------------------------------------------------------
*  -- /meth/ = 3 --
*  Regular upwind differencing.
   30 continue
      do 33 iy = j4, j8
         do 32 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            trax(ix2,iy) = upwind(flox(ix2,iy),
     .                                phi(ix,iy), phi(ix1,iy)) -
     .                       difx(ix2,iy)*(phi(ix1,iy)-phi(ix,iy))
   32    continue
   33 continue
      goto 100

*  ------------------------------------------------------------------
*  -- /meth/ = 4 --
*  Simple hybrid scheme
   40 continue
      do 43 iy = j4, j8
         do 42 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            tpv = dim(difx(ix2,iy), abs(flox(ix2,iy))/2.)
            trax(ix2,iy) = upwind(flox(ix2,iy),
     .                                phi(ix,iy), phi(ix1,iy)) -
     .                         tpv * (phi(ix1,iy)-phi(ix,iy))
   42    continue
   43 continue
      goto 100

*  ------------------------------------------------------------------
*  -- /meth/ = 5 --
*  Fifth power scheme.
   50 continue
      do 53 iy = j4, j8
         do 52 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            tpv = difx(ix2,iy) * (1 - abs(flox(ix2,iy))/
     .        max(10.*difx(ix2,iy),abs(flox(ix2,iy)),cutlo))**5
            trax(ix2,iy) = upwind(flox(ix2,iy),
     .                                phi(ix,iy), phi(ix1,iy)) -
     .                         tpv * (phi(ix1,iy)-phi(ix,iy))
   52    continue
   53 continue
      goto 100

*  ------------------------------------------------------------------
*  -- /meth/ = 6 --
*  Regular upwind differencing. Can be used for methg=66 in nonorthogonal diff.
   60 continue
      do 63 iy = j4, j8
         do 62 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            trax(ix2,iy) = upwind(flox(ix2,iy),
     .                                phi(ix,iy), phi(ix1,iy)) -
     .                       difx(ix2,iy)*(phi(ix1,iy)-phi(ix,iy))
   62    continue
   63 continue
      goto 100

*  ------------------------------------------------------------------
*  -- /meth/ = 7 --
*  Regular upwind differencing. Can be used for methg=77 in nonorthogonal diff.
   70 continue
      do 73 iy = j4, j8
         do 72 ix = i1, i5
            ix1 = ixp1(ix,iy)
            ix2= ix*(1 - posx) + ix1*posx
            trax(ix2,iy) = upwind(flox(ix2,iy),
     .                                phi(ix,iy), phi(ix1,iy)) -
     .                       difx(ix2,iy)*(phi(ix1,iy)-phi(ix,iy))
   72    continue
   73 continue
      goto 100

*=======================================================================
  100 continue
      if (isnonog .eq. 1) goto 200
*=======================================================================

*  ------------------------------------------------------------------
*  Compute the radial transport for orthogonal grid (isnonog=0)
*  ------------------------------------------------------------------

*  -- simulate a CASE statement --
      goto (101, 110, 120, 130, 140, 150), abs(methy) + 1
*  -- if code drops through this goto, improper value of methy
        call xerrab('*** methy has improper value in fd2tra ***')

*  ------------------------------------------------------------------
*  -- /meth/ = 0 --
*  Assume flo = 0 and use central differencing.
  101 continue
      do 103 iy = j1, j5-posy
         do 102 ix = i4, i8
            tray(ix,iy+posy) = -dify(ix,iy+posy)*
     .                          (phi(ix,iy+1)-phi(ix,iy))
  102    continue
  103 continue
      return

*  ---------------------------------------------------------------------
*  -- /meth/ = 1 --
*  Assume dif = 0 and use upwind differencing.
  110 continue
      do 114 iy = j1, j5-posy
         do 112 ix = i4, i8
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy),
     .                                phi(ix,iy),phi(ix,iy+1))
  112    continue
  114 continue
      return

*  ---------------------------------------------------------------------
*  -- /meth/ = 2 --
*  Regular central differencing.
  120 continue
      do 123 iy = j1, j5-posy
         do 122 ix = i4, i8
            tray(ix,iy+posy) = floy(ix,iy+posy) *
     .                        (phi(ix,iy+1)+phi(ix,iy))/2. -
     .                       dify(ix,iy+posy)*(phi(ix,iy+1)-phi(ix,iy))
  122    continue
  123 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 3 --
*  Regular upwind differencing.
  130 continue
      do 133 iy = j1, j5-posy
         do 132 ix = i4, i8
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy),
     .                                phi(ix,iy), phi(ix,iy+1)) -
     .                       dify(ix,iy+posy)*(phi(ix,iy+1)-phi(ix,iy))
  132    continue
  133 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 4 --
*  Simple hybrid scheme.
  140 continue
      do 143 iy = j1, j5-posy
         do 142 ix = i4, i8
            tpv = dim(dify(ix,iy+posy), abs(floy(ix,iy+posy))/2.)
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy),
     .                                phi(ix,iy), phi(ix,iy+1)) -
     .                         tpv * (phi(ix,iy+1)-phi(ix,iy))
  142    continue
  143 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 5 --
*  Fifth power scheme.
  150 continue
      do 153 iy = j1, j5-posy
         do 152 ix = i4, i8
            tpv = dify(ix,iy+posy) * (1 - abs(floy(ix,iy+posy))/
     .        max(10.*dify(ix,iy+posy),abs(floy(ix,iy+posy)),cutlo))**5
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy),
     .                                phi(ix,iy), phi(ix,iy+1)) -
     .                         tpv * (phi(ix,iy+1)-phi(ix,iy))
  152    continue
  153 continue
      return

*====================================================================
  200 continue
*====================================================================
*  ------------------------------------------------------------------
*  Compute the radial transport for nonorthogonal grid; isnonog=1
*  ------------------------------------------------------------------

*  -- simulate a CASE statement --
      goto (201, 210, 220, 230, 240, 250, 260, 270, 280), abs(methy) + 1
*  -- if code drops through this goto, improper value of methy
        call xerrab('** methy(isnonog=1) has improper value in fd2tra **')

*  ------------------------------------------------------------------
*  -- /meth/ = 0 --
*  Assume flo = 0 and use central differencing.
  201 continue
      do 203 iy = j1, j5-posy
         do 202 ix = i4, i8
            py0 = fxm (ix,iy,0)*phi(ixm1(ix,iy)  ,iy  ) + 
     .            fx0 (ix,iy,0)*phi(ix           ,iy  ) +
     .            fxp (ix,iy,0)*phi(ixp1(ix,iy)  ,iy  ) +
     .            fxmy(ix,iy,0)*phi(ixm1(ix,iy+1),iy+1) +
     .            fxpy(ix,iy,0)*phi(ixp1(ix,iy+1),iy+1)
            py1 = fxm (ix,iy,1)*phi(ixm1(ix,iy+1),iy+1) + 
     .            fx0 (ix,iy,1)*phi(ix           ,iy+1) +
     .            fxp (ix,iy,1)*phi(ixp1(ix,iy+1),iy+1) +
     .            fxmy(ix,iy,1)*phi(ixm1(ix,iy)  ,iy  ) +
     .            fxpy(ix,iy,1)*phi(ixp1(ix,iy)  ,iy  )
            tray(ix,iy+posy) = -dify(ix,iy+posy) * (py1-py0)
  202    continue
  203 continue
      return

*  ---------------------------------------------------------------------
*  -- /meth/ = 1 --
*  Assume dif = 0 and use upwind differencing.
  210 continue
      do 214 iy = j1, j5-posy
         do 212 ix = i4, i8
            py0 = fxm (ix,iy,0)*phi(ixm1(ix,iy)  ,iy  ) + 
     .            fx0 (ix,iy,0)*phi(ix           ,iy  ) +
     .            fxp (ix,iy,0)*phi(ixp1(ix,iy)  ,iy  ) +
     .            fxmy(ix,iy,0)*phi(ixm1(ix,iy+1),iy+1) +
     .            fxpy(ix,iy,0)*phi(ixp1(ix,iy+1),iy+1)
            py1 = fxm (ix,iy,1)*phi(ixm1(ix,iy+1),iy+1) + 
     .            fx0 (ix,iy,1)*phi(ix           ,iy+1) +
     .            fxp (ix,iy,1)*phi(ixp1(ix,iy+1),iy+1) +
     .            fxmy(ix,iy,1)*phi(ixm1(ix,iy)  ,iy  ) +
     .            fxpy(ix,iy,1)*phi(ixp1(ix,iy)  ,iy  )
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy), py0, py1)
  212    continue
  214 continue
      return

*  ---------------------------------------------------------------------
*  -- /meth/ = 2 --
*  Regular central differencing. 
  220 continue
      do 223 iy = j1, j5-posy
         do 222 ix = i4, i8
            py0 = fxm (ix,iy,0)*phi(ixm1(ix,iy)  ,iy  ) + 
     .            fx0 (ix,iy,0)*phi(ix           ,iy  ) +
     .            fxp (ix,iy,0)*phi(ixp1(ix,iy)  ,iy  ) +
     .            fxmy(ix,iy,0)*phi(ixm1(ix,iy+1),iy+1) +
     .            fxpy(ix,iy,0)*phi(ixp1(ix,iy+1),iy+1)
            py1 = fxm (ix,iy,1)*phi(ixm1(ix,iy+1),iy+1) + 
     .            fx0 (ix,iy,1)*phi(ix           ,iy+1) +
     .            fxp (ix,iy,1)*phi(ixp1(ix,iy+1),iy+1) +
     .            fxmy(ix,iy,1)*phi(ixm1(ix,iy)  ,iy  ) +
     .            fxpy(ix,iy,1)*phi(ixp1(ix,iy)  ,iy  )
            tray(ix,iy+posy) = floy(ix,iy+posy) * (py1+py0)/2. -
     .                              dify(ix,iy+posy)*(py1-py0)
  222    continue
  223 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 3 --
*  Regular upwind differencing.
  230 continue
      do 233 iy = j1, j5-posy
         do 232 ix = i4, i8
            py0 = fxm (ix,iy,0)*phi(ixm1(ix,iy)  ,iy  ) + 
     .            fx0 (ix,iy,0)*phi(ix           ,iy  ) +
     .            fxp (ix,iy,0)*phi(ixp1(ix,iy)  ,iy  ) +
     .            fxmy(ix,iy,0)*phi(ixm1(ix,iy+1),iy+1) +
     .            fxpy(ix,iy,0)*phi(ixp1(ix,iy+1),iy+1)
            py1 = fxm (ix,iy,1)*phi(ixm1(ix,iy+1),iy+1) + 
     .            fx0 (ix,iy,1)*phi(ix           ,iy+1) +
     .            fxp (ix,iy,1)*phi(ixp1(ix,iy+1),iy+1) +
     .            fxmy(ix,iy,1)*phi(ixm1(ix,iy)  ,iy  ) +
     .            fxpy(ix,iy,1)*phi(ixp1(ix,iy)  ,iy  )
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy), py0, py1) -
     .                                 dify(ix,iy+posy)*(py1-py0)
  232    continue
  233 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 4 --
*  Simple hybrid scheme.
  240 continue
      do 243 iy = j1, j5-posy
         do 242 ix = i4, i8
            py0 = fxm (ix,iy,0)*phi(ixm1(ix,iy)  ,iy  ) + 
     .            fx0 (ix,iy,0)*phi(ix           ,iy  ) +
     .            fxp (ix,iy,0)*phi(ixp1(ix,iy)  ,iy  ) +
     .            fxmy(ix,iy,0)*phi(ixm1(ix,iy+1),iy+1) +
     .            fxpy(ix,iy,0)*phi(ixp1(ix,iy+1),iy+1)
            py1 = fxm (ix,iy,1)*phi(ixm1(ix,iy+1),iy+1) + 
     .            fx0 (ix,iy,1)*phi(ix           ,iy+1) +
     .            fxp (ix,iy,1)*phi(ixp1(ix,iy+1),iy+1) +
     .            fxmy(ix,iy,1)*phi(ixm1(ix,iy)  ,iy  ) +
     .            fxpy(ix,iy,1)*phi(ixp1(ix,iy)  ,iy  )
            tpv = dim(dify(ix,iy+posy), abs(floy(ix,iy+posy))/2.)
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy), py0, py1) -
     .                                              tpv * (py1-py0)
  242    continue
  243 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 5 --
*  Fifth power scheme.
  250 continue
      do 253 iy = j1, j5-posy
         do 252 ix = i4, i8
            py0 = fxm (ix,iy,0)*phi(ixm1(ix,iy)  ,iy  ) + 
     .            fx0 (ix,iy,0)*phi(ix           ,iy  ) +
     .            fxp (ix,iy,0)*phi(ixp1(ix,iy)  ,iy  ) +
     .            fxmy(ix,iy,0)*phi(ixm1(ix,iy+1),iy+1) +
     .            fxpy(ix,iy,0)*phi(ixp1(ix,iy+1),iy+1)
            py1 = fxm (ix,iy,1)*phi(ixm1(ix,iy+1),iy+1) + 
     .            fx0 (ix,iy,1)*phi(ix           ,iy+1) +
     .            fxp (ix,iy,1)*phi(ixp1(ix,iy+1),iy+1) +
     .            fxmy(ix,iy,1)*phi(ixm1(ix,iy)  ,iy  ) +
     .            fxpy(ix,iy,1)*phi(ixp1(ix,iy)  ,iy  )
            tpv = dify(ix,iy+posy) * (1 - abs(floy(ix,iy+posy))/
     .        max(10.*dify(ix,iy+posy),abs(floy(ix,iy+posy)),cutlo))**5
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy), py0, py1) -
     .                                              tpv * (py1-py0)
  252    continue
  253 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 6 --
*  Regular upwind differencing with log(phi) interpolation
  260 continue
      do 263 iy = j1, j5-posy
         do 262 ix = i4, i8
            py0 = exp( fxm (ix,iy,0)*log(phi(ixm1(ix,iy)  ,iy  )) + 
     .                 fx0 (ix,iy,0)*log(phi(ix           ,iy  )) +
     .                 fxp (ix,iy,0)*log(phi(ixp1(ix,iy)  ,iy  )) +
     .                 fxmy(ix,iy,0)*log(phi(ixm1(ix,iy+1),iy+1)) +
     .                 fxpy(ix,iy,0)*log(phi(ixp1(ix,iy+1),iy+1)) ) 
            py1 = exp( fxm (ix,iy,1)*log(phi(ixm1(ix,iy+1),iy+1)) + 
     .                 fx0 (ix,iy,1)*log(phi(ix           ,iy+1)) +
     .                 fxp (ix,iy,1)*log(phi(ixp1(ix,iy+1),iy+1)) +
     .                 fxmy(ix,iy,1)*log(phi(ixm1(ix,iy)  ,iy  )) +
     .                 fxpy(ix,iy,1)*log(phi(ixp1(ix,iy)  ,iy  )) ) 
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy), py0, py1) -
     .                                 dify(ix,iy+posy)*(py1-py0)
  262    continue
  263 continue
      return

*  ------------------------------------------------------------------
*  -- /meth/ = 7 --
*  Regular upwind differencing with inverse (1/phi) interpolation
  270 continue
      do 273 iy = j1, j5-posy
         do 272 ix = i4, i8
            py0 = 1/( fxm (ix,iy,0)/phi(ixm1(ix,iy)  ,iy  ) + 
     .                fx0 (ix,iy,0)/phi(ix           ,iy  ) +
     .                fxp (ix,iy,0)/phi(ixp1(ix,iy)  ,iy  ) +
     .                fxmy(ix,iy,0)/phi(ixm1(ix,iy+1),iy+1) +
     .                fxpy(ix,iy,0)/phi(ixp1(ix,iy+1),iy+1) )
            py1 = 1/( fxm (ix,iy,1)/phi(ixm1(ix,iy+1),iy+1) + 
     .                fx0 (ix,iy,1)/phi(ix           ,iy+1) +
     .                fxp (ix,iy,1)/phi(ixp1(ix,iy+1),iy+1) +
     .                fxmy(ix,iy,1)/phi(ixm1(ix,iy)  ,iy  ) +
     .                fxpy(ix,iy,1)/phi(ixp1(ix,iy)  ,iy  ) )
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy), py0, py1) -
     .                                 dify(ix,iy+posy)*(py1-py0)
  272    continue
  273 continue

*  ------------------------------------------------------------------
*  -- /meth/ = 8 --
*  Regular upwind differencing, on an approximate staggered stencil
*  (for velocities)
  280 continue
      do 283 iy = j1, j5-posy
         do 282 ix = i4, i8
            ix1=ixp1(ix,iy)
            py0= (
     .            fxmv (ix,iy,0)*phi(ixm1(ix,iy)  ,iy  )+ 
     .            fx0v (ix,iy,0)*phi(ix           ,iy  )+
     .            fxpv (ix,iy,0)*phi(ixp1(ix,iy)  ,iy  )+
     .            fxmyv(ix,iy,0)*phi(ixm1(ix,iy+1),iy+1)+
     .            fxpyv(ix,iy,0)*phi(ixp1(ix,iy+1),iy+1) )
            py1= (
     .            fxmv (ix,iy,1)*phi(ixm1(ix,iy+1),iy+1)+ 
     .            fx0v (ix,iy,1)*phi(ix           ,iy+1)+
     .            fxpv (ix,iy,1)*phi(ixp1(ix,iy+1),iy+1)+
     .            fxmyv(ix,iy,1)*phi(ixm1(ix,iy)  ,iy  )+
     .            fxpyv(ix,iy,1)*phi(ixp1(ix,iy)  ,iy  ) )
            tray(ix,iy+posy) = upwind(floy(ix,iy+posy), py0, py1) -
     .                                 dify(ix,iy+posy)*(py1-py0)
  282    continue
  283 continue
      if (TimingPandfOn.gt.0) TotTimefd2tra=TotTimefd2tra+tock(Timefd2tra)
      return
      end
c**  End of subroutine fd2tra ******************
c-----------------------------------------------------------------------
      subroutine pandf (xc, yc, neq, time, yl, yldot)

c... Calculates matrix A and the right-hand side depending on the values
c... of xc, yc.
c  Definitions for argument list
c
c  Input variables:
c    xc is poloidal index of perturbed variablefor Jacobian calc, 
c       or =-1 for full RHS evaluation
c    yc is radial index for perturbed variable for Jacobian calc, 
c       or =-1 for full RHS evaluation
c    neq is the total number of variables
c    time is the present physical time; useable by VODPK but not NKSOL
c    yl is the vector of unknowns
c  Output variables:
c    yldot is the RHS of ODE solver or RHS=0 for Newton solver (NKSOL)

      implicit none

*  -- input arguments
      integer xc, yc, neq      
      real time, yl(*),yldot(*)

*  -- set local array dimension
      integer nigmx
      parameter (nigmx=100)

*  -- slocal variables
      integer ifld, jfld, zn, k, k1, k2, jx, ixt, ixt1, ixr, ixr1, iixt,
     .        ixt0
      real fxet, fxit, qr, vt0, vt1, vtn, vtn2, pradold, eeliold, 
     .     erlizold, erlrcold, psorrgold(nigmx), psorcxgold(nigmx),
     .     nuizold(nigmx), nucxold(nigmx), nurcold(nigmx), nuixold(nigmx),
     .     psorgold(nigmx), tsfe, tsjf, niavex, niavey, teave, tiave, tgavex,
     .     zeffave, noavex, noavey, tiavey, tgavey, psordisold, 
     .     nucxiold(nigmx), nueliold(nigmx), nuelgold(nigmx), rrfac, visxtmp,
     .     vttn, vttp, neavex, pwrebkgold, pwribkgold, feexflr, feixflr,
     .     naavex,naavey,nuelmolx,nuelmoly,fniycboave, corecells, sycore
      real fqpo, fqpom, friceo, friceom, upeo, upeom, fricio(100), 
     .     friciom(100), upio(100), upiom(100), uupo(100), uupom(100)
      real nevol, ngvol, kionz, krecz, kcxrz, kionm, krecm, kcxrm, nzbg,
     .     niz_floor, hflux, zflux, psorv, kionz0, pscx0, pxri, kcxrzig,
     .     nizm_floor, argx, massfac, ae, geyym, geyy0, geyyp, dgeyy0,
     .     dgeyy1, te_diss, wallfac, z1fac, bpolmin, rt2nus, epstmp, tv2
      real awoll,awll
      integer izch, ihyd, iimp, jg, jz, nsm1, ifld_fcs, ifld_lcs
      real uuv, ne_sgvi, nbarx, argth, fac_rad, ffyi, ffyo
      real grdnv, qflx, qfly, cshx, cshy, qshx, qshy, lxtec, lxtic
      real lmfpn, lmfppar, lmfpperp
      real temp1, temp2, temp3, temp4, cutlo3, lambd_ci, lambd_ce
      real upxavep1,upxave0,upxavem1,upf0,upfm1
      real teev,nexface,loglmcc,nmxface
      logical xccuts, xcturb
      integer iy1, ixmp2, iyp1, iyp2, iym1, ixs, ixf, iys, iyf,
     .        methnx, methny, iy2, i2pwr, i5pwr, j2pwr, j5pwr,
     .        iysepu, ixpt1u, ixpt2u
      integer iy0, jy, jylo, jyhi, iy_min, iy_max
      real diffustotal
      real difnimix, kyemix, kyimix, v2dia, bscalf, bpfac
      real rdum,rdumaray(1),afqp,ltmax,lmfpe,lmfpi,flxlimf
      real rdumx,rdumy,dr1,dr2,qion
      real b_ctr,dbds_m,dbds_p,eta_h0,eta_hm,eta_hp,drag_1,drag_2,
     .     drag_3,mf_path,nu_ii,frac_col,fniy_recy
      real thetacc,dupdx,dupdy
      real dndym1,dndy0,dndyp1,d2ndy20,d2ndy2p1,d3ndy3
      real dtdym1,dtdy0,dtdyp1,d2tdy20,d2tdy2p1,d3tdy3,nhi_nha
      integer idum, idumaray(1)
      real(Size4) sec4, gettime, tsimpfe, tsimp, tsnpg, ueb
      integer impflag
      # former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real tv,t0,t1,t2,a,t1old,t1new,t2old,t2new
cnxg      data igs/1/

      Use(Dim)      # nx,ny,nhsp,nusp,nzspt,nzsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
                               # iysptrx
      Use(Math_problem_size)   # neqmx,numvar
      Use(Timing)   # istimingon,ttotfe,ttotjf,ttimpfe,ttimpjf,ttnpg
      Use(Share)    # nxpt,nxomit,geometry,nxc,cutlo,islimon,ix_lim,iy_lims
                    # istabon,isudsym
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methe,methu,methn,methi,methg,concap,lnlam,
                    # convis,icnuiz,icnucx,cnuiz,cnucx,isrecmon,
                    # ngbackg,ingb,eion,ediss,afix,coef,ce,ci,
                    # dp1,qfl,csh,qsh,mfl,msh,cs,fxe,ctaue,fxi,ctaui,
                    # zcoef,coef1,nurlxu,isphiofft,isintlog,isgxvon
                    # isgpye,frnnuiz,nzbackg,inzb,nlimix,nlimiy,
                    # isofric,isteaven
                    # isnionxy,isuponxy,isteonxy,istionxy,isngonxy,isphionxy,
                    # isupon,isteon,istion,isphion
      Use(Aux)      # ixmp
      Use(Coefeq)
      Use(Bcond)    # albedoo,albedoi,isfixlb,isfixrb
                    # xcnearlb,xcnearrb,openbox
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Fixsrc)
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,xlinc,xrinc,ixs1,
                    # j1,j1p,j2,j2p,j3,j4,j5,j6,j5p,j6p,j7,j8,j5m
      Use(Comgeo)   # isxptx,isxpty
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx,fxm,fx0,fxp,fxmy,fxpy
      Use(Compla)   # znucl,rtaux,rtauy,rtau,rt_scal
      Use(Comflo)
      Use(Conduc)   # lmfplim
      Use(Rhsides)
      Use(Save_terms)   # psorold,psorxrold
      Use(Indexes)  
      Use(Ynorm)    # isflxvar,nnorm,ennorm,fnorm,n0,n0g
      Use(Poten)    # bcee,bcei,cthe,cthi
      Use(Comtra)   # parvis,travis,difni,difnit,difpr,difni2,
                    # difpr2,vcony,flalfe,flalfi,flgam,kxe,kxecore,
                    # kye,kyet,kxi,kxicore,kxn,kyi,kyit,kyn,feqp,
                    # flalfgx,flalfgy,flalfv,flalfvgx,flalfvgy,
                    # flalftgx,flalftgy,alfeqp,facbni,facbni2,facbee,
                    # facbei,isbohmcalc,diffusivity,diffusivwrk,
                    # diffusivloc,isdifxg_aug,isdifyg_aug,sigvi_floor,
                    # facbup
      Use(Bfield)   # btot,rbfbt,rbfbt2
      Use(Locflux)
      Use(Wkspace)
      Use(Gradients)
      Use(Imprad)   # isimpon,nzloc,impradloc,prad,pradz,na,ntau,nratio,
                    # afrac,atau,ismctab,rcxighg,pwrze

      Use(Volsrc)   # pwrsore,pwrsori,volpsor
      Use(Model_choice)   # iondenseqn
      Use(Time_dep_nwt)   # ylodt
      Use(Cfric)          # frice,frici
      Use(Turbulence)     # isturbnloc,isturbcons,diffuslimit,diffuswgts
      Use(Turbulence_diagnostics)   # chinorml,chinormh
      Use(MCN_dim)
      Use(MCN_sources)	  # uesor_ni,uesor_up,uesor_ti,uesor_te
      Use(Ext_neutrals)          # isextneuton, extneutopt
      Use(PNC_params)            # dtneut, dtold
      Use(PNC_data)              # ni_pnc, etc.
      Use(Reduced_ion_interface) # misotope,natomic
      Use(Indices_domain_dcl)    # ixmxbcl
      Use(Indices_domain_dcg)    # ndomain
      Use(Npes_mpi)              # mype
      Use(RZ_grid_info)  		 # bpol
      Use(Interp)				 # ngs, tgs 
      Use(ParallelEval)          # ParallelJac,ParallelPandf1
      Use(PandfTiming)

*  -- procedures for atomic and molecular rates --
      integer zmax,znuc
      real dene,denz(0:1),radz(0:1)
      real rsa, rra, rqa, rcx, emissbs, erl1, erl2, radneq, radimpmc
      real radmc, svdiss, vyiy0, vyiym1, v2ix0, v2ixm1
      real tgupyface, ngupyface, n1upyface, ng2upyface
      external rsa, rra, rqa, rcx, emissbs, erl1, erl2, radneq, radimpmc
      external radmc, svdiss
      real tick,tock
      external tick,tock
	  
ccc      save

*  -- procedures --
      real ave, etaper, interp_yf, interp_xf
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)
      etaper(ix,iy) = 3.234e-9*loglambda(ix,iy)/(max(te(ix,iy),temin*ev)
     .                                          /(1000.*ev))**(1.5)
      interp_yf(ix,iy,t0,t1) = (t0*gy(ix,iy) + t1*gy(ix,iy+1)) /
     .                                       (gy(ix,iy)+gy(ix,iy+1))
      interp_xf(ix,iy,t0,t1) =(t0*gx(ix,iy) + t1*gx(ixp1(ix,iy),iy)) /
     .                                (gx(ix,iy)+gx(ixp1(ix,iy),iy))
      cutlo3 = cutlo**0.3

c  Check array sizes
      if (ngsp > nigmx .or. nisp > nigmx) then
         call xerrab("***PANDF in oderhs.m: increase nigmx, recompile")
      endif
c... Timing of pandf components (added by J. Guterl)
        if (TimingPandf.gt.0
     . .and. yl(neq+1) .lt. 0 .and. ParallelPandf1.eq.0) then
        TimingPandfOn=1
        else
        TimingPandfOn=0
        endif
        if (TimingPandfOn.gt.0) TimePandf=tick()
c... Roadblockers for  call to pandf through openmp structures (added by J.Guterl)
      if ((isimpon.gt.0 .and. isimpon.ne.6) .and. (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)) then
      call xerrab('Only isimpon=0 or 6 is validated with openmp.
     .Contact the UEDGE team to use other  options with openmp.')
      endif

      if ((ismcnon.gt.0) .and. (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)) then
      call xerrab('Only ismcnon=0 is validated with openmp.
     .Contact the UEDGE team to use other options with openmp.')
      endif


************************************************************************
*  -- initialization --
************************************************************************
c     Set ix index for outer midplane turbulence
      if (isudsym==1) then
         ixmp2 = nxc + 2
      elseif (geometry=='dnull' .or. geometry(1:9)=="snowflake" .or.
     .        geometry=="dnXtarget" .or. geometry=="isoleg") then
         ixmp2 = ixmdp(2) + 1
      else
         ixmp2 = ixpt1(1) + nxomit + 3*(ixpt2(1)-ixpt1(1))/4
      endif

c     Set switches for neutrals-related source terms in plasma equations
c     (MER 1996/10/28)
c     (IJ  2015/04/06) add ismcnon>=3 for external call to run_neutrals 
      if (ismcnon .eq. 1) then        # use MC sources only:
         if (svrpkg .eq. "cvode") then
           call xerrab('*** ismcnon=1 not allowed for cvode ***')
         endif
         cfneut=0.
         if (isupgon(1) .eq. 1) then
            cfvgpx(iigsp)=0.
            cfvgpy(iigsp)=0.
            cfvcsx(iigsp)=0.
            cfvcsy(iigsp)=0.
         endif
         cmneut=1.
      else if (ismcnon .eq. 2) then    # switch between two models:
         if (yl(neq+1) .gt. 0) then   # use fluid model for Jacobian
            cfneut=1.
            if (isupgon(1) .eq. 1) then
               cfvgpx(iigsp)=1.
               cfvgpy(iigsp)=1.
               cfvcsx(iigsp)=1.
               cfvcsy(iigsp)=1.
            endif
            cmneut=0.
         elseif (yl(neq+1) .lt. 0) then     # use MC model for evaluation
            cfneut=0.
            if (isupgon(1) .eq. 1) then
               cfvgpx(iigsp)=0.
               cfvgpy(iigsp)=0.
               cfvcsx(iigsp)=0.
               cfvcsy(iigsp)=0.
            endif
            cmneut=1.
         else
            call xerrab('*** PANDF: ismcnon=2 & yl(neq+1)=0 ???')
         endif
      else if (ismcnon .eq. 3) then         # switch between two neutral models internally
         if (yl(neq+1) .gt. 0) then         # use fluid model for preconditioner
            if (extneutmeth .eq. 1) then				#fluid source & implicit MC flux
               cfneut=1.     #turn on  fluid sources
               cfneutdiv=0.  #turn off fluid div fluxes  
               cmneut=0.     #turn off MC sources
               cmneutdiv=1.  #turn on  MC div fluxes    
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif         
            else         								#fluid source & fluid flux
               cfneut=1.     #turn on fluid sources
               cfneutdiv=1.  #turn on fluid div fluxes  
               cmneut=0.     #turn off MC sources
               cmneutdiv=0.  #turn off MC div fluxes    
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif
            endif 
         elseif (yl(neq+1) .lt. 0) then     # use MC model for RHS evaluation
            if (extneutmeth .eq. 1) then				#fluid source & implicit MC flux
               cfneut=1.     #turn on  fluid sources
               cfneutdiv=0.  #turn off fluid div fluxes  
               cmneut=0.     #turn off MC sources
               cmneutdiv=1.  #turn on  MC div fluxes             
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif
            else         								#MC source & fluid flux
               cfneut=0.     #turn off fluid sources
               cfneutdiv=1.  #turn on  fluid div fluxes  
               cmneut=1.     #turn on  MC sources
               cmneutdiv=0.  #turn off MC div fluxes    
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=0.
                  cfvgpy(iigsp)=0.
                  cfvcsx(iigsp)=0.
                  cfvcsy(iigsp)=0.
               endif
            endif 

         else
            call xerrab('*** PANDF: ismcnon=3 & yl(neq+1)=0 ???')
         endif
         if (yl(neq+1) .lt. 0) then       # RHS eval (not precon eval)
			if (isextneuton .ne. 0) then  # implicit use of external neutrals inside exmain
                #Neutral step
                dtold=dtreal
                dtreal=dtneut
                call store_neutrals
                call run_neutrals		  # extneutopt sets choice of model
                call update_neutrals
                dtreal=dtold
            endif
        endif
      else if (ismcnon .eq. 4) then # test a different fluid model in the preconditioner
         if (yl(neq+1) .gt. 0) then   # Precon eval
            parvis=parvis*pnc_cfparvis
            travis=travis*pnc_cftravis
            do ifld=1,nisp
              ni(:,:,ifld)=ni(:,:,ifld)*pnc_cfni(ifld)
              up(:,:,ifld)=up(:,:,ifld)*pnc_cfup(ifld)
            enddo
c            write(*,*) 'ismcnon=4'
c            write(*,*) parvis
         endif
      end if #ismcnon

************************************************************************
*   This section is to use in the calculation of the jacobian locally.
************************************************************************

c ... Get initial value of system cpu timer.
      if(xc .lt. 0) then
         tsfe = gettime(sec4)
      else
         tsjf = gettime(sec4)
      endif

      if ( (xc .lt. 0) .or. 
     .     ((0<=yc).and.(yc-yinc<=0)).and.isjaccorall==1 ) then  
                                              # use full ix range near yc=0
                                              # with integrated core flux BC
         i1 = 0  # 1-ixmnbcl
         i2 = 1
         i2p = 1
         i3 = 0  # 1-ixmnbcl
         i4 = 0  # 1-ixmnbcl
         i5 = nx
         i5m = nx-1
         i6 = nx+1  # nx+ixmxbcl
         i7 = nx+1  # nx+ixmxbcl
         i8 = nx+1  # nx+ixmxbcl
      else
         i1 = max(0, xc-xlinc-1)
         i2 = max(1, xc-xlinc)
         i2p = max(1, xc-xrinc-1)
         i3 = xc-xlinc     # not used as loop indice (can be < 0)
         i4 = max(0, xc-xlinc)
         i5 = min(nx, xc+xrinc)
         i5m = min(nx-1, xc+xrinc)
         i6 = min(nx+1, xc+xrinc+1)
         i7 = xc+xrinc     # not used as loop indice (can be > nx)
         i8 = min(nx+1, xc+xrinc)
      endif
      if (yc .lt. 0) then
         j1 = 0 
         j1p = 0
         j2 = 1
         j2p = 1
         j3 = 0 
         j4 = 0 
         j5 = ny
         j5m = ny-1
         j6 = ny+1 
         j5p = ny
         j6p = ny+1
         j7 = ny+1 
         j8 = ny+1 
      else
         j1 = max(0, yc-yinc-1)
         j2 = max(1, yc-yinc)
         j1p = max(0, yc-yinc-2)
         j2p = max(1, yc-yinc-1)
         j3 = yc-yinc
         j4 = max(0, yc-yinc)
         j5 = min(ny, yc+yinc)
         j5m = min(ny-1, yc+yinc)
         j6 = min(ny+1, yc+yinc)
         j5p = min(ny, yc+yinc+1)
         j6p = min(ny+1, yc+yinc+1)
c         j6 = min(ny+1, yc+yinc+1)
         j7 = yc+yinc
         j8 = min(ny+1, yc+yinc)
      endif

c...  We will expand the range of possible responses when perturbing the
c...  plasma in a cell near one of the cuts.
      xccuts = .false.
      do jx = 1, nxpt
        if ( (xc-xlinc<=ixpt1(jx)+1) .and. (xc+xrinc+1>=ixpt1(jx)) .and.
     .       (yc-yinc<=iysptrx1(jx)) .and. (iysptrx1(jx)>0) ) xccuts=.true.
        if ( (xc-xlinc<=ixpt2(jx)+1) .and. (xc+xrinc+1>=ixpt2(jx)) .and.
     .       (yc-yinc<=iysptrx2(jx)) .and. (iysptrx2(jx)>0) ) xccuts=.true.
      enddo

c...  We must expand the range of ix in the vicinity of cells on
c...  which turbulent transport depends.
      xcturb = .false.
      do jx = 1, nxpt
         xcturb = xcturb .or. (xc.eq.ixlb(jx).and.ixmnbcl==1) .or.
     .                        (xc.eq.(ixrb(jx)+1).and.ixmxbcl==1)
      enddo
      xcturb = xcturb .or. (xc.eq.ixmp2)
      xcturb = xcturb .and. (isturbnloc.eq.1)
c...  NOTE: For a full double-null configuration, if there are 2 separatrices
c...  we use the innermost one (at iysptrx) to define the radial boundary
c...  of the turbulence.
      if (isturbcons .eq. 1) then
         xcturb = xcturb .and. yc .eq. iysptrx+1
      elseif (isturbcons .eq. 2) then
         xcturb = xcturb .and. yc .ge. iysptrx+1-diffuslimit
      else
         xcturb = xcturb .and. yc .ge. iysptrx+1
      endif

      if (xccuts .or. xcturb) then
         i1 = 0  
         i2 = 1
         i3 = 0  
         i4 = 0  
         i5 = nx
         i6 = nx+1 
         i7 = nx+1 
         i8 = nx+1 
      endif

c...  Define range for source terms to minimize calls to adpak-based routines
            ixs = i2
            ixf = i5
            iys = j2
            iyf = j5
            ixs1 = i1
            ixf6 = i6
            iys1 = j1
            iyf6 = j6
c...  Reset ioniz. and rad. indices to a point if this is a Jacobian calc.
         if (xc .ge.0 .and. yc .ge. 0) then
            ixs = xc
            ixf = xc
            iys = yc
            iyf = yc
            ixs1 = xc
            ixf6 = xc
            if (xrinc .ge. 20) then
               ixs1 = 0  
               ixf6 = nx+1 
            endif
            iys1 = yc
            iyf6 = yc
            if (yinc .ge. 20) then
               iys1 = 0 
               iyf6 = ny+1 
            endif
         endif

c...  Set flag that indicates wide open Jac'n "box" for subroutine bouncon.
      if (xc .lt. 0) then
         openbox = .true.
      elseif (xccuts .or. xcturb) then
         openbox = .true.
      elseif ( (0<=yc).and.(yc<=yinc) ) then # for integrated core flux BC
         openbox = .true.
      else
         openbox = .false.
      endif

c...  Set flags that indicate when the Jac'n "box" overlaps left or right
c...  boundary cells of a mesh region.  Used in subroutine bouncon.
      xcnearlb = .false.
      do jx = 1, nxpt
         xcnearlb = xcnearlb .or.
     .       ( (xc-xlinc.le.ixlb(jx)) .and. (xc+xrinc.ge.ixlb(jx)) )
      enddo
      if (xc .lt. 0) xcnearlb = .true.
      xcnearrb = .false.
      do jx = 1, nxpt
         xcnearrb = xcnearrb .or.
     .      ( (xc-xlinc.le.ixrb(jx)+1) .and. (xc+xrinc.ge.ixrb(jx)) )
      enddo
      if (xc .lt. 0) xcnearrb = .true.

************************************************************************
c... First, we convert from the 1-D vector yl to the plasma variables.
************************************************************************
         if (TimingPandfOn.gt.0) TimeConvert0=tick()
         call convsr_vo (xc, yc, yl)  # pre 9/30/97 was one call to convsr
         if (TimingPandfOn.gt.0) TotTimeConvert0=TotTimeConvert0+tock(TimeConvert0)
         if (TimingPandfOn.gt.0) TimeConvert1=tick()
         call convsr_aux (xc, yc)
         if (TimingPandfOn.gt.0) TotTimeConvert1=TotTimeConvert1+tock(TimeConvert1)

c ... Set variable controlling upper limit of species loops that
c     involve ion-density sources, fluxes, and/or velocities.

      nfsp = nisp
      if (isimpon .eq. 3 .or. isimpon .eq. 4) nfsp = nhsp

c ... Calculate the Bohm diffusion rates (units are m**2/s)
      do ifld = 1, nisp
       if (facbni+facbup+facbee+facbei>0 .and. isbohmcalc>0) then
         do iy = j1, j6
            iyp1 = min(ny+1, iy+1)
            do ix = i1, i6
               ix1 = ixp1(ix,iy)
               kybohm(ix,iy) = (te(ix,iy)+te(ix,iyp1)) /
     .                        (16*ev*(btot(ix,iy)+btot(ix,iyp1)))
               kxbohm(ix,iy) = (te(ix,iy)+te(ix1,iy)) /
     .                        (16*ev*(btot(ix,iy)+btot(ix1,iy)))
               dif_use(ix,iy,ifld)  = facbni*kybohm(ix,iy)
               dif2_use(ix,iy,ifld) = facbni2*kxbohm(ix,iy)
               tray_use(ix,iy,ifld)  = facbup*kybohm(ix,iy)
               kye_use(ix,iy)  = facbee*kybohm(ix,iy)
               kyi_use(ix,iy)  = facbei*kybohm(ix,iy)
	       dutm_use(ix,iy,ifld) = kybohm(ix,iy)
            enddo
         enddo
         if (isbohmcalc.eq.2) then  # calc. recip. average with const D
           fcdif = 0.   # used to zero constant diff. if recip. ave used
           do iy = j1, j6
             do ix = i1, i6
               dif_use(ix,iy,ifld)  = 0.5*ave(difni(ifld),  dif_use(ix,iy,ifld))
               dif2_use(ix,iy,ifld) = 0.5*ave(difni2(ifld), dif2_use(ix,iy,ifld))
               tray_use(ix,iy,ifld)  = 0.5*ave(travis(ifld), tray_use(ix,iy,ifld))
               kye_use(ix,iy)  = 0.5*ave(kye, kye_use(ix,iy))
               kyi_use(ix,iy)  = 0.5*ave(kyi, kyi_use(ix,iy))
               dutm_use(ix,iy,ifld)  = 0.5*ave(difutm(1), dutm_use(ix,iy,ifld))
             enddo
           enddo
         endif
       endif

c ... If isbohmcalc=3, then give (B0/B)**inbdif profile to diff
       if (isbohmcalc==3) then  # use inbtdif, inbpdif for btot, bpol scaling
         bpolmin = bpol(ixpt2(1)-ixbpmin,iysptrx,0)
         do iy = j1, j6
            do ix = i1, i6
              ix1 = ixp1(ix,iy)
	      bscalf=((.5*(btot(ixmp,iysptrx)/btot(ix,iy)) +
     .               .5*(btot(ixmp,iysptrx)/btot(ix1,iy)))**inbtdif)
     .       * ((bpol(ixmp,iysptrx,3)+bpol(ixmp,iysptrx,4))/
     .          (bpol(ix,iy,3)+bpol(ix,iy,4)+bpolmin))**inbpdif
	      dif_use(ix,iy,ifld)  = difniv(iy,ifld)*bscalf
	      difp_use(ix,iy,ifld) = difprv(iy,ifld)*bscalf          
              dif2_use(ix,iy,ifld) = difniv2(iy,ifld)*bscalf
              tray_use(ix,iy,ifld)  = travisv(iy,ifld)*bscalf
              kye_use(ix,iy)  = kyev(iy)*bscalf
              kyi_use(ix,iy)  = kyiv(iy)*bscalf
              dutm_use(ix,iy,ifld) = difutmv(iy,ifld)*bscalf
	      vy_use(ix,iy,ifld) = vconyv(iy,ifld)*bscalf
            enddo
         enddo
       endif
      enddo  # loop over species lfld

c ,,, Add diffusion propto betap**iexpbp and (B0/B)**inbdif (as for isbohmcalc=3)
      if (isdifbetap == 1) then
       do ifld = 1, nisp
         if(zi(ifld) > 0.) then
           bpolmin = bpol(ixpt2(1)-ixbpmin,iysptrx,0)
           do iy = j1, j6
             do ix = i1, i6
               ix1 = ixp1(ix,iy)
               betap(ix,iy) = 8.*pi*1.e-7*pr(ix,iy)/bpol(ix,iy,0)**2
               bpfac = betap(ix,iy)**iexpbp
 	       bscalf = ((.5*(btot(ixmp,iysptrx)/btot(ix,iy)) +
     .                   .5*(btot(ixmp,iysptrx)/btot(ix1,iy)))**inbtdif)
     .                  * ((bpol(ixmp,iysptrx,3)+bpol(ixmp,iysptrx,4))/
     .                   (bpol(ix,iy,3)+bpol(ix,iy,4)+bpolmin))**inbpdif
	       dif_use(ix,iy,ifld)  = difniv(iy,ifld)*bscalf +
     .                                                    dfacbp*bpfac
	       difp_use(ix,iy,ifld) = difprv(iy,ifld)*bscalf          
               dif2_use(ix,iy,ifld) = difniv2(iy,ifld)*bscalf + 
     .                                                    dfacbp*bpfac
               tray_use(ix,iy,ifld)  = travisv(iy,ifld)*bscalf + 
     .                                                   trfacbp*bpfac
               trax_use(ix,iy,ifld) = trfacbp*bpfac
               kye_use(ix,iy)  = kyev(iy)*bscalf +  kefacbp*bpfac
               kyi_use(ix,iy)  = kyiv(iy)*bscalf + kifacbp*bpfac
               kxe_use(ix,iy) = kefacbp*bpfac
               kxi_use(ix,iy) = kifacbp*bpfac
               dutm_use(ix,iy,ifld) = difutmv(iy,ifld)*bscalf
	       vy_use(ix,iy,ifld) = vconyv(iy,ifld)*bscalf
             enddo
           enddo
         endif
       enddo
      endif   # test on isdifbetap


************************************************************************
*  Transverse Drifts in y-direction and in 2-direction 
*  (normal to B and y)
************************************************************************
*  ---------------------------------------------------------------------
*  compute drifts
*  ---------------------------------------------------------------------

c ... Compute log_lambda
      do iy = j1, j6
        do ix = i1, i6
          ix1 = ixp1(ix,iy)
          teev = 0.5*(te(ix,iy)+te(ix1,iy))/ev
          nexface = 0.5*(ne(ix,iy)+ne(ix1,iy))
          if (islnlamcon == 1) then
            loglambda(ix,iy) = lnlam  # set to constant
          elseif (teev < 50.) then    # Braginskii formula, teev<50
            loglambda(ix,iy) = 23.4-1.15*log10(1.e-6*nexface)+
     .                              3.45*log10(teev)
          else                        #teev > 50
            loglambda(ix,iy) = 25.3-1.15*log10(1.e-6*nexface)+
     .              2.33167537087122D+00*log10(teev)
          endif
          ctaui(ix,iy,1) = 2.1e13*sqrt(mi(1)/mp)/ loglambda(ix,iy)
          ctaue(ix,iy,1) = 3.5e11/loglambda(ix,iy) #both for zi=1
        enddo
      enddo

c ... Calculate collis. factors eta1 and rtaue for the simple Braginski model
      do iy = j1, j6
        do ix = i1, i6
           eta1(ix,iy) = cfeta1*0.3*nm(ix,iy,1)*ti(ix,iy)*
     .                   (1/(qe*btot(ix,iy))) / omgci_taui
           rtaue(ix,iy) = cfrtaue*(1/(qe*btot(ix,iy))) / omgce_taue
           dclass_i(ix,iy) = cfcl_i*eta1(ix,iy)/(0.3*nm(ix,iy,1))
           dclass_e(ix,iy) = cfcl_e*te(ix,iy)*rtaue(ix,iy)
        enddo
      enddo

*  -- loop over species number --

      do 100 ifld = 1, nfsp
c --- If this is the neutral species (zi(ifld).eq.0)) we dont want velocities
        if(zi(ifld) > 1.e-10) then  # if not, skip to end of 100 loop
         qion = zi(ifld)*qe
         do 18 iy = j1, j5
            iyp1 = min(iy+1,ny+1)
            iym1 = max(iy-1,0)
            do 17 ix = i1, i6
              ix3 = ixm1(ix,iy)
              ix4 = ixm1(ix,iy+1)
              temp1 = 
     .                ( ex(ix ,iy) + ex(ix ,iy+1) +
     .                  ex(ix3,iy) + ex(ix4,iy+1) )
c... sknam: grad P from priv, prev
              temp1 = (-4.0)*(phiv(ix,iy) - phiv(ix3,iy))*gxc(ix,iy)
              temp2 = 4.0*(priv(ix,iy,ifld) - priv(ix3,iy,ifld))*gxc(ix,iy) 
              temp3 = 4.0*(prev(ix,iy) - prev(ix3,iy))*gxc(ix,iy) 
 
c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
              if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 ) then
                temp1 = (-4.0)*(phiv(ix,iy) - phiv(ix3,iy))*gxc(ix,iy)
                temp2 = 4.0*(priv(ix,iy,ifld) - priv(ix3,iy,ifld))*gxc(ix,iy) 
                temp3 = 4.0*(prev(ix,iy) - prev(ix3,iy))*gxc(ix,iy) 
              endif
c...    Calc collisionality factors nu_s/(1 + nu_s) = 1/(1 + lambda_s)
              lambd_ci = 1e16*(ti(ix,iy)/ev)**2/nit(ix,iy)  # approx
              lambd_ce = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)   # approx
              coll_fi(ix,iy) = cfnus_i/(cfnus_i + (lambd_ci/(lconi(ix,iy))))
              coll_fe(ix,iy) = cfnus_e/(cfnus_e + (lambd_ce/(lcone(ix,iy))))
              vyce(ix,iy,ifld) = 0.125 * temp1
     .                          * ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) )
              vycb(ix,iy,ifld) = ( cfcurv*( 0.5*(ti(ix,iy)+ti(ix,iyp1)) +
     .                mi(ifld)*(0.25*(up(ix,iy,ifld)+up(ix,iyp1,ifld)+
     .                         up(ix3,iy,ifld)+up(ix4,iyp1,ifld)))**2 )*
     .                                       curvrby(ix,iy)/qion +
     .                   cfgradb*0.5*( ti(ix,iy)+ti(ix,iyp1) )*
     .                               gradby(ix,iy)/qion )*coll_fi(ix,iy)
              veycb(ix,iy) = ( -cfcurv*0.5*(te(ix,iy)+te(ix,iyp1))*
     .                                               curvrby(ix,iy)/qe -
     .                   cfgradb*0.5*( te(ix,iy)+te(ix,iyp1) )*
     .                                  gradby(ix,iy)/qe )*coll_fe(ix,iy)

              vycp(ix,iy,ifld) = -0.25 * temp2
     .               * (rbfbt2(ix,iy)+rbfbt2(ix,iy+1)) /
     .            (qion*(niy0(ix,iy,ifld)+niy1(ix,iy,ifld)))
              veycp(ix,iy) =  0.25 * temp3
     .               * (rbfbt2(ix,iy)+rbfbt2(ix,iy+1)) /
     .                    (qe*(ney0(ix,iy)+ney1(ix,iy)))
c...   zero the vy-diamagnetic velocity on the y guard-cell faces
              vycp(ix,0,ifld) = 0.
              vycp(ix,ny,ifld) = 0.
              veycp(ix,0) = 0.
              veycp(ix,ny) = 0.

c...  Precompute radial velocities from fixed BOUT turbulence fluxes
              vy_cft(ix,iy,ifld) = 2*fniyos_use(ix,iy,ifld)/
     .                            (niy0(ix,iy,ifld)+niy1(ix,iy,ifld))
              vyte_cft(ix,iy) = 2*feeyosn_use(ix,iy)/
     .                            (tey0(ix,iy)+tey1(ix,iy))
              vyti_cft(ix,iy) = 2*feiyosn_use(ix,iy)/
     .                            (tiy0(ix,iy)+tiy1(ix,iy))
              vyrd(ix,iy,ifld) = - 2. * gpry(ix,iy) /
     .                       (btot(ix,iy)**2/etaper(ix,iy) +
     .                        btot(ix,iy+1)**2/etaper(ix,iy+1) )
              vydd(ix,iy,ifld) = vcony(ifld) + vy_use(ix,iy,ifld) +
     .                                         vy_cft(ix,iy,ifld) -
     .                         (difpr(ifld) + difp_use(ix,iy,ifld)) * 
     .                    ( 2*gpry(ix,iy)/(pr(ix,iy+1) + pr(ix,iy)) -
     .                      3.0*gtey(ix,iy)/(tey1(ix,iy)+tey0(ix,iy)) )
c ...   Note that the density grad. term for vydd added below
           if (cfrtaue.ne.0.) then  #special classical mom. transfer term
              vycr(ix,iy) = -0.5*(rtaue(ix,iy)+rtaue(ix,iyp1)) * (
     .                        (gpiy(ix,iy,1) + gpey(ix,iy))/
     .                        (0.5*(niy1(ix,iy,1)+niy0(ix,iy,1))) -
     .                          1.5*gtey(ix,iy) )
           endif
           if (cfeta1.ne.0. .and. iy.le.ny-1 .and. iy.gt.0) then  
                                             #special classical vis. term
              geyym = 2*gpiy(ix,iym1,1)/(ney1(ix,iym1)+ney0(ix,iym1))-
     .                qe*ey(ix,iym1)
              geyy0 = 2*gpiy(ix,iy,1)/(ney1(ix,iy)+ney0(ix,iy)) -
     .                qe*ey(ix,iy)
              geyyp = 2*gpiy(ix,iyp1,1)/(ney1(ix,iyp1)+ney0(ix,iyp1))-
     .                qe*ey(ix,iyp1)
              dgeyy0 = (geyy0-geyym)*eta1(ix,iy)*gy(ix,iy)
              dgeyy1 = (geyyp-geyy0)*eta1(ix,iyp1)*gy(ix,iyp1)
              vycf(ix,iy) = 2*(dgeyy1-dgeyy0)*gy(ix,iy) /
     .                                 ( (ney1(ix,iy)+ney0(ix,iy))*
     .                     (qe*0.5*(btot(ix,iy)+btot(ix,iym1)))**2 )
           endif

              diffusivwrk(ix,iy)=fcdif*difni(ifld)+dif_use(ix,iy,ifld)
  17       continue
  18     continue

c
c ... Compute diffusive part of radial velocity.
c .. Needs further cleaning; no turbulence model used now TDR 9/1/15
         do iy = j1, j5
            do ix = i1, i6
              difnimix = diffusivwrk(ix,iy)

c ... Alter diffusivity in the SOL by mixing fixed diffusivity
c     with anomalous diffusivity computed in subroutine turb_diffus but
c     reduced by the factor difnit(ifld).  The mixing ratio is given by
c     cdifnit.  Diffusivity is unaltered if difnit(ifld) = 0.
c...MER NOTE: For a full double-null configuration, the SOL is defined to
c...  be the region outside the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
cc              if (difnit(ifld) .gt. 1.e-20 .and. zi(ifld) .eq. 1.
cc     .                                 .and. iy .gt. iysptrx) then
cc                 difnimix = (1. - cdifnit) * 
cc     .                      (fcdif*difni(ifld) + dif_use(ix,iy,ifld)) +
cc     .                               cdifnit * difnit(ifld) * difnimix
cc              endif

              vydd(ix,iy,ifld) = vydd(ix,iy,ifld) 
     .           -1. * difnimix * (
     .            2*(1-isvylog)*( (niy1(ix,iy,ifld) - niy0(ix,iy,ifld)) /
     .              dynog(ix,iy) ) / (niy1(ix,iy,ifld)+niy0(ix,iy,ifld))+
     .              isvylog*(log(niy1(ix,iy,ifld)) - 
     .                            log(niy0(ix,iy,ifld))) /dynog(ix,iy) )

c ... Compute total radial velocity.
              vy(ix,iy,ifld) = cfydd *bfacyrozh(ix,iy) *
     .                                 vycp(ix,iy,ifld) + 
     .                         cfrd  * vyrd(ix,iy,ifld) +
     .                                 vydd(ix,iy,ifld) + 
     .                         cfyef * vyce(ix,iy,ifld) +
     .                         cfybf * vycb(ix,iy,ifld) +
     .                        cfvycf * vycf(ix,iy) + 
     .                        cfvycr * vycr(ix,iy)
c ... Compute radial vel v_grad_P eng eqn terms;cfydd+cfybf=1 or 0
              vygp(ix,iy,ifld) = (cfydd+cfybf)*bfacyrozh(ix,iy) *
     .                                         vycp(ix,iy,ifld) + 
     .                                 cfrd  * vyrd(ix,iy,ifld) +
     .                                         vydd(ix,iy,ifld) + 
     .                                 cfyef * vyce(ix,iy,ifld) +
     .                                cfvycf * vycf(ix,iy) + 
     .                                cfvycr * vycr(ix,iy)
              if (isybdrywd == 1) then  #make vy diffusive in wall cells
                 if (iy==0 .and. matwalli(ix) > 0) then
                    vy(ix,iy,ifld) = vydd(ix,iy,ifld)
                 elseif (iy==ny .and. matwallo(ix) > 0) then
                    vy(ix,iy,ifld) = vydd(ix,iy,ifld)
                 endif
              endif
            enddo  #loop over iy
         enddo     #loop over ix

	 do 20 iy = j1, j6
	    do 19 ix = i1, i6
	      iy1 = max(0,iy-1)            # does iy=0 properly
              iy2 = min(ny+1,iy+1) # use ex*fqx since phi(0,) may be large 
	      ix2 = ixp1(ix,iy)
	      ix4 = ixp1(ix,iy1)
              ix6 = ixp1(ix,iy2)
              do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))*gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                 elseif (ix==ixrb(jx) .and. ixmxbcl==1) then
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))*gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                 else  # not a boundary
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*
     .                                                       gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))* gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                endif
              enddo  #vis end do-loop over nxpt mesh regions
c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
              if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 ) then
                 temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                 temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*
     .                                                       gyc(ix,iy)
              endif

              v2ce(ix,iy,ifld) = - 0.5 * temp1
     .             / ( btot(ix,iy) + btot(ix2,iy) )
              v2cb(ix,iy,ifld) =(cfcurv*( 0.5*(tiv(ix,iy)+tiv(ix,iy1)) +
     .                 mi(ifld)*up(ix,iy,ifld)**2 )*curvrb2(ix,iy) +
     .                     cfgradb*0.5*( tiv(ix,iy)+tiv(ix,iy1) )*
     .                            gradb2(ix,iy))/qion
              ve2cb(ix,iy) = -(cfcurv*0.5*(tev(ix,iy)+tev(ix,iy1))*
     .                                                 curvrb2(ix,iy) +
     .                          cfgradb*0.5*(tev(ix,iy)+tev(ix,iy1))*
     .                            gradb2(ix,iy))/qe
              v2cd(ix,iy,ifld) = temp2
     .                    / ((btot(ix,iy)+btot(ix2,iy))*qion*
     .                           (ni(ix,iy,ifld)+ni(ix2,iy,ifld)))
              ve2cd(ix,iy,1) = -temp3
     .                    / ((btot(ix,iy)+btot(ix2,iy))*qe*
     .                           (ni(ix,iy,ifld)+ni(ix2,iy,ifld)))
              q2cd(ix,iy,ifld) = (priv(ix,iy,ifld)+priv(ix,iy1,ifld))*temp4
     .                    / ( (btot(ix,iy)+btot(ix2,iy))*qion )

c...  Calculate plate electr diamag flux used to find sheath potential
              do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                 # use ix=ixlb+1 values to avoid BC variations
                    v2dia = -0.5*( gpey(ixlb(jx)+1,iy)+gpey(ixlb(jx)+1,iy1) ) /
     .                          ( btot(ixlb(jx)+1,iy)*qe*ne(ixlb(jx)+1,iy) )
                    fdiaxlb(iy,jx) = ne(ixlb(jx)+1,iy) * sx(ixlb(jx),iy) *
     .                              v2dia * rbfbt(ixlb(jx)+1,iy)
                 endif
                 if (ix==ixrb(jx) .and. ixmxbcl==1) then
                    v2dia = -0.5*( gpey(ixrb(jx),iy)+gpey(ixrb(jx),iy1) ) /
     .                          ( btot(ixrb(jx),iy)*qe*ne(ixrb(jx),iy) )
                    fdiaxrb(iy,jx) = ne(ixrb(jx),iy) * sx(ixrb(jx),iy) *
     .                              v2dia * rbfbt(ixrb(jx),iy)
                 endif
              enddo  # end do-loop over nxpt mesh regions

              v2rd(ix,iy,ifld) = - 2. * gprx(ix,iy) /
     .           ( btot(ix,iy)/(etaper(ix,iy)*rbfbt2(ix,iy)) +
     .             btot(ix2,iy)/(etaper(ix2,iy)*rbfbt2(ix2,iy)) )
              v2dd(ix,iy,ifld) = - 2. * difpr2(ifld) * gprx(ix,iy) /
     .                                ( pr(ix2,iy)/rbfbt(ix2,iy) +
     .                                  pr(ix,iy)/rbfbt(ix,iy) ) -
     .            2. * (fcdif*difni2(ifld) + dif2_use(ix,iy,ifld)) * 
     .                               (ni(ix2,iy,ifld)-ni(ix,iy,ifld)) /
     .                      (ni(ix2,iy,ifld)/(rbfbt(ix2,iy)*gx(ix2,iy))+
     .                          ni(ix,iy,ifld)/(rbfbt(ix,iy)*gx(ix,iy)))
              v2(ix,iy,ifld) = cf2dd * bfacxrozh(ix,iy) *
     .                                 v2cd(ix,iy,ifld) + 
     .                         cfrd  * v2rd(ix,iy,ifld) +
     .                                 v2dd(ix,iy,ifld) + 
     .                         cf2ef * v2ce(ix,iy,ifld) +
     .                         cf2bf * v2cb(ix,iy,ifld)
c ...         Compute v2 for v2x_gradx_P eng terms; cf2dd+cf2bf=1 or 0
              v2xgp(ix,iy,ifld) =  0.5*(rbfbt(ix,iy)+rbfbt(ix2,iy)) * (
     .                 (cf2dd+cf2bf) * bfacxrozh(ix,iy) *
     .                                 v2cd(ix,iy,ifld) + 
     .                         cfrd  * v2rd(ix,iy,ifld) +
     .                                 v2dd(ix,iy,ifld) + 
     .                         cf2ef * v2ce(ix,iy,ifld) )
         if (isnonog.eq.1 .and. iy.le.ny) then
c            grdnv = ( 1/( fym (ix,iy,1)/ni(ix2,iy1,ifld) +
c     .                    fy0 (ix,iy,1)/ni(ix2,iy ,ifld) +
c     .                    fyp (ix,iy,1)/ni(ix2,iy2,ifld) +
c     .                    fymx(ix,iy,1)/ni(ix ,iy1,ifld) +
c     .                    fypx(ix,iy,1)/ni(ix, iy2,ifld) ) -
c     .                1/( fym (ix,iy,0)/ni(ix ,iy1,ifld) +
c     .                    fy0 (ix,iy,0)/ni(ix ,iy ,ifld) +
c     .                    fyp (ix,iy,0)/ni(ix ,iy2,ifld) + 
c     .                    fymx(ix,iy,0)/ni(ix4,iy1,ifld) +
c     .                    fypx(ix,iy,0)/ni(ix6,iy2,ifld) ) )
c     .                                                 / dxnog(ix,iy)
cc            grdnv = ( exp( fym (ix,iy,1)*log(ni(ix2,iy1,ifld)) + 
cc     .                     fy0 (ix,iy,1)*log(ni(ix2,iy ,ifld)) +
cc     .                     fyp (ix,iy,1)*log(ni(ix2,iy2,ifld)) +
cc     .                     fymx(ix,iy,1)*log(ni(ix ,iy1,ifld)) +
cc     .                     fypx(ix,iy,1)*log(ni(ix, iy2,ifld)) ) 
cc     .               -exp( fym (ix,iy,0)*log(ni(ix ,iy1,ifld)) + 
cc     .                     fy0 (ix,iy,0)*log(ni(ix ,iy ,ifld)) +
cc     .                     fyp (ix,iy,0)*log(ni(ix ,iy2,ifld)) +
cc     .                     fymx(ix,iy,0)*log(ni(ix4,iy1,ifld)) +
cc     .                     fypx(ix,iy,0)*log(ni(ix6,iy2,ifld)) ) ) /
cc     .                                                    dxnog(ix,iy)
            grdnv = (    ( fym (ix,iy,1)*log(ni(ix2,iy1,ifld)) + 
     .                     fy0 (ix,iy,1)*log(ni(ix2,iy ,ifld)) +
     .                     fyp (ix,iy,1)*log(ni(ix2,iy2,ifld)) +
     .                     fymx(ix,iy,1)*log(ni(ix ,iy1,ifld)) +
     .                     fypx(ix,iy,1)*log(ni(ix, iy2,ifld)) ) 
     .                  -( fym (ix,iy,0)*log(ni(ix ,iy1,ifld)) + 
     .                     fy0 (ix,iy,0)*log(ni(ix ,iy ,ifld)) +
     .                     fyp (ix,iy,0)*log(ni(ix ,iy2,ifld)) +
     .                     fymx(ix,iy,0)*log(ni(ix4,iy1,ifld)) +
     .                     fypx(ix,iy,0)*log(ni(ix6,iy2,ifld)) ) ) /
     .                                                      dxnog(ix,iy)
            vytan(ix,iy,ifld)=(fcdif*difni(ifld) + dif_use(ix,iy,ifld)) *
     .                                      (grdnv/cos(angfx(ix,iy)) - 
     .                       (log(ni(ix2,iy,ifld)) - log(ni(ix,iy,ifld)))
     .                                                 * gxf(ix,iy) )
            if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
              vytan(ix,iy,ifld) = 0.
            endif
            if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
c             non-physical interface between upper target plates for dnull
              vytan(ix,iy,ifld) = 0.
            endif
         endif

 19        continue
 20     continue

         do 21 ix = i1, i6
            vy(ix,ny+1,ifld) = 0.0   
 21      continue
        else    # test on zi > 1.e-10 to skip whole loop
        endif
  100 continue  # Giant loop over ifld (species)

c ... Save values returned by Hirshman mombal for Jacobian calc. to
c ... minimize calls - restore the "m" or ix-1 values at the end of pandf
c ... The Jacobian ix loop can then be reduced to only include ix-1 and ix
c ... Suffix "o" refers to "old" value at ix, and suffix "om" means "old" 
c ... value at ix-1.

      if (xc.ge.0 .and. yc.ge.0) then
         ix1 = ixm1(xc,yc)
         fqpom = fqp(ix1,yc)
         friceom = frice(ix1,yc)
         upeom = upe(ix1,yc)
         fqpo = fqp(xc,yc)
         friceo = frice(xc,yc)
         upeo = upe(xc,yc)
         do ifld = 1, nfsp
            friciom(ifld) = frici(ix1,yc,ifld)    # dimension req. nfsp<101
            upiom(ifld) = upi(ix1,yc,ifld)
            uupom(ifld) = uup(ix1,yc,ifld)
            fricio(ifld) = frici(xc,yc,ifld)
            upio(ifld) = upi(xc,yc,ifld)
            uupo(ifld) = uup(xc,yc,ifld)
         enddo
      endif


c ... Need to calculate new currents (fqp) after saving old & before frice,i
      if(isphion+isphiofft .eq. 1)  call calc_currents

c ... Add anomalous perp vis vy using calc_currents result - awkward,change
      if (cfvyavis > 0.) then
        do ifld = 1, 1  # nfsp  # only good for ifld=1
          do iy = max(j1,2), min(j5,ny-1)
            do ix = max(i1,2), min(i6,nx-1)
              vyavis(ix,iy,ifld) = fqya(ix,iy)*2/(
     .                  qe*(niy1(ix,iy,1)+niy0(ix,iy,1))*sy(ix,iy) )
              vy(ix,iy,ifld) = vy(ix,iy,ifld) + cfvyavis*vyavis(ix,iy,ifld)
            enddo
          enddo
        enddo
      endif          

c ... Calc friction forces from Braginskii; no individ chg-states;isimpon < 5.

      if (isimpon < 5) then
         do iy = j1, j6    #iys1, iyf6
            do ix = i1, i6
               ix2 = ixp1(ix,iy)
               nbarx = 0.5*(ne(ix,iy)+ne(ix2,iy))
               ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
               lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
               flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
               frice(ix,iy) = -cthe*flxlimf*nbarx*rrv(ix,iy)*gtex(ix,iy) +
     .                  cfnetap*qe*netap(ix,iy)*fqp(ix,iy)/sx(ix,iy)
               frici(ix,iy,1) = - frice(ix,iy)
               if (fac2sp .gt. 1.1 .and. nusp .eq. 2) then
                  frici(ix,iy,1) = - frice(ix,iy)/fac2sp
                  frici(ix,iy,2) = - frice(ix,iy)/fac2sp
               endif
            enddo
         enddo
      endif

c ... For use within subroutine mombal, the poloidal electric field is
c     calculated from || Ohms law if isphion = 0.  This field is not
c     intended for use in computing cross-field drifts, so a test of
c     cfyef is included. Both isphiofft=0 or 1 cases included in one loop

      if (isphion .eq. 0) then   # ex calc here assumes no parallel current
         do iy = iys1, iyf6
            do ix = i1, i6
               ix1 = ix
               do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                   ix1 = ixlb(jx) + 1
                 elseif (ix==ixrb(jx) .and. ixmxbcl==1) then
                   ix1 = ixrb(jx) - 1 
                 endif
               enddo
               ix2 = ixp1(ix1,iy)
               ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
               lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
               flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
               nexface = 0.5*(ne(ix2,iy)+ne(ix1,iy))
               ex(ix,iy) = (1 - isphiofft) * (
     .                                     -( gpex(ix1,iy)/nexface +
     .                                cthe*flxlimf*gtex(ix1,iy) )/qe -
     .                                              gpondpotx(ix,iy) +
     .                 pondomfpare_use(ix,iy)/(qe*rrv(ix,iy)*nexface) ) + 
     .                      isphiofft * (
     .                            (phi(ix1,iy)-phi(ix2,iy))*gxf(ix1,iy) )
            enddo
         enddo
      endif

c ... Loop over cells (restricted to poloidal slice of box if doing
c     Jacobian), calling mombal if it is providing parallel flow, and
c     taking poloidal projection of parallel flow to get poloidal flow.
c     Unperturbed values of the parallel-flow contribution to uu are
c     saved here so they can be restored below.



      do iy = iys1, iyf6
         if (xc .gt. 0) then
            ix = xc
            ix1 = ixm1(ix,iy)
            if (isimpon .eq. 5) then   # Hirshmans reduced-ion approx.
               if (istimingon .eq. 1) tsimp = gettime(sec4)
               call mombal (ix1,ix,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            elseif(isimpon .eq. 6 .or. isimpon .eq. 7) then # Force balance without inertia
               if (istimingon .eq. 1) tsimp = gettime(sec4)
               call mombalni (ix1,ix,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            endif
            do ifld = 1, nfsp
               if (ifld .le. nusp) then
                 upi(ix1,iy,ifld) = up(ix1,iy,ifld)
               else
                 do jx = 1, nxpt
                    if ( (ix1==ixlb(jx).and.ixmnbcl==1) .or.
     .                   (ix1==ixrb(jx).and.ixmxbcl==1) ) then
                       # constrain boundary velocity
                       if (zi(ifld) .gt. 1.e-10) then
                          argx = abs((2-2*upi(ix1,iy,ifld)/
     .                                     (upi(ix1,iy,1)+cutlo3))**3)
                          argx = min(20., argx)
                          upi(ix1,iy,ifld) = upi(ix1,iy,1) + 
     .                       (upi(ix1,iy,ifld)-upi(ix1,iy,1))*exp(-argx)
                       endif  # end if-test on zi
                    endif  # end if-test on ix
                 enddo  # end do-loop over nxpt mesh regions
               endif
               uup(ix1,iy,ifld) = rrv(ix1,iy)*upi(ix1,iy,ifld)
            enddo
         endif
         do ix = ixs1, min(ixf6, nx+1-ixmxbcl)
            ix2 = ixp1(ix,iy)
            if (isimpon .eq. 5) then
               if (istimingon .eq. 1) tsimp = gettime(sec4)
               call mombal (ix,ix2,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            elseif(isimpon .eq. 6 .or. isimpon .eq. 7) then # Force balance without inertia
               if (istimingon .eq. 1) tsimp = gettime(sec4)
               call mombalni (ix,ix2,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            endif
            do ifld = 1, nfsp
               if (ifld .le. nusp) then
                 upi(ix,iy,ifld) = up(ix,iy,ifld)
               else
                 do jx = 1, nxpt
                    if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                   (ix==ixrb(jx).and.ixmxbcl==1) ) then
                       # constrain boundary velocity
                       if (zi(ifld) .gt. 1.e-10) then
                          argx = abs((2-2*upi(ix,iy,ifld)/
     .                                     (upi(ix,iy,1)+cutlo3))**3)
                          argx = min(20., argx)
                          upi(ix,iy,ifld) = upi(ix,iy,1) + 
     .                         (upi(ix,iy,ifld)-upi(ix,iy,1))*exp(-argx)
                       endif  # end if-test on zi
                    endif  # end if-test on ix
                 enddo  # end do-loop over nxpt mesh regions
               endif
               uup(ix,iy,ifld) = rrv(ix,iy)*upi(ix,iy,ifld)
            enddo
         enddo
      enddo

c ... Add contributions to poloidal velocity from cross-field drifts
c     to those from parallel flow.

      do ifld = 1, nfsp
         do iy = j1, j6
            if (i1 .gt. 0) then  # il is initial ix; here for uu(ixm1(i1-1,,)
               ix = i1
               ix1 = ixm1(ix,iy)
               uu(ix1,iy,ifld) = uup(ix1,iy,ifld) +
     .                           0.5 * (rbfbt(ix,iy) + rbfbt(ix1,iy)) *
     .                           v2(ix1,iy,ifld) - vytan(ix1,iy,ifld) -
     .                         difax(ifld) * 0.5 * ( ( 0.5*(
     .                         ni(ix1,iy,ifld)/ni(ix,iy,ifld) +
     .                         ni(ix,iy,ifld)/ni(ix1,iy,ifld)) -1)**2 ) *
     .                        (ni(ix,iy,ifld)-ni(ix1,iy,ifld))*gxf(ix1,iy)
     .                       /(ni(ix,iy,ifld)+ni(ix1,iy,ifld))
               uz(ix1,iy,ifld) = -uup(ix1,iy,ifld)/rrv(ix1,iy)*
     .          0.5*(rbfbt(ix,iy) + rbfbt(ix1,iy)) + sign(1.,b02d(ix,iy))*
     .               (cftef*v2ce(ix1,iy,ifld)+cftdd*v2cd(ix1,iy,ifld))*
     .                                                       rrv(ix1,iy)
            endif
            do ix = i1, i6    #now the remainder of the uu(ix,,)
               ix2 = ixp1(ix,iy)
               uu(ix,iy,ifld) = uup(ix,iy,ifld) +
     .                          0.5 * (rbfbt(ix,iy) + rbfbt(ix2,iy)) *
     .                          v2(ix,iy,ifld) - vytan(ix,iy,ifld) -
     .                         difax(ifld) * 0.5 * ( ( 0.5*(
     .                         ni(ix,iy,ifld)/ni(ix2,iy,ifld) +
     .                         ni(ix2,iy,ifld)/ni(ix,iy,ifld)) -1)**2 ) *
     .                        (ni(ix2,iy,ifld)-ni(ix,iy,ifld))*gxf(ix,iy)
     .                       /(ni(ix2,iy,ifld)+ni(ix,iy,ifld))
               uz(ix,iy,ifld) = -uup(ix,iy,ifld)/rrv(ix,iy)*
     .          0.5*(rbfbt(ix,iy) + rbfbt(ix2,iy)) + sign(1.,b02d(ix,iy))*
     .               (cftef*v2ce(ix,iy,ifld)+cftdd*v2cd(ix,iy,ifld))*
     .                                                       rrv(ix,iy)
            enddo
         enddo
      enddo

c...  If upi not from full ||mom eq (e.g.,isimpon=6), set impurity
c...  uu(ixrb,,) & upi(ixrb,,) via generalized Bohm cond.
      if(isimpon > 0) then
        do jx = 1, nxpt
	  ixt0 = ixlb(jx)
          ixt = ixrb(jx)+1 
          ixt1 = ixrb(jx)
	  do ifld = nhsp+1, nfsp
            if(ifld > nusp) then  #species without full mom eqn
	      do iy = j1, j6
c ..          first left plate(s)
                if(isfixlb(jx) == 0) then # set upi for left plate
                  cs = csfacrb(ifld,jx)*sqrt( (te(ixt0,iy) +
     .                            csfacti*ti(ixt0,iy))/mi(ifld) )
                  ueb = cfueb*( cf2ef*v2ce(ixt0,iy,ifld)*rbfbt(ixt0,iy) -
     .                            vytan(ixt0,iy,ifld) )/rrv(ixt0,iy)
	          uu(ixt0,iy,ifld) = -rrv(ixt0,iy)*cs
                  upi(ixt0,iy,ifld) = -(cs - ueb)
                endif
c ..          switch to right plate(s)
                if(isfixrb(jx) == 0) then
                  cs = csfacrb(ifld,jx)*sqrt( (te(ixt1,iy) +
     .                            csfacti*ti(ixt1,iy))/mi(ifld) )
                  ueb = cfueb*( cf2ef*v2ce(ixt1,iy,ifld)*rbfbt(ixt,iy) -
     .                            vytan(ixt1,iy,ifld) )/rrv(ixt1,iy)
	          uu(ixt1,iy,ifld) = rrv(ixt1,iy)*cs
	          uu(ixt,iy,ifld) = uu(ixt1,iy,ifld)
                  upi(ixt1,iy,ifld) = cs - ueb
                  upi(ixt,iy,ifld) = upi(ixt1,iy,ifld)
                endif
              enddo
            endif   #checks if ifld > nusp
          enddo
        enddo
      endif         # checks if isimpon > 0

************************************************************************
*     Calculate the currents fqx, fqy, fq2 and fqp, if isphion = 1
*     or if isphiofft = 1.
************************************************************************
ccc      if(isphion+isphiofft .eq. 1)  call calc_currents

************************************************************************
*     Calculate the electron velocities, vex, upe, ve2, vey
************************************************************************
      do 25 iy = j1, j6
	 do 24 ix = i1, i6
	    vex(ix,iy) = 0.
	    vey(ix,iy) = 0.
   24    continue
   25 continue

      if (isimpon.eq.5) goto 29    # have upe from mombal

      do iy = j1, j6    #iys1, iyf6
         do ix = i1, i6
            upe(ix,iy) = 0.
         enddo
      enddo

      do 27 ifld = 1, nfsp
         do iy = j1, j6    #iys1, iyf6
	    do ix = i1, i6
               ix1 = ixp1(ix,iy)
	       upe(ix,iy) = upe(ix,iy) + upi(ix,iy,ifld)*zi(ifld)*0.5*
     .                      ( ni(ix,iy,ifld)+ni(ix1,iy,ifld) )
            enddo
         enddo
   27 continue
      afqp = 1.
      if (isimpon.eq.6 .or. isimpon.eq.7) afqp = fupe_cur  #allows gradual fix for old cases
      do iy = j1, j6    #iys1, iyf6
         do ix = i1, i6
            ix1 = ixp1(ix,iy)
	    upe(ix,iy) = (upe(ix,iy) -afqp*fqp(ix,iy)/
     .                               (rrv(ix,iy)*sx(ix,iy)*qe))/
     .                             (0.5*( ne(ix,iy)+ne(ix1,iy) ))
         enddo
      enddo

  29  continue

      do 731 iy = j1, j6   # ExB same all species;if cf2dd=1, no imp yet
	 do 730 ix = i1, i6
            ix1 = ixp1(ix,iy)
            vex(ix,iy) = upe(ix,iy)*rrv(ix,iy) + 
     .                   (cf2ef*v2ce(ix,iy,1) + cf2bf*ve2cb(ix,iy) + 
     .                         cf2dd*bfacxrozh(ix,iy)*ve2cd(ix,iy,1) ) *
     .                           0.5*(rbfbt(ix,iy) + rbfbt(ix1,iy)) -
     .                                               vytan(ix,iy,1) 
     
  730    continue
  731 continue

      do 734 ifld = 1, nfsp
	 do 733 iy = j1, j5
	    do 732 ix = i1, i6   # grad_B will be ok as next fqy is subtr.
	       vey(ix,iy) = vey(ix,iy) + vy(ix,iy,ifld)*zi(ifld)*0.5*
     .                      ( niy0(ix,iy,ifld)+niy1(ix,iy,ifld) )
  732       continue
  733    continue
  734 continue

      do 36 iy = j1, j5
	 do 35 ix = i1, i6
	    vey(ix,iy) = (vey(ix,iy)-cfjve*fqy(ix,iy)/(sy(ix,iy)*qe))/
     .                    (0.5*( ney0(ix,iy)+ney1(ix,iy) ))
   35    continue
   36 continue
	 
c ... if isnewpot=0, vey(,0) needs to be redone since fqy(,0)=0
      if (isnewpot==1) then
        do ix = i1, i6  # ExB vyce same all species
          vey(ix,0) = cfybf*veycb(ix,0) + vydd(ix,0,1) +
     .                cfyef*vyce(ix,0,1)
        enddo
      endif

c ... If isybdrywd = 1, make vey diffusive, just like vy
      if (isybdrywd == 1) then  #make vy diffusive in wall cells
        do ix = i1, i6
          if (matwalli(ix) > 0) vey(ix,0)  = vydd(ix,0,1)
          if (matwallo(ix) > 0) vey(ix,ny) = vydd(ix,ny,1)
        enddo
      endif
      

************************************************************************
*   We Calculate the source terms now.
************************************************************************
*  ---------------------------------------------------------------------
*  Coefficients for the source terms.
*  ---------------------------------------------------------------------

      do 702 iy = j2, j5
         do 701 ix = i2, i5
            do ifld = 1, nfsp
               snic(ix,iy,ifld) = 0.0
               sniv(ix,iy,ifld) = 0.0
               psori(ix,iy,ifld) = 0.0
            enddo
            do ifld = 1, nusp
               smoc(ix,iy,ifld) = 0.0
               smov(ix,iy,ifld) = 0.0
            enddo
            seec(ix,iy) = 0.0
            seev(ix,iy) = 0.0
            seic(ix,iy) = 0.0
            seiv(ix,iy) = 0.0
	    psorbgz(ix,iy) = 0.    # diagnostic only
  701    continue
  702 continue

************************************************************************
*  -- steady sources
************************************************************************
*  ---------------------------------------------------------------------
*  volume sources. (old PHYSRC)
*  ---------------------------------------------------------------------

c ... Calculate effective Lyman-alpha optical depth factors used if
c ... istabon=14 or 15 for hydr. rate look-ups set rtauxfac<=0 to bypass
c ----------------------------------------------------------------------
      if (rtauxfac .gt. 0.) then

CC .. FIRST GET THE POLOIDAL OPTICAL DEPTH FACTOR
c ------------------------------------------------
         do iy = 0, ny+1

c ... get optical-depth to left (ix=0) boundary
            rdumx = 0.
            do ix = 0, nx+1
               dr1  = 0.5*dx(ix,iy)
               dr2  = 0.5*dx(ix,iy)
               if(ix==(nx+1)/2) rdumx=0.  #make two plate regions independent
               rdumx = rdumx + rtauxfac*dr1*ng(ix,iy,1)
               rtaux(ix,iy) = rdumx*rt_scal
               rdumx = rdumx + rtauxfac*dr2*ng(ix,iy,1)
            enddo

c ... get optical-depth to right (ix=nx+1) bdry; initial selection of min rtaux
            rdumx = 0.
            do ix = nx+1, 0, -1
               dr1  = 0.5*dx(ix,iy)
               dr2  = 0.5*dx(ix,iy)
               if(ix==(nx+1)/2) rdumx=0.  #make two plate regions independent
               rdumx = rdumx + rtauxfac*dr1*ng(ix,iy,1)
               rtaux(ix,iy) = min(rdumx*rt_scal, rtaux(ix,iy))
               rdumx = rdumx + rtauxfac*dr2*ng(ix,iy,1)
            enddo
         enddo   # iy loop

CC .. NOW GET THE RADIAL OPTICAL DEPTH FACTOR
c --------------------------------------------
         do ix = 0, nx+1

c ... get optical-depth to inner (iy=0) bdry
            rdumy = 0.
            do iy = 0, ny+1
               dr1  = 0.5*dy(ix,iy)
               dr2  = 0.5*dy(ix,iy)
               rdumy = rdumy + rtauyfac*dr1*ng(ix,iy,1)
               rtauy(ix,iy) = rdumy*rt_scal
               rdumy = rdumy + rtauyfac*dr2*ng(ix,iy,1)
            enddo

c ... get optical-depth to outer (iy=ny+1) bdry; selection of min rtau
            rdumy = 0.
            do iy = ny+1, 0, -1
               dr1  = 0.5*dy(ix,iy)
               dr2  = 0.5*dy(ix,iy)
               rdumy = rdumy + rtauyfac*dr1*ng(ix,iy,1)
               rtauy(ix,iy) = min(rdumy*rt_scal, rtauy(ix,iy))
               rtau(ix,iy) = min(rtaux(ix,iy), rtauy(ix,iy))
               rdumy = rdumy + rtauyfac*dr2*ng(ix,iy,1)
            enddo
         enddo   # ix loop

      endif     # test on rtauxfac, skip if rtauxfac is negative

*     The following is a temporary recycling model.

*  -- recalculate particle source psor if ifixpsor=0 --

c...  Initialize save-variables if this is a Jacobian (xc,yc > -1)
         if (xc .ge. 0 .and. yc .ge. 0) then
            psordisold = psordis(xc,yc)
cc            write(*,*) 'Just after psordisold; xc,yc=',xc,yc
            do ifld = 1, nfsp
               psorold(ifld) = psorc(xc,yc,ifld)
               psorxrold(ifld) = psorxr(xc,yc,ifld)
               msorold(ifld) = msor(xc,yc,ifld)
               msorxrold(ifld) = msorxr(xc,yc,ifld)
               nucxiold(ifld) = nucxi(xc,yc,ifld)
               nueliold(ifld) = nueli(xc,yc,ifld)
            enddo
            do igsp = 1, ngsp
               nucxold(igsp) = nucx(xc,yc,igsp)
               nurcold(igsp) = nurc(xc,yc,igsp)
               nuizold(igsp) = nuiz(xc,yc,igsp)
               nuixold(igsp) = nuix(xc,yc,igsp)
               nuelgold(igsp) = nuelg(xc,yc,igsp)
               psorgold(igsp) = psorgc(xc,yc,igsp)
               psorrgold(igsp) = psorrgc(xc,yc,igsp)
               psorcxgold(igsp) = psorcxgc(xc,yc,igsp)
            enddo
         endif

c...  The particle source can be frozen if ifixpsor.ne.0
      if(ifixpsor .eq. 0) then
            
        igsp = 0
        do ifld = 1, nhsp  # Hydrogen-only loop
         if (zi(ifld) > 0.) then  #calc only for hydrogen ions
          igsp = igsp + 1
          do iy = iys1, iyf6
            do ix = ixs1, ixf6

c     Ionization of neutral hydrogen by electrons and recombination--
               if (icnuiz .eq. 0) then
                  ne_sgvi = ne(ix,iy)
                  if (ifxnsgi.eq.1) ne_sgvi = cne_sgvi  # fix density dependence
                  nuiz(ix,iy,igsp) = chioniz *  ne(ix,iy) * (
     .                           rsa(te(ix,iy),ne_sgvi,rtau(ix,iy),0)
     .                         + sigvi_floor )
                  if (xc .ge. 0) then        # limit Jacobian element
                     nuiz(ix,iy,igsp) = fnnuiz*nuiz(ix,iy,igsp) +
     .                               (1-fnnuiz)*nuizold(igsp)
                  endif
               elseif (icnuiz .eq. 1) then
                  nuiz(ix,iy,igsp) = cnuiz
               endif
               if (isrecmon == 1) then
                  nurc(ix,iy,igsp) = cfrecom * ne(ix,iy) 
     .                         * rra(te(ix,iy),ne(ix,iy),rtau(ix,iy),1)
                  if (xc .ge. 0) then        # limit Jacobian element
                     nurc(ix,iy,igsp) = fnnuiz*nurc(ix,iy,igsp) +
     .                             (1-fnnuiz)*nurcold(igsp)
                  endif
               else
                   nurc(ix,iy,igsp) = 0.
               endif
               psorbgg(ix,iy,igsp) = ngbackg(igsp)*( (0.9 + 0.1*
     .                            (ngbackg(igsp)/ng(ix,iy,igsp))**ingb) ) * 
     .                             nuiz(ix,iy,igsp) * vol(ix,iy)
               psorgc(ix,iy,igsp) = -ng(ix,iy,igsp)*nuiz(ix,iy,igsp)*vol(ix,iy) +
     .                              psorbgg(ix,iy,igsp)
               psorc(ix,iy,ifld) = - psorgc(ix,iy,igsp)
               psordis(ix,iy) = cfdiss*psorc(ix,iy,1)  # overwritten below if ishymol=1
               psorxrc(ix,iy,ifld) = -ni(ix,iy,ifld)*nurc(ix,iy,igsp)*vol(ix,iy)
               psorrgc(ix,iy,igsp) = -psorxrc(ix,iy,ifld)
               msor(ix,iy,ifld) = 0.
               msorxr(ix,iy,ifld) = 0.


c     Charge exchange on neutral hydrogen --
              if (icnucx .eq. 0) then
	         t0 = max(ti(ix,iy),temin*ev)
ccc   we omit the weak velocity dependence as it brings in ni(ix+1) in Jac
                 t1 = t0/(mi(ifld)/mp)
                 nucx(ix,iy,igsp) = ni(ix,iy,ifld) * rcx(t1,ni(ix,iy,ifld),1)
              elseif (icnucx .eq. 1) then
                 nucx(ix,iy,igsp) = cnucx
              elseif (icnucx == 2) then
	         t0 = max(ti(ix,iy),temin*ev)
                 nucx(ix,iy,igsp) = sqrt(t0/mi(ifld))*
     .                         sigcx*(ni(ix,iy,ifld)+rnn2cx*ng(ix,iy,igsp))
              endif
              nuix(ix,iy,igsp) = fnuizx*nuiz(ix,iy,igsp) + 
     .                           fnucxx*nucx(ix,iy,igsp)
                              #dont use neutral-neutral collisions here
c
c   neutral particle source/sink for isupgon=1 (reset below if multispecies
c   models are on [isimpon = 5 or 6 or 7])
              if(isupgon(igsp) .eq. 1)then # inertia gas species is ifld+1
                 psorc(ix,iy,ifld+1)= -psorc(ix,iy,ifld)
                 psorxrc(ix,iy,ifld+1)= -psorxrc(ix,iy,ifld)
                 msor(ix,iy,ifld+1)= 0.
                 msorxr(ix,iy,ifld+1)= 0.
              endif
c
            enddo   #end loop over ix
          enddo     #end loop over iy
         endif      #omit whole loop if zi(ifld) = 0. (neutrals)
        enddo       #end loop over hydrogen species (ifld)

c*****************************************************************
c ... Average psorgc and psorc over cell vol with simple 5pt ave
c*****************************************************************
        igsp = 0
        do ifld = 1, nhsp  # Hydrogen-only loop
         if (zi(ifld) > 0.) then  #calc only for hydr ions, not neuts
           igsp = igsp + 1
           if (ispsorave.eq.0.) then  #use only single-cell value
             do iy = iys1, iyf6
               do ix = ixs1, ixf6
                 psorg(ix,iy,igsp) = psorgc(ix,iy,igsp)
                 psor(ix,iy,ifld) =  psorc(ix,iy,ifld)
                 psorxr(ix,iy,ifld) = psorxrc(ix,iy,ifld)
                 psorrg(ix,iy,igsp) = psorrgc(ix,iy,igsp)
                 if(isupgon(igsp) .eq. 1) then # inertia gas is ifld+1
                   psor(ix,iy,ifld+1)= -psor(ix,iy,ifld)
                   psorxr(ix,iy,ifld+1)= -psorxr(ix,iy,ifld)
                 endif
               enddo
             enddo

           elseif (ispsorave > 0.) # use 5pt ave; first divide by vol

             if (xc < 0) then  #full RHS eval
               j2pwr = j2
               j5pwr = j5
             else  # Jacobian eval
               j2pwr = max(1, yc-1)
               j5pwr = min(ny, yc+1)
             endif 
             do iy = j2pwr, j5pwr
               if (xc < 0) then #full RHS eval
                 i2pwr = i2
                 i5pwr = i5
               else  #Jacobian eval
                 i2pwr = max(1,ixm1(xc,yc))
                 i5pwr = min(nx, ixp1(xc,yc))
               endif
               do ix = i2pwr, i5pwr
                 ix1 = ixm1(ix,iy)
                 ix2 = ixp1(ix,iy)
                 psorg(ix,iy,igsp) = (1.-ispsorave*0.5)*
     .                                  psorgc(ix,iy,igsp)+ 
     .                               0.125*ispsorave*vol(ix,iy)*
     .                          ( psorgc(ix,iy-1,igsp)/vol(ix,iy-1) + 
     .                            psorgc(ix,iy+1,igsp)/vol(ix,iy+1) +
     .                            psorgc(ix1,iy,igsp)/vol(ix1,iy)   + 
     .                            psorgc(ix2,iy,igsp)/vol(ix2,iy) )
                 psorxr(ix,iy,ifld) = (1.-ispsorave*0.5)*
     .                                 psorxrc(ix,iy,ifld) + 
     .                                  0.125*ispsorave*vol(ix,iy)*
     .                           ( psorxrc(ix,iy-1,ifld)/vol(ix,iy-1) + 
     .                             psorxrc(ix,iy+1,ifld)/vol(ix,iy+1) +
     .                             psorxrc(ix1,iy,ifld)/vol(ix1,iy)   + 
     .                             psorxrc(ix2,iy,ifld)/vol(ix2,iy) )
                 psor(ix,iy,ifld) = -psorg(ix,iy,igsp)
                 psorrg(ix,iy,igsp) = -psorxr(ix,iy,ifld)
                 if(isupgon(igsp) .eq. 1) then # inertia gas is ifld+1
                   psor(ix,iy,ifld+1)= -psor(ix,iy,ifld)
                   psorxr(ix,iy,ifld+1)= -psorxr(ix,iy,ifld)
                 endif

               enddo   #end loop over ix
             enddo     #end loop over iy
           endif       #if-loop on ipsorave
         endif         #omit whole loop if zi(ifld) = 0. (neutrals)
        enddo          #end loop over hydrogen species (ifld)

c ... Can now calc current from nucx since it is updated
      if (cfqyn .gt. 0.) call calc_curr_cx

c ... Ionization and recombination of impurities.
c     (Note that charge-exchange recombination is implemented for
c     impurities but the corresponding terms need to be added to
c     equations for hydrogenic-ion and gas continuity and electron
c     and ion energy.) ???
c     The total source is written as psor+psorxr, where
c     psor = n_(z-1) ne K^i_(z-1) - n_z ne K^i_z     # ionization gain/loss
c     psorxr = -n_z[ne K^r_z+ng K^cx_z]              # cx/r loss from z to z-1
c              +n_(z+1)[ne K^r_(z+1)+ng K^cx_(z+1)]  # cx/r gain to z from z+1

        if (isimpon .ge. 5 .and. nzspt .gt. 0) then
          do 604 iy = iys1, iyf6
             do 603 ix = ixs1, ixf6

                  if (istimingon .eq. 1) tsimp = gettime(sec4)
                  nevol = ne(ix,iy) * vol(ix,iy)
                  ngvol = ng(ix,iy,1) * vol(ix,iy)

                  jg = nhgsp
                  ifld_lcs = nhsp
                  do jz = 1, ngspmx-1           # for all impurity isotopes
                     if (nzsp(jz)==0) break
                     ifld_fcs = ifld_lcs + 1
                     ifld_lcs = ifld_fcs + nzsp(jz) - 1
                     if (ngsp .gt. nhgsp) then  # impurity gas is present
                         jg = jg + 1
                         if (ismctab .eq. 1) then  # rates for Z=0 gas
                            call imprates(te(ix,iy), 0, nzsp(jz), kionz0,
     .                                    krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy),te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   0, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz0, krecz, kcxrz)
                         endif
                         kionz0 = kionz0 + sigvi_floor
			 psorbgg(ix,iy,jg)= ngbackg(jg)*
     .                     (0.9+0.1*(ngbackg(jg)/ng(ix,iy,jg))**ingb) * 
     .                                                      nevol*kionz0
                         psorg(ix,iy,jg) = -ng(ix,iy,jg)*nevol*kionz0 +
     .                                      psorbgg(ix,iy,jg)
                         psor(ix,iy,ifld_fcs) = - psorg(ix,iy,jg)
                         msor(ix,iy,ifld_fcs)= 0.  # zero gas mom. assumed
                         if (ismctab .eq. 1) then  # rates for Z=1 ions
                            call imprates(te(ix,iy), 1, nzsp(jz), kionz,
     .                                    krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy), te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   1, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz, krecz, kcxrz)
                            kcxrz = cfcximp1*kcxrz  # rescale if desired
                         endif
                         kionz = kionz + sigvi_floor # only to set kionm below
                         kcxrzig = rcxighg(jg)*kcxrz  # K_cx of ng(jg)+ni(1)->
                         niz_floor = nzbackg(ifld_fcs) * (0.9 + 0.1*
     .                          (nzbackg(ifld_fcs)/ni(ix,iy,ifld_fcs))**inzb)
                         pscx0 = ngvol*(ni(ix,iy,ifld_fcs)-niz_floor)*kcxrz - 
     .                           ng(ix,iy,jg)*ni(ix,iy,1)*vol(ix,iy)*
     .                                                        kcxrzig
                         psorcxg(ix,iy,jg) = pscx0
                         psorcxg(ix,iy,1) = -pscx0
                         psorrg(ix,iy,jg) = nevol*(ni(ix,iy,ifld_fcs)-
     .                                                   niz_floor)*krecz
                         psorxr(ix,iy,ifld_fcs)= -psorrg(ix,iy,jg) - pscx0 
                         psorxr(ix,iy,1) = psorxr(ix,iy,1) + pscx0
cc                    Note: summed over ion/neutrals here backgrd source=0
                         massfac = cfmassfac*16*mi(1)/(3*(mg(jg)+mi(1)))
                         nuiz(ix,iy,jg) = kionz0*ne(ix,iy)
                         nuix(ix,iy,jg) = fnuizx*nuiz(ix,iy,jg) + 
     .                                   kcxrzig*ni(ix,iy,1) +
     .                         massfac*( kelighi(jg)*ni(ix,iy,1) +
     .                                   kelighg(jg)*ng(ix,iy,1) )
                         nucxi(ix,iy,ifld_fcs) = sigcxms(ifld_fcs,jg)*
     .                          sqrt(ti(ix,iy)/mi(ifld_fcs))*ng(ix,iy,jg)
                         nucx(ix,iy,jg) = sigcxms(ifld_fcs,jg)*
     .                        sqrt(ti(ix,iy)/mi(ifld_fcs))*ni(ix,iy,ifld_fcs)
                         massfac = cfmassfac*16*mi(ifld_fcs)/
     .                                               (3*(mg(jg)+mi(ifld_fcs)))
                         nueli(ix,iy,ifld_fcs) = massfac*( keligii(jg)*
     .                                                    ng(ix,iy,jg) )
                         nuelg(ix,iy,jg) = massfac*( keligii(jg)*
     .                                                    ni(ix,iy,ifld_fcs) )

                     else                       # no impurity gas present
                         izch = nint(zi(ifld_fcs))
                         if (ismctab .eq. 1) then
                            call imprates(te(ix,iy), izch, nzsp(jz), 
     .                                   kionz, krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy), te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   izch, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz, krecz, kcxrz)
                         endif
                         kionz = kionz + sigvi_floor
                         psor(ix,iy,ifld_fcs) = 0.   # in case nzspt=1
                         psorxr(ix,iy,ifld_fcs) = 0.
                         msor(ix,iy,ifld_fcs) = 0.   # in case nzspt=1
                         msorxr(ix,iy,ifld_fcs) = 0.
                         krecz = 0.
                         kcxrz = 0.
                     endif   # end if-branches for Z=1 with/wo impurity gas
                     if (znucl(ifld_fcs)==1) then
                         # hydrogenic impurity: include cx on ifld_fcs in nuix
                         nuix(ix,iy,jg) = nuix(ix,iy,jg) + 
     .                                   kcxrz*ni(ix,iy,ifld_fcs)
                         nuix(ix,iy,1) = nuix(ix,iy,1) + 
     .                                   kcxrzig*ni(ix,iy,ifld_fcs)
                     endif

                     kionm = kionz     # set values as previous charge state
                     krecm = krecz
                     kcxrm = kcxrz

                     do ifld = ifld_fcs + 1, ifld_lcs  # for charge states Z > 1
                        izch = nint(zi(ifld))
                        if (ismctab .eq. 1) then
                           call imprates(te(ix,iy), izch, nzsp(jz), 
     .                                  kionz, krecz, kcxrz)
                        elseif (ismctab .eq. 2) then
                           call mcrates(ne(ix,iy), te(ix,iy),
     .                                  ti(ix,iy)*mp/mi(1),
     .                                  izch, nzsp(jz), znucl(ifld),
     .                                  kionz, krecz, kcxrz)
                            kcxrz = cfcximp2*kcxrz   #rescale if desired
                        endif
                        kionz = kionz + sigvi_floor
			if (ifld==ifld_lcs) kionz = 0. #ensure no lcs ioniz
                        pxri = 0.    # gets reset if ifld.eq.ifld_fcs+1
                        z1fac = 1.   # gets reset = 0 if ifld.eq.ifld_fcs+1

                        if (ifld .eq. ifld_fcs+1) then #for 2nd charge-state
                           nizm_floor = nzbackg(ifld_fcs) * (0.9 + 0.1*
     .                         (nzbackg(ifld_fcs)/ni(ix,iy,ifld_fcs))**inzb)
                           psor(ix,iy,ifld_fcs) = psor(ix,iy,ifld_fcs)-
     .                            nevol*(ni(ix,iy,ifld_fcs)-nizm_floor)*
     .                                                            kionm
			   psorbgz(ix,iy) = psorbgz(ix,iy)-nevol*nizm_floor*
     .                                                            kionm
	                   msor(ix,iy,ifld_fcs) = msor(ix,iy,ifld_fcs)-
     .                            nevol*(ni(ix,iy,ifld_fcs))*
     .                            kionm*mi(ifld_fcs)*up(ix,iy,ifld_fcs)    
                           pxri = psorxr(ix,iy,ifld_fcs) #set in Z=1 loop
                           z1fac = 0.
                        endif

                        niz_floor = nzbackg(ifld) * (0.9 + 0.1*
     .                            (nzbackg(ifld)/ni(ix,iy,ifld))**inzb)
                        psor(ix,iy,ifld) = nevol *
     .                                     ( ni(ix,iy,ifld-1) * kionm -
     .                             (ni(ix,iy,ifld)-niz_floor) * kionz )
			psorbgz(ix,iy) = psorbgz(ix,iy)-nevol*niz_floor*kionz
                        msor(ix,iy,ifld) = nevol *
     .                                  ( ni(ix,iy,ifld-1) * kionm *
     .                                   mi(ifld-1) * up(ix,iy,ifld-1) -
     .                            (ni(ix,iy,ifld)) * kionz *
     .                                     mi(ifld) * up(ix,iy,ifld) )
                        psorxr(ix,iy,ifld-1) = pxri-z1fac*(nevol*krecm +
     .                                                     ngvol*kcxrm) *
     .                                    (ni(ix,iy,ifld-1)-nizm_floor) +
     .                                   (nevol*krecz + ngvol*kcxrz) *
     .                                                   ni(ix,iy,ifld)
			psorbgz(ix,iy) = psorbgz(ix,iy) + z1fac*
     .                                   (nevol*krecm + ngvol*kcxrm) *
     .                                   nizm_floor
                        msorxr(ix,iy,ifld-1) = 0. - (nevol*krecm +
     .                                               ngvol*kcxrm) *
     .                                    (ni(ix,iy,ifld-1)) *
     .                                      mi(ifld-1)*up(ix,iy,ifld-1) +
     .                        (nevol*krecz + ngvol*kcxrz)*ni(ix,iy,ifld)*
     .                                            mi(ifld)*up(ix,iy,ifld)
                        psorxr(ix,iy,1) = psorxr(ix,iy,1) + ngvol*
     .                                               ni(ix,iy,ifld)*kcxrz
                        psorcxg(ix,iy,1) = psorcxg(ix,iy,1) - ngvol*
     .                                               ni(ix,iy,ifld)*kcxrz
                        nucxi(ix,iy,ifld) = sigcxms(ifld,jg)*
     .                              sqrt(ti(ix,iy)/mi(ifld))*ng(ix,iy,jg)
                        nucx(ix,iy,jg) = nucx(ix,iy,jg) + sigcxms(ifld,jg)*
     .                          sqrt(ti(ix,iy)/mi(ifld))*ni(ix,iy,ifld)
                        massfac = cfmassfac*16*mi(ifld)/
     .                                              (3*(mg(jg)+mi(ifld)))
                        nueli(ix,iy,ifld) = massfac*( keligii(jg)*
     .                                                   ng(ix,iy,jg) )
                        nuelg(ix,iy,jg) = nuelg(ix,iy,jg) + massfac*
     .                                  ( keligii(jg)*ni(ix,iy,ifld) )

                        kionm = kionz
                        krecm = krecz
                        kcxrm = kcxrz
                        nizm_floor = niz_floor
                        if (ifld .eq. ifld_lcs) then  # last charge-state
                          psorxr(ix,iy,ifld) = -(nevol * krecz +
     .                                           ngvol * kcxrz) *
     .                                          (ni(ix,iy,ifld)-niz_floor) 
			  psorbgz(ix,iy) = psorbgz(ix,iy) + niz_floor *
     .                                      (nevol*krecz + ngvol*kcxrz)
                          msorxr(ix,iy,ifld) = -(nevol * krecz +
     .                                           ngvol * kcxrz) *
     .                                          (ni(ix,iy,ifld))*
     .                                             mi(ifld)*up(ix,iy,ifld)
                          nuix(ix,iy,jg) = nuix(ix,iy,jg) + nucx(ix,iy,jg) +
     .                                     nuelg(ix,iy,jg)
                        endif


                     enddo   # end do-loop on charge states for isotope jz
                  enddo  # end do-loop on impurity isotopes

                  if (istimingon .eq. 1) call timimpfj (tsimp, xc)
c
c   neutral particle source/sink for isupgon=1; make consistent with impurity
c   contributions just calculated for multispecies
              if (isupgon(1) .eq. 1) then #should be generalized to D & T
                 psor(ix,iy,iigsp)= -psor(ix,iy,1)
                 psorxr(ix,iy,iigsp)= -psorxr(ix,iy,1)
              endif
c
  603        continue
  604     continue
        endif            # end of if (isimpon .ge. 5 .and. nzspt .gt. 0)

c ... Add volume particle loss terms for quasi 1-D radial model
      if (l_parloss .le. 1e9) then
        do iy = iys1, iyf6  # core region has no loss
          do ix = ixs1, ixf6
            do ifld = 1, nfsp
              if (iy .le. iysptrx) then # inside the LCFS when nxpt>1
                                        # (see definition of iysptrx in nphygeo)
                nuvl(ix,iy,ifld) = 0.
              else
                nuvl(ix,iy,ifld) = ( cfvlh*
     .                             sqrt((te(ix,iy)+ti(ix,iy))/mi(1))+
     .                             cfvli(ifld)*
     .              sqrt((zi(ifld)*te(ix,iy)+ti(ix,iy))/mi(ifld)) ) /
     .                                                      l_parloss
              endif
            enddo
          enddo
        enddo
      endif

c ... Set up nuiz & sources for hydrogen molecular gas
      if (ishymol .eq. 1) then
         if (nhgsp .eq. 1) then
          call xerrab('*** nhgsp must exceed 1 for ishymol=1 ***')
        endif
        do iy = iys1, iyf6
         do ix = ixs1, ixf6
           nuiz(ix,iy,2) = ne(ix,iy) * (  #.. Tom 
     .                          svdiss( te(ix,iy) )
     .                        + cfizmol*rsa(te(ix,iy),ne_sgvi,rtau(ix,iy),0)
     .                        + sigvi_floor )
           massfac = 16*mi(1)/(3*(mg(2)+mi(1)))
           nuix(ix,iy,2)= fnuizx*nuiz(ix,iy,2) + 
     .                           massfac*( kelighi(2)*ni(ix,iy,1)+
     .                                     kelighg(2)*ng(ix,iy,1) )
c ...  molecule-molecule collisions would enter viscosity, not nuix
           psorbgg(ix,iy,2) = ngbackg(2)* 
     .                     (0.9+0.1*(ngbackg(2)/ng(ix,iy,2))**ingb ) * 
     .                                        nuiz(ix,iy,2) * vol(ix,iy)
           psorgc(ix,iy,2) = - ng(ix,iy,2)*nuiz(ix,iy,2)*vol(ix,iy) +
     .                        psorbgg(ix,iy,2)
           psorg(ix,iy,2) = psorgc(ix,iy,2)  # no mol sor averaging
           psordis(ix,iy) = -2*psorgc(ix,iy,2)  # 2 atoms per molecule
           if(isupgon(1) .eq. 1) then
             psor(ix,iy,iigsp) = psor(ix,iy,iigsp) + psordis(ix,iy)
           endif
         enddo
        enddo 
      endif  # end of loop for ishymol=1 (hydrogen molecules on)


c  *** Now integrate sources over cell volume if ishosor=1 & yl(neq+1)=-1,
c  *** where the last condition means this is only a full RHS eval, not
c  *** a Jacobian calculation

         if (ishosor.eq.1) then  #full RHS eval

           if (svrpkg.eq."cvode") then    # cannot access yl(neq+1)
            call xerrab('*** svrpkg=cvode not allowed for ishosor=1 **')
           endif 

          if (yl(neq+1).lt.0) then  #full RHS eval

c ...    integ. sources over cells (but not for Jac) for higher-order accuracy

             do ifld = 1, nfsp  # loop over ions
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                 ixm1(0:nx+1,0:ny+1),fsprd, vol(0:nx+1,0:ny+1), 
     .                 psor_tmpov(0:nx+1,0:ny+1), psor(0:nx+1,0:ny+1,ifld))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                 ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                 psor_tmpov(0:nx+1,0:ny+1), psorxr(0:nx+1,0:ny+1,ifld))
             enddo

c *** Now do the gas
             do igsp = 1, ngsp  # now loop over gas species
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                psor_tmpov(0:nx+1,0:ny+1), psorg(0:nx+1,0:ny+1,igsp))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                psor_tmpov(0:nx+1,0:ny+1), psorrg(0:nx+1,0:ny+1,igsp))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                psor_tmpov(0:nx+1,0:ny+1), psorcxg(0:nx+1,0:ny+1,igsp))
             enddo
       
          endif   # end of if (yl(neq+1).lt.0) test
         endif    # end of integrating over sources and ishosor test
      
      endif              # end of big loop starting if (ifixpsor .eq. 0) 

*-----------------------------------------------------------------------
*  -- Calculates the fixed source if it is on
*-----------------------------------------------------------------------

      if (ifixsrc .ne. 0) then
         do 920 iy = j2, j5
            do 910 ix = i2, i5
               snic(ix,iy,1) = snic(ix,iy,1) + vol(ix,iy) * a1n *
     .                          exp(-b1n*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1n*(yyc(iy)-yysrc)**2)
               seic(ix,iy) = seic(ix,iy) + vol(ix,iy) * a1i * ev *
     .                          exp(-b1i*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1i*(yyc(iy)-yysrc)**2)
               seec(ix,iy) = seec(ix,iy) + vol(ix,iy) * a1e * ev *
     .                          exp(-b1e*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1e*(yyc(iy)-yysrc)**2)
 910        continue
 920     continue
      endif



*****************************************************************
c In the case of neutral parallel mom, call neudif to get
c flux fngy, vy and uu, now that we have evaluated nuix etc.
*****************************************************************
      if (TimingPandfOn.gt.0) TimeNeudif=tick()
ccc      if (isupgon .eq. 1 .and. zi(ifld) .eq. 0.0) call neudif
      if (ineudif .eq. 1) then
         call neudif
      elseif (ineudif .eq. 2) then
c ..Timing
      if(istimingon==1) tsnpg=gettime(sec4)
         call neudifpg
c ..Timing
      if(istimingon==1) ttnpg=ttnpg+(gettime(sec4)-tsnpg)
      elseif (ineudif .eq. 3) then
         call neudifl
      else
         call neudifo
      endif
      if (TimingPandfOn.gt.0) TotTimeNeudif=TotTimeNeudif+tock(TimeNeudif)
*****************************************************************
*  Other volume sources calculated in old SRCMOD
*****************************************************************
*  ---------------------------------------------------------------------
*  electron-ion transfer terms and an
*  approximation to the ion-ion thermal force.
*  cfw: For the neutral momentum eq there is also a v_gas*grad(p_gas)
*       term which is evaluated using ng, ti and gpiy
*  ---------------------------------------------------------------------

c...  Force fluxes and gradients on cuts to be zero for half-space problems
      if (isfixlb(1).eq.2.or. isfixrb(1).eq.2) then
         if (isfixlb(1).eq.2) then
            ix = ixpt2(1)
         else
            ix = ixpt1(1)
         endif
         if (ix.ge.i2 .and. ix.le.i5+1 .and. iysptrx1(1) > 0) then
            do iy = 0, iysptrx1(1)
               gpex(ix,iy) = 0.
               frice(ix,iy) = 0.
               ex(ix,iy) = 0.
               upe(ix,iy) = 0.
               do ifld = 1, nfsp
                  gpix(ix,iy,ifld) = 0.
                  frici(ix,iy,ifld) = 0.
                  uu(ix,iy,ifld) = 0.
                  upi(ix,iy,ifld) = 0.
               enddo
            enddo
         endif
      endif

*  -- Set up electron parallel contribution to seec & smoc
      do iy = j2, j5
         do ix = i2, i5
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            nexface = 0.5*(ne(ix2,iy)+ne(ix,iy))
            t1old =.5*cvgp*(upe(ix,iy)*rrv(ix,iy)*
     .          ave(gx(ix,iy),gx(ix2,iy))*gpex(ix,iy)/gxf(ix,iy) +
     .          upe(ix1,iy)*rrv(ix1,iy)*
     .          ave(gx(ix,iy),gx(ix1,iy))*gpex(ix1,iy)/gxf(ix1,iy) )
            t2old = 1.e-20* 0.25*(fqp(ix,iy)+fqp(ix1,iy))*
     .          (ex(ix,iy)+ex(ix1,iy))/gx(ix,iy)
            iyp1 = min(iy+1,ny+1)
            iym1 = max(iy-1,0)
            t1new = .5*cvgp*( vex(ix,iy)*
     .          ave(gx(ix,iy),gx(ix2,iy))*gpex(ix,iy)/gxf(ix,iy) +
     .          vex(ix1,iy)*
     .          ave(gx(ix,iy),gx(ix1,iy))*gpex(ix1,iy)/gxf(ix1,iy) )
            t2new = .5*cvgp*( vey(ix,iy)*
     .          ave(gy(ix,iy),gy(ix,iyp1))*gpey(ix,iy)/gyf(ix,iy) +
     .          vey(ix,iy)*
     .          ave(gy(ix,iy),gy(ix,iym1))*gpey(ix,iym1)/gyf(ix,iym1)
     .                                                            )
            seec(ix,iy) = seec(ix,iy)
     .          + (t1old*vol(ix,iy) - t2old)*oldseec
     .          + ((t1new+t2new)*vol(ix,iy))*(1-oldseec)
            if (nusp-isupgon(1).eq.1) smoc(ix,iy,1)=(( -cpgx*gpex(ix,iy)-
     .                   qe*nexface*gpondpotx(ix,iy) )*rrv(ix,iy)  +
     .                     pondomfpare_use(ix,iy) )*sx(ix,iy)/gxf(ix,iy)
         enddo 
      enddo

*  -- Now loop over all ion species for seec, seic, and smoc --

      do 101 ifld = 1, nusp  #not nfsp; up only for ifld<=nusp
* ------ *
        if(zi(ifld) > 1.e-20) then  #only ions here; atoms follow
* ------ *
*     -- coupling in the x-direction --
*     -- (note sign change in pondomfpari_use term starting 031722)
           do 31 iy = j2, j5
             do 30 ix = i2, i5
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               tv = gpix(ix ,iy,ifld)/gxf(ix,iy)
               t1 = gpix(ix1,iy,ifld)/gxf(ix1,iy)
               t1 = .5*cvgp*( up(ix,iy,ifld)*rrv(ix,iy)*
     .                                 ave(gx(ix2,iy),gx(ix,iy))*tv
     .                      + up(ix1,iy,ifld)*rrv(ix1,iy)*
     .                                 ave(gx(ix,iy),gx(ix1,iy))*t1 )
               seic(ix,iy) = seic(ix,iy) + cfvgpx(ifld)*t1*vol(ix,iy)
               t0 = - cpiup(ifld)*( gpix(ix,iy,ifld)*rrv(ix,iy) -
     .                                  pondomfpari_use(ix,iy,ifld) )*
     .                                            sx(ix,iy)/gxf(ix,iy)
               if (nusp-isupgon(1) .eq. 1) then  # single ion mom. eq.
                  smoc(ix,iy,1) = smoc(ix,iy,1) + cpgx*t0
               else                # multiple mom. eq., so nusp=nisp
                  t0 = t0 +( qe*zi(ifld)*0.5*( ni(ix2,iy,ifld)+
     .                       ni(ix,iy,ifld) )*ex(ix,iy)*rrv(ix,iy) +
     .                       frici(ix,iy,ifld) )* sx(ix,iy)/gxf(ix,iy)
                  if (ifld <= nusp) smoc(ix,iy,ifld) = 
     .                                      smoc(ix,iy,ifld) + cpgx*t0 
               endif
c...  Add friction part of Q_e here
               tv = 0.25*(frice(ix,iy)+frice(ix1,iy))*
     .              ( upe(ix,iy)     + upe(ix1,iy) -
     .              upi(ix,iy,ifld) - upi(ix1,iy,ifld) )
               seec(ix,iy) = seec(ix,iy) - zi(ifld)**2*ni(ix,iy,ifld)*
     .                                        tv*vol(ix,iy)/nz2(ix,iy)
               
   30        continue
   31      continue

*     -- coupling in the x & y-directions --
           do 34 iy = j2, j5
            do 33 ix = i2, i5
             if (isgpye == 0) then
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               vyiy0 = fracvgpgp*vygp(ix,iy,ifld) 
     .                    +(1.-fracvgpgp)*vycb(ix,iy,ifld)
               vyiym1 = fracvgpgp*vygp(ix,iy-1,ifld) 
     .                    +(1.-fracvgpgp)*vycb(ix,iy-1,ifld)
               v2ix0 = fracvgpgp*v2xgp(ix,iy,ifld) 
     .                    +(1.-fracvgpgp)*v2cb(ix,iy,ifld)
               v2ixm1 = fracvgpgp*v2xgp(ix1,iy,ifld) 
     .                    +(1.-fracvgpgp)*v2cb(ix1,iy,ifld)
               t1 =.5*cvgp*( vygp(ix,iy  ,ifld)*gpiy(ix,iy  ,ifld) + 
     .                       vygp(ix,iy-1,ifld)*gpiy(ix,iy-1,ifld) +
     .                    v2xgp(ix ,iy,ifld)*ave(gx(ix,iy),gx(ix2,iy))*
     .                               gpix(ix ,iy,ifld)/gxf(ix ,iy) +    
     .                    v2xgp(ix1,iy,ifld)*ave(gx(ix,iy),gx(ix1,iy))*
     .                               gpix(ix1,iy,ifld)/gxf(ix1,iy) )     
               t2 = t1
             elseif (isgpye == 1) then    # Old B2 model with Jperp=0
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
               t2 = t1
             elseif (isgpye == 2) then    # Knoll expression
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
             endif
             seec(ix,iy) = seec(ix,iy) - fluxfacy*t1 * vol(ix,iy)
             seic(ix,iy) = seic(ix,iy) + fluxfacy*cfvgpy(ifld)*t2*
     .                                                     vol(ix,iy)
   33       continue
   34      continue

        endif  #test on zi(ifld) > 0, so only ion terms
  101 continue  #ifld loop over ion species

c ... Now include seic contribution from hydrogen atoms if isupgon=1
c ... Then "ion" species 2 (redundant as gas species 1) is hydr atom

      if(isupgon(1)==1 .and. zi(2)<1.e-20) then # .and. istgon(1)==0) then 
        if(cfvgpx(2) > 0.) then
          do iy = j2, j5
            do ix = i2, i5
              ix1 = ixm1(ix,iy)
              ix2 = ixp1(ix,iy)
              iy1 = max(0,iy-1)
              seic(ix,iy) = seic(ix,iy) + cftiexclg
     .                                   *0.5*cfvgpx(2)*( 
     .                       uuxg(ix, iy,1)*gpix(ix,iy,2) +
     .                       uuxg(ix1,iy,1)*gpix(ix1,iy,2) )*vol(ix,iy)
              seic(ix,iy) = seic(ix,iy) + cftiexclg
     .                                   *0.5*cfvgpy(2)*( 
     .                        vyg(ix, iy,1)*gpiy(ix,iy,2) +
     .                        vyg(ix,iy1,1)*gpiy(ix,iy1,2) )*vol(ix,iy)
            enddo
          enddo
        else  # Here if cfvgpx(2)=0, old vpar_g*grad_Pg only => ifld=2
          do iy = j2, j5
            do ix = i2, i5
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               tv = gpix(ix ,iy,2)/gxf(ix,iy)
               t1 = gpix(ix1,iy,2)/gxf(ix1,iy)
               t1 = .5*cvgp*( up(ix,iy,2)*rrv(ix,iy)*
     .                               ave(gx(ix2,iy),gx(ix,iy))*tv
     .                    + up(ix1,iy,2)*rrv(ix1,iy)*
     .                               ave(gx(ix,iy),gx(ix1,iy))*t1 )
               seic(ix,iy) = seic(ix,iy) + cftiexclg*t1*vol(ix,iy)
            enddo
          enddo
        endif  #test on cfvgpx(2) > 0 or = 0
      endif   #test for inertial neutrals

*****************************************************************
*  Other physics coefficients. (old PHYVIS)
*****************************************************************

* -- loop over species number --

      do 102 ifld = 1, nfsp

c
c     neutral viscosity for isupgon=1
c
c  Logic here is specialized; only ifld=2 (igsp=1) has viscosity calc.
c  If more neutral species have full parallel mom eqn, need to redo loops
         if(isupgon(1) .eq. 1 .and. zi(ifld) .eq. 0)then
            do 936 iy = j1,j6
               iyp1 = min(iy+1,ny+1)
               do 937 ix = i1,i6
c
                  ix1 = ixm1(ix,iy)
                  vtn = sqrt(max(tg(ix,iy,1),tgmin*ev)/mi(ifld))
 		  qfl = flalfvgxa(ix)*nm(ix,iy,ifld)*vtn**2
                  if(isvisxn_old == 1) then
                    lmfpn = 1./(sigcx * 
     .                          (ni(ix,iy,1) + rnn2cx*ni(ix,iy,ifld)))
                  elseif(isvisxn_old==0 .and. ishymol==0) then
                    lmfppar = vtn/(kelhihg*ni(ix,iy,1) +
     .                                         kelhghg*ni(ix,iy,ifld))
                    lmfpperp = vtn/( vtn*sigcx*ni(ix,iy,1) + 
     .                      kelhihg*ni(ix,iy,1)+kelhghg*ni(ix,iy,ifld) )
                    rrfac = rr(ix,iy)*rr(ix,iy)
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  else   # (isvisxn_old=0 .and. ishymol=1) then #with mols
                    lmfppar = vtn/(kelhihg*ni(ix,iy,1) +
     .                     kelhghg*ni(ix,iy,ifld) + kelhmhg*ng(ix,iy,2))
                    lmfpperp = vtn/( vtn*sigcx*ni(ix,iy,1) + 
     .                     kelhihg*ni(ix,iy,1) +kelhghg*ni(ix,iy,ifld) +
     .                     kelhmhg*ng(ix,iy,2) )
                    rrfac = rr(ix,iy)*rr(ix,iy)
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  endif
                  csh = lmfpn*nm(ix,iy,ifld)*vtn*
     .                                      lgvmax/(lgvmax + lmfpn)  
                  if (isgxvon .eq. 0) then 
                    qsh = csh * (up(ix1,iy,ifld)-up(ix,iy,ifld))
     .                                         *gx(ix,iy)
                  elseif (isgxvon .eq. 1) then
                    qsh = csh * (up(ix1,iy,ifld)-up(ix,iy,ifld))
     .                *2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
                  endif
                  visx(ix,iy,ifld)= cfvisxn*csh/ 
     .               (1 + (abs(qsh/(qfl+cutlo))**flgamvg))**(1./flgamvg)
     .               + cfanomvisxg*travis(ifld)*nm(ix,iy,ifld)

c    Now do y-direction; use ni on up y-face
                  ix2 = ixp1(ix,iy)
                  ix3 = ixp1(ix,iyp1)
                  tgupyface = 0.25*( tg(ix,iy,1)+
     .                         tg(ix,iyp1,1)+tg(ix2,iy,1)+
     .                         tg(ix3,iyp1,1) )
                  vtn = sqrt(max(tgupyface,tgmin*ev)/mi(ifld))
                  nmxface = 0.5*(nm(ix,iy,ifld)+nm(ix2,iy,ifld))
                  ngupyface = 0.25*( ni(ix,iy,ifld)+
     .                         ni(ix,iyp1,ifld)+ni(ix2,iy,ifld)+
     .                         ni(ix3,iyp1,ifld) )
                  n1upyface = 0.25*( ni(ix,iy,1)+
     .                         ni(ix,iyp1,1)+ni(ix2,iy,1)+
     .                         ni(ix3,iyp1,1) )
                  if(ishymol == 0) then
                    lmfppar = vtn/(kelhihg*n1upyface +
     .                                         kelhghg*ngupyface)
                    lmfpperp = vtn/( vtn*sigcx*n1upyface + 
     .                      kelhihg*n1upyface+kelhghg*ngupyface )
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  else  # ishymol=1
                    ng2upyface = 0.25*( ng(ix,iy,2)+
     .                         ng(ix,iyp1,2)+ng(ix2,iy,2)+
     .                         ng(ix3,iyp1,2) )
                    lmfppar = vtn/(kelhihg*n1upyface +
     .                     kelhghg*n1upyface + kelhmhg*ng2upyface)
                    lmfpperp = vtn/( vtn*sigcx*n1upyface + 
     .                     kelhihg*n1upyface +kelhghg*ngupyface +
     .                     kelhmhg*ng2upyface )
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  endif
        
                  csh = lmfpn*ngupyface*mi(ifld)*vtn*
     .                                      lgvmax/(lgvmax + lmfpn)  

 		  qfl = flalfvgya(iy)*ngupyface*mi(ifld)*vtn**2
                  qsh = csh * (up(ix,iy,ifld)-up(ix,iyp1,ifld)) * 
     .                                                        gyf(ix,iy)
                  visy(ix,iy,ifld)= cfvisyn*csh / 
     .               (1 + (abs(qsh/(qfl+cutlo))**flgamvg))**(1./flgamvg)
     .               + cfanomvisyg*travis(ifld)*nmxface
c
 937           continue
 936        continue
         endif
c
c
       if(zi(ifld) > 1.e-20) then
         do 39 iy = j1, j6
            do 38 ix = i1, i6
               w(ix,iy) = 0.0e0
 38         continue
 39      continue

         do 42 jfld = 1, nisp
            tv = zi(jfld)**2 / sqrt((mi(ifld)+mi(jfld))/(2*mp))
            do 41 iy = j1, j6
               do 40 ix = i1, i6
                  w(ix,iy) = w(ix,iy) + tv*ni(ix,iy,jfld)
   40          continue
   41       continue
   42    continue

         do 44 iy = j1, j6
            do 43 ix = i1, i6
	       ctaui(ix,iy,ifld) = 2.1e13/(loglambda(ix,iy)*zi(ifld)**2) # mass fac?
               tv2 = ctaui(ix,iy,ifld)/(ev*sqrt(ev))
               if (convis .eq. 0) then
                  a = max (ti(ix,iy), temin*ev)
               else
                  a = afix*ev
               endif
	       epstmp = max(epsneo(ix,iy), 1.e-50)
               visxtmp = tv2 * coef * rr(ix,iy) * rr(ix,iy) *
     .                            a*a*sqrt(a) * ni(ix,iy,ifld)/w(ix,iy)
               visx(ix,iy,ifld) = parvis(ifld)*visxtmp +
     .                            trax_use(ix,iy,ifld)*nm(ix,iy,ifld)
               nuii(ix,iy,ifld) = w(ix,iy)/(tv2*a*sqrt(a))
               nuiistar(ix,iy,ifld) = ( lconneo(ix,iy)*nuii(ix,iy,ifld)/
     .                                 (epstmp**1.5*(2*a/mi(ifld))**0.5)
     .                                           + 1.e-50 )
               visxneo(ix,iy,ifld) = visxtmp*
     .                 (1./(1.+epstmp**(-1.5)/nuiistar(ix,iy,ifld)))*
     .                               (1./(1.+1./nuiistar(ix,iy,ifld)))
               rt2nus = 1.414*nuiistar(ix,iy,ifld)
               ktneo(ix,iy,ifld) = (-0.17 + 1.05*rt2nus**.5 + 
     .                     2.7*rt2nus**2*epstmp**3) / ( 1.+
     .                0.7*rt2nus**.5 + rt2nus**2*epstmp**3 )
               alfneo(ix,iy,ifld) = (8./15.)*(ktneo(ix,iy,ifld) - 1.)*
     .                 (1./(1.+epstmp**(-1.5)/nuiistar(ix,iy,ifld)))*
     .                            ( 1./(1.+1./nuiistar(ix,iy,ifld)) )
               k2neo(ix,iy,ifld) =(.66 + 1.88*epstmp**.5 - 1.54*epstmp)/ 
     .                            (1. + 1.03*rt2nus**.5 + 0.31*rt2nus) +
     .             1.17*epstmp**3*rt2nus/(1. + 0.74*epstmp**1.5*rt2nus)
c...  flux limit the viscosity; beware of using visx(0,iy) and 
c...  visx(nx+1,iy) as they are meaningless when flux limited
               ix1 = ixm1(ix,iy)
               t0 = max (ti(ix,iy), temin*ev)
               vtn = sqrt(t0/mi(ifld))
               mfl = flalfv * nm(ix,iy,ifld) * rr(ix,iy) *
     .               vol(ix,iy) * gx(ix,iy) * (t0/mi(ifld)) 
ccc  Distance between veloc. cell centers:
               if (isgxvon .eq. 0) then     # dx(ix)=1/gx(ix)
                 csh = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .                            * gx(ix,iy)
               elseif (isgxvon .eq. 1) then # dx(ix)=.5/gxf(ix-1)+.5/gxf(ix)
                 csh = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .               * 2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
               endif
ccc 
               msh = abs( csh*(upi(ix1,iy,ifld) - upi(ix,iy,ifld)) )
               visx(ix,iy,ifld) = visx(ix,iy,ifld)
     .               / (1 + (msh/(mfl+1.e-20*msh))**flgamv )**(1/flgamv)
               visy(ix,iy,ifld)=(fcdif*travis(ifld)+ tray_use(ix,iy,ifld))*
     .                              nm(ix,iy,ifld) +  4*eta1(ix,iy)
   43       continue
   44    continue
       endif      # test if zi(ifld) > 1.e-20
  102 continue    # large loop for ifld = 1, nfsp

*****************************************************************
*****************************************************************
*  Heat Conduction. (old PHYTHC)
*****************************************************************
*  ---------------------------------------------------------------------
*  compute conductivities on cell faces
*  ---------------------------------------------------------------------

*  -- initialize to 0 --

      do 706 iy = j1, j6
         do 705 ix = i1, i6
            hcxe(ix,iy) = 0.0e0
            hcxi(ix,iy) = 0.0e0
            hcxineo(ix,iy) = 0.0e0
            hcye(ix,iy) = 0.0e0
            hcyi(ix,iy) = 0.0e0
            do ifld = 1, nisp
               hcxij(ix,iy,ifld) = 0.0e0
               hcyij(ix,iy,ifld) = 0.0e0
            enddo
  705    continue
  706 continue

*  -- loop over species number --

      do 103 ifld = 1, nisp
c -- Skip this if these are the neutrals (zi(ifld).eq.0)
         if (zi(ifld) .eq. 0.0e0) goto 103

c...  Initialize w1 and w2 for each species
         do 49 iy = j1, j6
            do 48 ix = i1, i6
               w1(ix,iy) = 0.0e0
               w2(ix,iy) = 0.0e0
 48         continue
 49      continue

*     -- conductivities --
*        The poloidal conductivities  are initially computed without
*        the factor rr**2 * tv**2.5)

         do 52 jfld = 1, nisp
            tv = zi(jfld)**2
            a = zi(jfld)**2 *
     .            sqrt(2*mi(ifld)*mi(jfld)/(mi(ifld)+mi(jfld)))
            do 51 iy = j1, j6
               do 50 ix = i1, i6
                  ix1 = ixp1(ix,iy)
                  w1(ix,iy) = w1(ix,iy) + tv*(ni(ix,iy,jfld)*gx(ix,iy) +
     .                                     ni(ix1,iy,jfld)*gx(ix1,iy)) / 
     .                                        (gx(ix,iy) + gx(ix1,iy))
                  w2(ix,iy) = w2(ix,iy) + a*(ni(ix,iy,jfld)*gx(ix,iy) +
     .                                     ni(ix1,iy,jfld)*gx(ix1,iy)) / 
     .                                        (gx(ix,iy) + gx(ix1,iy))
   50          continue       
   51       continue
   52    continue

         do 59 iy = j1, j6
            do 58 ix = i1, i6
               ix1 = ixp1(ix,iy)
               iyp1 = min(ny+1, iy+1)
               ctaue(ix,iy,ifld) = 3.5e11*zi(ifld)/loglambda(ix,iy)
               ctaui(ix,iy,ifld) =2.1e13/(loglambda(ix,iy)*zi(ifld)**2)
               fxe = kxe * ce * ctaue(ix,iy,ifld) / (me*ev*sqrt(ev))
               fxi = kxi * ci * ctaui(ix,iy,ifld) / (ev*sqrt(ev*mp))
               fxet = fxe
               fxit = fxi
               do jx = 1, nxpt  #reduce kxe inside sep by rkxecore fac
                  if ( (iy.le.iysptrx) .and. 
     .                    ix.gt.ixpt1(jx) .and. ix.le.ixpt2(jx) ) then
                     fxet = fxe/( 1. + (rkxecore-1.)*
     .                              (yyf(iy)/(yyf(0)+4.e-50))**inkxc )
                     fxit = kxicore * fxi
                  endif
               enddo
               niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                    ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                     (gx(ix,iy) + gx(ix1,iy))
               niavey = ( niy0(ix,iy,ifld)*gy(ix,iy) +
     .                    niy1(ix,iy,ifld)*gy(ix,iyp1)) /
     .                     (gy(ix,iy) + gy(ix,iyp1))
               hcxe(ix,iy) = hcxe(ix,iy)+fxet*niavex/w1(ix,iy)

c ... Use fixed diffusivity inside the separatrix, anomalous outside,
c     if anomalous-diffusivity multiplier is nonzero.
               kyemix = fcdif*kye + kye_use(ix,iy)
cccMER NOTE: when there are multiple x-points, as in 'dnull' configuration,
ccc          iysptrx is the last closed flux surface (see S.R. nphygeo)
               if(kyet .gt. 1.e-20 .and. iy .gt. iysptrx) then
                  kyemix = (1. - ckyet) * kyemix +
     .               ckyet * kyet * diffusivwrk(ix,iy)
               endif
               hcye(ix,iy) = hcye(ix,iy) + ( kyemix +
     .                       2.33*(dclass_e(ix,iy)+dclass_e(ix,iyp1)) )*
     .                                               zi(ifld) * niavey

               hcxij(ix,iy,ifld) = fxit*niavex/w2(ix,iy)
               kyimix = fcdif*kyi + kyi_use(ix,iy)
cccMER NOTE: when there are multiple x-points, as in 'dnull' configuration,
ccc          iysptrx is the last closed flux surface (see S.R. nphygeo)
               if(kyit .gt. 1.e-20 .and. iy .gt. iysptrx) then
                  kyimix = (1. - ckyit) * kyimix +
     .               ckyit * kyit * diffusivwrk(ix,iy)
               endif
               hcyij(ix,iy,ifld) = hcyij(ix,iy,ifld) + ( kyimix +
     .                           (dclass_i(ix,iy)+dclass_i(ix,iyp1)) )*
     .                                                         niavey
   58       continue
   59    continue

  103 continue
  
c ... Add ion temp. dep. for pol. terms, flux limit, & build total ion hcx,yi
      do 595 ifld = 1, nisp
         if (zi(ifld) .eq. 0.e0) goto 595
         do iy = j1, j6
            do ix = i1, i6
               ix1 = ixp1(ix,iy)
               if (concap .eq. 0) then
                  tiave = (ti(ix,iy)*gx(ix,iy) + ti(ix1,iy)*gx(ix1,iy)) /
     .                                         (gx(ix,iy) + gx(ix1,iy))
                  do jx = 1, nxpt
                    if (ix==ixlb(jx).and.ixmnbcl==1) tiave=ti(ixlb(jx)+1,iy)
                    if (ix==ixrb(jx).and.ixmxbcl==1) tiave=ti(ixrb(jx),iy)
                  enddo
                  a = max (tiave, temin*ev)
               else
                  a = afix*ev
               endif
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld)*rrv(ix,iy)*
     .                                            rrv(ix,iy)*a*a*sqrt(a)
c...  reduce hcxij if ti very flat; prevents large conduction for high ti
c...  or if lmfpi exceeds a mean-free path limit, lmfplim
               lmfpi = 1.e16*(tiave/ev)**2/ni(ix,iy,1) # i-mean-free-path
               niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                    ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                                     (gx(ix,iy) + gx(ix1,iy))
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld)/(1.+lmfpi/lmfplim)
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld) *
     .                         ( cutlo + (ti(ix,iy)-ti(ix1,iy))**2 ) /
     .                         ( cutlo + (ti(ix,iy)-ti(ix1,iy))**2 +
     .                      (0.5*alfkxi*(ti(ix,iy)+ti(ix1,iy)))**2 ) +
     .                         kxi_use(ix,iy)*niavex

c ... Flux limit individ. hcxij in poloidal direction if isflxldi=2
               if (isflxldi .eq. 2) then
                  niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                       ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                                        (gx(ix,iy) + gx(ix1,iy))
                  wallfac = 1.
                  do jx = 1, nxpt
                     if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                      (ix==ixrb(jx).and.ixmxbcl==1) )
     .                    .and. (isplflxl==0) ) wallfac = flalfipl/flalfi
                  enddo
                  qflx = wallfac*flalfi * rrv(ix,iy) *
     .                                    sqrt(a/mi(ifld)) * niavex * a
                  cshx = hcxij(ix,iy,ifld)
	          lxtic = 0.5*(ti(ix,iy)+ti(ix1,iy)) /
     .               (abs(ti(ix,iy)-ti(ix1,iy))*gxf(ix,iy) + 100.*cutlo)
	          qshx = cshx * (ti(ix,iy)-ti(ix1,iy)) * gxf(ix,iy) *
     .                                              (1. + lxtic/lxtimax)
                  hcxij(ix,iy,ifld) = cshx  / (1 + abs(qshx/qflx))
               endif
               hcxi(ix,iy) = hcxi(ix,iy) + hcxij(ix,iy,ifld)
               hcxineo(ix,iy) = hcxineo(ix,iy) + hcxij(ix,iy,ifld)*
     .                      1.5676*epsneo(ix,iy)**1.5/k2neo(ix,iy,ifld)
               hcyi(ix,iy) = hcyi(ix,iy) + hcyij(ix,iy,ifld)
               qipar(ix,iy,ifld) = hcxij(ix,iy,ifld)*gxf(ix,iy)*
     .                                (ti(ix,iy)-ti(ix1,iy))/rrv(ix,iy)
            enddo
         enddo
 595  continue 

c...  Now include elec. temp and other dep. in poloidal terms + diff. neut.
      do 61 iy = j1, j6
         do 60 ix = i1, i6
            ix1 = ixp1(ix,iy)
            iyp1 = min(ny+1, iy+1)
            if (concap .eq. 0) then
               teave = (te(ix,iy)*gx(ix,iy) + te(ix1,iy)*gx(ix1,iy)) /
     .                                       (gx(ix,iy) + gx(ix1,iy))
               do jx = 1, nxpt
                 if(ix==ixlb(jx).and.ixmnbcl==1) teave=te(ixlb(jx)+1,iy)
                 if(ix==ixrb(jx).and.ixmxbcl==1) teave=te(ixrb(jx),iy)
               enddo
               a = max (teave, temin*ev)
            else
               a = afix*ev
            endif
            zeffave = (zeff(ix,iy)*gx(ix,iy) + zeff(ix1,iy)*gx(ix1,iy)) /
     .                                         (gx(ix,iy) + gx(ix1,iy))
            zcoef = 0.308 + 0.767*zeffave - 0.075*zeffave**2
            hcxe(ix,iy) = hcxe(ix,iy)*rrv(ix,iy)*rrv(ix,iy)*a*a*sqrt(a)
     .                    *zcoef

c...  reduce hcxe if te very flat; prevents very large conduction for high te
c...  or if the mean-free path exceeds lmfplim
            lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)  #mfp for elec [m]
            neavex = ( ne(ix ,iy)*gx(ix ,iy) +
     .                 ne(ix1,iy)*gx(ix1,iy))/(gx(ix,iy) + gx(ix1,iy))
            hcxe(ix,iy) = hcxe(ix,iy) *
     .                         ( cutlo + (te(ix,iy)-te(ix1,iy))**2 ) /
     .                         ( cutlo + (te(ix,iy)-te(ix1,iy))**2 +
     .                     (0.5*alfkxe*(te(ix,iy)+te(ix1,iy)))**2 ) +
     .                        kxe_use(ix,iy)*neavex
            hcxe(ix,iy) = hcxe(ix,iy) /( (1. + lmfpe/lmfplim) *
     .                  (1+hcxe(ix,iy)*gx(ix,iy)**2*tdiflim/ne(ix,iy)) )
            if (isupgon(1).eq.0) then   # add diff. gas cx contrib. to hcxi
               hcxn(ix,iy) = 0.
               hcyn(ix,iy) = 0.
c..1dn0802
c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add
               hcxi(ix,iy) = hcxi(ix,iy)
     .           + cftiexclg*cfneut*cfneutsor_ei*kxn*( ng(ix ,iy,1)*ti(ix ,iy)
     .                              +ng(ix1,iy,1)*ti(ix1,iy) ) /
     .                  (mi(1)*(nucx(ix,iy,1) + nucx(ix1,iy,1)))
               hcyi(ix,iy) = hcyi(ix,iy)
     .           + cftiexclg*cfneut*cfneutsor_ei*kyn*( ngy0(ix,iy,1)*tiy0(ix,iy)
     .                              +ngy1(ix,iy,1)*tiy1(ix,iy) ) /
     .                  (mi(1)*(nucx(ix,iy,1) + nucx(ix,iyp1,1)))
            endif

   60    continue
 61   continue
c
c
      if (isupgon(1).eq.1) then
c
c ----- Section for the inertial neutral fluid; we need to do different
c ----- things than for the ions. Note third index=iigsp is neutral species
c ----- The inertial neutrals coeff. are flux-limited and add to total here
         do 62 iy = j1, j6
            iy1 = min(iy,ny)   #dont use j5 because hcx also in loop (not imp.)
            do 63 ix = i1, i6
               ix1 = ixp1(ix,iy)
               tgavex = max(0.5*(tg(ix,iy,1) + tg(ix1,iy,1)), temin*ev)
               tgavey= max(0.5*(tgy0(ix,iy,1)+tgy1(ix,iy,1)), temin*ev)
               niavex = 0.5*(ni(ix,iy,1) + ni(ix1,iy,1)) #only for coll. term
               niavey = 0.5*(niy0(ix,iy1,1) + niy1(ix,iy1,1)) #only coll. term
               noavex = ( ni(ix ,iy,iigsp)*gx(ix ,iy) +
     .                    ni(ix1,iy,iigsp)*gx(ix1,iy)) /
     .                     (gx(ix,iy) + gx(ix1,iy))
               noavey = 0.5*(niy0(ix,iy1,iigsp) + niy1(ix,iy1,iigsp))

c          Set up flux-limit variables (no rrv here) 
c          First limit the poloidal coeff, then radial
c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add
               qflx = flalftgxa(ix) * sqrt(tgavex/mi(iigsp)) * noavex *
     .                                                     tgavex
               lmfpn = 1./(sigcx * (niavex + rnn2cx*noavex))
               cshx = lmfpn*sqrt(tgavex/mi(iigsp))*noavex * 
     .                         lgtmax(iigsp)/(lmfpn + lgtmax(iigsp))
               qshx = cshx * (tg(ix,iy,1)-tg(ix1,iy,1)) * gxf(ix,iy)
	       hcxn(ix,iy) = cshx  / 
     .                      (1 + (abs(qshx/qflx))**flgamtg)**(1./flgamtg)
               hcxi(ix,iy) = hcxi(ix,iy) + 
     .                          cftiexclg*cfneut*cfneutsor_ei*hcxn(ix,iy)
c          Now for the radial flux limit - good for nonorthog grid too
               qfly = flalftgya(iy) * sqrt(tgavey/mi(iigsp)) * noavey *
     .                                                     tgavey
               lmfpn = 1./(sigcx * (niavey + rnn2cx*noavey))
               cshy = lmfpn*sqrt(tgavey/mi(iigsp))*noavey * 
     .                         lgtmax(iigsp)/(lmfpn + lgtmax(iigsp))
               qshy = cshy * (tgy0(ix,iy1,1)-tgy1(ix,iy1,1))/dynog(ix,iy)
               hcyn(ix,iy) = cshy / 
     .                      (1 + (abs(qshy/qfly))**flgamtg)**(1./flgamtg)
               hcyi(ix,iy) = hcyi(ix,iy) + 
     .                          cftiexclg*cfneut*cfneutsor_ei*hcyn(ix,iy)
c     
  63        continue
  62     continue
      endif

*  Equipartition (old PHYEQP)
*****************************************************************
*  ---------------------------------------------------------------------
*  compute equipartition.
*  ---------------------------------------------------------------------

*     -- initialize w3 --
      do 69 iy = j1, j6
         do 68 ix = i1, i6
            w3(ix,iy) = 0.0e0
 68      continue
 69   continue

      do 72 ifld = 1, nisp
         tv = zi(ifld)**2/mi(ifld)
         do 71 iy = j2, j5
            do 70 ix = i2, i5
               w3(ix,iy) = w3(ix,iy) + tv*ni(ix,iy,ifld)
   70       continue
   71    continue
   72 continue

*  -- compute equipartition --
ccc In detail, coef1 = qe**4*sqrt(me)*lnlam / ((2*pi)**1.5*eps0**2)
      do 74 iy = j2, j5
         do 73 ix = i2, i5
            ix2 = ixm1(ix,iy)
            a = max (te(ix,iy), temin*ev)
            loglmcc = 0.5*(loglambda(ix,iy)+loglambda(ix2,iy))
            coef1 = feqp*4.8e-15*loglmcc*sqrt(ev)*ev*mp
            eqp(ix,iy) = coef1 * w3(ix,iy) * ne(ix,iy) / (a*sqrt(a))
c...       reduce eqp when (te-ti)/(te+ti) << 1
            eqp(ix,iy) = eqp(ix,iy) * (a-ti(ix,iy))**2 / ( cutlo +
     .                (a-ti(ix,iy))**2 + (alfeqp*(a+ti(ix,iy)))**2 )
   73    continue
   74 continue

*********************************************************
c ... Gas thermal coefficients, initially for molecules *
*********************************************************
*
c ... Gas thermal conductivity coeffs - from self-collisions
*****************************************************************

      if (nisp >= 2) then  # uses ni(,,2), so must have atoms
       do igsp = 1, ngsp
        do iy = j1, j6
        iy1 = min(iy,ny)
          do ix = i1, i6
            ix1 = ixp1(ix,iy)
            tgavex = max( (tg(ix,iy,igsp)*gx(ix,iy) +
     .                              tg(ix1,iy,igsp)*gx(ix1,iy)) /
     .                             (gx(ix,iy) + gx(ix1,iy)), temin*ev )
            tgavey=max(0.5*(tgy0(ix,iy,igsp)+tgy1(ix,iy,igsp)),temin*ev)
            noavex = ( ng(ix,iy,igsp)*gx(ix,iy) +
     .                                   ng(ix1,iy,igsp)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            niavex = ( ni(ix,iy,1)*gx(ix,iy) +
     .                                   ni(ix1,iy,1)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            naavex = ( ni(ix,iy,2)*gx(ix,iy) +
     .                                   ni(ix1,iy,2)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            noavey = 0.5*(ngy0(ix,iy1,igsp) + ngy1(ix,iy1,igsp))
            niavey = 0.5*(niy0(ix,iy1,1) + niy1(ix,iy1,1))
            naavey = 0.5*(niy0(ix,iy1,2) + niy1(ix,iy1,2))
            nuelmolx = noavex*kelhmhm + niavex*kelhmhg + 
     .                 naavex*kelhmhg
            qflx = flalftmx*sqrt(tgavex/mg(igsp))*noavex*tgavex
            cshx = cftgcond*noavex*tgavex/(mg(igsp)*nuelmolx)  #assume K not fcn Tg
            qshx = cshx * (tg(ix,iy,igsp)-tg(ix1,iy,igsp)) * gxf(ix,iy)
            hcxg(ix,iy,igsp) = cshx / 
     .                     (1.+ (abs(qshx/qflx))**flgamtg)**(1./flgamtg)
            hcxg(ix,iy,igsp)=(1.-cfhcxgc(igsp))*hcxg(ix,iy,igsp)+
     .                     cfhcxgc(igsp)*noavex*kxg_use(ix,iy,igsp)
c..   Now radial direction
            nuelmoly = noavey*kelhmhm + niavey*kelhmhg + 
     .                 naavey*kelhmhg
            qfly = flalftmy*sqrt(tgavey/mg(igsp))*noavey*tgavey
            cshy = cftgcond*noavey*tgavey/(mg(igsp)*nuelmoly)  #assume Kel_s not fcn Tg
            qshy = cshy*(tgy0(ix,iy1,igsp)-tgy1(ix,iy1,igsp))/
     .                                                  dynog(ix,iy)
            hcyg(ix,iy,igsp) = cshy / 
     .                     (1 + (abs(qshy/qfly))**flgamtg)**(1./flgamtg)
            hcyg(ix,iy,igsp)=(1-cfhcygc(igsp))*hcyg(ix,iy,igsp)+
     .                     cfhcygc(igsp)*noavey*kyg_use(ix,iy,igsp)
          enddo
        enddo
        if (igsp.eq.1 .and. isupgon(igsp).eq.1) then 
          hcxg(:,:,igsp) = hcxn(:,:)
          hcyg(:,:,igsp) = hcyn(:,:)
        endif
       enddo
      endif

c ... Gas molecule thermal equipartition with hydrogen ions and atoms
*****************************************************************
      if (nisp >= 2) then   # uses ni(,,2), so must have atoms
       do igsp = 1, ngsp
        do iy = j1, j6
          do ix = i1, i6
	    nhi_nha = ni(ix,iy,1)+ni(ix,iy,2)
#            eqpg(ix,iy,igsp) = cftgeqp*ng(ix,iy,igsp)*nhi_nha*
#     .                                            keligig(igsp)
            eqpg(ix,iy,igsp) = cftgeqp*ng(ix,iy,igsp)*
     .                   (ni(ix,iy,1)+cftiexclg*ni(ix,iy,2))*keligig(igsp)
          enddo
        enddo
       enddo
      endif

c ... Call routine to evaluate gas energy fluxes
****************************************************************
      call engbalg


*****************************************************************
*****************************************************************
*  Here starts the old routine PARBAL
*****************************************************************
      do 104 ifld = 1, nfsp
*  ---------------------------------------------------------------------
*     compute flux, residual
*     The residual is: res := snic + sniv * ni - outflow(ni).
*  ---------------------------------------------------------------------

*  -- compute fnix --

         methnx = mod(methn, 10)
         methny = methn/10
         do 81 iy = j4, j8
            do 80 ix = i1, i5
              if ( zi(ifld).eq.0. .and. ineudif.ne.0 .and.
     .                                   1.-rrv(ix,iy) > 1.e-4 ) then
                 fnix(ix,iy,ifld) = fngx(ix,iy,1)
              else
               ix2 = ixp1(ix,iy)

               if (methnx .eq. 2) then   # central differencing
                  t2 = ( ni(ix, iy,ifld) + ni(ix2,iy,ifld) ) / 2

               elseif (methnx .eq. 3) then   # upwind differencing

                  if( uu(ix,iy,ifld) .ge. 0.) then
                     t2 = ni(ix,iy,ifld)
                  else
                     t2 = ni(ix2,iy,ifld)
                  endif

               elseif (methnx .eq. 6) then   # log central differencing
                  t2 = exp(0.5*
     .                ( log(ni(ix,iy,ifld)) + log(ni(ix2,iy,ifld)) ))

               else   # interp. ave or harmonic ave depending on wind*grad

                  t0 = ( ni(ix, iy,ifld)*gx(ix, iy) +
     .                   ni(ix2,iy,ifld)*gx(ix2,iy) ) / 
     .                                      ( gx(ix,iy)+gx(ix2,iy) )
                  t1 = ( gx(ix,iy)+gx(ix2,iy) ) * ni(ix,iy,ifld) *
     .                   ni(ix2,iy,ifld) / ( cutlo + ni(ix,iy,ifld)*
     .                   gx(ix2,iy) + ni(ix2,iy,ifld)*gx(ix,iy) )
                  if( uu(ix,iy,ifld)*(ni(ix,iy,ifld)-ni(ix2,iy,ifld)) 
     .                                                     .ge. 0.) then
                     t2 = t0
                  else
                     t2 = t1
                  endif

               endif

               fnix(ix,iy,ifld) = cnfx*uu(ix,iy,ifld) * sx(ix,iy) * t2
               fnixcb(ix,iy,ifld)=cnfx*sx(ix,iy) * t2 * 0.5*
     .                 (rbfbt(ix,iy) + rbfbt(ix2,iy))*v2cb(ix,iy,ifld)
                  fnix(ix,iy,ifld) = fnix(ix,iy,ifld)/sqrt( 1 +
     .              (nlimix(ifld)*ni(ix ,iy,ifld)/ni(ix2,iy,ifld))**2 +
     .              (nlimix(ifld)*ni(ix2,iy,ifld)/ni(ix ,iy,ifld))**2 )
              endif
   80      continue
           if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
              fnix(nxc-1,iy,ifld)=0.
              fnix(nxc,  iy,ifld)=0.
              fnix(nxc+1,iy,ifld)=0.
              uu(nxc-1,iy,ifld) = 0.
              uu(nxc  ,iy,ifld) = 0.
              uu(nxc+1,iy,ifld) = 0.
              vytan(nxc-1,iy,ifld) = 0.
              vytan(nxc  ,iy,ifld) = 0.
              vytan(nxc+1,iy,ifld) = 0.
 
           endif
           if (islimon.ne.0 .and. iy.ge.iy_lims) fnix(ix_lim,iy,ifld)=0.
           if (nxpt==2 .and. ixmxbcl==1) fnix(ixrb(1)+1,iy,ifld)=0.
   81    continue


*  -- compute fniy  --

         do 83 iy = j1, j5
            do 82 ix = i4, i8
               if (zi(ifld).eq.0.) then #inertial gas must follow ion index
                  fniy(ix,iy,ifld) = fngy(ix,iy,ifld-1)
               else
                  if (methny .eq. 2) then   # central differencing
                     t2 = ( niy0(ix,iy,ifld) + niy1(ix,iy,ifld) ) / 2

                  elseif (methny .eq. 3) then   # upwind differencing

                     if( vy(ix,iy,ifld) .ge. 0.) then
                        t2 = niy0(ix,iy,ifld)
                     else
                        t2 = niy1(ix,iy,ifld)
                     endif

                  elseif (methny .eq. 6) then   # log central differencing
                     t2 = exp( 0.5*
     .                   (log(niy0(ix,iy,ifld))+log(niy1(ix,iy,ifld))) )

                  else    # interp. ave or harmonic ave depending on wind*grad

                     t0 = ( niy0(ix,iy,ifld)*gy(ix,iy  ) +
     .                    niy1(ix,iy,ifld)*gy(ix,iy+1) ) / 
     .                    ( gy(ix,iy)+gy(ix,iy+1) )
                     t1 = ( gy(ix,iy)+gy(ix,iy+1) ) * niy0(ix,iy,ifld)*
     .                    niy1(ix,iy,ifld) / ( cutlo + niy0(ix,iy,ifld)*
     .                    gy(ix,iy+1) + niy1(ix,iy,ifld)*gy(ix,iy) )
                     if( (niy0(ix,iy,ifld)-niy1(ix,iy,ifld))*
     .                    vy(ix,iy,ifld) .ge. 0.) then
                        t2 = t0
                     else
                        t2 = t1
                     endif
                  
                  endif
               
                  fniy(ix,iy,ifld) = cnfy*vy(ix,iy,ifld)*sy(ix,iy)*t2
                  fniycb(ix,iy,ifld) = cnfy*vycb(ix,iy,ifld)*sy(ix,iy)*t2
                  if (vy(ix,iy,ifld)*(ni(ix,iy,ifld)-ni(ix,iy+1,ifld))
     .                                                      .lt. 0.) then
                     fniy(ix,iy,ifld) = fniy(ix,iy,ifld)/( 1 +
     .                               (nlimiy(ifld)/ni(ix,iy+1,ifld))**2 +
     .                               (nlimiy(ifld)/ni(ix,iy  ,ifld))**2 )
                  endif
c...  Note: nonorthogonality comes in through calc. of vy
               endif
 82         continue
 83      continue

c ... cosmetic setting of fniy - not used         
         do ix = i4, i8
            fniy(ix,ny+1,ifld) = 0.0e0
         enddo

 104  continue

c ... Add rad flux of 4th order diff operator; damp grid-scale oscillations
      do ifld = 1, nfsp
        if (abs(dif4order(ifld)) > 1.e-50) then
          do iy = j2p, j5m   #limits to range iy=1:ny-1 for fniy4ord
            iym1 = max(iy-1,0)
            iyp1 = min(iy+1,ny+1)
            iyp2 = min(iy+2,ny+1)
            do ix = i4, i8
              dndym1 = (ni(ix,iy,ifld)-ni(ix,iym1,ifld))*gyf(ix,iym1)
              dndy0 = (ni(ix,iyp1,ifld)-ni(ix,iy,ifld))*gyf(ix,iy)
              dndyp1 = (ni(ix,iyp2,ifld)-ni(ix,iyp1,ifld))*gyf(ix,iyp1)
              d2ndy20 = (dndy0 - dndym1)*gy(ix,iy)
              d2ndy2p1 = (dndyp1 - dndy0)*gy(ix,iyp1)
              d3ndy3 = (d2ndy2p1 - d2ndy20)*gyf(ix,iy)
              fniy4ord(ix,iy,ifld) = dif4order(ifld)*d3ndy3*sy(ix,iy)/
     .                                                  gyf(ix,iy)**2
              fniy(ix,iy,ifld) = fniy(ix,iy,ifld) + fniy4ord(ix,iy,ifld)
            enddo
          enddo
        endif
      enddo

c ... Setup a correction to surface-flux for grad_B and grad_P effects at iy=0
      do ifld = 1, nfsp
        do ix = i4, i8
           fniycbo(ix,ifld) = 0.0
        enddo
         do ix = i4, i8
            fniycbo(ix,ifld) = ( ni(ix,0,ifld)*sy(ix,0) ) *
     .                         ( (1-cfniybbo)*cfybf*vycb(ix,0,ifld) -
     .                            cfniydbo*(1-cfydd)*vycp(ix,0,ifld) )
         enddo
      enddo

c ... Normalize core flux to zero to avoid introducing artifical core source/sink
      do ifld = 1, nfsp
          if (isfniycbozero(ifld) .gt. 0) then 
            fniycboave = 0
            sycore = 0
            do ix = ixpt1(1)+1, ixpt2(1)
              fniycboave = fniycboave + fniycbo(ix, ifld)
            end do
            corecells =  (ixpt2(1) - ixpt1(1))
            if (geometry == 'dnull') then
                do ix = ixpt1(2)+1, ixpt2(2)
                  fniycboave = fniycboave + fniycbo(ix, ifld)
                end do
                corecells =  corecells + (ixpt2(2) - ixpt1(2))
            end if 
c ... TODO: Add double-null fix here (now only does one half-mesh...)            
            fniycboave = fniycboave / corecells
            do ix = ixpt1(1)+1, ixpt2(1)
              fniycbo(ix, ifld) = fniycbo(ix, ifld) - isfniycbozero(ifld)*fniycboave
            end do
          else if (isfniycbozero(ifld) .lt. 0) then 
            do ix = ixpt1(1)+1, ixpt2(1)
              fniycbo(ix, ifld) = 0
            end do
            if (geometry == 'dnull') then
                do ix = ixpt1(2)+1, ixpt2(2)
                  fniycbo(ix, ifld) = 0
                end do
            end if
          end if
      end do
             


c----------------------------------------------------------------------c
c          SCALE SOURCE TERMS FROM MONTE-CARLO-NEUTRALS MODEL
c
c     These sources are used in the residuals (resco,resmo,resee,resei)
c     so the call to scale_mcn must occur BEFORE these residuals are
c     evaluated.  Since they scale with fnix at the divertor plates,
c     the call to scale_mcn must occur AFTER fnix has been calculated.


      if (ismcnon .ne. 0) then
c 	     write(*,*) 'TEST ISMCNON START: ismcnon=',ismcnon
c         call scale_mcn
         call scale_mcnsor
      endif
c----------------------------------------------------------------------c

*  -- compute the residual if isnion = 1 --

      do ifld = 1, nfsp
       do 86 iy = j2, j5
         do 85 ix = i2, i5
	   if(isnionxy(ix,iy,ifld) == 1) then
              resco(ix,iy,ifld) = 
     .           snic(ix,iy,ifld)+sniv(ix,iy,ifld)*ni(ix,iy,ifld) +
     .           volpsor(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psor(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psorxr(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psori(ix,iy,ifld) -
     .           nuvl(ix,iy,ifld)*vol(ix,iy)*ni(ix,iy,ifld) +
     .           voljcsor(ix,iy)/qe
           endif
c           if (ifld .ne. iigsp) then
	       if(zi(ifld) .ne. 0) then # IJ 2016 skip if neutral zi(ifld)=0
              resco(ix,iy,ifld) = resco(ix,iy,ifld) + cmneut * uesor_ni(ix,iy,ifld)
           else # IJ 2016 zi==0, assume neutral->ifld and ion->ifld-1
              resco(ix,iy,ifld) = resco(ix,iy,ifld) - cmneut * uesor_ni(ix,iy,ifld-1)
           endif
   85    continue
   86  continue

       do 302 iy = j2, j5
         do 301 ix = i2, i5
	       if(isnionxy(ix,iy,ifld) == 1) then
              ix1 = ixm1(ix,iy)
	        if(zi(ifld) .ne. 0) then # IJ 2016 skip if neutral  zi(ifld)=0
                 resco(ix,iy,ifld) = resco(ix,iy,ifld)
     .                                   - ( (fnix(ix,iy,ifld) - fnix(ix1,iy, ifld))
     .                            + fluxfacy*(fniy(ix,iy,ifld) - fniy(ix,iy-1,ifld)) )
              else ## zi==0
                 resco(ix,iy,ifld) = resco(ix,iy,ifld)
     .              - cfneutdiv*cfneutdiv_fng*( (fnix(ix,iy,ifld)-fnix(ix1,iy, ifld))
     .                               + fluxfacy*(fniy(ix,iy,ifld)-fniy(ix,iy-1,ifld)) )
           
c ... IJ 2016/10/19 add MC neutral flux if flags set
                 if (get_neutral_moments .and. cmneutdiv_fng .ne. 0.0) then 
                    jfld=1  ## assume main ions in ifld=1
                    sng_ue(ix,iy,jfld) = - ( (fngx_ue(ix,iy,jfld) - fngx_ue(ix1,iy, jfld))
     .                        +   fluxfacy*(fngy_ue(ix,iy,jfld) - fngy_ue(ix,iy-1,jfld)) )
     .                        *( (ng(ix,iy,jfld)*ti(ix,iy))/(ng(ix,iy,jfld)*ti(ix,iy)) )
c                   if (ix .eq. 1 .and. iy .eq. 1) write(*,*) 'sng_ue', ifld, jfld
                    resco(ix,iy,ifld) = resco(ix,iy,ifld) + 
     .                                  cmneutdiv*cmneutdiv_fng*sng_ue(ix,iy,jfld)
                 endif
              endif
           endif
 301     continue
 302   continue
       enddo       # end of ifld loop

*********************************************************************
c  Here we do the neutral gas diffusion model
c  The diffusion is flux limited using the thermal flux
**********************************************************************

ccc         if(isngon .eq. 1) call neudif

*****************************************************************
*****************************************************************
*  Here starts the old MOMBAL_B2
*****************************************************************


*  ---------------------------------------------------------------------
*  loop over all species.
*  ---------------------------------------------------------------------

      do 105 ifld = 1, nusp
      if(isupon(ifld) .eq. 0) goto 105
*     ------------------------------------------------------------------
*     compute the residual.
*     ------------------------------------------------------------------

*  -- evaluate flox and conx --

         do 91 iy = j4, j8
            flox(0,iy) = 0.0e0
            conx(0,iy) = 0.0e0
            do 90 ix = i2, i6
               ix1 = ixm1(ix,iy)
               if (isimpon.ge.5 .and. ifld.eq.1) then
                   #up(,,1) is total mass vel, whereas uu(,,i) for each ion
                  uuv =0.5*( (up(ix1,iy,ifld)*rrv(ix1,iy)+
     .                        up(ix,iy,ifld)*rrv(ix,iy)) +
     .                     (v2(ix1,iy,ifld)+v2(ix,iy,ifld))*rbfbt(ix,iy)-
     .                     (vytan(ix1,iy,ifld)+vytan(ix,iy,ifld)) )
               else
                  uuv = 0.5 * (uu(ix1,iy,ifld)+uu(ix,iy,ifld))
               endif
               flox(ix,iy) = cmfx * nm(ix,iy,ifld) * uuv *
     .                          vol(ix,iy) * gx(ix,iy)
ccc Distance between veloc. cell centers:
               if (isgxvon .eq. 0) then     # dx(ix)=1/gx(ix)
                 conx(ix,iy) = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .                            * gx(ix,iy)
               elseif (isgxvon .eq. 1) then # dx(ix)=.5/gxf(ix-1) + .5/gxf(ix)
                 conx(ix,iy) = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .               * 2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
               endif
   90       continue
   91    continue

*  -- evaluate floy and cony without averaging over two ix cells --

         do 93 iy = j1, j5
            if (nxpt == 1 .or. iy <= iysptrx1(1)) then
              iysepu = iysptrx1(1)
              if (ndomain > 1) iysepu = iysptrxg(mype+1)  # and ixpt1,2u??
              ixpt1u = ixpt1(1)
              ixpt2u = ixpt2(1)
            else  # nxpt=2 and iy > iysptrx1(1), use second separatrix
              iysepu = iysptrx1(2)
              ixpt1u = ixpt1(2)
              ixpt2u = ixpt2(2)
            endif
            do 92 ix = i4, i8
               ix2 = ixp1(ix,iy)
               ix4 = ixp1(ix,iy+1)
               cony(ix,iy) = .5 * syv(ix,iy) *
     .                       (ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                           visy(ix,iy+1,ifld)*gy(ix,iy+1)) +
     .                        ave(visy(ix2,iy,ifld)*gy(ix2,iy),
     .                            visy(ix4,iy+1,ifld) * gy(ix4,iy+1)))
               if (iy==iysepu .and. (ix==ixpt1u .or. ix==ixpt2u)) then
                 cony(ix,iy) = syv(ix,iy) *
     .                       ( ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                             visy(ix,iy+1,ifld)*gy(ix,iy+1)) )
                 floy(ix,iy) = (cmfy/2) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) )*
     .                                         vy(ix,iy,ifld)
                 if (ifld==1) then  # add user-specified convection
                   floy(ix,iy) = floy(ix,iy)+(cmfy/2)*syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) )*
     .                                                 vyup_use(ix,iy)
                 endif
               elseif(isugfm1side == 1 .and. zi(ifld) == 0.) then
                 floy(ix,iy) = (cmfy/4) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vy(ix,iy,ifld) + vy(ix,iy,ifld) )
               else
                 floy(ix,iy) = (cmfy/4) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vy(ix,iy,ifld) + vy(ix2,iy,ifld) )
                 if (ifld==1) then  # add user-specified convection
                    floy(ix,iy) = floy(ix,iy)+(cmfy/4)*syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vyup_use(ix,iy) + vyup_use(ix2,iy) )
                 endif
               endif
	     if(ishavisy == 1) then
               cony(ix,iy) = .5 * syv(ix,iy) *
     .                       (ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                           visy(ix,iy+1,ifld)*gy(ix,iy+1)) +
     .                        ave(visy(ix2,iy,ifld)*gy(ix2,iy),
     .                            visy(ix4,iy+1,ifld) * gy(ix4,iy+1)))
             else
               cony(ix,iy) = .25 * cfaccony*syv(ix,iy) *
     .                       ( visy(ix,iy,ifld)*gy(ix,iy) +
     .                         visy(ix,iy+1,ifld)*gy(ix,iy+1) +
     .                         visy(ix2,iy,ifld)*gy(ix2,iy) +
     .                         visy(ix4,iy+1,ifld)*gy(ix4,iy+1) )
             endif

   92       continue
   93    continue

*  -- compute the momentum transport --

         call fd2tra (nx,ny,flox,floy,conx,cony,
     .                up(0:nx+1,0:ny+1,ifld),fmix(0:nx+1,0:ny+1,ifld),
     .                fmiy(0:nx+1,0:ny+1,ifld),1, methu)

      if (isnonog .eq. 1) then

c     Compute y-component fmixy of nonorthogonal diffusive momentum flux.
c     The convective component is already already added through uu(ix,iy).
c     Average fym, etc in ix to get staggered velocity-grid values fymv, etc.
c     The density-stencil dxnog has to be averaged as well.
         do 96 iy = j2, j5
            iy1 = max(iy-1,0)
            do 97 ix = i2, i5+1    # ixp1(i5,iy)
               ix1 = ixm1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               grdnv = ( 
     .                  fymv (ix,iy,1)*up(ix ,iy1 ,ifld)+
     .                  fy0v (ix,iy,1)*up(ix ,iy  ,ifld)+
     .                  fypv (ix,iy,1)*up(ix ,iy+1,ifld)+
     .                  fymxv(ix,iy,1)*up(ix3,iy1 ,ifld)+
     .                  fypxv(ix,iy,1)*up(ix5,iy+1,ifld)-
     .                  fymv (ix,iy,0)*up(ix3,iy1 ,ifld)-
     .                  fy0v (ix,iy,0)*up(ix1,iy  ,ifld)-
     .                  fypv (ix,iy,0)*up(ix5,iy+1,ifld)-
     .                  fymxv(ix,iy,0)*up(ix ,iy1 ,ifld)-
     .                  fypxv(ix,iy,0)*up(ix ,iy+1,ifld) )*2/
     .                              (dxnog(ix,iy)+dxnog(ix1,iy))
               if (isgxvon .eq. 0) then
                  fmixy(ix,iy,ifld) = cfvisxy(ifld)*visy(ix,iy,ifld) *
     .              ( grdnv/cos(0.5*(angfx(ix1,iy)+angfx(ix,iy))) - 
     .               (up(ix,iy,ifld) - up(ix1,iy,ifld))*gx(ix,iy) ) *
     .              0.5*(sx(ix1,iy)+sx(ix,iy))
               elseif (isgxvon .eq. 1) then
                  fmixy(ix,iy,ifld) = cfvisxy(ifld)*visy(ix,iy,ifld) *
     .              ( grdnv/cos(0.5*(angfx(ix1,iy)+angfx(ix,iy))) - 
     .               (up(ix,iy,ifld) - up(ix1,iy,ifld))*
     .                     ( 2*gxf(ix,iy)*gxf(ix1,iy) /
     .                        (gxf(ix,iy)+gxf(ix1,iy)) ) ) *
     .                     0.5*(sx(ix1,iy)+sx(ix,iy))
               endif
c...  Now flux limit with flalfvgxy if ifld=2
               if (ifld==2) then
                 t0 = max(tg(ix,iy,1),tgmin*ev)
                 vtn = sqrt(t0/mg(1))
                 qfl = flalfvgxya(ix)*0.5*(sx(ix,iy)+sx(ix1,iy))*vtn**2*
     .                                        nm(ix,iy,ifld) + cutlo
                 fmixy(ix,iy,ifld) = fmixy(ix,iy,ifld) /
     .                             sqrt(1+(fmixy(ix,iy,ifld)/qfl)**2)
               endif

 97         continue
 96      continue
      endif

c...  Compute viscous drag from nonuniform B-field, then add to smoc
      if (isupdrag .eq. 1 .and. ifld .eq. 1) then
        do iy = j2, j5
          do ix = i2, i5
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
c ...   First, the short mfp drag 
            b_ctr = 0.5*(btot(ix,iy)+btot(ix2,iy)) 
                    # derviatives dbds_m and dbds_p are one-sided
            dbds_m = (btot(ix,iy) - btot(ix1,iy))*
     .                                      gxf(ix1,iy)*rrv(ix1,iy)
            dbds_p = (btot(ix2,iy) - btot(ix,iy))*
     .                                        gxf(ix,iy)*rrv(ix,iy)
            eta_h0 = visx(ix,iy,1)/b_ctr**2.5
            eta_hm = 0.5*(visx(ix,iy,1)+visx(ix1,iy,1))/
     .                                              btot(ix,iy)**2.5
            eta_hp = 0.5*(visx(ix2,iy,1)+visx(ix,iy,1))/
     .                                              btot(ix2,iy)**2.5
            drag_1 = -2*eta_h0*0.5*(up(ix2,iy,1)-up(ix1,iy,1))*
     .                            (btot(ix2,iy)-btot(ix,iy))*
     .                                   (gxf(ix,iy)*rrv(ix,iy))**2
            drag_2 = up(ix,iy,1)*(eta_hp*dbds_p - eta_hm*dbds_m)*
     .                                        gxf(ix,iy)*rrv(ix,iy)
c ...   now for the trapped particle drag (sloppy)
            nu_ii = ni(ix,iy,1)*(2*mp/mi(1))**0.5/
     .                                   (3e12*(ti(ix,iy)*ev)**1.5)
            drag_3 = -mi(1)*ni(ix,iy,1)*up(ix,iy,1)*nu_ii*frac_pt
            mf_path = (2*ti(ix,iy)/mi(1))**0.5 / nu_ii
            frac_col = 1 / (1 + (mf_path/con_leng)**2)
            smoc(ix,iy,1) = smoc(ix,iy,1) + (0.6666667*frac_col*
     .                      b_ctr**2.5*(drag_1 + drag_2) +
     .                      (1 - frac_col)*drag_3) * volv(ix,iy)
          enddo
        enddo
      endif
      
c...  Compute total viscosity for nonuniform B-field; put in visvol_v,q
      if (cfvisxneov+cfvisxneoq > 0.) call upvisneo

c...  Now fix the fluxes touching the x-point(s):


         do k = 1, nxpt   # loop over all x-points
           k1 = k      # region argument of ixpt1 that touches this x-point
           k2 = k-1    # region argument of ixpt2 that touches this x-point
           if (k==1) k2 = nxpt

           if (nxpt==2) then      # set ghxpt,gvxpt,sxyxpt for full double null
              if (k==1) then      # this is the lower x-point
                 ghxpt = ghxpt_lower
                 gvxpt = gvxpt_lower
                 sxyxpt = sxyxpt_lower
              elseif (k==2) then  # this is the upper x-point
                 ghxpt = ghxpt_upper
                 gvxpt = gvxpt_upper
                 sxyxpt = sxyxpt_upper
              endif
           endif

	   if( ((2*(yc-iysptrx1(k1))-1)/4 .le. 1) .or. j1 == 0 ) then
           if( ((2*(xc-ixpt1(k1))-1)/4)*((2*(xc-ixpt2(k2))-1)/4).eq.0 .or. 
     .                                                        i1.eq.0 ) then
           if(isnfmiy .eq. 1) then

           fmiy(ixpt1(k1),iysptrx1(k1),ifld) = 0.
           fmiy(ixpt2(k2),iysptrx2(k2),ifld) = 0.
           nixpt(ifld,k1) = 0.125 * ( 
     .           ni(ixpt1(k1),iysptrx1(k1)  ,ifld) + ni(ixpt1(k1)+1,iysptrx1(k1)  ,ifld)
     .         + ni(ixpt1(k1),iysptrx1(k1)+1,ifld) + ni(ixpt1(k1)+1,iysptrx1(k1)+1,ifld)
     .         + ni(ixpt2(k2),iysptrx2(k2)  ,ifld) + ni(ixpt2(k2)+1,iysptrx2(k2)  ,ifld)
     .         + ni(ixpt2(k2),iysptrx2(k2)+1,ifld) + ni(ixpt2(k2)+1,iysptrx2(k2)+1,ifld) )
           visyxpt(ifld,k1) = 0.125 * ( 
     .        visy(ixpt1(k1),iysptrx1(k1)  ,ifld) + visy(ixpt1(k1)+1,iysptrx1(k1)  ,ifld)
     .      + visy(ixpt1(k1),iysptrx1(k1)+1,ifld) + visy(ixpt1(k1)+1,iysptrx1(k1)+1,ifld)
     .      + visy(ixpt2(k2),iysptrx2(k2)  ,ifld) + visy(ixpt2(k2)+1,iysptrx2(k2)  ,ifld)
     .      + visy(ixpt2(k2),iysptrx2(k2)+1,ifld) + visy(ixpt2(k2)+1,iysptrx2(k2)+1,ifld) )
           upxpt(ifld,k1) = 0.25 * (
     .           up(ixpt1(k1),iysptrx1(k1)  ,ifld) + up(ixpt2(k2),iysptrx2(k2)  ,ifld)
     .         + up(ixpt1(k1),iysptrx1(k1)+1,ifld) + up(ixpt2(k2),iysptrx2(k2)+1,ifld) )
           vyvxpt(ifld,k1) = (0.707*0.25) * (
     .           vy(ixpt1(k1)  ,iysptrx1(k1),ifld) - vy(ixpt2(k2)  ,iysptrx2(k2),ifld)
     .         - vy(ixpt1(k1)+1,iysptrx1(k1),ifld) + vy(ixpt2(k2)+1,iysptrx2(k2),ifld) )
           vyhxpt(ifld,k1) = (0.707*0.25) * (
     .         - vy(ixpt1(k1)  ,iysptrx1(k1),ifld) + vy(ixpt2(k2)  ,iysptrx2(k2),ifld)
     .         - vy(ixpt1(k1)+1,iysptrx1(k1),ifld) + vy(ixpt2(k2)+1,iysptrx2(k2),ifld) )
cccMER The convective contributions to fmihxpt and fmivxpt seem to have an
cccMER erroneous multiplicative factor -1/2 (from original code) ???
           fmihxpt(ifld,k1) = cfnfmiy*( - cmfy*(mi(ifld)/2)*sxyxpt*nixpt(ifld,k1)*
     .                     vyhxpt(ifld,k1)*upxpt(ifld,k1)
     .                   - sxyxpt*visyxpt(ifld,k1)*(up(ixpt2(k2),iysptrx2(k2)+1,ifld)
     .                     - up(ixpt1(k1),iysptrx1(k1)+1,ifld))*ghxpt )
           fmivxpt(ifld,k1) = cfnfmiy*(- cmfy*(mi(ifld)/2)*sxyxpt*nixpt(ifld,k1)*
     .                     vyvxpt(ifld,k1)*upxpt(ifld,k1)
     .                   - sxyxpt*visyxpt(ifld,k1)*(up(ixpt2(k2),iysptrx2(k2),ifld)
     .                     - up(ixpt1(k1),iysptrx1(k1),ifld))*gvxpt )
           smoc(ixpt1(k1),iysptrx1(k1)+1,ifld) = smoc(ixpt1(k1),iysptrx1(k1)+1,ifld) 
     .                                - fmihxpt(ifld,k1)
           smoc(ixpt2(k2),iysptrx2(k2)+1,ifld) = smoc(ixpt2(k2),iysptrx2(k2)+1,ifld) 
     .                                + fmihxpt(ifld,k1)
           smoc(ixpt1(k1),iysptrx1(k1)  ,ifld) = smoc(ixpt1(k1),iysptrx1(k1)  ,ifld) 
     .                                - fmivxpt(ifld,k1)
           smoc(ixpt2(k2),iysptrx2(k2)  ,ifld) = smoc(ixpt2(k2),iysptrx2(k2)  ,ifld) 
     .                                + fmivxpt(ifld,k1)
         
           endif # end if-test on isnfmiy
           endif # end if-test on xc
           endif # end if-test on yc

         enddo # end do-loop over nxpt x-points

*  -- source term and pressure gradient --

         do 99 iy = j2, j5
            do 98 ix = i2, i5
               ix2 = ixp1(ix,iy)
               if (zi(ifld) .ne. 0) then  # additions only for charged ions
                  dp1 =  cngmom(ifld)*(1/fac2sp)*
     .                      ( ng(ix2,iy,1)*tg(ix2,iy,1)-
     .                        ng(ix ,iy,1)*tg(ix ,iy,1) ) 
                  resmo(ix,iy,ifld) = 0.
                  resmo(ix,iy,ifld) = 
     .                  smoc(ix,iy,ifld)
     .                + smov(ix,iy,ifld) * up(ix,iy,ifld)
     .                - cfneut * cfneutsor_mi * sx(ix,iy) * rrv(ix,iy) * dp1
     .                - cfneut * cfneutsor_mi * cmwall(ifld)*0.5*(ng(ix,iy,1)+ng(ix2,iy,1))
     .                      * mi(ifld)*up(ix,iy,ifld)*0.5
     .                      *(nucx(ix,iy,1)+nucx(ix2,iy,1))*volv(ix,iy)
     .                + cmneut * cmneutsor_mi * uesor_up(ix,iy,ifld)
     .                + cfmsor*(msor(ix,iy,ifld) + msorxr(ix,iy,ifld)) #### IJ 2017: needs *cfneut for multi-charge state ions & MC neutrals?
     .                + volmsor(ix,iy,ifld)
     .                + cfvisxneov*visvol_v(ix,iy,ifld)
     .                + cfvisxneoq*visvol_q(ix,iy,ifld)
c  Add drag with cold, stationary impurity neutrals
                  if (ifld > nhsp) then
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld) - cfupimpg*
     .                    0.25*mi(ifld)*(ni(ix,iy,ifld)+ni(ix2,iy,ifld))*
     .                     ( nucxi(ix,iy,ifld)+nucxi(ix2,iy,ifld)+
     .                       nueli(ix,iy,ifld)+nueli(ix2,iy,ifld) )*
     .                       up(ix,iy,ifld)*volv(ix,iy)
                  endif
               endif

               if (isupgon(1) .eq. 1) then

c     If we are solving the parallel neutral mom eq. we need different/addtnl
c     source terms. Beware that cngmom and cmwall should be zero so that the
c     main ions do not get coupled twice to the neutrals!
c     The CX, ionization, recomb. friction for the parallel momentum eqs
c     for the neutrals and main ions are included here.
c     Assumes the neutrals have index iigsp and corresponding ions index 1

                  if (ifld .eq. 1) then
c     The main ions, momentum coupling:
                     resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                   + cfneut * cfneutsor_mi * cfupcx*0.25*volv(ix,iy)*
     .                       (nucx(ix,iy,1)+nucx(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                       (up(ix,iy,iigsp)-up(ix,iy,1))
     .                   + cfneut * cfneutsor_mi * 0.25*volv(ix,iy)*
     .                       (  (nuiz(ix,iy,1)+nuiz(ix2,iy,1))*
     .                          (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                          up(ix,iy,iigsp)
     .                        - (nurc(ix,iy,1)+nurc(ix2,iy,1))*
     .                          (nm(ix,iy,1)+nm(ix2,iy,1))*up(ix,iy,1)
     .                       )
                  elseif ((isupgon(1) .eq. 1) .and. ifld .eq. iigsp) then
c     The neutral species, momentum coupling AND other source terms:
                      resmo(ix,iy,iigsp) =   # TR resmo(ix,iy,ifld) #IJ 2016
     .                    - cmneut * cmneutsor_mi * uesor_up(ix,iy,1) 
     .                    -sx(ix,iy) * rrv(ix,iy) * 
     .                       cpgx*( cftiexclg*(ni(ix2,iy,iigsp)*ti(ix2,iy)-
     .                                ni(ix,iy,iigsp)*ti(ix,iy))+
     .                              (1.0-cftiexclg)*
     .                               (ni(ix2,iy,iigsp)*tg(ix2,iy,1)-
     .                                ni(ix,iy,iigsp)*tg(ix,iy,1)) ) 
     .                    -cfupcx*0.25*volv(ix,iy)*
     .                       (nucx(ix,iy,1)+nucx(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                       (up(ix,iy,iigsp)-up(ix,iy,1))
     .                    -0.25*volv(ix,iy)*(
     .                       (nuiz(ix,iy,1)+nuiz(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                                            up(ix,iy,iigsp)
     .                       -(nurc(ix,iy,1)+nurc(ix2,iy,1))*
     .                       (nm(ix,iy,1)+nm(ix2,iy,1))*up(ix,iy,1) )
                  endif
               endif
 98         continue
 99      continue

*  -- divergence of momentum flow --

         if (isnonog.eq.1) then
            do 3051 iy = j2, j5
               do 3061 ix = i2, i5
                  ix2 = ixp1(ix,iy)
c ... IJ 2016/10/10 use cfneutdiv_fmg multiplier for neutrals 
c                 if (ifld .ne. iigsp) then
	          if(zi(ifld) > 1.e-20) then # IJ 2016; depends if ion or neut
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                 + (fmixy(ix2,iy,ifld) - fmixy(ix,iy,ifld))
                  else
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                 + cfneutdiv*cfneutdiv_fmg*(fmixy(ix2,iy,ifld) - fmixy(ix,iy,ifld))
c***	IJ 2017/09/21: Need to add similar fmgxy calculation for MC neutrals on nonorthogonal mesh *** 
                  endif
 3061          continue
 3051       continue
         endif

         do 305 iy = j2, j5
            do 306 ix = i2, i5
               ix2 = ixp1(ix,iy)
c IJ 2016/10/10 add cfneutdiv_fmg multiplier for neutrals to control fraction of momentum to add 
c               if (ifld .ne. iigsp) then
c ... IJ 2016 resmo contrib changes if ion or neut
	       if(zi(ifld) > 1.e-20) then
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                                 - (fmix(ix2,iy,ifld) - fmix(ix,iy  ,ifld)
     .                        + fluxfacy*(fmiy(ix ,iy,ifld) - fmiy(ix,iy-1,ifld)) )
               else
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                      - cfneutdiv*cfneutdiv_fmg*(fmix(ix2,iy,ifld) - fmix(ix,iy  ,ifld)
     .                        + fluxfacy*(fmiy(ix ,iy,ifld) - fmiy(ix,iy-1,ifld)) )
                 if(cmneutdiv_fmg .ne. 0.0) then 
                    jfld=1
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                    - cmneutdiv*cmneutdiv_fmg*( (fmgx_ue(ix2,iy,jfld) - fmgx_ue(ix,iy  ,jfld))
     .                        + fluxfacy*(fmgy_ue(ix ,iy,jfld) - fmgy_ue(ix,iy-1,jfld)) )
     .                        * (ni(ix,iy,ifld)*ti(ix,iy))/(ni(ix,iy,ifld)*ti(ix,iy))
                 endif
              endif
  306       continue
  305    continue

c  -- Include frictional drag in parallel flow here if isofric=1; otherwise
c  -- it is included in frici from mombal or mombalni

        if (isofric.eq.1 .and. nusp .gt.1) then
*  -- w0 now accumulates friction coefficient --
*     -- set w2 = vol*ti**(-1.5) --
         do iy = j1, j6
           do ix = i1, i6
             fricnrl(ix,iy,ifld) = 0.  #diagnostic ~ ni*mi*nu*(up1-up2)
             w0(ix,iy) = 0.0e0
             w2(ix,iy) = vol(ix,iy) / (ti(ix,iy)*sqrt(ti(ix,iy)))
           enddo
         enddo

*  -- consider all other species --

         do jfld = 1, nusp
           if (jfld .ne. ifld) then
*     -- common factor in collision frequency --
             awoll = zi(ifld)**2 * zi(jfld)**2 *
     .            (qe**4/(12*pi**2*eps0**2)) *
     .            sqrt (2*pi*mi(ifld)*mi(jfld)/(mi(ifld)+mi(jfld)))

*     -- frictional coupling --
             do iy = j1, j6
               do ix = i1, i5
                 ix2 = ixp1(ix,iy)
                 t0 = ni(ix,iy,ifld) * ni(ix,iy,jfld) * w2(ix,iy)
                 t1 = ni(ix2,iy,ifld)*ni(ix2,iy,jfld)*w2(ix2,iy)
                 awll = awoll*loglambda(ix,iy)
                 tv  = awll*(t0+t1)/2
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld) +
     .                        tv * (up(ix,iy,jfld)-up(ix,iy,ifld))
                 fricnrl(ix,iy,ifld) = fricnrl(ix,iy,ifld) +
     .                tv*(up(ix,iy,jfld)-up(ix,iy,ifld))/vol(ix,iy)
                 w0(ix,iy) = w0(ix,iy) - tv
               enddo
             enddo
           endif
         enddo
        endif   # if test on isofric.eq.1

 105  continue


*****************************************************************
*****************************************************************
*  Here starts the old ENEBAL
*****************************************************************
*  ---------------------------------------------------------------------
*  compute temperature conductances.
*  ---------------------------------------------------------------------
*  -- initialize to 0 --

      do 708 iy = j1, j6
         do 707 ix = i1, i6
            floxe(ix,iy) = 0.0e0
            floxi(ix,iy) = 0.0e0
            floye(ix,iy) = 0.0e0
            floyi(ix,iy) = 0.0e0
            feiycbo(ix) = 0.0e0
            feeycbo(ix) = 0.0e0
            w0(ix,iy) = 0.0e0
            w1(ix,iy) = 0.0e0
  707    continue
  708 continue
      wvh = 0.0e0    #wvh is 3D vector (ix,iy,nusp)

*  -- compute conxe and conxi --

*     (The computation of conxe involves a flux limit)

      do 121 iy = j4, j8
         do 120 ix = i1, i5
            ix2 = ixp1(ix,iy)
            t0 = max (te(ix,iy), temin*ev)
            t1 = max (te(ix2,iy), temin*ev)
            vt0 = sqrt(t0/me)
            vt1 = sqrt(t1/me)
            wallfac = 1.
            do jx = 1, nxpt
               if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                (ix==ixrb(jx).and.ixmxbcl==1) )
     .              .and. (isplflxl==0) ) wallfac = flalfepl/flalfe
            enddo
            qfl = wallfac*flalfe * sx(ix,iy) * rrv(ix,iy) *
     .          (ne(ix,iy)*vt0*t0 + ne(ix2,iy)*vt1*t1) / 2
            csh = sx(ix,iy) * hcxe(ix,iy) * gxf(ix,iy)
	    lxtec = 0.5*(te(ix,iy)+te(ix2,iy)) /
     .             ( abs(te(ix,iy)-te(ix2,iy))*gxf(ix,iy) + 100.*cutlo )
	    qsh = csh * (te(ix,iy)-te(ix2,iy)) * (1. + lxtec/lxtemax)
            qr = (1-isflxlde)*abs(qsh/qfl)
            conxe(ix,iy) = (1-isflxlde)*csh / (1 + qr)**2 +
     .                isflxlde*csh /(1 + abs(qsh/qfl)**flgam)**(1/flgam)
            floxe(ix,iy) = floxe(ix,iy) + (sign(qr*qr,qsh)/(1 + qr)**2)
     .                  *flalfea(ix) * sx(ix,iy) * ( ne(ix,iy)*
     .                   rr(ix,iy)*vt0 + ne(ix2,iy)*rr(ix2,iy)*vt1 ) / 2
c.... Now do the ions (hcxi is flux-limited previously when it is built)
          if (isflxldi .ne. 2) then    # Else flux limit done on hcxij
            t0 = max (ti(ix,iy), temin*ev)
            t1 = max (ti(ix2,iy), temin*ev)
            vt0 = sqrt(t0/mi(1))
            vt1 = sqrt(t1/mi(1))
            wallfac = 1.
            do jx = 1, nxpt
               if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                (ix==ixrb(jx).and.ixmxbcl==1) )
     .              .and. (isplflxl==0) ) wallfac = flalfipl/flalfi
            enddo
            qfl = wallfac*flalfia(ix) * sx(ix,iy) * rrv(ix,iy) *
     .                    (ne(ix,iy)*vt0*t0 + ne(ix2,iy)*vt1*t1) / 2
            csh = sx(ix,iy) * hcxi(ix,iy) * gxf(ix,iy)
            lxtic = 0.5*(ti(ix,iy)+ti(ix2,iy)) /
     .             ( abs(ti(ix,iy)-ti(ix2,iy))*gxf(ix,iy) + 100.*cutlo)
	    qsh = csh * (ti(ix,iy)-ti(ix2,iy)) * (1. + lxtic/lxtimax)
            qipar(ix,iy,1) = qsh/(rrv(ix,iy)*sx(ix,iy))
            qr = (1-isflxldi)*abs(qsh/qfl)
            conxi(ix,iy) = (1-isflxldi)*csh / (1 + qr)**2 +
     .                isflxldi*csh / (1 + abs(qsh/qfl)**flgam)**(1/flgam)
            floxi(ix,iy) = floxi(ix,iy) + (sign(qr*qr,qsh)/(1 + qr)**2)
     .                  *flalfia(ix) * sx(ix,iy) *( ne(ix,iy)*rr(ix,iy)*
     .                   vt0 + ne(ix2,iy)*rr(ix2,iy)*vt1 ) / 2 
          else
            conxi(ix,iy) = sx(ix,iy) * hcxi(ix,iy) * gxf(ix,iy)
          endif
  120    continue
         conxe(nx+1,iy) = 0
         conxi(nx+1,iy) = 0
  121 continue

*  -- compute conye and conyi --

      do 123 iy = j1, j5
         do 122 ix = i4, i8
            conye(ix,iy) = sy(ix,iy) * hcye(ix,iy) / dynog(ix,iy)
            conyi(ix,iy) = sy(ix,iy) * hcyi(ix,iy) / dynog(ix,iy)
  122    continue
  123 continue

      do 124 ix = i1, i6
         conye(ix,ny+1) = 0.0e0
         conyi(ix,ny+1) = 0.0e0
  124 continue

*  ---------------------------------------------------------------------
*  compute the strength of convection for te and ti.
*  ---------------------------------------------------------------------
*  includes a correction to account for e-velocity .ne. i-velocity
*  (ccn term), and also a conduction contribution to electron heat flux
*  due to friction (cthe term).  (M.E. Rensink, 07/20/90)
*  ---------------------------------------------------------------------
*     floxe, floxi  contain the cross-derivative terms now
*                        JLM      5/3/90
*  ---------------------------------------------------------------------

      do 126 iy = j4, j8
         do 125 ix = i1, i5  
            ix1 = ixp1(ix,iy)
            ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
            lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
            flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
            floxe(ix,iy) = floxe(ix,iy) + cfcvte*1.25*
     .                  (ne(ix,iy)+ne(ix1,iy))*vex(ix,iy)*sx(ix,iy)
     .                   - cthe*flxlimf*cfjhf*fqp(ix,iy)/ev
  125    continue
         floxe(nx+1,iy) = 0.0e0
  126 continue

c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add 
      do 729 ifld = 1, nfsp
         if ((isupgon(1) .eq. 1) .and. (ifld .eq. iigsp)) then  #neutrals
            do 726 iy = j4, j8
               do 725 ix = i1, i5
                  floxi(ix,iy) = floxi(ix,iy) +
     .                 cftiexclg*cfcvti*2.5*cfneut*cfneutsor_ei*fnix(ix,iy,ifld) 
 725           continue   # next correct for incoming neut pwr = 0
               do jx = 1, nxpt  #if at plate, sub (1-cfloxiplt)*neut-contrib
                 if(ixmnbcl==1) then  #real plate-need for parallel UEDGE
                   iixt = ixlb(jx) #left plate
                   if(fnix(iixt,iy,ifld) > 0.) then
                     floxi(iixt,iy) = floxi(iixt,iy) - (1.-cfloxiplt)*
     .                 cftiexclg*cfcvti*2.5*cfneut*cfneutsor_ei*fnix(iixt,iy,ifld)
                   endif
                 endif
                 if(ixmxbcl==1) then #real plate-need for parallel UEDGE
                   iixt = ixrb(jx) # right plate
                   if(fnix(iixt,iy,ifld) < 0.) then
                     floxi(iixt,iy) = floxi(iixt,iy) - (1.-cfloxiplt)*
     .                 cftiexclg*cfcvti*2.5*cfneut*cfneutsor_ei*fnix(iixt,iy,ifld)
                   endif
                   floxi(ixrb(jx)+1,iy) = 0.0e0  #cosmetic
                 endif
               enddo
 726        continue
         else  #ions
            do 728 iy = j4, j8
               do 727 ix = i1, i5
                  floxi(ix,iy) = floxi(ix,iy) +
     .                           cfcvti*2.5*fnix(ix,iy,ifld)
 727           continue
               floxi(nx+1,iy) = 0.0e0
 728        continue

         endif
 729  continue 

*  -- compute floye and floyi --

      do 129 iy = j1, j5    # note: cfloye usually = 2.5 or 1.5 (ExB turb)
         do 128 ix = i4, i8
            floye(ix,iy) = floye(ix,iy) + (cfloye/2.)*
     .                    (ney0(ix,iy)+ney1(ix,iy))*vey(ix,iy)*sy(ix,iy)
     .                + (vyte_use(ix,iy)+vyte_cft(ix,iy))*0.5*sy(ix,iy)*
     .                     (ney0(ix,iy)+ney1(ix,iy))
            if (iy == 0) then
               feeycbo(ix) =  cfloye*
     .                          ( ne(ix,0)*te(ix,0)*sy(ix,0) ) *
     .                         ( (1-cfeeybbo)*cfybf*veycb(ix,0) -
     .                             cfeeydbo*(1-cfydd)*veycp(ix,0) )
            endif
  128    continue
  129 continue

      do 629 ifld = 1, nfsp
         if ((isupgon(1) .eq. 1) .and. (ifld .eq. iigsp)) then
            do iy = j1, j5
               do ix = i4, i8
                  floyi(ix,iy) = floyi(ix,iy)
     .                 + cftiexclg*cfneut * cfneutsor_ei * 2.5 * fniy(ix,iy,ifld)
               enddo
            enddo
c ...       Make correction at walls to prevent recyc neutrals injecting pwr
            do ix = i4, i8
              if (matwallo(ix) > 0 .and. recycwot(ix,1)>0.) then
                fniy_recy = max(recycwot(ix,1)*fac2sp*fniy(ix,ny,1), 0.)
                floyi(ix,ny) = floyi(ix,ny) + 
     .                 cftiexclg*cfneut*cfneutsor_ei*2.5*(1.-cfloygwall)*fniy_recy
              endif
              if (matwalli(ix) > 0 .and. recycwit(ix,1,1)>0.) then
                fniy_recy = min(recycwit(ix,1,1)*fac2sp*fniy(ix,0,1), 0.)
                floyi(ix,0) = floyi(ix,0) +
     .                 cftiexclg*cfneut*cfneutsor_ei*2.5*(1.-cfloygwall)*fniy_recy
              endif
            enddo 

         else
            do 628 iy = j1, j5 # note: cfloyi usually = 2.5 or 1.5 (ExB turb)
               do 627 ix = i4, i8
                  floyi(ix,iy) = floyi(ix,iy)
     .                            + cfloyi * fniy(ix,iy,ifld)
     .                            + (vyti_use(ix,iy)+vyti_cft(ix,iy))*
     .                                                  0.5*sy(ix,iy)*
     .                              (niy0(ix,iy,ifld)+niy1(ix,iy,ifld))
                  if (iy == 0) then
                     feiycbo(ix) = feiycbo(ix) + cfloyi*fniycbo(ix,ifld)*
     .                                                  ti(ix,0)
                  endif
 627           continue
 628        continue
         endif
 629  continue

c...  Next B x grad(T), first for the ions
      if(abs(cfbgt) .gt. 0 .or. cfeexdbo+cfeixdbo > 0.) then

      do 133 ifld = 1, nfsp
        do 132 iy = j4, j8
           do 131 ix = i1, i5
	     iy1 = max(0,iy-1)
             ix1 = ixp1(ix,iy)
             ix5 = ixp1(ix,iy1)
	     if (iy.eq.0 .or. iy.eq.ny+1) goto 131
c... sknam: grad T from tiv
             temp1 = 4.0*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 )
     .                 temp1 = 4.0*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)

             if (zi(ifld) > 1.e-10) then
                floxibgt(ix,iy,ifld)=(5*sx(ix,iy)/(32*qe*zi(ifld) )) *
     .                           ( ni(ix,iy,ifld) + ni(ix1,iy,ifld) ) *
     .                           ( rbfbt2(ix,iy) + rbfbt2(ix1,iy) ) *
     .                            temp1
             endif
             floxi(ix,iy) = floxi(ix,iy) + cfbgt*floxibgt(ix,iy,ifld)
  131      continue
  132  continue
  133 continue
 
      do 136 ifld = 1, nfsp
        do 135 iy = j1, j5
           do 134 ix = i4, i8
             ix3 = ixm1(ix,iy)
             ix4 = ixm1(ix,iy+1)
             do jx = 1, nxpt
                if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .               (ix==ixrb(jx)+1.and.ixmxbcl==1) ) goto 134
             enddo
c... sknam: grad T from tiv
             temp1 = 4.0*(tiv(ix,iy) - tiv(ix3,iy))*gxc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 )
     .              temp1 = 4.0*(tiv(ix,iy) - tiv(ix3,iy))*gxc(ix,iy)
             if (zi(ifld) > 1.e-10) then
                floyi(ix,iy) = floyi(ix,iy) - 
     .                      cfbgt*( 5*sy(ix,iy) / (32*qe*zi(ifld) )) * 
     .                          ( ni(ix,iy,ifld) + ni(ix,iy+1,ifld) ) *
     .                          ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) ) *
     .                           temp1
             endif
 134      continue
 135   continue
 136  continue

c...  Now B x grad(T) for the electrons

      do 138 iy = j4, j8
         do 137 ix = i1, i5
	     iy1 = max(0,iy-1)
             ix1 = ixp1(ix,iy)
             ix5 = ixp1(ix,iy1)
	     if (iy.eq.0 .or. iy.eq.ny+1) goto 137
c             temp1 = (gtey(ix,iy ) + gtey(ix1,iy ) +
c     .                gtey(ix,iy1) + gtey(ix5,iy1))
c... sknam: grad T from tev
             temp1 = 4.0*(tev(ix,iy) - tev(ix,iy1))*gyc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 )
     .                    temp1 = 4.0*(tev(ix,iy) - tev(ix,iy1))*gyc(ix,iy)
             floxebgt(ix,iy) = ( 5*sx(ix,iy) / (32*qe) ) *
     .                           ( ne(ix,iy) + ne(ix1,iy) ) *
     .                          ( rbfbt2(ix,iy) + rbfbt2(ix1,iy) ) *
     .                           temp1
             floxe(ix,iy) = floxe(ix,iy) - cfbgt*floxebgt(ix,iy)
  137    continue
  138 continue
 
      do 140 iy = j1, j5
	 do 139 ix = i4, i8
	    ix3 = ixm1(ix,iy)
	    ix4 = ixm1(ix,iy+1)
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx)+1.and.ixmxbcl==1) ) goto 139
            enddo
c... sknam: grad T from tev
            temp1 = 4.0*(tev(ix,iy) - tev(ix3,iy))*gxc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
            if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 )
     .               temp1 = 4.0*(tev(ix,iy) - tev(ix3,iy))*gxc(ix,iy)
            floye(ix,iy) = floye(ix,iy) + 
     .                       cfbgt*( 5*sy(ix,iy) / (32*qe) ) * 
     .                         ( ne(ix,iy) + ne(ix,iy+1) ) *
     .                         ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) ) *
     .                          temp1

  139    continue
  140 continue

      endif

c...Add the charge-exhange neutral contributions to ion+neutral temp eq.


         do 141 iy = j4, j8
            do 142 ix = i1, i5
               floxi(ix,iy) = floxi(ix,iy) +
     .          cftiexclg*cfneut*cfneutsor_ei*cngtgx(1)*cfcvti*2.5*fngx(ix,iy,1)
 142        continue
            floxi(nx+1,iy) = 0.0e0
 141        continue

*  --Adds to floyi --

         do 145 iy = j1, j5
            do 144 ix = i4, i8
               floyi(ix,iy) = floyi(ix,iy)
     .             + cftiexclg*cfneut*cfneutsor_ei*cngtgy(1)*2.5*fngy(ix,iy,1)
 144        continue
 145     continue

*  ---------------------------------------------------------------------
*  compute the electron and the ion energy flow.
*  ---------------------------------------------------------------------

*  -- compute the electron energy flow --
      if(isteon .eq. 1) call fd2tra (nx,ny,floxe,floye,conxe,conye,
     .                               te,feex,feey,0,methe)

*  -- compute the ion thermal energy flow --
      if(istion .eq. 1) call fd2tra (nx,ny,floxi,floyi,conxi,conyi,
     .                               ti,feix,feiy,0,methi)

c  -- Add rad flux of 4th order diff operator; damp grid-scale oscillations
      if (abs(kye4order)>1.e-50 .or. abs(kyi4order)>1.e-50) then
        do iy = j2p, j5m   #range to iy=1:ny-1 for feey4ord,feiy4ord
          iym1 = max(iy-1,0)
          iyp1 = min(iy+1,ny+1)
          iyp2 = min(iy+2,ny+1)
          do ix = i4, i8
            dtdym1 = (te(ix,iy)-te(ix,iym1))*gyf(ix,iym1)
            dtdy0 = (te(ix,iyp1)-te(ix,iy))*gyf(ix,iy)
            dtdyp1 = (te(ix,iyp2)-te(ix,iyp1))*gyf(ix,iyp1)
            d2tdy20 = (dtdy0 - dtdym1)*gy(ix,iy)
            d2tdy2p1 = (dtdyp1 - dtdy0)*gy(ix,iyp1)
            d3tdy3 = (d2tdy2p1 - d2tdy20)*gyf(ix,iy)
            feey4ord(ix,iy) = kye4order*d3tdy3*ney1(ix,iy)*
     .                                     sy(ix,iy)/gyf(ix,iy)**2
            feey(ix,iy) = feey(ix,iy) + feey4ord(ix,iy)
          enddo
          do ix = i4, i8
            dtdym1 = (ti(ix,iy)-ti(ix,iym1))*gyf(ix,iym1)
            dtdy0 = (ti(ix,iyp1)-ti(ix,iy))*gyf(ix,iy)
            dtdyp1 = (ti(ix,iyp2)-ti(ix,iyp1))*gyf(ix,iyp1)
            d2tdy20 = (dtdy0 - dtdym1)*gy(ix,iy)
            d2tdy2p1 = (dtdyp1 - dtdy0)*gy(ix,iyp1)
            d3tdy3 = (d2tdy2p1 - d2tdy20)*gyf(ix,iy)
            feiy4ord(ix,iy) = kyi4order*d3tdy3*niy1(ix,iy,1)*
     .                                     sy(ix,iy)/gyf(ix,iy)**2
            feiy(ix,iy) = feiy(ix,iy) + feiy4ord(ix,iy)
          enddo
        enddo
      endif

*  ---------------------------------------------------------------------
*  compute the energy residuals.
*  ---------------------------------------------------------------------

*  -- source terms --

      do 150 iy = j2, j5
         do 149 ix = i2, i5
            resee(ix,iy) = 
     .             seec(ix,iy) + seev(ix,iy) * te(ix,iy)
     .           + pwrsore(ix,iy)
     .           + cmneut * cmneutsor_ee * uesor_te(ix,iy)
     .           - nuvl(ix,iy,1)*vol(ix,iy)*bcee*ne(ix,iy)*te(ix,iy) 
            resei(ix,iy) = 
     .             seic(ix,iy) + seiv(ix,iy) * ti(ix,iy)
     .           + pwrsori(ix,iy)
     .           + cmneut * cmneutsor_ei * uesor_ti(ix,iy)
     .           - nuvl(ix,iy,1)*vol(ix,iy)*bcei*ne(ix,iy)*ti(ix,iy) 
  149    continue
  150 continue

*  -- divergence of electron and ion energy flows --

c...  Add y-component of nonorthogonal diffusive flux; convective component 
c...  already added to uu(ix,iy)
      if (isnonog .eq. 1) then

         do 3094 iy = j1, j6
            if (iy .gt. ny) goto 3094
            iy1 = max(iy-1,0)
            do 3093 ix = i1, i6
c...  First do the Te equation
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1)
                 grdnv=(    ( fym (ix,iy,1)*log(te(ix2,iy1 )) +  
     .                        fy0 (ix,iy,1)*log(te(ix2,iy  )) +
     .                        fyp (ix,iy,1)*log(te(ix2,iy+1)) +  
     .                        fymx(ix,iy,1)*log(te(ix ,iy1 )) +
     .                        fypx(ix,iy,1)*log(te(ix ,iy+1)) ) 
     .                     -( fym (ix,iy,0)*log(te(ix ,iy1 )) +
     .                        fy0 (ix,iy,0)*log(te(ix ,iy  )) +
     .                        fyp (ix,iy,0)*log(te(ix ,iy+1)) +
     .                        fymx(ix,iy,0)*log(te(ix4,iy1 )) +  
     .                        fypx(ix,iy,0)*log(te(ix6,iy+1)) ) ) / 
     .                                                   dxnog(ix,iy)  
               feexy(ix,iy) = exp( 0.5*
     .                         (log(te(ix2,iy)) + log(te(ix,iy))) )* 
     .                               (fcdif*kye+kye_use(ix,iy))*0.5*
     .                                       (ne(ix2,iy)+ne(ix,iy))*
     .                                     (grdnv/cos(angfx(ix,iy)) - 
     .                         (log(te(ix2,iy)) - log(te(ix,iy)))* 
     .                                         gxf(ix,iy))*sx(ix,iy)

c...  Now do the Ti equation.
c --- If we are using the parallel neutral momentum equation, we automatically
c --- change to a combined neutral+ion energy equation. We thus need to
c --- include the neutral heat conductivity. Since it is is isotropic
c --- we could use hcxn though we take the radial derivative; but this is
c --- only true if we dont flux limit.  Thus, we use 4-pt average of hcyn.
c --- Note: this four-point average results in not getting the full Jac. for
c --- a nonorthogonal mesh because of niy1,0 - see def. of hcyn

                 grdnv =(    ( fym (ix,iy,1)*log(ti(ix2,iy1 )) +  
     .                         fy0 (ix,iy,1)*log(ti(ix2,iy  )) +
     .                         fyp (ix,iy,1)*log(ti(ix2,iy+1)) +  
     .                         fymx(ix,iy,1)*log(ti(ix ,iy1 )) +
     .                         fypx(ix,iy,1)*log(ti(ix ,iy+1)) ) 
     .                      -( fym (ix,iy,0)*log(ti(ix ,iy1 )) +
     .                         fy0 (ix,iy,0)*log(ti(ix ,iy  )) +
     .                         fyp (ix,iy,0)*log(ti(ix ,iy+1)) +
     .                         fymx(ix,iy,0)*log(ti(ix4,iy1 )) +  
     .                         fypx(ix,iy,0)*log(ti(ix6,iy+1)) ) ) / 
     .                                                   dxnog(ix,iy)  
               feixy(ix,iy) = exp( 0.5*
     .                       (log(ti(ix2,iy)) + log(ti(ix,iy))) )*
     .                           ( (fcdif*kyi+kyi_use(ix,iy))*0.5*
     .                                     (nit(ix2,iy)+nit(ix,iy))
     .          + cftiexclg*cfneut*cfneutsor_ei*0.25*(hcyn(ix ,iy)+hcyn(ix ,iy1)
     .                              +hcyn(ix2,iy)+hcyn(ix4,iy1)) ) *
     .                                 (  grdnv/cos(angfx(ix,iy))
     .                         - (log(ti(ix2,iy)) - log(ti(ix,iy)))*
     .                                        gxf(ix,iy) )*sx(ix,iy)
c...  Flux limit with flalftxt even though hcys have parallel FL built in
               t0 = max(ti(ix,iy),temin*ev)
               t1 = max(ti(ix2,iy),temin*ev)
               vttn = t0*sqrt( t0/mi(1) )
               vttp = t1*sqrt( t1/mi(1) )
               qfl = flalftxy * (cftiexclg*0.125+(1.-cftiexclg)*0.25) * sx(ix,iy) * (vttn + vttp) * 
     .               (ni(ix,iy,1)+cftiexclg*ng(ix,iy,1)+ni(ix2,iy,1)+cftiexclg*ng(ix2,iy,1))
               feixy(ix,iy) = feixy(ix,iy) /
     .                              sqrt(1. + (feixy(ix,iy)/qfl)**2)

 3093      continue
 3094    continue

c...  Fix the fluxes with the same indice range as in fd2tra
         do iy = j4, j8
            do ix = i1, i5
               feex(ix,iy) = feex(ix,iy) - feexy(ix,iy)
               feix(ix,iy) = feix(ix,iy) - feixy(ix,iy)
            enddo
         enddo

      endif
c...  Finished with nonorthogonal mesh part

c ... Demand that net feex cannot be out of the plates
      if (isfeexpl0 == 1) then
        do iy = j4, j8
          do jx = 1, nxpt
            if(feex(ixlb(jx),iy) > 0. .and. ixmnbcl==1) then
              feexflr = ni(ixlb(jx),iy,1)*1.e4*ev*sx(ixlb(jx),iy)
              feex(ixlb(jx),iy) = feex(ixlb(jx),iy)/
     .                (1.+ (feex(ixlb(jx),iy)/feexflr)**4)
            endif
            if(feex(ixrb(jx),iy) < 0. .and. ixmxbcl==1) then
              feexflr = ni(ixrb(jx),iy,1)*1.e4*ev*sx(ixrb(jx),iy)
              feex(ixrb(jx),iy) = feex(ixrb(jx),iy)/
     .                (1.+ (feex(ixrb(jx),iy)/feexflr)**4)
            endif
          enddo
        enddo
      endif

      if (isfeixpl0 == 1) then
        do iy = j4, j8
          do jx = 1, nxpt
            if(feix(ixlb(jx),iy) > 0.) then
              feixflr = ni(ixlb(jx),iy,1)*1.e3*ev*sx(ixlb(jx),iy)
              feix(ixlb(jx),iy) = feix(ixlb(jx),iy)/
     .                (1.+ (feix(ixlb(jx),iy)/feixflr)**4)
            endif
            if(feix(ixrb(jx),iy) < 0.) then
              feixflr = ni(ixrb(jx),iy,1)*1.e3*ev*sx(ixrb(jx),iy)
              feix(ixrb(jx),iy) = feix(ixrb(jx),iy)/
     .                (1.+ (feix(ixrb(jx),iy)/feixflr)**4)
            endif
          enddo
        enddo
      endif

      do 310 iy = j2, j5
	    if((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
            feex(nxc-1,iy) = 0.
            feix(nxc-1,iy) = 0.
            feex(nxc,iy) = 0.
            feix(nxc,iy) = 0.
            feex(nxc+1,iy) = 0.
            feix(nxc+1,iy) = 0.
         endif
         if (islimon .ne. 0 .and. iy .ge. iy_lims) then
            feex(ix_lim,iy) = 0.
            feix(ix_lim,iy) = 0.
         endif
         if (nxpt==2 .and. ixmxbcl==1) then
            feex(ixrb(1)+1,iy) = 0.
            feix(ixrb(1)+1,iy) = 0.
         endif
         do 309 ix = i2, i5
            ix1 = ixm1(ix,iy)
            resee(ix,iy) = resee(ix,iy)
     .                  - ( feex(ix,iy) - feex(ix1,iy)
     .          + fluxfacy*(feey(ix,iy) - feey(ix,iy-1)) )
c ... ## IJ 2017 cfneutsor_ei flags above control neutral contrib.
            resei(ix,iy) = resei(ix,iy)					
     .                  - ( feix(ix,iy) - feix(ix1,iy)
     .          + fluxfacy*(feiy(ix,iy) - feiy(ix,iy-1)) )

c ... ## IJ 2016/10/19 add MC neutral flux
           if(get_neutral_moments .and. cmneutdiv_feg .ne. 0.0) then   
              jfld=1
              seg_ue(ix,iy,jfld)=-( (fegx_ue(ix,iy,jfld)-fegx_ue(ix1,iy,  jfld))
     .                   + fluxfacy*(fegy_ue(ix,iy,jfld)-fegy_ue(ix, iy-1,jfld)) )
     .                  *( (ni(ix,iy,jfld)*ti(ix,iy))/(ni(ix,iy,jfld)*ti(ix,iy)) )
              resei(ix,iy) = resei(ix,iy) +
     .                    cftiexclg*cmneutdiv*cmneutdiv_feg*seg_ue(ix,iy,jfld)
              reseg(ix,iy,1) = reseg(ix,iy,1) +
     .                             cmneutdiv*cmneutdiv_feg*seg_ue(ix,iy,jfld)
            endif
  309    continue
  310 continue


*  -- total energy residual and equipartition --

c...  Electron radiation loss -- ionization and recombination
            do 316 iy = iys1, iyf6  #iys,iyf
               do 315 ix = ixs1, ixf6  #iys, iyf
                  erlizold = erliz(ix,iy)
                  erlrcold = erlrc(ix,iy)
                  ne_sgvi = ne(ix,iy)
                  if (ifxnsgi.eq.1) ne_sgvi = cne_sgvi  # fix density dependence
                  if (istabon==16) then      # compute from b2frates data file
                     zmax=1
                     znuc=1
                     denz(0)=ng(ix,iy,1)     # use ngbackg as below ?
                     denz(1)=ni(ix,iy,1)     # use fac2sp  as below ?
                     dene=ne_sgvi
                     rdum=radmc(zmax,znuc,te(ix,iy),dene,denz,radz)
                     erliz(ix,iy)=chradi*radz(0)*vol(ix,iy)
                     if (isrecmon .ne. 0) erlrc(ix,iy)=chradr*radz(1)*vol(ix,iy)
                  else                       # compute from other data files
                     erliz(ix,iy) = chradi *
     .                           erl1(te(ix,iy),ne_sgvi,rtau(ix,iy))
     .                                  * (ng(ix,iy,1)-ngbackg(1)* 
     .                    (0.9+0.1*(ngbackg(1)/ng(ix,iy,1))**ingb) ) * 
     .                                                       vol(ix,iy)
                     if (isrecmon .ne. 0) erlrc(ix,iy) = chradr *
     .                               erl2(te(ix,iy),ne_sgvi,rtau(ix,iy))
     .                             * fac2sp*ni(ix,iy,1) * vol(ix,iy)
                  endif
                  eeliold = eeli(ix,iy)
                  if (icnuiz.le.1 .and. psor(ix,iy,1).ne.0.) 
     .                                           eeli(ix,iy) = 13.6*ev + 
     .                               erliz(ix,iy)/(fac2sp*psor(ix,iy,1))

                  pradhyd(ix,iy)= ( (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1)+
     .                                         erlrc(ix,iy) )/vol(ix,iy)
 315           continue
 316        continue

      do iy = iys1, iyf6  #j2, j5
        do ix = ixs1, ixf6  #i2, i5
          vsoreec(ix,iy) =
     .          - cfneut*cfneutsor_ee*cnsor*13.6*ev*fac2sp*psorc(ix,iy,1)
     .          + cfneut*cfneutsor_ee*cnsor*13.6*ev*fac2sp*psorrgc(ix,iy,1)
     .          - cfneut*cfneutsor_ee*cnsor*erliz(ix,iy)
     .          - cfneut*cfneutsor_ee*cnsor*erlrc(ix,iy)
     .          - cfneut*cfneutsor_ee*cnsor*ediss*ev*(0.5*psordis(ix,iy))
        enddo
      enddo

ccc         if (ishosor.eq.1) then  #full RHS eval
ccc
ccc           if (svrpkg.eq."cvode") then    # cannot access yl(neq+1)
ccc            call xerrab('*** svrpkg=cvode not allowed for ishosor=1 **')
ccc           endif 
ccc
ccc           if (yl(neq+1).lt.0) then  #full RHS eval
ccc
cccc ...    integ. source over cells (but not for Jac) for higher-order accuracy
ccc
ccc             call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
ccc     .                         fsprd, vol(0,0), psor_tmpov(0,0), vsoree)
ccc       
ccc           endif   # end of if (yl(neq+1).lt.0) test
ccc         endif    # end of integrating over sources and ishosor test

c*************************************************************
c   Perform 5pt average of source terms as volume integral
c*************************************************************
         if (iseesorave == 0.) then  #use only single-cell value
           do iy = iys1, iyf6
             do ix = ixs1, ixf6
               vsoree(ix,iy) = vsoreec(ix,iy)
             enddo
           enddo

         elseif (iseesorave > 0.) 

            if (xc < 0) then  #full RHS eval
              j2pwr = j2
              j5pwr = j5
            else  # Jacobian
              j2pwr = max(1, yc-1)
              j5pwr = min(ny, yc+1)
            endif 
            do iy = j2pwr, j5pwr
              if (xc < 0) then #full RHS eval
                i2pwr = i2
                i5pwr = i5
              else  #Jacobian eval
                i2pwr = max(1,ixm1(xc,yc))
                i5pwr = min(nx, ixp1(xc,yc))
              endif
              do ix = i2pwr, i5pwr
                ix1 = ixm1(ix,iy)
                ix2 = ixp1(ix,iy)
                vsoree(ix,iy) = (1.-iseesorave*0.5)*
     .                                  vsoreec(ix,iy)+ 
     .                               0.125*iseesorave*vol(ix,iy)*
     .                           ( vsoreec(ix,iy-1)/vol(ix,iy-1) + 
     .                             vsoreec(ix,iy+1)/vol(ix,iy+1) +
     .                             vsoreec(ix1,iy)/vol(ix1,iy)   + 
     .                             vsoreec(ix2,iy)/vol(ix2,iy) )
              enddo
            enddo

         endif    # end of integrating over sources and iseesorave test


      do 152 iy = j2, j5
         do 151 ix = i2, i5
            ix1 = ixm1(ix,iy)
            w0(ix,iy) = vol(ix,iy) * eqp(ix,iy) * (te(ix,iy)-ti(ix,iy))
            resee(ix,iy) = resee(ix,iy) - w0(ix,iy) + vsoree(ix,iy)
            if (isupgon(1).eq.1) then
c These terms include electron-ion equipartition as well as terms due
c to the friction force between neutrals and ions
               t1 = 0.5*(up(ix,iy,1)+up(ix1,iy,1))
               t2 = 0.5*(up(ix,iy,iigsp)+up(ix1,iy,iigsp))
               temp3 = cfnidhgy*0.25*(vy(ix,iy,iigsp)+vy(ix1,iy,iigsp))
     .                              *(vy(ix,iy,iigsp)+vy(ix1,iy,iigsp))
               temp4 = cfnidhg2*0.25*(v2(ix,iy,iigsp)+v2(ix1,iy,iigsp))
     .                              *(v2(ix,iy,iigsp)+v2(ix1,iy,iigsp))
               tv = cfticx*nucx(ix,iy,1)*ng(ix,iy,1)*vol(ix,iy)
               t0 = 1.5*( tg(ix,iy,1)* (psor(ix,iy,1)+tv)
     .                     -ti(ix,iy) * (psorrg(ix,iy,1)+tv) )
               resei(ix,iy) = resei(ix,iy) + w0(ix,iy)
     .             + cfneut * cfneutsor_ei * cfnidh * 0.5*mi(1) * 
     .                          ( (t1-t2)*(t1-t2)+temp3+temp4 ) * 
     .                    (  psor(ix,iy,1) + cftiexclg*psorrg(ix,iy,1)
     .              + tv + cftiexclg * tv  )
     .              + (1.0-cftiexclg) * t0
     .             + cftiexclg * cfneut * cfneutsor_ei * cnsor
     .               *( eion*ev+cfnidhdis*
     .                  0.5*mg(1)*(t2*t2+temp3+temp4) )*psordis(ix,iy) 
     .             + cfnidh2* 
     .                       ( -mi(1)*t1*t2*(psor(ix,iy,1)+tv)
     .                         +0.5*mi(1)*t1*t1*
     .                          (psor(ix,iy,1)+psorrg(ix,iy,1)+2*tv) )
               reseg(ix,iy,1) = reseg(ix,iy,1)
     .                            - t0+0.5*mg(1) * ( (t1-t2)*(t1-t2)
     .                                              +temp3+temp4 )
     .                            * (psorrg(ix,iy,1)+tv)
     .                            + ( eion*ev + cfnidh*cfnidhdis*
     .                   0.5*mg(1)*(t2*t2+temp3+temp4) )*psordis(ix,iy)
     .                     + cfnidh2* 
     .                       ( -mg(1)*t1*t2*(psorrg(ix,iy,1)+tv)
     .                         +0.5*mg(1)*(t2*t2+temp3+temp4)*
     .                          (psor(ix,iy,1)+psorrg(ix,iy,1)+2*tv) )
            else
               resei(ix,iy) = resei(ix,iy) + w0(ix,iy)
     .             + cfneut * cfneutsor_ei * ctsor*1.25e-1*mi(1)*
     .                    (upi(ix,iy,1)+upi(ix1,iy,1))**2*
     .                    fac2sp*psor(ix,iy,1)
     .             + cfneut * cfneutsor_ei * ceisor*cnsor* eion*ev*psordis(ix,iy)
     .             - cfneut * cfneutsor_ei * ccoldsor*ng(ix,iy,1)*nucx(ix,iy,1)*
     .                    (  1.5*ti(ix,iy)
     .                     - 0.125*mi(1)*(upi(ix,iy,1)+upi(ix1,iy,1))**2
     .                     - eion*ev  ) * vol(ix,iy)
            endif
  151    continue
  152 continue


c ... If molecules are present as gas species 2, add ion/atom cooling
      # energy transfer between ions and molecueles due to 
      # ion/molecule elastic collisions have been moved in
      # engbalg subroutine, so comment the following lines...
#      if(ishymol == 1) then
#        do iy = j2, j5
#          do ix = i2, i5
#            resei(ix,iy) = resei(ix,iy) - vol(ix,iy)*eqpg(ix,iy,2)*
#     .                                     (ti(ix,iy)-tg(ix,iy,2))
#          enddo
#        enddo
#      endif

*  -- Energy transfer to impurity neutrals at tg(,,igsp)
      if (ngsp >= 2) then   # for now, specialized to igsp=2 only
        do ifld = nhsp+1, nisp
          do iy = j2, j5    # iys,iyf limits dont seem to work(?)
            do ix = i2, i5
              resei(ix,iy) =resei(ix,iy) -cftiimpg*1.5*ni(ix,iy,ifld)*
     .                      (nucxi(ix,iy,ifld)+nueli(ix,iy,ifld))*
     .                      (ti(ix,iy) - tg(ix,iy,2))*vol(ix,iy)
            enddo
          enddo
        enddo
      endif

*  -- impurity radiation --

      if (isimpon .ge. 2) then
         if (istimingon .eq. 1) tsimp = gettime(sec4)

         do 533 iy = iys1, iyf6  #if Jacobian, only 1 cell done - local sor
            do 532 ix = ixs1, ixf6
               ntau(ix,iy) = atau(ix,iy) * ne(ix,iy)
               nratio(ix,iy) = ng(ix,iy,1)/ne(ix,iy)
               pradold = pwrzec(ix,iy)
               if (isimpon .eq. 2) then   # fixed-fraction model
                  na(ix,iy) = afrac(ix,iy) * ne(ix,iy)
                  pradcff(ix,iy) = na(ix,iy) * ne(ix,iy) *
     .               emissbs (te(ix,iy), nratio(ix,iy), ntau(ix,iy))
                  pradc(ix,iy) = pradcff(ix,iy)
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (isimpon .eq. 3) then   # average-ion model
                  na(ix,iy) = ni(ix,iy,nhsp+1)
                  pradc(ix,iy) = na(ix,iy) * ne(ix,iy) *
     .               radneq (te(ix,iy), nratio(ix,iy))
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (nzspt .eq. 0) then
                  # Multicharge group not allocated, so avoid radimpmc
                  pradc(ix,iy) = 0.
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (isimpon .ge. 4) then  # multi-charge model
                  pradc(ix,iy) = 0.
                  pwrzec(ix,iy) = 0.
                  nsm1 = nhsp 
                  do igsp = nhgsp+1, ngsp  # loop over diff imp species
                     jz = igsp - nhgsp
                     zn = znucl(nsm1+nzsp(jz))
                     if (ngsp .ge. nhgsp) nzloc(0) = ng(ix,iy,igsp)
                     do ifld = 1, nzsp(jz)
                        nzloc(ifld) = ni(ix,iy,nsm1+ifld)
                     enddo
                     nsm1 = nsm1 + nzsp(jz)   # setup for next igsp
                     argth = (te(ix,iy)-1.*ev)/(del_te_ro*ev)
                     fac_rad = 1.
                     if(del_te_ro.lt. 100.) fac_rad=0.5*(1+tanh(argth))
                     if (ismctab .eq. 1) then
                        pwrzec(ix,iy)= pwrzec(ix,iy) + fac_rad*
     .                                radimpmc (nzsp(jz), te(ix,iy), 
     .                                    ne(ix,iy), nzloc, impradloc)
                     elseif (ismctab .eq. 2) then
                        pwrzec(ix,iy)= pwrzec(ix,iy) + fac_rad*
     .                                radmc(nzsp(jz), zn, te(ix,iy), 
     .                                   ne(ix,iy), nzloc, impradloc)
                     endif
                     
                     do ifld = 0, nzsp(jz)
                        pradzc(ix,iy,ifld,jz) = impradloc(ifld)
                        pradc(ix,iy) = pradc(ix,iy)+impradloc(ifld)
                     enddo

                  enddo

		  if (isimpon .eq. 7) then  # add fixed-fraction contrib
                     na(ix,iy) = afrac(ix,iy) * ne(ix,iy)
                     pradcff(ix,iy) = na(ix,iy)* ne(ix,iy)*
     .                     emissbs(te(ix,iy), nratio(ix,iy), ntau(ix,iy))
                     pradc(ix,iy) = pradc(ix,iy) + pradcff(ix,iy)
                     pwrzec(ix,iy) = pwrzec(ix,iy) + pradcff(ix,iy)
                  endif

               endif
 532        continue
 533     continue

c*************************************************************
c   Perform 5pt average of source terms as volume integral
c*************************************************************
cc         if (ishosor.eq.0) then  #use only single-cell value
         if (iseesorave.eq.0.) then  #use only single-cell value
           do iy = iys1, iyf6
             do ix = ixs1, ixf6
               pwrze(ix,iy) = pwrzec(ix,iy)
               prad(ix,iy) = pradc(ix,iy)
               do igsp = nhgsp+1, ngsp
                 jz = igsp - nhgsp
                 do ifld = 0, nzsp(jz)
                   pradz(ix,iy,ifld,jz) = pradzc(ix,iy,ifld,jz)
                 enddo
               enddo
             enddo
           enddo

cc         elseif (ishosor .ne. 0) 
         elseif (iseesorave > 0.) 

            if (xc < 0) then  #full RHS eval
              j2pwr = j2
              j5pwr = j5
            else  # Jacobian
              j2pwr = max(1, yc-1)
              j5pwr = min(ny, yc+1)
            endif 
            do iy = j2pwr, j5pwr
              if (xc < 0) then #full RHS eval
                i2pwr = i2
                i5pwr = i5
              else  #Jacobian eval
                i2pwr = max(1,ixm1(xc,yc))
                i5pwr = min(nx, ixp1(xc,yc))
              endif
              do ix = i2pwr, i5pwr
                ix1 = ixm1(ix,iy)
                ix2 = ixp1(ix,iy)
                pwrze(ix,iy) = (1.-iseesorave*0.5)*pwrzec(ix,iy) + 
     .                                            0.125*iseesorave*
     .                         ( pwrzec(ix,iy-1)+ pwrzec(ix,iy+1)+
     .                           pwrzec(ix1,iy) + pwrzec(ix2,iy) )
                if (isimpon < 4) prad(ix,iy) = pwrze(ix,iy)
                if (isimpon >= 4) then  #prad, pradz only diagnostic
                  prad(ix,iy) = (1.-iseesorave*0.5)*pradc(ix,iy) + 
     .                                         0.125*iseesorave*
     .                          ( pradc(ix,iy-1)+ pradc(ix,iy+1)+
     .                            pradc(ix1,iy) + pradc(ix2,iy) )
                  do igsp = nhgsp+1, ngsp
                    jz = igsp - nhgsp
                    do ifld = 0, nzsp(jz)
                      pradz(ix,iy,ifld,jz) = (1.-iseesorave*0.5)*
     .                                        pradzc(ix,iy,ifld,jz) + 
     .                                            0.125*iseesorave*
     .              ( pradzc(ix,iy-1,ifld,jz)+ pradzc(ix,iy+1,ifld,jz)+
     .                pradzc(ix1,iy,ifld,jz) + pradzc(ix2,iy,ifld,jz) )
                    enddo
                  enddo
                endif
              enddo
            enddo

         endif    # end of integrating over sources and iseesorave test

c*******************************************************************
c ... Define a background elec energy source to prevent very low Te
c******************************************************************
      do iy = iys, iyf  #j2, j5
        do ix = ixs, ixf  #i2, i5
          pwrebkgold = pwrebkg(ix,iy)
          if (isimpon == 0) then
            pwrebkg(ix,iy) = (tebg*ev/te(ix,iy))**iteb*pwrbkg_c
          else  #add impurity rad loss
            pwrebkg(ix,iy) = (tebg*ev/te(ix,iy))**iteb*pwrbkg_c
          endif
        enddo
      enddo

c******************************************************************
c...  Update resee over whole "box" because initially set to zero 
c******************************************************************
         do 536 iy = j2, j5
            do 535 ix = i2, i5               
               resee(ix,iy) = resee(ix,iy) -
     .                            cnimp*pwrze(ix,iy)*vol(ix,iy) +
     .                                pwrebkg(ix,iy)*vol(ix,iy)
 535        continue
 536     continue

         if (istimingon .eq. 1) call timimpfj (tsimp, xc)
      endif  #loop for isimpon==2
  
*  -- joule heating --

      if (jhswitch > 0) then  # relies on div(J)=0, so omit iy=1 & ny
         if (isnewpot .eq. 1) then
            iy_min = 2
            iy_max = ny-1
         else
            iy_min = 1
            iy_max = ny
         endif
         if (jhswitch == 1) then   # div(J)=0 gives -grad(phi).J=-div(phi.J)
           do iy = max(iy_min, j2), min(iy_max, j5)
             do ix = i2, i5   
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               wjdote(ix,iy) = 
     .                          - 0.5*(fqp(ix,iy)+fq2(ix,iy))*
     .                                (phi(ix2,iy)+phi(ix,iy))
     .                          + 0.5*(fqp(ix1,iy)+fq2(ix1,iy))*
     .                                (phi(ix,iy)+phi(ix1,iy))
     .                          - 0.5*fqygp(ix,iy)*
     .                                (phi(ix,iy+1)+phi(ix,iy))
     .                          + 0.5*fqygp(ix,iy-1)*
     .                                (phi(ix,iy)+phi(ix,iy-1))
               resee(ix,iy) = resee(ix,iy) + wjdote(ix,iy) / ( 1. +
     .                             cfwjdotelim*(te(ix,iy)/tebg)**iteb )
             enddo
           enddo
         else  # for jhswitch > 1
           do iy = max(iy_min, j2), min(iy_max, j5)
             do ix = i2, i5    # use ex*fqx since phi(0,) may be large 
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               wjdote(ix,iy) = 
     .                       0.5*( ex(ix1,iy) *fqx(ix1,iy) +
     .                             ex(ix, iy) *fqx(ix, iy) )/gx(ix,iy)
     .                     + 0.5*( ey(ix, iy) *fqy(ix, iy) +
     .                             ey(ix,iy-1)*fqy(ix,iy-1) )/gy(ix,iy) 
               resee(ix,iy) = resee(ix,iy) + wjdote(ix,iy)
             enddo
           enddo
         endif
      endif

*  -- Now we introduce the viscous heating; one-side derviatives are used
*  -- on either side of the x-point where isxpty = 0

      do 157 iy = j2, j5
         do 156 ix = i2, i5
            do 155 ifld = 1, nusp  # if nusp --> nfsp, problems from y-term
               ix1 = ixm1(ix,iy)
               ix2 = ixm1(ix,iy+1)
               ix3 = ixm1(ix,iy-1)
	       thetacc = 0.5*(angfx(ix1,iy) + angfx(ix,iy))
	       dupdx = gx(ix,iy)*(upi(ix,iy,ifld)-upi(ix1,iy,ifld))
               wvh(ix,iy,ifld) = cfvcsx(ifld)*cfvisx*cos(thetacc)*
     .                                    visx(ix,iy,ifld)*dupdx**2
               if ( isxpty(ix,iy)==0 ) then  #1-sided deriv down in y
                 dupdy = 0.5*( upi(ix,iy,  ifld)+upi(ix1,iy  ,ifld) -
     .                         upi(ix,iy-1,ifld)-upi(ix3,iy-1,ifld) )*
     .                                                    gyf(ix,iy-1)
               elseif (isxpty(ix,iy)== -1) then #1-sided up in y
                 dupdy = 0.5*( upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld) -
     .                         upi(ix,iy  ,ifld)-upi(ix1,iy  ,ifld) )*
     .                                                    gyf(ix,iy)
               elseif (isxpty(ix,iy)==1.and.isvhyha==1) then
                                 #use harm y-ave for up face-values
                                 #take abs() to avoid near-zero denomin;
                                 #small err in wvh because up then small
                 upxavep1 = 0.5*(upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld))
                 upxave0 =  0.5*(upi(ix,iy  ,ifld)+upi(ix1,iy  ,ifld))
                 upxavem1 = 0.5*(upi(ix,iy-1,ifld)+upi(ix3,iy-1,ifld))
                 upf0  = 2.*upxavep1*upxave0*(upxavep1+upxave0) /
     .                           ( (upxavep1+upxave0)**2 + upvhflr**2 )
                 upfm1 = 2.*upxave0*upxavem1*(upxave0+upxavem1) /
     .                           ( (upxave0+upxavem1)**2 + upvhflr**2 )
                 dupdy = (upf0 - upfm1)*gy(ix,iy)
               else	#V7.08.04 option - linear ave in y-direction
		 dupdy = 0.25*( (upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld) - 
     .                           upi(ix,iy  ,ifld)-upi(ix1,iy  ,ifld))*
     .                                                     gyf(ix,iy) +
     .                          (upi(ix,iy  ,ifld)+upi(ix1,iy  ,ifld) -
     .                           upi(ix,iy-1,ifld)-upi(ix3,iy-1,ifld))*
     .                                                     gyf(ix,iy-1) )
               endif
               wvh(ix,iy,ifld) = wvh(ix,iy,ifld) + cfvcsy(ifld)*cfvisy*
     .                                   visy(ix,iy,ifld)*dupdy**2
	       wvh(ix,iy,ifld) = wvh(ix,iy,ifld) -
     .                             sin(thetacc)*cfvcsy(ifld)*cfvisy*
     .                                   visy(ix,iy,ifld)*dupdx*dupdy
            if (zi(ifld)==0.0 .and. ifld.eq.iigsp) then 
              resei(ix,iy) = resei(ix,iy) + cftiexclg*wvh(ix,iy,ifld)*vol(ix,iy)
              reseg(ix,iy,1) = reseg(ix,iy,1) + wvh(ix,iy,ifld)*vol(ix,iy)
            else
              resei(ix,iy) = resei(ix,iy) + wvh(ix,iy,ifld)*vol(ix,iy)
            endif
  155       continue   # loop over up species ifld
  156    continue
 157  continue


c*******************************************************************
c ... Define a background ion energy source to prevent very low Ti
c******************************************************************
      do iy = iys, iyf  #j2, j5
        do ix = ixs, ixf  #i2, i5
          pwribkgold = pwribkg(ix,iy)
          pwribkg(ix,iy) = (tibg*ev/ti(ix,iy))**iteb*pwribkg_c
        enddo
      enddo
 
      do iy = j2, j5
        do ix = i2, i5
          resei(ix,iy) = resei(ix,iy) + pwribkg(ix,iy)*vol(ix,iy)
        enddo
      enddo


**********************************************************************
*  --  Equations to be solved --
**********************************************************************
      do 270 iy = j2, j5
         do 260 ix = i2, i5
            do 254 ifld = 1, nisp
	       if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  yldot(iv) = (1-iseqalg(iv)) *
     .                        resco(ix,iy,ifld)/(vol(ix,iy)*n0(ifld))
               endif
 254        continue
            do 255 ifld = 1, nusp
	       if(isuponxy(ix,iy,ifld) .eq. 1) then
                  iv = idxu(ix,iy,ifld)
                  yldot(iv) = (1-iseqalg(iv)) *
     .                        resmo(ix,iy,ifld)/(volv(ix,iy)*fnorm(ifld))
                  do jx = 1, nxpt
                     if (ix.eq.ixrb(jx) .and. ixmxbcl.eq.1) yldot(iv) =
     .                        resmo(ix,iy,ifld)/(volv(ix,iy)*fnorm(ifld))
                  enddo
               endif
 255        continue
            if(isteonxy(ix,iy) == 1) then
              iv =  idxte(ix,iy)
	      yldot(iv) = (1-iseqalg(iv)) *
     .                                 resee(ix,iy)/(vol(ix,iy)*ennorm)
            endif
            if(istionxy(ix,iy) == 1) then
              iv1 = idxti(ix,iy)
	      yldot(iv1) = (1-iseqalg(iv1)) *
     .                                 resei(ix,iy)/(vol(ix,iy)*ennorm)
            endif
            do 256 igsp = 1, ngsp
	      if(isngonxy(ix,iy,igsp).eq.1) then
                iv2 = idxg(ix,iy,igsp)
                yldot(iv2) = (1-iseqalg(iv2)) *
     .                      resng(ix,iy,igsp)/(vol(ix,iy)*n0g(igsp))
              endif
	      if(istgonxy(ix,iy,igsp).eq.1) then
                iv2 = idxtg(ix,iy,igsp)
                yldot(iv2) = (1-iseqalg(iv2)) *
     .                      reseg(ix,iy,igsp)/(vol(ix,iy)*ennorm)
              endif
 256        continue
 260     continue
 270  continue
c ... The factor (1-iseqalg(iv)) above forces yldot=0 for algebraic
c ... equations, except up(nx,,); these yldot are subsequently set in
c ... subroutine bouncon.


c  POTEN calculates the electrostatic potential, and BOUNCON calculates the 
c  equations for the boundaries. For the vodpk solver, the B.C. are ODEs 
c  in time (rate equations).  Both bouncon and poten must be called before
c  the perturbed variables are reset below to get Jacobian correct

      if (isphion.eq.1) call poteneq (neq, yl, yldot)

      call bouncon (neq, yl, yldot)

c...  Finally, reset some source terms if this is a Jacobian evaluation
         if (xc .ge. 0 .and. yc .ge. 0) then
            ix1 = ixm1(xc,yc)
            if(isimpon.gt.0) pwrzec(xc,yc) = pradold
            pwrebkg(xc,yc) = pwrebkgold
            pwribkg(xc,yc) = pwribkgold
            erliz(xc,yc) = erlizold
            erlrc(xc,yc) = erlrcold
            eeli(xc,yc) = eeliold
            fqp(ix1,yc) = fqpom
            fqp(xc,yc) = fqpo
            frice(ix1,yc) = friceom
            frice(xc,yc) = friceo
            upe(ix1,yc) = upeom
            upe(xc,yc) = upeo
            psordis(xc,yc) = psordisold
            do ifld = 1, nfsp
               psorc(xc,yc,ifld) = psorold(ifld)
               psorxr(xc,yc,ifld) = psorxrold(ifld)
               frici(ix1,yc,ifld) = friciom(ifld)
               frici(xc,yc,ifld) = fricio(ifld)
               upi(ix1,yc,ifld) = upiom(ifld)
               upi(xc,yc,ifld) = upio(ifld)
               uup(ix1,yc,ifld) = uupom(ifld)
               uup(xc,yc,ifld) = uupo(ifld)
               nucxi(xc,yc,ifld) = nucxiold(ifld)
               nueli(xc,yc,ifld) = nueliold(ifld)
            enddo
            do igsp = 1, ngsp
               nucx(xc,yc,igsp) = nucxold(igsp)
               nurc(xc,yc,igsp) = nurcold(igsp)
               nuiz(xc,yc,igsp) = nuizold(igsp)
               nuelg(xc,yc,igsp) = nuelgold(igsp)
               nuix(xc,yc,igsp) = nuixold(igsp)
               psorgc(xc,yc,igsp) = psorgold(igsp)
               psorrgc(xc,yc,igsp) = psorrgold(igsp)
               psorcxgc(xc,yc,igsp) = psorcxgold(igsp)
            enddo
         endif

      if (ismcnon .eq. 4) then # test a different fluid model in the preconditioner
         if (yl(neq+1) .gt. 0) then   # Precon eval
            parvis=parvis/pnc_cfparvis
            travis=travis/pnc_cftravis
            do ifld=1,nisp
              ni(:,:,ifld)=ni(:,:,ifld)/pnc_cfni(ifld)
              up(:,:,ifld)=up(:,:,ifld)/pnc_cfup(ifld)
            enddo
         endif
      end if #ismcnon

c ... Accumulate cpu time spent here.
      if(xc .lt. 0) then
         ttotfe = ttotfe + gettime(sec4) - tsfe
      else
         ttotjf = ttotjf + gettime(sec4) - tsjf
      endif
      if (TimingPandfOn.gt.0) TotTimePandf=TotTimePandf+tock(TimePandf)
      return
      end
c****** end of subroutine pandf ************
c-----------------------------------------------------------------------
      subroutine mombal0 (nisp, nhsp, nzsp, minu, ziin,
     .                                      misotope, natomic, nchstate)
c ... Compute 'misotope', 'nchstate', and 'natomic', and allocate memory
c     for arrays used in subroutine mombal.
 
      implicit none

c ... Input arguments:
      integer nisp      # total number of ion species
      integer nhsp      # total number of hydrogenic ion species
      integer nzsp(ngspmx-1)   # number of charge states for each imp isotope
      real minu(nisp)   # mass (in amu) of ion species
      real ziin(nisp)   # charge (in units of e) of ion species

c ... Output arguments:
      integer misotope  # total number of isotopes (including electrons)
      integer natomic(*)   # maximum charge state of each isotope
      integer nchstate     # maximum charge state among all isotopes

c ... Local variables:
      integer misa, ifld, jz
	  
c ... Loop over ion species, looking for change to a new isotope, and
c     finding maximum charge state.
      natomic(1) = 1   # electrons are "isotope 1"
      nchstate = 0
      misa = 2
      do ifld = 1, nhsp
         natomic(misa) = max(nint(ziin(ifld)), 1)   # must be .ge. 1
         nchstate = max(nchstate, natomic(misa)) 
         if (ifld .eq. nhsp) go to 50
         if (minu(ifld+1) .ne. minu(ifld)) misa = misa + 1
      enddo
 50   misotope = misa
      do jz = 1, ngspmx-1
         if (nzsp(jz)==0) break
         misotope = misotope + 1
	 if (misotope .gt. MXMISO) then
           call remark("subroutine mombal0 error: ")
           call remark("To avoid write out-of-bounds for array natomic")
           call remark("increase the value of MXMISO and recompile.")
	       call xerrab("")
         endif
         natomic(misotope) = nzsp(jz)
         nchstate = max(nchstate, natomic(misotope)) 
      enddo

c ... Allocate memory for arrays used in subroutine mombal.
      call gallot("Reduced_ion_interface", 0)

      return
      end
c****** end of subroutine mombal0 ************
c-----------------------------------------------------------------------
      subroutine mombal (ix,ix1,iy)
c ... Prepare information needed to call Steve Hirshman reduced-ion
c     momentum-balance routine, and distribute results from it.
c     Note that results are computed at a poloidal density-cell face.
c     The parallel flow velocities are returned in arrays up and upe.

      implicit none

c ... Input arguments:
      integer ix    # poloidal index of density cell to left of face
      integer ix1   # poloidal index of density cell to right of face
      integer iy    # radial index of density cell

c ... Common blocks:
      Use(Dim)          # nx,ny,nisp
      Use(Selec)        # isupgoon
      Use(Comgeo)       # sx,rrv
      Use(UEpar)        # lnlam,isupgon
      Use(Cfric)        # frice,frici
      Use(Compla)       # ni,up,upe,upi,te,ti,ng,ne
      Use(Gradients)    # gpix,gtix,ex,gtex,gpex
      Use(Comflo)       # fqp
      Use(UEint)        # minu,ziin
      Use(Reduced_ion_interface) # misotope,nchstate,natomic,
                                 # amu,tempa,qneut,uneut,den,gradp,
                                 # gradt,friction,nuion,nurec,qcond,ucond
                                 # friccomp
      Use(Imprad)       # ismctab
      Use(Phyvar)       # ev
      Use(Comtra)       # fricflf,sigvi_floor
      Use(Share)        # cutlo

c ... External functions:
      real rra, rsa

c ... Local variables:
      integer ifld, misa, nz
      integer ldir
      real rdum, dloglam, epar, parcurrent, umass, lmfpe, lmfpi,
     .     ltmax, tif, flxlimf, umassni, massni

c ... Determine a flux-limit factor for all input terms by finding the min
c ... scale length. First consider the electrons
      den(1,1) = 0.5 * (ne(ix,iy) + ne(ix1,iy))
      tempa(1) = 0.5 * (te(ix,iy) + te(ix1,iy))
      tif = 0.5 * (ti(ix,iy) + ti(ix1,iy))
      ltmax = min( abs(tempa(1)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .             abs(tif/(rrv(ix,iy)*gtix(ix,iy) + cutlo)),
     .             abs(den(1,1)*tempa(1)/
     .                          (rrv(ix,iy)*gpex(ix,iy) + cutlo)) )
      lmfpe = 1.e16*((tempa(1)/ev)**2/den(1,1)) #Approx Coulomb e-mean-free-path
      lmfpi = 1.e16*(tif/ev)**2/den(1,1)  # Approx. Coulomb i-mean-free-path
c ... Now check ion pressure gradient scale-lengths
      do nz = 1, 1   # previously nisp - only do majority species now
         if (zi(nz) .gt. 1.e-10) then  # omits neutrals species
            ltmax = min( ltmax, abs( 0.5*(pri(ix,iy,nz)+pri(ix1,iy,nz))/
     .                           (rrv(ix,iy)*gpix(ix,iy,nz) + cutlo) ) )
         endif
      enddo
      flxlimf = 1 / (1 + fricflf*((lmfpe+lmfpi)/ltmax)**2)

c ... Store electron density and temperature, gradients of electron
c     pressure and temperature, (hydrogenic) gas density, ionization rate
c     of gas, and a few other (physically meaningless) quantities for which
c     there are array locations.
      amu(1) = 5.45e-4
      gradp(1,1) = flxlimf*rrv(ix,iy) * gpex(ix,iy)
      gradt(1,1) = flxlimf*rrv(ix,iy) * gtex(ix,iy)
      den(2,0) = 0.5 * (ng(ix,iy,1) + ng(ix1,iy,1))
      nuion(2,0) = den(1,1) * (rsa(tempa(1), den(1,1), 0., 0)
     .                         + sigvi_floor)
      den(1,0) = 0.
      qneut(1) = 0.
      uneut(1) = up(ix,iy,1)  # netural velocity. Cant be zero?; was umass
      nuion(1,0) = 0.
      nuion(1,1) = 0.
      nurec(1,1) = 0.

c ... Loop over isotopes.
      ifld = 1
      do misa = 2, misotope
         amu(misa) = minu(ifld)   # Store mass of this isotope

c ... Store ionization rate of neutral if this isotope is an impurity.
         if (natomic(misa) .gt. 1) then
            den(misa,0) = 0.              # impurity neutral density
            if (ismctab .eq. 1) then
               call imprates(tempa(1), 0, natomic(misa),
     .                       nuion(misa,0), rdum, rdum) 
               nuion(misa,0) = den(1,1)*nuion(misa,0) #sigv-->ne*sigv
            elseif (ismctab .eq. 2) then
               call mcrates(den(1,1), tempa(1), 0., 0, natomic(misa),
     .                      znucl(ifld), nuion(misa,0), rdum, rdum)
               nuion(misa,0) = den(1,1)*nuion(misa,0) #sigv-->ne*sigv
            endif
         endif

c ...    Loop over charged states, storing ni and parallel gradients
c        of pressure and Ti at poloidal density-cell faces.
         do nz = 1, natomic(misa)
            den(misa,nz) = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            gradp(misa,nz) = flxlimf*rrv(ix,iy) * gpix(ix,iy,ifld)
            tempa(misa) = 0.5 * (ti(ix,iy) + ti(ix1,iy))
            gradt(misa,nz) = flxlimf*rrv(ix,iy) * gtix(ix,iy) 
c .......   Get ionization and recombination rates.
c           Note that nuion has no meaning for the fully-stripped state,
c           but space is available to store zero returned by imprates.
            if (isupimpap .eq. 0) then  # omit atomic physics coupling in up
               nuion(misa,nz) = 0.
               nurec(misa,nz) = 0.
            elseif (natomic(misa) .eq. 1) then   # hydrogenic isotope
               nuion(misa,1) = 0.
               nurec(misa,1) = den(1,1) *
     .                         rra(tempa(1), den(1,1), 0., 1)
            else                             # impurity isotope
               if (ismctab .eq. 1) then
                  call imprates(tempa(1), nz, natomic(misa),
     .                          nuion(misa,nz), nurec(misa,nz), rdum)
                  nuion(misa,nz) = den(1,1)*nuion(misa,nz) #sigv-->ne*sigv
                  nurec(misa,nz) = den(1,1)*nurec(misa,nz)
               elseif (ismctab .eq. 2) then
                  call mcrates(den(1,1), tempa(1), 0., nz, natomic(misa),
     .                znucl(ifld), nuion(misa,nz), nurec(misa,nz), rdum)
                  nuion(misa,nz) = den(1,1)*nuion(misa,nz) #sigv-->ne*sigv
                  nurec(misa,nz) = den(1,1)*nurec(misa,nz)
               endif
            endif
            ifld = ifld + 1
            if (natomic(misa).eq.1.and.isupgon(ifld-1).eq.1) ifld=ifld+1
         enddo    # end of loop over charge states with index nz
      enddo    # end of loop over isotopes with index misa

c ... Set up other inputs for fmombal, including flow velocities.
 50   epar = flxlimf*rrv(ix,iy) * ex(ix,iy)
      umassni = 0.
      massni = 0.
      do ifld = 1, nusp
        if(zi(ifld) .gt. 1.e-20) then
          umassni = umassni + mi(ifld)*ni(ix,iy,ifld)*up(ix,iy,ifld)
          massni = massni + mi(ifld)*ni(ix,iy,ifld)
        endif
      enddo
      umass = umassni/massni
      do misa = 2, misotope
         qneut(misa) = 0.
         uneut(misa) = up(ix,iy,1)  # use hydr ion as default (orig)
      enddo
      ldir = 2
      dloglam = loglambda(ix,iy)

c ... Call Steve Hirshman reduced-ion momentum-balance routine.
      call fmombal(amu,den,dloglam,epar,friction,gradp,gradt,
     >	 nuion,nurec,qcond,qneut,ucond,uneut,umass,
     >   parcurrent,tempa,natomic,misotope,nchstate,ldir,friccomp)

c ... Distribute results into arrays used in pandf.  Note that we do
c     nothing with qcond.
      fqp(ix,iy) = cfparcur*parcurrent * rrv(ix,iy)*sx(ix,iy)
      frice(ix,iy) = friction(1,1)
      upe(ix,iy) = ucond(1,1)
      ifld = 0
      do misa = 2, misotope
         do nz = 1, natomic(misa)
 60         ifld = ifld + 1
            if (ziin(ifld) .lt. 1.e-10) goto 60    #omits neutral species
            frici(ix,iy,ifld) = friction(misa,nz)
            upi(ix,iy,ifld) = ucond(misa,nz)
            upifmb(ix,iy,ifld) = ucond(misa,nz)    #diagnostic only
         enddo
      enddo

      return
      end
c****** end of subroutine mombal ************
c-----------------------------------------------------------------------
      subroutine mombalni (ix,ix1,iy)
c ... Use force balance for the impurity momentum equation, neglecting
c     inertia, viscosity, and atomic-physics coupling ala Knoll, Campbell.
c     We also use Keilhacker, et al., Nucl. Fusion., Vol. 31, 537 (1991)
c     which differs somewhat from Campbell: alfi gets divided by zeffv
c     but we retain the extra term 0.6*... in betai from Campbell.
c     Note that results are computed at a poloidal density-cell face.
c     The parallel flow velocities are returned in arrays up and upe.

      implicit none

c ... Input arguments:
      integer ix    # poloidal index of density cell to left of face
      integer ix1   # poloidal index of density cell to right of face
      integer iy    # radial index of density cell

c ... Common blocks:
      Use(Dim)          # nx,ny,nisp,nhsp,nusp
      Use(Selec)        # isupgon
      Use(Comgeo)       # sx,rrv,vol
      Use(UEpar)        # lnlam,isupgon,isofric,is_z0_imp_const,z0_imp_const
      Use(Cfric)        # frice,frici,cfgte,cfgti,cftaud,alfe,betai
      Use(Compla)       # ni,up,upe,upi,te,ti,ng,ne,zeff,netap
      Use(Gradients)    # gpix,gtix,ex,gtex,gpex
      Use(Comflo)       # fqp
      Use(UEint)        # minu,ziin
      Use(Reduced_ion_interface) # misotope,nchstate,natomic,dztot,
                                 # amu,tempa,qneut,uneut,den,gradp,
                                 # gradt,friction,nuion,nurec,qcond,ucond
      Use(Imprad)       # ismctab
      Use(Phyvar)       # ev,qe
      Use(Comtra)       # fricflf,fupe_cur
      Use(Share)        # cutlo
      Use(Coefeq)       # cfnetap
      Use(Volsrc)       # volmsor
      Use(Npes_mpi)
      Use(Conduc)       # pondomfpari_use


c ... Local variables:
      integer ifld, misa, nz
      integer ldir
      real rdum, tdum, dloglam, epar, parcurrent, umass, lmfpe, lmfpi,
     .     zeffv, z0, taud, taudeff, ltmax, tif, flxlimf,
     .     dzj, dzz2tot


c ... Determine a flux-limit factor for all input terms by finding the min
c ... scale length. First consider the electrons
      den(1,1) = 0.5 * (ne(ix,iy) + ne(ix1,iy))
      tempa(1) = 0.5 * (te(ix,iy) + te(ix1,iy))
      tif = 0.5 * (ti(ix,iy) + ti(ix1,iy))
      ltmax = min( abs(tempa(1)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .             abs(tif/(rrv(ix,iy)*gtix(ix,iy) + cutlo)),
     .             abs(den(1,1)*tempa(1)/
     .                          (rrv(ix,iy)*gpex(ix,iy) + cutlo)) )
      lmfpe = 1.e16*((tempa(1)/ev)**2/den(1,1)) #Approx Coulomb e-mean-free-path
      lmfpi = 1.e16*((tif/ev)**2/den(1,1))  # Approx. Coulomb i-mean-free-path
c ... Now check ion pressure gradient scale-lengths
      do nz = 1, 1      # previously nisp - only do majority species now
         if (zi(nz) .gt. 1.e-10) then  # omits neutrals species
            ltmax = min( ltmax, abs( 0.5*(pri(ix,iy,nz)+pri(ix1,iy,nz))/
     .                           (rrv(ix,iy)*gpix(ix,iy,nz) + cutlo) ) )
         endif
      enddo
      flxlimf = 1 / (1 + fricflf*((lmfpe+lmfpi)/ltmax)**2)

c ... Store electron density and temperature, gradients of electron
c     pressure and temperature, (hydrogenic) gas density, ionization rate
c     of gas, and a few other (physically meaningless) quantities for which
c     there are array locations.
      gradp(1,1) = rrv(ix,iy) * gpex(ix,iy)
      gradt(1,1) = rrv(ix,iy) * gtex(ix,iy)
      frice(ix,iy) = -0.71*flxlimf*den(1,1)*gradt(1,1) +
     .                cfnetap*qe*netap(ix,iy)*fqp(ix,iy)/sx(ix,iy)

c ... Loop over charge states to get total impurity density
      ifld = nhsp
      dztot = 0.
      dzz2tot = 0.
      do misa = 3, misotope
         do nz = 1, natomic(misa)
            ifld = ifld + 1
            dzj = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            dztot = dztot + dzj
            dzz2tot = dzz2tot + dzj*zi(ifld)**2
         enddo
      enddo

c ... Set the hydrogen values based on electrons (trace limit)
ccc      frici(ix,iy,1) = - frice(ix,iy)   # needed for hydrogen
ccc   For arbitrary impurity concentration use the following frici:
      frici(ix,iy,1) = - frice(ix,iy) *
     .                   ni(ix,iy,1)*zi(1)**2/(ne(ix,iy)*zeff(ix,iy))
      upe(ix,iy) = up(ix,iy,1) - fupe_cur*fqp(ix,iy)/( sx(ix,iy)*qe*
     .                         rrv(ix,iy)*0.5*(ne(ix,iy)+ne(ix1,iy)) )
      upi(ix,iy,1) = up(ix,iy,1)
      den(2,1) = 0.5 * (ni(ix,iy,1) + ni(ix1,iy,1))
           
c ... Loop over isotopes for friction coefficients
      ifld = 1
      do misa = 3, misotope      # only executed if impurities are present

c ... Loop over charged states, storing ni and parallel gradients
c ... of pressure and Ti at poloidal density-cell faces.
         do nz = 1, natomic(misa)  # note: hydrogen friction set above
            ifld = ifld + 1
            if (ziin(ifld) .lt. 1.e-10) ifld = ifld+1 #skip gas index
            den(misa,nz) = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            gradp(misa,nz) = rrv(ix,iy) * gpix(ix,iy,ifld) -
     .                                      pondomfpari_use(ix,iy,ifld)
            tempa(misa) = 0.5 * (ti(ix,iy) + ti(ix1,iy))
            gradt(misa,nz) = rrv(ix,iy) * gtix(ix,iy) 
            zeffv = 0.5*(zeff(ix,iy)+zeff(ix1,iy))
            if (is_z0_imp_const == 0) then
              z0 = den(1,1)*zeffv/den(2,1) - 1.
            else # fixed user input
              z0 = z0_imp_const
            endif
            if (isbetaicalc(ifld) == 1) then
              betai(ifld)=cfgti*1.56*zi(ifld)**2*(1+1.414*z0)*(1+.52*z0)/
     .                       ( (1+2.65*z0)*(1+.285*z0)*( z0 + 
     .                         sqrt( 0.5*(mi(1)+mi(ifld))/mi(ifld)) ) )
     .                       + 0.6*(zi(ifld)**2*den(misa,nz)/dzz2tot - 1.)
            endif
            if (isalfecalc(ifld) == 1) then
              alfe(ifld) = cfgte*2.2*zi(ifld)**2*(1+.52*zeffv) /
     .                      ( (1+2.65*zeffv)*(1+.285*zeffv)*zeffv )
            endif
c... NOTE:next coefficient 12*pi*sqrt(pi/2)*epsilon**2/e**4 = 5.624e54 in mks
            taud =cftaud*5.624e54*sqrt(mi(1))*mi(ifld)*tempa(misa)**1.5 /
     .             ( loglambda(ix,iy)*den(misa,nz)*zi(ifld)**2*
     .                                                 (mi(1)+mi(ifld)) )
            taudeff = flxlimf*taud*den(misa,nz)*(1+2.65*z0)*(1+.285*z0) /
     .                         ( den(1,1)*(1+.24*z0)*(1+.93*z0) )
            upi(ix,iy,ifld) = up(ix,iy,1) + (taudeff/mi(1)) * (
     .                         - gradp(misa,nz)/den(misa,nz)
     .                         + alfe(ifld)*gradt(1,1) 
     .                         + betai(ifld)*gradt(misa,nz) 
     .                         + qe*zi(ifld)*rrv(ix,iy)*ex(ix,iy) 
     .                         + volmsor(ix,iy,ifld)/
     .                                       (den(misa,nz)*vol(ix,iy)) )
c ...       For force balance, frici just balances E-field and pressure
c ...       No flxlimf for 1st option; it only enhances (1/taudeff)
            if (nusp-isupgon(1) .eq. 1) then  #only frici(,,1) used here
              frici(ix,iy,ifld) =-qe*zi(ifld)*den(misa,nz)*
     .                           rrv(ix,iy)*ex(ix,iy) + gradp(misa,nz)  
            else # multi ion mom eqns; drag calc elsewhere (w0) if isofric=1
              frici(ix,iy,ifld) =flxlimf*den(misa,nz)*( 
     .                                          alfe(ifld)*gradt(1,1) + 
     .                                     betai(ifld)*gradt(misa,nz) + 
     .                                               (1-isofric)*mi(1)*
     .                           (up(ix,iy,1)-up(ix,iy,ifld))/ taudeff )
            endif
         enddo
      enddo

      return
      end
c****** end of subroutine mombalni ************
c-----------------------------------------------------------------------
      subroutine timimpfj (tsimp, xc)
      real(Size4) tsimp
      integer xc
      Use(Timing)   # ttimpfe,ttimpjf
      real(Size4) sec4, gettime, dtimp

      dtimp = gettime(sec4) - tsimp
      if (xc .lt. 0) then
         ttimpfe = ttimpfe + dtimp
      else
         ttimpjf = ttimpjf + dtimp
      endif

      return
      end
c---- end of subroutine timimpfj ---------------------------------------
c-----------------------------------------------------------------------
      subroutine turbdif (ix, iy, ixmp3, iyp1, ifld)
c ... For a grid cell outside the separatrix, calculate anomalous
c     diffusivity due to turbulence and return it via common array
c     diffusivwrk.  In addition, this subroutine computes values for
c     arrays chinorml and chinormh.
      implicit none

c ... Input arguments:
      integer ix, iy, ixmp3, iyp1
      integer ifld

c ... Common blocks:
      Use(Dim)              # nx,ny,nisp
      Use(Xpoint_indices)   # iysptrx
      Use(Comgeo)           # gyf,linelen
      Use(Compla)           # ney0,ney1,nity0,nity1,
                            # tiy0,tiy1,tey0,tey1,priy0,priy1,mi
      Use(Gradients)        # gtey,gpiy
      Use(Bfield)           # btot
      Use(Comtra)           # diffusivwrk
      Use(Turbulence)       # lambdan,lambdat,isturbnloc
      Use(Turbulence_diagnostics)   # chinorml,chinormh

c ... Local variables:
      integer ix0
      real drdr0, glte, lte, glpi, lpi, teyf, tiyf
      real neyf, nityf, priyf, btotyf, zavg
      real ted, tid, densd

      if (iy .gt. iysptrx) then
c...  For a full double-null configuration, iysptrx refers to the last
c...  closed flux surface (see definition in subroutine nphygeo).

c ... Select local or midplane value of ix.
         if (isturbnloc .ne. 1) then
           ix0 = ix
         else
           ix0 = ixmp3
         endif

c ... Compute mean Z at the center of the y-face.
         neyf = 0.5 * (ney1(ix0,iy) + ney0(ix0,iy))
         nityf = 0.5 * (nity1(ix0,iy) + nity0(ix0,iy))
         zavg = neyf / nityf

c ... Compute radial scale lengths of Te and ion pressure based on
c     y-face values and midplane grid spacing.
         drdr0 = gyf(ixmp3,iy) / gyf(ix0,iy)
         tiyf = 0.5 * (tiy1(ix0,iy) + tiy0(ix0,iy))
         teyf = 0.5 * (tey1(ix0,iy) + tey0(ix0,iy))
         glte = abs(gtey(ix0,iy)) * drdr0 / teyf
         glte = max(glte, 1.)
         lte = 1. / glte
         priyf = 0.5 * (priy1(ix0,iy,ifld) + priy0(ix0,iy,ifld))
         glpi = abs(gpiy(ix0,iy,ifld)) * drdr0 / priyf
         glpi = max(glpi, 1.)
         lpi = 1. / glpi

c ... Reduce local temperatures to values approximating those at divertor
c     plates.  Similarly increase the local density to approximate plate
c     value.  (This procedure is preliminary to a more rigorous one of
c     using variables at plates.)
c ... Compute values at divertor plates in either of two ways:
c     for local calculation of anomalous diffusivity, use input variables
c     lambdan and lambdat to get approximations to plate values ...
         if (isturbnloc .ne. 1) then
           ted = teyf / lambdat
           tid = tiyf / lambdat
           densd = lambdan * neyf
c     or for nonlocal calculation of diffusivity, use average of values
c     at the two plates.
         else
           ted = 0.25 * (tey1(0   ,iy) + tey0(0   ,iy) +
     .                   tey1(nx+1,iy) + tey0(nx+1,iy))
           tid = 0.25 * (tiy1(0   ,iy) + tiy0(0   ,iy) +
     .                   tiy1(nx+1,iy) + tiy0(nx+1,iy))
           densd = 0.25 * (ney1(0   ,iy) + ney0(0   ,iy) +
     .                     ney1(nx+1,iy) + ney0(nx+1,iy))
         endif

c ... Compute anomalous diffusivity.
         btotyf = 0.5 * (btot(ix0,iy) + btot(ix0,iyp1))
         call turb_diffus (btotyf, lte, lpi, teyf, tiyf, neyf,
     .                     ted, tid, densd, mi(ifld), zavg, linelen,
     .                     diffusivwrk(ix,iy),
     .                     chinorml(ix,iy), chinormh(ix,iy))
      endif

      return
      end
c---- end of subroutine turbdif ----------------------------------------
c-----------------------------------------------------------------------

c ======================================================================
c  Here we do the neutral gas diffusion model, if isngon=1.
c  The diffusion is flux limited using the thermal flux.
c
c  If we are solving for the neutral parallel momentum (isupgon=1)
c  we only want the perpendicular flux of neutrals
c  and the corresponding perp velocity.
c  Franck-Condon neutrals complicate the picture with parallel momentum.
c  Not yet resolved in the coding.
c ======================================================================
c
      subroutine neudif

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,nu1,nu2,ffyi,ffyo
      integer iy1, methgx, methgy, iy2, jx
      logical isxyfl
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1
      Use(Share)    # geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methg,
                    # qfl,csh,qsh,cs,
                    # isupgon,iigsp,nlimgx,nlimgy,nlimiy,rld2dx,rld2dy
      Use(Coefeq)   # cngfx,cngfy,cngmom,cmwall,cdifg,rld2dxg,rld2dyg
      Use(Bcond)    # albedoo,albedoi
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    # xlinc,xrinc,yinc,ixm1,ixp1,stretcx
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # fngx,fngy,fngxy,fnix,fniy
      Use(Conduc)   # nuiz, nucx, nuix
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
	  
*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10
	  
      do 895 igsp = 1, ngsp

c.... First the flux in the x-direction

      do 888 iy = j4, j8
         do 887 ix = i1, i5
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
            qfl = flalfgxa(ix,igsp) * sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
            csh = (1-isgasdc) * cdifg(igsp) * sx(ix,iy) * gxf(ix,iy) *
     .                ave( stretcx(ix,iy)*vtn**2/nu1,
     .                     stretcx(ix2,iy)*vtnp**2/nu2 ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng +
     .                        rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                          0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
            qtgf = cngfx(igsp) * fgtdx(ix) * sx(ix,iy) *
     .            ave( stretcx(ix,iy)*gx(ix,iy)/nu1 ,
     .                 stretcx(ix2,iy)*gx(ix2,iy)/nu2 )
     .                     * (vtn**2 - vtnp**2)
            if (isupgon(igsp).eq.1) then # only 2->pol project.;up has rest
               csh = csh*(1 - rrv(ix,iy)**2)
               qtgf = qtgf*(1 - rrv(ix,iy)**2)
            endif
            vygtan(ix,iy,igsp) = 0.  # vygtan is grad(T) rad vel on x-face
            if (isnonog .eq. 1 .and. iy .le. ny) then
                  grdnv =(    ( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                       -( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               vygtan(ix,iy,igsp) = exp( 0.5*
     .                     (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                      ( cngfx(igsp) / (mg(igsp)*0.5*(nu1+nu2)) ) *
     .                                     ( grdnv/cos(angfx(ix,iy)) - 
     .                       (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))
     .                                                 * gxf(ix,iy) ) 
             if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
               vygtan(ix,iy,igsp) = 0.
             endif
             if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
               vygtan(ix,iy,igsp) = 0.
             endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - vygtan(ix,iy,igsp)*sx(ix,iy)
            if (isupgon(igsp) .eq. 1) then
              qtgf = qtgf + rrv(ix,iy)*up(ix,iy,iigsp)*sx(ix,iy)
            endif
            nconv = 2.0*(ng(ix,iy,igsp)*ng(ix2,iy,igsp)) /
     .                  (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
c...   Use upwind for "convective" grad T term if methgx .ne. 2
            if(methgx.ne.2) nconv =
     .                         ng(ix ,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ng(ix2,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ng(ix,iy,igsp)-ng(ix2,iy,igsp)) + qtgf * nconv
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  bndry face - should not matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floxg
            floxg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = 
     .              floxg(ix,iy) + cngflox(1)*sx(ix,iy)*uu(ix,iy,1)

  887    continue
         conxg(nx+1,iy) = 0
  888 continue

c.... Now the flux in the y-direction

      do 890 iy = j1, j5
         do 889 ix = i4, i8
	    t0 = max(tg(ix,iy,igsp),temin*ev)
	    t1 = max(tg(ix,iy+1,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix,iy+1,igsp) + vtnp/lgmax(igsp)
            qfl = flalfgya(iy,igsp) * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .                              ( ngy0(ix,iy,igsp)*gy(ix,iy) + 
     .                                ngy1(ix,iy,igsp)*gy(ix,iy+1) ) / 
     .                                     (8*(gy(ix,iy)+gy(ix,iy+1)))
             csh = (1-isgasdc) * cdifg(igsp) *sy(ix,iy)/(dynog(ix,iy)) *
     .                                  ave(vtn**2/nu1, vtnp**2/nu2) +
     .            isgasdc * sy(ix,iy) * gyf(ix,iy) * difcng +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*(1/gyf(ix,iy))*
     .                       0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))
            qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                            ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2)
     .                                      * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0 ) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) )
              endif
              qtgf = cngfy(igsp) * fgtdy(iy)* sy(ix,iy) * 
     .                           ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ngy0(ix,iy,igsp)-ngy1(ix,iy,igsp)) + qtgf*nconv
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floyg
	    floyg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) =  
     .             floyg(ix,iy)+cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)

  889    continue
  890 continue

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             ng(0:nx+1,0:ny+1,igsp),fngx(0:nx+1,0:ny+1,igsp),
     .             fngy(0:nx+1,0:ny+1,igsp),0,methg)

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then

         do 8905 iy = j1, min(j6, ny)
            iy1 = max(iy-1,0)
            do 8904 ix = i1, min(i6, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
	       t0 = max(tg(ix ,iy,igsp),temin*ev)
	       t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn =  sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               nu1 = nuix(ix ,iy,igsp) + vtn/lgmax(igsp)
               nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
ccc            MER: Set flag to apply xy flux limit except at target plates
               isxyfl = .true.
               do jx = 1, nxpt
                  if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                 (ix==ixrb(jx).and.ixmxbcl==1) ) isxyfl = .false.
               enddo
               if (methgx .eq. 6) then  # log interpolation
               grdnv =(   ( fym (ix,iy,1)*log(ng(ix2,iy1 ,igsp)) +
     .                      fy0 (ix,iy,1)*log(ng(ix2,iy  ,igsp)) +
     .                      fyp (ix,iy,1)*log(ng(ix2,iy+1,igsp)) + 
     .                      fymx(ix,iy,1)*log(ng(ix ,iy1 ,igsp)) +
     .                      fypx(ix,iy,1)*log(ng(ix, iy+1,igsp)) )
     .                   -( fym (ix,iy,0)*log(ng(ix ,iy1 ,igsp)) +
     .                      fy0 (ix,iy,0)*log(ng(ix ,iy  ,igsp)) +
     .                      fyp (ix,iy,0)*log(ng(ix ,iy+1,igsp)) +
     .                      fymx(ix,iy,0)*log(ng(ix4,iy1 ,igsp)) + 
     .                      fypx(ix,iy,0)*log(ng(ix6,iy+1,igsp)) ) )/ 
     .                                                  dxnog(ix,iy)
               elseif (methgx .eq. 7) then  # inverse interpolation
               grdnv =( 1/(fym (ix,iy,1)/ng(ix2,iy1 ,igsp) + 
     .                     fy0 (ix,iy,1)/ng(ix2,iy  ,igsp) +
     .                     fyp (ix,iy,1)/ng(ix2,iy+1,igsp) + 
     .                     fymx(ix,iy,1)/ng(ix ,iy1 ,igsp) +
     .                     fypx(ix,iy,1)/ng(ix, iy+1,igsp))
     .                - 1/(fym (ix,iy,0)/ng(ix ,iy1 ,igsp) +
     .                     fy0 (ix,iy,0)/ng(ix ,iy  ,igsp) +
     .                     fyp (ix,iy,0)/ng(ix ,iy+1,igsp) +
     .                     fymx(ix,iy,0)/ng(ix4,iy1 ,igsp) + 
     .                     fypx(ix,iy,0)/ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               else                   # linear interpolation
               grdnv =( (fym (ix,iy,1)*ng(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*ng(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*ng(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*ng(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*ng(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*ng(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*ng(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*ng(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*ng(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               endif
               difgx2 = ave( tg(ix ,iy,igsp)/nu1,
     .                       tg(ix2,iy,igsp)/nu2 )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
              if (methgx .eq. 6) then
               fngxy(ix,iy,igsp) =  exp( 0.5*
     .                     (log(ng(ix2,iy,igsp))+log(ng(ix,iy,igsp))) )*
     .                               difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                     (log(ng(ix2,iy,igsp)) - log(ng(ix,iy,igsp)))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              else
               fngxy(ix,iy,igsp) = difgx2*( grdnv/cos(angfx(ix,iy)) -
     .                             (ng(ix2,iy,igsp) - ng(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              endif
c...  Now flux limit with flalfgxy
               t0 = max(tg(ix,iy,igsp),temin*ev)
               t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn = sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               qfl = flalfgxya(ix,igsp)*sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
ccc   MER NOTE:  no xy flux limit for ix=0 or ix=nx in original code
               if (isxyfl) then
                  fngxy(ix,iy,igsp) = fngxy(ix,iy,igsp) /
     .                           sqrt( 1 + (fngxy(ix,iy,igsp)/qfl)**2 )
               endif
 8904       continue
 8905    continue

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Flux-limit the total poloidal flux here
            do iy = j4, j8
               do ix = i1, i5
                  ix2 = ixp1(ix,iy)
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)-fngxy(ix,iy,igsp)
                  t0 = max(tg(ix,iy,igsp),temin*ev)
                  t1 = max(tg(ix2,iy,igsp),temin*ev)
                  vtn = sqrt( t0/mg(igsp) )
                  vtnp = sqrt( t1/mg(igsp) )
                  qfl = flalfgnx * sx(ix,iy)*(vtn + vtnp)*rt8opi*
     .                             (ng(ix,iy,igsp)+ng(ix2,iy,igsp))/16
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/
     .                              sqrt(1 + (fngx(ix,iy,igsp)/qfl)**2)
c ...          adjust fluxes to prevent pumpout
                    fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/( 1 - 2*nlimgx +
     .                          nlimgx*(ng(ix2,iy,igsp)/ng(ix,iy,igsp) +
     .                                  ng(ix,iy,igsp)/ng(ix2,iy,igsp)) )
               enddo
            enddo
c ...   adjust y-fluxes to prevent pumpout
            do iy = j1, j5    # same loop ranges as for fngy in fd2tra
               do ix = i4, i8
                     fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/( 1 - 2*nlimgy +
     .                         nlimgy*(ng(ix,iy+1,igsp)/ng(ix,iy,igsp)+
     .                                 ng(ix,iy,igsp)/ng(ix,iy+1,igsp)) )
	            t0 = max(tg(ix,iy,igsp),temin*ev)
	            t1 = max(tg(ix,iy+1,igsp),temin*ev)
                    vtn = sqrt( t0/mg(igsp) )
                    vtnp = sqrt( t1/mg(igsp) )
                    qfl = flalfgny*sy(ix,iy)*(vtn+vtnp)*rt8opi*
     .                        (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp))/16
                    fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/
     .                            sqrt(1 + (fngy(ix,iy,igsp)/qfl)**2)
               enddo
            enddo     

      endif
c...  Finished with nonorthogonal mesh part

c ... Calculate the neutral flow velocity from v = flux/ng; these are
c ... diagnostic if isupgon=0, but used for uu and vy of the inertial
c ... gas if isupgon=1.  However, even when isupgon=1, the particle
c ... fluxes are fngx -> fnix and fngy -> fniy in pandf, i.e., methg
c ... determines the differencing for the inertial particle fluxes, not 
c ... methn
      do iy = j1, j5
         do ix = i1,i5
            ix1 = ixp1(ix,iy)
            if (1.-rrv(ix,iy) > 1.e-4) then #combine binormal & par components
              uug(ix,iy,igsp) = fngx(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                      *sx(ix,iy) )
            else   # binormal component negligable small
               uug(ix,iy,igsp) = up(ix,iy,iigsp)
            endif
            vyg(ix,iy,igsp) = fngy(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
     .                                                      *sy(ix,iy) )
c --------------- transfer inertial gas velocities to neutral ion species
            if (isupgon(igsp).eq.1) then
               vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
cc               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            end if
         enddo
cfw ----   If doing only the outer half we want this boundary condition:
         if (iy.le.iysptrx2(1) .and. isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp) = 0
      enddo

c **- loop for uu just as in the previous version - needed for correct Jac?
      if (isupgon(igsp) .eq. 1) then
         do iy = j4, j6
            do ix = i1, i6
               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            enddo
         enddo
      endif

c.... Calculate the residual for the gas equation for diffusive neutral case

      if (isupgon(igsp).eq.0) then
         do 892 iy = j2, j5
	      if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
	      fngx(nxc-1,iy,igsp) = 0.
	      fngx(nxc,  iy,igsp) = 0.
	      fngx(nxc+1,iy,igsp) = 0.
	    endif
            if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
            if (nxpt==2) fngx(ixrb(1)+1,iy,igsp)=0.
       	    do 891 ix = i2, i5
               ix1 = ixm1(ix,iy)
               resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .                         psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .                       + volpsorg(ix,iy,igsp)                    
     .                       - fngx(ix,iy,igsp) + fngx(ix1,iy  ,igsp)
     .             - fluxfacy*(fngy(ix,iy,igsp) - fngy(ix ,iy-1,igsp))
     .                       + psgov_use(ix,iy,igsp)*vol(ix,iy)
               if (igsp.eq.1 .and. ishymol.eq.1) resng(ix,iy,igsp) =
     .                                  resng(ix,iy,igsp)+psordis(ix,iy)
 891        continue
 892     continue
      endif

 895  continue   # end of igsp loop from the beginning of subroutine

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                         nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif

      return

      end
c
c --------------------------------------------------------------------------
c END subroutine neudif
c --------------------------------------------------------------------------

c ======================================================================
c  Below neudifpg similar to neudif, except that the neutral pressure, ng*tg, 
c  is the dependent variable differenced, rather than ng and tg separately
c  Here we do the neutral gas diffusion model, if isngon=1.
c  The diffusion is flux limited using the thermal flux.
c
c  If we are solving for the neutral parallel momentum (isupgon=1)
c  we only want the perpendicular flux of neutrals
c  and the corresponding perp velocity.
c  Franck-Condon neutrals complicate the picture with parallel momentum.
c  Not yet resolved in the coding.
c ======================================================================
c
      subroutine neudifpg

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,nu1,nu2,tgf,ffyi,ffyo
      real tsngxlog,tsngylog,tsngfd2,tsngfxy
      real dndym1,dndy0,dndyp1,d2ndy20,d2ndy2p1,d3ndy3
      real dndxm1,dndx0,dndxp1,d2ndx20,d2ndx2p1,d3ndx3
      real flalfgx_adj, flalfgy_adj, flalfgxy_adj, ngxface, ngyface
      integer iy1, methgx, methgy, iy2, jx, jfld, ifld
      integer iym1,iyp1,iyp2,ixm1b,ixp1b,ixp2b
      logical isxyfl
      real(Size4) sec4, gettime
      # Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a

      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1
      Use(Share)    # geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methg,
                    # qfl,csh,qsh,cs,
                    # isupgon,iigsp,nlimgx,nlimgy,nlimiy,rld2dx,rld2dy
      Use(Coefeq)   # cngfx,cngfy,cngmom,cmwall,cdifg,rld2dxg,rld2dyg
      Use(Bcond)    # albedoo,albedoi
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    # xlinc,xrinc,yinc,ixm1,ixp1
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # fngx,fngy,fngxy,fnix,fniy
      Use(Conduc)   # nuiz, nucx, nuix
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
      Use(Timing)   # ttngxlog,ttngylog,ttngfd2,ttngfxy

      Use(Ext_neutrals) # get_neutral_moments, ...
      Use(MCN_dim)      # ngsp, ...
      Use(MCN_sources)  # cfneut_sng, cfneutdiv_fng, ... mcfngx, mcfngy, ...
      Use(Interp)		# ngs, tgs 
      Use(Bfield)   # rbfbt 

*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10
	  
c      write (*,*) "neudifpg"
      do 895 igsp = 1, ngsp

c *********************************************
c.... First the flux in the x-direction
c *********************************************

c ..Timing;initialize 
      if(istimingon==1) tsngxlog = gettime(sec4)

      do 888 iy = j4, j8
         do 887 ix = i1, i5
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            ngxface = 0.5*(ng(ix,iy,igsp)+ng(ix2,iy,igsp))
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
	    tgf = 0.5*(tg(ix,iy,igsp)+tg(ix2,iy,igsp))
            flalfgx_adj = flalfgxa(ix,igsp)*( 1. +
     .                    (cflbg*ngbackg(igsp)/ngxface)**inflbg )
            qfl = flalfgx_adj * sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
            csh = (1-isgasdc) * cdifg(igsp) * sx(ix,iy) * gxf(ix,iy) *
     .                (1/mg(igsp))* ave( 1./nu1,1./nu2 ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng / tgf +
     .                        rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                    0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))/tgf
            qtgf = alftng * fgtdx(ix) * sx(ix,iy) *
     .            ave( gx(ix,iy)/nu1 ,gx(ix2,iy)/nu2 )
     .                     * (vtn**2 - vtnp**2)
            if (isupgon(igsp).eq.1) then # only 2->pol project.;up has rest
               csh = csh*(1 - rrv(ix,iy)**2)
               qtgf = qtgf*(1 - rrv(ix,iy)**2)
            endif
            vygtan(ix,iy,igsp) = 0.  # vygtan is grad(T) rad vel on x-face
                                     # vygtan) only from thermal force
            if (isnonog .eq. 1 .and. iy .le. ny) then
                  grdnv =(    ( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                       -( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               vygtan(ix,iy,igsp) = exp( 0.5*
     .                     (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                      ( alftng / (mg(igsp)*0.5*(nu1+nu2)) ) *
     .                                     ( grdnv/cos(angfx(ix,iy)) - 
     .                       (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))
     .                                                 * gxf(ix,iy) ) 
             if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
               vygtan(ix,iy,igsp) = 0.
             endif
             if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
               vygtan(ix,iy,igsp) = 0.
             endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - vygtan(ix,iy,igsp)*sx(ix,iy)
            if (isupgon(igsp) .eq. 1) then
              qtgf = qtgf + rrv(ix,iy)*up(ix,iy,iigsp)*sx(ix,iy)
            endif
            nconv = 2.0*(ng(ix,iy,igsp)*ng(ix2,iy,igsp)) /
     .                  (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
c...   Use upwind for "convective" grad T term if methgx .ne. 2
            if(methgx.ne.2) nconv =
     .                         ng(ix ,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ng(ix2,iy,igsp)*0.5*(1-sign(1.,qtgf))

            qsh = csh * (pg(ix,iy,igsp)-pg(ix2,iy,igsp)) + qtgf * nconv
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  boundary face - shouldnot matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  themal force temperature gradient term is included in floxg
	    floxg(ix,iy) = (qtgf/tgf) / (1 + qr**flgamg)**(1/flgamg)
c...  now add the ion convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = floxg(ix,iy) +
     .                            cngflox(1)*sx(ix,iy)*uu(ix,iy,1)/tgf
c...  For one impurity, add convect vel for elastic scattering with ions
         if(igsp == 2 .and. ishymol == 0) then #Caution; need to weight uu
           do ifld = nhsp+1, nisp
             floxg(ix,iy) = floxg(ix,iy) +
     .                cngniflox(ifld,igsp)*sx(ix,iy)*uu(ix,iy,ifld)/tgf
           enddo
         endif  

  887    continue
         conxg(nx+1,iy) = 0
  888 continue

c ..Timing; add info if timing is on
      if(istimingon==1) ttngxlog=ttngxlog+(gettime(sec4)-tsngxlog)

c *******************************************************
c.... Now the flux in the y-direction
c *******************************************************

c ..Timing; initiate time for y-direction calc
      if(istimingon==1) tsngylog = gettime(sec4)

      do 890 iy = j1, j5
         do 889 ix = i4, i8
            ngyface = 0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
	    t0 = max(tg(ix,iy,igsp),tgmin*ev)
	    t1 = max(tg(ix,iy+1,igsp),tgmin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix,iy+1,igsp) + vtnp/lgmax(igsp)
	    tgf = 0.5*(tg(ix,iy,igsp)+tg(ix,iy+1,igsp))
            flalfgy_adj = flalfgya(iy,igsp)*( 1. +
     .                   (cflbg*ngbackg(igsp)/ngyface)**inflbg )
            qfl = flalfgy_adj * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .                              ( ngy0(ix,iy,igsp)*gy(ix,iy) + 
     .                                ngy1(ix,iy,igsp)*gy(ix,iy+1) ) / 
     .                                     (8*(gy(ix,iy)+gy(ix,iy+1)))
            if (iy==0) then  #at bdry, ng ave to avoid ng->0 prob
              qfl = flalfgy_adj * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .              (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) / 8.
            endif
            csh = (1-isgasdc) * (cdifg(igsp) *sy(ix,iy)/dynog(ix,iy)) *
     .                          (1/mg(igsp))* ave(1./nu1, 1./nu2) +
     .            isgasdc * sy(ix,iy) * difcng /(dynog(ix,iy)*tgf) +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*dynog(ix,iy)*
     .                     0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))/tgf

            qtgf = alftng * fgtdy(iy) * sy(ix,iy) * 
     .                            ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2)
     .                                      * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0 ) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) )
              endif
              qtgf = alftng * fgtdy(iy)* sy(ix,iy) * 
     .                           ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (pgy0(ix,iy,igsp)-pgy1(ix,iy,igsp)) + qtgf*nconv
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  thermal force temperature gradient term is included in floyg
	    floyg(ix,iy) = (qtgf/tgf) / (1 + qr**flgamg)**(1/flgamg)
c...  now add ion convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) = floyg(ix,iy) +
     .                            cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)/tgf
c...  For impurities, add convect vel for elastic scattering with ions
         if(igsp == 2 .and. ishymol == 0) then #Caution; need to weight uu 
           do ifld = nhsp+1, nisp
             floyg(ix,iy) = floyg(ix,iy) +
     .                cngnifloy(ifld,igsp)*sy(ix,iy)*vy(ix,iy,ifld)/tgf
           enddo
         endif

  889    continue
  890 continue

c ..Timing; add increment if timing is on
      if(istimingon==1) ttngylog=ttngylog+(gettime(sec4)-tsngylog)

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------
c ..Timing
      if(istimingon==1) tsngfd2 = gettime(sec4)

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             pg(0:nx+1,0:ny+1,igsp),fngx(0:nx+1,0:ny+1,igsp),
     .             fngy(0:nx+1,0:ny+1,igsp),0,methg)
c ..Timing
      if(istimingon==1) ttngfd2=ttngfd2+(gettime(sec4)-tsngfd2)

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then
c ..Timing
      if(istimingon==1) tsngfxy = gettime(sec4)

         do 8905 iy = j1, min(j6, ny)
            iy1 = max(iy-1,0)
            do 8904 ix = i1, min(i6, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
	       t0 = max(tg(ix ,iy,igsp),temin*ev)
	       t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn =  sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               nu1 = nuix(ix ,iy,igsp) + vtn/lgmax(igsp)
               nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
ccc            MER: Set flag to apply xy flux limit except at target plates
               isxyfl = .true.
               do jx = 1, nxpt
                 if( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .               (ix==ixrb(jx).and.ixmxbcl==1) ) isxyfl = .false.
               enddo
               if (methgx .eq. 6) then  # log interpolation
               grdnv =(   ( fym (ix,iy,1)*log(pg(ix2,iy1 ,igsp)) + 
     .                      fy0 (ix,iy,1)*log(pg(ix2,iy  ,igsp)) +
     .                      fyp (ix,iy,1)*log(pg(ix2,iy+1,igsp)) + 
     .                      fymx(ix,iy,1)*log(pg(ix ,iy1 ,igsp)) +
     .                      fypx(ix,iy,1)*log(pg(ix, iy+1,igsp)) )
     .                   -( fym (ix,iy,0)*log(pg(ix ,iy1 ,igsp)) +
     .                      fy0 (ix,iy,0)*log(pg(ix ,iy  ,igsp)) +
     .                      fyp (ix,iy,0)*log(pg(ix ,iy+1,igsp)) +
     .                      fymx(ix,iy,0)*log(pg(ix4,iy1 ,igsp)) + 
     .                      fypx(ix,iy,0)*log(pg(ix6,iy+1,igsp)) ) )/ 
     .                                                  dxnog(ix,iy)
               elseif (methgx .eq. 7) then # inverse interpolation
               grdnv =( 1/(fym (ix,iy,1)/pg(ix2,iy1 ,igsp) + 
     .                     fy0 (ix,iy,1)/pg(ix2,iy  ,igsp) +
     .                     fyp (ix,iy,1)/pg(ix2,iy+1,igsp) + 
     .                     fymx(ix,iy,1)/pg(ix ,iy1 ,igsp) +
     .                     fypx(ix,iy,1)/pg(ix, iy+1,igsp))
     .                - 1/(fym (ix,iy,0)/pg(ix ,iy1 ,igsp) +
     .                     fy0 (ix,iy,0)/pg(ix ,iy  ,igsp) +
     .                     fyp (ix,iy,0)/pg(ix ,iy+1,igsp) +
     .                     fymx(ix,iy,0)/pg(ix4,iy1 ,igsp) + 
     .                     fypx(ix,iy,0)/pg(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               else                   # linear interpolation
               grdnv =( (fym (ix,iy,1)*pg(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*pg(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*pg(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*pg(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*pg(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*pg(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*pg(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*pg(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*pg(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*pg(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               endif
               difgx2 = ave( 1./nu1,
     .                       1./nu2 )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
              if (methgx .eq. 6) then
               fngxy(ix,iy,igsp) =  exp( 0.5*
     .                     (log(pg(ix2,iy,igsp))+log(pg(ix,iy,igsp))) )*
     .                               difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                     (log(pg(ix2,iy,igsp)) - log(pg(ix,iy,igsp)))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              else
               fngxy(ix,iy,igsp) = difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                             (pg(ix2,iy,igsp) - pg(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              endif
c...  Now flux limit with flalfgxy
               ngxface = 0.5*(ng(ix,iy,igsp)+ng(ix2,iy,igsp))
               t0 = max(tg(ix,iy,igsp),temin*ev)
               t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn = sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               flalfgxy_adj = flalfgxya(ix,igsp)*( 1. +
     .                     (cflbg*ngbackg(igsp)/ngxface)**inflbg )
               qfl = flalfgxy_adj*sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
ccc   MER NOTE:  no xy flux limit for ix=0 or ix=nx in original code
               if (isxyfl) then
                  fngxy(ix,iy,igsp) = fngxy(ix,iy,igsp) /
     .                           sqrt( 1 + (fngxy(ix,iy,igsp)/qfl)**2 )
               endif
 8904       continue
 8905    continue
c ..Timing
      if(istimingon==1) ttngfxy=ttngfxy+(gettime(sec4)-tsngfxy)

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Flux-limit the total poloidal flux here
            do iy = j4, j8
               do ix = i1, i5
                  ix2 = ixp1(ix,iy)
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)-fngxy(ix,iy,igsp)
                  t0 = max(tg(ix,iy,igsp),temin*ev)
                  t1 = max(tg(ix2,iy,igsp),temin*ev)
                  vtn = sqrt( t0/mg(igsp) )
                  vtnp = sqrt( t1/mg(igsp) )
                  qfl = flalfgnx * sx(ix,iy)*(vtn + vtnp)*rt8opi*
     .                             (ng(ix,iy,igsp)+ng(ix2,iy,igsp))/16
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/
     .                              sqrt(1 + (fngx(ix,iy,igsp)/qfl)**2)
c ...          adjust fluxes to prevent pumpout
                    fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/( 1 - 2*nlimgx +
     .                          nlimgx*(ng(ix2,iy,igsp)/ng(ix,iy,igsp) +
     .                                  ng(ix,iy,igsp)/ng(ix2,iy,igsp)) )
               enddo
            enddo
c ...   adjust y-fluxes to prevent pumpout
            do iy = j1, j5    # same loop ranges as for fngy in fd2tra
               do ix = i4, i8
                     fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/( 1-2*nlimgy +
     .                         nlimgy*(ng(ix,iy+1,igsp)/ng(ix,iy,igsp)+
     .                                ng(ix,iy,igsp)/ng(ix,iy+1,igsp)) )
	            t0 = max(tg(ix,iy,igsp),temin*ev)
	            t1 = max(tg(ix,iy+1,igsp),temin*ev)
                    vtn = sqrt( t0/mg(igsp) )
                    vtnp = sqrt( t1/mg(igsp) )
                    qfl = flalfgny*sy(ix,iy)*(vtn+vtnp)*rt8opi*
     .                        (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp))/16
                    fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/
     .                            sqrt(1 + (fngy(ix,iy,igsp)/qfl)**2)
               enddo
            enddo     

      endif
c...  Finished with nonorthogonal mesh part

c...  Add 4th order radial diffusion op; damp grid-scale oscill
        if (abs(difgy4order(igsp)) > 1.e-50 .and. isngon(igsp)==1) then
          do iy = j2p, j5m   #limits to range iy=1:ny-1 for fngy4ord
            iym1 = max(iy-1,0)
            iyp1 = min(iy+1,ny+1)
            iyp2 = min(iy+2,ny+1)
            do ix = i4, i8
              dndym1 = (ng(ix,iy,igsp)-ng(ix,iym1,igsp))*gyf(ix,iym1)
              dndy0 = (ng(ix,iyp1,igsp)-ng(ix,iy,igsp))*gyf(ix,iy)
              dndyp1 = (ng(ix,iyp2,igsp)-ng(ix,iyp1,igsp))*gyf(ix,iyp1)
              d2ndy20 = (dndy0 - dndym1)*gy(ix,iy)
              d2ndy2p1 = (dndyp1 - dndy0)*gy(ix,iyp1)
              d3ndy3 = (d2ndy2p1 - d2ndy20)*gyf(ix,iy)
              fngy4ord(ix,iy,igsp) = difgy4order(igsp)*d3ndy3*sy(ix,iy)/
     .                                                  gyf(ix,iy)**2
              fngy(ix,iy,igsp) = fngy(ix,iy,igsp) + fngy4ord(ix,iy,igsp)
            enddo
          enddo
        endif

c...  Add 4th order poloidal diffusion op; damp grid-scale oscill
C...  NOTE: PRESENTLY ONLY CODED FOR SIMPLY-CONNECTED DOMAIN
        if (abs(difgx4order(igsp)) > 1.e-50 .and. isngon(igsp)==1) then
          do ix = i2p, i5m   #limits to range ix=1:nx-1 for fngx4ord
            ixm1b = max(ix-1,0)
            ixp1b = min(ix+1,nx+1)
            ixp2b = min(ix+2,nx+1)
            do iy = j4, j8
              dndxm1 = (ng(ix,iy,igsp)-ng(ixm1b,iy,igsp))*gxf(ixm1b,iy)
              dndx0 = (ng(ixp1b,iy,igsp)-ng(ix,iy,igsp))*gxf(ix,iy)
              dndxp1 = (ng(ixp2b,iy,igsp)-ng(ixp1b,iy,igsp))*gxf(ixp1b,iy)
              d2ndx20 = (dndx0 - dndxm1)*gx(ix,iy)
              d2ndx2p1 = (dndxp1 - dndx0)*gx(ixp1b,iy)
              d3ndx3 = (d2ndx2p1 - d2ndx20)*gxf(ix,iy)
              fngx4ord(ix,iy,igsp) = difgx4order(igsp)*d3ndx3*sx(ix,iy)/
     .                                                  gxf(ix,iy)**2
              fngx(ix,iy,igsp) = fngx(ix,iy,igsp) + fngx4ord(ix,iy,igsp)
            enddo
          enddo
        endif

c ... Calculate the neutral flow velocity from v = flux/ng; these are
c ... diagnostic if isupgon=0, but used for uu and vy of the inertial
c ... gas if isupgon=1.  However, even when isupgon=1, the particle
c ... fluxes are fngx -> fnix and fngy -> fngy in pandf, i.e., methg
c ... determines the differencing for the inertial particle fluxes, not 
c ... methn
      do iy = j1, j5
        do ix = i1,i5
          ix1 = ixp1(ix,iy)
          if (1.-rrv(ix,iy) > 1.e-4 .or. isupgon(igsp)==0) then 
                           #combine binormal/par comps or x only diffusive
             uug(ix,iy,igsp) = fngx(ix,iy,igsp) / (
     .                      0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                    *sx(ix,iy) )
c ...    remove nonorthogonal vy component from uug for seic contrib
             uuxg(ix,iy,igsp) = (fngx(ix,iy,igsp)+fngxy(ix,iy,igsp))/ (
     .                      0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                    *sx(ix,iy) )
          else   # binormal component negligable small
             uug(ix,iy,igsp) = up(ix,iy,iigsp)
          endif
          vyg(ix,iy,igsp) = fngy(ix,iy,igsp) / (
     .                      0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
     .                                                      *sy(ix,iy) )
c --------------- transfer inertial gas velocities to neutral ion species
          if (isupgon(igsp).eq.1) then
             vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
          endif
        enddo
cfw ----   If doing only the outer half we want this boundary condition:
        if(iy.le.iysptrx2(1).and.isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp)=0
      enddo

c **- loop for uu just as in the previous version - needed for correct Jac?
      if (isupgon(igsp) .eq. 1) then
         do iy = j4, j6
            do ix = i1, i6
               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
               v2(ix,iy,iigsp) = ( uuxg(ix,iy,igsp) 
     .                            - up(ix,iy,iigsp)*rrv(ix,iy) )
     .                       /(rbfbt(ix,iy) + rbfbt(ixp1(ix,iy),iy))*2.
            enddo
         enddo
      endif

c.... Calculate the residual for the gas equation for diffusive neutral case

      if (isupgon(igsp).eq.0) then
        do 892 iy = j2, j5
          if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then 
	     fngx(nxc-1,iy,igsp) = 0.
	     fngx(nxc,  iy,igsp) = 0.
	     fngx(nxc+1,iy,igsp) = 0.
	  endif
          if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
          if (nxpt==2.and.ixmxbcl==1) fngx(ixrb(1)+1,iy,igsp)=0.
          do 891 ix = i2, i5
            ix1 = ixm1(ix,iy)

c ... 2016/09/16 IJ: coding to blend MC neutral flux !!! here ***
c ... is it correct to use ng instead of ni??? i.e. will ng enter jacobian?
            resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .               psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .             + volpsorg(ix,iy,igsp)
     .             + psgov_use(ix,iy,igsp)*vol(ix,iy)
            if (igsp.eq.1 .and. ishymol.eq.1)
     .            resng(ix,iy,igsp) = resng(ix,iy,igsp) + psordis(ix,iy)
            resng(ix,iy,igsp) = resng(ix,iy,igsp) - cfneutdiv*
     .          cfneutdiv_fng*((fngx(ix,iy,igsp) - fngx(ix1,iy, igsp)) +
     .          fluxfacy*(fngy(ix,iy,igsp) - fngy(ix,iy-1,igsp)) )

c ... IJ 2016/10/19 add MC neut flux if flags set
             if (get_neutral_moments .and. cmneutdiv_fng .ne. 0.0) then  
               jfld=1  
               sng_ue(ix,iy,jfld) = - ( 
     .                     (fngx_ue(ix,iy,jfld)-fngx_ue(ix1,iy, jfld))
     .           +fluxfacy*(fngy_ue(ix,iy,jfld)-fngy_ue(ix,iy-1,jfld)) )
     .           *( (ng(ix,iy,jfld)*ti(ix,iy))/
     .                                  (ng(ix,iy,jfld)*ti(ix,iy)) )
               resng(ix,iy,igsp) = resng(ix,iy,igsp) + 
     .                     cmneutdiv*cmneutdiv_fng*sng_ue(ix,iy,igsp)
             endif

 891      continue
 892    continue
      endif

 895  continue   # end of igsp loop from the beginning of subroutine

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
            lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                         nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
            lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
               lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
               lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif

      return

      end
c
c --------------------------------------------------------------------------
c END subroutine neudifpg
c --------------------------------------------------------------------------


c --------------------------------------------------------------------------
c   Below subroutine neudifl is just like subroutine neudif, except that the
c   log of the gas density is used, and then converted back to give the 
c   physically meaningful gas variables (flux, velocity, etc)
c --------------------------------------------------------------------------

      subroutine neudifl

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,ffyo,ffyi
      integer iy1, methgx, methgy, iy2, jx
      logical isxyfl
      # Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1
      Use(Share)    # geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methg,
                    # qfl,csh,qsh,cs,
                    # isupgon,iigsp,nlimgx,nlimgy,nlimiy,rld2dx,rld2dy
      Use(Coefeq)   # cngfx,cngfy,cngmom,cmwall,cdifg,rld2dxg,rld2dyg
      Use(Bcond)    # albedoo,albedoi
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    # xlinc,xrinc,yinc,ixm1,ixp1,stretcx
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # fngx,fngy,fngxy,fnix,fniy
      Use(Conduc)   # nuiz, nucx, nuix
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
	  
*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10
	  
      do 895 igsp = 1, ngsp

c.... First the flux in the x-direction

      do 888 iy = j4, j8
         do 887 ix = i1, i5
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgxa(ix,igsp) * sx(ix,iy) * (vtn + vtnp)*rt8opi/8
            csh = (1-isgasdc) * cdifg(igsp) * sx(ix,iy) * gxf(ix,iy) *
     .                ave( stretcx(ix,iy)*vtn**2/nuix(ix,iy,igsp),
     .                     stretcx(ix2,iy)*vtnp**2/nuix(ix2,iy,igsp) ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng +
     .                        rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                          0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
            qtgf = cngfx(igsp) * fgtdx(ix) * sx(ix,iy) *
     .            ave( stretcx(ix,iy)*gx(ix,iy)/nuix(ix,iy,igsp) ,
     .                 stretcx(ix2,iy)*gx(ix2,iy)/nuix(ix2,iy,igsp) )
     .                     * (vtn**2 - vtnp**2)
            if (isupgon(igsp).eq.1) then # only 2->pol project.;up has rest
               csh = csh*(1 - rrv(ix,iy)**2)
               qtgf = qtgf*(1 - rrv(ix,iy)**2)
            endif
            vygtan(ix,iy,igsp) = 0.  # vygtan is grad(T) rad vel on x-face
            if (isnonog .eq. 1 .and. iy .le. ny) then
                  grdnv =(    ( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                       -( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               vygtan(ix,iy,igsp) = exp( 0.5*
     .                     (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                                  ( cngfx(igsp) / (mg(igsp)*0.5*
     .                         (nuix(ix,iy,igsp)+nuix(ix2,iy,igsp))) ) *
     .                             ( grdnv/cos(angfx(ix,iy)) - 
     .                       (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))
     .                                                 * gxf(ix,iy) ) 
               if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
                  vygtan(ix,iy,igsp) = 0.
               endif
               if (nxpt==2 .and. ix==ixrb(1)+1) then
                  vygtan(ix,iy,igsp) = 0.
               endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - vygtan(ix,iy,igsp)*sx(ix,iy)
            if (isupgon(igsp) .eq. 1) then
              qtgf = qtgf + rrv(ix,iy)*up(ix,iy,iigsp)*sx(ix,iy)
            endif
            qsh = csh * (lng(ix,iy,igsp)-lng(ix2,iy,igsp)) + qtgf
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  boundary face - shouldnt matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floxg
            floxg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = 
     .              floxg(ix,iy) + cngflox(1)*sx(ix,iy)*uu(ix,iy,1)
         floxg(ix,iy) = floxg(ix,iy)*2/(lng(ix,iy,igsp)+lng(ix2,iy,igsp))

  887    continue
         conxg(nx+1,iy) = 0
  888 continue

c.... Now the flux in the y-direction

      do 890 iy = j1, j5
         do 889 ix = i4, i8
	    t0 = max(tg(ix,iy,igsp),temin*ev)
	    t1 = max(tg(ix,iy+1,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgya(iy,igsp) * sy(ix,iy) * (vtn + vtnp)*rt8opi/8
            csh = (1-isgasdc) * (cdifg(igsp) *sy(ix,iy)/dynog(ix,iy)) *
     .                            ave( vtn**2/nuix(ix,iy,igsp) ,
     .                                 vtnp**2/nuix(ix,iy+1,igsp) ) +
     .            isgasdc * sy(ix,iy) * difcng / dynog(ix,iy) +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*dynog(ix,iy)*
     .                       0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))
c               csh = sy(ix,iy) * gyf(ix,iy) * ( (vtn**2+vtnp**2)/
c     .                 (nuix(ix,iy,igsp)+nuix(ix,iy+1,igsp)) )
            qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                     ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                          gy(ix,iy+1)/nuix(ix,iy+1,igsp) )
     .                    * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0 ) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp )) )
              endif
              qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                      ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                           gy(ix,iy+1)/nuix(ix,iy+1,igsp) ) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (lng(ix,iy,igsp)-lng(ix,iy+1,igsp)) + qtgf
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floyg
	    floyg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) =  
     .             floyg(ix,iy)+cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)

         floyg(ix,iy)=floyg(ix,iy)*2/(lng(ix,iy,igsp)+lng(ix,iy+1,igsp))
  889    continue
  890 continue

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             lng(0:nx+1,0:ny+1,igsp),flngx(0:nx+1,0:ny+1,igsp),
     .             flngy(0:nx+1,0:ny+1,igsp),0,methg)

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then

         do 8905 iy = j1, min(j6, ny)
            iy1 = max(iy-1,0)
            do 8904 ix = i1, min(i6, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
ccc            MER: Set flag to apply xy flux limit except at target plates
               isxyfl = .true.
               do jx = 1, nxpt
                  if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                 (ix==ixrb(jx).and.ixmxbcl==1) )isxyfl = .false.
               enddo
               grdnv =( (fym (ix,iy,1)*lng(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*lng(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*lng(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*lng(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*lng(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*lng(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*lng(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*lng(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*lng(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*lng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)

               difgx2 = ave( tg(ix ,iy,igsp)/nuix(ix ,iy,igsp),
     .                       tg(ix2,iy,igsp)/nuix(ix2,iy,igsp) )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))

               flngxy(ix,iy,igsp) = difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                             (lng(ix2,iy,igsp) - lng(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)

c...  Now flux limit with flalfgxy
               t0 = max(tg(ix,iy,igsp),temin*ev)
               t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn = sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               qfl = flalfgxya(ix,igsp)*sx(ix,iy)* (vtn + vtnp)*rt8opi/8
ccc   MER NOTE:  no xy flux limit for ix=0 or ix=nx in original code
               if (isxyfl) then
                  fngxy(ix,iy,igsp) = fngxy(ix,iy,igsp) /
     .                           sqrt( 1 + (fngxy(ix,iy,igsp)/qfl)**2 )
               endif
 8904       continue
 8905    continue

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Flux-limit the total poloidal flux here
            do iy = j4, j8
              do ix = i1, i5
                flngx(ix,iy,igsp)= flngx(ix,iy,igsp)-flngxy(ix,iy,igsp)
              enddo
            enddo

      endif
c...  Finished with nonorthogonal mesh part

c ... Calculate the neutral flow velocity from v = flux/ng; these are
c ... diagnostic if isupgon=0, but used for uu and vy of the inertial
c ... gas if isupgon=1.  However, even when isupgon=1, the particle
c ... fluxes are fngx -> fnix and fngy -> fngy in pandf, i.e., methg
c ... determines the differencing for the inertial particle fluxes, not 
c ... methn
      do iy = j1, j5
         do ix = i1,i5
            ix1 = ixp1(ix,iy)
            uug(ix,iy,igsp) = flngx(ix,iy,igsp) / sx(ix,iy)
            vyg(ix,iy,igsp) = flngy(ix,iy,igsp) / sy(ix,iy)
c --------------- transfer inertial gas velocities to neutral ion species
            if (isupgon(igsp).eq.1) then
               vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
cc               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            end if
         enddo
cfw ----   If doing only the outer half we want this boundary condition:
         if (iy.le.iysptrx1(1).and.isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp) = 0
      enddo

c **- loop for uu just as in the previous version - needed for correct Jac?
      if (isupgon(igsp) .eq. 1) then
         do iy = j4, j6
            do ix = i1, i6
               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            enddo
         enddo
      endif

c.... Calculate the particle flux, fnix,y, from flux of lng, i.e., flngx,y
      do iy = j4, j8
         do ix = i1, i5
            ix2 = ixp1(ix,iy)
            fngx(ix,iy,igsp) = flngx(ix,iy,igsp) *
     .                     exp(0.5*(lng(ix,iy,igsp)+lng(ix2,iy,igsp)))
         enddo
      enddo
c ...   now do fniy
      do iy = j1, j5    # same loop ranges as for fngy in fd2tra
         do ix = i4, i8
            fngy(ix,iy,igsp) = flngy(ix,iy,igsp) *
     .                     exp(0.5*(lng(ix,iy,igsp)+lng(ix,iy+1,igsp)))
         enddo
      enddo

c.... Calculate the residual for the gas equation for diffusive neutral case

      if (isupgon(igsp).eq.0) then
         do 892 iy = j2, j5
            if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
	      fngx(nxc-1,iy,igsp) = 0.
	      fngx(nxc,  iy,igsp) = 0.
	      fngx(nxc+1,iy,igsp) = 0.
	    endif
            if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
            if (nxpt==2.and.ixmxbcl==1) fngx(ixrb(1)+1,iy,igsp)=0.
                   do 891 ix = i2, i5
               ix1 = ixm1(ix,iy)
               resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .                         psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .                       + volpsorg(ix,iy,igsp)                     
     .                       - fngx(ix,iy,igsp) + fngx(ix1,iy  ,igsp)
     .             - fluxfacy*(fngy(ix,iy,igsp) - fngy(ix ,iy-1,igsp))
     .                       + psgov_use(ix,iy,igsp)*vol(ix,iy)

               if (igsp.eq.1 .and. ishymol.eq.1) resng(ix,iy,igsp) =
     .                                  resng(ix,iy,igsp)+psordis(ix,iy)
 891        continue
 892     continue
      endif

 895  continue   # end of igsp loop from the beginning of subroutine

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
	   lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                         nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
	   lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
	      lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
	      lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif

      return

      end
c
c --------------------------------------------------------------------------
c END subroutine neudifl
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------
      subroutine neudifo

c ..  Older version of neudif where the gas velocities are deduced from
c ..  the gas fluxes and then used to form fnix if isupgon(igsp)=1

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1
      integer iy1, methgx, methgy, iy2, jx

      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx
      Use(Share)    # nxpt,geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methg,
                    # qfl,csh,qsh,cs,
                    # isupgon,iigsp,nlimgx,nlimgy
       # Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
      Use(Coefeq)   # cngfx,cngfy,cngmom,cmwall,cdifg,rld2dxg,rld2dyg
      Use(Bcond)    # albedoo,albedoi
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    # xlinc,xrinc,yinc,ixm1,ixp1,stretcx
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # fngx,fngy,fngxy,fnix,fniy
      Use(Conduc)   # nuiz, nucx, nuix
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
	  
*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10

      do 895 igsp = 1, ngsp

c.... First the flux in the x-direction

      do 888 iy = j4, j8
         do 887 ix = i1, i5
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgxa(ix,igsp) * sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
c            csh = sx(ix,iy) * gxf(ix,iy) *
c     .                 (stretcx(ix,iy)*vtn**2+stretcx(ix2,iy)*vtnp**2)/
c     .                 (nuix(ix,iy,igsp)+nuix(ix2,iy,igsp))
            csh = (1-isgasdc) * cdifg(igsp)* sx(ix,iy) * gxf(ix,iy) *
     .                ave( stretcx(ix,iy)*vtn**2/nuix(ix,iy,igsp),
     .                     stretcx(ix2,iy)*vtnp**2/nuix(ix2,iy,igsp) ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng +
     .                       rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                          0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
            qtgf = cngfx(igsp) * fgtdx(ix) * sx(ix,iy) *
     .            ave( stretcx(ix,iy)*gx(ix,iy)/nuix(ix,iy,igsp) ,
     .                 stretcx(ix2,iy)*gx(ix2,iy)/nuix(ix2,iy,igsp) )
     .                     * (vtn**2 - vtnp**2)
            vygtan(ix,iy,igsp) = 0.
            if (isnonog .eq. 1 .and. iy .le. ny) then
               if (isintlog .eq. 0) then
                  grdnv =( fym (ix,iy,1)*tg(ix2,iy1,igsp) +
     .                     fy0 (ix,iy,1)*tg(ix2,iy ,igsp) +
     .                     fyp (ix,iy,1)*tg(ix2,iy2,igsp) +
     .                     fymx(ix,iy,1)*tg(ix ,iy1,igsp) +
     .                     fypx(ix,iy,1)*tg(ix, iy2,igsp) -
     .                     fym (ix,iy,0)*tg(ix ,iy1,igsp) -
     .                     fy0 (ix,iy,0)*tg(ix ,iy ,igsp) -
     .                     fyp (ix,iy,0)*tg(ix ,iy2,igsp) - 
     .                     fymx(ix,iy,0)*tg(ix4,iy1,igsp) -
     .                     fypx(ix,iy,0)*tg(ix6,iy2,igsp) )/dxnog(ix,iy)
               elseif (isintlog .eq. 1) then
                  grdnv =( exp( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                    -exp( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               endif
               vygtan(ix,iy,igsp) = ( cngfx(igsp) / (mg(igsp)*0.5*
     .                         (nuix(ix,iy,igsp)+nuix(ix2,iy,igsp))) ) *
     .                             ( grdnv/cos(angfx(ix,iy)) - 
     .                             (tg(ix2,iy,igsp) - tg(ix,iy,igsp))
     .                                                 * gxf(ix,iy) ) 
               if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
                  vygtan(ix,iy,igsp) = 0.
               endif
               if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
                  vygtan(ix,iy,igsp) = 0.
               endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - (1-isupgon(igsp))*vygtan(ix,iy,igsp)*sx(ix,iy)
            nconv = 2.0*(ng(ix,iy,igsp)*ng(ix2,iy,igsp)) /
     .                  (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
c...   Use upwind for "convective" grad T term if methgx .ne. 2
            if(methgx.ne.2) nconv =
     .                         ng(ix ,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ng(ix2,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ng(ix,iy,igsp)-ng(ix2,iy,igsp)) + qtgf * nconv
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  boundary face - shouldnt matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .             (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floxg
            floxg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = 
     .             floxg(ix,iy) + cngflox(1)*sx(ix,iy)*uu(ix,iy,1)

  887    continue
         conxg(nx+1,iy) = 0
  888 continue

c.... Now the flux in the y-direction

      do 890 iy = j1, j5
         do 889 ix = i4, i8
	    t0 = max(tg(ix,iy,igsp),temin*ev)
	    t1 = max(tg(ix,iy+1,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgya(iy,igsp) * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .                              ( ngy0(ix,iy,igsp)*gy(ix,iy) + 
     .                                ngy1(ix,iy,igsp)*gy(ix,iy+1) ) / 
     .                                     (8*(gy(ix,iy)+gy(ix,iy+1)))
            csh = (1-isgasdc) * (cdifg(igsp) *sy(ix,iy) /dynog(ix,iy)) *
     .                            ave( vtn**2/nuix(ix,iy,igsp) ,
     .                                 vtnp**2/nuix(ix,iy+1,igsp) ) +
     .            isgasdc * sy(ix,iy) * difcng / dynog(ix,iy) +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*dynog(ix,iy)*
     .                       0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))
c               csh = sy(ix,iy) * ( ((vtn**2+vtnp**2)/ dynog(ix,iy)) /
c     .                 (nuix(ix,iy,igsp)+nuix(ix,iy+1,igsp)) )
            qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                     ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                          gy(ix,iy+1)/nuix(ix,iy+1,igsp) )
     .                    * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) )
              endif
              qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                      ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                           gy(ix,iy+1)/nuix(ix,iy+1,igsp) ) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ngy0(ix,iy,igsp)-ngy1(ix,iy,igsp)) + qtgf*nconv
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floyg
	    floyg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)  
c...  now add the convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) = 
     .               floyg(ix,iy)+cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)

  889    continue
  890 continue

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             ng(0:nx+1,0:ny+1,igsp),fngx(0:nx+1,0:ny+1,igsp),
     .             fngy(0:nx+1,0:ny+1,igsp),0,methg)

c ... Calculate the neutral flow velocity from v = flux/ng
      do 8903 iy = j1, j5
         do 8902 ix = i1,i5
            ix1 = ixp1(ix,iy)
            uug(ix,iy,igsp) = fngx(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                      *sx(ix,iy) )
            vyg(ix,iy,igsp) = fngy(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
     .                                                      *sy(ix,iy) )
            if (isupgon(igsp).eq.1) then
c --------------- We need to transfer the diffusive radial neutral
c --------------- velocity to the "ion" species containing the neutrals
               vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
            end if
 8902    continue
cfw ----   If doing only the outer half we want this boundary condition:
         if (iy.le.iysptrx1(1).and.isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp) = 0
 8903 continue

c ... For nonorthogonal mesh and diffusive neutrals, limit fngy for pump out
      if (isnonog.eq.1 .and. isupgon(igsp).eq.0) then
         do iy = j1, j5
            do ix = i4, i8
               if (nlimgy*fngy(ix,iy,igsp)*
     .                 (ng(ix,iy,igsp)-ng(ix,iy+1,igsp)) .lt. 0) then
                  fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/( 1 - 2*nlimgy +
     .                      nlimgy*(ng(ix,iy+1,igsp)/ng(ix,iy,igsp)+
     .                              ng(ix,iy,igsp)/ng(ix,iy+1,igsp)) )
               endif
            enddo
         enddo     
      endif            

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then

         do 8905 iy = j1, min(j6, ny)
            iy1 = max(iy-1,0)
            do 8904 ix = i1, min(i6, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
               if (methgx .eq. 6) then  # log interpolation
               grdnv =( exp(fym (ix,iy,1)*log(ng(ix2,iy1 ,igsp)) + 
     .                      fy0 (ix,iy,1)*log(ng(ix2,iy  ,igsp)) +
     .                      fyp (ix,iy,1)*log(ng(ix2,iy+1,igsp)) + 
     .                      fymx(ix,iy,1)*log(ng(ix ,iy1 ,igsp)) +
     .                      fypx(ix,iy,1)*log(ng(ix, iy+1,igsp)))
     .                - exp(fym (ix,iy,0)*log(ng(ix ,iy1 ,igsp)) +
     .                      fy0 (ix,iy,0)*log(ng(ix ,iy  ,igsp)) +
     .                      fyp (ix,iy,0)*log(ng(ix ,iy+1,igsp)) +
     .                      fymx(ix,iy,0)*log(ng(ix4,iy1 ,igsp)) + 
     .                      fypx(ix,iy,0)*log(ng(ix6,iy+1,igsp))) ) / 
     .                                                  dxnog(ix,iy)
               elseif (methgx .eq. 7) then  # inverse interpolation
               grdnv =( 1/(fym (ix,iy,1)/ng(ix2,iy1 ,igsp) + 
     .                     fy0 (ix,iy,1)/ng(ix2,iy  ,igsp) +
     .                     fyp (ix,iy,1)/ng(ix2,iy+1,igsp) + 
     .                     fymx(ix,iy,1)/ng(ix ,iy1 ,igsp) +
     .                     fypx(ix,iy,1)/ng(ix, iy+1,igsp))
     .                - 1/(fym (ix,iy,0)/ng(ix ,iy1 ,igsp) +
     .                     fy0 (ix,iy,0)/ng(ix ,iy  ,igsp) +
     .                     fyp (ix,iy,0)/ng(ix ,iy+1,igsp) +
     .                     fymx(ix,iy,0)/ng(ix4,iy1 ,igsp) + 
     .                     fypx(ix,iy,0)/ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               else                   # linear interpolation
               grdnv =( (fym (ix,iy,1)*ng(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*ng(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*ng(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*ng(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*ng(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*ng(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*ng(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*ng(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*ng(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               endif
               difgx2 = ave( tg(ix ,iy,igsp)/nuix(ix ,iy,igsp),
     .                       tg(ix2,iy,igsp)/nuix(ix2,iy,igsp) )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
               fngxy(ix,iy,igsp) = difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                             (ng(ix2,iy,igsp) - ng(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
 8904       continue
 8905    continue

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Dont bother if  solving the neutral momentum equation
         if (isupgon(igsp).eq.0) then
            do iy = j4, j8
               do ix = i1, i5
                  ix2 = ixp1(ix,iy)
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp) - fngxy(ix,iy,igsp)
c ...          adjust fluxes to prevent pump out
                  if (nlimgx*fngx(ix,iy,igsp)
     .                   *(ng(ix,iy,igsp)-ng(ix2,iy,igsp)) .lt. 0.) then
                    fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/( 1 - 2*nlimgx +
     .                          nlimgx*(ng(ix2,iy,igsp)/ng(ix,iy,igsp) +
     .                                  ng(ix,iy,igsp)/ng(ix2,iy,igsp)) )
                  endif
               enddo
            enddo
         endif

      endif
c...  Finished with nonorthogonal mesh part

c.... Calculate the residual or right-hand-side for the gas equation

      if (isupgon(igsp).eq.0) then
         do 892 iy = j2, j5
            if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
	      fngx(nxc-1,iy,igsp) = 0.
	      fngx(nxc,  iy,igsp) = 0.
	      fngx(nxc+1,iy,igsp) = 0.
	    endif
            if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
            if (nxpt==2.and.ixmxbcl==1) fngx(ixrb(1)+1,iy,igsp)=0.
                   do 891 ix = i2, i5
               ix1 = ixm1(ix,iy)
               resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .                         psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .                       + volpsorg(ix,iy,igsp)                     
     .                       - fngx(ix,iy,igsp) + fngx(ix1,iy  ,igsp)
     .                       - fngy(ix,iy,igsp) + fngy(ix ,iy-1,igsp)
     .                       + psgov_use(ix,iy,igsp)*vol(ix,iy)
               if (igsp.eq.1 .and. ishymol.eq.1) resng(ix,iy,igsp) =
     .                                  resng(ix,iy,igsp)+psordis(ix,iy)
 891        continue
 892     continue
      else if (isupgon(igsp).eq.1) then

c --- Form the poloidal velocity uu(,,2) from
c --- a) the projection of the neutral parallel velocity up(,,2),
c --- b) the projection of the 2-direction gas velocity and
c --- c) vygtan, the grad(T) part of the non-orth diffusive radial velocity.
c --- d) fngxy, the grad(n) part of the non-orth diffusive radial velocity.
c --- The fngxy contribution could have been kept separate and added to
c --- fnix in PARBAL, but we include it here so that it automatically gets
c --- taken into account in PARBAL and MOMBAL_B2.
c --- By multiplying uu(,,2) with sx*ng in PARBAL we get the 
c --- TOTAL neutral particle flux out of the poloidal face.
c --- By multiplying uu(,,2) with sx*ng*mi(1)*up(,,iigsp) in MOMBAL_B2
c --- we get the TOTAL parallel momentum flux out of the poloidal face.
c --- Note that rrv=Bpol/B is defined at a velocity point
c --- In order to get Bt/B at a VELOCITY point we cannot use rbfbt,
c --- use (Bt/B)**2=1-rrv**2 instead.
         do 23 iy = j4, j6
            do 22 ix = i1, i6
               ix2 = ixp1(ix,iy)
               uu(ix,iy,iigsp) = rrv(ix,iy)*up(ix,iy,iigsp) +
     .              (1.-rrv(ix,iy)*rrv(ix,iy))*uug(ix,iy,igsp)
               if (isnonog .eq. 1) uu(ix,iy,iigsp) =
     .              uu(ix,iy,iigsp) - vygtan(ix,iy,igsp) -
     .              fngxy(ix,iy,igsp) / (sx(ix,iy)*
     .                    0.5*(ng(ix,iy,igsp)+ng(ix2,iy,igsp)))
 22         continue
 23      continue
      end if

 895  continue

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                          nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif

      return

      end
c
c --------------------------------------------------------------------------
c END subroutine neudifo - old version of neudif
c --------------------------------------------------------------------------


c --------------------------------------------------------------------------
c SUBROUTINE TO SET UP ENERGY EQUATION FOR NEUTRAL GAS
c --------------------------------------------------------------------------

      subroutine engbalg

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,nu1,nu2,ffyi,ffyo
      real vt0,vt1,wallfac,lxtgc,dupdx,dupdy,fniy_recy,thetacc
      real vttn,vttp
      integer ifld,iixt,iy1, methgx, methgy, iy2, jx
      #Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t0,t1,t2,tv,a
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1
      Use(Share)    # geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methg,
                    # qfl,csh,qsh,cs,
                    # isupgon,iigsp,nlimgx,nlimgy,nlimiy,rld2dx,rld2dy

      Use(Coefeq)   # cngfx,cngfy,cngmom,cmwall,cdifg,rld2dxg,rld2dyg
      Use(Bcond)    # albedoo,albedoi
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    # xlinc,xrinc,yinc,ixm1,ixp1,stretcx
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # fngx,fngy,fngxy,fnix,fniy
      Use(Conduc)   # nuiz, nucx, nuix
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
      Use(Wkspace)  # w0,w1,etc
      Use(MCN_dim)  #
      Use(MCN_sources)   # cfneutsor_ei

*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

*  -- Compute sources terms v_grad_Pg; first initialize array over range
c  -- This v_grad_Pg term first added by MZhao

      do igsp = 1, ngsp
        do iy = j2, j5
          do ix = i2, i5
            segc(ix,iy,igsp) = 0.0
          enddo
        enddo
      enddo

      do igsp = 1, ngsp
        if(istgon(igsp) == 1) then 
          do iy = j2, j5
            do ix = i2, i5
              ix1 = ixm1(ix,iy)
              ix2 = ixp1(ix,iy)
              iy1 = max(0,iy-1)
              tv = (pg(ix2,iy,igsp) - pg(ix ,iy,igsp))
              t1 = (pg(ix ,iy,igsp) - pg(ix1,iy,igsp))
              segc(ix,iy,igsp) = 0.5*cvgpg*( 
     ,                 uuxg(ix, iy,igsp)*ave(gx(ix2,iy),gx(ix, iy))*tv +
     .                 uuxg(ix1,iy,igsp)*ave(gx(ix ,iy),gx(ix1,iy))*t1 )*
     .                                                      vol(ix,iy)
              t2 = cvgpg*0.5*( vyg(ix, iy,igsp)*dynog(ix, iy)*
     .                       (pgy1(ix, iy,igsp)-pgy0(ix, iy,igsp)) +
     .                   vyg(ix,iy1,igsp)*dynog(ix,iy1)*
     .                       (pgy1(ix,iy1,igsp)-pgy0(ix,iy1,igsp)) )
              segc(ix,iy,igsp)=segc(ix,iy,igsp) + cvgpg*t2*vol(ix,iy)
            enddo
          enddo
        endif
      enddo



*  -- Compute flux terms; initialize some arrays to 0 --

      do igsp = 1, ngsp
        do iy = j1, j6
          do ix = i1, i6
            floxge(ix,iy,igsp) = 0.0e0
            floyge(ix,iy,igsp) = 0.0e0
            conxge(ix,iy,igsp) = 0.0e0
            conyge(ix,iy,igsp) = 0.0e0
          enddo
        enddo
      enddo

*  ---------------------------------------------------------------------
*  Compute thermal conductances
*  ---------------------------------------------------------------------
c ... Compute poloidal conduction
      do igsp = 1,ngsp
        do iy = j4, j8
          do ix = i1, i5
            ix2 = ixp1(ix,iy)

            t0 = max (tg(ix,iy,igsp), temin*ev)
            t1 = max (tg(ix2,iy,igsp), temin*ev)
            vt0 = sqrt(t0/mg(igsp))
            vt1 = sqrt(t1/mg(igsp))
c... flux-limit occurs in building hcxg - do not flux-limit 2nd time
            conxge(ix,iy,igsp) = sx(ix,iy) * hcxg(ix,iy,igsp) * gxf(ix,iy) 
          enddo
          conxge(nx+1,iy,igsp) = 0
        enddo
      enddo

*  -- compute radial conduction conyge
      do igsp = 1, ngsp
        do iy = j1, j5
          do ix = i4, i8
            conyge(ix,iy,igsp) = sy(ix,iy)*hcyg(ix,iy,igsp)/dynog(ix,iy)
          enddo
        enddo
      enddo

      do igsp = 1, ngsp
        do ix = i1, i6
          conyge(ix,ny+1,igsp) = 0.0e0
        enddo
      enddo
*  ---------------------------------------------------------------------
*  compute convective flow of tg
*  ---------------------------------------------------------------------

*  -- compute floxge --

      do igsp = 1, ngsp
        do iy = j4, j8
          do ix = i1, i5
            floxge(ix,iy,igsp) = cfcvtg*2.5*fngx(ix,iy,igsp)
          enddo
          floxge(nx+1,iy,igsp) = 0.
        enddo
      enddo

*  -- Correct bdry:remove any inward power from plates; ok in parallel
      do igsp = 1, ngsp
        do iy = j4, j8
          do jx = 1, nxpt  #if at plate, sub (1-cfloxiplt)*neut-contrib
            if(ixmnbcl==1) then  #real div plt -need for parallel UEDGE
              iixt = ixlb(jx) #left plate
              if(fngx(iixt,iy,igsp) > 0.) then
                floxge(iixt,iy,igsp) = floxge(iixt,iy,igsp) -
     .                    (1.-cfloxiplt)*cfcvti*2.5*fngx(iixt,iy,igsp)
              endif
            endif
            if(ixmxbcl==1) then #real div plt -need for parallel UEDGE
              iixt = ixrb(jx) # right plate
              if(fngx(iixt,iy,igsp) < 0.) then
                floxge(iixt,iy,igsp) = floxge(iixt,iy,igsp) -
     .                   (1.-cfloxiplt)*cfcvti*2.5*fngx(iixt,iy,igsp)
              endif
              floxge(ixrb(jx)+1,iy,igsp) = 0.0e0 #cosmetic
            endif
          enddo
        enddo
      enddo

*  -- compute floyge --

      do igsp = 1, ngsp
        do iy = j1, j5
          do ix = i4, i8
            floyge(ix,iy,igsp) = cfcvtg*2.5*fngy(ix,iy,igsp)
          enddo
        enddo
      enddo
 
*  -- Correct bdry:remove any inward power from walls; ok in parallel
      do igsp = 1, ngsp
        do ix = i4, i8
          do jx = 1, nxpt  #if on PF wall, sub (1-cfloygw)*neut-contrib
            if(iymnbcl==1) then  #real PFw-need for parallel UEDGE?
              if(ix <= ixpt1(jx) .or. ix > ixpt2(jx)) then
                if(fngy(ix,0,igsp) > 0.) then
                  floyge(ix,0,igsp) = floyge(ix,0,igsp) -
     .                       (1.-cfloygwall)*cfcvtg*2.5*fngy(ix,0,igsp)
                endif
              endif
            endif
          enddo
          if(iymxbcl==1) then #real outer-w-need for parallel UEDGE?
            if(fngy(ix,ny,igsp) < 0.) then
              floyge(ix,ny,igsp) = floyge(ix,ny,igsp) -
     .                    (1.-cfloygwall)*cfcvtg*2.5*fngy(ix,ny,igsp)
            endif
          endif
          floyge(ix,ny+1,igsp) = 0.0e0 #cosmetic
        enddo   #ix loop
      enddo

*  -- Combine conduction/convection to compute thermal energy flow --
      do igsp = 1,ngsp
        if(istgon(igsp) == 1) then
          call fd2tra (nx,ny,floxge(0:nx+1,0:ny+1,igsp), 
     .          floyge(0:nx+1,0:ny+1,igsp), conxge(0:nx+1,0:ny+1,igsp),
     .          conyge(0:nx+1,0:ny+1,igsp),tg(0:nx+1,0:ny+1,igsp),
     .          fegx(0:nx+1,0:ny+1,igsp),fegy(0:nx+1,0:ny+1,igsp),
     .          0,methi)
        endif
      enddo

c...  Add y-component of nonorthogonal diffusive flux; convective component 
c...  already added to uug(ix,iy,igsp)
      if (isnonog == 1) then
        do igsp = 1, ngsp
          if(istgon(igsp) == 0) cycle
          do iy = j1, j6
            if (iy .gt. ny) cycle
            iy1 = max(iy-1,0)
            do ix = i1, i6
              ix1 = ixm1(ix,iy)
              ix2 = ixp1(ix,iy)
              ix3 = ixm1(ix,iy1)
              ix4 = ixp1(ix,iy1)
              ix5 = ixm1(ix,iy+1)
              ix6 = ixp1(ix,iy+1)
              t0 = max(tg(ix,iy,igsp),tgmin*ev) 
              t1 = max(tg(ix2,iy,igsp),tgmin*ev)
              vtn = sqrt( t0/mg(igsp) )
              vtnp = sqrt( t1/mg(igsp) )
              nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
              nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)

c --- Note: this four-point average results in not getting the full Jac. for
c --- a nonorthogonal mesh because of ngy1,0 - see def. of hcyn

                grdnv =( ( fym (ix,iy,1)*log(tg(ix2,iy1 ,igsp)) +  
     .                     fy0 (ix,iy,1)*log(tg(ix2,iy  ,igsp)) +
     .                     fyp (ix,iy,1)*log(tg(ix2,iy+1,igsp)) +  
     .                     fymx(ix,iy,1)*log(tg(ix ,iy1 ,igsp)) +
     .                     fypx(ix,iy,1)*log(tg(ix ,iy+1,igsp)) ) 
     .                  -( fym (ix,iy,0)*log(tg(ix ,iy1 ,igsp)) +
     .                     fy0 (ix,iy,0)*log(tg(ix ,iy  ,igsp)) +
     .                     fyp (ix,iy,0)*log(tg(ix ,iy+1,igsp)) +
     .                     fymx(ix,iy,0)*log(tg(ix4,iy1 ,igsp)) +  
     .                     fypx(ix,iy,0)*log(tg(ix6,iy+1,igsp)) ) ) / 
     .                                                  dxnog(ix,iy)  
               difgx2 = ave( tg(ix ,iy,igsp)/nu1,
     .                       tg(ix2,iy,igsp)/nu2 )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))

               fegxy(ix,iy,igsp) = cfegxy*exp( 0.5*
     .                 (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                      difgx2*ave(ng(ix2,iy,igsp),ng(ix,iy,igsp))*
     .                                 ( grdnv/cos(angfx(ix,iy))
     .                  - (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))*
     .                                        gxf(ix,iy) )*sx(ix,iy)
c...  Flux limit with flalftxt even though hcys have parallel FL built in
               t0 = max(tg(ix,iy,igsp),tgmin*ev)
               t1 = max(tg(ix2,iy,igsp),tgmin*ev)
               vttn = t0*sqrt( t0/mg(igsp) )
               vttp = t1*sqrt( t1/mg(igsp) )
       if(isfegxyqflave == 0) then
               qfl = flalftgxy(igsp)*0.25*sx(ix,iy) * (vttn+vttp) * 
     .                              (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
       else  #use harmonic average of T*vt and ng to face
               qfl = flalftgxy(igsp) * sx(ix,iy) * ave(vttn,vttp) * 
     .                            ave(ng(ix,iy,igsp),ng(ix2,iy,igsp))
       endif
               fegxy(ix,iy,igsp) = fegxy(ix,iy,igsp) /
     .                          sqrt(1. + (fegxy(ix,iy,igsp)/qfl)**2)
               fegx(ix,iy,igsp) = fegx(ix,iy,igsp) - fegxy(ix,iy,igsp)
             enddo  #loop over ix
           enddo    #loop over iy
        enddo       #loop over igsp
      endif         #if-test on isnonog

*  ---------------------------------------------------------------------
*  compute the energy residuals.
*  ---------------------------------------------------------------------

*  -- total energy residual and equipartition --

      do igsp = 1, ngsp
        do iy = j2, j5
          iy1 = max(0,iy-1)
          do ix = i2, i5
            ix1 = ixm1(ix,iy)
            reseg(ix,iy,igsp)= -( fegx(ix,iy,igsp)-fegx(ix1,iy,  igsp)+
     .                            fegy(ix,iy,igsp)-fegy(ix, iy1,igsp) )
     .                                                + segc(ix,iy,igsp)
            reseg(ix,iy,igsp)= reseg(ix,iy,igsp) + vol(ix,iy)* 
     .                      eqpg(ix,iy,igsp)*(ti(ix,iy)-tg(ix,iy,igsp))#+
#     .                   cftgdiss(igsp)*psorg(ix,iy,igsp)*tg(ix,iy,igsp)
            if (igsp.eq.1) then  #..for D0, we should include D+ and D0 in Ti
              seic(ix,iy) = seic(ix,iy)- vol(ix,iy)*(1.0-cftiexclg)*
     .                                               eqpg(ix,iy,igsp)*
     .                                      (ti(ix,iy)-tg(ix,iy,igsp))
            else
              seic(ix,iy) = seic(ix,iy)- vol(ix,iy)*
     .                                               eqpg(ix,iy,igsp)*
     .                                      (ti(ix,iy)-tg(ix,iy,igsp))
              reseg(ix,iy,igsp) = reseg(ix,iy,igsp)
     .                                 + cftgeqp*ng(ix,iy,igsp)*
     .                           (1.0-cftiexclg)*ng(ix,iy,1)*kelighg(igsp)*
     .                        (tg(ix,iy,1)-tg(ix,iy,igsp))*vol(ix,iy)
              reseg(ix,iy,1) = reseg(ix,iy,1) - cftgeqp*ng(ix,iy,igsp)*
     .                                           ng(ix,iy,1)*kelighg(igsp)*
     .                        (tg(ix,iy,1)-tg(ix,iy,igsp))*vol(ix,iy)
              if (ishymol.eq.1 .and. igsp.eq.2) then  #..D2 dissociation
                 reseg(ix,iy,igsp) =
     .                             reseg(ix,iy,igsp)+psorg(ix,iy,igsp)
     .                                               *1.5*tg(ix,iy,igsp)
                 t0 = cfnidhmol*0.25*(uuxg(ix,iy,igsp)+uuxg(ix1,iy,igsp))
     .                              *(uuxg(ix,iy,igsp)+uuxg(ix1,iy,igsp))
                 t1 = cfnidhmol*0.25*(vyg(ix,iy,igsp)+vyg(ix1,iy,igsp))
     .                              *(vyg(ix,iy,igsp)+vyg(ix1,iy,igsp))
                 t2 = 0. #.. molecule v in the tol direction, it seems to be assumed as 0 in neudifpg?
                 reseg(ix,iy,1) = reseg(ix,iy,1) + cfnidhdis*
     .                            0.5*mg(1)*(t0+t1+t2)*psordis(ix,iy)
                 seic(ix,iy) = seic(ix,iy) + cftiexclg*cfnidhdis*
     .                            0.5*mg(1)*(t0+t1+t2)*psordis(ix,iy)
                 t0 = cfnidhmol*0.25*(uuxg(ix,iy,igsp)+uuxg(ix1,iy,igsp))
     .                              *(uuxg(ix,iy,1)+uuxg(ix1,iy,1))
                 t1 = cfnidhmol*0.25*(vyg(ix,iy,igsp)+vyg(ix1,iy,igsp))
     .                              *(vyg(ix,iy,1)+vyg(ix1,iy,1))
                 t2 = 0.
                 reseg(ix,iy,1) = reseg(ix,iy,1) - cfnidhdis*
     .                                mg(1)*(t0+t1+t2)*psordis(ix,iy)
                 seic(ix,iy) = seic(ix,iy) - cftiexclg*cfnidhdis*
     .                                mg(1)*(t0+t1+t2)*psordis(ix,iy)
              endif
            endif
	    #..zml place holder for neutral-neutral collision,
	    #..    not included above?
          enddo
        enddo
      enddo

      if((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
        do igsp = 1, ngsp
          do iy = j2, j5
            fegx(nxc-1,iy,igsp) = 0. 
            fegx(nxc  ,iy,igsp) = 0. 
            fegx(nxc+1,iy,igsp) = 0. 
          enddo
        enddo
      endif

*  -- Energy transfer to impurity neutrals at tg(,,igsp)
      if (ngsp >= 2) then   # for now, specialized to igsp=2 only
        do ifld = nhsp+1, nisp
          do iy = j2, j5    # iys,iyf limits dont seem to work(?)
            do ix = i2, i5
              #      possible bugs here? resei is replaced with seic
	      #      since resei is not defined before this subroutine called
	      #      more needs to be done here.. e.g. for segc
              seic(ix,iy) =seic(ix,iy) -cftiimpg*1.5*ni(ix,iy,ifld)*
     .                      (nucxi(ix,iy,ifld)+nueli(ix,iy,ifld))*
     .                      (ti(ix,iy) - tg(ix,iy,2))*vol(ix,iy)
              #..zml place holder for things remain to be done for segc
	      #..    to include neutral-impurity ion collisions
            enddo
          enddo
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
c  END subroutine engbalg - THE NEUTRAL GAS ENERGY EQUATION
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
      subroutine pandf1(xc, yc, ieq, neq, time, yl, yldot)

c ... Calculates matrix A and the right-hand side depending on the 
c     values of xc, yc.
c  Definitions for argument list
c
c  Input variables:
c    xc is poloidal index of perturbed variablefor Jacobian calc, 
c       or =-1 for full RHS evaluation
c    yc is radial index for perturbed variable for Jacobian calc, 
c       or =-1 for full RHS evaluation
c    ieq is the eqn number for Jacobian eval; not presently used
c    neq is the total number of variables
c    time is the present physical time; useable by VODPK but not NKSOL
c    yl is the vector of unknowns
c  Output variables:
c    yldot is the RHS of ODE solver or RHS=0 for Newton solver (NKSOL)

      implicit none
      Use(Dim)     # nusp,nisp,ngsp
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      Use(UEpar)   # svrpkg,isbcwdt,isnionxy,isuponxy,isteonxy,istionxy,
                   # isngonxy,isphionxy
cc      Use(Selec)   # i2,i5,j2,j5
      Use(Time_dep_nwt)   # nufak,dtreal,ylodt,dtuse
      Use(Indexes) # idxn,idxg,idxu,dxti,idxte,idxphi
      Use(Ynorm)   # isflxvar,isrscalf
      Use(Share)    # geometry,nxc,isnonog,cutlo
      Use(Indices_domain_dcl) # ixmnbcl,ixmxbcl,iymnbcl,iymxbcl
      Use(Compla)  # zi
      Use(Xpoint_indices)      # ixpt1,ixpt2,iysptrx

*  -- arguments
      integer xc, yc, ieq, neq     # ieq is the equation index for Jac. calc
      real time, yl(neqmx),yldot(neq)

*  -- local variables
      integer ix,iy,igsp,iv,iv1,ifld,j2l,j5l,i2l,i5l
      character*80 msgjm
      integer nrcv, ierrjm, ijmgetmr

ccc      save

c
c     Check if "k" or "kaboom" has been typed to jump back to the parser
c
      if (((svrpkg.eq.'nksol') .or. (svrpkg.eq.'petsc')) .and. iskaboom.eq.1) then
                              #can only call once - preserves 's' in vodpk
        ierrjm = ijmgetmr(msgjm,80,1,nrcv)
        if (ierrjm .eq. 0) then
          if (msgjm(1:nrcv).eq.'kaboom' .or. msgjm(1:nrcv).eq.'k')then
            call xerrab("")
          endif
        endif
      endif

c     check if a "ctrl-c" has been type to interrupt - from basis
      call ruthere

c
c  PANDF calculates the equations in the interior of the grid, plus calls
c  bouncon for B.C. and poten for potential
c
      call pandf (xc, yc, neq, time, yl, yldot)
c
c...  If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and the ODEs need
c...  to be modified as original equations are for d(nv)/dt, etc
c...  If isflxvar=2, variables are ni,v,nTe,nTi,ng. Boundary equations and
c...  potential equations are not reordered.

      if(isflxvar.ne.1 .and. isrscalf.eq.1) call rscalf(yl,yldot)
c
c ... Now add psuedo or real timestep for nksol method, but not both
      if (nufak.gt.1.e5 .and. dtreal.lt.1.e-5) then
         call xerrab('***Both 1/nufak and dtreal < 1.e5 - illegal***')
      endif

c...  Add a real timestep, dtreal, to the nksol equations 
c...  NOTE!! condition yl(neq+1).lt.0 means a call from nksol, not jac_calc

      if(dtreal < 1.e15) then
       if((svrpkg=='nksol' .and. yl(neq+1)<0) .or. svrpkg == 'petsc') then
         if (isbcwdt .eq. 0) then  # omit b.c. eqns
cccMER   NOTE: what about internal guard cells (for dnbot,dnull,limiter) ???
            j2l = 1
            j5l = ny
            i2l = 1
            i5l = nx
         else                      # include b.c. eqns
            j2l = (1-iymnbcl)
            j5l = ny+1-(1-iymxbcl)
            i2l = (1-ixmnbcl)
            i5l = nx+1-(1-ixmxbcl)
         endif           
         do iy = j2l, j5l    # if j2l=j2, etc., omit the boundary equations
            do ix = i2l, i5l
              do ifld = 1, nisp
                if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  yldot(iv) = (1.-fdtnixy(ix,iy,ifld))*yldot(iv)
                  if(zi(ifld).eq.0. .and. ineudif.eq.3) then
                    yldot(iv) = yldot(iv) - (1/n0(ifld))*
     .                          (exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                  else
                    yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
                endif
              enddo
               if(ix.ne.nx+2*isbcwdt) then  
                              # nx test - for algebr. eq. unless isbcwdt=1
                  do ifld = 1, nusp
                    if(isuponxy(ix,iy,ifld).eq.1) then
                      iv = idxu(ix,iy,ifld)
                      yldot(iv) = (1.-fdtupxy(ix,iy,ifld))*yldot(iv)
                      yldot(iv) = yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                    endif
                  enddo
               endif
               if (isteonxy(ix,iy) == 1) then
                 iv =  idxte(ix,iy)
                 yldot(iv) = (1.-fdttexy(ix,iy))*yldot(iv)
                 yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv) 
               endif
               if (istionxy(ix,iy) == 1) then
                 iv1 = idxti(ix,iy)
                 yldot(iv1) = (1.-fdttixy(ix,iy))*yldot(iv1)
                 yldot(iv1)=yldot(iv1) - (yl(iv1)-ylodt(iv1))/dtuse(iv1)
               endif
               do igsp = 1, ngsp
                  if(isngonxy(ix,iy,igsp).eq.1) then
                     iv = idxg(ix,iy,igsp)
                     yldot(iv) = (1.-fdtngxy(ix,iy,igsp))*yldot(iv)
                     if(ineudif.eq.3) then
                       yldot(iv) = yldot(iv) - (1/n0g(igsp))*
     .                            (exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                     else
                       yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                     endif
                  endif
               enddo
               do igsp = 1, ngsp
                  if(istgonxy(ix,iy,igsp).eq.1) then
                     iv = idxtg(ix,iy,igsp)
                     yldot(iv) = (1.-fdttgxy(ix,iy,igsp))*yldot(iv)
                     yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
               enddo
               if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
                  iv = idxphi(ix,iy)
                  yldot(iv) = (1.-fdtphixy(ix,iy))*yldot(iv)
                  yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv)
               endif

            enddo
         enddo
      
C...  Now do an additional relaxation of the potential equations with
c...  timestep dtphi
        if (dtphi < 1e10) then
          do iy = 0, ny+1
            do ix = 0, nx+1
              if (isphionxy(ix,iy) == 1) then
                iv = idxphi(ix,iy)
                yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtphi
              endif
            enddo
          enddo
        endif

       endif   #if-test on svrpkg and yl(neq+1)
      endif    #if-test on dtreal

      return
      end
c****** end of subroutine pandf1 ************
c-----------------------------------------------------------------------
      subroutine rscalf(yl, yldot)

c...  This routine reorders the equation for yldot depending on what
c...  variables are used, i.e., n,nv,nT, or n,v,T, or n,v,nT

      implicit none
      Use(Dim)       # nx,ny,nusp,nisp
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      Use(Ynorm)     # isflxvar,nnorm,n0
      Use(Compla)
      Use(Indexes)
      Use(Selec)     # i2,i5,j2,j5,ixp1

      Use(Coefeq)    # cngtgx
      Use(UEpar)     # isnion,isupon,isteon,istion,isngon,isnionxy,isuponxy,
                     # isteonxy,istionxy,isngonxy,isphionxy
      Use(Rhsides)   # resco
      Use(Comgeo)    # vol

*  -- Input parameters
      real yl(*), yldot(*)

*  -- Local variables
      integer ifld
      #Former Aux module variables
      integer ix,iy,iv,iv1,iv2,ix1,igsp
      real nbv, nbvdot, nbidot, nbedot, nbgdot, yldot_np1, nbg2dot(ngsp)

c...  If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and the ODEs need
c...  to be modified as original equations are for d(nv)/dt, etc
c...  If isflxvar=2, variables are ni,v,nTe,nTi,ng.  Boundary eqns and
c...  potential are not reordered (loops from i2-i5 and j2-j5).

*****************************************************************
*  --  ODE Equations to be solved - rescaling
*****************************************************************

      do 270 iy = j2, j5
         do 260 ix = i2, i5
            nbedot = 0.
            nbidot = 0.
            nbgdot = 0.
            nbvdot = 0.
ccc            if (isngonxy(ix,iy,1) .eq. 1) nbidot = cngtgx(1)*yldot(idxg(ix,iy,1))
            do 255 ifld = 1, nisp
	       if (isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  #.. separate ions and atoms
                  #nbidot = nbidot + yldot(iv)*n0(ifld)
                  if (isupgon(1)==1 .and. zi(ifld)==0) then  #neutral hyd
                    nbgdot = yldot(iv)*n0(ifld)
                  else
                    nbidot = nbidot + yldot(iv)*n0(ifld)
                  endif
                  nbedot = nbedot + zi(ifld)*yldot(iv)*n0(ifld)
               endif
  255       continue
            do igsp = 1, ngsp
               nbg2dot(igsp) = 0.
               if(isngonxy(ix,iy,igsp) == 1) then
                 iv = idxg(ix,iy,igsp)
                 nbg2dot(igsp) = yldot(iv)*n0g(igsp)
               endif
            enddo
c...    Omit cases where iseqalg=1 indicating algebraic b.c. and ix=nx for up.
c...    Could cause trouble for nisp.ne.nusp, but majority species should
c...    be the first ones in the nisp series, so probably ok.
            do 257 ifld = 1, nusp
	       if (isuponxy(ix,iy,ifld) .eq. 1) then
                  ix1 = ixp1(ix,iy)
                  iv2 = idxu(ix,iy,ifld)
		if (iseqalg(iv2)== 0.and.isnionxy(ix,iy,ifld)==1) then 
                  iv = idxn(ix,iy,ifld)
                  iv1 = idxn(ix1,iy,ifld)
                     yldot_np1 = resco(ix1,iy,ifld)/(vol(ix1,iy)*n0(ifld))
# need to use resco rather than yldot if dtreal is added; recursive prob.
c ....            Fix limiter case with algebraic eqns, not ODEs          
                     if (iseqalg(iv).eq.1) then
                        if (isnupdot1sd == 0) then
                           nbvdot = yldot_np1*n0(ifld)
                        else
                           nbvdot = yldot(iv1)*n0(ifld)
                        endif
                        nbv = ni(ix1,iy,ifld)
                     elseif (iseqalg(iv1).eq.1) then
                        nbvdot = yldot(iv)*n0(ifld)
                        nbv = ni(ix,iy,ifld)
                     else
                        if (isnupdot1sd == 0) then
                          nbvdot = 0.5*(yldot(iv) + yldot_np1)*n0(ifld)
                        else
                          nbvdot = yldot(iv)*n0(ifld)
                        endif
                         nbv = 0.5*(ni(ix,iy,ifld) + ni(ix1,iy,ifld))
                     endif
                     yldot(iv2) = (yldot(iv2)*n0(ifld) - yl(iv2)*nbvdot)
     .                                                              /nbv
                endif
               endif
  257       continue
            if (isflxvar .eq. 0) then
               if(isteonxy(ix,iy) == 1) then
                 iv =  idxte(ix,iy)
	         if(iseqalg(iv) == 0) then
                   yldot(iv) = ( yldot(iv)*nnorm - 
     .                                       yl(iv)*nbedot ) / ne(ix,iy)
                 endif
               endif
	       if(istionxy(ix,iy) == 1) then
                 iv1 = idxti(ix,iy)
                 if (iseqalg(iv1) == 0) then
                   if(isupgon(1)==1) then
                     yldot(iv1) = ( yldot(iv1)*nnorm - 
     .                              yl(iv1)*(nbidot + cftiexclg*nbgdot) ) /
     .                               (nit(ix,iy) + cftiexclg*ni(ix,iy,2) )
                   else
                     yldot(iv1) = ( yldot(iv1)*nnorm - yl(iv1)*
     .                             ( nbidot + cngtgx(1)*nbg2dot(1) ) ) /
     .                            (nit(ix,iy) + cngtgx(1)*ng(ix,iy,1))
                   endif
                 endif
               endif
               do igsp = 1, ngsp 
                 if (igsp == 1 .and. isupgon(1) == 1) then
                   if (istgonxy(ix,iy,igsp) == 1) then
                     iv = idxtg(ix,iy,igsp)
                     if (iseqalg(iv).eq.0) then
#                       yldot(iv) = ( yldot(iv)*n0g(igsp) -
#     .                              yl(iv)*yldot(idxn(ix,iy,iigsp))*n0(ifld) )
#     .                               /ni(ix,iy,iigsp)
                       #..the above one causes problem when turning off ni(:,:,2)
                        yldot(iv) = ( yldot(iv)*n0g(igsp) -
     .                               yl(iv)*nbgdot ) / ni(ix,iy,iigsp)
                     endif
                   endif
                 else
                   if (istgonxy(ix,iy,igsp) == 1) then
                     iv = idxtg(ix,iy,igsp)
                     if (iseqalg(iv) == 0) then
#                       yldot(iv) = ( yldot(iv)*n0g(igsp) -
#     .                         yl(iv)*yldot(idxg(ix,iy,igsp))*n0g(igsp) )
#     .                               /ng(ix,iy,igsp)
                       #..the above one causes problem when turning off ng(:,:,2)
                        yldot(iv) = ( yldot(iv)*n0g(igsp) -
     .                               yl(iv)*nbg2dot(igsp) ) / ng(ix,iy,igsp)
                     endif
                   endif
                 endif
              enddo
            endif   #end if-test on isflxvar
  260    continue
 270  continue

      return
      end

c   -------------------------------------------------------------------------
      subroutine volsor

*     VOLSOR defines volume particle and power sources for the ions and 
*     electrons. Input is total current, power, and shape of 2-D Gaussians

      implicit none

      Use(Dim)      # nx,nxm,ny,nisp
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      Use(Comgeo)   # vol
      Use(RZ_grid_info)   # rm,zm
      Use(Volsrc)   # pwrsore,pwrsori,volpsor,pvole,pvoli,ivolcur,
                    # z0pe,z0pi,r0pe,r0pi,zwpe,zwpi,rwpe,rwpi,
                    # z0ni,r0ni,zwni,rwni,voljcsor,jcvsor,
                    # ix_sjcsor, ix_ejcsor, iy_sjcsor, iy_ejcsor,
                    # thetarot,rcutmin,zcutmin,effvng,
      Use(Phyvar)   # ev
      Use(Bcond)    # islimsor,rlimiter
      Use(Parallv)  # nxg,nyg
      Use(Xpoint_indices)  # ixpt1,ixpt2,iysptrx
      Use(Share)    # nxomit

*  -- local scalars --
      real effvni, effvup, effvpe, effvpi, effvjel, zc, rc, ivolcurt,
     &     ivolcurgt, mvolcurt
      real argr, argz
      integer isxjcsor, iexjcsor, isyjcsor, ieyjcsor, ifld, nj
      #Former Aux module variables
      integer ix,iy,igsp


c...  Initialize values and arrays
      nj = nxomit
      effvni = 0.   
      effvup = 0.
      effvpe = 0.
      effvpi = 0.
      effvjel = 0.
      ivolcurt = 0.
      mvolcurt = 0.
      ivolcurgt = 0.
      do igsp = 1, ngsp
        effvng(igsp) = 0.
      enddo
      call s2fill (nx+2, ny+2, 0., pwrsori, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., pwrsore, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., voljcsor, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., pondpot, 1, nx+2)
      do ifld = 1, nisp
        call s2fill (nx+2, ny+2, 0., volpsor(0:nx+1,0:ny+1,ifld), 1, nx+2)
        call s2fill (nx+2, ny+2, 0., volmsor(0:nx+1,0:ny+1,ifld), 1, nx+2)
        ivolcurt = ivolcurt + ivolcur(ifld)
        mvolcurt = mvolcurt + mvolcur(ifld)
      enddo
      do igsp = 1, ngsp
	call s2fill (nx+2, ny+2, 0., volpsorg(0:nx+1,0:ny+1,igsp), 1, nx+2)
        ivolcurgt = ivolcurgt + ivolcurg(igsp)
      enddo        
cccMER NOTE: generalize the following for multiple x-points
cc Define index ranges for a localized ion-loss sink; crude & temporary
      if (ix_sjcsor .gt. 0) then
        isxjcsor = ix_sjcsor
      else
        isxjcsor = (ixpt1(1) + ixpt2(1))/2
      endif
      if (ix_ejcsor .gt. 0) then
        iexjcsor = ix_ejcsor
      else
        iexjcsor = ixpt2(1)
      endif
      if (iy_sjcsor > 0) then
        isyjcsor = iy_sjcsor
      else
        isyjcsor = 1
      endif
      if (iy_ejcsor > 0) then
        ieyjcsor = iy_ejcsor
      else
        ieyjcsor = iysptrx
      endif

      if ( abs( pvoli+pvole+ivolcurt+mvolcurt+
     .                        ivolcurgt+jcvsor+ponderompot )
     .                      > 0. ) then   # skip, big jump to end

      do 20 iy = 0, ny+1
         if (rm(0+nj,iy,0).lt.rlimiter .or. islimsor.eq.1) then
         do 10 ix = 0, nx+1
            zc = z0ni - (rm(ix+nj,iy,0)-r0ni)*sin(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*cos(thetarot)
            rc = r0ni + (rm(ix+nj,iy,0)-r0ni)*cos(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*sin(thetarot)
            if (zc.lt.zcutmin .or. rc.lt.rcutmin) goto 10
             argz = min(25., ((zc-z0ni)/zwni)**2)
             argr = min(25., ((rc-r0ni)/rwni)**2)            
            effvni = effvni + vol(ix,iy) * exp(-argz -argr)
             argz = min(25., ((zc-z0up)/zwup)**2)
             argr = min(25., ((rc-r0up)/rwup)**2)            
            effvup = effvup + vol(ix,iy) * exp(-argz -argr)
             argz = min(25., ((zc-z0pe)/zwpe)**2)
             argr = min(25., ((rc-r0pe)/rwpe)**2)            
            effvpe = effvpe + vol(ix,iy) * exp(-argz -argr)
             argz = min(25., ((zc-z0pi)/zwpi)**2)
             argr = min(25., ((rc-r0pi)/rwpi)**2)            
            effvpi = effvpi + vol(ix,iy) * exp(-argz -argr)
            do igsp = 1, ngsp
               argz = min(25., ((zc-z0ng(igsp))/zwng(igsp))**2)
               argr = min(25., ((rc-r0ng(igsp))/rwng(igsp))**2)            
              effvng(igsp) = effvng(igsp) + vol(ix,iy)*exp(-argz-argr)
            enddo
cccMER For full double-null configuration, iysptrx is last closed flux surface.
cc  Temporary localized current source (for prompt ion loss)
            if (iy >= isyjcsor .and. iy <= ieyjcsor) then
              if (ix .ge. isxjcsor .and. ix .le. iexjcsor) then
                effvjel = effvjel + vol(ix,iy)
              endif
            endif
 10      continue 
         endif
 20   continue

      do 40 iy = 0, ny+1
        if (rm(0+nj,iy,0).lt.rlimiter .or. islimsor.eq.1) then
         do 30 ix = 0, nx+1
            zc = z0ni - (rm(ix+nj,iy,0)-r0ni)*sin(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*cos(thetarot)
            rc = r0ni + (rm(ix+nj,iy,0)-r0ni)*cos(thetarot) +
     .                  (zm(ix+nj,iy,0)-z0ni)*sin(thetarot)
            if (zc < zcutmin .or. rc < rcutmin) goto 30
               argz = min(25., ((zc-z0pe)/zwpe)**2)
               argr = min(25., ((rc-r0pe)/rwpe)**2)            
              pwrsore(ix,iy) =pvole*vol(ix,iy)*exp(-argz-argr)/effvpe
               argz = min(25., ((zc-z0pi)/zwpi)**2)
               argr = min(25., ((rc-r0pi)/rwpi)**2)            
              pwrsori(ix,iy) =pvoli*vol(ix,iy)*exp(-argz-argr)/effvpi
               argz = min(25., ((zc-z0pondp)/zwpondp)**2)
               argr = min(25., ((rc-r0pondp)/rwpondp)**2)            
              pondpot(ix,iy) = ponderompot*exp(-argz-argr)
              do ifld = 1, nisp
                 argz = min(25., ((zc-z0ni)/zwni)**2)
                 argr = min(25., ((rc-r0ni)/rwni)**2)            
                volpsor(ix,iy,ifld) = ivolcur(ifld) * vol(ix,iy) *
     .                                   exp(-argz-argr)/(effvni*ev)
                 argz = min(25., ((zc-z0up)/zwup)**2)
                 argr = min(25., ((rc-r0up)/rwup)**2)            
                volmsor(ix,iy,ifld) = mvolcur(ifld) * vol(ix,iy) *
     .                                   exp(-argz-argr)/(effvup*ev)
              enddo
              do igsp = 1, ngsp
                 argz = min(25., ((zc-z0ng(igsp))/zwng(igsp))**2)
                 argr = min(25., ((rc-r0ng(igsp))/rwng(igsp))**2)            
                volpsorg(ix,iy,igsp) = ivolcurg(igsp) * vol(ix,iy) *
     .                                 exp(-argz-argr)/(effvng(igsp)*ev)
              enddo
cccMER For full double-null configuration, iysptrx is last closed flux surface
cc  Temporary localized current source (for prompt ion loss)
              if (iy >= isyjcsor .and. iy <= ieyjcsor) then
                if (ix .ge. isxjcsor .and. ix .le. iexjcsor) then
                  voljcsor(ix,iy) = jcvsor*vol(ix,iy) / effvjel
                endif
              endif
cccs            endif  # test on zcutmin and rcutmin
 30      continue
        endif   # test on rlimiter and islimsor
 40   continue

      endif  # end of large if beginning if (abs(pvoli+pvole+...

      return
      end
c-----------------------------------------------------------------------
      subroutine resid (t, yl, yldot, cj, delta, ires, rpar, ipar)

c ... Calculate residuals for differential-algebraic solution of the
c     UEDGE equations.
c
c     Note that the residuals can be expressed in terms of f(i) =
c     right-hand sides of rate equations used in differential
c     solution of the UEDGE equations.
c
c     The residuals have the form
c        delta(i) = f(i)
c     for the algebraic equations, and
c        delta(i) = f(i) - yldot(i)
c     for the differential equations.

      implicit none

c ... Input arguments:
      real t          # physical time
      real yl(*)      # most recent iterate of solution vector
      real yldot(*)   # most recent iterate of time-derivative of yl
      real cj         # proportional to 1/del_t; can rescale alg. constraint eq.
      real rpar(*)    # real parameter communication (not used)
      integer ipar(*) # integer parameter communication

c ... In-out argument:
      integer ires    # 0 on entry, set to -1 or -2 if errors occur

c ... Output argument:
      real delta(*)   # residuals

c ... Function:
      integer res_algeb

c ... Local variables:
      integer neq, i, ifail

c ... Get total number of equations (all grid points).
      neq = ipar(1)

c ... Calculate f(i), storing them in delta(i).
      call rhsdpk (neq, t, yl, delta, ifail)
      if (ifail .ne. 0) then
         ires = -1
         return
      endif

c ... Loop through all i, skipping those that correspond to algebraic
c     equations, and subtract yldot(i) from delta(i).
      do i = 1, neq
         if (res_algeb (i) .eq. 0) delta(i) = delta(i) - yldot(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ffun (neq, t, yl, yldot)

c ... Calculate the right-hand sides of the UEDGE rate equations.

      implicit none

      Use(Math_problem_size)   # neqmx
      Use(Constraints)         # icflag,rlx,icnstr,ylprevc,ylchng

c ... Input arguments used for all entry points:
      integer neq
      real yl(*)   # most recent iterate of solution vector

c ... Additional input argument required for ffun and rhsvd:
      real t         # physical time

c ... Additional input arguments required for rhsvd:
      real rpar(*)
      integer ipar(*)

c ... Output argument used for all entry points:
      real yldot(neq)   # right-hand sides

c ... Output argument for rhsvd and rhsdpk:
      integer ifail

c ... Local variables:
      integer i, ivar
      real tloc, tau, rlxl
      data tau /1.e0/, rlxl /1.e20/   # dummy argument here for cnstrt

c ... Beginning of execution for call ffun (by newton or lsode).
      goto 8

c ... Beginning of execution for call rhsvd (by vodpk), check constraints
      entry rhsvd (neq, t, yl, yldot, rpar, ipar, ifail)
      
      if (icflag .gt. 0) then     
         if (icflag .eq. 2) rlxl = rlx
         do 5 i = 1, neq     
            ylchng(i) = yl(i) - ylprevc(i)
 5       continue
         call cnstrt (neq,ylprevc,ylchng,icnstr,tau,rlxl,ifail,ivar) 
         if (ifail .ne. 0) then
            call remark ('***Constraint failure in VODPK, dt reduced***')
            write (*,*) 'variable index = ',ivar,'   time = ',t
            goto 20
         endif
         call scopy (neq, yl, 1, ylprevc, 1)  #put yl into ylprevc
      else
         ifail = 0
      endif

      go to 8

c ... Beginning of execution for call rhsdpk (by daspk), check constraints
      entry rhsdpk (neq, t, yl, yldot, ifail)
      
      if (icflag .gt. 0 .and. t .gt. 0.) then     
         if (icflag .eq. 2) rlxl = rlx
         do 6 i = 1, neq     
            ylchng(i) = yl(i) - ylprevc(i)
 6       continue
         call cnstrt (neq,ylprevc,ylchng,icnstr,tau,rlxl,ifail,ivar) 
         if (ifail .ne. 0) then
            call remark ('***Constraint failure in DASPK, dt reduced***')
            write (*,*) 'variable index = ',ivar,'   time = ',t
            goto 20
         endif
      else
         ifail = 0
      endif
      call scopy (neq, yl, 1, ylprevc, 1)  #put yl into ylprevc 

 8    tloc = t
      go to 10

c ... Beginning of execution for call rhsnk (by nksol).
      entry rhsnk (neq, yl, yldot)
      tloc = 0.

c ... Calculate right-hand sides for interior and boundary points.
ccc 10   call convsr_vo (-1,-1, yl)  # test new convsr placement
ccc      call convsr_aux (-1,-1, yl) # test new convsr placement
 10   call pandf1rhs_interface ( neq, tloc, yl, yldot)

 20   continue
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_calc (neq, t, yl, yldot00, ml, mu, wk,
     .                     nnzmx, jac, ja, ia)

c ... Calculate Jacobian matrix (derivatives with respect to each
c     dependent variable of the right-hand side of each rate equation).
c     Lower and upper bandwidths are used to select for computation
c     only those Jacobian elements that may be nonzero.
c     Estimates of Jacobian elements are computed by finite differences.
c     The Jacobian is stored in compressed sparse row format.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real t           # physical time
      real yl(*)       # dependent variables
      real yldot00(neq+2) # right-hand sides evaluated at yl
      integer ml, mu   # lower and upper bandwidths
      integer nnzmx    # maximum number of nonzeros in Jacobian

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... Output arguments:
      real jac(nnzmx)     # nonzero Jacobian elements
      integer ja(nnzmx)   # col indices of nonzero Jacobian elements
      integer ia(neq+1)   # pointers to beginning of each row in jac,ja

c ... Common blocks:
      Use(Dim)                     # nx,ny,
                                   # nusp[for fnorm not used here]
      Use(Timing)                  # istimingon,ttjstor,ttotjf,ttimpjf
      Use(Math_problem_size)       # neqmx,numvar
      Use(Grid)                    # ngrid,ig,ijac,ijactot
      Use(Indexes)                 # igyl,iseqalg
      Use(Variable_perturbation)   # del,dylconst
      Use(Jacobian_clipping)       # jaccliplim,istopjac,irstop,icstop
      Use(Jacobian_csc)            # rcsc,jcsc,icsc,yldot_pert
      Use(Ynorm)                   # suscal,sfscal
      Use(UEpar)                   # isphion,isnewpot,svrpkg,isbcwdt
      Use(Model_choice)            # iondenseqn
      Use(Imprad)                  # isimpon
      Use(Bcond)                   # isextrnpf,isextrtpf,isextrngc,
                                   # isextrnw,isextrtw
      Use(Parallv)                 # nxg,nyg
      Use(Time_dep_nwt)            # nufak,dtreal,ylodt,dtuse
      Use(Selec)                   # yinc
      Use(Jacaux)                  # ExtendedJacPhi
      Use(Flags)                  # ExtendedJacPhi

c ... Functions:
      logical tstguardc
      real(Size4) gettime
cc      real(Size4) ranf

c ... Local variables:
      integer nnz, ii, iv, ii1, ii2, xc, yc, ix, iy
      real yold, dyl, jacelem
      real(Size4) sec4, tsjstor, tsimpjf, dtimpjf

ccc      save

c ... Get initial value of system cpu timer.
      if (istimingon .eq. 1) tsjstor = gettime(sec4)

c ... Pause from BASIS if a ctrl_c is typed
      call ruthere

      ijac(ig) = ijac(ig) + 1

      if ((svrpkg.eq.'nksol') .and. (iprint .ne. 0)) write(*,*) ' Updating Jacobian, npe =  ', 
     .                                                          ijac(ig)

c ... Set up diagnostic arrays for debugging
      do iv = 1, neq
        yldot_unpt(iv) = yldot00(iv)  # for diagnostic only
        yldot_pert(iv) = 0.
      enddo

c############################################
c ... Begin loop over dependent variables.
c############################################
      nnz = 1
      do iv = 1, neq

ccc ... Only perturb variables that are being solved for (for Daspk option)
ccc      if (iseqon(iv) .eq. 0) goto 18

c ... Set beginning and ending indices of right-hand sides that might be
c     perturbed.
         ii1 = max(iv-mu, 1)
         ii2 = min(iv+ml, neq)
c... >>>>> WARNING (added by J.Guterl)
c... >>>>> WARNING -> Do not forget to keep consistency between this subroutine and localjacbuilder in parallel.F90 when modifiying jacobian eval.
c... >>>>> WARNING
c ... Reset range if this is a potential perturbation with isnewpot=1
ccc         if (isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
cc  Comment out;storage for Jac inconsistent if mu & ml above not used
cc  Reported by R. Smirnov Feb. 2020
cc         if (isphion*isnewpot.eq.1 .and. mfnkso < 4) then
cc            ii1 = max(iv-4*numvar*nx, 1)      # 3*nx may be excessive
cc            ii2 = min(iv+4*numvar*nx, neq)    # 3*nx may be excessive
cc         endif

cc J.Guterl: I disagree with the statement above about inconsistent storage. Storage is consistent for abs(mfnksol)<4 with ilut preconditionner
cc         : as the reconstruction of the jacobian is not based on ubw and lbw.
cc         : For the banded preconditioner, out of bound checks have been added to csrbnd
cc         : where lbw abd ubw are simply given to estimate the max dimension of the jacobian.
cc         : The banded method can be then used even without fully consistent ubw and lbw.
cc         : Further validation of this statement is needed though.
cc         : Also, stencil analysis can be performed with the UEDGEToolBox extension.

c... Option added to keep extended Jacobian when phi eq. is on (Added by J.Guterl). See comments above
        if (ExtendedJacPhi.gt.0) then
        if(isphion*isnewpot.eq.1 .and. mod(iv,numvar).eq.0) then
            ii1 = max(iv-4*numvar*nx, 1)      # 3*nx may be excessive
            ii2 = min(iv+4*numvar*nx, neq)    # 3*nx may be excessive
        endif
       endif

c ... Reset range if extrapolation boundary conditions are used
cc  This reset of ii1,2 may also cause storage prob.; see just above
         if (isextrnpf+isextrtpf+isextrngc+isextrnw+isextrtw.gt.0) then
            ii1 = max(iv-2*numvar*(nx+3), 1)      # guess to include extrap. bc
            ii2 = min(iv+2*numvar*(nx+3), neq)    # guess to include extrap. bc
         endif

c ... Initialize all of those right-hand sides to their unperturbed
c     values.
         do ii = ii1, ii2   # below wk is reset, but only over limited range
            wk(ii) = yldot00(ii)
         enddo

c ... Set spatial-location indices for this dependent variable.
         xc = igyl(iv,1)
         yc = igyl(iv,2)

c ... Save value of dependent variable, then perturb it.
c     The perturbation to the variable is proportional to parameter
c     del and to a measure of the size of the variable.  That measure
c     increases with the absolute value of the variable if it exceeds
c     the typical size given by dylconst/suscal but can never be less
c     than that typical size.
         yold = yl(iv)
         dyl = del * (abs(yold) + dylconst / suscal(iv))
         yl(iv) = yold + dyl

c ... Calculate right-hand sides near location of perturbed variable.
         call pandf1 (xc, yc, iv, neq, t, yl, wk)

c ... Calculate possibly nonzero Jacobian elements for this variable,
c     and store nonzero elements in compressed sparse column format.
         jcsc(iv) = nnz      # sets index for first Jac. elem. of var. iv
         do ii = ii1, ii2
            jacelem = (wk(ii) - yldot00(ii)) / dyl
ccc            jacelem = (wk(ii) - yldot0(ii)) / (2*dyl)  # for 2nd order Jac
c ...  Add diagonal 1/dt for nksol
            if (((svrpkg.eq."nksol") .or. (svrpkg.eq."petsc")) .and. iv.eq.ii) then
              if (iseqalg(iv)*(1-isbcwdt).eq.0) then
                jacelem = jacelem - 1/dtuse(iv)
              endif
              ix = igyl(iv,1)
              iy = igyl(iv,2)
              if (idxphi(ix,iy)==iv .and. dtphi<1e10) then #selects phi eqn
                jacelem = jacelem - 1/dtphi
              endif
            endif

c ...  Add a pseudo timestep to the diagonal ## if eqn is not algebraic
            if (svrpkg .ne. "cvode" .and. nufak .gt. 0) then
               if (iv.eq.ii .and. yl(neq+1).eq.1) 
     .             jacelem = jacelem - nufak  #omit .and. iseqalg(iv).eq.0)
            endif
            if (abs(jacelem*sfscal(iv)) .gt. jaccliplim) then
               if (nnz .gt. nnzmx) then
                  write(STDOUT,*)
     .             '*** jac_calc -- More storage needed for Jacobian.',
     .             ' Storage exceeded at (i,j) = (',ii,',',iv,').',
     .             ' Increase lenpfac.'
	          call xerrab("")
               endif
cc               if (rdoff.ne.0.e0) jacelem=jacelem*(1.0e0+ranf()*rdoff)
               rcsc(nnz) = jacelem
               icsc(nnz) = ii
               nnz = nnz + 1
            endif
           
            if (istopjac.gt.0 .and. ii.eq.irstop .and. iv.eq.icstop) then
               yldot_pert(ii) = wk(ii)      # for diagnostic only
               if (istopjac == 2) then
                 yl(iv) = yold
                 call pandf1 (xc, yc, iv, neq, t, yl, wk)
               endif
               call remark("***** non-zero jac_elem at irstop,icstop")
               write(*,*) 'irstop = ', irstop, ', icstop = ', icstop
               call xerrab("")
            endif

         enddo  # end of ii loop over equations

c ... Restore dependent variable yl & assoicated plasma vars near perturbation
         yl(iv) = yold
         call pandf1 (xc, yc, iv, neq, t, yl, wk)

c...  If this is the last variable before jumping to new cell, reset pandf 
ccc  Call not needed because goto 18 svrpkg=daspk option disabled above
ccc         if (mod(iv,numvar).eq.0 .and. isjacreset.ge.1) then
ccc            call pandf1 (xc, yc, iv, neq, t, yl, wk)
ccc         endif
   
c ... End loop over dependent variables and finish Jacobian storage.
c##############################################################      
      enddo             # end of main iv-loop over yl variables
c##############################################################

      jcsc(neq+1) = nnz

c ... Convert Jacobian from compressed sparse column to compressed
c     sparse row format.
      call csrcsc (neq, 1, 1, rcsc, icsc, jcsc, jac, ja, ia)

c ... Count Jacobian evaluations, both for total and for this case
      ijactot = ijactot + 1   #note: ijactot set 0 in exmain if icntnunk=0

c ... Accumulate cpu time spent here.
      if (istimingon .eq. 1) ttjstor = ttjstor + gettime(sec4) - tsjstor
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_lu_decomp (neq, jac, ja, ia, wp, iwp)

c ... Compute LU decomposition of Jacobian and return it in any one
c     of the storage formats.

      implicit none

c ... Input arguments:
      integer neq        # total number of equations (all grid points)
      real jac(*)        # nonzero Jacobian elements
      integer ja(*)      # col indices of nonzero Jacobian elements
      integer ia(neq+1)  # pointers to beginning of each row in jac,ja

c ... Output arguments:
      real wp(*)         # matrix elements of LU
      integer iwp(*)     # sizes and array indices for elements of LU
      
c ... Common blocks:
      Use(Timing)                  # ttmatfac
      Use(Decomp)                  # lbw,ubw
      Use(Grid)                    # ngrid,ig,ijac
      Use(Condition_number)        # rconds
      Use(Preconditioning)         # premeth,lenplumx
      Use(Ilutv)                    # lfililut,tolilut
      Use(Nonzero_diagonals)       # lfilinel,ndiagmx,ndiag,ndiagm,
                       # adiag(neq,ndiagmx),siginel(neq),fmuinel(neq)
                       # iwkd1(2*neq-1),iwkd2(ndiagmx),rwkd(ndiagmx)
      Use(Temporary_work_arrays)   # rwk1,rwk2,iwk1,iwk2,iwk3

c ... Function:
      real(Size4) gettime

c ... Local variables:
      integer lowd, ierr, i, idum(1)
      real rcond, dum(1)
      real(Size4) sec4
      real tsmatfac

c ... Convert compressed sparse row to banded format and use exact
c     factorization routine sgbco from Linpack/SLATEC.
      if (premeth .eq. 'banded') then
         lowd = 2 * lbw + ubw + 1
         call csrbnd (neq, jac, ja, ia, 0, wp, lowd, lowd,
     .                lbw, ubw, ierr)
         if (ierr .ne. 0) then
            write(STDOUT,*)
     .         '*** jac_lu_decomp -- csrbnd returned ierr =', ierr
            call xerrab("")
         endif
         tsmatfac = gettime(sec4)
         call sgbco (wp, lowd, neq, lbw, ubw, iwp(4), rcond, rwk1)
         iwp(1) = lowd
         iwp(2) = lbw
         iwp(3) = ubw

c ... Save condition number.
         i = ijac(ig)
         if (i .le. 300) rconds(i,ig) = rcond
         go to 99
      endif

c ... If sparse Jacobian matrix is in compressed sparse row storage
c     format, ...
      if (premeth .eq. 'ilut') then

c ... Reorder Jacobian rows and columns, if desired.
         call jac_reorder (neq, jac, ja, ia, wp, iwp(neq+2), iwp)

c ... Use incomplete factorization routine ilut from SparsKit.
         tsmatfac = gettime(sec4)
         call ilut (neq,jac,ja,ia,lfililut,tolilut,wp,iwp(neq+1),
     .              iwp,lenplumx,rwk1,rwk2,iwk1,
     .              iwk2,iwk3,ierr) 
         if (ierr .ne. 0) then
            write(STDOUT,*) ' Error return from ilut:  ierr = ',ierr
            write(STDOUT,9000)
 9000       format(
     ./'    ierr >  0   --> Zero pivot encountered at step number ierr.'
     ./'    ierr = -1   --> Error. input matrix may be wrong.'
     ./'                     (The elimination process has generated a'
     ./'                     row in L or U with length > n.)'
     ./'    ierr = -2   --> Matrix L overflows.'
     ./'    ierr = -3   --> Matrix U overflows.'
     ./'    ierr = -4   --> Illegal value for lfililut.'
     ./'    ierr = -5   --> Zero row encountered.'
     ./'    '
     ./'    For ierr = -2 or -3, increase the value of lenplufac or'
     ./'    decrease the value of lfililut if lenplufac cannot be'
     ./'    increased.'
     .)
            call xerrab("")
         endif

c ... Use incomplete factorization routine precond5 from INEL.
c     SparsKit routines are used in preliminary steps to convert to
c     diagonal format.
      elseif (premeth .eq. 'inel') then

c ... Get the number of nonzero diagonals and the maximum in LU.
         call infdia (neq, ja, ia, iwkd1, ndiag)
         if (ndiag .gt. ndiagmx) then
            call remark('More storage for diagonals of the Jacobian')
            call remark('is needed.  Increase value of ndiagmx.')
            call xerrab("")
         endif
         ndiagm = min(lfilinel+ndiag, ndiagmx)
         iwp(1) = ndiag
         iwp(2) = ndiagm

c ... Convert to diagonal format.
         call csrdia (neq, ndiag, 10, jac, ja, ia, neq, adiag,
     .                iwp(3), dum, idum, idum, iwkd1)

c ... Reorder rows to be in increasing column order.
         call cdiagsrt (neq, adiag, ndiag, iwp(3), iwkd1, iwkd2,
     .                  rwkd)

c ... Finally, calculate approximate LU decomposition.
         tsmatfac = gettime(sec4)
         call precond5 (neq, ndiag, ndiagm, adiag, wp, rwk2, rwk1,
     .                  iwk3, iwk2, siginel, fmuinel, iwp(3))
      endif

c ... Accumulate cpu time spent here.
 99   ttmatfac = ttmatfac + (gettime(sec4) - tsmatfac)
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_reorder (neq, jac, ja, ia, awk, jwk, iwk)

c ... If desired, reorder the Jacobian matrix.

      implicit none

c ... Input arguments:
      integer neq        # total number of equations (all grid points)
      real jac(*)        # nonzero Jacobian elements
      integer ja(*)      # column indices of nonzero Jacobian elements
      integer ia(neq+1)  # indices of 1st nonzero element in each row

c ... Work-array arguments:
      real awk(*)
      integer jwk(*)
      integer iwk(neq+1)

c ... Common blocks:
      Use(Preconditioning) # premeth
      Use(Timing)          # ttjreorder
      Use(Jacreorder)      # perm,qperm,levels,nlev,mask,maskval,ireorder

c ... Function:
      real(Size4) gettime

c ... Local variables:
      real(Size4) sec4
      real tsjreorder
      integer i, nfirst

c ... Get initial value of system cpu timer.
      tsjreorder = gettime(sec4)

      if ((ireorder .eq. 1) .and. (premeth .eq. 'ilut')) then

c ... Copy jac, ja, and ia to awk, jwk, and iwk.
         call atob (neq, jac, ja, ia, awk, jwk, iwk)

c ... Perform a Cuthill-McKee reordering of the Jacobian.
         nfirst = 1
         perm(1) = 0
         do i = 1, neq
            mask(i) = 1
         enddo
         maskval = 1
         qperm(1) = 1
         call bfs (neq,jwk,iwk,nfirst,perm,mask,maskval,qperm,levels,
     .             nlev)

c ... Reverse the permutation to obtain the reverse Cuthill-McKee
c     reordering.
         call reversp (neq,qperm)

c ... Calculate the inverse of qperm and put it in perm.
         do i = 1, neq
            perm(qperm(i)) = i
         enddo

c ... Permute rows and columns of Jacobian using perm.
         call dperm (neq,awk,jwk,iwk,jac,ja,ia,perm,perm,1)

c ... End of If block
      endif

c ... Accumulate cpu time spent here.
      ttjreorder = ttjreorder + (gettime(sec4) - tsjreorder)
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_norm_rows (neq, jac, ja, ia)

c ... If desired, normalize each row using one of three types of norm.

      implicit none

c ... Input arguments:
      integer neq        # total number of equations (all grid points)
      real jac(*)        # nonzero Jacobian elements
      integer ja(*)      # column indices of nonzero Jacobian elements
      integer ia(neq+1)  # indices of 1st nonzero element in each row

c ... Common blocks:
      Use(Timing)   # ttjrnorm
      Use(Math_problem_size)   # neqmx(for arrays in Jacaux not used here)
      Use(Jacaux)   # isrnorm,normtype,fnormnw

c ... Function:
      real(Size4) gettime

c ... Local variables:
      real(Size4) sec4
      real tsjrnorm

c ... Get initial value of system cpu timer.
      tsjrnorm = gettime(sec4)

      if (isrnorm .eq. 1) call roscal (neq, 0, normtype, jac, ja, ia,
     .                               fnormnw, jac, ja, ia)

c ... Accumulate cpu time spent here.
      ttjrnorm = ttjrnorm + (gettime(sec4) - tsjrnorm)
      return
      end
c-----------------------------------------------------------------------
      subroutine jac_sub_cj (neq, jac, ja, ia, cj)

c ... Loop through all Jacobian elements, skipping off-diagonal
c     elements and those that correspond to algebraic equations, and
c     subtract cj.

      implicit none

c ... Input arguments:
      integer neq        # total number of equations (all grid points)
      integer ja(*)      # column indices of nonzero Jacobian elements
      integer ia(neq+1)  # indices of 1st nonzero element in each row
      real cj            # scalar to be subtracted

c ... In-out argument:
      real jac(*)    # nonzero Jacobian elements

c ... Function:
      integer res_algeb

c ... Local variables:
      integer i   # index to row (i.e., equation)
      integer k   # index to nonzero element and its column index

      do i = 1, neq
         do k = ia(i), ia(i+1)-1
            if (ja(k) .ne. i) go to 10
            if (res_algeb (i) .eq. 1) go to 10
            jac(k) = jac(k) - cj
 10         continue
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function res_algeb (i)
      implicit none
c ... Return 1 if the equation index corresponds to a potential
c     equation or a boundary point, otherwise return 0.

c ... Input argument:
      integer i   # equation index

c ... Common blocks:
      Use(Dim)       # nx,ny,nxpt
      Use(Xpoint_indices)      # ixlb,ixrb
      Use(Math_problem_size)   # neqmx(for arrays in Indexes not used here)
      Use(Indexes)   # igyl,idxphi
      Use(UEpar)     # isphionxy

c ... Local variables:
      integer ix, iy, jx

      ix = igyl(i,1)   # get spatial location for this equation
      iy = igyl(i,2)
      if (isphionxy(ix,iy).eq. 1) then
         if (i .eq. idxphi(ix,iy)) then
            res_algeb = 1
            return
         endif
      endif
      if (ix .eq. 0 .or. ix .eq. nx+1 .or. iy .eq. 0 .or. iy .eq. ny+1)
     .   then
         res_algeb = 1
      else
         res_algeb = 0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine jacd1 (t, yl, yldot, pd, cj, rpar, ipar)

c ... Calculate matrix A in banded format for daspk solver.
c     Matrix A has the form
c        A(i,j) = dG(i)/dYL(j) + CJ*dG(i)/dYLDOT(j)
c     where G is the residual described in the comments for resid.
c     Therefore, for boundary points and for the potential equation,
c        A(i,j) = df(i)/dYL(j) = Jacobian
c     Otherwise,
c        A(i,j) = df(i)/dYL(j) - CJ*(Kronecker delta ij)

      implicit none

c ... Input arguments:
      real t          # physical time
      real yl(*)      # most recent iterate of solution vector
      real yldot(*)   # most recent iterate of time-derivative of yl
                      # (not used)
      real rpar(*)    # real parameter communication (not used)
      integer ipar(3) # integer parameter communication
      real cj         # scalar provided by daspk

c ... Output argument:
      real pd(*)      # matrix A in (Linpack) banded format

c ... Common blocks:
      Use(Math_problem_size)   # neqmx
      Use(Jacaux)     # yldot0,yldot1
c     Temporary_work_arrays cannot be used here because neq is not
c     an argument.
      Use(Jacobian)   # nnzmx,jac,jacj,jaci

c ... Local variables:
      integer neq, lbw, ubw
      integer lowd, ierr

c ... Get total number of equations (all grid points) and
c     lower and upper bandwidths.
      neq = ipar(1)
      lbw = ipar(2)
      ubw = ipar(3)

c ... Calculate right-hand sides at unperturbed values of variables for
c ... Jacobian calculation.  
      call ffun (neq, t, yl, yldot0)

c ... Calculate Jacobian of right-hand sides of UEDGE equations.
      call jac_calc (neq, t, yl, yldot0, lbw, ubw, yldot1,
     .               nnzmx, jac, jacj, jaci)

c ... Subtract cj from appropriate Jacobian elements.
      call jac_sub_cj (neq, jac, jacj, jaci, cj)

c ... Convert Jacobian from compressed sparse row to (Linpack) banded
c     format.
      lowd = 2 * lbw + ubw + 1
      call csrbnd (neq, jac, jacj, jaci, 0, pd, lowd, lowd,
     .             lbw, ubw, ierr)
      if (ierr .ne. 0) then
         write(STDOUT,*) '*** jacd1 -- ierr =', ierr
         call xerrab("")
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine jacd2 (resid, ires, neq, t, yl, yldot, rewt, savr, wk,
     .                  h, cj, wp, iwp, ierr, rpar, ipar)

c ... Calculate preconditioning matrix that approximates A (see
c     subroutine jacd1) for daspk solver.

      implicit none

c ... Input arguments:
      external resid   # function that evaluates residuals (not used)
      integer ires     # error flag from resid (not used)
      integer neq      # total number of equations (all grid points)
      real t           # physical time
      real yl(*)       # most recent iterate of solution vector
      real yldot(neq)  # most recent iterate of time-derivative of yl
                       # (not used)
      real rewt(neq)   # reciprocal error weights 
      real savr(neq)   # residual values G(t,yl,yldot) (not used)
      real h           # step size (not used)
      real cj          # scalar provided by daspk
      real rpar(*)     # real parameter communication (not used)
      integer ipar(3)  # integer parameter communication

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... Output arguments:
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # array indices for elements of LU
      integer ierr     # error flag (0 means success, else failure)

c ... Common blocks:
      Use(Math_problem_size)   # neqmx
      Use(Jacaux)     # yldot0,jscalcol
      Use(Jacobian)   # nnzmx,jac,jacj,jaci

c ... Local variables:
      integer lbw, ubw, ii
      real sqrtn

c ... Unpack parameters.
      lbw = ipar(2)
      ubw = ipar(3)

c ... Calculate right-hand sides at unperturbed values of variables.
c     (These values could be obtained with more coding but less
c     computation by taking savr and adding yldot for differential
c     equations.)
      call ffun (neq, t, yl, yldot0)

c ... Calculate Jacobian.
      call jac_calc (neq, t, yl, yldot0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)

c ... Subtract cj from appropriate Jacobian elements.
      call jac_sub_cj (neq, jac, jacj, jaci, cj)

C ... Multiply Jacobian columns by inverse of scaling vector REWT.
C     In PSOLILU, the WGHT array equals REWT/SQRT(NEQ), so we must
C     be consistent here.  Copied from coding by Peter Brown (3/10/97)
      if (jscalcol .eq. 1) then
         sqrtn = sqrt(real(neq))
         do 10 ii = 1, neq
            wk(ii) = sqrtn / rewt(ii)
 10      continue
         call amudia (neq, 0, jac, jacj, jaci, wk, jac, jacj, jaci)
      endif

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine jacnw (neq, yl, f0, dt, wk, wp, iwp)

c ... Calculate LU decomposition of the Jacobian matrix
c     for use as a preconditioner by the newton solver.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real yl(*)     # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)
      real dt          # false timestep to improve condition number

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... Output arguments:
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # array indices for elements of LU

c ... Common blocks:
      Use(Decomp)     # lbw,ubw
      Use(Jacobian)   # nnzmx,jac,jacj,jaci

c ... Local variables:
      real tp

c ... Flag these calls to RHS (pandf) as for the Jacobian
      yl(neq+1) = 1.

c ... Call pandf to set terms with yl(neq+1) flag on for Jacobian
      tp = 0.
      call ffun (neq, tp, yl, f0)

c ... Calculate Jacobian matrix.
      call jac_calc (neq, tp, yl, f0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)

      yl(neq+1) = -1.        # Turn-off Jacobian flag for pandf

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

      return
      end
c-----------------------------------------------------------------------
      subroutine jacvd (f, neq, tp, yl, ylsv, rewt, fty, wk, hrl1,
     .                  wp, iwp, ierr, rpar, ipar)

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix), and return
c     exact or approximate LU decomposition for the vodpk solver.

      implicit none

c ... Input arguments:
      external f       # function that evaluates f(tp,yl)
      integer neq      # total number of equations (all grid points)
      real tp          # physical time
      real yl(*)     # most recent iterate of solution vector
      real ylsv(neq)   # a copy of yl (not used)
      real rewt(neq)   # reciprocol error weights (not used)
      real fty(neq)    # function values f(tv,yl)
      real hrl1        # scalar provided by vodpk
      real rpar(*)     # real parameter communication (not used)
      integer ipar(2)  # integer parameter communication

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... Output arguments:
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # array indices for elements of LU
      integer ierr     # error flag (0 means success, else failure)

c ... Common block:
      Use(Jacobian)   # nnzmx,jac,jacj,jaci

c ... Local variables:
      integer lbw, ubw, nnz

c ... Unpack parameters.
      lbw = ipar(1)
      ubw = ipar(2)

c ... Calculate Jacobian.
      call jac_calc (neq, tp, yl, fty, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)

c ... Jacobian elements could be saved here for reuse if change in hrl1
c     is more significant than changes in yl.

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix).
      nnz = jaci(neq+1) - 1
      call sscal (nnz, -hrl1, jac, 1)
      call aplsca (neq, jac, jacj, jaci, 1., iwp)

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine FPRECO (neq, tp, yl, f0, jok, jcur, hrl1, rewt, h,
     .                   uround, nfe, v1, v2, v3, ierr)

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix), and return
c     exact or approximate LU decomposition for the pvode solver.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real tp          # physical time
      real yl(neq)     # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)
      integer jok      # Jacobian ok used by pvode
      real hrl1        # scalar provided by pvode
      real rewt(neq)   # reciprocol error weights (not used)
      real h           # step-size from pvode, don't change
      real uround      # from pvode, don't change
      integer nfe      # input/output for number of RHS evaluations

c ... Work-array argument:
      real v1(neq)     # work space available to this subroutine
      real v2(neq)     # work space available to this subroutine
      real v3(neq)     # work space available to this subroutine

c ... Output arguments:
      integer jcur     # =1 if FPRECO has updated Jacobian
      integer ierr     # error flag (0 means success, else failure)

c ... Common block:
      Use(Jacobian)   # nnzmx,jac,jacj,jaci
      Use(Decomp)     # lbw,ubw
      Use(Jac_work_arrays) # iwwp, wwp (use instead of iwp and wp for cvode cases)

c ... Local variables:
      integer nnz

c ... Calculate Jacobian.

      call jac_calc (neq, tp, yl, f0, lbw, ubw, v1, nnzmx, 
     .               jac, jacj, jaci)

c ... Jacobian elements could be saved here for reuse if change in hrl1
c     is more significant than changes in yl.

c ... Calculate (identity matrix) - hrl1 * (Jacobian matrix).
      nnz = jaci(neq+1) - 1
      call sscal (nnz, -hrl1, jac, 1)
      call aplsca (neq, jac, jacj, jaci, 1., iwwp)

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wwp, iwwp)

cxqx Alan suggestion
      jcur = 1

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine jacvnk (neq, yl, f0, v, z, wp, iwp)

c ... This subroutine is present merely to satisfy the requirements of

c     the option mdif=1 is not working in nksol and will be removed in
c     the future.  The Jacobian for nksol is calculated through psetnk.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real yl(*)       # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)
      real v(neq)      # arbitrary vector
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # array indices for elements of LU

c ... Output argument:
      real z(neq)      # (most recently calculated Jacobian) * v

      return
      end
c-----------------------------------------------------------------------
      subroutine psetnk (neq, yl, f0, su, sf, wk, f, wp, iwp, ierr)

c ... Calculate LU decomposition of the Jacobian matrix
c     for use as a preconditioner by the nksol solver.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real yl(*)       # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)
      real su(neq)     # scale factors for yl
      real sf(neq)     # scale factors for function values f(yl)
      external f       # function that evaluates residuals f(yl)

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... Output arguments:
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # array indices for elements of LU
      integer ierr     # error flag (0 means success, else failure)

c ... Common blocks:
      Use(Decomp)         # lbw,ubw
      Use(Jacobian)       # nnzmx,jac,jacj,jaci
      Use(Math_problem_size)   # neqmx
      Use(Dim)            # nx,ny
      Use(Time_dep_nwt)   # nufak,ydt_max,ydt_max0,alfnuf,expnuf,nufak0
                          # inufaknk,dtoptx,dtoptv
      Use(Indexes)        # idxn,idxu,idxte,idxti,idxng,idxphi
      Use(UEpar)          # isnion,isupon,isteon,istion,isngon,isphion,isnionxy,
                          # isuponxy,isteonxy,istionxy,isngonxy,isphionxy
      Use(Jac_work_arrays) # iwwp,wwp,liwp,lwp  # diagnostic arrays in this sub

c ... Local variables:
      real tp
      integer i

c ... Calculate maximum of f0*sf to control yl(neq+2) = nufak
      ydt_max = 1.e-100
      do i = 1, neq    # need to avoid neq+1 and neq+2
         if (abs(f0(i)*sf(i)) .gt. ydt_max) 
     .                           ydt_max = abs(f0(i)*sf(i))
      enddo
      if (ydt_max0 .eq. 0) ydt_max0 = ydt_max
      nufak = min(nufak*alfnuf*(ydt_max/ydt_max0)**expnuf, nufak0)
      if (inufaknk .eq. 1) then   # deter. if nufak is used in Krylov step
         yl(neq+2) = nufak
      else
         yl(neq+2) = 0.
      endif
      if (expnuf.ne.0.) write(*,*) ' nufak = ', nufak
      ydt_max0 = ydt_max
 
c ... Flag these calls to RHS (pandf) as for the Jacobian
      yl(neq+1) = 1.

c ... Call pandf to set terms with yl(neq+1) flag on for Jacobian
      call rhsnk (neq, yl, f0)

c ... Calculate Jacobian matrix. 
      tp = 0.
      call jac_calc_interface (neq, tp, yl, f0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)
      yl(neq+1) = -1.             # Turn-off Jacobian flag for pandf
      call rhsnk (neq, yl, f0)    # Reset f0 with nufak off

c ... Multiply Jacobian columns by inverse of scaling vector su.
      do i = 1, neq
         wk(i) = 1. / su(i)
      enddo
      call amudia (neq, 0, jac, jacj, jaci, wk, jac, jacj, jaci)

c ... Multiply Jacobian rows by scaling vector sf.  Scaling the
c     columns and rows allows nksol to work with quantities of
c     order unity regardless of the size of f0 and yl.
      call diamua (neq, 0, jac, jacj, jaci, sf, jac, jacj, jaci)

c ... Normalize Jacobian rows, if desired.
      call jac_norm_rows (neq, jac, jacj, jaci)

c ... Initialize iwwp and wwp
      do i = 1, liwp
        iwwp(i) = 0
      enddo
      do i = 1, lwp
        wwp(i) = 0.
      enddo

c ... Compute LU decomposition.
      call jac_lu_decomp (neq, jac, jacj, jaci, wp, iwp)

c ... Copy LU decomp into common arrays for diagostic
      do i = 1, liwp
        iwwp(i) = iwp(i)
      enddo
      do i = 1, lwp
        wwp(i) = wp(i)
      enddo

      ierr = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine psolnk (neq, yl, f0, su, sf, f, jvsol, wk, wp, iwp, bl,
     .                    ierr)

c ... Interface between linear-system solver nksol and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real yl(neq)     # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)
      real su(neq)     # scale factors for yl
      real sf(neq)     # scale factors for function values f(yl)
      external f       # function that evaluates residuals f(yl)
      external jvsol   # function that evaluates jac(u)*v
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # dimensions and array indices for elements of LU

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... In-out argument:
      real bl(neq)     # on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     # error flag (0 OK, >0 try again, <0 failure)

c ... Local variable:
      logical usingsu

      usingsu = .true.
      call psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine psold (neq, t, yl, yldot, f0, wk, cj, wght, wp, iwp,
     .                  bl, eplin, ierr, rpar, ipar)

c ... Interface between linear-system solver daspk and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real t           # physical time (not used here)
      real yl(neq)     # most recent iterate of solution vector
      real yldot(neq)  # most recent iterate of time-derivative of yl
      real f0(neq)     # G(t,yl,yldot) (not used)
      real cj          # scalar provided by daspk (not used)
      real wght(neq)   # error weights 
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # dimensions and array indices for elements of LU
      real eplin       # bound on solution error (not used)
      real rpar(*)     # real parameter communication (not used)
      integer ipar(*)  # integer parameter communication (not used)

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... In-out argument:
      real bl(neq)     # on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     # error flag (0 OK, >0 try again, <0 failure)

c ... Local variables:
      logical usingsu
##      real su(1)       # array that is not used when usingsu is false

      usingsu = .false.# this may get reset in psolbody if jscalcol=1
##      su(1) = 0.       # set unused array so flint will be happy
      call psolbody (neq, usingsu, wght, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine psolnw (neq, yl, wk, wp, iwp, bl, ierr)

c ... Interface between linear-system solver newton and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real yl(neq)     # most recent iterate of solution vector
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # dimensions and array indices for elements of LU

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... In-out argument:
      real bl(neq)     # on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     # error flag (0 OK, >0 try again, <0 failure)

c ... Local variables:
      logical usingsu
      real su(1)       # array that is not used when usingsu is false

      usingsu = .false.
      su(1) = 0.       # set unused array so flint will be happy
      call psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine psolvd (neq, t, yl, f0, wk, hrl1, wp, iwp, bl, lr,
     .                   ierr, rpar, ipar)

c ... Interface between linear-system solver vodpk and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real t           # physical time (not used here)
      real yl(neq)     # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)
      real hrl1        # scalar (not used for lr = 1)
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # dimensions and array indices for elements of LU
      integer lr       # type of calc. to be done (only 1 set up here)
      real rpar(*)     # real parameter communication (not used)
      integer ipar(*)  # integer parameter communication (not used)

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... In-out argument:
      real bl(neq)     # on input, c of P*x=c; on output, x

c ... Output argument:
      integer ierr     # error flag (0 OK, >0 try again, <0 failure)

c ... Local variables:
      logical usingsu
      real su(1)       # array that is not used when usingsu is false

      usingsu = .false.
      su(1) = 0.       # set unused array so flint will be happy
      call psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine FPSOL (neq, t, yl, f0, wk, hrl1, rewt, delta, nfe, bl, 
     .                  lr, zl, ierr)

c ... Interface between linear-system solver pvode and subroutine
c     psolbody.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real t           # physical time (not used here)
      real yl(neq)     # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)
      real hrl1        # scalar (not used for lr = 1)
      real rewt(*)     # reciprocal of error vector
      real delta       # input for iterative method (not used)
      integer nfe      # input/output for number of RHS evaluations
      integer lr       # type of calc. to be done (only 1 set up here)
      real bl(neq)     # RHS of P*zl=bl; P is preconditioning matrix

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... Output argument:
      real zl(neq)     # solution of P*zl=bl
      integer ierr     # error flag (0 OK, >0 try again, <0 failure)

c ... Common block:
      Use(Jac_work_arrays) # iwwp, wwp

c ... Local variables:
      logical usingsu
      real su(1)       # array that is not used when usingsu is false

      usingsu = .false.
      su(1) = 0.       # set unused array so flint will be happy
      call scopy (neq, bl, 1, zl, 1)
      call psolbody (neq, usingsu, su, wk, wwp, iwwp, zl, ierr)

      return
      end
c  **  End of subroutine FPSOL *********
c-----------------------------------------------------------------------
      subroutine psolbody (neq, usingsu, su, wk, wp, iwp, bl, ierr)

c ... Solve the linear system P*x=c, using elements of P loaded into
c     array wp.  Loading was done by the following subroutines:
c
c       caller   loading by
c       ------   ----------
c       psold      jacd2
c       psolnk     psetnk
c       psolnw     jacnw
c       psolvd     jacvd
c       fpsol      fpreco

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      logical usingsu  # .true. if su is used (svrpkg = "nksol" only)
      real su(neq)     # scale factors for yl
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # dimensions and array indices for elements of LU

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... In-out argument:
      real bl(neq)     # on input, c of P*x=c; on output, x

c ... Output arguments:
      integer ierr     # error flag (0 OK, >0 try again, <0 failure)

c ... Common blocks:
      Use(Timing)   # ttmatsol
      Use(Math_problem_size)   # neqmx(for arrays in Jacaux not used here)
      Use(Jacaux)            # isrnorm,fnormnw,jscalcol
      Use(Preconditioning)   # premeth
      Use(Jacreorder)        # perm,qperm,ireorder
      Use(Dim)               # nisp,ngsp
      Use(UEpar)             # svrpkg (used to reset usingsu for daspk)

c ... Function:
      real(Size4) gettime

c ... Local variables:
      integer i
      integer lowd, lbw, ubw
      integer ndiag, ndiagm
      real(Size4) sec4
      real tsmatsol

c ... Get initial value of system cpu timer.
      tsmatsol = gettime(sec4)

c ... Scale c by multiplying by row-normalization factors, if used,
c     and by column scaling vector su (for svrpkg='nksol' only).
      if (isrnorm .eq. 1) then
         do i = 1, neq
            bl(i) = bl(i) * fnormnw(i)
         enddo
      endif
      if (usingsu) then    # omit for daspk column scaling - correct?
         do i = 1, neq
            bl(i) = bl(i) * su(i)
         enddo
      endif

c ... Solve P*x=c for a preconditioner stored as a banded matrix.
      if (premeth .eq. 'banded') then
         lowd = iwp(1)
         lbw = iwp(2)
         ubw = iwp(3)
         call sgbsl (wp, lowd, neq, lbw, ubw, iwp(4), bl, 0)
         call scopy (neq, bl, 1, wk, 1)

c ... Solve P*x=c for a preconditioner stored as a sparse matrix in
c     compressed sparse row format.
c     If rows and columns of P were reordered (permuted), permute c,
c     then use inverse permutation on x.
      elseif (premeth .eq. 'ilut') then
	 if (ireorder .eq. 1) call dvperm (neq, bl, perm)
         call lusol0 (neq, bl, wk, wp, iwp(neq+1), iwp)
	 if (ireorder .eq. 1) call dvperm (neq, wk, qperm)

c ... Solve P*x=c for a preconditioner stored as a sparse matrix in
c     diagonal storage format.
      else
         ndiag = iwp(1)
         ndiagm = iwp(2)
         call minvmul (neq, ndiag, ndiagm, wp, iwp(3), wk, bl)
      endif

c ... Divide solution x by column-scaling factors (for svrpkg='nksol').
      if (usingsu) then
         do i = 1, neq
            bl(i) = wk(i) / su(i)
         enddo
      elseif (svrpkg.eq."daspk" .and. jscalcol.eq.1) then
         do i = 1, neq
            bl(i) = wk(i) / su(i)
         enddo
      elseif (premeth .ne. 'banded') then
         call scopy (neq, wk, 1, bl, 1)
      endif

c ... Accumulate cpu time spent here.
      ierr = 0
      ttmatsol = ttmatsol + (gettime(sec4) - tsmatsol)
      return
      end
c-----------------------------------------------------------------------
      subroutine sfsetnk (neq, yl, su, sf)

c ... Calculate the scaling vector of function values for the nksol
c     routine.  The Jacobian at yl is calculated, and sf(i) is set
c     to the reciprocal of the norm of row i.

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real yl(*)       # dependent variables
      real su(neq)     # scale factors for yl

c ... Output arguments:
      real sf(neq)     # scaling vector (after use as work array in
                       # jac_calc and here)

c ... Common blocks:
      Use(Decomp)              # lbw,ubw
      Use(Jacobian)            # nnzmx,jac,jacj,jaci
      Use(Math_problem_size)   # neqmx
      Use(Jacaux)              # yldot0,normtype
      Use(Dim)                 # nx,ny
      Use(Time_dep_nwt)        # ydt_max0,nufak0
      Use(Share)               # cutlo

c ... Local variables:
      real tp
      integer i
      #Add working array to avoid recursive variable assignment
      real wk(neq)
      wk(1:neq)=0.0
c ... Flag these calls to RHS (pandf) as Jacobian calculations
      yl(neq+1) = 1.  # with ccc, dont include nufak in scaling Jacobian

c ... Calculate right-hand sides at unperturbed values of variables.
      call rhsnk (neq, yl, yldot0)

c ... Calculate Jacobian matrix.
      tp = 0.
      #Working array wk in place of sf in the call to jac_calc (J.Guterl)
      call jac_calc_interface (neq, tp, yl, yldot0, lbw, ubw, wk,
     .               nnzmx, jac, jacj, jaci)

      yl(neq+1) = -1.      # Turn-off Jacobian flag for pandf
c ... Compute inverse of column-scaling vector, and perform column
c     scaling.  
      do i = 1, neq
         sf(i) = 1. / su(i)
      enddo
      call amudia (neq, 0, jac, jacj, jaci, sf, jac, jacj, jaci)

c ... Calculate one of three types of norm for each row.
c ... Also find initial maximum of yldot*sf = ydt_max0 for scaling nufak
      call rnrms (neq, normtype, jac, jacj, jaci, sf)
      nufak0 = nufak                    # record initial nufak value
      ydt_max0 = cutlo
      do i = 1, neq
         if (abs(sf(i)) .lt. 1e20*cutlo) then
            write(*,*) '*** Error: Jacobian row = 0 for eqn iv =', i
            call xerrab("")
         else
            sf(i) = 1./sf(i)
         endif
         if (abs(yldot0(i)*sf(i)) .gt. ydt_max0)
     .                            ydt_max0 = abs(yldot0(i)*sf(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_dt(neq, yl, f0)

c ... Calculates the time step to be used with svrpkg="nksol" based
c ... on dtreal and various cases of (yl/yldot)

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real yl(*)       # most recent iterate of solution vector
      real f0(neq)     # function values f(yl)

c ... Common blocks:
      Use(Math_problem_size)   # neqmx
      Use(Dim)            # nx,ny
      Use(Time_dep_nwt)   # nufak,ydt_max,ydt_max0,alfnuf,expnuf,nufak0
                          # inufaknk,dtoptx,dtoptv
      Use(Indexes)        # idxn,idxu,idxte,idxti,idxng,idxphi,iseqalg
      Use(UEpar)          # isbcwdt,isnionxy,isuponxy,isteonxy,istionxy,
                          # isngonxy,isphionxy
      Use(Share)          # geometry,nxc,isnonog,cutlo
      Use(Selec)          # ixm1,ixp1
      Use(Indices_domain_dcl)    # ixmnbcl,ixmxbcl,iymnbcl,iymxbcl

c ... Local variables:
      integer iv,iv1,iv2,iv3,iv4,ifld,igsp,ix,iy,iym1,iyp1,ixm1u,ixp1u
      real up_5ca


      call rhsnk (neq, yl, f0)    # Reset f0 with nufak off

c ... special new section for adjustable timesteps
      do iy = 1-iymnbcl, ny+iymxbcl
         iym1 = max(1-iymnbcl,iy-1)
	 iyp1 = min(ny+iymxbcl,iy+1)
         do ix = 1-ixmnbcl, nx+ixmxbcl
            do ifld = 1, nisp
               if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
                  dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
                  if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                  endif
               endif
            enddo

            if(ix.ne.nx+2*isbcwdt) then  
                           # nx test - for algebr. eq. unless isbcwdt=1
               do ifld = 1, nusp
                  if(isuponxy(ix,iy,ifld).eq.1) then
                    ixm1u = max(1-ixmnbcl, ixm1(ix,iy))
                    ixp1u = min(nx+ixmxbcl, ixp1(ix,iy))
                    iv = idxu(ix,iy,ifld)
                    iv1 = idxu(ixm1u,iy,ifld)
                    iv2 = idxu(ixp1u,iy,ifld)
                    iv3 = idxu(ix,iyp1,ifld)
                    iv4 = idxu(ix,iym1,ifld)
                    up_5ca = ( abs(ylodt(iv)) + abs(ylodt(iv1)) +
     .                       abs(ylodt(iv2))+ abs(ylodt(iv3)) +
     .                       abs(ylodt(iv4)) )/5
                    if (abs(f0(iv)).gt.cutlo) dtoptv(iv) = 
     .                                     deldt*abs(up_5ca/(f0(iv)))
                    if (model_dt .eq. 0) then
                      dtuse(iv) = dtreal
                    elseif (model_dt .eq. 1) then
                      dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                    elseif (model_dt .eq. 2) then
                      dtuse(iv) = dtoptv(iv)
                    elseif (model_dt .eq. 3) then
                      dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                    endif
                  endif
               enddo
            endif

            if(isteonxy(ix,iy).eq.1) then
               iv =  idxte(ix,iy)
               dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
               dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
               if (model_dt .eq. 0) then
                 dtuse(iv) = dtreal
               elseif (model_dt .eq. 1) then
                 dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
               elseif (model_dt .eq. 2) then
                 dtuse(iv) = dtoptv(iv)
               elseif (model_dt .eq. 3) then
                 dtuse(iv) = sqrt(dtreal*dtoptv(iv))
               endif
            endif

            if(istionxy(ix,iy).eq.1) then
               iv = idxti(ix,iy)
               dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
               dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
               if (model_dt .eq. 0) then
                 dtuse(iv) = dtreal
               elseif (model_dt .eq. 1) then
                 dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
               elseif (model_dt .eq. 2) then
                 dtuse(iv) = dtoptv(iv)
               elseif (model_dt .eq. 3) then
                 dtuse(iv) = sqrt(dtreal*dtoptv(iv))
               endif
            endif

            do igsp = 1, ngsp
               if(isngonxy(ix,iy,igsp).eq.1) then
                  iv = idxg(ix,iy,igsp)
                  dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
                  dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
                  if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                  endif
               endif
            enddo

            do igsp = 1, ngsp
               if(istgonxy(ix,iy,igsp).eq.1) then
                  iv = idxtg(ix,iy,igsp)
                  dtoptv(iv) = deldt*abs(ylodt(iv)/(f0(iv)+cutlo))
                  dtoptx(ix,iy) = min(dtoptv(iv), dtoptx(ix,iy))
                  if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
                  endif
               endif
            enddo

            if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
               iv = idxphi(ix,iy)
               dtoptv(iv) = dtoptv(idxte(ix,iy))  # same as for Te
               if (model_dt .eq. 0) then
                    dtuse(iv) = dtreal
                  elseif (model_dt .eq. 1) then
                    dtuse(iv) = dtreal*dtoptv(iv)/(dtreal+dtoptv(iv))
                  elseif (model_dt .eq. 2) then
                    dtuse(iv) = dtoptv(iv)
                  elseif (model_dt .eq. 3) then
                    dtuse(iv) = sqrt(dtreal*dtoptv(iv))
               endif
            endif
         enddo
      enddo
ccc
c ... If model_dt < 4, then jump over this to 23; otherwise use this
c ... section to define time-step based on cell minimum-step, dtoptx
ccc
      if (model_dt .lt. 4) goto 23
      do iy = 1-iymnbcl, ny+iymxbcl 
         do ix = 1-ixmnbcl, nx+ixmxbcl
            do ifld = 1, nisp
               if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  if (model_dt .eq. 4) then
                    dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                          (dtreal+dtoptx(ix,iy))
                  elseif (model_dt .eq. 5) then
                    dtuse(iv) = dtoptx(ix,iy)
                  elseif (model_dt .eq. 6) then
                    dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
                  endif
               endif
            enddo
            if(ix.ne.nx+2*isbcwdt) then  
                           # nx test - for algebr. eq. unless isbcwdt=1
               do ifld = 1, nusp
                  if(isuponxy(ix,iy,ifld).eq.1) then

                  iv = idxu(ix,iy,ifld)
                  if (model_dt .eq. 4) then
                    dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                          (dtreal+dtoptx(ix,iy))
                  elseif (model_dt .eq. 5) then
                    dtuse(iv) = dtoptx(ix,iy)
                  elseif (model_dt .eq. 6) then
                    dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
                  endif

                  endif
               enddo
            endif
            if(isteonxy(ix,iy).eq.1) then
               iv =  idxte(ix,iy)
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            endif
            if(istionxy(ix,iy).eq.1) then
            endif
            do igsp = 1, ngsp
               if(isngonxy(ix,iy,igsp).eq.1) then
                  iv = idxg(ix,iy,igsp)
               endif
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            enddo
            do igsp = 1, ngsp
               if(istgonxy(ix,iy,igsp).eq.1) then
                  iv = idxtg(ix,iy,igsp)
               endif
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            enddo
            if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
               iv = idxphi(ix,iy)
               if (model_dt .eq. 4) then
                 dtuse(iv) = dtreal*dtoptx(ix,iy)/
     .                       (dtreal+dtoptx(ix,iy))
               elseif (model_dt .eq. 5) then
                 dtuse(iv) = dtoptx(ix,iy)
               elseif (model_dt .eq. 6) then
                 dtuse(iv) = sqrt(dtreal*dtoptx(ix,iy))
               endif
            endif
         enddo
      enddo
 23   continue

c ... Be sure that algebraic equations (esp. limiter) have large dt
      if (isbcwdt .eq. 0) then
        do iv = 1, neq
          if (iseqalg(iv).eq.1) dtuse(iv) = 1.e20
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------

      subroutine jacmap

      implicit none
c_mpi      include 'mpif.h'

c ... Output a map of the Jacobian matrix.

c ... All information is passed through common for convenience in
c     calling this routine from the UEDGE> prompt.

c ... Common blocks:
      Use(Math_problem_size)   # neqmx(for arrays in Lsode not used here)
      Use(Lsode)           # neq
      Use(Jacobian)        # jac,jacj,jaci
      Use(Jacobian_full)   # jacfull
      Use(Jacreorder)     # ireorder

c ... Local variables:
      integer ierr
      integer us
c_mpi      integer my_pe
      character*24 filename

c ... Allocate full Jacobian for jacmap; warning of size
      call remark("*** CAUTION: allocating large jacfull(neq,neq)***")
      call gallot("Jacobian_full",0)
      write (STDOUT,*) '*** Full Jacobian size is neq**2 = ', neq*neq

c ... Issue a warning if reordering is on which may rearrange Jacobian
      if (ireorder .eq. 1) then
         write (STDOUT,*) '***ireorder=1, Jacobian may be rearranged***'
      endif

c ... Convert Jacobian matrix to full storage format.
      call csrdns (neq, neq, jac, jacj, jaci, jacfull, neq, ierr)
      if (ierr .ne. 0) then
         write (STDOUT,*)
     .      '*** jacmap got error return ierr =', ierr,
     .      ' from csrdns.'
         call xerrab("")
      endif

c ... Open a file, and output the map.
      call freeus (us)
      filename = 'Jacobian_map.dat'
c_mpi      if(MY_PE().eq.0) then
c_mpi        us = 59
c_mpi        filename = 'Jacobian_map.dat0'
c_mpi      else
c_mpi        us = 69
c_mpi        filename = 'Jacobian_map.dat1'
c_mpi      endif
      open(unit=us, file=filename, status='unknown')
      call jmap (neq, jacfull, us)

c ... Close file, and report file name.
      close(us)
      write (STDOUT,*) ' Jacobian map in data file:  ', filename

      return
      end
c-----------------------------------------------------------------------
      subroutine jacout

      implicit none

c ... Output Jacobian matrix and right-hand side in Boeing-Harwell
c     format.

c ... All information is passed through common for convenience in
c     calling this routine from the UEDGE> prompt.

c ... Common blocks:
      Use(Dim)        # nusp(for array fnorm in Ynorm not used here)
      Use(Math_problem_size)   # neqmx
      Use(Lsode)      # neq,yldot
      Use(Ynorm)      # sfscal
      Use(Jacobian)   # jac,jacj,jaci
      Use(UEpar)      # svrpkg

c ... Local variables:
      integer i, us, ifmt
      character*24 filename
      character*72 title

c ... For the nksol solver, scale the right-hand side.
      if ((svrpkg .eq. 'nksol') .or. (svrpkg.eq.'petsc')) then
         do i = 1, neq
            yldot(i) = yldot(i) * sfscal(i)
         enddo
      endif

c ... Open a file, and output data.
      call freeus (us)
      filename = 'Uedge_Test_Matrix.dat'
      open(unit=us, file=filename, status='unknown')
      title = ' UEDGE Test Matrix '
      ifmt = 15
      call prtmt (neq,neq,jac,jacj,jaci,yldot,
     .   'NN',title,'SPARSKIT','RUA',ifmt,3,us)

c ... Close file, and report file name.
      close(us)
      write (STDOUT,*) ' Jacobian matrix in data file:  ', filename

      return
      end
******* end of subroutine jacout ******* 
c-----------------------------------------------------------------------
      subroutine engbal(pwrin)

      implicit none

c ... Calculates various components of the 2-D energy flow and the 
c ... ionization and radiation for use in the postprocessing file
c ... balancee to determine energy balance; these 2-D loops become 
c ... expensive when done from the parser.

c ... Input arguments:
      real pwrin            #total input power for normalization

c ... Common blocks:
      Use(Dim)              #nx,ny,nzspt,nzsp,ngsp,nfsp,nxpt
      Use(Xpoint_indices)   #ixlb,ixrb
      Use(Comgeo)           #gx,gy,vol,sx,sy
      Use(Selec)            #ixm1,ixp1
      Use(Comflo)           #feex,feix,feey,feiy,fqx,fqy,fnix,fniy
      Use(Compla)           #up,mi,ng,v2,vy
      Use(UEpar)            #ediss,eion,ebind,iigsp,ishosor,fsprd,ishymol
      Use(Phyvar)           #ev,qe
      Use(Share)            #cutlo
      Use(Rhsides)          #psor,psorg,psorxr,erliz,erlrc,psor_tmpov
      Use(Coefeq)           #cfvisx,cfvisy,cnsor,cngmom
      Use(Bcond)            #ckinfl
      Use(Parallv)          # nxg,nyg
      Use(Conduc)           #visx,visy,eeli,nuiz,nucx
      Use(Imprad)           #prad,pradz
      Use(Postproc)         #complete group
      Use(Indices_domain_dcl)  # ixmnbcl,ixmxbcl
      Use(Noggeo)           # angfx

c ... Local variables:
      integer ix,iy,ix1,ix2,ix3,ix4,iimp,id,igsp,jz,jx,ixt,ixt1,ixr,ixr1
      real ave,a1,a2,t1,t2,thetaix,thetaix2,eta_dup2dy

c ... Implicit function:
      ave(a1,a2) = a1 * a2 / (a1 + a2 + cutlo)
   

      if (ishosor .eq. 1) then  # averge power terms is ishosor=1
         call volavenv(nx, ny, 1, ny, 1, nx, ixp1(0:nx+1,0:ny+1), ixm1(0:nx+1,0:ny+1),
     .                               fsprd, psor_tmpov(0:nx+1,0:ny+1), prad)
         do igsp = nhgsp+1, ngsp
           jz = igsp - nhgsp
           do iimp = 0, nzsp(jz)       
             call volavenv(nx, ny, 1, ny, 1, nx, ixp1(0:nx+1,0:ny+1), ixm1(0:nx+1,0:ny+1),
     .                     fsprd, psor_tmpov(0:nx+1,0:ny+1), pradz(0:nx+1,0:ny+1,iimp,jz))
           enddo
        enddo
      endif

# Set arrays to check energy conserv; add ion parallel drift and visc heat

      do jx = 1, nxpt
        do ix=ixlb(jx),ixrb(jx)
          do iy=0,ny
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            fety(ix,iy) = 0.
            fetx(ix,iy) = 0.
            do id = 1, nusp
	       thetaix =  0.5*(angfx(ix1,iy) + angfx(ix,iy))
               thetaix2 = 0.5*(angfx(ix,iy) + angfx(ix2,iy))
               eta_dup2dy = 0.25*( visy(ix,iy+1,id)*
     .                       (upi(ix1,iy+1,id)+upi(ix,iy+1,id))**2 -
     .                             visy(ix,iy  ,id)*
     .                       (upi(ix1,iy  ,id)+upi(ix,iy  ,id))**2 )
               fety(ix,iy) = fety(ix,iy) + (mi(id)/32)*( upi(ix1,iy,id)+
     .                         upi(ix,iy,id)+upi(ix1,iy+1,id)+
     .                         upi(ix,iy+1,id) )**2*fniy(ix,iy,id) -
     .                         cfvisy*0.5*sy(ix,iy)*gyf(ix,iy)*eta_dup2dy
               fetx(ix,iy) = fetx(ix,iy) + 0.5*mi(id)*upi(ix,iy,id)**2*
     .                          fnix(ix,iy,id) - cfvisx*0.25*sx(ix,iy)*(
     .                         visx(ix ,iy,id)*gx(ix ,iy)*cos(thetaix)* 
     .                          ( upi(ix,iy,id)**2 - upi(ix1,iy,id)**2 ) +
     .                        visx(ix2,iy,id)*gx(ix2,iy)*cos(thetaix2)* 
     .                          ( upi(ix2,iy,id)**2 - upi(ix,iy,id)**2 ) )
               fetx(ix,iy) = fetx(ix,iy) - upi(ix,iy,id)*fmixy(ix,iy,id)
            enddo
            fety(ix,iy) = fety(ix,iy) + feey(ix,iy) + feiy(ix,iy)
            fetx(ix,iy) = fetx(ix,iy) + feex(ix,iy) + feix(ix,iy)
          enddo
        enddo
      enddo

# Now correct the boundary x-fluxes if non-unity ckinfl
      if (abs(ckinfl-1.) > 1.e-10) then
       do jx = 1, nxpt
         ixt  = ixlb(jx)
         ixt1 = ixt + 1
         ixr  = ixrb(jx)
         ixr1 = ixr - 1
         do 15 iy = 0, ny
            fetx(ixt,iy) = 0.
            fetx(ixr,iy) = 0.
            do id = 1, nfsp
               fetx(ixt,iy) = fetx(ixt,iy) +
     .                        0.5*mi(id)*up(ixt,iy,id)**2*fnix(ixt,iy,id) -
     .                   ckinfl*0.5*sx(ixt,iy)*visx(ixt1,iy,id)*gx(ixt1,iy)*
     .                          ( up(ixt1,iy,id)**2 - up(ixt,iy,id)**2 )
               fetx(ixr,iy) = fetx(ixr,iy) +
     .                        0.5*mi(id)*up(ixr,iy,id)**2*fnix(ixr,iy,id) - 
     .                   ckinfl*0.5*sx(ixr,iy)*visx(ixr,iy,id)*gx(ixr,iy)*
     .                          ( up(ixr,iy,id)**2 - up(ixr1,iy,id)**2 )
            enddo
            fetx(ixt,iy) = fetx(ixt,iy) + feex(ixt,iy) + feix(ixt,iy)
            fetx(ixr,iy) = fetx(ixr,iy) + feex(ixr,iy) + feix(ixr,iy)
 15      continue
       enddo  # end do-loop over nxpt mesh regions
      endif   # test on ckinfl-1

      pvmomcx = 0.e0
      ptjdote = 0.e0 

      do jx = 1, nxpt
        do ix=ixlb(jx)+1,ixrb(jx)
          do iy=1,ny
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
# Here peirad includes sum of electron and ion energy losses; note that binding
# energy is included in eeli term, but it is carried by the ions.
# Note also that eion and ediss generally balance in the next line
# because should have ediss=2*eion - transfer from electron to ion energy
            peirad(ix,iy) = cnsor*( erliz(ix,iy) + erlrc(ix,iy) +
     .                              ebind*ev*psor(ix,iy,1) -
     .                              ebind*ev*psorrg(ix,iy,1) +
     .                              ediss*ev*(0.5*psordis(ix,iy)) -
     .                        ceisor*eion*ev*(psordis(ix,iy)) )

# other energy diagnostics are given below

cc            jdote(ix,iy) = -   # this energy is included in resee, not lost
cc     .                  0.5 * fqx(ix ,iy)*(phi(ix2,iy  )+phi(ix ,iy)) +
cc     .                  0.5 * fqx(ix1,iy)*(phi(ix ,iy  )+phi(ix1,iy)) -
cc     .                  0.5 * fqy(ix ,iy)*(phi(ix ,iy+1)+phi(ix ,iy)) +
cc     .                  0.5 * fqy(ix,iy-1)*(phi(ix,iy)+phi(ix,iy-1))
cc            ptjdote = ptjdote + jdote(ix,iy)
            ptjdote = ptjdote + wjdote(ix,iy)

            if (isupgon(1) .eq. 0) then
               pmomv(ix,iy)=cngmom(1)*up(ix,iy,1)*sx(ix,iy)*rrv(ix,iy)*
     .                           ( ng(ix2,iy,1)*tg(ix2,iy,1)- 
     .                             ng(ix ,iy,1)*tg(ix ,iy,1) ) +
     .             cmwall(1)*0.125*mi(1)*(up(ix,iy,1)+up(ix1,iy,1))**2*
     .                ng(ix,iy,1)*nucx(ix,iy,1)*vol(ix,iy)
            elseif (isupgon(1) .eq. 1) then    # inertial neutrals
               pmomv(ix,iy) = 0.  # coupled back to therm eng for inertial neut
            endif
            pvmomcx = pvmomcx + pmomv(ix,iy)

            engerr(ix,iy) = ( fetx(ix1,iy  )-fetx(ix,iy)+
     .                        fety(ix ,iy-1)-fety(ix,iy)-
     .                        peirad(ix,iy)-png2ni(ix,iy) ) /
     .                         abs(pwrin)
            if (isimpon.ne.0) then  # prad allocated only if isimpon.ne.0
               engerr(ix,iy) = engerr(ix,iy) - prad(ix,iy)*vol(ix,iy)/
     .                                                     abs(pwrin)
            endif
          enddo
        enddo
      enddo

# ionization and background sources

      iion_tot = 0.e0
      irecomb_tot = 0.e0
      icxgas_tot = 0.e0
      pradrc = 0.e0
      pradiz = 0.e0
      pradht = 0.e0
      prdiss = 0.e0
      pibirth = 0.e0
      pbinde = 0.e0
      pbindrc = 0.e0
      pradzbind = 0.e0
      do igsp = 1, ngsp
         iion(igsp) = 0.e0
         irecomb(igsp) = 0.e0
         icxgas(igsp) = 0.e0
      enddo
      do igsp = 1, max(1, ngsp-nhgsp)
         pradimpt(igsp) = 0.
         if (nzsp(igsp) .ne. 0) then
            do iimp = 0, nzsp(igsp)
               pradimp(iimp,igsp) = 0.
            enddo
         endif
      enddo    
      pradfft = 0.
      
      do jx = 1, nxpt
        do ix = ixlb(jx)+1,ixrb(jx)
          do iy = 1, ny
            do igsp = 1, ngsp
              if (ishymol.eq.0 .or. igsp.ne.2) then
               iion(igsp) = iion(igsp) - cnsor*qe*psorg(ix,iy,igsp)
               irecomb(igsp) = irecomb(igsp) -cnsor*qe*psorrg(ix,iy,igsp)
               icxgas(igsp) = icxgas(igsp) - qe*psorcxg(ix,iy,igsp)
              endif
            enddo
          enddo
        enddo
      enddo

      do igsp = 1, ngsp
         iion_tot = iion_tot + iion(igsp)
         irecomb_tot = irecomb_tot + irecomb(igsp)
         icxgas_tot = icxgas_tot + icxgas(igsp)
      enddo

      do iy = 1, ny
        do jx = 1, nxpt
          do ix = ixlb(jx)+1, ixrb(jx)
            ix1 = ixm1(ix,iy)
            pradrc = pradrc + cnsor*erlrc(ix,iy)
            pradiz = pradiz + (eeli(ix,iy)-ebind*ev) * psor(ix,iy,1) 
            pbinde = pbinde + ebind*ev * psor(ix,iy,1)
            pbindrc = pbindrc + ebind*ev*psorrg(ix,iy,1)
            prdiss = prdiss + ediss * ev * (0.5*psordis(ix,iy))
            pibirth = pibirth + ceisor* eion*ev * (psordis(ix,iy)) -
     .                ccoldsor*ng(ix,iy,1)*(1.5*ti(ix,iy)-eion*ev)*
     .                                          nucx(ix,iy,1)*vol(ix,iy) 
          enddo
        enddo
      enddo
      pradht = pradiz + pradrc
   
      if (isimpon .eq. 2 .or. isimpon .eq. 7) then #fixed fraction model
         do iy = 1, ny
           do jx = 1, nxpt
             do ix = ixlb(jx)+1, ixrb(jx)
               pradfft = pradfft + pradcff(ix,iy)*vol(ix,iy)
             enddo
           enddo
         enddo
      endif

      if (isimpon .gt. 2) then  #separate impurity species
         do iy = 1, ny
           do jx = 1, nxpt
             do ix = ixlb(jx)+1, ixrb(jx)
               do igsp = nhgsp+1, ngsp
                  jz = igsp - nhgsp
                  do iimp = 0, nzsp(jz)
                     pradimp(iimp,jz) = pradimp(iimp,jz) +
     .                             pradz(ix,iy,iimp,jz)*vol(ix,iy)
                  enddo
               enddo
               pradzbind = pradzbind + (pwrze(ix,iy)-prad(ix,iy))*
     .                                vol(ix,iy) # only total pradzbind calc
             enddo
           enddo
         enddo
         do igsp = nhgsp+1, ngsp
            jz = igsp - nhgsp
            do iimp = 0, nzsp(jz)
               pradimpt(jz) = pradimpt(jz) + pradimp(iimp,jz)
            enddo
         enddo
      endif
  
      return
      end          

******* end of subroutine engbal *******

c----------------------------------------------------------------------c
      subroutine pradpltwl

      Implicit none

c ... Calc radiation power to divertor and outer wall surfaces
c ... Use as a diagnostic to call from the BASIS parser

      Use(Dim)            # nx,ny
      Use(Share)          # nxomit
      Use(Imprad)         # isimpon, prad
      Use(Comgeo)         # sx, sy, vol
      Use(Noggeo)         # angfx
      Use(RZ_grid_info)   # rm, zm
      Use(Phyvar)         # pi, ev
      Use(Xpoint_indices) # ixlb, ixrb
      Use(Conduc)         # eeli
      Use(Rhsides)        # psor
      Use(UEpar)          # ebind
      Use(Postproc)	  #pwr_pltz,pwr_plth,pwr_wallz,pwr_wallh
                          #pwr_pfwallz, pwr_pfwallh

c ... Local variables
      real prdu(0:nx+1,0:ny+1)   
      real theta_ray1, theta_ray2, dthgy, dthgx, sxo, frth
      integer ixv, iyv, nj, ix, iy, ip, iodd

# Initialize arrays
      call sfill ((ny+2)*2*nxpt, 0., pwr_pltz, 1)   
      call sfill ((ny+2)*2*nxpt, 0., pwr_plth, 1)   
      call sfill ((nx+2), 0., pwr_wallz, 1)   
      call sfill ((nx+2), 0., pwr_wallh, 1) 
      call sfill ((nx+2)*nxpt, 0., pwr_pfwallz, 1)   
      call sfill ((nx+2)*nxpt, 0., pwr_pfwallh, 1) 
  
      if (isimpon > 0) then   # use prdu since prad might not be dimensioned
         call s2copy (nx+2,ny+2, prad,1,nx+2, prdu,1,nx+2) #prad --> prdu
      else
         call s2fill (nx+2, ny+2, 0., prdu, 1, nx+2)
      endif

      nj = nxomit

c ... First do the divertor plate surfaces
      do ip = 1, 2*nxpt
        iodd = (ip+1)/2 - ip/2  # =1 if ip odd, =0 if ip even
	if (iodd==1) then
          ixv = ixlb(ip/2+1)    # viewing ix position
	else
	  ixv = ixrb(ip/2) + 1
	endif
        do iyv = 1, ny		# loop over viewing ix
          do iy = 1, ny	        # loop over source iy
            do ix = 1, nx	# loop over source ix
              theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0), 
     .                            rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
              theta_ray2 = atan2( zm(ixv+nj,iyv,3)-zm(ix+nj,iy,0), 
     .                            rm(ixv+nj,iyv,3)-rm(ix+nj,iy,0) )
              dthgy = abs(theta_ray1-theta_ray2)
              frth = min(dthgy, 2*pi-dthgy)/(2*pi)  # frac.; need angle < pi
              sxo = sx(ixv,iyv)/(cos(angfx(ixv,iyv)))
              pwr_pltz(iyv,ip) = pwr_pltz(iyv,ip) + 
     .                                 prdu(ix,iy)*vol(ix,iy)*frth/sxo
              pwr_plth(iyv,ip) = pwr_plth(iyv,ip) + (
     .                           (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1)
     .                                        + erlrc(ix,iy))*frth/sxo
            enddo
          enddo
        enddo
c ... Set corner values
      pwr_pltz(0,ip)    = pwr_pltz(1,ip)	
      pwr_pltz(ny+1,ip) = pwr_pltz(ny,ip)
      pwr_plth(0,ip)    = pwr_plth(1,ip)	
      pwr_plth(ny+1,ip) = pwr_plth(ny,ip)

      enddo             # end of ip loop over divertor plates


c ... Now do the "outer" wall surface, i.e., iy=ny+1
      iyv = ny+1                # viewing iy position
      do ixv = 1, nx		# loop over viewing ix
        do iy = 1, ny	        # loop over source iy
          do ix = 1, nx	        # loop over source ix
            theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
            theta_ray2 = atan2( zm(ixv+nj,iyv,2)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,2)-rm(ix+nj,iy,0) )
            dthgx = abs(theta_ray1-theta_ray2)             
            frth = min(dthgx, 2*pi-dthgx)/(2*pi)  # frac; need angle < pi
            pwr_wallz(ixv) = pwr_wallz(ixv) + 
     .                           prdu(ix,iy)*vol(ix,iy)*frth/sy(ixv,iyv) 
            pwr_wallh(ixv) = pwr_wallh(ixv) + (
     .                       (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1) +
     .                          erlrc(ix,iy) )*frth/sy(ixv,iyv) 
          enddo
        enddo
      enddo
      pwr_wallz(0) = pwr_wallz(1)	# Because prad(0,) is not calculated
      pwr_wallz(nx+1) = pwr_wallz(nx)
      pwr_wallh(0) = pwr_wallh(1)
      pwr_wallh(nx+1) = pwr_wallh(nx)

c ... Finally do the private-flux wall surfaces, i.e., iy=0
      iyv = 0                        # viewing iy position
      do ip = 1, nxpt                # loop over number of x-points
        do ixv = 1, nx   	     # loop over viewing ix
          do iy = 1, ny	             # loop over source iy
            do ix = 1, nx	     # loop over source ix
              theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
              theta_ray2 = atan2( zm(ixv+nj,iyv,2)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,2)-rm(ix+nj,iy,0) )
              dthgx = abs(theta_ray1-theta_ray2)             
              frth = min(dthgx, 2*pi-dthgx)/(2*pi)  # frac; need angle < pi
              pwr_pfwallz(ixv,ip) = pwr_pfwallz(ixv,ip) + 
     .                           prdu(ix,iy)*vol(ix,iy)*frth/sy(ixv,iyv) 
              pwr_pfwallh(ixv,ip) = pwr_pfwallh(ixv,ip) + (
     .                       (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1) +
     .                          erlrc(ix,iy) )*frth/sy(ixv,iyv) 
            enddo
          enddo
          if(ixv>ixpt1(ip) .and. ixv <ixpt2(ip)+1) then  # 0 in core region
            pwr_pfwallh(ixv,ip) = 0.
            pwr_pfwallh(ixv,ip) = 0.
            pwr_pfwallz(ixv,ip) = 0.
            pwr_pfwallz(ixv,ip) = 0.
          endif
        enddo
        pwr_pfwallz(0,ip) = pwr_pfwallz(1,ip)  # prad(0,) is not calculated
        pwr_pfwallz(nx+1,ip) = pwr_pfwallz(nx,ip)
        pwr_pfwallh(0,ip) = pwr_pfwallh(1,ip)
        pwr_pfwallh(nx+1,ip) = pwr_pfwallh(nx,ip)
      enddo

      return
      end
******* end of subroutine pradpltwl *******
c-------------------------------------------------------------------------c

      subroutine plateflux

      Implicit none

c ... Calc particle & heat fluxes to divertor plates
c ... Use as diagnostic called from BASIS/Python parser (alt. to balancee)

      Use(Dim)            # nx,ny
      Use(Comgeo)         # sx,vol,gx,gxf
      Use(Noggeo)         # angfx
      Use(Phyvar)         # pi, ev
      Use(Postproc)       # pwr_plth,Pwr_pltz,sdel,rb;sdil,rb; sdbindl,rb,
                          # sdtl,rb;gdil,rd
      Use(Xpoint_indices) # ixlb, ixrb
      Use(Compla)         # mi,ni,up,te,ti
      Use(Comflo)         # fnix,feex,feix
      Use(Conduc)         # visx,hcxn
      Use(UEpar)          # ebind
      Use(Bcond)          # ckinfl
      Use(Parallv)        # nxg,nyg
      Use(Poten)          # phi0l,phi0r

c  Local variables
      integer iu,jx,id,ixi,ixo 
      #Former Aux module variables
      integer ix,iy
      real tempvo,tempvi
      real pdivilb(1:nxpt),pdivirb(1:nxpt),pdivelb(1:nxpt),
     .     pdiverb(1:nxpt)
      real pbindlb(1:nxpt),pbindrb(1:nxpt),pdivnlb(1:nxpt),
     .     pdivnrb(1:nxpt)
      real sdilbd(0:ny+1,1:nfsp,1:nxpt),sdirbd(0:ny+1,1:nfsp,1:nxpt)
      real sdnlb(0:ny+1,1:nxpt),sdnrb(0:ny+1,1:nxpt)
      real sxi(0:ny+1,1:nxpt),sxo(0:ny+1,1:nxpt)

########################################
# First do the particle flux
##############################################################
      do jx=1,nxpt
        ixi=ixlb(jx)	# ixi=0
        ixo=ixrb(jx)	# ixo=nx
        do iy=1,ny+1
          sxo(iy,jx) = sx(ixo,iy)/(cos(angfx(ixo,iy)))
          sxi(iy,jx) = sx(ixi,iy)/(cos(angfx(ixi,iy)))
          do id = 1, nfsp
	    gdilb(iy,id,jx) = -fnix(ixi,iy,id)/sxi(iy,jx)
	    gdirb(iy,id,jx) =  fnix(ixo,iy,id)/sxo(iy,jx)
            engilb(iy,id,jx) = (2.*ti(ixi,iy)/ev + zi(id)*
     .                          (phi(ixi,iy)-phi0l(iy,jx)) )
            engirb(iy,id,jx) = 2.*ti(ixo+1,iy)/ev + zi(id)*
     .                          (phi(ixo+1,iy)-phi0r(iy,jx))
          enddo
        enddo
      enddo

# Fix corner boundary values
      do jx = 1, nxpt
        do id = 1, nfsp
          gdilb(0,id,jx) =     gdilb(1,id,jx)
          gdilb(ny+1,id,jx) =  gdilb(ny,id,jx)
          gdirb(0,id,jx) =     gdirb(1,id,jx)
          gdirb(ny+1,id,jx) =  gdirb(ny,id,jx)
          engilb(0,id,jx) =    engilb(1,id,jx)
          engilb(ny+1,id,jx) = engilb(ny,id,jx)
          engirb(0,id,jx) =    engirb(1,id,jx)
          engirb(ny+1,id,jx) = engirb(ny,id,jx)
        enddo
      enddo

#######################################
# Now do the heat flux
##############################################################
#

# Get radiation power in pwr_plth and pwr_plt
      call pradpltwl

# Fill radiation heat flux arrays and initialize tot powers
      do jx = 1, nxpt
        pdivirb(jx) = 0.
        pdiverb(jx) = 0. 
        pdivilb(jx) = 0. 
        pdivelb(jx) = 0. 
        pbindrb(jx) = 0. 
        pbindlb(jx) = 0. 
        do iy = 0, ny+1
          iu = 2*(jx/2)+1  # iu/iu+1 gives "odd/even" plates (in/out)
          sdrlb(iy,jx) = pwr_plth(iy,iu)+pwr_pltz(iy,iu)
          sdrrb(iy,jx) = pwr_plth(iy,iu+1)+pwr_pltz(iy,iu+1)
        enddo
      enddo

# here the sds are ion and electron poloidal power fluxes in W/m**2
      do jx=1,nxpt
        ixi=ixlb(jx)	# ixi=0
        ixo=ixrb(jx)	# ixo=nx
        do iy=1,ny+1
          sdirb(iy,jx) = 0.
          sdilb(iy,jx) = 0.
          do id = 1, nfsp
           if (zi(id) .gt. 0) then
	    sdirbd(iy,id,jx) = ( 0.5*mi(id)*upi(ixo,iy,id)**2*
     .                           fnix(ixo,iy,id) )/sxo(iy,jx)
	    sdilbd(iy,id,jx) = (-0.5*mi(id)*upi(ixi,iy,id)**2*
     .                           fnix(ixi,iy,id) )/sxi(iy,jx)
           else    # note: upi=0 for neutral species; use up instead
	    sdirbd(iy,id,jx) = ( 0.5*mi(id)*up(ixo,iy,id)**2*
     .                           fnix(ixo,iy,id) )/sxo(iy,jx)
	    sdilbd(iy,id,jx) = (-0.5*mi(id)*up(ixi,iy,id)**2*
     .                           fnix(ixi,iy,id) )/sxi(iy,jx)
           endif
           sdirb(iy,jx) = sdirb(iy,jx) + sdirbd(iy,id,jx)
           sdilb(iy,jx) = sdilb(iy,jx) + sdilbd(iy,id,jx)
         enddo
         do id = 1, nusp      # note: up neutral species in nonzero
           tempvo =  - ckinfl*0.5*sx(ixo,iy)*visx(ixo,iy,id)*gx(ixo,iy)*
     .               ( up(ixo,iy,id)**2 - up(ixo-1,iy,id)**2 ) /sxo(iy,jx)
	   tempvi =  + ckinfl*0.5*sx(ixi,iy)*visx(ixi+1,iy,id)*gx(ixi+1,iy)*
     .                ( up(ixi+1,iy,id)**2 - up(ixi,iy,id)**2 ) / sxi(iy,jx)
	   sdirbd(iy,id,jx) = sdirbd(iy,id,jx) + tempvo
	   sdirb(iy,jx) = sdirb(iy,jx) + tempvo
	   sdilbd(iy,id,jx) = sdilbd(iy,id,jx) + tempvi
	   sdilb(iy,jx) = sdilb(iy,jx) + tempvi
         enddo
         sdirb(iy,jx) = sdirb(iy,jx) + feix(ixo,iy)/sxo(iy,jx)
         sbindrb(iy,jx) =  fnix(ixo,iy,1) * ebind*ev / sxo(iy,jx)
         sderb(iy,jx) =  ( feex(ixo,iy)+fqx(ixo,iy)*
     .                     (phi(ixo+1,iy)-phi0r(iy,jx)) )/sxo(iy,jx)
         sdtrb(iy,jx) = sderb(iy,jx) + sdirb(iy,jx) + sbindrb(iy,jx)
         sdilb(iy,jx) = sdilb(iy,jx) - feix(ixi,iy)/sxi(iy,jx) 
         sbindlb(iy,jx) = -fnix(ixi,iy,1) * ebind*ev / sxi(iy,jx)
         sdelb(iy,jx) = -( feex(ixi,iy)+fqx(ixi,iy)*
     .                     (phi(ixi  ,iy)-phi0l(iy,jx)) )/sxi(iy,jx)
         sdtlb(iy,jx) = sdelb(iy,jx) + sdilb(iy,jx) + sbindlb(iy,jx)
         pdivirb(jx) = pdivirb(jx) + sdirb(iy,jx)*sxo(iy,jx)
         pdiverb(jx) = pdiverb(jx) + sderb(iy,jx)*sxo(iy,jx)
         pdivilb(jx) = pdivilb(jx) + sdilb(iy,jx)*sxi(iy,jx)
         pdivelb(jx) = pdivelb(jx) + sdelb(iy,jx)*sxi(iy,jx)
         pbindrb(jx) = pbindrb(jx) + sbindrb(iy,jx)*sxo(iy,jx)
         pbindlb(jx) = pbindlb(jx) + sbindlb(iy,jx)*sxi(iy,jx)
         if (isupgon(1).eq.1) then    # Approx. neutral energy flux
	   sdnrb(iy,jx)=sdirbd(iy,jx,2) + ( sx(ixo,iy)*hcxn(ixo,iy)*
     .                    gxf(ixo,iy)*(ti(ixo,iy)-ti(ixo+1,iy)) + 
     .   		  2.5*fnix(ixo,iy,2)*ti(ixo+1,iy) ) / sxo(iy,jx)
	   sdnlb(iy,jx)=sdilbd(iy,jx,2) - ( sx(ixi,iy)*hcxn(ixi,iy)*
     .                    gxf(ixi,iy)*(ti(ixi,iy)-ti(ixi+1,iy)) + 
     .                    2.5*fnix(ixi,iy,2)*ti(ixi  ,iy) ) / sxi(iy,jx)
	   pdivnrb(jx) = pdivnrb(jx) + sdnrb(iy,jx)*sxo(iy,jx)
	   pdivnlb(jx) = pdivnlb(jx) + sdnlb(iy,jx)*sxi(iy,jx)
        endif
        enddo  # end do-loop for iy=1,ny+1
      enddo  # end do-loop for jx=1,nxpt

# Fix corner boundary values
      do jx = 1, nxpt
         sdtlb(0,jx) = sdtlb(1,jx)
         sdtlb(ny+1,jx) = sdtlb(ny,jx)
         sdelb(0,jx) = sdelb(1,jx)
         sdelb(ny+1,jx) = sdelb(ny,jx)
         sdilb(0,jx) = sdilb(1,jx)
         sdilb(ny+1,jx) = sdilb(ny,jx)
         sdrlb(ny+1,jx) = sdrlb(ny,jx)
         sdtrb(0,jx) = sdtrb(1,jx)
         sdtrb(ny+1,jx) = sdtrb(ny,jx)
         sderb(0,jx) = sderb(1,jx)
         sderb(ny+1,jx) = sderb(ny,jx)
         sdirb(0,jx) = sdirb(1,jx)
         sdirb(ny+1,jx) = sdirb(ny,jx)
         sdrrb(ny+1,jx) = sdrrb(ny,jx)
      enddo

      return
      end
******* end of subroutine plateflux *******
c-------------------------------------------------------------------------c

      subroutine wallflux

      Implicit none

c ... Calc particle & heat fluxes to outer wall surfaces; only PF rad flux?
c ... Use as diagnostic called from BASIS/Python parser (alt. to balancee)


      Use(Dim)            # nx,ny,nxpt
      Use(Comgeo)         # sy
      Use(Phyvar)         # pi, ev
      Use(Postproc)       # pwr_wallh, Pwr_wallz
                          # swallr, swalli, swalle, swbind, swallt
                          # spfwallr
      Use(Compla)         # mi,ni,up,te,ti
      Use(Comflo)         # fniy,feey,feiy
      Use(UEpar)          # ebind

c  Local variables
      integer jx,id,ixi,ixo,ip
      #Former Aux module variables
      integer ix,iy
      real pwallr,pwalli,pwalle,pwbind
      real swallid(0:nx+1,nfsp)

########################################
# First do the wall particle flux
##############################################################
      do ix = 1, nx
        do id = 1, nfsp
	  gwalli(ix,id) = fniy(ix,ny,id)/sy(ix,ny)
          engwalli(ix,id) = 2.*ti(ix,ny+1)/ev + zi(id)*phi(ix,ny+1)
        enddo
      enddo

# Fix corner boundary values
      do id = 1, nfsp
        gwalli(0,id) = gwalli(1,id)
        gwalli(nx+1,id) = gwalli(nx,id)
        engwalli(0,id) = engwalli(1,id)
        engwalli(nx+1,id) = engwalli(nx,id)
      enddo

#######################################
# Now do the wall heat flux
##############################################################
#
# Initialize total powers (local use only)
      pwallr = 0.
      pwalli = 0.
      pwalle = 0.
      pwbind = 0.

# Get radiation power in pwr_plth and pwr_plt
      call pradpltwl

# Fill radiation heat flux arrays
      do ix = 0, nx+1
	swallr(ix) = pwr_wallh(ix) + pwr_wallz(ix)
        do ip = 1, nxpt
          spfwallr(ix,ip) = pwr_pfwallz(ix,ip)+pwr_pfwallh(ix,ip)
        enddo
      enddo

# swalls are ion and electron radial power fluxes in W/m**2
      do ix = 1, nx
        swalli(ix) = 0.
        do id = 1, nfsp
          if (zi(id) .gt. 0) then
            swallid(ix,id) = ( 0.5*mi(id)*upi(ix,ny,id)**2*
     .                        fniy(ix,ny,id) )/sy(ix,ny)
          else    # note: upi=0 for neutral species; use up instead
	    swallid(ix,id) = ( 0.5*mi(id)*up(ix,ny,id)**2*
     .                         fniy(ix,ny,id) )/sy(ix,ny)
          endif
          swalli(ix) = swalli(ix) + swallid(ix,id)
        enddo
        swalli(ix) = swalli(ix) + feiy(ix,ny)/sy(ix,ny)
        swbind(ix) = fniy(ix,ny,1) * ebind*ev / sy(ix,ny)
        swalle(ix) = (feey(ix,ny)+fqy(ix,ny)*phi(ix,ny+1) )/sy(ix,ny)
        swallt(ix) = swalle(ix) + swalli(ix) + swbind(ix) + swallr(ix)
        pwallr = pwallr + swallr(ix)*sy(ix,ny)
        pwalle = pwalle + swalle(ix)*sy(ix,ny)
        pwalli = pwalli + swalli(ix)*sy(ix,ny)
        pwbind = pwbind + swbind(ix)*sy(ix,ny)
      enddo  # end do-loop for ix=1,nx

# Fix corner boundary values
      swallr(0) = swallr(1)
      swalli(0) = swalli(1)
      swalle(0) = swalle(1)
      swbind(0) = swbind(1)
      swallt(0) = swallt(1)
      swallr(nx+1) = swallr(nx)
      swalli(nx+1) = swalli(nx)
      swalle(nx+1) = swalle(nx)
      swbind(nx+1) = swbind(nx)
      swallt(nx+1) = swallt(nx)

      return
      end
******* end of subroutine wallflux *******
c-------------------------------------------------------------------------c

      subroutine bbb2wdf
      implicit none
Use(Dim)		# nx,ny
Use(Share)		# nycore,nysol,nxcore,nxleg,igrid,geometry
Use(RZ_grid_info)	# br,bz,bpol,bphi,b
Use(Bfield)             # rbfbt
Use(Bcond)		# fngysi,fngyso
Use(Parallv)            # nxg,nyg
Use(Selec)		# ixm1
Use(Comflo)		# fnix
Use(Compla)		# ni,uu,up,v2,vy,te,ti,ne
Use(Linkbbb)		# 

c     Local variables --
      integer ix,iy,nunit
      real uuc,upc,vyc,v2c

c     Compute data for output to wdf package

      nxbbb=nx
      nybbb=ny
      nycorebbb=nycore(igrid)
      nysolbbb=nysol(igrid)
      nxleg1bbb=nxleg(igrid,1)
      nxcore1bbb=nxcore(igrid,1)
      nxleg2bbb=nxleg(igrid,2)
      nxcore2bbb=nxcore(igrid,2)
      geometrybbb=geometry

      do ix=0,nx+1
         do iy=0,ny+1
            nibbb(ix,iy)=ni(ix,iy,1)
            tibbb(ix,iy)=ti(ix,iy)
            nebbb(ix,iy)=ne(ix,iy)
            tebbb(ix,iy)=te(ix,iy)
            fnixbbb(ix,iy)=fnix(ix,iy,1)
         enddo
         fngysibbb(ix)=fngysi(ix,1)
         fngysobbb(ix)=fngyso(ix,1)
      enddo

c     Compute ion velocity at cell centers :
      do ix=1,nx
	 do iy=1,ny
	    uuc=0.5*(uu(ix,iy,1)+uu(ixm1(ix,iy),iy,1))
	    vyc=0.5*(vy(ix,iy,1)+vy(ix,iy-1,1))
	    vflowxbbb(ix,iy) =   uuc*br(ix,iy,0)/bpol(ix,iy,0)
     &                       - vyc*bz(ix,iy,0)/bpol(ix,iy,0)
	    vflowzbbb(ix,iy) =   uuc*bz(ix,iy,0)/bpol(ix,iy,0)
     &                       + vyc*br(ix,iy,0)/bpol(ix,iy,0)

	    v2c=0.5*(v2(ix,iy,1)+v2(ixm1(ix,iy),iy,1))
	    upc=0.5*(up(ix,iy,1)+up(ixm1(ix,iy),iy,1))
	    vflowybbb(ix,iy) = - v2c*bpol(ix,iy,0)/b(ix,iy,0)
     &                       + upc*rbfbt(ix,iy)
	 enddo
      enddo

c ********* Write the file to be read by the wdf package *************

      call freeus (nunit)
      open (nunit, file='bbb-wdf', form='unformatted', status='unknown')
      write(nunit)nxbbb,nybbb,nycorebbb,nysolbbb,nxleg1bbb,nxcore1bbb,
     &            nxleg2bbb,nxcore2bbb
      write(nunit)nibbb,tibbb,nebbb,tebbb,vflowxbbb,vflowybbb,vflowzbbb,
     &            fnixbbb,fngysibbb,fngysobbb
      write(nunit)geometrybbb
      close (nunit)

      return
      end

c----------------------------------------------------------------------c

      subroutine scale_mcn
      implicit none
Use(Dim)	# nx,ny,nisp,nusp
Use(UEpar)      # isupgon,iigsp
Use(Comflo)	# fnix
Use(MCN_dim)
Use(MCN_sources)

      integer ix,iy,ifld,istra

      external remark

c     This subroutine scales plasma source terms obtained from the 
c     Monte-Carlo-Neutrals code.  These sources are assumed to scale
c     with the neutral source currents at the divertor plates.
c     The flag that controls the scaling is ismcnvar :
c         ismcnvar=0  -->  sources are fixed
c         ismcnvar=1  -->  sources scale with current
c
c     In the terminology of the EIRENE code, a stratum is a surface
c     or volume element where neutral particles originate.  Here
c     we consider only two possible strata: the inboard and outboard
c     divertor plates.  We neglect strata associated with gas puffing
c     or recombination.
c     
c     Plasma source terms from EIRENE are :
c         mcnsor_ni
c         mcnsor_up
c         mcnsor_te
c         mcnsor_ti
c     where mcnsor_ni(ix,iy,ifld,istra) is the particle source
c     for ion fluid 'ifld' at cell (ix,iy) due to neutral source
c     stratum 'istra'.
c
c     mcncurr(istra) is the incident ion current that
c     normalizes the plasma sources due to neutrals
c     that originate at stratum 'istra'.
c
c     The scaled plasma source terms for UEDGE are :
c         uesor_ni * cmneutsor_ni
c         uesor_up * cmneutsor_mi
c         uesor_te * cmneutsor_ee
c         uesor_ti * cmneutsor_ei
c     where uesor_ni is obtained by summing the scaled mcnsor_ni
c     over all strata, and similarly for up,te,ti.

ccc=================================================================
ccc     NOTE: For the cmod-box problem there is only ONE strata !!!
ccc     (in general, we need to devise a test to identify strata)
ccc=================================================================

c     Compute scale factors for each strata :
      do istra=1,nstra
         strascal(istra)=1.	# default is no scaling
      enddo
      if (ismcnvar==1) then	# sources scale with current
         do istra=1,nstra
            uecurr(istra)=0.
            if (istra==1) then	# for east target plate only
               do iy=1,ny
                  uecurr(istra)=uecurr(istra)+fnix(nx,iy,1)
               enddo
            else
               call remark("***")
               call remark("***    subroutine scale_mcn    ***")
               call remark("***  not defined for nstra > 1  ***")
               call remark("***")
            endif
            if(mcncurr(istra) > 0) then
               strascal(istra)=uecurr(istra)/mcncurr(istra)
            endif
         enddo
      endif

c     Scale source terms and sum over strata :
      do iy=0,ny+1
         do ix=0,nx+1
            do ifld=1,nisp
               uesor_ni(ix,iy,ifld)=0.
            enddo
            do ifld=1,nusp
               uesor_up(ix,iy,ifld)=0.
            enddo
            uesor_te(ix,iy)=0.
            uesor_ti(ix,iy)=0.
         enddo
      enddo
      do istra=1,nstra
         do iy=0,ny+1
            do ix=0,nx+1
               do ifld=1,nisp
                  uesor_ni(ix,iy,ifld)=uesor_ni(ix,iy,ifld)
     &                    +mcnsor_ni(ix,iy,ifld,istra)*strascal(istra)
               enddo
               do ifld=1,nusp
                  uesor_up(ix,iy,ifld)=uesor_up(ix,iy,ifld)
     &                    +mcnsor_up(ix,iy,ifld,istra)*strascal(istra)
               enddo
               uesor_te(ix,iy)=uesor_te(ix,iy)
     &                    +mcnsor_te(ix,iy,istra)*strascal(istra)
               uesor_ti(ix,iy)=uesor_ti(ix,iy)
     &                    +mcnsor_ti(ix,iy,istra)*strascal(istra)
            enddo
         enddo
      enddo

      return
      end

#----------------------------------------------------------------------#

      subroutine write30 (fname, runid)
      implicit none
      character*(*) fname, runid

Use(Dim)            # nx,ny
Use(UEint)          # mhdgeo
Use(Share)          # geometry,nxomit
Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx
Use(RZ_grid_info)   # rm,zm

      integer nunit, ix, iy, i, iv(4)
      integer nncut,nniso
      integer nxcut1(2),nxcut2(2),nycut1(2),nycut2(2)
      integer nxiso1(2),nxiso2(2),nyiso1(2),nyiso2(2)

      external remark

      data nncut/0/,nniso/0/
      data nxcut1/2*0/,nxcut2/2*0/,nycut1/2*0/,nycut2/2*0/
      data nxiso1/2*0/,nxiso2/2*0/,nyiso1/2*0/,nyiso2/2*0/

# This subroutine writes the geometry data file that is read by the
# EIRENE code, using the B2 code variable definitions.

cccMER NOTE: generalize the following for 2 or more x-points
# Convert UEDGE variables to B2/EIRENE data -
      if (mhdgeo .eq. 1) then
         nncut=2
         nxcut1(1)=ixpt1(1)
         nxcut2(1)=ixpt2(1)+1
         nycut1(1)=0
         nycut2(1)=iysptrx1(1)
         nxcut1(2)=ixpt2(1)
         nxcut2(2)=ixpt1(1)+1
         nycut1(2)=0
         nycut2(2)=iysptrx1(1)
      endif
# UEDGE cell vertex iv(i) corresponds to B2 cell vertex i :
      iv(1)=2
      iv(2)=4
      iv(3)=3
      iv(4)=1

# Write the data -
      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')
      write (nunit,*) runid
      write (nunit,*) ' '
      write (nunit,*) nx,ny,nncut
      write (nunit,*) 
     &              (nxcut1(i),nxcut2(i),nycut1(i),nycut2(i),i=1,nncut)
      if (nncut .gt. 2) then
         write (nunit,*) nniso
         write (nunit,*) 
     &              (nxiso1(i),nxiso2(i),nyiso1(i),nyiso2(i),i=1,nniso)
      endif
      write (nunit,*) ' '
      do ix=1,nx
         do iy=1,ny
            if (mhdgeo .eq. 1) then
               write (nunit,'(4e15.7)') (rm(ix,iy,iv(i)),i=1,4)
               write (nunit,'(4e15.7)') (zm(ix,iy,iv(i)),i=1,4)
            else
               write (nunit,'(4e15.7)') (zm(ix,iy,iv(i)),i=1,4)
               write (nunit,'(4e15.7)') (rm(ix,iy,iv(i)),i=1,4)
            endif
         enddo
      enddo

      close (nunit)
      call remark(" *** geometry file written for EIRENE ***")

      return
      end

#----------------------------------------------------------------------#


      subroutine write31 (fname, runid)
      implicit none
      character*(*) fname, runid

Use(Dim)            # nx,ny,nxm,nisp
Use(UEpar)          # isupgon
Use(Compla)         # ni,uu,up,vy,te,ti,pr,zi
Use(Comflo)         # fnix,fniy,feix,feiy,feex,feey
Use(RZ_grid_info)   # b
Use(Comgeo)         # vol,rr

      integer nunit,nxd,nyd,ifld
      external remark

# This subroutine writes the plasma data file that is read by the
# EIRENE code, using the B2 code variable definitions.  This file only
# includes data for charged species; inertial neutrals are excluded.

      nxd = nx
      nyd = ny

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0) 
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,ni(0:nxd+1,0:nyd+1,ifld))
      enddo
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0) 
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,uu(0:nxd+1,0:nyd+1,ifld))
      enddo
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0) 
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,vy(0:nxd+1,0:nyd+1,ifld))
      enddo
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,te(0:nxd+1,0:nyd+1))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,ti(0:nxd+1,0:nyd+1))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,pr(0:nxd+1,0:nyd+1))
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0) 
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,upi(0:nxd+1,0:nyd+1,ifld))
      enddo
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,rr(0:nxd+1,0:nyd+1))
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0) 
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,fnix(0:nxd+1,0:nyd+1,ifld))
      enddo
      do ifld = 1, nisp
         if (zi(ifld) .gt. 0) 
     .      call gfsub3 (nunit,nx,ny,nxd,nyd,1,fniy(0:nxd+1,0:nyd+1,ifld))
      enddo
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feix(0:nxd+1,0:nyd+1))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feiy(0:nxd+1,0:nyd+1))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feex(0:nxd+1,0:nyd+1))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,feey(0:nxd+1,0:nyd+1))
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,vol(0:nxd+1,0:nyd+1))

      nxd = nxm		# note dimension nxm rather than nx for b(,,)
      call gfsub3 (nunit,nx,ny,nxd,nyd,1,b(0:nxd+1,0:nyd+1,0))

      close (nunit)
      call remark(" *** background plasma file written for EIRENE ***")

      return
      end

#----------------------------------------------------------------------#

      subroutine gfsub3 (nunit,nx,ny,ndimx,ndimy,ndims,dummy)
      implicit none
      integer nunit,nx,ny,ndimx,ndimy,ndims
      real dummy(0:ndimx+1,0:ndimy+1,1:ndims)

      integer nd1,lim,is,iy,ix,iii
      real eps
      data eps/1.0e-90/

# This subroutine writes an array 'dummy' into an output file.
# (Copied from the B2 code by MER 95/12/12; modified 97/03/11)
# The additive constant, eps, ensures that extremely small non-zero
# values of array elements will not corrupt the e-format in the output file.

      nd1 = nx + 2
      lim = (nd1/5)*5 - 4
      do is=1,ndims
         do iy=0,ny+1
            do ix=1,lim,5
               write(nunit,910) (dummy(ix-2+iii,iy,is)+eps,iii=1,5)
            enddo
            if ( (lim+4) .lt. nd1 ) # write partial line
     .         write(nunit,910) (dummy(ix-1,iy,is)+eps,ix=lim+5,nd1)
         enddo
      enddo

 910  format(5(e16.8))

      return
      end

#----------------------------------------------------------------------#

      subroutine write_eirene
      implicit none

# This subroutine writes the geometry (fort.30) and plasma (fort.31)
# files that supply information to the EIRENE Monte Carlo neutrals code.


      call write30("fort.30","UEDGE geometry data")
      call write31("fort.31","UEDGE plasma data")

      return
      end

#----------------------------------------------------------------------#

      subroutine read32(fname)
      implicit none
      character*(*) fname

Use(Dim)		# nx,ny
Use(MCN_dim)
Use(MCN_sources)
      integer nunit,ix,iy,ifl,istra

c     Read data from EIRENE code output data file 'fort.32':

c     NOTE: Dimensions nxf and nyf should be read from
c            EIRENE code output data file 'fort.44'
c     Here, we assume nxf=nx and nyf=ny elements in ix and iy.
c     Also, nstra and nfl need to be set properly.

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

      do istra=1,nstra

c     Normalization constants:
         read (nunit,*) wsor(istra), esor(istra)

c     Normalized plasma source terms:
         do ifl=1,nfl
            read (nunit,*) ((sni(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smo(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
         enddo
         read (nunit,*) ((see(ix,iy,istra), ix=1,nx), iy=1,ny)
         read (nunit,*) ((sei(ix,iy,istra), ix=1,nx), iy=1,ny)

      enddo

      close (nunit)
      call remark(" *** plasma sources read from file fort.32 ***")

      return
      end

c----------------------------------------------------------------------c

      subroutine writemcnfile(fname, runidtag)

      implicit none
      character*(*) fname, runidtag

c     This subroutine writes the UEDGE mesh and plasma data file that
c     is read by the DEGAS2 code.  Plasma velocities are converted to
c     Cartesian or cylindrical components at density cell centers.

Use(Dim)             # nx,ny,nxpt,nisp
Use(Xpoint_indices)  # ixlb,ixpt1,ixmdp,ixpt2,ixrb,iysptrx1,iysptrx2
Use(RZ_grid_info)    # rm,zm,br,bz,bphi
Use(Ext_neutrals) 		 # ext_verbose
      integer nunit,ix,iy,n,jx

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

c write geometry data
      write (nunit,*) runidtag
      write (nunit,*) nx,ny,nxpt
      do jx=1,nxpt
         write (nunit,*) iysptrx1(jx),iysptrx2(jx)
         write (nunit,*) ixlb(jx),ixpt1(jx),ixmdp(jx),ixpt2(jx),ixrb(jx)
      enddo
      write (nunit,*) (((rm(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((zm(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((br(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((bz(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)
      write (nunit,*) (((bphi(ix,iy,n),ix=1,nx),iy=1,ny),n=0,4)

c write plasma data
      call gchange("MCN_bkgd",0)
      call writemcnbkgd(nunit)

      close (nunit)
      if (ext_verbose)
     . call remark(" *** data file "//fname//" written for DEGAS2 ***")

      return
      end

c----------------------------------------------------------------------c

      subroutine writemcnbkgd(nunit)
      implicit none
      integer nunit

Use(Dim)             # nx,ny,nxpt,nisp
Use(Xpoint_indices)  # ixlb,ixrb
Use(RZ_grid_info)    # rm,zm,br,bz,bphi
Use(Bfield)          # rbfbt
Use(Comgeo)          # rr
Use(Compla)          # zi,ni,upi,uu,vy,v2,te,ti
Use(Comflo)          # fnix
Use(MCN_bkgd)        # various velocities
Use(UEint)           # mhdgeo
Use(Selec)           # ixm1,ixp1

c local variables --
      integer ix,iy,ifld,jx,ixt,ixt1
      real sinb1,cosb1


c compute UEDGE flow velocity at cell centers:
      do ifld=1,nisp
         do ix=1,nx
            do iy=1,ny
               # 3 orthogonal components of the velocity are:
               v2c(ix,iy,ifld)=(v2(ixm1(ix,iy),iy,ifld)+v2(ix,iy,ifld))/2.
               vyc(ix,iy,ifld)=(vy(ix,iy-1,ifld)+vy(ix,iy,ifld))/2.
               upc(ix,iy,ifld)=(upi(ixm1(ix,iy),iy,ifld)+upi(ix,iy,ifld))/2.
               # from these we construct the poloidal and toroidal components:
               uuc(ix,iy,ifld)=upc(ix,iy,ifld)*rr(ix,iy)
     .                         +v2c(ix,iy,ifld)*rbfbt(ix,iy)
               utc(ix,iy,ifld)=upc(ix,iy,ifld)*rbfbt(ix,iy)
     .                         -v2c(ix,iy,ifld)*rr(ix,iy)
            enddo
         enddo
      enddo
      if (mhdgeo .eq. 1) then		# convert to cylindrical coordinates:
         do ifld=1,nisp
            do ix=1,nx
               do iy=1,ny
	          cosb1=br(ix,iy,0)/bpol(ix,iy,0)
	          sinb1=bz(ix,iy,0)/bpol(ix,iy,0)
	          vr(ix,iy,ifld)=uuc(ix,iy,ifld)*cosb1-vyc(ix,iy,ifld)*sinb1
	          vz(ix,iy,ifld)=uuc(ix,iy,ifld)*sinb1+vyc(ix,iy,ifld)*cosb1
	          vphi(ix,iy,ifld)=utc(ix,iy,ifld)
               enddo
            enddo
         enddo
      elseif (mhdgeo .eq. -1) then	# cartesian slab coordinates:
         do ifld=1,nisp
            do ix=1,nx
               do iy=1,ny
                  vr(ix,iy,ifld)=vyc(ix,iy,ifld)
	          vz(ix,iy,ifld)=uuc(ix,iy,ifld)
	          vphi(ix,iy,ifld)=utc(ix,iy,ifld)
               enddo
            enddo
         enddo
      endif

c write plasma data (charged species only)
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((ni(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((vr(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((vz(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) ((vphi(ix,iy,ifld),ix=1,nx),iy=1,ny)
         endif
      enddo
      write (nunit,*) ((ti(ix,iy),ix=1,nx),iy=1,ny)
      write (nunit,*) ((te(ix,iy),ix=1,nx),iy=1,ny)

c loop over nxpt mesh regions
      do jx=1,nxpt
c compute UEDGE flow velocity at left target plate:
      do ifld=1,nisp
         do iy=1,ny
            ixt=ixlb(jx)       # analog of ix=0
            ixt1=ixp1(ixt,iy)  # analog of ix=1
            # 3 orthogonal components of the velocity are:
            v2tg1(iy,ifld)=v2(ixt,iy,ifld)
            vytg1(iy,ifld)=(vy(ixt1,iy-1,ifld)+vy(ixt1,iy,ifld))/2. # <- NOTE
            uptg1(iy,ifld)=upi(ixt,iy,ifld)
            # from these we construct the poloidal and toroidal components:
            uutg1(iy,ifld)=uptg1(iy,ifld)*rr(ixt,iy)
     .                      +v2tg1(iy,ifld)*rbfbt(ixt,iy)
            uttg1(iy,ifld)=uptg1(iy,ifld)*rbfbt(ixt,iy)
     .                      -v2tg1(iy,ifld)*rr(ixt,iy)
         enddo
      enddo
      if (mhdgeo .eq. 1) then		# convert to cylindrical coordinates:
         do ifld=1,nisp
            do iy=1,ny
               cosb1=br(ixt,iy,0)/bpol(ixt,iy,0)
               sinb1=bz(ixt,iy,0)/bpol(ixt,iy,0)
               vrtg1(iy,ifld)=uutg1(iy,ifld)*cosb1-vytg1(iy,ifld)*sinb1
               vztg1(iy,ifld)=uutg1(iy,ifld)*sinb1+vytg1(iy,ifld)*cosb1
               vphitg1(iy,ifld)=uttg1(iy,ifld)
            enddo
         enddo
      elseif (mhdgeo .eq. -1) then	# cartesian slab coordinates:
         do ifld=1,nisp
            do iy=1,ny
               vrtg1(iy,ifld)=vytg1(iy,ifld)
               vztg1(iy,ifld)=uutg1(iy,ifld)
               vphitg1(iy,ifld)=uttg1(iy,ifld)
            enddo
         enddo
      endif

c write inner target data
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (ni(ixt,iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vrtg1(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vztg1(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vphitg1(iy,ifld),iy=1,ny)
         endif
      enddo
      write (nunit,*) (ti(ixt,iy),iy=1,ny)
      write (nunit,*) (te(ixt,iy),iy=1,ny)
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (-fnix(ixt,iy,ifld),iy=1,ny)
         endif
      enddo

c compute UEDGE flow velocity at right target plate:
      do ifld=1,nisp
         do iy=1,ny
            ixt=ixrb(jx)+1     # analog of nx+1
            ixt1=ixm1(ixt,iy)  # analog of nx
            # 3 orthogonal components of the velocity are:
            v2tg2(iy,ifld)=v2(ixt1,iy,ifld)
            vytg2(iy,ifld)=(vy(ixt1,iy-1,ifld)+vy(ixt1,iy,ifld))/2. # <- NOTE
            uptg2(iy,ifld)=upi(ixt1,iy,ifld)
            # from these we construct the poloidal and toroidal components:
            uutg2(iy,ifld)=uptg2(iy,ifld)*rr(ixt,iy)
     .                      +v2tg2(iy,ifld)*rbfbt(ixt,iy)
            uttg2(iy,ifld)=uptg2(iy,ifld)*rbfbt(ixt,iy)
     .                      -v2tg2(iy,ifld)*rr(ixt,iy)
         enddo
      enddo
      if (mhdgeo .eq. 1) then		# convert to cylindrical coordinates:
         do ifld=1,nisp
            do iy=1,ny
               cosb1=br(ixt,iy,0)/bpol(ixt,iy,0)
               sinb1=bz(ixt,iy,0)/bpol(ixt,iy,0)
               vrtg2(iy,ifld)=uutg2(iy,ifld)*cosb1-vytg2(iy,ifld)*sinb1
               vztg2(iy,ifld)=uutg2(iy,ifld)*sinb1+vytg2(iy,ifld)*cosb1
               vphitg2(iy,ifld)=uttg2(iy,ifld)
            enddo
         enddo
      elseif (mhdgeo .eq. -1) then	# cartesian slab coordinates:
         do ifld=1,nisp
            do iy=1,ny
               vrtg2(iy,ifld)=vytg2(iy,ifld)
               vztg2(iy,ifld)=uutg2(iy,ifld)
               vphitg2(iy,ifld)=uttg2(iy,ifld)
            enddo
         enddo
      endif

c write outer target data
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (ni(ixt,iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vrtg2(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vztg2(iy,ifld),iy=1,ny)
         endif
      enddo
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (vphitg2(iy,ifld),iy=1,ny)
         endif
      enddo
      write (nunit,*) (ti(ixt,iy),iy=1,ny)
      write (nunit,*) (te(ixt,iy),iy=1,ny)
      do ifld=1,nisp
         if (zi(ifld) .ne. 0) then
            write (nunit,*) (fnix(ixt1,iy,ifld),iy=1,ny)
         endif
      enddo

      enddo   # end do-loop over nxpt mesh regions

      return
      end

c----------------------------------------------------------------------c

      subroutine readmcnsor(fname)
      implicit none
      character*(*) fname

Use(Dim)		# nx,ny
Use(MCN_dim)
Use(MCN_sources)
Use(Ext_neutrals) 		 # ext_verbose
      integer nunit,ix,iy,ifl,istra

c     Read data from DEGAS2 code output data file 'sources.out':

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

      do istra=1,nstra

c     Integrated ion particle source:
         read (nunit,*) wsor(istra)

c     Plasma source terms:
         do ifl=1,nfl
            read (nunit,*) ((sni(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smor(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smophi(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
            read (nunit,*) ((smoz(ix,iy,ifl,istra), ix=1,nx), iy=1,ny)
         enddo
         read (nunit,*) ((see(ix,iy,istra), ix=1,nx), iy=1,ny)
         read (nunit,*) ((sei(ix,iy,istra), ix=1,nx), iy=1,ny)

      enddo

      close (nunit)
      if (ext_verbose)
     . call remark(" *** plasma sources read from DEGAS2 file "
     .                                                  //fname//" ***")

      return
      end

#----------------------------------------------------------------------#

      subroutine readmcntest(fname)
      implicit none
      character*(*) fname

Use(Dim)
Use(MCN_dim)

      integer nunit
c     Read data from DEGAS2 code output file 'testdata.out':

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

c     Dimensions:
      read (nunit,*) nxf, nyf, nmcsp

      if (nmcsp .gt. nmcmx) then
         call remark('***')
         call remark('*** READMCNTEST: nmcsp > nmcmx')
         call remark('                 re-compile with larger nmcmx')
         call remark('***')
         call xerrab("")
      endif

      call gchange("MCN_test",0)

      call readmcntesta(nunit)


      close(nunit)
      call remark(" *** neutral diagnostics read from DEGAS2 file "
     .                                                  //fname//" ***")

      return
      end

#----------------------------------------------------------------------#

      subroutine readmcntesta(nunit)
      implicit none
      integer nunit

Use(Dim)
Use(MCN_dim)
Use(MCN_test)

      integer ix,iy,id

c Species names:
      read (nunit,*) (labelmc(id), id=1,nmcsp)

c Densities and 'temperatures':
      do id=1,nmcsp
         read (nunit,*) ((nmc(ix,iy,id),ix=1,nxf),iy=1,nyf)
         read (nunit,*) ((tmc(ix,iy,id),ix=1,nxf),iy=1,nyf)
      enddo

c Particle fluxes:
      do id=1,nmcsp
         read (nunit,*) ((fnmcx(ix,iy,id),ix=1,nxf),iy=1,nyf)
         read (nunit,*) ((fnmcy(ix,iy,id),ix=1,nxf),iy=1,nyf)
         read (nunit,*) ((fnmcz(ix,iy,id),ix=1,nxf),iy=1,nyf)
      enddo

      return
      end

#----------------------------------------------------------------------#

      subroutine read44(fname)
      implicit none
      character*(*) fname

Use(Dim)
Use(MCN_dim)
Use(MCN_sources)

      integer nunit

c     Read data from EIRENE code output data file 'fort.44':

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

c     Dimensions:
      read (nunit,*) nxf, nyf
      read (nunit,*) natmi, nmoli, nioni

c     NOTE: 
c     Character arrays (labela,labelm,labeli) in MCN_sources have fixed
c     dimension; Dynamic allocation not allowed on C90 (MER 97/03/11).
      if ( (natmi .gt. nmcmx) .or. (nmoli .gt. nmcmx) 
     .                         .or. (nioni .gt. nmcmx) ) then
         call remark('***')
         call remark('*** READ44: natmi or nmoli or nioni > nmcmx')
         call remark('            re-compile with larger nmcmx')
         call remark('***')
         call xerrab("")
      endif

      call gchange("MCN_sources",0)

      call read44a(nunit)

      close (nunit)
      call remark(' *** neutral diagnostics read from file fort.44 ***')

      return
      end

#----------------------------------------------------------------------#

      subroutine read44a(nunit)
      implicit none
      integer nunit

Use(Dim)
Use(MCN_dim)
Use(MCN_sources)

      integer ix,iy,id

c     Species names:
      read (nunit,*) (labela(id), id=1,natmi)
      read (nunit,*) (labelm(id), id=1,nmoli)
      read (nunit,*) (labeli(id), id=1,nioni)

c     Densities and 'temperatures' of atoms, molecules, test ions:
      read (nunit,*) (((naf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((taf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((nmf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)
      read (nunit,*) (((tmf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)
      read (nunit,*) (((ntf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nioni)
      read (nunit,*) (((ttf(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nioni)

c     Radial particle fluxes:
      read (nunit,*) (((fnay(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((fnmy(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     'parallel' particle flux:
      read (nunit,*) (((fnax(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((fnmx(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     Radial energy flux:
      read (nunit,*) (((feay(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((femy(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     'parallel' energy flux:
      read (nunit,*) (((feax(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,natmi)
      read (nunit,*) (((femx(ix,iy,id),ix=1,nxf),iy=1,nyf),id=1,nmoli)

c     H-alpha emissivity from atoms and molecules:
      read (nunit,*) ((hatm(ix,iy),ix=1,nxf),iy=1,nyf)
      read (nunit,*) ((hmol(ix,iy),ix=1,nxf),iy=1,nyf)

      return
      end

#----------------------------------------------------------------------#

      subroutine volave(n1, n2, j2, j5, i2, i5, ixp, ixm, fsprd,
     .                                             vol, ps_tmp, ps)
c ... This subroutine does a volume integral, or average, of cell source-like
c ... quantities by including interpolation based on adjacent-cell quantities

c *** Input:
c        n1 is nx for one of the two array dimensions
c        n2 is ny for the second array dimension
c        j2 and j5 are lower and upper limits of iy loop
c        i2 and i5 are lower and upper limits of ix loop
c        ixp is the 2-D index array ixp1
c        ixm is the 2-D index array ixm1
c        fsprd is the fraction spread to each of 4 adjacent cells
c        vol is the 2-D real array of cell volume
c        pstmp is an (0:n1+1,0:n2+1) array used as a local work array
c        ps is an (0:n1+1,0:n2+1) array containing the input source-like array

c *** Output:
c        ps is an (n1,n2) array containing the result of volume-averaging ps

      implicit none

c *** Input and output variables
      integer n1, n2, j2, j5, i2, i5
      integer ixm(0:n1+1,0:n2+1), ixp(0:n1+1,0:n2+1)
      real vol(0:n1+1,0:n2+1), ps_tmp(0:n1+1,0:n2+1), ps(0:n1+1,0:n2+1)
      real fsprd

c *** Local variables
      integer iy, ix, ix1, ix2, iy1, iy2
      real fs0, signps
      
      fs0 = 1. - 4*fsprd    # fraction to central cell

            do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               if (abs( ps(ix,iy)*ps(ix1,iy)*ps(ix2,iy)*
     .                  ps(ix,iy1)*ps(ix,iy2) ) .gt. 1.e-200) then
                 signps = 1.
                  ps_tmp(ix,iy) = fs0*ps(ix,iy) + fsprd*(ps(ix1,iy)+
     .                                ps(ix2,iy)+ps(ix,iy1)+ps(ix,iy2))
               endif
             enddo
           enddo
           do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               ps(ix,iy) = ps_tmp(ix,iy)
             enddo
           enddo

      return
      end
#----------------------------------------------------------------------#

      subroutine volavenv(n1, n2, j2, j5, i2, i5, ixp, ixm, fsprd, 
     .                    ps_tmp, ps)
c ... This subroutine does a volume integral, or average, of cell source-like
c ... quantities by including interpolation based on adjacent-cell quantities
c ... like volave, except here per m**3 quantities are averaged directly,
c ... so volumes are not involved (note absence of vol(ix,iy))

c *** Input:
c        n1 is nx for one of the two array dimensions
c        n2 is ny for the second  array dimension
c        j2 and j5 are lower and upper limits of iy loop
c        i2 and i5 are lower and upper limits of ix loop
c        ixp is the 2-D index array ixp1
c        ixm is the 2-D index array ixm1
c        fsprd is the fraction spread to each of 4 adjacent cells
c        pstmp is an (0:n1+1,0:n2+1) array used as a local work array
c        ps is an (0:n1+1,0:n2+1) array containing the input source-like array
c        pstmp is an (0:n1+1,0:n2+1) array used as a local work array
c        ps is an (0:n1+1,0:n2+1) array containing the input source-like array

c *** Output:
c        ps is an (n1,n2) array containing the result of volume-averaging ps

      implicit none

c *** Input and output variables
      integer n1, n2, j2, j5, i2, i5
      integer ixm(0:n1+1,0:n2+1), ixp(0:n1+1,0:n2+1)
      real ps_tmp(0:n1+1,0:n2+1), ps(0:n1+1,0:n2+1)
      real fsprd

c *** Local variables
      integer iy, ix, ix1, ix2, iy1, iy2
      real fs0, signps
       
      fs0 = 1. - 4*fsprd    # fraction to central cell

            do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               if (abs( ps(ix,iy)*ps(ix1,iy)*ps(ix2,iy)*
     .                  ps(ix,iy1)*ps(ix,iy2) ) .gt. 1.e-200) then
                 signps = 1.
                 if (ps(ix,iy) .lt. 0.) signps = -1.
                 ps_tmp(ix,iy) =   fs0*  log(abs(ps(ix ,iy)))
     .                           + fsprd*log(abs(ps(ix1,iy)))
     .                           + fsprd*log(abs(ps(ix2,iy)))
     .                           + fsprd*log(abs(ps(ix ,iy1)))
     .                           + fsprd*log(abs(ps(ix ,iy2)))
               endif
             enddo
           enddo
           do iy = j2, j5
             do ix = i2, i5
               ix1 = ixm(ix,iy)
               ix2 = ixp(ix,iy)
               iy1 = max(0,iy-1)
               iy2 = min(n2+1,iy+1)
               if (abs( ps(ix,iy)*ps(ix1,iy)*ps(ix2,iy)*
     .                  ps(ix,iy1)*ps(ix,iy2) ) .gt. 1.e-200) then
                 signps = 1.
                 if (ps(ix,iy) .lt. 0.) signps = -1.
                 ps(ix,iy) = sign( exp(ps_tmp(ix,iy)), signps)
               endif
             enddo
           enddo

      return
      end
c ======================================================================
c
      subroutine coneq

c ... This subroutine calculates the fluxes needed for the ion continuity
c ... equations

      implicit none

*  -- local variables
      integer methnx,methny,ifld
      #Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp
      Use(Share)    # nxpt,geometry,nxc,cutlo
      Use(Xpoint_indices)      # ixpt1,ixpt2,iysptrx
      Use(UEpar)    # methn,nlimix,nlimiy,nlimiy
      Use(Coefeq)   # cnfx,cnfy
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    # xlinc,xrinc,yinc,ixm1,ixp1
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # flnix,flniy
      Use(Indices_domain_dcl)   # ixmxbcl
	  
c ------------------

      do 104 ifld = 1, nfsp
*  ---------------------------------------------------------------------
*     compute flux, residual
*     The residual is: res := snic + sniv * ni - outflow(ni).
*  ---------------------------------------------------------------------

*  -- compute flnix, flox, conx --

         methnx = mod(methn, 10)
         methny = methn/10
         do 81 iy = j4, j8
            do 80 ix = i1, i5
              if (zi(ifld).eq.0. .and. ineudif <= 2) then
                 flnix(ix,iy,ifld) = fngx(ix,iy,1)
              else
               ix2 = ixp1(ix,iy)

               if (methnx .eq. 2) then   # central differencing

                  t2 = ( lni(ix, iy,ifld) + lni(ix2,iy,ifld) ) / 2

               elseif (methnx .eq. 3) then   # upwind differencing

                  if( uu(ix,iy,ifld) .ge. 0.) then
                     t2 = lni(ix,iy,ifld)
                  else
                     t2 = lni(ix2,iy,ifld)
                  endif

               else   # interp. ave or harmonic ave depending on wind*grad

                  t0 = ( lni(ix, iy,ifld)*gx(ix, iy) +
     .                   lni(ix2,iy,ifld)*gx(ix2,iy) ) / 
     .                                      ( gx(ix,iy)+gx(ix2,iy) )
                  t1 = ( gx(ix,iy)+gx(ix2,iy) ) * lni(ix,iy,ifld) *
     .                   lni(ix2,iy,ifld) / ( cutlo + lni(ix,iy,ifld)*
     .                   gx(ix2,iy) + lni(ix2,iy,ifld)*gx(ix,iy) )
                  if( uu(ix,iy,ifld)*(lni(ix,iy,ifld)-lni(ix2,iy,ifld)) 
     .                                                     .ge. 0.) then
                     t2 = t0
                  else
                     t2 = t1
                  endif

               endif

               flnix(ix,iy,ifld) = cnfx*uu(ix,iy,ifld) * sx(ix,iy) * t2
                  flnix(ix,iy,ifld) = flnix(ix,iy,ifld)/( 1 +
     .                             (nlimix(ifld)/lni(ix2,iy,ifld))**2 +
     .                             (nlimix(ifld)/lni(ix ,iy,ifld))**2 )
               endif
   80      continue
           if ((isudsym==1.or.geometry.eq.'dnXtarget') 
     .                               .and. nxc > 1) flnix(nxc,iy,ifld)=0.
           if (islimon.ne.0 .and. iy.ge.iy_lims) flnix(ix_lim,iy,ifld)=0.
           if (nxpt==2.and.ixmxbcl==1) flnix(ixrb(1)+1,iy,ifld)=0.
   81    continue


*  -- compute flniy, floy, cony --

         do 83 iy = j1, j5
            do 82 ix = i4, i8
               if (zi(ifld).eq.0.) then
                  flniy(ix,iy,ifld) = fngy(ix,iy,1)
               else
                  if (methny .eq. 2) then   # central differencing

                     t2 = ( niy0(ix,iy,ifld) + niy1(ix,iy,ifld) ) / 2

                  elseif (methny .eq. 3) then    # upwind differencing

                     if( vy(ix,iy,ifld) .ge. 0.) then
                        t2 = niy0(ix,iy,ifld)
                     else
                        t2 = niy1(ix,iy,ifld)
                     endif

                  else    # interp. ave or harmonic ave depending on wind*grad

                     t0 = ( niy0(ix,iy,ifld)*gy(ix,iy  ) +
     .                    niy1(ix,iy,ifld)*gy(ix,iy+1) ) / 
     .                    ( gy(ix,iy)+gy(ix,iy+1) )
                     t1 = ( gy(ix,iy)+gy(ix,iy+1) ) * niy0(ix,iy,ifld)*
     .                    niy1(ix,iy,ifld) / ( cutlo + niy0(ix,iy,ifld)*
     .                    gy(ix,iy+1) + niy1(ix,iy,ifld)*gy(ix,iy) )
                     if( (niy0(ix,iy,ifld)-niy1(ix,iy,ifld))*
     .                    vy(ix,iy,ifld) .ge. 0.) then
                        t2 = t0
                     else
                        t2 = t1
                     endif
                  
                  endif
               
                  flniy(ix,iy,ifld) = cnfy*vy(ix,iy,ifld)*sy(ix,iy)*t2
                  if (vy(ix,iy,ifld)*(lni(ix,iy,ifld)-lni(ix,iy+1,ifld))
     .                                                      .lt. 0.) then
                     flniy(ix,iy,ifld) = flniy(ix,iy,ifld)/( 1 +
     .                               (nlimiy(ifld)/lni(ix,iy+1,ifld))**2 +
     .                               (nlimiy(ifld)/lni(ix,iy  ,ifld))**2 )
                  endif
               endif
 82         continue
 83      continue
         
         do ix = i4, i8
            flniy(ix,ny+1,ifld) = 0.0e0
         enddo

 104  continue

      return
      end
c ***** End of subroutine coneq **********
c ======================================================================
c
      subroutine upvisneo

c ... This subroutine calculates the total up ion viscosity term with 
c ... neoclassical effects

      implicit none

*  -- local variables
      integer ifld
      real tempp,tempm,diffp,diffm
      #Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nusp
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
                               # iysptrx
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(RZ_grid_info)  # bpol,b12,b32,bsqr
      Use(Bfield)   #rbfbt
      Use(Share)    # nxpt,geometry,nxc,cutlo
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,xlinc,xrinc,ixs1,
                    # j1,j1p,j2,j2p,j3,j4,j5,j6,j5p,j6p,j7,j8,
      Use(Compla)   # ni,up,te,ti,ng,phi,v2cd,
      Use(Comflo)   # qipar
      Use(Conduc)   # visxneo,nuii,alfneo,visvol_v,visvol_q

      do ifld = 1, nfsp
       if(zi(ifld) > 0.) then
        do iy = j2, j5
          do ix = i2, i5
            ix1 = ixp1(ix,iy)
            ix2 = ixm1(ix,iy)
c     First do the viscosity driven by particle flux            
            tempp = b12(ix1,iy)*( up(ix1,iy,ifld) + 
     .                (v2cd(ix1,iy,ifld)+v2ce(ix1,iy,ifld))*
     .                                       rbfbt(ix1,iy)/rrv(ix1,iy) )
            tempm = b12(ix,iy)*( up(ix,iy,ifld) + 
     .                (v2cd(ix,iy,ifld)+v2ce(ix,iy,ifld))*
     .                                       rbfbt(ix,iy)/rrv(ix,iy) )
            if(ix < nx) then
              diffp = (visxneo(ix1,iy,ifld)/rr(ix1,iy))*(tempp-tempm)*
     .                                          gx(ix1,iy)/bsqr(ix1,iy)
            else
              diffp = (visxneo(ix1,iy,ifld)/rr(ix1,iy))*(tempp-tempm)*
     .                                       2.*gx(nx,iy)/bsqr(ix1,iy)
            endif
            tempp = tempm
            tempm = b12(ix2,iy)*( up(ix2,iy,ifld) + 
     .               (v2cd(ix2,iy,ifld)+v2ce(ix2,iy,ifld))*
     .                                       rbfbt(ix2,iy)/rrv(ix2,iy) )
            diffm = (visxneo(ix,iy,ifld)/rr(ix,iy))*(tempp-tempm)*
     .                                             gx(ix,iy)/bsqr(ix,iy)
            visvol_v(ix,iy,ifld)=(4./3.)*rrv(ix,iy)*
     .                                 b32(ix,iy)*(diffp-diffm)*
     .                                   gxf(ix,iy)*volv(ix,iy)
c     Now do the viscosity driven by heat flux            
            tempp = b12(ix1,iy)*( qipar(ix1,iy,ifld) + 
     .                   q2cd(ix1,iy,ifld)*rbfbt(ix1,iy)/rrv(ix1,iy) )
            tempm = b12(ix,iy)*( qipar(ix,iy,ifld) + 
     .                   q2cd(ix,iy,ifld)*rbfbt(ix,iy)/rrv(ix,iy) )
            if(ix < nx) then
              diffp = alfneo(ix1,iy,ifld)*rr(ix1,iy)*(tempp-tempm)*
     .                  gx(ix1,iy) / (nuii(ix1,iy,ifld)*bsqr(ix1,iy))
            else
              diffp = alfneo(ix1,iy,ifld)*rr(ix1,iy)*(tempp-tempm)*
     .                 2.*gx(nx,iy) / (nuii(ix1,iy,ifld)*bsqr(ix1,iy))
            endif
            tempp = tempm
            tempm = b12(ix2,iy)*( qipar(ix2,iy,ifld) + 
     .                   q2cd(ix2,iy,ifld)*rbfbt(ix2,iy)/rrv(ix2,iy) )
            diffm = alfneo(ix,iy,ifld)*rr(ix,iy)*(tempp-tempm)*
     .                   gx(ix,iy) / (nuii(ix,iy,ifld)*bsqr(ix,iy))
            visvol_q(ix,iy,ifld) = rrv(ix,iy)*b32(ix,iy)*
     .                              (diffp-diffm)*gxf(ix,iy)*volv(ix,iy)
          enddo
        enddo
       endif  # if-test on zi(ifld)
      enddo
      return
      end
c ***** End of subroutine upvisneo **********
c ======================================================================
c
      subroutine jvisneo

c ... This subroutine calculates the plasma current from
c ... neoclassical viscosity effects

      implicit none

*  -- local variables
      integer ifld
      #Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real tv,t0,t1,t2,a
      real tempp,tempm,temp0,diffp,diffm,gradx_vpiv,gradx_vpiq

      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nusp
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
                               # iysptrx
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(RZ_grid_info)  # bpol,b12,b32,bsqr,b12ctr
      Use(Bfield)   #rbfbt
      Use(Share)    # nxpt,geometry,nxc,cutlo
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,xlinc,xrinc,ixs1,
                    # j1,j1p,j2,j2p,j3,j4,j5,j6,j5p,j6p,j7,j8,
      Use(Compla)   # ni,up,te,ti,ng,phi,v2cd,
      Use(Comflo)   # qipar
      Use(Conduc)   # visxneo,nuii,alfneo,visvol_v,visvol_q

      do ifld = 1, 1  # limit to main species for now
       if(zi(ifld) > 0.) then #place-holder for when ifld > 1 used
        do iy = j2, j5
          do ix = i2, i5
            ix1 = ixp1(ix,iy)
            ix2 = ixm1(ix,iy)
c     First do the current driven by neoclassical particle flux            
            tempp = up(ix1,iy,ifld) + 
     .                        (v2cd(ix1,iy,ifld)+v2ce(ix1,iy,ifld))*
     .                                     rbfbt(ix1,iy)/rrv(ix1,iy)
            temp0 = up(ix,iy,ifld) + 
     .                        (v2cd(ix2,iy,ifld)+v2ce(ix2,iy,ifld))*
     .                                     rbfbt(ix2,iy)/rrv(ix2,iy)
            tempm = up(ix2,iy,ifld) + 
     .                        (v2cd(ix,iy,ifld)+v2ce(ix,iy,ifld))*
     .                                     rbfbt(ix,iy)/rrv(ix,iy)
            gradx_vpiv = (visxneo(ix,iy,ifld)*b12(ix,iy)*rrv(ix,iy)/3.)
     .                    *( (tempp+temp0)*b12ctr(ix1,iy)-
     .                       (temp0+tempm)*b12ctr(ix,iy) )*0.5*gxf(ix,iy)
            fq2pneo(ix,iy) = -gradx_vpiv*dbm2dy(ix1,iy)
            fqypneo(ix,iy) = 0.5*(rbfbt(ix,iy)+rbfbt(ix1,iy))*
     .                            gradx_vpiv*dbm2dx(ix1,iy)

c     Now do the current driven by neoclassical heat flux            
            tempp = qipar(ix1,iy,ifld) + 
     .                   q2cd(ix1,iy,ifld)*rbfbt(ix1,iy)/rrv(ix1,iy)
            temp0 = qipar(ix,iy,ifld) + 
     .                   q2cd(ix,iy,ifld)*rbfbt(ix,iy)/rrv(ix,iy)
            tempm = qipar(ix2,iy,ifld) + 
     .                   q2cd(ix2,iy,ifld)*rbfbt(ix2,iy)/rrv(ix2,iy)
            gradx_vpiq =( alfneo(ix,iy,ifld)*0.24*
     .                     b12(ix,iy)*rrv(ix,iy)/nuii(ix,iy,ifld) )*
     .                     ( (tempp+temp0)*b12ctr(ix1,iy)-
     .                       (temp0+tempm)*b12ctr(ix,iy) )*0.5*gxf(ix,iy)
            fq2qneo(ix,iy) = -gradx_vpiq*dbm2dy(ix,iy)
            fqyqneo(ix,iy) = 0.5*(rbfbt(ix,iy)+rbfbt(ix1,iy))*
     .                            gradx_vpiq*dbm2dx(ix1,iy)

          enddo
        enddo
       endif  # if-test on zi(ifld)
      enddo
      return
      end
c ***** End of subroutine jvisneo **********
c
c-----------------------------------------------------------------------
c--  Used by ANL for PETSc development -------------------
c-------------------------------------------------------------
      subroutine jacwrite(n, jac, jacj, jaci)

c  This function serves to output the jacobian for viewing purposes
c  The output is the file jacwrite.txt

      integer n,j,k
      real jac(*)
      integer jacj(*), jaci(n+1)

      open(UNIT=88,FILE="jacwrite.txt",STATUS='REPLACE')
 77   format(/)

      write(88,*)"This is the jacobian after some scaling"
      do j=1,n
        do k=jaci(j),jaci(j+1)-1
          write(88,*)j,"  ",jacj(k),"  ",jac(k)
        end do
      end do

      close(88)
      write(*,*)"Jacobian written successfully to jacwrite.txt"
      end
c ***** End of subroutine jacwrite **********

      subroutine jac_calc_interface(neq, t, yl, yldot00, ml, mu, wk,
     .                     nnzmx, jac, ja, ia)

c ... Interface for Jacobian matrix calculation for nksol only(added by. J.Guterl)

      implicit none

c ... Input arguments:
      integer neq      # total number of equations (all grid points)
      real t           # physical time
      real yl(*)       # dependent variables
      real yldot00(neq+2) # right-hand sides evaluated at yl
      integer ml, mu   # lower and upper bandwidths
      integer nnzmx    # maximum number of nonzeros in Jacobian

c ... Work-array argument:
      real wk(neq)     # work space available to this subroutine

c ... Output arguments:
      real jac(nnzmx)     # nonzero Jacobian elements
      integer ja(nnzmx)   # col indices of nonzero Jacobian elements
      integer ia(neq+1)   # pointers to beginning of each row in jac,ja

      Use(ParallelEval)   # ParallelJac

c!omp  if (ParallelJac.eq.1) then
c!omp    call jac_calc_parallel (neq, t, yl, yldot00, ml, mu, wk,
c!omp.             nnzmx, jac, ja, ia)
c!omp  else
         call jac_calc (neq, t, yl, yldot00, ml, mu, wk,
     .               nnzmx, jac, ja, ia)
c!omp  endif
       end subroutine jac_calc_interface

      subroutine Pandf1rhs_interface(neq, time, yl, yldot)
c ... Interface for pandf1 rhs calculation for nksol only (added by. J.Guterl)
          implicit none
          Use(Math_problem_size) # neqmx
          Use(ParallelEval)      # ParallelPandf1
          integer neq
          real time, yl(neqmx),yldot(neq)

c!omp     if (ParallelPandf1.gt.0) then
c!omp         call OMPPandf1Rhs(neq, time, yl, yldot)
c!omp     else
              call pandf1(-1, -1, 0, neq, time, yl, yldot)
c!omp     endif

      end subroutine Pandf1rhs_interface

        subroutine PrintTimingPandf()
            Use(PandfTiming)
            write(*,*) '----- Timing Pandf as eval rhs ----'
            write(*,*) ' - TimePandf:',TotTimePandf
            if (TotTimePandf.gt.0) then
            write(*,*) ' - Convert0:', TotTimeConvert0,TotTimeConvert0/TotTimePandf
            write(*,*) ' - Convert1:', TotTimeConvert1,TotTimeConvert1/TotTimePandf
            write(*,*) ' - Neudif:', TotTimeNeudif,TotTimeNeudif/TotTimePandf
            write(*,*) ' - fd2tra:', TotTimefd2tra,TotTimefd2tra/TotTimePandf
            write(*,*) '-----------------------------------'
            endif
        end subroutine PrintTimingPandf

        real function tick()
        implicit none
            integer :: now, clock_rate
            call system_clock(now,clock_rate)
            tick=real(now)/real(clock_rate)
        end function tick

        real function tock(t)
         implicit none
            real, intent(in) :: t
            integer :: now, clock_rate
            call system_clock(now,clock_rate)

            tock = real(now)/real(clock_rate)-t
        end function tock

