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
cc      real(Size4) ranf

c ... Local variables:
      integer nnz, ii, iv, ii1, ii2, xc, yc, ix, iy
      real yold, dyl, jacelem
      real tick, tock, tsjstor, tsimpjf, dtimpjf
      external tick, tock

ccc      save

c ... Get initial value of system cpu timer.
      if (istimingon .eq. 1) tsjstor = tick()

c ... Pause from BASIS if a ctrl_c is typed
      call ruthere

      ijac(ig) = ijac(ig) + 1

      if ((svrpkg.eq.'nksol') .and. (iprint .ne. 0)) 
     .      write(*,*) ' Updating Jacobian, npe =  ', ijac(ig)

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
c ... and store nonzero elements in compressed sparse column format.
c ... Note that iv is the row (eqn) index and ii is column (var) index

         jcsc(iv) = nnz      # sets index for first Jac. elem. of eqn iv
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
            if (abs(jacelem*sfscal(ii)) .gt. jaccliplim) then
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
      if (istimingon .eq. 1) ttjstor = ttjstor + tock(tsjstor)
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

c ... Local variables:
      integer lowd, ierr, i, idum(1)
      real rcond, dum(1)
      real tsmatfac, tick, tock
      external tick, tock

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
         tsmatfac = tick()
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
         tsmatfac = tick()
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
         tsmatfac = tick()
         call precond5 (neq, ndiag, ndiagm, adiag, wp, rwk2, rwk1,
     .                  iwk3, iwk2, siginel, fmuinel, iwp(3))
      endif

c ... Accumulate cpu time spent here.
 99   ttmatfac = ttmatfac + tock(tsmatfac)
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

c ... Local variables:
      real tick, tock, tsjreorder
      integer i, nfirst
      external tick, tock

c ... Get initial value of system cpu timer.
      tsjreorder = tick()

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
      ttjreorder = ttjreorder + tock(tsjreorder)
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

c ... Local variables:
      real tsjrnorm, tick, tock
      external tick, tock

c ... Get initial value of system cpu timer.
      tsjrnorm = tick()

      if (isrnorm .eq. 1) call roscal (neq, 0, normtype, jac, ja, ia,
     .                               fnormnw, jac, ja, ia)

c ... Accumulate cpu time spent here.
      ttjrnorm = ttjrnorm + tock(tsjrnorm)
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

c ... Local variables:
      integer i
      integer lowd, lbw, ubw
      integer ndiag, ndiagm
      real tsmatsol, tick, tock
      external tick, tock

c ... Get initial value of system cpu timer.
      tsmatsol = tick()

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
      ttmatsol = ttmatsol + tock(tsmatsol)
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
      Use(UEpar)            #ediss,eion,ebind,iigsp,ishosor,fsprd,ishymol,ismolcrm
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
      real ave,a1,a2,t1,t2,thetaix,thetaix2,eta_dup2dy,sv_crumpet
      external sv_crumpet

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
                pmloss(ix,iy) =(1-ismolcrm)*cnsor*(ediss*ev*(0.5*psordis(ix,iy,2))+
     .                      ceisor*eion*ev*(psordis(ix,iy,2)) ) + 
     .                      ismolcrm*cnsor*(
     .                            cmesori*(emolia(ix,iy,1) + emolia(ix,iy,2))
     .                          + cmesore*edisse(ix,iy)
     .                      )
                pmpot(ix,iy) = ismolcrm*ng(ix,iy,2)*vol(ix,iy)*
     .                          sv_crumpet(te(ix,iy), ne(ix,iy), 22)
                pmrada(ix,iy) = ismolcrm*ng(ix,iy,2)*vol(ix,iy)*
     .                          sv_crumpet(te(ix,iy), ne(ix,iy), 23)
                pmradm(ix,iy) = ismolcrm*ng(ix,iy,2)*vol(ix,iy)*
     .                          sv_crumpet(te(ix,iy), ne(ix,iy), 24)
            peirad(ix,iy) = cnsor*( erliz(ix,iy) + erlrc(ix,iy) +
     .                              ebind*ev*psor(ix,iy,1) -
     .                              ebind*ev*psorrg(ix,iy,1) +
     .                              pmloss(ix,iy))

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
            prdiss = prdiss+(1-ismolcrm)*(ediss*ev*(0.5*psordis(ix,iy,2)))+
     .                              ismolcrm*cmesore*edisse(ix,iy)
            pibirth = pibirth+(1-ismolcrm)*(ceisor* eion*ev * (psordis(ix,iy,2)) -
     .                ccoldsor*ng(ix,iy,1)*(1.5*ti(ix,iy)-eion*ev)*
     .                                          nucx(ix,iy,1)*vol(ix,iy) )+
     .                  ismolcrm*( ceisor*cmesore*(emolia(ix,iy,1) + emolia(ix,iy,2))
     .                   - ccoldsor*ng(ix,iy,1)*(1.5*ti(ix,iy)-eion*ev)*
     .                                          nucx(ix,iy,1)*vol(ix,iy) )
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

c!omp if (ParallelJac.eq.1) then
c!omp    call jac_calc_parallel (neq, t, yl, yldot00, ml, mu, wk,
c!omp.             nnzmx, jac, ja, ia)
c!omp else
         call jac_calc (neq, t, yl, yldot00, ml, mu, wk,
     .               nnzmx, jac, ja, ia)
c!omp  endif
       end subroutine jac_calc_interface


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

