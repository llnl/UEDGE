c-----------------------------------------------------------------------
      logical function tstguardc (ix, iy)
c ... Return true if cell (ix,iy) is a guard cell, otherwise false.
      implicit none

      Use(Dim)   # nx,ny
 
      integer ix,iy

      if (ix .lt. 1 .or. ix .gt. nx .or. iy .lt. 1 .or. iy .gt. ny) then
         tstguardc = .true.
      else
         tstguardc = .false.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine guardc
      implicit none
      Use(Dim)            # nx,nxm,ny
      Use(RZ_grid_info)   # rm,zm,psi,br,bz,bpol,bphi,b,b12,b32,b12ctr,bsqr
      Use(Share)          # epslon
      Use(Cdv)

*  -- local variables
      integer ix,iy, i0, i1

*     Define guard cells at the edges by extrapolating the grid

*  -- for ix = nxomit :
      if (nxomit .ge. 0) then      # nxomit > 0 means omitting from left bdry
         i0 = nxomit
         i1 = nxomit+1
      elseif (nxomit .lt. 0) then  # nxomit < 0 means omitting from right bdry
         i0 = 0
         i1 = 1
      endif
      do 13 iy = 1, ny
         rm(i0,iy,1) = rm(i1,iy,1)
     .                - epslon * (rm(i1,iy,2)-rm(i1,iy,1))
         rm(i0,iy,2) = rm(i1,iy,1)
         rm(i0,iy,3) = rm(i1,iy,3)
     .                - epslon * (rm(i1,iy,4)-rm(i1,iy,3))
         rm(i0,iy,4) = rm(i1,iy,3)
         zm(i0,iy,1) = zm(i1,iy,1)
     .                - epslon * (zm(i1,iy,2)-zm(i1,iy,1))
         zm(i0,iy,2) = zm(i1,iy,1)
         zm(i0,iy,3) = zm(i1,iy,3)
     .                - epslon * (zm(i1,iy,4)-zm(i1,iy,3))
         zm(i0,iy,4) = zm(i1,iy,3)
         psi(i0,iy,1) = psi(i1,iy,1)
     .                 - epslon * (psi(i1,iy,2)-psi(i1,iy,1))
         psi(i0,iy,2) = psi(i1,iy,1)
         psi(i0,iy,3) = psi(i1,iy,3)
     .                 - epslon * (psi(i1,iy,4)-psi(i1,iy,3))
         psi(i0,iy,4) = psi(i1,iy,3)
         br(i0,iy,1) = br(i1,iy,1)
     .                - epslon * (br(i1,iy,2)-br(i1,iy,1))
         br(i0,iy,2) = br(i1,iy,1)
         br(i0,iy,3) = br(i1,iy,3)
     .                - epslon * (br(i1,iy,4)-br(i1,iy,3))
         br(i0,iy,4) = br(i1,iy,3)
         bz(i0,iy,1) = bz(i1,iy,1)
     .                - epslon * (bz(i1,iy,2)-bz(i1,iy,1))
         bz(i0,iy,2) = bz(i1,iy,1)
         bz(i0,iy,3) = bz(i1,iy,3)
     .                - epslon * (bz(i1,iy,4)-bz(i1,iy,3))
         bz(i0,iy,4) = bz(i1,iy,3)
         bpol(i0,iy,1) = bpol(i1,iy,1)
     .                  - epslon * (bpol(i1,iy,2)-bpol(i1,iy,1))
         bpol(i0,iy,2) = bpol(i1,iy,1)
         bpol(i0,iy,3) = bpol(i1,iy,3)
     .                  - epslon * (bpol(i1,iy,4)-bpol(i1,iy,3))
         bpol(i0,iy,4) = bpol(i1,iy,3)
         bphi(i0,iy,1) = bphi(i1,iy,1)
     .                  - epslon * (bphi(i1,iy,2)-bphi(i1,iy,1))
         bphi(i0,iy,2) = bphi(i1,iy,1)
         bphi(i0,iy,3) = bphi(i1,iy,3)
     .                  - epslon * (bphi(i1,iy,4)-bphi(i1,iy,3))
         bphi(i0,iy,4) = bphi(i1,iy,3)
         b(i0,iy,1) = b(i1,iy,1)
     .               - epslon * (b(i1,iy,2)-b(i1,iy,1))
         b(i0,iy,2) = b(i1,iy,1)
         b(i0,iy,3) = b(i1,iy,3)
     .               - epslon * (b(i1,iy,4)-b(i1,iy,3))
         b(i0,iy,4) = b(i1,iy,3)
         rm(i0,iy,0) = .25*(rm(i0,iy,1)+rm(i0,iy,2)+rm(i0,iy,3)
     .                     +rm(i0,iy,4))
         zm(i0,iy,0) = .25*(zm(i0,iy,1)+zm(i0,iy,2)+zm(i0,iy,3)
     .                     +zm(i0,iy,4))
         psi(i0,iy,0) = .25*(psi(i0,iy,1)+psi(i0,iy,2)
     .        +psi(i0,iy,3)+psi(i0,iy,4))
         br(i0,iy,0) = .25*(br(i0,iy,1)+br(i0,iy,2)+br(i0,iy,3)
     .                     +br(i0,iy,4))
         bz(i0,iy,0) = .25*(bz(i0,iy,1)+bz(i0,iy,2)+bz(i0,iy,3)
     .                     +bz(i0,iy,4))
         bpol(i0,iy,0) = .25*(bpol(i0,iy,1)+bpol(i0,iy,2)
     .                      +bpol(i0,iy,3)+bpol(i0,iy,4))
         bphi(i0,iy,0) = .25*(bphi(i0,iy,1)+bphi(i0,iy,2)
     .                      +bphi(i0,iy,3)+bphi(i0,iy,4))
         b(i0,iy,0) = .25*(b(i0,iy,1)+b(i0,iy,2)+b(i0,iy,3)+b(i0,iy,4))
 13   continue

*  -- for ix = nxm + 1 :
         if (nxomit .ge. 0) then
            i0 = nxm
            i1 = nxm+1
         elseif (nxomit .lt. 0) then
            i0 = nxm + nxomit
            i1 = nxm + nxomit + 1
         endif
      do 14 iy = 0, ny
         rm(i1,iy,1) = rm(i0,iy,2)
         rm(i1,iy,2) = rm(i0,iy,2)
     .                + epslon * (rm(i0,iy,2)-rm(i0,iy,1))
         rm(i1,iy,3) = rm(i0,iy,4)
         rm(i1,iy,4) = rm(i0,iy,4)
     .                + epslon * (rm(i0,iy,4)-rm(i0,iy,3))
         zm(i1,iy,1) = zm(i0,iy,2)
         zm(i1,iy,2) = zm(i0,iy,2)
     .                + epslon * (zm(i0,iy,2)-zm(i0,iy,1))
         zm(i1,iy,3) = zm(i0,iy,4)
         zm(i1,iy,4) = zm(i0,iy,4)
     .                + epslon * (zm(i0,iy,4)-zm(i0,iy,3))
         psi(i1,iy,1) = psi(i0,iy,2)
         psi(i1,iy,2) = psi(i0,iy,2)
     .                + epslon * (psi(i0,iy,2)-psi(i0,iy,1))
         psi(i1,iy,3) = psi(i0,iy,4)
         psi(i1,iy,4) = psi(i0,iy,4)
     .                + epslon * (psi(i0,iy,4)-psi(i0,iy,3))
         br(i1,iy,1) = br(i0,iy,2)
         br(i1,iy,2) = br(i0,iy,2)
     .                + epslon * (br(i0,iy,2)-br(i0,iy,1))
         br(i1,iy,3) = br(i0,iy,4)
         br(i1,iy,4) = br(i0,iy,4)
     .                + epslon * (br(i0,iy,4)-br(i0,iy,3))
         bz(i1,iy,1) = bz(i0,iy,2)
         bz(i1,iy,2) = bz(i0,iy,2)
     .                + epslon * (bz(i0,iy,2)-bz(i0,iy,1))
         bz(i1,iy,3) = bz(i0,iy,4)
         bz(i1,iy,4) = bz(i0,iy,4)
     .                + epslon * (bz(i0,iy,4)-bz(i0,iy,3))
         bpol(i1,iy,1) = bpol(i0,iy,2)
         bpol(i1,iy,2) = bpol(i0,iy,2)
     .                + epslon * (bpol(i0,iy,2)-bpol(i0,iy,1))
         bpol(i1,iy,3) = bpol(i0,iy,4)
         bpol(i1,iy,4) = bpol(i0,iy,4)
     .                + epslon * (bpol(i0,iy,4)-bpol(i0,iy,3))
         bphi(i1,iy,1) = bphi(i0,iy,2)
         bphi(i1,iy,2) = bphi(i0,iy,2)
     .                + epslon * (bphi(i0,iy,2)-bphi(i0,iy,1))
         bphi(i1,iy,3) = bphi(i0,iy,4)
         bphi(i1,iy,4) = bphi(i0,iy,4)
     .                + epslon * (bphi(i0,iy,4)-bphi(i0,iy,3))
         b(i1,iy,1) = b(i0,iy,2)
         b(i1,iy,2) = b(i0,iy,2)
     .                + epslon * (b(i0,iy,2)-b(i0,iy,1))
         b(i1,iy,3) = b(i0,iy,4)
         b(i1,iy,4) = b(i0,iy,4)
     .                + epslon * (b(i0,iy,4)-b(i0,iy,3))
         rm(i1,iy,0) = .25*(rm(i1,iy,1)+rm(i1,iy,2)
     .                       +rm(i1,iy,3)+rm(i1,iy,4))
         zm(i1,iy,0) = .25*(zm(i1,iy,1)+zm(i1,iy,2)
     .                       +zm(i1,iy,3)+zm(i1,iy,4))
         psi(i1,iy,0) = .25*(psi(i1,iy,1)+psi(i1,iy,2)
     .                        +psi(i1,iy,3)+psi(i1,iy,4))
         br(i1,iy,0) = .25*(br(i1,iy,1)+br(i1,iy,2)
     .                       +br(i1,iy,3)+br(i1,iy,4))
         bz(i1,iy,0) = .25*(bz(i1,iy,1)+bz(i1,iy,2)
     .                       +bz(i1,iy,3)+bz(i1,iy,4))
         bpol(i1,iy,0) = .25*(bpol(i1,iy,1)+bpol(i1,iy,2)
     .                         +bpol(i1,iy,3)+bpol(i1,iy,4))
         bphi(i1,iy,0) = .25*(bphi(i1,iy,1)+bphi(i1,iy,2)
     .                         +bphi(i1,iy,3)+bphi(i1,iy,4))
         b(i1,iy,0) = .25*(b(i1,iy,1)+b(i1,iy,2)
     .                      +b(i1,iy,3)+b(i1,iy,4))
 14   continue

*  -- for iy = 0 :
      do 15 ix = 0, nxm+1
         rm(ix,0,1) = rm(ix,1,1)
     .                - epslon * (rm(ix,1,3)-rm(ix,1,1))
         rm(ix,0,2) = rm(ix,1,2)
     .                - epslon * (rm(ix,1,4)-rm(ix,1,2))
         rm(ix,0,3) = rm(ix,1,1)
         rm(ix,0,4) = rm(ix,1,2)
         zm(ix,0,1) = zm(ix,1,1)
     .                - epslon * (zm(ix,1,3)-zm(ix,1,1))
         zm(ix,0,2) = zm(ix,1,2)
     .                - epslon * (zm(ix,1,4)-zm(ix,1,2))
         zm(ix,0,3) = zm(ix,1,1)
         zm(ix,0,4) = zm(ix,1,2)
         psi(ix,0,1) = psi(ix,1,1)
     .                - epslon * (psi(ix,1,3)-psi(ix,1,1))
         psi(ix,0,2) = psi(ix,1,2)
     .                - epslon * (psi(ix,1,4)-psi(ix,1,2))
         psi(ix,0,3) = psi(ix,1,1)
         psi(ix,0,4) = psi(ix,1,2)
         br(ix,0,1) = br(ix,1,1)
     .                - epslon * (br(ix,1,3)-br(ix,1,1))
         br(ix,0,2) = br(ix,1,2)
     .                - epslon * (br(ix,1,4)-br(ix,1,2))
         br(ix,0,3) = br(ix,1,1)
         br(ix,0,4) = br(ix,1,2)
         bz(ix,0,1) = bz(ix,1,1)
     .                - epslon * (bz(ix,1,3)-bz(ix,1,1))
         bz(ix,0,2) = bz(ix,1,2)
     .                - epslon * (bz(ix,1,4)-bz(ix,1,2))
         bz(ix,0,3) = bz(ix,1,1)
         bz(ix,0,4) = bz(ix,1,2)
         bpol(ix,0,1) = bpol(ix,1,1)
     .                - epslon * (bpol(ix,1,3)-bpol(ix,1,1))
         bpol(ix,0,2) = bpol(ix,1,2)
     .                - epslon * (bpol(ix,1,4)-bpol(ix,1,2))
         bpol(ix,0,3) = bpol(ix,1,1)
         bpol(ix,0,4) = bpol(ix,1,2)
         bphi(ix,0,1) = bphi(ix,1,1)
     .                - epslon * (bphi(ix,1,3)-bphi(ix,1,1))
         bphi(ix,0,2) = bphi(ix,1,2)
     .                - epslon * (bphi(ix,1,4)-bphi(ix,1,2))
         bphi(ix,0,3) = bphi(ix,1,1)
         bphi(ix,0,4) = bphi(ix,1,2)
         b(ix,0,1) = b(ix,1,1)
     .                - epslon * (b(ix,1,3)-b(ix,1,1))
         b(ix,0,2) = b(ix,1,2)
     .                - epslon * (b(ix,1,4)-b(ix,1,2))
         b(ix,0,3) = b(ix,1,1)
         b(ix,0,4) = b(ix,1,2)
         rm(ix,0,0) = .25*(rm(ix,0,1)+rm(ix,0,2)+rm(ix,0,3)+rm(ix,0,4))
         zm(ix,0,0) = .25*(zm(ix,0,1)+zm(ix,0,2)+zm(ix,0,3)+zm(ix,0,4))
         psi(ix,0,0) = .25*(psi(ix,0,1)+psi(ix,0,2)
     .        +psi(ix,0,3)+psi(ix,0,4))
         br(ix,0,0) = .25*(br(ix,0,1)+br(ix,0,2)+br(ix,0,3)+br(ix,0,4))
         bz(ix,0,0) = .25*(bz(ix,0,1)+bz(ix,0,2)+bz(ix,0,3)+bz(ix,0,4))
         bpol(ix,0,0) = .25*(bpol(ix,0,1)+bpol(ix,0,2)
     .                      +bpol(ix,0,3)+bpol(ix,0,4))
         bphi(ix,0,0) = .25*(bphi(ix,0,1)+bphi(ix,0,2)
     .                      +bphi(ix,0,3)+bphi(ix,0,4))
         b(ix,0,0) = .25*(b(ix,0,1)+b(ix,0,2)+b(ix,0,3)+b(ix,0,4))
 15   continue

*  -- for iy = ny + 1 :
      do 16 ix = 0, nxm+1
         rm(ix,ny+1,1) = rm(ix,ny,3)
         rm(ix,ny+1,2) = rm(ix,ny,4)
         rm(ix,ny+1,3) = rm(ix,ny,3)
     .               + epslon * (rm(ix,ny,3)-rm(ix,ny,1))
         rm(ix,ny+1,4) = rm(ix,ny,4)
     .               + epslon * (rm(ix,ny,4)-rm(ix,ny,2))
         zm(ix,ny+1,1) = zm(ix,ny,3)
         zm(ix,ny+1,2) = zm(ix,ny,4)
         zm(ix,ny+1,3) = zm(ix,ny,3)
     .               + epslon * (zm(ix,ny,3)-zm(ix,ny,1))
         zm(ix,ny+1,4) = zm(ix,ny,4)
     .               + epslon * (zm(ix,ny,4)-zm(ix,ny,2))
         psi(ix,ny+1,1) = psi(ix,ny,3)
         psi(ix,ny+1,2) = psi(ix,ny,4)
         psi(ix,ny+1,3) = psi(ix,ny,3)
     .               + epslon * (psi(ix,ny,3)-psi(ix,ny,1))
         psi(ix,ny+1,4) = psi(ix,ny,4)
     .               + epslon * (psi(ix,ny,4)-psi(ix,ny,2))
         br(ix,ny+1,1) = br(ix,ny,3)
         br(ix,ny+1,2) = br(ix,ny,4)
         br(ix,ny+1,3) = br(ix,ny,3)
     .               + epslon * (br(ix,ny,3)-br(ix,ny,1))
         br(ix,ny+1,4) = br(ix,ny,4)
     .               + epslon * (br(ix,ny,4)-br(ix,ny,2))
         bz(ix,ny+1,1) = bz(ix,ny,3)
         bz(ix,ny+1,2) = bz(ix,ny,4)
         bz(ix,ny+1,3) = bz(ix,ny,3)
     .               + epslon * (bz(ix,ny,3)-bz(ix,ny,1))
         bz(ix,ny+1,4) = bz(ix,ny,4)
     .               + epslon * (bz(ix,ny,4)-bz(ix,ny,2))
         bpol(ix,ny+1,1) = bpol(ix,ny,3)
         bpol(ix,ny+1,2) = bpol(ix,ny,4)
         bpol(ix,ny+1,3) = bpol(ix,ny,3)
     .               + epslon * (bpol(ix,ny,3)-bpol(ix,ny,1))
         bpol(ix,ny+1,4) = bpol(ix,ny,4)
     .               + epslon * (bpol(ix,ny,4)-bpol(ix,ny,2))
         bphi(ix,ny+1,1) = bphi(ix,ny,3)
         bphi(ix,ny+1,2) = bphi(ix,ny,4)
         bphi(ix,ny+1,3) = bphi(ix,ny,3)
     .               + epslon * (bphi(ix,ny,3)-bphi(ix,ny,1))
         bphi(ix,ny+1,4) = bphi(ix,ny,4)
     .               + epslon * (bphi(ix,ny,4)-bphi(ix,ny,2))
         b(ix,ny+1,1) = b(ix,ny,3)
         b(ix,ny+1,2) = b(ix,ny,4)
         b(ix,ny+1,3) = b(ix,ny,3)
     .               + epslon * (b(ix,ny,3)-b(ix,ny,1))
         b(ix,ny+1,4) = b(ix,ny,4)
     .               + epslon * (b(ix,ny,4)-b(ix,ny,2))
         rm(ix,ny+1,0) = .25*(rm(ix,ny+1,1)+rm(ix,ny+1,2)
     .                       +rm(ix,ny+1,3)+rm(ix,ny+1,4))
         zm(ix,ny+1,0) = .25*(zm(ix,ny+1,1)+zm(ix,ny+1,2)
     .                       +zm(ix,ny+1,3)+zm(ix,ny+1,4))
         psi(ix,ny+1,0) = .25*(psi(ix,ny+1,1)+psi(ix,ny+1,2)
     .                        +psi(ix,ny+1,3)+psi(ix,ny+1,4))
         br(ix,ny+1,0) = .25*(br(ix,ny+1,1)+br(ix,ny+1,2)
     .                       +br(ix,ny+1,3)+br(ix,ny+1,4))
         bz(ix,ny+1,0) = .25*(bz(ix,ny+1,1)+bz(ix,ny+1,2)
     .                       +bz(ix,ny+1,3)+bz(ix,ny+1,4))
         bpol(ix,ny+1,0) = .25*(bpol(ix,ny+1,1)+bpol(ix,ny+1,2)
     .                         +bpol(ix,ny+1,3)+bpol(ix,ny+1,4))
         bphi(ix,ny+1,0) = .25*(bphi(ix,ny+1,1)+bphi(ix,ny+1,2)
     .                         +bphi(ix,ny+1,3)+bphi(ix,ny+1,4))
         b(ix,ny+1,0) = .25*(b(ix,ny+1,1)+b(ix,ny+1,2)
     .                      +b(ix,ny+1,3)+b(ix,ny+1,4))
 16   continue

      return
      end

c----------------------------------------------------------------------c

      subroutine mpguardc
      implicit none
      Use(Dim)            # ny
      Use(RZ_grid_info)   # rm,zm,psi,br,bz,bpol,bphi,b
      Use(Share)          # nxc

c     local variables
      integer iy

c     Define guard cells at inboard midplane of double-null via symmetry
      do iy = 0, ny+1
         rm(nxc,iy,1) = rm(nxc-1,iy,2)
         zm(nxc,iy,1) = zm(nxc-1,iy,2)
         psi(nxc,iy,1) = psi(nxc-1,iy,2)
         br(nxc,iy,1) = -br(nxc-1,iy,2)
         bz(nxc,iy,1) = bz(nxc-1,iy,2)
         bpol(nxc,iy,1) = bpol(nxc-1,iy,2)
         bphi(nxc,iy,1) = bphi(nxc-1,iy,2)
         b(nxc,iy,1) = b(nxc-1,iy,2)

         rm(nxc,iy,2) = rm(nxc-1,iy,1)
         zm(nxc,iy,2) = 2*zm(nxc-1,iy,2) - zm(nxc-1,iy,1)
         psi(nxc,iy,2) = psi(nxc-1,iy,1)
         br(nxc,iy,2) = -br(nxc-1,iy,1)
         bz(nxc,iy,2) = bz(nxc-1,iy,1)
         bpol(nxc,iy,2) = bpol(nxc-1,iy,1)
         bphi(nxc,iy,2) = bphi(nxc-1,iy,1)
         b(nxc,iy,2) = b(nxc-1,iy,1)

         rm(nxc,iy,3) = rm(nxc-1,iy,4)
         zm(nxc,iy,3) = zm(nxc-1,iy,4)
         psi(nxc,iy,3) = psi(nxc-1,iy,4)
         br(nxc,iy,3) = -br(nxc-1,iy,4)
         bz(nxc,iy,3) = bz(nxc-1,iy,4)
         bpol(nxc,iy,3) = bpol(nxc-1,iy,4)
         bphi(nxc,iy,3) = bphi(nxc-1,iy,4)
         b(nxc,iy,3) = b(nxc-1,iy,4)

         rm(nxc,iy,4) = rm(nxc-1,iy,3)
         zm(nxc,iy,4) = 2*zm(nxc-1,iy,4) - zm(nxc-1,iy,3)
         psi(nxc,iy,4) = psi(nxc-1,iy,3)
         br(nxc,iy,4) = -br(nxc-1,iy,3)
         bz(nxc,iy,4) = bz(nxc-1,iy,3)
         bpol(nxc,iy,4) = bpol(nxc-1,iy,3)
         bphi(nxc,iy,4) = bphi(nxc-1,iy,3)
         b(nxc,iy,4) = b(nxc-1,iy,3)

         rm(nxc,iy,0) = .25*(rm(nxc,iy,1)+rm(nxc,iy,2)
     .                      +rm(nxc,iy,3)+rm(nxc,iy,4))
         zm(nxc,iy,0) = .25*(zm(nxc,iy,1)+zm(nxc,iy,2)
     .                      +zm(nxc,iy,3)+zm(nxc,iy,4))
         psi(nxc,iy,0) = .25*(psi(nxc,iy,1)+psi(nxc,iy,2)
     .                      +psi(nxc,iy,3)+psi(nxc,iy,4))
         br(nxc,iy,0) = .25*(br(nxc,iy,1)+br(nxc,iy,2)
     .                      +br(nxc,iy,3)+br(nxc,iy,4))
         bz(nxc,iy,0) = .25*(bz(nxc,iy,1)+bz(nxc,iy,2)
     .                      +bz(nxc,iy,3)+bz(nxc,iy,4))
         bpol(nxc,iy,0) = .25*(bpol(nxc,iy,1)+bpol(nxc,iy,2)
     .                      +bpol(nxc,iy,3)+bpol(nxc,iy,4))
         bphi(nxc,iy,0) = .25*(bphi(nxc,iy,1)+bphi(nxc,iy,2)
     .                      +bphi(nxc,iy,3)+bphi(nxc,iy,4))
         b(nxc,iy,0) = .25*(b(nxc,iy,1)+b(nxc,iy,2)
     .                      +b(nxc,iy,3)+b(nxc,iy,4))

      enddo

*     Define guard cells at outboard midplane of double-null via symmetry
      do iy = 0, ny+1
         rm(nxc+1,iy,2) = rm(nxc+2,iy,1)
         zm(nxc+1,iy,2) = zm(nxc+2,iy,1)
         psi(nxc+1,iy,2) = psi(nxc+2,iy,1)
         br(nxc+1,iy,2) = -br(nxc+2,iy,1)
         bz(nxc+1,iy,2) = bz(nxc+2,iy,1)
         bpol(nxc+1,iy,2) = bpol(nxc+2,iy,1)
         bphi(nxc+1,iy,2) = bphi(nxc+2,iy,1)
         b(nxc+1,iy,2) = b(nxc+2,iy,1)

         rm(nxc+1,iy,1) = rm(nxc+2,iy,2)
         zm(nxc+1,iy,1) = 2*zm(nxc+2,iy,1) - zm(nxc+2,iy,2)
         psi(nxc+1,iy,1) = psi(nxc+2,iy,2)
         br(nxc+1,iy,1) = -br(nxc+2,iy,2)
         bz(nxc+1,iy,1) = bz(nxc+2,iy,2)
         bpol(nxc+1,iy,1) = bpol(nxc+2,iy,2)
         bphi(nxc+1,iy,1) = bphi(nxc+2,iy,2)
         b(nxc+1,iy,1) = b(nxc+2,iy,2)

         rm(nxc+1,iy,4) = rm(nxc+2,iy,3)
         zm(nxc+1,iy,4) = zm(nxc+2,iy,3)
         psi(nxc+1,iy,4) = psi(nxc+2,iy,3)
         br(nxc+1,iy,4) = -br(nxc+2,iy,3)
         bz(nxc+1,iy,4) = bz(nxc+2,iy,3)
         bpol(nxc+1,iy,4) = bpol(nxc+2,iy,3)
         bphi(nxc+1,iy,4) = bphi(nxc+2,iy,3)
         b(nxc+1,iy,4) = b(nxc+2,iy,3)

         rm(nxc+1,iy,3) = rm(nxc+2,iy,4)
         zm(nxc+1,iy,3) = 2*zm(nxc+2,iy,3) - zm(nxc+2,iy,4)
         psi(nxc+1,iy,3) = psi(nxc+2,iy,4)
         br(nxc+1,iy,3) = -br(nxc+2,iy,4)
         bz(nxc+1,iy,3) = bz(nxc+2,iy,4)
         bpol(nxc+1,iy,3) = bpol(nxc+2,iy,4)
         bphi(nxc+1,iy,3) = bphi(nxc+2,iy,4)
         b(nxc+1,iy,3) = b(nxc+2,iy,4)

         rm(nxc+1,iy,0) = .25*(rm(nxc+1,iy,1)+rm(nxc+1,iy,2)
     .                      +rm(nxc+1,iy,3)+rm(nxc+1,iy,4))
         zm(nxc+1,iy,0) = .25*(zm(nxc+1,iy,1)+zm(nxc+1,iy,2)
     .                      +zm(nxc+1,iy,3)+zm(nxc+1,iy,4))
         psi(nxc+1,iy,0) = .25*(psi(nxc+1,iy,1)+psi(nxc+1,iy,2)
     .                      +psi(nxc+1,iy,3)+psi(nxc+1,iy,4))
         br(nxc+1,iy,0) = .25*(br(nxc+1,iy,1)+br(nxc+1,iy,2)
     .                      +br(nxc+1,iy,3)+br(nxc+1,iy,4))
         bz(nxc+1,iy,0) = .25*(bz(nxc+1,iy,1)+bz(nxc+1,iy,2)
     .                      +bz(nxc+1,iy,3)+bz(nxc+1,iy,4))
         bpol(nxc+1,iy,0) = .25*(bpol(nxc+1,iy,1)+bpol(nxc+1,iy,2)
     .                      +bpol(nxc+1,iy,3)+bpol(nxc+1,iy,4))
         bphi(nxc+1,iy,0) = .25*(bphi(nxc+1,iy,1)+bphi(nxc+1,iy,2)
     .                      +bphi(nxc+1,iy,3)+bphi(nxc+1,iy,4))
         b(nxc+1,iy,0) = .25*(b(nxc+1,iy,1)+b(nxc+1,iy,2)
     .                      +b(nxc+1,iy,3)+b(nxc+1,iy,4))

      enddo

      return
      end

c  ---------------------------------------------------------------------
      subroutine globalmesh

*     globalmesh defines the (R,Z) coordinates and magnetic field values
*     for the global mesh when domain decomposition is used (similar to first
*     section of nphygeo)

      implicit none
      Use(Dim)            # nx,ny
      Use(Math_problem_size) # neqmx
      Use(Share)          # nxomit,nxc,nxleg,nxcore,geometry,ismpsym
      Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2,iysptrx
      Use(RZ_grid_info)   # rm,zm,psi,br,bz,bpol,bphi,b
      Use(RZ_grid_global) # rmg,zmg,psig,brg,bzg,bpolg,bphig,bg
      Use(UEint)          # mhdgeo,gengrid
      Use(UEpar)          # thetar
      Use(Phyvar)         # pi
      Use(Bfield)         # b0old
      Use(Flags)          # iprint
      Use(Comgeo)         # area_core
*  -- local scalars --
      integer nj, ij, ix, iy, jx
      real str, ctr, rm0, zm0, dxc, dz, dacore
      character*200 fname
      character*60 runid

*=======================================================================
*//computation//

      if ((geometry.eq."dnbot").or.((geometry(1:9).eq."snowflake")
     .     .and.isudsym==1)) then
        isudsym=1
      else
        isudsym=0
      endif

c ... Generate grid information or read it from file 'gridue'.

      if(mhdgeo .eq. 1) then
        if (ismmon.ne.0 .and. isnonog.eq.0) then
           call remark('*** WARNING: ismmon.ne.0 BUT isnonog=0 ****')
        endif
        if(gengrid .eq. 1) then
c            write(6,*) "Calling flxrun in globalmesh."
            call flxrun
            call grdrun
         else
            fname = trim(GridFileName)
            call readgrid(fname, runid)
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) 'Read file "', fname, '" with runid:  ', runid
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*)
         endif
      elseif (mhdgeo .eq. 2) then
         if (gengrid == 1) then
            call torangrd
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo=2: Circ toroidal annulus generated *****'
         else
            fname = trim(GridFileName)
            call readgrid(fname, runid)
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) 'Read file "', fname, '" with runid:  ', runid
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*)
         endif
      elseif (mhdgeo .eq. 0) then
         call idealgrd
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo=0: cylindrical grid generated *****'
      elseif (mhdgeo .eq. -1) then
         call idealgrd
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo=-1: cartesian grid generated *****'
      elseif (mhdgeo .eq. -2) then
         call mirrorgrd
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo=-2: mag mirror grid generated *****'
      else
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo < -1: reading grid from file *****'
         fname = trim(GridFileName)
         call readgrid(fname, runid)
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) 'Read file "', fname, '" with runid:  ', runid
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*)
      endif

c...  Reset separatrix index if nyomitmx > 0
      if (nyomitmx >= nysol+nyout) then
        do jx = 1, nxpt
          iysptrx1(jx) = min( ny, iysptrx1(jx) )
          iysptrx2(jx) = min( ny, iysptrx2(jx) )
        enddo
      endif

c.... Set iy index of last closed flux surface
      iysptrx = nym
      do jx=1,nxpt
         iysptrx = min( iysptrx, iysptrx1(jx) )
         iysptrx = min( iysptrx, iysptrx2(jx) )
      enddo

c...  If new grid is generated, initialize b0old=1
      b0old = 1.0

c...  Here we reset ixpt1,2 because it is different in the grid package
c...  for nxomit > 0, and the grid package value is read from file gridue.
      nxc = nxleg(1) + nxcore(1)
      nj = 0
      if (nxomit .gt. 0) then
        do jx = 1, nxpt
          ixpt1(jx) = max(-1, ixpt1(jx) - nxomit)  # was max(0,
          ixpt2(jx) = ixpt2(jx) - nxomit
        enddo
	  nxc = max(nxc - nxomit,0)  # was ,1
        nj = nxomit
      endif

c.... Set indirect addressing arrays for x-direction
      call set_indirect_address(1)  #should be (1)

*     Rotate geometry by angle thetar (degrees) for diagnostic test
      if (thetar .ne. 0.) then
         str = sin(thetar*2*pi/360)
         ctr = cos(thetar*2*pi/360)
         do iy = 0, ny+1
            do ix = 0, nx+1
               do ij = 0, 4
                  rm0 = rm(ix+nj,iy,ij)
                  zm0 = zm(ix+nj,iy,ij)
                  rm(ix+nj,iy,ij) = rm0*ctr - zm0*str
                  zm(ix+nj,iy,ij) = rm0*str + zm0*ctr
               enddo
            enddo
         enddo 
      endif           

c     Re-define midplane guard cells for symmetric double-null
      if (isudsym==1 .and. ismpsym==1) call mpguardc

*     Define guard cells around the edge of the mesh --
      call guardc

      call gallot("RZ_grid_global",0)

*     Copy grid informtion into global arrays as local rm, etc will be
*     reallocated
      do ij = 0, 4
        call s2copy (nx+2,ny+2,rm(0,0,ij),1,nx+2,rmg(0,0,ij),1,nx+2)
        call s2copy (nx+2,ny+2,zm(0,0,ij),1,nx+2,zmg(0,0,ij),1,nx+2)
        call s2copy (nx+2,ny+2,br(0,0,ij),1,nx+2,brg(0,0,ij),1,nx+2)
        call s2copy (nx+2,ny+2,bz(0,0,ij),1,nx+2,bzg(0,0,ij),1,nx+2)
        call s2copy (nx+2,ny+2,bpol(0,0,ij),1,nx+2,bpolg(0,0,ij),1,nx+2)
        call s2copy (nx+2,ny+2,bphi(0,0,ij),1,nx+2,bphig(0,0,ij),1,nx+2)
        call s2copy (nx+2,ny+2,b(0,0,ij),1,nx+2,bg(0,0,ij),1,nx+2)
        call s2copy (nx+2,ny+2,psi(0,0,ij),1,nx+2,psig(0,0,ij),1,nx+2)
      enddo

c...  Compute surface area of core boundary
      nj = nxomit
      area_core = 0.
      do jx = 1, nxpt
        do ix = ixpt1(jx)+1, ixpt2(jx)
           dxc = sqrt((rm(ix+nj,1,2)-rm(ix+nj,1,1))**2
     .              + (zm(ix+nj,1,2)-zm(ix+nj,1,1))**2)
           dz = 2 * pi * 0.5*(rm(ix+nj,1,2)+rm(ix+nj,1,1))
           if (mhdgeo .le. -1) dz = 1.
             dacore = dxc * dz
          area_core = area_core + dacore
        enddo				  
      enddo

      return
      end
c ***  End of subroutine globalmesh *****
     
c-----------------------------------------------------------------------
      subroutine nphygeo

*     NPHYGEO defines the geometry for the edge plasma code by
*     reading grid information from a file.

      implicit none
      Use(Dim)            # nx,ny,nxm,nym,nxpt
      Use(Share)          # nxleg,nxcore,nxomit,isnonog,ismmon,
                          # geometry,ismpsym,simagxs,sibdrys
      Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2,iysptrx
      Use(Cut_indices)    # ixcut1,ixcut2,ixcut3,ixcut4
      Use(Aux)            # ixmp
      Use(Phyvar)         # pi
      Use(Selec)          # ixp1,ixm1
      Use(Comgeo)         # isxptx,isxpty
      Use(Noggeo)         # vtag,angfx
      Use(Bfield)         # btot,rbfbt,rbfbt2,b0old,rbpol,dbm2dx,dbm2dy
      Use(RZ_grid_info)   # rm,zm,rv,zv,bphi,bpol,b,b12,b12ctr,b32,bsqr
      Use(RZ_cell_info)   # rcn,zcn,rfn,zfn,rcv,zcv,rfv,zfv
      Use(UEint)          # lots ...
      Use(UEpar)          # is1D_gbx,ixgb,xgbx,r0slab,thetar,issvyxpt0
                          # isxmpog
      Use(Bcond)          # isfixlb,isfixrb
      Use(Math_problem_size)  # neqmx
      Use(Flags)          # iprint

*  -- local scalars --
      integer nj, iu, ik, ij, jx, iysi, iyso, iyp1, ix_last_core_cell,
     .        ixpt2_1temp, nxpt_temp
      #Former Aux module variables
      integer ix,iy,ix1,ix3,ix2
      data nj/0/
      real dxc, dyc, dz, str, ctr, rm0, zm0, cossr, cossp, s_bphi, rmmax,
     .     lcon_wk1, lcon_wk2
      real bsqrvoltot,voltot,bsqrave,b_yface,b_xface
      character*200 fname
      character*60 runid

*=======================================================================
*//computation//

c ... Generate grid information or read it from file 'gridue'.


        if (manualgrid == 0) then

      if(mhdgeo .eq. 1) then
         if (ismmon.ne.0 .and. isnonog.eq.0) then
           call remark('*** WARNING: ismmon.ne.0 BUT isnonog=0 ****')
         endif
         ixpt2_1temp = ixpt2(1)
         nxpt_temp = nxpt
         if(gengrid .eq. 1) then
            if(geometry=="isoleg") then
              geometry = "dnbot" #Used only to computer mesh, the set back
              ixpt2(1) = nxleg(1)+nxcore(1)+nxcore(2)
              call flxrun
              call grdrun
              geometry = "isoleg"
            else
              call flxrun
              call grdrun
            endif
c ...   Need to fix ixpt2(1) modified by grdrun is geometry=isoleg
            if(geometry=="isoleg") then
              ixpt2(1) = ixpt2_1temp
              nxpt = nxpt_temp
            endif
         else
            fname = trim(GridFileName)
            call readgrid(fname, runid)
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) 'Read file "', fname, '" with runid:  ', runid
c ...  now that the grid is read in, we can manipulate dnull for nxomit>0
            if (geometry=="dnull" .and. nxpt==2) then
              if (nxomit >= ixlb(2)) then
                call remark("*** nxomit>0: do outer quad as single-null")
                geometry = "snull"
                nxpt = 1
	        write(*,*) "ixpt1, ixpt2 = ", ixpt1, ixpt2
                write(*,*)
                ixpt1(1) = ixpt1(1)   # ixpt1,2 are shifted by nxomit later
                ixpt2(1) = ixpt2(2)
                ixlb(1) = 0
              endif
            endif
         endif
      elseif (mhdgeo .eq. 2) then
         if (gengrid == 1) then
            call torangrd
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '*** mhdgeo=2: Circ toroidal annulus ***'
         else
            fname = trim(GridFileName)
            call readgrid(fname, runid)
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) 'Read file "', fname, '" with runid:  ', runid
            if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*)
         endif
      elseif (mhdgeo .eq. 0) then
         call idealgrd
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo=0: cylindrical grid generated *****'
      elseif (mhdgeo .eq. -1) then
         call idealgrd
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo=-1: cartesian grid generated *****'
      elseif (mhdgeo .eq. -2) then
         call mirrorgrd
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo=-2: mag mirror grid generated *****'
      else
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) '**** mhdgeo < -1: reading grid from file *****'
         fname = trim(GridFileName)
         call readgrid(fname, runid)
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*) 'Read file "', fname, '" with runid:  ', runid
         if (iprint*(1-isgriduehdf5) .ne. 0) write(*,*)
      endif

        endif # End test on manualgrid

c...  Reset separatrix index if nyomitmx > 0
      if (nyomitmx >= nysol+nyout) then
        do jx = 1, nxpt
          iysptrx1(jx) = min( ny, iysptrx1(jx) )
          iysptrx2(jx) = min( ny, iysptrx2(jx) )
        enddo
      endif

c.... Set iy index of last closed flux surface
      iysptrx = nym
      do jx=1,nxpt
         iysptrx = min( iysptrx, iysptrx1(jx) )
         iysptrx = min( iysptrx, iysptrx2(jx) )
      enddo

c...  If new grid is generated, initialize b0old=1
      b0old = 1.0

c...  Here we reset ixpt1,2 because it is different in the grid package
c...  for nxomit > 0, and the grid package value is read from file gridue.
      if ((gengrid==0) .and.                # MER 28 Apr 2015
     .    (geometry=="dnXtarget")) then     # nxc is read from gridue file
         continue
      else
         nxc = nxleg(1) + nxcore(1) + 2*nxxpt
      endif
      nj = 0
      if (nxomit .gt. 0) then
        do jx = 1, nxpt
          ixpt1(jx) = max(-1, ixpt1(jx) - nxomit)  # was max(0,
          ixpt2(jx) = ixpt2(jx) - nxomit
          ixrb(1) = nxm - abs(nxomit)
        enddo
        nxc = max(nxc - nxomit,0)  # was,1)
        nj = nxomit
      endif

c.... Set indirect addressing arrays for x-direction
      call set_indirect_address(0)

*     Rotate geometry by angle thetar (degrees) for diagnostic test
      if (thetar .ne. 0.) then
         str = sin(thetar*2*pi/360)
         ctr = cos(thetar*2*pi/360)
         do iy = 0, ny+1
            do ix = 0, nx+1
               do ij = 0, 4
                  rm0 = rm(ix+nj,iy,ij)
                  zm0 = zm(ix+nj,iy,ij)
                  rm(ix+nj,iy,ij) = rm0*ctr - zm0*str
                  zm(ix+nj,iy,ij) = rm0*str + zm0*ctr
               enddo
            enddo
         enddo 
      endif           

c     Re-define midplane guard cells for symmetric double-null
      if (isudsym==1.and.ismpsym==1.and.nxc>1) call mpguardc

c...  For 1D poloidal flux-tube geometry (nysol=1), reset cell vertex
c...  locations near ix=ixpt2 on iy=1 surface if nxpt2msh or nxpt2psh > 1
      if (nxpt2msh+nxpt2psh > 0 .or. nxpt1msh+nxpt1psh > 0) call reset1dmeshpt

*     Define guard cells around the edge of the mesh --
      call guardc

*----------------------------------------------------------------------
*  -- Define the density and velocity cell data in (R,Z) coordinates --
*----------------------------------------------------------------------

*  -- Define (R,Z) coordinates at velocity cell corners.

      do 172 iy = 0, ny+1
         do 171 ix = 0, nx+1
            rv(ix,iy) = .5*(rm(ix+nj,iy,3)+rm(ix+nj,iy,4))
            zv(ix,iy) = .5*(zm(ix+nj,iy,3)+zm(ix+nj,iy,4))
 171     continue
         rv(nx+2,iy) = rm(nxm+1,iy,4)+.5*(rm(nxm+1,iy,4)-rm(nxm+1,iy,3))
         zv(nx+2,iy) = zm(nxm+1,iy,4)+.5*(zm(nxm+1,iy,4)-zm(nxm+1,iy,3))
 172  continue
      do 173 ix = 0, nx+1
         rv(ix,-1) = .5*(rm(ix+nj,0,1)+rm(ix+nj,0,2))
         zv(ix,-1) = .5*(zm(ix+nj,0,1)+zm(ix+nj,0,2))
 173  continue
      rv(nx+2,-1) = rm(nxm+1,0,2) + .5*(rm(nxm+1,0,2)-rm(nxm+1,0,1))
      zv(nx+2,-1) = zm(nxm+1,0,2) + .5*(zm(nxm+1,0,2)-zm(nxm+1,0,1))

*  -- Define (R,Z) coordinates for density cell faces and vel. cell centers.
      do 175 iy = 0, ny+1
         rfn(-1,iy) = .5*(rm(nj,iy,1)+rm(nj,iy,3))
         zfn(-1,iy) = .5*(zm(nj,iy,1)+zm(nj,iy,3))
         do 174 ix = 0, nx+1
            rfn(ix,iy) = .5*(rm(ix+nj,iy,2)+rm(ix+nj,iy,4))
            zfn(ix,iy) = .5*(zm(ix+nj,iy,2)+zm(ix+nj,iy,4))
            rcv(ix,iy) = rfn(ix,iy)
            zcv(ix,iy) = zfn(ix,iy)
 174     continue
 175  continue

*  -- Define (R,Z) coordinates for vel. cell faces and density cell centers.
      do 177 iy = 0, ny+1
         do 176 ix = 0, nx+1
            rcn(ix,iy) = .5*(rv(ix,iy)+rv(ix,iy-1))
            zcn(ix,iy) = .5*(zv(ix,iy)+zv(ix,iy-1))
            rfv(ix,iy) = rcn(ix,iy)
            zfv(ix,iy) = zcn(ix,iy)
 176     continue
         rfv(nx+2,iy) = .5*(rv(nx+2,iy)+rv(nx+2,iy-1))
         zfv(nx+2,iy) = .5*(zv(nx+2,iy)+zv(nx+2,iy-1))
 177  continue

*  -- Calculate angle the mesh lines make at the upper right vertex
*  -- of a cell, vtag; in radians (0 is orthogonal); angfx is angle on x-face
         call s2fill (nx+2, ny+2, 0., vtag, 1, nx+2)
         call s2fill (nx+2, ny+2, 0., angfx, 1, nx+2)
         call s2fill (nx+2, ny+2, 0., cosangfx, 1, nx+2)
         do iu = 0, 1
            call s2fill (nx+2, ny+2, 1., fx0(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxm(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxp(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 1., fy0(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fym(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fyp(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxmy(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxpy(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fymx(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fypx(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 1., fx0v(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxmv(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxpv(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 1., fy0v(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fymv(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fypv(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxmyv(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxpyv(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fymxv(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fypxv(0:nx+1,0:ny+1,iu), 1, nx+2)
         enddo
      if (isnonog .ge. 1) then
         call nonorthg   # Sets dist btwn interp pts as dxnog, dynog
                         # Also sets nonorth stencils:fx0,fxm,fy0,fym etc
         do ix = 0, nx+1    # Bdry value set: matters little except vygtan
ccc            gxfn(ix,0) = gxfn(ix,1)
            dxnog(ix,0) = dxnog(ix,1)
	    dxnog(ix,ny+1) = dxnog(ix,ny)
         enddo
         do iy = 0, ny+1  # likewise, dxnog(nx+1,) not relevant, but avoid 0
           dxnog(nx+1,iy) = 0.1*dxnog(nx,iy)
         enddo
      endif

*----------------------------------------------------------------------
*  -- Define density cell data --
*  -- define gx, gy, vol, rr, btot for density cells --
*  -- define sx at density faces in x, density centers in y --
*  -- define sy at density faces in y, density centers in x --
*----------------------------------------------------------------------
*  -- Define density cell data --
      do 302 iy = 0, ny+1
         do 301 ix = 0, nx+1

*  -- define gx, gy, vol, rr, btot for density cells --
            dx(ix,iy) = sqrt(
     .           (0.5*(rm(ix+nj,iy,4)+rm(ix+nj,iy,2)-rm(ix+nj,iy,3)
     .                -rm(ix+nj,iy,1)))**2
     .           + (0.5*(zm(ix+nj,iy,4)+zm(ix+nj,iy,2)-zm(ix+nj,iy,3)
     .                  -zm(ix+nj,iy,1)))**2)
            dy(ix,iy) = sqrt(
     .           (0.5*(rm(ix+nj,iy,3)+rm(ix+nj,iy,4)-rm(ix+nj,iy,1)
     .                -rm(ix+nj,iy,2)))**2
     .           + (0.5*(zm(ix+nj,iy,3)+zm(ix+nj,iy,4)-zm(ix+nj,iy,1)
     .                  -zm(ix+nj,iy,2)))**2)
cc      limit dx to avoid large Jacobian terms; test impact by changing
             dx(ix,iy) = max(dx(ix,iy), dxmin)
cc      for nonorthogonal mesh, dy gets reduced by cosine(average_angle)
            if ((islimon .ne. 0) .and. (ix .eq. ix_lim+1)) then
               dy(ix,iy) = dy(ix,iy)*cosangfx(ix,iy) 
            elseif (nxpt==2 .and. ix==ixlb(nxpt)) then  # full double null
c              Effectively, use angfx(ix-1,iy)=angfx(ix,iy) at left boundaries
               dy(ix,iy) = dy(ix,iy)*cosangfx(ix,iy) 
            else
               dy(ix,iy) = dy(ix,iy)*cos( 0.5*
     .                          (angfx(ix,iy) + angfx(ixm1(ix,iy),iy)) )
            endif
            dz = 2 * pi * 0.25 * (rm(ix+nj,iy,1)+rm(ix+nj,iy,2)
     .                           +rm(ix+nj,iy,3)+rm(ix+nj,iy,4))
            if (mhdgeo == -1) dz = 1.
            gx(ix,iy) = 1/dx(ix,iy)
            gy(ix,iy) = 1/dy(ix,iy)
            vol(ix,iy) = dx(ix,iy)*dy(ix,iy)*dz
            rr(ix,iy) = 0.25 * rr_fac *
     .           (bpol(ix+nj,iy,1)/b(ix+nj,iy,1)+bpol(ix+nj,iy,2)
     .                                                   /b(ix+nj,iy,2)
     .           +bpol(ix+nj,iy,3)/b(ix+nj,iy,3)+bpol(ix+nj,iy,4)
     .                                                   /b(ix+nj,iy,4))
            rbpol(ix,iy) = 0.25 *
     .           ((rm(ix+nj,iy,1)+r0slab)*bpol(ix+nj,iy,1)
     .           +(rm(ix+nj,iy,2)+r0slab)*bpol(ix+nj,iy,2)
     .           +(rm(ix+nj,iy,3)+r0slab)*bpol(ix+nj,iy,3)
     .           +(rm(ix+nj,iy,4)+r0slab)*bpol(ix+nj,iy,4))
            btot(ix,iy) = 0.25 * (b(ix+nj,iy,1)+b(ix+nj,iy,2)
     .                           +b(ix+nj,iy,3)+b(ix+nj,iy,4))
            b32(ix,iy) = ( 0.5*(b(ix,iy,2) + b(ix,iy,4)) )**1.5
            b12(ix,iy) = ( 0.5*(b(ix,iy,2) + b(ix,iy,4)) )**0.5
            b12ctr(ix,iy) = b(ix,iy,0)**0.5
            bsqr(ix,iy) = b(ix,iy,0)**2
            rbfbt(ix,iy) = 0.25 * ( 
     .              bphi(ix+nj,iy,1) + bphi(ix+nj,iy,2) +
     .              bphi(ix+nj,iy,3) + bphi(ix+nj,iy,4) ) / btot(ix,iy)
            rbfbt2(ix,iy) = rbfbt(ix,iy)/btot(ix,iy)
*     -- define sx at density faces in x, density centers in y --
*     -- 1/gyc is length*cos(angfx) between vertex (ix,iy) and (ix,iy-1)
            dyc = sqrt((rm(ix+nj,iy,4)-rm(ix+nj,iy,2))**2
     .              + (zm(ix+nj,iy,4)-zm(ix+nj,iy,2))**2) *
     .               cosangfx(ix,iy)
            gyc(ix,iy) = 1/dyc
            dz = 2 * pi * 0.5*(rm(ix+nj,iy,4)+rm(ix+nj,iy,2))
            if (mhdgeo == -1) dz = 1.
            sx(ix,iy) = dyc * dz
	    sxnp(ix,iy) = sx(ix,iy)/(cosangfx(ix,iy))
ctdr .. Special case - xgbx is used as reduction factor - not related to G.B.
ccc            if (ix.eq.ixpt1(1) .and. iy.le.iysptrx1(1)) then
ccc              sx(ix,iy) = xgbx*sx(ix,iy)
ccc            endif

*  -- define sy at density faces in y, density centers in x --
*  -- 1/gxc is the length between vertex (ix,iy) and (ix-1,iy)
            dxc = sqrt((rm(ix+nj,iy,4)-rm(ix+nj,iy,3))**2
     .              + (zm(ix+nj,iy,4)-zm(ix+nj,iy,3))**2)
            gxc(ix,iy) = 1/dxc
            dz = 2 * pi * 0.5*(rm(ix+nj,iy,4)+rm(ix+nj,iy,3))
            if (mhdgeo == -1) dz = 1.
            sy(ix,iy) = dxc * dz

 301     continue
 302  continue

c...  Compute area_core
        area_core = 0.
        do jx = 1, nxpt
          do ix = ixpt1(jx)+1, ixpt2(jx)
            area_core = area_core + sy(ix,0)
          enddo				  
        enddo

c...  Readjust rr in corner cells near x-point for the guard cells if
c...  this case is for core only (nyomitmx.ne.0)
      if (nyomitmx >= nysol+nyout) then
         rr(0,ny+1) = 0.5*(rr(0,ny)+rr(1,ny+1))
         rr(nx+1,ny+1) = 0.5*(rr(nx+1,ny)+rr(nx,ny+1))
      endif

c...  Define horizontal and vertical distance between velocity centers that
c...  straddle the x-point; 1/ghxpt is the horizontal distance and 1/gvxpt is 
c...  the vertical distance; sxyxpt is the average surface area
c
      if (geometry=="dnull" .or. geometry(1:9)=="snowflake" .or.
     .    geometry=="dnXtarget" .or. geometry=="isoleg") then
c         at lower x-point:
          ghxpt_lower = 2. / sqrt( 
     .         (rm(ixpt2(2)+nj,iysptrx2(2)+1,4) - rm(ixpt1(1)+nj,iysptrx1(1)+1,4))**2
     .       + (zm(ixpt2(2)+nj,iysptrx2(2)+1,4) - zm(ixpt1(1)+nj,iysptrx1(1)+1,4))**2 )
          gvxpt_lower = 2. / sqrt(     
     .         (rm(ixpt2(2)+nj,iysptrx2(2),2) - rm(ixpt1(1)+nj,iysptrx1(1),2))**2
     .       + (zm(ixpt2(2)+nj,iysptrx2(2),2) - zm(ixpt1(1)+nj,iysptrx1(1),2))**2 )
          sxyxpt_lower = 0.25 * (
     .          sy(ixpt1(1)+nj,iysptrx1(1)) + sy(ixpt1(1)+1+nj,iysptrx1(1)) +
     .          sy(ixpt2(2)+nj,iysptrx2(2)) + sy(ixpt2(2)+1+nj,iysptrx2(2)) ) 
c         at upper x-point:
          ghxpt_upper = 2. / sqrt( 
     .         (rm(ixpt2(1)+nj,iysptrx2(1)+1,4) - rm(ixpt1(2)+nj,iysptrx1(2)+1,4))**2
     .       + (zm(ixpt2(1)+nj,iysptrx2(1)+1,4) - zm(ixpt1(2)+nj,iysptrx1(2)+1,4))**2 )
          gvxpt_upper = 2. / sqrt(     
     .         (rm(ixpt2(1)+nj,iysptrx2(1),2) - rm(ixpt1(2)+nj,iysptrx1(2),2))**2
     .       + (zm(ixpt2(1)+nj,iysptrx2(1),2) - zm(ixpt1(2)+nj,iysptrx1(2),2))**2 )
          sxyxpt_upper = 0.25 * (
     .          sy(ixpt1(2)+nj,iysptrx1(2)) + sy(ixpt1(2)+1+nj,iysptrx1(2)) +
     .          sy(ixpt2(1)+nj,iysptrx2(1)) + sy(ixpt2(1)+1+nj,iysptrx2(1)) ) 
c
      else  # there is only one x-point in the simulation domain
        if (nyomitmx < nysol+nyout) then
        if (ixpt1(1)+nj.ge.0 .and. ixpt2(1)+nj.ge.0 .and. ixpt2(1)+nj.le.nx) then
          ghxpt = 2. / sqrt( 
     .         (rm(ixpt2(1)+nj,iysptrx2(1)+1,4) - rm(ixpt1(1)+nj,iysptrx1(1)+1,4))**2
     .       + (zm(ixpt2(1)+nj,iysptrx2(1)+1,4) - zm(ixpt1(1)+nj,iysptrx1(1)+1,4))**2 )
          gvxpt = 2. / sqrt(     
     .         (rm(ixpt2(1)+nj,iysptrx2(1),2) - rm(ixpt1(1)+nj,iysptrx1(1),2))**2
     .       + (zm(ixpt2(1)+nj,iysptrx2(1),2) - zm(ixpt1(1)+nj,iysptrx1(1),2))**2 )
          sxyxpt = 0.25 * (
     .          sy(ixpt1(1)+nj,iysptrx1(1)) + sy(ixpt1(1)+1+nj,iysptrx1(1)) +
     .          sy(ixpt2(1)+nj,iysptrx2(1)) + sy(ixpt2(1)+1+nj,iysptrx2(1)) )
        endif  # end if-test on ixpt1 and ixpt2
        endif  # end if-test on nyomitmx
      endif  # end if-test on geometry

*----------------------------------------------------------------------
*  -- Define velocity cell data --
*  -- define hxv, volv, rrv for velocity cells --
*  -- define syv for velocity cells --
*----------------------------------------------------------------------
      do 402 iy = 0, ny+1
         do 401 ix = 0, nx
            ix3 = ixp1(ix,iy)
            volv(ix,iy) = 0.5*(vol(ix3,iy) + vol(ix,iy))
            if (isrrvave .eq. 0) then
               rrv(ix,iy) = 0.5 * rr_fac *
     .                            (bpol(ix+nj,iy,2)/b(ix+nj,iy,2)+
     .                             bpol(ix+nj,iy,4)/b(ix+nj,iy,4))
            elseif (isrrvave .eq. 1) then  # avoids very small rrv at x-pt
               rrv(ix,iy) = 0.5*(rr(ix3,iy) + rr(ix,iy))
            elseif (isrrvave .eq. 2) then  # most accurate linear interp.
               rrv(ix,iy) = 0.25 * rr_fac *
     .                             (bpol(ix+nj,iy,2)/b(ix+nj,iy,2)+
     .                              bpol(ix+nj,iy,4)/b(ix+nj,iy,4)) +
     .                              0.25 * (rr(ix3,iy) + rr(ix,iy))
            endif
            syv(ix,iy) = 0.5*(sy(ix3,iy) + sy(ix,iy))
            hxv(ix,iy) = 
     &              2.0*(gx(ix3,iy)*gx(ix,iy)/(gx(ix3,iy)+gx(ix,iy)))
            dxvf(ix,iy) = dx(ix,iy)   # below ix=nx and nx+1 redefined
 401     continue
c        special cases --
            ix = nx+1
            ix3 = nx+2
            volv(ix,iy) = sx(ix,iy)*dx(ix,iy)
            rrv(ix,iy) = 0.5 * rr_fac *
     .                            (bpol(ix+nj,iy,2)/b(ix+nj,iy,2)+
     .                             bpol(ix+nj,iy,4)/b(ix+nj,iy,4))
            dz = 2 * pi * rm(ix+nj,iy,4)
            if (mhdgeo == -1) dz = 1.
            syv(ix,iy) = dx(ix,iy) * dz
            hxv(ix,iy) = 1/dx(ix,iy)
            dxvf(nx,iy) = 0.75*dx(nx,iy)
            dxvf(nx+1,iy) = 0.25*dx(nx,iy)
 402     continue

c ... radial velocity flow at x-pt is ambiguous; zero out surface area if..
      if (issyvxpt0==1) then
        do jx = 1, nxpt
        if (iysptrx1(jx).ge.0) then
          if (ixpt1(jx).gt.0 .and. ixpt1(jx).le.nx)
     .                              syv(ixpt1(jx),iysptrx1(jx)) = 0.
        endif
        if (iysptrx2(jx).ge.0) then
          if (ixpt2(jx).gt.0 .and. ixpt2(jx).le.nx)
     .                              syv(ixpt2(jx),iysptrx2(jx)) = 0.
        endif
        enddo
      endif

*---------------------------------------------------------------------
*     -- Define inverse of distance btwn density cell ctrs
*---------------------------------------------------------------------
      do 421 iy = 0,  ny+1
         do 420 ix = 0, nx
            ix1 = ixp1(ix,iy)
            gxf(ix,iy) = 2*gx(ix,iy)*gx(ix1,iy) / (gx(ix,iy)+gx(ix1,iy))
 420     continue
 421  continue

c...  Compute gyf for both isnonog=0,1; nonog mesh changes gy & gyf
      do 424 iy = 0, ny
         do 423 ix = 0, nx
            ix1 = ixp1(ix,iy)
            gyf(ix,iy) = 2*gy(ix,iy)*gy(ix,iy+1)/(gy(ix,iy)+gy(ix,iy+1))
 423     continue
         gyf(nx+1,iy) = 2*gy(nx+1,iy)*gy(nx+1,iy+1)/
     .                	 (gy(nx+1,iy)+gy(nx+1,iy+1))
 424  continue

c...  Set dxnog for orthog mesh; if isnonog=1, call to nonorthg above sets
      if (isnonog==0) then
        do iy = 0, ny+1
          do ix = 0, nx+1
             dxnog(ix,iy) = 1./(gxf(ix,iy) + 1.e-100)
             dynog(ix,iy) = 1./(gyf(ix,iy) + 1.e-100)
          enddo
        enddo
      endif

c...  Calculate a normalization constant for the iy=0 cells
c...  of the core boundary region
      if ((geometry(1:9)=="snowflake" .and. geometry(10:11).ne."15")
     &       .or. geometry=="dnXtarget") then
        ix_last_core_cell = ixpt2(1)
      else
        ix_last_core_cell = ixpt2(nxpt)
      endif
      ix = max(1, ixpt1(1)+1)
      ix = min(ix, nx+1)
      sygytotc = sy(ix,0)*gyf(ix,0)
      do ix = ixpt1(1)+1, ix_last_core_cell
         if ( .not. (isudsym==1 .and. ((ix==nxc) .or. (ix==nxc+1))) ) then
           sygytotc = sygytotc + sy(ix,0)*gyf(ix,0)
         endif
      enddo

c...  Setup the isixcore(ix) array: =1 if ix on iy=0 core bdry; =0 if not
      do ix = 0, nx+1
        if ((geometry(1:9)=="snowflake" .and. geometry(10:11).ne."15")
     &       .or. geometry=="dnXtarget") then
          if( ix > ixpt1(1) .and. ix <= ixpt2(1)) then
            isixcore(ix) = 1
          else
            isixcore(ix) = 0
          endif
        else  # all other geometries
          if(nxpt == 1) then
            isixcore(ix) = 0
            if(ix > ixpt1(1) .and. ix <= ixpt2(1)) isixcore(ix) = 1
          elseif(nxpt == 2) then
            isixcore(ix) = 0
            if( (ix > ixpt1(1) .and. ix <= ixpt2(1)).or.
     .        (ix > ixpt1(2) .and. ix <= ixpt2(2)) ) isixcore(ix) = 1
          endif
        endif
      enddo
          
*----------------------------------------------------------------------
*  -- Calculate geometrical factors needed for curvature and grad_B drifts

      if (mhdgeo .ge. 0) then   #skip the slab case, or add r0slab
      do iy = 0, ny+1
        do ix = 0, nx+1
          s_bphi = sign(1.,bphi(1,1,1))
          cossr = -s_bphi*(rm(ix+nj,iy,4)-rm(ix+nj,iy,3)) / (
     .                    (rm(ix+nj,iy,4)-rm(ix+nj,iy,3))**2 + 
     .                    (zm(ix+nj,iy,4)-zm(ix+nj,iy,3))**2 )**0.5
          curvrby(ix,iy) = -2*cossr/( rm(ix+nj,iy,4)*b(ix+nj,iy,4) +
     .                               rm(ix+nj,iy,3)*b(ix+nj,iy,3) )
          cossp = (rm(ix+nj,iy,4)-rm(ix+nj,iy,2)) / (
     .            (rm(ix+nj,iy,4)-rm(ix+nj,iy,2))**2 + 
     .            (zm(ix+nj,iy,4)-zm(ix+nj,iy,2))**2 )**0.5
          curvrb2(ix,iy) = -2*cossp/( rm(ix+nj,iy,4)*b(ix+nj,iy,4) +
     .                               rm(ix+nj,iy,2)*b(ix+nj,iy,2) )
          gradby(ix,iy) = -0.5*( bphi(ix+nj,iy,4)/b(ix+nj,iy,4)**3 +
     .                           bphi(ix+nj,iy,3)/b(ix+nj,iy,3)**3 )*
     .                    (b(ix+nj,iy,4)-b(ix+nj,iy,3))*gxc(ix,iy)
          gradb2(ix,iy) = 0.5*( 1/b(ix+nj,iy,4)**2 +
     .                          1/b(ix+nj,iy,2)**2 )*
     .                    (b(ix+nj,iy,4)-b(ix+nj,iy,2))*gyc(ix,iy)*
     .                        cosangfx(ix,iy) # cos corrects gyc
ccc          if (iy .ge. ny-1) then
ccc            curvrby(ix,iy) = 0.
ccc            gradby(ix,iy) = 0.
ccc          endif
        enddo
      enddo
      endif  #end test for mhdgeo

*  -- Calc Rozhansky factor [1-B**2/(B**2_ave)] for possible use
      bsqrvoltot = 0.
      voltot = 0.
      do iy = 1, iysptrx1(1)
        do ix = ixpt1(1)+1, ixpt2(1)
          bsqrvoltot = bsqrvoltot + btot(ix,iy)**2*vol(ix,iy)*dx(ix,iy)
          voltot = voltot + vol(ix,iy)*dx(ix,iy)
        enddo
      enddo
      if (voltot > 0.) bsqrave = bsqrvoltot/voltot

      do iy = 0, ny+1
        do ix = 0, nx+1
          b_yface = 0.5*(b(ix,iy,3) + b(ix,iy,4))
          b_xface = 0.5*(b(ix,iy,2) + b(ix,iy,4))
          bfacyrozh(ix,iy) = 1.
          bfacxrozh(ix,iy) = 1.
          if (isrozhfac == 1) then
            bfacyrozh(ix,iy) = (1. - b_yface**2/bsqrave)
            bfacxrozh(ix,iy) = (1. - b_xface**2/bsqrave)
          endif
        enddo
      enddo

*----------------------------------------------------------------
*  -- define some magnetic field arrays used for neoclassical calc
      do iy = 0, ny+1
        iyp1 = min(iy+1, ny+1)
        do ix = 0, nx+1
          dbm2dy(ix,iy) = ( 1./btot(ix,iyp1)**2 -
     .                      1./btot(ix,iy)**2 )*gyf(ix,iy)
          dbm2dx(ix,iy) = ( 1./btot(ixp1(ix,iy),iy)**2 -
     .                      1./btot(ix,iy)**2 )*gxf(ix,iy)
        enddo
      enddo

*----------------------------------------------------------------------
*  -- define (xcs,ycn) coordinate values --
*     x is measured poloidally from the inboard plate
*     y is measured radially outward from the separatrix

*  -- define x on density faces -- along the separatrix, inner, & outer wall
c     MER NOTE:
c     The general case has multiple separatrices so xfs,xcs,xfv,xcv
c     are not uniquely defined; below we choose the outermost separatrix.
      iyso = max(iysptrx1(1),iysptrx2(1))
      xfs(0) = 1./gx(0,iyso+1)
      xfwi(0) = 1./gx(0,0)
      xfwo(0) = 1./gx(0,ny)
      do ix = 1, nx+1
         xfwo(ix) = xfwo(ix-1) + 1./gx(ix,ny)
         xfs(ix) = xfs(ix-1) + 1./gx(ix,iyso+1)
         xfwi(ix) = xfwi(ix-1) + 1./gx(ix,0)
      enddo
*  -- define x on private flux walls --
c ... The p.f. walls may extend over adjoining mesh regions for a double
c ... null configuration; to ensure continuity we define the two halves
c ... sequentially. (MER 06 Nov 1999)
c ... First, between the left boundary and x-point cut:
      do jx = 1, nxpt
        if(ixpt1(jx) >= 0) then
          xfpf(ixlb(jx)) = 1./gx(ixlb(jx),0)
          do ix = ixlb(jx)+1, ixpt1(jx)
            xfpf(ix) = xfpf(ix-1) + 1./gx(ix,0)
          enddo
	  xfpf(ixpt1(jx)+1) = 1./gx(ixpt1(jx)+1,0)
	  do ix = ixpt1(jx)+2, ixpt2(jx)
	    xfpf(ix) = xfpf(ix-1) + 1./gx(ix,0)
          enddo
        endif
      enddo
c ... Next do the core boundary; later xfpf=0 for core - just convention
      do jx = 1, nxpt
        if(ixpt1(jx) >= 0) then
	  if (jx == 1) then
	    xfpf(ixpt1(jx)+1) = 1./gx(ixpt1(jx)+1,0)
          else # jx=2 for second (outer) core region if nxpt=2 only
	    xfpf(ixpt1(jx)+1) = xfpf(ixpt2(1)) + 1./gx(ixpt1(jx)+1,0)
          endif
	  do ix = ixpt1(jx)+2, ixpt2(jx)
	    xfpf(ix) = xfpf(ix-1) + 1./gx(ix,0)
          enddo
        endif
      enddo
c ... Then, between the x-point cut and right boundary:
      do jx = 1, nxpt
         do ix = max(1, ixpt2(jx)+1), ixrb(jx)+1
            xfpf(ix) = xfpf(ixm1(ix,0)) + 1./gx(ix,0)
         enddo
      enddo

*  -- define x on density centers -- midway between density faces
      do jx=1, nxpt
         if (jx==1) then
            xcs(ixlb(jx))  = 0.5*xfs(ixlb(jx))
            xcwo(ixlb(jx)) = 0.5*xfwo(ixlb(jx))
         else
	    xcs(ixlb(jx))  = 0.5 * ( xfs(ixlb(jx)-1) + xfs(ixlb(jx)) )
	    xcwo(ixlb(jx)) = 0.5 * ( xfwo(ixlb(jx)-1) + xfwo(ixlb(jx)) )
         endif
         xcpf(ixlb(jx)) = 0.5*xfpf(ixlb(jx))

         do ix = ixlb(jx)+1, ixrb(jx)+1
            ix1 = ixm1(ix,0)
            ix2 = ixm1(ix,ny)
            xcs(ix) = 0.5*(xfs(ix-1)+xfs(ix))
            xcwo(ix) = 0.5 * ( xfwo(ix2) + xfwo(ix) )
            xcpf(ix) = 0.5 * ( xfpf(ix1) + xfpf(ix) )
            xcwi(ix) = xcpf(ix)
         enddo
      enddo
      do jx = 1, nxpt
        if(ixpt1(jx) >= 0) then
          xcwi(ixpt1(jx)+1) = 0.5*xfpf(ixpt1(jx)+1)  # first cell in core
        endif
      enddo

c ... Fix the core boundary; just a convention
      do jx = 1, nxpt
	do ix  = max(0, ixpt1(jx)+1), ixpt2(jx)
	   xfpf(ix) = 0.           # just a convention
	   xcpf(ix) = 0.
        enddo
      enddo

*  -- define x on velocity faces -- coincident with density centers
      do 22 ix = 0, nx+1
         xfv(ix) = xcs(ix)
 22   continue

*  -- define x on velocity centers -- midway between velocity faces
      do 23 ix = 0, nx
         xcv(ix) = 0.5*(xfv(ix)+xfv(ix+1))
 23   continue
*     xcv(nx+1) should not matter
      xcv(nx+1) = xfv(nx+1)

*  -- define y on density faces -- at outboard midplane (?)
      if (ixpt2(1) > 0 .and. (isudsym.ne.1)) then
         rmmax = rm(nxleg(1)+nxcore(1)+1,ny,0)
         do ix = nxleg(1)+nxcore(1)+1, ixpt2(1)
           if (rm(ix,ny,0) >= rmmax) then
              rmmax = rm(ix,ny,0)
              ixmp = ix
           endif
         enddo
      endif
      if (geometry.eq.'dnbot') ixmp = nxc+1
      if (geometry.eq.'dnull' .or. geometry=='snowflake15' .or.
     .    geometry=='snowflake45' .or. geometry=='snowflake75' .or.
     .    geometry=='dnXtarget' .or. geometry=='isoleg') ixmp = ixmdp(2)
 
*  -- redefine ixmp if it is outside ix=0,nx domain
      if(ixmp.lt.0 .or. ixmp.gt.nx) then  #search for max rm
	 rmmax = rm(nxomit,0,0)
         do ix = nxomit, nx
            if(rm(ix,0,0) > rmmax) then
               rmmax = rm(ix,0,0)
               ixmp = ix
            endif
         enddo
      endif

c     MER NOTE:
c     The general case has multiple separatrices so yyf=0 is not
c     uniquely defined; below we choose the innermost separatrix.
      iysi = min(iysptrx1(1),iysptrx2(1))
      yyf(iysi) = 0.
      psinormf(iysi) = 1.                # normalized pol. flux
      do 24 iy = iysi+1, ny+1
         yyf(iy) = yyf(iy-1) + 1./gy(ixmp,iy)
 24   continue
      do 241 iy = iysi-1, 0, -1
         yyf(iy) = yyf(iy+1) - 1./gy(ixmp,iy+1)
 241  continue

*  -- define y on density centers -- midway btwn density faces
      yyc(0) =yyf(0) - 0.5/gy(ixmp,0)
      do 25 iy = 1, ny+1
         yyc(iy) = 0.5*(yyf(iy-1)+yyf(iy))
 25   continue

c...  Calculate normalized poloidal flux for mhdgeo=1 cases (tor. equil)
      if(mhdgeo == 1 .and. gengrid == 1) then  # otherwise sibdrys etc.=0
	do iy = 0, ny+1
	  psinormf(iy) = ( 0.5*(psi(ixmp,iy,3)+psi(ixmp,iy,4)) -
     .                             simagxs ) / (sibdrys - simagxs)
	  psinormc(iy) = (psi(ixmp,iy,0)-simagxs) / (sibdrys-simagxs)
        enddo
      endif

c --- compute normalized mesh for interpolation of solution
c --- Here the initial xnrm, ynrm, etc. are calculated
         do jx=1,nxpt
            call grdnrm (nx,ny,ixlb(jx),ixpt1(jx),ixpt2(jx),ixrb(jx),
     .                             iysptrx,isgindx,gx,gxf,gy,gyf,gyc,
     .                                         xnrm,xvnrm,ynrm,yvnrm)
         enddo

      if (isnonog .ge. 1) then
c...  Reset fxm, fx0, fxp around the x-point - only orthogonal coupling
c...  Also reset gyf, so we need the gy calc. after call to nonorthg
       if (isfixlb(1) .eq. 0 .and. isfixrb(1) .eq. 0) then 
        do ik = 0, 1
          do ij = 0, 1
            do jx = 1, nxpt
              fxm(ixpt1(jx)+ij,iysptrx1(jx),ik)  = 0.
              fx0(ixpt1(jx)+ij,iysptrx1(jx),ik)  = 1.
              fxp(ixpt1(jx)+ij,iysptrx1(jx),ik)  = 0.
              fxmy(ixpt1(jx)+ij,iysptrx1(jx),ik) = 0.
              fxpy(ixpt1(jx)+ij,iysptrx1(jx),ik) = 0.
              fxm(ixpt2(jx)+ij,iysptrx2(jx),ik)  = 0.
              fx0(ixpt2(jx)+ij,iysptrx2(jx),ik)  = 1.
              fxp(ixpt2(jx)+ij,iysptrx2(jx),ik)  = 0.
              fxmy(ixpt2(jx)+ij,iysptrx2(jx),ik) = 0.
              fxpy(ixpt2(jx)+ij,iysptrx2(jx),ik) = 0.
              gyf(ixpt1(jx)+ij,iysptrx1(jx)) = 
     .                                2*gy(ixpt1(jx)+ij,iysptrx1(jx)) *
     .                              gy(ixpt1(jx)+ij,iysptrx1(jx)+1) / (
     .                                 gy(ixpt1(jx)+ij,iysptrx1(jx)) + 
     .                                 gy(ixpt1(jx)+ij,iysptrx1(jx)+1) ) 
              gyf(ixpt2(jx)+ij,iysptrx2(jx)) = 
     .                                2*gy(ixpt2(jx)+ij,iysptrx2(jx)) *
     .                              gy(ixpt2(jx)+ij,iysptrx2(jx)+1) / (
     .                                 gy(ixpt2(jx)+ij,iysptrx2(jx)) + 
     .                                 gy(ixpt2(jx)+ij,iysptrx2(jx)+1) )
              dynog(ixpt1(jx)+ij,iysptrx1(jx)) = 
     .                                 1./gyf(ixpt1(jx)+ij,iysptrx1(jx))
              dynog(ixpt2(jx)+ij,iysptrx2(jx)) = 
     .                                 1./gyf(ixpt2(jx)+ij,iysptrx2(jx))
            enddo # end loop on jx 
          enddo # end loop on ij
        enddo # end loop on ik

c...  Reset fym, fy0, fyp around the x-point - only orthogonal coupling
        do ik = 0, 1
          do ij = 0, 1  #was -1,1; changed 12/13/19
            do jx = 1, nxpt
              fym(ixpt1(jx),iysptrx1(jx)+ij,ik)  = 0.
              fy0(ixpt1(jx),iysptrx1(jx)+ij,ik)  = 1.
              fyp(ixpt1(jx),iysptrx1(jx)+ij,ik)  = 0.
              fypx(ixpt1(jx),iysptrx1(jx)+ij,ik) = 0.
              fymx(ixpt1(jx),iysptrx1(jx)+ij,ik) = 0.
              fym(ixpt2(jx),iysptrx2(jx)+ij,ik)  = 0.
              fy0(ixpt2(jx),iysptrx2(jx)+ij,ik)  = 1.
              fyp(ixpt2(jx),iysptrx2(jx)+ij,ik)  = 0.
              fypx(ixpt2(jx),iysptrx2(jx)+ij,ik) = 0.
              fymx(ixpt2(jx),iysptrx2(jx)+ij,ik) = 0.
              dxnog(ixpt1(jx),iysptrx2(jx)+ij) = 
     .                               1./gxf(ixpt1(jx),iysptrx1(jx)+ij)
              dxnog(ixpt2(jx),iysptrx2(jx)+ij) = 
     .                               1./gxf(ixpt2(jx),iysptrx2(jx)+ij)
            enddo # end loop on jx
          enddo # end loop on ij
        enddo # end loop on ik
       endif

c...  Reset fxm, fx0, fxp at ixpt1,2 if half-space problem with no flux
        if (isfixlb(1).eq.2 .or. isfixrb(1).eq.2) then
           ix = ixpt2(1)
           if (isfixrb(1).eq.2) ix = ixpt1(1)
           do ik = 0, 1
             do iy = 0, iysptrx1(1)+1
               do ij = 0, 1
                  fxm(ix+ij,iy,ik)  = 0.
                  fx0(ix+ij,iy,ik)  = 1.
                  fxp(ix+ij,iy,ik)  = 0.
                  fxmy(ix+ij,iy,ik) = 0.
                  fxpy(ix+ij,iy,ik) = 0.
                  gyf(ix+ij,iy) = 2*gy(ix+ij,iy) * gy(ix+ij,iy+1) / (
     .                           gy(ix+ij,iy) + gy(ix+ij,iy+1) ) 
                enddo
             enddo
           enddo
        endif

c...  Reset fxm, fx0, fxp at iy=ny temporarily (6/13/96) to see if it fixes
c...  Porter problem with low density carbon
      if (isybdryog .eq. 1) then
      do iy = 0, ny, ny
        do ik = 0, 1
           do ix = 0, nx+1
              fxm(ix,iy,ik)  = 0.
              fx0(ix,iy,ik)  = 1.
              fxp(ix,iy,ik)  = 0.
              fxmy(ix,iy,ik) = 0.
              fxpy(ix,iy,ik) = 0.
              gyf(ix,iy) = 2*gy(ix,iy) * gy(ix,iy+1) / (
     .                           gy(ix,iy) + gy(ix,iy+1) ) 
           enddo
        enddo
      enddo
      endif

c...  Reset fym, fy0, fyp at ix=nxc-1 and ix=nxc+1 to impose zero-flux
c...  boundary conditions at midplane for geometry='dnbot'
      if (isudsym==1 .and. isxmpog==1) then
        do ix = nxc-1, nxc+1, 2
          do ik = 0, 1
            do iy = 0, ny+1
               fym(ix,iy,ik)  = 0.
               fy0(ix,iy,ik)  = 1.
               fyp(ix,iy,ik)  = 0.
               fymx(ix,iy,ik) = 0.
               fypx(ix,iy,ik) = 0.
               dxnog(ix,iy) = ( gx(ix,iy) + gx(ix+1,iy) )/
     .                        (2*gx(ix,iy) * gx(ix+1,iy))
            enddo
          enddo
        enddo
      endif

c...  Reset fym, fy0, fyp at ixpt1,2+0,1 if half-space problem with no flux
        if (isfixlb(1).eq.2 .or. isfixrb(1).eq.2) then
           ix = ixpt2(1)
           ix2 = 0
           if (isfixrb(1).eq.2) then
              ix = ixpt1(1)
              ix2 = nx
           endif
           do ik = 0, 1
              do iy = 0, iysptrx1(1)+1
                 fym(ix,iy,ik)  = 0.
                 fy0(ix,iy,ik)  = 1.
                 fyp(ix,iy,ik)  = 0.
                 fypx(ix,iy,ik) = 0.
                 fymx(ix,iy,ik) = 0.
                 dxnog(ix,iy) = 1./gxf(ix,iy)
              enddo
              do iy = 0, ny+1
                 fym(ix2,iy,ik)  = 0.
                 fy0(ix2,iy,ik)  = 1.
                 fyp(ix2,iy,ik)  = 0.
                 fypx(ix2,iy,ik) = 0.
                 fymx(ix2,iy,ik) = 0.
                 dxnog(ix2,iy) = 1./gxf(ix,iy)
              enddo
           enddo
        endif

      endif   # end of if (isnonog .ge. 1)


c...  As a test, if isnonog.ge.2, reset fym etc. to 5-point stencil
      if (isnonog .ge. 2) then
         do iu = 0, 1
            call s2fill (nx+2, ny+2, 1., fx0(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxm(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxp(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 1., fy0(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fym(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fyp(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxmy(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fxpy(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fymx(0:nx+1,0:ny+1,iu), 1, nx+2)
            call s2fill (nx+2, ny+2, 0., fypx(0:nx+1,0:ny+1,iu), 1, nx+2)
         enddo
      endif

c...  Set the nonorthogonal stencil for the velocity cells
      if (isnonog .ge. 1) then

        do iu = 0, 1
          do iy = 0, ny+1
            do ix = 0, nx+1
              ix1 = ixp1(ix,iy)
              ix2 = ixm1(ix,iy)
              fxmv (ix,iy,iu) = 0.5*(fxm (ix,iy,iu)+fxm (ix1,iy,iu))
              fx0v (ix,iy,iu) = 0.5*(fx0 (ix,iy,iu)+fx0 (ix1,iy,iu))
              fxpv (ix,iy,iu) = 0.5*(fxp (ix,iy,iu)+fxp (ix1,iy,iu))
              fxmyv(ix,iy,iu) = 0.5*(fxmy(ix,iy,iu)+fxmy(ix1,iy,iu))
              fxpyv(ix,iy,iu) = 0.5*(fxpy(ix,iy,iu)+fxpy(ix1,iy,iu))
              fymv (ix,iy,iu) = 0.5*(fym (ix,iy,iu)+fym (ix2,iy,iu))
              fy0v (ix,iy,iu) = 0.5*(fy0 (ix,iy,iu)+fy0 (ix2,iy,iu))
              fypv (ix,iy,iu) = 0.5*(fyp (ix,iy,iu)+fyp (ix2,iy,iu))
              fymxv(ix,iy,iu) = 0.5*(fymx(ix,iy,iu)+fymx(ix2,iy,iu))
              fypxv(ix,iy,iu) = 0.5*(fypx(ix,iy,iu)+fypx(ix2,iy,iu))
            enddo
          enddo
        enddo

c...  Special settings for velocity boundary cells at left and right edge
c...  of mesh:
ccccc** Shouldn't matter as up=cs BC, but also fx0v should =1
        do iu = 0, 1
          do iy = 0, ny+1
            do jx = 1, nxpt
              fxmv(ixlb(jx),iy,iu) = 0.
              fxmyv(ixlb(jx),iy,iu) = 0.
              fxpv(ixrb(jx),iy,iu) = 0.
              fxpyv(ixrb(jx),iy,iu) = 0.
            enddo
          enddo
        enddo

c   ... If a limiter is present, reset velocity stencil
        if (islimon .ne. 0) then 
          do iu = 0, 1
            do iy = 0, ny+1 # fyov... for iy>iy_lims not used, set anyway
               fy0v(ix_lim+iu,iy,1-iu) = 1.
               fymv(ix_lim+iu,iy,1-iu) = 0.
               fypv(ix_lim+iu,iy,1-iu) = 0.
               fymxv(ix_lim+iu,iy,1-iu) = 0.
               fypxv(ix_lim+iu,iy,1-iu) = 0.
             enddo
           enddo
        endif          

      endif   # end of if (isnonog .ge. 1)


c ... Find poloidal index for gas box in 1-D model
      if (is1D_gbx .eq. 1) then
         ixgb = 0
         do ix = 0, nx+1
            if (zm(ix,1,2) .lt. xgbx) ixgb = ix
         enddo
      endif 

c ... Set arrays that flag cells near x-points where one-sided derivatives
c ... may be required:
      do jx = 1, nxpt
        do ix = ixlb(jx), ixrb(jx)+1
          do iy = 0, ny+1
           isxptx(ix,iy)=1
           isxpty(ix,iy)=1
c ... isxptx=0 for cells whose x-face touches an x-point; else isxptx=1
            isxptx(ix,iy)=1
            if ( (ix==ixpt1(jx)) .and.
     .           (iy==iysptrx1(jx) .or. iy==iysptrx1(jx)+1) ) isxptx(ix,iy)=0
            if ( (ix==ixpt2(jx)) .and.
     .           (iy==iysptrx2(jx) .or. iy==iysptrx2(jx)+1) ) isxptx(ix,iy)=0
c ... isxpty=0 if y-face touches Xpt from below; =-1 if touches from above
c ... one-side y-deriv from below if isxpty=0 & from above if isxpty=-1
c ... May be needed for parallel domain decomp cases
            isxpty(ix,iy)=1
            if ( (iy==iysptrx1(jx)) .and. nyomitmx < nysol+nyout .and.
     .           (ix==ixpt1(jx) .or. ix==ixpt1(jx)+1) ) isxpty(ix,iy)=0
            if ( (iy==iysptrx2(jx)) .and. nyomitmx < nysol+nyout .and.
     .           (ix==ixpt2(jx) .or. ix==ixpt2(jx)+1) ) isxpty(ix,iy)=0
             if ( (iy==iysptrx1(jx)+1) .and. nyomitmx < nysol .and.
     .             ix==ixpt1(jx)+1 ) isxpty(ix,iy)=-1
             if ( (iy==iysptrx2(jx)+1) .and. nyomitmx < nysol .and.
     .             ix==ixpt2(jx)+1 ) isxpty(ix,iy)=-1
          enddo
        enddo
      enddo  # end do-loop over nxpt x-points
 
      return
      end
c ***** end of subroutine nphygeo ****
c ************************************
c-----------------------------------------------------------------------
      subroutine conlen
*     Calculate the ion and electron connection lengths for a poloidal
*     circuit or between end plates.  These lengths are interpolated 
*     within a banana width of the separatrix to account for finite
*     banana-width effects

      implicit none
      Use(Dim)            # nx,ny
      Use(Comgeo)         # lcon,lconi,lcone,rr,gx,yyc
      Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2,iysptrx
      Use(Aux)            # ixmp
      Use(RZ_grid_info)   # bpol
      Use(Phyvar)         # me,ev
      Use(Compla)         # mi
      Use(Comtra)         # tibsep,tebsep,cfelecbwd
      Use(UEint)          # mhdgeo

*  -- local scalars --
      integer iyso,iybwmni,iybwmne,iybwmxi,iybwmxe,ixref
      #Former Aux module variables
      integer ix,iy
      real lcon_wk1,lcon_wk2,omegcips,banwidi,banwide,dlcon,dyyc
      real eps_wk1,eps_wk2,rfac

ccc New computation of lconneo (2012) for single-null only
cc ------------------------------------------------------------------
      iyso = iysptrx1(1)
      do iy = 0, ny+1  #do core and SOL above x-point
        lcon_wk1 = 0.
        eps_wk2 = 0.
        eps_wk1 = 0.
        do ix = ixpt1(1)+1, ixpt2(1) #only use integral btwn x-points
          lcon_wk1 = lcon_wk1 + 1/(rr(ix,iy)*gx(ix,iy))
          eps_wk2 = eps_wk2 + 1/gx(ix,iy)
          eps_wk1 = eps_wk1 + 2.*pi*rm(ix,iy,0)/gx(ix,iy)
        enddo
        do ix = 0, nx+1
          lconneo(ix,iy) = lcon_wk1/(2.*pi)
          if (mhdgeo > 0) then
            epsneo(ix,iy) = eps_wk2**2/eps_wk1
          else  #cylinder or slab
            epsneo(ix,iy) = 1.e-100
          endif
        enddo
      enddo
ccc  Radially smooth, very ad hoc correction for PF region
      do iy = 0, iyso
        rfac = (float(iy+1)/float(iyso+1))**2
        do ix = 0, ixpt1(1)
          lconneo(ix,iy) = lconneo(ix,iyso+1)*rfac
          if (mhdgeo > 0) then
            epsneo(ix,iy) = epsneo(ix,iyso+1)*rfac
          else  #cylinder or slab
            epsneo(ix,iy) = 1.e-100
          endif
        enddo
        do ix = ixpt2(1)+1, nx+1
          lconneo(ix,iy) = lconneo(ix,iyso+1)*rfac
          if (mhdgeo > 0) then
            epsneo(ix,iy) = epsneo(ix,iyso+1)*rfac
          else  #cylinder or slab
            epsneo(ix,iy) = 1.e-100
          endif
        enddo
      enddo
cc ------------------------------------------------------------------

ccc  Computation of lconi & lcone for old (~2005) attempt at kinetics; 
      iyso = max(iysptrx1(1),iysptrx2(1))  # choose outer separatrix if 2
      do iy = 0, ny+1
         lcon_wk1 = 0.         # Place-holder work var for SOL or tot
         lcon_wk2 = 0.         # Place-holder work var for PF only
	 if(iy > iyso) then    # SOL, so only 1 region; lcon flux func
            do ix = 1, nx
	      lcon_wk1 = lcon_wk1 + 1/(rr(ix,iy)*gx(ix,iy))
            enddo
	    do ix = 0, nx+1
	      lcon(ix,iy) = lcon_wk1   #might limit by dist to wall
              lconi(ix,iy) = lcon(ix,iy)
              lcone(ix,iy) = lcon(ix,iy)
            enddo
	 else                  # Core and PF done together - carefully 
            do ix = 1, nx      # Do sum of core+SOL - then separate
              lcon_wk1 = lcon_wk1 + 1/(rr(ix,iy)*gx(ix,iy))
            enddo
            do ix = 1, ixpt1(1)  # Do PF only in next two wk2 loops
              lcon_wk2 = lcon_wk2 + 1/(rr(ix,iy)*gx(ix,iy))
            enddo
            do ix = max(1, ixpt2(1)+1), nx
              lcon_wk2 = lcon_wk2 + 1/(rr(ix,iy)*gx(ix,iy))
            enddo
            do ix = 0, nx+1    # Now separate into PF and core
              if (ix <= ixpt1(1) .or. ix > ixpt2(1)) then  # PF region
                lcon(ix,iy) = lcon_wk2
                lconi(ix,iy) = lcon(ix,iy)
                lcone(ix,iy) = lcon(ix,iy)
              else                                         # Core region
                lcon(ix,iy) = lcon_wk1 - lcon_wk2
                lconi(ix,iy) = lcon(ix,iy)
                lcone(ix,iy) = lcon(ix,iy)
              endif
            enddo
         endif
      enddo

c...  Correct ion/elec connect lengths near separatrix - banana widths
      omegcips = (qe*0.5*(bpol(ixmp,iysptrx,0) + bpol(ixmp,iysptrx+1,0)))
      banwidi = sqrt(2*tibsep*ev/mi(1))/omegcips
      banwide = cfelecbwd*sqrt(2*tebsep*ev/me)/((mi(1)/me)*omegcips)
      iybwmni = 0
      iybwmne = 0
      iybwmxi = 0
      iybwmxe = 0
      do iy = 0, ny+1   # Find i,e iy indices within banana width of sep
        if (yyc(iy) + 0.5*banwidi > 0. .and. iybwmni == 0) then
          iybwmni = iy
        elseif (yyc(iy) - 0.5*banwidi > 0. .and. iybwmxi == 0) then
          iybwmxi = iy-1
        endif
        if (yyc(iy) + 0.5*banwide > 0. .and. iybwmne == 0) then
          iybwmne = iy
        elseif (yyc(iy) - 0.5*banwide > 0. .and. iybwmxe == 0) then
          iybwmxe = iy-1
        endif
      enddo
      iybwmxi = max(iybwmni+1, iybwmxi)  # have at least one cell about sep
      iybwmxe = max(iybwmne+1, iybwmxe)

c...  First do ions
      do iy = 0, ny+1 
         dlcon = (lcon(ixmp,iybwmxi)-lcon(ixmp,iybwmni))
         dyyc = (yyc(iybwmxi)-yyc(iybwmni))
	 if(iy > iyso) then    # SOL, so only 1 region; lcon flux func
	    do ix = 0, nx+1  
              if (iy <= iybwmxi) then  # lin interp over banana width
                 lconi(ix,iy) = lcon(ixmp,iybwmni) + 
     .                            dlcon*(yyc(iy)-yyc(iybwmni))/dyyc
              else  
                 lconi(ix,iy) = lcon(ix,iy)
              endif
            enddo
	 else                  # PF and core
            do ix = 0, nx+1    # Now separate into PF and core
              if (ix <= ixpt1(1) .or. ix > ixpt2(1)) then  # PF region
                if(ix <= ixpt1(1)) then
                  ixref = ixpt1(1)
                else # ix >= ixpt2(1)
                  ixref = ixpt2(1)+1
                endif
                if (iy >= iybwmni) then  # lin interp over banana width
                     # PF region is approx, could be improved
                  lconi(ix,iy) = lcon(ixref,iybwmni) + 
     .                            dlcon*(yyc(iy)-yyc(iybwmni))/dyyc
                else
                  lconi(ix,iy) = lcon(ix,iy)
                endif
              else                                         # Core region
                if (iy >= iybwmni) then  # lin interp over banana width
                  lconi(ix,iy) = lcon(ixmp,iybwmni) + 
     .                            dlcon*(yyc(iy)-yyc(iybwmni))/dyyc
                else
                  lconi(ix,iy) = lcon(ix,iy)
                endif
              endif
            enddo
         endif
      enddo

c...  Now do the electrons
      do iy = 0, ny+1 
         dlcon = (lcon(ixmp,iybwmxe)-lcon(ixmp,iybwmne))
         dyyc = (yyc(iybwmxe)-yyc(iybwmne))
	 if(iy > iyso) then    # SOL, so only 1 region; lcon flux func
	    do ix = 0, nx+1  
              if (iy <= iybwmxe) then  # lin interp over banana width
                 lcone(ix,iy) = lcon(ixmp,iybwmne) + 
     .                            dlcon*(yyc(iy)-yyc(iybwmne))/dyyc
              else  
                 lcone(ix,iy) = lcon(ix,iy)
              endif
            enddo
	 else                  # PF and core
            do ix = 0, nx+1    # Now separate into PF and core
              if (ix <= ixpt1(1) .or. ix > ixpt2(1)) then  # PF region
                if(ix <= ixpt1(1)) then
                  ixref = ixpt1(1)
                else # ix >= ixpt2(1)
                  ixref = ixpt2(1)+1
                endif
                if (iy >= iybwmne) then  # lin interp over banana width
                     # PF region is approx, could be improved
                  lcone(ix,iy) = lcon(ixref,iybwmne) + 
     .                            dlcon*(yyc(iy)-yyc(iybwmne))/dyyc
                else
                  lcone(ix,iy) = lcon(ix,iy)
                endif
              else                                         # Core region
                if (iy >= iybwmne) then  # lin interp over banana width
                  lcone(ix,iy) = lcon(ixmp,iybwmne) + 
     .                            dlcon*(yyc(iy)-yyc(iybwmne))/dyyc
                else
                  lcone(ix,iy) = lcon(ix,iy)
                endif
              endif
            enddo
         endif
      enddo

c     Finally, fix-up the nonfunctional ix guard-cells if mhdgeo=2 (annulus)
      if (mhdgeo == 2) then
        do iy = 0, ny+1
	  lconi(0,iy)    = 1.e50
          lcone(0,iy)    = 1.e50
          lconi(nx+1,iy) = 1.e50
          lcone(nx+1,iy) = 1.e50
        enddo
      endif

      return
      end
c ***** end of subroutine conlen ****
c ************************************
c-----------------------------------------------------------------------
      subroutine nonorthg

*     NORTHGD calculates geometrical quantities need for the nonorthogonal
*     grid. The angle at the upper right vertex of the cell, called vtag,
*     where vtag=0 corresponds to orthogonal. Also calculate the weighting
*     functions for linear interpolation of the poloidal variables to give
*     the proper radial derivatives and averages.

      implicit none
      Use(Dim)            # nx,ny,nxpt
      Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2
      Use(Selec)          # ixp1,ixm1
      Use(Comgeo)         # gyf,gxf,dxnog
      Use(Noggeo)         # vtag,angfx,fxm,fx0,fxp,fym,fy0,fyp
      Use(RZ_grid_info)   # rm,zm
      Use(Share)          # geometry,nxc,nxomit,isoldgrid
                          # islimon,ix_lim,iy_lims,reset_core_og
      Use(Bcond)          # isfixlb,isfixrb
      Use(Phyvar)         # pi


*  -- local scalars --
      integer ixu1, iyu1,ixu2, iyu2, ishy, ishx, iym1, isht, jx
      real aa, bb, cc, cos1, cos4, ang1, ang4
      real x1, x2, x3, f1, f2, f3, tantx, tanty
      real slp1, slp0, zmid, rmid, zint, rint, d1, d2, d3
      real slpfs, angfs1, angfs2, denomf
      real rints(0:1), zints(0:1), dyf, dxf, errlim, bigslp, eps
      real z4, r4, delrm, delzm, thetax, thetay
      integer ifxfail, ifyfail, nj, itry, ik
      #Former Aux module variables
      integer ix,iy,ix1,ix3
      data errlim/1.e-10/, bigslp/1.e20/, eps/1e-3/
	  
*=======================================================================
*//computation//

      nj = max(0, nxomit)

c...  Calculate angle at upper-right vertex using the law of cosines
c              c**2 = a**2 + b**2 - 2*a*b*cos(C)
c     where a,b,c are the sides opposite the angles A,B,C of a plane
c     triangle.  In vtag(ix,iy), use the average of the angle at vertex 4
c     of cell (ix,iy) and vertex 1 of cell (ix+1,iy+1).
      do 11 iy = 0, ny
         do 10 ix = 0, nx    # ix=0 and ix=nx are done specially below
c     At vertex 4 of the (ix+nj,iy) cell:
               aa = sqrt( (rm(ix+nj,iy,3)-rm(ix+nj,iy,4))**2 +
     .                    (zm(ix+nj,iy,3)-zm(ix+nj,iy,4))**2 )
               bb = sqrt( (rm(ix+nj,iy,2)-rm(ix+nj,iy,4))**2 +
     .                    (zm(ix+nj,iy,2)-zm(ix+nj,iy,4))**2 )
               cc = sqrt( (rm(ix+nj,iy,3)-rm(ix+nj,iy,2))**2 +
     .                    (zm(ix+nj,iy,3)-zm(ix+nj,iy,2))**2 )
               cos4 = (aa**2 + bb**2 - cc**2)/(2 * aa * bb)
               ang4 = acos(cos4)
c     At vertex 1 of the (ix+1,iy+1) cell:
               aa = sqrt( (rm(ixp1(ix,iy+1)+nj,iy+1,3)
     .                    -rm(ixp1(ix,iy+1)+nj,iy+1,1))**2 +
     .                    (zm(ixp1(ix,iy+1)+nj,iy+1,3)
     .                    -zm(ixp1(ix,iy+1)+nj,iy+1,1))**2 )
               bb = sqrt( (rm(ixp1(ix,iy+1)+nj,iy+1,2)
     .                    -rm(ixp1(ix,iy+1)+nj,iy+1,1))**2 +
     .                    (zm(ixp1(ix,iy+1)+nj,iy+1,2)
     .                    -zm(ixp1(ix,iy+1)+nj,iy+1,1))**2 )
               cc = sqrt( (rm(ixp1(ix,iy+1)+nj,iy+1,3)
     .                    -rm(ixp1(ix,iy+1)+nj,iy+1,2))**2 +
     .                    (zm(ixp1(ix,iy+1)+nj,iy+1,3)
     .                    -zm(ixp1(ix,iy+1)+nj,iy+1,2))**2 )
               cos1 = (aa**2 + bb**2 - cc**2)/(2 * aa * bb)
               ang1 = acos(cos1)
               vtag(ix,iy) = 0.5*pi - 0.5*(ang4 + ang1)
 10      continue
         vtag(nx+1,iy) = vtag(nx,iy)
 11   continue

c...  Now calculate the angles at the x-faces, angfx
      do iy = 0, ny+1
         iym1 = max(0, iy-1)
         do ix = 0, nx+1
            angfx(ix,iy) = 0.5*(vtag(ix,iy)+vtag(ix,iym1))
            if (islimon .eq. 1) then
               if (iy.eq.iy_lims) then
                  angfx(ix,iy) = vtag(ix,iy)
               endif
               if (iy .eq. iy_lims-1) then
                  angfx(ix,iy) = vtag(ix,iym1)
               endif
            endif
         enddo
      enddo

c...  Roundoff problems may cause the ix=0 and ix=nx values of vtag to
c...  be in error.  Thus, linearly extrapolate anfgx at ix=0 and nx with
c...  adjacent values assuming a uniform mesh (gx not defined yet).
c...  Also redo vtag.  A similar situation exists near limiter guard
c...  cells at ix_lim and ix_lim+1. 
      if(redopltvtag == 1) then
       do iy = 0, ny+1
         angfx(0,iy) = 2*angfx(1,iy) - angfx(2,iy)
         vtag(0,iy) = 2*vtag(1,iy) - vtag(2,iy)
         if (nxomit .gt. 0) then
            angfx(0,iy) = angfx(1,iy)
            vtag(0,iy) = vtag(1,iy)
         endif
         if (islimon .ne. 0) then
            angfx(ix_lim-1,iy) = 2*angfx(ix_lim-2,iy)-angfx(ix_lim-3,iy)
            vtag(ix_lim-1,iy) = 2*vtag(ix_lim-2,iy)-vtag(ix_lim-3,iy)
            angfx(ix_lim,iy) = angfx(ix_lim-1,iy)
            vtag(ix_lim,iy) = vtag(ix_lim-1,iy)
            angfx(ix_lim+1,iy) = 2*angfx(ix_lim+2,iy)-angfx(ix_lim+3,iy)
            vtag(ix_lim+1,iy) = 2*vtag(ix_lim+2,iy)-vtag(ix_lim+3,iy)
         endif
         if (nxpt==2) then                # full double null
           angfx(ixrb(1),iy) = 2*angfx(ixrb(1)-1,iy) - angfx(ixrb(1)-2,iy)
           vtag(ixrb(1),iy) = 2*vtag(ixrb(1)-1,iy) - vtag(ixrb(1)-2,iy)
           angfx(ixrb(1)+1,iy) = angfx(ixrb(1),iy)
           vtag(ixrb(1)+1,iy) = vtag(ixrb(1),iy)
           angfx(ixlb(2),iy) = 2*angfx(ixlb(2)+1,iy) - angfx(ixlb(2)+2,iy)
           vtag(ixlb(2),iy) = 2*vtag(ixlb(2)+1,iy) - vtag(ixlb(2)+2,iy)
         endif
         angfx(nx,iy) = 2*angfx(nx-1,iy) - angfx(nx-2,iy)
         vtag(nx,iy) = 2*vtag(nx-1,iy) - vtag(nx-2,iy)
         angfx(nx+1,iy) = angfx(nx,iy)    # not really used
         vtag(nx+1,iy) = vtag(nx,iy)      # not really used
       enddo
      endif

c...  Set angfx and vtag in the y-guard cells to adjacent values
      do ix = 0, nx+1
         angfx(ix,0) = angfx(ix,1)
         angfx(ix,ny+1) = angfx(ix,ny)
         vtag(ix,0) = vtag(ix,1)
         vtag(ix,ny+1) = vtag(ix,ny)
      enddo


c...  reset values around x-point - make orthogonal for now
      do jx = 1, nxpt
         if (ixpt1(jx) .gt. 0) then
            angfx(ixpt1(jx),iysptrx1(jx)-1) = 0.
            angfx(ixpt1(jx),iysptrx1(jx)  ) = 0.
            angfx(ixpt1(jx),iysptrx1(jx)+1) = 0.
         endif
         if (ixpt2(jx) .gt. 0) then
            angfx(ixpt2(jx),iysptrx2(jx)-1) = 0.
            angfx(ixpt2(jx),iysptrx2(jx)  ) = 0.
            angfx(ixpt2(jx),iysptrx2(jx)+1) = 0.
         endif
      enddo

c...  reset values next to cut at ixpt1 or ixpt2 if half-space problem
      if (isfixlb(1).eq.2 .and. ixpt2(1).gt.0) then
         do iy = 0, iysptrx1(1)
            angfx(ixpt2(1),iy) = 0.
          enddo
       endif
       if (isfixrb(1).eq.2 .and. ixpt1(1).gt.0) then
         do iy = 0, iysptrx1(1)
            angfx(ixpt1(1),iy) = 0.
          enddo
       endif
        

c...  reset values for midplane guard cells of double-null configuration
      if ((isudsym==1.or.(geometry.eq.'dnXtarget')) .and. nxc.gt.0) then
         do iy = 0, ny+1
            angfx(nxc-1,iy) = 0.
            angfx(nxc  ,iy) = 0.
            angfx(nxc+1,iy) = 0.
            vtag(nxc-1,iy) = 0.
            vtag(nxc  ,iy) = 0.
            vtag(nxc+1,iy) = 0.
         enddo
      endif      

      do iy = 0, ny+1
        do ix = 0, nx+1
          cosangfx(ix,iy) = cos(angfx(ix,iy))
        end do
      end do

c...  Calculate the fraction-stencil for variables to use at ix-1, ix, ix+1
c...  when forming y-derivatives and averages for the nonorthogonal grid
C...  This is a somewhat complicated loop; we search for the crossing to the
c...  left of center point (ishx=0,ishy=1) using ave of ix+/-1 cells. If fails
c...  look right of the center point (ishx=1). Repeat for ishy=0 (lower iy
c...  intersection).  Consult separate diagrams give more graphical detail of
c...  the y-face normal & cells used to construct 2nd line for interpolation.
c...  We also calculate the stencil for the variables at iy (ishy=0) and 
c...  those at iy+1 (ishy=1)

      ifxfail = 0
      do 21 iy = 0, ny
         do 20 ix = 1, nx
c...        Skip non-physical interior guard cells for dnbot and limiter
c...        and upper divertor plates of a full double null.
c...        Note nxomit added to dnbot cases for outer DNB when nxc=0.
            if (isudsym==1 .and. ix+nxomit.eq.nxc) goto 20
            if (isudsym==1 .and. ix+nxomit.eq.nxc+1) goto 20
            if ((islimon .ne. 0) .and. (ix .eq. ix_lim  )) goto 20
            if ((islimon .ne. 0) .and. (ix .eq. ix_lim+1)) goto 20
            if ((nxpt==2) .and. (ix .eq. ixrb(1)+1)) goto 20
            if ((nxpt==2) .and. (ix .eq. ixlb(nxpt)  )) goto 20
c...  fix possible divide-by-zero
            if ( rm(ix+nj,iy,4) .eq. rm(ix+nj,iy,3) ) then
               slp1 = bigslp
            elseif ( zm(ix+nj,iy,4) .eq. zm(ix+nj,iy,3) ) then
               slp1 = 1/bigslp
            else
               slp1 =(zm(ix+nj,iy,4)-zm(ix+nj,iy,3))/
     .                                   (rm(ix+nj,iy,4)-rm(ix+nj,iy,3))
            endif

c...  Center y-face & ixu1=ix; iyu1 is below face (ishy=0) or above (ishy=1)
c...  iyu2 is opposite iyu1: below if ishy=0 or above if ishy=1
            zmid = 0.5*(zm(ix+nj,iy,4)+zm(ix+nj,iy,3))
            rmid = 0.5*(rm(ix+nj,iy,4)+rm(ix+nj,iy,3))
            ishy = 0
            ixu1 = ix
 13         iyu1 = iy + ishy
            iyu2 = iy + 1 - ishy

c           setup most like diagonal line for proper intersection
            if ( (vtag(ix,iy)+vtag(ixm1(ix,iy),iy))*(1-2*ishy) 
     .                                                  .ge. 0 ) then
               ishx = 1    # ixu1=ix for both ishx=0,1; unlike prev vers.
               ixu2 = ixp1(ix,iyu2)
            else
               ishx = 0
               ixu2 = ixm1(ix,iyu2)
            endif

c ...       Search for normal intersection btwn (ixu1,iyu1) & 2nd pt
c ...       given by ave of (ixu2,iyu2) & (ixu2,iyu2+ishy); see sub lindis
            isht = 1 - ishy
 15         call lindis(ixu1,iyu1,ixu2,iyu2,3,isht,rmid,zmid,slp1,
     .                                               rint,zint,d1,d2,d3)
ccc      if(ix==ixpt1(1)+1 .and. iy==1) then
ccc        write(*,*) "ixu1,iyu1,ixu2,iyu2 = ",ixu1,iyu1,ixu2,iyu2
ccc        write(*,*) "isht,rmid,zmid,slp1 = ",isht,rmid,zmid,slp1
ccc        write(*,*) "rint,zint = ",rint,zint
ccc        write(*,*) "d1,d2,d3 = ",d1,d2,d3
ccc      endif

c...    Confirm intersection, else retry for other diagonal
            if (d1.le.d3*1.0001 .and. d2.le.d3*1.0001) then 
               rints(ishy) = rint 
               zints(ishy) = zint
               fx0(ix,iy,ishy) = d2/d3
               fxm(ix,iy,ishy) =  (1-ishx)*0.5*d1/d3
               fxp(ix,iy,ishy) = ishx*0.5*d1/d3
               fxmy(ix,iy,ishy) = (1-ishx)*0.5*d1/d3
               fxpy(ix,iy,ishy) = ishx*0.5*d1/d3 
               itry = 1
            elseif (itry.eq.1) then  # switch to second diag for intersection
               if (ishx.eq.1) then
                  ishx = 0
                  ixu2 = ixm1(ix,iyu2)
               else
                  ishx = 1
                  ixu2 = ixp1(ix,iyu2)
               endif
               itry = 2
               goto 15
            else
               if (iy.ne.0 .and. iy.ne.ny) then
                  write(*,*) ' fx search failure at ix,iy = ',ix,iy
               endif
               if (iy.eq.0 .and. ishy.eq.1) then
                  write(*,*) ' fx search failure at ix,iy = ',ix,iy
               endif
               if (iy.eq.ny .and. ishy.eq.0) then
                  write(*,*) ' fx search failure at ix,iy = ',ix,iy
               endif
            endif

            if (ishy .eq. 0) then      # shift iy to do iy+1 point
               ishy = 1
               goto 13
            endif
c...  Recalculate the distance between y-points normal to y-face
c...  For isnonog=1, 1/gyf .ne. dynog
            dynog(ix,iy) = sqrt( (rints(1)-rints(0))**2 +
     .                             (zints(1)-zints(0))**2 )
 20      continue
 21   continue

c...  Approx. the gyf in guard cells around boundary; these should not
c...  matter except possibly gyf(0,0) for half-space problem & iflcore=1
      do iy = 0, ny
         dynog(0,iy) = dynog(1,iy)
         dynog(nx+1,iy) = dynog(nx,iy)
      enddo

c...  Approx. gyf in midplane guard cells for geometry=dnbot, only for 
c...  corner boundary condition
      if ((isudsym==1.or.(geometry.eq.'dnXtarget')) .and. nxc.gt.0) then
         do iy = 0, ny
            dynog(nxc,iy) = dynog(nxc-1,iy)
            dynog(nxc+1,iy) = dynog(nxc+2,iy)
         enddo
      endif
c...  Also set dynog(,ny+1) to nonzero; not use - avoid possible 1/dynog
      do ix = 0, nx+1
        dynog(ix,ny+1) = 0.1*dynog(ix,ny)
      enddo

c...  Similarly approximate gyf for limiter guard cells
      if (islimon .ne. 0) then
         do iy = 0, ny
            dynog(ix_lim,iy) = dynog(ix_lim-1,iy)
            dynog(ix_lim+1,iy) = dynog(ix_lim+2,iy)
         enddo
      endif
c...  Similarly approximate gyf for upper target plate guard cells
      if (nxpt==2) then
         do iy = 0, ny
            dynog(ixrb(1)+1,iy) = dynog(ixrb(1),iy)
            dynog(ixlb(2),iy) = dynog(ixlb(2)+1,iy)
         enddo
      endif

      if (ifxfail .gt. 0) then
         write(*,*) 'Total failures to find fxm,fx0,fxp = ', ifxfail
      endif

c...  Calculate the fraction-stencil for variables to use at iy-1, iy, iy+1
c...  when forming derivatives normal to the tilted x-face.
C...  This is a somewhat complicated loop; we search for the crossing to the
c...  above central point (ishy=0) using the central cell and ave of the
c...  2 left cells(ishx=0) or 2 right cells (ishx=1). Repeat process for 
c...  intersection below central cell center. Consult separate diagrams give 
c...  more graphical detail of the y-face normal & cells used to construct 
c...  2nd line for interpolation.
c...  Also calculate the stencil for the variables at ix (ishx=0) and 
c...  those at ixp1 (ishx=1)
c...  We now (8/25/94) include the bend in the flux surface through slpfs,
c...  the slope of the flux surface which can be different on each side
c...  of the x-face.

      ifyfail = 0
      do 41 ix = 0, nx 
ccc	 if (isudsym==1 .and. ix==nxc) goto 41
ccc         # skip non-physical interface between inboard and outboard midplane
         do 40 iy = 1, ny
c...  Fix possible divide-by-zero
           if ( isoldgrid .eq. 1) then #only use to retrieve pre-1/96 solution
            if ( rm(ix+nj,iy,4) .eq. rm(ix+nj,iy,2) .or.
     .              abs(zm(ix+nj,iy,4)-zm(ix+nj,iy,2)) .lt. 1.e-9 ) then
               dxnog(ix,iy) =
     .                sqrt( (rm(ixp1(ix,iy)+nj,iy,0)-rm(ix+nj,iy,0))**2 +
     .                      (zm(ixp1(ix,iy)+nj,iy,0)-zm(ix+nj,iy,0))**2 )
               if (isudsym==1 .and. ix==nxc) dxnog(ix,iy)= 1.e-20  #just small
               goto 40
            endif
           endif

c...  fix possible divide-by-zero; normally compute (R,Z) slope of x-face
           if ( rm(ix+nj,iy,4) .eq. rm(ix+nj,iy,2) ) then
              slp1 = bigslp
           elseif ( zm(ix+nj,iy,4) .eq. zm(ix+nj,iy,2) ) then
              slp1 = 1/bigslp
           else
              slp1 =(zm(ix+nj,iy,4)-zm(ix+nj,iy,2))/
     .                                  (rm(ix+nj,iy,4)-rm(ix+nj,iy,2))
           endif

cc_new            if ( rm(ix+nj,iy,4) .eq. rm(ix+nj,iy,2) ) then
cc_new                slp0 = bigslp
cc_new             elseif ( zm(ix+nj,iy,4) .eq. zm(ix+nj,iy,2) ) then
cc_new                slp0 = 1/bigslp
cc_new             else
cc_new                slp0 =(zm(ix+nj,iy,4)-zm(ix+nj,iy,2))/
cc_new      .                                   (rm(ix+nj,iy,4)-rm(ix+nj,iy,2))
cc_new             endif

c...  Set (R,Z) coordinate at midpoint of x-face
            zmid = 0.5*(zm(ix+nj,iy,4)+zm(ix+nj,iy,2))
            rmid = 0.5*(rm(ix+nj,iy,4)+rm(ix+nj,iy,2))
            ishx = 0

c...  First find intersection to normal to x-face to its left & set tilt param
cc_new            if(angfx(ix,iy) > 0) then
cc_new              ishy = 1  # search upward from (ix,iy) cell
cc_new            else
cc_new              ishy = 0  # search downward
cc_new            endif
cc_new
cc_new            ishx = 0
cc_new            if(angfx(ix,iy) > 0) then
cc_new              ishy = 1  # search upward from (ix,iy) cell
cc_new            else
cc_new              ishy = 0  # search downward
cc_new            endif
cc_new            call lindis2(ix,iy,rmid,zmid,slpl,ishx,ishy)

c...   Setup most likely diagonal for intersection with x-face normal
            iyu1 = iy
 33         ixu1 = (1-ishx)*ix + ishx*ixp1(ix,iyu1)
            if ( angfx(ix,iy)*(1-2*ishx) .ge. 0) then #setup most likely diag
               ishy = 1
               iyu2 = iy+1
               ixu2 = ishx*ix + (1-ishx)*ixp1(ix,iyu2)
            else
               ishy = 0
               iyu2 = iy-1
               ixu2 = ishx*ix + (1-ishx)*ixp1(ix,iyu2)
            endif

c ...       Search for intersection btwn (ixu1,iyu1) & 2nd pt given by
c ...       mid-point of (ixu2,iyu2) & (ixu2+ishx,iyu2) - see sub lindis
            isht = 1 - ishx
 35         call lindis(ixu1,iyu1,ixu2,iyu2,4,isht,rmid,zmid,slp1,
     .                                               rint,zint,d1,d2,d3)
            if (d1.le.d3*1.0001 .and. d2.le.d3*1.0001) then
               rints(ishx) = rint
               zints(ishx) = zint
               fy0(ix,iy,ishx) = d2/d3
               fym(ix,iy,ishx) = (1-ishy)*0.5*d1/d3
               fyp(ix,iy,ishx) = ishy*0.5*d1/d3
               fymx(ix,iy,ishx) = (1-ishy)*0.5*d1/d3
               fypx(ix,iy,ishx) = ishy*0.5*d1/d3
               itry = 1
            elseif (itry.eq.1) then  # switch to second diag for intersection
               if (ishy.eq.1) then
                  ishy = 0
                  iyu2 = iy-1
                  ixu2 = ishx*ix + (1-ishx)*ixp1(ix,iyu2)
               else
                  ishy = 1
                  iyu2 = iy+1
                  ixu2 = ishx*ix + (1-ishx)*ixp1(ix,iyu2)
               endif
               itry = 2
               goto 35
            else        # failed to find proper pt.
              do jx = 1, nxpt
                if (ix>ixlb(jx) .and. ix<ixrb(jx)) then
                  write(*,*) ' fy search failure at ix,iy = ',ix,iy
                endif
                if (ix==ixlb(jx) .and. ishx.eq.1) then
                  write(*,*) ' fy search failure at ix,iy = ',ix,iy
                endif
                if (ix==ixrb(jx) .and. ishx.eq.0) then
                  write(*,*) ' fy search failure at ix,iy = ',ix,iy
                endif
              enddo
            endif

            if (ishx .eq. 0) then    # shift ix to do ix+1 point
               ishx = 1
               goto 33
            endif
c...  Calculate the distance between interpolated pts normal to x-face
c...  Include bend in flux surface through angles angfs1,2; here cmt out.
c...  For isnonog=1, 1/gxf .ne. dxnog as interp pts not along cell ctrs


            dxnog(ix,iy) = sqrt( (rints(1)-rints(0))**2 + 
     .                           (zints(1)-zints(0))**2 )
ccc     .                  / abs( cos(0.5*(angfs1-angfs2)) )
ccc            gxfn(ix,iy) = 1/dxf
            
 40      continue
 41   continue

      if (ifyfail .gt. 0) then
         write(*,*) 'Total failures to find fym,fy0,fyp = ', ifyfail
      endif

c...  if isupstreamx=1, cells above x-point should be orthogonal; reset
      if (reset_core_og == 1) then 
        do jx = 1, nxpt
          do iy = 0, ny+1
            do ix = ixpt1(jx)+1, ixpt2(jx)
              do ik = 0, 1
                fxm(ix,iy,ik)  = 0.
                fx0(ix,iy,ik)  = 1.
                fxp(ix,iy,ik)  = 0.
                fxmy(ix,iy,ik) = 0.
                fxpy(ix,iy,ik) = 0.
                fxm(ix,iy,ik)  = 0.
                fx0(ix,iy,ik)  = 1.
                fxp(ix,iy,ik)  = 0.
                fxmy(ix,iy,ik) = 0.
                fxpy(ix,iy,ik) = 0.
              enddo
            enddo
          enddo
        enddo
      endif

c...  Reset plate guard-cell fy to orthogonal stencil
      do jx = 1, nxpt
        do iy = 0, ny+1
           ix = ixlb(jx) 
              fym(ix,iy,0)  = 0.
              fy0(ix,iy,0)  = 1.
              fyp(ix,iy,0)  = 0.
              fymx(ix,iy,0) = 0.
              fypx(ix,iy,0) = 0.
              fym(ix,iy,0)  = 0.
              fy0(ix,iy,0)  = 1.
              fyp(ix,iy,0)  = 0.
              fymx(ix,iy,0) = 0.
              fypx(ix,iy,0) = 0.
           ix = ixrb(jx)
              fym(ix,iy,1)  = 0.
              fy0(ix,iy,1)  = 1.
              fyp(ix,iy,1)  = 0.
              fymx(ix,iy,1) = 0.
              fypx(ix,iy,1) = 0.
              fym(ix,iy,1)  = 0.
              fy0(ix,iy,1)  = 1.
              fyp(ix,iy,1)  = 0.
              fymx(ix,iy,1) = 0.
              fypx(ix,iy,1) = 0.
        enddo
      enddo

      return
      end
c ***** end of subroutine nonorthg ****************
c ---------------------------------------------------------------------
      subroutine lindis (i1,j1,i2,j2,ipos2,ish,r0,z0,slp1,rx,zx,d1,d2,d3)

******************************************************************************
*     LINDIS finds the intersection of two lines, one defined by the
*     perpendicular to a face and obeying the equation z=z0-(r-z0)/slp1,
*     where (r0,z0) gives the center of the face and slp1 is the slope of the
*     face; thus, -1/slp1 is the slope of the perpendicular. The second line
*     connects two other points on the mesh given by [rm(i1,j1),zm(i1,j1)]
*     and a second point defined by the switch ipos2:
*        ipos2 = 0 2nd point is the cell center [rm(i2,j2),zm(i2,j2)] 
*        ipos2 = 1 2nd point is y-face for [rm(i1,j1),zm(i1,j1)]
*        ipos2 = 2 2nd point is x-face for [rm(i1,j1),zm(i1,j1)]
*        ipos2 = 3 2nd point is y-face for [rm(i2,j2),zm(i2,j2)]
*        ipos2 = 4 2nd point is x-face for [rm(i2,j2),zm(i1,j2)]
*     If ish=0, take face to increasing ix or iy; if ish=1, take face in
*     decreasing direction. LINDIS returns the values at the intersection
*     (rx,zx), and the distances d1, d2, and d3:
*
*        d1 = distance between (rx,zx) and (r1,z1), 
*        d2 = distance between (rx,zx) and (r2,z2) 
*        d3 = distance between (r1,z1) and (r2,z2), 
****************************************************************************
      implicit none
      Use(Dim)            # nx,ny
      Use(RZ_grid_info)   # rm,zm
      Use(Share)          # nxomit

*  -- Input scalars --
      integer i1,i2,j1,j2,     # indices of pts (i1,j1) & (i2,j2) & face switch
     .        ipos2, ish
      real slp1,r0,z0           # slope of face, and perp. midpoint (r0,z0)

*  -- Output scalars --          
      real rx,zx,d1,d2,d3    #intersection (rx,zx) and distances (see above)

*  -- Local scalars --
      real slp2,r1,z1,r2,z2,r2cc,z2cc,dface_cc
      integer nj
 
      nj = max(0, nxomit)

      r1 = rm(i1+nj,j1,0)
      z1 = zm(i1+nj,j1,0)
      r2cc = rm(i2+nj,j2,0)            # cell-center r for pt 2
      z2cc = zm(i2+nj,j2,0)            # cell-center z for pt 2
      if (ipos2 .eq. 0) then        # 2nd point is cell center
         r2 = rm(i2+nj,j2,0)
         z2 = zm(i2+nj,j2,0)
      elseif (ipos2 .eq. 1) then    # 2nd point is y-face of 1st point
         r2 = 0.5*(rm(i1+nj,j1,3-2*ish) + rm(i1+nj,j1,4-2*ish))
         z2 = 0.5*(zm(i1+nj,j1,3-2*ish) + zm(i1+nj,j1,4-2*ish))
      elseif (ipos2 .eq. 2) then    # 2nd point is x-face of 1st point
         r2 = 0.5*(rm(i1+nj,j1,2-ish) + rm(i1+nj,j1,4-ish))
         z2 = 0.5*(zm(i1+nj,j1,2-ish) + zm(i1+nj,j1,4-ish))
      elseif (ipos2 .eq. 3) then    # 2nd point is y-face of 2nd point
         r2 = 0.5*(rm(i2+nj,j2,3-2*ish) + rm(i2+nj,j2,4-2*ish))
         z2 = 0.5*(zm(i2+nj,j2,3-2*ish) + zm(i2+nj,j2,4-2*ish))
      elseif (ipos2 .eq. 4) then    # 2nd point is x-face of 2nd point
         r2 = 0.5*(rm(i2+nj,j2,2-ish) + rm(i2+nj,j2,4-ish))
         z2 = 0.5*(zm(i2+nj,j2,2-ish) + zm(i2+nj,j2,4-ish))
      endif

      slp2 = (z1-z2) / (r1-r2+1.e-20)
      if (abs(slp1) .gt. 1.e-9) then   # use expansion for small slp1
         rx = (z0 +r0/slp1 -z1 +slp2*r1) / (slp2 + 1/slp1)
         zx = z0 - (rx-r0) / slp1
      else
         rx = r0
         zx = z1 -slp2*(r1-r0)
         if (abs(slp2).gt.1.e9) then
            write(*,*) '* Warning: slp2*slp1 indeterminante in lindis *'
         endif
      endif
      d1 = sqrt((rx-r1)**2 + (zx-z1)**2)
      d2 = sqrt((rx-r2)**2 + (zx-z2)**2)
      d3 = sqrt((r1-r2)**2 + (z1-z2)**2)

      if (ipos2.eq.1 .or. ipos2.eq.2) then # add face-to-center distance
         dface_cc =  sqrt((r2-r2cc)**2 + (z2-z2cc)**2)
         d2 = d2 + dface_cc
         d3 = d3 + dface_cc
      endif

      return
      end
c ***** end of subroutine lindis *****************
c *****
      subroutine lindis2 (p1,q1,rmf,zmf,slpmf,isx,isy,idone)
C  THIS IS WORK-IN-PROGRESS FOR NEW INTERPOLATOR: NOT CALLED
c
******************************************************************************
*     LINDIS2 finds the intersection of two lines, one defined by the
*     perpendicular to a face obeying the equation z=zmf-(r-zmf)/slpmf,
*     where (rmf,zmf) gives the center of the face and slpmf is the face slope
*     face; thus, -1/slpmf is the slope of the perpendicular. The second line
*     connects two other points on the mesh as is defined below - search along
*     line segments between cell centers and midpoints of faces to follow
*     along the mesh lines even if the cells form a curved surface. The index
*     isx (ishx) specifies search to left of face (isx=0) or to the right 
*     (isx=1).  Also, isy is determined by angfx > 0 or < 0, which determines
*     the direction of the search (upward if isy=1, downward if isy=0).
*     See related figures of the stencil patterns for a visual picture.
*        d1 = distance between (rx,zx) and (r1,z1), 
*        d2 = distance between (rx,zx) and (r2,z2) 
*        d3 = distance between (r1,z1) and (r2,z2), 
****************************************************************************
      implicit none
      Use(Dim)            # nx,ny
      Use(RZ_grid_info)   # rm,zm
      Use(Share)          # nxomit
      Use(Noggeo)         # fy0,fym,fyp etc
      Use(Comgeo)         # dxnog
      Use(Selec)          # ixp1 

*  -- Input scalars --
      integer p1,q1,isx,isy          # indices of pts (p1,q1) & (i2,j2) & face switch
      real rmf,zmf,slpmf       # slope of face, and midpoint (rmf,zmf)

*  -- Output scalars --          
      real rx,zx,d12,d1x,d2x #intersection (rx,zx) and distances (see above)
      integer idone

*  -- Local scalars --
      real r1,z1,r2,z2	     #coordinates of two pts defining 2nd line
                             #connecting alternate cell ctrs & cell faces
      real slpn,slp2         #slope of normal to x-face (=-11/slpmf)
      real dcf
      integer nj,isg,ixc,iyc
 
      idone = 0
      nj = max(0, nxomit)

c... ************************************************************************
c...  Search up (isy=0) or down (isy=1) in iy for fix ix, then ix at fixed iy
c... ************************************************************************
      do isg = 1, 3   #consider 3 center-to-face vert sgmts for intersection
	iyc = q1 + sign(isg/2,-isy)
        if(isy==0) then
          ixc = p1 + nj
        else
          ixc = ixp1(p1+nj,iyc)
        endif
        if (isg/2 == (isg+1)/2) then   #isg even, r1,z1 on y-face
          r1 = 0.5*(rm(ixc,iyc,2+2*isy) + rm(ixc,iyc,1+2*isy)) # y-face
          z1 = 0.5*(zm(ixc,iyc,2+2*isy) + zm(ixc,iyc,1+2*isy)) # y-face
        else    #isg odd, so r1,z1 is cell ctr
          r1 = rm(ixc,iyc,0)	# cell ctr
          z1 = zm(ixc,iyc,0)	# cell ctr
        endif
c...  First line segment btw (r1,z1) y-face below (iys=0) or above (iys=1)
c...  Find slope of second line, and solve for intersection (rx,zx)
        if(isg == 1) then  #r2,z2 only set for isg=1; then use last r1,z1
          r2 =0.5*(rm(ixc,iyc,2+2*isy) + rm(ixc,iyc,1+2*isy)) # y-face
          z2 =0.5*(zm(ixc,iyc,2+2*isy) + zm(ixc,iyc,1+2*isy)) # y-face
        endif
        slp2 = (z1-z2) / (r1-r2+1.e-20)	 #slope of r1,z1-r2,z2 line segment
        slpn = -1./(slpmf+1.e-20)        #slope of normal to x-face
c...  Solution of 2 linear eqns (lines) for intersection rx,zx
        zx = (slpn*zmf - slp2*z2 + r2 - rmf)/(slpn-slp2+1.e-20)
        rx = slpn*(zx-zmf) + rmf

c...  Compute distances btwn pts 1,2, and btwn intersect rx,zx & pts 1,2
        d12 = sqrt((r1-r2)**2 + (z1-z2)**2)
        d1x = sqrt((rx-r1)**2 + (zx-z1)**2)
        d2x = sqrt((rx-r2)**2 + (zx-r2)**2)

c...  Set dist btwn cell ctr & adj face just beyond line segmt 1,2
        if(isg==1) then      #here r2,z2 is lower (upper) iy-face
          dcf = sqrt((r2-rm(ixc,iyc-1+2*isy,0))**2 +
     .               (z2-zm(ixc,iyc-1+2*isy,0))**2)
        elseif(isg==2) then  #here r1,z1 is upper (lower) iy-face
          dcf = sqrt((r1-rm(ixc,iyc-2,0))**2 +
     .               (z1-zm(ixc,iyc-2,0))**2)
        elseif(isg==3) then  #here r2,z2 is upper (lower) iy-face
          dcf = sqrt((r2-rm(ixc,iyc-1+isy,0))**2 +
     .               (z2-zm(ixc,iyc-1+isy,0))**2)
        endif
c...  Check if d1x and d2x < d12; if so, record stencil  values & exit
        if (d12 > d1x .and. d12 > d2x) then  #note then d12=d1x+d2x
          dxnog(p1+nj,q1) = dxnog(p1+nj,q1)+sqrt((rmf-rx)**2+(zmf-zx)**2)
          if(isy==0) then
            if(isg==1) then
              fym(p1+nj,q1,0) = d1x/(dcf+d12)
              fy0(p1+nj,q1,0) = (dcf+d2x)/(dcf+d12)
            elseif(isg==2) then
              fy0(p1+nj,q1,0) = (dcf+d1x)/(dcf+d12)
              fyp(p1+nj,q1,0) = d2x/(dcf+d12)
            elseif(isg==3) then
              fy0(p1+nj,q1,0) = d1x/(dcf+d12)
              fyp(p1+nj,q1,0) = (dcf+d2x)/(dcf+d12)
            endif
          elseif(isy==1) then
            if(isg==1) then
              fyp(p1+nj,q1,1) = d1x/(dcf+d12)
              fy0(p1+nj,q1,1) = (dcf+d2x)/(dcf+d12)
            elseif(isg==2) then
              fy0(p1+nj,q1,1) = (dcf+d1x)/(dcf+d12)
              fym(p1+nj,q1,1) = d2x/(dcf+d12)
            elseif(isg==3) then
              fy0(p1+nj,q1,1) = d1x/(dcf+d12)
              fym(p1+nj,q1,1) = (dcf+d2x)/(dcf+d12)
            endif
          endif
          idone = 1
          break
        endif     #testing if valid intersection is found

        r2 = r1
        z2 = z1
      enddo

c... ************************************************************************
c...  Similar coding to search right/left if no intersect in iy (idone=0)
c... ************************************************************************
      if(idone == 1) return
      do isg = 1,2   #consider 2 center-to-face horiz sgmts for intersection
        ixc = p1 + nj #sign(1,-isy)*isg/2
	iyc = q1 + 1 - 2*isy
        if (isg == 2) then   #isg even, r1,z1 at cell ctr
          if(isy==0) ixc = ixp1(ixc,iyc)
          r1 = rm(ixc,iyc,0)	# cell ctr
          z1 = zm(ixc,iyc,0)	# cell ctr
        else    #isg odd (=1), so r1,z1 on cell face
          if(isy==1) ixc = ixp1(ixc,iyc)
          r1 =0.5*(rm(ixc,iyc,2-isy) + rm(ixc,iyc,4-isy)) # x-face
          z1 =0.5*(zm(ixc,iyc,2-isy) + zm(ixc,iyc,4-isy)) # x-face
        endif
c...  First line segment btw (r1,z1) x-face to right (iys=0) or left (iys=1)
c...  Find slope of second line, and solve for intersection (rx,zx)
        if(isg == 1) then  #r2,z2 only set for isg=1; then use last r1,z1
          r2 = rm(ixc,iyc,0)            # cell ctr
          z2 = zm(ixc,iyc,0)            # cell ctr
        endif
        slp2 = (z1-z2) / (r1-r2+1.e-20)	 #slope of line segment
        slpn = -1./(slpmf+1.e-20)        #slope of normal to x-face
        zx = (slpn*zmf - slp2*z2 + r2 - rmf)/(slpn-slp2+1.e-20)
        rx = slpn*(zx-zmf) + rmf

c...  Compute distances btwn pts 1,2, and btwn intersect rx,zx & pts 1,2
        d12 = sqrt((r1-r2)**2 + (z1-z2)**2)
        d1x = sqrt((rx-r1)**2 + (zx-z1)**2)
        d2x = sqrt((rx-r2)**2 + (zx-r2)**2)

c...  If isg=1, first set dist btwn right-most(left-most)cell ctr & adj face
c...  Here for isg=1, r1,z1 is cell x-face, so rm,zm is shifted cell ctr
        if(isg==1) then
          if(isy==0) ixc=ixp1(p1+nj,iyc)
          if(isy==1) ixc=p1+nj
          dcf = sqrt((r1-rm(ixc,iyc,0))**2 + (z1-zm(ixc,iyc,0))**2)
        elseif(isg==2) then
          if(isy==0) ixc=p1+nj
          if(isy==1) ixc=ixp1(p1+nj,iyc)
          dcf = sqrt((r2-rm(ixc,iyc,0))**2 + (z2-zm(ixc,iyc,0))**2)
        endif

c...  Check if d1x and d2x < d12; if so, record stencil values & exit
        if (d12 > d1x .and. d12 > d2x) then  #note then d12=d1x+d2x
          dxnog(p1+nj,q1) = dxnog(p1+nj,q1)+sqrt((rmf-rx)**2+(zmf-zx)**2)
          if(isy==0) then
            if(isg==1) then
              fyp(p1+nj,q1,0) = d2x/(dcf+d12)
              fypx(p1+nj,q1,0) = (dcf+d1x)/(dcf+d12)
            elseif(isg==2) then
              fyp(p1+nj,q1,0) = d1x/(dcf+d12)
              fypx(p1+nj,q1,0) = (dcf+d2x)/(dcf+d12)
            endif
          elseif(isy==1) then
            if(isg==1) then
              fym(p1+nj,q1,1) = d2x/(dcf+d12)
              fymx(p1+nj,q1,1) = (dcf+d1x)/(dcf+d12)
            elseif(isg==2) then
              fym(p1+nj,q1,1) = d1x/(dcf+d12)
              fymx(p1+nj,q1,1) = (dcf+d2x)/(dcf+d12)           
            endif
          endif
          idone = 1
          break
        endif     #testing if valid intersection is found

        r2 = r1
        z2 = z1
      enddo

      return
      end
c ***** end of subroutine lindis2 *****************
c ---------------------------------------------------------------------

      subroutine grdnrm (nx,ny,ixlb,ixpt1,ixpt2,ixrb,iysptrx,isgindx,
     .                                  gx,gxf,gy,gyf,gyc,xn,xvn,yn,yvn)

c...  This subroutine constructs the normalized poloidal and radial grid
c...  point locations xn, xvn, yn, & yvn used for interpolating to a different 
c...  mesh size. Three poloidal regions are used for iy .le. iysep, and one
c...  for iy .gt. iysep. If isgindx=1, a uniform (index) mesh is assumed.

      implicit none

c  -- Input variables - geometrical distances
      integer nx, ny, ixlb, ixpt1, ixpt2, ixrb, iysptrx, isgindx
      real gx(0:nx+1,0:ny+1), gxf(0:nx+1,0:ny+1), gy(0:nx+1,0:ny+1),
     .     gyf(0:nx+1,0:ny+1), gyc(0:nx+1,0:ny+1)

c  -- Output varibles - normalized cell locations
      real xn(0:nx+1,0:ny+1), xvn(0:nx+1,0:ny+1), 
     .     yn(0:nx+1,0:ny+1), yvn(0:nx+1,0:ny+1)

c  -- local variables
      integer ix, iy, ixpt1t, ixpt2t
      real ytot, yvtot, xtot, xvtot, xns, xvns, yns, yvns, delyn, delxn
      integer iyst

      do 40 ix = 0, nx+1

      if (isgindx .eq. 0) then    # use true spatial distances for interp

         yn(ix,0) = 0.0e0
         yvn(ix,0) = 0.0e0
         do 10 iy = 1, ny+1
            yn(ix,iy) = yn(ix,iy-1) + 1.0e0/gyf(ix,iy-1)
            yvn(ix,iy) = yvn(ix,iy-1) + 0.5e0*(gyc(ix,iy-1)+gyc(ix,iy)) /
     .                                         (gyc(ix,iy-1)*gyc(ix,iy))
 10      continue

c...  normalize 0 -> 1, in two intervals: core-pf and sol

         yns = yn(ix,0)
         yvns = yvn(ix,0)
         ytot = 0.5e0/gy(ix,iysptrx) + yn(ix,iysptrx) - yns
         yvtot = 0.5e0/gyc(ix,iysptrx) + yvn(ix,iysptrx) - yvns
         do 20 iy = 1, iysptrx
            yn(ix,iy) = (yn(ix,iy) - yns) / ytot
            yvn(ix,iy) = (yvn(ix,iy) - yvns) / yvtot
 20      continue
         yns = yn(ix,iysptrx+1) - 0.5e0/gy(ix,iysptrx+1)
         yvns = yvn(ix,iysptrx+1) -0.5e0/gyc(ix,iysptrx+1)
         ytot = yn(ix,ny+1) - yns
         yvtot = yvn(ix,ny+1) - yvns
         do 30 iy = iysptrx+1, ny+1
            yn(ix,iy) = (yn(ix,iy) - yns) / ytot
            yvn(ix,iy) = (yvn(ix,iy) - yvns) / yvtot
 30      continue

      else     # large isgindx if-test      
c...  use index-based mesh if isgindx .ne. 0

         yn(ix,0) = 0.0e0
         yvn(ix,0) = 0.0e0
         if (iysptrx .gt. 0) then
            delyn = 1.0e0 / float(iysptrx)
            yn(ix,1) = 0.5e0 * delyn
            yvn(ix,1) = 0.5e0 * delyn
            do 34 iy = 2, iysptrx
               yn(ix,iy) = yn(ix,iy-1) + delyn
               yvn(ix,iy) = yn(ix,iy)
 34         continue
         endif

         if (iysptrx .lt. ny) then
            delyn = 1.0e0 / float(ny-iysptrx)
            yn(ix,iysptrx+1) = 0.5e0 * delyn
            yvn(ix,iysptrx+1) = 0.5e0 * delyn
            do 36 iy = iysptrx+2, ny
               yn(ix,iy) = yn(ix,iy-1) + delyn
               yvn(ix,iy) = yn(ix,iy)
 36         continue
         endif
         yn(ix,ny+1) = 1.e0
         yvn(ix,ny+1) = 1.e0


      endif    # large isgindx if-test
         
 40   continue

c...  For iy<iysep, normalize the poloidal or x-grid in three regions: left 
c...  leg, core, and right leg. We don't worry about using indirect addressing
c...  near the cuts as we normalize each region and any offset error 
c...  resulting from using gx(ix-1,iy) instead of gx(ixm1,iy) is unimportant
c...  For iy>iysep, only one poloidal region is used.

      if (isgindx .eq. 0) then     # use true spatial distances for interp 
                                   # large if-test on isgindx
      do 55 iy = 0, ny+1
         xn(0,iy) = 0.
         xvn(0,iy) = 0.
         do 50 ix = 1, nx+1
            xn(ix,iy) = xn(ix-1,iy) + 1.0e0/gxf(ix-1,iy)
            xvn(ix,iy) = xvn(ix-1,iy) + 1.0e0/gx(ix,iy)
 50      continue
 55   continue

c...  Normalize from 0 -> 1 in three regions: 0-ixpt1,ixpt1+1-ixpt2, and
c...  ixpt2+1-nx+1 if iy .le. iysptrx

      if (iysptrx .eq. 0) goto 95
      do 90 iy = 0, iysptrx
         xns = xn(0,iy)
         xvns = xvn(0,iy)
         if (ixpt1 .le.0) goto 65
         xtot = xn(ixpt1,iy) + 0.5e0/gx(ixpt1,iy) - xns
         xvtot = xvn(ixpt1,iy) + 0.5e0/gx(ixpt1,iy) - xvns
         do 60 ix = 1, ixpt1
            xn(ix,iy) = (xn(ix,iy) - xns) / xtot
            xvn(ix,iy) = (xvn(ix,iy) - xvns) / xvtot
 60      continue
         xns = xn(ixpt1+1,iy) - 0.5e0/gx(ixpt1+1,iy)
         xvns = xvn(ixpt1+1,iy) - 0.5e0/gx(ixpt1+1,iy)
 65      ixpt2t = min(ixpt2, nx+1)
         xtot = xn(ixpt2t,iy) + 0.5e0/gx(ixpt2t,iy) - xns
         xvtot = xvn(ixpt2t,iy) + 0.5e0/gx(ixpt2t,iy) - xvns
         do 70 ix = max(ixpt1+1, 0), ixpt2t
            xn(ix,iy) = (xn(ix,iy) - xns) / xtot
            xvn(ix,iy) = (xvn(ix,iy) - xvns) / xvtot
 70      continue
         if (ixpt2 .ge. nx) goto 90
         xns = xn(ixpt2+1,iy) - 0.5e0/gx(ixpt2+1,iy)
         xvns = xvn(ixpt2+1,iy) - 0.5e0/gx(ixpt2+1,iy)
         xtot = xn(nx+1,iy) - xns
         xvtot = xvn(nx+1,iy) - xvns
         do 80 ix = max(ixpt2+1, 0), nx+1
            xn(ix,iy) = (xn(ix,iy) - xns) / xtot
            xvn(ix,iy) = (xvn(ix,iy) - xvns) / xvtot
 80      continue   
 90   continue
 95   continue

c...  For iy .gt. iysptrx, normalize the single poloidal region
      iyst = iysptrx+1
      if (iysptrx .eq. 0) iyst = 0
      do 130 iy = iyst, ny+1
         xns = xn(0,iy)
         xvns = xvn(0,iy)
         xtot = xn(nx+1,iy) - xns
         xvtot = xvn(nx+1,iy) - xvns
         do 120 ix = 0, nx+1
            xn(ix,iy) = (xn(ix,iy) - xns) / xtot
            xvn(ix,iy) = (xvn(ix,iy) - xvns) / xvtot
 120     continue
 130  continue

      else    # treat iy<iysptrx & iy>iysptrx the same
              # large if-test on isgindx  

      do iy = 0, ny+1
         if (ixpt1 .le. ixlb) goto 315
         delxn = 1.e0 / (float(ixpt1-ixlb) + 1.e-20)
         xn(ixlb,iy) = 0.e0
         xn(ixlb+1,iy) = 0.5e0 * delxn
         xvn(ixlb,iy) = 0.25e0 * delxn
         xvn(ixlb+1,iy) = 1.e0 * delxn
         do ix = ixlb+2, ixpt1
            xn(ix,iy) = xn(ix-1,iy) + delxn
            xvn(ix,iy) = xvn(ix-1,iy) + delxn
         enddo
 315     continue
         if (ixpt2 .le. 0) goto 325
         ixpt1t = max(ixpt1, ixlb)
         ixpt2t = min(ixpt2, ixrb)
         delxn = 1.e0 / (float(ixpt2t-ixpt1t) + 1.e-20)
         xn(ixpt1t+1,iy) = 0.5e0 * delxn
         xvn(ixpt1t+1,iy) = 0.5e0 * delxn
         do ix = ixpt1t+2, ixpt2t
            xn(ix,iy) = xn(ix-1,iy) + delxn
            xvn(ix,iy) = xvn(ix-1,iy) + delxn
         enddo
 325     continue
         if (ixpt2 .ge. ixrb) goto 332
         ixpt2t = max(ixpt2, ixlb)
         delxn = 1.e0 / (float(ixrb-ixpt2t) + 1.e-20)
         xn(ixpt2t+1,iy) = 0.5e0 * delxn
         xvn(ixpt2t+1,iy) = 0.5e0 * delxn
         do ix = ixpt2t+2, ixrb
            xn(ix,iy) = xn(ix-1,iy) + delxn
            xvn(ix,iy) = xvn(ix-1,iy) + delxn
         enddo
 332     xn(ixrb+1,iy) = 1.e0 
         xvn(ixrb,iy) = 0.999999e0
         xvn(ixrb+1,iy) = 1.e0
      enddo

      endif     # large if-test on isgindx        

      return
      end   
c **** end of subroutine grdnrm ****
c **********************************
c-----------------------------------------------------------------------
c...  Subroutine to shift cell vertices near ix=ixpt2 for 1D flux tube ...
      subroutine reset1dmeshpt

      implicit none
      Use(Dim)            # nx,ny
      Use(RZ_grid_info)   # rm,zm
      Use(Xpoint_indices) # ixpt2
      Use(Share)          # nxomit,nxpt2msh,nxpt2psh,rxpt2msh,zxpt2psh

c...  local scalars
      integer ix,ix2,ix1 
      real lenginc,theta

      if (nxomit .ge. 0) then
        ix2 = nxomit+ixpt2(1)    		#consider outer flux-surface only
      else
	ix2 = ixpt2(1)
      endif

c...  First fix iy=1 cell vertices above X-point (ixpt2)
      theta = rxpt2msh/zxpt2msh   #defines curvature of flux surface
      do ix = ix2, ix2-nxpt2msh+1, -1
        rm(ix,1,3) = rm(ix,1,4)            # + zxpt2msh*theta
        zm(ix,1,3) = zm(ix,1,4) + zxpt2msh #+ zxpt2msh*(1-0.5*theta**2)
        rm(ix-1,1,4) = rm(ix,1,3)
        zm(ix-1,1,4) = zm(ix,1,3)
      enddo

c...  Second, fix iy=1 cell vertices below X-point
      theta = rxpt2psh/zxpt2psh   #defines curvature of flux surface
      do ix = ix2+1, ix2+nxpt2psh
        rm(ix,1,4) = rm(ix,1,3)             #  + zxpt2psh*theta
        zm(ix,1,4) = zm(ix,1,3) - zxpt2psh  #zxpt2psh*(1.-0.5*theta**2)
        rm(ix+1,1,3) = rm(ix,1,4)
        zm(ix+1,1,3) = zm(ix,1,4)
      enddo

#..   Do the inner X-pt
#..
#..
      if (nxomit .ge. 0) then
        ix1 = nxomit+ixpt1(1)+1
      else
        ix1 = ixpt1(1)+1
      endif
c...  First fix iy=1 cell vertices above X-point (ixpt2)
      theta = rxpt1msh/zxpt1msh   #defines curvature of flux surface
      do ix = ix1, ix1+nxpt1msh-1
        rm(ix,1,4) = rm(ix,1,3)
        zm(ix,1,4) = zm(ix,1,3) + zxpt1msh
        rm(ix-1,1,3) = rm(ix,1,4)
        zm(ix-1,1,3) = zm(ix,1,4)
      enddo

c...  Second, fix iy=1 cell vertices below X-point
      theta = rxpt1psh/zxpt1psh   #defines curvature of flux surface
      do ix = ix1-1, ix1-nxpt1psh, -1
        rm(ix,1,3) = rm(ix,1,4)
        zm(ix,1,3) = zm(ix,1,4) - zxpt1psh
        rm(ix+1,1,4) = rm(ix,1,3)
        zm(ix+1,1,4) = zm(ix,1,3)
      enddo

      return
      end
c***  End of subroutine reset1dmeshpt ****
c*****************************************

    
