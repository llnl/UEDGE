c-----------------------------------------------------------------------
      subroutine gridseq

*     GRIDSEQ provides saving of variables for grid interp. and restart

      implicit none

      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)   # ixlb,ixpt1,ixpt2,ixrb,iysptrx
      Use(Interp)   # ixlbo,ixpt1o,ixpt2o,ixrbo,iysptrxo
                    # xnrmo,xvnrmo,ynrmo,yvnrmo,
                    # nis,tes,tis,tgs,ups,phis,ngs,isimesh,afracs,
                    # ixst,ixend,ixsto,ixendo
      Use(Imprad)   # afrac,isimpon
      Use(Compla)   # ni,up,ng,te,ti,phi
      Use(Comgeo)   # xnrm,xvnrm,ynrm,yvnrm
      Use(Share)    # nyomitmx
      Use(Npes_mpi) # ismpion
      Use(Cut_indices)	# ixcut1,iycut1,ixcut2,iycut2,ixcut3,iycut3
                        # ixcut4,iycut4

      integer ifld,jx,ir
      #Former Aux module variables
      integer igsp,iy

      if (ismpion == 0) isimesh = 0  #turn-off switch; at least one mesh calc

      iysptrxo = iysptrx
      do jx=1,nxpt
         ixlbo(jx) = ixlb(jx)
         ixpt1o(jx) = ixpt1(jx)
         ixpt2o(jx) = ixpt2(jx)
         ixrbo(jx) = ixrb(jx)
      enddo

c...  Order poloidal regions from old mesh (as new mesh in ueinit)
      ixsto(1) = ixlb(1)
      ixendo(1) = ixcut1
      if (ixlb(1) == 0 .and. ixcut1 == 0) then  # no inner leg
        ixsto(2) = 0
      else
        ixsto(2) = max(ixlb(1), ixcut1+1)
      endif
      ixendo(2) = ixcut2
      if (nyomitmx >= nysol(1)) then   # no inner/outer leg region
         ixsto(2) = 0
         ixendo(2) = nx+1
      endif
      if (nx == 1) then  #special case: 1D in radial direction
         ixendo(2) = 2
      endif
c..   Now need to check if ixrb is > or < ixcut3
      ixsto(3) = ixcut2+1
      if (ixcut3 > ixrb(1) .or. nxpt==1) then  #3 regions in first domain
        ixendo(3) = ixrb(1)+1
      else    # 4 regions in 1st domain, end on ixcut3
        ixendo(3) = ixcut3
      endif

c..   Continue ordering if double null or snowflake
      if (nxpt == 2) then
        if (ixcut3 > ixrb(1)) then  # do 3-region 2nd domain
          ixsto(4) = ixlb(2)
          ixendo(4) = ixcut3
          ixsto(5) = ixcut3+1
        else # 4 regions in 1st domain, compl & do 2-region 2nd domain
          ixsto(4) = ixcut3+1
          ixendo(4) = ixrb(1)+1
          ixsto(5) = ixlb(2)
        endif  # remain indices are the same
          ixendo(5) = ixcut4
          ixsto(6) = ixcut4+1
          ixendo(6) = ixrb(2)+1
      endif  # if-test on nxpt

      call s2copy (nx+2, ny+2, xnrm, 1, nx+2, xnrmo, 1, nx+2)
      call s2copy (nx+2, ny+2, xvnrm, 1, nx+2, xvnrmo, 1, nx+2)
      call s2copy (nx+2, ny+2, ynrm, 1, nx+2, ynrmo, 1, nx+2)
      call s2copy (nx+2, ny+2, yvnrm, 1, nx+2, yvnrmo, 1, nx+2)

      do 705 ifld = 1, nisp
         if (nyomitmx >= nysol(1)+nyout(1)) then
c           # fill dead guard cells with adjacent values
            do iy = 0, ny+1
               ni(0,iy,ifld) = ni(1,iy,ifld)
               ni(nx+1,iy,ifld) = ni(nx,iy,ifld)
            enddo
         endif
         call s2copy (nx+2, ny+2, ni(0:nx+1,0:ny+1,ifld), 1, nx+2,
     .            nis(0:nx+1,0:ny+1,ifld), 1, nx+2)
  705 continue
      do 706 ifld = 1, nusp
         if (nyomitmx >= nysol(1)+nyout(1)) then
c           # fill dead guard cells with adjacent values
            do iy = 0, ny+1
               up(0,iy,ifld) = up(1,iy,ifld)
               up(nx+1,iy,ifld) = up(nx,iy,ifld)
            enddo
         endif
         call s2copy (nx+2, ny+2, up(0:nx+1,0:ny+1,ifld), 1, nx+2,
     .            ups(0:nx+1,0:ny+1,ifld), 1, nx+2)
  706 continue

      do 707 igsp = 1, ngsp
         if (nyomitmx >= nysol(1)+nyout(1)) then
c           # fill dead guard cells with adjacent values
            do iy = 0, ny+1
               ng(0,iy,igsp) = ng(1,iy,igsp)
               ng(nx+1,iy,igsp) = ng(nx,iy,igsp)
            enddo
         endif
         call s2copy (nx+2, ny+2, ng(0:nx+1,0:ny+1,igsp), 1, nx+2,
     .            ngs(0:nx+1,0:ny+1,igsp), 1, nx+2)
         call s2copy (nx+2, ny+2, tg(0:nx+1,0:ny+1,igsp), 1, nx+2,
     .            tgs(0:nx+1,0:ny+1,igsp), 1, nx+2)
  707 continue

         if (nyomitmx >=nysol(1)+nyout(1)) then
c           # fill dead guard cells with adjacent values
            do iy = 0, ny+1
               te(0,iy) = te(1,iy)
               te(nx+1,iy) = te(nx,iy)
               ti(0,iy) = ti(1,iy)
               ti(nx+1,iy) = ti(nx,iy)
               phi(0,iy) = phi(1,iy)
               phi(nx+1,iy) = phi(nx,iy)
            enddo
         endif
         call s2copy (nx+2, ny+2, te, 1, nx+2, tes, 1, nx+2)
         call s2copy (nx+2, ny+2, ti, 1, nx+2, tis, 1, nx+2)
         call s2copy (nx+2, ny+2, tg, 1, nx+2, tgs, 1, nx+2)
         call s2copy (nx+2, ny+2, phi, 1, nx+2, phis, 1, nx+2)
          if (isimpon.gt.0)
     .         call s2copy (nx+2, ny+2, afrac, 1, nx+2, afracs, 1, nx+2)

      return
      end
c****** end of subroutine gridseq ***
c************************************
c-----------------------------------------------------------------------
      subroutine refpla

*     REFPLA interpolates the plasma variables onto a grid which is
*     twice as fine as the preceding one

      implicit none

      Use(Share)    # geometry,nxc
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp
      Use(Compla)
      Use(Interp)   # nxold,nyold,nis,tes,tis,tgs,phis,ups,ngs
      Use(Selec)    # ixm1,ixp1


c...  local variables
      integer ifld
      integer ix2p, ix2m, ip1, im1, ixn, iy2, iy2p
      #Former Aux module variables
      integer ix,iy,igsp,ix2
c.... note: ix2 is already defined in Aux, but we can use it since it is
c.... meant to be a work variable, i.e., its value can change


c.... This is for doubling the grid in both x and y
         do 704 iy = 1, nyold-1
            iy2 = 2*iy
            iy2p = iy2 + 1
            do 703 ix = 1, nxold
               ix2 = 2*ix
               ix2p = ix2 + 1
               ix2m = max(ix2 - 1,0)
               ip1 = (ixp1(ix2,iy2) + 1) / 2
               im1 = ixm1(ix2m,iy2) / 2
               ixn = ixm1(ix2p,iy2) / 2
               do ifld = 1, nisp
               ni(ix2,iy2,ifld) = 9*nis(ix,iy,ifld)/16 +
     .                       3*(nis(ip1,iy,ifld)+nis(ix,iy+1,ifld))/16 +
     .                           nis(ip1,iy+1,ifld)/16
               ni(ix2p,iy2,ifld) = 9*nis(ix+1,iy,ifld)/16 +
     .                       3*(nis(ixn,iy,ifld)+nis(ix+1,iy+1,ifld))/16+
     .                           nis(ixn,iy+1,ifld)/16
               ni(ix2,iy2p,ifld) = 9*nis(ix,iy+1,ifld)/16 +
     .                       3*(nis(ix,iy,ifld)+nis(ip1,iy+1,ifld))/16 +
     .                             nis(ip1,iy,ifld)/16
               ni(ix2p,iy2p,ifld) = 9*nis(ix+1,iy+1,ifld)/16 +
     .                       3*(nis(ixn,iy+1,ifld)+nis(ix+1,iy,ifld))/16+
     .                             nis(ixn,iy,ifld)/16
               enddo

               do ifld = 1, nusp
               up(ix2,iy2,ifld) =  .75*ups(ix,iy,ifld) +
     .                             .25*ups(ix,iy+1,ifld)
               up(ix2,iy2p,ifld) = .75*ups(ix,iy+1,ifld) +
     .                             .25*ups(ix,iy,ifld)
               up(ix2m,iy2,ifld) = 3*(ups(ix,iy,ifld)+ups(im1,iy,ifld))/8
     .                         + (ups(ix,iy+1,ifld)+ups(im1,iy+1,ifld))/8
               up(ix2m,iy2p,ifld) = 3*(ups(ix,iy+1,ifld) +
     .                                             ups(im1,iy+1,ifld))/8
     .                            + (ups(ix,iy,ifld)+ups(im1,iy,ifld))/8
               enddo

               te(ix2,iy2) = 9*tes(ix,iy)/16 +
     .                        3*(tes(ip1,iy)+tes(ix,iy+1))/16 +
     .                           tes(ip1,iy+1)/16
               te(ix2p,iy2) = 9*tes(ix+1,iy)/16 +
     .                        3*(tes(ixn,iy)+tes(ix+1,iy+1))/16 +
     .                           tes(ixn,iy+1)/16
               te(ix2,iy2p) = 9*tes(ix,iy+1)/16 +
     .                        3*(tes(ix,iy)+tes(ip1,iy+1))/16 +
     .                             tes(ip1,iy)/16
               te(ix2p,iy2p) = 9*tes(ix+1,iy+1)/16 +
     .                        3*(tes(ixn,iy+1)+tes(ix+1,iy))/16 +
     .                             tes(ixn,iy)/16

               ti(ix2,iy2) = 9*tis(ix,iy)/16 +
     .                        3*(tis(ip1,iy)+tis(ix,iy+1))/16 +
     .                           tis(ip1,iy+1)/16
               ti(ix2p,iy2) = 9*tis(ix+1,iy)/16 +
     .                        3*(tis(ixn,iy)+tis(ix+1,iy+1))/16 +
     .                           tis(ixn,iy+1)/16
               ti(ix2,iy2p) = 9*tis(ix,iy+1)/16 +
     .                        3*(tis(ix,iy)+tis(ip1,iy+1))/16 +
     .                             tis(ip1,iy)/16
               ti(ix2p,iy2p) = 9*tis(ix+1,iy+1)/16 +
     .                        3*(tis(ixn,iy+1)+tis(ix+1,iy))/16 +
     .                             tis(ixn,iy)/16

               phi(ix2,iy2) = 9*phis(ix,iy)/16 +
     .                        3*(phis(ip1,iy)+phis(ix,iy+1))/16 +
     .                           phis(ip1,iy+1)/16
               phi(ix2p,iy2) = 9*phis(ix+1,iy)/16 +
     .                        3*(phis(ixn,iy)+phis(ix+1,iy+1))/16 +
     .                           phis(ixn,iy+1)/16
               phi(ix2,iy2p) = 9*phis(ix,iy+1)/16 +
     .                        3*(phis(ix,iy)+phis(ip1,iy+1))/16 +
     .                             phis(ip1,iy)/16
               phi(ix2p,iy2p) = 9*phis(ix+1,iy+1)/16 +
     .                        3*(phis(ixn,iy+1)+phis(ix+1,iy))/16 +
     .                             phis(ixn,iy)/16

            do igsp = 1, ngsp

               ng(ix2,iy2,igsp) = 9*ngs(ix,iy,igsp)/16 +
     .                        3*(ngs(ip1,iy,igsp)+ngs(ix,iy+1,igsp))/16 +
     .                           ngs(ip1,iy+1,igsp)/16
               ng(ix2p,iy2,igsp) = 9*ngs(ix+1,iy,igsp)/16 +
     .                       3*(ngs(ixn,iy,igsp)+ngs(ix+1,iy+1,igsp))/16+
     .                           ngs(ixn,iy+1,igsp)/16
               ng(ix2,iy2p,igsp) = 9*ngs(ix,iy+1,igsp)/16 +
     .                        3*(ngs(ix,iy,igsp)+ngs(ip1,iy+1,igsp))/16 +
     .                             ngs(ip1,iy,igsp)/16
               ng(ix2p,iy2p,igsp) = 9*ngs(ix+1,iy+1,igsp)/16 +
     .                       3*(ngs(ixn,iy+1,igsp)+ngs(ix+1,iy,igsp))/16+
     .                             ngs(ixn,iy,igsp)/16
            enddo

            do igsp = 1, ngsp

               tg(ix2,iy2,igsp) = 9*tgs(ix,iy,igsp)/16 +
     .                        3*(tgs(ip1,iy,igsp)+tgs(ix,iy+1,igsp))/16 +
     .                           tgs(ip1,iy+1,igsp)/16
               tg(ix2p,iy2,igsp) = 9*tgs(ix+1,iy,igsp)/16 +
     .                       3*(tgs(ixn,iy,igsp)+tgs(ix+1,iy+1,igsp))/16+
     .                           tgs(ixn,iy+1,igsp)/16
               tg(ix2,iy2p,igsp) = 9*tgs(ix,iy+1,igsp)/16 +
     .                        3*(tgs(ix,iy,igsp)+tgs(ip1,iy+1,igsp))/16 +
     .                             tgs(ip1,iy,igsp)/16
               tg(ix2p,iy2p,igsp) = 9*tgs(ix+1,iy+1,igsp)/16 +
     .                       3*(tgs(ixn,iy+1,igsp)+tgs(ix+1,iy,igsp))/16+
     .                             tgs(ixn,iy,igsp)/16
            enddo


  703       continue

            do ifld = 1, nisp
            ni(0,iy2,ifld) = .75*nis(0,iy,ifld)+.25*nis(0,iy+1,ifld)
            ni(1,iy2,ifld) = 3*(nis(0,iy  ,ifld)+nis(1,iy  ,ifld))/8 +
     .                         (nis(0,iy+1,ifld)+nis(1,iy+1,ifld))/8
            ni(nx+1,iy2,ifld) = .75*nis(nxold+1,iy,ifld)+
     .                          .25*nis(nxold+1,iy+1,ifld)
            ni(nx,iy2,ifld) = 3*(nis(nxold,iy  ,ifld)+
     .                                    nis(nxold+1,iy  ,ifld))/8 +
     .                          (nis(nxold,iy+1,ifld)+
     .                                    nis(nxold+1,iy+1,ifld))/8
            enddo
            do ifld = 1, nusp
            up(0,iy2,ifld) = .75*ups(0,iy,ifld)+.25*ups(0,iy+1,ifld)
            up(1,iy2,ifld) = 3*(ups(0,iy  ,ifld)+ups(1,iy  ,ifld))/8 +
     .                         (ups(0,iy+1,ifld)+ups(1,iy+1,ifld))/8
            up(nx+1,iy2,ifld) = .75*ups(nxold+1,iy,ifld)+
     .                          .25*ups(nxold+1,iy+1,ifld)
            up(nx,iy2,ifld) = .75*ups(nxold+1,iy,ifld)+
     .                          .25*ups(nxold+1,iy+1,ifld)
            enddo
            te(0,iy2) = .75*tes(0,iy)+.25*tes(0,iy+1)
            te(1,iy2) = 3*(tes(0,iy  ) + tes(1,iy  ))/8 +
     .                    (tes(0,iy+1) + tes(1,iy+1))/8
            te(nx+1,iy2) = .75*tes(nxold+1,iy)+.25*tes(nxold+1,iy+1)
            te(nx,iy2) = 3*(tes(nxold,iy  ) + tes(nxold+1,iy  ))/8 +
     .                     (tes(nxold,iy+1) + tes(nxold+1,iy+1))/8
            ti(0,iy2) = .75* tis(0,iy)+.25*tis(0,iy+1)
            ti(1,iy2) = 3*(tis(0,iy  ) + tis(1,iy  ))/8 +
     .                    (tis(0,iy+1) + tis(1,iy+1))/8
            ti(nx+1,iy2) = .75*tis(nxold+1,iy)+.25*tis(nxold+1,iy+1)
            ti(nx,iy2) = 3*(tis(nxold,iy  ) + tis(nxold+1,iy  ))/8 +
     .                     (tis(nxold,iy+1) + tis(nxold+1,iy+1))/8
            phi(0,iy2) = .75* phis(0,iy)+.25*phis(0,iy+1)
            phi(1,iy2) = 3*(phis(0,iy  ) + phis(1,iy  ))/8 +
     .                    (phis(0,iy+1) + phis(1,iy+1))/8
            phi(nx+1,iy2) = .75*phis(nxold+1,iy)+.25*phis(nxold+1,iy+1)
            phi(nx,iy2) = 3*(phis(nxold,iy  ) + phis(nxold+1,iy  ))/8 +
     .                     (phis(nxold,iy+1) + phis(nxold+1,iy+1))/8

         do igsp = 1, ngsp
            ng(0,iy2,igsp) = .75*ngs(0,iy,igsp)+.25*ngs(0,iy+1,igsp)
            ng(1,iy2,igsp) = 3*(ngs(0,iy  ,igsp)+ngs(1,iy,  igsp))/8 +
     .                      (ngs(0,iy+1,igsp)+ngs(1,iy+1,igsp))/8
            ng(nx+1,iy2,igsp) = .75*ngs(nxold+1,iy  ,igsp)+
     .                          .25*ngs(nxold+1,iy+1,igsp)
            ng(nx,iy2,igsp) = 3*(ngs(nxold  ,iy  ,igsp)+
     .                           ngs(nxold+1,iy  ,igsp))/8 +
     .                  (ngs(nxold,iy+1,igsp)+ngs(nxold+1,iy+1,igsp))/8
         enddo

         do igsp = 1, ngsp
            tg(0,iy2,igsp) = .75*tgs(0,iy,igsp)+.25*tgs(0,iy+1,igsp)
            tg(1,iy2,igsp) = 3*(tgs(0,iy  ,igsp)+tgs(1,iy,  igsp))/8 +
     .                      (tgs(0,iy+1,igsp)+tgs(1,iy+1,igsp))/8
            tg(nx+1,iy2,igsp) = .75*tgs(nxold+1,iy  ,igsp)+
     .                          .25*tgs(nxold+1,iy+1,igsp)
            tg(nx,iy2,igsp) = 3*(tgs(nxold  ,iy  ,igsp)+
     .                           tgs(nxold+1,iy  ,igsp))/8 +
     .                  (tgs(nxold,iy+1,igsp)+tgs(nxold+1,iy+1,igsp))/8
         enddo

            do ifld = 1, nisp
            ni(0,iy2p,ifld) = .25*nis(0,iy,ifld)+.75*nis(0,iy+1,ifld)
            ni(1,iy2p,ifld) = (nis(0,iy  ,ifld)+nis(1,iy  ,ifld))/8 +
     .                      3*(nis(0,iy+1,ifld)+nis(1,iy+1,ifld))/8
            ni(nx+1,iy2p,ifld) = .25*nis(nxold+1,iy,ifld)+
     .                          .75*nis(nxold+1,iy+1,ifld)
            ni(nx,iy2p,ifld) = (nis(nxold,iy  ,ifld)+
     .                                    nis(nxold+1,iy  ,ifld))/8 +
     .                       3*(nis(nxold,iy+1,ifld)+
     .                                    nis(nxold+1,iy+1,ifld))/8
            enddo
            do ifld = 1, nusp
            up(0,iy2p,ifld) = .25*ups(0,iy,ifld)+.75*ups(0,iy+1,ifld)
            up(1,iy2p,ifld) = (ups(0,iy  ,ifld)+ups(1,iy  ,ifld))/8 +
     .                      3*(ups(0,iy+1,ifld)+ups(1,iy+1,ifld))/8
            up(nx+1,iy2p,ifld) = .25*ups(nxold+1,iy  ,ifld) +
     .                           .75*ups(nxold+1,iy+1,ifld)
            up(nx,iy2p,ifld) = .25*ups(nxold+1,iy  ,ifld) +
     .                           .75*ups(nxold+1,iy+1,ifld)
            enddo
            te(0,iy2p) = .25*tes(0,iy)+.75*tes(0,iy+1)
            te(1,iy2p) = (tes(0,iy  ) + tes(1,iy  ))/8 +
     .                 3*(tes(0,iy+1) + tes(1,iy+1))/8
            te(nx+1,iy2p) = .25*tes(nxold+1,iy)+.75*tes(nxold+1,iy+1)
            te(nx,iy2p) = (tes(nxold,iy  ) + tes(nxold+1,iy  ))/8 +
     .                  3*(tes(nxold,iy+1) + tes(nxold+1,iy+1))/8
            ti(0,iy2p) = .25* tis(0,iy)+.75*tis(0,iy+1)
            ti(1,iy2p) = (tis(0,iy  ) + tis(1,iy  ))/8 +
     .                 3*(tis(0,iy+1) + tis(1,iy+1))/8
            ti(nx+1,iy2p) = .25*tis(nxold+1,iy)+.75*tis(nxold+1,iy+1)
            ti(nx,iy2p) = (tis(nxold,iy  ) + tis(nxold+1,iy  ))/8 +
     .                  3*(tis(nxold,iy+1) + tis(nxold+1,iy+1))/8
            phi(0,iy2p) = .25* phis(0,iy)+.75*phis(0,iy+1)
            phi(1,iy2p) = (phis(0,iy  ) + phis(1,iy  ))/8 +
     .                 3*(phis(0,iy+1) + phis(1,iy+1))/8
            phi(nx+1,iy2p) = .25*phis(nxold+1,iy)+.75*phis(nxold+1,iy+1)
            phi(nx,iy2p) = (phis(nxold,iy  ) + phis(nxold+1,iy  ))/8 +
     .                  3*(phis(nxold,iy+1) + phis(nxold+1,iy+1))/8

         do igsp = 1, ngsp
            ng(0,iy2p,igsp) = .25*ngs(0,iy,igsp)+.75*ngs(0,iy+1,igsp)
            ng(1,iy2p,igsp) = (ngs(0,iy  ,igsp)+ngs(1,iy  ,igsp))/8 +
     .                   3*(ngs(0,iy+1,igsp)+ngs(1,iy+1,igsp))/8
            ng(nx+1,iy2p,igsp) = .25*ngs(nxold+1,iy  ,igsp)+
     .                           .75*ngs(nxold+1,iy+1,igsp)
            ng(nx,iy2p,igsp) = (ngs(nxold  ,iy  ,igsp)+
     .                          ngs(nxold+1,iy  ,igsp))/8 +
     .                       3*(ngs(nxold  ,iy+1,igsp)+
     .                          ngs(nxold+1,iy+1,igsp))/8
         enddo
         do igsp = 1, ngsp
            tg(0,iy2p,igsp) = .25*tgs(0,iy,igsp)+.75*tgs(0,iy+1,igsp)
            tg(1,iy2p,igsp) = (tgs(0,iy  ,igsp)+tgs(1,iy  ,igsp))/8 +
     .                   3*(tgs(0,iy+1,igsp)+tgs(1,iy+1,igsp))/8
            tg(nx+1,iy2p,igsp) = .25*tgs(nxold+1,iy  ,igsp)+
     .                           .75*tgs(nxold+1,iy+1,igsp)
            tg(nx,iy2p,igsp) = (tgs(nxold  ,iy  ,igsp)+
     .                          tgs(nxold+1,iy  ,igsp))/8 +
     .                       3*(tgs(nxold  ,iy+1,igsp)+
     .                          tgs(nxold+1,iy+1,igsp))/8
         enddo

  704    continue

c.... Do the iy=0,1 and iy=ny+,ny1 points
         do 706 ix = 1 , nxold-1
         ix2 = 2*ix
         ix2p = ix2 + 1
         ix2m = max(ix2 - 1,0)
         ip1 = (ixp1(ix2,0) + 1) / 2
         im1 = ixm1(ix2m,0) / 2
         ixn = ixm1(ix2p,0) / 2
c....    the special ip1, im1, ixn are not done for ny+1 as the cuts are
c....    assumed not to extend there
         do ifld = 1, nisp
         ni(ix2,0,ifld) = .75*nis(ix,0,ifld)+.25*nis(ip1,0,ifld)
         ni(ix2,1,ifld) = 3*(nis(ix ,0,ifld)+nis(ix ,1,ifld))/8 +
     .                      (nis(ip1,0,ifld)+nis(ip1,1,ifld))/8
         ni(ix2,ny+1,ifld) = .75*nis(ix,nyold+1,ifld) +
     .                       .25*nis(ix+1,nyold+1,ifld)
         ni(ix2,ny,ifld) = 3*(nis(ix  ,nyold+1,ifld)+
     .                                     nis(ix  ,nyold,ifld))/8+
     .                       (nis(ix+1,nyold+1,ifld)+
     .                                     nis(ix+1,nyold,ifld))/8
         enddo
         do ifld = 1, nusp
         up(ix2,0,ifld) = ups(ix,0,ifld)
         up(ix2,1,ifld) = 0.5*(ups(ix,0,ifld)+ups(ix,1,ifld))
         up(ix2,ny+1,ifld) = ups(ix,nyold+1,ifld)
         up(ix2,ny,ifld) = 0.5*(ups(ix,nyold+1,ifld)+
     .                            ups(ix,nyold  ,ifld))
         enddo
         te(ix2,0) = .75*tes(ix,0) + .25*tes(ip1,0)
         te(ix2,1) = 3*(tes(ix ,0)+tes(ix ,1))/8 +
     .                 (tes(ip1,0)+tes(ip1,1))/8
         te(ix2,ny+1) = .75*tes(ix,nyold+1) + .25*tes(ix+1,nyold+1)
         te(ix2,ny) = 3*(tes(ix  ,nyold+1)+tes(ix  ,nyold))/8 +
     .                  (tes(ix+1,nyold+1)+tes(ix+1,nyold))/8
         ti(ix2,0) = .75*tis(ix,0) + .25*tis(ip1,0)
         ti(ix2,1) = 3*(tis(ix ,0)+tis(ix ,1))/8 +
     .                 (tis(ip1,0)+tis(ip1,1))/8
         ti(ix2,ny+1) = .75*tis(ix,nyold+1) + .25*tis(ix+1,nyold+1)
         ti(ix2,ny) = 3*(tis(ix  ,nyold+1)+tis(ix  ,nyold))/8 +
     .                  (tis(ix+1,nyold+1)+tis(ix+1,nyold))/8
         phi(ix2,0) = .75*phis(ix,0) + .25*phis(ip1,0)
         phi(ix2,1) = 3*(phis(ix ,0)+phis(ix ,1))/8 +
     .                 (phis(ip1,0)+phis(ip1,1))/8
         phi(ix2,ny+1) = .75*phis(ix,nyold+1) + .25*phis(ix+1,nyold+1)
         phi(ix2,ny) = 3*(phis(ix  ,nyold+1)+phis(ix  ,nyold))/8 +
     .                  (phis(ix+1,nyold+1)+phis(ix+1,nyold))/8

      do igsp = 1, ngsp
         ng(ix2,0,igsp) = .75*ngs(ix,0,igsp)+.25*ngs(ip1,0,igsp)
         ng(ix2,1,igsp) = 3*(ngs(ix ,0,igsp)+ngs(ix ,1,igsp))/8 +
     .                   (ngs(ip1,0,igsp)+ngs(ip1,1,igsp))/8
         ng(ix2,ny+1,igsp) = .75*ngs(ix,nyold+1,igsp)+
     .                       .25*ngs(ix+1,nyold+1,igsp)
         ng(ix2,ny,igsp) = 3*(ngs(ix  ,nyold+1,igsp)+
     .                        ngs(ix  ,nyold  ,igsp))/8+
     .                       (ngs(ix+1,nyold+1,igsp)+
     .                        ngs(ix+1,nyold  ,igsp))/8
      enddo

      do igsp = 1, ngsp
         tg(ix2,0,igsp) = .75*tgs(ix,0,igsp)+.25*tgs(ip1,0,igsp)
         tg(ix2,1,igsp) = 3*(tgs(ix ,0,igsp)+tgs(ix ,1,igsp))/8 +
     .                   (tgs(ip1,0,igsp)+tgs(ip1,1,igsp))/8
         tg(ix2,ny+1,igsp) = .75*tgs(ix,nyold+1,igsp)+
     .                       .25*tgs(ix+1,nyold+1,igsp)
         tg(ix2,ny,igsp) = 3*(tgs(ix  ,nyold+1,igsp)+
     .                        tgs(ix  ,nyold  ,igsp))/8+
     .                       (tgs(ix+1,nyold+1,igsp)+
     .                        tgs(ix+1,nyold  ,igsp))/8
      enddo

         do ifld = 1, nisp
         ni(ix2p,0,ifld) = .25*nis(ixn,0,ifld)+.75*nis(ix+1,0,ifld)
         ni(ix2p,1,ifld) = (nis(ixn ,0,ifld)+nis(ixn ,1,ifld))/8 +
     .                   3*(nis(ix+1,0,ifld)+nis(ix+1,1,ifld))/8
         ni(ix2p,ny+1,ifld) = .25*nis(ixn,nyold+1,ifld) +
     .                        .75*nis(ix+1,nyold+1,ifld)
         ni(ix2p,ny,ifld) = (nis(ixn ,nyold+1,ifld)+
     .                                nis(ixn ,nyold,ifld))/8+
     .                    3*(nis(ix+1,nyold+1,ifld)+
     .                                nis(ix+1,nyold,ifld))/8
         enddo
         do ifld = 1, nusp
         up(ix2m,0,ifld) = .5*(ups(ix,0,ifld)+ups(im1,0,ifld))
         up(ix2m,1,ifld) = .25*(ups(ix,0,ifld)+ups(im1,0,ifld)+
     .                          ups(ix,1,ifld)+ups(im1,1,ifld))
         up(ix2m,ny+1,ifld) = .5*(ups(ix,nyold+1,ifld) +
     .                            ups(ix-1,nyold+1,ifld))
         up(ix2m,ny,ifld) = .25*(ups(ix  ,nyold+1,ifld)+
     .                           ups(ix  ,nyold  ,ifld)+
     .                           ups(ix-1,nyold+1,ifld)+
     .                           ups(ix-1,nyold  ,ifld))
         enddo
         te(ix2p,0) = .25*tes(ixn,0) + .75*tes(ix+1,0)
         te(ix2p,1) = (tes(ixn ,0)+tes(ixn ,1))/8 +
     .              3*(tes(ix+1,0)+tes(ix+1,1))/8
         te(ix2p,ny+1) = .25*tes(ixn,nyold+1) + .75*tes(ix+1,nyold+1)
         te(ix2p,ny) = (tes(ixn ,nyold+1)+tes(ixn ,nyold))/8 +
     .               3*(tes(ix+1,nyold+1)+tes(ix+1,nyold))/8
         ti(ix2p,0) = .25*tis(ixn,0) + .75*tis(ix+1,0)
         ti(ix2p,1) = (tis(ixn ,0)+tis(ixn ,1))/8 +
     .              3*(tis(ix+1,0)+tis(ix+1,1))/8
         ti(ix2p,ny+1) = .25*tis(ixn,nyold+1) + .75*tis(ix+1,nyold+1)
         ti(ix2p,ny) = (tis(ixn ,nyold+1)+tis(ixn ,nyold))/8 +
     .               3*(tis(ix+1,nyold+1)+tis(ix+1,nyold))/8
         phi(ix2p,0) = .25*phis(ixn,0) + .75*phis(ix+1,0)
         phi(ix2p,1) = (phis(ixn ,0)+phis(ixn ,1))/8 +
     .              3*(phis(ix+1,0)+phis(ix+1,1))/8
         phi(ix2p,ny+1) = .25*phis(ixn,nyold+1) + .75*phis(ix+1,nyold+1)
         phi(ix2p,ny) = (phis(ixn ,nyold+1)+phis(ixn ,nyold))/8 +
     .               3*(phis(ix+1,nyold+1)+phis(ix+1,nyold))/8

      do igsp = 1, ngsp
         ng(ix2p,0,igsp) = .25*ngs(ixn,0,igsp)+.75*ngs(ix+1,0,igsp)
         ng(ix2p,1,igsp) = (ngs(ixn ,0,igsp)+ngs(ixn ,1,igsp))/8 +
     .                3*(ngs(ix+1,0,igsp)+ngs(ix+1,1,igsp))/8
         ng(ix2p,ny+1,igsp) = .25*ngs(ixn ,nyold+1,igsp)+
     .                        .75*ngs(ix+1,nyold+1,igsp)
         ng(ix2p,ny,igsp) = (ngs(ixn ,nyold+1,igsp)+
     .                       ngs(ixn ,nyold  ,igsp))/8+
     .                    3*(ngs(ix+1,nyold+1,igsp)+
     .                       ngs(ix+1,nyold  ,igsp))/8
      enddo

      do igsp = 1, ngsp
         tg(ix2p,0,igsp) = .25*tgs(ixn,0,igsp)+.75*tgs(ix+1,0,igsp)
         tg(ix2p,1,igsp) = (tgs(ixn ,0,igsp)+tgs(ixn ,1,igsp))/8 +
     .                3*(tgs(ix+1,0,igsp)+tgs(ix+1,1,igsp))/8
         tg(ix2p,ny+1,igsp) = .25*tgs(ixn ,nyold+1,igsp)+
     .                        .75*tgs(ix+1,nyold+1,igsp)
         tg(ix2p,ny,igsp) = (tgs(ixn ,nyold+1,igsp)+
     .                       tgs(ixn ,nyold  ,igsp))/8+
     .                    3*(tgs(ix+1,nyold+1,igsp)+
     .                       tgs(ix+1,nyold  ,igsp))/8
      enddo

  706    continue

c.... finally, reset the corners

      do ifld = 1, nisp
      ni(0,0,ifld) = nis(0,0,ifld)
      ni(0,ny+1,ifld) = nis(0,nyold+1,ifld)
      ni(nx+1,0,ifld) = nis(nxold+1,0,ifld)
      ni(nx+1,ny+1,ifld) = nis(nxold+1,nyold+1,ifld)
      ni(1,1,ifld) = 0.25*(nis(0,0,ifld)+nis(1,0,ifld)+
     .                     nis(0,1,ifld)+nis(1,1,ifld))
      ni(1,ny,ifld) = 0.25*(nis(0,nyold+1,ifld)+nis(1,nyold+1,ifld)+
     .                      nis(0,nyold  ,ifld)+nis(1,nyold  ,ifld))
      ni(nx,1,ifld) = 0.25*(nis(nxold+1,0,ifld)+nis(nxold,0,ifld)+
     .                      nis(nxold+1,1,ifld)+nis(nxold,1,ifld))
      ni(nx,ny,ifld) = 0.25*(nis(nxold+1,nyold+1,ifld)+
     .                                         nis(nxold,nyold+1,ifld)+
     .                       nis(nxold+1,nyold  ,ifld)+
     .                                         nis(nxold,nyold  ,ifld))
      ni(1,0,ifld) = 0.5*(nis(0,0,ifld)+nis(1,0,ifld))
      ni(0,1,ifld) = 0.5*(nis(0,0,ifld)+nis(0,1,ifld))
      ni(nx,0,ifld) = 0.5*(nis(nxold+1,0,ifld)+nis(nxold,0,ifld))
      ni(nx+1,1,ifld) = 0.5*(nis(nxold+1,0,ifld)+nis(nxold+1,1,ifld))
      ni(0,ny,ifld) = 0.5*(nis(0,nyold+1,ifld)+nis(0,nyold,ifld))
      ni(1,ny+1,ifld) = 0.5*(nis(0,nyold+1,ifld)+nis(1,nyold+1,ifld))
      ni(nx,ny+1,ifld) = 0.5*(nis(nxold+1,nyold+1,ifld)+
     .                                       nis(nxold,nyold+1,ifld))
      ni(nx+1,ny,ifld) = 0.5*(nis(nxold+1,nyold+1,ifld)+
     .                                       nis(nxold+1,nyold,ifld))
      enddo
      do ifld = 1, nusp
      up(0,0,ifld) = ups(0,0,ifld)
      up(0,ny+1,ifld) = ups(0,nyold+1,ifld)
      up(nx+1,0,ifld) = ups(nxold+1,0,ifld)
      up(nx+1,ny+1,ifld) = ups(nxold+1,nyold+1,ifld)
      up(1,1,ifld) = 0.25*(ups(0,0,ifld)+ups(1,0,ifld)+
     .                     ups(0,1,ifld)+ups(1,1,ifld))
      up(1,ny,ifld) = 0.25*(ups(0,nyold+1,ifld)+ups(1,nyold+1,ifld)+
     .                      ups(0,nyold  ,ifld)+ups(1,nyold  ,ifld))
      up(nx,1,ifld) = .5*(ups(nxold,0,ifld)+ups(nxold,1,ifld))
      up(nx,ny,ifld) = .5*(ups(nxold,nyold+1,ifld)+ups(nxold,nyold,ifld))
      up(1,0,ifld) = 0.5*(ups(0,0,ifld)+ups(1,0,ifld))
      up(0,1,ifld) = 0.5*(ups(0,0,ifld)+ups(0,1,ifld))
      up(nx,0,ifld) = 0.5*(ups(nxold+1,0,ifld)+ups(nxold,0,ifld))
      up(nx+1,1,ifld) = 0.5*(ups(nxold+1,0,ifld)+ups(nxold+1,1,ifld))
      up(0,ny,ifld) = 0.5*(ups(0,nyold+1,ifld)+ups(0,nyold,ifld))
      up(1,ny+1,ifld) = 0.5*(ups(0,nyold+1,ifld)+ups(1,nyold+1,ifld))
      up(nx,ny+1,ifld) = 0.5*(ups(nxold+1,nyold+1,ifld)+
     .                                       ups(nxold,nyold+1,ifld))
      up(nx+1,ny,ifld) = 0.5*(ups(nxold+1,nyold+1,ifld)+
     .                                       ups(nxold+1,nyold,ifld))
      up(nx-1,0,ifld) = 0.5*(ups(nxold,0,ifld)+ups(nxold-1,0,ifld))
      up(nx-1,1,ifld) = 0.25*(ups(nxold,0,ifld)+ups(nxold-1,0,ifld)+
     .                        ups(nxold,1,ifld)+ups(nxold-1,1,ifld))
      up(nx-1,ny+1,ifld) = .5*(ups(nxold  ,nyold+1,ifld) +
     .                         ups(nxold-1,nyold+1,ifld))
      up(nx-1,ny,ifld) = .25*(ups(nxold  ,nyold+1,ifld) +
     .                        ups(nxold-1,nyold+1,ifld) +
     .                        ups(nxold  ,nyold  ,ifld) +
     .                        ups(nxold-1,nyold  ,ifld))
      enddo
      te(0,0) = tes(0,0)
      te(0,ny+1) = tes(0,nyold+1)
      te(nx+1,0) = tes(nxold+1,0)
      te(nx+1,ny+1) = tes(nxold+1,nyold+1)
      te(1,1) = 0.25*(tes(0,0)+tes(1,0)+tes(0,1)+tes(1,1))
      te(1,ny) = 0.25*(tes(0,nyold+1)+tes(1,nyold+1)+
     .                 tes(0,nyold  )+tes(1,nyold  ))
      te(nx,1) = 0.25*(tes(nxold+1,0)+tes(nxold,0)+
     .                 tes(nxold+1,1)+tes(nxold,1))
      te(nx,ny) = 0.25*(tes(nxold+1,nyold+1)+tes(nxold,nyold+1)+
     .                  tes(nxold+1,nyold  )+tes(nxold,nyold  ))
      te(1,0) = 0.5*(tes(0,0)+tes(1,0))
      te(0,1) = 0.5*(tes(0,0)+tes(0,1))
      te(nx,0) =   0.5*(tes(nxold+1,0)+tes(nxold,0))
      te(nx+1,1) = 0.5*(tes(nxold+1,0)+tes(nxold+1,1))
      te(0,ny) =   0.5*(tes(0,nyold+1)+tes(0,nyold))
      te(1,ny+1) =  0.5*(tes(0,nyold+1)+tes(1,nyold+1))
      te(nx,ny+1) = 0.5*(tes(nxold+1,nyold+1)+tes(nxold,nyold+1))
      te(nx+1,ny) = 0.5*(tes(nxold+1,nyold+1)+tes(nxold+1,nyold))
      ti(0,0) = tis(0,0)
      ti(0,ny+1) = tis(0,nyold+1)
      ti(nx+1,0) = tis(nxold+1,0)
      ti(nx+1,ny+1) = tis(nxold+1,nyold+1)
      ti(1,1) = 0.25*(tis(0,0)+tis(1,0)+tis(0,1)+tis(1,1))
      ti(1,ny) = 0.25*(tis(0,nyold+1)+tis(1,nyold+1)+
     .                 tis(0,nyold  )+tis(1,nyold  ))
      ti(nx,1) = 0.25*(tis(nxold+1,0)+tis(nxold,0)+
     .                 tis(nxold+1,1)+tis(nxold,1))
      ti(nx,ny) = 0.25*(tis(nxold+1,nyold+1)+tis(nxold,nyold+1)+
     .                  tis(nxold+1,nyold  )+tis(nxold,nyold  ))
      ti(1,0) = 0.5*(tis(0,0)+tis(1,0))
      ti(0,1) = 0.5*(tis(0,0)+tis(0,1))
      ti(nx,0) =   0.5*(tis(nxold+1,0)+tis(nxold,0))
      ti(nx+1,1) = 0.5*(tis(nxold+1,0)+tis(nxold+1,1))
      ti(0,ny) =   0.5*(tis(0,nyold+1)+tis(0,nyold))
      ti(1,ny+1) =  0.5*(tis(0,nyold+1)+tis(1,nyold+1))
      ti(nx,ny+1) = 0.5*(tis(nxold+1,nyold+1)+tis(nxold,nyold+1))
      ti(nx+1,ny) = 0.5*(tis(nxold+1,nyold+1)+tis(nxold+1,nyold))
      phi(0,0) = phis(0,0)
      phi(0,ny+1) = phis(0,nyold+1)
      phi(nx+1,0) = phis(nxold+1,0)
      phi(nx+1,ny+1) = phis(nxold+1,nyold+1)
      phi(1,1) = 0.25*(phis(0,0)+phis(1,0)+phis(0,1)+phis(1,1))
      phi(1,ny) = 0.25*(phis(0,nyold+1)+phis(1,nyold+1)+
     .                 phis(0,nyold  )+phis(1,nyold  ))
      phi(nx,1) = 0.25*(phis(nxold+1,0)+phis(nxold,0)+
     .                 phis(nxold+1,1)+phis(nxold,1))
      phi(nx,ny) = 0.25*(phis(nxold+1,nyold+1)+phis(nxold,nyold+1)+
     .                  phis(nxold+1,nyold  )+phis(nxold,nyold  ))
      phi(1,0) = 0.5*(phis(0,0)+phis(1,0))
      phi(0,1) = 0.5*(phis(0,0)+phis(0,1))
      phi(nx,0) =   0.5*(phis(nxold+1,0)+phis(nxold,0))
      phi(nx+1,1) = 0.5*(phis(nxold+1,0)+phis(nxold+1,1))
      phi(0,ny) =   0.5*(phis(0,nyold+1)+phis(0,nyold))
      phi(1,ny+1) =  0.5*(phis(0,nyold+1)+phis(1,nyold+1))
      phi(nx,ny+1) = 0.5*(phis(nxold+1,nyold+1)+phis(nxold,nyold+1))
      phi(nx+1,ny) = 0.5*(phis(nxold+1,nyold+1)+phis(nxold+1,nyold))

      do igsp = 1, ngsp
       ng(0,0,igsp) = ngs(0,0,igsp)
       ng(0,ny+1,igsp) = ngs(0,nyold+1,igsp)
       ng(nx+1,0,igsp) = ngs(nxold+1,0,igsp)
       ng(nx+1,ny+1,igsp) = ngs(nxold+1,nyold+1,igsp)
       ng(1,1,igsp) = 0.25*(ngs(0,0,igsp)+ngs(1,0,igsp)+
     .                  ngs(0,1,igsp)+ngs(1,1,igsp))
       ng(1,ny,igsp) = 0.25*(ngs(0,nyold+1,igsp)+ngs(1,nyold+1,igsp)+
     .                   ngs(0,nyold  ,igsp)+ngs(1,nyold  ,igsp))
       ng(nx,1,igsp) = 0.25*(ngs(nxold+1,0,igsp)+ngs(nxold,0,igsp)+
     .                   ngs(nxold+1,1,igsp)+ngs(nxold,1,igsp))
       ng(nx,ny,igsp) = 0.25*(ngs(nxold+1,nyold+1,igsp)+
     .                                   ngs(nxold,nyold+1,igsp)+
     .                    ngs(nxold+1,nyold  ,igsp)+
     .                                   ngs(nxold,nyold  ,igsp))
       ng(1,0,igsp) = 0.5*(ngs(0,0,igsp)+ngs(1,0,igsp))
       ng(0,1,igsp) = 0.5*(ngs(0,0,igsp)+ngs(0,1,igsp))
       ng(nx,0,igsp) = 0.5*(ngs(nxold+1,0,igsp)+ngs(nxold,0,igsp))
       ng(nx+1,1,igsp) = 0.5*(ngs(nxold+1,0,igsp)+ngs(nxold+1,1,igsp))
       ng(0,ny,igsp) = 0.5*(ngs(0,nyold+1,igsp)+ngs(0,nyold,igsp))
       ng(1,ny+1,igsp) = 0.5*(ngs(0,nyold+1,igsp)+ngs(1,nyold+1,igsp))
       ng(nx,ny+1,igsp) = 0.5*(ngs(nxold+1,nyold+1,igsp)+
     .                                       ngs(nxold,nyold+1,igsp))
       ng(nx+1,ny,igsp) = 0.5*(ngs(nxold+1,nyold+1,igsp)+
     .                                       ngs(nxold+1,nyold,igsp))
      enddo

      do igsp = 1, ngsp
        tg(0,0,igsp) = tgs(0,0,igsp)
        tg(0,ny+1,igsp) = tgs(0,nyold+1,igsp)
        tg(nx+1,0,igsp) = tgs(nxold+1,0,igsp)
        tg(nx+1,ny+1,igsp) = tgs(nxold+1,nyold+1,igsp)
        tg(1,1,igsp) = 0.25*(tgs(0,0,igsp)+tgs(1,0,igsp)+
     .                  tgs(0,1,igsp)+tgs(1,1,igsp))
        tg(1,ny,igsp) = 0.25*(tgs(0,nyold+1,igsp)+tgs(1,nyold+1,igsp)+
     .                   tgs(0,nyold  ,igsp)+tgs(1,nyold  ,igsp))
        tg(nx,1,igsp) = 0.25*(tgs(nxold+1,0,igsp)+tgs(nxold,0,igsp)+
     .                   tgs(nxold+1,1,igsp)+tgs(nxold,1,igsp))
        tg(nx,ny,igsp) = 0.25*(tgs(nxold+1,nyold+1,igsp)+
     .                                   tgs(nxold,nyold+1,igsp)+
     .                    tgs(nxold+1,nyold  ,igsp)+
     .                                   tgs(nxold,nyold  ,igsp))
        tg(1,0,igsp) = 0.5*(tgs(0,0,igsp)+tgs(1,0,igsp))
        tg(0,1,igsp) = 0.5*(tgs(0,0,igsp)+tgs(0,1,igsp))
        tg(nx,0,igsp) = 0.5*(tgs(nxold+1,0,igsp)+tgs(nxold,0,igsp))
        tg(nx+1,1,igsp) = 0.5*(tgs(nxold+1,0,igsp)+tgs(nxold+1,1,igsp))
        tg(0,ny,igsp) = 0.5*(tgs(0,nyold+1,igsp)+tgs(0,nyold,igsp))
        tg(1,ny+1,igsp) = 0.5*(tgs(0,nyold+1,igsp)+tgs(1,nyold+1,igsp))
        tg(nx,ny+1,igsp) = 0.5*(tgs(nxold+1,nyold+1,igsp)+
     .                                       tgs(nxold,nyold+1,igsp))
        tg(nx+1,ny,igsp) = 0.5*(tgs(nxold+1,nyold+1,igsp)+
     .                                       tgs(nxold+1,nyold,igsp))
      enddo

c...  If this is a double-null case, we must adjust the values across
c...  the cut at nxc

      if (((geometry.eq.'dnbot').or.(geometry.eq.'dnXtarget'))
     &                                           .and. nxc.gt.0) then


      ix2 = nxc
      ix2p = ix2 + 1
      ix2m = max(ix2-1,0)
      do 910 iy = 0, ny+1
         do 902 ifld = 1, nisp
            ni(ix2 ,iy,ifld) = ni(ix2m  ,iy,ifld)
            ni(ix2p,iy,ifld) = ni(ix2p+1,iy,ifld)
 902     continue
         do 903 ifld = 1, nusp
            up(ix2m,iy,ifld) = 0.
            up(ix2 ,iy,ifld) = 0.
            up(ix2p,iy,ifld) = 0.
 903     continue

         te(ix2 ,iy) = te(ix2m  ,iy)
         te(ix2p,iy) = te(ix2p+1,iy)
         ti(ix2 ,iy) = ti(ix2m  ,iy)
         ti(ix2p,iy) = ti(ix2p+1,iy)
         phi(ix2 ,iy) = phi(ix2m  ,iy)
         phi(ix2p,iy) = phi(ix2p+1,iy)

         do 904 igsp = 1, ngsp
            ng(ix2,iy,igsp) = ng(ix2m,iy,igsp)
            ng(ix2p,iy,igsp) = ng(ix2p+1,iy,igsp)
 904     continue
         do igsp = 1, ngsp
            tg(ix2,iy,igsp) = tg(ix2m,iy,igsp)
            tg(ix2p,iy,igsp) = tg(ix2p+1,iy,igsp)
         enddo
 910  continue

      endif

      return
      end
c****** end of subroutine refpla ***
c***********************************
c -----------------------------------------------------------------------
      subroutine grdintpy (ixs,ixf,ixos,ixof,iys,iyf,iyos,iyof,nx,ny,
     .                     nxold,nyold,xn,yn,xno,yno,xny,yny,indxx,indxy)

c...  This subroutine calculates the "intermediate" mesh that has nxold
c...  poloidal points and ny radial points. Linear interpolation is used
c...  to find the intersections of the old mesh radial lines with the
c...  new mesh poloidal lines. This intermediate mesh is then used to
c...  interpolate variables in the radial direction to the new mesh, but at
c...  the old poloidal points. Subsequently, a poloidal interpolation it done.

c...  The input variables:
c...  ixs    -- initial poloidal index for new mesh
c...  ixf    -- final poloidal index for new mesh
c...  ixos   -- initial poloidal index for old mesh
c...  ixof   -- final poloidal index for old mesh
c...  iys    -- initial radial index for new mesh
c...  iyf    -- final radial index for new mesh
c...  iyos   -- initial radial index for old mesh
c...  iyof   -- final radial index for old mesh
c...  nx     -- total poloidal mesh size for new mesh
c...  ny     -- total radial mesh size for new mesh
c...  nxold  -- total poloidal mesh size for old mesh
c...  nyold  -- total radial mesh size for old mesh
c...  xn     -- normalized poloidal location for new mesh (nx,ny)
c...  yn     -- normalized radial location for new mesh (nx,ny)
c...  xno    -- normalized poloidal location for old mesh (nxold,nyold)
c...  yno    -- normalized radial location for old mesh (nxold,nyold)

c...  The output variables:
c...  xny    -- calculated poloidal location for intermediate mesh (nxold,ny)
c...  yny    -- calculated radial location for intermediate mesh (nxold,ny)
c...  indxx  -- poloidal index ixo used in constructing intermediate mesh
c...  indxy  -- radial index iy used in constructing intermediate mesh

      implicit none

c --  Input variables
      integer ixs,ixf,ixos,ixof,iys,iyf,iyos,iyof,nx,ny,nxold,nyold
      real xn(0:nx+1,0:ny+1), yn(0:nx+1,0:ny+1),
     .     xno(0:nxold+1,0:nyold+1), yno(0:nxold+1,0:nyold+1)

c --  Output variables
      real xny(0:nxold+1,0:ny+1), yny(0:nxold+1,0:ny+1)
      integer indxx(0:nxold+1,0:ny+1), indxy(0:nxold+1,0:ny+1)

c --  Local variables
ccc   integer ix,iyo
ccc   real d2
      integer iy,ixo,ixm,ixmp,iyom,iyomp,icount
      real spn,spo,xxn,yyn,d2min,smalln,almost1,delerr,ferr
      data smalln /1.e-07/, almost1 /9.9999e-01/, delerr /1.e-02/

c...  This routine searchs for the (ix,iyo) indice pair the gives a point
c...  on the new grid (xn,yn) that is closest to the old grid point
c...  (xno,yno) for a give (ixo,iy), i.e., the intermediate grid indices

      iyom = iyos
      do 90 iy = iys, iyf
         ixm = ixs
         do 80 ixo = ixos, ixof
            d2min = 1.e20

c...  Search for the closest point of old and new mesh for specific ixo,iy
c...  After the closest point, we must search for the points that straddle
cc            do 70 iyo = iyos, iyof
cc               do 60 ix = ixs, ixf
cc                  d2 = (xno(ixo,iyo) - xn(ix,iy))**2 +
cc     .                 (yno(ixo,iyo) - yn(ix,iy))**2
cc                  if (d2 .lt. d2min) then
cc                     d2min = d2
cc                     ixm = ix
cc                     iyom = iyo
cc                  endif
cc  60           continue
cc  70        continue
            if (ixm .eq. ixs) then
               ixmp = ixm + 1
            elseif (ixm .eq. ixf) then
               ixmp = ixm
               ixm = ixm - 1
            elseif (xn(ixm,iy) .ge. xno(ixo,iyom)) then
               ixmp = ixm
               ixm = ixm - 1
            else
               ixmp = ixm + 1
            endif
            if (iyom .eq. iyos) then
               iyomp = iyom + 1
            elseif (iyom .eq. iyof) then
               iyomp = iyom
               iyom = iyom - 1
            elseif (yno(ixo,iyom) .ge. yn(ixm,iy)) then
               iyomp = iyom
               iyom = iyom - 1
            else
               iyomp = iyom + 1
            endif

c...  Special case for xno = 0,1; yyn could be more general if xn.ne.0,1
            icount = 0
  72        continue
            icount = icount + 1
            if (icount .gt. 500) then
               call remark('***** grdinty cannot find straddling grid
     .points, check vel. grd at nx')
               write(*,*) 'ixo,iy,ixm,iyom = ',ixo,iy,ixm,iyom,
     .                    '  ixf,ixof,iyf,iyof = ',ixf,ixof,iyf,iyof
               call xerrab("")
            endif
            if (abs(xno(ixo,iyomp)-xno(ixo,iyom)) .lt. smalln) then
               xxn = xno(ixo,iyom)
               yyn = ( yn(ixm ,iy)*( xn(ixmp,iy  )-xno(ixo,iyom)) +
     .                 yn(ixmp,iy)*(xno(ixo ,iyom)- xn(ixm,iy  )) ) /
     .                            (xn(ixmp,iy)-xn(ixm,iy))
c...  Special if test prevents prob for vel. grid for nx pt.
               if (xxn.gt.almost1 .and. ixm.eq.ixf-1) yyn = yn(ixmp,iy)
            else
               spo = (yno(ixo,iyomp) - yno(ixo,iyom)) /
     .               (xno(ixo,iyomp) - xno(ixo,iyom))
               spn = (yn(ixmp,iy) - yn(ixm,iy)) /
     .               (xn(ixmp,iy) - xn(ixm,iy))
               xxn = (spn*xn(ixm,iy) - spo*xno(ixo,iyom) +
     .                yno(ixo,iyom) - yn(ixm,iy)) / (spn - spo +1e-200)
               yyn = yno(ixo,iyom) + spo*(xxn - xno(ixo,iyom))
            endif
c...  Verify that old grid pts nearly straddle new grd pt, or try again
            ferr = delerr*( yno(ixo,iyomp)-yno(ixo,iyom) )
            if (yyn.ge.yno(ixo,iyom)-ferr .and.
     .                              yyn.le.yno(ixo,iyomp)+ferr) goto 74
            if (yyn .lt. yno(ixo,iyom)) then
               if (iyom .eq. iyos) goto 74
               iyom = iyom - 1
               iyomp = iyom + 1
            else
               if (iyom .eq. iyof-1) goto 74
               iyom = iyom + 1
               iyomp = iyom + 1
            endif
            goto 72
  74        continue
            ferr = delerr*( xn(ixmp,iy)-xn(ixm,iy) )
            if (xxn.ge.xn(ixm,iy)-ferr .and.
     .                                xxn.le.xn(ixmp,iy)+ferr) goto 76
c...  Special if test prevents prob for vel. grid for nx pt.
            if (xxn.gt.almost1 .and. ixm.eq.ixf-1) goto 76
            if (xxn .lt. xn(ixm,iy)) then
               if (ixm .eq. ixs) goto 76
               ixm = ixm - 1
               ixmp = ixm + 1
            else
               if (ixm .eq. ixf-1) goto 76
               ixm = ixm + 1
               ixmp = ixm + 1
            endif
            goto 72
  76        continue

            xny(ixo,iy) = xxn
            yny(ixo,iy) = yyn
            indxx(ixo,iy) = ixm
            indxy(ixo,iy) = iyom

  80     continue
  90  continue

      return
      end
c **** end of subroutine grdintpy ****
c ************************************
c -----------------------------------------------------------------------
      subroutine intpvar (varo, varn, ivel, nnxold, nnyold)

c...  This subroutine does the combined radial and poloidal interpolation
c...  of the variable varo

      implicit none

c  -- Common blocks
      Use(Dim)              # nx,ny,nisp,ngsp(for arrays in Interp not used)]
      Use(Share)            # geometry,nxc,nyomitmx
      Use(Comgeo)           # xnrm,xvnrm
      Use(Xpoint_indices)   # ixlb,ixpt1,ixpt2,ixrb,iysptrx
      Use(Interp)           # ixlbo,ixpt1o,ixpt2o,ixrbo,iysptrxo,
                            # xnrmox,xvnrmox,ynrmox,yvnrmox,ynrmo,yvnrmo,
                            # xnrmnx,xvnrmnx,ynrmnx,yvnrmnx,wrkint,wrkint2

c  -- Input variables - vel. grid switch, mesh sizes, variable on old mesh
      integer ivel, nnxold, nnyold
      real varo(0:nnxold+1,0:nnyold+1)

c  -- Output variables - variable on new mesh
      real varn(0:nx+1,0:ny+1)

c  -- Local variables
      integer iy, ix2, ix2p, ix2m, ixxst, ixxsto, iyst, iysto, ixxend,
     .        ixxendo, iyend, iyendo, jx, ir

c...  Set double null indices
      ix2 = nxc
      ix2p = ix2 + 1
      ix2m = max(ix2-1,0)

c...  If ivel=0, use xnrm grid; if ivel=1, use xvnrm grid
c##########################################################
      if (ivel.eq.0) then    # Do the density mesh; very large ivel if-test
c##########################################################

c...  First we do the radial linear interpolation for two regions for
c...  first intermediate mesh (xnrmox,ynrmox)
         if (iysptrx .gt. 0) then
            iyend = iysptrx
            iyendo = iysptrxo
            if (nyomitmx >= nysol(1)+nyout(1)) then
               iyend = iysptrx+1
               iyendo = iysptrxo+1
            endif
            call radintp (0,iyend,0,iyendo,0,nnxold+1,
     .                 nx,ny,nnxold,nnyold,ynrmox,ynrmo,varo,wrkint)
         endif
         iyst = iysptrx + 1
         iysto = iysptrxo + 1
         if (iysptrx .eq. 0) iyst = 0
         if (iysptrxo .eq. 0) iysto = 0
         if (nyomitmx < nysol(1)+nyout(1)) then
            call radintp (iyst,ny+1,iysto,nnyold+1,0,nnxold+1,
     .                 nx,ny,nnxold,nnyold,ynrmox,ynrmo,varo,wrkint)
         endif

c...  Now do poloidal interpolation using radial result in wrkint from
c...  first intermediate mesh to second intermediate mesh (xnrmnx,ynrmnx)
      do ir = 1, 3*nxpt
         call polintp (ixst(ir),ixend(ir),ixsto(ir),ixendo(ir),0,ny+1,
     .                 nx,ny,nnxold,nnyold,xnrmnx,xnrmox,wrkint,wrkint2)
      enddo

c...  Next do radial interpolation using wrkint2 on second intermediate
c...  mesh to the final mesh (xnrm,ynrm); result is varn
         if (iysptrx .gt. 0) then
            iyend = iysptrx
            if (nyomitmx >= nysol(1)+nyout(1)) then
               iyend = iysptrx+1
            endif
            call radintp (0,iyend,0,iyend,0,nx+1,
     .                 nx,ny,nx,ny,ynrm,ynrmnx,wrkint2,varn)
         endif
         iyst = iysptrx + 1
         if (iysptrx .eq. 0) iyst = 0
         if (nyomitmx < nysol(1)+nyout(1)) then
            call radintp (iyst,ny+1,iyst,ny+1,0,nx+1,
     .                    nx,ny,nx,ny,ynrm,ynrmnx,wrkint2,varn)
         endif

c...  Last, fix midway boundaries for double null case
         if (((geometry.eq.'dnbot').or.(geometry.eq.'dnXtarget'))
     &                                           .and. nxc.gt.1) then
            do 20 iy = 0, ny+1
               varn(ix2 ,iy) = varn(ix2m  ,iy)
               varn(ix2p,iy) = varn(ix2p+1,iy)
  20        continue
         endif
#     End of the density mesh calculation

c...  This next case is for ivel=1, i.e., the parallel velocity equation
c##################################################
	    else   # now if ivel=1, branch of large ivel if-test
c##################################################
c...   First do the radial interpolation
         if (iysptrx .gt. 0) then
            iyend = iysptrx
            iyendo = iysptrxo
            if (nyomitmx >= nysol(1)+nyout(1)) then
               iyend = iysptrx+1
               iyendo = iysptrxo+1
            endif
            call radintp (0,iyend,0,iyendo,0,nnxold+1,
     .                 nx,ny,nnxold,nnyold,yvnrmox,yvnrmo,varo,wrkint)
         endif
         iyst = iysptrx + 1
         iysto = iysptrxo + 1
         if (iysptrx .eq. 0) iyst = 0
         if (iysptrxo .eq. 0) iysto = 0
         if (nyomitmx < nysol(1)+nyout(1)) then
            call radintp (iyst,ny+1,iysto,nnyold+1,0,nnxold+1,
     .                 nx,ny,nnxold,nnyold,yvnrmox,yvnrmo,varo,wrkint)
         endif

c...  poloidal interpolation for velocity grid
      do ir = 1, 3*nxpt
         call polintp (ixst(ir),ixend(ir),ixsto(ir),ixendo(ir),0,ny+1,
     .               nx,ny,nnxold,nnyold,xvnrmnx,xvnrmox,wrkint,wrkint2)
      enddo

c...  Next do radial interpolation using wrkint2 on second intermediate
c...  mesh to the final mesh (xvnrm,yvnrm); result is varn
         if (iysptrx .gt. 0) then
            iyend = iysptrx
            if (nyomitmx >= nysol(1)+nyout(1)) then
               iyend = iysptrx+1
            endif
            call radintp (0,iyend,0,iyend,0,nx+1,
     .                 nx,ny,nx,ny,yvnrm,yvnrmnx,wrkint2,varn)
         endif
         iyst = iysptrx + 1
         if (iysptrx .eq. 0) iyst = 0
         if (nyomitmx < nysol(1)+nyout(1)) then
            call radintp (iyst,ny+1,iyst,ny+1,0,nx+1,
     .                 nx,ny,nx,ny,yvnrm,yvnrmnx,wrkint2,varn)
         endif

c...  Last, fix midway boundaries for double null case
         if (((geometry.eq.'dnbot').or.(geometry.eq.'dnXtarget'))
     &                                           .and. nxc.gt.0) then
            do 40 iy = 0, ny+1
               varn(ix2m,iy) = 0.
               varn(ix2 ,iy) = 0.
               varn(ix2p,iy) = 0.
  40        continue
         endif

c###############################################
      endif   # end of large ivel if-test
c###############################################
      return
      end
c ***** end of subroutine intpvar ****
c ************************************
c ----------------------------------------------------------------------
      subroutine radintp (iys,iyf,iyos,iyof,ixs,ixf,nx,ny,nxold,nyold,
     .                                             yn,yo,varo,varn)

cccMER Note that argument nx is not used in subroutine radintp

c...  This subroutine does a linear interpolation in the radial direction
c...  but with the old poloidal mesh structure

      implicit none

c --  Input variables
      integer iys,iyf,iyos,iyof,ixs,ixf,nx,ny,nxold,nyold
      real yn(0:nxold+1,0:ny+1), yo(0:nxold+1,0:nyold+1),
     .     varo(0:nxold+1,0:nyold+1)

c --  Output variables
      real varn(0:nxold+1,0:ny+1)

c --  Local variables
      integer ix, iy, iyo, iyl
      real chng, avn, avo
      data chng/1.7/

      do 40 iy = iys, iyf
         do 30 ix = ixs, ixf
            iyl = iyos
            do 20 iyo = iyos, iyof
               if(yo(ix,iyo).gt.yn(ix,iy) .or. iyo.eq.iyof) goto 25
               iyl = iyo
  20        continue
  25        continue
            varn(ix,iy) = (varo(ix,iyl  )*(yo(ix,iyl+1)-yn(ix,iy )) +
     .                     varo(ix,iyl+1)*(yn(ix,iy   )-yo(ix,iyl))) /
     .                                    (yo(ix,iyl+1)-yo(ix,iyl))
c ...    check for extrapolation and limit change to chng
            if(yn(ix,iy) .lt. yo(ix,iyl)) then
               avn = abs(varn(ix,iy))
               avo = abs(varo(ix,iyl))
               if( avn .lt. avo ) then
                  avn = max(avn, avo/chng)
               else
                  avn = min(avn, avo*chng)
               endif
               varn(ix,iy) = sign (avn, varo(ix,iyl))
            endif
            if(yn(ix,iy) .gt. yo(ix,iyl+1)) then
               avn = abs(varn(ix,iy))
               avo = abs(varo(ix,iyl+1))
               if( avn .lt. avo ) then
                  avn = max(avn, avo/chng)
               else
                  avn = min(avn, avo*chng)
               endif
               varn(ix,iy) = sign (avn, varo(ix,iyl+1))
            endif
  30     continue
  40  continue

      return
      end
c **** end of subroutine radintp ****
c **********************************
c -------------------------------------------------------------------
      subroutine polintp (ixs,ixf,ixos,ixof,iys,iyf,nx,ny,nxold,nyold,
     .                                             xn,xo,varo,varn)

c...  This subroutine does the poloidal interpolation using the result
c...  of the radial interpolation on the new radial mesh

      implicit none

c --  Input variables
      integer ixs,ixf,ixos,ixof,iys,iyf,nx,ny,nxold,nyold
      real xn(0:nx+1,0:ny+1), xo(0:nxold+1,0:ny+1),
     .     varo(0:nxold+1,0:ny+1)

c --  Output variables
      real varn(0:nx+1,0:ny+1)

c --  Local variables
      integer ix, iy, ixo, ixl
      real chng, avn, avo
      data chng/1.7/

      do 40 ix = ixs, ixf
         do 30 iy = iys, iyf
            ixl = ixos
            do 30 ixo = ixos, ixof
               if(xo(ixo,iy).gt.xn(ix,iy) .or. ixo.eq.ixof) goto 25
               ixl = ixo
  20        continue
  25        continue
            varn(ix,iy) = (varo(ixl  ,iy)*(xo(ixl+1,iy)-xn(ix ,iy)) +
     .                     varo(ixl+1,iy)*(xn(ix   ,iy)-xo(ixl,iy))) /
     .                                    (xo(ixl+1,iy)-xo(ixl,iy))
c ...    check for extrapolation and limit change to chng
            if(xn(ix,iy) .lt. xo(ixl,iy)) then
               avn = abs(varn(ix,iy))
               avo = abs(varo(ixl,iy))
               if( avn .lt. avo ) then
                  avn = max(avn, avo/chng)
               else
                  avn = min(avn, avo*chng)
               endif
               varn(ix,iy) = sign (avn, varo(ixl,iy))
            endif
            if(xn(ix,iy) .gt. xo(ixl+1,iy)) then
               avn = abs(varn(ix,iy))
               avo = abs(varo(ixl+1,iy))
               if( avn .lt. avo ) then
                  avn = max(avn, avo/chng)
               else
                  avn = min(avn, avo*chng)
               endif
               varn(ix,iy) = sign (avn, varo(ixl+1,iy))
            endif
  30     continue
  40  continue

      return
      end
c ***** end of subroutine polintp ****

c__________________________________________________________
      subroutine filldead_guardcells
      
c...  For core-only simulations, the poloidal guard-cells are not included
c...  in the simulation but must have finite values, which is done here

      Use(Dim)     # nx,ny,nisp,ngsp
      Use(Compla)  # ni,up,ng,te,ti,phi
      Use(Share)   # nyomitmx
      
c...  Local variables
      integer ifld, igsp, ix, iy
      do ifld = 1, nisp
         do iy = 0, ny+1
            ni(0,iy,ifld) = ni(1,iy,ifld)
            ni(nx+1,iy,ifld) = ni(nx,iy,ifld)
         enddo
      enddo
         
      do ifld = 1, nusp
         do iy = 0, ny+1
            up(0,iy,ifld) = up(1,iy,ifld)
            up(nx+1,iy,ifld) = up(nx,iy,ifld)
         enddo
      enddo

      do igsp = 1, ngsp
         do iy = 0, ny+1
            ng(0,iy,igsp) = ng(1,iy,igsp)
            ng(nx+1,iy,igsp) = ng(nx,iy,igsp)
         enddo
      enddo

      do iy = 0, ny+1
         te(0,iy) = te(1,iy)
         te(nx+1,iy) = te(nx,iy)
         ti(0,iy) = ti(1,iy)
         ti(nx+1,iy) = ti(nx,iy)
         phi(0,iy) = phi(1,iy)
         phi(nx+1,iy) = phi(nx,iy)
      enddo

      return
      end
