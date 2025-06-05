c     ------------------------------------------------------------------
      subroutine flxfin
      implicit none

Use(Share)	#nycore,nysol
Use(Dimflxgrd)	#jdim,noregs,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Aeqflxgrd)	#
Use(Dim)
Use(Inpf)

      integer nunit
#
# ********* Write the file to be read by the grd package *************

      data nunit /66/
      open (nunit, file='flx-grd', form='unformatted', status='unknown')

      write(nunit) jdim,noregs,npts,nycore,nysol
      write(nunit) jmin,jmax,jsptrx,jaxis
      write(nunit) npoint,xcurve,ycurve
      write(nunit) nxefit,nyefit
      write(nunit) xold,yold,fold
      write(nunit) fpol
      write(nunit) bcentr,rcentr,rmagx,zmagx,simagx,sibdry,
     &              rgrid1,xdim,zdim,zmid
      write(nunit) nlim
      write(nunit) xlim,ylim
      write(nunit) eshot,etime,rseps,zseps,
     &              rvsin,zvsin,rvsout,zvsout
      write(nunit) xlbnd,xubnd,ylbnd,yubnd
      write(nunit) runid
      write(nunit) geometry

      close (nunit)

ccc   call remark("***** subroutine flxfin completed")

      return
      end

