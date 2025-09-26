c-----------------------------------------------------------------------
      real function erl1 (te, ne, tau)
      implicit none
      real te, ne, tau
Use(Dim)
Use(Share)                # istabon
Use(Physical_constants)   # ev
Use(Data_input)
Use(Rtdata)
Use(Rtdegas)

c     local variables --
      integer je,jd,jr
      real zloge,zlogd,zlogr,rle,rld,fje,fjd,fjr,
     .     erl111,erl112,erl121,erl122,erl11,erl12

c     procedures --
      real rqa, rsa
      external rqa, rsa

c     Compute electron radiation loss rate per neutral H atom due to
c     "ionization" processes - D. Stotler's "coupling to the ground state"
c     te [J]          = electron temperature
c     ne [/m**3]      = electron density
c     erl1 [J/sec]    = radiation rate

c----------------------------------------------------------------------c
      if (istabon .le. 7) then   # various older models

         erl1 = (rqa(te,ne,0)-13.6*ev_aph*rsa(te,ne,0.,0))*ne

c----------------------------------------------------------------------c
      elseif (istabon .eq. 8 .or. istabon .eq. 9) then # linear interpolation

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     radiation rate --
         erl111=welms1(je,jd,jr)
         erl112=welms1(je,jd+1,jr)
         erl121=welms1(je+1,jd,jr)
         erl122=welms1(je+1,jd+1,jr)
         erl11=erl111+fjd*(erl112-erl111)
         erl12=erl121+fjd*(erl122-erl121)
         erl1 = erl11 + fje*(erl12-erl11)

c----------------------------------------------------------------------c
      elseif ((istabon>9 .and. istabon<14) .or. istabon .eq. 17) then 
c     log-log interp on Stotler DEGAS2 tables

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     radiation rate --
         erl111=log( welms1(je,jd,jr) )
         erl112=log( welms1(je,jd+1,jr) )
         erl121=log( welms1(je+1,jd,jr) )
         erl122=log( welms1(je+1,jd+1,jr) )
         erl11=erl111+fjd*(erl112-erl111)
         erl12=erl121+fjd*(erl122-erl121)
         erl1 = exp( erl11 + fje*(erl12-erl11) )

c----------------------------------------------------------------------c
      elseif (istabon>13 .and. istabon<16) then # spatially-dependent

         if (tau.le.taumin) then
            jr  = 1
            fjr = tau/taumin
         else
            zlogr = log10(tau/taumin)/deltau + 2.
            zlogr = min(zlogr,real(mpr-1))
            jr  = int(zlogr)
            fjr = zlogr - jr
         endif

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     radiation rate --
         if (istabon .eq. 14) then
            erl111=(1.-fjr)*welms1(je,jd,jr)     +
     .                 fjr *welms1(je,jd,jr+1)
            erl112=(1.-fjr)*welms1(je,jd+1,jr)   +
     .                 fjr *welms1(je,jd+1,jr+1)
            erl121=(1.-fjr)*welms1(je+1,jd,jr)   +
     .                 fjr *welms1(je+1,jd,jr+1)
            erl122=(1.-fjr)*welms1(je+1,jd+1,jr) +
     .                 fjr *welms1(je+1,jd+1,jr+1)
            erl11=erl111+fjd*(erl112-erl111)
            erl12=erl121+fjd*(erl122-erl121)
            erl1 = erl11 + fje*(erl12-erl11)
         elseif (istabon .eq. 15) then
            erl111=(1.-fjr)*log( welms1(je,jd,jr) )     +
     .                 fjr *log( welms1(je,jd,jr+1) )
            erl112=(1.-fjr)*log( welms1(je,jd+1,jr) )   +
     .                 fjr *log( welms1(je,jd+1,jr+1) )
            erl121=(1.-fjr)*log( welms1(je+1,jd,jr) )   +
     .                 fjr *log( welms1(je+1,jd,jr+1) )
            erl122=(1.-fjr)*log( welms1(je+1,jd+1,jr) ) +
     .                 fjr *log( welms1(je+1,jd+1,jr+1) )
            erl11=erl111+fjd*(erl112-erl111)
            erl12=erl121+fjd*(erl122-erl121)
            erl1 = exp( erl11 + fje*(erl12-erl11) )
         endif

c----------------------------------------------------------------------c
      else				# write error message
         call xerrab('function erl1 not defined for istabon > 17')
      endif

      return
      end

c-----------------------------------------------------------------------
      real function erl2 (te, ne, tau)
      implicit none
      real te, ne, tau
Use(Dim)
Use(Share)                # istabon
Use(Physical_constants)   # ev
Use(Data_input)
Use(Rtdata)
Use(Rtdegas)

c     local variables --
      integer je,jd,jr
      real zloge,zlogd,zlogr,rle,rld,fje,fjd,fjr,
     .     erl211,erl212,erl221,erl222,erl21,erl22

c     procedures --
      real rra
      external rra

c     Compute electron radiation loss rate per H ion due to
c     "recombination" processes - D. Stotler's "coupling to the continuum"
c     te [J]          = electron temperature
c     ne [/m**3]      = electron density
c     erl2 [J/sec]    = radiation rate

c----------------------------------------------------------------------c
      if (istabon .le. 7) then   # various older models

	 erl2 = (13.6*ev_aph+1.5*te)*ne*rra(te,ne,0.,1)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 8 .or. istabon .eq. 9) then # linear interpolation

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     radiation rate --
         erl211=welms2(je,jd,jr)
         erl212=welms2(je,jd+1,jr)
         erl221=welms2(je+1,jd,jr)
         erl222=welms2(je+1,jd+1,jr)
         erl21=erl211+fjd*(erl212-erl211)
         erl22=erl221+fjd*(erl222-erl221)
         erl2 = erl21 + fje*(erl22-erl21)

c----------------------------------------------------------------------c
      elseif ((istabon>9 .and. istabon<14) .or. istabon .eq. 17) then 
c     log-log interp on Stotler DEGAS2 tables

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     radiation rate --  now logarithm of rate
         erl211=log( welms2(je,jd,jr) )
         erl212=log( welms2(je,jd+1,jr) )
         erl221=log( welms2(je+1,jd,jr) )
         erl222=log( welms2(je+1,jd+1,jr) )
         erl21=erl211+fjd*(erl212-erl211)
         erl22=erl221+fjd*(erl222-erl221)
         erl2 = exp( erl21 + fje*(erl22-erl21) )

c----------------------------------------------------------------------c
      elseif (istabon>13 .and. istabon<16) then # spatially-dependent

         if (tau.le.taumin) then
            jr  = 1
            fjr = tau/taumin
         else
            zlogr = log10(tau/taumin)/deltau + 2.
            zlogr = min(zlogr,real(mpr-1))
            jr  = int(zlogr)
            fjr = zlogr - jr
         endif

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     radiation rate --
         if (istabon .eq. 14) then
            erl211=(1.-fjr)*welms2(je,jd,jr)     +
     .                 fjr *welms2(je,jd,jr+1)
            erl212=(1.-fjr)*welms2(je,jd+1,jr)   +
     .                 fjr *welms2(je,jd+1,jr+1)
            erl221=(1.-fjr)*welms2(je+1,jd,jr)   +
     .                 fjr *welms2(je+1,jd,jr+1)
            erl222=(1.-fjr)*welms2(je+1,jd+1,jr) +
     .                 fjr *welms2(je+1,jd+1,jr+1)
            erl21=erl211+fjd*(erl212-erl211)
            erl22=erl221+fjd*(erl222-erl221)
            erl2 = erl21 + fje*(erl22-erl21)
        elseif (istabon .eq. 15) then
            erl211=(1.-fjr)*log( welms2(je,jd,jr) )     +
     .                 fjr *log( welms2(je,jd,jr+1) )
            erl212=(1.-fjr)*log( welms2(je,jd+1,jr) )   +
     .                 fjr *log( welms2(je,jd+1,jr+1) )
            erl221=(1.-fjr)*log( welms2(je+1,jd,jr) )   +
     .                 fjr *log( welms2(je+1,jd,jr+1) )
            erl222=(1.-fjr)*log( welms2(je+1,jd+1,jr) ) +
     .                 fjr *log( welms2(je+1,jd+1,jr+1) )
            erl21=erl211+fjd*(erl212-erl211)
            erl22=erl221+fjd*(erl222-erl221)
            erl2 = exp( erl21 + fje*(erl22-erl21) )
        endif

c-----------------------------------------------------------------------
      else				# write error message
         call xerrab('function erl2 not defined for istabon > 17')
      endif

      return
      end

c-----------------------------------------------------------------------
      real function rcx (t0, n0, k)
      implicit none
      integer k
      real t0, n0
Use(Dim)
Use(Share)                # istabon
Use(Physical_constants)   # ev,m_prot
Use(Data_input)           # issgvcxc,sgvcxc
Use(Rtdata)
Use(Rtdegas)

c     local variables --
      real a
      integer ini,iti
      real rlni,rlti,fxni,fxti,a0,a1
      integer je,j0
      real zloge,rle,fje,rcx1,rcx2
      real kdum
      integer zn,za,zamax
      external mcrates
      real rcxcopy

c     Compute rate parameter for k--->k-1 charge exchange on neutral hydrogen
c     k               = initial charge state
c     t0 [J]          = effective hydrogen temperature (per AMU)
c     n0 [/m**3]      = density (not used)
c     rcx [m**3/sec]  = <sigma v>

c----------------------------------------------------------------------c
      if ((istabon .eq. 1) .or. (istabon .eq. 2)) then
c			use ADPAK/STRAHL table look-up

c     indices for interpolation --
         iti = 0
         ini = 0
c     compute abscissae --
         rlti = log(t0/ev_aph)
         rlni = max(htln(0),min(htln(htnn),log(n0)))
c     find iti --
 51      if (iti .lt. htnt-1) then
            if (htlt(iti+1) .le. rlti) then
               iti = iti + 1
               goto 51
            endif
         endif
 52      if (0 .lt. iti) then
            if (rlti .lt. htlt(iti)) then
               iti = iti - 1
               goto 52
            endif
         endif
c     find ini --
 53      if (ini .lt. htnn-1) then
            if (htln(ini+1) .le. rlni) then
               ini = ini + 1
               goto 53
            endif
         endif
 54      if (0 .lt. ini) then
            if (rlni .lt. htln(ini)) then
               ini = ini - 1
               goto 54
            endif
         endif
c     compute coefficients for linear interpolation --
         fxni = (rlni-htln(ini))/(htln(ini+1)-htln(ini))
         fxti = (rlti-htlt(iti))/(htlt(iti+1)-htlt(iti))
c     compute charge exchange rate parameter for k --> k-1 process --
         a0 = (1-fxni)*htlcx(iti,ini,k) + fxni*htlcx(iti,ini+1,k)
         a1 = (1-fxni)*htlcx(iti+1,ini,k) + fxni*htlcx(iti+1,ini+1,k)
         rcx = exp((1-fxti)*a0+fxti*a1)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 3) then
c			use DEGAS table look-up

c     compute abscissa --
         zloge=log(t0/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
c     table index for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
c     fractional part of interval (je,je+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
c     charge exchange rate --
         j0=je		# assume neutral temperature is same as ions
         rcx1=svphcx(je,j0)
         rcx2=svphcx(je+1,j0)
         rcx = rcx1 + fje*(rcx2-rcx1)

c----------------------------------------------------------------------c
      elseif (istabon==16) then  # For data in b2frates format (B2,SOLPS)
         zn=1     # nuclear charge for hydrogenic species
         zamax=1  # maximum atomic charge for hydrogenic species
         za=1     # compute c-x rate for this charge state
         call mcrates(n0,t0,t0,za,zamax,zn,kdum,kdum,rcxcopy)
         rcx=rcxcopy

c----------------------------------------------------------------------c
      else     #     use analytic model (hydrogen) for all other istabon

         a = 3*t0 / (10*ev_aph)
         rcx = 1.7e-14 * a**0.333
         if (issgvcxc.eq.1) rcx = sgvcxc # use fixed sig-v 
         if (issgvcxc.eq.2) rcx = sgvcxc*sqrt(t0/m_prot) # fixed sig

c----------------------------------------------------------------------c
      endif

      return
      end
c-----------------------------------------------------------------------
      real function rqa (te, ne, k)
      implicit none
      integer k
      real te, ne
Use(Dim)
Use(Share)                # istabon
Use(Physical_constants)   # ev
Use(Data_input)
Use(Ionization_energy)    # erad
Use(Rtdata)
Use(Rtdegas)
Use(Aphwrk)
      Use(Timespl)
      real(Size4) sec4, gettime, tsval

c     external procedures --
      real rsa, svradp, B2VAhL
      external rsa, svradp, B2VAhL, gettime

c     local variables --
      real a
      integer ine,ite
      real rlne,rlte,fxne,fxte,t0,t1
      integer je,jd,jr
      real zloge,zlogd,rle,rld,fje,fjd,w11,w12,w21,w22,w1,w2,w,vlogw
      integer nxcoef,nycoef
      real xuse,yuse

c     Compute electron energy loss rate parameter for processes
c     starting from charge state k --
c     k                = initial charge state
c     te [J]           = electron temperature
c     ne [/m**3]       = electron density
c     rqa [J*m**3/sec] = <sigma*v>*dE where dE is electron energy loss

c----------------------------------------------------------------------c
      if (istabon .eq. 0) then   # use analytic model (hydrogen) with
c                                  # constant energy loss per ionization
         a = te / (10*ev_aph)
         rqa = erad * ev_aph * 3.0e-14 * a*a / (3.0 + a*a)

c----------------------------------------------------------------------c
      elseif ((istabon .eq. 1) .or. (istabon .eq. 2)) then
c			use ADPAK table look-up

c     indices for interpolation --
         ite = 0
         ine = 0
c     compute abscissae --
         rlte = log(te/ev_aph)
         rlne = max(htln(0),min(htln(htnn),log(ne)))
c     find ite --
 51      if (ite .lt. htnt-1) then
            if (htlt(ite+1) .le. rlte) then
               ite = ite + 1
               goto 51
            endif
         endif
 52      if (0 .lt. ite) then
            if (rlte .lt. htlt(ite)) then
               ite = ite - 1
               goto 52
            endif
         endif
c     find ine --
 53      if (ine .lt. htnn-1) then
            if (htln(ine+1) .le. rlne) then
               ine = ine + 1
               goto 53
            endif
         endif
 54      if (0 .lt. ine) then
            if (rlne .lt. htln(ine)) then
               ine = ine - 1
               goto 54
            endif
         endif
c     compute coefficients for linear interpolation --
         fxne = (rlne-htln(ine))/(htln(ine+1)-htln(ine))
         fxte = (rlte-htlt(ite))/(htlt(ite+1)-htlt(ite))
c     compute electron energy loss rate parameter --
         t0 = (1-fxne)*htlqa(ite,ine,k) + fxne*htlqa(ite,ine+1,k)
         t1 = (1-fxne)*htlqa(ite+1,ine,k) + fxne*htlqa(ite+1,ine+1,k)
         rqa = ev_aph*exp((1-fxte)*t0+fxte*t1)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 3) then
c			use DEGAS table look-up

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     electron energy loss per ionization --
         w11=welms(je,jd)
         w12=welms(je,jd+1)
         w21=welms(je+1,jd)
         w22=welms(je+1,jd+1)
         w1=w11+fjd*(w12-w11)
         w2=w21+fjd*(w22-w21)
         w = w1 + fje*(w2-w1)
c     electron energy loss rate parameter --
         rqa = ev_aph * w * rsa(te,ne,0.,k)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 4) then
c			use POST93 table look-up

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     hydrogen line radiation rate --
         w11=wlemiss(je,jd)
         w12=wlemiss(je,jd+1)
         w21=wlemiss(je+1,jd)
         w22=wlemiss(je+1,jd+1)
         w1=w11+fjd*(w12-w11)
         w2=w21+fjd*(w22-w21)
         w = w1 + fje*(w2-w1)
c     electron energy loss rate parameter --
         rqa = w + 13.6 * ev_aph * rsa(te,ne,0.,k)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 5) then
c			use spline fit to POST93 table data
c      xuse=min(max(xdata_aph(1),log(te/ev_aph)),xdata_aph(nxdata_aph))
c      yuse=min(max(ydata_aph(1),log10(ne)),ydata_aph(nydata_aph))

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))

      xuse=rle
      yuse=rld

      nxcoef=nxdata_aph
      nycoef=nydata_aph

      vlogw = B2VAhL(xuse, yuse, 0, 0, xknots_aph, yknots_aph, nxcoef,
     .          nycoef, kxords_aph, kyords_aph, rqacoef, ldf_aph, workh, iflag_aph)
      w=10**vlogw
c     electron energy loss rate parameter --
         rqa = w + 13.6 * ev_aph * rsa(te,ne,0.,k)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 6) then
c			use spline fit to POST93 table data
c      xuse=min(max(xdata_aph(1),log(te/ev_aph)),xdata_aph(nxdata_aph))
c      yuse=min(max(ydata_aph(1),log10(ne)),ydata_aph(nydata_aph))

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))

      xuse=rle
      yuse=rld

      nxcoef=nxdata_aph
      nycoef=nydata_aph

      tsval = gettime(sec4)
      w = B2VAhL(xuse, yuse, 0, 0, xknots_aph, yknots_aph, nxcoef, nycoef,
     .          kxords_aph, kyords_aph, rqacoef, ldf_aph, workh, iflag_aph)
      totb2val = totb2val + gettime(sec4) - tsval
c     electron energy loss rate parameter --
         rqa = w + 13.6 * ev_aph * rsa(te,ne,0.,k)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 7) then
c                       use polynomial fit from Bob Campbell -  8/93
c     Note that the 13.6 * ev_aph * rsa(te,ne,k) is omitted here as Campbell
c     has already added it in
         rqa = svradp(te/ev_aph,ne)

c----------------------------------------------------------------------c
      elseif (istabon .gt. 7) then	# write error message 
         call xerrab('**** function rqa is not defined for istabon > 7')

c----------------------------------------------------------------------c
      endif

      return
      end
c-----------------------------------------------------------------------
      real function rra (te, ne, tau, k)
      implicit none
      integer k
      real te, ne, tau
Use(Dim)
Use(Share)                # istabon
Use(Physical_constants)   # ev
Use(Data_input)
Use(Rtdata)
Use(Rtdegas)
Use(Aphwrk)
      Use(Timespl)
      real(Size4) sec4, gettime, tsval

c     external procedures --
      real srecf, B2VAhL
      external srecf, B2VAhL, gettime

c     local variables --
      integer ine,ite
      real rlne,rlte,fxne,fxte,t0,t1
      integer je,jd,jr
      real zloge,zlogd,zlogr,rle,rld,fje,fjd,fjr,
     .     rra11,rra12,rra21,rra22,rra1,rra2
      integer nxcoef,nycoef
      real xuse,yuse,vlog10rra
      real kdum
      integer zn,za,zamax
      real rracopy

c     Compute rate parameter for k--->k-1 recombination
c     k               = initial charge state
c     te [J]          = electron temperature
c     ne [/m**3]      = electron density
c     rra [m**3/sec]  = <sigma v>

c----------------------------------------------------------------------c
      if (istabon .eq. 0) then   # use analytic model
         rra = 0.

c----------------------------------------------------------------------c
      elseif ((istabon .eq. 1) .or. (istabon .eq. 2)) then
c			use ADPAK table look-up

c     indices for interpolation --
         ite = 0
         ine = 0
c     compute abscissae --
         rlte = log(te/ev_aph)
         rlne = max(htln(0),min(htln(htnn),log(ne)))
c     find ite --
 51      if (ite .lt. htnt-1) then
            if (htlt(ite+1) .le. rlte) then
               ite = ite + 1
               goto 51
            endif
         endif
 52      if (0 .lt. ite) then
            if (rlte .lt. htlt(ite)) then
               ite = ite - 1
               goto 52
            endif
         endif
c     find ine --
 53      if (ine .lt. htnn-1) then
            if (htln(ine+1) .le. rlne) then
               ine = ine + 1
               goto 53
            endif
         endif
 54      if (0 .lt. ine) then
            if (rlne .lt. htln(ine)) then
               ine = ine - 1
               goto 54
            endif
         endif
c     compute coefficients for linear interpolation --
         fxne = (rlne-htln(ine))/(htln(ine+1)-htln(ine))
         fxte = (rlte-htlt(ite))/(htlt(ite+1)-htlt(ite))
c     compute recombination rate parameter for k --> k-1 process  --
         t0 = (1-fxne)*htlra(ite,ine,k) + fxne*htlra(ite,ine+1,k)
         t1 = (1-fxne)*htlra(ite+1,ine,k) + fxne*htlra(ite+1,ine+1,k)
         rra = exp((1-fxte)*t0+fxte*t1)

c----------------------------------------------------------------------c
      elseif ((istabon .eq. 3) .or. (istabon .eq. 4)
     .   .or. (istabon .eq. 8) .or. (istabon .eq. 9)) then
c             use DEGAS or POST93 or DEGAS93 or Stotler95 table look-up

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     recombination rate parameter --
         rra11=wsveh0(je,jd,jr)
         rra12=wsveh0(je,jd+1,jr)
         rra21=wsveh0(je+1,jd,jr)
         rra22=wsveh0(je+1,jd+1,jr)
         rra1=rra11+fjd*(rra12-rra11)
         rra2=rra21+fjd*(rra22-rra21)
         rra = rra1 + fje*(rra2-rra1)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 5) then
c			use spline fit to POST93 table data
c      xuse=min(max(xdata_aph(1),log(te/ev_aph)),xdata_aph(nxdata_aph))
c      yuse=min(max(ydata_aph(1),log10(ne)),ydata_aph(nydata_aph))

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))

      xuse=rle
      yuse=rld

      nxcoef=nxdata_aph
      nycoef=nydata_aph

      vlog10rra = B2VAhL(xuse, yuse, 0, 0, xknots_aph, yknots_aph, nxcoef,
     .          nycoef, kxords_aph, kyords_aph, rracoef, ldf_aph, workh, iflag_aph)
      rra=10**vlog10rra

c----------------------------------------------------------------------c
      elseif (istabon .eq. 6) then
c			use spline fit to POST93 table data
c      xuse=min(max(xdata_aph(1),log(te/ev_aph)),xdata_aph(nxdata_aph))
c      yuse=min(max(ydata_aph(1),log10(ne)),ydata_aph(nydata_aph))

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))

      xuse=rle
      yuse=rld

      nxcoef=nxdata_aph
      nycoef=nydata_aph

      tsval = gettime(sec4)
      rra = B2VAhL(xuse, yuse, 0, 0, xknots_aph, yknots_aph, nxcoef, nycoef,
     .          kxords_aph, kyords_aph, rracoef, ldf_aph, workh, iflag_aph)
      totb2val = totb2val + gettime(sec4) - tsval
c----------------------------------------------------------------------c
      elseif (istabon .eq. 7) then
c                       use polynomial fit from Bob Campbell -  8/93
         rra = srecf(te/ev_aph,ne)

c----------------------------------------------------------------------c
      elseif ((istabon>9 .and. istabon<14) .or. istabon .eq. 17) then 
c     log-log interp on Stotler DEGAS2 tables
         jr = 1
c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     recombination rate parameter --
         rra11=log( wsveh0(je,jd,jr) )
         rra12=log( wsveh0(je,jd+1,jr) )
         rra21=log( wsveh0(je+1,jd,jr) )
         rra22=log( wsveh0(je+1,jd+1,jr) )
         rra1=rra11+fjd*(rra12-rra11)
         rra2=rra21+fjd*(rra22-rra21)
         rra = exp( rra1 + fje*(rra2-rra1) )

c----------------------------------------------------------------------c
      elseif (istabon>13 .and. istabon<16) then # spatially-dep opt-depth

         if (tau.le.taumin) then
            jr  = 1
            fjr = tau/taumin
         else
            zlogr = log10(tau/taumin)/deltau + 2.
            zlogr = min(zlogr,real(mpr-1))
            jr  = int(zlogr)
            fjr = zlogr - jr
         endif

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     recombination rate parameter --
         if (istabon .eq. 14) then
            rra11=(1.-fjr)*wsveh0(je,jd,jr)     +
     .                fjr *wsveh0(je,jd,jr+1)
            rra12=(1.-fjr)*wsveh0(je,jd+1,jr)   +
     .                fjr *wsveh0(je,jd+1,jr+1)
            rra21=(1.-fjr)*wsveh0(je+1,jd,jr)   +
     .                fjr *wsveh0(je+1,jd,jr+1)
            rra22=(1.-fjr)*wsveh0(je+1,jd+1,jr) +
     .                fjr *wsveh0(je+1,jd+1,jr+1)
            rra1=rra11+fjd*(rra12-rra11)
            rra2=rra21+fjd*(rra22-rra21)
            rra =rra1 + fje*(rra2-rra1)
         elseif (istabon .eq. 15) then
            rra11=(1.-fjr)*log( wsveh0(je,jd,jr) )     +
     .                fjr *log( wsveh0(je,jd,jr+1) )
            rra12=(1.-fjr)*log( wsveh0(je,jd+1,jr) )   +
     .                fjr *log( wsveh0(je,jd+1,jr+1) )
            rra21=(1.-fjr)*log( wsveh0(je+1,jd,jr) )   +
     .                fjr *log( wsveh0(je+1,jd,jr+1) )
            rra22=(1.-fjr)*log( wsveh0(je+1,jd+1,jr) ) +
     .                fjr *log( wsveh0(je+1,jd+1,jr+1) )
            rra1=rra11+fjd*(rra12-rra11)
            rra2=rra21+fjd*(rra22-rra21)
            rra = exp( rra1 + fje*(rra2-rra1) )
         endif

c----------------------------------------------------------------------c
         elseif (istabon==16) then
            zn=1     # nuclear charge for hydrogenic species
            zamax=1  # maximum atomic charge for hydrogenic species
            za=1     # compute recombination rate for this charge state
            call mcrates(ne,te,te,za,zamax,zn,kdum,rracopy,kdum)
            rra=rracopy

c----------------------------------------------------------------------c
      endif

      return
      end
c-----------------------------------------------------------------------
      real function rsa (te, ne, tau, k)
      implicit none
      integer k
      real te, ne, tau
Use(Dim)
Use(Share)                # istabon
Use(Physical_constants)   # ev
Use(Data_input)
Use(Rtdata)
Use(Rtdegas)
Use(Aphwrk)
      Use(Timespl)
      real(Size4) sec4, gettime, tsval

c     external procedures --
      real sionf, B2VAhL
      external sionf, B2VAhL, gettime

c     local variables --
      real a
      integer ine,ite
      real rlne,rlte,fxne,fxte,t0,t1
      integer je,jd,jr
      real zloge,zlogd,zlogr,rle,rld,fje,fjd,fjr,
     .     rsa11,rsa12,rsa21,rsa22,rsa1,rsa2
      integer nxcoef,nycoef
      real xuse,yuse,vlog10rsa
      real kdum
      integer zn,za,zamax
      external mcrates
      real rsacopy

c     Compute rate parameter for k--->k+1 ionization by electrons
c     k               = initial charge state
c     te [J]          = electron temperature
c     ne [/m**3]      = electron density
c     rsa [m**3/sec]  = <sigma v>

c----------------------------------------------------------------------c
      if (istabon .eq. 0) then   # use analytic model (hydrogen)
         a = te / (10*ev_aph)
         rsa = 3.0e-14 * a*a / (3.0 + a*a)

c----------------------------------------------------------------------c
      elseif ((istabon .eq. 1) .or. (istabon .eq. 2)) then
c			use ADPAK table look-up

c     indices for interpolation --
         ite = 0
         ine = 0
c     compute abscissae --
         rlte = log(te/ev_aph)
         rlne = max(htln(0),min(htln(htnn),log(ne)))
c     find ite --
 51      if (ite .lt. htnt-1) then
            if (htlt(ite+1) .le. rlte) then
               ite = ite + 1
               goto 51
            endif
         endif
 52      if (0 .lt. ite) then
            if (rlte .lt. htlt(ite)) then
               ite = ite - 1
               goto 52
            endif
         endif
c     find ine --
 53      if (ine .lt. htnn-1) then
            if (htln(ine+1) .le. rlne) then
               ine = ine + 1
               goto 53
            endif
         endif
 54      if (0 .lt. ine) then
            if (rlne .lt. htln(ine)) then
               ine = ine - 1
               goto 54
            endif
         endif
c     compute coefficients for linear interpolation --
         fxne = (rlne-htln(ine))/(htln(ine+1)-htln(ine))
         fxte = (rlte-htlt(ite))/(htlt(ite+1)-htlt(ite))
c     compute ionization rate parameter for k --> k+1 process --
         t0 = (1-fxne)*htlsa(ite,ine,k) + fxne*htlsa(ite,ine+1,k)
         t1 = (1-fxne)*htlsa(ite+1,ine,k) + fxne*htlsa(ite+1,ine+1,k)
         rsa = exp((1-fxte)*t0+fxte*t1)

c----------------------------------------------------------------------c
      elseif ((istabon .eq. 3) .or. (istabon .eq. 4)
     .   .or. (istabon .eq. 8) .or. (istabon .eq. 9)) then
c       	use DEGAS or POST93 or DEGAS93 or Stotler95 table look-up

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     ionization rate parameter --
         rsa11=wsveh(je,jd,jr)
         rsa12=wsveh(je,jd+1,jr)
         rsa21=wsveh(je+1,jd,jr)
         rsa22=wsveh(je+1,jd+1,jr)
         rsa1=rsa11+fjd*(rsa12-rsa11)
         rsa2=rsa21+fjd*(rsa22-rsa21)
         rsa = rsa1 + fje*(rsa2-rsa1)

c----------------------------------------------------------------------c
      elseif (istabon .eq. 5) then
c			use spline fit to POST93 table data
c      xuse=min(max(xdata_aph(1),log(te/ev_aph)),xdata_aph(nxdata_aph))
c      yuse=min(max(ydata_aph(1),log10(ne)),ydata_aph(nydata_aph))

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))

      xuse=rle
      yuse=rld

      nxcoef=nxdata_aph
      nycoef=nydata_aph

      vlog10rsa = B2VAhL(xuse, yuse, 0, 0, xknots_aph, yknots_aph, nxcoef,
     .          nycoef, kxords_aph, kyords_aph, rsacoef, ldf_aph, workh, iflag_aph)
      rsa=10**vlog10rsa

c----------------------------------------------------------------------c
      elseif (istabon .eq. 6) then
c			use spline fit to POST93 table data
c      xuse=min(max(xdata_aph(1),log(te/ev_aph)),xdata_aph(nxdata_aph))
c      yuse=min(max(ydata_aph(1),log10(ne)),ydata_aph(nydata_aph))

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))

      xuse=rle
      yuse=rld

      nxcoef=nxdata_aph
      nycoef=nydata_aph

      tsval = gettime(sec4)
      rsa = B2VAhL(xuse, yuse, 0, 0, xknots_aph, yknots_aph, nxcoef, nycoef,
     .          kxords_aph, kyords_aph, rsacoef, ldf_aph, workh, iflag_aph)
      totb2val = totb2val + gettime(sec4) - tsval      
c----------------------------------------------------------------------c
      elseif (istabon .eq. 7) then
c                       use polynomial fit from Bob Campbell -  8/93
         rsa = sionf(te/ev_aph,ne)

c----------------------------------------------------------------------c
      elseif ((istabon>9 .and. istabon<14) .or. istabon .eq. 17) then 
c     log-log interp on Stotler DEGAS2 tables

         jr = 1

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     ionization rate parameter -- now logarithm of rates
         rsa11=log( wsveh(je,jd,jr) )
         rsa12=log( wsveh(je,jd+1,jr) )
         rsa21=log( wsveh(je+1,jd,jr) )
         rsa22=log( wsveh(je+1,jd+1,jr) )
         rsa1=rsa11 + fjd*(rsa12-rsa11) 
         rsa2=rsa21 + fjd*(rsa22-rsa21) 
         rsa = exp( rsa1 + fje*(rsa2-rsa1) )

c----------------------------------------------------------------------c
      elseif (istabon>13 .and. istabon<16) then # spatially-dep opt-depth

         if (tau.le.taumin) then
            jr  = 1
            fjr = tau/taumin
         else
            zlogr = log10(tau/taumin)/deltau + 2.
            zlogr = min(zlogr,real(mpr-1))
            jr  = int(zlogr)
            fjr = zlogr - jr
         endif

c     compute abscissae --
         zloge=log(te/ev_aph)
         rle=max(rlemin, min(zloge,rlemax))
         zlogd=log10(ne)
         rld=max(rldmin, min(zlogd,rldmax))
c     table indicies for interpolation --
         je=int((rle-rlemin)/delekpt) + 1
         je=min(je,mpe-1)
         jd=int((rld-rldmin)/deldkpt) + 1
         jd=min(jd,mpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-ekpt(je))/(ekpt(je+1)-ekpt(je))
         fjd=(rld-dkpt(jd))/(dkpt(jd+1)-dkpt(jd))
c     ionization rate parameter
         if (istabon .eq. 14) then
            rsa11=(1.-fjr)*wsveh(je,jd,jr)     +
     .                fjr *wsveh(je,jd,jr+1)
            rsa12=(1.-fjr)*wsveh(je,jd+1,jr)   +
     .                fjr *wsveh(je,jd+1,jr+1)
            rsa21=(1.-fjr)*wsveh(je+1,jd,jr)   +
     .                fjr *wsveh(je+1,jd,jr+1)
            rsa22=(1.-fjr)*wsveh(je+1,jd+1,jr) +
     .                fjr *wsveh(je+1,jd+1,jr+1)
            rsa1=rsa11+fjd*(rsa12-rsa11)
            rsa2=rsa21+fjd*(rsa22-rsa21)
            rsa =rsa1 + fje*(rsa2-rsa1)
         elseif (istabon .eq. 15) then
            rsa11=(1.-fjr)*log( wsveh(je,jd,jr) )     +
     .                fjr *log( wsveh(je,jd,jr+1) )
            rsa12=(1.-fjr)*log( wsveh(je,jd+1,jr) )   +
     .                fjr *log( wsveh(je,jd+1,jr+1) )
            rsa21=(1.-fjr)*log( wsveh(je+1,jd,jr) )   +
     .                fjr *log( wsveh(je+1,jd,jr+1) )
            rsa22=(1.-fjr)*log( wsveh(je+1,jd+1,jr) ) +
     .                fjr *log( wsveh(je+1,jd+1,jr+1) )
            rsa1=rsa11+fjd*(rsa12-rsa11)
            rsa2=rsa21+fjd*(rsa22-rsa21)
            rsa = exp( rsa1 + fje*(rsa2-rsa1) )
         endif

c----------------------------------------------------------------------c
         elseif (istabon==16) then
            zn=1     # nuclear charge for hydrogenic species
            zamax=1  # maximum atomic charge for hydrogenic species
            za=0     # compute ionization rate for this charge state
            call mcrates(ne,te,te,za,zamax,zn,rsacopy,kdum,kdum)
            rsa=rsacopy

c----------------------------------------------------------------------c
      endif

      return
      end

c----------------------------------------------------------------------c
      real function sionf(temp,den)
c
c  this function is the new H ionization curve fit from
c  R.B. Campbell 1/94
c
c    0.1 eV < te < 500 eV
c
c    1.0e18 1/m3 < ne < 1.022 1/m3
c
c    y = log10(te(eV))
c    x = log10(ne(1/m3))
c caution: other version may reverse x,y and then correct later
c 
c Both ionization and recombination rates are in m3/sec.
c Etai is in eV.
c
      implicit none

      real temp,den,x,y,ain,bin,cin,din,ein,gin,hin,riin

c Fit for Ionization Rate

      ain(x) = -49.05905 + 2.51313783 * x - 0.049159714*x*x
      bin(x) = 41.1855162  - 2.3298672 * x + 4.24769144e-2*x*x
      cin(x) = -32.798921+1.72102919*x-0.038692357*x*x
      din(x) = 27.370466-1.6824361*x+0.0462317894*x*x
      ein(x) = -7.9990454+0.127573157*x-6.3586911e-3*x*x
      gin(x) = -4.5832951+0.776264783*x-1.8866089e-2*x*x
      hin(x) = 3.08056833-0.39114789*x+9.86833304e-3*x*x
      riin(x) = -0.4648639+0.0551428018*x-1.404213e-3*x*x
c
c *************************************************************
c
      x = min(22.e0,log10(den))
      y = log10(temp)
c
c
      sionf = 10**( ain(x)+bin(x)*y+cin(x)*y*y+din(x)*y*y*y
     .             + ein(x)*y*y*y*y + gin(x)*y*y*y*y*y
     .             + hin(x)*y*y*y*y*y*y + riin(x)*y*y*y*y*y*y*y )

c
      return
      end
c
c----------------------------------------------------------------------c

      real function srecf(temp,den)
c
c  this function is the H recombination curve fit from
c  R.B. Campbell 8/93
c
c    0.1 eV < te < 500 eV
c
c    1.0e18 1/m3 < ne < 1.022 1/m3
c
c    y = log10(te(eV))
c    x = log10(ne(1/m3))
c caution: other version may reverse x,y and then correct later
c 
c Both ionization and recombination rates are in m3/sec.
c Etai is in eV.
c
c
      implicit none

      real temp,den,x,y,ar,br,cr,dr,er,gr
c
c Fit for Recombination Rate
c
      ar(x) = -0.4575652 - 2.144012 * x + 6.7072142e-2 * x * x
     .        -1.391667e-4 * x * x * x
      br(x) = -121.8401 + 18.001822 * x -0.8679488 * x * x
     .        + 1.33165e-2 * x * x * x
      cr(x) = 80.897256 -13.29602 * x + 0.71881414 * x * x
     .        -0.0126549 * x * x * x
      dr(x) = 56.406823 - 7.301996 * x + 0.29339793 * x * x
     .       -3.50898e-3 * x * x * x
      er(x) = -55.73559 + 7.9634283 * x - 0.370274 * x * x
     .       + 5.567961e-3 * x * x * x
      gr(x) = 10.866692 - 1.584193 * x + 0.07563791 * x * x
     .       -1.177562e-3 * x * x * x
c
c *************************************************************
c
      x = min(22.e0,log10(den))
      y = log10(temp)
c
c
      srecf = 10**( ar(x)+br(x)*y+cr(x)*y*y+dr(x)*y*y*y
     .                 + er(x)*y*y*y*y + gr(x)*y*y*y*y*y )
c
      return
      end

c----------------------------------------------------------------------c

      real function svradp(temp,den)
c
c  this function is the H radiation curve fit from
c  R.B. Campbell 8/93
c
c    0.1 eV < te < 500 eV
c
c    1.0e18 1/m3 < ne < 1.022 1/m3
c
c    y = log10(te(eV))
c    x = log10(ne(1/m3))
c caution: other version may reverse x,y and then correct later
c 
c Both ionization and recombination rates are in m3/sec.
c Etai is in eV.
c

      implicit none

      real temp,den,x,y,ai,bi,ci,di,ei,gi,ae,be,ce,de,ee,ge,sionfl,etai

c Fit for Ionization Rate

      ai(x) = -275.845 + 37.010817 * x -1.788045 * x * x
     .      +  0.029078333 * x * x * x
      bi(x) = 2200.9478 - 326.1153 * x + 16.148655 * x * x
     .       - 0.2660702 * x * x * x
      ci(x) = -2.935221e3 + 4.3757698e2 * x -21.73964 * x * x
     .       + 0.358962 * x * x * x
      di(x) = 1604.1466 - 239.6959 * x + 11.923707 * x * x
     .       - 0.1970501 * x * x * x
      ei(x) = -390.8635 + 58.474495 * x -2.910997 * x * x
     .       + 0.048133829 * x * x * x
      gi(x) = 35.012574 - 5.24202 * x + 0.26109962 * x * x
     .       -4.319238e-3 * x * x * x
c
c Fit for (Erad-13.6)*svion Net hydrogenic  energy loss / ionization

      ae(x) = 2860.4173 - 610.2452 * x + 48.275821 * x * x
     .   -1.687994 * x * x * x + 0.02201375 * x * x * x * x
      be(x) = 10612.067 - 2046.397 * x + 147.73914 * x * x
     .       - 4.729973 * x * x * x + 0.056671796 * x * x * x * x
      ce(x) = -4.231708e4 + 8494.6102 * x - 639.0226 * x * x
     . + 21.350311 * x * x * x -0.2673466 * x * x * x * x
      de(x) = -8.385144e3 + 1887.6244 * x -157.8502 * x * x
     .  + 5.820501 * x * x * x - 0.07992837 * x * x * x * x
      ee(x) = 3.938282e4 -8.131339e3 * x + 628.8119 * x * x
     . -21.58636 * x * x * x + 0.27756029 * x * x * x * x
      ge(x) = -1.038281e4 + 2.1349333e3 * x -164.4201 * x * x
     .  + 5.6210487 * x * x * x - 0.07197622 * x * x * x * x
c
      sionfl(x,y) = 10**( ai(x)+bi(x)*y+ci(x)*y*y+di(x)*y*y*y
     .                 + ei(x)*y*y*y*y + gi(x)*y*y*y*y*y )
c
      etai(x,y) =( 10**( ae(x)+be(x)*y+ce(x)*y*y+de(x)*y*y*y
     .                 + ee(x)*y*y*y*y + ge(x)*y*y*y*y*y ) )
     .           / sionfl(x,y)
c
c *************************************************************
c
      x = min(22.e0,log10(den))
      y = log10(temp)
c
c
c Above 100eV, etai is constant at the 100eV value
c
      svradp = max( 0.e0,(13.6e0+etai(x,min(2.e0,y))) ) * 
     .         1.602e-19 * sionfl(x,y)
c
      return
      end

c----------------------------------------------------------------------c

      real function svdiss (te)
      implicit none
      real te

Use(Physical_constants)   # ev

c     Compute rate parameter for molecular dissociation by electrons
c        input:  te [J]   =  electron temperature (J)
c        output: svdiss   = <sigma v> (m**3/s)

c     Polynomial fit from HYDHEL data file of EIRENE Monte Carlo code,
c     reaction number 2.2.5, attributed to preprint by Janev.

c     local variables --
      real logt, logsv
      real b0,b1,b2,b3,b4,b5,b6,b7,b8

      b0 = -2.787217511174e+01
      b1 =  1.052252660075e+01
      b2 = -4.973212347860e+00   
      b3 =  1.451198183114e+00
      b4 = -3.062790554644e-01
      b5 =  4.433379509258e-02   
      b6 = -4.096344172875e-03
      b7 =  2.159670289222e-04
      b8 = -4.928545325189e-06   

      logt = log (te/ev_aph)
      logsv = b0+logt*(b1+logt*(b2+logt*(b3+logt*(b4
     .          +logt*(b5+logt*(b6+logt*(b7+logt*b8)))))))
      svdiss = (1e-6)*exp(logsv)
      return
      end

c----------------------------------------------------------------------c
      real function sv_crumpet (te, ne, rate)
      implicit none
      real te, ne
      integer rate

Use(Dim)
Use(UEpar)
Use(Share)                # istabon
Use(Physical_constants)   # ev
Use(Data_input)
Use(Rtdata)
Use(Rtcrumpet)

c     local variables --
      integer je,jd
      real zloge,zlogd,rle,rld,fje,fjd,mins,maxs,epsilon,
     .     svd_crm11,svd_crm12,svd_crm21,svd_crm22,svd_crm1,svd_crm2


c     Compute molecular dissociation rate - A. Holm's CRUMPET
c     te [J]          = electron temperature
c     ne [/m**3]      = electron density
c     rate            = rate to access
c                       10: mol. diss, X=[]
c                       11: atom source from mol dis, X=[]
c                       20: el. energy change from mol interactions, X=[J]
c                       21: i/a energy change from mol interactions, X=[J]
c                       22: Epot (binding E) change from mol interactions, X=[J]
c                       23: Atom radiation source from mol interactionsm X=[J]
c                       24: Mol. radiation source from mol interactionsm X=[J]
c     svma_crumpet [X/s]    = radiation rate

c----------------------------------------------------------------------c
c     Logarithmic interpolation as in Stotler-95



      if (ishymol .eq. 0) then
            sv_crumpet=0
      elseif (ismolcrm .eq. 0) then
            sv_crumpet=0
      else
c     compute abscissae --
         zloge=log(te/ev)
         rle=max(crlemin, min(zloge,crlemax))
         zlogd=log10(ne)
         rld=max(crldmin, min(zlogd,crldmax))
c     table indicies for interpolation --
         je=int((rle-crlemin)/cdelekpt) + 1
         je=min(je,cmpe-1)
         jd=int((rld-crldmin)/cdeldkpt) + 1
         jd=min(jd,cmpd-1)
c     fractional parts of intervals (je,je+1) and (jd,jd+1) --
         fje=(rle-cekpt(je))/(cekpt(je+1)-cekpt(je))
         fjd=(rld-cdkpt(jd))/(cdkpt(jd+1)-cdkpt(jd))

c     pick the appropriate matrix to interpolate from
         if (rate .eq. 10) then
              svd_crm11=crmdiss(je,jd) 
              svd_crm12=crmdiss(je,jd+1)
              svd_crm21=crmdiss(je+1,jd) 
              svd_crm22=crmdiss(je+1,jd+1)
         elseif (rate .eq. 11) then
              svd_crm11=crmarate(je,jd) 
              svd_crm12=crmarate(je,jd+1)
              svd_crm21=crmarate(je+1,jd) 
              svd_crm22=crmarate(je+1,jd+1)
         elseif (rate .eq. 20) then
              svd_crm11=crmselm(je,jd) 
              svd_crm12=crmselm(je,jd+1)
              svd_crm21=crmselm(je+1,jd) 
              svd_crm22=crmselm(je+1,jd+1)
         elseif (rate .eq. 21) then
              svd_crm11=crmsiam(je,jd) 
              svd_crm12=crmsiam(je,jd+1)
              svd_crm21=crmsiam(je+1,jd) 
              svd_crm22=crmsiam(je+1,jd+1)
         elseif (rate .eq. 22) then
              svd_crm11=crmspotm(je,jd) 
              svd_crm12=crmspotm(je,jd+1)
              svd_crm21=crmspotm(je+1,jd) 
              svd_crm22=crmspotm(je+1,jd+1)
         elseif (rate .eq. 23) then
              svd_crm11=crmsrada(je,jd) 
              svd_crm12=crmsrada(je,jd+1)
              svd_crm21=crmsrada(je+1,jd) 
              svd_crm22=crmsrada(je+1,jd+1)
         elseif (rate .eq. 24) then
              svd_crm11=crmsradm(je,jd) 
              svd_crm12=crmsradm(je,jd+1)
              svd_crm21=crmsradm(je+1,jd) 
              svd_crm22=crmsradm(je+1,jd+1)
         else    
              call remark("Rate option not recognized!")
         endif
              
         
c ... Check whether we have all positive, all negative, or zero crossing
            mins=min(svd_crm11,svd_crm12,svd_crm21,svd_crm22)
            maxs=max(svd_crm11,svd_crm12,svd_crm21,svd_crm22)
            if (mins .gt. 0) then
c ... Minimum is positive, normal interpolation
                svd_crm11=log(svd_crm11)
                svd_crm12=log(svd_crm12)
                svd_crm21=log(svd_crm21)
                svd_crm22=log(svd_crm22)
                svd_crm1=svd_crm11+fjd*(svd_crm12-svd_crm11)
                svd_crm2=svd_crm21+fjd*(svd_crm22-svd_crm21)
                sv_crumpet = exp( svd_crm1 + fje*(svd_crm2-svd_crm1) )

            elseif (maxs .lt. 0) then
c ... Maximum is negative, interpolate of negative values
                svd_crm11=log(-svd_crm11)
                svd_crm12=log(-svd_crm12)
                svd_crm21=log(-svd_crm21)
                svd_crm22=log(-svd_crm22)
                svd_crm1=svd_crm11+fjd*(svd_crm12-svd_crm11)
                svd_crm2=svd_crm21+fjd*(svd_crm22-svd_crm21)
                sv_crumpet = -exp( svd_crm1 + fje*(svd_crm2-svd_crm1) )

            else
c ... Zero-crossing – linear interpolation   
                jd=jd-1
                je=je-1
                fjd=((ne/1e6)-(10**(10+0.5*jd)))/ (10**(10+0.5*jd)*(10**0.5 - 1) )
                fje=(te/ev - (10**(-1.2+0.1*je)))/( (10**(-1.2+0.1*je))*(10**0.1-1))

                svd_crm1=svd_crm11+fjd*(svd_crm12-svd_crm11)
                svd_crm2=svd_crm21+fjd*(svd_crm22-svd_crm21)
                sv_crumpet = svd_crm1 + fje*(svd_crm2-svd_crm1)
                

            endif

        endif
        

c----------------------------------------------------------------------c
      return
      end












