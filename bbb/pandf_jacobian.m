c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"


      SUBROUTINE jacobian_store_momentum(xc, yc)
      IMPLICIT NONE
      Use(Selec)
      Use(Jacobian_restore)
      Use(Comflo)
      Use(Cfric)
      Use(Compla)
      Use(Dim)
      integer xc, yc
      integer ix1, ifld

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


      END SUBROUTINE jacobian_store_momentum

      SUBROUTINE jacobian_reset(xc, yc)
      IMPLICIT NONE
      Use(Selec)
      Use(Imprad)
      Use(Rhsides)
      Use(Jacobian_restore)
      Use(Conduc)
      Use(Comflo)
      Use(Cfric)
      Use(Compla)
      Use(Dim)
      Use(MCN_sources)
      Use(Lsode)
      Use(Comtra)
      Use(PNC_params)
      integer ix1, ifld, igsp, xc, yc
c...  Finally, reset some source terms if this is a Jacobian evaluation
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
c ...       TODO: Make psordisold ifld-dependent?
            psordis(xc,yc,ifld) = psordisold
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

c ...   TODO: Don't know what this block is doing, and doing here...
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


      END SUBROUTINE jacobian_reset


