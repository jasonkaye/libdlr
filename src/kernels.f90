      !
      !
      ! This file contains functions for evaluating Lehmann kernels 
      !
      !


      real *8 function kfunf(t,om)

      ! Fermionic kernel K(tau,omega), tau in absolute format

      implicit none
      real *8 t,om

      if (om.ge.0.0d0) then
        kfunf = exp(-t*om)/(1.0d0+exp(-om))
      else
        kfunf = exp((1.0d0-t)*om)/(1.0d0+exp(om))
      endif

      end function kfunf


      real *8 function kfunf_rel(t,om)

      ! Fermionic kernel K(tau,omega), tau in relative format

      implicit none
      real *8 t,om

      real *8, external :: kfunf

        if (t.ge.0.0d0) then
          kfunf_rel = kfunf(t,om)
        else
          kfunf_rel = kfunf(-t,-om)
        endif

      end function kfunf_rel


      complex *16 function kfunf_mf(n,om)

      ! Fermionic kernel K(i omega_n,omega)

      implicit none
      integer n
      real *8 om

      real *8 pi
      complex *16 eye

      pi = 4*atan(1.0d0)
      eye = (0.0d0,1.0d0)

      kfunf_mf = 1.0d0/((2*n+1)*pi*eye+om)

      end function kfunf_mf


      real *16 function qkfunf(t,om)

      ! Fermionic kernel K(tau,omega), tau in absolute format
      ! Quad precision

      implicit none
      real *16 t,om

      if (om.ge.0.0q0) then
        qkfunf = exp(-t*om)/(1.0q0+exp(-om))
      else
        qkfunf = exp((1.0q0-t)*om)/(1.0q0+exp(om))
      endif

      end function qkfunf


      real *16 function qkfunf_rel(t,om)

      ! Fermionic kernel K(tau,omega), tau in relative format
      ! Quad precision

      implicit none
      real *16 t,om

      real *16, external :: qkfunf

        if (t.ge.0.0q0) then
          qkfunf_rel = qkfunf(t,om)
        else
          qkfunf_rel = qkfunf(-t,-om)
        endif

      end function qkfunf_rel




!      real *8 function kfunb(t,om)
!
!      implicit none
!      real *8 t,om
!
!      if (om.ge.0.0d0) then
!        call evalexpfun(om,kfunb)
!        kfunb = exp(-t*om)/kfunb
!      else
!        call evalexpfun(-om,kfunb)
!        kfunb = exp((1.0d0-t)*om)/kfunb
!      endif
!
!      end function kfunb
