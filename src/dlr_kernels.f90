      !
      !
      ! This file contains functions for evaluating Lehmann kernels 
      !
      !




      !> Evaluate Lehmann kernel with imaginary time point given in
      !! relative format.
      !!
      !! @param[in] t   imaginary time point, in relative format
      !! @param[in] om  real frequency point

      function kfunf_rel(t,om)

      implicit none
      real *8 kfunf_rel,t,om

      real *8, external :: kfunf

        if (t.ge.0.0d0) then
          kfunf_rel = kfunf(t,om)
        else
          kfunf_rel = kfunf(-t,-om)
        endif

      end function kfunf_rel





      !> Evaluate Lehmann kernel with imaginary time point given in
      !! absolute format.
      !!
      !! Note: the result will not be accurate if t is very close to 1.
      !! To maintain full accuracy, represent t in relative format and
      !! use the function kfunf_rel.
      !!
      !! @param[in] t   imaginary time point, in absolute format
      !! @param[in] om  real frequency point

      function kfunf(t,om)

      implicit none
      real *8 kfunf,t,om

      if (om.ge.0.0d0) then
        kfunf = exp(-t*om)/(1.0d0+exp(-om))
      else
        kfunf = exp((1.0d0-t)*om)/(1.0d0+exp(om))
      endif

      end function kfunf




      !> Evaluate Lehmann kernel at a Matsubara frequency point.
      !!
      !! @param[in] n   index of Matsubara frequency point i*pi*n
      !! @param[in] om  real frequency point

      function kfunmf(n,om)

      implicit none
      integer n
      real *8 om
      complex *16 kfunmf

      real *8 pi
      complex *16 eye

      pi = 4*atan(1.0d0)
      eye = (0.0d0,1.0d0)

      kfunmf = 1.0d0/(n*pi*eye+om)

      end function kfunmf




!      real *8 function kfunb(t,om)
!     
!      ! A proposed bosonic kernel -- no longer used
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
