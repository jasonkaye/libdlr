      ! -------------------------------------------------------------
      !
      ! This file contains functions for kernel evaluation
      !
      ! -------------------------------------------------------------
      !
      ! Copyright (C) 2021 The Simons Foundation
      ! 
      ! Author: Jason Kaye
      ! 
      ! -------------------------------------------------------------
      ! 
      ! libdlr is licensed under the Apache License, Version 2.0 (the
      ! "License"); you may not use this file except in compliance with
      ! the License.  You may obtain a copy of the License at
      ! 
      !     http://www.apache.org/licenses/LICENSE-2.0
      ! 
      ! Unless required by applicable law or agreed to in writing,
      ! software distributed under the License is distributed on an "AS
      ! IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
      ! express or implied.  See the License for the specific language
      ! governing permissions and limitations under the License.
      ! 
      ! -------------------------------------------------------------



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
        kfunf = -exp(-t*om)/(1.0d0+exp(-om))
      else
        kfunf = -exp((1.0d0-t)*om)/(1.0d0+exp(om))
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

      kfunmf = 1.0d0/(n*pi*eye-om)

      end function kfunmf





      !> Evaluate Lehmann kernel with imaginary time point given in
      !! relative format, quadruple precision.
      !!
      !! @param[in] t   imaginary time point, in relative format
      !! @param[in] om  real frequency point

      real *16 function qkfunf_rel(t,om)

      implicit none
      real *16 t,om

      real *16, external :: qkfunf

      if (t.ge.0.0q0) then
        qkfunf_rel = qkfunf(t,om)
      else
        qkfunf_rel = qkfunf(-t,-om)
      endif

      end function qkfunf_rel





      !> Evaluate Lehmann kernel with imaginary time point given in
      !! absolute format, quadruple precision.
      !!
      !! Note: the result will not be accurate if t is very close to 1.
      !! To maintain full accuracy, represent t in relative format and
      !! use the function qkfunf_rel.
      !!
      !! @param[in] t   imaginary time point, in absolute format
      !! @param[in] om  real frequency point

      real *16 function qkfunf(t,om)

      implicit none
      real *16 t,om

      if (om.ge.0.0q0) then
        qkfunf = -exp(-t*om)/(1.0q0+exp(-om))
      else
        qkfunf = -exp((1.0q0-t)*om)/(1.0q0+exp(om))
      endif

      end function qkfunf
