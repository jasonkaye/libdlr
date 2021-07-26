      !
      !
      ! This file contains various utility routines to work with the
      ! discrete Lehmann representation
      !
      !




      !> Get equispaced points on [0,1] in relative format
      !!
      !! For the definitions of absolute and relative formats, see the
      !! readme.
      !!
      !! Endpoints of [0,1] are included
      !!
      !! @param[in]   n number of points on [0,1]
      !! @param[out]  t equispaced points on [0,1], in relative format

      subroutine eqpts_rel(n,t)

      implicit none
      integer n
      real *8 t(n)

      integer i

      do i=1,n-1

        if (i.le.n/2) then
          t(i) = (i-1)*1.0d0/(n-1)
        else
          t(i) = -(n-i)*1.0d0/(n-1)
        endif

      enddo

      t(n) = 1.0d0

      end subroutine eqpts_rel





      !> Convert points on [0,1] from relative format to absolute format
      !!
      !! For the definitions of absolute and relative formats, see the
      !! readme.
      !!
      !! Note: converting a point from relative to absolute format will,
      !! in general, constitute a loss of relative accuracy in the
      !! location of the point if the point is close to t = 1. For
      !! example, in three-digit arithmetic, the point t = 0.999111 could
      !! be stored as t* = -0.889e-3 in the relative format, but only as
      !! t = 0.999 in the absolute format.
      !!
      !! @param[in]   n     number of points to convert
      !! @param[in]   t     array of points in relative format
      !! @param[out]  tabs  array of points in absolute format

      subroutine rel2abs(n,t,tabs)

      implicit none
      integer n
      real *8 t(n),tabs(n)

      integer i

      do i=1,n

        if (t(i).lt.0.0d0) then
          tabs(i) = t(i)+1.0d0
        else
          tabs(i) = t(i)
        endif

      enddo

      end subroutine rel2abs





      !> Convert a point on [0,1] from absolute format to relative
      !! format
      !!
      !! For the definitions of absolute and relative formats, see the
      !! readme.
      !!
      !! Note: In order to specify points -- for example points at which
      !! to sample or evaluate a DLR -- in absolute format, those points
      !! must first be converted to relative format using this subroutine
      !! before using them as inputs into any other subroutines. Of
      !! course, in order to maintain full relative precision in all
      !! calculations, it is necessary to specify points in relative
      !! format from the beginning, but in most cases at most a mild loss
      !! of accuracy will result from using the absolute format.
      !!
      !! @param[in]   n     number of points to convert
      !! @param[in]   tabs  array of points in absolute format
      !! @param[out]  t     array of points in relative format


      subroutine abs2rel(n,tabs,t)

      implicit none
      integer n
      real *8 t(n),tabs(n)

      integer i

      do i=1,n

        if (t(i).gt.0.5d0.and.t(i).lt.1.0d0) then
          t(i) = tabs(i)-1.0d0
        else
          t(i) = tabs(i)
        endif

      enddo

      end subroutine abs2rel
