cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the
c        start of the actual quadrature routines.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine adapgaus_quad(ier,a,b,fun,par1,par2,m,eps,
     1      rint,maxrec,numint)
        implicit real *16 (a-h,o-z)
cccc        save
        dimension t(100),w(100),stack(400),vals(200),
     1      par1(1),par2(1)
        external fun
cccc        data m7/-2341034/
c 
c       this subroutine uses the adaptive chebychev quadrature
c       to evaluate the integral of the user-supplied function
c       fun on the user-specified interval [a,b]
c 
c                       input parameters:
c 
c  a,b - the ends of the interval on which the integral is to
c       be evaluated
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be
c 
c        fun(x,par1,par2).                            (1)
c 
c        in (1), x is a point on the interval [a,b] where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be
c        variables or arrays, real or integer, as desired.
c        fun is assumed to be real *8.
c  par1, par2 - partameters to be used by the user-supplied
c       function fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=8 means that at some point, one subinterval in the
c                subdivision was smaller than (b-a)/2**200. this
c                is a fatal error.
c          ier=16 means that the total number of subintervals in the
c                adaptive subdivision of [a,b] turned out to be greater
c                than 100000.  this is a fatal error.
c 
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the totla number of intervals in the subdivision. can not
c         be greater than 100000,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c 
c 
c          . . . check if the gauss quadarture has been
c                initialized at a preceeding call to this routine
c 
cccc        if(m .eq. m7) goto 1200
cccc        call gauswhts(m,t,w)
c 
        ifwhts=1
        call legewhts_quad(m,t,w,ifwhts)
c 
cccc        m7=m
cccc 1200 continue
c 
c        integrate the user-supplied function using the
c        adaptive gaussian quadratures
c 
        nnmax=100000
        maxdepth=200
c 
        call adinrec_quad(ier, stack, a, b, fun, par1, par2, t,
     1      w, m, vals, nnmax, eps, rint, maxdepth, maxrec, numint)
        return
        end
c 
c 
c 
c 
c 
        subroutine adinrec_quad(ier,stack,a,b,fun,
     1      par1,par2,t,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint)
        implicit real *16 (a-h,o-z)
cccc        save
        dimension stack(2,1),t(1),w(1),vals(1),par1(1),par2(1)
        external fun
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        call oneint_quad(a,b,fun,par1,par2,t,w,m,vals(1))
c 
c       recursively integrate the thing
c 
        j=1
        rint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
ccc        call prinf('j=*',j,1)
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
        call oneint_quad(stack(1,j),c,fun,
     1      par1,par2,t,w,m,value2)
c 
        call oneint_quad(c,stack(2,j),fun,
     1      par1,par2,t,w,m,value3)
c 
        dd=abs(value2+value3-vals(j))
cccc         call prin2('in adinrec, dd=*',dd,1)
        ifdone=0
        if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if(ifdone  .eq. 0) goto 2000
c 
        rint=rint+value2+value3
        j=j-1
c 
c        if the whole thing has been integrated - return
c 
        if(j .eq. 0) return
        goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(j+1)=value2
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
        vals(j)=value3
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .le. maxdepth) goto 3000
        ier=8
        return
 3000 continue
        ier=16
        return
        end
c 
c 
c 
c 
c 
        subroutine oneint_quad(a,b,fun,par1,par2,t,w,m,rint)
        implicit real *16 (a-h,o-z)
        external fun
        dimension t(1),w(1),par1(1),par2(1)
c 
c       integrate the function fun on the interval [a,b]
c 
        rint=0
        u=(b-a)/2
        v=(b+a)/2
c
        do 1200 i=1,m
            tt=u*t(i)+v
            rint=rint+fun(tt,par1,par2)*w(i)
 1200 continue
c
        rint=rint*u
c
        return
        end
c
c
c
c
c
        subroutine legewhts_quad(n,ts,whts,ifwhts)
        implicit real *16 (a-h,o-z)
        dimension ts(1),whts(1),ws2(1000),rats(1000)
c 
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on
c        the interval [-1,1]
c 
c                input parameters:
c 
c  n - the number of nodes in the quadrature
c 
c                output parameters:
c 
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c 
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c 
        eps=1.0q-28
        ZERO=0
        DONE=1
        pi=atan(done)*4
        h=pi/(2*n)
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=cos(t)
1200  CONTINUE
c 
c         use newton to find all roots of the legendre polynomial
c 
        ts(n/2+1)=0
        do 2000 i=1,n/2
c 
        xk=ts(i)
        ifout=0
        deltold=1
        do 1400 k=1,10
        call legepol_sum_quad(xk,n,pol,der,sum)
        delta=-pol/der
        xk=xk+delta
        if(abs(delta) .lt. eps) ifout=ifout+1
c 
        if(ifout .eq. 3) goto 1600
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c 
c        construct the weights via the orthogonality relation
c 
        if(ifwhts .eq. 0) return
c 
        do 2400 i=1,(n+1)/2
        call legepol_sum_quad(ts(i),n,pol,der,sum)
        whts(i)=1/sum
        whts(n-i+1)=whts(i)
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legepol_sum_quad(x,n,pol,der,sum)
        implicit real *16 (a-h,o-z)
c 
        done=1
        sum=0
c 
        pkm1=1
        pk=x
        sum=sum+pkm1**2 /2
        sum=sum+pk**2 *(1+done/2)
c 
        pk=1
        pkp1=x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
  
        sum=0
c 
        pol=1
        der=0
        sum=sum+pol**2 /2
        if(n .eq. 0) return
c 
        pol=x
        der=1
        sum=sum+pol**2*(1+done/2)
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        sum=sum+pkp1**2*(k+1+done/2)
 2000 continue
c 
c        calculate the derivative
c 
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
