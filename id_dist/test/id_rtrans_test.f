c
c
c       dependencies: prini, id_rand
c
c
        implicit real *8 (a-h,o-z)
c 
        call prini(6,13)
c 
c        set all parameters
c 
        print *, 'Enter n:'
        read *,n
        call prinf('n=*',n,1)
c
        call lams_test(n)
        call lams_test_complex(n)
c
        stop
        end
c
c 
c 
c 
c 
        subroutine lams_test_complex(n)
        implicit real *8 (a-h,o-z)
        save 
        complex *16 w(1000 000),x(1000 000),y(1000 000),z(1000 000)
c
c       tests the complex-valued random transformations.
c
c       input:
c       n -- size of transform to be tested
c
c        construct the input vector
c
        do 1200 i=1,n
c
        x(i)=i
 1200 continue
c
        nsteps=10
c
c        apply the operator
c
        call idz_random_transf_init(nsteps,n,w,keep)
        call idz_random_transf(x,y,w)
c
        call prin2('x=*',x,n*2)        
        call prin2('y=*',y,n*2)        
c
c        apply inverse operator
c
        call idz_random_transf_inverse(y,z,w)
        call prin2('and z=*',z,n*2)        
c
        return
        end
c
c 
c 
c 
c 
        subroutine lams_test(n)
        implicit real *8 (a-h,o-z)
        save 
        dimension w(1000 000),x(1000 000),y(1000 000),z(1000 000)
c
c       tests the real-valued random transformations.
c
c       input:
c       n -- size of transform to be tested
c
c        construct the input vector
c
        do 1200 i=1,n
c
        x(i)=i
 1200 continue
c
        nsteps=10
c
c        apply the operator
c
        call idd_random_transf_init(nsteps,n,w,keep)
        call idd_random_transf(x,y,w)
c
        call prin2('x=*',x,n)        
        call prin2('y=*',y,n)        
c
c        apply inverse operator
c
        call idd_random_transf_inverse(y,z,w)
        call prin2('and z=*',z,n)        
c
        return
        end
c
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
