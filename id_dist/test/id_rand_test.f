c
c
c       dependencies: prini
c
c
        implicit none
c
        integer len
        parameter(len = 5 000 000)
c
        integer n,k,m,ind(len)
        real*8 r(len),s(55),temp,r2(len),diff
c
        data s/
     1  0.2793574644042651d0, 0.1882566493961346d0,
     2  0.5202478134503912d0, 0.7568505373052146d0,
     3  0.5682465992936152d0, 0.5153148754383294d0,
     4  0.7806554095454596d0, 1.982474428974643d-2,
     5  0.2520464262278498d0, 0.6423784715775962d0,
     6  0.5802024387972178d0, 0.3784471040388249d0,
     7  7.839919528229308d-2, 0.6334519212594525d0,
     8  3.387627157788001d-2, 0.1709066283884670d0,
     9  0.4801610983518325d0, 0.8983424668099422d0,
     *  5.358948687598758d-2, 0.1265377231771848d0,
     1  0.8979988627693677d0, 0.6470084038238917d0,
     2  0.3031709395541237d0, 0.6674702804438126d0,
     3  0.6318240977112699d0, 0.2235229633873050d0,
     4  0.2784629939177633d0, 0.2365462014457445d0,
     5  0.7226213454977284d0, 0.8986523045307989d0,
     6  0.5488233229247885d0, 0.3924605412141200d0,
     7  0.6288356378374988d0, 0.6370664115760445d0,
     8  0.5925600062791174d0, 0.4322113919396362d0,
     9  0.9766098520360393d0, 0.5168619893947437d0,
     *  0.6799970440779681d0, 0.4196004604766881d0,
     1  0.2324473089903044d0, 0.1439046416143282d0,
     2  0.4670307948601256d0, 0.7076498261128343d0,
     3  0.9458030397562582d0, 0.4557892460080424d0,
     4  0.3905930854589403d0, 0.3361770064397268d0,
     5  0.8303274937900278d0, 0.3041110304032945d0,
     6  0.5752684022049654d0, 7.985703137991175d-2,
     7  0.5522643936454465d0, 1.956754937251801d-2,
     8  0.9920272858340107d0/
c
c
        call prini(6,13)
c
c
        print *,'Enter n:'
        read *,n
        call prinf('n = *',n,1)
c
c
c       Generate n random numbers uniformly drawn from [0,1].
c
        call id_frand(n,r)
        call prin2('r = *',r,n)
c
c       Generate n more random numbers uniformly drawn from [0,1].
c
        call id_frand(n,r)
        call prin2('r = *',r,n)
c
c       Initialize the seed values in id_frand
c       to their original values.
c
        call id_frando()
c
c       Generate n more random numbers uniformly drawn from [0,1].
c
        call id_frand(n,r)
        call prin2('r = *',r,n)
c
c
c       Print the percentiles of r.
c
        m = 10
        call histogram(n,r,m)
c
c
c       Reverse the order of the seed values in s.
c
        do k = 1,55/2
          temp = s(k)
          s(k) = s(55-k+1)
          s(55-k+1) = temp
        enddo ! k
c
c
c       Generate r2 using id_srand so that it should match r
c       generated using id_frand.
c
        call id_srandi(s)
        call id_srand(n,r2)
c
c
c       Compute and print the difference between r and r2.
c
        diff = 0
c
        do k = 1,n
          diff = diff+abs(r(k)-r2(k))
        enddo ! k
c
        call prin2('diff = *',diff,1)
c
c
c       Generate and display a random permutation.
c
        call id_randperm(n,ind)
        call prinf('ind = *',ind,n)
c
c
        stop
        end
c
c
c
c
        subroutine histogram(n,r,m)
c
c       counts and prints the number of entries of r falling
c       into m equally wide bins partitioning [0,1].
c
c       input:
c       n -- length of r
c       r -- array to be binned
c       m -- number of bins
c
        implicit none
        integer m,n,nbin,j,k,iarr(2)
        real*8 r(n),width,r1
c
        r1 = 1
c
c
        width = r1/m
c
c
        do j = 1,m
c
          nbin = 0
c
          do k = 1,n
c
            if(r(k) .gt. (j-1)*width .and. r(k) .le. j*width)
     1       nbin = nbin+1
c
          enddo ! k
c
          iarr(1) = j
          iarr(2) = nbin
          call prinf('(j,nbin) = *',iarr,2)
c
        enddo ! j
c
c
        return
        end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
