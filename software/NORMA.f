      PROGRAM NORMA
c
c     find a minimum for URO and NMA
c
c     Wed Sep 21 10:53:37 CEST 2005
c     adapted to URO cases for Gif
c
c     Fri Sep 23 11:49:40 CEST 2005
c     add external parameterization using NAMELIST
c
c     Tue Oct  4 13:35:53 CEST 2005
c     add reading of initial perturbation from NAMELIST
c
c     Fri Oct 21 22:51:10 CEST 2005
c     adapt to g77
c
c     Sat Oct 22 21:16:27 CEST 2005
c     add printing of simplex edges at every step
c
      IMPLICIT NONE 

      INTEGER maxdim
      PARAMETER (maxdim=100)

c     .. loop counter ..
      INTEGER i, j, k, l

c     .. variables for amebsa ..
      INTEGER niter, nanneal, iter,mp,ndim,np,NMAX, nround
      REAL tol, ftol,tt,temptr,yb,
     +     p(maxdim+1,maxdim),pb(maxdim),y(maxdim+1),funk,
     +     bestyb, bestpb(maxdim)
      REAL dy
      REAL ptest(maxdim)
      INTEGER mode(maxdim)
      REAL seed
      LOGICAL lrand

c     .. verbosity level ..
      INTEGER NVERB

      EXTERNAL funk
      COMMON /NDIM/ ndim, mode, nverb

c     .. NAMELIST for parameterization ..
      NAMELIST /param/ ndim, dy, niter, temptr, nanneal, nround,
     +                 lrand, tol, mode, pb, seed, nverb

c     .. initial NMA perturbation (overwrite using Namelist) ..
      do j = 1, maxdim
        pb(j) = 0.0
      end do

c     .. number of modes to use
      ndim = 1

c     .. expected DQ range in which minimum shall be searched
      dy = +100.0

c     .. tolerance ..
      tol = 1e-3

c     .. number of iterations at fixed temperature
      niter = 50

c     .. starting temperature
      temptr = 100.0

c     .. number of simulated annealing steps (0=no annealing) ..
      nanneal = 0

c     .. number of rounds ..
      nround = 1

c     .. add a random perturbation onto the original coordinates
c        so that not all procs do the same search (multiprocessor)
      lrand = .false.

c     .. seed for random number generator ..
      seed = 42.4711

c     .. initialize mode numbers ..
      do i=1, maxdim 
        mode(i) = i + 6
      end do

c     .. verbosity level (0=sparse, 1=normal, 2=much)
      nverb = 1

c     .. read NAMELIST 
      print *, "opening NORMA.inp for NAMELIST input"
      open(42, file="NORMA.inp", status="old", err=888)
      print *, "reading NAMELIST param"
      read(42, param, err=889)
      goto 887
888   print *, "WARNING: failed to open NORMA.inp"
889   print *, "WARNING: failed to read NAMELIST param"
      print *, "         default parameters will be used"
887   continue

c     .. set unused modes to zero ..
      do i=ndim + 1, maxdim
        mode(i) = 0
      end do

c     .. print parameters ..
      print *, "dumping NAMELIST param:"
      write(*, param)

      np=ndim
      mp=ndim+1

c     .. start random number generator ..
      write(*,*) "SEED = ", seed 
      call srand(int(seed*1000))

      if (lrand) then
        print *, 'initial values will be randomized'
      else
        print *, 'initial values will NOT be randomized'
      end if

c     .. outer loop for repeated simulated annealing
c     .. after each round we initialize around the minimu
c        using a smaller perturbation dy ..

      do l = nround, 1, -1

        print *, "Round ", nround - l + 1, " out of ", nround
c       .. initialize ..
        do i = 1, mp
          do j = 1, np

      if (lrand) then
            p(i,j) = pb(j) - 1./(2.*ndim) * ((rand(0)+0.5)*l/nround)*dy
      else
            p(i,j) = pb(j) - 1./(2.*ndim) * (1.*l/nround)*dy
      end if

            if ( j.eq.i ) then
      if (lrand) then
              p(i,j) = p(i,j) + ((rand(0)+0.5)*l/nround)*dy
      else
              p(i,j) = p(i,j) + (1.*l/nround)*dy
      end if
            end if
          enddo
        enddo

c       .. search in the reverse direction on next round
        dy = -dy 

        print *, "computing initial values"
        do i = 1, mp
          do j = 1, np
            ptest(j) = p(i,j)
          enddo
          y(i) = funk(ptest)
          if (i.eq.1) yb = y(i)
          if (y(i).lt.yb) then
            do j = 1, np
              pb(j) = p(i,j)
              yb = y(i)
            enddo
          end if
        enddo

ccc     .. annealing / minimization steps ..
        do k = nanneal, 0, -1

          iter = niter
          ftol = tol
          tt = 1.0*k
          if (nanneal.gt.0) tt = tt / nanneal

          print *, "amebsa minimization/annealing step ", nanneal-k+1,
     +             " out of ", nanneal+1
          print *, "temp =       ", tt*temptr

          call amebsa(p,y,maxdim+1,maxdim,ndim,pb,yb,ftol,funk,iter,
     +              tt*temptr)

          print *, "minimum =    ", yb
          print *, "iter =       ", iter


        enddo

        print *, "*** result for this round ***"
        do j = 1, np
          print *, j, pb(j)
        enddo 
        print *, "minimum = ", yb

c       .. keep the best result between rounds ..
        if (l .eq. nround) then
          bestyb = yb
          do j = 1, np
            bestpb(j) = pb(j)
          end do
        else if (yb < bestyb) then
          bestyb = yb
          do j = 1, np
            bestpb(j) = pb(j)
          end do
        else
          yb= bestyb
          do j = 1, np
            pb(j) = bestpb(j)
          end do
        end if

c     .. outer loop for repeated simulated annealing
      end do

      print *, "computing final solution"

      yb = funk ( pb )

      print *, "*** Final result ***"
      do j = 1, np
        print *, j, bestpb(j)
      enddo 
      print *, "Final minimum = ", bestyb

      END
      REAL FUNCTION funk (dq)
c     .. function to minimize ..

      IMPLICIT NONE
      INTEGER ndim
      INTEGER maxdim
      PARAMETER (maxdim=100)
      REAL dq(maxdim)
      INTEGER mode(maxdim)

c     .. to execute an external program using IFPOSIX ..
      INTEGER ierror, i
      INTEGER nverb
      COMMON /NDIM/ ndim, mode, nverb

c     .. write func.in ..
      OPEN(77,NAME="func.in")

      WRITE(77,"( 99(I3,F14.6) )")
     +           ( mode(i),dq(i), i=1, ndim )

      CLOSE(77)

      call system ("./func.sh")

c     .. read func.out ..
      OPEN(78,NAME="func.out", ERR=999)
      READ(78,*, ERR=999, END=999) funk
      CLOSE(78)
      GOTO 998

c     .. if the external program fails, stop ..
999   CONTINUE
      PRINT *, "ERROR: failed to read func.out"
      STOP

C     .. print some information ..
998   CONTINUE
      IF (nverb > 0) THEN
        PRINT '(A,F14.6,A,99(I3,F14.6))', "monitor=", funk, ", p=",
     +           (mode(i), dq(i), i=1, ndim)
      END IF

      RETURN 
      END

      REAL FUNCTION ran1 (idum)
      INTEGER idum
      ran1 = rand(0)
      RETURN
      END

      SUBROUTINE amebsa(p,y,mp,np,ndim,pb,yb,ftol,funk,iter,temptr)
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
      INTEGER iter,mp,ndim,np,NMAX
      REAL ftol,temptr,yb,p(mp,np),pb(np),y(mp),funk
      PARAMETER (NMAX=200)
      EXTERNAL funk
CU    USES amotsa,funk,ran1
      INTEGER i,idum,ihi,ilo,inhi,j,m,n
      REAL rtol,sum,swap,tt,yhi,ylo,ynhi,ysave,yt,ytry,psum(NMAX),
     *amotsa,ran1
      COMMON /ambsa/ tt,idum
      tt=-temptr
1     do 12 n=1,ndim
        sum=0.
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      inhi=1
      ihi=2
      ylo=y(1)+tt*log(ran1(idum))
      ynhi=ylo
      yhi=y(2)+tt*log(ran1(idum))
      if (ylo.gt.yhi) then
        ihi=1
        inhi=2
        ilo=2
        ynhi=yhi
        yhi=ylo
        ylo=ynhi
      endif
      do 13 i=3,ndim+1
        yt=y(i)+tt*log(ran1(idum))
        if(yt.le.ylo) then
          ilo=i
          ylo=yt
        endif
        if(yt.gt.yhi) then
          inhi=ihi
          ynhi=yhi
          ihi=i
          yhi=yt
        else if(yt.gt.ynhi) then
          inhi=i
          ynhi=yt
        endif
13    continue
      rtol=2.*abs(yhi-ylo)/(abs(yhi)+abs(ylo))
      if (rtol.lt.ftol.or.iter.lt.0) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
      iter=iter-2
      ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,-1.0)
      if (ytry.le.ylo) then
        ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,2.0)
      else if (ytry.ge.ynhi) then
        ysave=yhi
        ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,0.5)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum)
            endif
16        continue
          iter=iter-ndim
          goto 1
        endif
      else
        iter=iter+1
      endif
      goto 2
      END
      FUNCTION amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,fac)
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
      INTEGER ihi,mp,ndim,np,NMAX
      REAL amotsa,fac,yb,yhi,p(mp,np),pb(np),psum(np),y(mp),funk
      PARAMETER (NMAX=200)
      EXTERNAL funk
CU    USES funk,ran1
      INTEGER idum,j
      REAL fac1,fac2,tt,yflu,ytry,ptry(NMAX),ran1
      COMMON /ambsa/ tt,idum
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry)
      if (ytry.le.yb) then
        do 12 j=1,ndim
          pb(j)=ptry(j)
12      continue
        yb=ytry
      endif
      yflu=ytry-tt*log(ran1(idum))
      if (yflu.lt.yhi) then
        y(ihi)=ytry
        yhi=yflu
        do 13 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
13      continue
      endif
      amotsa=yflu
      return
      END
