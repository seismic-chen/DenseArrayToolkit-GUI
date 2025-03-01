******************************************************************************
c
      program iterdeconfd
c
******************************************************************************
c
c     Chuck Ammon - Saint Louis University
c     November 1996 / September 1998 
c     VERSION 1.04
c
c     Based on the Kikuchi and Kanamori (1981) iterative
c      deconvolution algorithm. The output of this code
c      is the deconvolution of the "denominator" from
c      the "numerator"
c
c     The final deconvolution is called "decon.out"
c       the observed time series (original convolved with Gaussian)
c        is called "observed"
c       the predicted is call "predicted"
c
c     Header values from the "numerator" are copied into the
c       decon.out file. The gwidth is stored in user0.
c
c******************************************************************************
c
c     if you choose verbose output
c       the deconvolutions are called d???
c       the predictions are called p???
c       the residuals are called r???
c
c      where ??? is corresponds to the number of "bumps"
c
******************************************************************************
c
c     Because the correlation is computed in the 
c     Frequency Domain, the maximum lag computed is npts / 2
c
c     That is, negative lags are not used - so make sure you wavelet starts
c       before the signal.
c
c     A low-pass Gaussian filter is applied to the signals before the
c      deconvolution is performed. The predicted will not match the
c      original signal, but the filtered signal, which is stored in the
c      file called "numerator".
c
******************************************************************************
c
      integer MAXPTS, MAXG
      parameter(MAXPTS = 8192, MAXG = 200)
      real f(MAXPTS),g(MAXPTS),p(MAXPTS),r(MAXPTS)
      real amps(MAXG)
      integer shifts(MAXG)
      character*256 numerator, denominator
      character*12 resfile, filename
      integer stdin,stdout,ounit,inunit,forward,inverse
      logical lpositive,verbose
      
      stdin = 5
      stdout = 6
      ounit = 9
      inunit = 10
      forward =  1
      inverse = -1
      gwidth = 2.5
c
      write(stdout,*) ' '
      write(stdout,*) 'Program iterdeconfd - Version 1.04, 1997-98'
      write(stdout,*) 'Chuck Ammon, Saint Louis University'
      write(stdout,*) ' '
c
c     read in the names of the input files
c      
      write(stdout,*)'What is the numerator file?'
      read(stdin,'(a)') numerator
      write(stdout,*)'What is the denominator file?'
      read(stdin,'(a)') denominator
      write(stdout,*)'What is the max number of iterations?'
      read(stdin,*) maxbumps
      if(maxbumps .gt. MAXG)then
	write(stdout,*) 'Maximum Number of bumps is ',MAXG
	maxbumps = 200
      end if
      write(stdout,*)'What is the phase shift (secs) for the output?'
      read(stdin,*) theshift
      write(stdout,*)'What is minimum percent error increase to accept?'
      read(stdin,*) tol
      write(stdout,*)'What is is the Gaussian filter width factor?'
      read(stdin,*) gwidth
      write(stdout,*)'Allow negative pulses? (1->y, 0->no)'
      read(stdin,*)idum
      if(idum .eq. 1) then
	lpositive = .false.
      else
	lpositive = .true.
      end if 
c
      write(stdout,*)'Minimal (0) or verbose output(1)?'
      read(stdin,*)idum
      if(idum .eq. 0) then
	verbose = .false.
      else
	verbose = .true.
      end if 
c
******************************************************************************
c     
      call rsac1(numerator,f,npts,beg,delta,MAXPTS,nerr)
      if(nerr .ne. 0)then
	write(stdout,*)'Problem reading the numerator file'
	stop
      end if
c
      call rsac1(denominator,g,nptsd,b,dt,MAXPTS,nerr)
      if(nerr .ne. 0)then
	write(stdout,*)'Problem reading the denominator file'
	stop
      end if      
c
******************************************************************************
c
c     Find the next power of two greater than the data
c       dimensions - use the numerator, zero pad
c
      n = 1
119   continue
      if(n .ge. npts) go to 120
      n = n * 2
      go to 119
c      
120   continue
      if(n .gt. MAXPTS)then
        write(stdout,*) 'Too many points needed.'
        write(stdout,*) 'n = ', n
        stop
      end if
c
******************************************************************************
c     zero-pad the data
c
      npts = n
c
******************************************************************************
c     FINISHED READING FILES
c      
c     Now begin the cross-correlation procedure
c
c      Put the filter in the signals
c
      call gfilter(f,gwidth,npts,dt)
      call gfilter(g,gwidth,npts,dt)
      call wsac1('numerator',f,npts,beg,dt,nerr)
      call wsac1('observed',f,npts,beg,dt,nerr)
      call wsac1('denominator',g,npts,beg,dt,nerr)
c
c     compute the power in the "numerator" for error scaling
c
      power = 0
      do 100 i = 1, npts
        power = power + f(i)*f(i)
100   continue
c
c     correlate the signals
c 
      call fcorrelate(f,g,npts,MAXPTS,dt)
c     call wsac1('ccor0',g,npts,beg,dt,nerr)
c
c     find the peak in the correlation
c
      maxlag = npts/2
      write(stdout,'(/,a27,f10.5)') 'The maximum spike delay is ', 
     &   real(maxlag) * dt
c
      if(lpositive) then
	call getmax(g,maxlag,amps(1),shifts(1))
      else
	call getabsmax(g,maxlag,amps(1),shifts(1))
      end if
      amps(1) = amps(1) / dt
c
      nshifts = 1
c
c     read in the signals again
c
      call zero(f,MAXPTS)
      call zero(g,MAXPTS)
      call rsac1(numerator,f,ndummy,beg,delta,MAXPTS,nerr)
      call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
c
c     compute the predicted deconvolution result
c
      call zero(p,MAXPTS)
      call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
      if(verbose) then
          call phs_shift(p,theshift,npts,dt)      
        call wsac1('d001',p,npts,-theshift,dt,nerr)
          call phs_shift(p,-theshift,npts,dt)      
      end if
c
c     convolve the prediction with the denominator signal
c      
      call convolve(p,g,npts,dt)
c
      if(verbose) then
        call wsac1('p001',p,npts,beg,dt,nerr)
      end if
c
c     filter the signals
c     
      call gfilter(f,gwidth,npts,dt)
      call gfilter(g,gwidth,npts,dt)
c
      if(verbose)then
        write(resfile,'(a1,i3.3)') 'r',0
        call wsac1(resfile,f,npts,beg,dt,nerr)
      end if
c      
c     compute the residual (initial error is 1.0)
c
      call getres(f,p,npts,r,sumsq_ip1)
c
      sumsq_i = 1.0
      sumsq_ip1 = sumsq_ip1 / power
      d_error = 100*(sumsq_i - sumsq_ip1) 
c
      write(resfile,'(a1,i3.3)') 'r',1
      if(verbose)then
        call wsac1(resfile,r,npts,beg,dt,nerr)
      end if
c     
      write(stdout,1000)
      write(stdout,1001)
     &  resfile, dt*amps(1),(shifts(1)-1)*dt,100*sumsq_ip1,
     &  d_error
1000  format(/,1x,'File',9x,
     & 'Spike amplitude   Spike delay   Misfit   Improvement')
1001  format(1x,a10,2x,e16.9,2x,f10.3,3x,f7.2,'%',3x,f9.4,'%')
c
******************************************************************************
c    
      do while(d_error .gt. tol .and. nshifts .lt. (maxbumps))
c
        nshifts = nshifts + 1
	sumsq_i = sumsq_ip1
c
        call zero(g,MAXPTS)
        call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
        call gfilter(g,gwidth,npts,dt)
        call fcorrelate(r,g,npts,MAXPTS,dt)
	if(lpositive)then
	 call getmax(g,maxlag,amps(nshifts),shifts(nshifts))
        else
         call getabsmax(g,maxlag,amps(nshifts),shifts(nshifts))
        end if
        amps(nshifts) = amps(nshifts) / dt
c
        call zero(p,MAXPTS)
        call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
	if(verbose)then
          write(filename,'(a1,i3.3)') 'd',nshifts
          call phs_shift(p,theshift,npts,dt)      
	  call wsac1(filename,p,npts,-theshift,dt,nerr)
          call phs_shift(p,-theshift,npts,dt)      
        end if
c        
	call zero(g,MAXPTS)
        call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
        call convolve(p,g,npts,dt)
	if(verbose)then
          write(filename,'(a1,i3.3)') 'p',nshifts
	  call wsac1(filename,p,npts,beg,dt,nerr)
        end if
c                
        call zero(f,MAXPTS)
        call rsac1(numerator,f,ndummy,beg,delta,MAXPTS,nerr)
        call gfilter(f,gwidth,npts,dt)
        call getres(f,p,npts,r,sumsq_ip1)
        
        sumsq_ip1 = sumsq_ip1/ power
        write(resfile,'(a1,i3.3)') 'r',nshifts
	if(verbose)then
          call wsac1(resfile,r,npts,beg,dt,nerr)
	end if
	d_error = 100*(sumsq_i - sumsq_ip1)
	
	write(stdout,1001)
     &   resfile,dt*amps(nshifts),(shifts(nshifts)-1)*dt,
     &   100*sumsq_ip1,d_error
c    
      enddo
c
******************************************************************************
c      
      write(stdout,1010) d_error
1010  format(/,1x,'Last Error Change = ',f9.4,'%',/)
c
c     if the last change made no difference, drop it
c      
      fit = 100 - 100*sumsq_ip1
c
      if(d_error .le. tol)then
         nshifts = nshifts - 1
         fit = 100 - 100*sumsq_i
         write(stdout,*)'Hit the min improvement tolerance - halting.'
      end if
c
      if(nbumps .ge. maxbumps)then
         write(stdout,*)'Hit the max number of bumps - halting.'
      end if
c
      write(stdout,*)'Number of bumps in final result: ', nshifts
      write(stdout,1011) fit
1011  format(1x,'The final deconvolution reproduces ',
     &    f6.1,'% of the signal.',/)
c
******************************************************************************
c
c     compute the final prediction
c
      call zero(p,MAXPTS)
      call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
      call zero(g,MAXPTS)
      call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
      call convolve(p,g,npts,dt)
      call wsac1('predicted',p,npts,beg,dt,nerr)
      call zero(g,MAXPTS)
c
c     write out the answer
c
      call zero(p,MAXPTS)
      call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
      call phs_shift(p,theshift,npts,dt)      
c
      call newhdr
      call rsac1(numerator,g,ndummy,b,dt,MAXPTS,nerr)
      call setnhv('NPTS',npts,nerr)
      call setfhv('B',-theshift,nerr)
      theend = -thshift + (npts-1)*dt
      call setfhv('E',theend,nerr)
      call setnhv('NZSEC',-12345,nerr)
      call setfhv('USER0',gwidth,nerr)
c     call setkhv('KUSER0','Rftn',nerr)
c     call setkhv('KUSER1','IT_DECON',nerr)
      call wsac0('decon.out',xdummy,p,nerr)
c
c     write out the gaussian filter
c
      if(verbose)then
        call newhdr
        call zero(p,MAXPTS)
        p(1) = 1 / dt
        call phs_shift(p,theshift,npts,dt)      
        call gfilter(p,gwidth,npts,dt)
        call wsac1('thefilter',p,npts,beg,dt,nerr)
      end if
      
      stop
      end
*
******************************************************************************
******************************************************************************
*
*
*
*
************************************************************************
*
*     correlation routine - correlates f and g and replaces the 
*       g with the cross-correlation the value is normalized
*       by the zero-lag autocorrelation of g
*
************************************************************************
*
      subroutine fcorrelate(f,g,n,MAXPTS,dt)
      real f(MAXPTS), g(MAXPTS), c(8192)
      real sum0, temp
      integer i,n,n2,n2o2
c
c     compute the zero-lag autocorrelation of g
c
      sum0 = 0
      do 1 i = 1, n
        sum0 = sum0 + g(i)*g(i)
1     continue
      sum0 = sum0 * dt
c
c     compute the next power of 2 greater than n
c
      n2 = 1
5     n2 = n2 * 2
      if(n2 .lt. n) go to 5
c     
6     continue
      n2o2 = n2 / 2
c
c     Use the Numerical Recipes routine to compute the cross correlation
c
      call correl(f,g,n2,c)
c 
      temp = dt / sum0
c
      do 20 i = 1,n2
        g(i) = c(i) * temp
20    continue
      
      return
      end   
*
************************************************************************
*
*     zero a real array
*
************************************************************************
*
      subroutine zero(x,n)
      real x(n)
      integer i,n
      
      do 1 i = 1,n
        x(i) = 0
1     continue
      return
      end
*
************************************************************************
*
*     get max value of array and its index
*
************************************************************************
*
      subroutine getmax(x,n,maxvalue,maxindex)
      real x(n), maxvalue
      integer i,n,maxindex
      
      maxvalue = x(1)
      maxindex = 1
      do 20 i = 2, n
	if(x(i) .gt. maxvalue) then
	   maxvalue = x(i)
	   maxindex = i
        end if
20    continue

      return
      end
*
************************************************************************
*
*     find max absolute value of array and its index
*
************************************************************************
*
      subroutine getabsmax(x,n,thevalue,maxindex)
      real x(n), maxvalue, thevalue
      integer i,n,maxindex
      
      maxvalue = abs(x(1))
      maxindex = 1
      thevalue = x(1)
      do 20 i = 2, n
	if(abs(x(i)) .gt. maxvalue) then
	   maxvalue = abs(x(i))
	   thevalue = x(i)
	   maxindex = i
        end if
20    continue

      return
      end
*
*
************************************************************************
*
*     getres
*
************************************************************************
*
      subroutine getres(x,y,n,r,sumsq)
      real x(n), y(n), r(n), sumsq
      integer i,n
      
      sumsq = 0 
      do 20 i = 1, n
       r(i) = x(i) - y(i)
       sumsq = sumsq + r(i)*r(i) 
20    continue

      return
      end
*
************************************************************************
*
*     compute the predicted time series from a set of
*       amplitudes and shifts
*
************************************************************************
*
      subroutine build_decon(amps,shifts,nshifts,p,n,gwidth,dt)
      real p(n), amps(nshifts)
      integer shifts(nshifts)
      integer i, n, nshifts
      
      call zero(p,n)
      do 1 i = 1, nshifts
        p(shifts(i)) = p(shifts(i)) + amps(i)
1     continue

      call gfilter(p,gwidth,n,dt)

      return
      end
*
************************************************************************
*
*     convolve a function with a unit-area Gaussian filter.
*
************************************************************************
*
      subroutine gfilter(x,gwidth_factor,n,dt)
      real x(n), pi, two_pi, gauss, d_omega, omega
      real gwidth, gwidth_factor, sum
      integer i, j, n, n2, halfpts
      integer forward, inverse
c      
      forward = 1
      inverse = -1
      pi = acos(-1.0)
      two_pi = 2 * pi
      sum = 0
c      
      n2 = 1
1     n2 = n2 * 2
      if(n2 .ge. n) goto 2
      go to 1
c      
2     continue
      halfpts = n2 / 2
c
      call realft(x,halfpts,forward)
c
      df = 1 / (float(n2) * dt)
      d_omega = two_pi * df
      gwidth = 4.0*gwidth_factor*gwidth_factor
c 
c     Handle the nyquist frequency
c
      omega = two_pi/(2.0*dt)
      gauss = exp(-omega*omega / gwidth)
      x(2) = x(2) * gauss
c 
      do 5 i = 2, halfpts
          j = i*2
          omega = (i-1) * d_omega
          gauss = exp(-omega*omega / gwidth)
          x(j-1) = x(j-1) * gauss
          x(j)   = x(j)   * gauss
5     continue 
c
      call realft(x,halfpts,inverse)
c      
      scalefactor = dt * (2 * df)
      do 10 i = 1, n
        x(i) = x(i) * scalefactor
10    continue    
c
      return
      end


*
************************************************************************
*
*     replace x with x convolved with y
*
************************************************************************
*
      subroutine convolve(x,y,n,dt)
      real x(n), y(n), dt, scale0, scale1
      integer i, j, n, n2, halfpts
      integer forward, inverse, stderr
c      
      forward = 1
      inverse = -1
      stderr = 6
c      
      if(mod(n,2) .ne. 0) then
        write(stderr,*) 'Error in convolve - n .ne. power of two.'
        stop
      endif
      n2 = n
      halfpts = n2 / 2
c
      call realft(x,halfpts,forward)
      call realft(y,halfpts,forward)
c
      df = 1 / (float(n2) * dt)
c 
c     Handle the zero & nyquist frequency
c
      x(1) = x(1) * y(1)
      x(2) = x(2) * y(2)
c 
      do 5 i = 2, halfpts
          j = i*2
          a = x(j-1)
          b = x(j)
          c = y(j-1)
          d = y(j)
          x(j-1) = a * c - b * d
          x(j)   = a * d + b * c
5     continue  
c
      call realft(x,halfpts,inverse)
      call realft(y,halfpts,inverse)
c      
      scale0 = dt * (2 * df)
      scale1 = dt * scale0
      do 10 i = 1, n
        y(i) = y(i) * scale0
        x(i) = x(i) * scale1
10    continue    
c
      return
      end
*
************************************************************************
*
*     phase shifts a signal
*      
*
************************************************************************
*
      subroutine phs_shift(x,theshift,n,dt)
      real x(n), pi, two_pi, theshift, d_omega, omega
      integer i, j, n, n2, halfpts
      integer forward, inverse
c      
      forward = 1
      inverse = -1
      pi = acos(-1.0)
      two_pi = 2 * pi
c      
      n2 = 1
1     n2 = n2 * 2
      if(n2 .ge. n) goto 2
      go to 1
c      
2     continue
      halfpts = n2 / 2
c
      call realft(x,halfpts,forward)
c
      df = 1 / (float(n2) * dt)
      d_omega = two_pi * df
c 
c     Handle the nyquist frequency
c
      omega = two_pi/(2.0*dt)
      x(2) = x(2) * cos(omega * tshift)
c 
      do 5 i = 2, halfpts
          j = i*2
          omega = (i-1) * d_omega
          a = x(j-1)
          b = x(j)
          c = cos(omega*theshift)
          d = sin(omega*theshift)
          x(j-1) = a*c-b*d
          x(j)   = a*d+b*c
5     continue 
c
      call realft(x,halfpts,inverse)
c      
      scalefactor = dt * (2 * df)
      do 10 i = 1, n
        x(i) = x(i) * scalefactor
10    continue    
c
      return
      end

