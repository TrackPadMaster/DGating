	program DGating

! Most recent as of 6/6/19 and finished
! No new edits besides commenting and documentation

! Welcome to the second iteration of the DGating program
! The previous one is based on old versions of code
! Rather than go through every single line fixing things,
!	it's probably easier to just do this whole thing over

! So this is going to be based on DLattice rather than DBiharmonic
! If questions you have aren't answered here or the README, go ahead
!	and check the DLattice documentation
! I've learned quite a few lessons recently that are talked about there

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaring Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	implicit none
! This is the part where I regret every variable I've made
! Integers needed for loops
	integer :: i,j
! Other integers
	integer :: atom_num,tsteps,idum
! Now every other variable
! We have our inputs
	real*8 :: GammaP_in,U_in,dt_in,driveratio,A,B,pqratio,phi
! Then calculated values from those inputs
	real*8 :: GammaP,U0,dt,w_v,w_1,w_2,state
! We've got some convenience variables
	real*8 :: ninth, cos2z,pi
! We've got our big three we take in our loops
	real*8 :: p,z,t,p0,z0,t0,p_f,z_f,t_f
! And some other variables needed in the loops
	real*8 :: kick,deltaP,rand,Dfsn,Fd,jumprate,gasdev
! And some things very unique to the gating ratchet
	real*8 :: U0_t,Gamma_t
! Finally this whole mess of other things we need for data analysis
	real*8 :: dispmean,vmean,dispsqrd,diff,v_sum,diff_sum
	real*8 :: vmean2,diff2,v_sum2,diff_sum2,v_var,diff_var
	real*8 :: v_std,diff_std,v_err,diff_err,peclet,v_avg,diff_avg

! And then we need our MPI and random number jargon
	integer :: ierr,myid,numprocs,n
	integer,allocatable :: seed(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Establish all the MPI stuff
! This allows parallel processing, though it adds some work for later
	include 'mpif.h'

! Set up the variables we'll need to for the MPI stuff

! Now start the MPI function, give the ID of each process
	call MPI_INIT (ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading in Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading in from text files outside of the program
! This way we don't need to recompile every time we change variables

! All input files are in the 10-series, output in 20-series
! Every input variable is in its own row, with a label after
! Then we need to put each "read" operation on its own row as well
	open(11,file='input.dat')
	read(11,*) atom_num
	read(11,*) tsteps
	read(11,*) U_in
	read(11,*) GammaP_in
	read(11,*) phi
	read(11,*) dt_in
	read(11,*) driveratio
	read(11,*) A
	read(11,*) B
	read(11,*) pqratio

! Remember that UNLIKE MultiFreq, A/B and pqratio are NOT related


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random Number Generator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we need a random number generator
! If we just feed the processors the same seed, they'll produce the
! same results for the simulation. This is bad.
! So intead, we are going to create seeds for the number generator that
! is based on a combination of time and processor ID. That way, the
! simulation I run tomorrow is different from today, and each process
! will be running its own simulation, rather than copying.

! Special note here!
! When you ask for the random_seed generator like this, it calls from the OS
! That means it's a little dependent on how good that RNG actually is
! Here, I'm trusting the OS to give me good random numbers
! If this is in doubt, it needs to be replaced by something else in a similar manner

	call random_seed(size=n)
	allocate(seed(n))
	call random_seed(get=seed)

	call random_seed(put=seed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fixing Values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We have now read in the values we need, let's fix them up
! Following the instructions of where to put the 1/2's

	GammaP = GammaP_in * 0.5d0
	U0 = U_in * 0.5d0
! Then we also scale the time step based on the scattering rate
	dt = dt_in/GammaP
! And the driving frequencies, in a round-a-bout kind of way
! First is the vibrational frequency, w_v = 2sqrt(U0)*w_r
! Where w_r is still 1/2, so 
	w_v = sqrt(U0)
! Now we pull w_1 as a ratio from w_v
	w_1 = w_v * driveratio
! And finally w_2 is pulled from w_1
	w_2 = w_1 * pqratio

! Some other details we should establish here as well

! We'll start the initial state as 1, this doesn't really matter
	state = 1.d0
! We acutally don't re-establish this each time because again,
!	it doesn't really matter
! Rather than do the same division each time, name it a variable
	ninth = 1.d0/9.d0

! And let's also set up a couple variables we'll use later
	vmean = 0.d0
	dispmean = 0.d0
	dispsqrd = 0.d0
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LOOPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! There is no need for a potentials loop because we only use one
! Then we start with the atom loop
	do i = 1,atom_num
! Set initial conditions for our new baby atom
! the gasdev actually doesn't give us realistic values for the start
! Not really a big deal as long as you give it long enough inside of
!	the lattice
	   p = gasdev(idum)
	   z = gasdev(idum)
	   t = 0.d0
	   p0 = p
	   z0 = z
	   t0 = t
! And then initial state is just whatever state we had last

! Now we go into our time loop
	   do j = 1,tsteps
! Now things get a bit more interesting
! Not only do we have this time-dependent driving force
! We ALSO have a time-dependent potential and gamma
	      U0_t=U0*(1+A*dcos(w_1*t))
	      Gamma_t=GammaP*(1+A*dcos(w_1*t))
! Now everything else looks pretty similar
! Less characters is easier, just call it a variable
	      cos2z = dcos(2.d0*z)
! The diffusion term, now including a time-dependence
	      Dfsn=(0.1d0*ninth*Gamma_t)*(35.d0+(state*7.d0*cos2z))
! The atom recieves a kick from this diffusion
	      kick=((2.d0*Dfsn*dt)**(0.5d0))*gasdev(idum)
! The driving force acting here
	      Fd = B * dcos(w_2*t + phi)
! And then put all of it together for a momentum change
	      deltaP=kick+(state*U0_t*dsin(2.d0*z)*dt) + Fd*dt

! Now let's determine if we're going to jump or not
! First, the probability of that happening
	      jumprate = ninth*dt*Gamma_t*(1.d0+(state*cos2z))
! And then the random number we compare it to
	      call random_number(rand)
! 	and compare them
	      if (rand .lt. jumprate) then
!	did we jump?
	         state = -state
	      end if
! And now we update our running variables
	      p = p + deltaP
	      z = z + p*dt
	      t = t + dt

	   end do ! end of the time loop

! Now we're not super concerned with the temperature
! Let's instead compare start and end states
	   z_f = z - z0
	   p_f = p - p0
	   t_f = t - t0
! I don't know that I needed this anyways, but we have them
! Now let's look at our running averages
! We have three we're interested in
! First is the final displacement of the atom, averaged together
	   dispmean = dispmean + z_f/atom_num
! Second is the overarching velocity, just total displacement over time 
	   vmean = vmean + z_f/(t_f*atom_num)
! Last is square of displacement, used for diffusion later
	   dispsqrd = dispsqrd + (z_f * z_f)/atom_num
	   
	end do ! end of the atom loop
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finalizing Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! So now we're outside of our computation, time to finalize everything
! First is to calculate diffusion
	diff = (dispsqrd - (dispmean**2.d0))/(2.d0*t_f)
! Now we have a bunch of these variables that need putting together
! Specifically, we're looking towards vmean and diff
	call MPI_REDUCE(vmean,v_sum,1,MPI_DOUBLE_PRECISION,
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
! Remember that this '0' says this value ONLY goes to processor 0
	call MPI_REDUCE(diff,diff_sum,1,MPI_DOUBLE_PRECISION,
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
! We also want to figure out the squares of these values
! This is for the error calculation later
! Start by squaring those two values on EACH processor
	vmean2 = vmean * vmean
	diff2 = diff * diff
	call MPI_REDUCE(vmean2,v_sum2,1,MPI_DOUBLE_PRECISION,
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(diff2,diff_sum2,1,MPI_DOUBLE_PRECISION,
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)

! And sure, we've summed those values, but they're not averages until..
	v_avg = v_sum/numprocs
	diff_avg = diff_sum/numprocs

! And then we can calculate standard error from that
! Only do this once, so:
	if (myid.eq.0) then
! This is the variance in velocity
	   v_var = (v_sum2-(numprocs*v_sum*v_sum))/(numprocs-1)
	   diff_var=(diff_sum2-(numprocs*diff_sum*diff_sum))
     & /(numprocs-1)
! So then standard deviation is one more step away from variance
	   v_std = dsqrt(v_var)
	   diff_std = dsqrt(diff_var)
! And lastly we have the standard error we're looking for
	   v_err = v_std/sqrt(dble(numprocs))
	   diff_err = diff_std/sqrt(dble(numprocs))

! And the Peclet number?
	   pi = 3.141592653589793d0
	   peclet = v_avg*pi/diff_avg
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing and Debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Going through, checking these values, they both seem reasonable
	!print *,p,z,t
	!print *,p_f,z_f,t_f

! Next are the values directly derived from those
! Again, not zero, not immediate red flags
	!print *,dispmean,vmean,dispsqrd

! Now we start looking into the "Finalizing Data" section
! This isn't true! --> wow these were all zero
! Remember that after the "reduce" function, only processor 0 has
! 	the value stored
! Fixed! --> But diff is still zero for some reason
	!print *,diff,v_sum,diff_sum
! These values are NEARLY zero if you only use one or two atoms
! Make sure to have just a small pool to make sure these averages work

! Step by step
! And these seem fine
	!print *,v_sum,diff_sum,v_sum2,diff_sum2

! Towards the end here now
! Oh no
	!print *,v_var,diff_var,v_std,diff_std,v_err,diff_err

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data Writing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First, reminder that all ID's for inputs are 10's, 20's for output
! Matching the outputs of DBiharmonic and DMultiFreq
! There's actually quite a bit of information here
! Not all of it is picked up by DataCollector
! Consider changing this in the future?

! We only want this written once, so only do one processor
	if (myid.eq.0) then
	   open(unit=21,file='output.dat',status='replace')
! First is a list of the input variables and time spent
	   write(21,*) 'dt = ',dt
	   write(21,*) 'U0/Er = ',U_in
	   write(21,*) 'time spent = ',t_f
	   write(21,*) 'w_1 = ',w_1
	   write(21,*) 'w_2 = ',w_2
	   write(21,*) 'phi = ',phi
	   write(21,*) 'Amplitude(A) = ',A
	   write(21,*) 'Force(B) = ',B
	   write(21,*) 'Processors: ',numprocs
! And next is the data we actually want
	   write(21,*) 'Output Data:'
	   write(21,*) '<v> = ',v_avg
	   write(21,*) 'v_err = ',v_err
	   write(21,*) 'Diff = ',diff_avg
	   write(21,*) 'Diff_err = ',diff_err
	   write(21,*) 'Peclet = ',peclet
	end if

	call MPI_FINALIZE(ierr)
	end program


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Algorithms to generate random numbers(Ref. W.H.Press et Al.         !
! "NUMERICAL RECIPES".Cambridge U.P. 1986)                            !
! Function GASDEV generates a Gaussian Distribution                   !
! of zero mean and variance unity.(Box-Mueller formula)               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION GASDEV(IDUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      REAL*8 GASDEV,X(2)
      SAVE ISET,GSET
      DATA ISET/0/
      IF(ISET.EQ.0) THEN
 1       CALL RANDOM_NUMBER(X)
         V1=2.0d0*X(1)-1.d0
         V2=2.0d0*X(2)-1.d0
         R=V1**2+V2**2
      IF(R.GE.1.0d0) GO TO 1
         FAC=DSQRT(-2.0d0*DLOG(R)/R)
         GSET=V1*FAC
         GASDEV=V2*FAC
         ISET=1
      ELSE
         GASDEV=GSET
         ISET=0
      ENDIF
      RETURN
      END

