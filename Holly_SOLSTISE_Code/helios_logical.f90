program helios_add_on2
   implicit none


  double precision omega, tcyc, theta, Elab, phi 
  double precision V, vper,vpar,t,x,y,z,a,c,d,f,g,l
  double precision pi,q,B,m,x1,x2,x3,y1,y2,y3,phi2,h
  double precision z1,z2,z3,E,zdif,Elab_2
  logical discrep,discrep2
  integer i,n,j,k,ios,ios2
  real,allocatable ::  arr(:)
  
     
  pi=3.141592653589793238462643383
  q=1.6E-19
  B=1.915
  m=1.6726219E-27

!.................Variables..................................

!omega..... angular frequency of the ejectile and is directly proportional to B

!tcyc...... cyclotron frequency of the ejectile (time for one full orbit)

!theta..... initial angle of the ejectile in the xy plane

!phi........initial angle of the ejectile in the z plane

!Elab.......initial energy of the ejectile in the lab frame, measured in MeV

!V..........initial velocity of the ejectile calculated from the intial energy Elab

!vper.......velocity perpendicular to the magnetic field 

!vpar.......velocity parallel to the magnteic field (constant)



!..............Files...........................................

!'out_eject_unsmeared_all.dat'.... VIKAR input files, has the ejectiles' theta and Elab values

!'all_particles_trajec.dat'.........write to this file the x,y,z trajectories without obstacle
!read from this file and use coordinates to calc trajectories that are blocked

!'E_z.dat'................ E vs z for all particles

!'all_particleswo_trajec.dat'...........write to this file the x,y,z trajectories with obstacle
!if these are different to those in helios_out_test3.dat then determines trajectories that are blocked

!'blocked_trajec.dat'...............write to this file the trajectories of blocked paths

!blocked_zproject.dat'..............final z position histogram for blocked particles

!'unblocked_trajec.dat'...........write to this file the tracjectories of unblocked paths

!unblocked_zproject.dat'...........final z position histogram for unblocked particles

!'z_phi.dat'..............write to this file the final z position and phi values for unblocked particles

!'xyz3.dat'...........write to this file the difference between the coordinates of blocked and 
!unblocked paths

!'phi.dat'..........write the phi values randomly generated for each row of VIKAR in order to be used
!in the DO loop with the obstacle without having to generate a new random value
write(6,*)''
write(6,*)''
write(6,*)''
write(6,*)' '
write(6,*)'########   ########  ##       ########  ########  ########  ########  ######## '
write(6,*)'##         ##    ##  ##       ##           ##        ##     ##        ##          '
write(6,*)'##         ##    ##  ##       ##           ##        ##     ##        ##          '
write(6,*)'########   ##    ##  ##       ########     ##        ##     ########  ########    '
write(6,*)'      ##   ##    ##  ##             ##     ##        ##           ##  ##          '
write(6,*)'      ##   ##    ##  ##             ##     ##        ##           ##  ##          '
write(6,*)'########   ########  #######  ########     ##      #######  ########  ########   '
write(6,*)'            Solenoid & Supersonic Target In Structure Experiments'                    
write(6,*)' '
write(6,*)''
write(6,*)'                          VIKAR add-on for HELIOS'
write(6,*)''
write(6,*)''


  
  open(unit=1,file='out_eject_unsmeared_all.dat',status='unknown') 
  open(unit=2,file='all_particles_trajec.dat',status='unknown')
  open(unit=3,file='all_particleswo_trajec.dat',status='unknown')
  open(unit=14,file='blocked_trajec.dat',status='unknown')
  open(unit=15,file='unblocked_trajec.dat',status='unknown') 
  open(unit=18,file='phi.dat',status='unknown')
  open(unit=46,file='zEphi_withoutob.dat',status='unknown')
  open(unit=47,file='zEphi_withob.dat',status='unknown')
 

!Ascertain how many VIKAR events the program should calculate trajectories for
! The VIKAR file has 10,000 events in total
  
 
 write(6,*) 'Enter the number of events to simulate (max 100,000)'
  read(5,*) n
  write(6,*) 'Un (or more) momento por favor...'
  call random_seed()
  
  omega=(q*B)/m
  tcyc=(2*pi)/omega

! Obtain the velocity and original theta value of the ejectiles from the VIKAR file
! Assign a random isotropic value of phi between -pi and pi


!..........................WITHOUT OBSTACLE.................................................. 

   outer: DO i=1,n
           read(1,*) theta, Elab
           theta=(theta/180)*pi
           allocate(arr(n))
           call random_number(arr(n))
           phi=360*arr(n)-180
           phi=(phi/180)*pi
           write(18,*) phi
            V=sqrt((2*Elab*1.6E-13)/m)
           vper=V*sin(theta)
           vpar=V*cos(theta)
           t=0

      
    inner1: DO while (t<=tcyc)
               x= (-(vper/omega)*cos((omega*t)+phi))+((vper/omega)*cos(phi))
               y= ((vper/omega)*sin((omega*t)+phi))-((vper/omega)*sin(phi))
               z=vpar*t
               t=t+1.0E-10
               if ( t>=tcyc) write(46,*) (vpar*t),Elab,phi
               write(2,*) x,y,z,Elab
            END DO inner1

    deallocate(arr)
   
    END DO outer
  close(18)
  close(1)
  close(46)
!.....................WITH OBSTACLE.............................................

   open(unit=25,file='out_eject_unsmeared_all.dat',status='unknown')
   open(unit=21,file='phi.dat',status='unknown')

   outer2: DO j=1,n
           read(25,*) theta, Elab
           theta=(theta/180)*pi
           read(21,*) phi
           V=sqrt((2*Elab*1.6E-13)/m)   !1.6E-13 converts from MeV to J
           vper=V*sin(theta)
           vpar=V*cos(theta)
           t=0

    inner2: DO while (t<=tcyc) 
               x= (-(vper/omega)*cos((omega*t)+phi))+((vper/omega)*cos(phi))
               y= ((vper/omega)*sin((omega*t)+phi))-((vper/omega)*sin(phi))
               z=vpar*t
               
!If x,y,z coordinates in range of the obstacle then sets vper to 0 so that they are different from now on from the first DO loop 
!without the obstacle present             
             
              t=t+1.0E-10
              if (x>=0.1 .and. x<=0.2 .and. y>=0.1 .and. y<=0.2 .and. z>=0.1 .and. z<=0.2) then
              vpar=0
              else if (t>=tcyc) then 
              write(47,*) (vpar*t),Elab,phi
              end if 
              write(3,*) x,y,z,t
            END DO inner2                 

    END DO outer2
  close(2)
  close(3)
  close(25)
  close(47)
!.......................THE TRAJECTORIES THAT NEVER HAPPENED.............................................

   open(unit=22,file='all_particles_trajec.dat',status='old')
   open(unit=23,file='all_particleswo_trajec.dat',status='old')
   open(unit=24,file='xyz3.dat',status='old')
   
  
!If the trajectories are different between the two DO loops then should write the coordinates of the tracjectory that should have happened
          outer3: DO
                      read(22,*,end=40)x1,y1,z1
                      read(23,*,end=50)x2,y2,z2                     
                      x3=abs(x1)-abs(x2)
                      y3=abs(y1)-abs(y2)
                      z3=abs(z1)-abs(z2)
                      discrep=(z3==0)
                      if (discrep .eqv. .true.) then 
                      write(15,*) x1,y1,z1 
                      else if (discrep .eqv. .false.) then
                      write(14,*) x1,y1,z1
                      end if  
                      END DO outer3
             40 close(22)
             50 close(23)           
          
              close(14)
              close(15)


!.......................Obtaining E vs z and z projection plots..........................................
 
  open(unit=48,file='zEphi_withoutob.dat',status='unknown')
  open(unit=51,file='zEphi_withob.dat',status='unknown')
  open(unit=49, file='E_z_no_obstacle.dat',status='unknown')
  open(unit=50,file='z_phi_no_obstacle.dat',status='unknown')
  open(unit=45,file='unblocked_zproject.dat',status='unknown')
  open(unit=36,file='blocked_zproject.dat',status='unknown')


  outer6: DO
                      read(48,*,end=65)z,Elab,phi  
                      read(51,*,end=66) z2,Elab_2,phi2
                      write(49,*) z,Elab
                      write(50,*) z,phi 
                      zdif=z-z2
                      discrep2=(zdif==0)
                      write(6,*) discrep2
                      if (discrep2 .eqv. .true.) then 
                      write(45,*) z,Elab
                      else if (discrep2 .eqv. .false.) then
                      write(36,*) z,Elab
                      end if
                      END DO outer6
            65 close(48)
            66 close(51)
               


write(6,*) 'Fin' 
end 



