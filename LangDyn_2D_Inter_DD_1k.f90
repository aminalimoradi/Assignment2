module globals
! Global variables
implicit none
integer :: n=1000   ! number of particles
real(8) :: L=1d0        ! box has dimensions L X L
integer :: nsec =4
double precision, parameter :: pi=2q0*asin(1q0) ! numerical constant
end module globals 

module Langevin
! Initialization and update rule for Langevin particles
use globals
implicit none
double precision :: dt,kT,g,m                   ! time step size and physical parameters
double precision :: pref1,pref2                 ! auxiliary parameters
double precision, allocatable, dimension(:) :: x,y,vx,vy,ax,ay,vhx,vhy
contains
subroutine set_parameters

! Set time step and physical parameters
dt=0.001d0
kT=1d0
g=1d0
m=1d0

! Set auxiliary parameters
pref1=g
pref2=sqrt(24d0*kT*g/dt)

end subroutine set_parameters
subroutine initialize_particles
integer :: i,inx,iny,m
double precision :: ran1,ran2,gr1,gr2,dx,dy
! Give particles initial position and velocity

m = ceiling(sqrt(real(n,8)))
dx = L/m
dy = L/m
do i=1,n
   inx=mod(i-1,m)
   iny=(i-1)/m
   x(i)=-L/2d0 + dx/2d0 + real(inx,8) * dx
   y(i)=-L/2d0 + dy/2d0 + real(iny,8) * dy
   ax(i)=0d0
   ay(i)=0d0
   call random_number(ran1)
   call random_number(ran2)
   gr1=sqrt(kT/(m))*sqrt(-2*log(ran1))*cos(2*pi*ran2) ! Box-Mueller transform
   gr2=sqrt(kT/(m))*sqrt(-2*log(ran1))*sin(2*pi*ran2)
   vx(i)=gr1
   vy(i)=gr2
end do

end subroutine initialize_particles
end module Langevin

module BC
  use globals
  implicit none
contains
  subroutine check_BC(x,y,vx,vy)
    integer :: i
    real(8) :: x(n),y(n),vx(n),vy(n)

    do i=1,n
       if(x(i) > L/2d0) then
          x(i) = L/2d0
          vx(i) = -vx(i)
       end if
       if(y(i) > L/2d0) then
          y(i) = L/2d0
          vy(i) = -vy(i)
       end if
       if(x(i) < -L/2d0) then
          x(i) = -L/2d0
          vx(i) = -vx(i)
       end if
       if(y(i) < -L/2d0) then
          y(i) = -L/2d0
          vy(i) = -vy(i)
       end if
    end do
    return
  end subroutine check_BC
end module BC

module Domain_decomposition
  use globals
  implicit none

contains
  include "sector_by_position.f90"
  include "particles_list_by_sector.f90"
  include "Sector_neighbours.f90"
  ! Here, we will include at least three functions: one the constructs a list of neighbours for each sector, one that find the dector index based on position and one that lists (or somehow specifies) all particles in a given sector.
end module Domain_decomposition

program main
use globals
use Langevin
use BC
use Domain_decomposition
use omp_lib
implicit none
integer :: i,j,s1,s2,p,q,p1,p2
integer, allocatable, dimension(:) :: S, p_list,q_list
integer, allocatable, dimension(:,:) :: neighbour_array
double precision, allocatable, dimension(:) :: ran1, ran2
double precision :: t,t_max,m1,m2
double precision :: wtime,begin,end,d,rc,sigma,eps,rx,ry,F

! Open files
open(11,file='trajectories')
open(12,file='means')

! Set interaction parameters
sigma=1d-3
eps=1d0
rc=sigma*2d0**(1d0/6d0)

! Allocate arrays
allocate(x(n),y(n),vx(n),vy(n),ax(n),ay(n),vhx(n),vhy(n),ran1(n),ran2(n),S(n),p_list(n+1),q_list(n+1),neighbour_array(0:nsec**2-1,0:9))
L = 1.0d0    ! L by L box centered at the origin.

call neighbour_array_generator(neighbour_array)
! Note that the last column of this array now holds the number of neighbours.

t=0d0
t_max=0.1d0


call set_parameters
call initialize_particles
call neighbour_array_generator(neighbour_array)
!begin = omp_get_wtime()
call cpu_time(begin)
!$omp parallel 
 do while(t.lt.t_max)
!$omp single
   vhx=vx+0.5*ax*dt
   vhy=vy+0.5*ay*dt
   x=x+vhx*dt
   y=y+vhy*dt

   call check_BC(x,y,vhx,vhy)
!$omp end single
      
   ax=0d0                   ! Add forces here if any
   ay=0d0                   ! Add forces here if any

   ! Before the computation of the total force, we need to compute the interaction forces:
! loop over sectors s1
   call neighbour_array_generator(neighbour_array)
!$omp do private(s1,s2,p1,p2,q,p,rx,ry,d,F)
! loop over sectors (s1)
   do s1=0,nsec**2-1
      call particle_list(x,y, s1, p_list)
      ! Note that the last column of p_list now holds the number of particles in the sector.
      ! loop over particles in s1 (p)   
      do p1=1,p_list(n+1)
         p=p_list(p1)
         ! loop over neighbours of s1 (including s1) (s2)
         do s2=neighbour_array(s1,0),neighbour_array(s1,neighbour_array(s1,9)-1)
            call particle_list(x,y, s2, q_list)
            ! loop over particles in neighbour s2 (q)
            do p2=1,q_list(n+1)
               q=q_list(p2)           
               if(p.ne.q) then
                  rx=x(q)-x(p)
                  ry=y(q)-y(p)
                  d=sqrt(rx**2 + ry**2)
                  if(d.lt.rc) then
                  	if (d .eq. 0d0) then
                  		d = 1*10E-5
                  	end if	
                     F=4d0*eps*( -12d0*sigma**12/d**13 + 6D0* sigma**6/d**7 )
                     ax(p)=ax(p)+F*rx/(d*m)
                     ay(p)=ay(p)+F*ry/(d*m)
                  end if
               end if
            end do
         end do
      end do
   end do
!$omp end do
   

!$omp single
   call random_number(ran1)
   ran1=ran1-0.5d0	
   call random_number(ran2)
   ran2=ran2-0.5d0
!$omp end single	
   
!$omp barrier
!$omp sections
   
!$omp section
   ax=ax-pref1*vhx+pref2*ran1
   ay=ay-pref1*vhy+pref2*ran2
      
   vx=vhx+0.5*ax*dt
   vy=vhy+0.5*ay*dt
 !$omp section
print*, real(n,8)
 write(12,*) t,sqrt(sum(x**2+y**2)/n)


   t=t+dt
!$omp end sections
   !Print *, 'hello3 from thread num=', omp_get_thread_num()
   do i=1,n
      write(11,*) x(i),y(i)
  end do
   ! write(11,*) ''
  
   
 end do
!$omp end parallel
!end = omp_get_wtime()
call cpu_time(end)

print *,'Wtime=',end-begin
 
! De-allocate arrays
deallocate(x,y,vx,vy,ax,ay,ran1,ran2,S,vhx,vhy,q_list,p_list)
! Close files
close(11)
close(12)

end program main
