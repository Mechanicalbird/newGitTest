Program CavityFlow02

!
!    Purpose:
!    To Solve 2-D Cavity Flow.
!
!    Record of revisions:
!        Date               Programer          Description of change
!        ====               =========          =====================
!      01/30/2013          M. S. Samara            Original code
!      02/15/2013          M. S. Samara        Olive Oile Cavity flow
!
!

Implicit none

Integer,PARAMETER :: nx = 300
Integer,PARAMETER :: ny = 300
Integer,PARAMETER :: nt = 400
Integer,PARAMETER :: nit = 1000


Integer :: iit
Integer :: it
Integer :: i
Integer :: j

Real, Dimension ( nx,ny ):: p
Real, Dimension ( nx,ny ):: b
Real, Dimension ( nx,ny ):: u
Real, Dimension ( nx,ny ):: un
Real, Dimension ( nx,ny ):: v
Real, Dimension ( nx,ny ):: vn
Real, Dimension ( nx,ny ):: pn
Real, Dimension ( nx ):: x
Real, Dimension ( ny ):: y

Real :: dx
Real :: dy
Real :: Xdis= 0.01
Real :: Ydis= 0.01
Real :: dt = 0.01
Real :: vis = 1.0125/(10**4)
Real :: rho = 800
real :: Re

dx = (Xdis)/(nx-1.0)
dy = (Ydis)/(ny-1.0)

open(Unit=10,file='OliveCavityFlowresults.PLT')


x(1) = 0.0
y(1) = 0.0

!connecting x with i  Values .......................................
Do i = 2, nx
  x(i) = x(i-1)+dx
End do
!...................................................................

!connecting y with j  Values .......................................
Do j = 2, ny
  y(j) = y(j-1)+dy
End do
!...................................................................

!initial condition..................................................
Do i = 1, nx
  Do j = 1, ny
  p(i,j)=0.0
  b(i,j)=0.0
  u(i,j)=0.0
  v(i,j)=0.0
  end do
end do
!...................................................................


!.............time...........................
Do it = 1, nt 
  Do i = 2, nx-1
   Do j = 2, ny-1
   b(i,j)=rho*((u(i+1,j)-u(i-1,j))/2/dx+(v(i,j+1)-v(i,j-1))/2/dy)/dt+((u(i+1,j)-u(i-1,j))/2/dx)**2&
   &+2*(u(i,j+1)-u(i,j-1))/2/dy*(v(i+1,j)-v(i,j-1))/2/dx+(((v(i,j+1)-v(i,j-1))/2/dy)**2)
   end Do
  end DO

Do iit = 1, nit
  pn = p
  Do i = 2, nx-1
    Do j = 2, ny-1
    p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy**2+(pn(i,j+1)+pn(i,j-1))*dx**2-b(i,j)*dx**2*dy**2)/(dx**2+dy**2)/2
    end Do
  end DO
!
Do j = 2, ny-1
p(1,j)= p(2,j)
end do
!
Do j = 2, ny-1
p(nx,j)= p(nx-1,j)
end do
!
Do i = 2, nx-1
p(i,1)= p(i,2)
end do
!
Do i = 2, nx-1
p(i,ny)= p(i,ny-1)
end do
!
end Do


un=u
vn=v
Do i = 2, nx-1
  Do j = 2, ny-1
    u(i,j)=un(i,j)-un(i,j)*dt/dx*(un(i,j)-un(i-1,j))-vn(i,j)*dt/dy*(un(i,j)-un(i,j-1))&
    &-1/rho*(p(i+1,j)-p(i-1,j))*dt/2/dx+vis*dt/dx**2*(un(i+1,j)-2*un(i,j)+un(i-1,j))&
    &+vis*dt/dy**2*(un(i,j+1)-2*un(i,j)+un(i,j-1))
    v(i,j)=vn(i,j)-un(i,j)*dt/dx*(vn(i,j)-vn(i-1,j))-vn(i,j)*dt/dy*(vn(i,j)&
    &-vn(i,j-1))-1/rho*(p(i,j+1)-p(i,j-1))*dt/2/dy&
    &+vis*dt/dx**2*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))&
    &+vis*dt/dy**2*(vn(i,j+1)-2*vn(i,j)+vn(i,j-1))
 end Do
end DO


!
Do j = 1, ny
u(1,j)=0.0
end do
!
Do j = 1, ny
u(nx,j)=0.0
end do
!
Do j = 1, ny
v(1,j)=0.0
end do
!
Do j = 1, ny
v(nx,j)=0.0
end do
!
Do i = 1, nx
u(i,1)=0.0
end do
!
Do i = 1, nx
u(i,ny)=1.0
end do
!
Do i = 1, nx
v(i,ny)=0.0
end do
!
Do i = 1, nx
v(i,1)=0.0
end do
!

  Write (10,*)'VARIABLES="X","Y","P","U","V"'
  Write(10,8800) it,it*dt,nx,ny
8800 format('ZONE T="STEP: ', I8,'", STRANDID=1, SOLUTIONTIME=', &
& E15.8,', I=',I5,',J=',I5,', F=POINT ')


Do i = 1, nx
 Do j = 1, ny
  Write (10,*) x(i),y(j),p(i,j),u(i,j),v(i,j)
 end Do
end DO
  Write (*,*)'time', it
end do
!..............end of the time loop........

Re = (u(1,ny)*Ydis)/vis
Write (*,*) 'Re  =  ' , Re 

Close (unit=10)

End program CavityFlow02
