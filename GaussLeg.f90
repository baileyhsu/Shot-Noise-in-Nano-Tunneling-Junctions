SUBROUTINE GaussLeg(x1,x2,n,x,w)
  Integer n
 Real*8 x1,x2,x(n),w(n)
 Real*8, PARAMETER :: EPS=3.d-14  !** EPS is the relative precision.

 Real*8 ::  p1,p2,p3,pp,xl,xm,z,z1,pi
  Integer :: i,j,m

  pi=acos(-1.d0)
                     !** High precision is a good idea for this routine.
  m=(n+1)/2          !** The roots are symmetric in the interval, so we
  xl=0.5*(x2-x1)     !** only have to do half of them.
  xm=0.5d0*(x2+x1)
  do i=1,m           !** Loop over the desired roots.
     z=cos(pi*(i-.25d0)/(n+.5d0))
     !** Starting with the above approximation to the ith root,
     !** we enter the main loop of refinement by Newton's method.
     z1 = z+10*EPS               !** Cheat to enter loop

     Do While (ABS(z-z1) > EPS)
        p1=1.0d0
        p2=0.0d0

        do j=1,n        !** Loop up the recurrence relation to get the Legendre
           p3=p2        !** polynomial evaluated at z.
           p2=p1
           p1=((2.0*j-1.d0)*z*p2-(j-1.0)*p3)/j
        end do

        !** p1 is now the desired Legendre polynomial. We next compute pp,
        !** its derivative, by a standard relation involving also p2,
        !** the polynomial of one lower order.
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp  !**Newton's method.
     End Do
     x(i)=xm-xl*z                   !** Scale the root to the desired interval,
     x(n+1-i)=xm+xl*z               !** and put in its symmetric counterpart.
     w(i)=2.0*xl/((1.0-z*z)*pp*pp)  !** Compute the weight
     w(n+1-i)=w(i)                  !**and its symmetric counterpart.
  end do
  return
End Subroutine  