CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 																C	 
C			   CYLINDER.FOR PROGRAMI							C
C																C
C	YARI-SILINDIR ÜZERINDE AKIS VE ISI TRANFERINI HESAPLAR		C		C
C             													C
C																C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       
     
     
      implicit real(a-h,o-z)
      real llo1,llo2,llo3,jac,alfa,beta,j2
	parameter(n1=80,n11=100,n12=0,m=40)
	parameter(n=n1+n11+n12)
	parameter(llo1=2.,hh=1.)
	parameter(dx=llo1/n1,dy=hh/m)
	parameter(ee=1.,r=.5)
	parameter(ik1=(ee-r)/dx,ik2=(llo1-ee-r)/dx,irm=2*r/dx)
  	parameter(errormax=0.01)
	parameter(eps=0.015,eps2=0.015,omega=1.8,nstepmax=10000)
	parameter(aa1=15.,cc1=5.)
	parameter(aa2=15.,cc2=5.)
	parameter(dt=0.01)
	parameter(re=100.,pr=0.71)
   	dimension x(n+1,m+1),y(n+1,m+1),xx(n+1,m+1),yy(n+1,m+1)
	*,a(n+1,m+1),b(n+1,m+1),
     *c(n+1,m+1),d(n+1,m+1),e(n+1,m+1),
     *eta(n+1,m+1),zeta(n+1,m+1),x_zeta(n+1,m+1),
     *x_eta(n+1,m+1),y_zeta(n+1,m+1),y_eta(n+1,m+1),x_zeta2(n+1,m+1),
     *y_zeta2(n+1,m+1),x_eta2(n+1,m+1),y_eta2(n+1,m+1),
     *x_etazeta(n+1,m+1),y_etazeta(n+1,m+1),jac(n+1,m+1),
     *psi(n+1,m+1),teta(100),psicon(n+1,m+1),error_psi(n+1,m+1),
	* error_xy(n+1,m+1),xcon(n+1,m+1),ycon(n+1,m+1),
	*p(n+1,m+1),q(n+1,m+1),om(n+1,m+1),u(n+1,m+1),v(n+1,m+1),
     *a_zeta(n+1),b_zeta(n+1),c_zeta(n+1),d_zeta(n+1),omj(n+1),
	*a_zeta2(n+1),b_zeta2(n+1),c_zeta2(n+1),d_zeta2(n+1),tempj(n+1),
	*omcon(n+1,m+1),error_om(n+1,m+1),omnew(n+1,m+1),tempnew(n+1,m+1),
	*a_eta(n+1),b_eta(n+1),c_eta(n+1),d_eta(n+1),omi(m+1),
	*a_eta2(n+1),b_eta2(n+1),c_eta2(n+1),d_eta2(n+1),tempi(m+1),
	*uc(n+1,m+1),vc(n+1,m+1),time(nstepmax),error_temp(n+1,m+1)
     *,ucon(nstepmax),ucon2(nstepmax),temp(n+1,m+1),tempcon(n+1,m+1)
	*,tcon(nstepmax),tcon2(nstepmax)

	goto 13
cccccccccc continuation  ccccccccccccccc
	open(781,file='re1000.dat')
    
      read(781,*)
      read(781,*)
    	  do J=1,m+1
	  do I=1,n+1
	read(781,*) x(i,j),y(i,j),psi(i,j),u(i,j),v(i,j),om(i,j),temp(i,j)
	  enddo
	  enddo
	read(781,*) nstep
     	print*, ' continuation ...psi_u_v.old  has been read.' 
	  close(781)
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

13	continue
	
c	boundary  conditions  for grids
		
	pi=acos(-1.)
	do j=1,m+1			 !
 	x(1,j)=0.0			 !left surface
	y(1,j)=(j-1)*dy 	 !
	enddo				 ! 		    

	do j= 1,m+1			 !
	x(n1+1,j)=llo1		 !right surface
	y(n1+1,j)=dy*(j-1)	 !
	enddo				 !
		 
	do i=2,n1			 !
	x(i,m+1)=(i-1)*dx	 !
	y(i,m+1)=hh			 !top surface
	enddo
     	dik1=(ee-r)/ik1					   !
 	do ii=1,ik1+1                      !
 	i=ii                               !
	x(i,1)=(ii-1)*dik1				   !
	y(i,1)=0.						   !
	enddo
								       !
	dteta=pi/irm					   !
	do ir=1,irm							!
     	i=ir+ik1+1							!  bottom surface
	teta(ir)=ir*dteta					!
	y(i,1)=r*sin(teta(ir))				!
	x(i,1)=ee-r*cos(teta(ir))			!
 	enddo
	dik2=(llo1-(ee+r))/ik2
  	do ii=1,ik2+1						 !
	i=ii+ik1+irm						 !
 	x(i,1)=ee+r+(ii-1)*dik2				 !
	y (i,1)=0.							 !
	enddo
  	
c	initial conditons  for grid distribution

	do i=2,n1
   	do j=2,m
      x(i,j)=x(i,1)
	y(i,j)=((y(i,m+1)-y(i,1))/m)*(j-1)+y(i,1)
	enddo
	enddo
 	  open(6,file='inigrid.dat')
 	  write(6,*)'VARIABLES = "X","Y"'
	  WRITE(6,*)'ZONE I=',n1+1,'J=',m+1,' F=POINT'
    	  do J=1,m+1
	  do I=1,n1+1
	  write(6,*) x(i,j),y(i,j)
	  enddo
	  enddo
	  close(6)

C	 GRÝD GENERATION

	deta=hh/m
	dzeta=llo1/n1
   	do i=1,n1+1
	do j=1,m+1
	zeta(i,j)=(i-1)*dzeta
	eta(i,j)=(j-1)*deta				           !for grid1
	enddo
	enddo

	do i=2,n1
 	do j=2,m
	q1=-aa1*(sign(1.0,(eta(i,j)-eta(i,1))))*
	+       exp(-cc1*abs(eta(i,j)-eta(i,1)))
 	q2=-aa2*(sign(1.0,(eta(i,j)-eta(i,m+1))))*
	+       exp(-cc2*abs(eta(i,j)-eta(i,m+1)))
	q(i,j)=q1+q2
 	enddo
	enddo
	p=0.
  	
  50	xcon=x
	ycon=y
 	do i=2,n1
      do j=2,m
 	x_zeta(i,j)=(x(i+1,j)-x(i-1,j))/(2*dzeta)
  	y_zeta(i,j)=(y(i+1,j)-y(i-1,j))/(2*dzeta)
	x_eta(i,j)=(x(i,j+1)-x(i,j-1))/(2*deta)
	y_eta(i,j)=(y(i,j+1)-y(i,j-1))/(2*deta)

 	  a(i,j)=(x_eta(i,j))**2+(y_eta(i,j))**2
	  c(i,j)=(x_zeta(i,j))**2.+(y_zeta(i,j))**2.
	  b(i,j)=x_zeta(i,j)*x_eta(i,j)+y_zeta(i,j)*y_eta(i,j)

	  jac(i,j)=1./(x_zeta(i,j)*y_eta(i,j)-y_zeta(i,j)*x_eta(i,j))

	flag4 = a(i,j)/dzeta/dzeta
 	flag5 = c(i,j)/deta/deta
	flag1 = p(i,j)/(2.*(jac(i,j)**2.)*dzeta)
	flag2 = q(i,j)/(2.*(jac(i,j)**2.)*deta)

	flagx3 = (x(i+1,j+1)-x(i+1,j-1)+x(i-1,j-1)-x(i-1,j+1))*
	+        b(i,j)/2.0/dzeta/deta
	
	x(i,j)= ((flag4+flag1)*x(i+1,j)+(flag4-flag1)*x(i-1,j)
     *	    +(flag5+flag2)*x(i,j+1)+(flag5-flag2)*x(i,j-1)
     *        -flagx3)/(2.*(flag4+flag5))

	flagy3 = (y(i+1,j+1)-y(i+1,j-1)+y(i-1,j-1)-y(i-1,j+1))*
 	+        b(i,j)/2.0/dzeta/deta
		
	 y(i,j)= ((flag4+flag1)*y(i+1,j)+(flag4-flag1)*y(i-1,j)
     *	    +(flag5+flag2)*y(i,j+1)+(flag5-flag2)*y(i,j-1)
     *        -flagy3)/(2.*(flag4+flag5))

       y(n1+1,j) = (4.0*y(n1,j) - y(n1-1,j))/3.0 !turev sinir sarti
 	 y(1,j) = (4.0*y(2,j) - y(3,j))/3.0
  	
	enddo
	enddo
 	error_xy=abs(x-xcon+y-ycon)

 	do i=1,n+1
	do j=1,m+1
	toperror_xy=toperror_xy+error_xy(i,j) 
	enddo
	enddo

	if(toperror_xy.gt.errormax)then
	toperror_xy=0.
	 
	goto 50
	else
	continue
	toperror_xy=0.
	endif

	dh2=x(2,1)-x(1,1)
	llo3=dh2*(((1.+eps)**n12)-1.)/eps

      do J=1,m+1
      do I=1,n1+1

	ii=i+n12
	xx(ii,j)=x(i,j)+llo3
 											!grid1
	yy(ii,j)=y(i,j)
	enddo
	enddo
	x=xx
	y=yy
	open(8,file='grid1.dat')
	  write(8,*)'VARIABLES = "X","Y"'
	  WRITE(8,*)'ZONE I=',n1+1,'J=',m+1,' F=POINT'
    	  do J=1,m+1
	  do I=n12+1,n1+n12+1
	  write(8,*) x(i,j),y(i,j)
	  enddo
	  enddo
	  close(8)

	dh=x(n1+n12+1,1)-x(n1+n12,1)
	llo2=dh*(((1.+eps)**n11)-1.)/eps

	do i=1,n11+1
	do j=1,m+1
	ii=i+n1+n12		                                    !grid 2
	x(ii,j)=llo3+llo1+dh*((1.+eps)**(i-1)-1.)/eps
	y(ii,j)=y(n1+n12+1,j)
	enddo
	enddo				                     ! 		    

	do i=1,n12+1
	do j=1,m+1
 	x(n12+2-i,j)=llo3-dh2*((1.+eps2)**(i-1)-1.)/eps2 !grid3
   	y(i,j)=y(n12+1,j)
	enddo
	enddo				                     ! 		    

         open(93,file='grid3.dat')
    	  write(93,*)'VARIABLES = "X","Y"'
	  WRITE(93,*)'ZONE I=',n12+1,'J=',m+1,' F=POINT'
   	  do J=1,m+1
	  do I=1,n12+1
	  write(93,*) x(i,j),y(i,j)
	  enddo
	  enddo
	  close(93)

       open(9,file='grid2.dat')
    	  write(9,*)'VARIABLES = "X","Y"'
	  WRITE(9,*)'ZONE I=',n11+1,'J=',m+1,' F=POINT'
   	  do J=1,m+1
     	   do I=n1+n12+1,n1+n11+n12+1
	  write(9,*) x(i,j),y(i,j)
	  enddo
	  enddo
	  close(9)

	  open(87,file='grid.dat')
	  write(87,*)'VARIABLES = "X","Y"'
	  WRITE(87,*)'ZONE I=',n+1,'J=',m+1,' F=POINT'
    	  do J=1,m+1
	  do I=1,n+1
	  write(87,*) x(i,j),y(i,j)
	  enddo
	  enddo
	  close(87)
	 
	do i=1,n+1
  	do j=1,m+1
	zeta(i,j)=(i-1)*dzeta
	eta(i,j)=(j-1)*deta				  !for all grid
	enddo
	enddo

c     COMPUTATIONAL DOMAIN

	open(18,file='compdomain.dat')
	  write(18,*)'VARIABLES = "zeta","eta"'
	  WRITE(18,*)'ZONE I=',n+1,'J=',m+1,' F=POINT'
    	  do J=1,m+1
	  do I=1,n+1
	  write(18,*) zeta(i,j),eta(i,j)
	  enddo
	  enddo
	  close(18)
	print*, 'grid is completed'
123	continue

c c 	 METRICS
c     interior domain

	do i=2,n
	do j=2,m
   	x_zeta(i,j)=(x(i+1,j)-x(i-1,j))/(2*dzeta)
	y_zeta(i,j)=(y(i+1,j)-y(i-1,j))/(2*dzeta)
	x_eta(i,j)=(x(i,j+1)-x(i,j-1))/(2*deta)
	y_eta(i,j)=(y(i,j+1)-y(i,j-1))/(2*deta)
   	x_etazeta(i,j)=(x(i+1,j+1)+x(i-1,j-1)-x(i+1,j-1)-x(i-1,j+1))
	*                /(4.*deta*dzeta)
   	y_etazeta(i,j)=(y(i+1,j+1)+y(i-1,j-1)-y(i+1,j-1)-y(i-1,j+1))
	*                /(4.*deta*dzeta)
   	x_zeta2(i,j)=(x(i+1,j)-2.*x(i,j)+x(i-1,j))/dzeta**2.
   	y_zeta2(i,j)=(y(i+1,j)-2.*y(i,j)+y(i-1,j))/dzeta**2.
    	x_eta2(i,j)=(x(i,j+1)-2.*x(i,j)+x(i,j-1))/deta**2.
      y_eta2(i,j)=(y(i,j+1)-2.*y(i,j)+y(i,j-1))/deta**2.
   	enddo
	enddo
   
c     METRICS ON BOUNDARIES

	do i=1,n
	x_zeta(i,1)=(x(i+1,1)-x(i,1))/dzeta
	y_zeta(i,1)=(y(i+1,1)-y(i,1))/dzeta
	x_eta(i,1)=(x(i,2)-x(i,1))/deta			     !BOTTOM
	y_eta(i,1)=(y(i,2)-y(i,1))/deta
	enddo

	do j=2,m+1
	x_zeta(1,j)=(x(2,j)-x(1,j))/dzeta
	y_zeta(1,j)=(y(2,j)-y(1,j))/dzeta
	x_eta(1,j)=(x(1,j)-x(1,j-1))/deta			 !left
	y_eta(1,j)=(y(1,j)-y(1,j-1))/deta
	enddo

	do i=2,n+1
	x_zeta(i,m+1)=(x(i,m+1)-x(i-1,m+1))/dzeta
	y_zeta(i,m+1)=(y(i,m+1)-y(i-1,m+1))/dzeta
	x_eta(i,m+1)=(x(i,m+1)-x(i,m))/deta			 !top
	y_eta(i,m+1)=(y(i,m+1)-y(i,m))/deta
	enddo

	do j=1,m
 	x_zeta(n+1,j)=(x(n+1,j)-x(n,j))/dzeta
 	y_zeta(n+1,j)=(y(n+1,j)-y(n,j))/dzeta
	x_eta(n+1,j)=(x(n+1,j+1)-x(n+1,j))/deta			 !right
	y_eta(n+1,j)=(y(n+1,j+1)-y(n+1,j))/deta
	enddo
   
s	secondary metrics on boundaries
   	do i=2,n
	x_zeta2(i,1)=(x(i+1,1)-2.*x(i,1)+x(i-1,1))/dzeta**2.
	y_zeta2(i,1)=(y(i+1,1)-2.*y(i,1)+y(i-1,1))/dzeta**2.
      enddo
     	do i=1,n+1  											     ! bottom
 	x_eta2(i,1)=(x(i,3)-2.*x(i,2)+x(i,1))/deta**2.
	y_eta2(i,1)=(y(i,3)-2.*y(i,2)+y(i,1))/deta**2.
	enddo

	do i=2,n
	x_zeta2(i,m+1)=(x(i+1,m+1)-2.*x(i,m+1)+x(i-1,m+1))/dzeta**2.
	y_zeta2(i,m+1)=(y(i+1,m+1)-2.*y(i,m+1)+y(i-1,m+1))/dzeta**2.
      enddo
     	do i=1,n+1  												  !top
 	x_eta2(i,m+1)=(x(i,m+1)-2.*x(i,m)+x(i,m-1))/deta**2.
	y_eta2(i,m+1)=(y(i,m+1)-2.*y(i,m)+y(i,m-1))/deta**2.
	enddo
   
	do j=1,m+1
	x_zeta2(1,j)=(x(3,j)-2.*x(2,j)+x(1,j))/dzeta**2.
	y_zeta2(1,j)=(y(3,j)-2.*y(2,j)+y(1,j))/dzeta**2.
	enddo
     	do j=2,m  												  !left
 	x_eta2(1,j)=(x(1,j+1)-2.*x(1,j)+x(1,j-1))/deta**2.
	y_eta2(1,j)=(y(1,j+1)-2.*y(1,j)+y(1,j-1))/deta**2.
	enddo
 
	do j=1,m+1
	x_zeta2(n+1,j)=(x(n+1,j)-2.*x(n,j)+x(n-1,j))/dzeta**2.
	y_zeta2(n+1,j)=(y(n+1,j)-2.*y(n,j)+y(n-1,j))/dzeta**2.
	enddo
     	do j=2,m  												  !right
 	x_eta2(n+1,j)=(x(n+1,j+1)-2.*x(n+1,j)+x(n+1,j-1))/deta**2.
	y_eta2(n+1,j)=(y(n+1,j+1)-2.*y(n+1,j)+y(n+1,j-1))/deta**2.
	enddo

c     secondary metrics on the corners
	x_zeta2(1,1)=(x(3,1)-2.*x(2,1)+x(1,1))/dzeta**2.
	y_zeta2(1,1)=(y(3,1)-2.*y(2,1)+y(1,1))/dzeta**2.
 	x_eta2(1,1)=(x(1,3)-2.*x(1,2)+x(1,1))/deta**2.
	y_eta2(1,1)=(y(1,3)-2.*y(1,2)+y(1,1))/deta**2.
 	x_zeta2(n+1,1)=(x(n+1,1)-2.*x(n,1)+x(n-1,1))/dzeta**2.
 	y_zeta2(n+1,1)=(y(n+1,1)-2.*y(n,1)+y(n-1,1))/dzeta**2.
 	x_eta2(n+1,1)=(x(n+1,3)-2.*x(n+1,2)+x(n+1,1))/deta**2.
	y_eta2(n+1,1)=(y(n+1,3)-2.*y(n+1,2)+y(n+1,1))/deta**2.
 	x_zeta2(n+1,m+1)=(x(n+1,m+1)-2.*x(n,m+1)+x(n-1,m+1))/dzeta**2.
 	y_zeta2(n+1,m+1)=(y(n+1,m+1)-2.*y(n,m+1)+y(n-1,m+1))/dzeta**2.
 	x_eta2(n+1,m+1)=(x(n+1,m+1)-2.*x(n+1,m)+x(n+1,m-1))/deta**2.
	y_eta2(n+1,m+1)=(y(n+1,m+1)-2.*y(n+1,m)+y(n+1,m-1))/deta**2.
 	x_zeta2(1,m+1)=(x(1,m+1)-2.*x(2,m+1)+x(3,m+1))/dzeta**2.
 	y_zeta2(1,m+1)=(y(1,m+1)-2.*y(2,m+1)+y(3,m+1))/dzeta**2.
 	x_eta2(1,m+1)=(x(1,m+1)-2.*x(1,m)+x(1,m-1))/deta**2.
	y_eta2(1,m+1)=(y(1,m+1)-2.*y(1,m)+y(1,m-1))/deta**2.

c	mixed metrics on boundaries
	do i=2,n
	x_etazeta(i,1)=(x(i+1,2)+x(i-1,1)-x(i-1,2)-x(i+1,1))
	*                /(2.*deta*dzeta)
     	y_etazeta(i,1)=(y(i+1,2)+y(i-1,1)-y(i-1,2)-y(i+1,1))
	*                /(2.*deta*dzeta)

 	x_etazeta(i,m+1)=(x(i+1,m+1)+x(i-1,m)-x(i-1,m+1)-x(i+1,m))
 	*                /(2.*deta*dzeta)
     	y_etazeta(i,m+1)=(y(i+1,m+1)+y(i-1,m)-y(i-1,m+1)-y(i+1,m))
	*                /(2.*deta*dzeta)
	enddo

	do j=2,m
	x_etazeta(1,j)=(x(2,j+1)+x(1,j-1)-x(1,j+1)-x(2,j-1))
	*                /(2.*deta*dzeta)
     	y_etazeta(1,j)=(y(2,j+1)+y(1,j-1)-y(1,j+1)-y(2,j-1))
	*                /(2.*deta*dzeta)
 	x_etazeta(n+1,j)=(x(n+1,j+1)+x(n,j-1)-x(n,j+1)-x(n+1,j-1))
 	*                /(2.*deta*dzeta)
     	y_etazeta(n+1,j)=(y(n+1,j+1)+y(n,j-1)-y(n,j+1)-y(n+1,j-1))
	*                /(2.*deta*dzeta)
	enddo 

c     mixed metrics on corners
	x_etazeta(1,1)=(x(2,2)+x(1,1)-x(1,2)-x(2,1))
	*                /(deta*dzeta)
     	y_etazeta(1,1)=(y(2,2)+y(1,1)-y(1,2)-y(2,1))
	*                /(deta*dzeta)
 	x_etazeta(n+1,1)=(x(n+1,2)+x(n,1)-x(n,2)-x(n+1,1))
 	*                /(deta*dzeta)
     	y_etazeta(n+1,1)=(y(n+1,2)+y(n,1)-y(n,2)-y(n+1,1))
	*                /(deta*dzeta)
	x_etazeta(1,m+1)=(x(2,m+1)+x(1,m)-x(1,m+1)-x(2,m))
	*                /(deta*dzeta)
     	y_etazeta(1,m+1)=(y(2,m+1)+y(1,m)-y(1,m+1)-y(2,m))
	*                /(deta*dzeta)
 	x_etazeta(n+1,m+1)=(x(n+1,m+1)+x(n,m)-x(n,m+1)-x(n+1,m))
 	*                /(deta*dzeta)
     	y_etazeta(n+1,m+1)=(y(n+1,m+1)+y(n,m)-y(n,m+1)-y(n+1,m))
	*                /(deta*dzeta)

c     jacobien and coeficients
	do i=1,n+1
	do j=1,m+1
	  a(i,j)=(x_eta(i,j))**2.+(y_eta(i,j))**2.
	  c(i,j)=(x_zeta(i,j))**2.+(y_zeta(i,j))**2.
	  b(i,j)=x_zeta(i,j)*x_eta(i,j)+y_zeta(i,j)*y_eta(i,j)

	  jac(i,j)=1./(x_zeta(i,j)*y_eta(i,j)-y_zeta(i,j)*x_eta(i,j))

	  alfa=a(i,j)*x_zeta2(i,j)-2.*b(i,j)*x_etazeta(i,j)+
	*           c(i,j)*x_eta2(i,j)
	  beta=a(i,j)*y_zeta2(i,j)-2.*b(i,j)*y_etazeta(i,j)+
	*           c(i,j)*y_eta2(i,j)

	  d(i,j)=jac(i,j)*(y_zeta(i,j)*alfa-x_zeta(i,j)*beta)
        e(i,j)=jac(i,j)*(x_eta(i,j)*beta-y_eta(i,j)*alfa)

	 
	enddo
	enddo

C	   SOLUTION OF STREAM FUNCTION EQ
c	BOUNDARY and INITIAL CONDÝTÝONS FOR psi
  	do i=2,n
	do j=2,m
cccccccccccccccc  ATTENTION!!! FOR CONTINUATION SET  CCCCC
	psi(i,j)=0.5
	enddo
	enddo

	do i=1,n+1  
	psi(i,m+1)=1.          !top and bottom surface
	psi(i,1)=0.
	enddo

	do j=1,m+1
	psi(1,j)=y(1,j)            !left surface
	enddo

	do j=1,m+1
	u(1,j)=1.				   !velocity boundaries
	v(1,j)=0
	v(n+1,j)=0 
	enddo

c     CALCULATION OF psi BY GAUSS SEÝDEL METHOD
300 	nstep=nstep+1		              !  main LOOP
      omcon=om
	tempcon=temp
100	psicon=psi                        !   AKIM FONK: LOOPU
	isayac = isayac + 1
	 
	  do i=2,n
	  do j=2,m
      flagpsi1=a(i,j)/(dzeta**2.)
	flagpsi2=e(i,j)/(2.*dzeta)
	flagpsi3=c(i,j)/(deta**2.)
	flagpsi4=d(i,j)/(2.*deta)
	flagpsi5=2.*b(i,j)/(4.*deta*dzeta)

	    psi(i,j)=
	*    (flagpsi1+flagpsi2)*psi(i+1,j)
	*   +(flagpsi1-flagpsi2)*psi(i-1,j)
	*   +(flagpsi3+flagpsi4)*psi(i,j+1)
     *   +(flagpsi3-flagpsi4)*psi(i,j-1)
	*   -flagpsi5*(psi(i+1,j+1)+psi(i-1,j-1)-psi(i+1,j-1)-psi(i-1,j+1))
     *            +om(i,j)/jac(i,j)**2
       psi(i,j)=psi(i,j)/(2.*flagpsi1+2.*flagpsi3)
       psi(N+1,j) = (4.0*psi(N,j) - psi(N-1,j))/3.0
	 
	enddo
	enddo

      error_psi=abs(psi-psicon)
	do i=1,n+1
	do j=1,m+1
	toperror_psi=toperror_psi+error_psi(i,j) 
	enddo
	enddo

	if(toperror_psi.gt.errormax)then
	toperror_psi=0.
	goto 100
	else
	print*, isayac
	isayac = 0.0
	continue
	endif
 
c	computation of velocities
  	do i=2,n
 	do j=2,m
 	flag1=(psi(i+1,j)-psi(i-1,j))
	flag2=(psi(i,j+1)-psi(i,j-1))
 	u(i,j)=jac(i,j)*(x_zeta(i,j)*flag2/(deta*2.)
	*                 -x_eta(i,j)*flag1/(2.*dzeta))
	v(i,j)=jac(i,j)*(y_zeta(i,j)*flag2/(deta*2.)
	*                 -y_eta(i,j)*flag1/(2.*dzeta))
 	uc(i,j)=jac(i,j)*(u(i,j)*y_eta(i,j)-v(i,j)*x_eta(i,j))
	vc(i,j)=jac(i,j)*(v(i,j)*x_zeta(i,j)-u(i,j)*y_zeta(i,j))
 	enddo
	enddo
 	do j=1,m+1
  	u(n+1,j)=u(n,j) 
      enddo	  
 
c	boundary conditions for vorticities
	do i=1,n+1
	om(i,1)=2*(jac(i,1)**2)*c(i,1)*(psi(i,1)-psi(i,2))/(deta**2.)	!bottom surface 
	enddo

	do i=1,n+1
    	om(i,m+1)=2*(jac(i,m+1)**2)*c(i,m+1)*(psi(i,m+1)-psi(i,m))
	*          /(deta**2.)                                           !top  surface	  
      enddo

	do j=2,m
    	om(1,j)=a(1,j)*2.*(psi(1,j)-psi(2,j))/(dzeta**2.)	                   !left inflow surface
	*		 +c(1,j)*(psi(1,j+1)-2.*psi(1,j)+psi(1,j-1))/(deta**2.)
     *		 +d(1,j)*(psi(1,j+1)-psi(1,j-1))/(2.*deta)
      om(1,j)=om(1,j)*(jac(1,j)**2.)
	enddo
     
     	do j=2,m							
    	om(n+1,j)=a(n+1,j)*2.*(psi(n+1,j)-psi(n,j))/(dzeta**2.)	                !right surface
	* -c(n+1,j)*(psi(n+1,j+1)-2.*psi(n+1,j)+psi(n+1,j-1))/(deta**2.)
     * -d(n+1,j)*(psi(n+1,j+1)-psi(n+1,j-1))/(2.*deta)
      om(n+1,j)=om(n+1,j)*(jac(n+1,j)**2.)
    	enddo

c     boundary conditions for temperature
    	do j=1,m+1
	temp(1,j)=1.	            !left
	enddo

	do i=1,n+1
	temp(i,1)=0                !top and bottom
	temp(i,m+1)=0
	enddo

	do ir=1,irm							!
     	i=ir+ik1+1							!  cylinder surface
	temp(i,1)=0.142				
	enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      EK BASLANGIC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	do j=1,m+1
	temp(n+1,j)=temp(n,j)          !right
	temp(n+1,j)=0.
	enddo

c	updating vorticities
c	derivation in zeta direction

	do j=2,m
	do i=2,n

	czeta=uc(i,j)*dt/dzeta
	ceta=vc(i,j)*dt/deta
	dizeta=dt/(re*(dzeta**2.))
	dieta=dt/(re*(deta**2.))
	ezeta=dt/(re*dzeta)
	eeta=dt/(re*deta)
	fetazeta=dt/(re*deta*dzeta)
	j2=jac(i,j)**2
      a_zeta(i)=-.5*(.5*czeta+j2*a(i,j)*dizeta-.5*j2*e(i,j)*ezeta)
      b_zeta(i)=1.+j2*a(i,j)*dizeta
      c_zeta(i)=.5*(.5*czeta-j2*a(i,j)*dizeta-.5*j2*e(i,j)*ezeta)

      d_zeta(i)=.5*(.5*ceta+j2*c(i,j) *dieta-
     *	           .5*j2*d(i,j)*eeta)*om(i,j-1)
     *          +(1.-j2*c(i,j)*dieta)*om(i,j)
     *	+.5*(-.5*ceta+j2*c(i,j)*dieta+.5*j2*d(i,j)*eeta)*om(i,j+1)
     *	-fetazeta*j2*b(i,j)*.25* (om(i+1,j+1)+om(i-1,j-1)-om(i-1,j+1)
     *	-om(i+1,j-1))

	dizeta2=dt/(re*pr*(dzeta**2.)) 
	dieta2=dt/(re*pr*(deta**2.))
	ezeta2=dt/(re*pr*dzeta)
	eeta2=dt/(re*pr*deta)
	fetazeta2=dt/(re*pr*deta*dzeta)

	a_zeta2(i)=-.5*(.5*czeta+j2*a(i,j)*dizeta2-.5*j2*e(i,j)*ezeta2)
	b_zeta2(i)=1.+j2*a(i,j)*dizeta2
      c_zeta2(i)=.5*(.5*czeta-j2*a(i,j)*dizeta2-.5*j2*e(i,j)+ezeta2)
	d_zeta2(i)=.5*(.5*ceta+j2*c(i,j)*dieta2-
     *	       .5*j2*d(i,j)*eeta2)*temp(i,j-1)
     *           +(1.-j2*c(i,j)*dieta2)*temp(i,j)
     *	+.5*(-.5*ceta+j2*c(i,j)*dieta2+.5*j2*d(i,j)*eeta2)*temp(i,j+1)
     *	-fetazeta2*j2*b(i,j)*.25*(temp(i+1,j+1)+temp(i-1,j-1)
     * -temp(i-1,j+1)-temp(i+1,j-1))
      enddo

c	appling boundaries

      d_zeta(2)=d_zeta(2) - a_zeta(2)*om(1,j)
      d_zeta (n)=d_zeta(n) - c_zeta(n)*om(n+1,j)
      d_zeta2(2)=d_zeta2(2)-a_zeta2(2)*temp(1,j)
      b_zeta2(n)=b_zeta2 (n)+c_zeta2(n)	                          !turev party
      d_zeta2(n)=d_zeta2(n)- c_zeta2(n)*temp(n+1,j)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call tridag (a_zeta,b_zeta,c_zeta,d_zeta,omj,n)

      call tridag (a_zeta2,b_zeta2,c_zeta2,d_zeta2,tempj,n)
    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	c

	do i=2,n
      omnew(i,j)=omj(i)
	tempnew(i,j)=tempj(i)
	enddo
	enddo

c      derivation in eta direction
      do i=2,n
	do j=2,m
      om(i,j) = omnew(i,j)
	temp(i,j)=tempnew(i,j)
	czeta=uc(i,j)*dt/dzeta
	ceta=vc(i,j)*dt/deta
	dizeta=dt/(re*(dzeta**2.))
	dieta=dt/(re*(deta**2.)) 
	ezeta=dt/(re*dzeta) 
	eeta=dt/(re*deta) 
	fetazeta=dt/(re*deta*dzeta)

	j2=jac(i,j)**2
      a_eta(j)=-0.5*(0.5*ceta+j2*c(i,j)*dieta-.5*j2*d(i,j)*eeta)
      b_eta(j)=1.+j2*c(i,j)*dieta
      c_eta(j)=0.5*(0.5*ceta-j2*c(i,j)*dieta-.5*j2*d(i,j)*eeta)
      d_eta(j)=.5*(.5*czeta+j2*a(i,j)*dizeta-
     *	.5*j2*e(i,j)+ezeta)*om(i-1,j)+(1.-j2*a(i,j)*dizeta)*om(i,j)
     *	+.5*(-.5*czeta+j2*a(i,j)*dizeta+.5*j2*e(i,j)*ezeta)*om(i+1,j)
     *	+fetazeta*j2*b(i,j)*.25*(om(i+1,j+1)+om(i-1,j-1)-om(i-1,j+1)
     *	-om(i+1,j-1))

	dizeta2=dt/(re*pr*(dzeta**2.))
	dieta2=dt/(re*pr*(deta**2.))
	ezeta2=dt/(re*pr*dzeta)
	eeta2=dt/(re*pr*deta)
	fetazeta2=dt/(re*pr*deta*dzeta)

	a_eta2(j)=0.5*(0.5*ceta+j2*c(i,j)*dieta2-.5*j2*d(i,j)*eeta2)
      b_eta2(j)=1.+j2*c(i,j)*dieta2
      c_eta2(j)=0.5*(0.5+ceta-j2*c(i,j)*dieta2-.5*j2*d(i,j)*eeta2)

	d_eta2(j)=.5*(.5*czeta+j2*a(i,j)*dizeta2-
     *.5*j2*e(i,j)*ezeta2)*temp(i-1,j)+(1.-j2*a(i,j)*dizeta2)*temp(i,j)
     *+.5*(-.5*czeta+j2*a(i,j)*dizeta2+.5*j2*e(i,j)*ezeta2)*temp(i+i,j)
     *+fetazeta2*j2*b(i,j)*.25*(temp(i+1,j+1)+temp(i-1,j-1)
     * -temp(i-1,j+1)-temp(i+1,j-1))
      enddo


c     appling boundaries

      d_eta(2)=d_eta(2)-a_eta(2)*om(i,1)
	d_eta(m)=d_eta(m)-c_eta(m)*om(i,m+1)

      d_eta2(2)=d_eta2(2)-a_eta2(2)*temp(i,1)
	d_eta2(m)=d_eta2(m)-c_eta2(m)*temp(i,m+1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       CALL tridag(a_eta,b_eta,c_eta,d_eta,omi,m)
	 CALL tridag(a_eta2,b_eta2,c_eta2,d_eta2,tempi,m)

CCCCCCCCCCCCCCCCGCCCCGCCCGCGCCCCCGCCCCCCCCCCCCCCCCCCCC

	do j=2,m
         om(i,j)=omi(j)
	   temp(i,j)=tempi(j)
	enddo
	enddo

      error_om=abs(om-omcon)

      error_temp=abs(temp-tempcon)

      do i=1,n+1
      do j=1,m+1
      toperror_om=toperror_om+error_om(i,j)
      toperror_temp=toperror_temp+error_temp(i,j)
      enddo
      enddo

c    if(toperror_om.gt.errormax)then

      print*, toperror_om,toperror_temp,nstep
      time(nstep)=nstep*dt
      ucon(nstep)=u(n1,m/2)
      ucon2(nstep) =u (n-20,m/4)
      tcon(nstep)=temp(n1,m/2)
      tcon2(nstep)=temp(n-20,m/4)

	

      if (nstep.lt.nstepmax) then 

      toper ror_om=0.
	toperror_temp=0.

      open (55,file='osilasyon_u-T.dat')
	do i=1,nstep
      write (55,*)time(i),ucon(i),ucon2(i),tcon(i)
	enddo
	close (55)


	open (7,file='ps i_u_v.dat')
      write (7,*) 'VARIABLES="X","Y","psi","u","v","om","T"'
	WRITE (7,*) 'ZONE I=',n+1,'J=',m+1,'F=POINT'
      do J=1,m+1
      do I=1,n+1
      write (7,*) x(i,j),y(i,j),psi(i,j),u(i,j),v(i,j),om(i,j),temp(i,j)
      enddo
      enddo
      write (7,*) nstep
	close (7)
	goto 300
	else
      continue
	endif
	stop
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridag (a,b,c,d,r,n)
	real beta(500),s(500),m(500)
	real a(500),b(500),c(500),d(500),r(500)
	beta(2)=b(2)
	s(2)=d(2)
      do  j=3,n
      m(j)=a(j)/beta(j-1)
      beta(j)=b(j)-(m(j)*c(j-1))
      if (beta(j).eq.0.) pause 'tridag failed'
      enddo
      r(n)=s(n)/beta(n)
      do  j=n-1,2,-1
      r(j)=(s(j)-c(j)*r(j+1))/beta(j)
      enddo
	return
      end
	
