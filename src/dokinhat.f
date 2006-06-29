C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C subroutine dokinhat originated from dokhat, adapted for inhomogeneous
C K function by Pingping Zheng
      
c-------------------------------------------------------------------------
      function plarea(xp,yp,np)
c-------------------------------------------------------------------------
c
c find the area of the polygon defined by points in xp,yp
c
      implicit real*8 (a-h,o-z)

      dimension xp(np+1),yp(np+1)

      totare=0

      do is=1,np

        x1=xp(is)
        y1=yp(is)

        if(is.eq.np)then
          x2=xp(1)
          y2=yp(1)
        else
          x2=xp(is+1)
          y2=yp(is+1)
        end if

c Find the area of the trapezium
        totare = totare + (x2-x1)*(y2+y1)/2.0

      end do

c return a positive value
      plarea = abs(totare)

      end

CCCCCCCCCCCCCCCCCCCCCC

      function iplace(s,ns,t)
c
c which of the variable width bins s is t in?
c
      implicit real*8 (a-h,o-z)

      dimension s(ns)

      do ib=1,ns
        if(s(ib).ge.t)then
          iplace=ib
          return
        end if
      end do
c
c if it is outside the range of s
c
      iplace=ns+1

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sort2(x,n)                                                  
c                                                                           
c     shellsort algorithm                                                   
c     n     : number of elements to be sorted                               
c     x     : on enter an array of dimension at least n containing          
c             real numbers                                                  
c             on output first n elements of x are sorted from smallest      
c             to largest                                                    
c                                                                           
      implicit real*8 (a-h,o-z)                                             
      dimension x(n)                                                        
      i=1                                                                   
    1 i=i+1                                                                 
      if (i.le.n) goto 1                                                    
      m=i-1                                                                 
    2 m=m/2                                                                 
      if (m.eq.0) return                                                    
      k=n-m                                                                 
      do 4 j=1,k                                                            
      kk=j                                                                  
    3 if (kk.lt.1) goto 4                                                   
      if (x(kk+m).ge.x(kk)) goto 4                                          
      w=x(kk+m)                                                             
      x(kk+m)=x(kk)                                                         
      x(kk)=w                                                               
      kk=kk-m                                                               
      goto 3                                                                
    4 continue                                                              
      goto 2                                                                
      end   
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function ipippa(x,y,xc,yc,nc)
c 
c point in polygon routine.
c
c returns 0 if point x,y not in the bound polygon defined by xc,yc
c
c fortran version of C routine by Ken McElvain
c

      implicit real*8 (a-h,o-z)
      common /bounds/area,iconvx
      real*8 area
      integer iconvx


      dimension xc(nc+1),yc(nc+1)

        iwind = 0
        xlastp = xc(nc)
        ylastp = yc(nc)
        ioldq = iquad(xlastp,ylastp,x,y)
        do i=1,nc 
c for each point in the polygon 
                xthisp=xc(i)
                ythisp=yc(i)
                inewq = iquad(xthisp,ythisp,x,y)
                if(ioldq.ne.inewq) then
                        if(mod(ioldq+1,4).eq.inewq) then
                          iwind=iwind+1
                        else if(mod(inewq+1,4).eq.ioldq) then
                          iwind = iwind - 1
                        else 
                          a = (ylastp-ythisp)*(x-xlastp)
                          b = xlastp-xthisp
                          a = a + ylastp * b
                          b=b*y
                             if (a.gt.b) then
                               iwind=iwind+2
                             else
                               iwind=iwind-2
                             end if
                        end if
                end if
                xlastp=xthisp
                ylastp=ythisp
                ioldq=inewq
      end do
c 
c quadrant winding is either -4,0,+4 so divide down and take abs.
c
      ipippa = abs(iwind/4)

      end

      function iquad(xp,yp,xo,yo)
c
c determine which quadrant xp,yp is in relative to xo,yo as origin
c
      implicit real*8 (a-h,o-z)

        if(xp.lt.xo)then
                if(yp.lt.yo) then
                   iquad=2
                else 
                   iquad=1
                end if
        else
                if(yp.lt.yo)then
                   iquad = 3
                else 
                   iquad = 0
                end if
        end if

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function cncvwt(x,y,r,xp,yp,np)
c
c compute the weight given to a point at x,y according to how much
c of a circle of radius r is inside the bounding polygon
c 
c
      implicit real*8 (a-h,o-z)
      common /bounds/area,iconvx
      real*8 area
      integer iconvx
      dimension xp(np+1),yp(np+1)
      parameter(pi=3.141592654d0)
c store circle/poly intersections here
      parameter(maxcrs=40)
      dimension cross(maxcrs+1)
      parameter(tiny=1.0e-7)
c set count of crossing points to zero
      ncross = 0

c first loop over the boundary and find the crossing points
      do ic=1,np
c work with the trial point at origin
        x1=xp(ic)-x
        y1=yp(ic)-y
        x2=xp(ic+1)-x
        y2=yp(ic+1)-y

        cx=x2-x1
        cy=y2-y1
 
c these are the coefficients of the quadratic giving the intercept of
c line and circle.
        a=cx*cx+cy*cy
        b=2*(x1*cx+y1*cy)
        c=x1*x1+y1*y1-r*r

c find out if real solutions exist...
        b2m4ac=b*b-4*a*c

c ... and if they do, find them.
        if (b2m4ac.ge.0) then
          t1=(-b+sqrt(b2m4ac))/(2*a)
          t2=(-b-sqrt(b2m4ac))/(2*a)

c see if the solutions lie in the line segments
          if ((t1.gt.tiny).and.(t1-1.0.le.tiny)) then
             ncross=ncross+1
c find the angle to this point on thecircle
             ctemp=atan2(y1+t1*cy,x1+t1*cx)
             if(ctemp.lt.0)ctemp=2*pi+ctemp
             cross(ncross)=ctemp
c check crossing of circle with vertex
           else if (abs(t1).le.tiny)then
c compare this polygon segment's direction with that of the previous one
             nprev = (mod((ic+ (np-2)),np)+1)
             x0 = xp(nprev) - x
             y0 = yp(nprev) - y
             idp1 = isig8((x2-x1)*x1+ (y2-y1)*y1,tiny)
             idp2 = isig8((x1-x0)*x1+ (y1-y0)*y1,tiny)
c see if the polygon passes through the circle here
             if ((idp1-idp2).ne.1 .and.
     +        abs(idp1+idp2).ne.2) then
               ncross = ncross + 1
               ctemp = atan2(y1+t1*cy,x1+t1*cx)
               if (ctemp.lt.0.0) ctemp = 2*pi + ctemp
               cross(ncross) = ctemp
             end if
          end if

          if ((t2.gt.tiny).and.(t2-1.0.lt.tiny)) then
             ncross=ncross+1
             ctemp=atan2(y1+t2*cy,x1+t2*cx)
             if(ctemp.lt.0)ctemp=2*pi+ctemp
             cross(ncross)=ctemp
c check crossing of circle with vertex
           else if (abs(t2).le.tiny)then
c compare this polygon segment's direction with that of the previous one
             nprev = (mod((ic+ (np-2)),np)+1)
             x0 = xp(nprev) - x
             y0 = yp(nprev) - y
             idp1 = isig8((x2-x1)*x1+ (y2-y1)*y1,tiny)
             idp2 = isig8((x1-x0)*x1+ (y1-y0)*y1,tiny)
c see if the polygon passes through the circle here
             if ((idp1-idp2).ne.1 .and.
     +        abs(idp1+idp2).ne.2) then
               ncross = ncross + 1
               ctemp = atan2(y1+t2*cy,x1+t2*cx)
               if (ctemp.lt.0.0) ctemp = 2*pi + ctemp
               cross(ncross) = ctemp
             end if
          end if
        end if
      end do

c now we have all the crossing point angles stored in cross(1:ncross)

c if ncross = 0 then the total angle within the poly is 2*pi unless the
c circle is large and spans the polygon. this should be checked
c beforehand so it's okay to assume 2*pi here.

      if (ncross.eq.0) then
        totang=2*pi
      else

c sort into ascending order
        call sort2(cross,ncross)

c fix the ncross+1'th element to be the first plus 2pi so that the
c   list is circular...
        cross(ncross+1)=cross(1)+2*pi

c check that the number of crossings is even - if not then error.
        if (mod(ncross,2).ne.0) then
          cncvwt=-1
          return
        end if
c now find a nice spot to do the point-in-poly search
        sepmax=0.0
        icm=0

        do ic=1,ncross
          if (cross(ic+1)-cross(ic).gt.sepmax) then
            sepmax=cross(ic+1)-cross(ic)
            icm=ic
          end if
        end do
  
c icm is now the index of the crossing with the largest gap between it
c and the next crossing point

c test for point in poly of the point on the circle between these points
        angtes=(cross(icm)+cross(icm+1))/2.

        xtest=x+r*cos(angtes)
        ytest=y+r*sin(angtes)

c find out if test point is in the polygon boundary
        linpol=ipippa(xtest,ytest,xp,yp,np)

c find the total angle between (odd-even) crossings (i.e. 1-2 + 3-4 + ...
        totang = 0.
        do ic=1,ncross-1,2
          totang = totang + (cross(ic+1)-cross(ic))
        end do

c If the point we tested for p-i-p was on an odd-even
c section and was in the poly, then totang is the amount of circle inside the
c polygon. if the point was outside the polygon, then we need to subtract
c totang from 2*pi radians to get the angle inside the polygon. conversely,
c if the point tested was between even-odd crossings and outside the polygon,
c then totang is the angle we want, and if inside the polygon then again we
c have to do 2*pi-totang

        if ( (((mod(icm,2).eq.1).and.(linpol.eq.0))  .or.
     &        ((mod(icm,2).eq.0).and.(linpol.eq.1)) ) ) then
          totang = 2*pi-totang
        end if

      end if
c now totang is the angle contained in the polygon

c weight is proportion of total angle in the poly
      cncvwt = (2*pi)/(totang)
      return
      end

      integer function isig8(value,tiny)
c return the sign (+1,0,-1) of a value
            real*8 tiny,value
          if (value.gt.tiny) then
            isig8 = 1
          else if (value.lt.-tiny) then
            isig8 = -1
          else
            isig8 = 0
          end if
        return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function weight(x,y,r,xp,yp,np)
c
c find the weight for the point at x,y, radius r
c
      implicit real*8 (a-h,o-z)

      common /bounds/area,iconvx
      real*8 area
      integer iconvx

      dimension xp(np+1),yp(np+1)

      weight=cncvwt(x,y,r,xp,yp,np)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCC

C Copyright Barry Rowlingson <b.rowlingson@lancaster.ac.uk> and 
C Peter Diggle (c) 1991-3; http://www.maths.lancs.ac.uk/~rowlings/Splancs/
C R port: copyright 1998-2002 by Roger S. Bivand
C x,y pts, z=lambda(x,y), xp, yp polygon

      subroutine trykh(x,y,n,z,xp,yp,np,s,ns,hkhat,icounts,hkhats,nptns)

      implicit real*8(a-h,o-z)
      real*8 wij, wji

      common /bounds/area,iconvx
      real*8 area
      integer iconvx

      dimension x(n),y(n),xp(np+1),yp(np+1),s(ns),hkhat(ns)
      dimension icounts(nptns), hkhats(nptns)

      area = plarea(xp,yp,np)

      pi=3.141592654d0

      tmax=(s(ns))**2  

      do i=1,ns
        hkhat(i)=0.0d0
      end do


      do i=2,n
        i1=i-1                                                                
        xi=x(i)                                                              
        yi=y(i)                                                              
        do  j=1,i1                                                           
          xj=xi-x(j)                                                           
          yj=yi-y(j)                                                           
          t=xj*xj+yj*yj                                                         
          if (t.lt.tmax) then                                                 

            t=dsqrt(t)
            it=iplace(s,ns,t) 
 
            if(it.le.ns) then
              wij=weight(xi,yi,t,xp,yp,np)
              wji=weight(x(j),y(j),t,xp,yp,np)

              hkhat(it)=hkhat(it)+wij+wji
              ipos=n*(it-1)
              hkhats(i+ipos)=hkhats(i+ipos) + wij
              hkhats(j+ipos)=hkhats(j+ipos) + wji
              icounts(i+ipos)=icounts(i+ipos) + 1
              icounts(j+ipos)=icounts(j+ipos) + 1
            end if
          end if
        end do
      end do

      do i=2,ns
        hkhat(i)=hkhat(i)+hkhat(i-1)
        do j=1,n
          jj=j+(n*(i-1))
          jj1=j+(n*(i-2))
          hkhats(jj)=hkhats(jj)+hkhats(jj1)
        end do
      end do

      dn=dfloat(n)*dfloat(n-1)
      adn=area/dn

      do i=1,ns                                              
        hkhat(i)=hkhat(i)*adn
        do j=1,n
          jj=j+(n*(i-1))
          hkhats(jj)=hkhats(jj)*adn
        end do
      end do

      return                                                                
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C x, y pts coordinates, z lambda at(x,y), xp, yp polygon vertices,
C khat calculated at s[ns]

      subroutine dokinhat(x,y,n,z,xp,yp,np,s,ns,hkhat)

      implicit real*8(a-h,o-z)
      double precision wij, wji

      common /bounds/area,iconvx
      real*8 area
      integer iconvx

      dimension x(n),y(n),z(n),xp(np+1),yp(np+1),s(ns),hkhat(ns)

      area = plarea(xp,yp,np)

      pi=3.141592654d0

      tmax=(s(ns))**2

      do i=1,ns
        hkhat(i)=0.0d0
      end do

      do i=2,n
        i1=i-1
        xi=x(i)
        yi=y(i)
        do  j=1,i1
          xj=xi-x(j)
          yj=yi-y(j)
          t=xj*xj+yj*yj
          if (t.lt.tmax) then

            t=dsqrt(t)
            it=iplace(s,ns,t)

            if(it.le.ns) then
              wij=weight(xi,yi,t,xp,yp,np)
              wji=weight(x(j),y(j),t,xp,yp,np)
              hkhat(it)=hkhat(it)+(wij+wji)/(z(i)*z(j))
            end if
          end if
        end do
      end do

      do i=2,ns
        hkhat(i)=hkhat(i)+hkhat(i-1)
      end do

      dn=dfloat(n)*dfloat(n-1)

      do i=1,ns                                              
        hkhat(i)=hkhat(i)/area
      end do

      return                                                                
      end

