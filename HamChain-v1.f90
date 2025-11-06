module param
!real(8)::bz,w
!complex(8),allocatable,dimension(:,:)::ham
!real(8),allocatable,dimension(:)::eiro
!integer nn
end module param

use param
Implicit none
real(8)::x,z,pr,d,c,v,b,yc1,yc2,t,plz,pi
real(8),allocatable,dimension(:)::y,dydx,yout
integer i,j,k,l,m,n,don
complex(8),allocatable,dimension(:)::WORK
real(8),allocatable,dimension(:)::rwork,ra,ev,evf,exn,az,ax,eiro,smas,smenos
complex(8),allocatable,dimension(:,:)::ham,hoz,hox
complex(8)::Im,cc,cc0,cfree
integer::idum,mmm,inicial,Nt,it,Nreal,i1,i2,infor,LWORK,n1,n2
real(8)::je,uno,cero,ep,halfa,dt,time,alfa1,alfa2,ran1,rui,fac
real(8)::ab,bx,h,cNt,cNt0,norm,t1,t2,tt,bmier,bz1,bz2,aux,sumhoz,sumhox
real(8)::sumhop,sumhom,auxm,auxp
real(8)::xssg,yssg,zssg,Bz,g
integer,allocatable,dimension(:,:)::basis
integer,allocatable,dimension(:)::tbas1,tbas2,tbas3,tbas4,tbas5

External derivs
character(len = 1)::a1,a2

!!!!!!!!! Inputs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Write(*,*)'tamaño de la cadena, primera cifra, segunda cifra'! 
        read(*,*)n,a1,a2!,alfa1,alfa2

       ! Write(*,*)'paso y Nt'
       ! read(*,*)h,Nt
        Write(*,*)'B_z y B_x'
        read(*,*)Bz,Bx
 

        print*,a1,a2

!!!!!!!!! Decalracion de Variables y Allocateos importantes !!!!!!!
        don=2**n
        je = real(1,kind(8))
        uno = real(1,kind(8))
        cero = real(0,kind(8))
        ep= 0.4d0
        Im=(0.d0,1.d0)
        LWORK=4*don
        Nt=0
        allocate(WORK(LWORK),RWORK(3*don-2),exn(n))!,H(icon,icon))
        allocate(ham(don,don),eiro(don),ra(don),ev(Nt),evf(Nt))!voct(n,n))
        allocate(basis(don,n+1),hoz(don,don))
        allocate(hox(don,don))
        allocate(az(n),ax(n),smas(n),smenos(n))
        allocate(tbas1(n+1),tbas2(n+1),tbas3(n+1),tbas4(n+1))

       !  allocate (y(don),dydx(don),yout(don))
 

       g=2.d0
       az=g*Bz
       ax=g*Bx



        ev=0
        evf=0
        open(8,file='Autoenrgias.dat')!!! Para guardar data
        


!!! Armemos la base

!! Creamos el indice de excitacines que nos da la dimenison de 
!! el espacio de "i" excitaciones

!basis(1,:)=0
!basis(don,:)=1

!!!!!!!!!!!! Armo la Base !!!!!!!!!!!!!!!!!!!!!!!

!!!!! |m_s^1,m_s^2,........,m_s^N>



open(9,file='guide'//a1//a2//'.in',status='old')

basis=0

do i=1,n
exn(i)=fac(n)/(fac(n-i)*fac(i))
end do

do i=1,don
   do j=1,n
     read(9,*)basis(i,j)
   end do
end do

do i=1,don
 do j=1,n
  basis(i,n+1)=basis(i,n+1)+basis(i,j)*2**(n-j)
  end do
end do  




do i=1,don
print*,'conf=',i,'->',basis(i,1:n)
print*,'bin',basis(i,n+1)
print*,'---------------------------------------'
end do

!!!!!!!!! Calculo de los elementos de matriz!!!!!!!!!!!!!!!!!!!!!!!!!!



do i=1,don !!! Loop elemento de la base
do j=1,don!!! Loop elemento de la base pero solo LU


!!!! Campo en z, el Hamiltiniano es Diagonal ahi

    if (i.eq.j) then

      aux=0.d0
      sumhoz=0.d0

      do k=1,n

        if (basis(j,k).eq.1) then
             aux=az(k)/2.d0
        else
             aux=-az(k)/2.d0 
        end if


      sumhoz=sumhoz+aux

      end do

  hoz(i,j)=sumhoz

   end if
!!!!!!!!!!!!!!!!!! Termino de campo en Z listo

!!!! Campo en x, el Hamiltiniano [ax(k)/2]*(Sk+ + Sk-)

!!! sum_k=1^N (Sk+|j>+Sk-|j>)


       aux=0.d0
       sumhop=0.d0

 do k=1,n
   tbas1=basis(j,:)
   tbas2=basis(j,:)

     if (basis(j,k).eq.1) then
       auxp=0.d0
       auxm=1.d0
       tbas1(k)=0
         do i1=1,n
           tbas1(n+1)=tbas1(n+1)+tbas1(i1)*2**(n-i1)
         end do

         if (tbas1(n+1).eq.basis(i,n+1))then
           smenos(k)=auxm
         else
           smenos(k)=0.d0
         end if
      end if
      if (basis(j,k).eq.0) then
        auxp=1.d0
        auxm=0.d0
        tbas2(k)=1

          do i1=1,n
            tbas2(n+1)=tbas2(n+1)+tbas2(i1)*2**(n-i1)
          end do

          if (tbas2(n+1).eq.basis(i,n+1))then
            smas(k)=auxp
          else
            smas(k)=0.d0
          end if
       end if

        hox(i,j)=hox(i,j)+(ax(k)/2.d0)*(smas(k)+smenos(k))
 end do

!!!!!!!!!!!!!!!!!!! Termino de campo en X listo       


!!!!!!!!!!!!!!!!!!!!!!!!!!!! Terminos interacruantes!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! hint=\sum_k [(J_k/2)*(S_{k}^{+}*S_{k+1}^{-} + 
!!!!               S_{k}^{-}*S_{k+1}^{+})+J_k*S_{k}^{z}*S_{k+1}^{z})
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do k=1,n !!!!!!!! Empieza la suma

   tbas1=basis(j,:)
   tbas2=basis(j,:)
   tbas3=basis(j,:)
   tbas4=basis(j,:)

! S_{k}^{z}*S_{k+1}^{z} ----------- Es diagonal

   if (i.eq.j)


        if (basis(j,k).eq.1) then
             aux1=1.d0/2.d0
        else
             aux1=-1.d0/2.d0
        end if
        if (basis(j,k+1).eq.1) then
             aux2=1.d0/2.d0
        else
             aux2=-1.d0/2.d0
        end if


      sumhintz=sumhintz+jex(k)*aux1*aux2
    end if


    !!!! Vamos con el termino de S+S-


!!!! Termino 1 segun las notas    
!!!!  S_{k}^{+}*S_{k+1}^{-}

      if ((basis(j,k+1).eq.1).and.(basis(j,k).eq.0))) then
       auxm=1.d0
       auxp=1.d0
       tbas1(k+1)=0
       tbas1(k)=1

         do i1=1,n
           tbas1(n+1)=tbas1(n+1)+tbas1(i1)*2**(n-i1)
         end do
       
         if (tbas1(n+1).eq.basis(i,n+1))then
           sint1(k)=auxm*auxp
         else
           sint1(k)=0.d0
         end if

      end if

!!!! Termino 2 segun las notas    
!!!!  S_{k}^{-}*S_{k+1}^{+}

      if ((basis(j,k).eq.1).and.(basis(j,k+1).eq.0))) then
       auxm=1.d0
       auxp=1.d0
       tbas2(k)=0
       tbas2(k+1)=1

         do i1=1,n
           tbas2(n+1)=tbas2(n+1)+tbas2(i1)*2**(n-i1)
         end do

         if (tbas2(n+1).eq.basis(i,n+1))then
           sint2(k)=auxm*auxp
         else
           sint2(k)=0.d0
         end if

      end if



        hox(i,j)=hox(i,j)+(jex(k)/2.d0)*(sint1(k)+sint2(k))


end do   !!!!!!!!!! Muere la suma k para interaccion

hiz(i,j)=sumhintz
himp(i,j)=



  end do ! End Loop elemento de la base pero solo LU
end do ! End Loop elemento de la base

ham=hoz+hox+



do i=1,don
   do j=1,don
   print*,i,j,real(ham(i,j))
   end do
end do

!   do i=1,exn(i)

      
        

!!!!!!!!!!!!!!!!!!!! Diagonalizo H_0

        infor=0
        call ZHEEV('V','U',don,ham,don,eiro,WORK,LWORK,RWORK,INFOR)

        do i=1,don
         print*,i,eiro(i)
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!! Dynamics!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Hacemos primero la Evolucion Libre
     !   open(7,file='Evolfree.dat')

!        time=0
!        do it=0,Nt
!          time=it*h
!          cfree=0
!          do i=1,n
!             cfree=cfree+exp(-Im*eiro(i)*time)*ham(n,i)*ham(1,i)
!          end do
!        ! write(7,*)time,abs(cfree)**2

!         ev(it)=ev(it)+abs(cfree)**2

!         end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! Ahora comenzamos la parte critica --> Evol con driving OCT

!open(8,file='Evoldriving.dat')!!! Para guardar data

!open(9,file='pulse.dat',status='old') !! Leemos pulso OCT de SSG

        ! Condicion inicial --> Excitacion en el Sitio 1
!        do i=1,n
!        y(i)=dreal(ham(1,i))
!        y(i+n)=dimag(ham(1,i))
!        end do


!            do it=1,nt !!!!! Evoluciono usando RK4 Valentina
!            read(9,*)t,bz,bx,bx

!         call derivs(t,y,dydx)  


!         call rk4(y,dydx,don,t,h,yout,derivs)

!         y=yout

!            cNt=0;cNt0=0;cc=0;cc0=0

!              do i=1,n
!                cc=cc+(y(i)+Im*y(i+n))*dconjg(ham(n,i)) !! Target(Pob sitio N)
!                cc0=cc0+(y(i)+Im*y(i+n))*ham(1,i)       !! Pob sitio 1
!              end do

!           cNt=abs(cc)**2
!           cNt0=abs(cc0)**2

           !!calculo de la norma por las dudas
!           norm=0.d0
!           do i=1,n
!             norm=norm+abs(y(i)+Im*y(i+n))**2
!           end do
!           !write(8,*)t,cNt,cNt0,norm
!           evf(it)=evf(it)+CNt
!
!           end do
!
!           close(9)
!
!
!       end do
!
!       time=0
!        do it=0,Nt
!          time=it*h
!   
!       write(8,*)time,evf(it)/Nreal,ev(it)/Nreal
!
!        end do



end

!        subroutine derivs(x,y,dy)
!        use param
!        real(8)::dy(2*nn),y(2*nn)
!        real(8)::x 
!        integer n2,i
!        complex(8)::Im
!        n2=nn
!        Im=(0.d0,1.d0)
!       ! allocate(dy(n2),y(n2))

!        do i=1,nn

!        dy(i)=y(i+nn)*eiro(i)
!        dy(i+nn)=-y(i)*eiro(i)
!
!do k=1,nn

!         dy(i)=dy(i)+bz*dimag((y(k)+Im*y(k+nn))* &
!         &dconjg((ham(2,k)*ham(1,i))+dconjg(ham(2,i))*ham(1,k)))
!
!         dy(i+nn)=dy(i+nn)-bz*real((y(k)+Im*y(k+nn))* &
!         &(dconjg(ham(2,k)*ham(1,i))+dconjg(ham(2,i))*ham(1,k)))

!        end do



!        end do


        
        !return

!dy(1)=-0.5d0*(e0+e*sin(w*x))*y(2)-0.5d0*delta*y(4)
!dy(2)=0.5d0*(e0+e*sin(w*x))*y(1)+0.5d0*delta*y(3)
!dy(3)=-0.5d0*delta*y(2)+0.5d0*(e0+e*sin(w*x))*y(4)
!dy(4)=0.5d0*delta*y(1)-0.5d0*(e0+e*sin(w*x))*y(3)

!        end

         


!        subroutine rk4(y,dydx,n,x,h,yout,derivs)
!	implicit none
!	integer n,nmax,i
!	double precision h,x
!        real(8)::dydx(n),y(n),yout(n)
!	external derivs
!	parameter(nmax=80)
!	double precision h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax)
!	hh=h*0.5d0
!	h6=h/6d0
!	xh=x+hh
!	do i=1,n
!	   yt(i)=y(i)+hh*dydx(i)
!	enddo
!	call derivs(xh,yt,dyt)
!	do i=1,n
!	   yt(i)=y(i)+hh*dyt(i)
!	enddo
!	call derivs(xh,yt,dym)
!	do i=1,n
!	   yt(i)=y(i)+h*dym(i)
!	   dym(i)=dyt(i)+dym(i)
!	enddo
!	call derivs(x+h,yt,dyt)
!	do i=1,n
!	   yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2d0*dym(i))
!	enddo
!	return
!	end

!        FUNCTION ran1(idum)
!        INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
!        REAL(8):: ran1,AM,EPS,RNMX
!        PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
!        PARAMETER (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!        INTEGER j,k,iv(NTAB),iy
!        SAVE iv,iy
!        DATA iv /NTAB*0/, iy /0/
!       if (idum.le.0.or.iy.eq.0) then !Initialize.
!       idum=max(-idum,1) !Be sure to prevent idum = 0.
!       do j=NTAB+8,1,-1 !Load the shuffle table (after 8 warm-ups).
!       k=idum/IQ
!       idum=IA*(idum-k*IQ)-IR*k
!       if (idum.lt.0) idum=idum+IM
!       if (j.le.NTAB) iv(j)=idum
!       enddo
!        iy=iv(1)
!        endif
!        k=idum/IQ
!        idum=IA*(idum-k*IQ)-IR*k
!        if (idum.lt.0) idum=idum+IM
!        j=1+iy/NDIV !Will be in the range 1:NTAB.
!        iy=iv(j) !Output previously stored value and refill the shuffle
!       iv(j)=idum! ble.
!       ran1=min(AM*iy,RNMX)! Because users don’t expect endpoint values.
!       return
!        END

         FUNCTION fac(N)

         integer N,i,j,k
         real(8)::x,y,z,a,fac

         a=1.d0
         do i=1,n
         a=a*dfloat(i)
         end do
          
         fac=a

         return
         end
 
