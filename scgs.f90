PROGRAM Lid_driven_SCGS_vanka
    Implicit None
    Integer::i,j,k,NI,NJ,Ihalo,Jhalo,MAXIT,IT,l,t1,t2,ITC
    Integer,Parameter::NX=150,NY=150
    Real::U(-2:NX,-2:NY),V(-2:NX,-2:NY),P(-2:NX,-2:NY),X(-2:NX,-2:NY),Y(-2:NX,-2:NY),APU(0:1), &
    APV(0:1),AEU(0:1),AEV(0:1),AWU(0:1),AWV(0:1),ANU(0:1),ANV(0:1),ASU(0:1),ASV(0:1),BU(0:1),BV(0:1),&
    DX,DY,REN,URFU,URFV,URFP,RESORU,RESORV,RESORM,SORMAX,RESU(-2:NX,-2:NY),RESV(-2:NX,-2:NY),RESM(-2:NX,-2:NY), &
    UC(-2:NX,-2:NY),VC(-2:NX,-2:NY),XC,YC,a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5,x1,x2,x3,x4,x5,nu,ts,te
    OPEN(Unit=10,File='result.dat')
    OPEN(Unit=11,File='iteration.dat')
    OPEN(Unit=13,File='gihacomparison.dat')
    OPEN(Unit=133,File='gihacomparison2.dat')
    OPen(unit=14,File='timereport.dat')
    Print*,'Enter max iteration= '
    Read(*,*)MAXIT
    Print*,'  Renolds number= '
    Read(*,*)REN
    Print*,' Enter Grid size '
    Read(*,*)NI,NJ
call cpu_time(ts)
    nu=1.0/REN
    SORMAX=0.001
    URFU=.3
    URFV=URFU
    Ihalo=2
    Jhalo=2
    DX=1./Float(NI)
    DY=1./Float(NJ)
    print*,DX,DY,nu
!Grid generation
    Do i=-Ihalo,NI+Ihalo
        Do j=-Jhalo,NJ+Jhalo
            X(i,j)=DX*i
            Y(i,j)=DY*j
        end do
    end do
! initial Data set up
    Do i=-Ihalo,NI+Ihalo
        Do j=-Jhalo,NJ+Jhalo
            U(i,j)=0.0
            V(i,j)=0.0
            P(i,j)=0.0
        end do
    end do

!SCGS start
ITC=1
 109   DO IT=ITC,MAXIT
      RESORU=0.
      RESORV=0.
      RESORM=0.

!block SOR
      DO i=1,NI
      DO j=1,NJ

! for u(i  ,j)
      AEU(1)=max(-.5*(U(i,j)+U(i+1,j))*DY,0.0)+(nu*(DY/DX))
      AWU(1)=max(.5*(U(i,j)+U(i-1,j))*DY,0.0)+(nu*(DY/DX))
      ANU(1)=max(-.5*(V(i,j)+V(i+1,j))*DX,0.0)+(nu*(DX/DY))
      ASU(1)=max(.5*(V(i,j-1)+V(i+1,j-1))*DX,0.0)+(nu*(DX/DY))
      APU(1)=max(.5*(U(i,j)+U(i+1,j))*DY,0.0)+max(-.5*(U(i,j)+U(i-1,j))*DY,0.0)+max(.5*(V(i,j)+V(i+1,j))*DX,0.0)+ &
      max(-.5*(V(i,j-1)+V(i+1,j-1))*DX,0.0)+2*(nu*(DY/DX))+2*(nu*(DX/DY))
       BU(1)=AEU(1)*U(i+1,j)+AWU(1)*U(i-1,j)+ANU(1)*U(i,j+1)+ASU(1)*U(i,j-1)+((P(i,j)-P(i+1,j))*DY)-APU(1)*U(i,j)

     ! pause
!for u(i-1,j)

      AEU(0)=max(-.5*(U(i-1,j)+U(i,j))*DY,0.0)+(nu*(DY/DX))
      AWU(0)=max(.5*(U(i-1,j)+U(i-2,j))*DY,0.0)+(nu*(DY/DX))
      ANU(0)=max(-.5*(V(i-1,j)+V(i,j))*DX,0.0)+(nu*(DX/DY))
      ASU(0)=max(.5*(V(i-1,j-1)+V(i,j-1))*DX,0.0)+(nu*(DX/DY))
      APU(0)=max(.5*(U(i-1,j)+U(i,j))*DY,0.0)+max(-.5*(U(i-1,j)+U(i-2,j))*DY,0.0)+max(.5*(V(i-1,j)+V(i,j))*DX,0.0)+ &
      max(-.5*(V(i-1,j-1)+V(i,j-1))*DX,0.0)+2*(nu*(DY/DX))+2*(nu*(DX/DY))
       BU(0)=AEU(1)*U(i,j)+AWU(1)*U(i-2,j)+ANU(1)*U(i-1,j+1)+ASU(1)*U(i-1,j-1)+((P(i-1,j)-P(i,j))*DY)-APU(0)*U(i-1,j)

! for v(i,j  )

      AEV(1)=max(-.5*(U(i,j)+U(i,j+1))*DY,0.0)+(nu*(DY/DX))
      AWV(1)=max(.5*(U(i-1,j+1)+U(i-1,j))*DY,0.0)+(nu*(DY/DX))
      ANV(1)=max(-.5*(V(i,j)+V(i,j+1))*DX,0.0)+(nu*(DX/DY))
      ASV(1)=max(.5*(V(i,j)+V(i,j-1))*DX,0.0)+(nu*(DX/DY))
      APV(1)=max(.5*(U(i,j)+U(i,j+1))*DY,0.0)+max(-.5*(U(i-1,j+1)+U(i-1,j))*DY,0.0)+max(.5*(V(i,j)+V(i,j+1))*DX,0.0)+ &
      max(-.5*(V(i,j)+V(i,j-1))*DX,0.0)+2*(nu*(DY/DX))+2*(nu*(DX/DY))
      BV(1)=AEV(1)*V(i+1,j)+AWV(1)*V(i-1,j)+ANV(1)*V(i,j+1)+ASV(1)*V(i,j-1)+((P(i,j)-P(i,j+1))*DX)-APV(1)*V(i,j)


! for v(i,j-1)

      AEV(0)=max(-.5*(U(i,j-1)+U(i,j))*DY,0.0)+(nu*(DY/DX))
      AWV(0)=max(.5*(U(i-1,j)+U(i-1,j-1))*DY,0.0)+(nu*(DY/DX))
      ANV(0)=max(-.5*(V(i,j-1)+V(i,j))*DX,0.0)+(nu*(DX/DY))
      ASV(0)=max(.5*(V(i,j-1)+V(i,j-2))*DX,0.0)+(nu*(DX/DY))
      APV(0)=max(.5*(U(i,j-1)+U(i,j))*DY,0.0)+max(-.5*(U(i-1,j)+U(i-1,j-1))*DY,0.0)+max(.5*(V(i,j-1)+V(i,j))*DX,0.0)+ &
      max(-.5*(V(i,j-1)+V(i,j-2))*DX,0.0)+2*(nu*(DY/DX))+2*(nu*(DX/DY))
       BV(0)=AEV(0)*V(i+1,j-1)+AWV(0)*V(i-1,j-1)+ANV(0)*V(i,j)+ASV(0)*V(i,j-2)+((P(i,j-1)-P(i,j))*DX)-APV(0)*V(i,j-1)


      a11=APU(0)
      a15=DY
      a22=APU(1)
      a25=-DY
      a33=APV(0)
      a35=DX
      a44=APV(1)
      a45=-DX
! under-relaxation
      a11=a11/URFU
      a22=a22/URFU
      a33=a33/URFV
      a44=a44/URFV

!for p(i,j)

      a51=-DY
      a52=DY
      a53=-DX
      a54=DX

      b1=BU(0)
      b2=BU(1)
      b3=BV(0)
      b4=BV(1)
      b5=-((U(i,j)-U(i-1,j))*DY+(V(i,j)-V(i,j-1))*DX)

! Boundary Conditions for u', v'
        If (i==1)then
          a51=0.0
          a11=1.0
          a15=0.0
          b1=0.0
          END if
        If (i==NI)then
          a52=0.0
          a22=1.0
          a25=0.0
          b2=0.0
          END if
        If (j==1)then
            a53=0.0
            a33=1.0
            a35=0.0
            b3=0.0
          END if
        If(j==NJ) Then
            a54=0.0
            a44=1.0
            a45=0.0
            b4=0.0
            end if
      CALL VANKA(a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5,x1,x2,x3,x4,x5)
! correcting U,V and P
      U(i-1,j)=U(i-1,j)+x1
      U(i,j)=U(i,j)+x2
      V(i,j-1)=V(i,j-1)+x3
      V(i,j)=V(i,j)+x4
      P(i,j)=P(i,j)+x5
! calculate residuals at centroid of P-CV
      RESORU=RESORU+(abs(b1)+abs(b2))/2.
      RESORV=RESORV+(abs(b3)+abs(b4))/2.
      RESORM=RESORM+abs(b5)

      END DO
      END DO

! specify BCs for U and V
    do k=1,NI-1
        U(k,NJ+1)=2.0-U(k,NJ)
        U(k,0)=-U(k,1)
        end do
        do k=1,NJ+1
        U(0,k)=0.0
        U(NI+1,k)=0.0
        end do
        do l=1,NJ
            V(l,0)=0.0
            V(l,NJ+1)=0.0
            end do
            do l=1,NJ
            V(0,l)=-V(1,l)
            V(NI+1,l)=-V(NI,l)
            end do
        PRINT*,IT,RESORU,RESORV,RESORM
        Write(11,*)IT,RESORU,RESORV,RESORM

      IF(MAX(RESORU,RESORV,RESORM).LE.SORMAX.AND.IT.GT.1) GO TO 1000
      IF(IT==MAXIT) Then
        Write(*,*)'Did not converge to run more place the number of iteration or 0 to abort'
        Read(*,*)MAXIT
        If (MAXIT==0) Go to 1000
        ITC=IT+1
        MAXIT=MAXIT+IT
        Go to 109
        END IF
      END DO

1000 Do i=1,NI
        Do j=1,NJ
            UC(i,j)=.5*(U(i,j)+U(i-1,j))
            VC(i,j)=.5*(V(i,j)+V(i-1,j))
            End Do
            End Do
   Write(10,*)'Variables= "X" , "Y" , "U" , "V", "P"'
   Write(10,*)'Zone F=Point, I=',NI,',J= ',NJ
   Do j=1,NJ
        Do i=1,NI
            XC=0.5*DX+(i-1)*DX
            YC=0.5*DY+(j-1)*DY
            Write(10,*)XC,YC,UC(i,j),VC(i,j),P(i,j)
           End do
        End do
         Write(13,*)'Variables= "V" , "X"'
        do j=0,NJ
            XC=j*DX
            Write(13,*)VC(j,NJ/2),XC
            end do
            Write(133,*)'Variables= "U" , "Y"'
             do i=0,NI
            YC=i*DY
            Write(133,*)UC(NI/2,i),YC
            end do
    Call cpu_time(te)
Write(14,98)NI,NJ,REN
98 Format(1X,'At grid size ',I4,' * ',I4,' And Reynolds number = 'f10.4 )
 Write(14,99)IT,abs(te-ts)
 99 Format(1X,'After 'I7,' total time required = ',F10.4,'s')
 END PROGRAM Lid_driven_SCGS_vanka

      Subroutine VANKA(a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5,x1,x2,x3,x4,x5)
      Implicit None
      Real::a11,a15,a22,a25,a33,a35,a44,a45,a51,a52,a53,a54,b1,b2,b3,b4,b5,x1,x2,x3,x4,x5,r1,r2,r3,r4,DEN
      r1=a51/a11
      r2=a52/a22
      r3=a53/a33
      r4=a54/a44
      DEN=r1*a15+r2*a25+r3*a35+r4*a45
      x5=(r1*b1+r2*b2+r3*b3+r4*b4-b5)/DEN
      x1=(b1-a15*x5)/a11
      x2=(b2-a25*x5)/a22
      x3=(b3-a35*x5)/a33
      x4=(b4-a45*x5)/a44
      RETURN
      END


