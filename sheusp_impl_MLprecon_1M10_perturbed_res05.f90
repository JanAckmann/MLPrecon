PROGRAM IS1SPSL
use implicit_functions_DP

implicit none

INTEGER, PARAMETER :: NLON=128/2 +2,NLAT=32
INTEGER, PARAMETER :: N=NLON,M=NLAT,NM=N*M
 character(len=150) :: EXP_NAME, Dp_depth_str

double precision :: U(N,M,0:1),&
               & V(N,M,0:1),  &
               & PD(N,M),     &
               & QX(N,M),     &
               & QY(N,M),     &
               & PT(N,M),     & 
               & P0(N,M),     &
               & HX(N,M),     &
               & HY(N,M),     &
               & DHX2Y(N,M),  &
               & S(N,M),      &
               & F1(N,M,0:1), &
               & F2(N,M,0:1), & 
               & E1(N,M,-1:0),&
               & E2(N,M,-1:0),&
               & UA(N,M),     &
               & VA(N,M+1),   &
               & PC(N,M),     &
               & X(N),        &
               & Y(M+1),        &
               & COR(N,M),    &
               & DIFF_U(N,M), &
               & DIFF_V(N,M), &
               & Velfx(N,M),  &
               & Velfy(N,M)

double precision :: MGH1IHX(M),  &
               & MGH2IHY(M),  &
               & AC(M),       &
               & BC(M),       &
               & AD(M),       &
               & BD(M), A, B, C, D
double precision ::U_52(N,M),&
                 & V_52(N,M),     &
                 & HX_52(N,M),     &
                 & HY_52(N,M),     &
                 & DHX2Y_52(N,M),  &
                 & S_52(N,M),      &
                 & X_52(N),        &
                 & Y_52(M+1),        &
                 & COR_52(N,M)

double precision :: TIME,&
               & DX, &
               & DY, &
               & DT, &
               & mue

double precision :: DX_52, &
               & DY_52, &
               & DT_52

INTEGER :: IP(N)
double precision :: sum_time, sum_lp_time

!!! store old values in low prec. to calc final tend.

double precision ::    PD_T(N,M), &
                  & QX_T(N,M), &
                  & QY_T(N,M)


                  
!!! the HIGH-PRECISION VARIABLES
double precision ::    PD_HP(N,M), &
                  & QX_HP(N,M), &
                  & QY_HP(N,M), &
                  & UA_HP(N,M), &
                  & VA_HP(N,M), &
                  & PT_HP(N,M), &
                  & P0_HP(N,M)
                  
                  
! CORIOLIS, GRAVITY AND EARTH RADIUS SPECIFICATION
INTEGER, PARAMETER  :: ICORIO=1
double precision :: F0
double precision :: G 
double precision :: R

! CHARACTERISTICS OF THE FLOW
double precision :: USCAL 
double precision :: H00   
double precision :: HMTS  

! Relaxation at the poles
double precision :: QX_REL(N,M),     &
               & QY_REL(N,M)
double precision :: Relaxation, QRelax_str
double precision :: Relax_M(M)
integer    :: DP_relax


double precision ::  F0_52
double precision :: G_52
double precision ::  R_52

! CHARACTERISTICS OF THE FLOW
double precision :: USCAL_52 
double precision ::  H00_52   
double precision ::  HMTS_52, GMM_52

Integer :: IPRINT


! CONTROL TESTS: ZONAL FLOW OR ROSSBY WAVE
INTEGER:: IRHW, &
        & DP_Depth


! ! CONTROL EVALUATION OF THE PRESSURE GRADIENTS: IPS=0 CENTERED DIFFERNCING
! CHOOSING IPS=1 GIVES WEIGHTED AVERAGE OF ONE-SIDED DERIVATIVES THAT
! CONVERGES TO ONE SIDED DERIVATIVE AT THE CROSSECTION WITH THE BOTTOM
! !JA IPS=1 is deactivated
INTEGER, PARAMETER  :: IPS=0

! CREATE TAPEWR
INTEGER, PARAMETER :: IANAL=0,IRST=0,IWRITE=0
INTEGER, PARAMETER :: NFIL=50
double precision :: DTFIL(NFIL),&
     & NTFIL(NFIL)

INTEGER:: NITER,  &
        & NITSM,  &
        & ICOUNT
double precision :: ERROR


double precision :: PI, &
     & PI2, &
     & PIH, &
     & PVEL, &
     & BETA, &
     & GI, &
     & EP
     
double precision :: PI_52, &
     & PI2_52, &
     & PIH_52, &
     & PVEL_52, &
     & BETA_52, &
     & GI_52

INTEGER :: IORD, &
         & ISOR, &
         & NONOS, &
         & IDIV, &
         & ISGNSPL

! grid creation
double precision :: GC1, &
     & GC2, &
     & GH1, &
     & GH2

double precision :: GC1_52, &
     & GC2_52, &
     & GH1_52, &
     & GH2_52
     
double precision :: D0, &
     & S_full         , &
     & Exit_Cond 

double precision :: D_Adv(N,M)
 
! the rest
double precision :: GMM, SUM1, SUM0
INTEGER :: I, J, KF, KT, NPLOT, NPRINT, NT, num_of_bits, ID_PREC
INTEGER :: stencil, ueber

logical :: codesignQ, codesignD, mountain, gcr52, save_time, QRelax

!!! Machine Learning
character(len=150) :: filename_NN, latitude_layer
double precision :: Layer1(32,177,1) , &
                    Layer2(32,6  ,1), &
                    max_val
INTEGER          :: Order_Layer1(32,176), max_ind, k



!! read in weights for L0N0, 5x5 preconditioner
do i=0,31
  write(latitude_layer,*) i
  filename_NN= 'weights/FortranWeights_model2x2_1x0_l1w5coeff_R_Lat'//trim(adjustl(latitude_layer))//'_layer0.txt' 
  call read_layer(filename_NN,Layer1(i+1,:,:),177,1)
enddo

Order_Layer1(:,:)=0


!! setup shallow water model
!!


mountain = .false.


  codesignQ = .FALSE.
  codesignD = .TRUE.
  gcr52     = .false.



  QRelax    =.false.     !! no polar filtering
  stencil=0              
  QRelax_str=0.0d0       !! filtering strength !4.6d0*10.0d0**(-2.0d0)
!0.00000001d0

do ID_PREC=5,5,-5  !! preconditioner, later overwritten in linear solver: "gcr_pre"
  do IRHW = 3,3   !! zonal flow with ETOPO5 topography

  
  DP_Depth=0          !! no varying precision
  do DP_relax=0,0     !! no polar filtering

   do num_of_bits=60,60,1  !! standard double precision: 52 bits
   write(Dp_depth_str,*) num_of_bits


   EXP_NAME= 'data/data_ADI_NOTgcr_EXP3_Dp_1M10_MLPrecon_iterINF_FULL_res05'
  ! EXP_NAME= 'data_ADI_Precon_init52'

   !EXP_NAME= 'data/PiotrV2_ADI_L2_1M3_R10_D'//trim(adjustl(Dp_depth_str))//'_NoQfil_dt200_res4'


    save_time=.false.
   
    If (num_of_bits>10) then
! if .not. init52


  write(*,*) 'Experiment' , IRHW
  write(*,*) 'Resolution' , NLON, NLAT
  write(*,*) 'Preconditioner', ID_PREC
  write(*,*) 'Dp Zonal bands', DP_Depth
  write(*,*) 'Experiment folder ',   EXP_NAME
  write(*,*) 'Relaxation ',  QRelax
  write(*,*) 'Relaxation strength ',  QRelax_str
  write(*,*) 'Relaxation depth ',  DP_relax
  write(*,*) 'Relax_Stencil', stencil

  sum_time=0.0d0
  sum_lp_time=0.0d0


 call rpenum_init(num_of_bits)

   write(*,*) 'codesign: Q | D ', codesignQ, codesignD, 'default bits: ', num_of_bits


if (IRHW==2) then
NT = 4*60*30 ! 345600/4/4/2 !322560 !60480 !40320 
NPRINT =60*24/4 ! 23040/4/4/2 !4320!2880
DT_52= 4.0d0* 60  !3.0d0*4.0d0*4.0d0*2.0d0 !20.0d0 

elseif (IRHW==3) then    !! zonal flow with ETOPO5
NT = 4*60*24*30/4 ! 345600/4/4/2 !322560 !60480 !40320 
NPRINT = 60/4 !hourly 60*24/4 !daily 23040/4/4/2 !4320!2880
DT_52= 4.0d0* 60  !3.0d0*4.0d0*4.0d0*2.0d0 !20.0d0 

else
NT = 4*60*30 ! 345600/4/4/2 !322560 !60480 !40320 
NPRINT = 60/4 !hourly 60*24/4 !daily 23040/4/4/2 !4320!2880
DT_52= 4.0d0* 60
endif
write(*,*) 'Timestep:', DT_52

 DT = DT_52

 
 
F0_52=1.4584E-4
F0=F0_52

 
G_52 =9.80616
G=G_52
GI_52=1.0d0/G_52
GI = GI_52

! end param test1
R_52=6371.22E+03
R=R_52


! CHARACTERISTICS OF THE FLOW
USCAL_52 = 5.0d0 
H00_52  = 8.E3
HMTS_52  = 1.E-6 


if (IRHW==2) then
!  
  mountain=.true.
  USCAL_52 = 20.0d0 
  H00_52   = 5960.0d0 
  HMTS_52  = 1.0d0 
!
elseif (IRHW==3) then

!  
  mountain=.true.
  USCAL_52 = 20.0d0 
  H00_52   = 5960.0d0 
  HMTS_52  = 1.0d0 
!
  endif

USCAL=USCAL_52
H00 = H00_52
HMTS= HMTS_52


DTFIL(:) = NFIL*40.0d0
NTFIL(:) = NFIL*2160


NITSM  = 0
ICOUNT = 0




!COMPUTE SOME RELEVANT CONSTANTS

PI_52=ACOS(-1.0d0)
PI = PI_52
PI2_52=2.0d0*ACOS(-1.0d0)
PI2 = PI2_52
PIH_52=0.5d0*ACOS(-1.0d0)
PIH = PIH_52
TIME=0.0d0
PVEL_52=USCAL_52/R_52
PVEL = PVEL_52
F0_52=ICORIO*F0_52
F0 = F0_52

! BETA IS AN ANGLE IN THE ZONAL FLOW TEST

BETA_52=0.0d0 
if (IRHW==2) then
BETA_52=0.0d0 
elseif (IRHW==3) then
BETA_52=0.0d0 
endif
BETA = BETA_52
H00_52=G_52*H00_52
H00 = H00_52


! PARAMETERS FOR MPDATA ADVECTION
IORD=2
ISOR=1
NONOS=1
IDIV=1               
ISGNSPL=0

! SMALL CONSTANT (SMALL COMPARED TO DEPTH OF THE ATMOSPHERE)
EP= 1.E-6   !! ersetzen durch tiny times factor oder was aehnliches

!COMPUTE GRID      


DX_52= 2.0d0*ACOS(-1.0d0)/FLOAT(N-2)
DX = DX_52
DY_52= ACOS(-1.0d0)/FLOAT(M)
DY = DY_52

GC1_52= DT_52/(2.0d0*ACOS(-1.0d0)/FLOAT(N-2))                                ! DT/DX
GC1 = GC1_52
GC2_52= DT_52/(ACOS(-1.0d0)/FLOAT(M))                                        ! DT/DY
GC2 = GC2_52

GH1_52= 0.5d0* DT_52/(2.0d0*ACOS(-1.0d0)/FLOAT(N-2))                         ! 0.5d0*GC1
GH1 = GH1_52
GH2_52= 0.5d0* DT_52/(ACOS(-1.0d0)/FLOAT(M))                                 ! 0.5d0*GC2
GH2 = GH2_52


! end param test2

DO J=1,M+1
  Y_52(J)=-PIH_52+(float(J)-0.5d0)*DY_52

  Y(J) = Y_52(J)

end do


DO I=2,N-1
  X_52(I)=(float(I)-1)*DX_52
  X(I) = X_52(I)
end do
X_52(1)=X_52(N-1)
X_52(N)=X_52(2)

X(1)=X(N-1)
X(N)=X(2)

DO I=1,N
  IP(I)=MOD(I+(N-2)/2-1,N-2)+1
end do
!COMPUTE METRIC TERMS FOR SPHERICAL COORDINATES

DO J=1,M
  DO I=1,N

    HX_52(I,J)=R_52*COS(Y_52(J))

    HX(I,J) = HX_52(I,J)
    HY_52(I,J)=R_52
    HY(I,J) = HY_52(I,J)

    S_52(I,J)=HX_52(I,J)*HY_52(I,J)
    S(I,J) = S_52(I,J)

  end do

end do

DO I=2,N-1
  DO J=2,M-1
    DHX2Y_52(I,J)= (HX_52(I,J+1)-HX_52(I,J-1))*GC2_52/S_52(I,J)*0.5d0
    DHX2Y(I,J) = DHX2Y_52(I,J) 
  end do

  DHX2Y_52(I,1)= (HX_52(I,  2)+HX_52(I,  1))*GC2_52/S_52(I,1)*0.5d0
  DHX2Y(I,1) = DHX2Y_52(I,1)
  DHX2Y_52(I,M)=-(HX_52(I,  M)+HX_52(I,M-1))*GC2_52/S_52(I,M)*0.5d0
  DHX2Y(I,M) = DHX2Y_52(I,M) 
end do


      
CALL XBC(DHX2Y,N,M)
CALL XBC_52(DHX2Y_52,N,M)




!CONDITIONS OF THE INITIAL STATE ***********************************


If(IRHW .NE. 3) then
CALL TOPOGR(P0_HP,X_52,Y_52,N,M, mountain)
else
CALL ETOPO5(P0_HP,X_52,Y_52,N,M, mountain)  !! load ETOPO5 dataset
endif
!call write_MLfields(P0_HP, 0, N, M, EXP_NAME, 'Topo')

IF(IRHW.EQ.0) CALL INITZON(U_52,V_52,PT_HP,COR_52,X_52,Y_52,N,M,F0_52,BETA_52,H00_52,R_52,PVEL_52)
IF(IRHW.EQ.1) CALL INITRHW(U_52,V_52,PT_HP,COR_52,X_52,Y_52,N,M,F0_52,R_52)
IF(IRHW.EQ.2) CALL INITZON(U_52,V_52,PT_HP,COR_52,X_52,Y_52,N,M,F0_52,BETA_52,H00_52,R_52,PVEL_52)
IF(IRHW.EQ.3) CALL INITZON(U_52,V_52,PT_HP,COR_52,X_52,Y_52,N,M,F0_52,BETA_52,H00_52,R_52,PVEL_52)


IF(IRST.EQ.0) THEN
! INITIATE PRIMARY VARIABLES (PD, QX, QY)
  DO J=1,M
    DO I=1,N
      COR(I,J)=COR_52(I,J)*DT_52 
      P0_HP(I,J)=  P0_HP(I,J)*HMTS_52
      PD_HP(I,J)= max(rpe_0, PT_HP(I,J)*GI_52-P0_HP(I,J))

      QX_HP(I,J)=PD_HP(I,J)*U_52(I,J)
      QY_HP(I,J)=PD_HP(I,J)*V_52(I,J)

    end do

  end do

!! call calculateAVGdepth(PD)
D0=0.0d0
S_full=0.0d0

  DO J=1,M
    DO I=2,N-1
      S_full=S_full+S_52(I,J)
    end do
  end do
  write(*,*) 'Area Earth', S_full
  DO J=1,M
    DO I=2,N-1
      D0=D0+PD_HP(I,J)*S_52(I,J)/S_full
    end do
  end do
  write(*,*) 'average Depth', D0  

If(QRelax) then

do J=1,M
 do I=1+stencil,N-stencil
  QX_REL(I,J)=SUM(QX_HP(I-stencil:I+stencil,J))/(1.0d0+2.0d0*stencil)
  QY_REL(I,J)=SUM(QY_HP(I-stencil:I+stencil,J))/(1.0d0+2.0d0*stencil)
 enddo
enddo
do J=1,M
 do I=2,stencil
  ueber=2-(I-stencil)
  QX_REL(I,J)=(SUM(QX_HP(N-ueber:N,J))+SUM(QX_HP(2+1:I+stencil,J)))/(1.0d0+2.0d0*stencil)
  QY_REL(I,J)=(SUM(QY_HP(N-ueber:N,J))+SUM(QY_HP(2+1:I+stencil,J)))/(1.0d0+2.0d0*stencil)
 enddo
enddo
do J=1,M
 do I=N-stencil+1,N-1
  ueber=(I+stencil)-(N-1)
  QX_REL(I,J)=(SUM(QX_HP(I-stencil:N-1,J))+SUM(QX_HP(2:1+ueber,J)))/(1.0d0+2.0d0*stencil)
  QY_REL(I,J)=(SUM(QY_HP(I-stencil:N-1,J))+SUM(QY_HP(2:1+ueber,J)))/(1.0d0+2.0d0*stencil)
 enddo
enddo
CALL XBC(QX_REL,N,M)
CALL XBC(QY_REL,N,M)


endif

If(QRelax) then


  Do J=1,M

    Relax_M(J)=max(0.0d0,&
                 & 1.0d0/COS(Y_52(J))-1.0d0/COS(Y_52(DP_relax+1))) &
                 & /(1.0d0/COS(Y_52(1))-1.0d0/COS(Y_52(DP_relax+1))) !scale back to 0,1 interval

  end do

Relax_M(:)=QRelax_str*Relax_M(:)

endif

!! Prepare linear Operator coefficients for subroutine gcr_pre
If(.not. QRelax) then
      DO J=1,M

       GMM=0.5d0*COR(1,J)
       A= rpe_1/(rpe_1+GMM*GMM)
       B=GMM*A
       C=GH1*HY(1,J)*0.5d0
       D=GH2*HX(1,J)*0.5d0

       MGH1IHX(J)=-0.5d0*GC1/HX(1,J)
       MGH2IHY(J)=-0.5d0*GC2/HY(1,J)
       AC(J)=A*C
       BC(J)=B*C
       AD(J)=A*D
       BD(J)=B*D
      enddo
else

      DO J=1+DP_relax,M-DP_relax
       GMM=0.5d0*COR(1,J)
       A= rpe_1/(rpe_1+GMM*GMM)
       B=GMM*A
       C=GH1*HY(1,J)*0.5d0
       D=GH2*HX(1,J)*0.5d0

       MGH1IHX(J)=-0.5d0*GC1/HX(1,J)
       MGH2IHY(J)=-0.5d0*GC2/HY(1,J)
       AC(J)=A*C
       BC(J)=B*C
       AD(J)=A*D
       BD(J)=B*D
      enddo

      DO J=1,DP_relax
        
        Relaxation=1.0d0+0.5d0*DT_52*Relax_M(J)
       GMM=0.5d0*COR(1,J)
       A= rpe_1/(Relaxation+GMM*GMM/(Relaxation))
       B=GMM*A/(Relaxation)
       C=GH1*HY(1,J)*0.5d0
       D=GH2*HX(1,J)*0.5d0

       MGH1IHX(J)=-0.5d0*GC1/HX(1,J)
       MGH2IHY(J)=-0.5d0*GC2/HY(1,J)
       AC(J)=A*C
       BC(J)=B*C
       AD(J)=A*D
       BD(J)=B*D
      enddo

      DO J=M+1-DP_relax,M

        Relaxation=1.0d0+0.5d0*DT_52*Relax_M(J)
       GMM=0.5d0*COR(1,J)
       A= rpe_1/(Relaxation+GMM*GMM/(Relaxation))
       B=GMM*A/(Relaxation)
       C=GH1*HY(1,J)*0.5d0
       D=GH2*HX(1,J)*0.5d0

       MGH1IHX(J)=-0.5d0*GC1/HX(1,J)
       MGH2IHY(J)=-0.5d0*GC2/HY(1,J)
       AC(J)=A*C
       BC(J)=B*C
       AD(J)=A*D
       BD(J)=B*D
      enddo

endif


  DO J=1,M
    DO I=1,N
      P0(I,J) = P0_HP(I,J)
      PT(I,J) = PT_HP(I,J) 
      PD(I,J) = PD_HP(I,J)
      QX(I,J) = QX_HP(I,J) 
      QY(I,J) = QY_HP(I,J) 
    end do
  end do

  DO J=1,M
    DO I=1,N
      U(I,J,0) =U_52(I,J)
      V(I,J,0) =V_52(I,J)
      U(I,J,1)=U(I,J,0)
      V(I,J,1)=V(I,J,0)
    end do
  end do
  DO J=1,M
    DO I=1,N
      PD_T(I,J)=rpe_0
      QX_T(I,J)=rpe_0
      QY_T(I,J)=rpe_0
    end do
  end do

!! prepare output
   call init_perf_markers(PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:), rpe_0, &
                   & codesignQ, codesignD, IRHW, X, Y, N, M, num_of_bits, ID_prec, EXP_NAME)
!! write initial fields
   call write_fields(PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:), rpe_0, &
& codesignQ, codesignD, IRHW, X, Y, N, M, num_of_bits, ID_prec, EXP_NAME)  !PD_HP(:,:)

  
  
! INITIATE CORRESPONDING FORCES R_x R_y
  CALL PRFC0(PT,F1(:,:,0),F2(:,:,0),PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)

  DO J=1,M
    DO I=1,N 

      E1(I,J,0)=-DHX2Y(I,J)*QX(I,J)*QY(I,J)/max(PD(I,J),EP)
      E2(I,J,0)= DHX2Y(I,J)*QX(I,J)*QX(I,J)/max(PD(I,J),EP)
      E1(I,J,-1)=E1(I,J,0)
      E2(I,J,-1)=E2(I,J,0)
    end do
  end do


  if(.not. QRelax) then
    DO J=1,M
      DO I=1,N

        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0)
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0)

      end do
    end do
  else

      DO J=1+DP_relax,M-DP_relax
        DO I=1,N

        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0)
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0)
        enddo
      enddo

      DO J=1,DP_relax
        
        DO I=1,N
        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0) -DT_52*Relax_M(J)*(QX(I,J)-QX_REL(I,J))
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0) -DT_52*Relax_M(J)*(QY(I,J)-QY_REL(I,J))
        enddo
      enddo

      DO J=M+1-DP_relax,M
       
        DO I=1,N
        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0) -DT_52*Relax_M(J)*(QX(I,J)-QX_REL(I,J))
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0) -DT_52*Relax_M(J)*(QY(I,J)-QY_REL(I,J))
        enddo
      enddo
  endif


ELSE
  DO KF=1,NFIL
    READ(10) PD,PT,QX,QY,U,V,F1,F2,E1,E2
    IF(IWRITE.EQ.1) WRITE(9) PD,PT,QX,QY,U,V,F1,F2,E1,E2
    TIME=TIME+DTFIL(KF)*NTFIL(KF)
  end do
ENDIF

!! model diagnostics
CALL DIAGNOS(QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:),PD_HP,&
                 & PT,HX,HY,IP,S,TIME,DX,DY,DT, SUM0,SUM1, &
                 & KT,N,M,0, NITER,NITSM,ICOUNT,ERROR, sum_time, sum_lp_time)

! CLOSE INITIAL CONDITIONS *************************************

! COMPUTE SOLUTION IN TIME *************************************

IF(IANAL.EQ.0) THEN
  DO KT=1,NT
    !write(*,*) KT
    IPRINT=0
    IF(KT/NPRINT*NPRINT.EQ.KT) IPRINT=1

    ! COMPUTE ADVECTIVE COURANT NUMBERS
    ! COMPUTE VELOCITY PREDICTOR
    CALL VELPRD(U,V,F1,F2,PD,HX,HY,IP,N,M,GC1,GC2,EP)
      
  !IF(IPRINT.EQ.1) then
  ! call divergence (U,V) for exit condition
    Call DIVER(D_Adv(:,:),U(:,:,1),V(:,:,1),HX,HY,S,N,M,IP,1)

    Exit_Cond=maxval(ABS(D_Adv(:,:)))
    Exit_Cond=2.0d0*Exit_Cond*D0*G_52*DT_52

  !endif 

    DO J=1,M
      DO I=1,N
        U(I,J,1)=U(I,J,1)*HY(I,J)
        V(I,J,1)=V(I,J,1)*HX(I,J)
      end do
    end do
    ! COMPUTE COURANT NUMBERS AT STAGGERED TIME/SPACE POSITIONS
    DO J=1,M
      DO I=2,N-1
        UA(I,J)=(U(I,J,1)+U(I-1,J,1))*GH1
      end do
    end do

    CALL XBC(UA,N,M)
    
    DO I=1,N
      DO J=2,M
        VA(I,J)=(V(I,J,1)+V(I,J-1,1))*GH2
      end do
      VA(I,  1)=rpe_0
      VA(I,M+1)=rpe_0
    end do

!   call compint(UA,VA,N,M)

! CLOSE ADVECTIVE COURANT NUMBERS

! COLLECT EXPLICIT PARTS OF CONTINUITY AND MOMENTUM EQUATIONS:

! COLLECT FROM CONTINUITY EQUATION
    DO J=1,M
      DO I=1,N
        U(I,J,1)=QX(I,J)*GC1
        V(I,J,1)=QY(I,J)*GC2
      end do
    end do

    CALL DIVER(F1(:,:,1),U(:,:,1),V(:,:,1),HX,HY,S,N,M,IP,1)

    DO J=1,M
      DO I=1,N
        F1(I,J,1)=PD(I,J)+P0(I,J)-rpe_05*F1(I,J,1)
      end do
    end do

! C--->                       ADVECTION

    DO J=1,M
      DO I=1,N
        QX(I,J)=QX(I,J)+rpe_05*F1(I,J,0)
        QY(I,J)=QY(I,J)+rpe_05*F2(I,J,0)
!       PC(I,1)=PD(I,1)
      end do
    end do

    
    IF (codesignQ) then
      DO J=1,M
        DO I=1,N
          QX_T(I,J)=QX_T(I,J)+rpe_05*F1(I,J,0)
          QY_T(I,J)=QY_T(I,J)+rpe_05*F2(I,J,0)

        end do
      end do
    end if


       !! MPDATA for explicit computation of Momenta
      CALL MPDATT(UA,VA,QX,S,N,M,IORD,ISOR,NONOS,IDIV,-1, IP, QX_T, codesignQ)
      CALL MPDATT(UA,VA,QY,S,N,M,IORD,ISOR,NONOS,IDIV,-1, IP, QY_T, codesignQ)


!--->                CORIOLIS AND METRIC FORCES
    DO J=1,M
      DO I=1,N
        UA(I,J)=QX(I,J)+rpe_05*(rpe_2*E1(I,J,0)-E1(I,J,-1))
        VA(I,J)=QY(I,J)+rpe_05*(rpe_2*E2(I,J,0)-E2(I,J,-1))
      end do   
   end do
   
    IF (codesignQ) then
      DO J=1,M
        DO I=1,N
          QX_T(I,J)=QX_T(I,J)+rpe_05*(rpe_2*E1(I,J,0)-E1(I,J,-1))
          QY_T(I,J)=QY_T(I,J)+rpe_05*(rpe_2*E2(I,J,0)-E2(I,J,-1))

        end do
      end do
    end if  
    
   
     IF (.not. codesignQ) then

  if(.not. QRelax) then
      DO J=1,M
        DO I=1,N
          GMM=rpe_05*COR(I,J)
          QX(I,J)=(UA(I,J)+GMM*VA(I,J))/(rpe_1+GMM**2)
          QY(I,J)=(VA(I,J)-GMM*UA(I,J))/(rpe_1+GMM**2)

        end do
      end do
  else
     !! bring in possible relaxation
      DO J=1,DP_relax
       
        DO I=1,N
        UA(I,J)=UA(I,J)+0.5d0*DT_52*Relax_M(J)*QX_REL(I,J)
        VA(I,J)=VA(I,J)+0.5d0*DT_52*Relax_M(J)*QY_REL(I,J)
      end do   
   end do

      DO J=M+1-DP_relax,M
      
        DO I=1,N
        UA(I,J)=UA(I,J)+0.5d0*DT_52*Relax_M(J)*QX_REL(I,J)
        VA(I,J)=VA(I,J)+0.5d0*DT_52*Relax_M(J)*QY_REL(I,J)
      end do   
   end do

      DO J=1+DP_relax,M-DP_relax
        DO I=1,N

          GMM=rpe_05*COR(I,J)
          QX(I,J)=(UA(I,J)+GMM*VA(I,J))/(rpe_1+GMM**2)
          QY(I,J)=(VA(I,J)-GMM*UA(I,J))/(rpe_1+GMM**2)

        enddo
      enddo

      DO J=1,DP_relax
      
        Relaxation=1.0d0+0.5d0*DT_52*Relax_M(J)
        DO I=1,N
          GMM=rpe_05*COR(I,J)
          
          QX(I,J)=(UA(I,J)+GMM/Relaxation*VA(I,J))/(Relaxation+GMM**2/Relaxation)
          QY(I,J)=(VA(I,J)-GMM/Relaxation*UA(I,J))/(Relaxation+GMM**2/Relaxation)

        enddo
      enddo

      DO J=M+1-DP_relax,M
      
        Relaxation=1.0d0+0.5d0*DT_52*Relax_M(J)
        DO I=1,N
          GMM=rpe_05*COR(I,J)
          QX(I,J)=(UA(I,J)+GMM/Relaxation*VA(I,J))/(Relaxation+GMM**2/Relaxation)
          QY(I,J)=(VA(I,J)-GMM/Relaxation*UA(I,J))/(Relaxation+GMM**2/Relaxation)
        enddo
      enddo
    endif


      DO J=1,M
        DO I=1,N

            QX_HP(I,J)= QX(I,J)
            QY_HP(I,J)= QY(I,J)
        end do
      end do
    
      DO J=1,M
        DO I=1,N

          QX(I,J)=QX(I,J)*GC1
          QY(I,J)=QY(I,J)*GC2

        end do

      end do
   
   
    elseif (codesignQ) then 
    
       DO J=1,M
        DO I=1,N

            UA_HP(I,J)= QX_HP(I,J)+ QX_T(I,J)
            VA_HP(I,J)= QY_HP(I,J)+ QY_T(I,J)

        end do
      end do
      
     DO J=1,M
       DO I=1,N

         GMM_52=0.5d0*COR_52(I,J)*DT_52 
         QX_HP(I,J)=(UA_HP(I,J)+GMM_52*VA_HP(I,J))/(rpe_1+GMM_52**2)
         QY_HP(I,J)=(VA_HP(I,J)-GMM_52*UA_HP(I,J))/(rpe_1+GMM_52**2)

       end do
     
     end do
      
  
      DO J=1,M
        DO I=1,N

            QX(I,J)= QX_HP(I,J)
            QY(I,J)= QY_HP(I,J)
            
            QX(I,J)= QX(I,J)*GC1
            QY(I,J)= QY(I,J)*GC2

            QX_T(I,J)= rpe_0
            QY_T(I,J)= rpe_0
        end do
      end do
    

    end if  
    
     DO J=1,M
       DO I=1,N
         E1(I,J,-1)=E1(I,J,0)
         E2(I,J,-1)=E2(I,J,0)
       end do 
     end do
    


  
    
!--->            DIVERGENCE OF ALL COLLECTED TERMS:  !! divergence of incoming mass from explicit part of new momenta QX, QY
    CALL DIVER(F2(:,:,1),QX,QY,HX,HY,S,N,M,IP,1)


    DO J=1,M
      DO I=1,N

        F1(I,J,1)=(F1(I,J,1)-rpe_05*F2(I,J,1))*G   !! (h+h0)*g pressure minus incoming mass of 0.5 old and 0.5 new momentum
        F2(I,J,1)=PT(I,J)                          !! old one without new contributions
        PD(I,J)=P0(I,J)*G
      end do
    end do
  
    CALL PRFORC(PT,E1(:,:,0),E2(:,:,0),PD,F2(:,:,1), &
            &   F1(:,:,0),F2(:,:,0),HX,HY,COR,N,M,IP,GC1,GC2,0,0)


! COMPUTE ELLIPTIC EQUATION WITH Conjugate Residual SCHEME

  !call write_MLfields(PT, KT, N, M, EXP_NAME, 'H_in')
  !call write_MLfields(F1(:,:,1), KT, N, M, EXP_NAME, 'RHS')


  CALL  GCR_PRE(PT,F1(:,:,0),F2(:,:,0),HX,HY,S,F1(:,:,1),F2(:,:,1), &
       &      PD(:,:),E1(:,:,0),E2(:,:,0),COR,IP, &
       &      U(:,:,0),U(:,:,1),V(:,:,0),V(:,:,1),N,M,GC1,GC2,   &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
  &      niter,nitsm,icount,error, PD_T, sum_time, sum_lp_time, ID_PREC,.FALSE., save_time,&
       &      TIME, codesignQ, IRHW, X, Y, Exit_Cond,               &
       &  EXP_NAME, iprint, num_of_bits, DP_Depth, KT, Layer1, Layer2)

 If(save_time) exit  ! if iteration counter reached 100 once, exit simulation
  
   

    
   IF ( codesignD) then
         
      DO J=1,M
        DO I=1,N
            PT_HP(I,J)= PT_HP(I,J)+ PD_T(I,J)

        end do
      end do 

      DO J=1,M
        DO I=1,N
            PT(I,J)= PT_HP(I,J)
        end do
      end do
    

     end if
   
  If(.not. QRelax) then     
    CALL PRFORC(PT,F1(:,:,0),F2(:,:,0),PD,F2(:,:,1), &
      &      E1(:,:,0),E2(:,:,0),HX,HY,COR,N,M,IP,GC1,GC2,1,1)
  else
    CALL PRFORC_FIN(PT,F1(:,:,0),F2(:,:,0),PD,F2(:,:,1), &
       &      E1(:,:,0),E2(:,:,0),HX,HY,COR,N,M,IP,GC1,GC2,Relax_M(:),DT_52,1,1)
  endif

       
   !!! update PD, because PRFORC needs an intermediate value of PD
    IF (.not. codesignD) then   
      
      DO J=1,M
        DO I=1,N
          PD(I,J)= max(EP, PT(I,J)*GI-P0(I,J))
   
        end do
      end do
   
        
      DO J=1,M
        DO I=1,N
            PD_HP(I,J)= PD(I,J)
   
        end do
      end do
    
    else
    
      
      DO J=1,M
        DO I=1,N
            PD_HP(I,J)= max(EP, PT_HP(I,J)*GI_52-P0_HP(I,J))

        end do
      end do
      

      
      DO J=1,M
        DO I=1,N
            PD(I,J)= PD_HP(I,J)
            PD_T(I,J)= rpe_0
        end do
      end do
    

     end if
      ! call write_MLfields(PT, KT, N, M, EXP_NAME, 'H_out')    
     
! COMPUTE SOLUTION'S UPDATE

    DO J=1,M
      DO I=1,N
        QX(I,J)=(QX(I,J)+F1(I,J,0)*GI)/GC1
        QY(I,J)=(QY(I,J)+F2(I,J,0)*GI)/GC2

      end do
    end do

    IF (codesignQ) then
      DO J=1,M
        DO I=1,N
          QX_T(I,J)=QX_T(I,J)+F1(I,J,0)*GI/GC1
          QY_T(I,J)=QY_T(I,J)+F2(I,J,0)*GI/GC2

        end do
      end do
    end if
   
 
If(QRelax) then
!! update relaxation

do J=1,M
 do I=1+stencil,N-stencil
  QX_REL(I,J)=SUM(QX(I-stencil:I+stencil,J))/(1.0d0+2.0d0*stencil)
  QY_REL(I,J)=SUM(QY(I-stencil:I+stencil,J))/(1.0d0+2.0d0*stencil)
 enddo
enddo
do J=1,M
 do I=2,stencil

  ueber=2-(I-stencil)
  QX_REL(I,J)=(SUM(QX(N-ueber:N,J))+SUM(QX(2+1:I+stencil,J)))/(1.0d0+2.0d0*stencil)
  QY_REL(I,J)=(SUM(QY(N-ueber:N,J))+SUM(QY(2+1:I+stencil,J)))/(1.0d0+2.0d0*stencil)
 enddo
enddo
do J=1,M
 do I=N-stencil+1,N-1

  ueber=(I+stencil)-(N-1)
  QX_REL(I,J)=(SUM(QX(I-stencil:N-1,J))+SUM(QX(2:1+ueber,J)))/(1.0d0+2.0d0*stencil)
  QY_REL(I,J)=(SUM(QY(I-stencil:N-1,J))+SUM(QY(2:1+ueber,J)))/(1.0d0+2.0d0*stencil)
 enddo
enddo
CALL XBC(QX_REL,N,M)
CALL XBC(QY_REL,N,M)

endif
    
! COMPUTE NEW FORCES



    CALL PRFC0(PT,F1(:,:,0),F2(:,:,0),PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)
    
    DO J=1,M
      DO I=1,N
        E1(I,J,0)=-DHX2Y(I,J)*QX(I,J)*QY(I,J)/PD(I,J)
        E2(I,J,0)= DHX2Y(I,J)*QX(I,J)*QX(I,J)/PD(I,J)

      end do
    end do

  if(.not. QRelax) then
    DO J=1,M
      DO I=1,N

        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0)
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0)

      end do
    end do
  else

      DO J=1+DP_relax,M-DP_relax
        DO I=1,N

        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0)
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0)

        enddo
      enddo

      DO J=1,DP_relax
      
        DO I=1,N
        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0) -DT_52*Relax_M(J)*(QX(I,J)-QX_REL(I,J))
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0) -DT_52*Relax_M(J)*(QY(I,J)-QY_REL(I,J))
        enddo
      enddo

      DO J=M+1-DP_relax,M
    
        DO I=1,N
        F1(I,J,0)=F1(I,J,0)+COR(I,J)*QY(I,J)+E1(I,J,0) -DT_52*Relax_M(J)*(QX(I,J)-QX_REL(I,J))
        F2(I,J,0)=F2(I,J,0)-COR(I,J)*QX(I,J)+E2(I,J,0) -DT_52*Relax_M(J)*(QY(I,J)-QY_REL(I,J))
        enddo
      enddo



  endif

    DO J=1,M
      DO I=1,N

        U(I,J,0)=QX(I,J)/PD(I,J)
        V(I,J,0)=QY(I,J)/PD(I,J)
      end do
    end do

!COMPUTE OUTPUTED FIELDS ****************************************

    IF(.not. (KT/NPRINT*NPRINT.NE.KT)) then
    
      
      IF (.not. codesignQ) then
      

    
        DO J=1,M
          DO I=1,N

            QX_HP(I,J)= QX(I,J)
            QY_HP(I,J)= QY(I,J)
          end do
        end do
    
    else
    
      DO J=1,M
        DO I=1,N

            QX_HP(I,J)= QX_HP(I,J)+ QX_T(I,J)
            QY_HP(I,J)= QY_HP(I,J)+ QY_T(I,J)

        end do
      end do
      
   
      DO J=1,M
        DO I=1,N

            QX_T(I,J)= rpe_0
            QY_T(I,J)= rpe_0
        end do
      end do
    

    end if  
    
    
      IF(IWRITE.EQ.1) WRITE(9) PD,PT,QX,QY,U,V,F1,F2,E1,E2
      CALL DIAGNOS(QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:),PD_HP,&
                 & PT,HX,HY,IP,S,TIME,DX,DY,DT, SUM0,SUM1, &
                 & KT,N,M,1, NITER,NITSM,ICOUNT,ERROR, sum_time, sum_lp_time)
      call write_perf_markers (PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:)&
                  &, TIME, codesignQ, codesignD, IRHW, X, Y, N, M, &
                  & num_of_bits,NITER,NITSM,ICOUNT, sum_time, sum_lp_time)
      call write_fields(PD_HP(:,:)+P0(:,:),QX_HP(:,:)/PD_HP(:,:),QY_HP(:,:)/PD_HP(:,:)&
                  &, TIME, codesignQ, codesignD, IRHW, X, Y, N, M, num_of_bits, ID_prec, EXP_NAME)
    end if

  end do
!!CLOSE TIME INTEGRATION

ENDIF

IF(IANAL.EQ.1) THEN
  
  TIME=0.0
  NPLOT=5
 
  DO KF=1,NFIL
    TIME=TIME+DTFIL(KF)*NTFIL(KF)/3600.0d0
    PRINT 300, TIME      
    300 FORMAT(14X,5HTIME=,F7.2,6H HOURS)
    READ(10) PD,PT,QX,QY,U,V,E1,E2,F1,F2
    IF(KF/NPLOT*NPLOT.EQ.KF) then
      ! plot not yet finished
      !CALL PLOT(PT, P0, U(:,:,0), V(:,:,0),HX,HY,N,M,IRHW)
    end if
  end do
ENDIF

endif !!
end do !! end loop over all number of bits

call close_perf_markers
  end do
 end do
enddo

END program


subroutine topogr(h0,x,y,n,m, mountain)
use implicit_functions_DP

implicit none

INTEGER :: n, m
double precision :: h0(n,m)
double precision :: x(n),y(m)
double precision :: hs0, Rad, x_c, y_c, dist, pi, sigma
integer :: i, j
logical :: mountain

pi = acos(-1.0d0)

hs0=1800.0d0 ! 2000.0d0
Rad=pi/9.0d0
x_c=1.0d0*pi/2.0d0
y_c=pi/8.0d0!pi/6.0d0


sigma=Rad/2.150d0

! mountain shape functions
!     dist(xln,ylt) =sqrt( (cos(ylt)*sin((xln-x0)/2))**2
!    .                      +       (sin((ylt-y0)/2))**2 )/xl
!     profm(rad)=.5*(1.+cos(pi*rad))
! special for the conical mountain on the pole ccccc

!dist(xln,ylt)=abs(ylt-y0)
!profm(rad)=amax1(0., 1.-gamm*rad/rnot)
!rnot = 2*acos(-1.)/128. * 10.
!gamm=1.
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

!pi = acos(-1.)
!xl=dx * 5.
!x0 = pi
!y0 = pi*.5

do j=1,m
  do i=1,n
    h0(i,j)=0.0d0
  end do
end do

If (mountain) then
  
  do j=1,m
    do i=1,n

      dist=min(Rad, sqrt( (x(i)-x_c)**2 + (y(j)-y_c)**2) )
      h0(i,j)=hs0*(rpe_1-dist/Rad)

    end do

  end do

end if

end subroutine

!!!!!!!!! INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITZON(U,V,PT,COR,X,Y,N,M,F0,BETA,H00,R,Q)
use implicit_functions_DP

implicit none

double precision :: PT(N,M)
double precision :: U(N,M),V(N,M),COR(N,M),X(N),Y(M), F0, beta, H00, R, Q, pert
INTEGER :: N, M, seed(N*M*2)


INTEGER :: I, J
 
call random_seed(put=seed)
DO J=1,M
  DO I=1,N
    COR(I,J)=F0*(-COS(X(I))*COS(Y(J))*SIN(BETA)+SIN(Y(J))*COS(BETA))
    PT(I,J)=H00-R**2.0d0*(F0+Q)*0.5d0*Q*  &
      &     (-COS(X(I))*COS(Y(J))*SIN(BETA)+SIN(Y(J))*COS(BETA))**2
    !CALL RANDOM_NUMBER(pert)
    !PT(I,J)=PT(I,J)+PT(I,J)*(pert-0.5d0)*0.01
  end do
end do

DO J=1,M
  DO I=1,N
    U(I,J)=Q*(COS(BETA)+TAN(Y(J))*COS(X(I))*SIN(BETA))*R*COS(Y(J))
    CALL RANDOM_NUMBER(pert)
    U(I,J)=U(I,J)+U(I,J)*(pert-0.5d0)*0.05

    V(I,J)=-Q*SIN(X(I))*SIN(BETA)*R
    CALL RANDOM_NUMBER(pert)
    V(I,J)=V(I,J)+V(I,J)*(pert-0.5d0)*0.05
  end do
end do
      write (6,*)  'initzon called'

END SUBROUTINE

SUBROUTINE INITRHW(U,V,PT,COR,X,Y,N,M,F0,A)
use implicit_functions_DP

implicit none
double precision :: PT(N,M)

double precision ::U(N,M),V(N,M),F0, A,COR(N,M),X(N),Y(M)
INTEGER :: N, M


double precision ::     ATH(M), BTH(M), CTH(M), TH
double precision ::     OM,K,PH0
INTEGER :: R, I, J

OM=7.848E-6
K=7.848E-6
R=4
PH0= 78.4E3

DO J=1,M
  TH=Y(J)
  ATH(J)=OM*0.5d0*(F0+OM)*(COS(TH))**2                          &
     &  +0.25d0*K**2*(COS(TH))**(2*R)*( (R+1)*(COS(TH))**2      &
     &   +FLOAT(2*R**2-R-2)-2.0d0*R**2/(COS(TH))**2 )
  BTH(J)=(F0+2.*OM)*K/FLOAT((R+1)*(R+2))*(COS(TH))**R         &
     &       *( FLOAT(R**2+2*R+2)-((R+1)*COS(TH))**2 )
  CTH(J)=0.25d0*K**2*(COS(TH))**(2*R)*( FLOAT(R+1)*(COS(TH))**2 &
     &       -FLOAT(R+2) )  
end do

DO J=1,M
  DO I=1,N
    COR(I,J)=F0*SIN(Y(J))
    U(I,J)=A*OM*COS(Y(J))+A*K*COS(R*X(I))                      &
       &   *(COS(Y(J)))**(R-1)*(R*(SIN(Y(J)))**2-(COS(Y(J)))**2)
    V(I,J)=-A*K*R*(COS(Y(J)))**(R-1)*SIN(Y(J))*SIN(R*X(I))
    PT(I,J)=PH0+A**2*ATH(J)+A**2*BTH(J)*COS(R*X(I))      &
       &   +A**2*CTH(J)*COS(2.0d0*R*X(I))
  end do
end do
  write (6,*)  'initrhw called'

END SUBROUTINE


SUBROUTINE PRFC0(A,F1,F2,PD,HX,HY,IP,IPS,GH1,GH2,EP,N,M)
use implicit_functions_DP

implicit none
double precision :: A(N,M),F1(N,M),F2(N,M),PD(N,M),HX(N,M),HY(N,M)
INTEGER :: IP(N)
double precision :: GH1, GH2, EP
INTEGER :: IPS, N, M


double precision :: GP, GN
INTEGER :: I, J

DO I=2,N-1

  DO J=1,M
    F1(I,J)=-GH1*PD(I,J)/HX(I,J)*                          &
       & (A(I+1,J)-A(I-1,J))
  end do    

  DO J=2,M-1
    F2(I,J)=-GH2*PD(I,J)/HY(I,J)*                          &
     &  (A(I,J+1)-A(I,J-1)) 
  end do  

  F2(I,1)=-GH2*PD(I,1)/HY(I,1)*                            &
     &  (A(I,2)-A(IP(I),1)) 

  F2(I,M)=-GH2*PD(I,M)/HY(I,M)*                            &
     &  (A(IP(I),M)-A(I,M-1)) 
end do
 
CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

END SUBROUTINE


SUBROUTINE PRFORC( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,            &  !! pressure
     &            N,M,IP,GC1,GC2,NOR,IRS)
use implicit_functions_DP

implicit none
double precision :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
double precision :: GC1, GC2
INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS


INTEGER :: I, J, NM

double precision :: GH1, GH2, UTILD, VTILD, GMM

GH1=rpe_05*GC1
GH2=rpe_05*GC2

DO J=2,M-1
  DO I=2,N-1
    F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
    F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
  end do
end do
   
DO I=2,N-1
  F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
  F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
  F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
  F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

IF(NOR.EQ.1) THEN
  DO J=1,M
    DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
    end do
  end do
CALL XBC(F1,N,M)
CALL XBC(F2,N,M)
ENDIF

END SUBROUTINE


SUBROUTINE PRFORC_FIN( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,       &   !! pressure force with relaxation
     &            N,M,IP,GC1,GC2,Relax_M, DT_52, NOR,IRS)
use implicit_functions_DP

implicit none
double precision :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
double precision :: GC1, GC2
double precision :: Relax_M(M), DT_52, Relaxation
INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS


INTEGER :: I, J, NM

double precision :: GH1, GH2, UTILD, VTILD, GMM

GH1=rpe_05*GC1
GH2=rpe_05*GC2

DO J=2,M-1
  DO I=2,N-1
    F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
    F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
  end do
end do
   
DO I=2,N-1
  F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
  F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
  F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
  F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

IF(NOR.EQ.1) THEN

  DO J=1,M
    DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      Relaxation=1.0d0+0.5d0*DT_52*Relax_M(J)
      F1(I,J)=(UTILD+GMM/Relaxation*VTILD)/(Relaxation+GMM**2/Relaxation)*GH1
      F2(I,J)=(VTILD-GMM/Relaxation*UTILD)/(Relaxation+GMM**2/Relaxation)*GH2
    end do
  end do
CALL XBC(F1,N,M)
CALL XBC(F2,N,M)
ENDIF

END SUBROUTINE



SUBROUTINE DIVER(F, U,  V, HX,HY,S, N, M,IP,IFLG)  !! divergence operator

use implicit_functions_DP

 implicit none
 
double precision :: F(N,M),U(N,M),V(N,M),HX(N,M),HY(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M, IFLG


INTEGER:: I, J

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
        &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
  end do
end do    

DO I=2,N-1
  F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
      &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
  F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
      &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
end do

DO J=1,M
  DO I=2,N-1
    F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
  end do
end do

DO J=1,M
  F(1,J)=F(N-1,J)
  F(N,J)=F(2  ,J)
end do

END SUBROUTINE


SUBROUTINE PRFORC_depth( P,F1,F2,PB,P0,E1,E2,HX,HY,COR,            &  !! pressure
     &            N,M,IP,GC1,GC2,NOR,IRS, num_of_bits, DP_Depth)      !! same as PRFORC
use implicit_functions_DP

implicit none
double precision :: P(N,M),F1(N,M),F2(N,M),PB(N,M),P0(N,M),E1(N,M),E2(N,M), &
     & HX(N,M),HY(N,M),COR(N,M)
double precision :: GC1, GC2
INTEGER  :: IP(N)
INTEGER  :: N, M, NOR, IRS, DP_Depth, num_of_bits


INTEGER :: I, J, NM

double precision :: GH1, GH2, UTILD, VTILD, GMM

GH1=rpe_05*GC1
GH2=rpe_05*GC2

If (DP_Depth<=0) then
 DO J=2,M-1
   DO I=2,N-1
     F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
     F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
   end do
 end do
 
 DO I=2,N-1
   F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
   F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
   F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
   F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
 end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
         F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
         F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
        enddo
      enddo

      DO J=2,DP_Depth
        DO I=2,N-1
         F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
         F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
         F1(I,J)=-GH1*(P(I+1,J)-P(I-1,J))/HX(I,J)
         F2(I,J)=-GH2*(P(I,J+1)-P(I,J-1))/HY(I,J)
        enddo
      enddo

      DO I=2,N-1
       F1(I,1)=-GH1*(P(I+1,1)-P(I-1,1))/HX(I,1)
       F1(I,M)=-GH1*(P(I+1,M)-P(I-1,M))/HX(I,M)
       F2(I,1)=-GH2*(P(I,2)-P(IP(I),1))/HY(I,1)
       F2(I,M)=-GH2*(P(IP(I),M)-P(I,M-1))/HY(I,M)
      end do

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)

endif



IF(NOR.EQ.1) THEN

  
      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
      UTILD=F1(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E1(I,J)
      VTILD=F2(I,J)*(P0(I,J)-PB(I,J))+(P(I,J)-IRS*P0(I,J))*E2(I,J)
      GMM=rpe_05*COR(I,J)
      F1(I,J)=(UTILD+GMM*VTILD)/(rpe_1+GMM**2)*GH1
      F2(I,J)=(VTILD-GMM*UTILD)/(rpe_1+GMM**2)*GH2
        enddo
      enddo

CALL XBC(F1,N,M)
CALL XBC(F2,N,M)



ENDIF

END SUBROUTINE


SUBROUTINE DIVER_depth(F, U,  V, HX,HY,S, N, M,IP,IFLG, num_of_bits, DP_Depth)  !! divergence operator
                                                                                !! same as DIVER
use implicit_functions_DP

 implicit none
 
double precision :: F(N,M),U(N,M),V(N,M),HX(N,M),HY(N,M),S(N,M)

INTEGER :: IP(N)
INTEGER :: N, M, IFLG, DP_Depth, num_of_bits


INTEGER:: I, J


If (DP_Depth<=0) then
DO J=2,M-1
  DO I=2,N-1
    F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
        &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
  end do
end do    


DO I=2,N-1
  F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
      &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
  F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
      &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
end do

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
          F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
          &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
        enddo
      enddo


      DO J=2,DP_Depth
        DO I=2,N-1
          F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
          &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
          F(I,J)= HY(I+1,J)*U(I+1,J)-HY(I-1,J)*U(I-1,J)        & 
          &  +HX(I,J+1)*V(I,J+1)-HX(I,J-1)*V(I,J-1)
        enddo
      enddo

      DO I=2,N-1
        F(I,1)= HY(I+1,1)*U(I+1,1)-HY(I-1,1)*U(I-1,1) &
         &  +(HX(I,2)*V(I,2)+HX(I,1)*V(I,1))
        F(I,M)= HY(I+1,M)*U(I+1,M)-HY(I-1,M)*U(I-1,M)  &
         &  -(HX(I,M)*V(I,M)+HX(I,M-1)*V(I,M-1))
      end do

endif


      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
          F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-1
          F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
          F(I,J)=rpe_05*IFLG*F(I,J)/S(I,J)
        enddo
      enddo
DO J=1,M
  F(1,J)=F(N-1,J)
  F(N,J)=F(2  ,J)
end do



END SUBROUTINE



SUBROUTINE LAP0(A11,A12,A21,A22,B11,B22,    &
     &          PB,P0,E1,E2,HX,HY,COR,N,M,GC1,GC2)
use implicit_functions_DP

implicit none
double precision :: A11(N,M),A12(N,M),A21(N,M),A22(N,M),B11(N,M),B22(N,M),  &
     &      PB(N,M),P0(N,M),E1(N,M),E2(N,M),HX(N,M),HY(N,M),COR(N,M)
double precision :: GC1, GC2
INTEGER :: N, M


double precision :: GH1, GH2, C1, C2, GMM, A, B, C, D
INTEGER :: I, J


GH1=rpe_05*GC1
GH2=rpe_05*GC2

DO J=1,M
  DO I=1,N
    C1=-GH1/HX(I,J)*(P0(I,J)-PB(I,J))
    C2=-GH2/HY(I,J)*(P0(I,J)-PB(I,J))
    GMM=rpe_05*COR(I,J)
    A= rpe_1/(rpe_1+GMM*GMM)
    B=GMM*A
    C=GH1*HY(I,J)*rpe_05
    D=GH2*HX(I,J)*rpe_05
    A11(I,J)=-C1*A*C
    A12(I,J)=-C2*B*C
    A21(I,J)=-C2*A*D
    A22(I,J)= C1*B*D
    B11(I,J)=-   A*E1(I,J)*C                  &
       &     -   B*E2(I,J)*C
    B22(I,J)=-   A*E2(I,J)*D                  &
       &     +   B*E1(I,J)*D
  end do
ENDDO
end subroutine

SUBROUTINE LAP0_depth(A11,A12,A21,A22,B11,B22,  &              !! coefficients of the linear operator
           &    PB,P0,E1,E2,HX,HY,COR,N,M,GC1,GC2, &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
           & DP_Depth)
use implicit_functions_DP

implicit none
double precision :: A11(N,M),A12(N,M),A21(N,M),A22(N,M),B11(N,M),B22(N,M),  &
     &      PB(N,M),P0(N,M),E1(N,M),E2(N,M),HX(N,M),HY(N,M),COR(N,M)

double precision :: MGH1IHX(M), MGH2IHY(M), AC(M), BC(M), AD(M), BD(M)
double precision :: GC1, GC2
INTEGER :: N, M, DP_Depth


double precision :: GH1, GH2, C1, C2, GMM, A, B, C, D
INTEGER :: I, J

DO J=1+DP_Depth,M-DP_Depth
  AC(J)=AC(J)
  BD(J)=BD(J)
  BC(J)=BC(J)
  AD(J)=AD(J)
  MGH1IHX(J)=MGH1IHX(J)
  MGH2IHY(J)=MGH2IHY(J)
end do



      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N

    C1=MGH1IHX(J)*(P0(I,J)-PB(I,J))
    A11(I,J)=-C1*AC(J)
    A22(I,J)= C1*BD(J)

    C2=MGH2IHY(J)*(P0(I,J)-PB(I,J))
    A12(I,J)=-C2*BC(J)
    A21(I,J)=-C2*AD(J)
        enddo
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
    B11(I,J)=-   E1(I,J)*AC(J)                  &
       &     -   E2(I,J)*BC(J)
    B22(I,J)=-   E2(I,J)*AD(J)                  &
       &     +   E1(I,J)*BD(J)
        enddo
      enddo



      DO J=1,DP_Depth
        DO I=1,N

    C1=MGH1IHX(J)*(P0(I,J)-PB(I,J))
    A11(I,J)=-C1*AC(J)
    A22(I,J)= C1*BD(J)

    C2=MGH2IHY(J)*(P0(I,J)-PB(I,J))
    A12(I,J)=-C2*BC(J)
    A21(I,J)=-C2*AD(J)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
    B11(I,J)=-   E1(I,J)*AC(J)                  &
       &     -   E2(I,J)*BC(J)
    B22(I,J)=-   E2(I,J)*AD(J)                  &
       &     +   E1(I,J)*BD(J)
        enddo
      enddo


      DO J=M+1-DP_Depth,M
        DO I=1,N

    C1=MGH1IHX(J)*(P0(I,J)-PB(I,J))
    A11(I,J)=-C1*AC(J)
    A22(I,J)= C1*BD(J)

    C2=MGH2IHY(J)*(P0(I,J)-PB(I,J))
    A12(I,J)=-C2*BC(J)
    A21(I,J)=-C2*AD(J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
    B11(I,J)=-   E1(I,J)*AC(J)                  &
       &     -   E2(I,J)*BC(J)
    B22(I,J)=-   E2(I,J)*AD(J)                  &
       &     +   E1(I,J)*BD(J)
        enddo
      enddo

END SUBROUTINE


SUBROUTINE LAPL_depth(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP, num_of_bits, DP_Depth)
use implicit_functions_DP                        !! apply linear Operator to field P

implicit none
double precision :: P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M, num_of_bits, DP_Depth

double precision :: UTIL, VTIL
INTEGER :: I, J


 

If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
  end do
end do
  
DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
ENDDO  
CALL XBC(U,N,M)
CALL XBC(V,N,M)
else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo


      DO J=2,DP_Depth
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
      end do
CALL XBC(U,N,M)
CALL XBC(V,N,M)


endif

  
      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J)) 
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo
CALL XBC(U,N,M)
CALL XBC(V,N,M)



If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
  end do
end do    
 
DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
ENDDO
CALL XBC(F,N,M)

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO J=2,DP_Depth
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
      end do
  CALL XBC(F,N,M)

endif


END SUBROUTINE


SUBROUTINE LAPLfirst_depth(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP, num_of_bits, DP_Depth)
use implicit_functions_DP              !! apply linear Operator to field P: first iteration

implicit none
double precision :: P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M, num_of_bits, DP_Depth

double precision :: UTIL, VTIL
INTEGER :: I, J



If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
  end do
end do
  
DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
ENDDO  
CALL XBC(U,N,M)
CALL XBC(V,N,M)
else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo




      DO J=2,DP_Depth
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
        enddo
      enddo

      DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
      end do
CALL XBC(U,N,M)
CALL XBC(V,N,M)


endif



      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P0(I,J))  
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P0(I,J))  
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P0(I,J)) 
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
        enddo
      enddo
CALL XBC(U,N,M)
CALL XBC(V,N,M)



If (DP_Depth<=0) then

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
  end do
end do    
 
DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
ENDDO
CALL XBC(F,N,M)

else

      DO J=2+(DP_Depth-1),M-1-(DP_Depth-1)
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo


      DO J=2,DP_Depth
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M-1
        DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
        enddo
      enddo

      DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
      end do
  CALL XBC(F,N,M)

endif

END SUBROUTINE



SUBROUTINE LAPL(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)   !! apply linear Operator to field P
use implicit_functions_DP

implicit none
double precision :: P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M

double precision :: UTIL, VTIL
INTEGER :: I, J


DO J=2,M-1
  DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
  end do
end do
  
DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
ENDDO   

CALL XBC(U,N,M)
CALL XBC(V,N,M)



DO J=1,M
  DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J)) 
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
  ENDDO
ENDDO

CALL XBC(U,N,M)
CALL XBC(V,N,M)



DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
  end do
end do    
 

DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
ENDDO

CALL XBC(F,N,M)


END SUBROUTINE

SUBROUTINE LAPLfirst(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP) !! apply linear Operator to field P
                                                                  !! first iteration
use implicit_functions_DP

implicit none
double precision :: P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M

double precision :: UTIL, VTIL
INTEGER :: I, J


DO J=2,M-1
  DO I=2,N-1
    U(I,J)=P(I+1,J)-P(I-1,J)
    V(I,J)=P(I,J+1)-P(I,J-1)
  end do
end do
  
DO I=2,N-1
  U(I,1)=P(I+1,1)-P(I-1,1)
  U(I,M)=P(I+1,M)-P(I-1,M)
  V(I,1)=P(I,2)-P(IP(I),1)
  V(I,M)=P(IP(I),M)-P(I,M-1)
ENDDO   

CALL XBC(U,N,M)
CALL XBC(V,N,M)



DO J=1,M
  DO I=1,N
    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J)-P0(I,J)) 
    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J)-P0(I,J))
    U(I,J)=UTIL
    V(I,J)=VTIL
  ENDDO
ENDDO

CALL XBC(U,N,M)
CALL XBC(V,N,M)



DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
  end do
end do    
 

DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
ENDDO

CALL XBC(F,N,M)


END SUBROUTINE

!! advective velocities

SUBROUTINE VELPRD(U,V,F,G,PD,HX,HY,IP,N,M,A,B,EP)
use implicit_functions_DP

implicit none
INTEGER :: N, M
double precision :: U(N,M,0:1),V(N,M,0:1),F(N,M,0:1),G(N,M,0:1),         &
     &          PD(N,M),HX(N,M),HY(N,M)
INTEGER :: IP(N)
double precision :: EP, A, B

double precision :: UU(N,M),VV(N,M)
double precision :: CF, C1, C2, C1H, C2H, ALFA, BETA, ALF1, BET1, ALFM, BETM
INTEGER :: IORT, I, J


 IORT=1
 CF=rpe_05
 C1=A*CF
 C2=B*CF
 C1H=C1*rpe_05
 C2H=C2*rpe_05

!  COMPUTE V+R/PHI*DT FIELD FOR LAGRANGIAN ESTIMATES
 
DO J=1,M
  DO I=1,N
    F(I,J,1)=U(I,J,0)+CF*F(I,J,0)/max(PD(I,J),EP)
    G(I,J,1)=V(I,J,0)+CF*G(I,J,0)/max(PD(I,J),EP)
  enddo
end do

! COMPUTE U AND V TO FIRST ORDER

DO J=2,M-1
  DO I=2,N-1
    ALFA=U(I,J,0)/HX(I,J)*C1
    BETA=V(I,J,0)/HY(I,J)*C2
    U(I,J,1)=F(I,J,1)-max(rpe_0,ALFA)*(F(I,J,1)-F(I-1,J,1))   &
      &              -min(rpe_0,ALFA)*(F(I+1,J,1)-F(I,J,1))   &
      &              -max(rpe_0,BETA)*(F(I,J,1)-F(I,J-1,1))   &
      &              -min(rpe_0,BETA)*(F(I,J+1,1)-F(I,J,1)) 
    V(I,J,1)=G(I,J,1)-max(rpe_0,ALFA)*(G(I,J,1)-G(I-1,J,1))   &
      &              -min(rpe_0,ALFA)*(G(I+1,J,1)-G(I,J,1))   &
      &              -max(rpe_0,BETA)*(G(I,J,1)-G(I,J-1,1))   &
      &              -min(rpe_0,BETA)*(G(I,J+1,1)-G(I,J,1))
  end do
end do

DO I=2,N-1
  ALF1=U(I,1,0)/HX(I,1)*C1
  BET1=V(I,1,0)/HY(I,1)*C2
  U(I,1,1)=F(I,1,1)-max(rpe_0,ALF1)*(F(I,1,1)-F(I-1,1,1))       & 
     &                 -min(rpe_0,ALF1)*(F(I+1,1,1)-F(I,1,1))   &
     &                 -max(rpe_0,BET1)*(F(I,1,1)-F(IP(I),1,1)) &
     &                 -min(rpe_0,BET1)*(F(I,2,1)-F(I,1,1))
  V(I,1,1)=G(I,1,1)-max(rpe_0,ALF1)*(G(I,1,1)-G(I-1,1,1))   &
     &                 -min(rpe_0,ALF1)*(G(I+1,1,1)-G(I,1,1))   & 
     &                 -max(rpe_0,BET1)*(G(I,1,1)+G(IP(I),1,1)) &
     &                 -min(rpe_0,BET1)*(G(I,2,1)-G(I,1,1))

  ALFM=U(I,M,0)/HX(I,M)*C1
  BETM=V(I,M,0)/HY(I,M)*C2
  U(I,M,1)=F(I,M,1)-max(rpe_0,ALFM)*(F(I,M,1)-F(I-1,M,1))        &
     &                 -min(rpe_0,ALFM)*(F(I+1,M,1)-F(I,M,1))    &
     &                 -max(rpe_0,BETM)*(F(I,M,1)-F(I,M-1,1))    &
     &                 -min(rpe_0,BETM)*(F(IP(I),M,1)-F(I,M,1)) 
  V(I,M,1)=G(I,M,1)-max(rpe_0,ALFM)*(G(I,M,1)-G(I-1,M,1))        &
     &                 -min(rpe_0,ALFM)*(G(I+1,M,1)-G(I,M,1))    &
     &                 -max(rpe_0,BETM)*(G(I,M,1)-G(I,M-1,1))    &
     &                 +min(rpe_0,BETM)*(G(IP(I),M,1)+G(I,M,1))
end do

CALL XBC(U(:,:,1),N,M)
CALL XBC(V(:,:,1),N,M)

IF(IORT.EQ.2) THEN
 ! COMPUTE U AND V TO SEMI-SECOND ORDER

  DO J=1,M
    DO I=1,N
      UU(I,J)=rpe_05*(U(I,J,0)+U(I,J,1))
      VV(I,J)=rpe_05*(V(I,J,0)+V(I,J,1))
    ENDDO
  ENDDO

  DO J=2,M-1
    DO I=2,N-1
      ALFA=UU(I,J)/HX(I,J)*C1
      BETA=VV(I,J)/HY(I,J)*C2
      U(I,J,1)=F(I,J,1)-max(rpe_0,ALFA)*(F(I,J,1)-F(I-1,J,1))     &
       &                 -min(rpe_0,ALFA)*(F(I+1,J,1)-F(I,J,1))   &
       &                 -max(rpe_0,BETA)*(F(I,J,1)-F(I,J-1,1))   &
       &                 -min(rpe_0,BETA)*(F(I,J+1,1)-F(I,J,1))
      V(I,J,1)=G(I,J,1)-max(rpe_0,ALFA)*(G(I,J,1)-G(I-1,J,1))     &
       &                 -min(rpe_0,ALFA)*(G(I+1,J,1)-G(I,J,1))   &
       &                 -max(rpe_0,BETA)*(G(I,J,1)-G(I,J-1,1))   &
       &                 -min(rpe_0,BETA)*(G(I,J+1,1)-G(I,J,1))
    end do
  end do

  DO I=2,N-1
    ALF1=UU(I,1)/HX(I,1)*C1
    BET1=VV(I,1)/HY(I,1)*C2
    U(I,1,1)=F(I,1,1)-max(rpe_0,ALF1)*(F(I,1,1)-F(I-1,1,1))     &
     &                 -min(rpe_0,ALF1)*(F(I+1,1,1)-F(I,1,1))   &
     &                 -max(rpe_0,BET1)*(F(I,1,1)-F(IP(I),1,1)) &
     &                 -min(rpe_0,BET1)*(F(I,2,1)-F(I,1,1)) 
    V(I,1,1)=G(I,1,1)-max(rpe_0,ALF1)*(G(I,1,1)-G(I-1,1,1))     &
      &                -min(rpe_0,ALF1)*(G(I+1,1,1)-G(I,1,1))   & 
      &                -max(rpe_0,BET1)*(G(I,1,1)+G(IP(I),1,1)) &
      &                -min(rpe_0,BET1)*(G(I,2,1)-G(I,1,1))
    ALFM=UU(I,M)/HX(I,M)*C1
    BETM=VV(I,M)/HY(I,M)*C2

    U(I,M,1)=F(I,M,1)-max(rpe_0,ALFM)*(F(I,M,1)-F(I-1,M,1))     &
     &                 -min(rpe_0,ALFM)*(F(I+1,M,1)-F(I,M,1))   &
     &                 -max(rpe_0,BETM)*(F(I,M,1)-F(I,M-1,1))   &
     &                 -min(rpe_0,BETM)*(F(IP(I),M,1)-F(I,M,1))
    V(I,M,1)=G(I,M,1)-max(rpe_0,ALFM)*(G(I,M,1)-G(I-1,M,1))     &
     &                 -min(rpe_0,ALFM)*(G(I+1,M,1)-G(I,M,1))   &
     &                 -max(rpe_0,BETM)*(G(I,M,1)-G(I,M-1,1))   &
     &                 +min(rpe_0,BETM)*(G(IP(I),M,1)+G(I,M,1))
  end do

  CALL XBC(U(:,:,1),N,M)
  CALL XBC(V(:,:,1),N,M)
ENDIF

END SUBROUTINE


!! MPDATA subroutine
SUBROUTINE MPDATT(U1,U2,X,H,N,M,IORD,ISOR,NONOS,IDIV,IBC, IP, X_T, codes)

use implicit_functions_DP
implicit none

  INTEGER, PARAMETER  :: LW=1,MP=1-LW

  INTEGER :: N, M

  double precision :: U1(N,M),U2(N,M+1),X(N,M),H(N,M), X_T(N,M)
  INTEGER :: IP(N)
  INTEGER  :: IORD, ISOR, NONOS, IDIV, IBC 
      

  double precision :: V1(N,M),V2(N,M+1),F1(N,M),F2(N,M+1)  &
      &      ,CP(N,M),CN(N,M)   
  double precision :: MX(N,M),MN(N,M)
  double precision :: EP, C1, C2, V1D, V2D, V2D1, V2DN
  INTEGER :: N1, N2, I, J, K
  LOGICAL :: codes

  N1=N
  N2=M

  EP= 1.E-10


  IF(ISOR.EQ.3) IORD=MAX(IORD,3)
  IF(LW.EQ.1) IORD=2

      !! take predicted advective velocities
  DO J=1,N2
    DO I=1,N1
     V1(I,J)=U1(I,J)
    end do
  end do
  
  DO J=1,N2+1
    DO I=1,N1
     V2(I,J)=U2(I,J)
    end do
  end do

  
  IF(NONOS.EQ.1) THEN           
    DO J=2,N2-1
      DO I=2,N1-1
        MX(I,J)=max(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
        MN(I,J)=min(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
      end do
    end do
    
    DO I=2,N1-1
      MX(I,1)=max(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MN(I,1)=min(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MX(I,N2)=max(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                &
         &                         IBC*X(IP(I),N2))
      MN(I,N2)=min(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                &
         &                         IBC*X(IP(I),N2))
    end do  
      
    CALL XBC(MX,N1,N2)
    CALL XBC(MN,N1,N2)
  ENDIF

  
  C1=rpe_1
  C2=rpe_0
  
  DO K=1,IORD   !! k-th order correction
   
   ! COMPUTE DONOR-CELL FLUXES
    DO J=1,N2
      DO I=2,N1-1
        F1(I,J)=DONOR(C1*X(I-1,J)+C2,C1*X(I,J)+C2,V1(I,J))
      end do
    end do
    
    DO J=2,N2
      DO I=2,N1-1
        F2(I,J)=DONOR(C1*X(I,J-1)+C2,C1*X(I,J)+C2,V2(I,J))
      end do
    end do
      
    DO I=2,N1-1
      F2(I,N2+1)=DONOR(C1*X(I,N2)+C2,C1*IBC*X(IP(I),N2)+C2,V2(I,N2+1))
      F2(I,1)=DONOR(C1*IBC*X(IP(I),1)+C2,C1*X(I,1)+C2,V2(I,1))
    end do
    
    CALL XBC(F1,N1,N2)
    CALL XBC(F2,N1,N2+1)
   ! COMPUTE NEW UPWIND-SOLUTION
    
    DO J=1,N2
      DO I=2,N1-1
        X(I,J)=X(I,J)-(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
      end do
    end do
    
    
    IF (codes) then
      DO J=1,N2
        DO I=2,N1-1
          X_T(I,J)=X_T(I,J) -(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
        end do
      end do
      
      CALL XBC(X_T,N1,N2)
    end if
    
    
    CALL XBC(X,N1,N2)

    IF(.not. (K.EQ.IORD)) then  ! GO TO 6
    
      C1=FLOAT(MP)
      C2=FLOAT(LW)
   ! CONVERT VELOCITIES TO LOCAL STORAGE
      DO J=1,N2
        DO I=1,N1
          F1(I,J)=V1(I,J)
        end do
      end do
      
      DO J=1,N2+1
        DO I=1,N1
          F2(I,J)=V2(I,J) 
        end do
      end do 
   ! CALCULATE PSEUDO VELOCITIES      
   ! COMPUTE FIRST DIRECTION    
      
      DO J=2,N2-1
        DO I=2,N1-1
          V1(I,J)=VDYF_D(X(I-1,J),X(I,J),F1(I,J),rpe_05*(H(I-1,J)+H(I,J)),MP, EP)     &
               & +VCORR_D(F1(I,J), F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),                 &
               &   X(I-1,J-1),X(I,J-1),X(I-1,J+1),X(I,J+1),                   &
               &       rpe_05*(H(I-1,J)+H(I,J)),MP, EP)
        end do
      end do
   ! COMPUTE B.C IN Y DIRECTION
      DO I=2,N1-1
        V1(I,1)=VDYF_D(X(I-1,1),X(I,1),F1(I,1),rpe_05*(H(I-1,1)+H(I,1)),MP, EP)         &
           & +VCORR_D(F1(I,1), F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),                 &
           & IBC*X(IP(I-1),1),IBC*X(IP(I),1),X(I-1,2),X(I,2),                       &
           &           rpe_05*(H(I-1,1)+H(I,1)),MP, EP)
        V1(I,N2)=VDYF_D(X(I-1,N2),X(I,N2),F1(I,N2),rpe_05*(H(I-1,N2)+H(I,N2)),MP, EP)     &
           & +VCORR_D(F1(I,N2), F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),              &
           &   X(I-1,N2-1),X(I,N2-1),IBC*X(IP(I-1),N2),IBC*X(IP(I),N2),              &
           &           rpe_05*(H(I-1,N2)+H(I,N2)),MP, EP)        
      end do
      
      IF(IDIV.EQ.1) THEN
 
 
    ! COMPUTE FLOW-DIVERGENCE CORRECTION
        DO J=1,N2
          DO I=2,N1-1
            V1D=-VDIV1(F1(I-1,J),F1(I,J),F1(I+1,J),rpe_05*(H(I-1,J)+H(I,J)))              &
              & -VDIV2(F1(I,J),F2(I-1,J+1),F2(I,J+1),F2(I-1,J),F2(I,J),               &
              &        rpe_05*(H(I-1,J)+H(I,J)))                            
            V1(I,J)=V1(I,J)+LW*(PP(V1D)*X(I-1,J)-PN(V1D)*X(I,J))+MP*V1D
          end do
        end do
        
      ENDIF   

      
    ! COMPUTE SECOND DIRECTION
      DO J=2,N2
        DO I=2,N1-1
          V2(I,J)=VDYF_D(X(I,J-1),X(I,J),F2(I,J),rpe_05*(H(I,J-1)+H(I,J)),MP, EP)         &
             & +VCORR_D(F2(I,J), F1(I,J-1)+F1(I,J)+F1(I+1,J)+F1(I+1,J-1),             &
             &          X(I-1,J-1),X(I-1,J),X(I+1,J-1),X(I+1,J),                      &
             &         rpe_05*(H(I,J-1)+H(I,J)),MP, EP)
        end do
      end do

    ! COMPUTE B.C IN Y-DIRECTION
      DO I=2,N1-1
        V2(I,1)=VDYF_D(IBC*X(IP(I),1),X(I,1),F2(I,1),rpe_05*(H(IP(I),1)+H(I,1)),MP, EP)   &
          & +VCORR_D(F2(I,1), F1(IP(I),1)+F1(I,1)+F1(I+1,1)+F1(IP(I+1),1),            &
          &    IBC*X(IP(I-1),1),X(I-1,1),X(I+1,1),IBC*X(IP(I+1),1),                   &
          &            rpe_05*(H(IP(I),1)+H(I,1)),MP, EP)
        V2(I,N2+1)=VDYF_D(X(I,N2),IBC*X(IP(I),N2),F2(I,N2+1),                         &
          &                   rpe_05*(H(I,N2)+H(IP(I),N2)),MP, EP)                        &
          & +VCORR_D(F2(I,N2+1),F1(I,N2)+F1(IP(I),N2)+F1(IP(I+1),N2)+F1(I+1,N2),      &
          & X(I-1,N2),IBC*X(IP(I-1),N2),X(I+1,N2),IBC*X(IP(I+1),N2),                  &
          &                   rpe_05*(H(I,N2)+H(IP(I),N2)),MP, EP)
      end do
      
      IF(IDIV.EQ.1) THEN
        DO J=2,N2
          DO I=2,N1-1
            V2D=-VDIV1(F2(I,J-1),F2(I,J),F2(I,J+1),rpe_05*(H(I,J-1)+H(I,J)))              &
              & -VDIV2(F2(I,J),F1(I+1,J-1),F1(I+1,J),F1(I,J-1),F1(I,J),               &
              &   rpe_05*(H(I,J-1)+H(I,J)))
            V2(I,J)=V2(I,J)+LW*(PP(V2D)*X(I,J-1)-PN(V2D)*X(I,J))+MP*V2D
          end do
        end do
        
        DO I=2,N1-1
          V2D1=-VDIV1(-F2(IP(I),1),F2(I,1),F2(I,2),rpe_05*(H(IP(I),1)+H(I,1)))            &
            & -VDIV2( F2(I,1),F1(IP(I+1),1),F1(I+1,1),F1(IP(I),1),F1(I,1),            &
            &      rpe_05*(H(IP(I),1)+H(I,1))) 
          V2(I,1)=V2(I,1)                                                             &
            & +LW*(PP(V2D1)*IBC*X(IP(I),1)-PN(V2D1)*X(I,1))+MP*V2D1
          V2DN=-VDIV1(F2(I,N2),F2(I,N2+1),-F2(IP(I),N2+1)                             &
            &              ,rpe_05*(H(I,N2)+H(IP(I),N2)))                                 &
            &  -VDIV2(F2(I,N2+1),F1(I+1,N2),F1(IP(I+1),N2),                           &
            &     F1(I,N2),F1(IP(I),N2),rpe_05*(H(I,N2)+H(IP(I),N2))) 
          V2(I,N2+1)=V2(I,N2+1)                                                       &
            &  +LW*(PP(V2DN)*X(I,N2)-PN(V2DN)*IBC*X(IP(I),N2))+MP*V2DN
        end do
       
      ENDIF
      
!
    ! THIRD ORDER CORRECTION
      
      IF(ISOR.EQ.3) THEN
      
    ! FIRST DIRECTION
        DO J=1,N2
          DO I=3,N1-1
            V1(I,J)=V1(I,J)     +VCOR31_D(F1(I,J),                                    &
              &  X(I-2,J),X(I-1,J),X(I,J),X(I+1,J),rpe_05*(H(I-1,J)+H(I,J)),MP, EP)
          end do
   
          V1(2,J)=V1(2,J)     +VCOR31_D(F1(2,J),                                      &
            &  X(N1-2,J),X(1,J),X(2,J),X(3,J),rpe_05*(H(1,J)+H(2,J)),MP, EP)
        end do
        
        DO J=2,N2-1
          DO I=2,N1-1               
            V1(I,J)=V1(I,J)                                                           &
              &  +VCOR32_D(F1(I,J),F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),           &
              &   X(I,J-1),X(I-1,J+1),X(I-1,J-1),X(I,J+1),                            &
              &  rpe_05*(H(I-1,J)+H(I,J)),MP, EP)
          end do
        end do
    ! C B.C. FOLLOW
        DO I=2,N1-1
          V1(I,1)=V1(I,1)                                                             &
            & +VCOR32_D(F1(I,1),F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),                  &
            & IBC*X(IP(I),1),X(I-1,2),X(I,2),IBC*X(IP(I-1),1),                        &
            & rpe_05*(H(I-1,1)+H(I,1)),MP, EP)
          V1(I,N2)=V1(I,N2)                                                           &
            & +VCOR32_D(F1(I,N2),F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),         &
            & X(I,N2-1),IBC*X(IP(I-1),N2),IBC*X(IP(I),N2),X(I-1,N2-1),                &
            & rpe_05*(H(I-1,N2)+H(I,N2)),MP, EP)
        end do
      
        DO J=1,N2
          V1(1,J)=V1(N1-1,J)
          V1(N1,J)=V1(2,J)
        end do

    !  SECOND DIRECTION
        
        DO J=3,N2-1
          DO I=2,N1-1
            V2(I,J)=V2(I,J)     +VCOR31_D(F2(I,J),                                    &
              &  X(I,J-2),X(I,J-1),X(I,J),X(I,J+1),rpe_05*(H(I,J-1)+H(I,J)),MP, EP)
          end do
        end do
        
        DO I=2,N1-1
          V2(I,1)=V2(I,1)     +VCOR31_D(F2(I,1),IBC*X(IP(I),2),                       &
            &  IBC*X(IP(I),1),X(I,1),X(I,2),rpe_05*(H(IP(I),1)+H(I,1)),MP, EP)  
          V2(I,2)=V2(I,2)     +VCOR31_D(F2(I,2),                                      &
            &  IBC*X(IP(I),1),X(I,1),X(I,2),X(I,3),rpe_05*(H(I,1)+H(I,2)),MP, EP)   
          V2(I,N2)=V2(I,N2)     +VCOR31_D(F2(I,N2),X(I,N2-2),                         &
            &  X(I,N2-1),X(I,N2),IBC*X(IP(I),N2),rpe_05*(H(I,N2-1)+H(I,N2)),MP, EP)  
          V2(I,N2+1)=V2(I,N2+1) +VCOR31_D(F2(I,N2+1), X(I,N2-1),X(I,N2),              &
            &  IBC*X(IP(I),N2),IBC*X(IP(I),N2-1),rpe_05*(H(I,N2)+H(IP(I),N2)),MP, EP)  
        end do
        
        DO J=2,N2
          DO I=2,N1-1
            V2(I,J)=V2(I,J)                                                           &
              & +VCOR32_D(F2(I,J),F1(I,J-1)+F1(I+1,J-1)+F1(I+1,J)+F1(I,J),            &
              &  X(I+1,J-1),X(I-1,J),X(I-1,J-1),X(I+1,J),                             &
              &         rpe_05*(H(I,J-1)+H(I,J)),MP, EP)
          end do
        end do
        
        
      ENDIF  !! end third order correction

     ! CALL B.C IN X DIRECTION
      CALL XBC(V1,N1,N2)
      CALL XBC(V2,N1,N2+1)

      
      
      IF(.not. (NONOS.EQ.0)) then  !! if non-osc is not turned off
     ! NON-OSCILLATORY OPTION
        DO J=2,N2-1
          DO I=2,N1-1
            MX(I,J)=max(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MX(I,J))
            MN(I,J)=min(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MN(I,J))
          end do
        end do
      
        DO I=2,N1-1
          MX(I,1)=max(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),                      &
             &                X(I,2),MX(I,1))
          MN(I,1)=min(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),                      &
             &                X(I,2),MN(I,1))
          MX(I,N2)=max(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                       &
             &                       IBC*X(IP(I),N2),MX(I,N2))
          MN(I,N2)=min(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),                       &
             &                      IBC*X(IP(I),N2),MN(I,N2))
        end do
        
        CALL XBC(MX,N1,N2)
        CALL XBC(MN,N1,N2)

        DO J=1,N2
          DO I=2,N1-1
            F1(I,J)=DONOR(C1*X(I-1,J)+C2,C1*X(I,J)+C2,V1(I,J))
          end do
        end do
      
        DO J=2,N2
          DO I=2,N1-1
            F2(I,J)=DONOR(C1*X(I,J-1)+C2,C1*X(I,J)+C2,V2(I,J))
          end do
        end do
      
        DO I=2,N1-1    
          F2(I,N2+1)=DONOR(C1*X(I,N2)+C2,C1*IBC*X(IP(I),N2)+C2,V2(I,N2+1))
          F2(I,1)=DONOR(C1*IBC*X(IP(I),1)+C2,C1*X(I,1)+C2,V2(I,1))
        end do
      
        CALL XBC(F1,N1,N2)
        CALL XBC(F2,N1,N2+1)

        DO J=1,N2
          DO I=2,N1-1
            CP(I,J)=(MX(I,J)-X(I,J))*H(I,J)/                                     &
               & (PN(F1(I+1,J))+PP(F1(I,J))+PN(F2(I,J+1))+PP(F2(I,J))+EP)
            CN(I,J)=(X(I,J)-MN(I,J))*H(I,J)/                                     &                        
               & (PP(F1(I+1,J))+PN(F1(I,J))+PP(F2(I,J+1))+PN(F2(I,J))+EP)
          end do
        end do
        
        CALL XBC(CP,N1,N2)
        CALL XBC(CN,N1,N2)

        IF(LW.EQ.0) THEN
          
          DO J=2,N2
            DO I=2,N1-1
              V1(I,J)=PP(V1(I,J))*                                              &
                & ( min(rpe_1,CP(I,J),CN(I-1,J))*PP(SIGN(rpe_1, X(I-1,J)))          &
                &  +min(rpe_1,CP(I-1,J),CN(I,J))*PP(SIGN(rpe_1,-X(I-1,J))) )         &
                & -PN(V1(I,J))*                                                 &
                & ( min(rpe_1,CP(I-1,J),CN(I,J))*PP(SIGN(rpe_1, X(I ,J )))          &
                &  +min(rpe_1,CP(I,J),CN(I-1,J))*PP(SIGN(rpe_1,-X(I ,J ))) )     
              V2(I,J)=PP(V2(I,J))*                                              &
                & ( min(rpe_1,CP(I,J),CN(I,J-1))*PP(SIGN(rpe_1, X(I,J-1)))          &
                &  +min(rpe_1,CP(I,J-1),CN(I,J))*PP(SIGN(rpe_1,-X(I,J-1))) )         &
                & -PN(V2(I,J))*                                                 &
                & ( min(rpe_1,CP(I,J-1),CN(I,J))*PP(SIGN(rpe_1, X(I ,J )))          &
                &  +min(rpe_1,CP(I,J),CN(I,J-1))*PP(SIGN(rpe_1,-X(I ,J ))) )     
             end do
           end do
       ! B.C. FOLLOW
       
          DO I=2,N1-1
            V2(I,1)=PP(V2(I,1))*                                              &
              & ( min(rpe_1,CP(I,1),CN(IP(I),1))*PP(SIGN(rpe_1, IBC*X(IP(I),1)))    &
              & +min(rpe_1,CP(IP(I),1),CN(I,1))*PP(SIGN(rpe_1,-IBC*X(IP(I),1))) )   &
              & -PN(V2(I,1))*                                                 &
              & ( min(rpe_1,CP(IP(I),1),CN(I,1))*PP(SIGN(rpe_1, X(I ,1 )))          &
              & +min(rpe_1,CP(I,1),CN(IP(I),1))*PP(SIGN(rpe_1,-X(I ,1 ))) )       
            V2(I,N2+1)=PP(V2(I,N2+1))*                                          &
              & ( min(rpe_1,CP(IP(I),N2),CN(I,N2))*PP(SIGN(rpe_1, X(I,N2)))         & 
              & +min(rpe_1,CP(I,N2),CN(IP(I),N2))*PP(SIGN(rpe_1,-X(I,N2))) )        &
              & -PN(V2(I,N2+1))*                                                &
              & ( min(rpe_1,CP(I,N2),CN(IP(I),N2))*PP(SIGN(rpe_1, IBC*X(IP(I),N2))) &
              & +min(rpe_1,CP(IP(I),N2),CN(I,N2))*PP(SIGN(rpe_1,-IBC*X(IP(I),N2)))) 
            V1(I,1)=PP(V1(I,1))*                                                & 
              & ( min(rpe_1,CP(I,1),CN(I-1,1))*PP(SIGN(rpe_1, X(I-1,1)))            &
              & +min(rpe_1,CP(I-1,1),CN(I,1))*PP(SIGN(rpe_1,-X(I-1,1))) )           &
              & -PN(V1(I,1))*                                                   &
              & ( min(rpe_1,CP(I-1,1),CN(I,1))*PP(SIGN(rpe_1, X(I ,1 )))            &
              & +min(rpe_1,CP(I,1),CN(I-1,1))*PP(SIGN(rpe_1,-X(I ,1 ))) )
          end do
          
        ELSE
     
          DO J=2,N2
            DO I=2,N1-1
              V1(I,J)=PP(V1(I,J))*min(rpe_1,CP(I,J),CN(I-1,J))                   &
                & -PN(V1(I,J))*min(rpe_1,CP(I-1,J),CN(I,J))
              V2(I,J)=PP(V2(I,J))*min(rpe_1,CP(I,J),CN(I,J-1))                   & 
                & -PN(V2(I,J))*min(rpe_1,CP(I,J-1),CN(I,J))
            end do
          end do
      ! B.C. FOLLOW
      
          DO I=2,N1-1
            V2(I,1)=PP(V2(I,1))*min(rpe_1,CP(I,1),CN(IP(I),1))                   &
              & -PN(V2(I,1))*min(rpe_1,CP(IP(I),1),CN(I,1))
            V2(I,N2+1)=PP(V2(I,N2+1))*min(rpe_1,CP(IP(I),N2),CN(I,N2))           &
              & -PN(V2(I,N2+1))*min(rpe_1,CP(I,N2),CN(IP(I),N2))
            V1(I,1)=PP(V1(I,1))*min(rpe_1,CP(I,1),CN(I-1,1))                     &
              & -PN(V1(I,1))*min(rpe_1,CP(I-1,1),CN(I,1))
          end do

        ENDIF
      
        CALL XBC(V1,N1,N2)
        CALL XBC(V2,N1,N2+1)
!
!         END OF NONOSCILLATORY OPTION
!
      end if
      
    end if ! IF last iteration exit loop
  
  end do
  
END SUBROUTINE   

  
!! periodic boundary condition
SUBROUTINE XBC(X,N,M)
use implicit_functions_DP

implicit none

INTEGER :: N, M
double precision :: X(N,M)

INTEGER :: J

DO J=1,M
  X(1,J)=X(N-1,J)
  X(N,J)=X(2,J)
end do

END subroutine

!! periodic boundary condition: same as XBC
SUBROUTINE XBC_52(X,N,M)
use implicit_functions_DP

implicit none

INTEGER :: N, M
double precision :: X(N,M)

INTEGER :: J

DO J=1,M
  X(1,J)=X(N-1,J)
  X(N,J)=X(2,J)
end do

END subroutine


SUBROUTINE DIAGNOS(U,V,PD,PT,HX,HY,IP,S,TIME,DX,DY,DT, SUM0,SUM1, &
         & KT,N,M, IFLG, NITER,NITSM,ICOUNT,ERROR, sum_time, sum_lp_time)

use implicit_functions_DP

implicit none

INTEGER :: N, M

double precision :: U(N,M),V(N,M),PD(N,M),PT(N,M),HX(N,M),HY(N,M),  &
     &        S(N,M)
double precision :: TIME,DX,DY,DT, SUM0,SUM1,ERROR
INTEGER  :: IP(N)
INTEGER  :: KT, IFLG, NITER,NITSM,ICOUNT
double precision :: avg_time, sum_time, avg_lp_time, sum_lp_time, NITAV
 double precision :: GC1, GC2, COUR1, COUR2, PDMX,PDMN,PDAV, SUMER, DLI
INTEGER ::  I, J


GC1=DT/DX
GC2=DT/DY

IF(IFLG.EQ.0) THEN

  KT=0
  NITER=0
  ERROR=0.0d0
  SUM0=0.0d0

  DO J=1,M
    DO I=2,N-1
      SUM0=SUM0+PD(I,J)*S(I,J)
    end do
  end do       

ENDIF

TIME=KT*DT/3600.0
      PRINT 299
  299 FORMAT(1X,1H )
      PRINT *, 'HTIME= ', TIME,'KT= ',KT
  300 FORMAT(14X,5HTIME=,F7.2,7H HOURS;,5H  KT=,I5)

  ! CHECK COURANT NUMBERS

 COUR1=0.0d0
 COUR2=0.0d0


DO J=1,M
  DO I=2,N
    DLI=SQRT((GC1/HX(I,J))**2+(GC2/HY(I,J))**2)
    COUR1=max(COUR1,DLI*SQRT(abs(PT(I,J))))
    COUR2=max(COUR2, GC1*ABS(U(I,J)/HX(I,J))                          &
      &                  +GC2*ABS(V(I,J)/HY(I,J)) )
   ! COUR2=max(COUR2, max(DT*ABS(U(I,J)/ (HX(I,J)/FLOAT(N-2)) ), & !! COURANT NUMBER
   !   &             DT*ABS(V(I,J)/ (HY(I,J)/FLOAT(M)  ) ) ) &
   !   &       )
  end do
end do

      PRINT *,'COUR1,COUR2: ', COUR1,COUR2
  301 FORMAT(4X,'COUR1,COUR2:',2E11.4)

SUM1=0.0d0
PDMX=-1.E30
PDMN= 1.E30
PDAV=0.0d0

DO J=1,M
  DO I=2,N-1
    PDMX=max(PDMX,PD(I,J))
    PDMN=min(PDMN,PD(I,J))
    PDAV=PDAV+PD(I,J)
    SUM1=SUM1+PD(I,J)*S(I,J)
  end do
end do

PDAV=PDAV/FLOAT(M*(N-2))
      write(*,*) 'SUM1, SUM0, MAX(SUM0,PDMN)', SUM1, SUM0, MAX(SUM0,PDMN)
SUMER=(SUM1-SUM0)/max(SUM0,PDMN)

      PRINT *, 'PDMX,PDMN,PDAV: ', PDMX,PDMN,PDAV,' SUMER: ',SUMER

avg_time=sum_time/MAX(ICOUNT,1)
avg_lp_time=sum_lp_time/MAX(ICOUNT,1)
NITAV=float(NITSM)/MAX(ICOUNT,1)

      write(*,*) 'ERROR:', ERROR,'NITER, NITAV (GCR ITERATIONS): ',NITER,NITAV
      write(*,*) 'Computing time per implicit solve, low precision:', avg_time, avg_lp_time
END subroutine



subroutine  init_perf_markers(H_rpe,U_rpe,V_rpe, TIME_rpe, codesignQ, codesignD, &
                    &   IRHW, X_rpe, Y_rpe, N, M, bits, ID_prec, EXP_NAME)
use implicit_functions_DP

double precision :: H_rpe(N, M), U_rpe(N, M), V_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)

integer :: N, M, IRHW, bits
double precision :: H(N, M), U(N, M), V(N, M), X(N), Y(M), TIME
logical :: codesignQ,codesignD, itsopen

 character(len=150) path, file_name, experiment, simtime, codesQ, bit_count, Precond, EXP_NAME, codesD
INTEGER I,J


H(:, :)=H_rpe(:, :)
U(:, :)=U_rpe(:, :)
V(:, :)=V_rpe(:, :)
TIME=TIME_rpe

path ='../'//trim(adjustl(EXP_NAME))//'/Precon'
write(experiment,*) IRHW
write(Precond,*) ID_prec
write(codesQ,*) codesignQ
write(codesD,*) codesignD
write(bit_count,*) bits

inquire(unit=324, opened=itsopen)
If(.not. itsopen) then 
  file_name = trim(adjustl(Precond))//'_EXP_'//trim(adjustl(experiment)) &
      &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'.txt'
  Open(unit=324, file=trim(path)//trim(file_name), status='replace', form = 'formatted')
 write(unit=324,fmt='(I5)') bits
else

 write(unit=324,fmt=*)
 write(unit=324,fmt='(I5)') bits
endif

end subroutine init_perf_markers

subroutine  close_perf_markers()
use implicit_functions_DP

 close(Unit=324)

end subroutine close_perf_markers




subroutine write_perf_markers(H_rpe,U_rpe,V_rpe, TIME_rpe,  codesignQ, codesignD, IRHW, X_rpe, Y_rpe, &
                         & N, M, bits,NITER,NITSM,ICOUNT, sum_time, sum_lp_time)
use implicit_functions_DP

double precision :: H_rpe(N, M), U_rpe(N, M), V_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)

integer :: N, M, IRHW, bits,NITER,NITSM,ICOUNT
double precision :: H(N, M), U(N, M), V(N, M), X(N), Y(M), TIME, sum_time, sum_lp_time
logical ::  codesignQ, codesignD

 character(len=150) path, file_name, experiment, simtime, codesD, codesQ, bit_count
INTEGER I,J


H(:, :)=H_rpe(:, :)
U(:, :)=U_rpe(:, :)
V(:, :)=V_rpe(:, :)
TIME=TIME_rpe
X(:)=X_rpe(:)
Y(:)=Y_rpe(:)

    write(unit=324,fmt='(I5, F6.2, F8.4, F8.4)') int(time), float(NITSM)/MAX(ICOUNT,1), &
                        & sum_time/MAX(ICOUNT,1), sum_lp_time/MAX(ICOUNT,1)

sum_lp_time=0.0d0
sum_time=0.0d0
NITAV=0
NITSM=0
ICOUNT=0

end subroutine write_perf_markers

subroutine  write_residual(R_rpe, exitcond, iteration, TIME_rpe,  codesignQ, codes, &
                        &  IRHW, X_rpe, Y_rpe, N, M, bits, ID_prec,EXP_NAME)
use implicit_functions_DP



integer :: N, M, IRHW, bits, iteration
double precision :: R_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)
double precision :: R(N, M), X(N), Y(M), TIME, exitcond
logical :: codesignQ, codes

 character(len=150) path, file_name, experiment, simtime, codesQ, codesD, bit_count, Precond, EXP_NAME, str_iter
INTEGER I,J

R(:,:) = R_rpe(:, :)
TIME=TIME_rpe
X(:)=X_rpe(:)
Y(:)=Y_rpe(:)


write(experiment,*) IRHW
write(simtime,*) INT(TIME)
write(codesQ,*) codesignQ
write(codesD,*) codes
write(bit_count,*) bits
write(Precond,*) ID_prec
write(str_iter,*) iteration

path ='../'//trim(adjustl(EXP_NAME))//& 
                        & '/Precon'//trim(adjustl(Precond))

file_name = '_R_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))//&
     & '_iter_'//trim(adjustl(str_iter))      &
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599,fmt=*) X(I), Y(J), R(I,J), exitcond
  end do
end do
 close(unit=599)


end subroutine



subroutine  write_fields(H_rpe,U_rpe,V_rpe, TIME_rpe,  codesignQ, codesignD, IRHW, X_rpe, Y_rpe, N, M, bits, ID_prec,EXP_NAME)
use implicit_functions_DP

double precision :: H_rpe(N, M), U_rpe(N, M), V_rpe(N, M), TIME_rpe, X_rpe(N), Y_rpe(M)

integer :: N, M, IRHW, bits
double precision :: H(N, M), U(N, M), V(N, M), X(N), Y(M), TIME
logical :: codesignQ, codesignD

 character(len=150) path, file_name, experiment, simtime, codesQ, codesD, bit_count, Precond, EXP_NAME
INTEGER I,J

H(:, :)=H_rpe(:, :)
U(:, :)=U_rpe(:, :)
V(:, :)=V_rpe(:, :)
TIME=TIME_rpe
X(:)=X_rpe(:)
Y(:)=Y_rpe(:)


write(experiment,*) IRHW
write(simtime,*) INT(TIME)
write(codesQ,*) codesignQ
write(codesD,*) codesignD
write(bit_count,*) bits
write(Precond,*) ID_prec

path ='../'//trim(adjustl(EXP_NAME))//& 
                        & '/Precon'//trim(adjustl(Precond))

file_name = '_H_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))&
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599,fmt=*) X(I), Y(J), H(I,J)
  end do
end do
 close(unit=599)

file_name = '_U_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))&
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599, fmt=*) X(I), Y(J), U(I,J)
  end do
end do

 close(unit=599)


file_name = '_V_exp'//trim(adjustl(experiment))//'_time'//trim(adjustl(simtime))&
     &  //'_codes_'//trim(adjustl(codesQ))//trim(adjustl(codesD))//'_bits'//trim(adjustl(bit_count))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599, fmt=*) X(I), Y(J), V(I,J)
  end do
end do

 close(unit=599)


end subroutine

subroutine  write_MLfields(H_rpe, TIMESTEP, N, M, EXP_NAME, field)
use implicit_functions_DP

double precision :: H_rpe(N, M)

integer :: N, M, TIMESTEP
double precision :: H(N, M)


 character(len=150) path, file_name, simtime, EXP_NAME
 character(len=*) field
INTEGER I,J

H(:, :)=H_rpe(:, :)


write(simtime,*) TIMESTEp

path ='../'//trim(adjustl(EXP_NAME))//& 
                        & '/'//trim(adjustl(field))//'/'


file_name = 'Timestep'//trim(adjustl(simtime))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599, fmt=*) H(I,J)
  end do
end do

 close(unit=599)


end subroutine

subroutine  write_MLfields_iter(H_rpe, TIMESTEP, N, M, EXP_NAME, field, iteration)
use implicit_functions_DP

double precision :: H_rpe(N, M)

integer :: N, M, TIMESTEP, iteration
double precision :: H(N, M)


 character(len=150) path, file_name, simtime, EXP_NAME, str_iteration
 character(len=*) field
INTEGER I,J

H(:, :)=H_rpe(:, :)


write(simtime,*) TIMESTEp
write(str_iteration,*) iteration

path ='../'//trim(adjustl(EXP_NAME))//& 
                        & '/'//trim(adjustl(field))//'/'


file_name = 'Timestep'//trim(adjustl(simtime))//'_iter'//trim(adjustl(str_iteration))//'.txt'
Open(unit=599, file=trim(path)//trim(file_name), status='replace', form = 'formatted')

do I=2,N-1
  do J=1,M
    write(unit=599, fmt=*) H(I,J)
  end do
end do

 close(unit=599)


end subroutine


subroutine read_layer(filename_NN,Layer_NN,input_neurons,output_neurons)
character(len=150) :: filename_NN
integer :: input_neurons, output_neurons
double precision :: Layer_NN(input_neurons, output_neurons)

Open(unit=599, file=filename_NN, status='old', form = 'formatted')

do I=1,input_neurons
  !do J=1,output_neurons
    read(unit=599, fmt=*) Layer_NN(I,:)
    !write(*,*) I, J, Layer_NN(I,J)
  !end do
end do

 close(unit=599)


end subroutine


!! preconditioned Linear Solver
subroutine GCR_PRE(p,pfx,pfy,hx,hy,s,b,p0,pb,e1,e2,cor,ip  &
              & ,d,q,r,ar,n1,n2,gc1,gc2, &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
           &  niter,nitsm,icount,error, p_T, sum_time,&
           &  sum_lp_time,ID_PREC, codes, save_time , &
           & TIME, codesQ, IRHW, DX_rpe, DY_rpe, Exit_cond, EXP_NAME&
           & , iprint, num_of_bits, DP_Depth, KT, Layer1, Layer2)
use implicit_functions_DP

implicit none
INTEGER, parameter :: kord=2, lord=kord-1

INTEGER :: n1, n2, KT

double precision :: p(n1,n2),pfx(n1,n2),pfy(n1,n2),hx(n1,n2),hy(n1,n2),s(n1,n2), &
     &   b(n1,n2),pb(n1,n2),p0(n1,n2),                         &
     &   e1(n1,n2),e2(n1,n2),cor(n1,n2),d(n1,n2),q(n1,n2),r(n1,n2),ar(n1,n2), &
     &   p_T(n1,n2), r_HP(n1,n2), r_true(n1,n2), r0_true(n1,n2), p_true(n1,n2), &
     &   p0_true(n1, n2), b_true(n1, n2), PMB(n1, n2), PMP0(n1, n2)
double precision :: MGH1IHX(n2), MGH2IHY(n2), AC(n2), BC(n2), AD(n2), BD(n2)
!! preconditioning
double precision :: qu(n1,n2), aqu(n1,n2),  A_c(n1,n2), B_c(n1,n2), C_c(n1,n2), ps(n1+1,n2), divi(n1,n2)
INTEGER :: ID_PREC
!! end preconditioning
INTEGER :: IP(n1)
double precision :: GC1, GC2,error, max_QX_QY, epa
double precision :: res_lats0(n2), res_lats(n2)
double precision :: start, finish, sum_time, sum_lp_time, startLP, endLP
integer :: num_of_bits

INTEGER :: niter,nitsm,icount
INTEGER :: iprint, DP_Depth
LOGICAL :: codes, save_time


double precision :: x(n1,n2,lord),ax(n1,n2,lord),ax2(lord),axaqu(lord),del(lord),  &
     & a11(n1,n2),a12(n1,n2),a21(n1,n2),a22(n1,n2),b11(n1,n2),b22(n1,n2)
double precision :: err0, rax, beta, errn, x2, y2, T_step
INTEGER :: itr, J, I, l, ll, i1, it
double precision :: eps, help1, help2, quotient, lowprectime, err_true, err0_true
double precision :: Exit_cond

double precision :: TIME, DX_rpe(n1), DY_rpe(n2)
 character(len=150) :: EXP_NAME
LOGICAL :: codesQ

INTEGER :: IRHW 

!ML preconditioner
double precision :: Layer1(32 ,177,1), &
                 &  Layer2(32 ,6,1) 
double precision :: NN_r_HP(n1-2,n2)
double precision :: NNa11(n1-2,n2),NNa12(n1-2,n2),NNa21(n1-2,n2)&
                 & ,NNa22(n1-2,n2),NNb11(n1-2,n2),NNb22(n1-2,n2), NN_Y(n2), MLx(n1,n2)
double precision :: inflation_factor0, inflation_factorit

ps(:,:)=0.0d0
divi(:,:)=0.0d0


lowprectime=0.0d0


eps=10.e-11 
itr=100


epa=1.e-30



 DO J=1,n2
   DO I=1,n1
    PMB(I, J)= p(I,J)-b(I,J)
    PMP0(I,J)= p(I,J)-p0(I,J)
   enddo
 enddo

p_true(:,:)=p(:,:)
p0_true(:,:)=p0(:,:)
b_true(:,:)=b(:,:)

DO J=1,n2
  DO I=1,n1
    r(I,J)=rpe_0
    ar(I,J)=rpe_0
  enddo
enddo

do l=1,lord
  DO J=1,n2
    DO I=1,n1
      x(I,J,l)=rpe_0
      ax(I,J,l)=rpe_0
    enddo
  enddo
enddo



call lap0_depth(a11,a12,a21,a22,b11,b22,                   &
     &          pb,p0,e1,e2,hx,hy,cor,n1,n2,gc1,gc2, &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
           & DP_Depth)


! calculate initial residual

call laplfirst_depth(p(:,:),r_HP(:,:), a11,a12,a21,a22,b11,b22,PMP0,&
     &                     pfx,pfy,S,n1,n2,IP, 52,DP_Depth)



call cpu_time(startLP)

 DO J=1,n2
   DO I=1,n1
 
    r_HP(I,J)=rpe_05*r_HP(I,J)-(PMB(I,J))
    r(I,J)  = r_HP(I,J)
   enddo
 enddo

!! calculate initial residual done



!DO J=1,n2
!  DO I=1,n1
!
!      err0=err0+r_HP(I,J)*r_HP(I,J)
!
!  enddo
!enddo


err0 =maxval(ABS(r_HP(:,:)))

!call write_MLfields_iter(r_HP(:,:), KT, n1, n2, EXP_NAME, 'R_iter', 0)
if (iprint==1) then
    call write_residual(r,eps*Exit_cond, 0, TIME, codesQ, codes, IRHW, DX_rpe, DY_rpe,&
                     & n1, n2, 52, 5 ,EXP_NAME)
endif

     
!write(*,*) 0, maxval(ABS(r_HP(:,:))), eps*Exit_cond

! Apply Full Machine-learned preconditioner

      NN_Y(:) = DY_rpe(:)/(3.5d0/2.0d0)
    do J=0,n2-1

      NN_r_HP(:,J+1)= (r_HP(2:n1-1,n2-J))/9.80616
      NNa11(:,J+1) = (a11(2:n1-1,n2-J)-1.17612126e+11 )/(8.54163233e+11 - 2.48437139e+10)
      NNa12(:,J+1) = a12(2:n1-1,n2-J)/(2.0d0*7.32606531e+08)
      NNa21(:,J+1) = (a21(2:n1-1,n2-J)-2.57373289e+10 )/(4.34271611e+10 - 1.28429346e+09) 
      NNa22(:,J+1) = a22(2:n1-1,n2-J)/(2.0d0*7.32606531e+08)
      NNb11(:,J+1) =  b11(2:n1-1,n2-J)/(2.0d0*5.26020440e+09)
      NNb22(:,J+1) = b22(2:n1-1,n2-J)/(2.0d0*2.63625191e+09 )
    enddo


call ML_precon(MLx(:,:),Layer1, Layer2,NN_Y,NN_r_HP,NNa11,&
          & NNa12,NNa21,NNa22,NNb11,NNb22,n1-2,n2)

do J=0,n2-1

  x(:,n2-J,1)=MLx(:,J+1)

enddo
    x(:,:,1)=x(:,:,1)*9.80616
!call precon(r,x(:,:,1),ax(:,:,1), T_step,  A_c, ps, divi,a11,a12,a21,a22,b11,b22,p0,  &
!                &   pfx,pfy,s,n1,n2,ip,ID_PREC, 52, DP_Depth)


! Apply linear operator
call lapl_depth(x(:,:,1),ax(:,:,1), A11,A12,A21,A22,B11,B22,P0,pfx,pfy,S,n1,n2,IP,52,DP_Depth)


  DO J=1,n2
    DO I=1,n1
      ax(I,J,1)=rpe_05*ax(I,J,1)-x(I,J,1)
    enddo
  enddo

! Apply linear operator done

call cpu_time(endLP)
lowprectime=lowprectime + endLP-startLP

! start solver iterations
do it=1,itr

 If (it==100) then
   save_time=.true.
   write(*,*) 'too many iterations'
   exit
 endif
  do l=1,lord
  

    ! calculate beta
    ax2(l)=rpe_0
    rax=rpe_0
    DO J=1,n2
      DO I=1,n1
        rax=rax+r(I,J)*ax(I,J,l)
        ax2(l)=ax2(l)+ax(I,J,l)*ax(I,J,l)
      enddo
    enddo

    ax2(l)=max(epa,ax2(l))
    beta=-rax/ax2(l) 
    errn=0.0d0

  if (iprint==1) then
     write(*,*) 'beta', it, beta
  endif


      ! update Phi and residual
      DO J=1,n2
        DO I=1,n1
         p_T(I,J)=p_T(I,J) +beta* x(I,J,l) 
         r_HP(I,J)  =r_HP(I,J)   +beta*ax(I,J,l) 
         r(I,J)  = r_HP(I,J)
         !errn=errn+r_HP(I,J)*r_HP(I,J)
        enddo
      enddo

  !call write_MLfields_iter(p_true(:,:)+p_T(:,:), KT, n1, n2, EXP_NAME, 'H_iter', it)
  !call write_MLfields_iter(r_HP(:,:), KT, n1, n2, EXP_NAME, 'R_iter', it)


 if (iprint==1) then
    call write_residual(r,eps*Exit_cond, it, TIME, codesQ, codes, IRHW, DX_rpe, DY_rpe, n1, n2, 52, 5 ,EXP_NAME)
endif


!write(*,*) it, maxval(ABS(r_HP(:,:))), eps*Exit_cond
!read(*,*)

!! check whether exit condition is satisfied
if(maxval(ABS(r_HP(:,:))) .lt. eps*err0) exit


! machine-learned preconditioner

    do J=0,n2-1

      NN_r_HP(:,J+1)= (r_HP(2:n1-1,n2-J))/9.80616

    enddo


call ML_precon(MLx(:,:),Layer1, Layer2,NN_Y,NN_r_HP,NNa11,&
          & NNa12,NNa21,NNa22,NNb11,NNb22,n1-2,n2)

    do J=0,n2-1

     qu(:,n2-J)=MLx(:,J+1)
    enddo
    qu(:,:)=qu(:,:)*9.80616

   ! call precon(r,qu, aqu , T_step,  A_c, ps, divi,a11,a12,a21,a22,b11,b22,p0,   &
   !             &   pfx,pfy,s,n1,n2,ip,ID_PREC, 52, DP_Depth)




! Apply linear operator

  call lapl_depth(qu(:,:),aqu(:,:), A11,A12,A21,A22,B11,B22,P0,pfx,pfy,S,n1,n2,IP , 52,DP_Depth)


      DO J=1,n2
        DO I=1,n1
      aqu(I,J)=rpe_05*aqu(I,J)-qu(I,J)
        enddo
      enddo


! calculate alpha
    do ll=1,l
      axaqu(ll)=rpe_0
      DO J=1,n2
        DO I=1,n1
          axaqu(ll)=axaqu(ll)+ax(I,J,ll)*aqu(I,J)
        enddo
      enddo
      del(ll)=-axaqu(ll)/ax2(ll)

      if (iprint==1) then
        write(*,*) 'alpha', it, del(ll)
      endif

    enddo


! new direction that is conjugated to previous
    if(l.lt.lord) then

      DO J=1,n2
        DO I=1,n1
          x(I,J,l+1)= qu(I,J)
          ax(I,J,l+1)=aqu(I,J)
        enddo
      enddo
      do ll=1,l
        DO J=1,n2
          DO I=1,n1
            x(I,J,l+1)= x(I,J,l+1)+del(ll)* x(I,J,ll)
            ax(I,J,l+1)=ax(I,J,l+1)+del(ll)*ax(I,J,ll)
          enddo
        enddo
      enddo

    else

 
      DO J=1,n2
        DO I=1,n1
          x(I,J,1)= qu(I,J)+del(1)* x(I,J,1)
          ax(I,J,1)=aqu(I,J)+del(1)*ax(I,J,1)
        enddo
      enddo

      do ll=2,l
        DO J=1,n2
          DO I=1,n1
            x(I,J,1 )= x(I,J,1)+del(ll)* x(I,J,ll)
            x(I,J,ll)=rpe_0
            ax(I,J,1 )=ax(I,J,1)+del(ll)*ax(I,J,ll)
            ax(I,J,ll)=rpe_0
          enddo
        enddo
      enddo

    endif

 
  enddo

  if(maxval(ABS(r_HP(:,:))) .lt. eps*err0) then !! to replace the go to 200
    
    niter=it
    exit
  else
    niter = itr
  end if

end do


if (iprint==1) then

call lap0_depth(a11,a12,a21,a22,b11,b22,                   &
     &          pb,p0,e1,e2,hx,hy,cor,n1,n2,gc1,gc2, &
           &    MGH1IHX, MGH2IHY, AC, BC, AD, BD,  &
           & 0)



! sanity check to screen
r0_true(:,:)=0.0d0
call laplfirst(p_true(:,:),r0_true(:,:),a11,a12,a21,a22,b11,b22, p0_true,   &
     &                           pfx,pfy,s,n1,n2,ip)

DO J=1,n2
  DO I=1,n1
    r0_true(I,J)=0.5d0*r0_true(I,J)-(p_true(I,J)-b_true(I,J))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo

r_true(:,:)=0.0d0
call laplfirst(p_true(:,:)+p_T(:,:),r_true(:,:),a11,a12,a21,a22,b11,b22, p0_true,   &
     &                           pfx,pfy,s,n1,n2,ip)

DO J=1,n2
  DO I=1,n1
    r_true(I,J)=0.5d0*r_true(I,J)-(p_true(I,J)+p_T(I,J)-b_true(I,J))
 ! write(*,*), i, J, P(i,J)
  enddo
enddo
err_true=0.0d0
err0_true=0.0d0

DO J=1,n2
  DO I=1,n1
    err0_true=err0_true+r0_true(I,J)*r0_true(I,J)
    err_true=err_true+r_true(I,J)*r_true(I,J)
  enddo
enddo
write(*,*) niter
 write(*,*) 'truth ACC',sqrt(err_true/err0_true),'max(rtrue)' ,&
           & maxval(ABS(r_true(:,:))),'max(r0true)', maxval(ABS(r0_true(:,:))), 'max(r)',maxval(ABS(r_HP(:,:))),&
           &'max(r0)',err0 , 'EXIT', eps*Exit_cond 
endif

call cpu_time(finish)

icount=icount+1

nitsm=nitsm+niter
sum_time=sum_time+(finish-start)
sum_lp_time=sum_lp_time+lowprectime
end subroutine


subroutine precon_prep(T_step,A, B, C, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)

use implicit_functions_DP

implicit none
double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
double precision ::  A(N,M), B(N,M), C(N,M)
double precision :: AQU(N,M), T_step, Delta_t, max_QX_QY
INTEGER :: IP(N), ID_PREC
INTEGER :: N, M, I, J

IF (ID_PREC==5) then !! ADI type preconditioner

! 1) get timestep length for preconditioner Delta_t from linear stability argument
!write(*,*) 'in precon ADI'
max_QX_QY=(abs(A11(1,1))+abs(A21(1,1)))/(rpe_2*S(1,1))
DO J=1,M

  DO I=1,N
     max_QX_QY=max(max_QX_QY,(abs(A11(I,J))+abs(A21(I,J)))/(rpe_2*S(I,J)) )     
  ENDDO

ENDDO
T_step=rpe_025/max_QX_QY !0.92d0!
Delta_t=rpe_1/T_step

  DO J=1,M
   DO I=2,N-1
     A(I,J)= -Delta_t*A11(I-1,J)/(rpe_2*S(I,J))
   end do
  end do  

  !DO J=1,M
  !   A(1,J)= A(N-1,J)
  !   A(0,J)= A(N-2,J)
  !end do  

  DO J=1,M
   DO I=2,N-1
     B(I,J)= rpe_1 +Delta_t+ Delta_t* (A11(I+1,J)+A11(I-1,J))/(rpe_2*S(I,J))
   end do
  end do  
 
  !DO J=1,M
  !   B(1,J)= B(N-1,J)
  !   B(0,J)= B(N-2,J)
  !end do 

  DO J=1,M
   DO I=2,N-1
     C(I,J)= -Delta_t*A11(I+1,J)/(rpe_2*S(I,J))
   end do
  end do  

  CALL XBC(A,N,M)
  CALL XBC(B,N,M)
  CALL XBC(C,N,M)

end if



end subroutine

subroutine precon_prep_depth(T_step,A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC, DP_Depth)

use implicit_functions_DP

implicit none
double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
double precision ::  A(N,M), B(N,M), C(N,M),ps(N+1,M), divi(N,M)
double precision :: AQU(N,M), T_step, Delta_t, max_QX_QY
INTEGER :: IP(N), ID_PREC, DP_Depth
INTEGER :: N, M, I, J

IF (ID_PREC==5) then !! ADI type preconditioner

! 1) get timestep length for preconditioner Delta_t from linear stability argument
!write(*,*) 'in precon ADI'
max_QX_QY=(abs(A11(1,1))+abs(A21(1,1)))/(rpe_2*S(1,1))
DO J=1,M

  DO I=1,N
     max_QX_QY=max(max_QX_QY,(abs(A11(I,J))+abs(A21(I,J)))/(rpe_2*S(I,J)) )     
  ENDDO

ENDDO
T_step=rpe_025/max_QX_QY !0.92d0!
Delta_t=rpe_1/T_step



      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N-1
     A(I,J)= -Delta_t*A11(I-1,J)/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N-1
     B(I,J)= rpe_1 +Delta_t+ Delta_t* (A11(I+1,J)+A11(I-1,J))/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N-1
     C(I,J)= -Delta_t*A11(I+1,J)/(rpe_2*S(I,J))
        enddo
      enddo


      DO J=1,DP_Depth
        DO I=1,N-1
     A(I,J)= -Delta_t*A11(I-1,J)/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N-1
     A(I,J)= -Delta_t*A11(I-1,J)/(rpe_2*S(I,J))
        enddo
      enddo


      DO J=1,DP_Depth
        DO I=1,N-1
     B(I,J)= rpe_1 +Delta_t+ Delta_t* (A11(I+1,J)+A11(I-1,J))/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N-1
     B(I,J)= rpe_1 +Delta_t+ Delta_t* (A11(I+1,J)+A11(I-1,J))/(rpe_2*S(I,J))
        enddo
      enddo


      DO J=1,DP_Depth
        DO I=1,N-1
     C(I,J)= -Delta_t*A11(I+1,J)/(rpe_2*S(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N-1
     C(I,J)= -Delta_t*A11(I+1,J)/(rpe_2*S(I,J))
        enddo
      enddo



  CALL XBC(A,N,M)
  CALL XBC(B,N,M)
  CALL XBC(C,N,M)

      do J=1,M
       ps(2,J)=0.0d0
       ps(3,J)=0.0d0 
      enddo

      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
     ps(I+2,J)= -C(I,J)/divi(I,J)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-1
     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
     ps(I+2,J)= -C(I,J)/divi(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
     ps(I+2,J)= -C(I,J)/divi(I,J)
        enddo
      enddo

end if



end subroutine

subroutine ML_precon (QU,Layer1, Layer2,NN_Y,NNr_HP,NN_a11,NN_a12,NN_a21,NN_a22,NN_b11,NN_b22, N, M)
use implicit_functions_DP
implicit none
double precision :: QU(N+2,M)
double precision :: Layer1(32 ,177,1), &
                 &  Layer2(32 ,6,1) 
double precision :: NNr_HP(N,M),NN_a11(N,M),NN_a12(N,M),NN_a21(N,M),NN_a22(N,M),NN_b11(N,M),NN_b22(N,M), NN_Y(M)
INTEGER :: N, M

double precision :: ML_input0(176),ML_input1(5), scaling
INTEGER :: I, J,k, lat_box,lon_box,  pos_lat, pos_lon
INTEGER :: merid_ext, zonal_ext, counter

merid_ext=2
zonal_ext=2

Do J=0,M-1

    ML_input0(1)=NN_Y(J+1)
  do I=0,N-1



   ! get stencil from fields NN_Y,r_HP,NN_a11,NN_a12,NN_a21,NN_a22,NN_b11,NN_b22
   
      scaling=0.0d0
     
      counter=2
      do lat_box = -merid_ext,merid_ext
        do lon_box = -zonal_ext,zonal_ext
           pos_lat= J+lat_box
           pos_lon=Modulo(I+lon_box,N)
           if (pos_lat<0) then
             pos_lat=abs(pos_lat)-1
             pos_lon=Modulo(I+int(N/2)+lon_box,N)
           elseif (pos_lat>=M) then
             pos_lat=pos_lat-(Mod(pos_lat,M)+1)
             pos_lon=Modulo(I+int(N/2)+lon_box,N)
           endif
           if (pos_lon<0) then
             !write(*,*) 'Exe1'
             !write(*,*) pos_lon
             pos_lon=N-pos_lon+1
             !write(*,*) pos_lon
           endif
           !write(*,*) counter, pos_lon, pos_lat
           !read(*,*)
           ML_input0(counter)=  NNr_HP(pos_lon+1, pos_lat+1)
           if(abs(NNr_HP(pos_lon+1, pos_lat+1))> scaling) then
             scaling=abs(NNr_HP(pos_lon+1, pos_lat+1))
           endif
           ML_input0(counter+1)=NN_a11(pos_lon+1, pos_lat+1)
           ML_input0(counter+2)=NN_a12(pos_lon+1, pos_lat+1)
           ML_input0(counter+3)=NN_a21(pos_lon+1, pos_lat+1)
           ML_input0(counter+4)=NN_a22(pos_lon+1, pos_lat+1)
           ML_input0(counter+5)=NN_b11(pos_lon+1, pos_lat+1)
           ML_input0(counter+6)=NN_b22(pos_lon+1, pos_lat+1)

           counter=counter+7  ! continue counting
       end do
     end do
     ML_input0(2:176:7)=ML_input0(2:176:7)/scaling
   !write(*,*) J, I
   !write(*,*) (ML_input0(k), k=1,counter-1)

   ! apply neural network
   !call apply_layer_ReLu(ML_input0,Layer1(J+1,:,:),ML_input1,177,5)

   call apply_layer_Lin(ML_input0,Layer1(J+1,:,:),QU((I+1)+1,J+1),177,1)
   QU((I+1)+1,J+1)=QU((I+1)+1,J+1)*scaling
   !write(*,*) 'ML_Output'
   !write(*,*) QU(I+1+1,J+1)
   !read(*,*)

 end do
end do

  CALL XBC(QU,N+2,M)

end subroutine

subroutine apply_layer_ReLu(Input,Weights,Output,In_dim,Out_dim)
double precision :: Input(176)
double precision :: Weights(177, 5)
double precision :: Output(5)
double precision :: summand
integer :: In_dim, Out_dim
integer :: I, J



do J=1,Out_dim

  summand=0.0d0

  do I=1,In_dim-1
    summand=summand+Input(I)*Weights(I,J)
  end do
  summand=summand+Weights(In_dim,J)   !! add the bias
  
  Output(J)=dmax1(0.0d0,summand)   !! RELU activation function, hardcoded for now
end do


end subroutine

subroutine apply_layer_Lin(Input,Weights,Output,In_dim,Out_dim)
double precision :: Input(176)
double precision :: Weights(177,1)
double precision :: Output(1)
double precision :: summand
integer :: In_dim, Out_dim
integer :: I, J



do J=1,Out_dim

  summand=0.0d0

  do I=1,In_dim-1
    summand=summand+Input(I)*Weights(I,J)
  end do
  summand=summand+Weights(In_dim,J)   !! add the bias
  
  Output(J)=summand   !! linear activation function, hardcoded for now
end do


end subroutine


subroutine precon(R,QU,AQU, T_step,A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC, num_of_bits,DP_Depth)

use implicit_functions_DP

implicit none
double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
double precision ::  A(N,M), B(N,M), C(N,M), ps(N+1,M), divi(N,M)
double precision :: AQU(N,M), T_step
INTEGER :: IP(N), ID_PREC, num_of_bits,DP_Depth
INTEGER :: N, M, I, J

IF (ID_PREC==5) then !! ADI type preconditioner



call precon_ADI(R,QU , T_step, A, ps, divi, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
 ! write(*,*) 'does he get out?'
!! regardless of preconditioners L(q^(n+1)) is needed
 




! orig
!  DO J=1,M
!    DO I=1,N
!      AQU(I,J)=rpe_05*AQU(I,J)-QU(I,J)
!    ENDDO
!  ENDDO
! end orig


  ! write(*,*) 'AQU complete'
elseif(ID_PREC==1) then  !! possible choice of several preconditioners
call precon_LAPL(R,QU ,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)
elseif(ID_PREC==2) then 
call precon_LAPL_MRes(R,QU ,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)
!! regardless of preconditioners L(q^(n+1)) is needed
  call lapl(QU,AQU, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

  DO J=1,M
    DO I=1,N
      AQU(I,J)=rpe_05*AQU(I,J)-QU(I,J)
    ENDDO
  ENDDO
elseif(ID_PREC==3) then

call precon_LAPL_MRes_opt(R,QU,AQU,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)

elseif(ID_PREC==4) then

call precon_LAPL_MRes2_opt(R,QU,AQU,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)

end if



end subroutine

!   implement ADI type preconditioner based on q_n+1 =q_n + dt{ Lz(q_n+1) + Lm(q_n) + H(q_n+1) -R}
!
SUBROUTINE precon_ADI(R,QU , T_step,  A, ps, divi,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC, num_of_bits,DP_Depth)
use implicit_functions_DP

implicit none

double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N), ID_PREC, num_of_bits,DP_Depth
INTEGER :: N, M

double precision ::  max_QX_QY, F(N,M), rhs(N,M), qs(N+1,M,0:4), ps(N+1,M), ws(N+1,M,0:4), A(N,M), B(N,M), C(N,M), divi(N,M)
double precision ::  aa(M,1:4), deti, det40, det41, det42, det43, det44, det3,   &
                           & d11, d12, d13, d14, d21, d22, d23, d24,       &
                           & d31, d32, d33, d34, d41, d42, d43, d44,       &
                           & s1, s2, s3, s4
double precision :: T_step, Delta_t
!double precision :: 
integer :: iter, max_iter, time_scale  !! number of richardson iterations
INTEGER :: I, J, iteration

max_iter=1
!initialize the inverse of R with 0
DO J=1,M
  DO I=1,N
    QU(I,J)=rpe_0
  end do
end do


Delta_t=rpe_1/T_step
!Delta_t=1.0d0/T_step
!write(*,*) Delta_t




! 2) loop of following ADi iterations
do iteration=1,max_iter
 ! 2.1 calculate new right hand side using old QU and R to get tilde{tilde(R)}
If(iteration==1) then

! 1) first iteration

      DO J=1+DP_Depth,M-DP_Depth
        DO I=1,N
      rhs(I,J)=Delta_t*( - R(I,J))
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=1,N
      rhs(I,J)=Delta_t*( - R(I,J))
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=1,N
      rhs(I,J)=Delta_t*( - R(I,J))
        enddo
      enddo


else
  DO J=2,M-1
    DO I=2,N-1
      V(I,J)=QU(I,J+1)-QU(I,J-1)
    end do
  end do
  
  DO I=2,N-1
    V(I,1)=QU(I,2)-QU(IP(I),1)
    V(I,M)=QU(IP(I),M)-QU(I,M-1)
  ENDDO   

  CALL XBC(V,N,M)

  DO J=1,M
    DO I=1,N

      V(I,J)=V(I,J)*A21(I,J) ! +B22(I,J)*(QU(I,J))

    ENDDO
  ENDDO

  CALL XBC(V,N,M)

  DO J=2,M-1
    DO I=2,N-1
      F(I,J)= (V(I,J+1)-V(I,J-1))/(rpe_2*S(I,J))
    end do
  end do    
 

  DO I=2,N-1
    F(I,1)= ((V(I,2)+V(I,1)))/(rpe_2*S(I,1))
    F(I,M)= -((V(I,M)+V(I,M-1)))/(rpe_2*S(I,M))
  ENDDO


  CALL XBC(F,N,M)


  DO J=1,M
    DO I=1,N    
 
      rhs(I,J)=QU(I,J)+Delta_t*(  F(I,J)   &
                          &   - R(I,J))
    END DO
  END DO
 CALL XBC(rhs,N,M)

end if
!DO J=1,M
!
!   rhs(0,J)= rhs(N-2,J)
!
!end do 
  !write(*,*) 'right hand side finished'
 ! 2.2 solve for new QU implicitly, aka solving a tridiagonal system with periodic BC 

  ! calculate and store coefficients ai bi and ci

  !DO J=1,M
  ! DO I=2,N-1
  !   write(*,*) A(I,J)%val, B(I,J)%val, C(I,J)%val
  ! end do
  ! read(*,*)
  !end do 
  
  !DO J=1,M
  !   C(1,J)= C(N-1,J)
  !   C(0,J)= C(N-2,J)
  !end do 
  ! write(*,*) 'A B C coeffficients finished'
  ! calculate p and q's, p is the same for all, 4 q's are needed , first with (bc1, bc2) =(0,0), second with (1,0), third with (0,1)

  ! initialize the 5 linear subsystems with boundary conditions (left boundary)
  DO J=1,M
!    ps(2,J)=rpe_0
!    ps(3,J)=rpe_0

    qs(2,J,0)  =rpe_0
    qs(3,J,0)  =rpe_0

    qs(2,J,1)  =rpe_1
    qs(3,J,1)  =rpe_0

    qs(2,J,2)  =rpe_0
    qs(3,J,2)  =rpe_1

    qs(2,J,3)  =rpe_0
    qs(3,J,3)  =rpe_0

    qs(2,J,4)  =rpe_0
    qs(3,J,4)  =rpe_0
  end do


  !! rewritten to change precision at poles as desired   

!      DO J=1+DP_Depth,M-DP_Depth
!        DO I=2,N-1
!     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
!     ps(I+2,J)= -C(I,J)/divi(I,J)
!        enddo
!      enddo

!      DO J=1,DP_Depth
!        DO I=2,N-1
!     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
!     ps(I+2,J)= -C(I,J)/divi(I,J)
!        enddo
!      enddo
!
!      DO J=M+1-DP_Depth,M
!        DO I=2,N-1
!     divi(I,J)=A(I,J)*ps(I,J)+B(I,J)
!     ps(I+2,J)= -C(I,J)/divi(I,J)
!        enddo
!      enddo

   !! end of rewrite


!  DO J=1,M
!   DO I=2,N-1
!
!     qs(I+2,J,0)= (rhs(I,J)-A(I,J)*qs(I,J,0))/divi(I,J)
!     qs(I+2,J,1)= (rpe_0   -A(I,J)*qs(I,J,1))/divi(I,J)
!     qs(I+2,J,2)= (rpe_0   -A(I,J)*qs(I,J,2))/divi(I,J)
!     qs(I+2,J,3)= (rpe_0   -A(I,J)*qs(I,J,3))/divi(I,J)
!     qs(I+2,J,4)= (rpe_0   -A(I,J)*qs(I,J,4))/divi(I,J)
!
!   end do
!  end do 
  !! rewritten to change precision at poles as desired    
      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
     qs(I+2,J,0)= (rhs(I,J)-A(I,J)*qs(I,J,0))/divi(I,J)
     qs(I+2,J,1)= (rpe_0   -A(I,J)*qs(I,J,1))/divi(I,J)
     qs(I+2,J,2)= (rpe_0   -A(I,J)*qs(I,J,2))/divi(I,J)
     qs(I+2,J,3)= (rpe_0   -A(I,J)*qs(I,J,3))/divi(I,J)
     qs(I+2,J,4)= (rpe_0   -A(I,J)*qs(I,J,4))/divi(I,J)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-1
     qs(I+2,J,0)= (rhs(I,J)-A(I,J)*qs(I,J,0))/divi(I,J)
     qs(I+2,J,1)= (rpe_0   -A(I,J)*qs(I,J,1))/divi(I,J)
     qs(I+2,J,2)= (rpe_0   -A(I,J)*qs(I,J,2))/divi(I,J)
     qs(I+2,J,3)= (rpe_0   -A(I,J)*qs(I,J,3))/divi(I,J)
     qs(I+2,J,4)= (rpe_0   -A(I,J)*qs(I,J,4))/divi(I,J)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
     qs(I+2,J,0)= (rhs(I,J)-A(I,J)*qs(I,J,0))/divi(I,J)
     qs(I+2,J,1)= (rpe_0   -A(I,J)*qs(I,J,1))/divi(I,J)
     qs(I+2,J,2)= (rpe_0   -A(I,J)*qs(I,J,2))/divi(I,J)
     qs(I+2,J,3)= (rpe_0   -A(I,J)*qs(I,J,3))/divi(I,J)
     qs(I+2,J,4)= (rpe_0   -A(I,J)*qs(I,J,4))/divi(I,J)
        enddo
      enddo



  !write(*,*) 'q finished'
  ! calculate the 5 linear subsystems with boundary conditions ( right boundary)
 
  DO J=1,M  
  ws(N  ,J,0)=rpe_0
  ws(N+1,J,0)=rpe_0

  ws(N  ,J,1)=rpe_0
  ws(N+1,J,1)=rpe_0

  ws(N  ,J,2)=rpe_0
  ws(N+1,J,2)=rpe_0

  ws(N  ,J,3)=rpe_1
  ws(N+1,J,3)=rpe_0

  ws(N  ,J,4)=rpe_0
  ws(N+1,J,4)=rpe_1
   end do
  
!  DO J=1,M
!   DO I=N-1,2,-1
!     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)
!     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
!     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
!     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
!     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
!    end do
!  end do 

  !! rewritten to change precision at poles as desired    
      DO J=1+DP_Depth,M-DP_Depth
        DO I=N-1,2,-1
     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)

     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=N-1,2,-1
     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)

     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=N-1,2,-1
     ws(I,J,0) = ps(I+2,J)*ws(I+2,J,0)+qs(I+2,J,0)

     ws(I,J,1) = ps(I+2,J)*ws(I+2,J,1)+qs(I+2,J,1)
     ws(I,J,2) = ps(I+2,J)*ws(I+2,J,2)+qs(I+2,J,2)
     ws(I,J,3) = ps(I+2,J)*ws(I+2,J,3)+qs(I+2,J,3)
     ws(I,J,4) = ps(I+2,J)*ws(I+2,J,4)+qs(I+2,J,4)
        enddo
      enddo

   !! end of rewrite

  ! write(*,*) 'w finished'
  ! solve the subsystems for the coefficients a,b,g,d for each latitude

  DO J=1,M
    if (J<=DP_depth .or. J>=M+1-DP_depth) then

    end if
      
       d21=ws(N-2,J,1)-rpe_1
       d22=ws(N-2,J,2)
       d23=ws(N-2,J,3)
       d24=ws(N-2,J,4)
       s2 =-ws(N-2,J,0)

       d11=ws(N-1,J,1)
       d12=ws(N-1,J,2)-rpe_1
       d13=ws(N-1,J,3)
       d14=ws(N-1,J,4)
       s1 =-ws(N-1,J,0)

       d31=ws(2,J,1)
       d32=ws(2,J,2)
       d33=ws(2,J,3)-rpe_1
       d34=ws(2,J,4)
       s3 =-ws(2,J,0)

       d41=ws(3,J,1)
       d42=ws(3,J,2)
       d43=ws(3,J,3)
       d44=ws(3,J,4)-rpe_1
       s4 =-ws(3,J,0)

      det40=d11*det3(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
        &  -d21*det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
        &  +d31*det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
        &  -d41*det3(d12,d13,d14,d22,d23,d24,d32,d33,d34) 
      deti=rpe_1/det40
      det41=s1 *det3(d22,d23,d24,d32,d33,d34,d42,d43,d44) &
        &  -s2 *det3(d12,d13,d14,d32,d33,d34,d42,d43,d44)  &
        &  +s3 *det3(d12,d13,d14,d22,d23,d24,d42,d43,d44)  &
        &  -s4 *det3(d12,d13,d14,d22,d23,d24,d32,d33,d34)  
      det42=d11*det3( s2,d23,d24, s3,d33,d34, s4,d43,d44) &
        &  -d21*det3( s1,d13,d14, s3,d33,d34, s4,d43,d44)  &
        &  +d31*det3( s1,d13,d14, s2,d23,d24, s4,d43,d44)  &
        &  -d41*det3( s1,d13,d14, s2,d23,d24, s3,d33,d34)  
      det43=d11*det3(d22, s2,d24,d32, s3,d34,d42, s4,d44) &
        &  -d21*det3(d12, s1,d14,d32, s3,d34,d42, s4,d44)  &
        &  +d31*det3(d12, s1,d14,d22, s2,d24,d42, s4,d44)  &
        &  -d41*det3(d12, s1,d14,d22, s2,d24,d32, s3,d34)
      det44=d11*det3(d22,d23, s2,d32,d33, s3,d42,d43, s4) &
        &  -d21*det3(d12,d13, s1,d32,d33, s3,d42,d43, s4)  &
        &  +d31*det3(d12,d13, s1,d22,d23, s2,d42,d43, s4)  &
        &  -d41*det3(d12,d13, s1,d22,d23, s2,d32,d33, s3)
       
      aa(J,4)=det44*deti
      aa(J,3)=det43*deti
      aa(J,2)=det42*deti
      aa(J,1)=det41*deti


  end do 
  !write(*,*) 'aa's finished'
!  DO J=1,M
!   DO I=2,N-1
!      QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
!                    &   +aa(J,2)*ws(I,J,2)  &
!                    &   +aa(J,3)*ws(I,J,3)  &
!                    &   +aa(J,4)*ws(I,J,4)
!  enddo
!   write(*,*) J, (QU(I,J), I=1,N)
!read(*,*)
!  enddo

  !! rewritten to change precision at poles as desired    
      DO J=1+DP_Depth,M-DP_Depth
        DO I=2,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

      DO J=1,DP_Depth
        DO I=2,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

      DO J=M+1-DP_Depth,M
        DO I=2,N-1
          QU(I,J)=ws(I,J,0) +aa(J,1)*ws(I,J,1) &
                    &   +aa(J,2)*ws(I,J,2)  &
                    &   +aa(J,3)*ws(I,J,3)  &
                    &   +aa(J,4)*ws(I,J,4)
        enddo
      enddo

  CALL XBC(QU,N,M)

   !! end of rewrite
   
  !write(*,*) 'QU finished'






end do
  !write(*,*) 'BC QU finished'
END SUBROUTINE

function det3(r11,r12,r13,r21,r22,r23,r31,r32,r33) 


implicit none
double precision :: r11,r12,r13,r21,r22,r23,r31,r32,r33, det3

        det3=    r11*r22*r33+r12*r23*r31+r13*r21*r32    &
            &     -r31*r22*r13-r32*r23*r11-r33*r21*r12

end function det3


!! implementation of Smolarkiewicz, Margolin "Variational methods for elliptic problems in Fluid 
!! Models" equation (42) on page 153
SUBROUTINE precon_LAPL(R,QU ,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)
use implicit_functions_DP

implicit none
double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N), ID_PREC
INTEGER :: N, M

double precision :: AQU(N,M), delta_t, check_PI
double precision :: QU_test(N,M), AQU_test(N,M), Save_QU(N,M), SaveInv_test, Inv_test, base_t, SaveDelta_t
integer :: rich_iter, max_rich, time_scale  !! number of richardson iterations
INTEGER :: I, J

max_rich=2
base_t=1.0**2

DO J=1,M
  DO I=1,N
    QU(I,J)=rpe_0
  end do
end do

!! calculate || r ||_2 as reference; Aim: || L(P^(-1)(r))-r ||_2 should be smaller than || r ||_2
SaveInv_test=rpe_0
DO J=1,M
  DO I=1,N
    SaveInv_test=SaveInv_test+(R(I,J))**2
  ENDDO
ENDDO
SaveInv_test=sqrt(SaveInv_test)
       !! end calculate || r ||_2 as reference

do rich_iter =1,max_rich
  !write(*,*) rich_iter
     !! calculate L(q^(mue))
  call lapl(QU,AQU, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

  DO J=1,M
    DO I=1,N
      AQU(I,J)=rpe_05*AQU(I,J)-QU(I,J)
    ENDDO
  ENDDO
      !! end calculate L(q^(mue))



  do time_scale=1,10  !! initialize richardson iterations with different timestep lengths

    delta_t=base_t/(2**time_scale)    !! scaling timestep lentgth by 2^n

    !! iterate new test q^(mue+1) from q^(mue) by richardson iteration
    DO J=1,M  
      DO I=1,N
        QU_test(I,J)=QU(I,J)+delta_t*(AQU(I,J)-R(I,J))
      ENDDO
    ENDDO

    DO J=1,M   !! zonal boundary conditions of preconditioner
      QU_test(1,J)=QU_test(N-1,J)
      QU_test(N,J)=QU_test(2  ,J)
    ENDDO 
    !! end iterate new test q^(mue+1)
    
    !! Checking: is it a good delta_t: check if || L(q_test^(mue)) -r ||_2
    call lapl(QU_test,AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

    DO J=1,M
      DO I=1,N
        AQU_test(I,J)=rpe_05*AQU_test(I,J)-QU_test(I,J)
      ENDDO
    ENDDO

     !! calculate || L(q_test^(mue)) -r ||_2
    Inv_test=rpe_0
    DO J=1,M
      DO I=1,N
        Inv_test=Inv_test+(AQU_test(I,J)-R(I,J))**2
      ENDDO
    ENDDO
    !write(*,*) 'convergence of inverse', time_scale, sqrt(Inv_test), SaveInv_test
    If(sqrt(Inv_test)<=SaveInv_test) then
      SaveInv_test= sqrt(Inv_test) !! new minimal norm 
      Save_QU = QU_test       !! new best q
    end if
    !! Checking end

  end do

    DO J=1,M
      DO I=1,N
        QU(I,J)=Save_QU(I,J)
      end do
    end do

end do

END SUBROUTINE


SUBROUTINE precon_LAPL_MRes(R,QU ,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)
use implicit_functions_DP

implicit none
double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N), ID_PREC
INTEGER :: N, M

double precision :: AQU(N,M), delta_t, Ares2, resAres
double precision :: QU_test(N,M), AQU_test(N,M), residuum(N,M) 
integer :: rich_iter, max_rich !! number of richardson iterations
INTEGER :: I, J
double precision:: epa
epa=1.e-30

max_rich=1


DO J=1,M
  DO I=1,N
    QU(I,J)=rpe_0
  end do
end do


do rich_iter =1,max_rich
  !write(*,*) rich_iter
     !! calculate L(q^(mue))
  call lapl(QU,AQU, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

  DO J=1,M
    DO I=1,N
      AQU(I,J)=rpe_05*AQU(I,J)-QU(I,J)
    ENDDO
  ENDDO
      !! end calculate L(q^(mue))



  DO J=1,M
    DO I=1,N
        residuum(I,J)=AQU(I,J)-R(I,J)
    ENDDO
  ENDDO


  call lapl(residuum,AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

  DO J=1,M
    DO I=1,N
      AQU_test(I,J)=rpe_05*AQU_test(I,J)-residuum(I,J)
    ENDDO
  ENDDO



    Ares2=rpe_0
    resAres=rpe_0
    DO J=1,M
      DO I=1,N
        resAres=resAres+residuum(I,J)*AQU_test(I,J)
        Ares2=Ares2+AQU_test(I,J)*AQU_test(I,J)
      enddo
    enddo

    Ares2=max(epa,Ares2)
    delta_t=-resAres/Ares2

    write(*,*) 'min_res', delta_t

    DO J=1,M  
      DO I=1,N
        QU(I,J)=QU(I,J)+delta_t*(residuum(I,J))
      ENDDO
    ENDDO

    DO J=1,M   !! zonal boundary conditions of preconditioner
      QU(1,J)=QU(N-1,J)
      QU(N,J)=QU(2  ,J)
    ENDDO 
    !! end iterate new test q^(mue+1)


end do

END SUBROUTINE

SUBROUTINE precon_LAPL_MRes_opt(R,QU,AQU,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)
use implicit_functions_DP

implicit none
double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N), ID_PREC
INTEGER :: N, M

double precision :: AQU(N,M), delta_t, Ares2, resAres
double precision :: QU_test(N,M), AQU_test(N,M), residuum(N,M)
integer :: rich_iter, max_rich !! number of richardson iterations
INTEGER :: I, J
double precision:: epa
epa=1.e-30

max_rich=1


QU(:,:)=rpe_0
AQU(:,:)=rpe_0  ! zero by default


  !write(*,*) rich_iter


!  DO J=1,M
!    DO I=1,N
!        residuum(I,J)=AQU(I,J)-R(I,J)
!    ENDDO
!  ENDDO


!  call lapl(residuum,AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)
call lapl(-R(:,:),AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

!call LAPL_tilde(-R(:,:),AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)
  DO J=1,M
    DO I=1,N
!      AQU_test(I,J)=rpe_05*AQU_test(I,J)-residuum(I,J)
       AQU_test(I,J)=rpe_05*AQU_test(I,J)+R(I,J)
    ENDDO
  ENDDO



    Ares2=rpe_0
    resAres=rpe_0
    DO J=1,M
      DO I=1,N
 !       resAres=resAres+residuum(I,J)*AQU_test(I,J)
        resAres=resAres-R(I,J)*AQU_test(I,J)
        Ares2=Ares2+AQU_test(I,J)*AQU_test(I,J)
      enddo
    enddo

    Ares2=max(epa,Ares2)
    delta_t=-resAres/Ares2

    DO J=1,M  
      DO I=1,N
       ! QU(I,J)=QU(I,J)+delta_t*(residuum(I,J))
         QU(I,J)=delta_t*(-R(I,J))
      ENDDO
    ENDDO

    DO J=1,M   !! zonal boundary conditions of preconditioner
      QU(1,J)=QU(N-1,J)
      QU(N,J)=QU(2  ,J)
    ENDDO 
    !! end iterate new test q^(mue+1)
    DO J=1,M  
      DO I=1,N
       ! QU(I,J)=QU(I,J)+delta_t*(residuum(I,J))
         AQU(I,J)=delta_t*(AQU_test(I,J))
      ENDDO
    ENDDO


END SUBROUTINE

SUBROUTINE precon_LAPL_MRes2_opt(R,QU,AQU,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP,ID_PREC)
use implicit_functions_DP

implicit none
double precision :: R(N,M),QU(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N), ID_PREC
INTEGER :: N, M

double precision :: AQU(N,M), delta_t, Ares2, resAres
double precision :: QU_test(N,M), AQU_test(N,M), residuum(N,M)
integer :: rich_iter, max_rich !! number of richardson iterations
INTEGER :: I, J
double precision:: epa
epa=1.e-30

max_rich=2


QU(:,:)=rpe_0
AQU(:,:)=rpe_0  ! zero by default


  !write(*,*) rich_iter


!  DO J=1,M
!    DO I=1,N
!        residuum(I,J)=AQU(I,J)-R(I,J)
!    ENDDO
!  ENDDO


!  call lapl(residuum,AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)
  call lapl(-R(:,:),AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

  DO J=1,M
    DO I=1,N
!      AQU_test(I,J)=rpe_05*AQU_test(I,J)-residuum(I,J)
       AQU_test(I,J)=rpe_05*AQU_test(I,J)+R(I,J)
    ENDDO
  ENDDO



    Ares2=rpe_0
    resAres=rpe_0
    DO J=1,M
      DO I=1,N
 !       resAres=resAres+residuum(I,J)*AQU_test(I,J)
        resAres=resAres-R(I,J)*AQU_test(I,J)
        Ares2=Ares2+AQU_test(I,J)*AQU_test(I,J)
      enddo
    enddo

    Ares2=max(epa,Ares2)
    delta_t=-resAres/Ares2

    DO J=1,M  
      DO I=1,N
       ! QU(I,J)=QU(I,J)+delta_t*(residuum(I,J))
         QU(I,J)=delta_t*(-R(I,J))
      ENDDO
    ENDDO

    DO J=1,M   !! zonal boundary conditions of preconditioner
      QU(1,J)=QU(N-1,J)
      QU(N,J)=QU(2  ,J)
    ENDDO 
    !! end iterate new test q^(mue+1)
    DO J=1,M  
      DO I=1,N
       ! QU(I,J)=QU(I,J)+delta_t*(residuum(I,J))
         AQU(I,J)=delta_t*(AQU_test(I,J))
      ENDDO
    ENDDO

do rich_iter =2,max_rich
  !write(*,*) rich_iter
     !! calculate L(q^(mue))
!  call lapl(QU,AQU, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

!  DO J=1,M
!    DO I=1,N
!      AQU(I,J)=rpe_05*AQU(I,J)-QU(I,J)
!    ENDDO
!  ENDDO
      !! end calculate L(q^(mue))



  DO J=1,M
    DO I=1,N
        residuum(I,J)=AQU(I,J)-R(I,J)
    ENDDO
  ENDDO


  call lapl(residuum,AQU_test, A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)

  DO J=1,M
    DO I=1,N
      AQU_test(I,J)=rpe_05*AQU_test(I,J)-residuum(I,J)
    ENDDO
  ENDDO



    Ares2=rpe_0
    resAres=rpe_0
    DO J=1,M
      DO I=1,N
        resAres=resAres+residuum(I,J)*AQU_test(I,J)
        Ares2=Ares2+AQU_test(I,J)*AQU_test(I,J)
      enddo
    enddo

    Ares2=max(epa,Ares2)
    delta_t=-resAres/Ares2

    DO J=1,M  
      DO I=1,N
        QU(I,J)=QU(I,J)+delta_t*(residuum(I,J))
      ENDDO
    ENDDO

    DO J=1,M   !! zonal boundary conditions of preconditioner
      QU(1,J)=QU(N-1,J)
      QU(N,J)=QU(2  ,J)
    ENDDO 
    !! end iterate new test q^(mue+1)
    DO J=1,M  
      DO I=1,N
        AQU(I,J)=AQU(I,J)+delta_t*(AQU_test(I,J))
      ENDDO
    ENDDO

end do

END SUBROUTINE



SUBROUTINE LAPL_tilde(P,F,A11,A12,A21,A22,B11,B22,P0,U,V,S,N,M,IP)
use implicit_functions_DP

implicit none
double precision :: P(N,M),F(N,M),A11(N,M),A12(N,M),A21(N,M),A22(N,M),  &
    &      B11(N,M),B22(N,M),P0(N,M),U(N,M),V(N,M),S(N,M)
INTEGER :: IP(N)
INTEGER :: N, M

double precision :: UTIL, VTIL
INTEGER :: I, J


!DO J=2,M-1
!  DO I=2,N-1
!    U(I,J)=P(I+1,J)-P(I-1,J)
!    V(I,J)=P(I,J+1)-P(I,J-1)
!  end do
!end do
  
!DO I=2,N-1
!  U(I,1)=P(I+1,1)-P(I-1,1)
!  U(I,M)=P(I+1,M)-P(I-1,M)
!  V(I,1)=P(I,2)-P(IP(I),1)
!  V(I,M)=P(IP(I),M)-P(I,M-1)
!ENDDO   

!CALL XBC(U,N,M)
!CALL XBC(V,N,M)

DO J=1,M
  DO I=1,N
!    UTIL=U(I,J)*A11(I,J)+A12(I,J)*V(I,J)+B11(I,J)*(P(I,J))  ! why only B11*P and not B11*(P-P0)
!    VTIL=V(I,J)*A21(I,J)+A22(I,J)*U(I,J)+B22(I,J)*(P(I,J))
!    U(I,J)=UTIL
!    V(I,J)=VTIL
     UTIL=U(I,J)*A11(I,J)+B11(I,J)*(P(I,J))
     VTIL=V(I,J)*A21(I,J)+B22(I,J)*(P(I,J))
     U(I,J)=UTIL
     V(I,J)=VTIL
  ENDDO
ENDDO

CALL XBC(U,N,M)
CALL XBC(V,N,M)

DO J=2,M-1
  DO I=2,N-1
    F(I,J)= (U(I+1,J)-U(I-1,J)+V(I,J+1)-V(I,J-1))/S(I,J)
  end do
end do    
 
DO I=2,N-1
  F(I,1)= (U(I+1,1)-U(I-1,1)+(V(I,2)+V(I,1)))/S(I,1)
  F(I,M)= (U(I+1,M)-U(I-1,M)-(V(I,M)+V(I,M-1)))/S(I,M)
ENDDO

CALL XBC(F,N,M)


END SUBROUTINE

subroutine ETOPO5(h0,x,y,n,m, mountain)
use implicit_functions_DP

implicit none

INTEGER :: n, m, counter
INTEGER*2 :: TOPO5(2160,4320)
integer :: trans_x(n),trans_y(m)
double precision :: h0(n,m)
double precision :: x(n),y(m), x_top(4320),y_top(2160)
double precision :: hs0, Rad, x_c, y_c, dist, pi, sigma, step, distance
integer :: i, j
logical :: mountain

pi = acos(-1.0d0)



do j=1,m
  do i=1,n
    h0(i,j)=0.0d0
  end do
end do

If (mountain) then
  open (13, file = '../Model/ETOPO5.DOS',form='unformatted', access = 'direct', recl = 4320*2)

 do i=1,2160
  read (13, rec=i) (TOPO5(i,j), j=1,4320 )
 enddo
  close(unit=13)
! write(*,*) 'value', (TOPO5(1,j), j=1,4320 )

step=(359.0d0*60.0d0+55.0d0)/(360.0d0*60.0d0)
step=step/4319
do j=1,4320
x_top(j)=(j-1)*step*2*pi
enddo

step=(179.0d0*60.0d0+55.0d0)/(180.0d0*60.0d0)
step=step/2159
do i=1,2160
y_top(i)=pi/2.0d0-(i-1)*step*pi
end do

do i=1,n
 trans_x(i)=1
 distance=sqrt( (x_top(1)-x(i))**2)
 do j=2,4320
  if (distance > sqrt( (x_top(j)-x(i))**2)) then
   trans_x(i)=j
   distance=sqrt( (x_top(j)-x(i))**2)
  endif
 enddo
enddo

do i=1,m
 trans_y(i)=1
 distance=sqrt( (y_top(1)-y(i))**2)
 do j=2,2160
  if (distance > sqrt( (y_top(j)-y(i))**2)) then
   trans_y(i)=j
   distance=sqrt( (y_top(j)-y(i))**2)
  endif
 enddo
enddo


  do j=1,m
    do i=n/2,n-1

      h0(i,j)=max(0.0d0, 0.5d0*real( TOPO5(trans_y(j),trans_x(i+1-n/2)) ) )

    end do

  end do
  do j=1,m
    do i=1,n/2-1
      h0(i,j)=max(0.0d0, 0.5d0*real( TOPO5(trans_y(j),trans_x(i+n/2)) ) )

   end do

  end do

CALL XBC(h0,n,m)
!  do j=1,m
!    do i=1,n
!      dist= sqrt( (x(i)-x_c)**2 + (y(j)-y_c)**2) 
!      h0(i,j)=0.4d0*hs0*(1.0d0/(sigma*sqrt(2*pi))*exp(-0.5d0*(dist/sigma)**2) )
!    end do
!  end do

end if

end subroutine

! subroutine initc(xlon, ylat, rho)
! 
! implicit none
! 
! double precision :: xlon, ylat, rho
! 
! double precision :: pi, anot, rnot, xcent, ycent, gamm, dist
! !
! pi = acos(-1.0)
! anot = 1.0
! 
! !      rnot = 2*pi*7./128.
! !      xcent = 3*pi/2.
! 
! rnot = 2*pi*10.0/128.
! xcent = 2*pi/2.
! ycent = pi/2. *1.0
! !     
! gamm=1.0
! rho =rpe_0
! !
! 
! dist = 2.*sqrt( (cos(ylat)*sin( (xlon-xcent)/2 ) )**2    &
!   &   + sin((ylat-ycent)/2)**2       )
! 
! if (dist.le.rnot) rho = anot                             &
!   &   *(1.0-gamm*dist/rnot)                               &
!   &   + rho
! !
! 
! end subroutine



! subroutine dep (xar, yar, xde, yde, rvel, beta, dt)
! 
! implicit none
! 
! double precision :: xar, yar, xde, yde, rvel, beta, dt
! 
! 
! double precision :: cosalf, sinalf, zz, yy, xx, yyp, xxp, zzp, ypm, xpm, pi, xpofs
! !     
! !     distance in radians moved per time step
!       xpofs = rvel * dt
! !     
! !     pole of velocity rotated by angle beta
! !     
!       cosalf = cos(beta)
!       sinalf = sin(beta)
! !     
!       zz = sin(yar)
!       yy = sin(xar)*cos(yar)
!       xx = cos(xar)*cos(yar)
! !     
!       yyp = yy
!       xxp = xx*cosalf + zz*sinalf
!       zzp = zz*cosalf - xx*sinalf
! !     
!       ypm = asin(zzp)
!       xpm = atan2(yyp,xxp)
! !     
!       xpm = xpm + xpofs
! !     
!       zzp = sin(ypm)
!       yyp = sin(xpm)*cos(ypm)
!       xxp = cos(xpm)*cos(ypm)
! !     
!       yy = yyp
!       xx = xxp*cosalf - zzp*sinalf
!       zz = zzp*cosalf + xxp*sinalf
! !     
!       yde = asin(zz)
!       xde = atan2(yy,xx)
! !     
!       pi = 4*atan(1.0)
! !     
!       if(yde.lt.-pi/2) then
!          yde = -pi-yde
!          xde = xde + pi/2.
!       write(6,*) '-pi/2',xde,yde,xar,yar
!       endif
!       if(yde.gt.pi/2)  then
!          yde = pi-yde
!          xde = xde + pi/2.
!       write(6,*) '+pi/2',xde,yde,xar,yar
!       endif
!       if(xde.lt.0.0) then
!          xde = xde + 2*pi
!       endif
!       if(xde.gt.2*pi) then
!          xde = xde - 2*pi
!       endif
! !
! end subroutine   

! 
!       SUBROUTINE PLOT(H,H0,U,V,HX,HY,N1,N2,IRHW)
!       PARAMETER(N=130,M=64,NM=N*M*2)
!       DIMENSION H(N1,N2),H0(N1,N2),U(N1,N2),V(N1,N2)
!      1         ,HX(N1,N2),HY(N1,N2),WORK(NM)
!       COMMON// UP(N,M),VP(N,M),F1(N,M)
!       CHARACTER*80 LHEAD
!       COMMON /SMOLAB/ ISWT,ILAB, IOFFM
!       DATA ISWT,ILAB, IOFFM /1,1,1/
!       COMMON /RECINT/ IRECMJ     ,IRECMN     ,IRECTX
!       COMMON/VPLT/ IFLV,IVU1,IVU2,IVRT
!       COMMON /VEC2/ BIG,IX0,IX1,INCX,IY0,IY1,INCY
!       COMMON /STR03/  INITA , INITB , AROWL , ITERP , ITERC , IGFLG     
!      1             ,  IMSG , UVMSG , ICYC , DISPL , DISPC , CSTOP       
!       COMMON/DISPL/ TIME,U00,DX,DY,H00,HMTS,DT
!       COMMON/TOPOG/ TP(6000),ITC,ZCT
!       DATA TP/6000*0./
!       REAL LM0,LM1
!       DATA IVCTPL,ISTRPL,IUPL,IVPL,ISPL/0,0,1,0,0/
!       DATA VECMN,VECMX,LENV/1.E-3,0.,0/
!       DATA IFLV,IVRT/1,0/
!       DATA ICOLOR/0/
! C
!       ITC=0
!       ZCT=1. 
! C
!       iswt=1
!       ilab=0
!       ioffm=1
!       incx=4
!       incy=4
!       inita=8
!       initb=8
! C
!       I1=INT(102.4+409.6)
!       IPT1=INT(192.8+819.2)-150
!       TIMED=TIME/24
! C
! C
! CONTOUR PLOTS FOLLOW
! C SURFACE HEIGHT
!       IF(ICOLOR.EQ.0) THEN
!       CALL SET(.1,.9,.3,.7,0.,2.,-0.5,0.5,1)
!       CALL LABMOD('(F4.1)','(F4.1)',4,4,2,2,20,20,0)
!       CALL PERIML(4,5,2,5)
!       X2=FLOAT(N1)
!       Y2=FLOAT(N2)
!       CALL SET(.1,.9,.3,.7,1.,X2,1.,Y2,1)
! c     ENCODE(45,50,LHEAD) TIMED
!       WRITE (UNIT=LHEAD, FMT=50) TIMED
!    50 FORMAT('GEOPOTENTIAL PERTURBATION AT T=',F8.2,'  DAYS')
!       CALL WTSTR(CPUX(522),CPUY(IPT1),LHEAD(1:45),2,0,0)
!       CALL WTSTR(CPUX(I1),CPUY(200),'x/Pi',3,0,0)
!       CALL WTSTR(CPUX(22),CPUY(I1),'y/Pi',3,90,0)
!       CNT=0.0025
!       IF(IRHW.EQ.1) CNT=0.025
!       DO 214 J=1,N2
!       DO 214 I=1,N1
!       F1(I,J)=(H(I,J)+0.*78.4e3-H00)/H00
!       IF(H(I,J)/9.81-H0(I,J).LE.0.1) F1(I,J)=0.
!   214 CONTINUE
!       CALL CONREC(F1,N1,N1,N2,CNT,1.5,CNT,1,-1,-1252)
!       DO 215 J=1,N2
!       DO 215 I=1,N1
!   215 F1(I,J)=-F1(I,J)
!       CALL CONREC(F1,N1,N1,N2,CNT,1.5,CNT,1,-1,1252)
!       SPEED=-1.E15
!       DO 216 J=1,N2
!       DO 216 I=1,N1
!       VP(I,J)=V(I,J)
!       UP(I,J)=U(I,J)
!       IF(H(I,J)/9.81-H0(I,J).LE.0.1) THEN
!       UP(I,J)=0.
!       VP(I,J)=0.
!       ENDIF
!       SPEED=AMAX1(SPEED, UP(I,J)**2+VP(I,J)**2)
!   216 CONTINUE
!       SPEED=SQRT(SPEED)
!       PRINT 2784, SPEED
!  2784 FORMAT(5X,'MAX SPEED=',E11.4,' M/S')
!       write(6,*) 'veklvct'
! c     IF(IVCTPL.EQ.1) CALL VELVCT(UP,N1,VP,N1,N1,N2,VECMN,VECMX,1,LENV)
!       IF(ISTRPL.EQ.1) CALL STRMLN(UP,VP,WORK,N1,N1,N2,1,IER)
!       DO 217 J=1,N2
!       DO 217 I=1,N1
!   217 F1(I,J)=H0(I,J)/HMTS
!       CALL SETUSV('LW',2000) 
!       write(6,*) 'dritter'
!       CALL CONREC(F1,N1,N1,N2,0.25,1.,-3.25,1,-1,-1252)
!       CALL SETUSV('LW',1000) 
!       CALL FRAME
!   
! C U VELOCITY
!       IF(IUPL.EQ.1) THEN
!       CALL SET(.1,.9,.3,.7,0.,2.,-0.5,0.5,1)
!       CALL LABMOD('(F4.1)','(F4.1)',4,4,2,2,20,20,0)
!       CALL PERIML(4,5,2,5)
!       X2=FLOAT(N1)
!       Y2=FLOAT(N2)
!       CALL SET(.1,.9,.3,.7,1.,X2,1.,Y2,1)
! c     ENCODE(30,51,LHEAD) TIME
!       WRITE (UNIT=LHEAD, FMT=51) TIME
!    51 FORMAT('U VELOCITY AT T=',F8.2,' HOURS')
!       CALL WTSTR(CPUX(512),CPUY(IPT1),LHEAD(1:30),2,0,0)
!       CALL WTSTR(CPUX(I1),CPUY(200),'x/Pi',3,0,0)
!       CALL WTSTR(CPUX(22),CPUY(I1),'y/Pi',3,90,0)
!       CNT=.5
!       IF(IRHW.EQ.1) CNT=20.
!       DO 314 J=1,N2
!       DO 314 I=1,N1
!   314 F1(I,J)=UP(I,J)
!       CALL CONREC(F1,N1,N1,N2,CNT,10.,CNT,1,-1,-1252)
!       DO 315 J=1,N2
!       DO 315 I=1,N1
!   315 F1(I,J)=-F1(I,J)
!       CALL CONREC(F1,N1,N1,N2,CNT,10.,CNT,1,-1,1252)
!       CALL FRAME
!       ENDIF
!   
! C V VELOCITY
!       IF(IVPL.EQ.1) THEN
!       CALL SET(.1,.9,.3,.7,0.,2.,-0.5,0.5,1)
!       CALL LABMOD('(F4.1)','(F4.1)',4,4,2,2,20,20,0)
!       CALL PERIML(4,5,2,5)
!       X2=FLOAT(N1)
!       Y2=FLOAT(N2)
!       CALL SET(.1,.9,.3,.7,1.,X2,1.,Y2,1)
! c     ENCODE(30,52,LHEAD) TIME
!       WRITE (UNIT=LHEAD, FMT=52) TIME
!    52 FORMAT('V VELOCITY AT T=',F8.2,' HOURS')
!       CALL WTSTR(CPUX(512),CPUY(IPT1),LHEAD(1:30),2,0,0)
!       CALL WTSTR(CPUX(I1),CPUY(200),'x/Pi',3,0,0)
!       CALL WTSTR(CPUX(22),CPUY(I1),'y/Pi',3,90,0)
!       CNT=.5
!       IF(IRHW.EQ.1) CNT=20.
!       DO 414 J=1,N2
!       DO 414 I=1,N1
!   414 F1(I,J)=VP(I,J)
!       CALL CONREC(F1,N1,N1,N2,CNT,10.,CNT,1,-1,-1252)
!       DO 415 J=1,N2
!       DO 415 I=1,N1
!   415 F1(I,J)=-F1(I,J)
!       CALL CONREC(F1,N1,N1,N2,CNT,10.,CNT,1,-1,1252)
!       CALL FRAME
!       ENDIF
! 
! C ISOTACHS
!       IF(ISPL.EQ.1) THEN
!       CALL SET(.1,.9,.3,.7,0.,2.,-0.5,0.5,1)
!       CALL LABMOD('(F4.1)','(F4.1)',4,4,2,2,20,20,0)
!       CALL PERIML(4,5,2,5)
!       X2=FLOAT(N1)
!       Y2=FLOAT(N2)
!       CALL SET(.1,.9,.3,.7,1.,X2,1.,Y2,1)
! c     ENCODE(28,53,LHEAD) TIME
!       WRITE (UNIT=LHEAD, FMT=53) TIME
!    53 FORMAT('ISOTACHS AT T=',F8.2,' HOURS')
!       CALL WTSTR(CPUX(512),CPUY(IPT1),LHEAD(1:30),2,0,0)
!       CALL WTSTR(CPUX(I1),CPUY(200),'x/Pi',3,0,0)
!       CALL WTSTR(CPUX(22),CPUY(I1),'y/Pi',3,90,0)
!       CNT=.5
!       IF(IRHW.EQ.1) CNT=20.
!       DO 514 J=1,N2
!       DO 514 I=1,N1
!   514 F1(I,J)=SQRT(VP(I,J)**2+UP(I,J)**2)
!       CALL CONREC(F1,N1,N1,N2,CNT,10.,CNT,1,-1,-1252)
!       CALL FRAME
!       ENDIF
!   
!       ELSE
! 
!       DO 814 J=1,N2
!       DO 814 I=1,N1
!   814 F1(I,J)= (H(I,J)-H00)/H00
!       X2=FLOAT(N1)
!       Y2=FLOAT(N2)
!       CALL SFLUSH
!       CALL GSTXCI(1)
!       CALL GSPLCI(14)
!       CALL COLORPL(F1,N1,N2)
!       CALL GETSET (XVPL,XVPR,YVPB,YVPT,XWDL,XWDR,YWDB,YWDT,LNLG)
!       CALL SET(XVPL,XVPR,YVPB,YVPT,0.,2.,-0.5,0.5,LNLG)
!       CALL LABMOD('(F4.1)','(F4.1)',4,4,2,2,20,20,0)
!       CALL PERIML(4,5,2,5)
!       CALL SET(XVPL,XVPR,YVPB,YVPT,1.,X2,1.,Y2,1)
! c     ENCODE(45,850,LHEAD) TIME
!       WRITE (UNIT=LHEAD, FMT=850) TIME
!   850 FORMAT('GEOPOTENTIAL PERTURBATION AT T=',F8.2,' HOURS')
!       CALL WTSTR(CPUX(512),CPUY(IPT1),LHEAD(1:45),2,0,0)
!       CALL WTSTR(CPUX(I1),CPUY(200),'x/Pi',3,0,0)
!       CALL WTSTR(CPUX(22),CPUY(I1),'y/Pi',3,90,0)
! 
!       SPEED=-1.E15
!       DO 816 J=1,N2
!       DO 816 I=1,N1
!       VP(I,J)=V(I,J)
!       UP(I,J)=U(I,J)
!       SPEED=AMAX1(SPEED, UP(I,J)**2+VP(I,J)**2)
!   816 CONTINUE
!       SPEED=SQRT(SPEED)
!       PRINT 2784, SPEED
!       CALL GSTXCI(1)
!       CALL GSPLCI(12)
!       IF(IVCTPL.EQ.1) CALL VELVCT(UP,N1,VP,N1,N1,N2,VECMN,VECMX,1,LENV)
!       IF(ISTRPL.EQ.1) CALL STRMLN(UP,VP,WORK,N1,N1,N2,1,IER)
!       DO 817 J=1,N2
!       DO 817 I=1,N1
!   817 F1(I,J)=H0(I,J)/HMTS
!       IRECMJ=13
!       IRECMN=13
!       CALL CONREC(F1,N1,N1,N2,0.25,1.,0.25,1,-1,-1252)
!       IRECMJ=1
!       IRECMN=1
!       CALL FRAME
!       ENDIF
! 
!       RETURN
!       END
