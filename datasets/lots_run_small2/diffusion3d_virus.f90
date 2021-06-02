

          PROGRAM diffusion3d_virus
	  IMPLICIT NONE
          REAL PIE, TOLER
          PARAMETER(PIE=3.141592654, TOLER=1.E-10) 
          integer, parameter :: sp = kind(1e0), dp = kind(1d0)
	  INTEGER NX,NY,NZ,NG,NONODS_NG, NLAY_IN,NLAY_OUT
!          parameter(NX=12,NY=12,NZ=3,NG=8,NONODS_NG=(NX-2)*(NY-2)*(NZ-2)*NG,NLAY_IN=NONODS_NG,NLAY_OUT=NONODS_NG) 
!          parameter(NX=42,NY=42,NZ=3,NG=8,NONODS_NG=(NX-2)*(NY-2)*(NZ-2)*NG,NLAY_IN=NONODS_NG,NLAY_OUT=NONODS_NG) 
          parameter(NX=12,NY=12,NZ=3,NG=8,NONODS_NG=(NX-2)*(NY-2)*(NZ-2)*NG,NLAY_IN=NONODS_NG,NLAY_OUT=NONODS_NG) 
! NITS_SOLV=no of non-linear iterations used for anything other than THETA=0.0 NITS_SOLV>1.
          REAL T_NEW(NX,NY,NZ,NG), &
	       T_OLD(NX,NY,NZ,NG), RHS(NX,NY,NZ,NG), A_MAT_T(NX,NY,NZ,NG)   
          REAL T_NEW_LONG(NONODS_NG), T_DECOD_LONG(NONODS_NG),T_EXACT_LONG(NONODS_NG)
          REAL N1(NX,NY,NZ),N2(NX,NY,NZ)
          REAL DT,V_N(NG),DX,DY,DZ,RELAX_SOLV
          REAL SIGMAT(NX,NY,NZ,NG),KK(NX,NY,NZ,NG),S(NX,NY,NZ,NG), SIGMA_F(NX,NY,NZ,NG,NG)
          REAL SIGMA_S_OFF(NX,NY,NZ,NG,NG)
! SIGMAT is diagonal sigma_t, sigma_f is total fission matrix, sigma_s_off is just offdiagonal scattering. 
          REAL A_MAT(NX,NY,NZ,NG,-1:1, -1:1, -1:1), B_MAT(NX,NY,NZ,NG,-1:1, -1:1,-1:1),&
               FIS_MAT(NX,NY,NZ,NG,-1:1, -1:1, -1:1)
          REAL V(NX,NY,NZ,NG) ! Temp vector used for normalizing fission as well as calculating fission source
          INTEGER MATERIAL(NX,NY,NZ) ! Material number
! ANN:
       REAL, ALLOCATABLE :: WEIGHT(:),NEURAL(:,:),IDEAL_INPUT(:,:),IDEAL_OUTPUT(:,:),IDEAL_MIDLE(:,:)
       REAL, ALLOCATABLE :: CT(:,:),ALPHA_NEW(:),ALPHA(:),NEURAL_ONE(:)
       REAL, ALLOCATABLE :: INTEGRAL_T_NEW(:,:), SAMPLE_T_NEW(:,:) 
       INTEGER :: ISAM_PT, JSAM_PT, KSAM_PT
       REAL :: W_EXP(1000)
       INTEGER :: NLAYER, NLAY(1000)
       INTEGER :: NCOLM,NEXAMPLE
! Local...
       INTEGER :: ILAY,IEXAMPLE,I,j,k,NOD,NITS_SOLV,NITS_SOLV_NG,&
		  NTIME,NTIME_ONE, NOREXP,NBATCH,NALPHA,IALPHA,ITEST,NTEST,ITS_ROM,NITS_ROM,NO_NUR
       INTEGER :: NCOLM2, NLAYER2, NLAY_IN2,NLAY_OUT2,NO_NUR2, NDEAL_MIDLE
       INTEGER :: II,JJ, ISTART,IFINI, JSTART,JFINI, IWIDTH,JWIDTH
       LOGICAL :: GOT_WEIGHT, FOUR_REGS
       REAL :: ERROR_TOLER,ERROR_SOLV_NG,NEU_ALPHA
       REAL :: SIG_RAN_NO,K_RAN_NO, DIAG_REG_TOLER(7)
       INTEGER :: SCALE_IN,SCALE_OUT, IRAN_NO, ISTART_NEU, NITS_EIG,ITS_EIG, NITS, NITS_GLOBAL, IG,JG
       INTEGER :: ITIME, ITS
       REAL :: RAN_NO(100000)
       logical :: RECALL, EIGEN_START, pickup_WEIGHT
       REAL :: MAX_SCALE_IN, MIN_SCALE_IN, MAX_SCALE_OUT, MIN_SCALE_OUT, KEFF, ERROR_SOLV,&
               ERROR_EIG,KEFF_EXACT
       real :: keff2, lambda, mean_keff_error, mean_t_error
       REAL :: ONE_DAY, ACCTIM, TIME_MAT
       REAL :: FORCE, H2M
       REAL :: DURATION_INFECT, INCUBATION_DURATION, VIRUS_SIGMA, VIRUS_GAMMA
       REAL :: XI, R0_HOME, R0_MOBILE, BETA, BETA2, MU, NU
       real :: VIRUS_N1_AIM, VIRUS_N2_AIM, length, LAMBDA_H_H, LAMBDA_M_M, R_EIGEN, RDAY
       real :: VIRUS_N1, VIRUS_N2
       real :: r_ratio, rswitch
       REAL :: SIGN_SMALL ! function
       LOGICAL :: EIGEN_VALUE, EIGEN_METHOD1, FOUND_FILE
       INTEGER :: NSAVE_TO_FILE,ISAVE_TO_FILE,IFILE_COUNT
       CHARACTER CHARA*240
       CHARACTER str_trim*240
! Local variables...
! Set up problem ********************************
         print *,'--hello'
         EIGEN_VALUE=.false.
         EIGEN_METHOD1=.true. ! (new method=9.1437710515488124) (old method1=9.1435851016370382)
         NSAVE_TO_FILE=4 ! Save the 1st file and then every N'the file unless =0 and then dont save.
         ISAVE_TO_FILE=0 ! for counting the number of time steps before outputing to disc
         IFILE_COUNT=1 ! Starting file number to generate. 
! original method1=0.99937582280063098; new method=0.99937582280063209
! ratio of indoor/outdoor = 1462/57 = 25.65
         DT=1.e+3 ! time in seconds.
         IF(EIGEN_VALUE) DT=1.e+10
         V_N=1.0 ! velocity
         RELAX_SOLV=1.0 ! RELAX_SOLV=1 (NOD RELAXATION), =0.5 TAKE 0.5 OF PREVIOUS VALUE TO CONVERGE BETTER IF PROBLEMS CONVERGING. 
         ONE_DAY=24.*3600.0 ! A DAY IN seconds. 

         DX = 100.0e+3/REAL(NX-2) ! domain 100km
         DY = 100.0e+3/REAL(NY-2)
         DZ = 1.0/REAL(NZ)
         NTIME_ONE=1

         if(EIGEN_VALUE) THEN
            NTIME=1
            NITS_SOLV=400
            NITS_SOLV_NG=100
            NITS_EIG=100
            EIGEN_START=.TRUE.
            NITS_GLOBAL=1
            ERROR_SOLV=1.0E-10
            ERROR_SOLV_NG=1.0E-10
         ELSE ! Time dep...
!            NTIME=1000
!            NTIME=2000
!            NTIME=23000
!            NTIME= 0.75 * one_day/dt
!            NTIME= 45.0 * one_day/dt
            NTIME= 45.5 * one_day/dt
!            NTIME= 9.0 * one_day/dt
!            NTIME= 10.25 * one_day/dt
!            ntime=1
!            NTIME=20
!            NITS_SOLV=40
!            NITS_SOLV_NG=10
            NITS_SOLV=400
            NITS_SOLV_NG=4
            ERROR_SOLV=1.0E-6!4
            ERROR_SOLV_NG=1.0E-4!4
            NITS_EIG=0
            EIGEN_START=.FALSE.
!            NITS_GLOBAL=10 ! no of non-linear iterations per time step.
            NITS_GLOBAL=2!4 ! no of non-linear iterations per time step.
         ENDIF
         ALLOCATE( INTEGRAL_T_NEW(NG,NTIME),SAMPLE_T_NEW(NG,NTIME)  ) 
         S(:,:,:,:) = 0.0 
         ERROR_EIG=1.0E-3

! Define the materials...
         MATERIAL=1
         IWIDTH=(NX-2)/5
         JWIDTH=(NY-2)/5

! Bottom middle...
         ISTART=2 + 2*IWIDTH; IFINI=1 + 3*IWIDTH
         JSTART=2 + 0*JWIDTH; JFINI=1 + 1*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=2
         ISAM_PT=ISTART; JSAM_PT=JSTART; KSAM_PT=2

! left ...
         ISTART=2 + 0*IWIDTH; IFINI=1 + 1*IWIDTH
         JSTART=2 + 2*JWIDTH; JFINI=1 + 3*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=3

! centre...
         ISTART=2 + 2*IWIDTH; IFINI=1 + 3*IWIDTH
         JSTART=2 + 2*JWIDTH; JFINI=1 + 3*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=4

! right...
         ISTART=2 + 4*IWIDTH; IFINI=1 + 5*IWIDTH
         JSTART=2 + 2*JWIDTH; JFINI=1 + 3*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=5

! Top middle...
         ISTART=2 + 2*IWIDTH; IFINI=1 + 3*IWIDTH
         JSTART=2 + 4*JWIDTH; JFINI=1 + 5*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=6

! either side of middle regions*********************
! Just above bottom middle...
         ISTART=2 + 2*IWIDTH; IFINI=1 + 3*IWIDTH
         JSTART=2 + 1*JWIDTH; JFINI=1 + 2*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=7

! next to left region...
         ISTART=2 + 1*IWIDTH; IFINI=1 + 2*IWIDTH
         JSTART=2 + 2*JWIDTH; JFINI=1 + 3*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=8

! next to right region...
         ISTART=2 + 3*IWIDTH; IFINI=1 + 4*IWIDTH
         JSTART=2 + 2*JWIDTH; JFINI=1 + 3*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=9

! next to top region...
         ISTART=2 + 2*IWIDTH; IFINI=1 + 3*IWIDTH
         JSTART=2 + 3*JWIDTH; JFINI=1 + 4*JWIDTH
         MATERIAL(ISTART:IFINI,JSTART:JFINI,:)=10

! Key variables...
         DURATION_INFECT=7.0*one_day ! 7 day infection duration. 
         INCUBATION_DURATION=4.5*one_day! 4.5 day inCUBATION duration.
         VIRUS_SIGMA=1.0/INCUBATION_DURATION
         VIRUS_GAMMA=1.0/DURATION_INFECT
! ξ (xi) is the rate which recovered individuals return to the susceptible statue due to loss of immunity.
         XI= 1./(365.0*one_day) 
         INQUIRE(FILE='r0-2values.csv', EXIST=FOUND_FILE) 
         IF(FOUND_FILE) THEN
            open(27, file='r0-2values.csv')
!            WRITE(27, *) '# R0_HOME, R0_MOBILE'
            READ(27, *) 
            READ(27, * ) R0_HOME, R0_MOBILE
         ELSE
            R0_HOME=0.2 !5.0 !0.2 ! No of people an infected person infects.
            R0_MOBILE=20.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 
         ENDIF
         print *,'FOUND_FILE, R0_HOME, R0_MOBILE:',FOUND_FILE, R0_HOME, R0_MOBILE
         MU=1.0*1./(60.*365.0*one_day) ! BIRTH RATES - HAVE THE SAME RATE as death - life time 60 years. 
         NU=1./(60.*365.0*one_day) ! DEATH RATES HAVE THE SAME RATE
         T_OLD=0.0
         T_NEW=0.0

! time stepping...
         ACCTIM=0.0
         DO ITIME=1,NTIME
         DO ITS=1,NITS_GLOBAL
! Material at which the material is evaluated. 
            IF(EIGEN_VALUE) THEN
               TIME_MAT=ONE_DAY/4.0
               R_EIGEN=1.0
            ELSE
               TIME_MAT=ACCTIM
               R_EIGEN=0.0
            ENDIF

            RDAY= 0.5*SIN( (TIME_MAT/ONE_DAY) * PIE * 2.0 ) + 0.5

          print *,'itime,ntime,TIME_MAT/max(toler,one_day):',itime,ntime,TIME_MAT/max(toler,one_day)
          WRITE(*, '(200E11.3)' ) (maxval(t_new(:,:,:,ig)),ig=1,ng)
!          print *,'ig,(maxval(t_new(:,:,:,ig)),ig=1,ng):',ig,(maxval(t_new(:,:,:,ig)),ig=1,ng)
          do ig=1,-ng
             print *,'ig,minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig)):',ig,minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig))
          end do
      if(.false.) then
            print *,'TIME_MAT,ONE_DAY, TIME_MAT/ONE_DAY, (TIME_MAT/ONE_DAY) * PIE * 2.0:', &
                     TIME_MAT,ONE_DAY, TIME_MAT/ONE_DAY, (TIME_MAT/ONE_DAY) * PIE * 2.0
            print *,'SIN( (TIME_MAT/ONE_DAY) * PIE * 2.0 ):',SIN( (TIME_MAT/ONE_DAY) * PIE * 2.0 )
      endif
!            LAMBDA_H_H= 0.0
!            beta=0.0
!            beta2=0.0
!            LAMBDA_M_M= 1.0/one_day
            
            KK=-1.0E+10
            SIGMAT=0.0
            SIGMA_S_OFF=0.0
            length=real(nx-2)*dx
            SIGMA_F=0.0 ! This is always true for time depdent problems
            DO K=2,NZ-1
            DO J=2,NY-1
            DO I=2,NX-1

               IF(MATERIAL(I,J,K)==2) THEN ! bottom middle region
                  VIRUS_N1_AIM=1000.0 * (1.0-RDAY) + 1000.0
                  VIRUS_N2_AIM=0.0
                  LAMBDA_H_H= 1000.0*1.0/one_day ! Push some people out of their homes on time scale of 0.1day. 
               ELSE
                  VIRUS_N1_AIM=0.0
                  VIRUS_N2_AIM=0.0
                  LAMBDA_H_H= 0.0/one_day ! Push some people out of their homes on time scale of 0.1day. 
                  if(EIGEN_VALUE) LAMBDA_H_H= 1000.0*1.0/one_day
               ENDIF


               IF((ITIME==1).AND.(.NOT.EIGEN_VALUE)) THEN
                  T_OLD(I,J,K,1)=VIRUS_N1_AIM*0.999 * 2./1.5
                  T_OLD(I,J,K,2)=VIRUS_N1_AIM*0.001 * 2./1.5! start off with 0.1% of people exposed and at home.
                  IF(ITS==1) T_NEW(I,J,K,:)=T_OLD(I,J,K,:) 
               ENDIF

               VIRUS_N1 = MAX(TOLER,SUM( T_NEW(I,J,K,1:4) ) )
               VIRUS_N2 = MAX(TOLER,SUM( T_NEW(I,J,K,5:8) ) )

!               print *,'i,j,k,MATERIAL(I,J,K):',i,j,k,MATERIAL(I,J,K),T_NEW(I,J,K,:)

               IF(MATERIAL(I,J,K)>=2) THEN ! CROSS REGION
!                  KK(I,J,K,:)=0.001*2.5*length**2/one_day
                  KK(I,J,K,:)=2.5*length**2/one_day
!                  KK(I,J,K,:)= max(KK(I,J,K,:) ,  (1.-rday**2)*10.0*length**2/one_day )
!                  IF(MATERIAL(I,J,K)==5) KK(I,J,K,:)=0.001*length**2/one_day
               ELSE
                  KK(I,J,K,:)=-1.0E+10
               ENDIF
               KK(I,J,K,1:4)=0.0 ! dont have diffusion in houses. 
               BETA=VIRUS_GAMMA*R0_HOME
               BETA2=VIRUS_GAMMA*R0_MOBILE
!               if(MATERIAL(I,J,K)==-5) BETA2=VIRUS_GAMMA*R0_MOBILE*1.0 ! area on right has good social distancing. 
               if(MATERIAL(I,J,K)==1) then
                  BETA=0.0 ! area on right has good social distancing. 
                  BETA2=0.0 ! area on right has good social distancing.
                  if(EIGEN_VALUE) LAMBDA_H_H= 100.0*1.0/one_day 
               endif
               LAMBDA_M_M=0.0
               if(EIGEN_VALUE) then
                  if(MATERIAL(I,J,K)==1) LAMBDA_M_M= 100.0*1.0/one_day 
               endif
!               if(MATERIAL(I,J,K)==5) BETA2=VIRUS_GAMMA*R0_HOME ! area on right has good social distancing.
!               if(MATERIAL(I,J,K)==9) BETA2=VIRUS_GAMMA*R0_HOME ! area on right has good social distancing. 
!               K(I,J,K,:)=0.0
!               K(I,J,K,:)=0.1*RDAY*1.0/(one_day*length*
!               IF(MATERIAL(I,J,K)==2) THEN
                  FORCE=(VIRUS_N1-VIRUS_N1_AIM)/MAX( TOLER, VIRUS_N1, VIRUS_N1_AIM ) 
!                  FORCE2=(VIRUS_N2-VIRUS_N2_AIM(I,J,K))/MAX(TOLER,VIRUS_N2)
                  H2M=0.5+0.5*SIGN(1.0, FORCE) 

                  SIGMAT(I,J,K,1)=-MU          + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M *0.01
                  SIGMAT(I,J,K,2)=+VIRUS_SIGMA + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMAT(I,J,K,3)=+VIRUS_GAMMA + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMAT(I,J,K,4)=+XI          + NU + (1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M *0.01

                  SIGMAT(I,J,K,5)=-MU          + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMAT(I,J,K,6)=+VIRUS_SIGMA + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMAT(I,J,K,7)=+VIRUS_GAMMA + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMAT(I,J,K,8)=+XI          + NU - (1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)

! Time dep...
                  SIGMA_S_OFF(I,J,K,1,5)= SIGMA_S_OFF(I,J,K,1,5) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMA_S_OFF(I,J,K,2,6)= SIGMA_S_OFF(I,J,K,2,6) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMA_S_OFF(I,J,K,3,7)= SIGMA_S_OFF(I,J,K,3,7) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)
                  SIGMA_S_OFF(I,J,K,4,8)= SIGMA_S_OFF(I,J,K,4,8) +(1.-R_EIGEN)*LAMBDA_H_H * FORCE * (1.0-H2M)

                  SIGMA_S_OFF(I,J,K,5,1)= SIGMA_S_OFF(I,J,K,5,1) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMA_S_OFF(I,J,K,6,2)= SIGMA_S_OFF(I,J,K,6,2) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMA_S_OFF(I,J,K,7,3)= SIGMA_S_OFF(I,J,K,7,3) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01
                  SIGMA_S_OFF(I,J,K,8,4)= SIGMA_S_OFF(I,J,K,8,4) -(1.-R_EIGEN)*LAMBDA_H_H * FORCE * H2M*0.01

                  SIGMA_S_OFF(I,J,K,1,3)=SIGMA_S_OFF(I,J,K,1,3) + (1.-R_EIGEN)*BETA *T_NEW(I,J,K,1)/MAX(VIRUS_N1,TOLER) 
!       if((i==6).and.(j==2)) then
!         print *,'i,j,k,(1.-R_EIGEN)*BETA *T_NEW(I,J,K,1)/MAX(VIRUS_N1,TOLER):', &
!                  i,j,k,(1.-R_EIGEN)*BETA *T_NEW(I,J,K,1)/MAX(VIRUS_N1,TOLER)
!         print *,'(1.-R_EIGEN),BETA ,T_NEW(I,J,K,1),MAX(VIRUS_N1,TOLER):',(1.-R_EIGEN),BETA ,T_NEW(I,J,K,1),MAX(VIRUS_N1,TOLER)
!         print *,'T_NEW(I,J,K,1:4):',T_NEW(I,J,K,1:4)
!         print *,'SIGMA_S_OFF(I,J,K,1,3):',SIGMA_S_OFF(I,J,K,1,3)
!       endif
                  SIGMA_S_OFF(I,J,K,2,3)=SIGMA_S_OFF(I,J,K,2,3) - (1.-R_EIGEN)*BETA *T_NEW(I,J,K,1)/MAX(VIRUS_N1,TOLER) 
                  SIGMA_S_OFF(I,J,K,5,7)=SIGMA_S_OFF(I,J,K,5,7) + (1.-R_EIGEN)*BETA2*T_NEW(I,J,K,5)/MAX(VIRUS_N2,TOLER)
                  SIGMA_S_OFF(I,J,K,6,7)=SIGMA_S_OFF(I,J,K,6,7) - (1.-R_EIGEN)*BETA2*T_NEW(I,J,K,5)/MAX(VIRUS_N2,TOLER) 
                  SIGMA_S_OFF(I,J,K,1,4) = SIGMA_S_OFF(I,J,K,1,4) - XI
                  SIGMA_S_OFF(I,J,K,5,8) = SIGMA_S_OFF(I,J,K,5,8) - XI

! exposed people resulting in infections...
                  SIGMA_S_OFF(I,J,K,3,2)=SIGMA_S_OFF(I,J,K,3,2) -VIRUS_SIGMA *(1.-R_EIGEN)
                  SIGMA_S_OFF(I,J,K,7,6)=SIGMA_S_OFF(I,J,K,7,6) -VIRUS_SIGMA *(1.-R_EIGEN)
!          print *,'virus_sigma:',virus_sigma
!         stop 393
! infected people that are recovered...
                  SIGMA_S_OFF(I,J,K,4,3)=SIGMA_S_OFF(I,J,K,4,3) -VIRUS_GAMMA 
                  SIGMA_S_OFF(I,J,K,8,7)=SIGMA_S_OFF(I,J,K,8,7) -VIRUS_GAMMA 
!                  SIGMA_S_OFF(I,J,K,4,3)=SIGMA_S_OFF(I,J,K,4,3) -VIRUS_GAMMA *(1.-R_EIGEN)
!                  SIGMA_S_OFF(I,J,K,8,7)=SIGMA_S_OFF(I,J,K,8,7) -VIRUS_GAMMA *(1.-R_EIGEN)

!                  S(I,J,K,1,5) = + (1.-R_EIGEN) * LAMBDA_H_H * VIRUS_N1_AIM(I,J,K,1) 
!                  S(I,J,K,2,6) = + (1.-R_EIGEN) * LAMBDA_H_H * VIRUS_N1_AIM(I,J,K,2)
!                  S(I,J,K,3,7) = + (1.-R_EIGEN) * LAMBDA_H_H * VIRUS_N1_AIM(I,J,K,3)
!                  S(I,J,K,4,8) = + (1.-R_EIGEN) * LAMBDA_H_H * VIRUS_N1_AIM(I,J,K,4)

!                  S(I,J,K,5,1) = + (1.-R_EIGEN) * LAMBDA_M_M * VIRUS_N2_AIM(I,J,K,5)
!                  S(I,J,K,6,2) = + (1.-R_EIGEN) * LAMBDA_M_M * VIRUS_N2_AIM(I,J,K,6)
!                  S(I,J,K,7,3) = + (1.-R_EIGEN) * LAMBDA_M_M * VIRUS_N2_AIM(I,J,K,7)
!                  S(I,J,K,8,4) = + (1.-R_EIGEN) * LAMBDA_M_M * VIRUS_N2_AIM(I,J,K,8)
                  ! Only for eigen-value problems...
!                  SIGMA_F(I,J,K,1,3)=+R_EIGEN*BETA
!                  SIGMA_F(I,J,K,5,7)=+R_EIGEN*BETA2
!                  SIGMA_F(I,J,K,2,3)=-R_EIGEN*BETA
!                  SIGMA_F(I,J,K,6,7)=-R_EIGEN*BETA2
!                  SIGMA_S_OFF(I,J,K,5,1)= SIGMA_S_OFF(I,J,K,5,1) -R_EIGEN*LAMBDA_H_H
!                  SIGMA_S_OFF(I,J,K,6,2)= SIGMA_S_OFF(I,J,K,6,2) -R_EIGEN*LAMBDA_H_H
!                  SIGMA_S_OFF(I,J,K,7,3)= SIGMA_S_OFF(I,J,K,7,3) -R_EIGEN*LAMBDA_H_H
!                  SIGMA_S_OFF(I,J,K,8,4)= SIGMA_S_OFF(I,J,K,8,4) -R_EIGEN*LAMBDA_H_H
          if(EIGEN_VALUE) then
! eigen-value...
               if(.false.) then ! original method...
                  KK(I,J,K,1:8)=0.00001 * KK(I,J,K,1:8)
                  SIGMAT(I,J,K,1)= SIGMAT(I,J,K,1) + 10000000000.0/one_day
                  SIGMAT(I,J,K,5)= SIGMAT(I,J,K,5) + 10000000000.0/one_day

                  SIGMAT(I,J,K,1)= SIGMAT(I,J,K,1) + R_EIGEN*LAMBDA_H_H
                  SIGMAT(I,J,K,2)= SIGMAT(I,J,K,2) + R_EIGEN*LAMBDA_H_H
                  SIGMAT(I,J,K,3)= SIGMAT(I,J,K,3) + R_EIGEN*LAMBDA_H_H
                  SIGMAT(I,J,K,4)= SIGMAT(I,J,K,4) + R_EIGEN*LAMBDA_H_H

                  SIGMAT(I,J,K,5)= SIGMAT(I,J,K,5) + R_EIGEN*LAMBDA_M_M
                  SIGMAT(I,J,K,6)= SIGMAT(I,J,K,6) + R_EIGEN*LAMBDA_M_M
                  SIGMAT(I,J,K,7)= SIGMAT(I,J,K,7) + R_EIGEN*LAMBDA_M_M
                  SIGMAT(I,J,K,8)= SIGMAT(I,J,K,8) + R_EIGEN*LAMBDA_M_M
               else ! consistent with the houses...
                  KK(I,J,K,1:8)=1.0*0.05 * KK(I,J,K,1:8)
                  r_ratio=25.65
                  rswitch=0.0
                  if(material(i,j,k)==2) rswitch=1.0 ! switch on where there are houses
!                  if(material(i,j,k)==3) rswitch=1.0 ! switch on where there are houses
                  LAMBDA_H_H=1.0*rswitch*1.0/one_day
                  LAMBDA_M_M=1.0*(1.0-rswitch)*10000.0/one_day

                  SIGMAT(I,J,K,1)= SIGMAT(I,J,K,1) + 10000000000000.0/one_day
                  SIGMAT(I,J,K,5)= SIGMAT(I,J,K,5) + 10000000000000.0/one_day
                  SIGMAT(I,J,K,4)= SIGMAT(I,J,K,4) + 10000000000000.0/one_day
                  SIGMAT(I,J,K,8)= SIGMAT(I,J,K,8) + 10000000000000.0/one_day

                  SIGMAT(I,J,K,1)= SIGMAT(I,J,K,1) + LAMBDA_H_H    + LAMBDA_M_M 
                  SIGMAT(I,J,K,2)= SIGMAT(I,J,K,2) + LAMBDA_H_H    + LAMBDA_M_M 
                  SIGMAT(I,J,K,3)= SIGMAT(I,J,K,3) + LAMBDA_H_H    + LAMBDA_M_M 
                  SIGMAT(I,J,K,4)= SIGMAT(I,J,K,4) + LAMBDA_H_H    + LAMBDA_M_M 

                  SIGMAT(I,J,K,5)= SIGMAT(I,J,K,5) + LAMBDA_H_H * r_ratio
                  SIGMAT(I,J,K,6)= SIGMAT(I,J,K,6) + LAMBDA_H_H * r_ratio
                  SIGMAT(I,J,K,7)= SIGMAT(I,J,K,7) + LAMBDA_H_H * r_ratio
                  SIGMAT(I,J,K,8)= SIGMAT(I,J,K,8) + LAMBDA_H_H * r_ratio

                  SIGMA_S_OFF(I,J,K,1,5)= SIGMA_S_OFF(I,J,K,1,5) - LAMBDA_H_H * r_ratio
                  SIGMA_S_OFF(I,J,K,2,6)= SIGMA_S_OFF(I,J,K,2,6) - LAMBDA_H_H * r_ratio
                  SIGMA_S_OFF(I,J,K,3,7)= SIGMA_S_OFF(I,J,K,3,7) - LAMBDA_H_H * r_ratio
                  SIGMA_S_OFF(I,J,K,4,8)= SIGMA_S_OFF(I,J,K,4,8) - LAMBDA_H_H * r_ratio

                  SIGMA_S_OFF(I,J,K,5,1)= SIGMA_S_OFF(I,J,K,5,1) - LAMBDA_H_H 
                  SIGMA_S_OFF(I,J,K,6,2)= SIGMA_S_OFF(I,J,K,6,2) - LAMBDA_H_H 
                  SIGMA_S_OFF(I,J,K,7,3)= SIGMA_S_OFF(I,J,K,7,3) - LAMBDA_H_H 
                  SIGMA_S_OFF(I,J,K,8,4)= SIGMA_S_OFF(I,J,K,8,4) - LAMBDA_H_H 
                  
               endif
! 
               if(EIGEN_METHOD1) then ! method that put the 1/keff over the I - infection eqn

                  SIGMA_F(I,J,K,3,2)=+VIRUS_SIGMA *R_EIGEN
                  SIGMA_F(I,J,K,7,6)=+VIRUS_SIGMA *R_EIGEN

                  SIGMA_S_OFF(I,J,K,1,3)=SIGMA_S_OFF(I,J,K,1,3)+R_EIGEN*BETA
                  SIGMA_S_OFF(I,J,K,5,7)=SIGMA_S_OFF(I,J,K,5,7)+R_EIGEN*BETA2
                  SIGMA_S_OFF(I,J,K,2,3)=SIGMA_S_OFF(I,J,K,2,3)-R_EIGEN*BETA
                  SIGMA_S_OFF(I,J,K,6,7)=SIGMA_S_OFF(I,J,K,6,7)-R_EIGEN*BETA2

!                  SIGMA_F(I,J,K,4,3)=+VIRUS_GAMMA *R_EIGEN
!                  SIGMA_F(I,J,K,8,7)=+VIRUS_GAMMA *R_EIGEN

!                  SIGMAT(I,J,K,8)= SIGMAT(I,J,K,8) + R_EIGEN*LAMBDA_H_H
!                  SIGMA_F(I,J,K,1,3)=-R_EIGEN*BETA
!                  SIGMA_F(I,J,K,5,7)=-R_EIGEN*BETA2
!                  SIGMA_F(I,J,K,2,3)=+R_EIGEN*BETA
!                  SIGMA_F(I,J,K,6,7)=+R_EIGEN*BETA2
               else
                  SIGMA_S_OFF(I,J,K,3,2)=SIGMA_S_OFF(I,J,K,3,2)-VIRUS_SIGMA *R_EIGEN
                  SIGMA_S_OFF(I,J,K,7,6)=SIGMA_S_OFF(I,J,K,7,6)-VIRUS_SIGMA *R_EIGEN

                  SIGMA_S_OFF(I,J,K,1,3)=SIGMA_S_OFF(I,J,K,1,3)+R_EIGEN*BETA
                  SIGMA_S_OFF(I,J,K,5,7)=SIGMA_S_OFF(I,J,K,5,7)+R_EIGEN*BETA2
!                  SIGMA_S_OFF(I,J,K,2,3)=SIGMA_S_OFF(I,J,K,2,3)-R_EIGEN*BETA
!                  SIGMA_S_OFF(I,J,K,6,7)=SIGMA_S_OFF(I,J,K,6,7)-R_EIGEN*BETA2
!                  SIGMA_F(I,J,K,1,3)=-R_EIGEN*BETA
!                  SIGMA_F(I,J,K,5,7)=-R_EIGEN*BETA2
                  SIGMA_F(I,J,K,2,3)=+R_EIGEN*BETA 
                  SIGMA_F(I,J,K,6,7)=+R_EIGEN*BETA2 
               endif
           endif ! if(EIGEN_VALUE) then
!                  K(I,J,K,5:8)=RDAY*1.0/(one_day*length**2)
!               ELSE IF(MATERIAL(I,J,K)==2) THEN
!                  SIGMAT=
!                  SIGMA_S_OFF=0.0
!                  SIGMA_F=1.0
!                  K=1.0
!               ENDIF
            END DO
            END DO
            END DO
!               stop 292

! change the sign as we have +ve on the rhs: 
            SIGMA_S_OFF=-SIGMA_S_OFF
        if(.false.) then
            print *,'minval(SIGMAT),maxval(SIGMAT):',minval(SIGMAT),maxval(SIGMAT)
            print *,'minval(SIGMA_S_OFF),maxval(SIGMA_S_OFF):',minval(SIGMA_S_OFF),maxval(SIGMA_S_OFF)
            print *,'minval(KK),maxval(KK):',minval(KK),maxval(KK)
          do ig=1,-ng
          do jg=1,ng
             print *,'ig,jg:',ig,jg
            print *,'minval(SIGMA_S_OFF(:,:,:,ig,jg)),maxval(SIGMA_S_OFF(:,:,:,ig,jg)):', &
                     minval(SIGMA_S_OFF(:,:,:,ig,jg)),maxval(SIGMA_S_OFF(:,:,:,ig,jg))
          end do
          end do
            print *,'minval(SIGMA_F),maxval(SIGMA_F):',minval(SIGMA_F),maxval(SIGMA_F)
        endif ! if(.false.) then

            CALL GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( A_MAT, B_MAT,  &
                         DX,DY,DZ,DT,V_N,SIGMAT,KK, NX,NY,NZ,NG) 

!      print *,'itime,its=',itime,its
!      print *,'minval(t_old),maxval(t_old):',minval(t_old),maxval(t_old)
!      print *,'minval(t_new),maxval(t_new):',minval(t_new),maxval(t_new)
            CALL EIG_OR_TIME_STEPING_DIFFUSION_CALC(T_NEW,KEFF, T_OLD,&
                 NTIME_ONE,NITS_SOLV,NITS_SOLV_NG,NITS_EIG,&
                 ERROR_SOLV,ERROR_SOLV_NG,ERROR_EIG, EIGEN_START, RELAX_SOLV,&
                 DX,DY,DZ,DT,V_N,SIGMAT,SIGMA_S_OFF,KK,SIGMA_F,S, NX,NY,NZ,NG )
!      print *,' '
!      print *,'minval(t_old),maxval(t_old):',minval(t_old),maxval(t_old)
!      print *,' '
!      print *,'minval(t_new),maxval(t_new):',minval(t_new),maxval(t_new)
!      print *,' '


         END DO ! DO ITS=1,NITS_GLOBAL
            ACCTIM=ACCTIM+DT
            DO IG=1,NG
               INTEGRAL_T_NEW(IG,ITIME) = SUM( T_NEW(2:NX-1, 2:NY-1, 2:NZ-1, IG) ) 
               SAMPLE_T_NEW(IG,ITIME) = T_NEW( ISAM_PT, JSAM_PT, KSAM_PT, IG)
            END DO
            T_OLD=T_NEW ! Prepare for next time step.

! ************* save to disc***************
         ISAVE_TO_FILE=ISAVE_TO_FILE+1

      IF((ISAVE_TO_FILE==NSAVE_TO_FILE) .or. ((NSAVE_TO_FILE>0).and. (ITIME==1))) THEN
         ISAVE_TO_FILE=0
!         IF(ITIME==1) THEN
!            open(27, file='group-output-time.csv', STATUS='new')
!            open(27, file='group-output-time.csv', status='replace')
! trim(chara4(ino))
            str_trim=chara(IFILE_COUNT); k=INDEX(str_trim,' ')-1
!            open(27, file='run/group-output-time'//str_trim(1:k)//'.csv', status='replace')
            open(27, file='group-output-time'//str_trim(1:k)//'.csv', status='replace')
!            open(27, file='group-output-time'//chara(IFILE_COUNT)//'.csv', status='replace')
            WRITE(27, *) '# DX,DY,DZ, NX,NY,NZ:', DX,DY,DZ, NX-2,NY-2,NZ-2
            WRITE(27, *) '# there are 8 groups: HOME-S, HOME-E, HOME-I, HOME-R, MOBILE-S, MOBILE-E, MOBILE-I, MOBILE-R'
            
            WRITE(27, *) '# DO K=1,NZ'
            WRITE(27, *) '# do j=1,NY,1,-1'
            WRITE(27, *) '#    WRITE(27, "(200E14.4)" ) (T_NEW( I, J, K, IG ) ,I=1,NX) ! in single procession for printing' 
            WRITE(27, *) '# end do'
            WRITE(27, *) '# end do'
!         ELSE
!            open(27, file='group-output-time.csv', Access = 'append')
!         ENDIF
         WRITE(27, *) '# TIME=', ACCTIM
         DO IG=1,NG
            WRITE(27, *) 'GROUP=',IG
            DO K=2,NZ-1
            do j=2,ny-1
               WRITE(27, '(200E14.4)' ) (SIGN_SMALL(T_NEW( I, J, K, IG )) ,I=2,NX-1) ! in single procession for printing
            end do
            end do
!           stop 382
         END DO
         close(27) 
!         CALL SYSTEM('gzip group-output.csv' )
         str_trim=chara(IFILE_COUNT); k=INDEX(str_trim,' ')-1
!         CALL SYSTEM('gzip run/group-output-time'//str_trim(1:k)//'.csv' ) ! factor of 12 reduction in disc space.
         CALL SYSTEM('rm group-output-time'//str_trim(1:k)//'.csv.gz' ) ! factor of 12 reduction in disc space.
         CALL SYSTEM('gzip group-output-time'//str_trim(1:k)//'.csv' ) ! factor of 12 reduction in disc space.
         IFILE_COUNT=IFILE_COUNT+1 ! next file.
       endif
! ************* save to disc***************

         END DO ! DO ITIME=1,NTIME

         CALL MAKE_LONG(T_NEW_LONG, T_NEW, NX,NY,NZ,NG, NONODS_NG) 

! SOLVE A PROBLEM USING ROM ********************

         IF(EIGEN_VALUE) PRINT *,'KEFF=',KEFF

! ******************************TEST THE INTERPOLATION*****************
         print *,'beta,beta2:',beta,beta2
         print *,'rday,1000.0 * (1.0-RDAY) + 1000.0:',rday,1000.0 * (1.0-RDAY) + 1000.0
         print *,'ntime, acctim, acctim/one_day, ntime*(one_day/max(toler,acctim)):', &
                  ntime, acctim, acctim/one_day, ntime*(one_day/max(toler,acctim))
      DO I=2,nx-1
      DO J=2,ny-1
      DO K=2,nz-1
         n1(i,j,k)=sum(t_new(i,j,k,1:4))
         n2(i,j,k)=sum(t_new(i,j,k,5:8))
      end do
      end do
      end do
      print *,'minval(n1),maxval(n1):',minval(n1),maxval(n1)
      print *,'minval(n2),maxval(n2):',minval(n2),maxval(n2)
      print *,'sum(n1),sum(n2),sum(n1)+sum(n2):',sum(n1),sum(n2),sum(n1)+sum(n2)
         DO IG=1,NG
            PRINT *,'GROUP=',IG
      print *,'minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig)):',minval(t_new(:,:,:,ig)),maxval(t_new(:,:,:,ig))
      print *,'sum(t_new(:,:,:,ig)):',sum(t_new(:,:,:,ig))
!       print *,'t_new(:,:,:,ig):',t_new(:,:,:,ig)
       IF(NX<55) THEN ! Only print small problems to screen'
            DO K=2,NZ-1
            do j=ny-1,2,-1
!               PRINT *,(T_NEW( I, J, K, IG ),I=2,NX-1)
!               PRINT *,(real(T_NEW( I, J, K, IG ), sp),I=2,NX-1) ! in single procession for printing
!               do i=2,nx-1
               WRITE(*, '(200E10.2)' ) (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
!               print *, (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
!               end do
!               WRITE(*, '(F10.6)' ) ( T_NEW( I, J, K, IG ),I=2,nx-1) ! in single procession for printing
!               i=2
!               WRITE(*, '(E12.6)' ) T_NEW( I, J, K, IG ) ! in single procession for printing
            end do
            end do
!           stop 382
      ENDIF ! IF(NX<55) THEN
         END DO
! 
         IF(EIGEN_VALUE) PRINT *,'KEFF=',KEFF

! output to disc the group data
    if(.true.) then
!         open(27, file='run/group-output-end.csv')
         open(27, file='group-output-end.csv')
            WRITE(27, *) '# DX,DY,DZ, NX,NY,NZ:', DX,DY,DZ, NX-2,NY-2,NZ-2
            WRITE(27, *) '# there are 8 groups: HOME-S, HOME-E, HOME-I, HOME-R, MOBILE-S, MOBILE-E, MOBILE-I, MOBILE-R'
            
            WRITE(27, *) '# DO K=1,NZ'
            WRITE(27, *) '# do j=1,NY,1,-1'
            WRITE(27, *) '#    WRITE(27, "(200E14.4)" ) (SIGN_SMALL(T_NEW( I, J, K, IG )) ,I=1,NX) ! in single procession for printing' 
            WRITE(27, *) '# end do'
            WRITE(27, *) '# end do'
!         ELSE
!            open(27, file='group-output-time.csv', Access = 'append')
!         ENDIF
         WRITE(27, *) '# TIME=', ACCTIM
         DO IG=1,NG
            WRITE(27, *) 'GROUP=',IG
            DO K=2,NZ-1
            do j=2,ny-1
!               PRINT *,(T_NEW( I, J, K, IG ),I=2,NX-1)
!               PRINT *,(real(T_NEW( I, J, K, IG ), sp),I=2,NX-1) ! in single procession for printing
!               do i=2,nx-1
!               WRITE(27, '(200E12.4)' ) (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
               WRITE(27, '(200E14.4)' ) (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
!               print *, (T_NEW( I, J, K, IG ) ,I=2,nx-1) ! in single procession for printing
!               end do
!               WRITE(*, '(F10.6)' ) ( T_NEW( I, J, K, IG ),I=2,nx-1) ! in single procession for printing
!               i=2
!               WRITE(*, '(E12.6)' ) T_NEW( I, J, K, IG ) ! in single procession for printing
            end do
            end do
!           stop 382
         END DO
         close(27) 
         CALL SYSTEM('rm group-output-end.csv.gz' )
         CALL SYSTEM('gzip group-output-end.csv' )

         open(27, file='group-time-integral')
         DO IG=1,NG
            ACCTIM=0.0
            DO ITIME=1,NTIME
               ACCTIM=ACCTIM+DT
               WRITE(27, * ) ACCTIM, INTEGRAL_T_NEW(IG,ITIME) 
            END DO
            WRITE(27, * ) 
         END DO
         close(27) 

         open(27, file='group-time-sample-pt')
         DO IG=1,NG
            ACCTIM=0.0
            DO ITIME=1,NTIME
               ACCTIM=ACCTIM+DT
               WRITE(27, * ) ACCTIM, SAMPLE_T_NEW(IG,ITIME) 
            END DO
            WRITE(27, * ) 
         END DO
         close(27) 

    endif
 

       STOP
       END PROGRAM diffusion3d_virus
! 
! 
       REAL FUNCTION SIGN_SMALL(T_NEW)
       REAL T_NEW
       SIGN_SMALL=T_NEW

       if(abs(T_NEW)<1.e-25) SIGN_SMALL=0.0
       END 
! 
! 
      CHARACTER*240 FUNCTION CHARA(CURPRO)
!     CURPRO must be the same in calling subroutine.
      INTEGER MXNLEN
      PARAMETER(MXNLEN=5)
      INTEGER DIGIT(MXNLEN),II,J,TENPOW,CURPRO,K
      CHARACTER STRI*240
      CHARACTER*1 CHADIG(MXNLEN)
      LOGICAL SWITCH
!     MXNLEN=max no of digits in integer I.
!     This function converts integer I to a string.
      II=CURPRO
      SWITCH=.FALSE.
      K=0
      STRI=' '
      do  J=MXNLEN,1,-1! Was loop 10
         TENPOW=10**(J-1)
         DIGIT(J)=II/TENPOW
         II=II-DIGIT(J)* TENPOW
!     ewrite(3,*)'DIGIT*J)=',DIGIT(J),' TENOW=',TENPOW
!     STRI=STRI//CHAR(DIGIT(J)+48)
         IF(DIGIT(J).EQ.0) CHADIG(J)='0'
         IF(DIGIT(J).EQ.1) CHADIG(J)='1'
         IF(DIGIT(J).EQ.2) CHADIG(J)='2'
         IF(DIGIT(J).EQ.3) CHADIG(J)='3'
         IF(DIGIT(J).EQ.4) CHADIG(J)='4'
         IF(DIGIT(J).EQ.5) CHADIG(J)='5'
         IF(DIGIT(J).EQ.6) CHADIG(J)='6'
         IF(DIGIT(J).EQ.7) CHADIG(J)='7'
         IF(DIGIT(J).EQ.8) CHADIG(J)='8'
         IF(DIGIT(J).EQ.9) CHADIG(J)='9'
!         Q=INDEX(STRI,' ')-1
!         IF(Q.LT.1) STRI=CHADIG
!         IF(Q.GE.1) STRI=STRI(1:Q)//CHADIG
!         IF(DIGIT(J).NE.0) SWITCH=.TRUE.
!         IF(SWITCH) K=K+1
      end do ! Was loop 10
       k=0
      do  J=MXNLEN,1,-1! Was loop 20
         IF(SWITCH) THEN 
           STRI=STRI(1:K)//CHADIG(J)
           K=K+1
         ENDIF
         IF((CHADIG(J).NE.'0').AND.(.NOT.SWITCH)) THEN 
           SWITCH=.TRUE.
           STRI=CHADIG(J)
           K=1
         ENDIF
      end do ! Was loop 20
       IF(K.EQ.0) THEN 
         STRI='0'
          K=1
        ENDIF
!       ewrite(3,*)'STRI=',STRI,' K=',K
!       ewrite(3,*)'STRI(MXNLEN-K+1:MXNLEN)=',STRI(MXNLEN-K+2:MXNLEN+1)
!       CHARA=STRI(MXNLEN-K+2:MXNLEN+1)
        CHARA=STRI(1:K)
       END
! In pythin this is:
! ALPHA_INTERP = SIMPLE_INTERPOLATION(MATERIAL_INTERP,   MATERIAL,ALPHA_COEFFS, NEXAMPLE,NALPHA,NMATERIAL) 
       SUBROUTINE SIMPLE_INTERPOLATION(ALPHA_INTERP, MATERIAL_INTERP,   MATERIAL,ALPHA_COEFFS, NEXAMPLE,NALPHA,NMATERIAL) 
! This sub calculates the interpolated ALPHA_INTERP coeffs at a point MATERIAL_INTERP in material space. 
! It does this by interpolating NEXAMPLE points in material space at positions defined by MATERIAL and 
! at these positions we have the ALPHA_COEFFS which we are interpolating. 
! NMATERIAL are the number of different materials. 
! This subroutine uses an inverse distance weighting to perform the interpolation. 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NEXAMPLE, NALPHA,NMATERIAL
       REAL, intent( out ) :: ALPHA_INTERP(NALPHA)
       REAL, intent( in ) :: MATERIAL(NMATERIAL,NEXAMPLE),ALPHA_COEFFS(NALPHA,NEXAMPLE), MATERIAL_INTERP(NMATERIAL)
! local variables...
       REAL TOLER,INFINY
       PARAMETER(TOLER=1.0E-10,INFINY=1.0E+15)
       REAL, ALLOCATABLE :: D2(:),D2_SMALL(:),D_SMALL(:),W(:)
       INTEGER, ALLOCATABLE :: EXAMP_INDEX(:)
       INTEGER IEXAMPLE,IMAT, IEXAMP_SMALL, IALPHA
       REAL D2_SMALLEST

       ALLOCATE(D2(NEXAMPLE),D2_SMALL(NMATERIAL+1),D_SMALL(NMATERIAL+1),W(NMATERIAL+1))
       ALLOCATE(EXAMP_INDEX(NMATERIAL+1))
       
       DO IEXAMPLE=1,NEXAMPLE
          D2(IEXAMPLE)=SUM((MATERIAL_INTERP(:)-MATERIAL(:,IEXAMPLE))**2)
       END DO 
       
! find NMATERIAL+1 smallest values...       
       DO IMAT=1,NMATERIAL+1
          D2_SMALLEST=INFINY 
          DO IEXAMPLE=1,NEXAMPLE
             IF(D2(IEXAMPLE).LT.D2_SMALLEST) THEN ! insert into list
                D2_SMALLEST=D2(IEXAMPLE) 
                IEXAMP_SMALL = IEXAMPLE
             ENDIF
          END DO 
          D2_SMALL(IMAT)=D2_SMALLEST
          EXAMP_INDEX(IMAT)=IEXAMP_SMALL
          D2(IEXAMP_SMALL)=INFINY ! Set to large no so not to include this in list again. 
       END DO 

       DO IMAT=1,NMATERIAL+1
          D_SMALL(IMAT)=MAX(SQRT(D2_SMALL(IMAT)),TOLER)
       END DO
       W(:)=(1./D_SMALL(:))/ SUM(1./D_SMALL(:))

       DO IALPHA=1,NALPHA
          ALPHA_INTERP(IALPHA)= SUM( W(:)*ALPHA_COEFFS(IALPHA,EXAMP_INDEX(:)) )
       END DO
       
       END SUBROUTINE SIMPLE_INTERPOLATION
! 
! 
! 
! 
       SUBROUTINE KEEP_RANDOM_NO(SIG_RAN_NO,IRAN_NO,RAN_NO,RECALL)
! GENERATE AND STORE RANDOM NUMBERS. 
       INTEGER IRAN_NO
       REAL SIG_RAN_NO, RAN_NO(*)
       LOGICAL RECALL
! 
       IRAN_NO=IRAN_NO+1
       IF(RECALL) THEN
          SIG_RAN_NO = RAN_NO(IRAN_NO) 
       ELSE
          CALL RANDOM_NUMBER(SIG_RAN_NO) 
          RAN_NO(IRAN_NO) = SIG_RAN_NO
       ENDIF
       END SUBROUTINE KEEP_RANDOM_NO
! 
! 
! 
! 
          SUBROUTINE GET_CT_T_NEW(WEIGHT,NO_NUR,NCOLM,W_EXP,NLAY,NLAYER, &
              ALPHA_NEW, NALPHA, T_NEW_LONG, CT, NONODS_NG, SCALE_OUT,MAX_SCALE_OUT,MIN_SCALE_OUT ) 
! Form CT and T_NEW
! for ann:
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! This sub calculates the neuron values. 
! If W_EXP=0.0 find exponent of NEURONS otherwise dont. 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NO_NUR,NCOLM, NLAYER
       REAL, intent( in ) :: WEIGHT(NCOLM),W_EXP(NLAYER)
       INTEGER , intent( in ):: NLAY(NLAYER)
! for this subroutine:
          INTEGER, intent( in ) :: NALPHA,NONODS_NG, SCALE_OUT
          REAL, intent( in ) :: ALPHA_NEW(NALPHA),MAX_SCALE_OUT,MIN_SCALE_OUT
          REAL, intent( out ) :: T_NEW_LONG(NONODS_NG), CT(NALPHA,NONODS_NG)
             
! LOCAL VARIABLES...
          real EPSILON_PERT
          parameter(EPSILON_PERT=1.e-4)
          INTEGER NLAY_RECALL, ISTART_NEU,ISTART_WEIGHT, NO_NUR_RECALL, NCOLM_RECALL,IALPHA
          REAL, ALLOCATABLE :: NEURAL(:), NEURAL_TEMP(:)

          ALLOCATE(NEURAL(NO_NUR),NEURAL_TEMP(NO_NUR))

!          PRINT *,'NEXAMPLE,NBATCH:',NEXAMPLE,NBATCH
          NLAY_RECALL=(NLAYER-1)/2+1
          ISTART_NEU= SUM(NLAY( 1:(NLAYER-1)/2 ))+1
          NO_NUR_RECALL = NO_NUR - ISTART_NEU +1
          NCOLM_RECALL = NCOLM/2
          ISTART_WEIGHT=NCOLM/2 +1

          NEURAL(ISTART_NEU:ISTART_NEU+NALPHA-1)=ALPHA_NEW(:)
          CALL GETNEUVALS_FAST(NEURAL(ISTART_NEU),WEIGHT(ISTART_WEIGHT),&
               NO_NUR_RECALL,NCOLM_RECALL,W_EXP((NLAYER-1)/2 +1),NLAY((NLAYER-1)/2 +1),NLAY_RECALL)
          T_NEW_LONG(1:NONODS_NG) = NEURAL(NO_NUR-NONODS_NG+1:NO_NUR)
          IF(SCALE_OUT.NE.0) T_NEW_LONG = MIN_SCALE_OUT  + T_NEW_LONG*(MAX_SCALE_OUT - MIN_SCALE_OUT)
! unscale... 

          DO IALPHA=1,NALPHA  
             NEURAL_TEMP=NEURAL
             NEURAL_TEMP(ISTART_NEU-1+IALPHA) = NEURAL_TEMP(ISTART_NEU-1+IALPHA) + EPSILON_PERT
             CALL GETNEUVALS_FAST(NEURAL_TEMP(ISTART_NEU),WEIGHT(ISTART_WEIGHT),&
                  NO_NUR_RECALL,NCOLM_RECALL,W_EXP((NLAYER-1)/2 +1),&
                  NLAY((NLAYER-1)/2 +1),NLAY_RECALL)
             CT(IALPHA,:) =  (NEURAL_TEMP(NO_NUR-NONODS_NG+1:NO_NUR) - NEURAL(NO_NUR-NONODS_NG+1:NO_NUR))/EPSILON_PERT
             IF(SCALE_OUT.NE.0) CT(IALPHA,:) =CT(IALPHA,:) * (MAX_SCALE_OUT - MIN_SCALE_OUT)
          END DO
          
          END SUBROUTINE GET_CT_T_NEW
! 
! 
! 
! In python:  
! ALPHA_NEW = ONE_AE_ITERATION( ALPHA, A_MAT,SIGMA_S_OFF,B_MAT,C, U_NEW, U_OLD, S, DIAG_REG_TOLER, NX,NY,NZ,NG, NALPHA,NONODS_NG) 
         SUBROUTINE ONE_AE_ITERATION(ALPHA_NEW,  ALPHA, A_MAT,SIGMA_S_OFF,B_MAT,CT, U_NEW, U_OLD, S,&
                                      DIAG_REG_TOLER, NX,NY,NZ,NG, NALPHA,NONODS_NG) 
	 IMPLICIT NONE
	 INTEGER, INTENT(IN) :: NX,NY,NZ,NG, NALPHA,NONODS_NG
! Perform matrix vector multiplication.
! B=A*X
! DIAG_REG_TOLER contains regularization terms.
         REAL, INTENT(OUT) :: ALPHA_NEW(NALPHA)
         REAL, INTENT(IN) :: U_NEW(NX,NY,NZ,NG), U_OLD(NX,NY,NZ,NG), S(NX,NY,NZ,NG), ALPHA(NALPHA), DIAG_REG_TOLER(7) 
         REAL, INTENT(IN) :: CT(NALPHA,NONODS_NG), A_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1), B_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1)
         REAL, INTENT(IN) :: SIGMA_S_OFF(NX,NY,NZ,NG,NG) 
! Local variables...
         REAL, ALLOCATABLE :: CTAC(:,:), D(:,:,:,:), A_U(:,:,:,:),P(:),C_P(:,:,:,:),A_C_P(:,:,:,:),B_U_OLD(:,:,:,:), C(:,:) 
         REAL, ALLOCATABLE :: RHS(:), S_LONG(:),A_U_LONG(:),A_C_P_LONG(:), D_LONG(:), B_U_OLD_LONG(:)
         REAL, ALLOCATABLE :: CT_A_U(:), CT_B_U_OLD(:), CT_A_C_P(:), C_P_LONG(:), DALPHA(:)
         REAL, ALLOCATABLE :: SIGMA_S_ZERO(:,:,:,:,:)
         INTEGER I,NOD,IALPHA

         allocate( CTAC(NALPHA,NALPHA),D(NX,NY,NZ,NG), A_U(NX,NY,NZ,NG),P(NALPHA),C_P(NX,NY,NZ,NG) ) 
         allocate( A_C_P(NX,NY,NZ,NG), B_U_OLD(NX,NY,NZ,NG),  RHS(NALPHA), C(NONODS_NG,NALPHA) ) 
         allocate( S_LONG(NONODS_NG),A_U_LONG(NONODS_NG),A_C_P_LONG(NONODS_NG), D_LONG(NONODS_NG), B_U_OLD_LONG(NONODS_NG)  ) 
         allocate( CT_A_U(NALPHA),CT_B_U_OLD(NALPHA),CT_A_C_P(NALPHA), C_P_LONG(NONODS_NG), DALPHA(NALPHA) )
         ALLOCATE( SIGMA_S_ZERO(NX,NY,NZ,NG,NG)) 
         SIGMA_S_ZERO=0.0

         DO NOD=1,NONODS_NG
            DO IALPHA=1,NALPHA
               C(NOD,IALPHA) = CT(IALPHA,NOD) 
            END DO
         END DO
         ! print *,'here1'
       
         CALL MATRIX_VEC_MULT_TIME_STEPING( A_U, A_MAT,SIGMA_S_OFF,U_NEW, NX,NY,NZ,NG) 
         CALL MAKE_LONG(A_U_LONG, A_U, NX,NY,NZ,NG,NONODS_NG) 
         CT_A_U=matmul(CT,A_U_LONG)

         CALL MATRIX_VEC_MULT_TIME_STEPING( B_U_OLD, B_MAT,SIGMA_S_ZERO,U_OLD, NX,NY,NZ,NG) 
         CALL MAKE_LONG(B_U_OLD_LONG, B_U_OLD, NX,NY,NZ,NG,NONODS_NG) 
         CT_B_U_OLD=matmul(CT,B_U_OLD_LONG)

         CALL MAKE_LONG(S_LONG, S,NX,NY,NZ,NG,NONODS_NG) 

         RHS=-CT_A_U + CT_B_U_OLD + matmul(CT,S_LONG)

! The matrix FORMATION:
         do i=1,nalpha
            p=0.0
            p(i)=1.0

            C_P_LONG=matmul(C,P)
            CALL MAKE_GRID( C_P, C_P_LONG, NX,NY,NZ,NG,NONODS_NG)

            CALL MATRIX_VEC_MULT_TIME_STEPING( A_C_P, A_MAT,SIGMA_S_OFF,C_P, NX,NY,NZ,NG) 
            CALL MAKE_LONG(A_C_P_LONG, A_C_P, NX,NY,NZ,NG,NONODS_NG)
            CT_A_C_P=matmul(CT,A_C_P_LONG)     
            CTAC(:,I)=CT_A_C_P(:)                  
         end do

! Regularize CTAC and RHS vector before solution:
         CALL REGULARIZE_CTAC_RHS(DIAG_REG_TOLER, CTAC, ALPHA, NALPHA) 

! Solve CTAC alpha = rhs:      
!         CALL SMLINN(CTAC,ALPHA_NEW,RHS,NALPHA,NALPHA) 
         CALL SMLINN(CTAC,DALPHA,RHS,NALPHA,NALPHA) 
         ALPHA_NEW =  DALPHA + ALPHA 

         END SUBROUTINE ONE_AE_ITERATION
! 
! 
! 
! 
         subroutine REGULARIZE_CTAC_RHS(DIAG_REG_TOLER, CTAC, ALPHA, NALPHA) 
! Regularize CTAC and RHS vector before solution:
! ALPHA is the lattest soln in compressed space. 
         implicit none
         INTEGER NALPHA
         REAL CTAC(NALPHA,NALPHA), ALPHA(NALPHA), DIAG_REG_TOLER(7) 
! local variables...
         REAL, ALLOCATABLE :: DIAG(:) 
         INTEGER IALPHA

         allocate( DIAG(NALPHA) )
! 
         DO IALPHA=1,NALPHA
            DIAG(IALPHA)= &
         + DIAG_REG_TOLER(1)                             + MAXVAL(ABS(CTAC(:,:)))*DIAG_REG_TOLER(2) &
         + MAXVAL(ABS(CTAC(IALPHA,:)))*DIAG_REG_TOLER(3) + ABS(CTAC(IALPHA,IALPHA))*DIAG_REG_TOLER(4) &
         + SUM(ABS(CTAC(:,:)))*DIAG_REG_TOLER(5)         + SUM(ABS(CTAC(IALPHA,:)))*DIAG_REG_TOLER(6) 
         END DO

         DO IALPHA=1,NALPHA
            CTAC(IALPHA,IALPHA) = CTAC(IALPHA,IALPHA) + DIAG(IALPHA)
!            RHS(IALPHA) = RHS(IALPHA) + DIAG(IALPHA)*ALPHA(IALPHA) 
         END DO     
         end subroutine REGULARIZE_CTAC_RHS
! 
! 
! 
! 
         subroutine MAKE_LONG(b_long, B_GRID, NX,NY,NZ,NG,NONODS_NG) 
         implicit none
         integer NX,NY,NZ,NG,NONODS_NG
         real B_GRID(NX,NY,NZ,NG)
         real b_long(NONODS_NG)
! Local variables...
         integer i,j,k,ig, nx2,ny2,nz2
         nx2=nx-2; ny2=ny-2; nz2=nz-2
         do ig=1,ng
            do k=1,nz2
               do j=1,ny2
                  do i=1,nx2
                     b_long(i+(j-1)*nx2 + (k-1)*nx2*ny2 + (ig-1)*nx2*ny2*nz2) = b_grid(i+1,j+1,k+1,ig) 
                  end do
               end do
            end do
         end do
         
         end subroutine MAKE_LONG
! 
! 
! 
! 
         subroutine MAKE_GRID(b_grid, B_LONG, NX,NY,NZ,NG,NONODS_NG) 
         implicit none
         integer NX,NY,NZ,NG,NONODS_NG
         real B_LONG(NONODS_NG)
         real b_grid(NX,NY,NZ,NG) 
! Local variables...
         integer i,j,k,ig, nx2,ny2,nz2
         nx2=nx-2; ny2=ny-2; nz2=nz-2
         b_grid=0.0
         do ig=1,ng
            do k=1,nz2
               do j=1,ny2
                  do i=1,nx2
                     b_grid(i+1,j+1,k+1,ig) = b_long(i+(j-1)*nx2+ (k-1)*nx2*ny2 + (ig-1)*nx2*ny2*nz2) 
                  end do
               end do
            end do
         end do
         
         end subroutine MAKE_GRID
! 
! 
!     
!	  
        SUBROUTINE SMLINN(A,X,B,NMX,N)
        IMPLICIT NONE
        INTEGER NMX,N
        REAL A(NMX,NMX),X(NMX),B(NMX)
        REAL R
        INTEGER K,I,J
!     Form X = A^{-1} B
!     Useful subroutine for inverse
!     This sub overwrites the matrix A. 
        DO K=1,N-1
           DO I=K+1,N
              A(I,K)=A(I,K)/A(K,K)
           END DO
           DO J=K+1,N
              DO I=K+1,N
                 A(I,J)=A(I,J) - A(I,K)*A(K,J)
              END DO
           END DO
        END DO
!     
!     Solve L_1 x=b
        DO I=1,N
           R=0.
           DO J=1,I-1
              R=R+A(I,J)*X(J)
           END DO
           X(I)=B(I)-R
        END DO
!     
!     Solve U x=y
        DO I=N,1,-1
           R=0.
           DO J=I+1,N
              R=R+A(I,J)*X(J)
           END DO
           X(I)=(X(I)-R)/A(I,I)
        END DO
        RETURN
        END
!     
! 
! In phython use:
! b = MATRIX_VEC_MULT_TIME_STEPING( A,SIGMA,X, NX,NY,NZ,NG) 
         SUBROUTINE MATRIX_VEC_MULT_TIME_STEPING( B, A,SIGMA_S_OFF,X, NX,NY,NZ,NG) 
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: NX,NY,NZ,NG
! Perform matrix vector multiplication.
! B=A*X
          REAL, INTENT(OUT) :: B(NX,NY,NZ,NG)
          REAL, INTENT(IN) :: X(NX,NY,NZ,NG), A(NX,NY,NZ,NG, -1:1, -1:1, -1:1),SIGMA_S_OFF(NX,NY,NZ,NG,NG)
! Local variables...
       INTEGER I,J,K,IG, II,JJ,KK, JG

       B=0.0
 
       DO IG=1,NG
       DO JG=1,NG
          B(:,:,:,IG)=B(:,:,:,IG) - SIGMA_S_OFF(:,:,:,IG,JG)*X(:,:,:,JG) 
       END DO
       END DO


       DO KK=-1,1
       DO JJ=-1,1
       DO II=-1,1

          DO IG = 1,NG
          DO K = 2,NZ-1
          DO J = 2,NY-1  
          DO I = 2,NX-1 

            B(I,J,K,IG) = B(I,J,K,IG) + A(I,J,K,IG,II,JJ,KK)*X(I+II,J+JJ,K+KK,IG) 
          END DO
          END DO 
          END DO 
          END DO 

       END DO
       END DO 
       END DO 
       END SUBROUTINE MATRIX_VEC_MULT_TIME_STEPING
! 
! 
! 
! 
       SUBROUTINE EIG_OR_TIME_STEPING_DIFFUSION_CALC(T_NEW,KEFF, T_INITIAL,&
                  NTIME,NITS_SOLV,NITS_SOLV_NG,NITS_EIG,ERROR_SOLV,ERROR_SOLV_NG,&
                  ERROR_EIG,EIGEN_START, RELAX_SOLV, DX,DY,DZ,DT,V_N,SIGMAT,&
                  SIGMA_S_OFF,K,SIGMA_F,S, NX,NY,NZ,NG )
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: NX,NY,NZ,NG
! NITS=no of non-linear iterations 
! NTIME=NO OF TIME STEPS
! DT=TIME STEP SIZE (CAN USE DT=1.E+15 for steady state cacls) 
          REAL, INTENT(OUT) :: T_NEW(NX,NY,NZ,NG),KEFF
          REAL, INTENT(IN) :: T_INITIAL(NX,NY,NZ,NG) 
          REAL, INTENT(IN) :: RELAX_SOLV,DX,DY,DZ,DT,V_N(NG), ERROR_SOLV,ERROR_SOLV_NG,ERROR_EIG
          REAL, INTENT(IN) :: SIGMAT(NX,NY,NZ,NG),K(NX,NY,NZ,NG),SIGMA_F(NX,NY,NZ,NG,NG),&
                              SIGMA_S_OFF(NX,NY,NZ,NG,NG),S(NX,NY,NZ,NG) 
	  INTEGER, INTENT(IN) :: NTIME,NITS_SOLV,NITS_SOLV_NG,NITS_EIG
          LOGICAL, INTENT(IN) :: EIGEN_START
! Local variables...
       REAL, ALLOCATABLE :: V(:,:,:,:),A_MAT_T(:,:,:,:),A_MAT(:,:,:,:,:,:,:), B_MAT(:,:,:,:,:,:,:)
       INTEGER ITIME,ITS,II,I,JJ,J,IG,JG,ITS_EIG
       REAL LAMBDA,KEFF2

       IF(NITS_EIG==0) THEN
          CALL TIME_STEPING_DIFFUSION_CALC(T_NEW, T_INITIAL, NTIME,NITS_SOLV,NITS_SOLV_NG,RELAX_SOLV,ERROR_SOLV,ERROR_SOLV_NG, &
                         DX,DY,DZ,DT,V_N,SIGMAT,SIGMA_S_OFF,K,S, NX,NY,NZ,NG)
       ELSE
          ALLOCATE(V(NX,NY,NZ,NG),A_MAT_T(NX,NY,NZ,NG), A_MAT(NX,NY,NZ,NG,-1:1,-1:1, -1:1),B_MAT(NX,NY,NZ,NG,-1:1,-1:1, -1:1) ) 
          KEFF=1.E+10
          T_NEW=T_INITIAL
          IF(EIGEN_START) T_NEW=1.0
          CALL GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( A_MAT, B_MAT,  &
                         DX,DY,DZ,DT,V_N,SIGMAT,K, NX,NY,NZ,NG) 
          DO ITS_EIG=1,NITS_EIG
!             V=SIGMA_F*T_NEW ! Matrix vector
             V=0.0
             DO IG=1,NG
                DO JG=1,NG
                   V(:,:,:,IG)=V(:,:,:,IG) + SIGMA_F(:,:,:,IG,JG)*T_NEW(:,:,:,JG) 
                END DO
             END DO
                  
             CALL TIME_STEPING_DIFFUSION_CALC(T_NEW, T_INITIAL, NTIME,NITS_SOLV,NITS_SOLV_NG,RELAX_SOLV,ERROR_SOLV,ERROR_SOLV_NG, &
                         DX,DY,DZ,DT,V_N,SIGMAT,SIGMA_S_OFF,K,V, NX,NY,NZ,NG)
! Normalize...
             V=0.0
             DO IG=1,NG
                DO JG=1,NG
                   V(:,:,:,IG)=V(:,:,:,IG) + SIGMA_F(:,:,:,IG,JG)*T_NEW(:,:,:,JG) 
                END DO
             END DO
             T_NEW = T_NEW / SUM(V) 
!             T_NEW = T_NEW / SUM(SIGMA_F*T_NEW) 
! Eigen-value...
             CALL MATRIX_VEC_MULT_TIME_STEPING( A_MAT_T, A_MAT,SIGMA_S_OFF,T_NEW, NX,NY,NZ,NG) 
             V=0.0
             DO IG=1,NG
                DO JG=1,NG
                   V(:,:,:,IG)=V(:,:,:,IG) + SIGMA_F(:,:,:,IG,JG)*T_NEW(:,:,:,JG) 
                END DO
             END DO
	     !PRINT *, SUM(A_MAT_T), SUM(V)
             LAMBDA = SUM(A_MAT_T) / SUM(V)
!             LAMBDA = SUM(A_MAT_T) / SUM(SIGMA_F*T_NEW)
             KEFF2=1.0/MAX(1.E-10, LAMBDA) 
!              print *,'its_eig,keff2=',its_eig,keff
              print *,'its_eig,NITS_EIG,keff2,ABS(KEFF2-KEFF),ERROR_EIG=', &
                       its_eig,NITS_EIG,keff, ABS(KEFF2-KEFF),ERROR_EIG
		!KEFF_ALL(ITS_EIG)=KEFF2
             IF(ABS(KEFF2-KEFF).LT.ERROR_EIG) THEN
                KEFF=KEFF2
                EXIT
             ELSE
                KEFF=KEFF2
             ENDIF
          END DO
       ENDIF
       END SUBROUTINE 

! 
! 
! 
! In phython use:
! A, B_MAT, RHS = GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( T_INITIAL, DX,DY,DZ,DT,V_N,SIGMAT,K,S, NX,NY,NZ,NG)
         SUBROUTINE GET_MATRIX_TIME_STEPING_DIFFUSION_CALC( A_MAT, B_MAT,   &
                         DX,DY,DZ,DT,V_N,SIGMAT,KDIFF, NX,NY,NZ,NG) 
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: NX,NY,NZ,NG
! NITS=no of non-linear iterations 
! NTIME=NO OF TIME STEPS
! DT=TIME STEP SIZE (CAN USE DT=1.E+15 for steady state cacls) 
          REAL, INTENT(OUT) :: A_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1), B_MAT(NX,NY,NZ,NG, -1:1, -1:1, -1:1)
          REAL, INTENT(IN) :: DX,DY,DZ,DT,V_N(NG)
          REAL, INTENT(IN) :: SIGMAT(NX,NY,NZ,NG),KDIFF(NX,NY,NZ,NG)
! Local variables...
       INTEGER I,J,K, IG
       REAL DIFF_I 
       REAL A_I_J_K, A_IP1_J_K, A_IM1_J_K, A_I_JP1_K, A_I_JM1_K, A_I_J_KP1, A_I_J_KM1
       A_MAT=0.0
       B_MAT=0.0

	DO IG = 1,NG
	DO K = 2,NZ-1
	DO J = 2,NY-1  
	DO I = 2,NX-1 

            DIFF_I = 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I-1,J,K,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I+1,J,K,IG),0.0))/(DX**2) &
                   + 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I,J-1,K,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I,J+1,K,IG),0.0))/(DY**2) &
                   + 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I,J,K-1,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I,J,K+1,IG),0.0))/(DZ**2) 

            A_I_J_K = DIFF_I + SIGMAT(I,J,K,IG) + (1./(DT*V_N(IG))) 
            A_IP1_J_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I+1,J,K,IG),0.0) /(DX**2) 
            A_IM1_J_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I-1,J,K,IG),0.0) /(DX**2)
            A_I_JP1_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J+1,K,IG),0.0) /(DY**2) 
            A_I_JM1_K =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J-1,K,IG),0.0) /(DY**2) 
            A_I_J_KP1 =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J,K+1,IG),0.0) /(DZ**2) 
            A_I_J_KM1 =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J,K-1,IG),0.0) /(DZ**2) 

            A_MAT(I,J,K,IG,0,0,0)   = A_I_J_K  
            A_MAT(I,J,K,IG,1,0,0) = A_IP1_J_K
            A_MAT(I,J,K,IG,-1,0,0) = A_IM1_J_K
            A_MAT(I,J,K,IG,0,1,0) = A_I_JP1_K
            A_MAT(I,J,K,IG,0,-1,0) = A_I_JM1_K

            A_MAT(I,J,K,IG,0,0,1)  = A_I_J_KP1
            A_MAT(I,J,K,IG,0,0,-1) = A_I_J_KM1

            B_MAT(I,J,K,IG,0,0,0) = 1./(V_N(IG)*DT) 

	END DO
	END DO 
	END DO 
	END DO 

	  END SUBROUTINE GET_MATRIX_TIME_STEPING_DIFFUSION_CALC
! 
! 
! 
! In phython use:
! T_NEW = TIME_STEPING_DIFFUSION_CALC(T_INITIAL,NTIME,NITS,RELAX,DX,DY,DZ,DT,V_N,SIGMAT,K,S, NX,NY,NZ,NG)
         SUBROUTINE TIME_STEPING_DIFFUSION_CALC(T_NEW,T_INITIAL, NTIME,NITS,NITS_SOLV_NG,RELAX,ERROR_SOLV,ERROR_SOLV_NG, &
                         DX,DY,DZ,DT,V_N,SIGMAT,SIGMA_S_OFF,KDIFF,S, NX,NY,NZ,NG) 
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: NX,NY,NZ,NG,NITS,NITS_SOLV_NG
! NITS=no of non-linear iterations 
! NTIME=NO OF TIME STEPS
! DT=TIME STEP SIZE (CAN USE DT=1.E+15 for steady state cacls) 
! NITS,NITS_NG are the max no of single group iterations and overall iterations. 
! ERROR_SOLV,ERROR_SOLV_NG are the error tolerances of the solver of single group iterations and overall iterations.
          REAL, INTENT(OUT) :: T_NEW(NX,NY,NZ,NG)
          REAL, INTENT(IN) :: T_INITIAL(NX,NY,NZ,NG) 
          REAL, INTENT(IN) :: RELAX,DX,DY,DZ,DT,V_N(NG), ERROR_SOLV, ERROR_SOLV_NG
          REAL, INTENT(IN) :: SIGMAT(NX,NY,NZ,NG),SIGMA_S_OFF(NX,NY,NZ,NG,NG),KDIFF(NX,NY,NZ,NG),S(NX,NY,NZ,NG) 
	  INTEGER, INTENT(IN) :: NTIME
! Local variables...
       REAL RELAX_NG
       PARAMETER(RELAX_NG=1.0) ! Multi-group relaxation term.
       INTEGER ISWITCH_FBGS
       PARAMETER(ISWITCH_FBGS=0) ! Switch on FBGS (ISWITCH_FBGS=1 switched on FBGS and ISWITCH_FBGS=0 used forward block Gauss Siedel.)
       REAL, ALLOCATABLE :: T_OLD(:,:,:,:), A_MAT_I_J_K(:,:,:,:), A_MAT_IP1_J_K(:,:,:,:)
       REAL, ALLOCATABLE :: A_MAT_IM1_J_K(:,:,:,:),A_MAT_I_JP1_K(:,:,:,:),A_MAT_I_JM1_K(:,:,:,:), A_MAT_I_J_KP1(:,:,:,:)
       REAL, ALLOCATABLE :: A_MAT_I_J_KM1(:,:,:,:)
       REAL, ALLOCATABLE :: RHS_SCAT(:,:,:), T_NEW_TEMP_NG(:,:,:)
       INTEGER ITIME,ITS,ITS_NG,II,I,JJ,J,KK,K,IG,JG
       REAL DIFF_I, R_T_NEW, R_T_MAX_DIF, R_T_MAX_DIF_NG, T_MAX
       REAL A_I_J_K, A_IP1_J_K, A_IM1_J_K, A_I_JP1_K, A_I_JM1_K, A_I_J_KP1, A_I_J_KM1
       INTEGER IIGG, IG_START, IG_FINISH, IG_STEP

       ALLOCATE(T_OLD(NX,NY,NZ,NG),A_MAT_I_J_K(NX,NY,NZ,NG),&
               A_MAT_IP1_J_K(NX,NY,NZ,NG),A_MAT_IM1_J_K(NX,NY,NZ,NG),&
               A_MAT_I_JP1_K(NX,NY,NZ,NG),A_MAT_I_JM1_K(NX,NY,NZ,NG),&
               A_MAT_I_J_KP1(NX,NY,NZ,NG),&
               A_MAT_I_J_KM1(NX,NY,NZ,NG),&
               RHS_SCAT(NX,NY,NZ),&
               T_NEW_TEMP_NG(NX,NY,NZ)) 

! Define matricies
	DO IG = 1,NG
	DO K = 2, NZ-1
	DO J = 2, NY-1 
	DO I = 2, NX-1 
            DIFF_I = 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I-1,J,K,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I+1,J,K,IG),0.0))/(DX**2) &
                   + 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I,J-1,K,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I,J+1,K,IG),0.0))/(DY**2) &
                   + 0.5*(max(KDIFF(I,J,K,IG)+KDIFF(I,J,K-1,IG),0.0)+ max(KDIFF(I,J,K,IG)+KDIFF(I,J,K+1,IG),0.0))/(DZ**2) 

            A_MAT_I_J_K(I,J,K,IG) = DIFF_I + SIGMAT(I,J,K,IG) + (1./(DT*V_N(IG))) 
            A_MAT_IP1_J_K(I,J,K,IG) =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I+1,J,K,IG),0.0) /(DX**2) 
            A_MAT_IM1_J_K(I,J,K,IG) =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I-1,J,K,IG),0.0) /(DX**2)
            A_MAT_I_JP1_K(I,J,K,IG) =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J+1,K,IG),0.0) /(DY**2) 
            A_MAT_I_JM1_K(I,J,K,IG) =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J-1,K,IG),0.0) /(DY**2) 
            A_MAT_I_J_KP1(I,J,K,IG) =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J,K+1,IG),0.0) /(DZ**2)
            A_MAT_I_J_KM1(I,J,K,IG) =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J,K-1,IG),0.0) /(DZ**2)
            
!		A_MAT_I_J_KM1(I,J,K,IG) =  - 0.5*max(KDIFF(I,J,K,IG)+KDIFF(I,J,K-1,IG),0.0) /(DZ**2) 


!            DIFF_I = 0.5*(max(KDIFF(I,J)+KDIFF(I-1,J),0.0)+ max(K(I,J)+K(I+1,J),0.0))/(DX**2) &
!                   + 0.5*(max(KDIFF(I,J)+KDIFF(I,J-1),0.0)+ max(K(I,J)+K(I,J+1),0.0))/(DY**2) 

!            A_MAT_I_J(I,J) = DIFF_I + SIGMAT(I,J) + (1./(DT*V_N)) 
!            A_MAT_IP1_J(I,J) =  - 0.5*max(K(I,J)+K(I+1,J),0.0) /(DX**2) 
!            A_MAT_IM1_J(I,J) =  - 0.5*max(K(I,J)+K(I-1,J),0.0) /(DX**2)
!            A_MAT_I_JP1(I,J) =  - 0.5*max(K(I,J)+K(I,J+1),0.0) /(DY**2) 
!            A_MAT_I_JM1(I,J) =  - 0.5*max(K(I,J)+K(I,J-1),0.0) /(DY**2) 
  
	END DO
	END DO 
	END DO 
	END DO 


         T_OLD = T_INITIAL

         DO ITIME=1,NTIME

         T_NEW = T_OLD
         DO ITS_NG=1,NITS_SOLV_NG
            R_T_MAX_DIF_NG = 0.0

! Forward backward Gauss - Seidel...
        DO IIGG=0,ISWITCH_FBGS ! Switch on FBGS
        IG_START =1*(1-IIGG)  + (NG-1)*IIGG
        IG_FINISH=NG*(1-IIGG) + 2*IIGG
        IG_STEP  =1*(1-IIGG)  -  1*IIGG
        DO IG=IG_START, IG_FINISH, IG_STEP

           RHS_SCAT=0.0 
           DO JG=1,NG
              RHS_SCAT(:,:,:) =  RHS_SCAT(:,:,:) + SIGMA_S_OFF(:,:,:,IG,JG) * T_NEW(:,:,:,JG) 
           END DO
        T_NEW_TEMP_NG(:,:,:)=T_NEW(:,:,:,IG)

        R_T_MAX_DIF = 0.0

        DO ITS=1,NITS

        DO KK=1,0,-1
	DO K = 2*KK + (NZ-1)*(1-KK), 2*(KK-1) + (NZ-1)*KK , 1*KK - 1*(KK-1)  
        DO JJ=1,0,-1
	DO J = 2*JJ + (NY-1)*(1-JJ), 2*(JJ-1) + (NY-1)*JJ , 1*JJ - 1*(JJ-1)  
        DO II=1,0,-1
	DO I = 2*II + (NX-1)*(1-II), 2*(II-1) + (NX-1)*II , 1*II - 1*(II-1)  

            A_I_J_K = A_MAT_I_J_K(I,J,K,IG)
            A_IP1_J_K =  A_MAT_IP1_J_K(I,J,K,IG) 
            A_IM1_J_K =  A_MAT_IM1_J_K(I,J,K,IG)
            A_I_JP1_K =  A_MAT_I_JP1_K(I,J,K,IG)
            A_I_JM1_K =  A_MAT_I_JM1_K(I,J,K,IG)

            A_I_J_KP1 =  A_MAT_I_J_KP1(I,J,K,IG)
            A_I_J_KM1=A_MAT_I_J_KM1(I,J,K,IG)

            R_T_NEW = ( - A_IP1_J_K*T_NEW(I+1,J,K,IG) - A_IM1_J_K*T_NEW(I-1,J,K,IG) - A_I_JP1_K*T_NEW(I,J+1,K,IG) &
                  - A_I_JM1_K*T_NEW(I,J-1,K,IG) - A_I_J_KP1*T_NEW(I,J,K+1,IG) - A_I_J_KM1*T_NEW(I,J,K-1,IG) &
                  + S(I,J,K,IG) + RHS_SCAT(I,J,K) + (1./(DT*V_N(IG)))*T_OLD(I,J,K,IG) ) / A_I_J_K
            R_T_MAX_DIF = MAX(ABS(R_T_NEW-T_NEW(I,J,K,IG) ), R_T_MAX_DIF)
            T_NEW(I,J,K,IG) = RELAX*R_T_NEW + (1.-RELAX) * T_NEW(I,J,K,IG)
							 ! Relax the soln.       
	END DO
	END DO 
	END DO
	END DO 
	END DO 
	END DO 
        IF(R_T_MAX_DIF.LT.ERROR_SOLV) EXIT
	!T_TOTAL(:,:,:,:,ITS)=T_NEW(:,:,:,:)
        END DO ! DO ITS=1,NITS
        T_NEW(:,:,:,IG) = RELAX_NG*T_NEW(:,:,:,IG) + (1.-RELAX_NG) * T_NEW_TEMP_NG(:,:,:) 
        R_T_MAX_DIF_NG = MAX( R_T_MAX_DIF_NG,  MAXVAL(ABS(T_NEW_TEMP_NG(:,:,:)-T_NEW(:,:,:,IG) ))  )
	END DO ! DO IG=IG_START, IG_FINISH, IG_STEP
	END DO ! DO IIGG=0,ISWITCH_FBGS ! Switch on FBGS
	!T_TOTAL(:,:,:,:,ITS_NG)=T_NEW(:,:,:,:)
!          PRINT *,'ITS,R_T_MAX_DIF,ERROR_SOLV:',ITS,R_T_MAX_DIF,ERROR_SOLV
         IF(R_T_MAX_DIF_NG.LT.ERROR_SOLV_NG) EXIT
         END DO ! ITS_G=1,NITS_SOLV_NG

         T_OLD=T_NEW ! Prepare for next time step. 

         END DO ! DO ITIME=1,NTIME

	  END SUBROUTINE TIME_STEPING_DIFFUSION_CALC
! 
! 
! 
! 
! *********************************************************************************************
! *********************************************************************************************
! **************************************** ANN SUBROUTINES ************************************
! *********************************************************************************************
! *********************************************************************************************


!
!
       SUBROUTINE NEURAL_NET_BATCHES(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,NONODS,&
                  NEXAMPLE,NBATCH,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY,NITS,&
                  ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER,&
                  IDEAL_MIDLE,NDEAL_MIDLE )
! Train (IF nits>0) OR SOLVE FOR THE NEURON VALUES OF AN ANN. 
! NITS= no of iterations of the ANN
! ERROR_TOLER = tolerance of the ANN. 
! This sub forms the neural network for a set of NEXAMPLE problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! If NOREXP=0 find exponent of OUTPUT NEURONS otherwise dont.

! THE ITERATION PARAMETERS: ALPHA,NTEST,TEST_BELOW,TEST_ABOVE,ANNEAL_FRAC
! Default values: NITS=10000,ERROR_TOLER=1.e-5, ALPHA=0.001,ANNEAL_FRAC=0.999
! ALPHA= initial value of the max value to add to the weights.
! ANNEAL_FRAC=what fraction to reduce ALPHA by if we are not converging well
! WEIGHT_DECAY is the weight magnitude penalty term
! 
       IMPLICIT NONE
       INTEGER, PARAMETER :: W_TOLER=0.5
       INTEGER, intent( inout ) :: NBATCH
       INTEGER, intent( in ) :: NCOLM,NONODS,NLAYER,NLAY_IN,NLAY_OUT,NEXAMPLE,nits,NDEAL_MIDLE
       INTEGER, intent( in ) :: NLAY(NLAYER),NOREXP
       REAL, intent( in ) :: IDEAL_INPUT(NLAY_IN,NEXAMPLE),IDEAL_OUTPUT(NLAY_OUT,NEXAMPLE),&
                             IDEAL_MIDLE(NDEAL_MIDLE,NEXAMPLE),WEIGHT_DECAY
       REAL, intent( inout ) :: WEIGHT(NCOLM)
       REAL, intent( inout )  :: NEURAL(NONODS,NEXAMPLE),FF
       REAL, intent( out )  :: ERROR_TOLER
       REAL, intent( inout ) :: ALPHA
       REAL, intent( in ) :: ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER
! local variables
       INTEGER :: ITS,NPRINT,IGOT_BETTER,I,J,IEXAMPLE,ITS2,NITS2
       REAL :: FF_KEEP, ACC_TOLER, RAN_NO, ff_best_q
       LOGICAL :: SWITCH, only_best
       REAL, ALLOCATABLE :: WEIGHT_KEEP(:),RAN_WEIGHT(:)
       REAL, ALLOCATABLE :: NEURAL2(:,:)
       REAL, ALLOCATABLE :: IDEAL_INPUT2(:,:),IDEAL_OUTPUT2(:,:),IDEAL_MIDLE2(:,:)
       LOGICAL, ALLOCATABLE :: IN_EXAMPLE(:)
       INTEGER, ALLOCATABLE :: BATCH_EXAMPLE(:)

       ALLOCATE(WEIGHT_KEEP(NCOLM),RAN_WEIGHT(NCOLM))
       ALLOCATE(IN_EXAMPLE(NEXAMPLE))


       only_best=.true. ! only keep a batch if it improves the iteration

          IF((NBATCH.NE.NEXAMPLE).AND.(NBATCH.NE.0)) THEN


       NITS2=INT(SQRT(REAL(NITS)))
       DO ITS2=1,NITS2
            PRINT *,'ITS2,INT(SQRT(REAL(NITS)))=',ITS2,NITS2
       ALLOCATE(NEURAL2(NONODS,NBATCH))
       ALLOCATE(IDEAL_INPUT2(NLAY_IN,NBATCH),IDEAL_OUTPUT2(NLAY_OUT,NBATCH),IDEAL_MIDLE2(NDEAL_MIDLE,NBATCH))
       ALLOCATE(BATCH_EXAMPLE(NBATCH))
!          IF(.true.) THEN
             if(.false.) then
             DO I=1,NBATCH
                IEXAMPLE=i
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             else
             DO I=1,NBATCH
                DO J=1,10000
!                    print *,'j=',j
                   CALL RANDOM_NUMBER(RAN_NO)
                   IEXAMPLE = MIN(NEXAMPLE,MAX(1, INT(RAN_NO*REAL(NEXAMPLE)+1.0)))
                   IF(.NOT.IN_EXAMPLE(IEXAMPLE)) THEN
                      IN_EXAMPLE(IEXAMPLE) = .TRUE.
                      EXIT
                   ENDIF
                END DO
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             endif
       if(only_best) then
       CALL NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,&
            NONODS,NEXAMPLE,NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,&
            WEIGHT_DECAY,1,ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,&
            ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
       ff_best_q=ff
       weight_keep = weight
       endif

       CALL NEURAL_NET(NEURAL2,FF,IDEAL_INPUT2,IDEAL_OUTPUT2,WEIGHT,NCOLM,NONODS,&
            NBATCH,NBATCH,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY,&
            INT(SQRT(REAL(NITS))),ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,&
            ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
             NEURAL(:,BATCH_EXAMPLE(:))=NEURAL2(:,:) ! Store the results
             IN_EXAMPLE(BATCH_EXAMPLE(:))=.FALSE. ! set everything back to .false.

       if(only_best) then
       CALL NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,NONODS,NEXAMPLE,&
            NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY,&
            1,ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
           print *,'ff,ff_best_q,nbatch:',ff,ff_best_q,nbatch
       if(ff.gt.ff_best_q) weight = weight_keep
       if(ff.gt.ff_best_q) print *,'****reset the weights'
          if(ff.gt.ff_best_q) then
             nbatch=nbatch+1
             if(nbatch.gt.nexample/2) nbatch=nexample
          else
             nbatch=nbatch-1
             if(nbatch.gt.nexample/2) nbatch=nexample/2
          endif
          print *,'new nbatch=',nbatch
       endif

       DEALLOCATE(NEURAL2, IDEAL_INPUT2, IDEAL_OUTPUT2,IDEAL_MIDLE2, BATCH_EXAMPLE )

       END DO ! ITS2

          ELSE

       CALL NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,&
            NONODS,NEXAMPLE,NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,NOREXP,WEIGHT_DECAY, &
            NITS,ERROR_TOLER, ALPHA, ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER, IDEAL_MIDLE,NDEAL_MIDLE)
          ENDIF

       RETURN
       END SUBROUTINE NEURAL_NET_BATCHES
!
!
!
!
!
!
       SUBROUTINE NEURAL_NET(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,WEIGHT,NCOLM,NONODS,&
                             &NEXAMPLE,NBATCH,NLAY,NLAYER,NLAY_IN,NLAY_OUT,&
                             &NOREXP,WEIGHT_DECAY,&
                             &NITS,ERROR_TOLER,ALPHA, ANNEAL_FRAC_SMALLER,&
                             &ANNEAL_FRAC_BIGGER,IDEAL_MIDLE,NDEAL_MIDLE)
! Train (IF nits>0) OR SOLVE FOR THE NEURON VALUES OF AN ANN. 
! NITS= no of iterations of the ANN
! ERROR_TOLER = tolerance of the ANN. 
! This sub forms the neural network for a set of NEXAMPLE problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! If NOREXP=0 find exponent of OUTPUT NEURONS otherwise dont.

! THE ITERATION PARAMETERS: ALPHA,NTEST,TEST_BELOW,TEST_ABOVE,ANNEAL_FRAC
! Default values: NITS=10000,ERROR_TOLER=1.e-5, ALPHA=0.001,ANNEAL_FRAC=0.999
! ALPHA= initial value of the max value to add to the weights.
! ANNEAL_FRAC=what fraction to reduce ALPHA by if we are not converging well
! WEIGHT_DECAY is the weight magnitude penalty term
! 
       IMPLICIT NONE
       INTEGER, PARAMETER :: W_TOLER=0.5
       INTEGER, intent( in ) :: NCOLM,NONODS,NLAYER,NLAY_IN,NLAY_OUT,NEXAMPLE,NBATCH,nits,NDEAL_MIDLE
       INTEGER, intent( in ) :: NLAY(NLAYER),NOREXP
       REAL, intent( in ) :: IDEAL_INPUT(NLAY_IN,NEXAMPLE),IDEAL_OUTPUT(NLAY_OUT,NEXAMPLE),IDEAL_MIDLE(NDEAL_MIDLE,NEXAMPLE)
       REAL, intent( in ) :: WEIGHT_DECAY
       REAL, intent( inout ) :: WEIGHT(NCOLM)
       REAL, intent( inout )  :: NEURAL(NONODS,NEXAMPLE),FF
       REAL, intent( out )  :: ERROR_TOLER
       REAL, intent( inout ) :: ALPHA
       REAL, intent( in ) :: ANNEAL_FRAC_SMALLER,ANNEAL_FRAC_BIGGER
! local variables
       INTEGER :: ITS,NPRINT,IGOT_BETTER,I,J,IEXAMPLE
       REAL :: FF_KEEP, ACC_TOLER, RAN_NO
       LOGICAL :: SWITCH
       REAL, ALLOCATABLE :: WEIGHT_KEEP(:),RAN_WEIGHT(:)
       REAL, ALLOCATABLE :: NEURAL2(:,:)
       REAL, ALLOCATABLE :: IDEAL_INPUT2(:,:),IDEAL_OUTPUT2(:,:),IDEAL_MIDLE2(:,:)
       LOGICAL, ALLOCATABLE :: IN_EXAMPLE(:)
       INTEGER, ALLOCATABLE :: BATCH_EXAMPLE(:)

       ALLOCATE(WEIGHT_KEEP(NCOLM),RAN_WEIGHT(NCOLM))
       ALLOCATE(NEURAL2(NONODS,NBATCH))
       ALLOCATE(IN_EXAMPLE(NEXAMPLE))
       ALLOCATE(IDEAL_INPUT2(NLAY_IN,NBATCH),IDEAL_OUTPUT2(NLAY_OUT,NBATCH),IDEAL_MIDLE2(NDEAL_MIDLE,NBATCH))
       ALLOCATE(BATCH_EXAMPLE(NBATCH))

       IN_EXAMPLE=.FALSE.

       CALL NEURAL_EXAMPLE(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,&
            WEIGHT,NCOLM,NONODS,NEXAMPLE,NLAY,NLAYER,NLAY_IN,&
            NLAY_OUT,NOREXP,WEIGHT_DECAY, IDEAL_MIDLE,NDEAL_MIDLE)
!            print *,'original ff=',ff

       ACC_TOLER=FF
       
       IGOT_BETTER=0
       SWITCH=.FALSE.
       DO ITS=1,NITS ! Train the ANN if NITS>0
          WEIGHT_KEEP=WEIGHT
          IF(SWITCH) THEN
             WEIGHT = WEIGHT - (RAN_WEIGHT-0.5)*ALPHA
          ELSE
             CALL RANDOM_NUMBER(RAN_WEIGHT)
!            alpha=1.e-5
             WEIGHT = WEIGHT + (RAN_WEIGHT-0.5)*ALPHA
          ENDIF
!          print *,'RAN_WEIGHT-0.5:',RAN_WEIGHT-0.5
!          stop 11

          FF_KEEP = FF
          IF((NBATCH.NE.NEXAMPLE).AND.(NBATCH.NE.0)) THEN
!          IF(.true.) THEN
             if(.false.) then
             DO I=1,NBATCH
                IEXAMPLE=i
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             else
             DO I=1,NBATCH
                DO J=1,10000
!                    print *,'j=',j
                   CALL RANDOM_NUMBER(RAN_NO)
                   IEXAMPLE = MIN(NEXAMPLE,MAX(1, INT(RAN_NO*REAL(NEXAMPLE)+1.0)))
                   IF(.NOT.IN_EXAMPLE(IEXAMPLE)) THEN
                      IN_EXAMPLE(IEXAMPLE) = .TRUE.
                      EXIT
                   ENDIF
                END DO
                BATCH_EXAMPLE(I)=IEXAMPLE
                IDEAL_INPUT2(1:NLAY_IN,I) = IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
                IDEAL_OUTPUT2(1:NLAY_OUT,I) = IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE) 
                IDEAL_MIDLE2(1:NDEAL_MIDLE,I) = IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE) 
                NEURAL2(:,I)=NEURAL(:,IEXAMPLE)
             END DO
             endif
!             print *,'Batch_example:',batch_example
!             print *,'NEURAL2:',NEURAL2
!             print *,'IDEAL_INPUT2:',IDEAL_INPUT2
!             print *,'IDEAL_OUTPUT2:',IDEAL_OUTPUT2
             CALL NEURAL_EXAMPLE(NEURAL2,FF,IDEAL_INPUT2,IDEAL_OUTPUT2,&
                  WEIGHT,NCOLM,NONODS,NBATCH,NLAY,NLAYER,NLAY_IN,&
                  NLAY_OUT,NOREXP,WEIGHT_DECAY, IDEAL_MIDLE2,NDEAL_MIDLE)
             NEURAL(:,BATCH_EXAMPLE(:))=NEURAL2(:,:) ! Store the results
             IN_EXAMPLE(BATCH_EXAMPLE(:))=.FALSE. ! set everything back to .false.
!             ff_keep = 0.9*ff_keep + 0.1*ff
!             print *,'ff=',ff
          ELSE
             CALL NEURAL_EXAMPLE(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,&
                  WEIGHT,NCOLM,NONODS,NEXAMPLE,NLAY,NLAYER,NLAY_IN,&
                  NLAY_OUT,NOREXP,WEIGHT_DECAY, IDEAL_MIDLE,NDEAL_MIDLE)
          ENDIF

          ACC_TOLER = (1.-W_TOLER) * ACC_TOLER + W_TOLER * ABS(FF-FF_KEEP) ! Make the tolerance change slowly.

          IF(ACC_TOLER<ERROR_TOLER) CYCLE
!          print *,'its,ff,ff_keep:',its,ff,ff_keep

          IF(FF<FF_KEEP) THEN
             IF(SWITCH) ALPHA=ALPHA/ANNEAL_FRAC_BIGGER  ! make bigger
             SWITCH=.FALSE.
             IGOT_BETTER=IGOT_BETTER+1
!              print *,'****'
          ELSE
             IF(SWITCH) ALPHA=ALPHA*ANNEAL_FRAC_SMALLER  ! make smaller
             SWITCH=.NOT.SWITCH
             FF=FF_KEEP
             WEIGHT=WEIGHT_KEEP
          ENDIF
!          print *,'alpha=',alpha

!          NPRINT=100000
!          NPRINT=100000
!          NPRINT=1000
          NPRINT=100
          IF(MOD(ITS,NPRINT)==0) THEN ! REAL(IGOT_BETTER)/REAL(NPRINT)=0.66 WHEN ALPHA=SMALL NO
             print *,'its,ff,ALPHA:',its,ff,ALPHA,REAL(IGOT_BETTER)/REAL(NPRINT)
             IGOT_BETTER=0
          ENDIF

       END DO ! DO ITS=1,NITS
!
       RETURN
       END SUBROUTINE NEURAL_NET
!
!
!
!
       SUBROUTINE NEURAL_EXAMPLE(NEURAL,FF,IDEAL_INPUT,IDEAL_OUTPUT,&
                  WEIGHT,NCOLM,NONODS,NEXAMPLE,NLAY,NLAYER,NLAY_IN,NLAY_OUT,&
                  NOREXP,WEIGHT_DECAY, IDEAL_MIDLE,NDEAL_MIDLE)
! This sub forms the neural network for a set of NEXAMPL problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! If NOREXP=0 find exponent of OUTPUT NEURONS otherwise dont.
! WEIGHT_DECAY is the weight magnitude penalty term
! 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NCOLM,NONODS,NLAYER,NLAY_IN,NLAY_OUT,NEXAMPLE,NDEAL_MIDLE
       INTEGER, intent( in ) :: NLAY(NLAYER),NOREXP
       REAL, intent( in ) :: IDEAL_INPUT(NLAY_IN,NEXAMPLE),IDEAL_OUTPUT(NLAY_OUT,NEXAMPLE),WEIGHT(NCOLM),&
                             WEIGHT_DECAY, IDEAL_MIDLE(NDEAL_MIDLE,NEXAMPLE)
       REAL, intent( inout )  :: NEURAL(NONODS,NEXAMPLE),FF
! Local varibales...
       INTEGER :: IEXAMPLE,ISTART
       REAL, ALLOCATABLE :: W_EXP(:)

       ALLOCATE(W_EXP(NLAYER)) 

       W_EXP(:)=0.0
       IF(NOREXP==0) W_EXP(NLAYER)=1.0 ! no exponential in final layer...
       IF(NOREXP==0) W_EXP(int(NLAYER/2) +1)=1.0 ! no exponential in final layer...

       ISTART = SUM(NLAY(1:NLAYER/2)) +1

       FF=0.0
       DO IEXAMPLE=1,NEXAMPLE
          NEURAL(1:NLAY_IN,IEXAMPLE)=IDEAL_INPUT(1:NLAY_IN,IEXAMPLE) 
          CALL GETNEUVALS_FAST(NEURAL(:,IEXAMPLE),WEIGHT,NONODS,NCOLM,W_EXP,NLAY,NLAYER)

          FF=FF + SUM( (NEURAL(NONODS-NLAY_OUT+1:NONODS,IEXAMPLE)-IDEAL_OUTPUT(1:NLAY_OUT,IEXAMPLE))**2 )
          FF=FF + 0.*SUM( (NEURAL(ISTART:ISTART+NDEAL_MIDLE-1,IEXAMPLE)-IDEAL_MIDLE(1:NDEAL_MIDLE,IEXAMPLE))**2 )
       END DO ! DO IEXAMPLE=1,NEXAMPLE
       FF=FF/REAL(NEXAMPLE)
       FF=FF+(WEIGHT_DECAY/REAL(NCOLM)) * SUM(WEIGHT(:)**2) 
!
       RETURN
       END SUBROUTINE NEURAL_EXAMPLE
!
!
!
!
       SUBROUTINE GETNEUVALS_FAST(NEUVAL,WEIGHT,NONODS,NCOLM,W_EXP,NLAY,NLAYER)
! This sub forms the neural network for a set of NEXAMPL problems
! The neural network neran values are in NEURAL and the functional used for training (output neuron data missmatch) is FF. 
! IDEAL_INPUT,IDEAL_OUTPUT are the ideal input and output neuron values for the example problems.
! WEIGHT are the weights of the neural network. 
! NLAY(ILAYER) = the number of neurons in layer ILAYER
! NLAYER = No of neuron layers including input and putout. 
! NLAY_IN,NLAY_OUT are the no of neurons in the input and output layers.
! NLAYER1,NLAFIN no of nodes in the 1st (input layer) and final (output layer)
! This sub calculates the neuron values. 
! If W_EXP=0.0 find exponent of NEURONS otherwise dont. 
       IMPLICIT NONE
       INTEGER, intent( in ) :: NONODS,NCOLM, NLAYER
       REAL, intent( in ) :: WEIGHT(NCOLM),W_EXP(NLAYER)
       REAL, intent( inout ) :: NEUVAL(NONODS)
       INTEGER , intent( in ):: NLAY(NLAYER)
! LOCAL VARIABLES...
       REAL :: SUMWEI_VAL
       INTEGER :: NOD,ILAY,ILAY1,ILAY2,NLAYACC_WEIT,NLAYACC_NOD,NLAYACC_NOD_PREV, IWEI_ADD
       
       NLAYACC_WEIT=0
       NLAYACC_NOD_PREV=0
       NLAYACC_NOD=NLAY(1)
       DO ILAY=2,NLAYER

          DO ILAY2=1,NLAY(ILAY) 
             SUMWEI_VAL=0.0
             IWEI_ADD = NLAYACC_WEIT + (ILAY2-1)*NLAY(ILAY-1)
             DO ILAY1=1,NLAY(ILAY-1)
                SUMWEI_VAL=SUMWEI_VAL+WEIGHT(IWEI_ADD + ILAY1)*NEUVAL(ILAY1 + NLAYACC_NOD_PREV)
             END DO
             NEUVAL(NLAYACC_NOD+ILAY2)= W_EXP(ILAY)*SUMWEI_VAL   +   (1.-W_EXP(ILAY))/(1.+EXP(-SUMWEI_VAL))
          END DO
          NLAYACC_NOD_PREV=NLAYACC_NOD_PREV+NLAY(ILAY-1)
          NLAYACC_NOD=NLAYACC_NOD+NLAY(ILAY)
          NLAYACC_WEIT=NLAYACC_WEIT + NLAY(ILAY) * NLAY(ILAY-1)

       END DO

       RETURN
       END SUBROUTINE GETNEUVALS_FAST
!
!

