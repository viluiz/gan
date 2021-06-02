

          PROGRAM MAIN
	  IMPLICIT NONE
          REAL PIE, TOLER
          PARAMETER(PIE=3.141592654, TOLER=1.E-10)
       real :: R0_HOME, R0_MOBILE, RAN_NO1, RAN_NO2
       INTEGER :: IRUNS,NRUNS, k
       CHARACTER CHARA*240
       CHARACTER str_trim*240
! Local variables...
! Set up problem ********************************

            
!            R0_HOME=0.2 !5.0 !0.2 ! No of people an infected person infects.
!            R0_MOBILE=20.0 !10.0 !20.0 !5.0 ! No of people an infected person infects. 

       NRUNS=1

       DO IRUNS=1,NRUNS

!         CALL SYSTEM('gzip run/group-output-time'//str_trim(1:k)//'.csv' ) ! factor of 12 reduction in disc space.
            str_trim=chara(IRUNS); k=INDEX(str_trim,' ')-1
            CALL SYSTEM('mkdir run'//str_trim(1:k) ) 

            CALL RANDOM_NUMBER(RAN_NO1)
            CALL RANDOM_NUMBER(RAN_NO2)
            R0_HOME = 0.2  + 19.8*RAN_NO1
            R0_MOBILE = 0.2  + 19.8*RAN_NO2
          print *,'iruns, R0_HOME,  R0_MOBILE:',iruns, R0_HOME,  R0_MOBILE

!            CALL SYSTEM('mkdir run'//str_trim(1:k) ) 

            open(27, file='run'//str_trim(1:k)//'/r0-2values.csv', status='replace')
            WRITE(27, *) '# R0_HOME, R0_MOBILE'
            WRITE(27, * ) R0_HOME, R0_MOBILE
            close(27) 

            CALL SYSTEM('cp diffusion3d_virus run'//str_trim(1:k)//'/diffusion3d_virus' ) 
            CALL SYSTEM('./run'//str_trim(1:k)//'/diffusion3d_virus' ) 
!            CALL SYSTEM('cd ./run1' ) 
!            CALL SYSTEM('chdir run1' ) 
!            CALL SYSTEM('cd run'//str_trim(1:k) ) 
!            CALL SYSTEM('.././diffusion3d_virus > /dev/null &')
!            CALL SYSTEM('./diffusion3d_virus')
!            CALL SYSTEM('more r0-2values.csv')
            CALL SYSTEM('pwd')
!            CALL SYSTEM('cd ..') 
       END DO 

       STOP
       END PROGRAM MAIN
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


!

