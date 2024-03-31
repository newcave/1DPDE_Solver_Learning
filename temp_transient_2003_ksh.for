!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Parabolic Partial Differenctial Equations Soliving 
! 1-Dimensional P-PDE 
! DIRICHLET BOUNDARY CONDITIONS
! Basic Schemes FTCS, DUFORT-FRANKEL, LAASONEN & CRANK-NICOLSON
! f77 to f90, Modified by KSH 2003.4.13 & 2024.04.26

!*******************************************************************
!*** MODULE BLOCKS
!*******************************************************************

MODULE G1
IMPLICIT NONE
INTEGER, SAVE :: JM 
END MODULE G1

MODULE G2
IMPLICIT NONE
SAVE
REAL (KIND=8) :: CONV, D 
END MODULE G2

	MODULE G3
	IMPLICIT NONE
	SAVE
	REAL (KIND=8) :: UC, US, NU, DLENGTH
	END MODULE G3

	MODULE G4
	IMPLICIT NONE
	SAVE
	INTEGER, PARAMETER :: JMAX=200, NMAX=700
 
!*******************************************************************
!*** PROGRAM MAIN 
!*******************************************************************

PROGRAM ODPPDE 

!*******************************************************************
!*** VARIABLE DECLARATION 변수 선언 
!*******************************************************************

USE G1
USE G2
USE G3
USE G4

IMPLICIT NONE

INTEGER :: J, OPT, INCASE, N, NM
INTEGER :: COUNT, ORDER, CHECK
REAL
REAL (KIND=8) :: DELY, DELT, TOTTIME, UP, PI
REAL (KIND=8) :: U(JMAX, NMAX), Y(JMAX), UANA(JMAX), U1(JMAX)
REAL (KIND=8) :: U1(JMAX) ,U2(JMAX), U3(JMAX), U4(JMAX), U5(JMAX)
REAL (KIND=8) :: U6(JMAX)
REAL (KIND=8) :: U1ANA(JMAX), U2ANA(JMAX), U3ANA(JMAX), U4ANA(JMAX)
REAL (KIND=8) :: U5ANA(JMAX), U6ANA(JMAX)


!*******************************************************************
!*** INPUT DATA/ INITIALIZATION
!*******************************************************************

PI = DACOS(-1.D0)
DLENGTH = 1.D0 		! 총 연장 
TOTTIME = 0.5D0		! 총 연산시간
UC = 50.D0			! 초기 온도 설정 
US = 0.D0			! 표면 온도 
NU = 0.03D0			! 확산계수 [L^2/T] dissusivity

PRINT*, ' INPUT DELY :'
READ*, DELY			! 공간 스텝 		
PRINT*, 'INPUT DELT'
READ*, DELT   		! 시간 스텝 

!*******************************************************************
!*** COMPUTE MAXIMUM NUMBER OF SPATIAL GRID POINTS (JM)
!*** NUMBER OF TIME STEPS (NM)
!*******************************************************************

JM = IDNINT(DLENGTH/DBLE(DELY))+1
NM = IDNINT(DBLE(TOTTIME/DELT))+1

CONV = 1.0D-3

!*******************************************************************
!*** COMPUTE DIFFUSION NUMBER
!*******************************************************************

D = NU*DELT/DELY/DELY

!*******************************************************************
!*** INITIAL AND BOUNDARY CONDITION
!*******************************************************************

DIST: DO J = 1,JM

	Y (J) = DELY*DBLE(J-1)

END DO DIST

T: DO N = 1, NM
		S: DO J = 1, JM
	IF (J == 1 .AND. J == JM) THEN
		U(J, N) = US
		ELSE
			UI(J) = UC*DSIN(PI*Y(J)/DLENGTH)
		ENDIF 	
		END DO S
END DO T

DO
	PRINT*, 'FTCS EXPLICIT 	:1'
	PRINT*, 'DUFORT-FRANKEL	:2'
	PRINT*, 'LAASONEN		:3'
	PRINT*, 'CRANK-NICOLSON :4'	
	READ*, OPT
	IF (OPT<1 .OR. OPT> 4) THEN
		PRINT*, "WRONG KEY ENTERED!"
		PRINT*, f"RE-ENTER THE KEY: 1-4"
		PRINT*, "  "
	ELSE
		EXIT
	END IF
END DO

CHECK = 0
COUNT = 0

TOTAL: DO N = 2,NM
	PRINT*, NM, N
	
	
!*******************************************************************
!*** DETERMINE TIME PRINTING TIME SEQUENCE
!*******************************************************************
	
	IF (DABS (DELT* DBLE(N-1)-0.1D0) <= 1.D-7) THEN
		ORDER = 1
	ELSEIF (DABS (DELT*DBLE(N-1)-0.2D0) <= 1. D-7) THEN
		ORDER = 1
	ELSEIF (DABS (DELT*DBLE(N-1)-0.3D0) <= 1. D-7) THEN
		ORDER = 1
	ELSEIF (DABS (DELT*DBLE(N-1)-0.4D0) <= 1. D-7) THEN
		ORDER = 1
	ELSEIF (DABS (DELT*DBLE(N-1)-0.5D0) <= 1. D-7) THEN
		ORDER = 1
	ELSE
		ORDER = 0
	END IF

	IF (ORDER == 1) CALL ANALYTIC1(Y, DELT*DBLE(N-1), UANA) 

	SCHEME: SELECT CASE (OPT)

	CASE (1)
		CALL FTCSEX(N, U) 
	CASE (2)
		CALL D_FEX(CHECK, N, U)
	CASE (3)
		CALL LAAIM(N, U) 
	CASE (4)
		CALL C_NIM (N, U)

	END SELECT SCHEME

!*******************************************************************
!*** READ THE RESULTS 
!*******************************************************************

	IF (ORDER == 1) THEN

		COUNT = COUNT+1	

		COLLECT: DO J = 1,JM

			IF (COUNT == 1) THEN
				U1 (J) = U (J, N)
				ULANA (J) = UANA (J)
			ELSEIF (COUNT == 2) THEN
				U2 (J) = U (J, N)
				ULANA (J) = UANA (J)
			ELSEIF (COUNT == 3) THEN
				U3 (J) = U (J, N)
				USANA (J) = UANA (J)
			ELSEIF (COUNT == 4) THEN
				U4 (J) = U (J, N)
				U4ANA (J) = UANA (J)
			ELSEIF (COUNT == 5) THEN
				US (J) = U (J, N)
				USANA (J) = UANA (J)
			ENDIF
		END DO COLLECT
	ENDIF

END DO TOTAL

FILES: SELECT CASE (OPT)
	CASE (1) 
		OPEN (UNIT=9, FILE=' FTCSA2.DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
		OPEN  (UNIT=29, FILE=' FTCSERRA2 DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
	CASE (2)
		OPEN (UNIT=9, FILE='D_FA2.DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
		OPEN  (UNIT=29, FILE='D_FA2ERRA2.DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
	CASE (3)
		OPEN (UNIT=9, FILE='LASSA2.DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
		OPEN  (UNIT=29, FILE='LAAERRA2.DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
	CASE (4)
		OPEN (UNIT=9, FILE='C_NA2.DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
		OPEN  (UNIT=29, FILE='C_NERRA2.DAT', STATUS='REPLACE', &
			ACTION= 'WRITE', IOSTAT=INCASE)
END SELECT FILES 

	OPEN (UNIT=19, FILE=' ANALYTICA2.DAT', STATUS=' REPLACE', &
			ACTION= 'WRITE', IOSTAT= INCASE)


!*******************************************************************
!*** WRITE OUTPUT TO FILES
!*******************************************************************

WRITE (9, *)
WRITE (9,10)
WRITE (9, *)
WRITE (19, *)
WRITE (19,10)
WRITE (19, *)

PRINTOUT: DO J = JM,1,-2

	IF (J == 1 .AND. J == JM) THEN
		UP = US
	ELSE
		UP = UI(J)
	ENDIF

	WRITE (9,20) Y(J), UP, U1(J), U2(J), U3(J), U4(J), U5(J)
	WRITE (19,20) Y(J), UP, UIANA(J), U2ANA(J), U3ANA (J), &
		U4ANA(J), U5ANA(J)
	WRITE (29, 30) Y(J), U1(J), UIANA(J), US(J), USANA(J)
END DO PRINTOUT

10 FORMAT (1X, 'Y', 3X, 'T =0.0, 3X, 'T=0.1', 3X, 'T=0.2, 3X, &
			'T= 0.3', 3X, 'T=0.4', 3X, 'T=0.5')
20 FORMAT (1X, F6.3, 3X, F6. 3, 5(3X, F6.3))
30 FORMAT (1X, F6.3,2 (3X, E16.9))
CLOSE (9)
CLOSE (19)

END PROGRAM ODPPDE

!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!*******************************************************************
!*** (1) THE FORWARD TIME/CENTRAL SPACE METHOD 
!*******************************************************************

SUBROUTINE FTCSEX (N, U)
USE G1
USE G2
USE G4

IMPLICIT NONE
INTEGER :: J, N

REAL (KIND=8) :: ERROR
REAL (KIND=8) :: U (JMAX, NMAX)

!*** THE FORWARD TIME/CENTRAL SPACE METHOD EQN

S: DO J = 2,JM-1
U(J,N) = U(J,N-1)+D*(U(J+1,N-1)-2.D0U(J, N-1)+U(J-1,N-1))
END DO S
RETURN
END SUBROUTINE FTCSEX


!*******************************************************************
!*** (2) THE DUFORT-FRANKEL METHOD
!*******************************************************************

SUBROUTINE D_FEX (CHECK, N, U)
USE G1
USE G2
USE G4
IMPLICIT NONE
INTEGER :: J, CHECK, N
REAL (KIND=8) :: ERROR
REAL (KIND=8) :: U(JMAX, NMAX), UNP1(JMAX), UNM1(JMAX)

!*** DUFORT-FRANKEL METHOD EQN

CHECK = CHECK+1

IF (CHECK == 1 .AND. D<1.D0) THEN
	CALL ETCSEX (N, U)
	GOTO 999

ELSEIF (CHECK == 1 .AND. D>1.D0) THEN
	CALL C NIM (N, U) 
	GOTO 999
END IF

NEW: DO J = 2, JM-1

	U(J, N) = ((1.D0-2.D0*D)*U(J,N-2)+2.D0*D*(U(J+1,N-1)+*(J-1,N-1)))/1.D0+2.D0*D)
END DO NEW


999 CONTINUE

RETURN
END SUBROUTINE D_FEX



!*******************************************************************
!*** (3) LAASONEN METHOD 
!*******************************************************************

SUBROUTINE LAAIM (N, U) 

USE G1
USE G2
USE G4

IMPLICIT NONE

INTEGER :: J, N

REAL (KIND=8) :: ERROR
REAL (KIND=8) :: U(JMAX, NMAX)
REAL (KIND=8), DIMENSION (JMAX) :: AA,BB,CC,DD

!*** LAASONEN EQN

DO J = 2, JM-1
	AA (J) = D
	BB (J) = -(2.D0*D+1.D0)
	CC (J) = D
    DD (J) = -U(J,N-1)
END DO

CALL TRID(N, AA, BB, CC DD, U)

RETURN
END SUBROUTINE LAAIM


!*******************************************************************
!*** SOLVER FOR TRIDIAGONAL MATRIX 
!*******************************************************************

CALL TRID(N, AA, BB, CC, DD, U)

USE G1
USE G4

IMPLICIT NONE

INTEGER :: J, N

REAL (KIND=8) , DIMENSION(JMAX) :: H, G, AA, BB, CC, DD
REAL (KIND=8) :: U (JMAX, NMAX)
REAL (KIND=8) :: ERROR 

G(1) = U(1, N-1)

HANDG: DO J = 2, M-1

	H(J)＝ CC(J) / BB(J) - AA(J)*H(J-1))
	
	G(J) = (DD(J)-AA(J)*G(J- 1)) / (BB(J)-AA(J)*H(J- 1))

END DO HANDG

CALU: DO J = JM-1,2,-1
		
	U(J, N) = -H(J)*U(J+1,N)+G(J)

END DO CALU


!*******************************************************************
!*** (4) CRANK-NICOLSON
!*******************************************************************

SUBROUTINE C_NIM (N, U) 

USE G1
USE G2
USE G4

IMPLICIT NONE

INTEGER :: J, N

REAL (KIND=8) :: ERROR
REAL (KIND=8) :: U(JMAX, NMAX)
REAL (KIND=8), DIMENSION (JMAX) :: AA, BB, CC, DD

!*** THE FORWARD TIME/CENTRAL SPACE METHOD EQN

DO J = 2,JM-1

	AA (J) = -D
	BB (J) = 2.D0*(D+1.D0)
	CC (J) = -D
	DD (J) = D*(U(J+1,N-1)-2.D0*U(J,N-1)+U(J-1,N-1))+2.D0*U(J,N-1)

END DO

CALLL TRID (N, AA, BB, CC, DD, U) 
RETURN
END SUBROUTINE C_NIM


!*******************************************************************
!*** (5) ANALYTICAL SOLUTION
!*******************************************************************

SUBROUTINE ANALYTIC1(Y, TT, T)

USE G3
USE G4

IMPLICIT NONE

INTEGER :: J
REAL (KIND=8) :: Y(JMAX), T(JMAX)
REAL (KIND=8) :: TT, L, PI

L = DLENGTH
PI = ACOS(-1.D0)

NODE: DO J = 1,JM
	T (J) = UC* DEXP(-NU*PI*PI/L/L*TT)*DSIN(PI*Y(J)/L)
END DO NODE

RETURN
END SUBROUTINE ANALYTIC1