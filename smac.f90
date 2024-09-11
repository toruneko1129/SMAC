!**********************************************************************
!     UNSTEADY FLOW IN CUBIC CAVITY   --- SMAC METHOD                 
!**********************************************************************
PROGRAM cubic_cavity_smac
IMPLICIT NONE

INTEGER, PARAMETER :: NX=300, NY=300
REAL :: U(NX,NY), V(NX,NY), P(NX,NY), Q(NX,NY), D(NX,NY), R(NX,NY)
REAL :: DT, RE, AFCT, DX, DY, XD, YD, XDH, YDH, TD, R1
INTEGER :: MX, MY, I21, J21, LM, KM, I, J, L, K
INTEGER :: I20, I19, J20, J19
REAL :: A1, A2, A3, U1, V2, UA, UB, VA, VB, DIVV, G2, G1, ULI, DIV2
REAL :: V1, UN, UV, VN, VV
  
! INPUT PARAMETERS
PRINT *, 'INPUT NUMBER OF MESH(<150) (20,20)'
MX = 130
MY = 130
I21 = MX + 1
J21 = MY + 1

PRINT *, 'INPUT TIME INCREMENT DT  (0.01)'
DT = 0.001

PRINT *, 'INPUT REYNOLDS NUMBER (40)'
READ *, RE

PRINT *, 'INPUT NUMBER OF TIME STEP (300)'
LM = 100 / DT

PRINT *, 'INPUT MAXIMUM NUMBER OF ITERATIONS FOR POISSON EQ. (50)'
KM = 100

! CALCULATE CONSTANTS
AFCT = 0.5
DX   = 1.0 / REAL(I21 - 3)
DY   = 1.0 / REAL(J21 - 3)
XD   = 1.0 / DX
YD   = 1.0 / DY
XDH  = 0.5 * XD
YDH  = 0.5 * YD
TD   = 1.0 / DT
R1   = 1.0 / RE
I20  = I21 - 1
I19  = I21 - 2
J20  = J21 - 1
J19  = J21 - 2
A1   = 0.5 * DY * DY / (DX * DX + DY * DY)
A2   = 0.5 * DX * DX / (DX * DX + DY * DY)
A3   = 0.5 * DY * DY / (1.0 + DY * DY / (DX * DX))

! INITIAL CONDITION
DO J = 1, J21
    DO I = 1, I21
        U(I, J) = 0.0
        V(I, J) = 0.0
        P(I, J) = 0.0
        Q(I, J) = 0.0
    END DO
END DO

! TIME MARCHING
DO L = 1, LM

    ! BOUNDARY CONDITION
    DO I = 1, I21
    	U(I, J20) = 1.0
		U(I, J21) = 2.0 - U(I, J19)
        V(I, J20) = 0.0
		V(I, J21) = -V(I, J19)
		U(I, 1)   = -U(I, 3)
		U(I, 2)   = 0.0
        V(I, 1)   = -V(I, 3)
        V(I, 2)   = 0.0
    END DO

    DO J = 1, J21
		U(I20, J) = 0.0
        U(I21, J) = -U(I19, J)
        V(I20, J) = 0.0
		V(I21, J) = -V(I19, J)
        U(1, J)   = -U(3, J)
        U(2, J)   = 0.0
        V(1, J)   = -V(3, J)
		V(2, J)   = 0.0
    END DO

	DO J = 2, J19
    	P(1,J)   = P(2,J)
    	P(I20,J) = P(I19,J)
	END DO

	DO I = 2, I19
    	P(I,1)   = P(I,2)
    	P(I,J20) = P(I,J19)
	END DO

	! TIME INTEGRATION FOR NAVIER-STOKES EQUATION
	DO J = 2, J19
    	DO I = 3, I19
        	V1 = 0.25 * (V(I,J) + V(I,J+1) + V(I-1,J+1) + V(I-1,J))
        	UN = U(I,J) * (U(I+1,J) - U(I-1,J)) * XDH &
            	+ V1 * (U(I,J+1) - U(I,J-1)) * YDH
        	UV = (U(I+1,J) - 2.0 * U(I,J) + U(I-1,J)) * XD * XD &
            	+ (U(I,J+1) - 2.0 * U(I,J) + U(I,J-1)) * YD * YD
        	U(I,J) = U(I,J) + DT * (-UN - (P(I,J) - P(I-1,J)) * XD + R1 * UV)
    	END DO
	END DO

	DO J = 3, J19
    	DO I = 2, I19
        	U1 = 0.25 * (U(I,J) + U(I+1,J) + U(I+1,J-1) + U(I,J-1))
        	VN = U1 * (V(I+1,J) - V(I-1,J)) * XDH &
            	+ V(I,J) * (V(I,J+1) - V(I,J-1)) * YDH
        	VV = (V(I+1,J) - 2.0 * V(I,J) + V(I-1,J)) * XD * XD &
            	+ (V(I,J+1) - 2.0 * V(I,J) + V(I,J-1)) * YD * YD
        	V(I,J) = V(I,J) + DT * (-VN - (P(I,J) - P(I,J-1)) * YD + R1 * VV)
    	END DO
	END DO

	! BOUNDARY CONDITION
    DO I = 1, I21
    	U(I, J20) = 1.0
		U(I, J21) = 2.0 - U(I, J19)
        V(I, J20) = 0.0
		V(I, J21) = V(I, J19)
		U(I, 1)   = -U(I, 3)
		U(I, 2)   = 0.0
        V(I, 1)   = V(I, 3)
        V(I, 2)   = 0.0
    END DO

    DO J = 1, J21
		U(I20, J) = 0.0
        U(I21, J) = U(I19, J)
        V(I20, J) = 0.0
		V(I21, J) = -V(I19, J)
        U(1, J)   = U(3, J)
        U(2, J)   = 0.0
        V(1, J)   = -V(3, J)
		V(2, J)   = 0.0
    END DO

	DO J = 2, J19
    	P(1,J)   = P(2,J)
    	P(I20,J) = P(I19,J)
	END DO

	DO I = 2, I19
    	P(I,1)   = P(I,2)
    	P(I,J20) = P(I,J19)
	END DO

     ! CALCULATE RIGHT HAND SIDE OF POISSON EQUATION
     DIVV = 0.0
     DO J = 2, J19
        DO I = 2, I19
           U1 = (U(I+1, J) - U(I, J)) * XD
           V2 = (V(I, J+1) - V(I, J)) * YD
           D(I, J) = U1 + V2
           DIVV = DIVV + ABS(U1 + V2)
           UA = 0.25 * (U(I, J) + U(I+1, J) + U(I+1, J+1) + U(I, J+1))
           UB = 0.25 * (U(I, J) + U(I+1, J) + U(I+1, J-1) + U(I, J-1))
           VA = 0.25 * (V(I, J) + V(I, J+1) + V(I+1, J+1) + V(I+1, J))
           VB = 0.25 * (V(I, J) + V(I, J+1) + V(I-1, J+1) + V(I-1, J))
           R(I, J) = -U1 * U1 - 2.0 * (UA - UB) * (VA - VB) * XD * YD - V2 * V2 + TD * (U1 + V2)
        END DO
     END DO

     ! SOLVING POISSON EQUATION FOR PRESSURE
     DO K = 1, KM
        G2 = 0.0
        ! PRESSURE BOUNDARY CONDITION
        DO J = 1, J21
           P(1, J) = P(2, J)
           P(I20, J) = P(I19, J)
        END DO
        DO I = 1, I21
           P(I, 1) = P(I, 2)
           P(I, J20) = P(I, J19)
        END DO

        ! SOR METHOD
        DO J = 2, J19
           DO I = 2, I19
              ULI = A1 * (P(I+1, J) + P(I-1, J)) + A2 * (P(I, J+1) + P(I, J-1)) - A3 * R(I, J) - P(I, J)
              G2 = G2 + ULI * ULI
              P(I, J) = ULI + P(I, J)
           END DO
        END DO
        IF (G2 <= 0.00001) EXIT
     END DO

	! *** SOLVING POISSON EQUATION FOR POTENTIAL ***
	DO K = 1, KM
    	G1 = 0.0
    	! BOUNDARY CONDITION FOR Q
    	DO J = 1, J21
        	Q(1,J)   = Q(2,J)
        	Q(I20,J) = Q(I19,J)
    	END DO
    	DO I = 1, I21
        	Q(I,1)   = Q(I,2)
        	Q(I,J20) = Q(I,J19)
    	END DO

    	! SOR METHOD FOR Q
    	DO J = 2, J19
        	DO I = 2, I19
            	ULI = A1 * (Q(I+1,J) + Q(I-1,J)) + A2 * (Q(I,J+1) + Q(I,J-1)) &
                	+ A3 * D(I,J) - Q(I,J)
            	G1 = G1 + ULI * ULI
            	Q(I,J) = ULI + Q(I,J)
        	END DO
    	END DO

    	IF (G1 <= 0.000001) EXIT
	END DO

	! CORRECT VELOCITIES U AND V
	DO J = 2, J19
    	DO I = 3, I19
        	U(I,J) = U(I,J) + (Q(I,J) - Q(I-1,J)) * XD * AFCT
    	END DO
	END DO

	DO J = 3, J19
    	DO I = 2, I19
        	V(I,J) = V(I,J) + (Q(I,J) - Q(I,J-1)) * YD * AFCT
    	END DO
	END DO

	! CHECK DIVERGENCE
	DIV2 = 0.0
	DO J = 3, J19
    	DO I = 3, I19
        	U1 = (U(I+1,J) - U(I,J)) * XD
        	V2 = (V(I,J+1) - V(I,J)) * YD
        	DIV2 = DIV2 + ABS(U1 + V2)
    	END DO
	END DO

	IF (MOD(L, 1000) == 0) PRINT *, L, K, G1, DIV2
END DO

PRINT *, '---U---'
DO J = 2, J20
    PRINT '(I3, 1X, F8.5)', J-1, U(MX/2+1, J)
END DO
PRINT *, '---V---'
DO I = 2, I20
    PRINT '(I3, 1X, F8.5)', I-1, V(I, MY/2+1)
END DO

END PROGRAM cubic_cavity_smac
