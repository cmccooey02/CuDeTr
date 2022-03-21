MODULE green_function

USE constants
USE linear_algebra

IMPLICIT NONE

PRIVATE

PUBLIC                   :: GREEN
PUBLIC                   :: GREEN_sur_lead
PUBLIC                   :: GREEN_lin_eq
CONTAINS

SUBROUTINE GREEN(BigG, n, go, vpotential)

  INTEGER                   :: n, r
  COMPLEX(dp)               :: go(1:n, 1:n)
  COMPLEX(dp)	              :: vpotential(1:n, 1:n)
  COMPLEX(dp)	              :: BigG(1:n, 1:n)
  COMPLEX(dp)               :: WORKA(1:n, 1:n)
  COMPLEX(dp)               :: WORKB(1:n, 1:n)
  COMPLEX(dp)               :: alpha, beta, ONE

  ONE =  DCMPLX(1.0_dp,0.0_dp)

  alpha = DCMPLX(-1.0_dp,0.0_dp)
  beta  = DCMPLX(0.0_dp,0.0_dp)

  CALL zmm(alpha,n,go,'N',vpotential, 'N', WORKA, beta)

  DO r = 1,n
  WORKA(r,r) = WORKA(r,r) + ONE
  ENDDO

  CALL inv_matrix(WORKB,WORKA, n)

  alpha = DCMPLX(1.0_dp,0.0_dp)
  beta  = DCMPLX(0.0_dp,0.0_dp)

  CALL zmm(alpha,n,WORKB,'N',go, 'N', BigG, beta)

END SUBROUTINE


SUBROUTINE GREEN_sur_lead(WORKA, Ndeep, gamm, E)

!!!! Calculates the initial surface greens function on the surface of the leads


  INTEGER                   :: m, n, p, Ndeep
  REAL(dp)                  :: ndone, gamm, abgamm, E, l, h
  COMPLEX(dp)               :: Gzero
  COMPLEX(dp)               :: WORKA(1:Ndeep, 1:Ndeep)

  abgamm = ABS(gamm)
  ndone = 1.0_dp/(REAL(Ndeep+1,dp))

  DO m = 1, Ndeep
      DO n = 1,Ndeep
          Gzero = DCMPLX(0.0_dp,0.0_dp)
          DO p = 1, Ndeep
             l = E + 2.0_dp*abgamm*COS((REAL(p,dp))*PI*ndone)
             h = l**2.0_dp - 4.0_dp*(abgamm*abgamm)
             IF ( h > 0.0_dp .AND. l < 0.0_dp) THEN
                Gzero = SIN((REAL(p*m)*PI*ndone))*SIN((REAL(p*n))*PI*ndone)*(l+SQRT(h+ZERO)) + Gzero
             ELSE
                Gzero = SIN((REAL(p*m))*PI*ndone)*SIN((REAL(p*n))*ndone*PI)*(l-SQRT(h+ZERO)) + Gzero
             END IF
          ENDDO
          WORKA(m,n) = ndone*Gzero/(abgamm*abgamm)
      ENDDO
  ENDDO

END SUBROUTINE


SUBROUTINE GREEN_lin_eq(go, vpotential, n)

!! Solves the Dyson equation using a set of linear equations
!! This method is similar to the one used in GREEN() but doesn't compute
!! the inverse of 1-gV.

INTEGER                   :: n, r
COMPLEX(dp)               :: go(1:n, 1:n)
COMPLEX(dp)	              :: vpotential(1:n, 1:n)
COMPLEX(dp)               :: WORKA(1:n, 1:n)
COMPLEX(dp)               :: alpha, beta, ONE

ONE =  DCMPLX(1.0_dp,0.0_dp)

alpha = DCMPLX(-1.0_dp,0.0_dp)
beta  = DCMPLX(0.0_dp,0.0_dp)

CALL zmm(alpha,n,go,'N',vpotential, 'N', WORKA, beta)

DO r = 1,n
WORKA(r,r) = WORKA(r,r) + ONE
ENDDO

CALL Ax_eq_b(WORKA, go, n )


END SUBROUTINE
END MODULE green_function
