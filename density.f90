MODULE density

USE constants
USE linear_algebra

IMPLICIT NONE

PRIVATE

PUBLIC                   :: DENSITYLEADS
PUBLIC                   :: DOSfunc
PUBLIC                   :: Current_gen
CONTAINS

SUBROUTINE DENSITYLEADS(DENSITY, Projector ,PVP, BigG, Gdag, n)

  INTEGER                   :: n, i, m, r
  COMPLEX(dp)               :: DENSITY(1:n, 1:n)
  COMPLEX(dp)	              :: Projector(1:n, 1:n)
  COMPLEX(dp)	              :: PVP(1:n, 1:n)
  COMPLEX(dp)	              :: BigG(1:n, 1:n)
  COMPLEX(dp)	              :: Gdag(1:n, 1:n)
  COMPLEX(dp)               :: WORKA(1:n, 1:n)
  COMPLEX(dp)               :: WORKB(1:n, 1:n)
  COMPLEX(dp)               :: alpha, beta
  EXTERNAL  MKL_ZOMATADD

!!! WORKA = PVP(Gdag)
  alpha = DCMPLX(1.0_dp,0.0_dp)
  beta =  DCMPLX(0.0_dp,0.0_dp)
  CALL zmm(alpha,n,PVP,'N',Gdag,'N',WORKA,beta)



!!! WORKB = G(WORKA)
  alpha = DCMPLX(1.0_dp,0.0_dp)
  beta =  DCMPLX(0.0_dp,0.0_dp)
  CALL zmm(alpha,n,BigG,'N',WORKA,'N',WORKB,beta)



!!! WORKA = P*Gdag
  alpha = DCMPLX(1.0_dp,0.0_dp)
  beta =  DCMPLX(0.0_dp,0.0_dp)
  CALL zmm(alpha,n,Projector,'N',Gdag,'N',WORKA,beta)



!!! WORKA = -G*P + WORKA
  alpha = CMPLX(-1.0_dp,0.0_dp)
  beta =  CMPLX(1.0_dp,0.0_dp)
  CALL zmm(alpha,n,BigG,'N',Projector,'N',WORKA,beta)



!!! D = WORKA + WORKB
  alpha = DCMPLX(0.0_dp,-0.5_dp/PI)
  beta =  DCMPLX(0.0_dp,0.5_dp/PI)


!  CALL mkl_zomatadd('C','N','N',n,n,alpha,WORKA,n,beta,WORKB,n,DENSITY,n)
  DO r = 1,n
     DO m = 1,n

       DENSITY(m,r) = WORKA(m,r) + WORKB(m,r)
       DENSITY(m,r) = alpha*DENSITY(m,r)
     ENDDO
  ENDDO

END SUBROUTINE


SUBROUTINE DOSfunc(DOSurf, BigG, n)

  INTEGER                   :: n, i
  REAL(dp)                  :: DOSurf(1:n,1)
  COMPLEX(dp)	              :: BigG(1:n, 1:n)

  DO i = 1, n
    DOSurf(i,1) = -DIMAG(BigG(i,i))/PI
  ENDDO
END SUBROUTINE

SUBROUTINE Current_gen(DOSurf, BigG, n)

  INTEGER                   :: n, i
  REAL(dp)                  :: DOSurf(1:n,1)
  COMPLEX(dp)	              :: BigG(1:n, 1:n)

  DO i = 1, n
    DOSurf(i,1) = -DIMAG(BigG(i,i))/PI
  ENDDO
END SUBROUTINE

END MODULE density
