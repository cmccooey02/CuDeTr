MODULE linear_algebra

USE constants

IMPLICIT NONE

PRIVATE

PUBLIC                           :: inv_matrix
PUBLIC                           :: zmm
PUBLIC                           :: Ax_eq_b

INTERFACE inv_matrix

   MODULE PROCEDURE dinv_matrix
   MODULE PROCEDURE zinv_matrix

END INTERFACE inv_matrix

CONTAINS


SUBROUTINE dinv_matrix(invmat, mat, n )
    INTEGER                   :: n
    REAL(dp)	                :: mat(1:n, 1:n)
    REAL(dp)                  :: invmat(1:n, 1:n)
    REAL(dp)                  :: work(1:n)!work array for LAPACK
    INTEGER 	                :: ipiv(1:n)   ! pivot indices
    INTEGER 	                :: info

    ! Store A in Ainv to prevent it from being overwritten by LAPACK

    invmat = mat

   ! ZGETRF computes an LU factorization of a general M-by-N matrix A
   ! using partial pivoting with row interchanges.

   CALL DGETRF(n, n, invmat, n, ipiv, info)

   if (info /= 0) then
      stop 'Matrix is numerically singular!'
   end if

   ! DGETRI computes the inverse of a matrix using the LU factorization
   ! computed by DGETRF.

   CALL DGETRI(n, invmat, n, ipiv, work, n, info)

   if (info /= 0) then
      stop 'Matrix inversion failed!'
   end if

END SUBROUTINE dinv_matrix


SUBROUTINE zinv_matrix( invmat, mat, n )
    INTEGER                   :: n
    COMPLEX(dp)	              :: mat(1:n, 1:n)
    COMPLEX(dp)               :: invmat(1:n, 1:n)
    COMPLEX(dp)               :: work(1:1) !work array for LAPACK
    COMPLEX(dp), ALLOCATABLE   :: workF(:) !work array for LAPACK
    INTEGER 	                :: ipiv(1:n)   ! pivot indices
    INTEGER 	                :: info, Lwork
    EXTERNAL ZGETRF
    EXTERNAL ZGETRI
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    invmat = mat

   ! ZGETRF computes an LU factorization of a general M-by-N matrix A
   ! using partial pivoting with row interchanges.

   CALL ZGETRF(n, n, invmat, n, ipiv, info)

   if (info /= 0) then
      stop 'Matrix is numerically singular!'
   end if

   ! ZGETRI computes the inverse of a matrix using the LU factorization
   ! computed by DGETRF.

   CALL ZGETRI(n, invmat, n, ipiv, work, -1, info)
   LWORK = INT(REAL(work(1)))
   ALLOCATE(workF(1:LWORK))
   CALL ZGETRI(n, invmat, n, ipiv, workF, LWORK, info)
   if (info /= 0) then
      stop 'Matrix inversion failed!'
   end if
   DEALLOCATE(workF)
END SUBROUTINE

SUBROUTINE zmm(alpha,n,A,opA, B, opB, C, Beta)
   INTEGER       :: n
   CHARACTER(LEN = 1)     :: opA, opB
   COMPLEX(dp)   :: Beta, alpha
   COMPLEX(dp)   :: A(1:n,1:n), B(1:n,1:n), C(1:n,1:n)

   CALL ZGEMM(opA,opB,n,n,n,alpha,A,n,B,n,Beta,C,n)

END SUBROUTINE



SUBROUTINE Ax_eq_b(mat, lilg, n )

  !!!! Routine to solve Dyson equation without inverting 1-gV directly
  !!!! Instead G is found by sovling a system of linear equations with
  !!! an LU factored square matrix.


    INTEGER                   :: n
    COMPLEX(dp)	              :: mat(1:n, 1:n)  ! matrix 1-gV
    COMPLEX(dp)               :: lilg(1:n, 1:n) ! g
    INTEGER 	                :: ipiv(1:n)   ! pivot indices
    INTEGER 	                :: info

    EXTERNAL ZGETRF
    EXTERNAL ZGETRS


    ! Store A in Ainv to prevent it from being overwritten by LAPACK
!    invmat = mat

   ! ZGETRF computes an LU factorization of a general M-by-N matrix A
   ! using partial pivoting with row interchanges.

   CALL ZGETRF(n, n, mat, n, ipiv, info)

   if (info /= 0) then
      stop 'Matrix is numerically singular!'
   end if




   ! Solves a system of linear equations with an LU-factored square coefficient matrix,
   ! with multiple right-hand sides.
   ! call zgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
   CALL ZGETRS('N', n, n, mat, n, ipiv, lilg, n, info)



   if (info /= 0) then
      stop 'Matrix inversion failed!'
   end if
END SUBROUTINE

END MODULE linear_algebra
