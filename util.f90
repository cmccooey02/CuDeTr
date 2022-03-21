!------------------------------------------------------------------------------!
!                                                                              !
!               Module Util                                                    !
!                                                                              !
!------------------------------------------------------------------------------!
!                                                                              !
! This module contains relativelly low level                                   !
! such as error handling and input/output                                      !
!                                                                              !
!------------------------------------------------------------------------------!
!                                                                              !
!------------------------------------------------------------------------------!



MODULE util

   USE constants

   IMPLICIT NONE
   
   PRIVATE

   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Public types                                                              !
   !                                                                           !
   !---------------------------------------------------------------------------!

   TYPE, PUBLIC :: CONTROL_TYPE
      
      LOGICAL              :: initialized
      INTEGER              :: output_file
      INTEGER              :: popinnerfile
      INTEGER              :: poptotalfile
      INTEGER              :: dipoleaccelerationfile
      INTEGER              :: dipolemomentfile
      INTEGER              :: harxcmomfile
      INTEGER              :: apulsefile
      INTEGER              :: epulsefile
      INTEGER              :: kefile
      INTEGER              :: iontempfile
      LOGICAL              :: graphics_output_desired
      LOGICAL              :: momentum_desired
      LOGICAL              :: output_prob_xy
      LOGICAL              :: output_prob_xz
      LOGICAL              :: output_prob_yz
      LOGICAL              :: output_psi
      LOGICAL              :: output_gridcoords
      LOGICAL              :: output_absorber
      LOGICAL              :: output_dipolemoment
      LOGICAL              :: output_harxcmom
      LOGICAL              :: screen_output_desired
      LOGICAL              :: checkpoint_desired
      LOGICAL              :: psi_reduced
      LOGICAL              :: calcionize
      INTEGER              :: min_xproc_for_psi_output
      INTEGER              :: max_xproc_for_psi_output
      INTEGER              :: min_yproc_for_psi_output
      INTEGER              :: max_yproc_for_psi_output
      INTEGER              :: min_zproc_for_psi_output
      INTEGER              :: max_zproc_for_psi_output
      INTEGER              :: output_level
      INTEGER              :: debug_level
      INTEGER              :: ran3_seed
      CHARACTER(LEN = 100) :: screen_file
      CHARACTER(LEN = 100) :: poptotal_file
      CHARACTER(LEN = 100) :: popinner_file
      CHARACTER(LEN = 100) :: dipoleacceleration_file
      CHARACTER(LEN = 100) :: dipolemoment_file
      CHARACTER(LEN = 100) :: harxcmom_file
      CHARACTER(LEN = 100) :: apulse_file
      CHARACTER(LEN = 100) :: epulse_file
      CHARACTER(LEN = 100) :: ke_file
      CHARACTER(LEN = 100) :: iontemp_file
      CHARACTER(LEN = 100) :: datadirectory
      CHARACTER(LEN = 20)  :: job_name
   
   END TYPE control_type

   TYPE, PUBLIC :: GENERAL_TYPE
   
      INTEGER             :: idevicetype
      CHARACTER(LEN = 20) :: cdevicetype
        
   END TYPE general_type

   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Public variables                                                          !
   !                                                                           !
   !---------------------------------------------------------------------------!

   TYPE(CONTROL_TYPE), PUBLIC, SAVE :: control_var
   TYPE(GENERAL_TYPE), PUBLIC, SAVE :: general

   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Public subroutines                                                        !
   !                                                                           !
   !---------------------------------------------------------------------------!

   PUBLIC :: leqi
   PUBLIC :: ran3  
   PUBLIC :: derf
   PUBLIC :: derfc
   PUBLIC :: linear_interpolate
   PUBLIC :: cubic_interpolate

   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Private variables                                                         !
   !                                                                           !
   !---------------------------------------------------------------------------!



   CONTAINS


      !------------------------------------------------------------------------!
      !                                                                        !
      ! Random number generation from Numerical Recipes                        !
      !                                                                        !
      !------------------------------------------------------------------------!



      FUNCTION ran3( )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'ran3'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp )                     :: ran3
      
      !--internal variables ---------------------------------------------------!  
      
      INTEGER, PARAMETER            :: MBIG  = 1000000000
      INTEGER, PARAMETER            :: MSEED = 161803398
      INTEGER, PARAMETER            :: MZ    = 0
      REAL(dp)                      :: FAC   = 1.0d0 / MBIG
      INTEGER                       :: i, ii, k
      INTEGER                       :: mj, mk
      INTEGER                       :: iff = 0
      INTEGER, SAVE                 :: inext
      INTEGER, SAVE                 :: inextp
      INTEGER, SAVE                 :: ma(55)
      
      !------------------------------------------------------------------------!
     
      
      
      IF (control_var%ran3_seed.LT.0 .OR. iff.EQ.0) THEN
	 
	 iff    = 1
         mj     = MSEED - ABS( control_var%ran3_seed )
         mj     = MOD( mj, MBIG )
         ma(55) = mj
         mk     = 1
	 
         DO i = 1, 54
	    
	    ii     = MOD( 21 * i, 55 )
            ma(ii) = mk
            mk     = mj - mk
	    
	    IF (mk.LT.MZ) mk = mk + MBIG
            
	    mj = ma(ii) 
	 
	 END DO
	 
	 DO k = 1, 4
            DO i = 1, 55
	       
	       ma(i) = ma(i) - ma(1 + MOD( i + 30, 55 ))
               
	       IF (ma(i).LT.MZ) ma(i) = ma(i) + MBIG
	    
	    ENDDO
         ENDDO
	 
	 inext                 = 0
         inextp                = 31
         control_var%ran3_seed = 1
      
      ENDIF
      
      inext = inext + 1
      
      IF (inext.EQ.56) inext = 1
      
      inextp = inextp + 1
      
      IF (inextp.EQ.56) inextp = 1
      
      mj = ma(inext) - ma(inextp)
      
      IF (mj.LT.MZ) mj = mj + mbig
      
      ma(inext) = mj
      ran3      = mj * fac


      
      END FUNCTION ran3


      
      !------------------------------------------------------------------------!
      !                                                                        !
      ! Fermi distribution function                                            !
      !                                                                        !
      !------------------------------------------------------------------------!



      FUNCTION fermi( temp, energy, mu )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'fermi'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp)                      :: fermi
      REAL(dp)                      :: temp
      REAL(dp)                      :: energy
      REAL(dp)                      :: mu
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp)                      :: expo
      
      !------------------------------------------------------------------------!



      expo = (energy - mu) / (temp * boltzmann_k)

      IF (expo.LT.-100.0_dp) THEN
	 
	 fermi = 1.0_dp
      
      ELSEIF (expo.GT.100.0_dp) THEN
	
	 fermi = 0.0_dp
      
      ELSE
	 
	 fermi = 1.0_dp / (EXP((energy - mu) / (temp * boltzmann_k)) + 1.0_dp)
      
      ENDIF



      END FUNCTION fermi



      !------------------------------------------------------------------------!
      !                                                                        !
      ! Differential of the Fermi distributionfunction wrt mu                  !
      !                                                                        !
      !------------------------------------------------------------------------!
      
      
      
      FUNCTION fermidiff( temp, energy, mu )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'fermidiff'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp) :: fermidiff
      REAL(dp) :: temp
      REAL(dp) :: energy
      REAL(dp) :: mu
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp) :: expo
      
      !------------------------------------------------------------------------!



      expo = (energy - mu) / (temp * boltzmann_k)

      IF (expo.LT.-100.0_dp) THEN
	 
	 fermidiff = 0.0_dp
      
      ELSEIF (expo.GT.100.0_dp) THEN
	 
	 fermidiff = 0.0_dp
      
      ELSE
	 
	 fermidiff = (EXP((energy - mu) / (temp * boltzmann_k)) /              &
	              (temp * boltzmann_k)) /                                  &
		     (EXP((energy - mu) / (temp * boltzmann_k)) + 1.0_dp) ** 2
      
      ENDIF



      END FUNCTION fermidiff



      !------------------------------------------------------------------------!
      !                                                                        !
      ! modulus of a 3D vector                                                 !
      !                                                                        !
      !------------------------------------------------------------------------!



      FUNCTION vecmod( x, y, z )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'vecmod'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp)                      :: vecmod
      REAL(dp)                      :: x, y, z
      
      !--internal variables ---------------------------------------------------!  
      
      !------------------------------------------------------------------------!



      vecmod = SQRT( x * x + y * y + z * z )



      END FUNCTION vecmod

      
      
      !------------------------------------------------------------------------!
      !                                                                        !
      ! Case-insensitive lexical equal-to comparison                           !
      !                                                                        !
      !------------------------------------------------------------------------!
      
      
      
      FUNCTION leqi( strng1, strng2 )
      
      IMPLICIT NONE
     
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'leqi'
      
      !--subroutine parameters ------------------------------------------------!
      
      LOGICAL                :: leqi
      CHARACTER(LEN = *)     :: strng1, strng2
      
      !--internal variables ---------------------------------------------------!  
      
      INTEGER                 :: len1, len2, lenc, i
      CHARACTER(LEN = 1)      :: s1, s2
      
      !------------------------------------------------------------------------!
      
      
      
      len1 = LEN( strng1 )
      len2 = LEN( strng2 )
      lenc = MIN( len1, len2 )
     
      leqi = .FALSE.
      
      DO i = 1, lenc
	 
	 s1 = strng1(i:i)
         s2 = strng2(i:i)
	 
	 CALL CHRCAP( s1, 1 )
         CALL CHRCAP( s2, 1 )
	 
	 IF (s1.NE.s2) RETURN
      
      ENDDO
      
      IF (len1.GT.lenc .AND. strng1((lenc + 1):len1).NE.' ') RETURN
      IF (len2.GT.lenc .AND. strng2((lenc + 1):len2).NE.' ') RETURN
      
      leqi = .TRUE.
      
      
      
      END FUNCTION leqi



      !------------------------------------------------------------------------!
      !                                                                        !
      ! gritoc                                                                 !
      !                                                                        !
      !------------------------------------------------------------------------!
     
     
     
      FUNCTION gritoc( int, str )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'gritoc'
      
      !--subroutine parameters ------------------------------------------------!
      
      INTEGER                 :: gritoc
      INTEGER                 :: int
      CHARACTER(LEN = *)      :: str
      
      !--internal variables ---------------------------------------------------!  
      
      CHARACTER(LEN = 10)     :: digits
      INTEGER                 :: d, i, intval, j, l
      CHARACTER               :: k

      !------------------------------------------------------------------------!
      
      
      
      DATA digits /'0123456789'/
      
      intval = ABS(int)
      i      = 0
      str    = ''
      
  10  CONTINUE
  
         i        = i + 1
         d        = 1 + MOD(intval, 10)
	 str(i:i) = digits(d:d)
         intval   = intval / 10
	 
	 IF (i.LT.LEN(str) .AND. intval.NE.0) GOTO 10
	 
      IF (int.LT.0 .AND. i.LT.LEN(str)) THEN
	 
	 i        = i + 1
         str(i:i) = '-'
      
      END IF
      
      gritoc = i
      l      = i / 2
      
      DO j = 1, l
	  
	 k        = str(i:i)
         str(i:i) = str(j:j)
         str(j:j) = k
         i        = i - 1
	  
      ENDDO
   
   
   
      END function gritoc



      !------------------------------------------------------------------------!
      !                                                                        !
      ! Inverse of a 3x3 matrix                                                !
      !                                                                        !
      !------------------------------------------------------------------------!



      SUBROUTINE minv( c, a )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'minv'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp), INTENT(IN)          :: a(3, 3) 
      REAL(dp), INTENT(OUT)         :: c(3, 3) 
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp)                      :: a11, a12, a13
      REAL(dp)                      :: a21, a22, a23
      REAL(dp)                      :: a31, a32, a33
      REAL(dp)                      :: denominator
      
      !------------------------------------------------------------------------!



      a11     = a(1, 1) 
      a12     = a(1, 2) 
      a13     = a(1, 3)
      a21     = a(2, 1) 
      a22     = a(2, 2) 
      a23     = a(2, 3)
      a31     = a(3, 1) 
      a32     = a(3, 2) 
      a33     = a(3, 3)
      
      denominator = (a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 -    &
		     a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31)
      denominator = 1.0_dp / denominator
      
      c(1, 1)     = (-(a23 * a32) + a22 * a33) * denominator
      c(1, 2)     = ( (a13 * a32) - a12 * a33) * denominator
      c(1, 3)     = (-(a13 * a22) + a12 * a23) * denominator
      c(2, 1)     = ( (a23 * a31) - a21 * a33) * denominator
      c(2, 2)     = (-(a13 * a31) + a11 * a33) * denominator
      c(2, 3)     = ( (a13 * a21) - a11 * a23) * denominator
      c(3, 1)     = (-(a22 * a31) + a21 * a32) * denominator
      c(3, 2)     = ( (a12 * a31) - a11 * a32) * denominator
      c(3, 3)     = (-(a12 * a21) + a11 * a22) * denominator




      END SUBROUTINE minv




      !------------------------------------------------------------------------!
      !                                                                        !
      ! scalar product of two 3 component vectors                              !
      !                                                                        !
      !------------------------------------------------------------------------!
      
      
      
      FUNCTION sprod( v1, v2 )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
     
      CHARACTER(LEN = *), PARAMETER :: myname = 'sprod'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp) :: sprod
      REAL(dp) :: v1(3), v2(3)
      
      !--internal variables ---------------------------------------------------!  
      
      !------------------------------------------------------------------------!



      sprod = SUM( v1 * v2 )



      END FUNCTION sprod



      !------------------------------------------------------------------------!
      !                                                                        !
      ! apply 3D transformation to vector                                      !
      !                                                                        !
      !------------------------------------------------------------------------!
      
      
      
      SUBROUTINE mvec( a, v )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'mvec'
      
      !--subroutine parameters ------------------------------------------------!
     
      REAL(dp), INTENT(IN)          :: a(3, 3) 
      REAL(dp), INTENT(INOUT)       :: v(3) 
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp)                      :: vt(3)
      
      !------------------------------------------------------------------------!



      vt(1) = a(1, 1) * v(1) + a(1, 2) * v(2) + a(1, 3) * v(3)
      vt(2) = a(2, 1) * v(1) + a(2, 2) * v(2) + a(2, 3) * v(3)
      vt(3) = a(3, 1) * v(1) + a(3, 2) * v(2) + a(3, 3) * v(3)

      v     = vt



      END SUBROUTINE mvec 
  
  
  
      !------------------------------------------------------------------------!
      !                                                                        !
      ! apply 3D transformation to vector                                      !
      !                                                                        !
      !------------------------------------------------------------------------!
      
      
      
      SUBROUTINE mvec3( a, v1, v2, v3 )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'mvec3'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp), INTENT(IN)          :: a(3, 3) 
      REAL(dp), INTENT(INOUT)       :: v1, v2, v3 
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp)                      :: vt1, vt2, vt3
      
      !------------------------------------------------------------------------!



      vt1 = a(1, 1) * v1 + a(1, 2) * v2 + a(1, 3) * v3
      vt2 = a(2, 1) * v1 + a(2, 2) * v2 + a(2, 3) * v3
      vt3 = a(3, 1) * v1 + a(3, 2) * v2 + a(3, 3) * v3

      v1  = vt1
      v2  = vt2
      v3  = vt3



      END SUBROUTINE mvec3 



      !------------------------------------------------------------------------!
      !                                                                        !
      ! build 3D rotation matrix from Euler angles                             !
      !                                                                        !
      !------------------------------------------------------------------------!
 


      SUBROUTINE rmatrix(a, phi, theta, psi)
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
     
      CHARACTER(LEN = *), PARAMETER :: myname = 'rmatrix'
     
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp), INTENT(OUT)         :: a(3, 3) 
      REAL(dp), INTENT(IN)          :: phi, theta, psi 
     
      !--internal variables ---------------------------------------------------!  
      
      !------------------------------------------------------------------------!



      a(1, 1) =  COS(phi)   * COS(psi)                -                            &
                 COS(theta) * SIN(phi)   * SIN(psi)
      a(1, 2) =  COS(psi)   * SIN(phi)                +                            &
                 COS(phi)   * COS(theta) * SIN(psi)
      a(1, 3) =  SIN(psi)   * SIN(theta)
      a(2, 1) = -COS(psi)   * COS(theta) * SIN(phi)   - COS(phi) * SIN(psi)
      a(2, 2) =  COS(phi)   * COS(psi)   * COS(theta) - SIN(phi) * SIN(psi)
      a(2, 3) =  COS(psi)   * SIN(theta)
      a(3, 1) =  SIN(phi)   * SIN(theta)
      a(3, 2) = -COS(phi)   * SIN(theta)
      a(3, 3) =  COS(theta)



      END SUBROUTINE rmatrix
   
   
   
      !------------------------------------------------------------------------!
      !                                                                        !
      ! calculates reciprocal lattice vectors. their                           !
      ! product with direct lattice vectors is 1 if                            !
      ! iopt = 0 or 2 * pi if iopt = 1                                         !
      !                                                                        !
      !------------------------------------------------------------------------!
   
   
   
      SUBROUTINE reclat( a, b, iopt )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'reclat'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp), INTENT(OUT)         :: b(3, 3) 
      REAL(dp), INTENT(IN)          :: a(3, 3)
      LOGICAL, INTENT(IN)           :: iopt
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp)                      :: c, ci
      INTEGER                       :: i
      
      !------------------------------------------------------------------------!



      b(1, 1) = a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3)
      b(2, 1) = a(3, 2) * a(1, 3) - a(1, 2) * a(3, 3)
      b(3, 1) = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3)
      b(1, 2) = a(2, 3) * a(3, 1) - a(3, 3) * a(2, 1)
      b(2, 2) = a(3, 3) * a(1, 1) - a(1, 3) * a(3, 1)
      b(3, 2) = a(1, 3) * a(2, 1) - a(2, 3) * a(1, 1)
      b(1, 3) = a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2)
      b(2, 3) = a(3, 1) * a(1, 2) - a(1, 1) * a(3, 2)
      b(3, 3) = a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2)

      c = 1.0_dp
      
      IF (iopt) c = 2.0_dp * pi

      DO i = 1, 3
	 
	 ci      = c / (a(1, i) * b(1, i) + a(2, i) * b(2, i) +                &
	                a(3, i) * b(3, i))
			
         b(1, i) = b(1, i) * ci
         b(2, i) = b(2, i) * ci
         b(3, i) = b(3, i) * ci
     
      END DO



      END SUBROUTINE reclat



      !------------------------------------------------------------------------!
      !                                                                        !
      ! Complementary error function from "numerical recipes"                  !
      !                                                                        !
      !------------------------------------------------------------------------!
     
     
     
      FUNCTION derfc( x )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'derfc'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp)                      :: derfc
      REAL(dp), INTENT(IN)          :: x 
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp)                      :: z, t
      
      !------------------------------------------------------------------------!



      z     = ABS(x)
      t     = 1.0_dp / (1.0_dp + 0.5_dp * z)
      
      derfc = t * EXP(-(z * z) - 1.26551223_dp +                               &
                      t * (1.00002368_dp + t * ( 0.37409196_dp +               &
                      t * (0.09678418_dp + t * (-0.18628806_dp +               &
                      t * (0.27886807_dp + t * (-1.13520398_dp +               &
                      t * (1.48851587_dp + t * (-0.82215223_dp +               &
		      t * 0.17087277_dp)))))))))
      
      IF (x.LT.0.0_dp) derfc = 2.0_dp - derfc



      END FUNCTION derfc

    
    
      !------------------------------------------------------------------------!
      !                                                                        !
      ! error function from "numerical recipes"                                !
      !                                                                        !
      !------------------------------------------------------------------------!



      FUNCTION derf( x )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'derf'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp)                      :: derf
      REAL(dp), INTENT(IN)          :: x 
      
      !--internal variables ---------------------------------------------------!  
      
      REAL(dp)                      :: z, t
      
      !------------------------------------------------------------------------!



      z    = ABS(x)
      t    = 1.0_dp / (1.0_dp + 0.5_dp * z)

      derf = t * EXP(-(z * z) - 1.26551223_dp +                                &
                     t * (1.00002368_dp + t * ( 0.37409196_dp +                &
                     t * (0.09678418_dp + t * (-0.18628806_dp +                &
                     t * (0.27886807_dp + t * (-1.13520398_dp +                &
                     t * (1.48851587_dp + t * (-0.82215223_dp +                &
		     t * (0.17087277_dp))))))))))

      IF (x.LT.0.0_dp) derf = 2.0_dp - derf

      derf = 1.d0 - derf



      END FUNCTION derf

    
    
      !------------------------------------------------------------------------!
      !                                                                        !
      ! Cubic interpolation of the function value from 4 points in             !
      ! neighbourhood                                                          !
      !                                                                        !
      !------------------------------------------------------------------------!



      FUNCTION cubic_interpolate( zpt, xpts, ypts, numpts, ierr )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'cubic_interpolate'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp)                      :: cubic_interpolate
      INTEGER, INTENT(INOUT)        :: ierr
      INTEGER, INTENT(IN)           :: numpts
      REAL(dp), INTENT(IN)          :: zpt
      REAL(dp), INTENT(IN)          :: xpts(1:numpts)
      REAL(dp), INTENT(IN)          :: ypts(1:numpts)
      
      !--internal variables ---------------------------------------------------!  
      
      INTEGER                       :: ipt
      REAL(dp)                      :: x1, x2, x3, x4, y1, y2, y3, y4
      REAL(dp)                      :: l1, l2, l3, l4
      
      !------------------------------------------------------------------------!


      
      ierr = 0
      
      DO ipt = 1, numpts - 1 
      
         IF (zpt.EQ.xpts(ipt)) THEN
         
	    cubic_interpolate = ypts(ipt)
	    
	    RETURN
         
         ELSEIF (zpt.EQ.xpts(ipt + 1)) THEN
         
	    cubic_interpolate = ypts(ipt + 1)
	    
	    RETURN
      
         ELSEIF (zpt.GT.xpts(ipt) .AND. zpt.LT.xpts(ipt + 1)) THEN
      
            IF (ipt.EQ.1) THEN
	    
	       x1 = xpts(1)
	       x2 = xpts(2)
	       x3 = xpts(3)
	       x4 = xpts(4)
	    
	       y1 = ypts(1)
	       y2 = ypts(2)
	       y3 = ypts(3)
	       y4 = ypts(4)
	 
	    ELSEIF (ipt.EQ.(numpts - 1)) THEN
	    
	       x1 = xpts(numpts - 3)
	       x2 = xpts(numpts - 2)
	       x3 = xpts(numpts - 1)
	       x4 = xpts(numpts)
	    
	       y1 = ypts(numpts - 3)
	       y2 = ypts(numpts - 2)
	       y3 = ypts(numpts - 1)
	       y4 = ypts(numpts)
	 
	    ELSE
	    
	       x1 = xpts(ipt - 1)
	       x2 = xpts(ipt)
	       x3 = xpts(ipt + 1)
	       x4 = xpts(ipt + 2)
	    
	       y1 = ypts(ipt - 1)
	       y2 = ypts(ipt)
	       y3 = ypts(ipt + 1)
	       y4 = ypts(ipt + 2)
	 
	    ENDIF
	 
	    l1 = (zpt - x2) * (zpt - x3) * (zpt - x4) /                        &
	         ((x1 - x2) * (x1 - x3) * (x1 - x4))
	    l2 = (zpt - x1) * (zpt - x3) * (zpt - x4) /                        &
	         ((x2 - x1) * (x2 - x3) * (x2 - x4))
	    l3 = (zpt - x1) * (zpt - x2) * (zpt - x4) /                        &
	         ((x3 - x1) * (x3 - x2) * (x3 - x4))
	    l4 = (zpt - x1) * (zpt - x2) * (zpt - x3) /                        &
	         ((x4 - x1) * (x4 - x2) * (x4 - x3))
         
	    cubic_interpolate = l1 * y1 + l2 * y2 + l3 * y3 + l4 * y4
	    
	    RETURN
      
         ENDIF
      
      ENDDO
      
      ! Cubic interpolation point does not lie within range of data
      
      ierr = 4


      END FUNCTION cubic_interpolate



      FUNCTION linear_interpolate( zpt, xpts, ypts, numpts, ierr )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'linear_interpolate'
      
      !--subroutine parameters ------------------------------------------------!
      
      REAL(dp)                      :: linear_interpolate
      INTEGER, INTENT(INOUT)        :: ierr
      INTEGER, INTENT(IN)           :: numpts
      REAL(dp), INTENT(IN)          :: zpt
      REAL(dp), INTENT(IN)          :: xpts(1:numpts)
      REAL(dp), INTENT(IN)          :: ypts(1:numpts)
      
      !--internal variables ---------------------------------------------------!  
      
      INTEGER                       :: ipt
      REAL(dp)                      :: x1, x2
      REAL(dp)                      :: y1, y2
      REAL(dp)                      :: l1, l2
      
      !------------------------------------------------------------------------!


      
      ierr = 0
      
      DO ipt = 1, numpts - 1 
      
         IF (zpt.EQ.xpts(ipt)) THEN
         
	    linear_interpolate = ypts(ipt)
	    
	    RETURN
         
         ELSEIF (zpt.EQ.xpts(ipt + 1)) THEN
         
	    linear_interpolate = ypts(ipt + 1)
	    
	    RETURN
      
         ELSEIF (zpt.GT.xpts(ipt) .AND. zpt.LT.xpts(ipt + 1)) THEN
      
	    x1 = xpts(ipt)
	    x2 = xpts(ipt + 1)
	    
	    y1 = ypts(ipt)
	    y2 = ypts(ipt + 1)
	 
	    l1 = (zpt - x2) / (x1 - x2)
	    l2 = (zpt - x1) / (x2 - x1)
         
	    linear_interpolate = l1 * y1 + l2 * y2
	    
	    RETURN
      
         ENDIF
      
      ENDDO
      
      ! Cubic interpolation point does not lie within range of data
      
      ierr = 4


      END FUNCTION linear_interpolate



END MODULE util



!------------------------------------------------------------------------------!
!                                                                              !
! End of Module Util                                                           !
!                                                                              !
!------------------------------------------------------------------------------!
