MODULE constants

IMPLICIT NONE

PRIVATE


INTEGER, PUBLIC                  :: NATOMS, Nlayers, Ndeep

INTEGER, PARAMETER, PUBLIC       :: dp       = KIND(1.0D0)
INTEGER, PARAMETER, PUBLIC       :: runkink = 1
INTEGER, PARAMETER, PUBLIC       :: runqdot = 2

REAL(dp), PARAMETER, PUBLIC	       :: SQRTTWO   = 1.4142135623730951_dp
COMPLEX(dp), PARAMETER, PUBLIC         :: ZIMAGONE  = DCMPLX(0.0_dp, 1.0_dp)
COMPLEX(dp), PARAMETER, PUBLIC         :: ZONE      = DCMPLX(1.0_dp, 0.0_dp)
COMPLEX(dp), PARAMETER, PUBLIC         :: ZERO      = DCMPLX(0.0_dp, 0.0_dp)

   ! physical and mathematical

   REAL(dp), PARAMETER, PUBLIC :: amu_to_internal       = 103.6426867_dp
   REAL(dp), PARAMETER, PUBLIC :: boltzmann_k           = 8.61734215E-5_dp         ! in eV / K
   REAL(dp), PARAMETER, PUBLIC :: ev_to_hartree         = 0.03674932600434263_dp
   REAL(dp), PARAMETER, PUBLIC :: hatree_to_ev          = 27.211383411_dp
   REAL(dp), PARAMETER, PUBLIC :: hbar                  = 0.65821188926_dp ! in eV fs
   REAL(dp), PARAMETER, PUBLIC :: pi                    = 3.141592653589793_dp
   REAL(dp), PARAMETER, PUBLIC :: bohr_to_ang           = 0.5291772083_dp
   REAL(dp), PARAMETER, PUBLIC :: efs_to_microamps      = 160.2176462_dp
   REAL(dp), PARAMETER, PUBLIC :: epsilon0              = 5.526349954E-3_dp



END MODULE constants
