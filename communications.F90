!------------------------------------------------------------------------------!
!		     		                                               !
!     Module Communications		                                       !
!		     		                                               !
!------------------------------------------------------------------------------!
!		     		                                               !
! Sets up communications packages for parallelization                          !
! Currently two protocols are coded for:                                       !
!		     		                                               !
!    1. Serial Communications		     		                       !
!    2. MPI-based communications	     		                       !
!		     		                                               !
! A user can switch between theese using the _COM_MPI compiler directive       !
!		     		                                               !
!------------------------------------------------------------------------------!
!		     		                                               !
! Written by Daniel Dundas (Belfast)                                           !
! September 2006		 	                                       !
!		     		                                               !
!------------------------------------------------------------------------------!



MODULE communications
   
   USE constants
   USE fdf
   USE util

#if _MPI

   USE mpi

#endif
  
   IMPLICIT NONE
  
  
  
   PRIVATE
  
  
  
   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Public subroutines                                                        !
   !                                                                           !
   !---------------------------------------------------------------------------!

   PUBLIC :: initialize_communications
   PUBLIC :: output_communications
   PUBLIC :: all_processor_barrier
   PUBLIC :: global_barrier
   PUBLIC :: check_processor_count
   PUBLIC :: finalize_communications
   PUBLIC :: stopping
   PUBLIC :: print_debug

      
      
   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Interfaces                                                                !
   !                                                                           !
   !---------------------------------------------------------------------------!
   
   
   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Public types                                                              !
   !                                                                           !
   !---------------------------------------------------------------------------!
   
   TYPE, PUBLIC :: COMMTYPE
   
      INTEGER                       :: iprocessor
      INTEGER                       :: numprocessorsglobal
      INTEGER                       :: numprocessors
      INTEGER                       :: maxprocessors
      CHARACTER(LEN = 3)            :: cprocessor
   
   END TYPE COMMTYPE
       
   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Public variables                                                          !
   !                                                                           !
   !---------------------------------------------------------------------------!

   TYPE(COMMTYPE), PUBLIC :: communicator
  
   !---------------------------------------------------------------------------!
   !                                                                           !
   ! Private variables                                                         !
   !                                                                           !
   !---------------------------------------------------------------------------!

  
  
   CONTAINS
  
  
  
      !------------------------------------------------------------------------!
      !		     		                                               !
      ! Initialises communications                                             !
      !		     		                                               !
      !------------------------------------------------------------------------!



      SUBROUTINE initialize_communications( )

      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'initialize_communications'
      
      !--subroutine parameters-------------------------------------------------!
      
      !--internal variables----------------------------------------------------!   

#if _COM_MPI

      INTEGER                       :: ierror

#endif
      
      !------------------------------------------------------------------------!
      
#if _COM_MPI

      ! Start up MPI

      CALL MPI_init( ierror )

      ! Find out number of processors.

      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, communicator%numprocessorsglobal,    &
                          ierror )

      ! Find out number of the processor we are working on.

      CALL MPI_COMM_RANK( MPI_COMM_WORLD, communicator%iprocessor, ierror )
      
      ! Initialize fdf
      
      CALL fdf_init( 'simulations.fdf', 'simulations.dat' )      
     
      communicator%numprocessors  = fdf_integer( 'NumProcessors', 1 )
      
      CALL fdf_shutdown( )

#else
   
      communicator%numprocessorsglobal = 1      
      communicator%numprocessors       = 1      
      communicator%iprocessor          = 0
      
#endif
      
      WRITE(communicator%cprocessor, '(I3.3)') communicator%iprocessor



      END SUBROUTINE initialize_communications

      
      
      SUBROUTINE output_communications( )

      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'output_communications'
      
      !--subroutine parameters-------------------------------------------------!
      
      !--internal variables----------------------------------------------------!   

      
      !------------------------------------------------------------------------!
      

      IF (communicator%IPROCESSOR.EQ.0) THEN
      
         WRITE(control_var%output_file, *) 
      
         WRITE(control_var%output_file, *)                                     &
            '--Communications Parameters-------------------------------------'
         
	 WRITE(control_var%output_file, *)                                     

         WRITE(control_var%output_file, '(A, I9)')                             &
            '    Number of Processors                    = ',                  &
	    communicator%numprocessors
         
         WRITE(control_var%output_file, *) 

         WRITE(control_var%output_file, *)                                     &
            '----------------------------------------------------------------'

         WRITE(control_var%output_file, *) 
      
      ENDIF


      END SUBROUTINE output_communications
      
      
      
      !------------------------------------------------------------------------!
      !		     		                                               !
      ! Sets up a barrier call for syncronizationn purposes.                   !
      ! Generally this is not needed for either serial communications          !
      ! or two-sided MPI communications. However, it would be required         !
      !	for one	sided MPI2 communications and it is also very handy for        !
      !	testing and debugging purposes.                                        !
      !		     		                                               !
      !------------------------------------------------------------------------!



      SUBROUTINE all_processor_barrier( )
      
#if _COM_MPI
   
      USE MPI
      
      IMPLICIT NONE
   
      INTEGER       :: ierror
      
      ! The MPI barrier routine.

      CALL MPI_barrier( communicator%COMMSIMULATION, ierror )
      
#endif      
      
      
      END SUBROUTINE all_processor_barrier



      SUBROUTINE global_barrier( )
      
#if _COM_MPI
   
      USE MPI
      
      IMPLICIT NONE
   
      INTEGER       :: ierror
      
      ! The MPI barrier routine.

      CALL MPI_barrier( MPI_COMM_WORLD, ierror )
      
#endif      
      
      
      END SUBROUTINE global_barrier
      
      
      
      !------------------------------------------------------------------------!
      !		     		                                               !
      ! Check that the processors requested in the input fdf files             !
      ! matches the number requested in the job submission stage.              !
      ! Stop if they do not match.                                             !
      !		     		                                               !
      !------------------------------------------------------------------------!



      SUBROUTINE check_processor_count( )
            
      IMPLICIT NONE
      
      INTEGER               :: ierror
      CHARACTER(LEN = 120)  :: outputstring
      CHARACTER(LEN = 4)    :: cprocused, cproc
      


      IF (communicator%numprocessors.NE.                                       &
          communicator%numprocessorsglobal) THEN
	 
	 WRITE(cprocused, '(I4.4)') communicator%numprocessorsglobal
	 WRITE(cproc, '(I4.4)')     communicator%numprocessors

	 outputstring = "Number of processors requested for job " //           &
	                "(N = " // cprocused // ") does not equal the " //     &
			"number specified in input file (N = " // cproc //     &
			") "
	 
	 CALL stopping( outputstring )
      
      ENDIF
      
      
      
      END SUBROUTINE check_processor_count
    
    
    
      !------------------------------------------------------------------------!
      !		     		                                               !
      ! Dinamo stop routine. Allows for graceful stopping of the code          !
      ! by ensuring final call to MPI_finalize if we are using MPI             !
      ! communications.                                                        !
      !		     		                                               !
      !------------------------------------------------------------------------!



      SUBROUTINE stopping( message )

      IMPLICIT NONE
      
      CHARACTER(LEN = *)      :: message
      
      
      
      IF (communicator%IPROCESSOR.EQ.0) WRITE(*, *) message

      CALL finalize_communications( )
         
	 
  
      END SUBROUTINE stopping
      
      
      
      !------------------------------------------------------------------------!
      !		     		                                               !
      ! Initialize the processor grid for matrix decomposition.                !
      !		     		                                               !
      !------------------------------------------------------------------------!






      !------------------------------------------------------------------------!
      !		     		                                               !
      ! Finalises communications                                               !
      !		     		                                               !
      !------------------------------------------------------------------------!
      
      
      
      SUBROUTINE finalize_communications ( )

#if _COM_MPI
   
      USE mpi

#endif

      IMPLICIT NONE

#if _COM_MPI
   
      INTEGER :: ierror
      
    
      CALL MPI_COMM_FREE( communicator%COMMSIMULATION, ierror )
      
      CALL global_barrier( )
      
#if _SCALAPACK
   
      CALL finalize_scalapack( communicator%scalapack%ICTXT )
      
#endif 
      
      CALL MPI_FINALIZE( ierror )
      
#endif 
      
      STOP
      
      

      END SUBROUTINE finalize_communications
   
      
      
      SUBROUTINE print_debug( message )
      
      IMPLICIT NONE
      
      CHARACTER(LEN = *)      :: message
      
      
      CALL all_processor_barrier( )
      
      IF (communicator%iprocessor.EQ.0) WRITE(*, *) message
         
	 
  
      END SUBROUTINE print_debug
   
   
      !------------------------------------------------------------------------!
      !                                                                        !
      ! Sum the elements of an integer array across processors.                !
      ! (Element 1 is summed across processors, element 2 is, and so on.)      !
      !                                                                        !
      !------------------------------------------------------------------------!






END MODULE communications
