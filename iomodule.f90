MODULE iomodule
   
   
   
   USE communications
   USE constants
   USE fdf
   USE util
   
   IMPLICIT NONE
   
   PRIVATE
   PUBLIC                         :: initialize_io
   PUBLIC                         :: read_simulation_parameters
   
   CHARACTER(LEN = 100), PUBLIC   :: iofilename
   
   CONTAINS
   
   
   
      SUBROUTINE initialize_io

      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'initialize_io'
      
      !--subroutine parameters ------------------------------------------------!
      
      !--internal variables ---------------------------------------------------!  
      
      CHARACTER(LEN = 200)          :: filename
      
      !------------------------------------------------------------------------!

      
      
      
      control_var%screen_output_desired   =                                    &
         fdf_boolean( 'ScreenOutputDesired', .TRUE. )

      control_var%checkpoint_desired =                                         &
         fdf_boolean( 'CheckpointDesired', .TRUE. )
      
      control_var%screen_file       = fdf_string( 'ScreenFile',                &
                                                  'echoscreen.dat' )
      
      control_var%output_file            = 6
      
      IF (communicator%iprocessor.EQ.0 .AND.                                   &
          .NOT.control_var%screen_output_desired) THEN
	 
         OPEN(UNIT = control_var%output_file, FORM = 'formatted',              &
	      FILE = control_var%screen_file)
	 
      ENDIF
      
      
      END SUBROUTINE initialize_io


      



      SUBROUTINE read_simulation_parameters( )
            
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER     :: myname = 'read_simulation_parameters'
      
      !--subroutine parameters-------------------------------------------------!

      !--internal variables----------------------------------------------------!   

      CHARACTER(LEN = 100) :: filenamein, filenameout
      CHARACTER(LEN = 4)   :: initialrun 
      
      !------------------------------------------------------------------------!
      
      
      filenamein  = communicator%cprocessor // '/conductor.fdf'
      filenameout = communicator%cprocessor // '/conductor.dat'
      
      ! Initialize fdf
      
      CALL fdf_init( filenamein, filenameout )
           
      END SUBROUTINE read_simulation_parameters



END MODULE iomodule
