program conductor

   USE constants
   USE communications
   USE fdf
   USE iomodule
   USE linear_algebra
   USE util

   implicit none
   INTEGER :: SNx, NE
   INTEGER :: i, k , m, n, p, j, nmesh
   INTEGER :: numcut, iblk, ik
   INTEGER, ALLOCATABLE :: cut(:)

   CHARACTER(LEN=3) :: cnum, cnumb
   CHARACTER(LEN=100) :: filename, fileid, path
   CHARACTER(8)  :: date
   CHARACTER(10) :: time

   REAL(dp):: ep, alattice, gamm, dE, eion
   COMPLEX(dp), allocatable  :: V(:,:)
   COMPLEX(dp), allocatable  :: PL(:,:), PR(:,:)
   COMPLEX(dp), allocatable  :: g(:,:), Gbig(:,:), gdag(:,:)
   COMPLEX(dp), allocatable  :: DLO(:,:), DRO(:,:), DL(:,:), DR(:,:), DOS(:,:)
   COMPLEX(dp), allocatable  :: INTMED(:,:), INVINTMED(:,:)
   REAL(dp),    allocatable  :: Current(:,:)
   COMPLEX(dp), allocatable  :: CEYE(:,:)
   REAL(dp), allocatable     :: E(:)
   REAL(dp) :: EHigh, eon, edif, ndone, l, h, Eo
   REAL(dp) :: arg, W
   COMPLEX(dp) :: Gzero, alpha, beta
   REAL(dp), allocatable  :: fx(:,:), fy(:,:), Fcurlz(:,:)
   REAL(dp), allocatable  :: ux(:,:), vy(:,:), RDOS(:,:)


! Initialize communications

   CALL initialize_communications( )
   CALL output_communications( )
   CALL check_processor_count( )
   CALL read_simulation_parameters( )

   nmesh = fdf_integer( 'NumMesh', 10 )
   numcut  = fdf_integer( 'NumCut', 16 )
   ALLOCATE(cut(1:numcut))

   IF (fdf_block( 'CutValues', iblk )) THEN
      DO ik = 1,  numcut
     	   READ(iblk, *) cut(ik)
      ENDDO
   ELSE
      CALL stopping( "No cut values" )
   ENDIF

   NE = 3
   SNx = 7*nmesh
   NLayers = (2 + SNx)
   Ndeep = 6*nmesh
   alattice = 5.00_dp/REAL(nmesh,dp)
   Natoms = Nlayers*Ndeep
   ep = 1.0E-20_dp
   gamm = -3.8099_dp/(alattice**2)

   IF (communicator%iprocessor.EQ.0) THEN

      filename = communicator%cprocessor // '/outputconstantscut.dat'
      OPEN(13, file = filename, form = 'unformatted')
      WRITE(13) SNx
      WRITE(13) Ndeep
      WRITE(13) NE
      CLOSE(unit = 13)

   ENDIF

   IF (communicator%iprocessor.EQ.0) THEN

      filename = communicator%cprocessor // '/cut' // communicator%cprocessor // '.dat'
      OPEN(11, file = filename, form = 'formatted')
      WRITE(11,*) 'SNx = ', SNx
      WRITE(11,*) 'Ndeep = ', Ndeep
      WRITE(11,*) 'allattice = ', alattice
      WRITE(11,*) 'gamm = ', gamm
   ENDIF

   allocate(g(1:Natoms,1:Natoms))
   allocate(gdag(1:Natoms,1:Natoms))
   allocate(V(1:Natoms,1:Natoms))
   allocate(PL(1:Natoms,1:Natoms))
   allocate(PR(1:Natoms,1:Natoms))
   allocate(E(1:NE))
   allocate(DLO(1:Natoms,1:Natoms))
   allocate(DRO(1:Natoms,1:Natoms))
   allocate(DL(1:Natoms,1:Natoms))
   allocate(DR(1:Natoms,1:Natoms))
   allocate(DOS(1:Natoms,1:Natoms))
   allocate(Gbig(1:Natoms,1:Natoms))
   allocate(INTMED(1:Natoms,1:Natoms))
   allocate(INVINTMED(1:Natoms,1:Natoms))
   allocate(CEYE(1:Natoms,1:Natoms))
   allocate(Current(1:Natoms,1:Natoms))
   allocate(fx(1:Ndeep,1:Nlayers))
   allocate(fy(1:Ndeep,1:Nlayers))
   allocate(Fcurlz(1:Ndeep,1:Nlayers))
   allocate(ux(1:Ndeep,1:Nlayers))
   allocate(vy(1:Ndeep,1:Nlayers))
   allocate(RDOS(1:NATOMS,1:NATOMS))

   V = ZERO
   PL = ZERO
   PR = ZERO
   g  = ZERO
   gdag = ZERO
   INTMED  = ZERO
   ux = 0.0_dp
   vy = 0.0_dp
   fx = 0.0_dp
   fy = 0.0_dp
   Fcurlz = 0.0_dp

! Setting up the Potential V
   DO i = 1, Natoms - Ndeep
      V(i,i+Ndeep) = gamm
      V(i + Ndeep , i) = gamm
   END DO

   DO i = 1,SNx
     DO j = 1,(Ndeep-1)
        V(i*Ndeep +j ,i*Ndeep+j+1 ) = gamm
        V(i*Ndeep+j+1,i*Ndeep +j ) = gamm
     END DO
   END DO

! LEFT AND RIGHT PROJECTOR OPERATORS

   DO i = 1,Ndeep
      PL(i,i) = 1.0_dp
      PR(Natoms+1-i,Natoms+1-i) =1.0_dp
   END DO

! CONSTRUCT COMPLEX IDENTITY

   DO i = 1,Natoms
      CEYE(i,i) = CMPLX(1.0_dp,0.0_dp)
   END DO


   eion = 10.0_dp
   Eo = 0.0_dp
   eon = -4.0_dp*ABS(gamm)

   dE = (eion + eon - Eo - eon) /(REAL(NE+1_dp, dp))
   edif  = eon - Eo

   DO i = 1, NE
      edif = edif + dE
      E(i) = edif
   ENDDO

   W = 1.0_dp
   EHigh = 2000.0
   ndone = 1.0_dp/(REAL(Ndeep+1,dp))

   IF (communicator%iprocessor.EQ.0) THEN
      WRITE(11,*) 'E0 = ', Eo
      WRITE(11,*) 'EHigh = ', EHigh
      CLOSE(UNIT = 11)
   ENDIF

   general%cdevicetype = fdf_string( 'DeviceType', 'ExtremeKink' )

   IF (leqi(general%cdevicetype, 'ExtremeKink' )) THEN
      general%idevicetype = runkink
   ELSEIF (leqi(general%cdevicetype, 'QuantumDot' )) THEN
      general%idevicetype = runqdot
   ELSE
      CALL stopping( "The requested Device Type is not implemented" )
   ENDIF

   DO k = 1,NE

      DO n  = Ndeep + 1, Natoms-Ndeep
         g(n,n) = 1.0_dp/(E(k) - Eo + CMPLX(0.0,ep))
      END DO

      IF (general%idevicetype.EQ.runkink) THEN
      ! Extreme Kink or cut

         DO j = 2, Nlayers - 1
            DO i = 1,5
               g(Ndeep*(j-1)+i,Ndeep*(j-1)+i) = 1.0_dp/(E(k) - (eion) + CMPLX(0.0,ep))
               g(Ndeep*j-i +1,Ndeep*j - i + 1) = 1.0_dp/(E(k) - (eion) + CMPLX(0.0,ep))
	          ENDDO
         ENDDO

         DO i = 1,SIZE(cut)
	          m = cut(i)
            DO j = 1,m
	             g(Ndeep*(1+nmesh*2+7+i)+5+j,Ndeep*(1+nmesh*2+7+i)+5+j) = 1.0_dp/(E(k) - (eion) +CMPLX(0.0,ep))
               g(Ndeep*(1+nmesh*2+8+i)-5-j,Ndeep*(1+nmesh*2+8+i)-5-j) = 1.0_dp/(E(k)- (eion) + CMPLX(0.0,ep))
            ENDDO
         ENDDO
      ELSEIF (general%idevicetype.EQ.runqdot) THEN
	    ! Settings soft wall

	       DO j = 2,Nlayers-1
	          DO i = 1,2*nmesh
	             g(Ndeep*(j-1)+i,Ndeep*(j-1)+i) = 1.0_dp/(E(k) - (eion) + CMPLX(0.0,ep))
               g(Ndeep*j-i+1,Ndeep*j-i+1) = 1.0_dp/(E(k) - (eion) + CMPLX(0.0,ep))
            ENDDO
         ENDDO

         DO j = 1,nmesh
            DO i = 1,nmesh
               g(i + nmesh*2 + Ndeep*j,i + nmesh*2 + Ndeep*j) = 1.0_dp/(E(k) - eion + CMPLX(0.0,ep))
               g(i + nmesh*4 + Ndeep*j,i + nmesh*4 + Ndeep*j) = 1.0_dp/(E(k) - eion + CMPLX(0.0,ep))
               g(i + 3*nmesh*Ndeep + 2*nmesh + Ndeep*j,i + 3*nmesh*Ndeep + 2*nmesh + Ndeep*j) = 1.0_dp/(E(k) - eion + CMPLX(0.0,ep))
               g(i + 3*nmesh*Ndeep + 4*nmesh + Ndeep*j,i + 3*nmesh*Ndeep + 4*nmesh + Ndeep*j) = 1.0_dp/(E(k) - eion + CMPLX(0.0,ep))
            ENDDO
        ENDDO

      ENDIF

   ! Constructing g for the leads
      DO m = 1, Ndeep
         DO n = 1,Ndeep
            Gzero = ZERO
            DO p = 1, Ndeep
               l = E(k) + 2.0_dp*ABS(gamm)*COS((REAL(p,dp))*PI*ndone)
               h = l**2.0_dp - 4.0_dp*(ABS(gamm)**2.0_dp)
               IF ( h > 0.0_dp .AND. l < 0.0_dp) THEN
                  Gzero = SIN((REAL(p*m)*PI*ndone))*SIN((REAL(p*n))*PI*ndone)*(l+SQRT(h+ZERO)) + Gzero
               ELSE
                  Gzero = SIN((REAL(p*m))*PI*ndone)*SIN((REAL(p*n))*ndone*PI)*(l-SQRT(h+ZERO)) + Gzero
               END IF
            END DO
         g(m,n) = ndone*Gzero/(gamm**2.0_dp)
         gdag(n,m) = CONJG(g(m,n))
         g(m + Natoms - Ndeep, n + Natoms - Ndeep) = ndone*Gzero/(gamm**2.0_dp)
         gdag(n + Natoms - Ndeep,m + Natoms - Ndeep) = CONJG(g(m + Natoms - Ndeep, n + Natoms - Ndeep))
         END DO
      END DO

    !  DLO   = ZERO
    !  alpha = CMPLX(1.0_dp,0.0_dp)
    !  beta  = ZERO
    !  CALL zmm(alpha, Natoms, PL,'N',g,'N', DLO,beta)
    !  alpha = CMPLX(0.0_dp,-0.5_dp/PI)
    !  beta  = CMPLX(0.0_dp,0.5_dp/PI)
    !  CALL zmm(alpha,Natoms,PL,'N',g,'C',DLO,beta)
    !  DLO  = CMPLX(0.0_dp,-0.5_dp/PI)*MATMUL(PL,(TRANSPOSE(CONJG(g))-g))


      CALL date_and_time( date, time )


      IF (communicator%IPROCESSOR.EQ.0) THEN
         WRITE(control_var%output_file, *)
         WRITE(control_var%output_file, *)                                     &
            'Compute DL0 and DRO',                                             &
            date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                 &
            time(1:2), ':', time(3:4), ':', time(5:6)
         WRITE(control_var%output_file, *)
      ENDIF

      DRO  = ZERO
      DLO  = ZERO

      DO i = 1,Ndeep
         DO j = 1,Ndeep
          DLO(i,j) = CMPLX(0.0_dp,-0.5/PI)*(g(i,j)-gdag(i,j))
          DRO(i,j) = CMPLX(0.0_dp,-0.5/PI)*(g(NATOMS-Ndeep+i,NATOMS-Ndeep+j)-gdag(NATOMS-Ndeep+i,NATOMS-Ndeep+j))
         ENDDO
      ENDDO

      CALL date_and_time( date, time )
            IF (communicator%IPROCESSOR.EQ.0) THEN
               WRITE(control_var%output_file, *)
               WRITE(control_var%output_file, *)                               &
                  'Compute DL0 and DRO start intmed           ',               &
                  date(7:8), '/', date(5:6), '/', date(1:4), ' at ',           &
                  time(1:2), ':', time(3:4), ':', time(5:6)
               WRITE(control_var%output_file, *)
    ENDIF



      CALL date_and_time( date, time )
      IF (communicator%IPROCESSOR.EQ.0) THEN
         WRITE(control_var%output_file, *)
         WRITE(control_var%output_file, *)                                     &
            'fin dr0 start intmed           ',                                 &
            date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                 &
            time(1:2), ':', time(3:4), ':', time(5:6)
         WRITE(control_var%output_file, *)
      ENDIF


      alpha = CMPLX(1.0_dp,0.0_dp)
      beta = CMPLX(0.0_dp,0.0_dp)

      CALL zmm(alpha, Natoms,g,'N',V,'N',INTMED,beta)
      INTMED = CEYE - INTMED
      CALL date_and_time( date, time )
      IF (communicator%IPROCESSOR.EQ.0) THEN
         WRITE(control_var%output_file, *)
	       WRITE(control_var%output_file, *)                                     &
         'start inversion     fin intmed                ',                     &
	       date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                    &
	       time(1:2), ':', time(3:4), ':', time(5:6)
         WRITE(control_var%output_file, *)
      ENDIF

      CALL inv_matrix(GBig,INTMED, NATOMS )

      CALL date_and_time( date, time )
      IF (communicator%IPROCESSOR.EQ.0) THEN
         WRITE(control_var%output_file, *)
	       WRITE(control_var%output_file, *)                                     &
         'fin inversion     Gbig started              ',                       &
	       date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                    &
	       time(1:2), ':', time(3:4), ':', time(5:6)
         WRITE(control_var%output_file, *)
      ENDIF


        alpha = CMPLX(1.0_dp,0.0_dp)
        beta =  CMPLX(0.0_dp,0.0_dp)
        CALL zmm(alpha,Natoms,Gbig,'N',g,'N',Gbig,beta)

        CALL date_and_time( date, time )
        IF (communicator%IPROCESSOR.EQ.0) THEN
           WRITE(control_var%output_file, *)
  	       WRITE(control_var%output_file, *)                                   &
           'Gbig fin     DL started             ',                             &
  	       date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                  &
  	       time(1:2), ':', time(3:4), ':', time(5:6)
           WRITE(control_var%output_file, *)
        ENDIF






        alpha = CMPLX(1.0_dp,0.0_dp)
        beta =  CMPLX(0.0_dp,0.0_dp)
        CALL zmm(alpha,Natoms,Gbig,'N',V,'N',DL,beta)
        DL = CEYE + DL
        alpha = CMPLX(1.0_dp,0.0_dp)
        beta =  CMPLX(0.0_dp,0.0_dp)
        CALL zmm(alpha,Natoms,DL,'N',DLO,'N',DL,beta)

        alpha = CMPLX(1.0_dp,0.0_dp)
        beta =  CMPLX(0.0_dp,0.0_dp)


        CALL zmm(alpha,Natoms,V,'N',Gbig,'C',DR,beta)
        CALL zmm(alpha,Natoms,V,'N',Gbig,'C',DR,beta)


        DL = MATMUL(DL, (CEYE + MATMUL(V, TRANSPOSE(CONJG(Gbig)))))

        CALL date_and_time( date, time )
        IF (communicator%IPROCESSOR.EQ.0) THEN
           WRITE(control_var%output_file, *)
  	       WRITE(control_var%output_file, *)                                   &
           'DL fin     DR started             ',                               &
  	       date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                  &
  	       time(1:2), ':', time(3:4), ':', time(5:6)
           WRITE(control_var%output_file, *)
        ENDIF
        DR = CEYE + MATMUL(Gbig, V)
        DR = MATMUL(DR, DRO)
        DR = MATMUL(DR, (CEYE + MATMUL(V, TRANSPOSE(CONJG(Gbig)))))


        DOS  = DR + DL


        Current = 2.0_dp*W*gamm*AIMAG(DL)
        RDOS = REAL(DOS,dp)

        CALL date_and_time( date, time )

        IF (communicator%IPROCESSOR.EQ.0) THEN

           WRITE(control_var%output_file, *)

           WRITE(control_var%output_file, *)                                     &
              'start calculating u, v fcurlz                     ',              &
              date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                 &
              time(1:2), ':', time(3:4), ':', time(5:6)

           WRITE(control_var%output_file, *)

        ENDIF





        DO i = 1, Ndeep
          ux(i,1) = Current(i,i+Ndeep)
          ux(i,Nlayers) = Current(Natoms- Ndeep - Ndeep +i, NATOMS - Ndeep+i)
          fx(i,1) = RDOS(i,i)*ux(i,1)
          fx(i,Nlayers) = RDOS(Natoms-Ndeep+i,Natoms-Ndeep+i)*ux(i,1)
        END DO


        DO j = 2,Nlayers-1
          DO i = 1,Ndeep
            ux(i,j) = Current( Ndeep*(j-2)+i,Ndeep*(j-1) + i ) + Current(Ndeep*(j-1)+i,Ndeep*(j) +i)
            fx(i,j) = RDOS( Ndeep*(j-1) +i, Ndeep*(j-1) + i)*ux(i,j)
          END DO
        END DO

        DO i = 1,NLayers
          vy(1,i) = 0.5*Current((i-1)*Ndeep+2,(i-1)*Ndeep+1)
          vy(Ndeep,i) = 0.5*Current(i*Ndeep,i*Ndeep-1)
          fy(1,i) = vy(1,i)*RDOS(1+(i-1)*Ndeep,1+(i-1)*Ndeep)
          fy(Ndeep,i) = vy(Ndeep,i)*RDOS(i*Ndeep,i*Ndeep)
        END DO


        DO i = 1,NLayers
            DO j = 2,Ndeep-1
            vy(j,i) = 0.5*(Current(j+(i-1)*Ndeep,j-1+(i-1)*Ndeep)+Current(j+1+(i-1)*Ndeep,j+(i-1)*Ndeep))
            fy(j,i) = vy(j,i)*RDOS(j+Ndeep*(i-1),j+Ndeep*(i-1))
            END DO
        END DO

        Do i = 2,NLayers-1
          Fcurlz(1,i) = 0.5_dp*(fy(1,i-1)-fy(1,i+1) -fx(2,i))/alattice
          Fcurlz(Ndeep,i) = 0.5_dp*(fy(Ndeep,i-1)-fy(Ndeep,i+1) +fx(Ndeep-1,i))/alattice
        END DO

        Do i = 2,Nlayers-1
          DO j = 2,Ndeep-1
            Fcurlz(j,i) = 0.5_dp*(fy(j,i-1)-fy(j,i+1) + fx(j-1,i)-fx(j+1,i))/alattice
          END DO
        END DO






      CALL date_and_time( date, time )

      IF (communicator%IPROCESSOR.EQ.0) THEN

         WRITE(control_var%output_file, *)

         WRITE(control_var%output_file, *)                                     &
            'finished calculating u, v fcurlz  writing out fcurls                   ',                             &
            date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                 &
            time(1:2), ':', time(3:4), ':', time(5:6)

         WRITE(control_var%output_file, *)

      ENDIF

      IF (communicator%iprocessor.EQ.0) THEN

         WRITE(cnum, '(I3.3)') k
         filename = communicator%cprocessor // '/outputcut.' // cnum // '.dat'

         OPEN(10, file = filename, form = 'unformatted')

        WRITE(10) REAL(ux,dp)
        WRITE(10) REAL(vy,dp)
        WRITE(10) REAL(Fcurlz,dp)
        WRITE(10) E(k)
        CLOSE(unit = 10)

      ENDIF

      CALL date_and_time( date, time )
      IF (communicator%IPROCESSOR.EQ.0) THEN
         WRITE(control_var%output_file, *)
         WRITE(control_var%output_file, *)                                     &
         'finished outputting file             ',                            &
         date(7:8), '/', date(5:6), '/', date(1:4), ' at ',                 &
         time(1:2), ':', time(3:4), ':', time(5:6)
         WRITE(control_var%output_file, *)
      ENDIF
STOP
     END DO


!     DO i = 1, Ndeep

!        WRITE(*, '(24(1X, F12.6))') (REAL(ux(i,j)), j = 1, Nlayers)

!     END DO

    deallocate(g)
    deallocate(V)
    deallocate(PL)
    deallocate(PR)
    deallocate(E)
    deallocate(DLO)
    deallocate(DRO)
    deallocate(DL)
    deallocate(DR)
    deallocate(Gbig)
    deallocate(INVINTMED)
    deallocate(INTMED)
    deallocate(CEYE)
    deallocate(Current)
    deallocate(DOS)
    deallocate(ux)
    deallocate(vy)
    deallocate(fx)
    deallocate(fy)
    deallocate(Fcurlz)
    deallocate(RDOS)
end program conductor
