
program conductor

   USE constants
   USE communications
   USE fdf
   USE iomodule
   USE linear_algebra
   USE util
   USE green_function
   USE density
   implicit none

!  Sample constants
   INTEGER :: SNx, SURN, SNatomsx, SNatomsy, nmesh
   INTEGER :: NE, NEproc

   INTEGER :: ierror
!  Indices
   INTEGER :: i, k , m, n, p, j
   INTEGER :: kstart
!  Cut Constants & Cut Arrays
   INTEGER :: iblk, ik, numcut, cutstart, inner_cutstart, inner_numcut, eppow, RandomOnsiteCutStart, RandomOnsiteNumCuts
   INTEGER :: botcountL, botcountR, topcountL, topcountR, innercountL, innercountR
   INTEGER, ALLOCATABLE :: topcut(:), botcut(:), innercut(:)  ! Cut array
   REAL(dp), ALLOCATABLE ::  RandomOnsiteCuts(:)
   REAL(dp) ::  inner_Energy, E_start, E_finish

! Files
   CHARACTER(LEN=3) :: cnum, cnumb
   CHARACTER(LEN=100) :: filename, fileid, path
   CHARACTER(LEN=2) :: SIMTYPE

!
   REAL(dp):: ep, alattice, gamm, dE, eion, gamm_2
   REAL(dp) :: start, finish, estart, efinish, gstart, gfinish
   REAL(dp), ALLOCATABLE     :: E(:), Energy_ev(:), DeBroglie(:), Conductance(:), red(:), blue(:)
   REAL(dp) :: EHigh, eon, edif, ndone, l, h, Eo
   REAL(dp) :: arg, W
   REAL(dp) :: Jin, Jout, workJ


   ! Test variable
   REAL(dp) :: remainder, rcount

   COMPLEX(dp)                :: Gzero, alpha, beta
   COMPLEX(dp), ALLOCATABLE   :: GLsurface(:,:,:),  gLedge(:,:), GBigLedge(:,:), GCON(:,:)
   COMPLEX(dp), ALLOCATABLE   :: WORKA(:,:), WORKB(:,:), WORKC(:,:), WORKD(:,:), WORKE(:,:)
   COMPLEX(dp), ALLOCATABLE   :: GRsurface(:,:,:),  gRedge(:,:), GBigRedge(:,:), GBigDagLedge(:,:)

   COMPLEX(dp), ALLOCATABLE   :: ident(:,:), identcount
   REAL(dp)                   :: conditionnumber
   COMPLEX(dp), ALLOCATABLE   :: VL(:,:), VR(:,:), VCON(:,:)
   COMPLEX(dp), ALLOCATABLE   :: PLVPR(:,:), PRVPL(:,:), PL(:,:), PR(:,:)
   COMPLEX(dp), ALLOCATABLE   :: BigG(:,:)
   COMPLEX(dp), ALLOCATABLE   :: DLedge(:,:), DRedge(:,:)
   REAL(dp),    ALLOCATABLE   :: Current(:,:,:), JX(:,:), JY(:,:), FX(:,:), FY(:,:), CURL(:,:), DOS(:,:)
   REAL(dp),    ALLOCATABLE   :: DOSL(:,:), JLength(:,:),  DSurface(:,:), DLSurface(:,:)
   COMPLEX(dp), ALLOCATABLE   :: CEYE(:,:), INTMED(:,:)


   CALL initialize_communications( )
   CALL output_communications( )
   CALL check_processor_count( )
   CALL read_simulation_parameters( )
   E_start = fdf_double('E_START',0.05_dp)
   E_finish = fdf_double('E_FINISH',0.25_dp)
   SNatomsx = fdf_integer( 'SNx', 3)
   SNatomsy = fdf_integer( 'SNy', 3)
   numcut  = fdf_integer( 'NUMCUT', 0 )
   cutstart = fdf_integer('CUTSTART', 10)
   NEproc  = fdf_integer('E_NUM', 10)
   NE  =  NEproc
   inner_numcut = fdf_integer('INNERNUMCUT', 0)
   RandomOnsiteCutStart = fdf_integer('RandomOnsiteCutStart', 1)
   RandomOnsiteNumCuts = fdf_integer('RandomOnsiteNumCuts', 0)
   inner_cutstart = fdf_integer('INNERCUTSTART',1)
   kstart = fdf_integer('KSTART',1)
   ep = fdf_double('ep',1.0E-10_dp)


   SIMTYPE = fdf_string('SIMTYPE', 'DH')
   EHigh = fdf_double('EHIGH_CUT', 2000.0_dp)





   CALL cpu_time(start)


  IF (numcut /= 0) THEN

     ALLOCATE(topcut(1:numcut))
     IF (fdf_block( 'TopCutValues', iblk )) THEN
        DO ik = 1,  numcut
            READ(iblk, *) topcut(ik)
        ENDDO
     ELSE
        CALL stopping( "No TOPcut values" )
     ENDIF


     ALLOCATE(botcut(1:numcut))
     IF (fdf_block( 'BotCutValues', iblk )) THEN
        DO ik = 1,  numcut
           READ(iblk, *) botcut(ik)
        ENDDO
     ELSE
        CALL stopping( "No BOTcut values" )
     ENDIF


     PRINT *, ' READING IN INNER CUTS'

     inner_Energy = fdf_double('EHIGH_INNER', 0.0_dp)


     ALLOCATE(innercut(1:2*inner_numcut))

      IF (fdf_block( 'InnerCutValues', iblk )) THEN
         DO ik = 1,  2*inner_numcut
            READ(iblk, *) innercut(ik)
         ENDDO
      ELSE
         CALL stopping( "No InnerCut Values" )
      ENDIF
      IF (communicator%iprocessor.EQ.0) THEN
         filename = communicator%cprocessor // '/processinginputs/topcut000.dat'
         OPEN(13, file = filename, form = 'unformatted')
         WRITE(13) topcut
         CLOSE(unit = 13)
      ENDIF
      IF (communicator%iprocessor.EQ.0) THEN
         filename = communicator%cprocessor // '/processinginputs/botcut000.dat'
         OPEN(14, file = filename, form = 'unformatted')
         WRITE(14) botcut
         CLOSE(unit = 14)
      ENDIF
  ENDIF

  IF (RandomOnsiteNumCuts /= 0) THEN

     PRINT *, ' READING IN RANDOM ONSITE CUTS'

     ALLOCATE(RandomOnsiteCuts(1:RandomOnsiteNumCuts))

      IF (fdf_block( 'RandomOnsiteCuts', iblk )) THEN
         DO ik = 1,  RandomOnsiteNumCuts
            READ(iblk, *) RandomOnsiteCuts(ik)
         ENDDO
      ELSE
         CALL stopping( "No Random OnSite Values" )
      ENDIF
      IF (communicator%iprocessor.EQ.0) THEN
         filename = communicator%cprocessor // '/processinginputs/randcuts.dat'
         OPEN(17, file = filename, form = 'unformatted')
         WRITE(17) RandomOnsiteNumCuts
         WRITE(17) RandomOnsiteCuts
         CLOSE(unit = 17)
      ENDIF
  ENDIF






   SNx = SNatomsx
   NLayers = (2 + SNx)
   Ndeep = SNatomsy
   Natoms = Nlayers*Ndeep
   SURN = 2*Ndeep
   PRINT *, 'epsilon', ep
   W = 1.0_dp






   IF (SIMTYPE.EQ.'TB') THEN
     alattice = fdf_double('TBALATTICE',4.9_dp)   ! Bohr Units
     gamm = fdf_double('TBGAMM',0.0_dp) ! Hartree
     PRINT *, ' TB'
   ENDIF
   IF (SIMTYPE.EQ.'DH') THEN
     nmesh = fdf_integer('NUMMESH',1)
     alattice = 7.5_dp/REAL(nmesh,dp)   ! Bohr Units
     gamm = -0.5_dp/(alattice**2) ! Hartree
     PRINT *, ' Discretised H'
   ENDIF







   IF (communicator%iprocessor.EQ.0) THEN

      filename = communicator%cprocessor // '/processinginputs/outputconstantscut.dat'
      OPEN(13, file = filename, form = 'unformatted')
      WRITE(13) SNx
      WRITE(13) Ndeep
      WRITE(13) NE
      WRITE(13) alattice
      WRITE(13) nmesh
      WRITE(13) SNatomsx
      WRITE(13) SNatomsy
      WRITE(13) E_start
      WRITE(13) E_finish
      CLOSE(unit = 13)
   ENDIF
   IF (communicator%iprocessor.EQ.0) THEN
      filename = communicator%cprocessor // '/processinginputs/cut' // communicator%cprocessor // '.dat'
      OPEN(15, file = filename, form = 'formatted')
      WRITE(15,*) 'SNx = ', SNx
      WRITE(15,*) 'Ndeep = ', Ndeep
      WRITE(15,*) 'allattice = ', alattice
      WRITE(15,*) 'gamm = ', gamm
      WRITE(15,*) 'SNx = ', SNatomsx
      WRITE(15,*) 'SNx = ', SNatomsy
      WRITE(15,*) 'Estart = ', E_start
      WRITE(15,*) 'Efinish = ', E_finish
      CLOSE(unit = 15)
   ENDIF

   ALLOCATE(E(1:NEproc))
   ALLOCATE(Energy_ev(1:NEproc))
   ALLOCATE(red(1:NEproc))
   ALLOCATE(blue(1:NEproc))
   ALLOCATE(Conductance(1:NEproc))
   ALLOCATE(DeBroglie(1:NEproc))
   ALLOCATE(GLsurface(1:Ndeep,1:Ndeep, 1:SNx+1))
   ALLOCATE(GRsurface(1:Ndeep,1:Ndeep, 1:SNx+1))
   ALLOCATE(WORKE(1:Ndeep,1:Ndeep))
   ALLOCATE(gLedge(1:SURN,1:SURN))
   ALLOCATE(gRedge(1:SURN,1:SURN))
   ALLOCATE(GBigLedge(1:SURN,1:SURN))
   ALLOCATE(GBigRedge(1:SURN,1:SURN))

   ALLOCATE(ident(1:SURN,1:SURN))

   ALLOCATE(GCON(1:SURN,1:SURN))
   ALLOCATE(GBigDagLedge(1:SURN,1:SURN))
   ALLOCATE(BigG(1:SURN,1:SURN))
   ALLOCATE(DSurface(1:SURN,1:Nlayers))

   ALLOCATE(DLSurface(1:SURN,1:Nlayers))

   ALLOCATE(Current(1:SURN,1:SURN,1:Nlayers))
   ALLOCATE(JX(1:Ndeep,1:NLayers))
   ALLOCATE(JY(1:Ndeep,1:NLayers))
   ALLOCATE(FX(1:Ndeep,1:NLayers))
   ALLOCATE(FY(1:Ndeep,1:NLayers))
   ALLOCATE(JLength(1:Ndeep,1:NLayers))
   ALLOCATE(CURL(1:Ndeep,1:NLayers))
   ALLOCATE(DOS(1:Ndeep,1:NLayers))
   ALLOCATE(DOSL(1:Ndeep,1:NLayers))
   ALLOCATE(DLedge(1:SURN,1:SURN))
   ALLOCATE(DRedge(1:SURN,1:SURN))
   ALLOCATE(WORKA(1:SURN,1:SURN))
   ALLOCATE(WORKB(1:SURN,1:SURN))
   ALLOCATE(WORKC(1:SURN,1:SURN))
   ALLOCATE(WORKD(1:SURN,1:SURN))
   ALLOCATE(VL(1:SURN,1:SURN))
   ALLOCATE(VR(1:SURN,1:SURN))
   ALLOCATE(VCON(1:SURN,1:SURN))
   ALLOCATE(PLVPR(1:SURN,1:SURN))
   ALLOCATE(PRVPL(1:SURN,1:SURN))
   ALLOCATE(PL(1:SURN,1:SURN))
   ALLOCATE(PR(1:SURN,1:SURN))
   ALLOCATE(CEYE(1:SURN,1:SURN))
   ALLOCATE(INTMED(1:SURN,1:SURN))


   Eo = 0.0_dp

   eon = -4.0_dp*ABS(gamm)
   edif  = eon + Eo + E_start

   dE = ABS(E_finish-E_start)/REAL(NE,dp)

   ndone = 1.0_dp/(REAL(Ndeep+1,dp))

   DO i = 1, NEproc
      E(i) = edif + i * dE
   ENDDO

   DO i = 1,NEproc
   DeBroglie(i) = PI*SQRT(2.0_dp/(ABS(E(i)- (Eo+eon))))
   ENDDO

   VL= ZERO
   VR= ZERO
   VCON = ZERO
   CEYE = ZERO

   DO i = 1,Ndeep
      VL(i,Ndeep+i) = gamm
      VL(Ndeep+i,i) = gamm
      VR(i,Ndeep+i) = gamm
      VR(Ndeep+i,i) = gamm
      VCON(i,Ndeep+i) = gamm
      VCON(Ndeep+i,i) = gamm
   ENDDO

   DO i = 1,Ndeep-1
      VL(Ndeep+i,Ndeep+i+1) = gamm
      VL(Ndeep+i+1,Ndeep+i) = gamm
      VR(i,i+1) = gamm
      VR(i+1,i) = gamm
   ENDDO

   DO i = 1,SURN
      CEYE(i,i) = DCMPLX(1.0_dp,0.0_dp)
   ENDDO

   PLVPR(:,:) = ZERO
   PRVPL(:,:) = ZERO
   PL(:,:) = ZERO
   PR(:,:) = ZERO

   DO i = 1,Ndeep
      PLVPR(Ndeep+i,i) = gamm
      PLVPR(i, Ndeep+i) = -gamm
      PRVPL(Ndeep+i,i) = -gamm
      PRVPL(i, Ndeep+i) = gamm
      PL(i,i) = DCMPLX(1.0_dp,0.0_dp)
      PR(i+Ndeep,i+Ndeep) = DCMPLX(1.0_dp,0.0_dp)
   ENDDO




   red(:) = 0.0_dp
   blue(:) = 0.0_dp
   IF (communicator%iprocessor.EQ.0) THEN
      filename = communicator%cprocessor // '/processinginputs/calcfactors' // communicator%cprocessor // '.dat'
      OPEN(16, file = filename, form = 'unformatted')
      WRITE(16) SNx
      WRITE(16) Ndeep
      WRITE(16) NE
      WRITE(16) alattice
      WRITE(16) nmesh
      WRITE(16) SNatomsx
      WRITE(16) SNatomsy
      WRITE(16) gamm
      WRITE(16) REAL(E,dp)

      CLOSE(unit = 16)
   ENDIF

   DO k = kstart, NEproc
      CALL cpu_time(estart)

!!! Setting up the surface greens function for the leads

      CALL GREEN_sur_lead(WORKE,Ndeep,gamm, E(k))


!!! Assigning WORKE to the left and right lead

      DO i = 1,Ndeep
         DO j = 1,Ndeep
            GLsurface(i,j,1) = WORKE(i,j)
            GRsurface(i,j,SNx+1) = WORKE(i,j)
         ENDDO
      ENDDO

      topcountL = 0
      topcountR = NumCut
      botcountL = 0
      botcountR = NumCut
      innercountL = 0
      innercountR = 2*inner_numcut

      DO i = 1,SNx
         gLedge(:,:) = ZERO
         gRedge(:,:) = ZERO

         DO n = 1,Ndeep
            DO m = 1,Ndeep
            gLedge(m,n) = GLsurface(m,n,i)
            gRedge(m+Ndeep,n+Ndeep) = GRsurface(m,n,Nlayers-i)
            ENDDO
         ENDDO

         DO m = 1,Ndeep
            gLedge(m+Ndeep,m+Ndeep) = 1.0_dp/(E(k) - Eo + DCMPLX(0.0_dp,ep))
            gRedge(m,m) = 1.0_dp/(E(k) - Eo + DCMPLX(0.0_dp,ep))
         ENDDO



         IF (numcut /= 0) THEN
!!!! Top CUT
           IF (i >= cutstart .AND. i <= (NumCut+cutstart-1)) THEN

                    topcountL = topcountL + 1
                    m = topcut(topcountL)
                    IF (m /= 0)THEN
                        DO n = 1,m
                         gLedge(n+Ndeep,n+Ndeep) = 1.0_dp/(E(k) - EHigh + DCMPLX(0.0_dp,ep))
                        ENDDO
                    ENDIF


           ENDIF
           IF (i > (SNx-cutstart-NumCut+1) .AND. i <= (SNx-cutstart+1)) THEN

                    m = topcut(topcountR)
                    topcountR = topcountR - 1
                    IF (m /= 0)THEN
                       DO n = 1,m
                          gRedge(n,n) = 1.0_dp/(E(k) - EHigh + DCMPLX(0.0_dp,ep))
                       ENDDO
                     ENDIF
           ENDIF
!!!!! BOTTOM CUT
           IF (i >= cutstart .AND. i <= (NumCut+cutstart-1)) THEN

                 botcountL = botcountL + 1
                 m = botcut(botcountL)
                 IF (m /= 0) THEN
                    DO n = 1,m
                       gLedge(SURN+1-n,SURN+1-n) = 1.0_dp/(E(k) - EHigh + DCMPLX(0.0_dp,ep))
                    ENDDO
                 ENDIF


           ENDIF
           IF (i > (SNx-cutstart-NumCut+1) .AND. i <= (SNx-cutstart+1)) THEN

                 m = botcut(botcountR)
                 botcountR = botcountR - 1
                 IF ( m /= 0) THEN
                   DO n = 1,m
                       gRedge(Ndeep -n + 1,Ndeep + 1-n) = 1.0_dp/(E(k) - EHigh + DCMPLX(0.0_dp,ep))
                   ENDDO
                 ENDIF

           ENDIF
         ENDIF


         !!!!!! Defects

         IF (inner_numcut /= 0) THEN
             IF (i >= inner_cutstart .AND. i <= (inner_numcut+inner_cutstart-1)) THEN
                      innercountL = innercountL + 1
                      m = InnerCut(innercountL)
                      innercountL = innercountL + 1
                      j = InnerCut(innercountL)
                      IF (m /= 0)THEN
                          DO n = m,j
                           gLedge(n+Ndeep,n+Ndeep) = 1.0_dp/(E(k) - inner_Energy + DCMPLX(0.0_dp,ep))
                          ENDDO
                      ENDIF
             ENDIF
             IF (i > (SNx-inner_cutstart-inner_numcut+1) .AND. i <= (SNx-inner_cutstart+1)) THEN
                      m = InnerCut(innercountR)
                      innercountR = innercountR - 1
                      j = InnerCut(innercountR)
                      innercountR = innercountR - 1
                      IF (m /= 0)THEN
                         DO n = j,m
                            gRedge(n,n) = 1.0_dp/(E(k) - inner_Energy + DCMPLX(0.0_dp,ep))
                         ENDDO
                     ENDIF
             ENDIF
         ENDIF

         !!!!!! Random Onsites

         IF (RandomOnsiteNumCuts /= 0) THEN

                 DO n = 1,Ndeep
                     gLedge(n+Ndeep,n+Ndeep) = 1.0_dp/(E(k) - RandomOnsiteCuts((i-1)*Ndeep+n) + DCMPLX(0.0_dp,ep))
                 ENDDO

                 DO n = 1,Ndeep
                            gRedge(n,n) = 1.0_dp/(E(k) - RandomOnsiteCuts(RandomOnsiteNumCuts-Ndeep*i +n) + DCMPLX(0.0_dp,ep))
                 ENDDO

         ENDIF




         CALL GREEN_lin_eq(gLedge, VL, SURN)



         DO n = 1,Ndeep
            DO m = 1,Ndeep
               GLsurface(m,n,i+1) = gLedge(m+Ndeep,n+Ndeep)
            ENDDO
         ENDDO

         CALL GREEN_lin_eq(gRedge, VR, SURN)



         DO n = 1,Ndeep
            DO m = 1,Ndeep
               GRsurface(m,n,Nlayers-i-1) = gRedge(m,n)
            ENDDO
         ENDDO

      ENDDO

!!! Connect both Leads together
      DO j = 1,SNx+1

         gLedge = ZERO

         DO n = 1,Ndeep
            DO m = 1, Ndeep
               gLedge(m,n) = GLsurface(m,n,j)
               gLedge(m+Ndeep,n+Ndeep) = GRsurface(m,n,j)
            ENDDO
         ENDDO



         CALL GREEN_lin_eq(gLedge, VCON, SURN)


         DO m = 1,SURN
           DO n = 1,SURN
              GBigDagLedge(m,n) = DCONJG(gLedge(n,m))
           ENDDO
         ENDDO

  !      IF (j .EQ. 196) THEN
!
!            filename = communicator%cprocessor // '/data/green.dat'
!            OPEN(14, file = filename, form = 'unformatted')
!            WRITE(14) REAL(gLedge,dp)
!            WRITE(14) DIMAG(gLedge)
!            WRITE(14) REAL(VCON,dp)
!            WRITE(14) DIMAG(VCON)
!            WRITE(14) REAL(ident,dp)
!            WRITE(14) REAL(GBigDagLedge,dp)
!            WRITE(14) DIMAG(GBigDagLedge)
!            WRITE(14) REAL(PL,dp)
!            WRITE(14) REAL(gLedge,dp)
!            WRITE(14) DIMAG(gLedge)
!            CLOSE(unit = 14)
!        ENDIF


!!!! Construct DL


         CALL DENSITYLEADS(DLedge, PL, PLVPR , gLedge, GBigDagLedge, SURN)
!!! Change back to plvpr


!!! Construct DR


         CALL DOSfunc(DSurface(:,j),gLedge,SURN)

         DO n = 1,SURN
            DO m = 1,SURN

               Current(m,n,j) = 2.0_dp*W*gamm*DIMAG(DLedge(m,n))
  !             DLSurface(m,j) = REAL(DLedge(m,m),dp)
            ENDDO
         ENDDO


!         IF (j .EQ. 50) THEN

!             filename = communicator%cprocessor // '/data/current196.dat'
!             OPEN(14, file = filename, form = 'unformatted')
!             WRITE(14) DIMAG(DLedge)
!             CLOSE(unit = 14)
!         ENDIF





      ENDDO






      DO i = 3, Nlayers-3
         remainder = 0.0
         rcount = 0.0
         DO n = 1, Ndeep
            DO m = 1, Ndeep
               IF (ABS(Current(m,n,i+1)) > 1.0E-12_dp) THEN
                   remainder = (ABS( Current(m+Ndeep,n+Ndeep,i)-Current(m,n,i+1) ))/ABS(Current(m+Ndeep,n+Ndeep,i)+Current(m,n,i+1)) + remainder
                   rcount = rcount + 1.0
               ENDIF
            ENDDO
         ENDDO
  !       PRINT *, 'bond current', i, remainder, rcount, remainder/rcount
      ENDDO


      !!! Construct x and y components of J and also the curl
      JX(:,:) = 0.0_dp
      JY(:,:) = 0.0_dp
      FX(:,:) = 0.0_dp
      FY(:,:) = 0.0_dp
      CURL(:,:) = 0.0_dp
      DOS(:,:) = 0.0_dp
      DOSL(:,:) = 0.0_dp
      Conductance(k) = 0.0_dp
      DO i = 1,Ndeep
        Conductance(k) = Current(i,Ndeep+i,4) + Conductance(k)
      ENDDO
      DO i = 2,Nlayers-1
         DO m = 1, Ndeep
            JX(m,i) = (Current(m,m+Ndeep,i-1)+Current(m,m+Ndeep,i))/2.0_dp
            FX(m,i) = DSurface(m,i)*JX(m,i)
            DOS(m,i) = REAL(DSurface(m,i),dp)
    !        DOSL(m,i) = DLSurface(m,i)
         ENDDO
      ENDDO

      DO i = 2,Nlayers-1
         DO m = 2, Ndeep-1
            JY(m,i) = (Current(m,m+1,i)+Current(m-1,m,i))/2.0_dp
            FY(m,i) = DSurface(m,i)*JY(m,i)
         ENDDO
      ENDDO
      DO i = 2,Nlayers-1
            JY(1,i) = (Current(1,2,i))/2.0_dp
            FY(1,i) = DSurface(1,i)*JY(1,i)
            JY(Ndeep,i) = (Current(Ndeep-1,Ndeep,i))/2.0_dp
            FY(Ndeep,i) = DSurface(Ndeep,i)*JY(Ndeep,i)
      ENDDO

      DO i = 3,Nlayers-2
        DO m = 2, Ndeep-1
           CURL(m,i) = 0.5_dp*(FY(m,i+1)-FY(m,i-1)+FX(m-1,i)-FX(m+1,i))/alattice

        ENDDO
      ENDDO
      DO i = 2,Nlayers-1
           CURL(1,i) = 0.5_dp*(FY(1,i+1)-FY(1,i-1)-FX(2,i))/alattice
           CURL(Ndeep,i) = 0.5_dp*(FY(Ndeep,i+1)-FY(Ndeep,i-1)+FX(Ndeep-1,i))/alattice
      ENDDO



      DO i = 2,Nlayers-1
         DO m = 1,Ndeep

           IF (CURL(m,i) > 0.0_dp) THEN
              red(k) = CURL(m,i) + red(k)
           ELSE
              blue(k) = CURL(m,i) + blue(k)
           ENDIF
         ENDDO
      ENDDO





      m = k
      WRITE(cnum, '(I3.3)') m
      filename =  './000/data/outputcut.' // cnum // '.dat'

      OPEN(10, file = filename, form = 'unformatted')
      WRITE(10) REAL(JX,dp)
      WRITE(10) REAL(JY,dp)
      WRITE(10) REAL(CURL,dp)
      WRITE(10) REAL(DOS,dp)
      WRITE(10) REAL(E(k),dp)
      WRITE(10) REAL(DeBroglie(k),dp)
      WRITE(10) REAL(Conductance(k),dp)
      WRITE(10) REAL(red(k),dp)
      WRITE(10) REAL(blue(k),dp)
  !    WRITE(10) REAL(DOSL,dp)
      CLOSE(unit = 10)



   ENDDO

   IF (communicator%iprocessor.EQ.0) THEN
      filename = communicator%cprocessor // '/processinginputs/postprocessingfactors.dat'
      OPEN(9, file = filename, form = 'unformatted')
      WRITE(9) CutStart
      WRITE(9) NumCut
      CLOSE(unit = 9)
   ENDIF

   CALL cpu_time(finish)

   finish = finish - start
   IF (communicator%iprocessor.EQ.0) THEN

      PRINT *, 'total time for 000', finish
      finish = finish/NEproc
      PRINT *, 'time per energy for 000', finish

    ENDIF

   !     DO i = 1, Ndeep

   !        WRITE(*, '(24(1X, F12.6))') (REAL(ux(i,j)), j = 1, Nlayers)

   !     END DO


   DEALLOCATE(E)
   DEALLOCATE(DeBroglie)
   DEALLOCATE(GLsurface)
   DEALLOCATE(GRsurface)
   DEALLOCATE(gLedge)
   DEALLOCATE(gRedge)
   DEALLOCATE(GBigLedge)
   DEALLOCATE(GBigRedge)
   DEALLOCATE(GBigDagLedge)
   DEALLOCATE(BigG)
   DEALLOCATE(DSurface)
   DEALLOCATE(DLedge)
   DEALLOCATE(WORKA)
   DEALLOCATE(WORKB)
   DEALLOCATE(WORKC)
   DEALLOCATE(WORKD)
   DEALLOCATE(VL)
   DEALLOCATE(VR)
   DEALLOCATE(VCON)
   DEALLOCATE(PLVPR)
   DEALLOCATE(PRVPL)
   DEALLOCATE(PL)
   DEALLOCATE(PR)
   DEALLOCATE(CEYE)
   DEALLOCATE(INTMED)
  CALL MPI_finalize( ierror )
end program conductor
