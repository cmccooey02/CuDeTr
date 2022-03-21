F90         =	mpif90 -openmp -I${MKLROOT}/include
F77         =	mpif90 -openmp -I${MKLROOT}/include

64BIT       =	
64BITAR     =	

OPT         =    -O3 -heap-arrays

F90FLAGS    =	$(OPT) $(64BIT)
FD90FLAGS   =	$(OPT) $(64BIT)
F90FLAGSFDF =	$(OPT) $(64BIT)

F77FLAGS    =	$(OPT) $(64BIT)
F77FLAGSFDF =	$(OPT) $(64BIT)

LINKFLAGS   =	$(OPT) $(64BIT)

LIBS	=	-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread \
		-lmkl_core -liomp5 -lpthread -lm -ldl

ARFLAGS     =   $(64BITAR) qv

DIRECTIVES  = -D_COM_MPI=1,-D_USE_BLAS=1

COMMTYPE    =	mpi



# -----------------------------------------------------------------------------#
#                                                                              #
# How to produce objects from source files                                     #
#                                                                              #
# -----------------------------------------------------------------------------#



# Fortran90

%.o:%.f90
	$(F90) $(F90FLAGS) -c $<



# Fortran77

%.o:%.f
	$(F77) $(F77FLAGS) -c $<



# Pre-Processed Fortran90

%.o:%.F90
	$(F90) $(FD90FLAGS) $(DIRECTIVES) -c $<



# -----------------------------------------------------------------------------#
#                                                                              #
# End of architecture Makefile                                                 #
#                                                                              #
# -----------------------------------------------------------------------------#
