include lobster.make


TARGET = conductor_$(COMMTYPE).exe

default: $(TARGET)

PROGRAM =		conductor.o

MODULES =		constants.o communications.o iomodule.o linalgebra.o   \
			util.o green.o density.o

OBJS =			$(PROGRAM) $(MODULES)

FDF =			lib/libfdf.a

# module dependencies


communications.o:	constants.o util.o
iodmoule.o:		communications.o constants.o  util.o
linalgebra.o:		constants.o
util.o:			constants.o
green.o:	        constants.o linalgebra.o
density.o:              constants.o linalgebra.o

$(MODULES) :		$(FDF)

$(PROGRAM) :		$(MODULES) $(FDF)

$(FDF):
	(cd fdf; $(MAKE) "F90=$(F90)" "F90FLAGS=$(F90FLAGSFDF)" "F77=$(F77)"   \
	 "F77FLAGS=$(F77FLAGSFDF)" "ARFLAGS=$(ARFLAGS)" module)

$(TARGET) :		$(MODULES) $(PROGRAM) $(FDF)
			$(F90) -o bin/$(TARGET) $(LINKFLAGS) $(OBJS) $(FDF) $(LIBS)

clean:
			rm -f $(TARGET) $(OBJS) $(FDF) *.o *.oo *.mod
			(cd fdf; make clean)

clearbin:
			rm -f bin/*.dat
run:
		        ./$(TARGET)
# DO NOT DELETE
