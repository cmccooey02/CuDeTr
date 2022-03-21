# 
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2003
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#  Makefile for fdf example and fdf module.
#
CFLAGS= 
LDFLAGS= $(CFLAGS)
RANLIB = ranlib
#
#
module:  libfdf.a
libfdf.a:	fdf.o fdf_mod.o parse.o  io.o
	ar $(ARFLAGS) libfdf.a fdf.o fdf_mod.o parse.o io.o
	$(RANLIB) libfdf.a
	cp libfdf.a ../lib/
#
# Note trick for Intel Fortran Compiler
#
fdf_mod.o: fdf_mod.f90 parse.o
	$(F90) $(F90FLAGS) -c fdf_mod.f90
	@cp -p *.*d ..
	@echo "fdf/work.pc" >> ../work.pcl
fdf.o:		fdf.h
#
#---------------------------------------------------
OBJS=	sample.o fdf.o io.o
#
sample: $(OBJS)
	$(FC) $(LDFLAGS) -o sample $(OBJS)
#
check:	sample sample.fdf
	sample < sample.fdf
#
io_sample: io_sample.o io.o
	$(FC) $(LDFLAGS) -o io_sample io_sample.o io.o
#
io_check: io_sample
	io_sample
clean:
	rm -f *.o libfdf.a *.mod fdf_mod.o io_sample.o sample io_sample
	rm -f OUTPUT FDF.debug file2 file3 sample.out
	rm -f *.pc *.pcl
#
sample.o:	fdfdefs.h

%.o:%.f90
	$(F90) $(F90FLAGS) -c $<

%.o:%.f  
	$(F77) $(F77FLAGS) -c $<



