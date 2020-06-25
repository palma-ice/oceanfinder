.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = .
bindir = .

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

FC = gfortran
INC_NC  = -I/opt/local/include
LIB_NC  = -L/opt/local/lib -lnetcdff -L/opt/local/lib -Wl,-headerpad_max_install_names -Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk -arch x86_64 -lnetcdf -lnetcdf -lm 

FFLAGS = -I$(objdir) -J$(objdir) $(INC_NC) -m64 -ffree-line-length-none 
LFLAGS = $(LIB_NC)
DFLAGS = -O2

###############################################
##
## Compilation of complete programs
##
###############################################

# Test programs that use yelmo-static
test : ncio.f90 test.f90 
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test.x ncio.f90 test.f90 $(LFLAGS)
		@echo " "
		@echo "    test.x is ready."
		@echo " "

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
