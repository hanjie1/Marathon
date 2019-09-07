# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib -lpacklib -lc -lm
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

BinCenter_objs = BinCenter.o f2glob.o ineft.o 
BinCenter_srcs = BinCenter.f f2glob.f ineft.f 

########################################

none: BinCenter

all:  BinCenter

BinCenter.o: BinCenter.f
		 $(F77) $(FFLAGS) -c $< -o $@

f2glob.o: f2glob.f
	$(F77) $(FFLAGS) -c $< -o $@

ineft.o: ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

BinCenter: $(BinCenter_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(BinCenter_objs) $(OTHERLIBS)


clean:
	rm -f *.o BinCenter
