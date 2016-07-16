project=main.f90
lib =markov_module.f90 
F90 =ifort
#opts= -check all -CB -warn all -r8 
opts=  -CB -r8  -fast
out =main
MKLROOT=/opt/intel/mkl
MKLPATH=${MKLROOT}/lib/intel64
INCLUDE=-I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include

#LIBRARY=  -lmkl_sequential  -lmkl_core

#LIBRARY= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
LIBRARY=-Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_sequential.a ${MKLPATH}/libmkl_core.a -Wl,--end-group -lpthread -lm


main: $(project) lib
#	$(F90) $(project) $(lib)  $(opts)   -o $(out)
	$(F90) $(project) $(lib)  $(opts) $(INCLUDE) $(LIBRARY)    -o $(out)
lib: $(lib)
	$(F90) -c $(lib)  $(opts) 
clean:
	rm -f *.dat fort.*

