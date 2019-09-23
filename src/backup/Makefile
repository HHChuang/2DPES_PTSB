FC = ifort
FCFLAGS = -traceback -debug all -O3
LIBPATH = -Wl,  /usr/local/lib/liblapack.a /usr/local/lib/libblas.a 

OverlapStruc:	OverlapStruc.f90
	$(FC) -o OverlapStruc $(FCFLAGS) OverlapStruc.f90 $(LIBPATH)

clean:	
	rm -rf OverlapStruc OverlapStruc.dSYM constants.mod 
