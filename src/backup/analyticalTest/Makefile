FC = ifort
FCFLAGS = -traceback -debug -O3
LIBPATH = -Wl,  /usr/local/lib/liblapack.a /usr/local/lib/libblas.a 

analytic2DPES:	analytic2DPES.f90
	$(FC) -o analytic2DPES $(FCFLAGS) analytic2DPES.f90 $(LIBPATH)

clean:	
	rm -rf analytic2DPES.dSYM analytic2DPES