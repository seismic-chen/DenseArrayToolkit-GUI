FFLAGS = -g -C
SACLIB = /usr/local/seis/sac/lib/libsacio.a
OBJECTS = iterdecon.o
#LIBS = /geo/NumRec/fnrlib.a $(SACLIB)
LIBS = recipes00.o $(SACLIB)

iterdeconfd: iterdeconfd.o $(LIBS)
	f77 $(FFLAGS) -o iterdecon iterdeconfd.o $(LIBS)

iterdeconfd.o: iterdeconfd.f
recipes00.o: recipes00.f

export:
	mkdir IterDecon
	cp iterdeconfd.f IterDecon
	cp recipes00.f IterDecon
	cp Makefile IterDecon
	cp -r TEST01 IterDecon
	tar cvf iterdecon.tar ./IterDecon
	compress iterdecon.tar
	csh install.it.csh
	/bin/rm -r IterDecon

clean:
	/bin/rm residual* d0 out.predict p1 numerator denominator
