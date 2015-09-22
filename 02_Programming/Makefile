#VEC = -no-vec

SRC = nbody.c

nbody: nbody.c
	icc -O3 -ansi-alias -openmp -openmp-report -vec-report=5 $(VEC) -std=gnu99 -debug inline-debug-info -o nbody $< -xHost

nbody.mic: nbody.c
	icc -O3 -ansi-alias -openmp -openmp-report -vec-report=5 $(VEC) -std=gnu99 -debug inline-debug-info -o nbody.mic $< -mmic

clean:
	-rm -rf *~ *.o nbody *.optrpt
