
CC   = gcc
COPT = -g -Wall -O3 -pg #-O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

INCLUDEDIR = -I../../libnum/ -I.

DEF = -DGCC #-DTRUE_ARRAY

CPROPEP_LIBNAME = libcpropep.a
THERMO_LIBNAME  = libthermo.a

THERMO_LIBOBJS  = load.o thermo.o
CPROPEP_LIBOBJS = equilibrium.o print.o performance.o derivative.o

.SUFFIXES: .c

all: $(CPROPEP_LIBNAME) $(THERMO_LIBNAME)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(CPROPEP_LIBNAME): $(CPROPEP_LIBOBJS)
	ar -r $@ $(CPROPEP_LIBOBJS)
	ranlib $@

$(THERMO_LIBNAME): $(THERMO_LIBOBJS)
	ar -r $@ $(THERMO_LIBOBJS)
	ranlib $@
	
clean:
	rm -f *.o *~

deep-clean: clean
	rm -f $(CPROPEP_LIBNAME) $(THERMO_LIBNAME)
