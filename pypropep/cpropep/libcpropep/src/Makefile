
CC   = gcc
COPT = -g -Wall -O3 

ROOT = ../..

INCLUDEDIR = -I$(ROOT)/libcompat/include \
             -I$(ROOT)/libthermo/include \
             -I$(ROOT)/libnum/include \
             -I../include/

DEF = -DGCC #-DTRUE_ARRAY

LIBNAME = libcpropep.a

LIBOBJS = equilibrium.o print.o performance.o derivative.o

all: $(LIBNAME)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(LIBNAME): $(LIBOBJS)
	ar -r $@ $(LIBOBJS)
	ranlib $@
	mv $(LIBNAME) ../lib

clean:
	rm -f *.o *~

deep-clean: clean
	rm -f ../lib/$(LIBNAME) 
