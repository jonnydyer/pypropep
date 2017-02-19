
CC   = gcc
COPT = -g -Wall -O3 # -pg

ROOT = ../..
INCLUDEDIR = -I$(ROOT)/libcpropep/include \
             -I$(ROOT)/libcompat/include \
             -I$(ROOT)/libnum/include \
             -I../include/

DEF = -DGCC 

LIBNAME  = libthermo.a

LIBOBJS  = load.o thermo.o

all: $(LIBNAME)

.c.o:
	$(CC) $(DEF) $(INCLUDEDIR) $(COPT) -c $*.c -o $*.o

$(LIBNAME): $(LIBOBJS)
	ar -r $@ $(LIBOBJS)
	ranlib $@
	mv $(LIBNAME) ../lib/
	
clean:
	rm -f *.o *~

deep-clean: clean
	rm -f ../lib/$(LIBNAME)
