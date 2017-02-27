
CC   = gcc
COPT = -g -Wall -O3 #-pg 

LIB  = -lnum -lm
LIBDIR = -L../lib/
INCDIR = -I../include/

DEF = -DLINUX

PROG = test
OBJS = test.o

LIBOBJS = lu.o rk4.o rkf.o general.o print.o sec.o newton.o ptfix.o\
          sysnewton.o trapeze.o simpson.o spline.o

LIBNUM = libnum.a

all: $(LIBNUM) $(PROG)

.c.o:
	$(CC) $(DEF) $(INCDIR) $(COPT) -c $*.c -o $*.o

$(LIBNUM): $(LIBOBJS)
	ar -r $@ $(LIBOBJS) 
	ranlib $@
	mv $(LIBNUM) ../lib/

$(PROG): $(LIBNUM) $(OBJS) 
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@


clean:
	rm -f $(PROG) *.o *~

deep-clean: clean
	rm -f ../lib/$(LIBNUM)

