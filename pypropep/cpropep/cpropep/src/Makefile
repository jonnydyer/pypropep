CC     = gcc
COPT   = -g -Wall -O3 #-pg 

LIB    = -lcpropep -lthermo -lnum -lm  
ROOT   = ../..
LIBDIR = -L$(ROOT)/libnum/lib \
         -L$(ROOT)/libthermo/lib \
         -L$(ROOT)/libcpropep/lib

INCDIR = -I$(ROOT)/libnum/include/ \
         -I$(ROOT)/libthermo/include/ \
         -I$(ROOT)/libcpropep/include/ \
         -I$(ROOT)/libcompat/include/

DEF    = -DGCC -DCONF_FILE=\"/etc/rocketworkbench/cpropep.conf\"
PROG   = cpropep
OBJS   = cpropep.o

all: $(PROG)

.c.o:
	$(CC) $(DEF) $(INCDIR) $(COPT) -c $*.c -o $*.o

$(PROG): $(OBJS)
	$(CC) $(COPT) $(OBJS) $(LIBDIR) $(LIB) -o $@

clean:
	rm -f *.o *~

deep-clean: clean
	rm -f $(PROG)
