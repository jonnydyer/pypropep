
CC   = gcc
COPT = -g -Wall -O3 #-pg -O6\
# -mpentium -ffast-math -funroll-loops -fnonnull-objects\
# -fno-exceptions -fforce-mem -fforce-addr -fcse-follow-jumps\
# -fexpensive-optimizations -march=pentium -fno-rtti #-fomit-frame-pointer

DEF = -DGCC #-DTRUE_ARRAY

all:
	make -C src all

clean:
	make -C src clean

deep-clean: clean
	make -C src deep-clean
