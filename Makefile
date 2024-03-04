CC = gcc
CC += -std=gnu11 -march=native
MYFLAGS=-Wall -pedantic -O2
MYFLAGS += -Wno-unused-result
MYLIBS=-lm

MYFFT?=-lfftw3

all: summary

summary: summary.c
	$(CC) $(MYFLAGS) -o summary summary.c $(MYFFT) $(MYLIBS)

clean:
	rm -f summary

distclean:
	rm -f *~ summary
