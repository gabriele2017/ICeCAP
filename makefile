CC      = gcc
CFLAGS  = -lm -g
CFLAGS3  = -lm -g -L/usr/X11R6/lib -lX11
CFLAGS2  = -lm -I/apps/well/gsl/2.3-gcc5.4.0/include -L/apps/well/gsl/2.3-gcc5.4.0/lib   -lgsl -lgslcblas

RM      = rm -f


default: all

all: ICeCap

ICeCap: ./reference_C/HiCinC.c
	$(CC) ./reference_C/HiCinC.c -o ./ICeCAP $(CFLAGS2) 

clean veryclean:
	$(RM) ICeCap

default: shiftrev

shiftrev: ./reference_C/shiftrev.c
	$(CC) $(CFLAGS) -o ./reference_C/shiftrev ./reference_C/shiftrev.c

default: plot
plot: ./reference_C/plot2.c
	$(CC) $(CFLAGS3) -o ./ICeCAP_plot ./reference_C/plot2.c


clean2 veryclean2:
	$(RM) shiftrev

BOWTIE2 := $(shell command -v bowtie2 2> /dev/null)

all:
ifndef BOWTIE2
    $(error "is not available. Please install bowtie2")
endif



