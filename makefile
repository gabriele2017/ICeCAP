CC      = cc
CFLAGS  =  -g -lm -w
CFLAGS2  = -lm -L/usr/X11R6/lib -lX11
CFLAGS3 = -lz
CFLAGS4 = -g -lm

RM      = rm -f


default: all

all: ICeCap

ICeCap: ./reference_C/new.c
	$(CC) ./reference_C/new.c -g -o ./ICeCAP $(CFLAGS)

clean veryclean:
	$(RM) ICeCap

default: shiftrev

shiftrev: ./reference_C/shiftrev.c
	$(CC) $(CFLAGS4) -g -o ./reference_C/shiftrev ./reference_C/shiftrev.c

default: gflash

gflash: ./reference_C/extend.c
	$(CC)  -o ./reference_C/gflash ./reference_C/complementReverse.c ./reference_C/combineReads.c ./reference_C/utilities2.c ./reference_C/extend.c $(CFLAGS3)

#default: plot
#plot: ./reference_C/plot.c
#	$(CC) $(CFLAGS2) -o ./ICeCap_plot ./reference_C/plot.c

clean2 veryclean2:
	$(RM) shiftrev

BOWTIE2 := $(shell command -v bowtie2 2> /dev/null)

all:
ifndef BOWTIE2
    $(error "is not available please install a suitable bowtie2!")
endif

