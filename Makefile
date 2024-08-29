# makefile for gaffer developed on Richard's Mac

#CFLAGS= -DMACOS -O3
CFLAGS= -DMACOS -g -Wall -target arm64-apple-macos11	# for debugging

ALL=svfind ONEview

DESTDIR=~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	\rm -r *.dSYM

### object files

UTILS_OBJS=array.o utils.o
UTILS_HEADERS=utils.h array.h
$(UTILS_OBJS): $(UTILS_HEADERS)

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

### programs

alncode.o: alncode.h align.h

alnseq.o: alnseq.h

seqio.o: seqio.h

svfind: svfind.c alncode.o alnseq.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lz

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

### end of file
