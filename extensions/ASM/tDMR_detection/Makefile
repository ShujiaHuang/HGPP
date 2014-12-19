#SHELL         = /bin/sh
PROG          = tdmr
DEBUG         = NO
PROFILE       = NO
PTHREADS      = YES
CC            = gcc
DEBUG_FLAGS   = -g3 -Wall -O2
PROFILE_FLAGS = -fprofile-arcs -ftest-coverage -pg 
RELEASE_FLAGS = -msse3 -O3 -funroll-loops -maccumulate-outgoing-args -fomit-frame-pointer 
STATIC_FLAGS  =
DFLAGS        = 
LIBS          = -lm -lz
OBJ           = fisher.o cytosine.o slide.o extend.o 
#TARBALL_EXCLUDE = "*.(o,gz,zip)"
#ZIP_EXCLUDE     = *.o *.gz *.zip

ifeq (YES, $(DEBUG))
        CFLAGS   = $(DEBUG_FLAGS)
        DFLAGS  += -DDEBUG
#        PTHREADS = NO
else
        CFLAGS   = $(RELEASE_FLAGS)
endif

ifeq (YES, $(PTHREADS))
        LIBS   +=  -lpthread
        DFLAGS +=  -DPTHREADS
endif

ifeq (YES, $(PROFILE))
       DFLAGS += $(PROFILE_FLAGS)
endif

.SUFFIX:
.SUFFIX: .c .o

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all: $(PROG)

$(PROG): $(OBJ) main.o
	$(CC) $(CFLAGS) $(DFLAGS) $(OBJ) main.o -o $@ $(LIBS)

fisher.o: fisher.h
cytosine.o: cytosine.h
slide.o: slide.h cytosine.h fisher.h
extend.o: extend.h cytosine.h fisher.h

clean:
	rm -f *.o 
