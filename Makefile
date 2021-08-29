CFLAGS=-I. -O3 -g
Q_AR = ar

CC=`which cc`

#JLSE
ifeq ($(JLSE),yes)
#CFLAGS+=-D_COBALT -fPIC
CFLAGS+=-fPIC
ifneq ($(NOMPI),yes)
CC=`which mpicc`
endif
ifneq ($(NOOMP),yes)
CFLAGS+=-fopenmp
endif
endif

#CRAY
ifeq ($(CRAY), yes)
CFLAGS+=-D_PMI -D_CRAY -fPIC
ifneq ($(NOOMP),yes)
CFLAGS+=-qopenmp
endif
endif

#BGQ
ifeq ($(BGQ),yes)
CC=mpixlc
ifneq ($(NOOMP),yes)
CFLAGS+= -qsmp=omp -qthreaded
endif
CFLAGS+=-D_BGQ -qpic
Q_AR = /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc64-bgq-linux-ar
endif

#other flags
ifeq ($(DEBUG),yes)
CFLAGS+=-D_DEBUG
endif
ifeq ($(TRACE),yes)
CFLAGS+=-D_TRACE
endif
ifeq ($(TIMER_OFF),yes)
CFLAGS+=-D_TIMER_OFF
endif
ifeq ($(BENCH),yes)
CFLAGS+=-D_BENCH
endif
ifeq ($(NOMPI),yes)
CFLAGS+=-D_NOMPI
endif
ifeq ($(NOOMP),yes)
CFLAGS+=-D_NOOMP
endif
ifeq ($(COBALT),yes)
CFLAGS+=-D_COBALT
endif

CFLAGS+=-I./include

LIB=

LIBDIR=lib
OBJDIR=bin

all: $(LIBDIR)/libpolimer.a $(LIBDIR)/libpolimer.so

OBJ = $(OBJDIR)/PoLiMEr.o $(OBJDIR)/PoLiLog.o $(OBJDIR)/msr-handler.o

ifeq ($(CRAY),yes)
OBJ+= $(OBJDIR)/cray_pm-handler.o
endif

ifeq ($(BGQ),yes)
OBJ+= $(OBJDIR)/bgq-handler.o
endif

$(OBJ): $(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LIB)

$(LIBDIR)/libpolimer.a: $(OBJ)
	$(Q_AR) rcs $@ $(OBJ)

$(LIBDIR)/libpolimer.so: $(OBJ)
	`which cc` -shared -o $@ $(OBJ)

clean:
	rm -f lib/*.a bin/*.o lib/*.so a.out
