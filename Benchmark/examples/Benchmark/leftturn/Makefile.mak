BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

BASENAME = smidgen
INSTALLDIR0 = $(BINDIR)
CPPSOURCES = $(BASENAME).C
TARGET0 = $(BASENAME)

LCPPINCS= -I$(BASEDIR)/../../../All/Benchmark/include
LCPPINCS+= $(CGALINCS)

include $(ROOT)/include/make/cgalrul.mak
