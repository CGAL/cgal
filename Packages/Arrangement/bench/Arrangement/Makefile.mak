BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

INSTALLDIR0 = $(BINDIR)
CPPSOURCES = benchPmwx.C

LOBJDIR =
ifeq ($(USE_LEDA_KERNEL), 1)
TARGET0 = benchPmwxLedaKernel$(EXEFILESUFFIX)
LCPPDEFS+= -DUSE_LEDA_KERNEL
LOBJDIR :=$(LOBJDIR)_leda_kernel
else
ifeq ($(USE_MY_KERNEL), 1)
TARGET0 = benchPmwxMyKernel$(EXEFILESUFFIX)
LCPPDEFS+= -DUSE_MY_KERNEL
LOBJDIR :=$(LOBJDIR)_my_kernel
else
TARGET0 = benchPmwx$(EXEFILESUFFIX)
endif
endif

LCPPINCS = -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)/../../../Benchmark/include
LCPPINCS+= -I$(BASEDIR)/../../../Planar_map/include
LCPPINCS+= -I$(BASEDIR)/../../../Trapezoidal_decomposition/include
LCPPINCS+= -I$(BASEDIR)/../../../Sweep_line_2/include
LCPPINCS+= -I$(BASEDIR)/../../../../src/LEDA_KERNEL/include
LCPPINCS+= $(CGALINCS)

include $(ROOT)/include/make/cgalrul.mak

use_leda_kernel:
	$(MAKEF) "USE_LEDA_KERNEL=1"

use_my_kernel:
	$(MAKEF) "USE_MY_KERNEL=1"
