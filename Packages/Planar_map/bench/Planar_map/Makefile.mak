BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

TARGET0 = benchPm$(EXEFILESUFFIX)
INSTALLDIR0 = $(BINDIR)
CPPSOURCES = benchPm.C

LCPPINCS = -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)/../../../Arrangement/include
LCPPINCS+= -I$(BASEDIR)/../../../Trapezoidal_decomposition/include
LCPPINCS+= -I$(BASEDIR)/../../../../src/LEDA_KERNEL/include
LCPPINCS+= $(CGALINCS)

LOBJDIR =
ifeq ($(USE_LEDA_KERNEL), 1)
LCPPDEFS+= -DUSE_LEDA_KERNEL
LOBJDIR :=$(LOBJDIR)_leda_kernel
else
ifeq ($(USE_MY_KERNEL), 1)
LCPPDEFS+= -DUSE_MY_KERNEL
LOBJDIR :=$(LOBJDIR)_my_kernel
endif
endif

include $(ROOT)/include/make/cgalrul.mak

use_leda_kernel:
	$(MAKEF) "USE_LEDA_KERNEL=1"

use_my_kernel:
	$(MAKEF) "USE_MY_KERNEL=1"
