BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

INSTALLDIR0 = $(BINDIR)
CPPSOURCES = benchSweep.C

LOBJDIR =
ifeq ($(USE_CONIC_TRAITS), 1)
TARGET0 = benchSweepConics$(EXEFILESUFFIX)
LCPPDEFS+= -DUSE_CONIC_TRAITS
LOBJDIR :=$(LOBJDIR)_conics
else
ifeq ($(USE_LEDA_KERNEL), 1)
TARGET0 = benchSweepLDSegments$(EXEFILESUFFIX)
LCPPDEFS+= -DUSE_LEDA_KERNEL
LOBJDIR :=$(LOBJDIR)_leda_kernel
else
ifeq ($(USE_MY_KERNEL), 1)
TARGET0 = benchSweepMyKernel$(EXEFILESUFFIX)
LCPPDEFS+= -DUSE_MY_KERNEL
LOBJDIR :=$(LOBJDIR)_my_kernel
else
TARGET0 = benchSweep$(EXEFILESUFFIX)
endif
endif
endif
LOBJDIR :=$(LOBJDIR)_$(COMPILER)$(COMPILER_VER)

ifeq ($(USE_INSERT_FAST), 1)
LCPPDEFS+= -DUSE_INSERT_FAST
LOBJDIR :=$(LOBJDIR)_fast
TARGET0 :=$(TARGET0)_fast
endif

LCPPINCS = -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)/../../../Benchmark/include
LCPPINCS+= -I$(BASEDIR)/../../../Planar_map/include
LCPPINCS+= -I$(BASEDIR)/../../../Trapezoidal_decomposition/include
LCPPINCS+= -I$(BASEDIR)/../../../Arrangement/include
LCPPINCS+= -I$(BASEDIR)/../../../Leda_rat_kernel/include
LCPPINCS+= -I$(BASEDIR)/../../../Arrangement/bench/Arrangement
LCPPINCS+= $(CGALINCS)

include $(ROOT)/include/make/cgalrul.mak

use_leda_kernel:
	$(MAKEF) "USE_LEDA_KERNEL=1"

use_my_kernel:
	$(MAKEF) "USE_MY_KERNEL=1"

use_conic_traits:
	$(MAKEF) "USE_CONIC_TRAITS=1"

use_insert_fast:
	$(MAKEF) "USE_INSERT_FAST=1"
