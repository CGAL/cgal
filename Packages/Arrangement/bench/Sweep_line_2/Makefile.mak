BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

INSTALLDIR0 = $(BINDIR)
CPPSOURCES = benchPmwxAggInsert.C

LOBJDIR =
ifeq ($(USE_CONIC_TRAITS), 1)
TARGET0 = benchPmwxAggInsertConics$(EXEFILESUFFIX)
LCPPDEFS+= -DUSE_CONIC_TRAITS
LOBJDIR :=$(LOBJDIR)_conics
else
TARGET0 = benchPmwxAggInsert$(EXEFILESUFFIX)
endif
LOBJDIR :=$(LOBJDIR)_$(COMPILER)$(COMPILER_VER)

ifeq ($(USE_LEDA_KERNEL), 1)
LCPPDEFS+= -DUSE_LEDA_KERNEL
endif

ifeq ($(USE_INSERT_FAST), 1)
LCPPDEFS+= -DUSE_INSERT_FAST
LOBJDIR :=$(LOBJDIR)_fast
TARGET0 :=$(TARGET0)_fast
endif

LCPPINCS = -I$(BASEDIR)/../Arrangement
LCPPINCS+= -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)/../../../Benchmark/include
LCPPINCS+= -I$(BASEDIR)/../../../Planar_map/include
LCPPINCS+= -I$(BASEDIR)/../../../Trapezoidal_decomposition/include
LCPPINCS+= -I$(BASEDIR)/../../../Sweep_line_2/include
LCPPINCS+= -I$(BASEDIR)/../../../Leda_rat_kernel/include
LCPPINCS+= $(CGALINCS)

include $(ROOT)/include/make/cgalrul.mak

use_insert:
	$(MAKEF) "USE_LEDA_KERNEL=1"

use_fast_insert:
	$(MAKEF) "USE_LEDA_KERNEL=1" "USE_INSERT_FAST=1"


