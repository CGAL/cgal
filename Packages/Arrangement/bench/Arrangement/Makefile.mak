BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

BASENAME = benchPmwx
INSTALLDIR0 = $(BINDIR)
CPPSOURCES = $(BASENAME).C
TARGET0 = $(BASENAME)

LOBJDIR =
ifeq ($(USE_CONIC_TRAITS), 1)
LCPPDEFS+= -DUSE_CONIC_TRAITS
TARGET0 := $(TARGET0)Conics
LOBJDIR :=$(LOBJDIR)_conics
else

ifeq ($(USE_LEDA_KERNEL), 1)
LCPPDEFS+= -DUSE_LEDA_KERNEL
TARGET0 := $(TARGET0)LedaKernel
LOBJDIR :=$(LOBJDIR)_leda_kernel
else
ifeq ($(USE_MY_KERNEL), 1)
LCPPDEFS+= -DUSE_MY_KERNEL
TARGET0 := $(TARGET0)MyKernel
LOBJDIR :=$(LOBJDIR)_my_kernel
else
ifeq ($(USE_LAZY_RAT), 1)
LCPPDEFS+= -DUSE_LAZY_RAT
TARGET0 := $(TARGET0)LazyRat
LOBJDIR :=$(LOBJDIR)_lazy_rat
else
ifeq ($(USE_LAZY_QUOTIENT), 1)
LCPPDEFS+= -DUSE_LAZY_QUOTIENT
TARGET0 := $(TARGET0)LazyQuotient
LOBJDIR :=$(LOBJDIR)_lazy_quotient
endif
endif
endif
endif
endif

ifeq ($(USE_CACHED_TRAITS), 1)
LCPPDEFS+= -DUSE_CACHED_TRAITS
TARGET0 := $(TARGET0)Cached
LOBJDIR :=$(LOBJDIR)_cached
endif

ifeq ($(USE_TIGHT_TRAITS), 1)
LCPPDEFS+= -DUSE_TIGHT_TRAITS
TARGET0 :=$(TARGET0)Tight
LOBJDIR :=$(LOBJDIR)_tight
endif

ifeq ($(USE_LEDA_SEGMENT_TRAITS), 1)
LCPPDEFS+= -DUSE_LEDA_SEGMENT_TRAITS
TARGET0 := $(TARGET0)LedaSegment
LOBJDIR :=$(LOBJDIR)_leda_segment
endif

TARGET0 := $(TARGET0)$(EXEFILESUFFIX)
LOBJDIR := $(LOBJDIR)_$(COMPILER)$(COMPILER_VER)

LCPPINCS = -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)/../../../Benchmark/include
LCPPINCS+= -I$(BASEDIR)/../../../Planar_map/include
LCPPINCS+= -I$(BASEDIR)/../../../Trapezoidal_decomposition/include
LCPPINCS+= -I$(BASEDIR)/../../../Sweep_line_2/include
LCPPINCS+= -I$(BASEDIR)/../../../Leda_rat_kernel/include
LCPPINCS+= $(CGALINCS)

include $(ROOT)/include/make/cgalrul.mak

use_leda_kernel:
	$(MAKEF) "USE_LEDA_KERNEL=1"

use_my_kernel:
	$(MAKEF) "USE_MY_KERNEL=1"

use_conic_traits:
	$(MAKEF) "USE_CONIC_TRAITS=1"

use_insert_fast:
	$(MAKEF) "USE_INSERT_TIGHT=1"
