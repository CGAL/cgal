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
ifeq ($(USE_MP_FLOAT), 1)
LCPPDEFS+= -DUSE_MP_FLOAT
TARGET0 := $(TARGET0)MPFloat
LOBJDIR :=$(LOBJDIR)_mp_float
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

cartesian:
	$(MAKEF)

leda_kernel:
	$(MAKEF) "USE_LEDA_KERNEL=1"

quotient_mp_float:
	$(MAKEF) "USE_MP_FLOAT=1"

lazy_rat:
	$(MAKEF) "USE_LAZY_RAT=1"

lazy_quotient_mp_float:
	$(MAKEF) "USE_LAZY_QUOTIENT=1"

# 
insert_fast:
	$(MAKEF) "USE_INSERT_TIGHT=1"

my_kernel:
	$(MAKEF) "USE_MY_KERNEL=1"

conic_traits:
	$(MAKEF) "USE_CONIC_TRAITS=1"

# Cached:

cached_traits:
	$(MAKEF) "USE_CACHED_TRAITS=1"

leda_kernel_cached_traits:
	$(MAKEF) "USE_LEDA_KERNEL=1" "USE_CACHED_TRAITS=1"

quotient_mp_float_cached_traits:
	$(MAKEF) "USE_MP_FLOAT=1" "USE_CACHED_TRAITS=1"

lazy_rat_cached_traits:
	$(MAKEF) "USE_LAZY_RAT=1" "USE_CACHED_TRAITS=1"

lazy_quotient_mp_float_cached_traits:
	$(MAKEF) "USE_LAZY_QUOTIENT=1" "USE_CACHED_TRAITS=1"

all_non_cached: cartesian \
	leda_kernel \
	quotient_mp_float \
	lazy_rat \
	lazy_quotient_mp_float

all_cached: cached_traits \
	leda_kernel_cached_traits \
	quotient_mp_float_cached_traits \
	lazy_rat_cached_traits \
	lazy_quotient_mp_float_cached_traits

# install:

cartesian_inst:
	$(MAKEF) install

leda_kernel_inst:
	$(MAKEF) "USE_LEDA_KERNEL=1" install

quotient_mp_float_inst:
	$(MAKEF) "USE_MP_FLOAT=1" install

lazy_rat_inst:
	$(MAKEF) "USE_LAZY_RAT=1" install

lazy_quotient_mp_float_inst:
	$(MAKEF) "USE_LAZY_QUOTIENT=1" install

# 
insert_fast_inst:
	$(MAKEF) "USE_INSERT_TIGHT=1" install

my_kernel_inst:
	$(MAKEF) "USE_MY_KERNEL=1" install

conic_traits_inst:
	$(MAKEF) "USE_CONIC_TRAITS=1" install

# Cached:

cached_traits_inst:
	$(MAKEF) "USE_CACHED_TRAITS=1" install

leda_kernel_cached_traits_inst:
	$(MAKEF) "USE_LEDA_KERNEL=1" "USE_CACHED_TRAITS=1" install

quotient_mp_float_cached_traits_inst:
	$(MAKEF) "USE_MP_FLOAT=1" "USE_CACHED_TRAITS=1" install

lazy_rat_cached_traits_inst:
	$(MAKEF) "USE_LAZY_RAT=1" "USE_CACHED_TRAITS=1" install

lazy_quotient_mp_float_cached_traits_inst:
	$(MAKEF) "USE_LAZY_QUOTIENT=1" "USE_CACHED_TRAITS=1" install

all_non_cached_inst: cartesian_inst \
	leda_kernel_inst \
	quotient_mp_float_inst \
	lazy_rat_inst \
	lazy_quotient_mp_float_inst

all_cached_inst: cached_traits_inst \
	leda_kernel_cached_traits_inst \
	quotient_mp_float_cached_traits_inst \
	lazy_rat_cached_traits_inst \
	lazy_quotient_mp_float_cached_traits_inst
