BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

BASENAME = benchPmwx
INSTALLDIR0 = $(BINDIR)
CPPSOURCES = $(BASENAME).C
TARGET0 = $(BASENAME)
LCPPDEFS+= -DCGAL_NO_PM_DEFAULT_POINT_LOCATION

# forced values:
ifeq ($(KERNEL), LEDA_KERNEL)
NT ?= RATIONAL_NT
endif

ifeq ($(KERNEL), MY_KERNEL)
NT ?= RATIONAL_NT
TRAITS ?= LEDA_SEGMENT_TRAITS
endif

ifeq ($(TRAITS), LEDA_SEGMENT_TRAITS)
NT ?= RATIONAL_NT
KERNEL ?= LEDA_KERNEL
endif

ifeq ($(TRAITS), CONIC_TRAITS)
NT ?= REAL_NT
endif

# default value:
KERNEL ?= CARTESIAN_KERNEL
NT ?= RATIONAL_NT
TRAITS ?= SEGMENT_TRAITS

# illegal combinations:
ifeq ($(KERNEL), LEDA_KERNEL)
ifneq ($(NT), RATIONAL_NT)
error "Leda kernel implies rational number type!"
endif
endif

ifeq ($(KERNEL), MY_KERNEL)
ifneq ($(NT), RATIONAL_NT)
error "My kernel implies rational number type!"
endif
ifneq ($(TRAITS), LEDA_SEGMENT_TRAITS)
error "My kernel implies leda segment traits!"
endif
endif

ifeq ($(TRAITS), CONIC_TRAITS)
ifneq ($(NT), REAL_NT)
error "Conic traits implies real number type!"
endif
ifneq ($(KERNEL), LEDA_KERNEL)
error "Conic traits implies non leda kernel!"
endif
ifneq ($(KERNEL), MY_KERNEL)
error "Conic traits implies non my kernel!"
endif
endif

ifeq ($(TRAITS), LEDA_SEGMENT_TRAITS)
ifneq ($(NT), RATIONAL_NT)
error "Leda segment traits implies rational number type!"
endif
ifneq ($(KERNEL), LEDA_KERNEL)
ifneq ($(KERNEL), MY_KERNEL)
error "Leda segment traits implies leda kernel or my kernel!"
endif
endif
endif

LOBJDIR =

ifeq ($(NT), RATIONAL_NT)
LCPPDEFS+= -DUSE_RATIONAL_NT
TARGET0 := $(TARGET0)Rat
LOBJDIR :=$(LOBJDIR)_rat
else
ifeq ($(NT), QUOTIENT_MP_FLOAT_NT)
LCPPDEFS+= -DUSE_QUOTIENT_MP_FLOAT_NT
TARGET0 := $(TARGET0)Quotient
LOBJDIR :=$(LOBJDIR)_quotient
else
ifeq ($(NT), GMPQ_NT)
LCPPDEFS+= -DUSE_GMPQ_NT
TARGET0 := $(TARGET0)Gmpq
LOBJDIR :=$(LOBJDIR)_gmpq
else
ifeq ($(NT), DOUBLE_NT)
LCPPDEFS+= -DUSE_DOUBLE_NT
TARGET0 := $(TARGET0)Double
LOBJDIR :=$(LOBJDIR)_double
else
ifeq ($(NT), REAL_NT)
LCPPDEFS+= -DUSE_REAL_NT
TARGET0 := $(TARGET0)Real
LOBJDIR :=$(LOBJDIR)_real
else
ifeq ($(NT), LAZY_RATIONAL_NT)
LCPPDEFS+= -DUSE_LAZY_RATIONAL_NT
TARGET0 := $(TARGET0)LazyRat
LOBJDIR :=$(LOBJDIR)_lazy_rat
else
ifeq ($(NT), LAZY_QUOTIENT_MP_FLOAT_NT)
LCPPDEFS+= -DUSE_LAZY_QUOTIENT_MP_FLOAT_NT
TARGET0 := $(TARGET0)LazyQuotient
LOBJDIR :=$(LOBJDIR)_lazy_quotient
else
ifeq ($(NT), LAZY_GMPQ_NT)
LCPPDEFS+= -DUSE_LAZY_GMPQ_NT
TARGET0 := $(TARGET0)LazyGmpq
LOBJDIR :=$(LOBJDIR)_lazy_gmpq
endif
endif
endif
endif
endif
endif
endif
endif

ifeq ($(KERNEL), CARTESIAN_KERNEL)
LCPPDEFS+= -DUSE_CARTESIAN_KERNEL
TARGET0 := $(TARGET0)Cartesian
LOBJDIR :=$(LOBJDIR)_cartesian
else
ifeq ($(KERNEL), SIMPLE_CARTESIAN_KERNEL)
LCPPDEFS+= -DUSE_SIMPLE_CARTESIAN_KERNEL
TARGET0 := $(TARGET0)SimpleCartesian
LOBJDIR :=$(LOBJDIR)_simple_cartesian
else
ifeq ($(KERNEL), LEDA_KERNEL)
LCPPDEFS+= -DUSE_LEDA_KERNEL
TARGET0 := $(TARGET0)Leda
LOBJDIR :=$(LOBJDIR)_leda
else
ifeq ($(KERNEL), MY_KERNEL)
LCPPDEFS+= -DUSE_MY_KERNEL
TARGET0 := $(TARGET0)My
LOBJDIR :=$(LOBJDIR)_my
endif
endif
endif
endif

ifeq ($(TRAITS), SEGMENT_TRAITS)
LCPPDEFS+= -DUSE_SEGMENT_TRAITS
TARGET0 := $(TARGET0)Segment
LOBJDIR :=$(LOBJDIR)_segment
else
ifeq ($(TRAITS), SEGMENT_CACHED_TRAITS)
LCPPDEFS+= -DUSE_SEGMENT_CACHED_TRAITS
TARGET0 := $(TARGET0)SegmentCached
LOBJDIR :=$(LOBJDIR)_segment_cached
else
ifeq ($(TRAITS), LEDA_SEGMENT_TRAITS)
LCPPDEFS+= -DUSE_LEDA_SEGMENT_TRAITS
TARGET0 := $(TARGET0)LedaSegment
LOBJDIR :=$(LOBJDIR)_leda_segment
else 
ifeq ($(TRAITS), CONIC_TRAITS)
LCPPDEFS+= -DUSE_CONIC_TRAITS
TARGET0 := $(TARGET0)Conic
LOBJDIR :=$(LOBJDIR)_conic
else
ifeq ($(TRAITS), POLYLINE_TRAITS)
LCPPDEFS+= -DUSE_POLYLINE_TRAITS
TARGET0 := $(TARGET0)Polyline
LOBJDIR :=$(LOBJDIR)_polyline
else
ifeq ($(TRAITS), LEDA_POLYLINE_TRAITS)
LCPPDEFS+= -DUSE_LEDA_POLYLINE_TRAITS
TARGET0 := $(TARGET0)LedaPolyline
LOBJDIR :=$(LOBJDIR)_leda_polyline
endif
endif
endif
endif
endif
endif

ifeq ($(USE_INSERT_OLD), 1)
LCPPDEFS+= -DUSE_INSERT_OLD
TARGET0 :=$(TARGET0)Old
LOBJDIR :=$(LOBJDIR)_old
endif

TARGET0 := $(TARGET0)$(EXEFILESUFFIX)
LOBJDIR := $(LOBJDIR)_$(COMPILER)$(COMPILER_VER)

LCPPINCS = -I$(BASEDIR)
LCPPINCS+= -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)/../../../Benchmark/include
LCPPINCS+= -I$(BASEDIR)/../../../Planar_map/include
LCPPINCS+= -I$(BASEDIR)/../../../Trapezoidal_decomposition/include
LCPPINCS+= -I$(BASEDIR)/../../../Sweep_line_2/include
LCPPINCS+= -I$(BASEDIR)/../../../Leda_rat_kernel/include
LCPPINCS+= -I$(BASEDIR)/../../../Cartesian_kernel/include
LCPPINCS+= $(CGALINCS)

include $(ROOT)/include/make/cgalrul.mak

cartesian:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

simple_cartesian:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_SIMPLE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

leda_kernel:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_LEDA_KERNEL=1" "USE_SEGMENT_TRAITS=1"

my_kernel:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_MY_KERNEL=1" "USE_LEDA_SEGMENT_TRAITS=1"

lazy_rat:
	$(MAKEF) "USE_LAZY_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

quotient_mp_float:
	$(MAKEF) "USE_QUOTIENT_MP_FLOAT_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

lazy_quotient_mp_float:
	$(MAKEF) "USE_LAZY_QUOTIENT_MP_FLOAT_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

gmpq:
	$(MAKEF) "USE_GMPQ_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

lazy_gmpq:
	$(MAKEF) "USE_LAZY_GMPQ_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

double:
	$(MAKEF) "USE_DOUBLE_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1"

# 
insert_old:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_TRAITS=1" "USE_INSERT_OLD=1"

conic_traits:
	$(MAKEF) "USE_REAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_CONIC_TRAITS=1"

polyline_traits:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_POLYLINE_TRAITS=1"

leda_polyline_traits:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_LEDA_KERNEL=1" "USE_LEDA_POLYLINE_TRAITS=1"

# Cached:

cached_traits:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

simple_cartesian_cached_traits:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_SIMPLE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

leda_kernel_cached_traits:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_LEDA_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

lazy_rat_cached_traits:
	$(MAKEF) "USE_LAZY_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

quotient_mp_float_cached_traits:
	$(MAKEF) "USE_QUOTIENT_MP_FLOAT_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

lazy_quotient_mp_float_cached_traits:
	$(MAKEF) "USE_LAZY_QUOTIENT_FLOAT_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

gmpq_cached_traits:
	$(MAKEF) "USE_GMPQ_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

lazy_gmpq_cached_traits:
	$(MAKEF) "USE_LAZY_GMPQ_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

double_cached_traits:
	$(MAKEF) "USE_DOUBLE_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1"

all_non_cached: cartesian \
        simple_cartesian \
	leda_kernel \
	lazy_rat \
	quotient_mp_float \
	lazy_quotient_mp_float \
	gmpq \
	lazy_gmpq \
	double

all_cached: cached_traits \
	leda_kernel_cached_traits \
	quotient_mp_float_cached_traits \
	lazy_rat_cached_traits \
	lazy_quotient_mp_float_cached_traits \
	gmpq_cached_traits \
	lazy_gmpq_cached_traits \
	double_cached_traits

# install:

cartesian_inst:
	$(MAKEF) install

simple_cartesian_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_SIMPLE_CARTESIAN_KERNEL=1" install

leda_kernel_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_LEDA_KERNEL=1" install

my_kernel_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_MY_KERNEL=1" install

lazy_rat_inst:
	$(MAKEF) "USE_LAZY_RATIONAL_NT=1" install

quotient_mp_float_inst:
	$(MAKEF) "USE_QUOTIENT_MP_FLOAT_NT=1" install

lazy_quotient_mp_float_inst:
	$(MAKEF) "USE_LAZY_QUOTIENT_MP_FLOAT_NT=1" install

gmpq_inst:
	$(MAKEF) "USE_GMPQ_NT=1" install

lazy_gmpq_inst:
	$(MAKEF) "USE_LAZY_GMPQ_NT=1" install

double_inst:
	$(MAKEF) "USE_DOUBLE_NT=1" install

# 
insert_old_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_INSERT_OLD=1" install

conic_traits_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_CONIC_TRAITS=1" install

polyline_traits_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=POLYLINE_TRAITS" install

leda_polyline_traits_inst:
	$(MAKEF) "NT=USE_RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=LEDA_POLYLINE_TRAITS" install

# Cached:

cached_traits_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

simple_cartesian_cached_traits_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_SIMPLE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

leda_kernel_cached_traits_inst:
	$(MAKEF) "USE_RATIONAL_NT=1" "USE_LEDA_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

lazy_rat_cached_traits_inst:
	$(MAKEF) "USE_LAZY_RATIONAL_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

quotient_mp_float_cached_traits_inst:
	$(MAKEF) "USE_QUOTIENT_MP_FLOAT_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

lazy_quotient_mp_float_cached_traits_inst:
	$(MAKEF) "USE_LAZY_QUOTIENT_MP_FLOAT_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

gmpq_cached_traits_inst:
	$(MAKEF) "USE_GMPQ_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

lazy_gmpq_cached_traits_inst:
	$(MAKEF) "USE_LAZY_GMPQ_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

double_cached_traits_inst:
	$(MAKEF) "USE_DOUBLE_NT=1" "USE_CARTESIAN_KERNEL=1" "USE_SEGMENT_CACHED_TRAITS=1" install

#
all_non_cached_inst: cartesian_inst \
        simple_cartesian_inst \
	leda_kernel_inst \
	quotient_mp_float_inst \
	lazy_rat_inst \
	lazy_quotient_mp_float_inst \
	gmpq_inst \
	lazy_gmpq_inst \
	double_inst

all_cached_inst: cached_traits_inst \
        simple_cartesian_cached_traits_inst \
	leda_kernel_cached_traits_inst \
	quotient_mp_float_cached_traits_inst \
	lazy_rat_cached_traits_inst \
	lazy_quotient_mp_float_cached_traits_inst \
	gmpq_cached_traits_inst \
	lazy_gmpq_cached_traits_inst \
	double_cached_traits_inst
