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
ifeq ($(KERNEL), LEDA_KERNEL)
error "Conic traits implies non leda kernel!"
endif
ifeq ($(KERNEL), MY_KERNEL)
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

GCPPOPTS = -O3

ifeq ($(USE_CGAL_WINDOW), 1)
LCPPDEFS+= -DUSE_CGAL_WINDOW
TARGET0 := $(TARGET0)CgalWindow
LOBJDIR :=$(LOBJDIR)_cgal_window
else
LLDLIBS+= -lCGALQt -lqt
LLDOPTS+= -L$(QTDIR)/lib
endif

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
ifeq ($(NT), QUOTIENT_GMPZ_NT)
LCPPDEFS+= -DUSE_QUOTIENT_GMPZ_NT
TARGET0 := $(TARGET0)QuotientGmpz
LOBJDIR :=$(LOBJDIR)_quotient_gmpz
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
ifeq ($(TRAITS), POLYLINE_CACHED_TRAITS)
LCPPDEFS+= -DUSE_POLYLINE_CACHED_TRAITS
TARGET0 := $(TARGET0)PolylineCached
LOBJDIR :=$(LOBJDIR)_polyline_cached
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

LCPPINCS = -I.
LCPPINCS+= -I$(BASEDIR)
LCPPINCS+= -I$(BASEDIR)/../../include
LCPPINCS+= -I$(BASEDIR)/../../../Benchmark/include
LCPPINCS+= -I$(BASEDIR)/../../../Planar_map/include
LCPPINCS+= -I$(BASEDIR)/../../../Arrangement/include
LCPPINCS+= -I$(BASEDIR)/../../../Trapezoidal_decomposition/include
LCPPINCS+= -I$(BASEDIR)/../../../Sweep_line_2/include
LCPPINCS+= -I$(BASEDIR)/../../../Leda_rat_kernel/include
LCPPINCS+= $(CGALINCS)

include $(ROOT)/include/make/cgalrul.mak

$(BASENAME).moc: $(BASENAME).C
	${QT_MOC} -o $@ $<

seg_cartesian:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

seg_simple_cartesian:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=SIMPLE_CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

seg_leda_kernel:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=SEGMENT_TRAITS"

seg_my_kernel:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=MY_KERNEL" "TRAITS=LEDA_SEGMENT_TRAITS"

lazy_rat_seg:
	$(MAKEF) "NT=LAZY_RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

quotient_mp_float_seg:
	$(MAKEF) "NT=QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

quotient_qmpz:
	$(MAKEF) "NT=QUOTIENT_GMPZ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

lazy_quotient_mp_float_seg:
	$(MAKEF) "NT=LAZY_QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

gmpq_seg:
	$(MAKEF) "NT=GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

lazy_gmpq_seg:
	$(MAKEF) "NT=LAZY_GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

double_seg:
	$(MAKEF) "NT=DOUBLE_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS"

# 
insert_old:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" "USE_INSERT_OLD=1"

conic:
	$(MAKEF) "NT=REAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=CONIC_TRAITS"

pol_cartesian:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=POLYLINE_TRAITS"

pol_leda_kernel:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=POLYLINE_TRAITS"

# Cached:

cached_seg_cartesian:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

cached_seg_simple_cartesian:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=SIMPLE_CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

cached_seg_leda_kernel:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

lazy_rat_cached_seg:
	$(MAKEF) "NT=LAZY_RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

quotient_mp_float_cached_seg:
	$(MAKEF) "NT=QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

quotient_gmpz_cached_seg:
	$(MAKEF) "NT=QUOTIENT_GMPZ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

lazy_quotient_mp_float_cached_seg:
	$(MAKEF) "NT=LAZY_QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

gmpq_cached_seg:
	$(MAKEF) "NT=GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

lazy_gmpq_cached_seg:
	$(MAKEF) "NT=LAZY_GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

double_cached_seg:
	$(MAKEF) "NT=DOUBLE_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS"

all_non_cached: seg_cartesian \
        seg_simple_cartesian \
	seg_leda_kernel \
	lazy_rat_seg \
	quotient_mp_float_seg \
	lazy_quotient_mp_float_seg \
	gmpq_seg \
	lazy_gmpq_seg \
	double_seg

all_cached: cached_seg_cartesian \
	cached_seg_leda_kernel \
	quotient_mp_float_cached_seg \
	lazy_rat_cached_seg \
	lazy_quotient_mp_float_cached_seg \
	gmpq_cached_seg \
	lazy_gmpq_cached_seg \
	double_cached_seg

# install:

seg_cartesian_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

seg_simple_cartesian_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=SIMPLE_CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

seg_leda_kernel_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=SEGMENT_TRAITS" install

seg_my_kernel_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=MY_KERNEL" "TRAITS=SEGMENT_TRAITS" install

lazy_rat_seg_inst:
	$(MAKEF) "NT=LAZY_RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

quotient_mp_float_seg_inst:
	$(MAKEF) "NT=QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

quotient_gmpz_seg_inst:
	$(MAKEF) "NT=QUOTIENT_GMPZ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

lazy_quotient_mp_float_seg_inst:
	$(MAKEF) "NT=LAZY_QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

gmpq_seg_inst:
	$(MAKEF) "NT=GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

lazy_gmpq_seg_inst:
	$(MAKEF) "NT=LAZY_GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

double_seg_inst:
	$(MAKEF) "NT=DOUBLE_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

# 
insert_old_inst:
	$(MAKEF) "USE_INSERT_OLD=1" "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_TRAITS" install

conic_inst:
	$(MAKEF) "NT=REAL_NT" "TRAITS=CONIC_TRAITS" "KERNEL=CARTESIAN_KERNEL" install

pol_cartesian_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=POLYLINE_TRAITS" install

pol_leda_kernel_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=POLYLINE_TRAITS" install

# Cached:

cached_seg_cartesian_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

cached_seg_simple_cartesian_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=SIMPLE_CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

cached_seg_leda_kernel_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

lazy_rat_cached_seg_inst:
	$(MAKEF) "NT=LAZY_RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

quotient_mp_float_cached_seg_inst:
	$(MAKEF) "NT=QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

quotient_gmpz_cached_seg_inst:
	$(MAKEF) "NT=QUOTIENT_GMPZ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

lazy_quotient_mp_float_cached_seg_inst:
	$(MAKEF) "NT=LAZY_QUOTIENT_MP_FLOAT_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

gmpq_cached_seg_inst:
	$(MAKEF) "NT=GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

lazy_gmpq_cached_seg_inst:
	$(MAKEF) "NT=LAZY_GMPQ_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

double_cached_seg_inst:
	$(MAKEF) "NT=DOUBLE_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=SEGMENT_CACHED_TRAITS" install

# Polyline:
cached_pol_cartesian_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=CARTESIAN_KERNEL" "TRAITS=POLYLINE_CACHED_TRAITS" install

cached_pol_leda_kernel_inst:
	$(MAKEF) "NT=RATIONAL_NT" "KERNEL=LEDA_KERNEL" "TRAITS=POLYLINE_CACHED_TRAITS" install

#
seg_std_inst: seg_cartesian_inst \
        seg_simple_cartesian_inst \
	quotient_mp_float_seg_inst \
	quotient_gmpz_seg_inst \
	lazy_rat_seg_inst \
	lazy_quotient_mp_float_seg_inst \
	gmpq_seg_inst \
	lazy_gmpq_seg_inst \
	double_seg_inst

#	seg_leda_kernel_inst \

seg_cached_inst: cached_seg_cartesian_inst \
        cached_seg_simple_cartesian_inst \
	quotient_mp_float_cached_seg_inst \
	quotient_gmpz_cached_seg_inst \
	lazy_rat_cached_seg_inst \
	lazy_quotient_mp_float_cached_seg_inst \
	gmpq_cached_seg_inst \
	lazy_gmpq_cached_seg_inst \
	double_cached_seg_inst

#	cached_seg_leda_kernel_inst \

pol_std_inst: pol_cartesian_inst \
        pol_simple_cartesian_inst \
	quotient_mp_float_pol_inst \
	quotient_gmpz_pol_inst \
	lazy_rat_pol_inst \
	lazy_quotient_mp_float_pol_inst \
	gmpq_pol_inst \
	lazy_gmpq_pol_inst \
	double_pol_inst

#	pol_leda_kernel_inst \

pol_cached_inst: cached_pol_cartesian_inst \
        cached_pol_simple_cartesian_inst \
	quotient_mp_float_cached_pol_inst \
	quotient_gmpz_cached_pol_inst \
	lazy_rat_cached_pol_inst \
	lazy_quotient_mp_float_cached_pol_inst \
	gmpq_cached_pol_inst \
	lazy_gmpq_cached_pol_inst \
	double_cached_pol_inst

#	cached_pol_leda_kernel_inst \

$(BASENAME).o: $(BASENAME).moc
