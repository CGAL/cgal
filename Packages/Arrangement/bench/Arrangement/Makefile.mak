BASEDIR =.

include $(ROOT)/include/make/cgaldef.mak

# Initialize:
CARTESIAN_KERNEL =              0
SIMPLE_CARTESIAN_KERNEL =       1
LEDA_KERNEL =                   2
MY_KERNEL =                     3

SEGMENT_TRAITS =                0
SEGMENT_CACHED_TRAITS =         1
LEDA_SEGMENT_TRAITS =           2
POLYLINE_TRAITS =               3
CACHED_POLYLINE_TRAITS =        4
CONIC_TRAITS =                  5
EXACUS_CONIC_TRAITS =           6
CK_CONIC_TRAITS =               7
CORE_CONIC_TRAITS =             8
CK_CIRCLE_TRAITS =              9

DOUBLE_NT =                     0
MP_FLOAT_NT =                   1
GMPZ_NT =                       2
LEDA_RAT_NT =                   3
QUOTIENT_MP_FLOAT_NT =          4
QUOTIENT_CGAL_GMPZ_NT =         5
GMPQ_NT =                       6
CGAL_GMPQ_NT =                  7
LAZY_LEDA_RAT_NT =              8
LAZY_CGAL_GMPQ_NT =             9
LAZY_QUOTIENT_MP_FLOAT_NT =     10
LEDA_REAL_NT =                  11
CORE_EXPR_NT =                  12
NIX_LEDA_FIELD_WITH_SQRT_NT =   13
NIX_CORE_FIELD_WITH_SQRT_NT =   14
LAZY_GMPZ_NT =                  15


# Force values:
ifeq ($(BENCH_KERNEL), $(LEDA_KERNEL))
BENCH_NT ?= $(LEDA_RAT_NT)
endif

ifeq ($(BENCH_KERNEL), $(MY_KERNEL))
BENCH_NT ?= $(LEDA_RAT_NT)
TRAITS ?= $(LEDA_SEGMENT_TRAITS)
endif

ifeq ($(BENCH_TRAITS), $(LEDA_SEGMENT_TRAITS))
BENCH_NT ?= $(LEDA_RAT_NT)
KERNEL ?= $(LEDA_KERNEL)
endif

ifeq ($(BENCH_TRAITS), $(CONIC_TRAITS))
BENCH_NT ?= $(LEDA_REAL_NT)
endif

ifeq ($(BENCH_TRAITS), $(CORE_CONIC_TRAITS))
BENCH_NT ?= $(CORE_EXPR_NT)
endif

# default value:
BENCH_KERNEL ?= $(CARTESIAN_KERNEL)
BENCH_NT ?= $(LEDA_RAT_NT)
BENCH_TRAITS ?= $(SEGMENT_TRAITS)

# illegal combinations:
ifeq ($(BENCH_KERNEL), $(LEDA_KERNEL))
ifneq ($(BENCH_NT), $(LEDA_RAT_NT))
error "Leda kernel implies rational number type!"
endif
endif

ifeq ($(BENCH_KERNEL), $(MY_KERNEL))
ifneq ($(BENCH_NT), $(LEDA_RAT_NT))
error "My kernel implies rational number type!"
endif
ifneq ($(BENCH_TRAITS), $(LEDA_SEGMENT_TRAITS))
error "My kernel implies leda segment traits!"
endif
endif

ifeq ($(BENCH_TRAITS), $(CONIC_TRAITS))
ifneq ($(BENCH_NT), $(LEDA_REAL_NT))
error "Conic traits implies real number type!"
endif
ifeq ($(BENCH_KERNEL), $(LEDA_KERNEL))
error "Conic traits implies non leda kernel!"
endif
ifeq ($(BENCH_KERNEL), $(MY_KERNEL))
error "Conic traits implies non my kernel!"
endif
endif

ifeq ($(BENCH_TRAITS), $(CORE_CONIC_TRAITS))
ifneq ($(BENCH_NT), $(CORE_EXPR_NT))
error "Conic traits implies Core Expr number type!"
endif
endif

ifeq ($(BENCH_TRAITS), $(LEDA_SEGMENT_TRAITS))
ifneq ($(BENCH_NT), $(LEDA_RAT_NT))
error "Leda segment traits implies rational number type!"
endif
ifneq ($(BENCH_KERNEL), $(LEDA_KERNEL))
ifneq ($(BENCH_KERNEL), $(MY_KERNEL))
error "Leda segment traits implies leda kernel or my kernel!"
endif
endif
endif

BASENAME = benchPmwx
INSTALLDIR0 = $(BINDIR)
CPPSOURCES = $(BASENAME).C
TARGET0 = $(BASENAME)
LCPPDEFS+= -DCGAL_NO_PM_DEFAULT_POINT_LOCATION
# LCPPDEFS+= -DVERBOSE

LOBJDIR =

ifeq ($(DEBUG),1)
GCPPOPTS = -g
# GCPPOPTS+= -DVERBOSE
else
GCPPOPTS = -O3

ifeq ($(BENCH_TRAITS), $(CORE_CONIC_TRAITS))
GCPPOPTS = -O2 -fno-strict-aliasing
endif

ifeq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
ifeq ($(BENCH_NT), $(NIX_CORE_FIELD_WITH_SQRT_NT))
GCPPOPTS = -O2 -fno-strict-aliasing
endif
endif

ifeq ($(BENCH_TRAITS), $(CK_CIRCLE_TRAITS))
ifeq ($(BENCH_NT), $(CORE_EXPR_NT))
GCPPOPTS = -O2 -fno-strict-aliasing
endif
endif

endif

ifeq ($(USE_CGAL_WINDOW), 1)
LCPPDEFS+= -DUSE_CGAL_WINDOW
TARGET0 := $(TARGET0)CgalWindow
LOBJDIR :=$(LOBJDIR)_cgal_window
endif

# Number Type:
LCPPDEFS+= -DBENCH_NT=$(BENCH_NT)

ifeq ($(BENCH_NT), $(DOUBLE_NT))
TARGET0 := $(TARGET0)Double
LOBJDIR :=$(LOBJDIR)_double
else
ifeq ($(BENCH_NT), $(MP_FLOAT_NT))
TARGET0 := $(TARGET0)MPFloat
LOBJDIR :=$(LOBJDIR)_mp_float
else
ifeq ($(BENCH_NT), $(GMPZ_NT))
TARGET0 := $(TARGET0)Gmpz
LOBJDIR :=$(LOBJDIR)_gmpz
else
ifeq ($(BENCH_NT), $(LEDA_RAT_NT))
TARGET0 := $(TARGET0)LedaRat
LOBJDIR :=$(LOBJDIR)_leda_rat
else
ifeq ($(BENCH_NT), $(QUOTIENT_MP_FLOAT_NT))
TARGET0 := $(TARGET0)Quotient
LOBJDIR :=$(LOBJDIR)_quotient
else
ifeq ($(BENCH_NT), $(QUOTIENT_CGAL_GMPZ_NT))
TARGET0 := $(TARGET0)QuotientCgalGmpz
LOBJDIR :=$(LOBJDIR)_quotient_cgal_gmpz
else
ifeq ($(BENCH_NT), $(GMPQ_NT))
TARGET0 := $(TARGET0)Gmpq
LOBJDIR :=$(LOBJDIR)_gmpq
else
ifeq ($(BENCH_NT), $(CGAL_GMPQ_NT))
TARGET0 := $(TARGET0)CgalGmpq
LOBJDIR :=$(LOBJDIR)_cgal_gmpq
else
ifeq ($(BENCH_NT), $(LAZY_LEDA_RAT_NT))
TARGET0 := $(TARGET0)LazyRat
LOBJDIR :=$(LOBJDIR)_lazy_rat
else
ifeq ($(BENCH_NT), $(LAZY_CGAL_GMPQ_NT))
TARGET0 := $(TARGET0)LazyCgalGmpq
LOBJDIR :=$(LOBJDIR)_lazy_cgal_gmpq
else
ifeq ($(BENCH_NT), $(LAZY_GMPZ_NT))
TARGET0 := $(TARGET0)LazyGmpz
LOBJDIR :=$(LOBJDIR)_lazy_gmpz
else
ifeq ($(BENCH_NT), $(LAZY_QUOTIENT_MP_FLOAT_NT))
TARGET0 := $(TARGET0)LazyQuotient
LOBJDIR :=$(LOBJDIR)_lazy_quotient
else
ifeq ($(BENCH_NT), $(LEDA_REAL_NT))
TARGET0 := $(TARGET0)LedaReal
LOBJDIR :=$(LOBJDIR)_leda_real
else
ifeq ($(BENCH_NT), $(CORE_EXPR_NT))
TARGET0 := $(TARGET0)CoreExpr
LOBJDIR :=$(LOBJDIR)_core_expr
else
ifeq ($(BENCH_NT), $(NIX_LEDA_FIELD_WITH_SQRT_NT))
TARGET0 := $(TARGET0)NixLeda
LOBJDIR :=$(LOBJDIR)_nix_leda
else
ifeq ($(BENCH_NT), $(NIX_CORE_FIELD_WITH_SQRT_NT))
TARGET0 := $(TARGET0)NixCore
LOBJDIR :=$(LOBJDIR)_nix_core
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

# Kernel:
LCPPDEFS+= -DBENCH_KERNEL=$(BENCH_KERNEL)

ifeq ($(BENCH_KERNEL), $(CARTESIAN_KERNEL))
TARGET0 := $(TARGET0)Cartesian
LOBJDIR :=$(LOBJDIR)_cartesian
else
ifeq ($(BENCH_KERNEL), $(SIMPLE_CARTESIAN_KERNEL))
TARGET0 := $(TARGET0)SimpleCartesian
LOBJDIR :=$(LOBJDIR)_simple_cartesian
else
ifeq ($(BENCH_KERNEL), $(LEDA_KERNEL))
TARGET0 := $(TARGET0)Leda
LOBJDIR :=$(LOBJDIR)_leda
else
ifeq ($(BENCH_KERNEL), $(MY_KERNEL))
TARGET0 := $(TARGET0)My
LOBJDIR :=$(LOBJDIR)_my
endif
endif
endif
endif

# Traits:
LCPPDEFS+= -DBENCH_TRAITS=$(BENCH_TRAITS)

ifeq ($(BENCH_TRAITS), $(SEGMENT_TRAITS))
TARGET0 := $(TARGET0)Segment
LOBJDIR :=$(LOBJDIR)_segment
else
ifeq ($(BENCH_TRAITS), $(SEGMENT_CACHED_TRAITS))
TARGET0 := $(TARGET0)SegmentCached
LOBJDIR :=$(LOBJDIR)_segment_cached
else
ifeq ($(BENCH_TRAITS), $(LEDA_SEGMENT_TRAITS))
TARGET0 := $(TARGET0)LedaSegment
LOBJDIR :=$(LOBJDIR)_leda_segment
else 
ifeq ($(BENCH_TRAITS), $(POLYLINE_TRAITS))
TARGET0 := $(TARGET0)Polyline
LOBJDIR :=$(LOBJDIR)_polyline
else
ifeq ($(BENCH_TRAITS), $(POLYLINE_CACHED_TRAITS))
TARGET0 := $(TARGET0)PolylineCached
LOBJDIR :=$(LOBJDIR)_polyline_cached
else 
ifeq ($(BENCH_TRAITS), $(CONIC_TRAITS))
TARGET0 := $(TARGET0)Conic
LOBJDIR :=$(LOBJDIR)_conic
else
ifeq ($(BENCH_TRAITS), $(CORE_CONIC_TRAITS))
TARGET0 := $(TARGET0)CoreConic
LOBJDIR :=$(LOBJDIR)_core_conic
else
ifeq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
TARGET0 := $(TARGET0)ExacusConic
LOBJDIR :=$(LOBJDIR)_exacus_conic
else 
ifeq ($(BENCH_TRAITS), $(CK_CIRCLE_TRAITS))
TARGET0 := $(TARGET0)CkCircle
LOBJDIR :=$(LOBJDIR)_ck_circle
else 
ifeq ($(BENCH_TRAITS), $(CK_CONIC_TRAITS))
TARGET0 := $(TARGET0)CkConic
LOBJDIR :=$(LOBJDIR)_ck_conic
endif
endif
endif
endif
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

# ifneq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
LCPPINCS+= -I$(COREROOT)/inc
# endif

ifeq ($(BENCH_KERNEL), $(LEDA_KERNEL))
LCPPINCS+= -I$(BASEDIR)/../../../Leda_rat_kernel/include
endif

ifeq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
LCPPINCS+= -I$(EXACUS_ROOT)/NumeriX/include
LCPPINCS+= -I$(EXACUS_ROOT)/Support/include
LCPPINCS+= -I$(EXACUS_ROOT)/SweepX/include
LCPPINCS+= -I$(EXACUS_ROOT)/ConiX/include
LCPPINCS+= -I/usr/local/boost
LCPPDEFS+= -DHAVE_CONFIG_H -DQT_NO_COMPAT -DQT_CLEAN_NAMESPACE
LCPPOPTS+= -ftemplate-depth-50 -Wno-deprecated
endif

ifeq ($(BENCH_TRAITS), $(CK_CIRCLE_TRAITS))
LCPPINCS+= -I$(CURVED_KERNEL_ROOT)/include
endif
ifeq ($(BENCH_TRAITS), $(CK_CONIC_TRAITS))
LCPPINCS+= -I$(CURVED_KERNEL_ROOT)/include
endif
LCPPINCS+= $(CGALINCS)

# ifneq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
LLDOPTS = -L$(CORE_LIB_DIR)
LLDLIBS = -lcore
# else
ifeq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
LLDLIBS+= $(EXACUS_ROOT)/ConiX/src/.libs/libCnX.so
LLDLIBS+= $(EXACUS_ROOT)/SweepX/src/.libs/libSoX.so
LLDLIBS+= $(EXACUS_ROOT)/NumeriX/src/.libs/libNiX.so
LLDLIBS+= $(EXACUS_ROOT)/Support/src/.libs/libLiS.so
endif
LLDLIBS+= $(CGALLIB) $(LEDALIBS) $(CGALQTLIB) $(QTLIB)
# ifeq ($(BENCH_TRAITS), $(EXACUS_CONIC_TRAITS))
# LLDLIBS+= $(CGALCORELIB)
# endif
LLDLIBS+= $(GMPLIBS)
# LLDLIBS+= -lX11
LLDLIBS+= -lm
LLDOPTS+= $(CGALLIBOPTS) $(CGALLIBDIRS)

include $(ROOT)/include/make/cgalrul.mak

$(BASENAME).moc: $(BASENAME).C
	${QT_MOC} -o $@ $<

leda_rat_cartesian_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

leda_rat_simple_cartesian_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

leda_kernel_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

my_kernel_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(MY_KERNEL)" "BENCH_TRAITS=$(LEDA_SEGMENT_TRAITS)"

lazy_rat_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

quotient_mp_float_seg:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

quotient_cgal_qmpz:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

lazy_quotient_mp_float_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

cgal_gmpq_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

lazy_cgal_gmpq_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

double_seg:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

# Cached:

leda_rat_cartesian_cached_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

leda_rat_simple_cartesian_cached_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

cached_seg_leda_kernel:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

lazy_rat_cached_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

quotient_mp_float_cached_seg:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

quotient_cgal_gmpz_cached_seg:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

lazy_quotient_mp_float_cached_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

cgal_gmpq_cached_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

lazy_cgal_gmpq_cached_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

double_cached_seg:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)"

all_non_cached: leda_rat_cartesian_seg \
        leda_rat_simple_cartesian_seg \
	leda_kernel_seg \
	lazy_rat_seg \
	quotient_mp_float_seg \
	lazy_quotient_mp_float_seg \
	cgal_gmpq_seg \
	lazy_cgal_gmpq_seg \
	double_seg

all_cached: leda_rat_cartesian_cached_seg \
	cached_seg_leda_kernel \
	quotient_mp_float_cached_seg \
	lazy_rat_cached_seg \
	lazy_quotient_mp_float_cached_seg \
	cgal_gmpq_cached_seg \
	lazy_cgal_gmpq_cached_seg \
	double_cached_seg

# install:

leda_rat_cartesian_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

leda_rat_simple_cartesian_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

leda_kernel_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

my_kernel_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(MY_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

lazy_rat_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

quotient_mp_float_seg_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

quotient_cgal_gmpz_seg_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

lazy_quotient_mp_float_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

cgal_gmpq_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

lazy_cgal_gmpq_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

double_seg_inst:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

# Cached:

leda_rat_cartesian_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

leda_rat_simple_cartesian_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

cached_seg_leda_kernel_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

lazy_rat_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

quotient_mp_float_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

quotient_cgal_gmpz_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

lazy_quotient_mp_float_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

cgal_gmpq_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

lazy_cgal_gmpq_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

double_cached_seg_inst:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_CACHED_TRAITS)" install

#
seg_std_inst: leda_rat_cartesian_seg_inst \
        leda_rat_simple_cartesian_seg_inst \
	quotient_mp_float_seg_inst \
	quotient_cgal_gmpz_seg_inst \
	lazy_rat_seg_inst \
	lazy_quotient_mp_float_seg_inst \
	cgal_gmpq_seg_inst \
	lazy_cgal_gmpq_seg_inst \
	double_seg_inst

#	leda_kernel_seg_inst \

seg_cached_inst: leda_rat_cartesian_cached_seg_inst \
        leda_rat_simple_cartesian_cached_seg_inst \
	quotient_mp_float_cached_seg_inst \
	quotient_cgal_gmpz_cached_seg_inst \
	lazy_rat_cached_seg_inst \
	lazy_quotient_mp_float_cached_seg_inst \
	cgal_gmpq_cached_seg_inst \
	lazy_cgal_gmpq_cached_seg_inst \
	double_cached_seg_inst

#	cached_seg_leda_kernel_inst \

# Polyline:
leda_rat_cartesian_pol:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)"

leda_kernel_pol:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)"

leda_rat_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

leda_kernel_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

leda_rat_cartesian_cached_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_CACHED_TRAITS)" install

leda_kernel_cached_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(POLYLINE_CACHED_TRAITS)" install

pol_std_inst: leda_rat_cartesian_pol_inst \
        leda_rat_simple_cartesian_pol_inst \
	quotient_mp_float_pol_inst \
	quotient_cgal_gmpz_pol_inst \
	lazy_rat_pol_inst \
	lazy_quotient_mp_float_pol_inst \
	cgal_gmpq_pol_inst \
	lazy_cgal_gmpq_pol_inst \
	double_pol_inst

#	leda_kernel_pol_inst \

pol_cached_inst: leda_rat_cartesian_cached_pol_inst \
        leda_rat_simple_cartesian_cached_pol_inst \
	quotient_mp_float_cached_pol_inst \
	quotient_cgal_gmpz_cached_pol_inst \
	lazy_rat_cached_pol_inst \
	lazy_quotient_mp_float_cached_pol_inst \
	cgal_gmpq_cached_pol_inst \
	lazy_cgal_gmpq_cached_pol_inst \
	double_cached_pol_inst

#	leda_kernel_cached_pol_inst \

# Conics:
cartesian_leda_conic:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)"

simple_cartesian_leda_conic:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)"

cartesian_core_conic:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "USE_CGAL_WINDOW=1"

simple_cartesian_core_conic:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "USE_CGAL_WINDOW=1"

leda_exacus_conic:
	$(MAKEF) "BENCH_NT=$(NIX_LEDA_FIELD_WITH_SQRT_NT)" "BENCH_TRAITS=$(EXACUS_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)"

core_exacus_conic:
	$(MAKEF) "BENCH_NT=$(NIX_CORE_FIELD_WITH_SQRT_NT)" "BENCH_TRAITS=$(EXACUS_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)"

lazy_cgal_gmpq_cartesian_ck_circle:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)"

lazy_cgal_gmpq_simple_cartesian_ck_circle:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)"

gmpq_cartesian_ck_conic:
	$(MAKEF) "BENCH_NT=$(GMPQ_NT)" "BENCH_TRAITS=$(CK_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)"

gmpq_simple_cartesian_ck_conic:
	$(MAKEF) "BENCH_NT=$(GMPQ_NT)" "BENCH_TRAITS=$(CK_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)"

conics: cartesian_leda_conic \
        simple_cartesian_leda_conic \
	cartesian_core_conic \
	simple_cartesian_core_conic \
	lazy_cgal_gmpq_cartesian_ck_circle \
	lazy_cgal_gmpq_simple_cartesian_ck_circle \
	gmpq_cartesian_ck_conic \
	gmpq_simple_cartesian_ck_conic

# With install:
cartesian_leda_conic_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

simple_cartesian_leda_conic_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

cartesian_core_conic_inst:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "USE_CGAL_WINDOW=1" install

simple_cartesian_core_conic_inst:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "USE_CGAL_WINDOW=1" install

cgal_conics_int: cartesian_leda_conic_inst \
        simple_cartesian_leda_conic_inst \
	cartesian_core_conic_inst \
	simple_cartesian_core_conic_inst \

leda_exacus_conic_inst:
	$(MAKEF) "BENCH_NT=$(NIX_LEDA_FIELD_WITH_SQRT_NT)" "BENCH_TRAITS=$(EXACUS_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

core_exacus_conic_inst:
	$(MAKEF) "BENCH_NT=$(NIX_CORE_FIELD_WITH_SQRT_NT)" "BENCH_TRAITS=$(EXACUS_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

exacus_conics_inst: leda_exacus_conic_inst \
	core_exacus_conic_inst

lazy_cgal_gmpq_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

lazy_cgal_gmpq_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

lazy_gmpz_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_GMPZ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

lazy_gmpz_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_GMPZ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

cgal_gmpq_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

cgal_gmpq_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

gmpz_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(GMPZ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

gmpz_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(GMPZ_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

leda_real_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

leda_real_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

leda_rat_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

leda_rat_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

quotient_mp_float_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

quotient_mp_float_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

mp_float_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(MP_FLOAT_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

mp_float_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(MP_FLOAT_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

core_expr_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

core_expr_simple_cartesian_ck_circle_inst:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CK_CIRCLE_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

ck_circles_inst: mp_float_cartesian_ck_circle_inst \
	mp_float_simple_cartesian_ck_circle_inst \
	cgal_gmpq_cartesian_ck_circle_inst \
	cgal_gmpq_simple_cartesian_ck_circle_inst \
	gmpz_cartesian_ck_circle_inst \
	gmpz_simple_cartesian_ck_circle_inst \
	quotient_mp_float_cartesian_ck_circle_inst \
	quotient_mp_float_simple_cartesian_ck_circle_inst \
	leda_real_cartesian_ck_circle_inst \
	leda_real_simple_cartesian_ck_circle_inst \
	core_expr_cartesian_ck_circle_inst \
	core_expr_simple_cartesian_ck_circle_inst \
	lazy_gmpz_cartesian_ck_circle_inst \
	lazy_gmpz_simple_cartesian_ck_circle_inst \
	lazy_cgal_gmpq_cartesian_ck_circle_inst \
	lazy_cgal_gmpq_simple_cartesian_ck_circle_inst

#	leda_rat_cartesian_ck_circle_inst
#	leda_rat_simple_cartesian_ck_circle_inst

gmpq_cartesian_ck_conic_inst:
	$(MAKEF) "BENCH_NT=$(GMPQ_NT)" "BENCH_TRAITS=$(CK_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

gmpq_simple_cartesian_ck_conic_inst:
	$(MAKEF) "BENCH_NT=$(GMPQ_NT)" "BENCH_TRAITS=$(CK_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

ck_conics_inst:	gmpq_cartesian_ck_conic_inst \
	gmpq_simple_cartesian_ck_conic_inst

conic_inst: cgal_conics_int \
	exacus_conics_inst \
	ck_circles_inst \
	ck_conics_inst

# Miscellaneous
insert_old:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" "USE_INSERT_OLD=$(1)"

insert_old_inst:
	$(MAKEF) "USE_INSERT_OLD=$(1)" "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

# Dependencies:
$(BASENAME).o: $(BASENAME).moc
