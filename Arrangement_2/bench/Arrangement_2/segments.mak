# Segemnets:
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

cgal_gmpq_exact_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

cgal_gmpq_cartesian_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

cgal_gmpq_simple_cartesian_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

cgal_gmpq_lazy_cartesian_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

cgal_gmpq_lazy_simple_cartesian_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

lazy_cgal_gmpq_cartesian_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

double_seg:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)"

# non caching:

leda_rat_cartesian_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

leda_rat_simple_cartesian_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

non_caching_leda_kernel_seg:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

lazy_rat_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

quotient_mp_float_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

quotient_cgal_gmpz_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

lazy_quotient_mp_float_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

cgal_gmpq_exact_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

cgal_gmpq_cartesian_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

cgal_gmpq_simple_cartesian_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

cgal_gmpq_lazy_cartesian_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

cgal_gmpq_lazy_simple_cartesian_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

lazy_cgal_gmpq_cartesian_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

double_non_caching_seg:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)"

all_non_simple: leda_rat_cartesian_seg \
        leda_rat_simple_cartesian_seg \
	lazy_rat_seg \
	quotient_mp_float_seg \
	lazy_quotient_mp_float_seg \
	cgal_gmpq_exact_seg \
	cgal_gmpq_cartesian_seg \
	cgal_gmpq_simple_cartesian_seg \
	cgal_gmpq_lazy_cartesian_seg \
	cgal_gmpq_lazy_simple_cartesian_seg \
	lazy_cgal_gmpq_cartesian_seg \
	double_seg

#	leda_kernel_seg

all_simple: leda_rat_cartesian_non_caching_seg \
	quotient_mp_float_non_caching_seg \
	lazy_rat_non_caching_seg \
	lazy_quotient_mp_float_non_caching_seg \
	cgal_gmpq_exact_non_caching_seg \
	cgal_gmpq_cartesian_non_caching_seg \
	cgal_gmpq_simple_cartesian_non_caching_seg \
	cgal_gmpq_lazy_cartesian_non_caching_seg \
	cgal_gmpq_lazy_simple_cartesian_non_caching_seg \
	lazy_cgal_gmpq_cartesian_non_caching_seg \
	double_non_caching_seg

#	non_caching_leda_kernel_seg

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

cgal_gmpq_exact_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

cgal_gmpq_cartesian_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

cgal_gmpq_simple_cartesian_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

cgal_gmpq_lazy_cartesian_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

cgal_gmpq_lazy_simple_cartesian_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

lazy_cgal_gmpq_cartesian_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

double_seg_inst:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(SEGMENT_TRAITS)" install

# non caching:

leda_rat_cartesian_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

leda_rat_simple_cartesian_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

non_caching_seg_leda_kernel_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

lazy_rat_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

quotient_mp_float_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

quotient_cgal_gmpz_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

lazy_quotient_mp_float_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

cgal_gmpq_exact_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

cgal_gmpq_cartesian_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

cgal_gmpq_simple_cartesian_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

cgal_gmpq_lazy_cartesian_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

cgal_gmpq_lazy_simple_cartesian_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

lazy_cgal_gmpq_cartesian_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

double_non_caching_seg_inst:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_SEGMENT_TRAITS)" install

#
seg_caching_inst: \
	quotient_mp_float_seg_inst \
	quotient_cgal_gmpz_seg_inst \
	lazy_quotient_mp_float_seg_inst \
	cgal_gmpq_cartesian_seg_inst \
	cgal_gmpq_simple_cartesian_seg_inst \
	lazy_cgal_gmpq_cartesian_seg_inst \
	double_seg_inst \
	cgal_gmpq_lazy_cartesian_seg_inst \
	cgal_gmpq_lazy_simple_cartesian_seg_inst \
	cgal_gmpq_exact_seg_inst \

# No more LEDA support:
#       leda_rat_cartesian_seg_inst
#       leda_rat_simple_cartesian_seg_inst
#	lazy_rat_seg_inst

# The leda kernel is missing the Min_vertex_2 and Max_vertex_2 constructors
#	leda_kernel_seg_inst

seg_non_caching_inst: \
	quotient_mp_float_non_caching_seg_inst \
	quotient_cgal_gmpz_non_caching_seg_inst \
	lazy_quotient_mp_float_non_caching_seg_inst \
	cgal_gmpq_cartesian_non_caching_seg_inst \
	cgal_gmpq_simple_cartesian_non_caching_seg_inst \
	lazy_cgal_gmpq_cartesian_non_caching_seg_inst \
	double_non_caching_seg_inst \
	cgal_gmpq_lazy_cartesian_non_caching_seg_inst \
	cgal_gmpq_lazy_simple_cartesian_non_caching_seg_inst \
	cgal_gmpq_exact_non_caching_seg_inst \

# No more LEDA support:
#       leda_rat_cartesian_non_caching_seg_inst
#       leda_rat_simple_cartesian_non_caching_seg_inst
#	lazy_rat_non_caching_seg_inst

# The leda kernel is missing the Min_vertex_2 and Max_vertex_2 constructors
#	non_caching_seg_leda_kernel_inst
