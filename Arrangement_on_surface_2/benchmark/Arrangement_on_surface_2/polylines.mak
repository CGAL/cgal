# Polyline:
leda_rat_cartesian_pol:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)"

leda_kernel_pol:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)"

# install - standard:

leda_rat_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

leda_rat_simple_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

leda_kernel_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

quotient_mp_float_pol_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

quotient_cgal_gmpz_pol_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

lazy_rat_pol_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

lazy_quotient_mp_float_pol_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

cgal_gmpq_exact_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

cgal_gmpq_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

cgal_gmpq_simple_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

cgal_gmpq_lazy_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

cgal_gmpq_lazy_simple_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

lazy_cgal_gmpq_cartesian_pol_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

double_pol_inst:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_TRAITS)" install

# INSTALL - non-caching:

leda_rat_cartesian_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

quotient_mp_float_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

leda_kernel_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

cgal_gmpq_exact_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

cgal_gmpq_cartesian_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

cgal_gmpq_simple_cartesian_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

cgal_gmpq_lazy_cartesian_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

cgal_gmpq_lazy_simple_cartesian_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(LAZY_SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

leda_rat_simple_cartesian_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

lazy_quotient_mp_float_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

quotient_cgal_gmpz_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_CGAL_GMPZ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

lazy_rat_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_LEDA_RAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

lazy_cgal_gmpq_cartesian_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(LAZY_CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

double_non_caching_pol_inst:
	$(MAKEF) "BENCH_NT=$(DOUBLE_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(NON_CACHING_POLYLINE_TRAITS)" install

pol_std_inst: quotient_mp_float_pol_inst \
	quotient_cgal_gmpz_pol_inst \
	lazy_quotient_mp_float_pol_inst \
	cgal_gmpq_exact_pol_inst \
	cgal_gmpq_cartesian_pol_inst \
	cgal_gmpq_simple_cartesian_pol_inst \
	cgal_gmpq_lazy_cartesian_pol_inst \
	cgal_gmpq_lazy_simple_cartesian_pol_inst \
	lazy_cgal_gmpq_cartesian_pol_inst \
	double_pol_inst

# No more LEDA support:
#	leda_kernel_pol_inst
#	leda_rat_cartesian_pol_inst
#	leda_rat_simple_cartesian_pol_inst
#	lazy_rat_pol_inst

pol_non_caching_inst: quotient_mp_float_non_caching_pol_inst \
	quotient_cgal_gmpz_non_caching_pol_inst \
	lazy_quotient_mp_float_non_caching_pol_inst \
	cgal_gmpq_exact_non_caching_pol_inst \
	cgal_gmpq_cartesian_non_caching_pol_inst \
	cgal_gmpq_simple_cartesian_non_caching_pol_inst \
	cgal_gmpq_lazy_cartesian_non_caching_pol_inst \
	cgal_gmpq_lazy_simple_cartesian_non_caching_pol_inst \
	lazy_cgal_gmpq_cartesian_non_caching_pol_inst \
	double_non_caching_pol_inst

# No more LEDA support:
#	leda_kernel_non_caching_pol_inst
#       leda_rat_cartesian_non_caching_pol_inst
#       leda_rat_simple_cartesian_non_caching_pol_inst
#	lazy_rat_non_caching_pol_inst
