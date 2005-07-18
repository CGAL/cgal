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

quotient_mp_float_cached_pol_inst:
	$(MAKEF) "BENCH_NT=$(QUOTIENT_MP_FLOAT_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_CACHED_TRAITS)" install

leda_kernel_cached_pol_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_RAT_NT)" "BENCH_KERNEL=$(LEDA_KERNEL)" "BENCH_TRAITS=$(POLYLINE_CACHED_TRAITS)" install

cgal_gmpq_cached_pol_inst:
	$(MAKEF) "BENCH_NT=$(CGAL_GMPQ_NT)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" "BENCH_TRAITS=$(POLYLINE_CACHED_TRAITS)" install

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

