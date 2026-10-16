# Conics:
cartesian_leda_conic:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(LEDA_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)"

simple_cartesian_leda_conic:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(LEDA_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)"

cartesian_core_conic:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)"

simple_cartesian_core_conic:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)"

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
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(LEDA_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

simple_cartesian_leda_conic_inst:
	$(MAKEF) "BENCH_NT=$(LEDA_REAL_NT)" "BENCH_TRAITS=$(LEDA_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

cartesian_core_conic_inst:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(CARTESIAN_KERNEL)" install

simple_cartesian_core_conic_inst:
	$(MAKEF) "BENCH_NT=$(CORE_EXPR_NT)" "BENCH_TRAITS=$(CORE_CONIC_TRAITS)" "BENCH_KERNEL=$(SIMPLE_CARTESIAN_KERNEL)" install

cgal_conics_inst: cartesian_core_conic_inst \
	simple_cartesian_core_conic_inst \

#	cartesian_leda_conic_inst
#        simple_cartesian_leda_conic_inst

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

conic_inst: cgal_conics_inst \
	exacus_conics_inst \
	ck_circles_inst \
	ck_conics_inst

