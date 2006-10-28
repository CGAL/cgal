# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -Iinclude \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                bench_interval$(EXE_EXT) \
                Coercion_traits$(EXE_EXT) \
                constant$(EXE_EXT) \
                CORE_BigFloat$(EXE_EXT) \
                CORE_BigInt$(EXE_EXT) \
                CORE_BigRat$(EXE_EXT) \
                CORE_Expr$(EXE_EXT) \
                Counted_number$(EXE_EXT) \
                double$(EXE_EXT) \
                doubletst$(EXE_EXT) \
                float$(EXE_EXT) \
                floattst$(EXE_EXT) \
                Gmpq$(EXE_EXT) \
                Gmpq_new$(EXE_EXT) \
                Gmpz$(EXE_EXT) \
                int$(EXE_EXT) \
                Interval_nt$(EXE_EXT) \
                Interval_nt_new$(EXE_EXT) \
                known_bit_size_integers$(EXE_EXT) \
                Lazy_exact_nt$(EXE_EXT) \
                Lazy_exact_nt_new$(EXE_EXT) \
                leda_bigfloat$(EXE_EXT) \
                leda_integer$(EXE_EXT) \
                leda_rational$(EXE_EXT) \
                leda_real$(EXE_EXT) \
                long_double$(EXE_EXT) \
                long_doubletst$(EXE_EXT) \
                long_long$(EXE_EXT) \
                MP_Float$(EXE_EXT) \
                MP_Float_new$(EXE_EXT) \
                mpq_class$(EXE_EXT) \
                mpz_class$(EXE_EXT) \
                NT_checker$(EXE_EXT) \
                Quotient$(EXE_EXT) \
                quotient_io$(EXE_EXT) \
                Quotient_new$(EXE_EXT) \
                root_of_2$(EXE_EXT) \
                rounding_modes$(EXE_EXT) \
                simplest_rational$(EXE_EXT) \
                Sqrt_extension$(EXE_EXT) \
                _test_valid_finite_double$(EXE_EXT) \
                _test_valid_finite_float$(EXE_EXT) \
                to_interval_test$(EXE_EXT) \
                utilities$(EXE_EXT) 

bench_interval$(EXE_EXT): bench_interval$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bench_interval bench_interval$(OBJ_EXT) $(LDFLAGS)

Coercion_traits$(EXE_EXT): Coercion_traits$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Coercion_traits Coercion_traits$(OBJ_EXT) $(LDFLAGS)

constant$(EXE_EXT): constant$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constant constant$(OBJ_EXT) $(LDFLAGS)

CORE_BigFloat$(EXE_EXT): CORE_BigFloat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CORE_BigFloat CORE_BigFloat$(OBJ_EXT) $(LDFLAGS)

CORE_BigInt$(EXE_EXT): CORE_BigInt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CORE_BigInt CORE_BigInt$(OBJ_EXT) $(LDFLAGS)

CORE_BigRat$(EXE_EXT): CORE_BigRat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CORE_BigRat CORE_BigRat$(OBJ_EXT) $(LDFLAGS)

CORE_Expr$(EXE_EXT): CORE_Expr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CORE_Expr CORE_Expr$(OBJ_EXT) $(LDFLAGS)

Counted_number$(EXE_EXT): Counted_number$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Counted_number Counted_number$(OBJ_EXT) $(LDFLAGS)

double$(EXE_EXT): double$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)double double$(OBJ_EXT) $(LDFLAGS)

doubletst$(EXE_EXT): doubletst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)doubletst doubletst$(OBJ_EXT) $(LDFLAGS)

float$(EXE_EXT): float$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)float float$(OBJ_EXT) $(LDFLAGS)

floattst$(EXE_EXT): floattst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)floattst floattst$(OBJ_EXT) $(LDFLAGS)

Gmpq$(EXE_EXT): Gmpq$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Gmpq Gmpq$(OBJ_EXT) $(LDFLAGS)

Gmpq_new$(EXE_EXT): Gmpq_new$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Gmpq_new Gmpq_new$(OBJ_EXT) $(LDFLAGS)

Gmpz$(EXE_EXT): Gmpz$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Gmpz Gmpz$(OBJ_EXT) $(LDFLAGS)

int$(EXE_EXT): int$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)int int$(OBJ_EXT) $(LDFLAGS)

Interval_nt$(EXE_EXT): Interval_nt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Interval_nt Interval_nt$(OBJ_EXT) $(LDFLAGS)

Interval_nt_new$(EXE_EXT): Interval_nt_new$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Interval_nt_new Interval_nt_new$(OBJ_EXT) $(LDFLAGS)

known_bit_size_integers$(EXE_EXT): known_bit_size_integers$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)known_bit_size_integers known_bit_size_integers$(OBJ_EXT) $(LDFLAGS)

Lazy_exact_nt$(EXE_EXT): Lazy_exact_nt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_exact_nt Lazy_exact_nt$(OBJ_EXT) $(LDFLAGS)

Lazy_exact_nt_new$(EXE_EXT): Lazy_exact_nt_new$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_exact_nt_new Lazy_exact_nt_new$(OBJ_EXT) $(LDFLAGS)

leda_bigfloat$(EXE_EXT): leda_bigfloat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)leda_bigfloat leda_bigfloat$(OBJ_EXT) $(LDFLAGS)

leda_integer$(EXE_EXT): leda_integer$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)leda_integer leda_integer$(OBJ_EXT) $(LDFLAGS)

leda_rational$(EXE_EXT): leda_rational$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)leda_rational leda_rational$(OBJ_EXT) $(LDFLAGS)

leda_real$(EXE_EXT): leda_real$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)leda_real leda_real$(OBJ_EXT) $(LDFLAGS)

long_double$(EXE_EXT): long_double$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)long_double long_double$(OBJ_EXT) $(LDFLAGS)

long_doubletst$(EXE_EXT): long_doubletst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)long_doubletst long_doubletst$(OBJ_EXT) $(LDFLAGS)

long_long$(EXE_EXT): long_long$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)long_long long_long$(OBJ_EXT) $(LDFLAGS)

MP_Float$(EXE_EXT): MP_Float$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MP_Float MP_Float$(OBJ_EXT) $(LDFLAGS)

MP_Float_new$(EXE_EXT): MP_Float_new$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MP_Float_new MP_Float_new$(OBJ_EXT) $(LDFLAGS)

mpq_class$(EXE_EXT): mpq_class$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mpq_class mpq_class$(OBJ_EXT) $(LDFLAGS)

mpz_class$(EXE_EXT): mpz_class$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mpz_class mpz_class$(OBJ_EXT) $(LDFLAGS)

NT_checker$(EXE_EXT): NT_checker$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)NT_checker NT_checker$(OBJ_EXT) $(LDFLAGS)

Quotient$(EXE_EXT): Quotient$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Quotient Quotient$(OBJ_EXT) $(LDFLAGS)

quotient_io$(EXE_EXT): quotient_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quotient_io quotient_io$(OBJ_EXT) $(LDFLAGS)

Quotient_new$(EXE_EXT): Quotient_new$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Quotient_new Quotient_new$(OBJ_EXT) $(LDFLAGS)

root_of_2$(EXE_EXT): root_of_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)root_of_2 root_of_2$(OBJ_EXT) $(LDFLAGS)

rounding_modes$(EXE_EXT): rounding_modes$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rounding_modes rounding_modes$(OBJ_EXT) $(LDFLAGS)

simplest_rational$(EXE_EXT): simplest_rational$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)simplest_rational simplest_rational$(OBJ_EXT) $(LDFLAGS)

Sqrt_extension$(EXE_EXT): Sqrt_extension$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Sqrt_extension Sqrt_extension$(OBJ_EXT) $(LDFLAGS)

_test_valid_finite_double$(EXE_EXT): _test_valid_finite_double$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)_test_valid_finite_double _test_valid_finite_double$(OBJ_EXT) $(LDFLAGS)

_test_valid_finite_float$(EXE_EXT): _test_valid_finite_float$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)_test_valid_finite_float _test_valid_finite_float$(OBJ_EXT) $(LDFLAGS)

to_interval_test$(EXE_EXT): to_interval_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)to_interval_test to_interval_test$(OBJ_EXT) $(LDFLAGS)

utilities$(EXE_EXT): utilities$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)utilities utilities$(OBJ_EXT) $(LDFLAGS)

clean: \
                   bench_interval.clean \
                   Coercion_traits.clean \
                   constant.clean \
                   CORE_BigFloat.clean \
                   CORE_BigInt.clean \
                   CORE_BigRat.clean \
                   CORE_Expr.clean \
                   Counted_number.clean \
                   double.clean \
                   doubletst.clean \
                   float.clean \
                   floattst.clean \
                   Gmpq.clean \
                   Gmpq_new.clean \
                   Gmpz.clean \
                   int.clean \
                   Interval_nt.clean \
                   Interval_nt_new.clean \
                   known_bit_size_integers.clean \
                   Lazy_exact_nt.clean \
                   Lazy_exact_nt_new.clean \
                   leda_bigfloat.clean \
                   leda_integer.clean \
                   leda_rational.clean \
                   leda_real.clean \
                   long_double.clean \
                   long_doubletst.clean \
                   long_long.clean \
                   MP_Float.clean \
                   MP_Float_new.clean \
                   mpq_class.clean \
                   mpz_class.clean \
                   NT_checker.clean \
                   Quotient.clean \
                   quotient_io.clean \
                   Quotient_new.clean \
                   root_of_2.clean \
                   rounding_modes.clean \
                   simplest_rational.clean \
                   Sqrt_extension.clean \
                   _test_valid_finite_double.clean \
                   _test_valid_finite_float.clean \
                   to_interval_test.clean \
                   utilities.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

