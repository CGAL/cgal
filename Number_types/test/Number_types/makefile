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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                bench_interval$(EXE_EXT) \
                constant$(EXE_EXT) \
                doubletst$(EXE_EXT) \
                floattst$(EXE_EXT) \
                Gmpq$(EXE_EXT) \
                Interval_nt$(EXE_EXT) \
                known_bit_size_integers$(EXE_EXT) \
                Lazy_exact_nt$(EXE_EXT) \
                long_doubletst$(EXE_EXT) \
                MP_Float$(EXE_EXT) \
                NT_checker$(EXE_EXT) \
                Quotient$(EXE_EXT) \
                quotient_io$(EXE_EXT) \
                root_of_2$(EXE_EXT) \
                rounding_modes$(EXE_EXT) \
                simplest_rational$(EXE_EXT) \
                _test_valid_finite_double$(EXE_EXT) \
                _test_valid_finite_float$(EXE_EXT) \
                to_interval_test$(EXE_EXT) \
                utilities$(EXE_EXT) 

bench_interval$(EXE_EXT): bench_interval$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bench_interval bench_interval$(OBJ_EXT) $(LDFLAGS)

constant$(EXE_EXT): constant$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constant constant$(OBJ_EXT) $(LDFLAGS)

doubletst$(EXE_EXT): doubletst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)doubletst doubletst$(OBJ_EXT) $(LDFLAGS)

floattst$(EXE_EXT): floattst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)floattst floattst$(OBJ_EXT) $(LDFLAGS)

Gmpq$(EXE_EXT): Gmpq$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Gmpq Gmpq$(OBJ_EXT) $(LDFLAGS)

Interval_nt$(EXE_EXT): Interval_nt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Interval_nt Interval_nt$(OBJ_EXT) $(LDFLAGS)

known_bit_size_integers$(EXE_EXT): known_bit_size_integers$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)known_bit_size_integers known_bit_size_integers$(OBJ_EXT) $(LDFLAGS)

Lazy_exact_nt$(EXE_EXT): Lazy_exact_nt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_exact_nt Lazy_exact_nt$(OBJ_EXT) $(LDFLAGS)

long_doubletst$(EXE_EXT): long_doubletst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)long_doubletst long_doubletst$(OBJ_EXT) $(LDFLAGS)

MP_Float$(EXE_EXT): MP_Float$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MP_Float MP_Float$(OBJ_EXT) $(LDFLAGS)

NT_checker$(EXE_EXT): NT_checker$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)NT_checker NT_checker$(OBJ_EXT) $(LDFLAGS)

Quotient$(EXE_EXT): Quotient$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Quotient Quotient$(OBJ_EXT) $(LDFLAGS)

quotient_io$(EXE_EXT): quotient_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quotient_io quotient_io$(OBJ_EXT) $(LDFLAGS)

root_of_2$(EXE_EXT): root_of_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)root_of_2 root_of_2$(OBJ_EXT) $(LDFLAGS)

rounding_modes$(EXE_EXT): rounding_modes$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rounding_modes rounding_modes$(OBJ_EXT) $(LDFLAGS)

simplest_rational$(EXE_EXT): simplest_rational$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)simplest_rational simplest_rational$(OBJ_EXT) $(LDFLAGS)

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
                   constant.clean \
                   doubletst.clean \
                   floattst.clean \
                   Gmpq.clean \
                   Interval_nt.clean \
                   known_bit_size_integers.clean \
                   Lazy_exact_nt.clean \
                   long_doubletst.clean \
                   MP_Float.clean \
                   NT_checker.clean \
                   Quotient.clean \
                   quotient_io.clean \
                   root_of_2.clean \
                   rounding_modes.clean \
                   simplest_rational.clean \
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

