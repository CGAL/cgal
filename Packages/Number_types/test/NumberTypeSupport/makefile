# Created by the script create_makefile
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
           -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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
                Gmpq$(EXE_EXT) \
                MP_Float$(EXE_EXT) \
                _test_to_interval$(EXE_EXT) \
                _test_valid_finite_double$(EXE_EXT) \
                _test_valid_finite_float$(EXE_EXT) \
                cln_io$(EXE_EXT) \
                doubletst$(EXE_EXT) \
                floattst$(EXE_EXT) \
                quotient_io$(EXE_EXT) \
                simplest_rational$(EXE_EXT) \
                utilities$(EXE_EXT) 

Gmpq$(EXE_EXT): Gmpq$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Gmpq Gmpq$(OBJ_EXT) $(LDFLAGS)

MP_Float$(EXE_EXT): MP_Float$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MP_Float MP_Float$(OBJ_EXT) $(LDFLAGS)

_test_to_interval$(EXE_EXT): _test_to_interval$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)_test_to_interval _test_to_interval$(OBJ_EXT) $(LDFLAGS)

_test_valid_finite_double$(EXE_EXT): _test_valid_finite_double$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)_test_valid_finite_double _test_valid_finite_double$(OBJ_EXT) $(LDFLAGS)

_test_valid_finite_float$(EXE_EXT): _test_valid_finite_float$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)_test_valid_finite_float _test_valid_finite_float$(OBJ_EXT) $(LDFLAGS)

cln_io$(EXE_EXT): cln_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)cln_io cln_io$(OBJ_EXT) $(LDFLAGS)

doubletst$(EXE_EXT): doubletst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)doubletst doubletst$(OBJ_EXT) $(LDFLAGS)

floattst$(EXE_EXT): floattst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)floattst floattst$(OBJ_EXT) $(LDFLAGS)

quotient_io$(EXE_EXT): quotient_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quotient_io quotient_io$(OBJ_EXT) $(LDFLAGS)

simplest_rational$(EXE_EXT): simplest_rational$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)simplest_rational simplest_rational$(OBJ_EXT) $(LDFLAGS)

utilities$(EXE_EXT): utilities$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)utilities utilities$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Gmpq.clean \
                   MP_Float.clean \
                   _test_to_interval.clean \
                   _test_valid_finite_double.clean \
                   _test_valid_finite_float.clean \
                   cln_io.clean \
                   doubletst.clean \
                   floattst.clean \
                   quotient_io.clean \
                   simplest_rational.clean \
                   utilities.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

