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
           -Iinclude \
           $(CGAL_CXXFLAGS) \
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
                MP_Float$(EXE_EXT) \
                _test_sign$(EXE_EXT) \
                _test_to_interval$(EXE_EXT) \
                _test_valid_finite_double$(EXE_EXT) \
                _test_valid_finite_float$(EXE_EXT) \
                cln_io$(EXE_EXT) \
                doubletst$(EXE_EXT) \
                floattst$(EXE_EXT) \
                quotient_io$(EXE_EXT) 

MP_Float$(EXE_EXT): MP_Float$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MP_Float MP_Float$(OBJ_EXT) $(LDFLAGS)

_test_sign$(EXE_EXT): _test_sign$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)_test_sign _test_sign$(OBJ_EXT) $(LDFLAGS)

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

clean: \
                   MP_Float.clean \
                   _test_sign.clean \
                   _test_to_interval.clean \
                   _test_valid_finite_double.clean \
                   _test_valid_finite_float.clean \
                   cln_io.clean \
                   doubletst.clean \
                   floattst.clean \
                   quotient_io.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

