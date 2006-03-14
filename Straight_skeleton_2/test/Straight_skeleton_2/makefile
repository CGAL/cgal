# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE=c:/cgal/CGAL-3.2-I-361/make/makefile_i686_CYGWINNT-5.1_CL.EXE-1310
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           -O2

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

all: \
  test_sls_traits$(EXE_EXT) \
  test_sls_builder$(EXE_EXT) \
  test_offset_builder$(EXE_EXT) \

test_sls_traits$(EXE_EXT): test_sls_traits$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) test_sls_traits test_sls_traits$(OBJ_EXT) $(LDFLAGS)

test_sls_builder$(EXE_EXT): test_sls_builder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) test_sls_builder test_sls_builder$(OBJ_EXT) $(LDFLAGS)

test_offset_builder$(EXE_EXT): test_offset_builder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) test_offset_builder test_offset_builder$(OBJ_EXT) $(LDFLAGS)

clean: \
        test_sls_traits.clean \
        test_sls_builder.clean \
        test_offset_builder.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

