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
                bench_IA$(EXE_EXT) \
                Filtered_kernel$(EXE_EXT) \
                gcc_3.0_bug$(EXE_EXT) \
                Interval_nt$(EXE_EXT) \
                Lazy_exact_nt$(EXE_EXT) \
                Static_filters$(EXE_EXT) \
                to_interval_test$(EXE_EXT) \
                tst0$(EXE_EXT) \
                tst10$(EXE_EXT) \
                tst5$(EXE_EXT) 

bench_IA$(EXE_EXT): bench_IA$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bench_IA bench_IA$(OBJ_EXT) $(LDFLAGS)

Filtered_kernel$(EXE_EXT): Filtered_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Filtered_kernel Filtered_kernel$(OBJ_EXT) $(LDFLAGS)

gcc_3.0_bug$(EXE_EXT): gcc_3.0_bug$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gcc_3.0_bug gcc_3.0_bug$(OBJ_EXT) $(LDFLAGS)

Interval_nt$(EXE_EXT): Interval_nt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Interval_nt Interval_nt$(OBJ_EXT) $(LDFLAGS)

Lazy_exact_nt$(EXE_EXT): Lazy_exact_nt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_exact_nt Lazy_exact_nt$(OBJ_EXT) $(LDFLAGS)

Static_filters$(EXE_EXT): Static_filters$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Static_filters Static_filters$(OBJ_EXT) $(LDFLAGS)

to_interval_test$(EXE_EXT): to_interval_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)to_interval_test to_interval_test$(OBJ_EXT) $(LDFLAGS)

tst0$(EXE_EXT): tst0$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst0 tst0$(OBJ_EXT) $(LDFLAGS)

tst10$(EXE_EXT): tst10$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst10 tst10$(OBJ_EXT) $(LDFLAGS)

tst5$(EXE_EXT): tst5$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst5 tst5$(OBJ_EXT) $(LDFLAGS)

clean: \
                   bench_IA.clean \
                   Filtered_kernel.clean \
                   gcc_3.0_bug.clean \
                   Interval_nt.clean \
                   Lazy_exact_nt.clean \
                   Static_filters.clean \
                   to_interval_test.clean \
                   tst0.clean \
                   tst10.clean \
                   tst5.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

