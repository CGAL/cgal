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
           $(LONG_NAME_PROBLEM_CXXFLAGS) -W

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
                Filtered_kernel$(EXE_EXT) \
                gcc_3.0_bug$(EXE_EXT) \
                Kernel_checker$(EXE_EXT) \
                Lazy_exact_nt$(EXE_EXT) \
                Static_filters$(EXE_EXT) \
                to_interval_test$(EXE_EXT) \
                tst0$(EXE_EXT) \
                tst10$(EXE_EXT) \
                tst1$(EXE_EXT) \
                tst2$(EXE_EXT) \
                tst3$(EXE_EXT) \
                tst4$(EXE_EXT) \
                tst5$(EXE_EXT) \
                tst6$(EXE_EXT) 

Filtered_kernel$(EXE_EXT): Filtered_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Filtered_kernel Filtered_kernel$(OBJ_EXT) $(LDFLAGS)

gcc_3.0_bug$(EXE_EXT): gcc_3.0_bug$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gcc_3.0_bug gcc_3.0_bug$(OBJ_EXT) $(LDFLAGS)

Kernel_checker$(EXE_EXT): Kernel_checker$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Kernel_checker Kernel_checker$(OBJ_EXT) $(LDFLAGS)

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

tst1$(EXE_EXT): tst1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst1 tst1$(OBJ_EXT) $(LDFLAGS)

tst2$(EXE_EXT): tst2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst2 tst2$(OBJ_EXT) $(LDFLAGS)

tst3$(EXE_EXT): tst3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst3 tst3$(OBJ_EXT) $(LDFLAGS)

tst4$(EXE_EXT): tst4$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst4 tst4$(OBJ_EXT) $(LDFLAGS)

tst5$(EXE_EXT): tst5$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst5 tst5$(OBJ_EXT) $(LDFLAGS)

tst6$(EXE_EXT): tst6$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst6 tst6$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Filtered_kernel.clean \
                   gcc_3.0_bug.clean \
                   Kernel_checker.clean \
                   Lazy_exact_nt.clean \
                   Static_filters.clean \
                   to_interval_test.clean \
                   tst0.clean \
                   tst10.clean \
                   tst1.clean \
                   tst2.clean \
                   tst3.clean \
                   tst4.clean \
                   tst5.clean \
                   tst6.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

