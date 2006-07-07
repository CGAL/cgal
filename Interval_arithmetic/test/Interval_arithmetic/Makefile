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
                Filtered_kernel$(EXE_EXT) \
                Lazy_exact_nt$(EXE_EXT) \
                Static_filters$(EXE_EXT) \
                tst10$(EXE_EXT) \
                tst5$(EXE_EXT) 

Filtered_kernel$(EXE_EXT): Filtered_kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Filtered_kernel Filtered_kernel$(OBJ_EXT) $(LDFLAGS)

Lazy_exact_nt$(EXE_EXT): Lazy_exact_nt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Lazy_exact_nt Lazy_exact_nt$(OBJ_EXT) $(LDFLAGS)

Static_filters$(EXE_EXT): Static_filters$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Static_filters Static_filters$(OBJ_EXT) $(LDFLAGS)

tst10$(EXE_EXT): tst10$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst10 tst10$(OBJ_EXT) $(LDFLAGS)

tst5$(EXE_EXT): tst5$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst5 tst5$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Filtered_kernel.clean \
                   Lazy_exact_nt.clean \
                   Static_filters.clean \
                   tst10.clean \
                   tst5.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

