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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
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
                dft \
                ist \
                unt 

dft$(EXE_EXT): dft$(OBJ_EXT) Bops_io1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)dft dft$(OBJ_EXT) Bops_io1$(OBJ_EXT) $(LDFLAGS)

ist$(EXE_EXT): ist$(OBJ_EXT) Bops_io2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ist ist$(OBJ_EXT) Bops_io2$(OBJ_EXT) $(LDFLAGS)

unt$(EXE_EXT): unt$(OBJ_EXT) Bops_io3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)unt unt$(OBJ_EXT) Bops_io3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   dft.clean Bops_io1.clean \
                   ist.clean Bops_io2.clean \
                   unt.clean Bops_io3.clean
	rm -f *$(OBJ_EXT)

clobber: clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(EXTRAFLAGS) $(CXXFLAGS) $(OBJ_OPT) $<

