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
           -I$(HOME)/CGAL/Arrangement_2/include \
	   -I$(HOME)/CGAL/Curved_kernel/include \
	   -I$(HOME)/CGAL/Algebraic_kernel/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(DEBUG_OPT)

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
		ex_arr$(EXE_EXT) \
                example1$(EXE_EXT) \
                example2$(EXE_EXT) \
                example3$(EXE_EXT) \
                example4$(EXE_EXT) \
                example5$(EXE_EXT) 

ex_arr$(EXE_EXT): ex_arr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ex_arr ex_arr$(OBJ_EXT) $(LDFLAGS)

example1$(EXE_EXT): example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example1 example1$(OBJ_EXT) $(LDFLAGS)

example2$(EXE_EXT): example2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example2 example2$(OBJ_EXT) $(LDFLAGS)

example3$(EXE_EXT): example3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example3 example3$(OBJ_EXT) $(LDFLAGS)

example4$(EXE_EXT): example4$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example4 example4$(OBJ_EXT) $(LDFLAGS)

example5$(EXE_EXT): example5$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example5 example5$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example1.clean \
                   example2.clean \
                   example3.clean \
                   example4.clean \
                   example5.clean \
		   ex_arr.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

