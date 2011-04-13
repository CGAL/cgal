# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
#CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE

!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

TESTSUITE_CXXFLAGS =

CXXFLAGS = \
	   $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                example1 \
                example2 \
                example3

example1: example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example1 example1$(OBJ_EXT) $(LDFLAGS)

example2: example2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example2 example2$(OBJ_EXT) $(LDFLAGS)

example3: example3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example3 example3$(OBJ_EXT) $(LDFLAGS)


clean:             example1.clean \
                   example2.clean \
                   example3.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c $<

