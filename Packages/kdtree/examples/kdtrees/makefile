# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

TESTSUITE_CXXFLAGS =

CXXFLAGS = -I../../include\
           $(TESTSUITE_CXXFLAGS) \
	   $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
	   $(LONG_NAME_PROBLEM_CXXFLAGS)

#long names flag for Sun+egcs platform

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS) \
	  $(LONG_NAME_PROBLEM_LDFLAGS)

#long names flag for Sun+egcs platform

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                example1 \
                example2 \
                example3
#                example4 \
#		example5 \
#		example6

example1: example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example1 example1$(OBJ_EXT) $(LDFLAGS)

example2: example2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example2 example2$(OBJ_EXT) $(LDFLAGS)

example3: example3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example3 example3$(OBJ_EXT) $(LDFLAGS)

#example4: example4$(OBJ_EXT)
#	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example4 example4$(OBJ_EXT) $(LDFLAGS)

#example5: example5$(OBJ_EXT)
#	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example5 example5$(OBJ_EXT) $(LDFLAGS)

#example6: example6$(OBJ_EXT)
#	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example6 example6$(OBJ_EXT) $(LDFLAGS)


clean:
	rm -f -r *$(OBJ_EXT) \
	 	   ii_files \
                   example1 \
                   example2 \
                   example3 \
                   example4 \
		   example5 \
		   example6 \
                   core \
		   *.exe *.tds

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c $<

