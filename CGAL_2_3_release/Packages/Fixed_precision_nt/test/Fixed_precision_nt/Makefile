# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.


include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS =   $(DEBUG_OPT) -I../../include $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(CGAL_LDFLAGS)

LIBPATH = $(TESTSUITE_LIBPATH) $(CGAL_LIBPATH)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

error.txt:
	cgal_test

test_integers: test_integers$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_integers \
	test_integers$(OBJ_EXT)   $(LDFLAGS)

test_not_integers: test_not_integers$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_not_integers \
	test_not_integers$(OBJ_EXT)   $(LDFLAGS)

clean: test_not_integers.clean test_integers.clean

clean_all:
	rm -f ProgramOutput.* error.txt

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#


.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

