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
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(EXTRA_DEBUG_OPT) \
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
                test01$(EXE_EXT) \
                test02$(EXE_EXT) \
                test03$(EXE_EXT) \
                test04$(EXE_EXT) \
                test05$(EXE_EXT) 

test01$(EXE_EXT): test01$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test01 test01$(OBJ_EXT) $(LDFLAGS)

test02$(EXE_EXT): test02$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test02 test02$(OBJ_EXT) $(LDFLAGS)

test03$(EXE_EXT): test03$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test03 test03$(OBJ_EXT) $(LDFLAGS)

test04$(EXE_EXT): test04$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test04 test04$(OBJ_EXT) $(LDFLAGS)

test05$(EXE_EXT): test05$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test05 test05$(OBJ_EXT) $(LDFLAGS)

clean:
	rm -f test??$(EXE_EXT) test??$(OBJ_EXT) core error.txt aresult result.*
	rm -f -r ii_files

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

