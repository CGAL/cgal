# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
           -I./include \
           $(TESTSUITE_CXXFLAGS) \
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
                test01 \
                test02 \
                test03 \
                test04 \
                test05

test01:	test01$(EXE_OPT)

test02:	test02$(EXE_OPT)

test03:	test03$(EXE_OPT)

test04:	test04$(EXE_OPT)

test05:	test05$(EXE_OPT)

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
        $(RM_FORCE) $(RMFLAGS) test??$(EXE_EXT) test??$(OBJ_EXT) aresult

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
$(OBJ_EXT)$(EXE_EXT):
	$(CGAL_CXX) $(EXE_OPT)$@ $< $(LDFLAGS)
.C$(EXE_EXT): 
	$(CGAL_CXX) $(EXE_OPT)$@ $(CXXFLAGS) $< $(LDFLAGS)
