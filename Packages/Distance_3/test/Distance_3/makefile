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
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                tst321 \
                tst322 \
                tst323 \
                tst324 \
                tst325 \
                tst326 \
                tst327 \
                tst328 \
                tst329 \
                tst330 \
                tst331 \
                tst332 \
                tst333 

tst321$(EXE_EXT): tst321$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst321 tst321$(OBJ_EXT) $(LDFLAGS)

tst322$(EXE_EXT): tst322$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst322 tst322$(OBJ_EXT) $(LDFLAGS)

tst323$(EXE_EXT): tst323$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst323 tst323$(OBJ_EXT) $(LDFLAGS)

tst324$(EXE_EXT): tst324$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst324 tst324$(OBJ_EXT) $(LDFLAGS)

tst325$(EXE_EXT): tst325$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst325 tst325$(OBJ_EXT) $(LDFLAGS)

tst326$(EXE_EXT): tst326$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst326 tst326$(OBJ_EXT) $(LDFLAGS)

tst327$(EXE_EXT): tst327$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst327 tst327$(OBJ_EXT) $(LDFLAGS)

tst328$(EXE_EXT): tst328$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst328 tst328$(OBJ_EXT) $(LDFLAGS)

tst329$(EXE_EXT): tst329$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst329 tst329$(OBJ_EXT) $(LDFLAGS)

tst330$(EXE_EXT): tst330$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst330 tst330$(OBJ_EXT) $(LDFLAGS)

tst331$(EXE_EXT): tst331$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst331 tst331$(OBJ_EXT) $(LDFLAGS)

tst332$(EXE_EXT): tst332$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst332 tst332$(OBJ_EXT) $(LDFLAGS)

tst333$(EXE_EXT): tst333$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst333 tst333$(OBJ_EXT) $(LDFLAGS)

clean: \
                   tst321.clean \
                   tst322.clean \
                   tst323.clean \
                   tst324.clean \
                   tst325.clean \
                   tst326.clean \
                   tst327.clean \
                   tst328.clean \
                   tst329.clean \
                   tst330.clean \
                   tst331.clean \
                   tst332.clean \
                   tst333.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

