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
                tst04 \
                tst05 \
                tst08 \
                tst09 \
                tst10 \
                tst11 \
                tst50 \
                tst51 \
                tst52 \
                tst53 \
                tst54 \
                tst56 \
                tst57 \
                tst58 \
                tst59 

tst04$(EXE_EXT): tst04$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst04 tst04$(OBJ_EXT) $(LDFLAGS)

tst05$(EXE_EXT): tst05$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst05 tst05$(OBJ_EXT) $(LDFLAGS)

tst08$(EXE_EXT): tst08$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst08 tst08$(OBJ_EXT) $(LDFLAGS)

tst09$(EXE_EXT): tst09$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst09 tst09$(OBJ_EXT) $(LDFLAGS)

tst10$(EXE_EXT): tst10$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst10 tst10$(OBJ_EXT) $(LDFLAGS)

tst11$(EXE_EXT): tst11$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst11 tst11$(OBJ_EXT) $(LDFLAGS)

tst50$(EXE_EXT): tst50$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst50 tst50$(OBJ_EXT) $(LDFLAGS)

tst51$(EXE_EXT): tst51$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst51 tst51$(OBJ_EXT) $(LDFLAGS)

tst52$(EXE_EXT): tst52$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst52 tst52$(OBJ_EXT) $(LDFLAGS)

tst53$(EXE_EXT): tst53$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst53 tst53$(OBJ_EXT) $(LDFLAGS)

tst54$(EXE_EXT): tst54$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst54 tst54$(OBJ_EXT) $(LDFLAGS)

tst56$(EXE_EXT): tst56$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst56 tst56$(OBJ_EXT) $(LDFLAGS)

tst57$(EXE_EXT): tst57$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst57 tst57$(OBJ_EXT) $(LDFLAGS)

tst58$(EXE_EXT): tst58$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst58 tst58$(OBJ_EXT) $(LDFLAGS)

tst59$(EXE_EXT): tst59$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst59 tst59$(OBJ_EXT) $(LDFLAGS)

clean: \
                   tst04.clean \
                   tst05.clean \
                   tst08.clean \
                   tst09.clean \
                   tst10.clean \
                   tst11.clean \
                   tst50.clean \
                   tst51.clean \
                   tst52.clean \
                   tst53.clean \
                   tst54.clean \
                   tst56.clean \
                   tst57.clean \
                   tst58.clean \
                   tst59.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

