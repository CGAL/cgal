# This is the makefile for compiling the Core library within CGAL.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
# CGAL_MAKEFILE = SHOULD_BE_SET_BY_COMMAND_LINE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# CXXFLAGS must be set on the make command-line
# CXXFLAGS = -I../../include  $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
CoreAux$(OBJ_EXT) \
extLong$(OBJ_EXT) \
CoreDefs$(OBJ_EXT) \
CoreIO$(OBJ_EXT) \
BigFloat$(OBJ_EXT) \
Real$(OBJ_EXT) \
Expr$(OBJ_EXT) \
GmpIO$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

# The shared_lib case is still buggy:
# - Currently MPFR is linked with it (through GMP_LDFLAGS), although it is
#   not needed, but anyway next Core release will need it...
# - When GMP or MPFR is a static library on Darwin, it fails because -fno-common is required...
# - The Darwin install_name is buggy (it's libCGAL)

shared_lib: $(CGAL_OBJECTS)
	$(CGAL_SHARED_LIB_CREATE)$(CGAL_CORE_SHARED_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_SHARED_LIB_LDFLAGS) $(GMP_LIBPATHFLAGS) $(MPFR_LIBPATHFLAGS) $(GMP_LDFLAGS)
	mv $(CGAL_CORE_SHARED_LIB) '$(CGAL_LIB_DESTINATION)'
	rm $(CGAL_OBJECTS)

static_lib: static_lib_no_install
	mv $(CGAL_CORE_LIB) $(CGAL_LIB_DESTINATION)

static_lib_no_install: $(CGAL_OBJECTS)
	$(CGAL_LIB_CREATE)$(CGAL_CORE_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_CORE_LIB)
	rm $(CGAL_OBJECTS)

clean:
	rm -f $(CGAL_CORE_LIB) $(CGAL_CORE_SHARED_LIB) $(CGAL_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.SUFFIXES: .cpp $(SUFFIXES)

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS) -c $<

