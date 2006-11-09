#
# $URL$
# $Id$
# 

# This is the makefile for compiling the CGAL GB/RS support.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# The following line should better be passed as command line argument.
# CGAL_MAKEFILE = SHOULD_BE_SET_BY_INSTALL_CGAL
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# CXXFLAGS needs to be passed as command line argument
# CXXFLAGS = $(CGAL_SHARED_LIB_CXXFLAGS)

#---------------------------------------------------------------------#
#                    object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
Gbrs_algebraic_1$(OBJ_EXT) \
Gbrs_polynomial_1$(OBJ_EXT) \
Gbrs_solve_1$(OBJ_EXT) \
Gbrs_polynomial_2$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

#--------------------------------------------------
# shared_lib: $(CGAL_OBJECTS)
# 	$(CGAL_SHARED_LIB_CREATE)$(CGAL_SHARED_LIB) \
# 	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
# 		$(CGAL_SHARED_LIB_LDFLAGS)
# 	mv $(CGAL_SHARED_LIB) '$(CGAL_LIB_DESTINATION)'
# 	rm $(CGAL_OBJECTS)
#-------------------------------------------------- 

static_lib: static_lib_no_install
	mv $(CGAL_GBRS_LIB) $(CGAL_LIB_DESTINATION)

static_lib_no_install: $(CGAL_OBJECTS)
	$(CGAL_LIB_CREATE)$(CGAL_GBRS_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_GBRS_LIB)
	rm $(CGAL_OBJECTS)

clean:
	rm -f $(CGAL_GBRS_LIB) $(CGAL_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.SUFFIXES: .cpp $(SUFFIXES)

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS)  -c $<

