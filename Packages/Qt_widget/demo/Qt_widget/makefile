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

CXXFLAGS = -I../../include \
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
          $(CGAL_QT_LDFLAGS)
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

DEMO_SUBDIRS= \
./Constrained_triangulation_2 \
./Convex_Hull_2 \
./Generator \
./Max_k-gon \
./Triangulation_2_constrained \
./Examples \
./basic/tutorial1 \
./basic/tutorial2 \
./basic/tutorial3 \
./basic/tutorial4 \
./basic/tutorial5 \
./basic/tutorial6 \
./basic/tutorial7 \
./basic/tutorial8

all:
	badexit=""; \
	for dir in $(DEMO_SUBDIRS); do \
	  (cd $$dir && make ) || badexit="yes";\
	done; \
        test -z "$$badexit"

$(DEMO_SUBDIRS):
	cd $@ && make

#	for dir in ${DEMO_SUBDIRS}; do make -C $$dir; done

clean:
	for dir in ${DEMO_SUBDIRS}; do make -C $$dir clean; done

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

