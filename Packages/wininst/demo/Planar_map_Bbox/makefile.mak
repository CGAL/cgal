# This is the makefile for compiling a CGAL application.
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

OBJ_EXT = .wrong_cgal_makefile

#CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
        $(EXTRA_FLAGS) \
	-I./include \
	-I../../include/ \
	$(CGAL_CXXFLAGS) \
        $(LONG_NAME_PROBLEM_CXXFLAGS) \
        $(DEBUG_OPT) \
        $(EXTRA_DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:	demo voronoi

demo$(EXE_EXT): demo$(OBJ_EXT) draw_map$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) draw_map$(OBJ_EXT) $(LDFLAGS)

voronoi$(EXE_EXT): voronoi$(OBJ_EXT) draw_map$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)voronoi voronoi$(OBJ_EXT) draw_map$(OBJ_EXT) $(LDFLAGS)

###########################################################################
#The following lines compile the demo with different options, we use
#them for experiments only.

demo.all:	\
		demo.float.rebuild\
		demo.rational.rebuild\
		demo.leda.rational.rebuild\
		demo.float.no-rebuild\
		demo.rational.no-rebuild\
		demo.leda.rational.no-rebuild\
		demo.float.walk\
		demo.rational.walk\
		demo.leda.rational.walk\
		demo.float.naive\
		demo.rational.naive\
		demo.leda.rational.naive

# executables ... default strategy
demo.float.rebuild$(EXE_EXT): demo.float.rebuild$(OBJ_EXT) draw_map.float.rebuild$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.float.rebuild demo.float.rebuild$(OBJ_EXT) draw_map.float.rebuild$(OBJ_EXT) $(LDFLAGS)

demo.rational.rebuild$(EXE_EXT): demo.rational.rebuild$(OBJ_EXT) draw_map.rational.rebuild$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.rational.rebuild demo.rational.rebuild$(OBJ_EXT) draw_map.rational.rebuild$(OBJ_EXT) $(LDFLAGS)

demo.leda.rational.rebuild$(EXE_EXT): demo.leda.rational.rebuild$(OBJ_EXT) draw_map.leda.rational.rebuild$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.leda.rational.rebuild demo.leda.rational.rebuild$(OBJ_EXT) draw_map.leda.rational.rebuild$(OBJ_EXT) $(LDFLAGS)

# executables ... default strategy without rebuilds
demo.float.no-rebuild$(EXE_EXT): demo.float.no-rebuild$(OBJ_EXT) draw_map.float.no-rebuild$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.float.no-rebuild demo.float.no-rebuild$(OBJ_EXT) draw_map.float.no-rebuild$(OBJ_EXT) $(LDFLAGS)

demo.rational.no-rebuild$(EXE_EXT): demo.rational.no-rebuild$(OBJ_EXT) draw_map.rational.no-rebuild$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.rational.no-rebuild demo.rational.no-rebuild$(OBJ_EXT) draw_map.rational.no-rebuild$(OBJ_EXT) $(LDFLAGS)

demo.leda.rational.no-rebuild$(EXE_EXT): demo.leda.rational.no-rebuild$(OBJ_EXT) draw_map.leda.rational.no-rebuild$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.leda.rational.no-rebuild demo.leda.rational.no-rebuild$(OBJ_EXT) draw_map.leda.rational.no-rebuild$(OBJ_EXT) $(LDFLAGS)

# executables ... walk strategy
demo.float.walk$(EXE_EXT): demo.float.walk$(OBJ_EXT) draw_map.float.walk$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.float.walk demo.float.walk$(OBJ_EXT) draw_map.float.walk$(OBJ_EXT) $(LDFLAGS)

demo.rational.walk$(EXE_EXT): demo.rational.walk$(OBJ_EXT) draw_map.rational.walk$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.rational.walk demo.rational.walk$(OBJ_EXT) draw_map.rational.walk$(OBJ_EXT) $(LDFLAGS)

demo.leda.rational.walk$(EXE_EXT): demo.leda.rational.walk$(OBJ_EXT) draw_map.leda.rational.walk$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.leda.rational.walk demo.leda.rational.walk$(OBJ_EXT) draw_map.leda.rational.walk$(OBJ_EXT) $(LDFLAGS)

# executables ... naive strategy 
demo.float.naive$(EXE_EXT): demo.float.naive$(OBJ_EXT) draw_map.float.naive$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.float.naive demo.float.naive$(OBJ_EXT) draw_map.float.naive$(OBJ_EXT) $(LDFLAGS)

demo.rational.naive$(EXE_EXT): demo.rational.naive$(OBJ_EXT) draw_map.rational.naive$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.rational.naive demo.rational.naive$(OBJ_EXT) draw_map.rational.naive$(OBJ_EXT) $(LDFLAGS)

demo.leda.rational.naive$(EXE_EXT): demo.leda.rational.naive$(OBJ_EXT) draw_map.leda.rational.naive$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo.leda.rational.naive demo.leda.rational.naive$(OBJ_EXT) draw_map.leda.rational.naive$(OBJ_EXT) $(LDFLAGS)

draw_map_h = draw_map.h ../../include/CGAL/Pm_straight_exact_traits.h

voronoi_h = voronoi.h ../../include/CGAL/Pm_dynamic_open_bounding_box.h ../../include/CGAL/Pm_dynamic_closed_bounding_box.h

# demo object files
demo.float.rebuild$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) $(EXE_OPT)demo.float.rebuild$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.rational.rebuild$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_RATIONAL $(EXE_OPT)demo.rational.rebuild$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.leda.rational.rebuild$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)demo.leda.rational.rebuild$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.float.no-rebuild$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_DEFAULT_WITHOUT_REBUILD $(EXE_OPT)demo.float.no-rebuild$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.rational.no-rebuild$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_DEFAULT_WITHOUT_REBUILD -DUSE_RATIONAL $(EXE_OPT)demo.rational.no-rebuild$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.leda.rational.no-rebuild$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_DEFAULT_WITHOUT_REBUILD -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)demo.leda.rational.no-rebuild$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.float.walk$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_WALK_POINT_LOCATION $(EXE_OPT)demo.float.walk$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.rational.walk$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_WALK_POINT_LOCATION -DUSE_RATIONAL $(EXE_OPT)demo.rational.walk$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.leda.rational.walk$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_WALK_POINT_LOCATION -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)demo.leda.rational.walk$(OBJ_EXT) $(OBJ_OPT) demo.C 


demo.float.naive$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_NAIVE_POINT_LOCATION $(EXE_OPT)demo.float.naive$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.rational.naive$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_NAIVE_POINT_LOCATION -DUSE_RATIONAL $(EXE_OPT)demo.rational.naive$(OBJ_EXT) $(OBJ_OPT) demo.C 

demo.leda.rational.naive$(OBJ_EXT): demo.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_NAIVE_POINT_LOCATION -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)demo.leda.rational.naive$(OBJ_EXT) $(OBJ_OPT) demo.C 


#draw_map object file
draw_map.float.rebuild$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) $(EXE_OPT)draw_map.float.rebuild$(OBJ_EXT) $(OBJ_OPT) draw_map.C

draw_map.rational.rebuild$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_RATIONAL $(EXE_OPT)draw_map.rational.rebuild$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

draw_map.leda.rational.rebuild$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)draw_map.leda.rational.rebuild$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

draw_map.float.no-rebuild$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_DEFAULT_WITHOUT_REBUILD $(EXE_OPT)draw_map.float.no-rebuild$(OBJ_EXT) $(OBJ_OPT) draw_map.C

draw_map.rational.no-rebuild$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_DEFAULT_WITHOUT_REBUILD -DUSE_RATIONAL $(EXE_OPT)draw_map.rational.no-rebuild$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

draw_map.leda.rational.no-rebuild$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_DEFAULT_WITHOUT_REBUILD -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)draw_map.leda.rational.no-rebuild$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

draw_map.float.walk$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_WALK_POINT_LOCATION $(EXE_OPT)draw_map.float.walk$(OBJ_EXT) $(OBJ_OPT) draw_map.C

draw_map.rational.walk$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_WALK_POINT_LOCATION -DUSE_RATIONAL $(EXE_OPT)draw_map.rational.walk$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

draw_map.leda.rational.walk$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_WALK_POINT_LOCATION -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)draw_map.leda.rational.walk$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

draw_map.float.naive$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_NAIVE_POINT_LOCATION $(EXE_OPT)draw_map.float.naive$(OBJ_EXT) $(OBJ_OPT) draw_map.C

draw_map.rational.naive$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_NAIVE_POINT_LOCATION -DUSE_RATIONAL $(EXE_OPT)draw_map.rational.naive$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

draw_map.leda.rational.naive$(OBJ_EXT): draw_map.C draw_map.h makefile 
	$(CGAL_CXX) $(CXXFLAGS) -DUSE_NAIVE_POINT_LOCATION -DUSE_LEDA_RAT_KERNEL $(EXE_OPT)draw_map.leda.rational.naive$(OBJ_EXT) $(OBJ_OPT) draw_map.C 

voronoi$(OBJ_EXT): voronoi.C $(voronoi_h) $(draw_map_h)
	$(CGAL_CXX) $(CXXFLAGS) $(EXE_OPT)voronoi$(OBJ_EXT) $(OBJ_OPT) voronoi.C 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

clean: voronoi.clean demo.clean draw_map.clean

cleanall:
	rm -f -r *$(OBJ_EXT) demo$(EXE_EXT) \
	demo.float.rebuild$(EXE_EXT) \
	demo.rational.rebuild$(EXE_EXT) \
	demo.leda.rational.rebuild$(EXE_EXT) \
	demo.float.no-rebuild$(EXE_EXT) \
	demo.rational.no-rebuild$(EXE_EXT) \
	demo.leda.rational.no-rebuild$(EXE_EXT) \
	demo.float.walk$(EXE_EXT) \
	demo.rational.walk$(EXE_EXT) \
	demo.leda.rational.walk$(EXE_EXT) \
	demo.float.naive$(EXE_EXT) \
	demo.rational.naive$(EXE_EXT) \
	demo.leda.rational.naive$(EXE_EXT) \
	voronoi$(EXE_EXT) \
	helputil$(EXE_EXT) \
	ii_files \
	core

# note - we use many compile flag dependent executables for optimization
