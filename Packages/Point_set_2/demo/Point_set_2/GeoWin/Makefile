include $(CGAL_MAKEFILE)

CC   = $(CGAL_CXX) $(CGAL_WINDOW_LIBPATH)

CGALCFLAGS =  $(CGAL_CXXFLAGS)


CGDL = $(CGAL_GEOWIN_LDFLAGS)

geowin_pset_2d$(EXE_EXT) :  geowin_pset_2d$(OBJ_EXT)
	$(CC) $(EXE_OPT)geowin_pset_2d geowin_pset_2d$(OBJ_EXT) $(CGDL)

geowin_k_nearest_neighbors$(EXE_EXT) : geowin_k_nearest_neighbors$(OBJ_EXT)
	$(CC) $(EXE_OPT)geowin_k_nearest_neighbors geowin_k_nearest_neighbors$(OBJ_EXT) $(CGDL)

	
clean :
	rm -f  tt *~ *.o *.obj *.exe
	rm -f  geowin_k_nearest_neighbors geowin_pset_2d
	

all: geowin_k_nearest_neighbors geowin_pset_2d
	
	
cleano :
	rm -f *~ *.o


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) -c $<

.c$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) -c $<

