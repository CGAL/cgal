include $(CGAL_MAKEFILE)

CC   = $(CGAL_CXX) $(CGAL_WINDOW_LIBPATH)

CGALCFLAGS =  $(CGAL_CXXFLAGS)


CGDL = $(CGAL_GEOWIN_LDFLAGS)

geowin_pset_2d_tb$(EXE_EXT) :  geowin_pset_2d_tb$(OBJ_EXT)
	$(CC) $(EXE_OPT)geowin_pset_2d_tb geowin_pset_2d_tb$(OBJ_EXT) $(CGDL)

geowin_k_nearest_neighbors_tb$(EXE_EXT) : geowin_k_nearest_neighbors_tb$(OBJ_EXT)
	$(CC) $(EXE_OPT)geowin_k_nearest_neighbors_tb geowin_k_nearest_neighbors_tb$(OBJ_EXT) $(CGDL)

	
clean :
	rm -f  tt *~ *.o *.obj *.exe
	rm -f  geowin_k_nearest_neighbors_tb geowin_pset_2d_tb
	

all: geowin_k_nearest_neighbors_tb geowin_pset_2d_tb
	
	
cleano :
	rm -f *~ *.o


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) -c $<

.c$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) -c $<

