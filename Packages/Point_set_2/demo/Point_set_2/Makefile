include $(CGAL_MAKEFILE)

CC   = $(CGAL_CXX) $(CGAL_WINDOW_LIBPATH)

CGALCFLAGS =  $(CGAL_CXXFLAGS) 

CGDL =  $(CGAL_WINDOW_LDFLAGS)


ps_test1_cgal$(EXE_EXT) : ps_test1_cgal$(OBJ_EXT)
	$(CC)  $(EXE_OPT)ps_test1_cgal ps_test1_cgal$(OBJ_EXT) $(CGDL)

nearest_neighbor$(EXE_EXT) :  nearest_neighbor$(OBJ_EXT)
	$(CC)  $(EXE_OPT)nearest_neighbor nearest_neighbor$(OBJ_EXT) $(CGDL)

clean :
	rm -f  tt *~ *.o *.obj *.exe core
	rm -f  nearest_neighbor ps_test1_cgal

all: nearest_neighbor ps_test1_cgal 

cleano :
	rm -f *~ *.o


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) $(OBJ_OPT) $<


