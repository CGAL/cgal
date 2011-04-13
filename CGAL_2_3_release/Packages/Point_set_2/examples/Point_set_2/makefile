include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(CGAL_CXXFLAGS)
  


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
                range_search_tr \
                nearest_nb1_tr \
		rs_example 

range_search_tr$(EXE_EXT) : range_search_tr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_search_tr range_search_tr$(OBJ_EXT) $(LDFLAGS)

rs_example$(EXE_EXT) : rs_example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rs_example rs_example$(OBJ_EXT) $(LDFLAGS)

nearest_nb1_tr$(EXE_EXT) : nearest_nb1_tr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_nb1_tr nearest_nb1_tr$(OBJ_EXT) $(LDFLAGS)
	
clean: \
	range_search_tr.clean rs_example.clean nearest_nb1_tr.clean

myclean: 
	rm -f range_search_tr rs_example nearest_nb1_tr range_search_tr_leda *$(OBJ_EXT) *.obj *.exe core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

