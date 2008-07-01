# Created by the script cgal_create_makefile
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
	-g \
	-I$(HOME)/EXACUS/Support/include\
	-I$(HOME)/EXACUS/NumeriX/include\
	-I$(HOME)/EXACUS/SweepX/include\
	-I$(HOME)/EXACUS/AlciX/include\
	-I$(HOME)/trunk/Installation/include/ \
	-I$(HOME)/trunk/Algebraic_foundations/include/ \
	-I$(HOME)/trunk/STL_Extension/include/ \
	-I$(HOME)/trunk/Circulator/include/ \
	-I$(HOME)/trunk/HalfedgeDS/include/ \
	-I$(HOME)/trunk/Number_types/include/ \
	-I$(HOME)/trunk/Polynomial/include \
	-I$(HOME)/trunk/Algebraic_kernel_d/include \
	-I$(HOME)/trunk/Curved_kernel_via_analysis_2/include \
	-I$(HOME)/trunk/Arrangement_on_surface_2/include/ \
	-I$(HOME)/trunk/Envelope_3/include/ \
        -Iinclude \
        -I../../include \
        -I../../../Number_types/test/Number_types/include \
	   $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
	$(LONG_NAME_PROBLEM_LDFLAGS) \
        $(CGAL_LDFLAGS) \
	/usr/local/gcc-4.1/ntl-5.4/lib/libntl.a \
	$(HOME)/EXACUS/NumeriX/src/.libs/libNiX.so 


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                Algebraic_curve_kernel_2$(EXE_EXT) \
                Curve_analysis_2$(EXE_EXT) \
                Curve_pair_analysis_2$(EXE_EXT) \
		Bitstream_descartes_bfs$(EXE_EXT) \
		alg_real_utils$(EXE_EXT) \
		Best_approximation_cache$(EXE_EXT)


Algebraic_curve_kernel_2$(EXE_EXT): Algebraic_curve_kernel_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Algebraic_curve_kernel_2 Algebraic_curve_kernel_2$(OBJ_EXT) $(LDFLAGS)

Curve_analysis_2$(EXE_EXT): Curve_analysis_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Curve_analysis_2 Curve_analysis_2$(OBJ_EXT) $(LDFLAGS)

Curve_pair_analysis_2$(EXE_EXT): Curve_pair_analysis_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Curve_pair_analysis_2 Curve_pair_analysis_2$(OBJ_EXT) $(LDFLAGS)

subresultants$(EXE_EXT): subresultants$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)subresultants subresultants$(OBJ_EXT) $(LDFLAGS)

sturm_habicht_sequence$(EXE_EXT): sturm_habicht_sequence$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)sturm_habicht_sequence sturm_habicht_sequence$(OBJ_EXT) $(LDFLAGS)

Bitstream_descartes_bfs$(EXE_EXT): Bitstream_descartes_bfs$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Bitstream_descartes_bfs Bitstream_descartes_bfs$(OBJ_EXT) $(LDFLAGS)

alg_real_utils$(EXE_EXT): alg_real_utils$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)alg_real_utils alg_real_utils$(OBJ_EXT) $(LDFLAGS)

Best_approximation_cache$(EXE_EXT): Best_approximation_cache$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Best_approximation_cache Best_approximation_cache$(OBJ_EXT) $(LDFLAGS)


clean: \
                   Algebraic_curve_kernel_2.clean \
                   Algebraic_kernel_1.clean \
                   Algebraic_real_pure.clean \
                   Bitstream_descartes.clean \
                   Descartes.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

