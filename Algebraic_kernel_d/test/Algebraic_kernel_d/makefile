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
           -Iinclude \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                Algebraic_curve_kernel_2$(EXE_EXT) \
                Algebraic_kernel_1$(EXE_EXT) \
                Algebraic_real_pure$(EXE_EXT) \
                Bitstream_descartes$(EXE_EXT) \
                Descartes$(EXE_EXT) 

Algebraic_curve_kernel_2$(EXE_EXT): Algebraic_curve_kernel_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Algebraic_curve_kernel_2 Algebraic_curve_kernel_2$(OBJ_EXT) $(LDFLAGS)

Algebraic_kernel_1$(EXE_EXT): Algebraic_kernel_1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Algebraic_kernel_1 Algebraic_kernel_1$(OBJ_EXT) $(LDFLAGS)

Algebraic_real_pure$(EXE_EXT): Algebraic_real_pure$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Algebraic_real_pure Algebraic_real_pure$(OBJ_EXT) $(LDFLAGS)

Bitstream_descartes$(EXE_EXT): Bitstream_descartes$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Bitstream_descartes Bitstream_descartes$(OBJ_EXT) $(LDFLAGS)

Descartes$(EXE_EXT): Descartes$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Descartes Descartes$(OBJ_EXT) $(LDFLAGS)

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

