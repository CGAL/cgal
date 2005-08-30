# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = /proj/geometrica/home/CGAL/CGAL-3.1-I-104/make/makefile_i686_Linux-2.4.28-pre3_g++-3.3.3
# penser a enlever io.h de Conic_2.h en repassant a CGAL-I
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
	   -I../../../Arrangement_2/include \
	   -I../../../Curved_kernel/include \
           -I../../../Algebraic_kernel/include \
           -I$(HOME)/synaps/include \
	   $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           -fmessage-length=0 \
	   $(DEBUG_OPT) \
	   -DHAVE_CONFIG_H -g #-DCGAL_CURVED_KERNEL_DEBUG
#	   -DCGAL_CURVED_KERNEL_DEBUG \
           -O2 -DNDEBUG \
           # -I$(HOME)/CGAL/Interval_arithmetic.RootOf_2/include \


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) \
	  -lgmpxx -lmpfr

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                demo$(EXE_EXT) \
                demo_elias$(EXE_EXT) \
		demo_line$(EXE_EXT) \
		demo_viewer$(EXE_EXT)\
		demo_line_circle$(EXE_EXT)
#                demo_conic$(EXE_EXT) 

demo.moc: demo.C
	${QT_MOC} -o demo.moc $<

demo_elias.moc: demo_elias.C
	${QT_MOC} -o demo_elias.moc $<

demo_line_circle.moc: demo_line_circle.C
	${QT_MOC} -o demo_line_circle.moc $<

demo_viewer.moc: demo_viewer.C
	${QT_MOC} -o demo_viewer.moc $<

#demo_conic.moc: demo_conic.C
#	${QT_MOC} -o demo_conic.moc $<

Qt_widget_get_arc.moc: Qt_widget_get_arc.h
	$(QT_MOC) -o Qt_widget_get_arc.moc $<

Qt_widget_get_segment.moc: Qt_widget_get_segment.h
	$(QT_MOC) -o Qt_widget_get_segment.moc $<

#Qt_widget_get_conic_2.moc: Qt_widget_get_conic_2.h
#	$(QT_MOC) -o Qt_widget_get_conic_2.moc $<

demo$(EXE_EXT): demo.moc demo$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

demo_elias$(EXE_EXT): demo_elias.moc demo_elias$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_elias demo_elias$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

demo_line$(EXE_EXT): demo.moc demo_line$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_line demo_line$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) $(LDFLAGS)

demo_line_circle$(EXE_EXT): demo_line_circle.moc demo_line_circle$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_line_circle demo_line_circle$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

demo_viewer$(EXE_EXT): demo_viewer.moc demo_viewer$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_viewer demo_viewer$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

#demo_conic$(EXE_EXT): demo_conic.moc demo_conic$(OBJ_EXT) Qt_widget_get_conic_2$(OBJ_EXT)
#	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_conic demo_conic$(OBJ_EXT) Qt_widget_get_conic_2$(OBJ_EXT) $(LDFLAGS)

#Qt_widget_get_conic_2$(OBJ_EXT): Qt_widget_get_conic_2.C Qt_widget_get_conic_2.moc
#	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_get_conic_2.C

Qt_widget_get_arc$(OBJ_EXT): Qt_widget_get_arc.C Qt_widget_get_arc.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_get_arc.C

Qt_widget_get_segment$(OBJ_EXT): Qt_widget_get_segment.C Qt_widget_get_segment.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_get_segment.C

clean: \
                   demo.clean demo_elias.clean demo_line.clean demo_line_circle.clean demo_viewer.clean\
                   Qt_widget_get_arc.clean Qt_widget_get_segment.clean #\
#                   Qt_widget_get_conic_2.clean \
#                   demo_conic.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

