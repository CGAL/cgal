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

MY_DIR = /home/hachenb/CGAL
HOME = /KM/guests/granados
POLY = $(MY_DIR)/Polynomial/include
LOCAL = $(MY_DIR)/local/include
NEF2 = $(MY_DIR)/Nef_2/include
NEF3 = $(MY_DIR)/Nef_3/include 
NEFS2 = $(MY_DIR)/Nef_S2/include
WALKER = $(HOME)/projects/SNC_walker
K3TREE = $(MY_DIR)/Nef_3/k3tree/include
BBOX = $(MY_DIR)/Box_intersection_d/include

# CUSTOM_CXXFLAGS = -I$(POLY) -I$(NEF3) -I$(NEF2) -I$(NEFS2) -I$(BBOX) 

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

CUSTOM_LDFLAGS = -lGLU

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                nef_3$(EXE_EXT) 

moc_Qt_widget_Nef_3.C:  Qt_widget_Nef_3.h
	$(QT_MOC) -o moc_Qt_widget_Nef_3.C Qt_widget_Nef_3.h

moc_Qt_widget_Nef_3$(OBJ_EXT): moc_Qt_widget_Nef_3.C

Qt_widget_Nef_3$(OBJ_EXT): Qt_widget_Nef_3.C

nef_3.clean:
	rm -f nef_3.o moc_Qt_widget_Nef_3.C moc_Qt_widget_Nef_3.o Qt_widget_Nef_3.o

nef_3$(EXE_EXT): nef_3$(OBJ_EXT) Qt_widget_Nef_3$(OBJ_EXT) moc_Qt_widget_Nef_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nef_3 nef_3$(OBJ_EXT) Qt_widget_Nef_3$(OBJ_EXT) moc_Qt_widget_Nef_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   nef_3.clean  

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

