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

# X_DIR = /opt/xfree86-4.0.1
# X_GLUTDIR = /usr

NEF3_DIR = $(MY_DIR)/Nef_3/include
NEFS2_DIR = $(MY_DIR)/Nef_S2/include
NEF2_DIR = $(MY_DIR)/Nef_2/include
LOCAL_DIR = $(MY_DIR)/local/include

# X_INCLUDE = -I$(X_DIR)/include -I$(X_GLUTDIR)/include

 CXXFLAGS = \
            -fprofile-arcs -ftest-coverage \
            -I$(NEF2_DIR) \
 	    -I$(NEF3_DIR) \
            -I$(NEFS2_DIR) \
 	    $(CGAL_CXXFLAGS) \
            $(LONG_NAME_PROBLEM_CXXFLAGS) \
            $(DEBUG_OPT) \
 	    -Wno-deprecated

#CXXFLAGS = \
#           -I$(NEF2_DIR) \
#	    -I$(NEF3_DIR) \
#           -I$(NEFS2_DIR) \
#	    $(CGAL_CXXFLAGS) \
#           $(LONG_NAME_PROBLEM_CXXFLAGS) \
#           $(DEBUG_OPT) \
#	    -Wno-deprecated

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

OGL_LIBS = \
        -Wl,-R$(X_DIR)/lib \
        -lGL -lGLU -lXext -lXmu -lX11 -lglut

CUSTOM_LDFLAGS = $(OGL_LIBS)

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH) 

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test_Nef_3$(EXE_EXT)

test_Nef_3$(EXE_EXT): test_Nef_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_Nef_3 test_Nef_3$(OBJ_EXT) $(LDFLAGS)

clean: \
		   test_Nef_3.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

