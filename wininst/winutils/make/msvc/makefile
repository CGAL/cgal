# This file contains CGAL makefile settings for the following platform:
# OS:		Win32
# COMPILER:	VC++ 12.00.*
# GMP:		supported

#---------------------------------------------------------------------#
#                    make and other utilities to use
#---------------------------------------------------------------------#

MAKE = nmake
RM = del
RM_FORCE = del
MV = move
CP = copy

#---------------------------------------------------------------------#
#                    include/lib directory settings
#---------------------------------------------------------------------#

CGAL_INCL_DIR = $(CGAL_ROOT)\include
CGAL_LIB_DIR = $(CGAL_ROOT)\lib

#---------------------------------------------------------------------#
#                    filename related flags
#---------------------------------------------------------------------#

OBJ_EXT = .obj
EXE_EXT =
LIB_SUFFIX = .lib

#---------------------------------------------------------------------#
#                    some specific compiler flags
#---------------------------------------------------------------------#

OBJ_OPT = -c
DEBUG_OPT = -Z7
START_LD_OPT = -link 
LIBPATH_OPT = -LIBPATH:
LIB_OPT = 

#---------------------------------------------------------------------#
#                    suffixes and rules
#---------------------------------------------------------------------#

.SUFFIXES:  .C .obj .exe .clean

#.C:
#	$(MAKE) $*.exe

.C.clean:
	if exist $*.exe $(RM_FORCE) $*.exe
	if exist $*.obj $(RM_FORCE) $*.obj
	if exist $*.pdb $(RM_FORCE) $*.pdb
	if exist vc60.pdb $(RM_FORCE) vc60.pdb
	if exist $*.ilk $(RM_FORCE) $*.ilk


#---------------------------------------------------------------------#
#                    compiler
#---------------------------------------------------------------------#

CGAL_CC  = cl
CGAL_CXX = cl

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# For more information about the compiler flags see the installation guide.

# compiler specific configuration dir
CGAL_INCL_CONF_DIR = \
	$(CGAL_INCL_DIR)\CGAL\config\msvc

# *** Fill in your GMP include directory ***
GMP_INCL_DIR = \
	$(CGAL_ROOT)\auxiliary\wingmp\gmp-2.0.2

# *** Fill in custom compiler flags ***
CUSTOM_CXXFLAGS = 

# *** to work around the long-name-problem with sun as/ld
LONG_NAME_PROBLEM_CXXFLAGS =

CGAL_CXXFLAGS = \
	$(CUSTOM1_CXXFLAGS) $(CUSTOM_CXXFLAGS) -TP -GR -GX -Zm900 -nologo \
	-I$(CGAL_ROOT)\stlport  \
	-DCGAL_USE_GMP \
	-I$(GMP_INCL_DIR) \
	-I$(CGAL_INCL_CONF_DIR) \
	-I$(CGAL_INCL_DIR) \
	$(LE_CXXFLAGS) \


CGAL_WINLIB_CXXFLAGS = \
	$(CGAL_LIB_CXXFLAGS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

EXE_OPT = -Fp

# *** Fill in your GMP lib directory ***
GMP_LIB_DIR = $(GMP_INCL_DIR)\msvc

# *** Fill in any additional linker flags you need ***
CUSTOM_LIBPATH = 

CUSTOM_LDFLAGS =

# *** to work around the long-name-problem with sun as/ld
LONG_NAME_PROBLEM_LDFLAGS =

CGAL_LIBPATH =

CGAL_LD_LIBPATH = \
	$(CUSTOM_LIBPATH) \
	-LIBPATH:$(CGAL_LIB_DIR)\$(CGAL_OS_COMPILER) $(LE_LIB_DIR) -LIBPATH:$(GMP_LIB_DIR) $(CGALWIN_LIB_DIR_PATH)

CGAL_LDFLAGS = $(CUSTOM1_LDFLAGS) \
	-link $(CGAL_LD_LIBPATH) -nologo \
	$(CUSTOM_LDFLAGS) \
	CGAL.lib  $(LE_LIBS_LIST) ws2_32.lib user32.lib gdi32.lib comdlg32.lib shell32.lib advapi32.lib gmp.lib

CGAL_WINDOW_LIBPATH = 

CGAL_LD_WINDOW_LIBPATH = \
	-LIBPATH:$(CGAL_LIB_DIR)\$(CGAL_OS_COMPILER) $(LE_LIB_DIR) -LIBPATH:$(GMP_LIB_DIR) $(CGALWIN_LIB_DIR_PATH)

CGAL_WINDOW_LDFLAGS = $(CUSTOM1_LDFLAGS) \
	-link $(CGAL_LD_WINDOW_LIBPATH) -nologo \
	$(CUSTOM_LDFLAGS) \
	CGAL.lib $(LE_LIBS_LIST) ws2_32.lib user32.lib gdi32.lib comdlg32.lib shell32.lib advapi32.lib gmp.lib  

CGAL_GEOWIN_LDFLAGS = \
	$(CGAL_WINDOW_LDFLAGS)



#---------------------------------------------------------------------#
#                    commands and flags for creating libraries
#---------------------------------------------------------------------#

CGAL_EXTRA_OBJECTS = workaround_4_ms.obj
CGAL_LIB = CGAL.lib
CGAL_WIN_LIB = CGALWin.lib
CGAL_WIN_EXTRA_DIR = mswin
CGAL_LIB_CXXFLAGS = \
	$(CGAL_CXXFLAGS)
CGAL_LIB_LDFLAGS = 

CGAL_LIB_CREATE = LIB /OUT:
CGAL_OBJ_PREFIX = 
CGAL_LIB_DESTINATION = $(CGAL_ROOT)\lib\$(CGAL_OS_COMPILER)


# *** only used for the testsuite
# *** program used to create stacktrace from coredump
PRINT_STACKTRACE_PROGRAM = 


#---------------------------------------------------------------------#
#                    EOF
#---------------------------------------------------------------------#


