# This file contains CGAL makefile settings for the following platform:
# OS:		Win32
# COMPILER:	Borland 5.[4,5]
# 
# Note - to run with MSVC make, 
# replace "MAKE = make -N" by "MAKE = nmake"
#
#---------------------------------------------------------------------#
#                    make and other utilities to use
#---------------------------------------------------------------------#

MAKE = make -N
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
#                    os/compiler description
#---------------------------------------------------------------------#

#CGAL_OS_COMPILER = i686_CYGWINNT-4.0-20.1_bcc32.exe-5.4_LEDA

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
DEBUG_OPT = -v
START_LD_OPT = 
LIBPATH_OPT = -L
LIB_OPT = 

#---------------------------------------------------------------------#
#                    suffixes and rules
#---------------------------------------------------------------------#

.SUFFIXES: .C .obj .exe .clean

#.C:
#	$(MAKE) $*.exe

.C.clean:
	if exist $*.exe $(RM_FORCE) $*.exe
	if exist $*.obj $(RM_FORCE) $*.obj
	if exist $*.tds $(RM_FORCE) $*.tds


#---------------------------------------------------------------------#
#                    compiler
#---------------------------------------------------------------------#

CGAL_CC  = bcc32
CGAL_CXX = bcc32

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

# For more information about the compiler flags see the installation guide.

# compiler specific configuration dir
CGAL_INCL_CONF_DIR = \
	$(CGAL_INCL_DIR)\CGAL\config\bcc

# *** Fill in any additional compiler flags you need ***
CUSTOM_CXXFLAGS = 

# *** Fill in your GMP include directory ***
GMP_INCL_DIR = \
	$(CGAL_ROOT)\auxiliary\wingmp\gmp-2.0.2

# *** to work around the long-name-problem with sun as/ld
LONG_NAME_PROBLEM_CXXFLAGS =

CGAL_CXXFLAGS = \
	$(CUSTOM1_CXXFLAGS) $(CUSTOM_CXXFLAGS) -P -vi- -w-inl \
	-DCGAL_USE_GMP \
	-I$(GMP_INCL_DIR) \
	-I$(CGAL_INCL_CONF_DIR) \
	$(LE_CXXFLAGS) \
	-I$(CGAL_INCL_DIR)

CGAL_WINLIB_CXXFLAGS = \
	$(CGAL_LIB_CXXFLAGS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

EXE_OPT = -e

# *** Fill in your GMP lib directory ***
# (e.g. /usr/local/lib)
GMP_LIB_DIR = $(GMP_INCL_DIR)\bcc

# *** Fill in any additional linker flags you need ***
CUSTOM_LIBPATH = 

CUSTOM_LDFLAGS = 

# *** to work around the long-name-problem with sun as/ld
LONG_NAME_PROBLEM_LDFLAGS =

CGAL_LD_LIBPATH =

CGAL_LIBPATH = $(CUSTOM1_LDFLAGS) \
	$(CUSTOM_LIBPATH)\
	-L$(CGAL_LIB_DIR)\$(CGAL_OS_COMPILER) $(LE_LIB_DIR) -L$(GMP_LIB_DIR) $(CGALWIN_LIB_DIR_PATH)

CGAL_LDFLAGS = \
	$(CGAL_LD_LIBPATH)  \
	$(CUSTOM_LDFLAGS) \
	CGAL.lib  $(LE_LIBS_LIST)  gmp.lib  

CGAL_WINDOW_LIBPATH = $(CUSTOM1_LDFLAGS) \
	-L$(CGAL_LIB_DIR)\$(CGAL_OS_COMPILER) $(LE_LIB_DIR) -L$(GMP_LIB_DIR) $(CGALWIN_LIB_DIR_PATH)

CGAL_WINDOW_LDFLAGS = \
	$(CGAL_WINDOW_LIBPATH)  \
	$(CUSTOM_LDFLAGS) \
	CGAL.lib $(LE_LIBS_LIST) gmp.lib 

CGAL_GEOWIN_LDFLAGS = \
	$(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    commands and flags for creating libraries
#---------------------------------------------------------------------#
CGAL_EXTRA_OBJECTS =
CGAL_LIB = CGAL.lib
CGAL_LIB_CXXFLAGS = \
	-P \
	$(CGAL_CXXFLAGS)
CGAL_LIB_LDFLAGS = 

#CGAL_LIB_CREATE = tlib /C /P512 ''
CGAL_LIB_CREATE = tlib /C /P512 
CGAL_OBJ_PREFIX = +
CGAL_LIB_DESTINATION = $(CGAL_ROOT)\lib\$(CGAL_OS_COMPILER)

# *** only used for the testsuite
# *** program used to create stacktrace from coredump
PRINT_STACKTRACE_PROGRAM = 


#---------------------------------------------------------------------#
#                    EOF
#---------------------------------------------------------------------#

