#-----------------------------------------------------------------------#
# This is a makefile that can be used for quickly installing the CGAL
# libraries. Use this method only if you are an experienced user! The
# following steps are required:
#
# 1) Set up the necessary configuration variables.
#
# 2) Enter 'make install' on the command line. This will create an
#    include makefile in the 'make' directory, by making a call to the
#    install_cgal script.
#
# 3) Enter 'make cgal_lib' or 'make cgal_shared' on the command line. This will
#    create the CGAL libraries.
# 
# N.B. This is NOT the recommended way to install CGAL!!!
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
#                         configuration
#-----------------------------------------------------------------------#

# DO NOT PUT SPACES BEFORE OR AFTER

# CGAL directory
CGAL_DIR=/users/jannes/CGAL-1.0           # enter full pathname

# CONFIGURATION
CGAL_CC=/bin/CC                           # full pathname to compiler

# COMPANION LIBRARIES

#-------------#
#    STL      #
#-------------#
# STL include directory (required if the compiler has no built-in STL!).
STL_DIR=/users/jannes/stl

#-------------#
#    LEDA     #
#-------------#
# Remove the leda entries from the install command if you don't use LEDA.
LEDA_DIR=/packages/LEDA-3.5

# For standard distributions of LEDA you should not have to change
# these variables.
LEDA_INCL_DIR=$(LEDA_DIR)/incl
LEDA_LIB_DIR=$(LEDA_DIR)/

#-------------#
#    GMP      #
#-------------#
# Remove the gmp entries from the install command if you don't use GMP.
GMP_DIR=/users/jannes/gmp-2.0.2
GMP_INCL_DIR=$(GMP_DIR)
GMP_LIB_DIR=$(GMP_DIR)

#-----------------------------------------------------------------------#
#                         targets
#-----------------------------------------------------------------------#

all: install cgal_lib cgal_sharedlib

install:
	./install_cgal \
	     -leda --LEDA_INCL_DIR $(LEDA_INCL_DIR) \
	           --LEDA_LIB_DIR $(LEDA_LIB_DIR) \
	     -gmp  --GMP_INCL_DIR $(GMP_INCL_DIR) \
 	           --GMP_LIB_DIR $(GMP_LIB_DIR) \
	           --STL_DIR $(STL_DIR) \
	     -ni $(CGAL_CC)

cgal_lib:
	cd src ; make -f makefile_lib

cgal_sharedlib:
	cd src ; make -f makefile_sharedlib

