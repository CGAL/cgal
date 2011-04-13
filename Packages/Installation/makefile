#-----------------------------------------------------------------------#
#                                                                       #
# This is a makefile that can be used for installing the CGAL libraries.# 
# Use this method only if you know what to do!                          #
#                                                                       #
# First-time installers (especially): use the interactive script        #
#                                                                       #
#             ./install_cgal -i                                         #
#                                                                       #
#-----------------------------------------------------------------------#
#                                                                       #
# The following steps are required:                                     #
#                                                                       #
# 1) Set up the necessary configuration variables.                      #
#                                                                       #
# 2) Enter 'make install' on the command line. This will create an      #
#    include makefile in the 'make' directory, by making a call to the  #
#    install_cgal script.                                               #
#                                                                       #
# 3) Enter 'make cgal_lib' or 'make cgal_shared' on the command line.   #
#    This will create the CGAL libraries.                               #
#                                                                       #
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
#                         configuration                                 #
#-----------------------------------------------------------------------#

# DO NOT PUT SPACES BEFORE OR AFTER

# CGAL directory
CGAL_DIR=/packages/CGAL                   # enter full pathname

# COMPILER
CGAL_CC=/bin/CC                           # full pathname to compiler

# COMPANION LIBRARIES

#-----------------------------------------------------------------------#
#                         LEDA                                          #
#-----------------------------------------------------------------------#
# Remove the leda entries from the install command if you don't use LEDA.
LEDA_DIR=/packages/LEDA

# For standard distributions of LEDA you should not have to change
# these variables.
LEDA_INCL_DIR=$(LEDA_DIR)/incl
LEDA_LIB_DIR=$(LEDA_DIR)/

#-----------------------------------------------------------------------#
#                         GMP                                           #
#-----------------------------------------------------------------------#
# Remove the gmp entries from the install command if you don't use GMP.
GMP_DIR=/packages/gmp
GMP_INCL_DIR=$(GMP_DIR)
GMP_LIB_DIR=$(GMP_DIR)

#-----------------------------------------------------------------------#
#                         CLN                                           #
#-----------------------------------------------------------------------#
# Remove the gmp entries from the install command if you don't use CLN.
CLN_DIR=/packages/cln
CLN_INCL_DIR=$(CLN_DIR)
CLN_LIB_DIR=$(CLN_DIR)

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
	     -cln  --CLN_INCL_DIR $(CLN_INCL_DIR) \
	           --CLN_LIB_DIR $(CLN_LIB_DIR) \
	           --STL_DIR $(STL_DIR) \
	     -ni $(CGAL_CC)

cgal_lib:
	cd src ; make -f makefile_lib

cgal_sharedlib:
	cd src ; make -f makefile_sharedlib

