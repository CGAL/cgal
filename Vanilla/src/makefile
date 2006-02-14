# This is the makefile for compiling the library libVanill.a 

#---------------------------------------------------------------------#
#                    the library name
#---------------------------------------------------------------------#
LIB_NAME = libVanill.a

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = 
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
#

CC = $(CGAL_CXX)

CXXFLAGS = \
           -I../include \
           $(CGAL_CXXFLAGS) \
           -g 


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

LIB_OBJECTS = flavor.o

all: $(LIB_OBJECTS)
	rm -f $(LIB_NAME)
	ar rcv $(LIB_NAME) $(LIB_OBJECTS) 
	mv $(LIB_NAME) ../lib

flavor.o: \
	../include/CEP/Vanilla/flavor.h 

clean:
	/bin/rm -f *.o \
        core

really_clean: clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CC) $(CXXFLAGS) -c  $<

