#CC      = g++ -Wall -O -DGEOWIN_USE_NAMESPACE -DGEOWIN_NAMESPACE_NAME=CGAL -DLEDA_STD_IO_HEADERS
CC      = g++ -Wall -O

#use -fPIC when working with non - optimizing compilation
CCLIB   = g++ -G
LD      = g++

RAN     = @touch

LIB     = $(DLIB)

SLIB    = libGeoW.a
DLIB    = libGeoW.so

CFLAGS  = 

DFLAGS  =  -I. -I/LEDA/BETA/incl  

o       = .o 

all: lib

.c.o:
	$(CC) -DLEDA_INSTALL  $(CFLAGS) $(DFLAGS) -c $*.c

.s.o:
	$(CC) -c $*.s
 
lib: $(LIB)
 
include Make.lst
 
$(SLIB): $(OBJECTS)
	ar cr $(SLIB) $?
	$(RAN) $(SLIB) > /dev/null

$(DLIB): $(SLIB)
	@rm -f -r libtmp; mkdir libtmp
	@cd libtmp; ar x ../$(SLIB); rm -f _sparc_add.o _sparc_sub.o
	$(CCLIB) -o $(DLIB) libtmp/*.o
	@rm -f -r libtmp
	cp lib*.* ../../src
 
clean:
	rm -f -r *.o *~ Templates.DB





