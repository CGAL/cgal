include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = $(CGAL_LIB_CXXFLAGS) -I.

#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

OBJECTS = \
  ImageIO$(OBJ_EXT) \
  analyze$(OBJ_EXT) \
  bmp$(OBJ_EXT) \
  bmpendian$(OBJ_EXT) \
  bmpread$(OBJ_EXT) \
  convert$(OBJ_EXT) \
  gif$(OBJ_EXT) \
  gis$(OBJ_EXT) \
  inr$(OBJ_EXT) \
  iris$(OBJ_EXT) \
  mincio$(OBJ_EXT) \
  pnm$(OBJ_EXT) \
  recbuffer$(OBJ_EXT) \
  recline$(OBJ_EXT) \
  reech4x4$(OBJ_EXT)

LIB_NAME=libImageIO.a

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: $(LIB_NAME)

$(LIB_NAME): $(OBJECTS)
	$(CGAL_LIB_CREATE)$(LIB_NAME) $(OBJECTS) $(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(LIB_NAME)

clean_obj: 
	rm -f $(OBJECTS)

clean: clean_obj
	rm -f $(LIB_NAME)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.c$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

