IMAGEIO_OBJECTS = \
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

#---------------------------------------------------------------------#
#                    build rules
#---------------------------------------------------------------------#

# Generated with the perl uniline:
# perl -ne '/([[:alnum:]]+).*OBJ_EXT/ && print "$1\$(OBJ_EXT): $1.c\n      \$(CGAL_CXX) \$(CXXFLAGS) \$(OBJ_OPT) \$(IMAGEIO_PATH)$1.c\n\n"' imageio.mk 

ImageIO$(OBJ_EXT): $(IMAGEIO_PATH)ImageIO.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)ImageIO.c

analyze$(OBJ_EXT): $(IMAGEIO_PATH)analyze.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)analyze.c

bmp$(OBJ_EXT): $(IMAGEIO_PATH)bmp.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)bmp.c

bmpendian$(OBJ_EXT): $(IMAGEIO_PATH)bmpendian.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)bmpendian.c

bmpread$(OBJ_EXT): $(IMAGEIO_PATH)bmpread.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)bmpread.c

convert$(OBJ_EXT): $(IMAGEIO_PATH)convert.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)convert.c

gif$(OBJ_EXT): $(IMAGEIO_PATH)gif.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)gif.c

gis$(OBJ_EXT): $(IMAGEIO_PATH)gis.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)gis.c

inr$(OBJ_EXT): $(IMAGEIO_PATH)inr.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)inr.c

iris$(OBJ_EXT): $(IMAGEIO_PATH)iris.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)iris.c

mincio$(OBJ_EXT): $(IMAGEIO_PATH)mincio.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)mincio.c

pnm$(OBJ_EXT): $(IMAGEIO_PATH)pnm.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)pnm.c

recbuffer$(OBJ_EXT): $(IMAGEIO_PATH)recbuffer.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)recbuffer.c

recline$(OBJ_EXT): $(IMAGEIO_PATH)recline.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)recline.c

reech4x4$(OBJ_EXT): $(IMAGEIO_PATH)reech4x4.c
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $(IMAGEIO_PATH)reech4x4.c


.PHONY: imageio-clean

imageio-clean:
	rm -f $(IMAGEIO_OBJECTS)

