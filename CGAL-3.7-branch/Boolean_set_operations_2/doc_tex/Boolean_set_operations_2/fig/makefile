include $(ROOT)/make/comdef.mak

GIFS = simple.gif
GIFS+= teaser.gif
# GIFS+= circles_rects.gif
# GIFS+= conic_arcs.gif
# GIFS+= disks.gif
# GIFS+= general_polygon.gif
# GIFS+= general_polygon_with_holes.gif
# GIFS+= non_strictly_simple.gif
# GIFS+= triangles.gif
# GIFS+= unique.gif

TARGET_GIFS = $(GIFS)
TARGET_GIFS+= non_strictly_simple.gif
TARGET_GIFS+= circles_rects.gif
TARGET_GIFS+= conic_arcs.gif
TARGET_GIFS+= disks.gif
TARGET_GIFS+= general_polygon.gif
TARGET_GIFS+= general_polygon_with_holes.gif
TARGET_GIFS+= triangles.gif
TARGET_GIFS+= unique.gif
TARGET_GIFS+= simpleDefsExamples.gif

TARGET_PDFS+= $(GIFS:.gif=.pdf)
TARGET_PDFS+= non_strictly_simple.pdf
TARGET_PDFS+= circles_rects.pdf
TARGET_PDFS+= conic_arcs.pdf
TARGET_PDFS+= disks.pdf
TARGET_PDFS+= general_polygon.pdf
TARGET_PDFS+= general_polygon_with_holes.pdf
TARGET_PDFS+= triangles.pdf
TARGET_PDFS+= unique.pdf
TARGET_PDFS+= simpleDefsExamples.pdf

TARGET_EPSS = $(TARGET_PDFS:.pdf=.eps)

TARGET0 = $(TARGET_GIFS)
TARGET1 = $(TARGET_PDFS)
TARGET2 = $(TARGET_EPSS)

default:: $(TARGET0) $(TARGET1) $(TARGET2)

include $(MAKEDIR)/comrul.mak

define toeps
sed -e "s/%%BoundingBox.*$$/`gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $(basename $@).ps 2>&1 | tr "\n" "\^"`/" $(basename $@).ps | tr "\^" "\n"
endef

%.dvi : %.tex
	cp $< tmp.tex
	latex fig.tex
	mv fig.dvi $@

%.bbox : %.ps
	gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $<
%.bbox2 : %.dvi
	dvips -x 2000 -y 2000 -D4000 -P cmz $< -o
	gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $(basename $<).ps

%.gif : %.dvi
	dvips  -T 17in,11in -x 2000 -y 2000 -D4000 -P cmz $< -o
	$(toeps) > $(basename $@).eps2 && rm $(basename $@).ps
	convert -interlace Line -transparent White $(basename $@).eps2 $@;

%.ps : %.dvi
	dvips -P cmz $< -o

%.eps : %.ps
	$(toeps) > $(basename $@).eps

%.pdf : %.eps
	epstopdf --hires $<

.SECONDARY: $(TARGET_PDFS:.pdf=.eps)

# non_strictly_simple.dvi: non_strictly_simple.tex
#  	latex $<

# circles_rects.dvi: circles_rects.tex
# 	latex $<

# conic_arcs.dvi: conic_arcs.tex
# 	latex $<

# disks.dvi: disks.tex
# 	latex $<

# general_polygon.dvi: general_polygon.tex
# 	latex $<

# general_polygon_with_holes.dvi: general_polygon_with_holes.tex
# 	latex $<

# triangles.dvi: triangles.tex
# 	latex $<

# unique.dvi: unique.tex
# 	latex $<

LDIRT = fig.dvi fig.aux fig.log $(TARGET_PDFS:.pdf=.ps) \
$(TARGET_GIFS:.gif=.eps2) $(TARGET_GIFS:.gif=.ps) $(TARGET_GIFS:.gif=.dvi) \
$(TARGET_GIFS:.gif=.log) $(TARGET_GIFS:.gif=.aux)
