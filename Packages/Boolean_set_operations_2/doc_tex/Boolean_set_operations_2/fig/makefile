include $(ROOT)/include/make/comdef.mak

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
TARGET_GIFS+= fig_non_strictly_simple.gif
TARGET_GIFS+= fig_circles_rects.gif
TARGET_GIFS+= fig_conic_arcs.gif
TARGET_GIFS+= fig_disks.gif
TARGET_GIFS+= fig_general_polygon.gif
TARGET_GIFS+= fig_general_polygon_with_holes.gif
TARGET_GIFS+= fig_triangles.gif
TARGET_GIFS+= fig_unique.gif

TARGET_PDFS+= $(GIFS:.gif=.pdf)
TARGET_PDFS+= non_strictly_simple.pdf
TARGET_PDFS+= circles_rects.pdf
TARGET_PDFS+= conic_arcs.pdf
TARGET_PDFS+= disks.pdf
TARGET_PDFS+= general_polygon.pdf
TARGET_PDFS+= general_polygon_with_holes.pdf
TARGET_PDFS+= triangles.pdf
TARGET_PDFS+= unique.pdf

TARGET_EPSS = $(TARGET_PDFS:.pdf=.eps)

TARGET0 = $(TARGET_GIFS)
TARGET1 = $(TARGET_PDFS)
TARGET2 = $(TARGET_EPSS)

default:: $(TARGET0) $(TARGET1) $(TARGET2)

include $(MAKEINCDIR)/comrul.mak

define toeps
sed -e "s/%%BoundingBox.*$$/`gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $(basename $@).ps 2>&1 | tr "\n" "\^"`/" $(basename $@).ps | tr "\^" "\n"
endef

%.dvi : %.tex
	cp $< tmp.tex
	latex fig.tex
	mv fig.dvi $@

%.bbox : %.ps
	gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $<
%.gif : %.dvi
	dvips -x 2000 -y 2000 -D4000 -P cmz $< -o
	$(toeps) > $(basename $@).eps2 && rm $(basename $@).ps
	convert -interlace Line -transparent White $(basename $@).eps2 $@;

%.ps : %.dvi
	dvips -P cmz $< -o

%.eps : %.ps
	$(toeps) > $(basename $@).eps

%.pdf : %.eps
	epstopdf --hires $<

.SECONDARY: $(TARGET_PDFS:.pdf=.eps)

fig_non_strictly_simple.dvi: fig_non_strictly_simple.tex
	latex $<

fig_circles_rects.dvi: fig_circles_rects.tex
	latex $<

fig_conic_arcs.dvi: fig_conic_arcs.tex
	latex $<

fig_disks.dvi: fig_disks.tex
	latex $<

fig_general_polygon.dvi: fig_general_polygon.tex
	latex $<

fig_general_polygon_with_holes.dvi: fig_general_polygon_with_holes.tex
	latex $<

fig_triangles.dvi: fig_triangles.tex
	latex $<

fig_unique.dvi: fig_unique.tex
	latex $<

LDIRT = fig.dvi fig.aux fig.log $(TARGET_PDFS:.pdf=.ps) \
$(TARGET_GIFS:.gif=.eps2) $(TARGET_GIFS:.gif=.ps) $(TARGET_GIFS:.gif=.dvi) \
$(TARGET_GIFS:.gif=.log) $(TARGET_GIFS:.gif=.aux)
