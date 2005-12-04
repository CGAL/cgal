include $(ROOT)/include/make/comdef.mak

TARGET0 =
# TARGET0+= circles_rects.gif
# TARGET0+= conic_arcs.gif
# TARGET0+= disks.gif
# TARGET0+= general_polygon.gif
# TARGET0+= general_polygon_with_holes.gif
# TARGET0+= non_strictly_simple.gif
TARGET0+= simple.gif
TARGET0+= teaser.gif
# TARGET0+= triangles.gif
# TARGET0+= unique.gif

TARGET1 = fig_non_strictly_simple.gif
TARGET1+= fig_circles_rects.gif
TARGET1+= fig_conic_arcs.gif
TARGET1+= fig_disks.gif
TARGET1+= fig_general_polygon.gif
TARGET1+= fig_general_polygon_with_holes.gif
TARGET1+= fig_triangles.gif
TARGET1+= fig_unique.gif

include $(MAKEINCDIR)/comrul.mak

%.dvi : %.tex
	cp $< tmp.tex
	latex fig.tex
	mv fig.dvi $@

%.gif : %.dvi
	dvips -x 2000 -y 2000 -D4000 -P cmz $< -o
	eps2eps $(basename $@).ps $(basename $@).eps;
	convert $(basename $@).eps $@;

%.eps : %.dvi
	dvips -P cmz $< -o
	eps2eps $(basename $@).ps $@

%.pdf : %.eps
	epstopdf --hires $<

.SECONDARY: $(TARGET1:.gif=.eps)

fig_non_strictly_simple.dvi: fig_non_strictly_simple.tex
	latex $<

fig_circles_rects.dvi: fig_circles_rects.tex
	latex $<

fig_conic_arcs.dvi: fig_conic_arcs.tex
	latex $<

fig_disks.dvi: fig_disks.tex
	latex $<

fig_general_polygon.vi: fig_general_polygon.tex
	latex $<

fig_general_polygon_with_holes.dvi: fig_general_polygon_with_holes.tex
	latex $<

fig_triangles.dvi: fig_triangles.tex
	latex $<

fig_unique.dvi: fig_unique.tex
	latex $<

LDIRT = fig.dvi fig.aux fig.log $(TARGET1:.pdf=.eps)
