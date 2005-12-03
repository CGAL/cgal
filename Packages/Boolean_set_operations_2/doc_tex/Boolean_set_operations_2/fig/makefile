include $(ROOT)/include/make/comdef.mak

TARGET0 = circles_rects.gif
TARGET0+= conic_arcs.gif
TARGET0+= disks.gif
TARGET0+= general_polygon.gif
TARGET0+= general_polygon_with_holes.gif
TARGET0+= non_strictly_simple.gif
TARGET0+= simple.gif
TARGET0+= teaser.gif
TARGET0+= triangles.gif
TARGET0+= unique.gif

include $(MAKEINCDIR)/comrul.mak

%.dvi : %.tex
	cp $< tmp.tex
	latex fig.tex
	mv fig.dvi $@

%.gif : %.dvi
	dvips -x 2000 -y 2000 -D4000 -P cmz $< -o
	eps2eps $(basename $@).ps $(basename $@).eps; rm $(basename $@).ps
	convert $(basename $@).eps $@; rm $(basename $@).eps

%.eps : %.dvi
	dvips -P cmz $< -o
	eps2eps $(basename $@).ps $@

%.pdf : %.eps
	epstopdf --hires $<

.SECONDARY: $(TARGET1:.pdf=.eps)

LDIRT = fig.dvi fig.aux fig.log $(TARGET1:.pdf=.eps)
