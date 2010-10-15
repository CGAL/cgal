include $(ROOT)/make/comdef.mak

TARGET_GIFS = compare_over_curve.gif
TARGET_GIFS+= compare_over_point.gif
TARGET_GIFS+= ex_triangles.gif
TARGET_GIFS+= ex_tri_le.gif
TARGET_GIFS+= ex_tri_ue.gif

TARGET_PDFS = $(TARGET_GIFS:.gif=.pdf)

TARGET0 = $(TARGET_GIFS)
TARGET1 = $(TARGET_PDFS)

include $(MAKEDIR)/comrul.mak

%.pstex_t : %.fig
	fig2dev -L pstex_t -p $(basename $@).pstex $< > $@

%.pstex : %.fig
	fig2dev -L pstex $< > $@

%.dvi : %.pstex_t %.pstex
	sed -e s/Envelope_3\\/fig\\/// $< > tmp.pstex_t
	latex tmp.tex
	mv tmp.dvi $@

define toeps
sed -e "s/%%BoundingBox.*$$/`gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $(basename $@).ps 2>&1 | tr "\n" "\^"`/" $(basename $@).ps | tr "\^" "\n"
endef

%.bbox : %.ps
	gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $<

%.gif : %.dvi
	dvips -x 1000 -y 1000 -D4000 -P cmz $< -o
	$(toeps) > $(basename $@).eps2
	convert -interlace Line -transparent White $(basename $@).eps2 $@;

%.ps : %.dvi
	dvips -E -P cmz $< -o

%.eps : %.ps
	$(toeps) > $(basename $@).eps

%.pdf : %.eps
	epstopdf --hires $<

insert.gif: insert.tex \
            insert_in_face.pstex_t \
	    insert_from_vertex.pstex_t \
	    insert_at_vertices.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// insert_in_face.pstex_t > insert1.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// insert_from_vertex.pstex_t > insert2.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// insert_at_vertices.pstex_t > insert3.pstex_t
	latex $<
	dvips -x 1800 -y 1800 -D4000 -P cmz $(basename $@).dvi -o
	eps2eps $(basename $@).ps $(basename $@).eps
	convert $(basename $@).eps $@

ex_3.gif: ex_3.tex \
            ex_3a.pstex_t ex_3b.pstex_t ex_3c.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// ex_3a.pstex_t > ex_3_1.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// ex_3b.pstex_t > ex_3_2.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// ex_3c.pstex_t > ex_3_3.pstex_t
	latex $<
	dvips -x 1800 -y 1800 -D4000 -P cmz $(basename $@).dvi -o
	eps2eps $(basename $@).ps $(basename $@).eps
	convert $(basename $@).eps $@

.PHONY :: gif pdf
gif :
	$(MAKE) $(TARGET_GIFS)

pdf :
	$(MAKE) $(TARGET_PDFS)

.SECONDARY: $(TARGET_PDFS:.pdf=.eps)

LDIRT = tmp.pstex_t tmp.dvi tmp.aux tmp.log $(TARGET_PDFS:.pdf=.eps2) $(TARGET_PDFS:.pdf=.ps)

# ex_triangles.eps.gif: fan_grids.eps
#	convert -resize 140x140% -interlace Line -transparent White $< $@
 
