include $(ROOT)/make/comdef.mak

TARGET_GIFS = triangle.gif
TARGET_GIFS+= connect_comp.gif
TARGET_GIFS+= pred_around_vertex.gif
TARGET_GIFS+= h_shape.gif

TARGET_PDFS = $(TARGET_GIFS:.gif=.pdf)

TARGET_PDFS+= arr_segs.pdf
TARGET_PDFS+= ex_1.pdf
TARGET_PDFS+= ex_2.pdf
TARGET_PDFS+= ex_4.pdf
TARGET_PDFS+= ex_5.pdf
TARGET_PDFS+= ex_8.pdf
TARGET_PDFS+= ex_10.pdf
TARGET_PDFS+= ex_12.pdf
TARGET_PDFS+= ex_13.pdf
TARGET_PDFS+= ex_14.pdf
TARGET_PDFS+= ex_16.pdf
TARGET_PDFS+= ex_17.pdf
TARGET_PDFS+= ex_18.pdf
TARGET_PDFS+= ex_19.pdf
TARGET_PDFS+= ex_20.pdf
TARGET_PDFS+= ex_22.pdf
TARGET_PDFS+= ex_24.pdf
TARGET_PDFS+= ex_25.pdf

TARGET_PDFS+= ex_3a.pdf
TARGET_PDFS+= ex_3b.pdf
TARGET_PDFS+= ex_3c.pdf
TARGET_PDFS+= insert_in_face.pdf
TARGET_PDFS+= insert_from_vertex.pdf
TARGET_PDFS+= insert_at_vertices.pdf
TARGET_PDFS+= unb_dcel.pdf

TARGET0 = $(TARGET_GIFS)
TARGET1 = $(TARGET_PDFS)

include $(MAKEDIR)/comrul.mak

%.pstex_t : %.fig
	fig2dev -L pstex_t -p $(basename $@).pstex $< > $@

%.pstex : %.fig
	fig2dev -L pstex $< > $@

%.dvi : %.pstex_t %.pstex
	sed -e s/Arrangement_on_surface_2\\/fig\\/// $< > tmp.pstex_t
	latex tmp.tex
	mv tmp.dvi $@

define toeps
sed -e "s/%%BoundingBox.*$$/`gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $(basename $@).ps 2>&1 | tr "\n" "\^"`/" $(basename $@).ps | tr "\^" "\n"
endef

%.bbox : %.ps
	gs -dSAFER -dCompatibilityLevel=1.2 -q -sDEVICE=bbox -dNOPAUSE -dBATCH $<

%.gif : %.dvi
	dvips -x 2000 -y 2000 -D4000 -P cmz $< -o
	$(toeps) > $(basename $@).eps2
	convert -interlace Line -transparent White $(basename $@).eps2 $@;

%.ps : %.dvi
	dvips -E -P cmz $< -o

%.eps : %.ps
	$(toeps) > $@

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

fan_grids.gif: fan_grids.eps
	convert -resize 137x137% -interlace Line -transparent White $< $@

Europe.gif: Europe.eps
	convert -resize 140x140% -interlace Line -transparent White $< $@
