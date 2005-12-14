include $(ROOT)/include/make/comdef.mak

TARGET0 = arr_segs.gif
TARGET0+= ex_1.gif
TARGET0+= ex_2.gif
TARGET0+= ex_4.gif
TARGET0+= ex_5.gif
TARGET0+= ex_8.gif
TARGET0+= ex_10.gif
TARGET0+= ex_12.gif
TARGET0+= ex_13.gif
TARGET0+= ex_14.gif
TARGET0+= ex_16.gif
TARGET0+= ex_17.gif
TARGET0+= ex_18.gif
TARGET0+= ex_19.gif
TARGET0+= ex_20.gif
TARGET0+= ex_22.gif
TARGET0+= ex_24.gif
TARGET0+= ex_25.gif

# wrapfigure:
TARGET0+= triangle.gif
TARGET0+= connect_comp.gif
TARGET0+= pred_around_vertex.gif
TARGET0+= h_shape.gif

TARGET1 = $(TARGET0:.gif=.pdf)
TARGET1+= ex_3a.pdf
TARGET1+= ex_3b.pdf
TARGET1+= ex_3c.pdf
TARGET1+= insert_in_face.pdf
TARGET1+= insert_from_vertex.pdf
TARGET1+= insert_at_vertices.pdf

include $(MAKEINCDIR)/comrul.mak

%.dvi : %.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// $< > tmp.pstex_t
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
	convert $(basename $@).eps2 $@;

%.ps : %.dvi
	dvips -P cmz $< -o

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

gif :
	$(MAKE) $(TARGET0)

.SECONDARY: $(TARGET1:.pdf=.eps) $(TARGET0:.pdf=.eps)

LDIRT = tmp.pstex_t tmp.dvi tmp.aux tmp.log $(TARGET1:.pdf=.eps) $(TARGET0:.gif=.eps2)
