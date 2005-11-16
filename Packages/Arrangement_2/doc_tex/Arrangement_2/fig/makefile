include $(ROOT)/include/make/comdef.mak

TARGET0 = arr_segs.gif
# TARGET0+= insert_in_face.gif
# TARGET0+= insert_from_vertex.gif
# TARGET0+= insert_at_vertices.gif
# TARGET0+= insert.gif
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

include $(MAKEINCDIR)/comrul.mak

%.gif : %.pstex_t
	sed -e s/Arrangement_2\\/fig\\/// $< > tmp.pstex_t
	latex tmp.tex
	dvips -x 2000 -y 2000 -D4000 -P cmz tmp.dvi -o
	eps2eps tmp.ps tmp.eps
	convert tmp.eps $@

#	epstopdf --hires tmp.eps

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

LDIRT = tmp.pstex_t tmp.dvi tmp.aux tmp.log tmp.eps tmp.ps
