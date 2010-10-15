include $(ROOT)/make/comdef.mak

TARGET_EPSS = approx_offset.eps
TARGET_EPSS+= sum_triangles.eps
TARGET_EPSS+= tight.eps

TARGET_PDFS = $(TARGET_EPSS:.eps=.pdf)
TARGET_GIFS = $(TARGET_PDFS:.pdf=.gif)
TARGET_GIFS+= onecyc_in.gif
TARGET_GIFS+= onecyc_out.gif
TARGET_GIFS+= ex_offset.gif
TARGET_GIFS+= convex_offset.gif
TARGET_GIFS+= offset_decomp.gif
TARGET_GIFS+= offset_conv.gif

TARGET0 = $(TARGET_GIFS)
TARGET1 = $(TARGET_PDFS)
TARGET2 = $(TARGET_EPSS)

default:: $(TARGET0) $(TARGET1) $(TARGET2)

include $(MAKEDIR)/comrul.mak

%.pstex_t : %.fig
	fig2dev -L pstex_t -p $(basename $@).pstex $< > $@

%.pstex : %.fig
	fig2dev -L pstex $< > $@

%.dvi : %.pstex_t %.pstex
	sed -e s/Minkowski_sum_2\\/fig\\/// $< > tmp.pstex_t
	latex tmp.tex
	mv tmp.dvi $@

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
	dvips -E -P cmz $< -o

%.eps : %.ps
	$(toeps) > $@

%.pdf : %.eps
	epstopdf --hires $<


onecyc_in.gif: onecyc_in.eps
	convert -resize 137x137% -interlace Line -transparent White $< $@

onecyc_out.gif: onecyc_out.eps
	convert -resize 137x137% -interlace Line -transparent White $< $@

ex_offset.gif: ex_offset.eps
	convert -resize 137x137% -interlace Line -transparent White $< $@

convex_offset.gif: convex_offset.eps
	convert -resize 137x137% -interlace Line -transparent White $< $@

offset_decomp.gif: offset_decomp.eps
	convert -resize 137x137% -interlace Line -transparent White $< $@

offset_conv.gif: offset_conv.eps
	convert -resize 137x137% -interlace Line -transparent White $< $@

.SECONDARY: $(TARGET_PDFS:.pdf=.eps)

LDIRT = fig.dvi fig.aux fig.log $(TARGET_PDFS:.pdf=.ps) \
$(TARGET_GIFS:.gif=.eps2) $(TARGET_GIFS:.gif=.ps) $(TARGET_GIFS:.gif=.dvi) \
$(TARGET_GIFS:.gif=.log) $(TARGET_GIFS:.gif=.aux)
