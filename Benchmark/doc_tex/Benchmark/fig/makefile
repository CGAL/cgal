include $(ROOT)/include/make/comdef.mak

TARGET_GIFS = classification.gif

TARGET_PDFS = $(TARGET_GIFS:.gif=.pdf)

TARGET0 = $(TARGET_GIFS)
TARGET1 = $(TARGET_PDFS)

include $(MAKEINCDIR)/comrul.mak

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
	dvips -P cmz $< -o

%.eps : %.ps
	$(toeps) > $(basename $@).eps

%.pdf : %.eps
	epstopdf --hires $<

.PHONY :: gif pdf
gif :
	$(MAKE) $(TARGET_GIFS)

pdf :
	$(MAKE) $(TARGET_PDFS)

.SECONDARY: $(TARGET_PDFS:.pdf=.eps)

LDIRT = $(TARGET_PDFS:.pdf=.eps2) $(TARGET_PDFS:.pdf=.ps)
