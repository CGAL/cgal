# latex ws.tex
# dvips -E ws.dvi -o
# epstool --bbox --copy ws.ps ws.eps
# ps2pdf -dEPSCrop ws.eps ws.pdf

TARGETS = simple_arr.pdf

.DEFAULT_GOAL = all
.PHONY: all clean

all: $(TARGETS)

%.dvi: %.tex
	latex $<

%.ps: %.dvi
	dvips -E $< -o

%.eps: %.ps
	epstool --bbox --copy $< $@

%.pdf: %.eps
	ps2pdf -dEPSCrop $< $@

clean:
	rm -f $(DOCUMENTS)
