.PHONY: dvi ps

all: ps
	cp CLN.ps ../../../doc_ps

ps: CLN.ps

dvi: wrapper.dvi

CLN.ps: wrapper.dvi
	dvips -o CLN.ps wrapper.dvi

wrapper.dvi: CLN.tex wrapper.tex cgal.bib
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
