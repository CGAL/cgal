.PHONY: dvi ps

all: ps

ps: CLN.ps

dvi: wrapper.dvi

CLN.ps: dvi
	dvips -o CLN.ps wrapper.dvi

wrapper.dvi: CLN.tex wrapper.tex cgal.bib
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
