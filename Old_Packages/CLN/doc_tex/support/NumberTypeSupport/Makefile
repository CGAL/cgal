.PHONY: dvi ps

all: ps

ps: ../../../doc_ps/CLN.ps

dvi: wrapper.dvi

../../../doc_ps/CLN.ps: wrapper.dvi
	dvips -o ../../../doc_ps/CLN.ps wrapper.dvi

wrapper.dvi: CLN.tex wrapper.tex cgal.bib
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
