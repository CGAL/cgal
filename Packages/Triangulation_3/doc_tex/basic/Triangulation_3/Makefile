.PHONY: dvi ps

ps: Triangulation3.ps

Triangulation3.ps: wrapper.dvi
	dvips -o Triangulation3.ps wrapper.dvi

dvi: wrapper.dvi

wrapper.dvi: *.tex *.eps
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

all: ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl *.toc
