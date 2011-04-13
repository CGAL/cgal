.PHONY: dvi ps

ps: Modifier.ps

Modifier.ps: dvi
	dvips -o Modifier.ps wrapper.dvi

dvi: wrapper.dvi

wrapper.dvi: main.tex modifier.tex wrapper.tex
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

all: Modifier.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
