.PHONY: dvi ps

ps: BooleanOperations.ps

BooleanOperations.ps: dvi
	dvips -o BooleanOperations.ps wrapper.dvi

dvi: wrapper.dvi

wrapper.dvi: main.tex wrapper.tex
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

all: BooleanOperations.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
