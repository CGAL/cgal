.PHONY: dvi ps

ps: bops.ps

bops.ps: dvi
	dvips -o bops.ps wrapper.dvi

dvi: wrapper.dvi

wrapper.dvi: main.tex wrapper.tex
	latex wrapper.tex
	latex wrapper.tex

all: bops.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
