.PHONY: dvi ps

ps: alpha.ps

alpha.ps: dvi
	dvips -o Alpha_shapes_2.ps wrapper.dvi

dvi: wrapper.dvi

wrapper.dvi: wrapper.tex
	latex wrapper.tex
	latex wrapper.tex

all: alpha.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
