.PHONY: dvi ps

ps: alpha.ps

alpha.ps: dvi
	dvips -o Alpha_shapes_2.ps wrapper.dvi

dvi: wrapper.dvi

wrapper.dvi: alpha2.tex
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

all: alpha.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl *.idx *.ind
