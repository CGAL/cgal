
Fixed_precision_nt.ps: Fixedprecision.tex wrapper.tex
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex
	dvips -o Fixed_precision_nt.ps wrapper.dvi

all: Fixed_precision_nt.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
