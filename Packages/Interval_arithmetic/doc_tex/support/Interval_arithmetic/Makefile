
Interval_arithmetic.ps: main.tex wrapper.tex
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex
	dvips -o Interval_arithmetic.ps wrapper.dvi

all: Interval_arithmetic.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
