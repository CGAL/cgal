.PHONY: dvi ps

ps: Interval_arithmetic.ps

Interval_arithmetic.ps: dvi
	dvips -o Interval_arithmetic.ps wrapper.dvi

dvi: wrapper.dvi

wrapper.dvi: Interval.tex wrapper.tex
	latex wrapper.tex
	bibtex wrapper
	latex wrapper.tex
	latex wrapper.tex

all: Interval_arithmetic.ps

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl
