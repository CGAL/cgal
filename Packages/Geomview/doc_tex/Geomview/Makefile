
all: Geomview.ps.gz

Geomview.ps.gz: wrapper.ps
	gzip -c wrapper.ps > Geomview.ps.gz

wrapper.ps : wrapper.dvi
	dvips wrapper.dvi -o

wrapper.dvi : main.tex wrapper.tex
	latex wrapper
	latex wrapper
	bibtex wrapper
	bibtex wrapper
	latex wrapper
	latex wrapper

clean:
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl *.toc *.idx Geomview.ps.gz
