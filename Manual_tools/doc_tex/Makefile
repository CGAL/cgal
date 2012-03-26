# WARNING: Use GNU/make or that will fail.
export TEXINPUTS := ${PWD}/../sty:${TEXINPUTS}

%.dvi: %.tex
	@latex $*
	@bibtex $* || echo "no bibtex needed"
	@makeindex $*
	@latex $*
	@latex $*

%.ps: %.dvi
	dvips -D 600 -t a4 $*.dvi -o $*.ps

%.ps.gz: %.ps
	rm -f $*.ps.gz
	gzip $*.ps

%.pdf:	%.tex
	@pdflatex $*
	@bibtex $* || echo "no bibtex needed"
	@makeindex $*
	@pdflatex $*
	@pdflatex $*

all: allpdf

allps: cc_manual.ps.gz cc_manual_index.ps.gz latex_to_html.ps.gz tidyup

allpdf: cc_manual.pdf cc_manual_index.pdf latex_to_html.pdf tidyup

tidyup:
	@rm -f *.log *.aux *.dvi *.bbl *.blg *.toc *.idx *.ilg *.ind

clean: tidyup
	@rm -f *.ps.gz *.ps *.pdf
