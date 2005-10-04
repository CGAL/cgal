%.dvi: %.tex
	../release_scripts/latexloop -bib -idx $*

%.ps: %.dvi
	dvips -D 600 -t a4 $*.dvi -o $*.ps

%.ps.gz: %.ps
	rm -f $*.ps.gz
	gzip $*.ps


all: cc_manual.ps.gz cc_manual_index.ps.gz cc_extract.ps.gz latex_to_html.ps.gz tidyup

tidyup:
	@rm -f *.log *.aux *.dvi *.bbl *.blg *.toc *.idx *.ilg *.ind

clean: tidyup
	@rm -f *.ps.gz *.ps
