
all: to_ps to_html

to_html:
	cc_manual_to_html -o ../doc_html manual.tex

to_ps:
	latex manual.tex
	bibtex manual
	makeindex manual
	index_fix manual.ind
	latex manual.tex
	dvips -o ../doc_ps/manual.ps manual.dvi
	gzip -f ../doc_ps/manual.ps

clean:
	\rm -f *.aux *.log *.dvi *.idx *.ilg *.ind *.ind.unfixed *.bbl *.blg
