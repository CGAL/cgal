
all: to_ps to_html

to_html:
	cc_manual_to_html -o html manual.tex

to_ps:
	latex manual.tex
	makeindex manual
	index_fix manual.ind
	latex manual.tex
	dvips -o manual.ps manual.dvi
	gzip -f manual.ps

clean:
	\rm -f *.aux *.log *.dvi *.idx *.ilg *.ind *.ind.unfixed 
