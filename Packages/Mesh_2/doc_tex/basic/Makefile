.PHONY: manual.ps
all: to_ps to_html

to_html:
	cc_manual_to_html -o /net/servers/www-sop/prisme/personnel/Laurent.Rineau/Mesh_2 manual.tex

manual.ps:
	latex manual.tex
	makeindex manual
	index_fix manual.ind
	latex manual.tex
	dvips -o manual.ps manual.dvi

manual.ps.gz: manual.ps
	gzip -f manual.ps

manual.pdf: manual.ps
	ps2pdf manual.ps manual.pdf

clean:
	\rm -f *.aux *.log *.dvi *.idx *.ilg *.ind *.ind.unfixed 
