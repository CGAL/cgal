# $Id$
#-----------------------------------------------------------

.PHONY:	tex

doc:
	latex wrapper ; latex wrapper
	dvips -o wrapper.ps wrapper
	cp iostream.tex Color.tex ../doc_tex/general
	gzip wrapper.ps
	mv wrapper.ps.gz ../doc_ps/iostream.ps.gz
	cc_manual_to_html -title iostreams -release 1.7 \
		-o ../doc_html/general/iostream/ \
		main.tex

.PHONY:	print

clean:
	- rm -rf *log *aux *dvi *toc *loc *bbl *blg
