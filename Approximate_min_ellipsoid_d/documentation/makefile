mel.ps:		makefile mel.dvi
		dvips -t a4 mel.dvi -o mel.ps

mel.pdf:	makefile mel.dvi
		dvips -P pdf -t a4 mel.dvi -o mel.ps
		ps2pdf mel.ps mel.pdf

mel.dvi:	makefile mel.tex
		latexmk mel.tex

clean:		
		rm -fr *~ tmp* *.dvi *.aux *.aux.bak *.toc *.ps *.gz \
		*.blg *.bbl *.pdf *.log
