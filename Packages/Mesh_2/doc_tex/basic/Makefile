.PHONY: manual.dvi to_html
all: to_ps to_html

to_html:
	mkdir -p html
	LATEX_CONV_INPUTS=.:./Mesh_2_ref:./Mesh_2:../../examples: cc_manual_to_html -o html manual.tex

manual.dvi:
	TEXINPUTS=${PWD}/Mesh_2_ref:./Mesh_2:../../examples:${TEXINPUTS} latex manual.tex
	makeindex manual
	index_fix manual.ind
	TEXINPUTS=${PWD}/Mesh_2_ref:./Mesh_2:../../examples:${TEXINPUTS} latex manual.tex

manual.ps: manual.dvi
	dvips -o manual.ps manual.dvi

manual.ps.gz: manual.ps
	gzip -f manual.ps

manual.pdf: manual.ps
	ps2pdf manual.ps manual.pdf

clean:
	\rm -f *.aux *.log *.dvi *.idx *.ilg *.ind *.ind.unfixed *.mtc*
