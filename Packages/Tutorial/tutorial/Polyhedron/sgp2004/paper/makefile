tex: *.tex
	latex paper
	bibtex paper
	latex paper
	latex paper
ps: 
	dvips -Ppdf -G0 -t letter -z -o paper.ps paper.dvi

pdf:
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true -dAntiAliasColorImages=true  -dColorImageResolution=1200 -dEncodeColorImages#true -dAutoFilterColorImages#false -dColorImageFilter#/FlateEncode paper.ps paper.pdf

all:
	make tex
	make ps
	make pdf

clean:
	\rm -f *.aux *.dvi *.log *.out *.blg *.bbl

superclean:
	\rm -f *.ps *.pdf *~
	make clean



