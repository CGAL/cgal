tex: *.tex
	latex tutorial
	bibtex tutorial
	latex tutorial
	latex tutorial
ps: 
	dvips -Ppdf -G0 -t letter -z -o tutorial.ps tutorial.dvi

pdf:
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true -dAntiAliasColorImages=true  -dColorImageResolution=1200 -dEncodeColorImages#true -dAutoFilterColorImages#false -dColorImageFilter#/FlateEncode tutorial.ps tutorial.pdf

all:
	make tex
	make ps
	make pdf

clean:
	\rm -f *.aux *.dvi *.log *.out *.blg *.bbl

superclean:
	\rm -f *.ps *.pdf *~
	make clean

abstract: *.tex
	latex abstract
	latex abstract
	latex abstract
	dvips -Ppdf -G0 -t letter -z -o abstract.ps abstract.dvi
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true -dAntiAliasColorImages=true  -dColorImageResolution=1200 -dEncodeColorImages#true -dAutoFilterColorImages#false -dColorImageFilter#/FlateEncode abstract.ps abstract.pdf
