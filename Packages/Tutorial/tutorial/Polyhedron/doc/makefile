tex:
	latex tutorial
	latex tutorial
	ls -l *.dvi
ps:
	dvips -Ppdf -G0 -t letter -z -o tutorial.ps tutorial.dvi
	ls -l *.ps

pdf:
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true -dAntiAliasColorImages=true  -dColorImageResolution=1200 -dEncodeColorImages#true -dAutoFilterColorImages#false -dColorImageFilter#/FlateEncode tutorial.ps tutorial.pdf

all:
	make tex
	make ps
	make pdf
	make clean
	\rm *.ps
	clear
	ls -l

clean:
	\rm *.aux *.dvi *.log

superclean:
	\rm *.ps *.pdf *~
	make clean



