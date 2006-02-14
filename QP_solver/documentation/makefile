UpperBounding.pdf:	UpperBounding.ps
			ps2pdf UpperBounding.ps UpperBounding.pdf

UpperBounding.ps:	UpperBounding.dvi
			dvips -P pdf UpperBounding.dvi -o UpperBounding.ps

UpperBounding.dvi:	UpperBounding.tex
			latexmk UpperBounding

Test_suite.pdf:		Test_suite.ps
			ps2pdf Test_suite.ps Test_suite.pdf

Test_suite.ps:		Test_suite.dvi
			dvips -P pdf Test_suite.dvi -o Test_suite.ps

Test_suite.dvi:		Test_suite.tex
			latexmk Test_suite
