all:		Min_sphere_of_spheres_d.ps

Min_sphere_of_spheres_d.ps:	Min_sphere_of_spheres_d.dvi
		dvips -t a4 $< -o $@

Min_sphere_of_spheres_d.dvi:	Min_sphere_of_spheres_d.tex
		latexmk $<

Min_sphere_of_spheres_d.tex:	Min_sphere_of_spheres_d.fw
		fw +t $<

Min_sphere_of_spheres_d.pdf:	Min_sphere_of_spheres_d.dvi
		dvips -t a4 -P pdf $< -o Min_sphere_of_spheres_d.tmp.ps
		ps2pdf Min_sphere_of_spheres_d.tmp.ps $@

clean:		
		rm -f *.aux *.bak *.dvi Min_sphere_of_spheres_d.tex \
		*~ *.log *.ps *.pdf *.toc *.lis *.tmp.ps *.bbl *.blg