wrapper.ps:	wrapper.dvi
		dvips -t a4 wrapper.dvi -o wrapper.ps

wrapper.dvi:	wrapper.tex \
		Optimisation_ref/main_Min_sphere_of_spheres_d.tex \
		Optimisation_ref/Min_sphere_of_spheres_d.tex \
		Optimisation_ref/Traits.tex \
		Optimisation_ref/Traits_model_d.tex \
		Optimisation_ref/Traits_model_2.tex \
		Optimisation_ref/Traits_model_3.tex
		latex wrapper.tex

clean:
		rm -f *.aux *.dvi *.toc *.idx *.ps *.log *~