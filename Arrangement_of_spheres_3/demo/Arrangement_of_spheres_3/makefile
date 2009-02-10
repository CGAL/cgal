all: sc sw ipa s2

makefiles: 
	qmake -o makefile.ipa interactive_point_insertion.pro
	qmake -o makefile.sc show_slice.pro
	qmake -o makefile.sw sweep.pro
	qmake -o makefile.s2 show_2d.pro



sc: makefiles
	$(MAKE) -f makefile.sc all

sw: makefiles
	$(MAKE) -f makefile.sw all

ipa: makefiles
	$(MAKE) -f makefile.ipa all

s2: makefiles
	$(MAKE) -f makefile.s2 all

clean: makefiles
	$(MAKE) -f makefile.sc clean
	$(MAKE) -f makefile.s2 clean
	$(MAKE) -f makefile.sw clean
	$(MAKE) -f makefile.ipa clean
