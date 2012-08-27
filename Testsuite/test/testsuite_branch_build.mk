# A GNU makefile which calls run_testsuite_with_cmake over all directories.

dirs:=$(wildcard ../../*/test/*/)
targets:=$(addsuffix pink_elephant,$(dirs))
cleans:=$(addsuffix green_elephant,$(dirs))

all: ${targets}

clean: ${cleans}

%/pink_elephant:
	@+./run_testsuite_with_cmake $*

%/green_elephant:
	@cd $* && $(MAKE) clean
