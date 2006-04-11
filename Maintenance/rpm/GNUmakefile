CGAL_URL := $(shell ./get_cgal_trunk_url)

RPM_BUILD_ARCH:=$(shell rpm --showrc | awk '/^build arch/ { print $$4}')

# SUDO variable
# If you are not root, put:
#  SUDO = sudo 
# of 
#  SUDO = su -c 
# or an equivalent command, to became root, in the local file my_makefile
# (which is not submitted to SVN).

SUDO:=

all:
	@printf "Specify:\n  make CGAL\n  make CGAL-manual-tools\n"
	@printf "or:\n  make CGAL.install\n  make CGAL-manual-tools.install\n"



#############################################################################
#
#                       CGAL-manual-tools package
#
#############################################################################
CGAL_MANUAL_REVISION:=$(shell awk '/Version:/ {print $$2}' SPECS/CGAL-manual-tools.spec)
CGAL_MANUAL_RELEASE_NUMBER:=$(shell awk '/Release:/ {print $$2}' SPECS/CGAL-manual-tools.spec)
CGAL_MANUAL_PACKAGE_VERSION=$(CGAL_MANUAL_REVISION)-$(CGAL_MANUAL_RELEASE_NUMBER).$(RPM_BUILD_ARCH)

CGAL_MANUAL_PACKAGE_PATH=RPMS/$(RPM_BUILD_ARCH)/CGAL-manual-tools-$(CGAL_MANUAL_PACKAGE_VERSION).rpm

CGAL_MANUAL_PACKAGE_INSTALLED_VERSION=$(shell rpm -q --qf '%{VERSION}-%{RELEASE}.%{ARCH}' CGAL-manual-tools)

CGAL_MANUAL_PACKAGE_SOURCES = \
        SOURCES/Manual_tools-$(CGAL_MANUAL_REVISION).tar.gz \
        SOURCES/Manual-$(CGAL_MANUAL_REVISION).tar.gz

CGAL_MANUAL_PACKAGE_DEPS = \
        $(CGAL_MANUAL_PACKAGE_SOURCES) \
	SPECS/CGAL-manual-tools.spec

.PHONY:: prepare-CGAL-manual-tools-sources

prepare-CGAL-manual-tools-sources: $(CGAL_MANUAL_PACKAGE_SOURCES)

SOURCES/Manual_tools-$(CGAL_MANUAL_REVISION).tar.gz:
	cd SOURCES && svn export -r $(CGAL_MANUAL_REVISION) --force \
	    $(CGAL_URL)/Manual_tools && \
	tar czf Manual_tools-$(CGAL_MANUAL_REVISION).tar.gz Manual_tools/
	rm -rf SOURCES/Manual_tools/

SOURCES/Manual-$(CGAL_MANUAL_REVISION).tar.gz:
	cd SOURCES && svn export -r $(CGAL_MANUAL_REVISION) --force \
	    $(CGAL_URL)/Manual && \
	tar czf Manual-$(CGAL_MANUAL_REVISION).tar.gz -C Manual/ .
	rm -rf SOURCES/Manual/

$(CGAL_MANUAL_PACKAGE_PATH): $(CGAL_MANUAL_PACKAGE_DEPS)
ifeq ($(CGAL_MANUAL_PACKAGE_VERSION), $(CGAL_MANUAL_PACKAGE_INSTALLED_VERSION))
	@echo CGAL-manual-tools.$(CGAL_MANUAL_PACKAGE_VERSION) already \
	  installed.
	@echo You should increase the Release: tag.
else
	./rpmbuild -bb SPECS/CGAL-manual-tools.spec
endif

.PHONY:: CGAL-manual-tools

CGAL-manual-tools: $(CGAL_MANUAL_PACKAGE_PATH)

CGAL-manual-tools.install: $(CGAL_MANUAL_PACKAGE_PATH)
ifeq ($(findstring not installed,$(CGAL_MANUAL_PACKAGE_INSTALLED_VERSION)),)
	$(SUDO) rpm -Uvh $(CGAL_MANUAL_PACKAGE_PATH)
else
	$(SUDO) rpm -ivh $(CGAL_MANUAL_PACKAGE_PATH)
endif

#############################################################################
#
#                       CGAL package
#
#############################################################################
CGAL_REVISION:=$(shell awk '/Version:/ {print $$2}' SPECS/CGAL.spec)
CGAL_RELEASE_NUMBER:=$(shell awk '/define internal_release/ {print $$3}' SPECS/CGAL.spec)
CGAL_PACKAGE_VERSION=$(CGAL_REVISION)-$(CGAL_RELEASE_NUMBER).$(RPM_BUILD_ARCH)

CGAL_PACKAGE_PATH=RPMS/$(RPM_BUILD_ARCH)/CGAL-$(CGAL_PACKAGE_VERSION).rpm
CGAL_DEVEL_PACKAGE_PATH=RPMS/$(RPM_BUILD_ARCH)/CGAL-devel-$(CGAL_PACKAGE_VERSION).rpm

CGAL_PACKAGE_INSTALLED_VERSION=$(shell rpm -q --qf '%{VERSION}-%{RELEASE}.%{ARCH}' CGAL)

CGAL_PACKAGE_SOURCES = \
        SOURCES/CGAL-$(CGAL_REVISION)-I-$(CGAL_RELEASE_NUMBER).tar.gz

CGAL_PACKAGE_DEPS = \
        $(CGAL_PACKAGE_SOURCES) \
	SPECS/CGAL.spec

$(CGAL_PACKAGE_PATH): $(CGAL_PACKAGE_DEPS)
ifeq ($(CGAL_PACKAGE_VERSION), $(CGAL_PACKAGE_INSTALLED_VERSION))
	@echo CGAL-manual-tools.$(CGAL_PACKAGE_VERSION) already \
	  installed.
	@echo You should increase the Release: tag.
else
	./rpmbuild -bb SPECS/CGAL.spec
endif

.PHONY:: CGAL

CGAL: $(CGAL_PACKAGE_PATH)

CGAL.install: $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH)
ifeq ($(findstring not installed,$(CGAL_PACKAGE_INSTALLED_VERSION)),)
	$(SUDO) rpm -Uvh $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH) 
else
	$(SUDO) rpm -ivh $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH)
endif

-include my_makefile
