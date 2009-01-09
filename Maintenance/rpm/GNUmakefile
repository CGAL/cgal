CGAL_URL := $(shell ./get_cgal_trunk_url)

RPM_BUILD_ARCH:=$(shell rpm --showrc | awk '/^build arch/ { print $$4}')

RPMBUILD:=./rpmbuild
EXTRA_RPM_OPTIONS:=
RPM_OPTIONS:=$(EXTRA_RPM_OPTIONS)

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

CGAL_MANUAL_REVISION:=$(shell LANG= LC_ALL= svn info . | awk '/Revision/ {print $$2}')
CGAL_MANUAL_RELEASE_NUMBER:=$(shell rpm -q --qf '%{release}\n' --specfile SPECS/CGAL-manual-tools.spec | head -1)

RPM_OPTIONS+= --define "cgal_manual_revision $(CGAL_MANUAL_REVISION)"

CGAL_MANUAL_PACKAGE_VERSION=$(CGAL_MANUAL_REVISION)-$(CGAL_MANUAL_RELEASE_NUMBER).$(RPM_BUILD_ARCH)

CGAL_MANUAL_PACKAGE_PATH=RPMS/$(RPM_BUILD_ARCH)/CGAL-manual-tools-$(CGAL_MANUAL_PACKAGE_VERSION).rpm
CGAL_MANUAL_SRC_PACKAGE_PATH=SRPMS/CGAL-manual-tools-$(CGAL_MANUAL_REVISION)-$(CGAL_MANUAL_RELEASE_NUMBER).src.rpm

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
	@echo "CGAL_MANUAL_REVISION=$(CGAL_MANUAL_REVISION)"
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
	$(RPMBUILD) $(RPM_OPTIONS) -bb SPECS/CGAL-manual-tools.spec
endif

$(CGAL_MANUAL_SRC_PACKAGE_PATH): $(CGAL_MANUAL_PACKAGE_DEPS)
	$(RPMBUILD) $(RPM_OPTIONS) -bs SPECS/CGAL-manual-tools.spec

CGAL-manual-tools.src: $(CGAL_MANUAL_SRC_PACKAGE_PATH)

.PHONY:: CGAL-manual-tools CGAL-manual-tools.install CGAL-manual-tools.rpmlint
.PHONY:: CGAL-manual-tools.src

CGAL-manual-tools: $(CGAL_MANUAL_PACKAGE_PATH) $(CGAL_MANUAL_SRC_PACKAGE_PATH)

CGAL-manual-tools.install: $(CGAL_MANUAL_PACKAGE_PATH)
ifeq ($(findstring not installed,$(CGAL_MANUAL_PACKAGE_INSTALLED_VERSION)),)
	$(SUDO) rpm -Uvh $(CGAL_MANUAL_PACKAGE_PATH)
else
	$(SUDO) rpm -ivh $(CGAL_MANUAL_PACKAGE_PATH)
endif

CGAL-manual-tools.rpmlint: $(CGAL_MANUAL_PACKAGE_PATH) $(CGAL_MANUAL_SRC_PACKAGE_PATH)
	rpmlint $(CGAL_MANUAL_PACKAGE_PATH) $(CGAL_MANUAL_SRC_PACKAGE_PATH)

#############################################################################
#
#                       CGAL package
#
#############################################################################

CGAL_VERSION:=$(shell rpm --specfile SPECS/CGAL.spec -q --qf '%{version}\n' | head -1)

RPM_OPTIONS += --define "cgal_version $(CGAL_VERSION)"

ifdef CGAL_INTERNAL_RELEASE
RPM_OPTIONS += --define "internal_release $(CGAL_INTERNAL_RELEASE)"
endif

CGAL_RELEASE_NUMBER:=$(shell rpm $(RPM_OPTIONS) --specfile SPECS/CGAL.spec -q --qf '%{release}\n' | head -1)

CGAL_PACKAGE_VERSION=$(CGAL_VERSION)-$(CGAL_RELEASE_NUMBER).$(RPM_BUILD_ARCH)

CGAL_PACKAGE_PATH=RPMS/$(RPM_BUILD_ARCH)/CGAL-$(CGAL_PACKAGE_VERSION).rpm
CGAL_DEVEL_PACKAGE_PATH=RPMS/$(RPM_BUILD_ARCH)/CGAL-devel-$(CGAL_PACKAGE_VERSION).rpm
CGAL_DEMO_PACKAGE_PATH=RPMS/$(RPM_BUILD_ARCH)/CGAL-demo-$(CGAL_PACKAGE_VERSION).rpm
CGAL_SRC_PACKAGE_PATH=SRPMS/CGAL-$(CGAL_VERSION)-$(CGAL_RELEASE_NUMBER).src.rpm

CGAL_PACKAGE_INSTALLED_VERSION=$(shell rpm -q --qf '%{VERSION}-%{RELEASE}.%{ARCH}' CGAL)

CGAL_PACKAGE_DEPS = \
	SPECS/CGAL.spec

$(CGAL_PACKAGE_PATH): $(CGAL_PACKAGE_DEPS)
ifeq ($(CGAL_PACKAGE_VERSION), $(CGAL_PACKAGE_INSTALLED_VERSION))
	@echo CGAL-manual-tools.$(CGAL_PACKAGE_VERSION) already \
	  installed.
	@echo You should increase the Release: tag.
else
	$(RPMBUILD) $(RPM_OPTIONS) -bb SPECS/CGAL.spec
endif

$(CGAL_SRC_PACKAGE_PATH): $(CGAL_PACKAGE_DEPS)
	$(RPMBUILD) $(RPM_OPTIONS) -bs SPECS/CGAL.spec

.PHONY:: CGAL.src CGAL CGAL.install CGAL.rpmlint

CGAL.src: $(CGAL_PACKAGE_DEPS)
	$(RPMBUILD) $(RPM_OPTIONS) -bs --nodeps SPECS/CGAL.spec

CGAL: $(CGAL_PACKAGE_PATH)

CGAL.install: $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH)
ifeq ($(findstring not installed,$(CGAL_PACKAGE_INSTALLED_VERSION)),)
	$(SUDO) rpm -Uvh $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH)
else
	$(SUDO) rpm -ivh $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH)
endif

CGAL.rpmlint: $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH) $(CGAL_SRC_PACKAGE_PATH)
	rpmlint $(CGAL_PACKAGE_PATH) $(CGAL_DEVEL_PACKAGE_PATH) $(CGAL_SRC_PACKAGE_PATH)


-include my_makefile
