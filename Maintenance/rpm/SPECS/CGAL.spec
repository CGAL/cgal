# ** This header may define:
# (*)   cgal_name:              the package name. 
#                               Default: CGAL
#                               (may be CGAL-I for internal releases)
# (*)   cgal_version:           the upstream CGAL version.
# (*)   internal_release:       version of the internal release. If 0, 
#                               the release is official.
#       release_number:         release number of this spec file.
#       boost_version:          minimal boost version required by this
#                               version of CGAL
# (*)   build_doc:              decides if the -doc subpackage is build.
# (*)   build_demo:             decides if the -demo subpackage is build.
# (*)   cgal_prefix:            Prefix directory, for the CGAL
#                               installation. If not defined, CGAL is
#                               installed in standard directories
# (*)   set_prefix:             If set_prefix==1, cgal_prefix is filled with a
#                               default value, if empty.
#       cgal_download:          The CGAL download URL.
#
# ** The "(*)" indicates that it can be overiden by the rpmbuild command line,
# this way:
#  rpmbuild --define 'build_doc 1' --define 'internal_release 417' CGAL.spec
#
# ** build_doc is always 0, when internal_release != 0
#
%{!?cgal_version:%define cgal_version 3.3}
%{!?cgal_name: %define cgal_name CGAL}
%{!?internal_release: %define internal_release 0}
%define release_number 13
%define boost_version 1.32
%{!?build_doc: %define build_doc 0}
%{!?build_demo: %define build_demo 1}
%{!?set_prefix: %define set_prefix 0}
%define cgal_download ftp://ftp.mpi-sb.mpg.de/pub/outgoing/CGAL/
# 
# ** End of the header.

# Specific handling of internal releases:
%if 0%{internal_release}
%define tarball_name CGAL-%{cgal_version}-I-%{internal_release}
%define release_value 0.%{internal_release}.%{release_number}
%else
%define tarball_name CGAL-%{cgal_version}
%define release_value %{release_number}
%endif

%if 0%{set_prefix}
%{?!cgal_prefix: %define cgal_prefix %{_libdir}/%{tarball_name}}
%endif

# Installation directories
%{?cgal_prefix: %define install_in_prefix_dir 1}
%{!?cgal_prefix: %define install_in_prefix_dir 0}

%if 0%{install_in_prefix_dir}
  %define cgal_scripts_dir %{cgal_prefix}/scripts
  %define cgal_headers_dir %{cgal_prefix}/include
  %define cgal_demo_src_dir %{cgal_prefix}/demo
  %define cgal_examples_src_dir %{cgal_prefix}/examples
  %define cgal_libs_dir %{cgal_prefix}/%{_lib}
  %define cgal_makefile_dir %{cgal_prefix}/make
%else
  %define cgal_scripts_dir %{_bindir}
  %define cgal_headers_dir %{_includedir}
  %define cgal_demo_src_dir %{_datadir}/CGAL/demo
  %define cgal_examples_src_dir %{_datadir}/CGAL/examples
  %define cgal_libs_dir %{_libdir}
  %define cgal_makefile_dir %{_datadir}/CGAL/make
%endif

# Disable automatic handling of Provides: if %{install_in_prefix_dir} == 1
# so that CGAL does not provides libCGAL.so.x, for example, if prefix is used.
%if 0%{install_in_prefix_dir} == 1
  %define __find_provides %{_builddir}/%{tarball_name}/find_provides.sh
  %define _use_internal_dependency_generator 0
%endif

# The documentation tarball has no license. The -doc subpackage cannot be
# build by the Fedora build system.
%{?fedora: %{!?force_build_doc:%define build_doc 0}}
# Macro force_build_doc, to force the build of doc on Fedora systems
%{?force_build_doc:%define build_doc 1}

# No documentation tarball anyway for internal releases.
%if 0%{internal_release}
%define build_doc 0
%endif

Name:           %{cgal_name}
Version:        %{cgal_version}
Release:        %{release_value}%{?dist}
Summary:        Computational Geometry Algorithms Library

Group:          System Environment/Libraries
License:        QPL/GPL
URL:            http://www.cgal.org/
Source0:        %{cgal_download}%{tarball_name}.tar.gz
%if 0%{build_doc}
Source1:        %{cgal_download}CGAL-%{version}-doc_pdf.tar.gz
Source2:        %{cgal_download}CGAL-%{version}-doc_html.tar.gz
%endif
Source10:       CGAL-README.Fedora
Patch1:         CGAL-install_cgal-SUPPORT_REQUIRED.patch
Patch2:         CGAL-build-library.dpatch

BuildRoot:      %{_tmppath}/CGAL-%{version}-%{release}-root-%(%{__id_u} -n)

# Required packages.
BuildRequires: gmp-devel
BuildRequires: boost-devel >= %boost_version
BuildRequires: qt-devel >= 3.0
BuildConflicts:qt-devel < 4
BuildRequires: zlib-devel

# Requires sub-packages
Requires: %{name}-libs
Requires: %{name}-devel
%if 0%{build_doc}
Requires: %{name}-doc
%endif
%if 0%{build_demo}
Requires: %{name}-demos-source
%endif

%description
Libraries for CGAL applications.
CGAL is a collaborative effort of several sites in Europe and
Israel. The goal is to make the most important of the solutions and
methods developed in computational geometry available to users in
industry and academia in a C++ library. The goal is to provide easy
access to useful, reliable geometric algorithms.

%package libs
Group:          System Environment/Libraries
Summary:        Computational Geometry Algorithms Library libraries
%description libs
Libraries for CGAL applications.
CGAL is a collaborative effort of several sites in Europe and
Israel. The goal is to make the most important of the solutions and
methods developed in computational geometry available to users in
industry and academia in a C++ library. The goal is to provide easy
access to useful, reliable geometric algorithms.

%package devel
Group:          Development/Libraries
Summary:        Development files and tools for CGAL applications
Requires:       %{name}-libs = %{version}-%{release}
Requires:       boost-devel >= %{boost_version}
Requires:       /etc/profile.d
%description devel
The %{name}-devel package provides the headers files and tools you may need to 
develop applications using CGAL.

%if 0%{build_doc}
%package doc
Group:          Documentation
Summary:        HTML and PDF documentation for developing with CGAL
Requires:       %{name}-libs  = %{version}-%{release}
%description doc
The %{name}-doc package provides the html and pdf documentation of CGAL.
%endif

%if 0%{build_demo}
%package demos-source
Group:          Documentation
Summary:        Examples and demos of CGAL algorithms
Requires:       %{name}-libs  = %{version}-%{release}
Obsoletes:      %{name}-demo < %{version}-%{release}
Provides:       %{name}-demo = %{version}-%{release}
%description demos-source
The %{name}-demos-source package provides the sources of examples and demos of
CGAL algorithms.
%endif

%prep
%setup -q -n %{tarball_name}

%patch1 -p0
%patch2 -p1 -b .debian.build-library.back

%if 0%{build_doc}
%setup -q -D -T -a 1
%setup -q -D -T -a 2
%endif

# fix end-of-lines of several files
sed -i 's/\r//' \
    examples/Surface_mesh_parameterization/data/mask_cone.off \
    examples/Boolean_set_operations_2/test.dxf

for f in demo/Straight_skeleton_2/data/vertex_event_9.poly \
         demo/Straight_skeleton_2/data/vertex_event_0.poly;
do
  [ -r $f ] && sed -i 's/\r//' $f;
done

# README.Fedora
install -m 644 %{SOURCE10} %{_builddir}/%{tarball_name}/README.Fedora

# Dummy find_provides
cat > %{_builddir}/%{tarball_name}/find_provides.sh <<EOF
#!/bin/sh

while read dummy; do 
  :
done;
EOF
chmod a+x %{_builddir}/%{tarball_name}/find_provides.sh

%build

source /etc/profile.d/qt.sh

./install_cgal -ni g++ --CUSTOM_CXXFLAGS "$RPM_OPT_FLAGS" \
               --without-autofind \
               --with-ZLIB \
               --with-BOOST \
               --with-BOOSTPROGRAMOPTIONS \
               --with-X11 \
               --with-GMP \
               --with-GMPXX \
               --with-MPFR \
               --with-CGALCORE \
               --with-QT3MT

%install
rm -rf %{buildroot}

# Install headers
mkdir -p %{buildroot}%{cgal_headers_dir}
cp -a include/* %{buildroot}%{cgal_headers_dir}
rm -rf %{buildroot}%{cgal_headers_dir}/CGAL/config/msvc7
mv %{buildroot}%{cgal_headers_dir}/CGAL/config/*/CGAL/compiler_config.h %{buildroot}%{cgal_headers_dir}/CGAL/
rm -rf %{buildroot}%{cgal_headers_dir}/CGAL/config

# Install scripts (only those prefixed with "cgal_").
mkdir -p %{buildroot}%{cgal_scripts_dir}
cp -a scripts/cgal_* %{buildroot}%{cgal_scripts_dir}

# Install libraries
mkdir -p %{buildroot}%{cgal_libs_dir}
cp -a lib/*/lib* %{buildroot}%{cgal_libs_dir}
ln -s libCGAL.so.2.0.0 %{buildroot}%{cgal_libs_dir}/libCGAL.so
ln -s libCGAL.so.2.0.0 %{buildroot}%{cgal_libs_dir}/libCGAL.so.2

# Install makefile:
mkdir -p %{buildroot}%{cgal_makefile_dir}
cp -p make/makefile_* %{buildroot}%{cgal_makefile_dir}/makefile

%if 0%{build_demo}
# Install demos and examples
cp -a demo %{buildroot}%{cgal_demo_src_dir}
cp -a examples %{buildroot}%{cgal_examples_src_dir}
%endif

# Modify makefile
cat > makefile.sed <<'EOF'
s,CGAL_INCL_DIR *=.*,CGAL_INCL_DIR = %{cgal_headers_dir},;
s,CGAL_LIB_DIR *=.*,CGAL_LIB_DIR = %{cgal_libs_dir},;
/CUSTOM_CXXFLAGS/ s/-O2 //;
/CUSTOM_CXXFLAGS/ s/-g //;
/CGAL_INCL_DIR/ s,/CGAL/config/.*,,;
s,/$(CGAL_OS_COMPILER),,g;
/-I.*CGAL_INCL_CONF_DIR/ d
EOF

sed -i -f makefile.sed %{buildroot}%{cgal_makefile_dir}/makefile

# check if the sed script above has worked:
grep -q %{_builddir} %{buildroot}%{cgal_makefile_dir}/makefile && false
grep -q %{buildroot} %{buildroot}%{cgal_makefile_dir}/makefile && false
grep -q CGAL/config %{buildroot}%{cgal_makefile_dir}/makefile && false
grep -q -E 'CUSTOM_CXXFLAGS.*(-O2|-g)' %{buildroot}%{cgal_makefile_dir}/makefile && false

# If CGAL is not installed in a prefix directory, remove -L and -R flags
# from the makefile
%if 0%{install_in_prefix_dir} == 0

cat > makefile-noprefix.sed <<'EOF'
/'-L$(CGAL_LIB_DIR)'/ d;
/-R$(CGAL_LIB_DIR)/ d;
/'-I$(CGAL_INCL_DIR)'/ d;
EOF

sed -i -f makefile-noprefix.sed  %{buildroot}%{cgal_makefile_dir}/makefile

# check that the sed script has worked
grep -q -E -- '-[LI]\$' %{buildroot}%{cgal_makefile_dir}/makefile && false
grep -q -E -- '-R' %{buildroot}%{cgal_makefile_dir}/makefile && false

%endif

# Create /etc/profile.d/ scripts
cd %{buildroot}
mkdir -p ./etc/profile.d
cat > ./etc/profile.d/cgal.sh <<EOF
if [ -z "\$CGAL_MAKEFILE" ] ; then
  CGAL_MAKEFILE="%{cgal_makefile_dir}/makefile"
fi
export CGAL_MAKEFILE
EOF

cat > ./etc/profile.d/cgal.csh <<EOF
if ( \$?CGAL_MAKEFILE ) then
  exit
endif
setenv CGAL_MAKEFILE "$MAKEFILE"
EOF
chmod 755 ./etc/profile.d/cgal.*sh

%clean
rm -rf %{buildroot}

%post libs -p /sbin/ldconfig

%postun libs -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%doc LICENSE* README.Fedora

%files libs
%defattr(-,root,root,-)
%if 0%{install_in_prefix_dir} == 1
%dir %{cgal_prefix}
%dir %{cgal_libs_dir}
%endif
%{cgal_libs_dir}/libCGAL.so.2
%{cgal_libs_dir}/libCGAL.so.2.0.0
%doc LICENSE*

%files devel
%defattr(-,root,root,-)
%if 0%{install_in_prefix_dir} == 1
%dir %{cgal_prefix}
%dir %{cgal_libs_dir}
%dir %{cgal_headers_dir}
%dir %{cgal_scripts_dir}
%endif
%{cgal_headers_dir}/CGAL
%{cgal_headers_dir}/OpenNL
%{cgal_headers_dir}/CORE
%{cgal_libs_dir}/libCGALQt.a
%{cgal_libs_dir}/libcore++.a
%{cgal_libs_dir}/libCGAL.so
%exclude %{cgal_libs_dir}/libCGAL.a
%if 0%{install_in_prefix_dir}
%dir %{_datadir}/CGAL
%endif
%dir %{cgal_makefile_dir}
%config(noreplace) %{cgal_makefile_dir}/makefile
%{cgal_scripts_dir}/*
%if 0%{install_in_prefix_dir}
/etc/profile.d/cgal.*
%else
%exclude /etc/profile.d/cgal.*
%endif
%doc LICENSE*

%if 0%{build_doc}
%files doc
%defattr(-,root,root,-)
%doc doc_html
%doc doc_pdf
%doc LICENSE*
%endif

%if 0%{build_demo}
%files demos-source
%defattr(-,root,root,-)
%doc LICENSE*
%{cgal_demo_src_dir}
%{cgal_examples_src_dir}
%endif

%changelog
* Mon Jul 17 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-13
- Remove unneeded  -R/-L/-I flags from %%{_datadir}/CGAL/make/makefile

* Mon Jul 17 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-12
- Fix %%{cgal_prefix} stuff!!
- Quote 'EOF', so that the lines are not expanded by the shell.

* Thu Jul 13 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-11
- soname is now libCGAL.so.2

* Tue Jul  4 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-10
- Fix makefile.sed so that %%{buildroot} does not appear in 
  %%{_datadir}/CGAL/make/makefile.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-9
- Remove Obsoletes: in the meta-package CGAL.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-8
- Fix the localisation of demo and examples.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-6
- Set Requires, in sub-packages.

* Sun Jul  2 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-5
- Sub-package "demo" is now named "demos-source" (Fedora guidelines).
- Fix some rpmlint warnings
- Added README.Fedora, to explain why the documentation is not shipped, and how CGAL is divided in sub-packages.


* Sat Jul  1 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-4
- Use %%{_datadir}/CGAL instead of %%{_datadir}/%%{name}-%%{version}
- Fix %%{_datadir}/CGAL/makefile, with a sed script.
- Added a new option %%set_prefix (see top of spec file).

* Sat Jul  1 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-3
- Use less "*" in %%files, to avoid futur surprises.
- Remove /etc/profile.d/cgal.* from %%files if %%cgal_prefix is not empty.
- Fix %%build_doc=0 when %%fedora is set. New option macro: %%force_build_doc.

* Fri Jun 30 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-2
- Fix some end-of-lines in %%prep, to please rpmlint.

* Mon May 22 2006 Laurent Rineau <laurent.rineau__fedora_extras@normalesup.org> - 3.3-1
- Remove README from %%doc file: it describes the tarball layout.
- Updated to CGAL-3.3.
- Added examples in the -demo subpackage.
- Cleaning up, to follow Fedora Guidelines.
- The -doc subpackage cannot be build on Fedora (no license).
- Add ldconfig back.
- No prefix.

* Fri Apr 28 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.447
- Update to CGAL-3.2-447.

* Fri Apr 21 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.440
- Updated to CGAL-3.2-I-440.

* Wed Apr 19 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.438
- Added a patch to install_cgal, to require support for BOOST, BOOST_PROGRAM_OPTIONS, X11, GMP, MPFR, GMPXX, CORE, ZLIB, and QT.
- Move scripts to %%{_bindir}
- %%{_libdir}/CGAL-I now belong to CGAL and CGAL-devel, so that it disappears when the packages are removed.

* Wed Apr 12 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org> - 3.2-0.431
- Updated to CGAL-3.2-I-431.
- Remove the use of ldconfig.
- Changed my email address.
- No longer need for patch0.
- Pass of rpmlint.
- Remove unneeded Requires: tags (rpm find them itself).
- Change the release tag.
- Added comments at the beginning of the file.
- Added custom ld flags, on 64bits archs (so that X11 is detected).

* Tue Apr 11 2006 Laurent Rineau <laurent.rineau__fc_extra@normalesup.org>
- Removed -g and -O2 from CUSTOM_CXXFLAGS, in the makefile only.
  They are kept during the compilation of libraries.
- Added zlib in dependencies.
- Added a patch to test_ZLIB.C, until it is merged upstream.

* Fri Mar 31 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- adding a test in the setup section.

* Mon Mar 13 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- delete the patch that fixes the perl path.
- add build_doc and build_demo flags.

* Fri Mar 10 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- Adding new sub-packages doc(pdf&html) and demo.
- Add internal_release flag. 

* Thu Mar 09 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- Cleanup a specfile.

