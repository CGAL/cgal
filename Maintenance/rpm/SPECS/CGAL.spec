%define internal_release 0
%define build_doc 0
%define build_demo 0

%if %{internal_release}
%define tarball_name CGAL-3.2-I-%{internal_release}
%define release %{internal_release}
%define CGAL_DIR %{_libdir}/CGAL-3.2-I
%else
%define tarball_name CGAL-3.2
%define release 1
%define CGAL_DIR %{_libdir}/CGAL-3.2
%endif

%define boost_version 1.32

Summary: Computational Geometry Algorithms Library
Name: CGAL
Version: 3.2
Release: %{release}
License: QPL/LGPL
URL: http://www.cgal.org/
Group: System Environment/Libraries
Source: %{tarball_name}.tar.gz
%if %{build_doc}
Source1: CGAL-3.2-doc_pdf.tar.gz
Source2: CGAL-3.2-doc_html.tar.gz
%endif
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root

#Prefix: %{CGAL_DIR}

Prereq: /sbin/ldconfig
Prereq: fileutils

# Required packages.
Requires: boost >= %boost_version
Requires: gmp
Requires: qt >= 3.0
Requires: libstdc++


BuildRequires: gmp-devel
BuildRequires: boost-devel >= %boost_version
BuildRequires: gcc-c++
BuildRequires: libstdc++-devel
BuildRequires: qt-devel >= 3.0



%description
Libraries for CGAL applications.
CGAL is a collaborative effort of several sites in Europe and
Israel. The goal is to make the most important of the solutions and
methods developed in computational geometry available to users in
industry and academia in a C++ library. The goal is to provide easy
access to useful, reliable geometric algorithms.

%package devel
Group: Development/Libraries
Summary: Development files and tools for %name applications
Requires: boost-devel >= %boost_version, gmp-devel
Requires: %{name}  = %{version}-%{release}

%description devel
The %{name}-devel package provides the headers files and tools you may need to 
develop applications using %name.

%if %{build_doc}
%package doc
Group: Documentation
Summary: HTML and PDF documentation for developing with %name

%description doc
The %{name}-doc package provides the html and pdf documentation of %name.

%endif

%if %{build_demo}
%package demo
Group: Development/Libraries
Summary: demo of %name algorithms.
Requires: %{name}-devel  = %{version}-%{release}

%description demo
The %{name}-demo package provides some demos of %name algorithms.(to be compiled)

%endif

%prep
%if%{internal_release}
%setup -n %{name}-%{version}-I-%{internal_release}
%else
%setup -n %{name}-%{version}
%endif

%if %{build_doc}
%setup -D -T -a 1
%setup -D -T -a 2
%endif
 
%build
rm -rf $RPM_BUILD_ROOT
./install_cgal -ni g++  --prefix ${RPM_BUILD_ROOT}%{CGAL_DIR} --CUSTOM_CXXFLAGS "$RPM_OPT_FLAGS"

%install
CGAL_OS=`./install_cgal -os g++`
MAKEFILE=makefile_${CGAL_OS}
ln -s $MAKEFILE ${RPM_BUILD_ROOT}%{CGAL_DIR}/make/makefile
sed -i "s,$RPM_BUILD_ROOT,,g; /CUSTOM_CXXFLAGS/ s/-O2 //; /CUSTOM_CXXFLAGS/ s/-g //" ${RPM_BUILD_ROOT}%{CGAL_DIR}/make/$MAKEFILE

MAKEFILE=%{CGAL_DIR}/make/makefile
cd $RPM_BUILD_ROOT
mkdir -p ./etc/profile.d
cat > ./etc/profile.d/cgal.sh <<EOF
if [ -z "\$CGAL_MAKEFILE" ] ; then
	CGAL_MAKEFILE="$MAKEFILE"
	CGAL_DIR="%{CGAL_DIR}"
fi
export CGAL_MAKEFILE CGAL_DIR
EOF

cat > ./etc/profile.d/cgal.csh <<EOF

if ( \$?CGAL_MAKEFILE ) then
         exit
endif
setenv CGAL_MAKEFILE "$MAKEFILE"
setenv CGAL_DIR "%{CGAL_DIR}"
EOF
chmod 755 ./etc/profile.d/cgal.*sh
%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%dir %{CGAL_DIR}
%dir %{CGAL_DIR}/lib
%dir %{CGAL_DIR}/lib/*
%{CGAL_DIR}/lib/*/libCGAL.so
%doc LICENSE*
%doc README

%files devel
%defattr(-,root,root,-)
%{CGAL_DIR}/include
%{CGAL_DIR}/lib/*/*.a
%{CGAL_DIR}/make
%{CGAL_DIR}/bin
/etc/profile.d/cgal.*
%doc LICENSE*
%doc README

%if %{build_doc}
%files doc
%defattr(-,root,root,-)
%doc doc_html
%doc doc_pdf
%doc LICENSE*
%endif

%if %{build_demo}
%files demo
%defattr(-,root,root,-)
%doc demo
%doc LICENSE*
%doc README
%endif

%changelog
* Tue Apr 11 2006 Laurent Rineau <laurent.rineau@ens.fr>
- Removed -g and -O2 from CUSTOM_CXXFLAGS, in the makefile only.
  They are kept during the compilation of libraries.

* Fri Mar 31 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- adding a test in the setup section.
* Mon Mar 13 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- delete the patch that fixes the perl path.
- add build_doc and build_demo flags.
* Fri Mar 10 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- adding new sub-packages doc(pdf&html) and demo.
- add internal_release flag. 
* Thu Mar 09 2006 Naceur MESKINI <nmeskini@sophia.inria.fr>
- cleanup a specfile.

