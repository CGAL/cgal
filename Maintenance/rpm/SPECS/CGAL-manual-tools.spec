Name:           CGAL-manual-tools
Version:        30080
Release:        1%{?dist}
Summary:        CC Manual Style and LaTeX Converter for HTML.

Group:          Development/Tools
License:        N/A
URL:            http://www.cgal.org/Members/Manual_tools/
Source0:        Manual_tools-%{version}.tar.gz
Source1:        Manual-%{version}.tar.gz
Patch0:         CGAL_manual_tools-config.patch
Patch1:         CGAL_manual_tools-rpm.patch
Patch2:         CGAL_manual_tools-perl.patch
Patch3:         CGAL_manual_tools-cgal_manual.patch
Patch4:         CGAL_manual_tools-cc_extract_html.patch
Patch5:         CGAL_manual_tools-latex_to_html.patch
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires: bison flex
BuildRequires:  /usr/bin/kpsewhich
BuildRequires:  sed >= 3.95
Requires:       tetex-latex tetex-dvips 
Requires:       ghostscript >= 6.0
Requires(post):	  /usr/bin/texhash
Requires(postun): /usr/bin/texhash

%description
Specification and Manual Writing Tools for C++ Reference Manuals

%prep
%setup -q -n Manual_tools
%setup -q -n Manual_tools -a 1
%patch0 -p0
%patch1 -p0
%patch2 -p0
%patch3 -p0
%patch4 -p0
%patch5 -p0

%build
source install.config
make -C src LATEX_CONV_INPUTS=$LATEX_CONV_INPUTS \
            CXXFLAGS="${CXXFLAGS:-%optflags}" || exit 1

%install
rm -rf $RPM_BUILD_ROOT
./install.sh
[ -d $RPM_BUILD_ROOT/usr/share/texmf/tex/latex/CGAL ] || mkdir -p $RPM_BUILD_ROOT/usr/share/texmf/tex/latex/CGAL
cp -r doc_tex/Manual $RPM_BUILD_ROOT/usr/share/texmf/tex/latex/CGAL
cp doc_tex/ipe.sty $RPM_BUILD_ROOT/usr/share/texmf/tex/latex/CGAL
[ -d $RPM_BUILD_ROOT/usr/share/texmf/bibtex/bib/CGAL ] || mkdir -p $RPM_BUILD_ROOT/usr/share/texmf/bibtex/bib/CGAL/Manual
mv $RPM_BUILD_ROOT/usr/share/texmf/tex/latex/CGAL/Manual/*.bib $RPM_BUILD_ROOT/usr/share/texmf/bibtex/bib/CGAL/Manual
cp developer_scripts/cgal_manual developer_scripts/bibmerge $RPM_BUILD_ROOT/usr/bin/ 

%post
texhash > /dev/null 2>&1 || :

%postun
texhash > /dev/null 2>&1 || :

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
/usr/bin/*
%dir /usr/share/texmf/tex/latex/CGAL
%dir /usr/share/texmf/bibtex/bib/CGAL
%dir /usr/share/texmf/bibtex/bib/CGAL/Manual
%dir /usr/share/texmf/tex/latex/CGAL/eps_tabs/
%dir /usr/share/CGAL_latex_conv_config
%dir /usr/share/CGAL_latex_conv_config/html
%dir /usr/share/CGAL_latex_conv_config/gif
/usr/share/texmf/tex/latex/CGAL
/usr/share/texmf/bibtex/bib/CGAL
/usr/share/texmf/tex/latex/CGAL/eps_tabs/
/usr/share/CGAL_latex_conv_config
/usr/share/CGAL_latex_conv_config/html
/usr/share/CGAL_latex_conv_config/gif

%changelog
* Thu Apr  6 2006 Laurent Rineau <laurent.rineau@ens.fr> - 30080
- Revision 30080.

* Wed Apr  5 2006 Laurent Rineau <laurent.rineau@ens.fr> - 30025
- Updated to revision 30025.
- Source0 and Source1 now have %version in their names.

* Tue Apr  4 2006 Laurent Rineau <laurent.rineau@ens.fr> - 29933
- Updated to revision 29933

* Mon Apr  3 2006 Laurent Rineau <laurent.rineau@ens.fr> - 29922
- Updated to revision 29922

* Fri Mar 31 2006 Laurent Rineau <laurent.rineau@ens.fr> - 29890
- Updated to 29890.
- New patch, for latex_to_html
- Updated the CGAL_manual_tools-config.patch, to add
  $LATEX_CONV_CONFIG/html to $LATEX_CONV_INPUTS
- New patch, for cc_extract_html.C. Will up submitted upstream soon.
- Added ipe.sty to files
- Move bib files to the bibtex installation.

* Thu Mar 30 2006 Laurent Rineau <laurent.rineau@ens.fr> - 29841
- Updated to rev. 29841. No longer need for patch2
- Specify the two scripts of developer_scripts/ that need to be installed,
to avoid things like cgal_manual.orig!

* Tue Mar  7 2006 Laurent Rineau <laurent.rineau@ens.fr> - 29111
- Updated to Manual_tools rev 29111., Manual rev 29112.
- Manual_tools are now pulled from SVN too.
- Some patches are now included in SVN.

* Fri Mar  3 2006 Laurent Rineau <laurent.rineau@ens.fr> - 3.12
- Added cgal_manual and dependancies
- Added a call to texhash

