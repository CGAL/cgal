// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_MATCHING_BUG_4.C
// revision      : 1.0
// revision_date : 5 August 2003
// author(s)     : Michael Hoffmann
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_MATCHING_BUG_4.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if a compiler cannot distinguish the signature
//| of overloaded function templates, which have arguments whose type
//| depends on the template parameter.
//| This bug apopears for example on Sunpro 5.3 and 5.4.

template <class K> void foo(typename K::P, typename K::L, K) {}

template <class K> void foo(typename K::L, typename K::P, K) {}

int main() { return 0; }

// EOF //
