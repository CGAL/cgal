// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_USING_NAMESPACE_BUG.C
// package       : Configuration (1.28)
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_USING_NAMESPACE_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler / STL implementation
// whether it supports the new standard headers (i.e. without the .h suffix)
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_USING_NAMESPACE_BUG is set, if a compiler does not 
//| not how to compile the following code.
//| Created to workaround a cl1300 bug

namespace CGAL {
  namespace CommonFunctors {
    template < class K >
    struct F {};
  }

  namespace CartesianFunctors {
    using namespace CommonFunctors;
  }

  template < class FT >
  struct Cartesian {
    typedef Cartesian<FT>  Self;
    typedef CartesianFunctors::F<Self> Func;
  };

} // end namespace CGAL

int main() {
  CGAL::Cartesian<double>::Func f;
  (void) f;
  return 0;
}