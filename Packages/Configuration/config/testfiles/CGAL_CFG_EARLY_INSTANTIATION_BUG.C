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
// file          : config/testfiles/CGAL_CFG_EARLY_INSTANTIATION_BUG.C
// package       : Configuration (1.28)
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_EARLY_INSTANTIATION_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler / STL implementation
// whether it supports the new standard headers (i.e. without the .h suffix)
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if a compiler instantiates a class template A<>
//| too early, where A<> is declared in the scope of another class
//| template B<>, such that B<> is the template parameter of A<>.
//| This bug shows up on VC7.0.

template < class Container > class Iterator;

template < class A >
struct Container {
  typedef Container<A>    Self;
  typedef Iterator<Self>  iterator;

  // VC++ 7.0 instantiates iterator too early, unless passed by reference.
  // Solution: void erase(iterator&){}
  void erase(iterator) {}
};

template < class Container >
class Iterator {
  typedef typename Container::Self pipo;
};

int main() {
  Container<int> C;
  typedef Container<int>::iterator It;
  It i;
  C.erase(i);
  return 0;
}
