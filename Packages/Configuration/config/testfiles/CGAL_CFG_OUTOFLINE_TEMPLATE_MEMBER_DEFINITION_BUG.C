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
// file          : config/testfiles/CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG.C
// package       : Configuration (2.62)
// maintainer    : Geert-Jan Giezeman <geert@cs.uu.nl>
// source        :
// revision      : 
// revision_date : 
// author(s)     : Radu Ursu
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG is set, if 
//| a compiler does not support the definition of the members templates 
//| out of line. The solution is to put the definition inside the class.
//| This is a feature of cl1200 and cl1300.

template<class A>
class B{
  template<class C>
    void fct(C *i);
};

template<class A>
template<class C>  //syntax error
void
B<A>::fct(C *i){
}

int main(){
  return 1;
}