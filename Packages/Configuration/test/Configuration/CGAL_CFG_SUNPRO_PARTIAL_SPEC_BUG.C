// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : config/testfiles/CGAL_CFG_SUNPRO_PARTIAL_SPEC_BUG.C
// package       : Configuration (2.3)
// author(s)     : Lutz & Sylvain
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_SUNPRO_PARTIAL_SPEC_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This is a test-case for a bug in SunPro 5.3 that occurs in the HDS.
//| When the bug is present, CGAL_CFG_SUNPRO_PARTIAL_SPEC_BUG is set.

template < class Refs, class D = int >
struct Halfedge_base;
 
template < class Refs >
struct Halfedge_base <Refs, int> {
    typedef Halfedge_base<Refs> Base;
    void           set_vertex( )  { }
};
 
struct HDS {
    typedef Halfedge_base<int>     Halfedge;
    typedef Halfedge::Base          HBase;
 
    void create_pair() {
        Halfedge h;
        h.HBase::set_vertex();
    }
};
 
int main() {
    HDS hds;
    hds.create_pair();
    return 0;
}
