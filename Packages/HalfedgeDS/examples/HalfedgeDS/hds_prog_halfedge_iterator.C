// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//
// examples/HalfedgeDS/hds_prog_halfedge_iterator.C
// ------------------------------------------------
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>

struct Traits { typedef int Point_2; };
typedef CGAL_HALFEDGEDS_DEFAULT<Traits> HDS;
typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;
typedef HDS::Halfedge_iterator          Iterator;

int main() {
    HDS hds;
    Decorator decorator(hds);
    decorator.create_loop();
    decorator.create_segment();
    CGAL_assertion( decorator.is_valid());
    int n = 0;
    for ( Iterator i = hds.halfedges_begin(); i != hds.halfedges_end(); ++i )
        ++n;
    CGAL_assertion( n == 4);  // == 2 edges
    return 0;
}
