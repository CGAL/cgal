// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//
// examples/HalfedgeDS/hds_prog_vector.C              
// -------------------------------------
#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/HalfedgeDS_decorator.h>

struct Traits { typedef int Point_2; };
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
typedef CGAL::HalfedgeDS_vector     < Traits, CGAL::HalfedgeDS_items_2> HDS;
#else
typedef CGAL::HalfedgeDS_vector::HDS< Traits, CGAL::HalfedgeDS_items_2> HDS;
#endif
typedef CGAL::HalfedgeDS_decorator<HDS>  Decorator;

int main() {
    HDS hds(1,2,2);
    Decorator decorator(hds);
    decorator.create_loop();
    CGAL_assertion( decorator.is_valid());
    return 0;
}
