// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//
// examples/Polyhedron/polyhedron_prog_color.C

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_3.h>

// A face type with a color member variable.
template <class Refs>
struct My_face : public CGAL::HalfedgeDS_face_base<Refs> {
    CGAL::Color color;
};

// An items type using my face.
struct My_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef My_face<Refs> Face;
    };
};

typedef CGAL::Cartesian<double>               Kernel;
typedef CGAL::Polyhedron_3<Kernel, My_items>  Polyhedron;
typedef Polyhedron::Halfedge_handle           Halfedge_handle;

int main() {
    Polyhedron P;
    Halfedge_handle h = P.make_tetrahedron();
    h->facet()->color = CGAL::RED;
    return 0;
}
