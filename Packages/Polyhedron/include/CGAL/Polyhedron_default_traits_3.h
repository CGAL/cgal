// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : Polyhedron_default_traits_3.h
// chapter       : $CGAL_Chapter: 3D-Polyhedral Surfaces $
// package       : $CGAL_Package: Polyhedron 2.9 (13 Sep 2000) $
// source        : polyhedron.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Default Traits for Polyhedral Surfaces in CGAL.
// ============================================================================

#ifndef CGAL_POLYHEDRON_DEFAULT_TRAITS_3_H
#define CGAL_POLYHEDRON_DEFAULT_TRAITS_3_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif
#ifndef CGAL_VECTOR_3_H
#include <CGAL/Vector_3.h>
#endif
#ifndef CGAL_PLANE_3_H
#include <CGAL/Plane_3.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class Rep >
class Polyhedron_default_traits_3 {
public:
    typedef Rep              R;
    typedef Point_3<Rep>     Point;
    typedef Vector_3<Rep>    Normal;
    typedef Plane_3<Rep>     Plane;
    void reverse_normal( Normal& normal) const { normal = - normal; }
    void reverse_plane( Plane& plane) const { plane  = plane.opposite(); }
};

CGAL_END_NAMESPACE
#endif // CGAL_POLYHEDRON_DEFAULT_TRAITS_3_H //
// EOF //
