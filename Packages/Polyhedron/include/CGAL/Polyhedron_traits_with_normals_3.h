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
// file          : include/CGAL/Polyhedron_traits_with_normals_3.h
// package       : Polyhedron 2.9 (13 Sep 2000)
// chapter       : 3D-Polyhedral Surfaces
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>)
// maintainer    :
// coordinator   : MPI Saarbruecken
//
// Default Traits for Polyhedral Surfaces in CGAL.
// ============================================================================

#ifndef CGAL_POLYHEDRON_TRAITS_WITH_NORMALS_3_H
#define CGAL_POLYHEDRON_TRAITS_WITH_NORMALS_3_H 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template < class Kernel_ >
class Polyhedron_traits_with_normals_3 {
public:
    typedef Kernel_                   Kernel;
    typedef typename Kernel::Point_3  Point_3;
    typedef typename Kernel::Vector_3 Plane_3;

    typedef typename Kernel::Construct_opposite_vector_3 
                                      Construct_opposite_plane_3;
private:
    Kernel m_kernel;

public:
    Polyhedron_traits_with_normals_3() {}
    Polyhedron_traits_with_normals_3( const Kernel& kernel)
        : m_kernel(kernel) {}

    Construct_opposite_plane_3 construct_opposite_plane_3_object() const {
        return m_kernel.construct_opposite_vector_3_object();
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_POLYHEDRON_TRAITS_WITH_NORMALS_3_H //
// EOF //
