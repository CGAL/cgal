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
// file          : Polyhedron_traits_3.h
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

#ifndef CGAL_POLYHEDRON_TRAITS_3_H
#define CGAL_POLYHEDRON_TRAITS_3_H 1

#include <CGAL/basic.h>
// MS Visual C++ 6.0 does not work with the new design.
#if defined( _MSC_VER) && (_MSC_VER <= 1200)
#ifndef CGAL_USE_POLYHEDRON_DESIGN_TWO
#define CGAL_USE_POLYHEDRON_DESIGN_ONE 1
#endif
#endif

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
#include <CGAL/Polyhedron_old/Polyhedron_default_traits_3.h>
#else // CGAL_USE_POLYHEDRON_DESIGN_ONE //
#define CGAL_USE_POLYHEDRON_DESIGN_TWO 1


CGAL_BEGIN_NAMESPACE

template < class Kernel_ >
class Polyhedron_traits_3 {
public:
    typedef Kernel_                   Kernel;
    typedef typename Kernel::Point_3  Point_3;
    typedef typename Kernel::Plane_3  Plane_3;

    typedef typename Kernel::Construct_opposite_plane_3 
                                      Construct_opposite_plane_3;
private:
    Kernel m_kernel;

public:
    Polyhedron_traits_3() {}
    Polyhedron_traits_3( const Kernel& kernel) : m_kernel(kernel) {}

    Construct_opposite_plane_3 construct_opposite_plane_3_object() const {
        return m_kernel.construct_opposite_plane_3_object();
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
#endif // CGAL_POLYHEDRON_TRAITS_3_H //
// EOF //
