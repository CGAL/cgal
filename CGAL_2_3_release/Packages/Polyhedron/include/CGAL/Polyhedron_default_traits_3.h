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

// Use renamed traits class and provide derived class for backwards
// compatibility.

#include <CGAL/Polyhedron_traits_3.h>

CGAL_BEGIN_NAMESPACE

template < class Kernel_ >
class Polyhedron_default_traits_3 : public Polyhedron_traits_3<Kernel_> {
public:
    typedef Kernel_ Kernel;
    Polyhedron_default_traits_3() {}
    Polyhedron_default_traits_3( const Kernel& kernel) 
        : Polyhedron_traits_3<Kernel>(kernel) {}
};

CGAL_END_NAMESPACE

#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //
#endif // CGAL_POLYHEDRON_DEFAULT_TRAITS_3_H //
// EOF //


