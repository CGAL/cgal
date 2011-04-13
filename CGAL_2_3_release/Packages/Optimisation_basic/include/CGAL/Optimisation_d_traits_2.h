// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Optimisation_d_traits_2.h
// package       : $CGAL_Package: Optimisation_basic $
// chapter       : Geometric Optimisation
//
// source        : web/Optimisation_d_traits.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Traits class (2D) for dD optimisation algorithms
// ============================================================================

#ifndef CGAL_OPTIMISATION_D_TRAITS_2_H
#define CGAL_OPTIMISATION_D_TRAITS_2_H

// includes
#ifndef CGAL_OPTIMISATION_ACCESS_DIMENSION_2_H
#  include <CGAL/Optimisation/Access_dimension_2.h>
#endif
#ifndef CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_2_H
#  include <CGAL/Optimisation/Access_coordinates_begin_2.h>
#endif
#ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_2_H
#  include <CGAL/Optimisation/Construct_point_2.h>
#endif

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class K_, class ET_ = CGAL_TYPENAME_MSVC_NULL K_::RT,
                     class NT_ = CGAL_TYPENAME_MSVC_NULL K_::RT >
class Optimisation_d_traits_2;

// Class interface
// ===============
template < class K_, class ET_, class NT_>
class Optimisation_d_traits_2 {
  public:
    // self
    typedef  K_                         K;
    typedef  ET_                        ET;
    typedef  NT_                        NT;
    typedef  Optimisation_d_traits_2<K,ET,NT>
                                        Self;

    // types
    typedef  typename K::Point_2        Point_d;

    typedef  typename K::Rep_tag        Rep_tag;

    typedef  typename K::RT             RT;
    typedef  typename K::FT             FT;

    typedef  Access_dimension_2<K>      Access_dimension_d;
    typedef  Access_coordinates_begin_2<K>
                                        Access_coordinates_begin_d;

    typedef  Construct_point_2<K>       Construct_point_d;

    // creation
    Optimisation_d_traits_2( ) { }
    Optimisation_d_traits_2( const Optimisation_d_traits_2<K_,ET_,NT_>&) {}

    // operations
    Access_dimension_d
    access_dimension_d_object( ) const
        { return Access_dimension_d(); }

    Access_coordinates_begin_d
    access_coordinates_begin_d_object( ) const
        { return Access_coordinates_begin_d(); }

    Construct_point_d
    construct_point_d_object( ) const
        { return Construct_point_d(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_OPTIMISATION_D_TRAITS_2_H

// ===== EOF ==================================================================
