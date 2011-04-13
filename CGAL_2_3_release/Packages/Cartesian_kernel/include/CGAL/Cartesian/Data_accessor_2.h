// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Data_accessor_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DATA_ACCESSOR_2_H
#define CGAL_CARTESIAN_DATA_ACCESSOR_2_H

#include <CGAL/Cartesian/redefine_names_2.h>

CGAL_BEGIN_NAMESPACE

// 2D Cartesian point data accessor
template < class R_ >
class Data_accessorC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
{
public:
    typedef R_                           R;
    typedef typename R::FT               FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
    typedef  typename R::Point_2          Point;
#else
    typedef  typename R::Point_2_base     Point;
#endif

    FT  get_x( const Point &p) const { return( p.x()); }
    FT  get_y( const Point &p) const { return( p.y()); }

    void get( const Point &p, FT &x, FT &y) const { x=get_x(p); y=get_y(p); }
    void set( Point& p, const FT &x, const FT &y) const { p=Point(x,y); }
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DATA_ACCESSOR_2_H
