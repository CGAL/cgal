// ======================================================================
//
// Copyright (c) 1999,2001 The CGAL Consortium
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
// file          : include/CGAL/Data_accessorH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_DATA_ACCESSORH2_H
#define CGAL_DATA_ACCESSORH2_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Data_accessorH2
{
public:
    typedef typename R::FT FT;
    typedef typename R::RT RT;
    typedef typename R::Kernel_base::Point_2      Point_2;

    RT  get_hx( const Point_2 & p) const { return( p.hx()); }
    RT  get_hy( const Point_2 & p) const { return( p.hy()); }
    RT  get_hw( const Point_2 & p) const { return( p.hw()); }

    void
    get( const Point_2 & p, RT& hx, RT& hy, RT& hw) const
    {
        hx = get_hx( p);
        hy = get_hy( p);
        hw = get_hw( p);
    }

    void
    set( Point_2& p, const RT & hx, const RT & hy, const RT & hw) const
    {
        p = Point_2( hx, hy, hw);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_DATA_ACCESSORH2_H
