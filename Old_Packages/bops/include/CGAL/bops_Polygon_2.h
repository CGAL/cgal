// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_Polygon_2.h
// package       : bops (2.2)
// source        : include/CGAL/bops_Polygon_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_POLYGON_2_H
#define CGAL_BOPS_POLYGON_2_H


#include <CGAL/bops_traits_2.h>

CGAL_BEGIN_NAMESPACE

template < class R, class Container >
bool do_intersect(
     const Polygon_2<Polygon_traits_2<R>, Container>& A,
     const Polygon_2<Polygon_traits_2<R>, Container>& B
);


template < class R, class Container, class OutputIterator >
OutputIterator
intersection( const Polygon_2<Polygon_traits_2<R>, Container>& A,
                   const Polygon_2<Polygon_traits_2<R>, Container>& B,
                   OutputIterator list_of_objects_it);

 
template < class R, class Container, class OutputIterator >
OutputIterator
UnionXS(        const Polygon_2<Polygon_traits_2<R>, Container>& A,
                   const Polygon_2<Polygon_traits_2<R>, Container>& B,
                   OutputIterator list_of_objects_it);
 
 
template < class R, class Container, class OutputIterator >
OutputIterator
difference(   const Polygon_2<Polygon_traits_2<R>, Container>& A,
                   const Polygon_2<Polygon_traits_2<R>, Container>& B,
                   OutputIterator list_of_objects_it);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_Polygon_2.C>
#endif

#endif // CGAL_BOPS_POLYGON_2_H
