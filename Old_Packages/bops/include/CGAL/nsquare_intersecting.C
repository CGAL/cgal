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
// file          : include/CGAL/nsquare_intersecting.C
// package       : bops (2.2)
// source        : include/CGAL/nsquare_intersecting.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Carl Van Geem <Carl.Van.Geem@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

CGAL_BEGIN_NAMESPACE

template < class R > 
Intersectionresult<R>::Intersectionresult() 
      {   
        _is_vertex_of_poly1 = false;
        _is_vertex_of_poly2 = false;
        _is_edge_of_poly1 = false;
        _is_edge_of_poly2 = false;
        _is_vertex_vertex_intersection = false;
        _is_vertex_edge_intersection = false;
        _is_edge_edge_intersection = false;
      }

template < class R > 
Intersectionresult<R>::Intersectionresult (const 
                                        Intersectionresult<R> &ires )
     {  
       _intersection_object = ires._intersection_object;
       _segments_poly1 = ires._segments_poly1;
       _segments_poly2 = ires._segments_poly2;
       _is_vertex_of_poly1 = ires._is_vertex_of_poly1;
       _is_vertex_of_poly2 = ires._is_vertex_of_poly2;
       _is_edge_of_poly1 = ires._is_edge_of_poly1;
       _is_edge_of_poly2 = ires._is_edge_of_poly2;
       _is_vertex_vertex_intersection = 
                 ires._is_vertex_vertex_intersection;
       _is_vertex_edge_intersection = ires._is_vertex_edge_intersection;
       _is_edge_edge_intersection = ires._is_edge_edge_intersection;
     }

template < class R > 
Intersectionresult<R>::Intersectionresult (
                           const Point_2<R> &iobj,
                           const Segment_2<R> &iseg1,
                           const Segment_2<R> &iseg2,
                           const bool &i_is_vertex_of_poly1,
                           const bool &i_is_vertex_of_poly2,
                           const bool &i_is_edge_of_poly1,
                           const bool &i_is_edge_of_poly2,
                           const bool &i_is_vertex_vertex_intersection,
                           const bool &i_is_vertex_edge_intersection,
                           const bool &i_is_edge_edge_intersection) 
     {
       _intersection_object = make_object(iobj);
       _segments_poly1.push_back(iseg1);
       _segments_poly2.push_back(iseg2);
       _is_vertex_of_poly1 = i_is_vertex_of_poly1;
       _is_vertex_of_poly2 = i_is_vertex_of_poly2;
       _is_edge_of_poly1 = i_is_edge_of_poly1;
       _is_edge_of_poly2 = i_is_edge_of_poly2;
       _is_vertex_vertex_intersection = i_is_vertex_vertex_intersection;
       _is_vertex_edge_intersection = i_is_vertex_edge_intersection;
       _is_edge_edge_intersection = i_is_edge_edge_intersection;
     }

template < class R > 
Intersectionresult<R>::Intersectionresult (
                           const Segment_2<R> &iobj,
                           const Segment_2<R> &iseg1,
                           const Segment_2<R> &iseg2,
                           const bool &i_is_vertex_of_poly1,
                           const bool &i_is_vertex_of_poly2,
                           const bool &i_is_edge_of_poly1,
                           const bool &i_is_edge_of_poly2,
                           const bool &i_is_vertex_vertex_intersection,
                           const bool &i_is_vertex_edge_intersection,
                           const bool &i_is_edge_edge_intersection) 
     {
       _intersection_object = make_object(iobj);
       _segments_poly1.push_back(iseg1);
       _segments_poly2.push_back(iseg2);
       _is_vertex_of_poly1 = i_is_vertex_of_poly1;
       _is_vertex_of_poly2 = i_is_vertex_of_poly2;
       _is_edge_of_poly1 = i_is_edge_of_poly1;
       _is_edge_of_poly2 = i_is_edge_of_poly2;
       _is_vertex_vertex_intersection = i_is_vertex_vertex_intersection;
       _is_vertex_edge_intersection = i_is_vertex_edge_intersection;
       _is_edge_edge_intersection = i_is_edge_edge_intersection;
     }


template < class R > 
Intersectionresult<R>& Intersectionresult<R>::operator=(
                                 const Intersectionresult<R>  &ires)  
     {
       _intersection_object = ires._intersection_object;
       _segments_poly1 = ires._segments_poly1;
       _segments_poly2 = ires._segments_poly2;
       _is_vertex_of_poly1 = ires._is_vertex_of_poly1;
       _is_vertex_of_poly2 = ires._is_vertex_of_poly2;
       _is_edge_of_poly1 = ires. _is_edge_of_poly1;
       _is_edge_of_poly2 = ires._is_edge_of_poly2;
       _is_vertex_vertex_intersection = ires._is_vertex_vertex_intersection;
       _is_vertex_edge_intersection = ires._is_vertex_edge_intersection;
       _is_edge_edge_intersection = ires._is_edge_edge_intersection;
       return *this;
     }

template < class R > 
Intersectionresult<R>::~Intersectionresult() 
     { 
     }



template < class R > 
List_of_intersections<R>::List_of_intersections() 
{
}
template < class R > 
List_of_intersections<R>::List_of_intersections(
                      const List_of_intersections<R> &ilist) 
{
       _list_of_intersections = ilist._list_of_intersections;
}
template < class R > 
List_of_intersections<R>& List_of_intersections<R>::operator=(
                              const List_of_intersections<R> &ilist) 
{
       _list_of_intersections = ilist._list_of_intersections;
       return *this;
}
template < class R > 
List_of_intersections<R>::~List_of_intersections() 
{
}

CGAL_END_NAMESPACE
