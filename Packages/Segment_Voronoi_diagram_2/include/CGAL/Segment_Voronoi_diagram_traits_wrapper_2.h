// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Segment_Voronoi_diagram_traits_wrapper_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_WRAPPER_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_WRAPPER_2_H

CGAL_BEGIN_NAMESPACE


template<class Gt_base>
class Segment_Voronoi_diagram_traits_wrapper_2 : public Gt_base
{
private:
  typedef Segment_Voronoi_diagram_traits_wrapper_2<Gt_base> Self;

public:
  struct Triangle_2 {};

  Segment_Voronoi_diagram_traits_wrapper_2() {}

  Segment_Voronoi_diagram_traits_wrapper_2
  (const Self&) {}

  Segment_Voronoi_diagram_traits_wrapper_2
  operator=(const Self&) {
    return (*this);
  }

  Segment_Voronoi_diagram_traits_wrapper_2(const Gt_base&) {}
};




CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_TRAITS_WRAPPER_2_H
