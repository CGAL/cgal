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
// release       : $CGAL_Revision: CGAL-2.4-I-84 $
// release_date  : $CGAL_Date: 2002/05/01 $
//
// file          : include/CGAL/Pm_segment_traits_2.h
// package       : Planar_map (5.109)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Eyal Flato
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//                 Shai Hirsch <shaihi@post.tau.ac.il>
//                 Efi Fogel <efifogel@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_BOPS_TRAITS_2_H
#define CGAL_BOPS_TRAITS_2_H

#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif

#ifndef CGAL_RAY_2_SEGMENT_2_INTERSECTION_H
#include <CGAL/Ray_2_Segment_2_intersection.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Kernel>
class Bops_traits_2
{
public:

  typedef Polygon_2<Kernel>          Polygon_2;
  
  typedef typename Polygon_2::Vertex_iterator        Vertex_iterator;
  typedef typename Polygon_2::Vertex_const_iterator  Vertex_const_iterator;
  typedef typename Polygon_2::Edge_const_iterator    Edge_const_iterator;
  
  typedef typename Polygon_2::Point_2                Point_2;
  typedef typename Polygon_2::Segment_2              Segment_2;
  
  typedef Ray_2<Kernel>                            Ray;
  
public:
  Bops_traits_2() {}
  
  Vertex_iterator vertices_begin(Polygon_2& polygon) 
  {
    return polygon.vertices_begin();
  }
  
  Vertex_iterator vertices_end(Polygon_2& polygon) 
  {
    return polygon.vertices_end();
  }
  
  Vertex_const_iterator vertices_begin(const Polygon_2& polygon) const
  {
    return polygon.vertices_begin();
  }
  
  Vertex_const_iterator vertices_end(const Polygon_2& polygon) const
  {
    return polygon.vertices_end();
  }
  
  Edge_const_iterator edges_begin(const Polygon_2& polygon) const
  {
    return polygon.edges_begin();
  }
  
  Edge_const_iterator edges_end(const Polygon_2& polygon) const
  {
    return polygon.edges_end();
  }  

  //---------------------- Ray interface
  
  bool do_intersect(const Segment_2& segment, const Ray& ray)
  {
    return do_intersect(segment, ray);
  }
  
  bool intersection(const Segment_2& segment, const Ray& ray, Point_2& p)
  {
    Object obj=CGAL::intersection(segment, ray);
    
    return assign(p, obj);
  }
  
};

CGAL_END_NAMESPACE

#endif // CGAL_BOPS_TRAITS_2_H
// EOF

