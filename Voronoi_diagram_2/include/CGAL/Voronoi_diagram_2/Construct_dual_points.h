// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_CONSTRUCT_DUAL_POINTS_H
#define CGAL_VORONOI_DIAGRAM_2_CONSTRUCT_DUAL_POINTS_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================

template<class AG2>
class Apollonius_graph_Voronoi_point_2
{
private:
  typedef typename AG2::Geom_traits      Geom_traits;

public:
  typedef typename Geom_traits::Point_2  result_type;
  typedef typename AG2::Face_handle      Face_handle;

  result_type operator()(const Face_handle& f) const {
    return Geom_traits().construct_Apollonius_vertex_2_object()
      (f->vertex(0)->site(), f->vertex(1)->site(), f->vertex(2)->site());
  }
};

//=========================================================================

template<class DT2>
class Delaunay_triangulation_Voronoi_point_2
{
private:
  typedef typename DT2::Geom_traits       Geom_traits;

public:
  typedef typename Geom_traits::Point_2   result_type;
  typedef typename DT2::Face_handle       Face_handle;

  result_type operator()(const Face_handle& f) const {
    return Geom_traits().construct_circumcenter_2_object()
      (f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
  }
};

//=========================================================================

template<class RT2>
class Regular_triangulation_Voronoi_point_2
{
private:
  typedef typename RT2::Geom_traits       Geom_traits;
public:
  typedef typename Geom_traits::Point_2   result_type;
  typedef typename RT2::Face_handle       Face_handle;

  result_type operator()(const Face_handle& f) const {
    return Geom_traits().construct_weighted_circumcenter_2_object()
      (f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
  }
};

//=========================================================================

template<class SVD2>
class Segment_Voronoi_diagram_Voronoi_point_2
{
private:
  typedef typename SVD2::Geom_traits   Geom_traits;

public:
  typedef typename SVD2::Point_2       result_type;
  typedef typename SVD2::Face_handle   Face_handle;

  result_type operator()(const Face_handle& f) const {
    return Geom_traits().construct_svd_vertex_2_object()
      (f->vertex(0)->site(), f->vertex(1)->site(), f->vertex(2)->site());
  }
};

template<class SDG2>
class Segment_Delaunay_graph_Voronoi_point_2
{
private:
  typedef typename SDG2::Geom_traits   Geom_traits;

public:
  typedef typename SDG2::Point_2       result_type;
  typedef typename SDG2::Face_handle   Face_handle;

  result_type operator()(const Face_handle& f) const {
    return Geom_traits().construct_svd_vertex_2_object()
      (f->vertex(0)->site(), f->vertex(1)->site(), f->vertex(2)->site());
  }
};

//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_CONSTRUCT_DUAL_POINTS_H
