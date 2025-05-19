// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_CONSTRUCT_DUAL_POINTS_H
#define CGAL_VORONOI_DIAGRAM_2_CONSTRUCT_DUAL_POINTS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


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

template<class DToS2>
struct DToS2_Voronoi_point_2
{
private:
  typedef typename DToS2::Geom_traits                      Geom_traits;
  typedef typename Geom_traits::Point_on_sphere_2          Point_on_sphere_2;

public:
  typedef typename DToS2::Face_handle                      Face_handle;
  typedef Point_on_sphere_2                                result_type;

  DToS2_Voronoi_point_2(const Geom_traits& gt) : gt(gt) { }

  result_type operator()(const Face_handle f) const
  {
    return gt.construct_circumcenter_on_sphere_2_object()(
             f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
  }

private:
  const Geom_traits& gt;
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
