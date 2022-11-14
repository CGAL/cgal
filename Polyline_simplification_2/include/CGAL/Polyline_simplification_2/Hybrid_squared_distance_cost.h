// Copyright (c) 2012 Geometry Factory. All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri, Fernando Cacciola
//
#ifndef CGAL_POLYLINE_SIMPLIFICATION_2_HYBRID_SQUARED_DISTANCE_COST_H
#define CGAL_POLYLINE_SIMPLIFICATION_2_HYBRID_SQUARED_DISTANCE_COST_H

#include <CGAL/license/Polyline_simplification_2.h>


#include <CGAL/algorithm.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

namespace CGAL {

namespace Polyline_simplification_2
{
/// \ingroup PkgPolylineSimplification2Classes

/// This class is a cost function which calculates the cost as the square of the distance between the original and simplified polylines,
/// possibly scaled based on a factor.
///
/// \cgalModels `PolylineSimplificationCostFunction`.
template<class FT>
class Hybrid_squared_distance_cost
{

public:

  /// Initializes the cost function with the specified `ratio`
  Hybrid_squared_distance_cost( FT ratio ) : mSquaredRatio(ratio*ratio) {}

  /// Compute the hybrid squared distance cost.
  ///
  /// Given a vertex in constraint iterator `vicq` computes `vicp=std::prev(vicq)` and `vicr=std::next(vicq)`,
  /// returns the maximal square distance between each point along the original subpolyline,
  /// between `vicp` and `vicr`,
  /// and the straight line segment from `*vicp->point() to *vicr->point()` divided by the smallest of
  /// - the square of the ratio given to the constructor of the cost function,
  /// - and the shortest squared distance between that segment and each of the vertices adjacent to `vicq`.
  ///
  /// \tparam CDT  must be `CGAL::Constrained_Delaunay_triangulation_2` with a vertex type that
  /// is model of  `PolylineSimplificationVertexBase_2`.
 template<class CDT>
 boost::optional<typename CDT::Geom_traits::FT>
    operator()(  const Constrained_triangulation_plus_2<CDT>& pct
                 , typename Constrained_triangulation_plus_2<CDT>::Vertices_in_constraint_iterator vicq) const
  {
    typedef typename Constrained_triangulation_plus_2<CDT>::Points_in_constraint_iterator Points_in_constraint_iterator;
    typedef typename Constrained_triangulation_plus_2<CDT>::Vertex_circulator Vertex_circulator;
    typedef typename Constrained_triangulation_plus_2<CDT>::Geom_traits Geom_traits ;
    typedef typename Geom_traits::Compute_squared_distance_2 Compute_squared_distance;
    typedef typename Geom_traits::Construct_segment_2        Construct_segment;
    typedef typename Geom_traits::Segment_2                  Segment;
    typedef typename Geom_traits::Point_2                    Point;

    Compute_squared_distance compute_squared_distance = pct.geom_traits().compute_squared_distance_2_object() ;
    Construct_segment        construct_segment        = pct.geom_traits().construct_segment_2_object() ;
    typedef typename Constrained_triangulation_plus_2<CDT>::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;

    Vertices_in_constraint_iterator vicp = boost::prior(vicq);
    Vertices_in_constraint_iterator vicr = boost::next(vicq);

    Point const& lP = (*vicp)->point();
    Point const& lR = (*vicr)->point();

    Segment lP_R = construct_segment(lP, lR) ;

    FT d1(0);
    Points_in_constraint_iterator pp(vicp), rr(vicr);
    ++pp;

    for ( ;pp != rr; ++pp )
      d1 = (std::max)(d1, compute_squared_distance( lP_R, *pp ) ) ;

    FT d2 = mSquaredRatio;

    Vertex_circulator vc = pct.incident_vertices(*vicq), done(vc);
    do {
      if((vc != pct.infinite_vertex()) && (vc != *vicp) && (vc != *vicr)){
        d2 = (std::min)(d2, compute_squared_distance(vc->point(), (*vicq)->point()));
      }
      ++vc;
    }while(vc != done);

    return d1 / d2;
  }

private:

  FT mSquaredRatio ;

};


} // namespace Polyline_simplification_2


} //namespace CGAL

#endif // CGAL_POLYLINE_SIMPLIFICATION_2_HYBRID_SQUARED_DISTANCE_COST_H


