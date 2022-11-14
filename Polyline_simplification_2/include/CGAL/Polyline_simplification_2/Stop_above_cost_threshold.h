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
#ifndef CGAL_POLYLINE_SIMPLIFICATION_2_STOP_ABOVE_COST_THRESHOLD_H
#define CGAL_POLYLINE_SIMPLIFICATION_2_STOP_ABOVE_COST_THRESHOLD_H

#include <CGAL/license/Polyline_simplification_2.h>


#include <CGAL/Constrained_triangulation_plus_2.h>

namespace CGAL {

namespace Polyline_simplification_2
{

/// \ingroup PkgPolylineSimplification2Classes


/// This class is a stop predicate returning `true` when the cost for
/// simplifying a vertex is greater than a certain threshold.
///
/// \cgalModels `PolylineSimplificationStopPredicate`.
class Stop_above_cost_threshold
{
public :

  /// Initializes it with the given threshold value.
  Stop_above_cost_threshold( double threshold ) : mThres(threshold) {}

  /// Returns `true` when `cost` is smaller or equal than the threshold.
  /// \tparam CDT  must be `CGAL::Constrained_Delaunay_triangulation_2` with a vertex type that
  /// is model of  `PolylineSimplificationVertexBase_2`.

  template<class CDT>
  bool operator()(const Constrained_triangulation_plus_2<CDT>&
                  , typename Constrained_triangulation_plus_2<CDT>::Vertex_handle
                  , typename CDT::Geom_traits::FT          cost
                 , std::size_t
                 , std::size_t
                 ) const
  {
    return cost >= mThres ;
  }

private:

  double mThres ;
};

} // namespace Polyline_simplification_2

} //namespace CGAL

#endif // CGAL_POLYLINE_SIMPLIFICATION_2_STOP_ABOVE_COST_THRESHOLD_H

