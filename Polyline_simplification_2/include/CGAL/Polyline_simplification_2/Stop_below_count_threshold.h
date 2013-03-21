// Copyright (c) 2012 Geometry Factory. All rights reserved.
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
// Author(s)     : Andreas Fabri, Fernando Cacciola
//
#ifndef CGAL_POLYLINE_SIMPLIFICATION_2_STOP_BELOW_COUNT_THRESHOLD_H
#define CGAL_POLYLINE_SIMPLIFICATION_2_STOP_BELOW_COUNT_THRESHOLD_H

namespace CGAL {

namespace Polyline_simplification_2
{

/// This class is a stop predicate returning the number of vertices is smaller than a certain threshold.
///
/// @heading Is Model for the Concepts: 'PolylineSimplificationStopPredicate'.
class Stop_below_count_threshold
{
public :
  
  /// Initializes it with the given threshold value
  Stop_below_count_threshold( std::size_t threshold ) : mThres(threshold) {}
  
  /// Returns true when "current_count" is smaller or equal than the threshold
  template<class ConstrainedDelaunayTriangulation, class VertexHandle>  
  bool operator()( ConstrainedDelaunayTriangulation const& cdt
                 , VertexHandle                     const& p
                 , double                                  cost
                 , std::size_t                             initial_count
                 , std::size_t                             current_count
                 ) const 
  {
    return current_count <= mThres ;
  }
  
private:
  
  std::size_t mThres ;
};    

} // namespace Polyline_simplification_2

} //namespace CGAL

#endif // CGAL_POLYLINE_SIMPLIFICATION_2_STOP_BELOW_COUNT_THRESHOLD_H
 
