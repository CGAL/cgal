// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>




#ifndef CGAL_SVD_IS_DEGENERATE_EDGE_2_H
#define CGAL_SVD_IS_DEGENERATE_EDGE_2_H

//#include <CGAL/predicates/Svd_basic_predicates_C2.h>
//#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>

#include <CGAL/predicates/Svd_are_same_points_C2.h>


CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------



template<class K, class Method_tag>
class Svd_is_degenerate_edge_C2
{
public:
  typedef typename K::Site_2      Site_2;

private:
  typedef CGAL::Svd_voronoi_vertex_2<K,Method_tag>  Voronoi_vertex_2;

  typedef Svd_are_same_points_C2<K>   Are_same_points_C2;

  bool is_endpoint(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );
    Are_same_points_C2 same_points;

    return
      same_points(p, s.source_site()) || same_points(p, s.target_site());
  }

public:
  typedef bool          result_type;
  typedef Site_2        argument_type;
  typedef Arity_tag<4>  Arity;

  bool operator()(const Site_2& p, const Site_2& q,
		  const Site_2& r, const Site_2& s) const
  {
    Voronoi_vertex_2 vpqr(p, q, r);
    if ( vpqr.incircle_no_easy(s) == POSITIVE ) { return false; }

    Voronoi_vertex_2 vqps(q, p, s);
    if ( vqps.incircle_no_easy(r) == POSITIVE ) { return false; }

    return true;
  }

};


//-----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_IS_DEGENERATE_EDGE_2_H
