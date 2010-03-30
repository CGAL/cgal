// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France) and
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_ORIENTATION_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_ORIENTATION_C2_H

#include <CGAL/Segment_Delaunay_graph_2/Basic_predicates_C2.h>

CGAL_BEGIN_NAMESPACE

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------



template<class K>
class Orientation_C2
  : private Basic_predicates_C2<K>
{
private:
  typedef Basic_predicates_C2<K>              Base;
  
public:
  typedef typename Base::Orientation          Orientation;

private:
  typedef typename Base::Point_2              Point_2;
  typedef typename Base::Segment_2            Segment_2;
  typedef typename Base::Site_2               Site_2;
  typedef typename Base::FT                   FT;
  typedef typename Base::RT                   RT;

  typedef typename Base::Comparison_result    Comparison_result;
  typedef typename Base::Oriented_side        Oriented_side;
  typedef typename Base::Sign                 Sign;

  typedef typename K::Kernel::Orientation_2   Orientation_2;

public:
  typedef Orientation                  result_type;
  typedef Site_2                       argument_type;

  inline
  Orientation operator()(const Point_2& p, const Point_2& q,
			 const Point_2& r) const
  {
    return Orientation_2()(p, q, r);
  }

  inline
  Orientation operator()(const Site_2& p, const Site_2& q,
			 const Site_2& r) const
  {
    CGAL_precondition( p.is_point() && q.is_point() && r.is_point() );

    return Orientation_2()(p.point(), q.point(), r.point());
  }
};


//-----------------------------------------------------------------------------

CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_ORIENTATION_C2_H
