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

#ifndef CGAL_SVD_ARE_PARALLEL_C2_H
#define CGAL_SVD_ARE_PARALLEL_C2_H


#include <CGAL/determinant.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------
//                           are parallel
//-----------------------------------------------------------------------

template< class K >
class Svd_are_parallel_C2
{

public:
  typedef typename K::Site_2       Site_2;
  typedef bool                     result_type;
  typedef Arity_tag<2>             Arity;

private:
  typedef typename K::Segment_2    Segment_2;
  typedef typename K::FT           FT;

private:
  bool predicate(const Site_2& p, const Site_2& q) const {
    CGAL_precondition( p.is_segment() && q.is_segment() );
    
    Segment_2 s1 = p.segment();
    Segment_2 s2 = q.segment();

    FT x1 = s1.source().x(),
      y1 = s1.source().y(),
      x2 = s1.target().x(),
      y2 = s1.target().y(),
      x3 = s2.source().x(),
      y3 = s2.source().y(),
      x4 = s2.target().x(),
      y4 = s2.target().y();

    FT det = det2x2_by_formula(x2 - x1, x4 - x3,
			       y2 - y1, y4 - y3);

    return ( CGAL::sign(det) == CGAL::ZERO );
  }

public:
  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    return predicate(p, q);
  }
};


CGAL_END_NAMESPACE


#endif // CGAL_SVD_ARE_PARALLEL_C2_H
