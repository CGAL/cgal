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


#include <CGAL/predicates/Svd_are_same_points_C2.h>


#ifndef CGAL_SVD_ARE_SAME_SEGMENTS_C2_H
#define CGAL_SVD_ARE_SAME_SEGMENTS_C2_H

CGAL_BEGIN_NAMESPACE

template<class K>
class Svd_are_same_segments_C2
{
private:
  typedef typename K::Site_2          Site_2;
  typedef Svd_are_same_points_C2<K>   Are_same_points_2;

private:
  Are_same_points_2  same_points;

public:
  bool operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() );

    return
      ( same_points(p.source_site(), q.source_site()) &&
        same_points(p.target_site(), q.target_site()) ) ||
      ( same_points(p.source_site(), q.target_site()) &&
        same_points(p.target_site(), q.source_site()) );
  }
};



CGAL_END_NAMESPACE

#endif // CGAL_SVD_ARE_SAME_SEGMENTS_C2_H
