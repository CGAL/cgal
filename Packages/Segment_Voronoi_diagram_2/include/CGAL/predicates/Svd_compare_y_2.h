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


#ifndef CGAL_SVD_COMPARE_Y_2_H
#define CGAL_SVD_COMPARE_Y_2_H


CGAL_BEGIN_NAMESPACE



//-----------------------------------------------------------------------
//                           compare y
//-----------------------------------------------------------------------

template< class K >
class Svd_compare_y_2
{
public:
  typedef typename K::Site_2         Site_2;
  typedef typename K::Point_2        Point_2;
  typedef Comparison_result          result_type;
  typedef Arity_tag<2>               Arity;

private:
  typedef typename K::Compare_y_2  compare_y_2;

public:

  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );
    return compare_y_2()( p.point(), q.point() );
  }
};


CGAL_END_NAMESPACE

#endif // CGAL_SVD_COMPARE_Y_2_H
