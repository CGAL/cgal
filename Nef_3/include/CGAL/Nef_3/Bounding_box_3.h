// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// Author(s)     : Miguel Granados <granados@mpi-sb.mpg.de>

#ifndef CGAL_BOUNDING_BOX_3_H
#define CGAL_BOUNDING_BOX_3_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/box_limits.h>

CGAL_BEGIN_NAMESPACE

template <typename Traits>
class Bounding_box_3 : 
public Box_intersection_d::Box_d< double, 3> {

  typedef Box_intersection_d::Box_d< double, 3>  Base;

  typedef typename Traits::Point_3             Point_3;

public:
  Bounding_box_3() : Base() {}
    
  void extend( const Point_3& p) {
    std::pair<double, double> q[3];
    q[0] = CGAL::to_interval( p.x() );
    q[1] = CGAL::to_interval( p.y() );
    q[2] = CGAL::to_interval( p.z() );
    Box_intersection_d::Box_d< double, 3 >::extend(q);
  }	
};

CGAL_END_NAMESPACE
#endif // CGAL_BOUNDING_BOX_3_H
