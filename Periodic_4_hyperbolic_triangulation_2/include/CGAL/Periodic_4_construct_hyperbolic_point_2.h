// Copyright (c) 1999-2004,2006-2009,2014-2015   INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pellé <Aymeric.Pelle@sophia.inria.fr>
#ifndef CGAL_PERIODIC_4_CONSTRUCT_HYPERBOLIC_POINT_2_H
#define CGAL_PERIODIC_4_CONSTRUCT_HYPERBOLIC_POINT_2_H

#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_group.h>

namespace CGAL
{
template < typename K, typename Construct_point_2_base>
class Periodic_4_construct_hyperbolic_point_2 : public Construct_point_2_base
{
  typedef K Kernel;

public:
  typedef typename Kernel::Point_2            Point;
  typedef HyperbolicOctagonGroup              Offset;
  typedef typename Kernel::Circle_2           Circle_2;

  typedef Point                               result_type;

  Periodic_4_construct_hyperbolic_point_2(const Circle_2 & dom) : _dom(dom) { }

  Point operator() ( const Point& p, const Offset& o ) const {
    pair<double, double> pt = o.apply(to_double(p.x()), to_double(p.y()));
    Point r = Point(pt.first, pt.second);

    // Do we need this? I (iiordano) think it's a good idea...
    assert( _dom.has_on_bounded_side(r) );
    
    return r;
  }

private:
  Circle_2 _dom;
};
}

#endif
