// Copyright (c) 2011 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_COREFINEMENT_PREDICATES_H
#define CGAL_INTERNAL_COREFINEMENT_PREDICATES_H

#include <CGAL/license/Polygon_mesh_processing.h>


namespace CGAL{

namespace Corefinement{

namespace OOP{
//Considering the plane with normal vector [O_prime,O] and containing O.
//We define the counterclockwise order around O when looking from
//the side of the plane into which the vector [O_prime,O] is pointing.
//We consider the portion of the plane defined by rotating a ray starting at O
//from the planar projection of P1 to the planar projection of P2 in
//counterclockwise order.
//The predicates indicates whether the planar projection of point Q lies in this
//portion of the plane.
//Preconditions:
//  O_prime,O,P1 are not collinear
//  O_prime,O,P2 are not collinear
//  O_prime,O,Q are not collinear
//  O_prime,O,P1,Q are not coplanar or coplanar_orientation(O,O_prime,P1,Q)==NEGATIVE
//  O_prime,O,P2,Q are not coplanar or coplanar_orientation(O,O_prime,P2,Q)==NEGATIVE
template <class Kernel>
bool  sorted_around_edge(
  const typename Kernel::Point_3& O_prime,const typename Kernel::Point_3& O,
  const typename Kernel::Point_3& P1,const typename Kernel::Point_3& P2,
  const typename Kernel::Point_3& Q)
{
  //guarantee to have non-flat triangles
  CGAL_precondition( !collinear(O_prime,O,P1) );
  CGAL_precondition( !collinear(O_prime,O,P2) );
  CGAL_precondition( !collinear(O_prime,O,Q)  );

  //no two triangles are coplanar and on the same side of their common edge
  CGAL_precondition( !coplanar(O_prime,O,P1,Q)
                     || coplanar_orientation(O,O_prime,P1,Q)==NEGATIVE );
  CGAL_precondition( !coplanar(O_prime,O,P2,Q)
                     || coplanar_orientation(O,O_prime,P2,Q)==NEGATIVE );

  Sign s0 = CGAL::sign( determinant(O-O_prime,P1-O,P2-O) );

  if ( s0==ZERO ) {
    //O, O_prime, P1 and P2 are coplanar
    Orientation o=orientation(O_prime,O,P1,Q);
    CGAL_precondition(o!=COPLANAR);
    return o==POSITIVE;
  }

  //O, O_prime, P1 and P2 are not coplanar
  Sign s1 = CGAL::sign( determinant(O-O_prime,P1-O,Q -O) );
  Sign s2 = CGAL::sign( determinant(O-O_prime,Q -O,P2-O) );

  if (s0 == POSITIVE) // the angle P1,O,P2 is smaller that Pi.
    return ( s1 == POSITIVE )
           && ( s2 ==POSITIVE ); //true if the angles P1,O,Q and Q,O,P2 are smaller than Pi
  else
    return ( s1 != NEGATIVE )
           || ( s2 !=
                NEGATIVE ); //true if the angle P1,O,Q or the angle Q,O,P2 is smaller than or equal to Pi
}

template <class PolyhedronPointPMap,class Nodes_vector, class Vertex_handle>
bool sorted_around_edge_filtered( int O_prime_index,
                                  int O_index,
                                  int P1_index,
                                  int P2_index,
                                  int Q_index,
                                  Vertex_handle P1,
                                  Vertex_handle P2,
                                  Vertex_handle Q,
                                  const Nodes_vector& nodes,
                                  PolyhedronPointPMap ppmap)
{
  typename Nodes_vector::Protector p;
  try {
    CGAL_USE(p);
    return sorted_around_edge<typename Nodes_vector::Ikernel>(
             nodes.interval_node(O_prime_index),
             nodes.interval_node(O_index),
             P1_index == -1 ? nodes.to_interval(get(ppmap,P1))
                            : nodes.interval_node(P1_index),
             P2_index == -1 ? nodes.to_interval(get(ppmap,P2))
                            : nodes.interval_node(P2_index),
             Q_index  == -1 ? nodes.to_interval(get(ppmap,Q))
                            : nodes.interval_node(Q_index )
           );
  } catch(Uncertain_conversion_exception&) {
    return sorted_around_edge<typename Nodes_vector::Exact_kernel>(
             nodes.exact_node(O_prime_index),
             nodes.exact_node(O_index),
             P1_index == -1 ? nodes.to_exact(get(ppmap,P1))
                            : nodes.exact_node(P1_index),
             P2_index == -1 ? nodes.to_exact(get(ppmap,P2))
                            : nodes.exact_node(P2_index),
             Q_index  == -1 ? nodes.to_exact(get(ppmap,Q))
                            : nodes.exact_node(Q_index )
           );
  }
}

}

}

} // end of namespace CGAL::Corefinement::OOP

#endif //CGAL_INTERNAL_COREFINEMENT_PREDICATES_H
