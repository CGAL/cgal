// Copyright (c) 2009 INRIA Sophia-Antipolis (France),
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
// Author(s)    : Samuel Hornus

/* RANDOM DESIGN IDEAS:
- Use a policy tag to choose for incremental with inserts only or
  incremental with  removals and inserts.
  In the first case: use Pure_complex for storage.
  In the second case: use Delaunay !
    In this second case, we must keeps the points that are inserted in the hull,
    as they may become part of the boundary later on, when some points are removed.
- Constructor with range argument uses quickhull.
*/

#ifndef CGAL_CONVEX_HULL_H
#define CGAL_CONVEX_HULL_H

namespace CGAL {

template <  class CHTraits, class PCDS_ = Default >
class Convex_hull
{
    typedef typename Ambient_dimension<typename PCTraits::Point_d>::type
                                                    Ambient_dimension_;
    typedef typename Default::Get<PCDS_, Pure_complex_data_structure
                    <   Ambient_dimension_,
                        Pure_complex_vertex<PCTraits>,
                        Pure_complex_simplex<PCTraits> >
                        >::type                     PCDS;
    typedef Pure_complex<PCTraits, PCDS_>           Self;

    typedef typename PCTraits::Coaffine_orientation_d
                                                    Coaffine_orientation_d;
    typedef typename PCTraits::Orientation_d        Orientation_d;

public:
};

} //namespace CGAL

#endif // CGAL_CONVEX_HULL_H