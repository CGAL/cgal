// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_MINKOWSKI_SUM_3_H
#define CGAL_MINKOWSKI_SUM_3_H

#include <CGAL/license/Minkowski_sum_3.h>


#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Minkowski_sum_3/bipartite_nary_union_sorted_combined.h> 
#include <CGAL/Is_extended_kernel.h>

/// \file minkowski_sum_3.h

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum3

The function `minkowski_sum_3()` computes the Minkowski sum of two 
given 3D Nef polyhedra \f$ N0\f$ and \f$ N1\f$. Note that the function runs in 
\f$ O(n^3m^3)\f$ time in the worst case, where \f$ n\f$ and 
\f$ m\f$ are the complexities of the two input polyhedra (the complexity of 
a `Nef_polyhedron_3` is the sum of its `Vertices`, 
`Halfedges` and `SHalfedges`). 

An input polyhedron may consist of: 
<OL> 
<LI>singular vertices 
<LI>singular edges 
<LI>singular convex facets without holes 
<LI>surfaces with convex facets that have no holes. 
<LI>three-dimensional features, whose coplanar facets have 
common selection marks (this includes open and closed solids) 
</OL> 

Taking a different viewpoint, the implementation is restricted as 
follows: 
<OL> 
<LI>The input polyhedra must be bounded (selected outer volume is ignored). 
<LI>All sets of coplanar facets of a full-dimensional 
feature must have the same selection mark (in case of different 
selection marks, unselected is assumed). 
<LI>All facets of lower-dimensional features need to be convex and 
must not have holes (non-convex facets and holes are ignored). 
</OL> 

\post If either of the input polyhedra is non-convex, it is modified during the computation, i.e., it is decomposed into convex pieces. 

\sa `CGAL::Nef_polyhedron_3<Traits>` 
\sa \link CGAL::convex_decomposition_3 `CGAL::convex_decomposition_3()`\endlink

*/
template<typename Nef_polyhedron_3>
Nef_polyhedron_3 
minkowski_sum_3(Nef_polyhedron_3& N0, Nef_polyhedron_3& N1) 
{
  typedef typename Nef_polyhedron_3::Kernel Kernel;
  typedef typename Is_extended_kernel<Kernel>::value_type Is_extended_kernel;
  if(check_tag(Is_extended_kernel())) {
    std::cerr << "extended kernel is not supported" << std::endl;
    return N0;
  }

  if(N0.volumes_begin()->mark()) {
    std::cerr << "first parameter is an infinite point set" << std::endl;
    return N0;
  }
   
  if(N1.volumes_begin()->mark()) {
    std::cerr << "second parameter is an infinite point set" << std::endl;
    return N1;
  }


  CGAL::convex_decomposition_3<Nef_polyhedron_3>(N0);
  CGAL::convex_decomposition_3<Nef_polyhedron_3>(N1);
  CGAL_assertion(N0.is_valid());
  CGAL_assertion(N1.is_valid());

  Nef_polyhedron_3 result =
    CGAL::bipartite_nary_union_sorted_combined(N0, N1);
  CGAL_assertion(result.is_valid());
  return result;
}

} //namespace CGAL
#endif // CGAL_MINKOWSKI_SUM_3_H
