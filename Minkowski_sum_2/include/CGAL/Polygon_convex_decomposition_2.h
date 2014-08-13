// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_POLYGON_CONVEX_DECOMPOSITION_H
#define CGAL_POLYGON_CONVEX_DECOMPOSITION_H

#include <CGAL/Minkowski_sum_2/Decomposition_strategy_adapter.h>
#include <vector>

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

\anchor mink_refopt_decomp 

The `Optimal_convex_decomposition_2` class provides an implementation of Greene's 
dynamic programming algorithm for optimal decomposition of a 
polygon into convex sub-polygons \cgalCite{g-dpcp-83}. Note that 
this algorithm requires \f$ O(n^4)\f$ time and \f$ O(n^3)\f$ space in 
the worst case, where \f$ n\f$ is the size of the input polygon. 


\tparam Kernel must be a geometric kernel that can be used for the polygon.
\tparam Container must be a container that can be used for the polygon. 
It is by default `std::vector<typename Kernel::Point_2>`. 

\cgalModels `PolygonConvexDecomposition_2`

\sa `CGAL::optimal_convex_partition_2()` 

*/
template <class Kernel_, 
          class Container_ = std::vector<typename Kernel_::Point_2> >
class Optimal_convex_decomposition_2 :
  public Polygon_decomposition_strategy_adapter<Kernel_, Container_,
						Tag_optimal_convex_parition>
{
public:

  typedef Kernel_                                  Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>      Polygon_2;
  typedef typename Kernel::Point_2                 Point_2;

};

/*!
\ingroup PkgMinkowskiSum2Classes

\anchor mink_refHM_decomp 

The `Hertel_Mehlhorn_convex_decomposition_2` class implements the approximation algorithm of Hertel 
and Mehlhorn for decomposing a polygon into convex 
sub-polygons \cgalCite{hm-ftsp-83}. This algorithm constructs a 
triangulation of the input polygon and proceeds by removing 
unnecessary triangulation edges. Given the triangulation, the 
algorithm requires \f$ O(n)\f$ time and space to construct a convex 
decomposition (where \f$ n\f$ is the size of the input polygon), whose 
size is guaranteed to be no more than four times the size of the 
optimal decomposition. 

\tparam Kernel must be a geometric kernel that can be used for the polygon.
\tparam Container must be a container that can be used for the polygon. 
It is by default `std::vector<typename Kernel::Point_2>`. 

\cgalModels `PolygonConvexDecomposition_2`

\sa `CGAL::approx_convex_partition_2()` 

*/
template <class Kernel_, 
          class Container_ = std::vector<typename Kernel_::Point_2> >
class Hertel_Mehlhorn_convex_decomposition_2 :
  public Polygon_decomposition_strategy_adapter<Kernel_, Container_,
						Tag_approx_convex_parition>
{
public:

  typedef Kernel_                                  Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>      Polygon_2;
  typedef typename Kernel::Point_2                 Point_2;

};

/*!
\ingroup PkgMinkowskiSum2

\anchor mink_refGreene_decomp 

The `Greene_convex_decomposition_2` class implements the approximation algorithm of 
Greene for the decomposition of an input polygon into convex 
sub-polygons \cgalCite{g-dpcp-83}. This algorithm takes \f$ O(n \log n)\f$ 
time and \f$ O(n)\f$ space, where \f$ n\f$ is the size of the input polygon, 
and outputs a decomposition whose size is guaranteed to be no more 
than four times the size of the optimal decomposition. 

\tparam Kernel must be a geometric kernel that can be used for the polygon.
\tparam Container must be a container that can be used for the polygon. 
It is by default `std::vector<typename Kernel::Point_2>`. 

\cgalModels `PolygonConvexDecomposition_2`

\sa `CGAL::greene_approx_convex_partition_2()` 

*/
template <class Kernel_, 
          class Container_ = std::vector<typename Kernel_::Point_2> >
class Greene_convex_decomposition_2 :
  public Polygon_decomposition_strategy_adapter<Kernel_, Container_,
						Tag_Greene_convex_parition>
{
public:

  typedef Kernel_                                  Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>      Polygon_2;
  typedef typename Kernel::Point_2                 Point_2;

};

} //namespace CGAL

#endif
