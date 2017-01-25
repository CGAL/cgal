// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_PARTITION_H
#define CGAL_PARTITION_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/Partition_2/partition_greene_approx_convex_2.h>
#include <CGAL/Partition_2/partition_optimal_convex_2.h>
#include <CGAL/Partition_2/partition_approx_convex_2.h>
#include <CGAL/Partition_2/partition_y_monotone_2.h>

namespace CGAL {

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator greene_approx_convex_partition_2(InputIterator first, 
                                                InputIterator beyond,
                                                OutputIterator result, 
                                                const Traits& traits)
{
   return partition_greene_approx_convex_2(first, beyond, result, traits);
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator greene_approx_convex_partition_2(InputIterator first, 
                                                InputIterator beyond,
                                                OutputIterator result) 
{
   return partition_greene_approx_convex_2(first, beyond, result);
}

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator optimal_convex_partition_2(InputIterator first, 
                                          InputIterator beyond,
                                          OutputIterator result, 
                                          const Traits& traits)
{
   return partition_optimal_convex_2(first, beyond, result, traits);
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator optimal_convex_partition_2(InputIterator first, 
                                          InputIterator beyond,
                                          OutputIterator result) 
{
   return partition_optimal_convex_2(first, beyond, result);
}

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator approx_convex_partition_2(InputIterator first, 
                                         InputIterator beyond,
                                         OutputIterator result, 
                                         const Traits& traits)
{
   return partition_approx_convex_2(first, beyond, result, traits);
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator approx_convex_partition_2(InputIterator first, 
                                         InputIterator beyond,
                                         OutputIterator result) 
{
   return partition_approx_convex_2(first, beyond, result);
}

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator y_monotone_partition_2(InputIterator first, 
                                      InputIterator beyond,
                                      OutputIterator result, 
                                      const Traits& traits)
{
   return partition_y_monotone_2(first, beyond, result, traits);
}

template <class InputIterator, class OutputIterator>
inline
OutputIterator y_monotone_partition_2(InputIterator first, 
                                      InputIterator beyond,
                                      OutputIterator result)
{
   return partition_y_monotone_2(first, beyond, result);
}

}

#endif // CGAL_PARTITION_H
