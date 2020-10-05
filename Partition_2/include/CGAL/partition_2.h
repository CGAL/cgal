// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
