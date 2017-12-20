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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Ron Wein   <wein_r@yahoo.com>
//             Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_MINKOWSKI_SUM_2_H
#define CGAL_MINKOWSKI_SUM_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Minkowski_sum_2/Hole_filter_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_by_reduced_convolution_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_conv_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_decomp_2.h>
#include <CGAL/Polygon_nop_decomposition_2.h>
#include <CGAL/Gps_segment_traits_2.h>

#include <list>

namespace CGAL {

/*!
 * Computes the Minkowski sum \f$ P \oplus Q\f$ of two given polygons.
 * The function computes the reduced convolution of the two polygons and
 * extracts those loops of the convolution which are part of the Minkowsi
 * sum. This method works very efficiently, regardless of whether `P` and
 * `Q` are convex or non-convex.
 * Note that as the input polygons may not be convex, their Minkowski
 * sum may not be a simple polygon. The result is therefore represented
 * as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \pre Both `P` and `Q` are simple, counterclockwise-oriented polygons.
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_reduced_convolution_2(const Polygon_2<Kernel_,
                                       Container_>& pgn1,
                                       const Polygon_2<Kernel_,
                                       Container_>& pgn2)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Minkowski_sum_by_reduced_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel, Container>                              sum_bound;
  std::list<Polygon_2<Kernel, Container> >                  sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum(pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else mink_sum(pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container>(sum_bound,
                                                 sum_holes.begin(),
                                                 sum_holes.end()));
}

/*!
 * Computes the Minkowski sum \f$ P \oplus Q\f$ of two given polygons with
 * holes.
 * The function computes the reduced convolution of the two polygons and
 * extracts those loops of the convolution which are part of the Minkowsi
 * sum. This method works very efficiently, regardless of whether `P` and
 * `Q` are convex or non-convex.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_reduced_convolution_2
(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
 const Polygon_with_holes_2<Kernel_, Container_>& pgn2)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Hole_filter_2<Kernel, Container> hole_filter;

  Polygon_with_holes_2<Kernel, Container> filtered_pgn1;
  Polygon_with_holes_2<Kernel, Container> filtered_pgn2;

  hole_filter(pgn1, pgn2, filtered_pgn1);
  hole_filter(pgn2, pgn1, filtered_pgn2);

  Minkowski_sum_by_reduced_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel, Container>                              sum_bound;
  std::list<Polygon_2<Kernel, Container> >                  sum_holes;

  mink_sum(filtered_pgn1, filtered_pgn2, sum_bound,
           std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container>(sum_bound,
                                                 sum_holes.begin(),
                                                 sum_holes.end()));
}

/*!
 * Computes the Minkowski sum \f$ P \oplus Q\f$ of a simple polygon and a
 * polygon with holes.
 * The function computes the reduced convolution of the two polygons and
 * extracts those loops of the convolution which are part of the Minkowsi
 * sum. This method works very efficiently, regardless of whether `P` and
 * `Q` are convex or non-convex.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The simple polygon.
 * \param[in] pgn2 The polygon with holes.
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_reduced_convolution_2
(const Polygon_2<Kernel_, Container_>& pgn1,
 const Polygon_with_holes_2<Kernel_, Container_>& pgn2)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel, Container> filtered_pgn2;
  hole_filter(pgn2, pgn1, filtered_pgn2);
  Minkowski_sum_by_reduced_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel, Container>                              sum_bound;
  std::list<Polygon_2<Kernel, Container> >                  sum_holes;
  mink_sum(pgn1, filtered_pgn2, sum_bound, std::back_inserter(sum_holes));
  return (Polygon_with_holes_2<Kernel,Container>(sum_bound,
                                                 sum_holes.begin(),
                                                 sum_holes.end()));
}

/*!
 * Computes the Minkowski sum \f$ P \oplus Q\f$ of a simple polygon and a
 * polygon with holes.
 * The function computes the reduced convolution of the two polygons and
 * extracts those loops of the convolution which are part of the Minkowsi
 * sum. This method works very efficiently, regardless of whether `P` and
 * `Q` are convex or non-convex.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The polygon with holes.
 * \param[in] pgn2 The simple polygon.
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_reduced_convolution_2
(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
 const Polygon_2<Kernel_, Container_>& pgn2)
{ return minkowski_sum_by_reduced_convolution_2(pgn2, pgn1); }

/*!
 * Compute the Minkowski sum of two simple polygons using the (full)
 * convolution method.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_full_convolution_2(const Polygon_2<Kernel_, Container_>& pgn1,
                                    const Polygon_2<Kernel_, Container_>& pgn2)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Minkowski_sum_by_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel, Container>                      sum_bound;
  std::list<Polygon_2<Kernel, Container> >          sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum(pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else mink_sum(pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));
  return (Polygon_with_holes_2<Kernel, Container>(sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

/*!
 * Compute the Minkowski sum of two simple polygons using the convolution
 * method. This function defaults to calling the reduced convolution method,
 * as it is more efficient in most cases.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \return The resulting polygon with holes, representing the sum.
 *
 * \sa `CGAL::minkowski_sum_by_reduced_convolution_2()`
 * \sa `CGAL::minkowski_sum_by_full_convolution_2()`
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2)
{ return minkowski_sum_by_reduced_convolution_2(pgn1, pgn2); }

/*!
 * Compute the Minkowski sum of two polygons with holes using the convolution
 * method. This function defaults to calling the reduced convolution method,
 * as it is more efficient in most cases.
 * Note that the result may not be a simple polygon. The result is therefore
 * represented as a polygon with holes.
 * \param[in] pgn1 The first polygon with holes.
 * \param[in] pgn2 The second polygon with holes.
 * \return The resulting polygon with holes, representing the sum.
 *
 * \sa `CGAL::minkowski_sum_by_reduced_convolution_2()`
 * \sa `CGAL::minkowski_sum_by_full_convolution_2()`
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2)
{ return minkowski_sum_by_reduced_convolution_2(pgn1, pgn2); }

/*!
 * Compute the Minkowski sum of a simple polygon and a polygon with holes
 * using the convolution method. This function defaults to calling the reduced
 * convolution method, as it is more efficient in most cases.
 * Note that the result may not be a simple polygon. The result is therefore
 * represented as a polygon with holes.
 * \param[in] pgn1 The simple polygon.
 * \param[in] pgn2 The polygon with holes.
 * \return The resulting polygon with holes, representing the sum.
 *
 * \sa `CGAL::minkowski_sum_by_reduced_convolution_2()`
 * \sa `CGAL::minkowski_sum_by_full_convolution_2()`
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2)
{ return minkowski_sum_by_reduced_convolution_2(pgn1, pgn2); }

/*!
 * Compute the Minkowski sum of a simple polygon and a polygon with holes
 * using the convolution method. This function defaults to calling the reduced
 * convolution method, as it is more efficient in most cases.
 * Note that the result may not be a simple polygon. The result is therefore
 * represented as a polygon with holes.
 * \param[in] pgn1 The polygon with holes.
 * \param[in] pgn2 The simple polygon.
 * \return The resulting polygon with holes, representing the sum.
 *
 * \sa `CGAL::minkowski_sum_by_reduced_convolution_2()`
 * \sa `CGAL::minkowski_sum_by_full_convolution_2()`
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2)
{ return minkowski_sum_by_reduced_convolution_2(pgn1, pgn2); }

/*!
 * \defgroup Minkowski sum by decomposition
 * @{
 */

/*! Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy2_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy2,
                         traits);
}

/*! Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 *
 * The type of the last argument, namely,
 *   Gps_segment_traits_2<Kernel_, Container_>
 * and the type
 *   const typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
 *                                                   DecompositionStrategy2_,
 *                                                   Container_>::Traits_2>
 * are exchangable except for in one case, where there is an ambiguity.
 * Thus, we use the former, even though it is less generic, as change to the
 * traits type in Minkowski_sum_by_decomposition_2 would require a similar
 * change here.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Container_                            Container;
  typedef DecompositionStrategy1_               Decomposition_strategy1;
  typedef DecompositionStrategy2_               Decomposition_strategy2;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy1,
                                   Decomposition_strategy2, Container>
    mink_sum(decomposition_strategy1, decomposition_strategy2, traits);
  return mink_sum(pgn1, pgn2);
}

/*!
 * Compute the Minkowski sum of two polygon with holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy2_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy2,
                         traits);
}

/*! Compute the Minkowski sum of two polygon with holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef DecompositionStrategy1_               Decomposition_strategy1;
  typedef DecompositionStrategy2_               Decomposition_strategy2;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy1,
                                   Decomposition_strategy2, Container>
    mink_sum(decomposition_strategy1, decomposition_strategy2, traits);
  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn1;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn2;
  hole_filter(pgn1, pgn2, filtered_pgn1);
  hole_filter(pgn2, pgn1, filtered_pgn2);
  return mink_sum(filtered_pgn1, filtered_pgn2);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy2_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy2,
                         traits);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The simple polygon.
 * \param[in] pgn2 The polygon with holes.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef DecompositionStrategy1_               Decomposition_strategy1;
  typedef DecompositionStrategy2_               Decomposition_strategy2;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy1,
                                   Decomposition_strategy2, Container>
    mink_sum(decomposition_strategy1, decomposition_strategy2, traits);
  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel, Container> filtered_pgn2;
  hole_filter(pgn2, pgn1, filtered_pgn2);
  return mink_sum(pgn1, filtered_pgn2);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The second polygon.
 * \param[in] pgn2 The first polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy2_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy2,
                         traits);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The polygon with holes.
 * \param[in] pgn2 The simple polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_, typename DecompositionStrategy2_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const DecompositionStrategy2_& decomposition_strategy2,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  return minkowski_sum_2(pgn2, pgn1,
                         decomposition_strategy2, decomposition_strategy1,
                         traits);
}

/*! Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy1_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy1,
                         traits);
}

/*! Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Container_                            Container;
  typedef DecompositionStrategy1_               Decomposition_strategy1;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy1,
                                   Decomposition_strategy1, Container>
    mink_sum(decomposition_strategy1, decomposition_strategy1, traits);
  return mink_sum(pgn1, pgn2);
}

/*!
 * Compute the Minkowski sum of two polygon with holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy1_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy1,
                         traits);
}

/*! Compute the Minkowski sum of two polygon with holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef DecompositionStrategy1_               Decomposition_strategy1;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy1,
                                   Decomposition_strategy1, Container>
    mink_sum(decomposition_strategy1, decomposition_strategy1, traits);
  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn1;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn2;
  hole_filter(pgn1, pgn2, filtered_pgn1);
  hole_filter(pgn2, pgn1, filtered_pgn2);
  return mink_sum(filtered_pgn1, filtered_pgn2);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy1_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy1,
                         traits);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The simple polygon.
 * \param[in] pgn2 The polygon with holes.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef DecompositionStrategy1_               Decomposition_strategy1;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy1,
                                   Decomposition_strategy1, Container>
    mink_sum(decomposition_strategy1, decomposition_strategy1, traits);
  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn2;
  hole_filter(pgn2, pgn1, filtered_pgn2);
  return mink_sum(pgn1, filtered_pgn2);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The second polygon.
 * \param[in] pgn2 The first polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy1_,
                                            DecompositionStrategy1_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2,
                         decomposition_strategy1, decomposition_strategy1,
                         traits);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The polygon with holes.
 * \param[in] pgn2 The simple polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy1_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy1_& decomposition_strategy1,
                const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  return minkowski_sum_2(pgn2, pgn1,
                         decomposition_strategy1, decomposition_strategy1,
                         traits);
}

/*! Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_, typename Decomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2(const Polygon_2<Kernel_, Container_>& pgn1,
                                 const Polygon_2<Kernel_, Container_>& pgn2,
                                 const Decomposition_& decomp)
{
  typename Minkowski_sum_by_decomposition_2<Decomposition_, Decomposition_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_by_decomposition_2(pgn1, pgn2, decomp, traits);
}

/*! Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_, typename Decomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2
(const Polygon_2<Kernel_, Container_>& pgn1,
 const Polygon_2<Kernel_, Container_>& pgn2,
 const Decomposition_& decomp,
 const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef Decomposition_                        Decomposition;
  typedef Polygon_nop_decomposition_2<Kernel>   Nop_decomposition;

  if (pgn1.is_convex()) {
    Nop_decomposition decomp_nop;
    if (pgn2.is_convex()) {
      Minkowski_sum_by_decomposition_2<Nop_decomposition, Nop_decomposition,
                                       Container>
        mink_sum(decomp_nop, decomp_nop, traits);
      return mink_sum(pgn1, pgn2);
    }
    Minkowski_sum_by_decomposition_2<Nop_decomposition, Decomposition,
                                     Container>
      mink_sum(decomp_nop, decomp, traits);
    return mink_sum(pgn1, pgn2);
  }

  if (pgn2.is_convex()) {
    Nop_decomposition decomp_nop;
    Minkowski_sum_by_decomposition_2<Decomposition, Nop_decomposition,
                                     Container>
      mink_sum(decomp, decomp_nop, traits);
    return mink_sum(pgn1, pgn2);
  }

  Minkowski_sum_by_decomposition_2<Decomposition, Decomposition, Container>
    mink_sum(decomp, decomp, traits);
  return mink_sum(pgn1, pgn2);
}

/*! Compute the Minkowski sum of two polygon with holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename NoHolesDecomposition_, typename WithHolesDecomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2
(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
 const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
 const NoHolesDecomposition_& decomp_no_holes,
 const WithHolesDecomposition_& decomp_with_holes)
{
  typename Minkowski_sum_by_decomposition_2<NoHolesDecomposition_,
                                            WithHolesDecomposition_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_by_decomposition_2(pgn1, pgn2,
                                          decomp_no_holes, decomp_with_holes,
                                          traits);
}

/*! Compute the Minkowski sum of two polygon with holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename NoHolesDecomposition_, typename WithHolesDecomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2
(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
 const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
 const NoHolesDecomposition_& decomp_no_holes,
 const WithHolesDecomposition_& decomp_with_holes,
 const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef NoHolesDecomposition_                 No_holes_decomposition;
  typedef WithHolesDecomposition_               With_holes_decomposition;
  typedef Polygon_nop_decomposition_2<Kernel>   Nop_decomposition;

  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel, Container> filtered_pgn1;
  Polygon_with_holes_2<Kernel, Container> filtered_pgn2;
  hole_filter(pgn1, pgn2, filtered_pgn1);
  hole_filter(pgn2, pgn1, filtered_pgn2);

  if (0 == filtered_pgn1.number_of_holes()) {
    const Polygon_2<Kernel, Container>& pnh1 = filtered_pgn1.outer_boundary();
    if (pnh1.is_convex()) {
      // pnh1 is convex
      Nop_decomposition decomp_nop;
      if (0 == filtered_pgn2.number_of_holes()) {
        const Polygon_2<Kernel, Container>& pnh2 =
          filtered_pgn2.outer_boundary();
        if (pnh2.is_convex()) {
          // pnh2 is convex
          Minkowski_sum_by_decomposition_2<Nop_decomposition, Nop_decomposition,
                                           Container>
            mink_sum(decomp_nop, decomp_nop, traits);
          return mink_sum(pnh1, pnh2);
        }
        // pnh2 is concave
        Minkowski_sum_by_decomposition_2<Nop_decomposition,
                                         No_holes_decomposition,
                                         Container>
          mink_sum(decomp_nop, decomp_no_holes, traits);
        return mink_sum(pnh1, pnh2);
      }
      // pnh2 has holes
      Minkowski_sum_by_decomposition_2<Nop_decomposition,
                                       With_holes_decomposition,
                                       Container>
        mink_sum(decomp_nop, decomp_with_holes, traits);
      return mink_sum(pnh1, filtered_pgn2);
    }
    // pnh1 is concave
    if (0 == filtered_pgn2.number_of_holes()) {
      const Polygon_2<Kernel, Container>& pnh2 =
        filtered_pgn2.outer_boundary();
      if (pnh2.is_convex()) {
        // pnh2 is convex
        Nop_decomposition decomp_nop;
        Minkowski_sum_by_decomposition_2<No_holes_decomposition,
                                         Nop_decomposition,
                                         Container>
          mink_sum(decomp_no_holes, decomp_nop, traits);
        return mink_sum(pnh1, pnh2);
      }
      // pnh2 is concave
      Minkowski_sum_by_decomposition_2<No_holes_decomposition,
                                       No_holes_decomposition,
                                       Container>
        mink_sum(decomp_no_holes, decomp_no_holes, traits);
      return mink_sum(pnh1, pnh2);
    }
    // pnh2 has holes
    Minkowski_sum_by_decomposition_2<No_holes_decomposition,
                                     With_holes_decomposition,
                                     Container>
      mink_sum(decomp_no_holes, decomp_with_holes, traits);
    return mink_sum(pnh1, filtered_pgn2);
  }

  // filtered_pgn1 has holes
  if (0 == filtered_pgn2.number_of_holes()) {
    const Polygon_2<Kernel, Container>& pnh2 =
      filtered_pgn2.outer_boundary();
    if (pnh2.is_convex()) {
      // pnh2 is convex
      Nop_decomposition decomp_nop;
      Minkowski_sum_by_decomposition_2<Nop_decomposition,
                                       With_holes_decomposition,
                                       Container>
        mink_sum(decomp_nop, decomp_with_holes, traits);
      return mink_sum(pnh2, filtered_pgn1);
    }
    // pnh2 is concave
    Minkowski_sum_by_decomposition_2<No_holes_decomposition,
                                     With_holes_decomposition,
                                     Container>
      mink_sum(decomp_no_holes, decomp_with_holes, traits);
    return mink_sum(pnh2, filtered_pgn1);
  }
  // pnh2 has holes
  Minkowski_sum_by_decomposition_2<With_holes_decomposition,
                                   With_holes_decomposition,
                                   Container>
    mink_sum(decomp_with_holes, decomp_with_holes, traits);
  return mink_sum(filtered_pgn1, filtered_pgn2);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename NoHolesDecomposition_, typename WithHolesDecomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2
(const Polygon_2<Kernel_, Container_>& pgn1,
 const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
 const NoHolesDecomposition_& decomp_no_holes,
 const WithHolesDecomposition_& decomp_with_holes)
{
  typename Minkowski_sum_by_decomposition_2<NoHolesDecomposition_,
                                            WithHolesDecomposition_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_by_decomposition_2(pgn1, pgn2,
                                          decomp_no_holes, decomp_with_holes,
                                          traits);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The simple polygon.
 * \param[in] pgn2 The polygon with holes.
 * \param[in] decomposition A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename NoHolesDecomposition_, typename WithHolesDecomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2
(const Polygon_2<Kernel_, Container_>& pgn1,
 const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
 const NoHolesDecomposition_& decomp_no_holes,
 const WithHolesDecomposition_& decomp_with_holes,
 const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef NoHolesDecomposition_                 No_holes_decomposition;
  typedef WithHolesDecomposition_               With_holes_decomposition;
  typedef Polygon_nop_decomposition_2<Kernel>   Nop_decomposition;

  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn2;
  hole_filter(pgn2, pgn1, filtered_pgn2);

  if (pgn1.is_convex()) {
    Nop_decomposition decomp_nop;
    if (0 == filtered_pgn2.number_of_holes()) {
      const Polygon_2<Kernel, Container>& pnh2 = filtered_pgn2.outer_boundary();
      if (pnh2.is_convex()) {
        Minkowski_sum_by_decomposition_2<Nop_decomposition, Nop_decomposition,
                                         Container>
          mink_sum(decomp_nop, decomp_nop, traits);
        return mink_sum(pgn1, pnh2);
      }
      // pnh2 is concave
      Minkowski_sum_by_decomposition_2<Nop_decomposition,
                                       No_holes_decomposition,
                                       Container>
        mink_sum(decomp_nop, decomp_no_holes, traits);
      return mink_sum(pgn1, pnh2);
    }
    // pnh2 has holes
    Minkowski_sum_by_decomposition_2<Nop_decomposition,
                                     With_holes_decomposition,
                                     Container>
      mink_sum(decomp_nop, decomp_with_holes, traits);
    return mink_sum(pgn1, filtered_pgn2);
  }
  // pgn1 is concave
  if (0 == filtered_pgn2.number_of_holes()) {
    const Polygon_2<Kernel, Container>& pnh2 = filtered_pgn2.outer_boundary();
    if (pnh2.is_convex()) {
      // pnh2 is convex
      Nop_decomposition decomp_nop;
      Minkowski_sum_by_decomposition_2<No_holes_decomposition,
                                       Nop_decomposition,
                                       Container>
        mink_sum(decomp_no_holes, decomp_nop, traits);
      return mink_sum(pgn1, pnh2);
    }
    // pnh2 is concave
    Minkowski_sum_by_decomposition_2<No_holes_decomposition,
                                     No_holes_decomposition,
                                     Container>
      mink_sum(decomp_no_holes, decomp_no_holes, traits);
    return mink_sum(pgn1, pnh2);
  }
  // pnh2 has holes
  Minkowski_sum_by_decomposition_2<No_holes_decomposition,
                                   With_holes_decomposition,
                                   Container>
    mink_sum(decomp_no_holes, decomp_with_holes, traits);
  return mink_sum(pgn1, filtered_pgn2);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The second polygon.
 * \param[in] pgn2 The first polygon.
 * \param[in] decomposition A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename NoHolesDecomposition_, typename WithHolesDecomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2
(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
 const Polygon_2<Kernel_, Container_>& pgn2,
 const NoHolesDecomposition_& decomp_no_holes,
 const WithHolesDecomposition_& decomp_with_holes)
{
  typename Minkowski_sum_by_decomposition_2<NoHolesDecomposition_,
                                            WithHolesDecomposition_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_by_decomposition_2(pgn1, pgn2,
                                          decomp_no_holes, decomp_with_holes,
                                          traits);
}

/*! Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The polygon with holes.
 * \param[in] pgn2 The simple polygon.
 * \param[in] decomposition A functor for decomposing polygons.
 * \param[in] traits The traits.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename NoHoleDecomposition_, typename WithHolesDecomposition_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_by_decomposition_2
(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
 const Polygon_2<Kernel_, Container_>& pgn2,
 const NoHoleDecomposition_& decomp_no_holes,
 const WithHolesDecomposition_& decomp_with_holes,
 const Gps_segment_traits_2<Kernel_, Container_>& traits)
{
  return minkowski_sum_by_decomposition_2(pgn2, pgn1,
                                          decomp_no_holes, decomp_with_holes,
                                          traits);
}

/*!@}*/

} //namespace CGAL

#endif
