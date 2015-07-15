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
// Author(s) : Ron Wein   <wein_r@yahoo.com>
//             Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_MINKOWSKI_SUM_2_H
#define CGAL_MINKOWSKI_SUM_2_H

#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <CGAL/Minkowski_sum_2/Hole_filter_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_by_reduced_convolution_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_conv_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_decomp_2.h>
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
 * Compute the Minkowski sum of two simple polygons by decomposing each
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
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2, decomposition_strategy, traits);
}

/*!
 * Compute the Minkowski sum of two simple polygons by decomposing each
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
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy,
                const typename
                Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                                 Container_>::Traits_2& traits)
{
  typedef Container_                            Container;
  typedef DecompositionStrategy_                Decomposition_strategy;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy, Container>
    mink_sum(decomposition_strategy, traits);
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
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2, decomposition_strategy, traits);
}

/*!
 * Compute the Minkowski sum of two polygon with holes by decomposing each
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
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy,
                const typename
                Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                                 Container_>::Traits_2& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef DecompositionStrategy_                Decomposition_strategy;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy, Container>
    mink_sum(decomposition_strategy, traits);
  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn1;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn2;
  hole_filter(pgn1, pgn2, filtered_pgn1);
  hole_filter(pgn2, pgn1, filtered_pgn2);
  return mink_sum(filtered_pgn1, filtered_pgn2);
}

/*!
 * Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The first polygon.
 * \param[in] pgn2 The second polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2, decomposition_strategy, traits);
}

/*!
 * Compute the Minkowski sum of a simple polygon and a polygon with holes
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
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_with_holes_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy,
                const typename
                Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                                 Container_>::Traits_2& traits)
{
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef DecompositionStrategy_                Decomposition_strategy;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy, Container>
    mink_sum(decomposition_strategy, traits);
  Hole_filter_2<Kernel, Container> hole_filter;
  Polygon_with_holes_2<Kernel,Container> filtered_pgn2;
  hole_filter(pgn2, pgn1, filtered_pgn2);
  return mink_sum(pgn1, filtered_pgn2);
}

/*!
 * Compute the Minkowski sum of a simple polygon and a polygon with holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param[in] pgn1 The second polygon.
 * \param[in] pgn2 The first polygon.
 * \param[in] decomposition_strategy A functor for decomposing polygons.
 * \return The resulting polygon with holes, representing the sum.
 */
template <typename Kernel_, typename Container_,
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy)
{
  typename Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                            Container_>::Traits_2 traits;
  return minkowski_sum_2(pgn1, pgn2, decomposition_strategy, traits);
}

/*!
 * Compute the Minkowski sum of a simple polygon and a polygon with holes
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
          typename DecompositionStrategy_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2,
                const DecompositionStrategy_& decomposition_strategy,
                const typename
                Minkowski_sum_by_decomposition_2<DecompositionStrategy_,
                                                 Container_>::Traits_2& traits)
{ return minkowski_sum_2(pgn2, pgn1, decomposition_strategy, traits); }

} //namespace CGAL

#endif
