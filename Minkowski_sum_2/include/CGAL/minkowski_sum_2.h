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

#include <CGAL/Minkowski_sum_2/Minkowski_sum_by_reduced_convolution_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_conv_2.h>
#include <CGAL/Minkowski_sum_2/Minkowski_sum_decomp_2.h>
#include <list>

namespace CGAL {

/*!
 * Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
 * The function computes the reduced convolution of the two polygons and
 * extracts those loops of the convolution which are part of the Minkowsi
 * sum. This method works very efficiently, regardless of whether `P` and
 * `Q` are convex or non-convex.
 * Note that as the input polygons may not be convex, their Minkowski
 * sum may not be a simple polygon. The result is therefore represented
 * as a polygon with holes.
 * \pre Both `P` and `Q` are simple, counterclockwise-oriented polygons.
*/

template <class Kernel_, class Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_reduced_convolution_2(const Polygon_2<Kernel_, Container_>& pgn1,
                                    const Polygon_2<Kernel_, Container_>& pgn2)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Minkowski_sum_by_reduced_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel,Container>                               sum_bound;
  std::list<Polygon_2<Kernel,Container> >                   sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum (pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else
    mink_sum (pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container> (sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

/*!
 * Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given
 * polygons-with-holes.
 * The function computes the reduced convolution of the two polygons and
 * extracts those loops of the convolution which are part of the Minkowsi
 * sum. This method works very efficiently, regardless of whether `P` and
 * `Q` are convex or non-convex.
 * The result is also represented as a polygon with holes.
*/
template <class Kernel_, class Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_reduced_convolution_2(const Polygon_with_holes_2<Kernel_, Container_>& pgn1,
                                    const Polygon_with_holes_2<Kernel_, Container_>& pgn2)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Minkowski_sum_by_reduced_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel,Container>                               sum_bound;
  std::list<Polygon_2<Kernel,Container> >                   sum_holes;

  mink_sum (pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));

  return (Polygon_with_holes_2<Kernel,Container> (sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

/*!
 * Compute the Minkowski sum of two simple polygons using the (full)
 * convolution method.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in) The second polygon.
 * \return The resulting polygon with holes, representing the sum.
 */
template <class Kernel_, class Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_full_convolution_2(const Polygon_2<Kernel_, Container_>& pgn1,
                                 const Polygon_2<Kernel_, Container_>& pgn2)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Minkowski_sum_by_convolution_2<Kernel, Container> mink_sum;
  Polygon_2<Kernel, Container>                      sum_bound;
  std::list<Polygon_2<Kernel, Container> >          sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum(pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else
    mink_sum(pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));
  return (Polygon_with_holes_2<Kernel, Container>(sum_bound,
                                                  sum_holes.begin(),
                                                  sum_holes.end()));
}

/*!
 * Compute the Minkowski sum of two simple polygons using the (full)
 * convolution method.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in) The second polygon.
 * \param kernel (in) The kernel.
 * \return The resulting polygon with holes, representing the sum.
 */
template <class Kernel_, class Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_full_convolution_2(const Polygon_2<Kernel_, Container_>& pgn1,
                                 const Polygon_2<Kernel_, Container_>& pgn2,
                                 const Kernel_& kernel)
{
  typedef Kernel_                                    Kernel;
  typedef Container_                                 Container;

  Minkowski_sum_by_convolution_2<Kernel, Container>  mink_sum(kernel);
  Polygon_2<Kernel, Container> sum_bound;
  std::list<Polygon_2<Kernel, Container> > sum_holes;

  if (pgn1.size() > pgn2.size())
    mink_sum(pgn1, pgn2, sum_bound, std::back_inserter(sum_holes));
  else
    mink_sum(pgn2, pgn1, sum_bound, std::back_inserter(sum_holes));
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
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in) The second polygon.
 * \return The resulting polygon with holes, representing the sum.
 *
 * \sa `CGAL::minkowski_sum_reduced_convolution_2()`
 * \sa `CGAL::minkowski_sum_full_convolution_2()`
 */
template <typename Kernel_, typename Container_>
Polygon_with_holes_2<Kernel_, Container_>
minkowski_sum_2(const Polygon_2<Kernel_, Container_>& pgn1,
                const Polygon_2<Kernel_, Container_>& pgn2)
{
  return minkowski_sum_reduced_convolution_2(pgn1, pgn2);
}

/*!
 * Compute the Minkowski sum of two simple polygons by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * Note that as the input polygons may not be convex, their Minkowski sum may
 * not be a simple polygon. The result is therefore represented as a polygon
 * with holes.
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in) The second polygon.
 * \param decomposition_strategy (in) A functor for decomposing polygons.
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
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in) The second polygon.
 * \param decomposition_strategy (in) A functor for decomposing polygons.
 * \param traits (in) traits The traits.
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
  typedef Kernel_                               Kernel;
  typedef Container_                            Container;
  typedef DecompositionStrategy_                Decomposition_strategy;

  Minkowski_sum_by_decomposition_2<Decomposition_strategy, Container>
    mink_sum(decomposition_strategy, traits);
  Polygon_with_holes_2<Kernel, Container> sum = mink_sum(pgn1, pgn2);
  return (sum);
}

/*!
 * Compute the Minkowski sum of two polygon-with-holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in) The second polygon.
 * \param decomposition_strategy (in) A functor for decomposing polygons.
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
 * Compute the Minkowski sum of two polygon-with-holes by decomposing each
 * polygon to convex sub-polygons and computing the union of the pairwise
 * Minkowski sums of the sub-polygons.
 * The result is also represented as a polygon with holes.
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in) The second polygon.
 * \param decomposition_strategy (in) A functor for decomposing polygons.
 * \param traits (in) The traits.
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
  Polygon_with_holes_2<Kernel, Container> sum = mink_sum(pgn1, pgn2);
  return sum;
}

/*!
 * Compute the Minkowski sum of one simple polygon and one polygon-with-holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in)The second polygon.
 * \param decomposition_strategy (in) A functor for decomposing polygons.
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
 * Compute the Minkowski sum of one simple polygon and one polygon-with-holes
 * by decomposing each polygon to convex sub-polygons and computing the union
 * of the pairwise Minkowski sums of the sub-polygons.  The result is also
 * represented as a polygon with holes.
 * \param pgn1 (in) The first polygon.
 * \param pgn2 (in)The second polygon.
 * \param decomposition_strategy (in) A functor for decomposing polygons.
 * \param traits (in) The traits.
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
  Polygon_with_holes_2<Kernel, Container> sum = mink_sum(pgn1, pgn2);
  return sum;
}

} //namespace CGAL

#endif
