// Copyright (c) 2016 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_SMS_2_PULLOUT_DIRECTIONS_SINGLE_MOLD_TRANSLATIONAL_CASTING_H
#define CGAL_SMS_2_PULLOUT_DIRECTIONS_SINGLE_MOLD_TRANSLATIONAL_CASTING_H

#include <CGAL/license/Set_movable_separability_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Set_movable_separability_2/internal/Utils.h>
#include <CGAL/Set_movable_separability_2/internal/Circle_arrangment.h>

namespace CGAL {
namespace Set_movable_separability_2 {
namespace Single_mold_translational_casting {

/*! Same as below with the additional traits argument.
 * \param[in] traits the traits to use.
 *
 *   algorithm:
 *           this function implements a very simple algorithm... it just keep at any stage the current
 *           intersection in  [firstClockwise,secondClockwise].
 *           When a new semicircle appear the possible cases are as such:
 *           (let f:=firstClockwise, s:=secondClockwise, a:=newSemicircleFirstClockwise , b:=newSemicircleSecondClockwise)
 *        REMEBER THAT THIS ARE SEGMENTS ON A CIRCLE! NOT ON A LINE!
 * 1. [f,s] contained in [a,b]
 *                f           s                    *      f                s   *          f        s  *     f        s
 *           a                b          *      a                b   *     a                b  *     a                b
 *           _________________               *      _________________  *     _________________*     _________________
 *                f           s                  *      f                s   *          f        s  *     f        s
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * 2.  a contained in (f,s] and b is not  /  or in other words / s in [a,b) and f is not in [a,b] (it is enough to ask if s is in [a,b] since fs+ab is less than 2*pi)
 *                    f                s          *         f        s
 *                   a                b *                   a        b
 *           _________________           *           _________________
 *                   f        s          *                   fs
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * 3. b contained in [f,s) and a is not /  or in other words / f in (a,b] and s is not in [a,b] (it is enough to ask if f is in [a,b] since fs is shorter the ab)
 *                        f            s          *                     f        s
 *       a                b          *   a                b
 *           _________________ *           _________________
 *              f        s          *                    fs
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * 4.         no intersection between [f,s] and [a,b] / case a: or in other words /  f,s are not in [a,b]
 *                f           s              *          f                s
 *           b                a   *        b                a
 *           _________________   *         _________________
 *           NO INTERSECTION!    *          NO INTERSECTION! (the only case in which this is possible is if (f,s) was not changes, and then (f,s) is an open arc)
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * 5. Illegal cases
 *          f                s *          f                s
 *                a           b          *                 b          a
 *           __________________*        __________________
 *           THIS CASE CANT HAPPEN!! [a,b] is an semicircle, and (f,s) is a semicircle or less
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& i,
 CGAL::Orientation orientation,
 CastingTraits_2& traits)
{
  CGAL_precondition(pgn.is_simple());
  CGAL_precondition(!internal::is_any_edge_colinear(pgn, traits));
  CGAL_precondition(pgn.edges_end()!=i);

  typedef CastingTraits_2 Casting_traits_2;
  //the returned range is [clock_first, clock_second]
  typename Casting_traits_2::Direction_2 clock_first, clock_second;


  auto segment_outer_circle =
    internal::get_segment_outer_circle<Casting_traits_2>(*i, orientation);
  clock_first = segment_outer_circle.first;
  clock_second = segment_outer_circle.second;
  //well theoretically, this is a bug since the current intersection is
  //currently (clock_first,clock_second) and not [clock_first,clock_second].. but
  //this edges will surly change since we are in a polygon

  bool is_range_smaller_than_semicircle(false);
  auto cc_in_between = traits.counterclockwise_in_between_2_object();

  for (auto e_it = pgn.edges_begin(); e_it != pgn.edges_end(); ++e_it) {
    if (e_it==i) continue;
    // std::cout << "f " << clock_first << " s " << clock_second << std::endl;
    auto segment_outer_circle =
      internal::get_segment_outer_circle<Casting_traits_2>(*e_it, orientation);
    // std::cout << "a "<< segment_outer_circle.second << " b "
    //           << segment_outer_circle.first<<std::endl;

    // notice that we are interested in the segment_inner_circle
    // (segment_outer_circle.second,segment_outer_circle.first)
    if (!is_range_smaller_than_semicircle) {
      if ((segment_outer_circle.first == clock_second) &&
          (segment_outer_circle.second == clock_first))
      {
        // std::cout<<"case 1b"<<std::endl<<std::endl;
        // the arc is the range case 1b
        continue;
      }
      if ((segment_outer_circle.first == clock_first) &&
          (segment_outer_circle.second == clock_second))
      {
        // std::cout<<"case 4b"<<std::endl<<std::endl;

        // the arc the opposite of the range case 4b
        return std::make_pair(false, std::make_pair(clock_first, clock_second));
      }
      is_range_smaller_than_semicircle = true;
    }
    bool f_between_ab = !cc_in_between(clock_first, segment_outer_circle.second,
                                     segment_outer_circle.first);
    //is true if segment_outer_circle \in [first,clock_first,clock_second]
    bool s_between_ab = !cc_in_between(clock_second,
                                     segment_outer_circle.second,
                                     segment_outer_circle.first);
    //is true if segment_outer_circle \in [first,clock_first,clock_second]
    if (f_between_ab && s_between_ab) {
      // std::cout<<"case 1"<<std::endl<<std::endl;
      // case 1 //surly not case 4b since [f,s] is less then a semicircle
      continue;
    }
    if (!f_between_ab && s_between_ab) {
      // std::cout<<"case 2"<<std::endl<<std::endl;
      // case 2 - return a,s
      clock_first = segment_outer_circle.second;
    }
    else if(f_between_ab && !s_between_ab) {
      // std::cout<<"case 3"<<std::endl<<std::endl;
      // case 3 - return f,b
      clock_second = segment_outer_circle.first;
    }
    else {
      //  std::cout<<"case 4a"<<std::endl<<std::endl;
      //case 4a
      return std::make_pair(false, std::make_pair(clock_first, clock_second));
    }
  }

  return std::make_pair(true, std::make_pair(clock_first, clock_second));
}

/*! Given a simple polygon and an edge of the polygon, this function determines
 * whether a cavity (of a mold in the plane) that has the shape of the polygon
 * can be used so that the polygon could be casted in the mold with the input
 * edge being the top edge and then pulled out of the mold without colliding
 * into the mold (but possibly sliding along the mold surface). If the polygon
 * is <em>castable</em> this way, the function computes the closed range of pull
 * directions.
 *
 * The type that substitutes the template parameter `%CastingTraits_2` must be
 * a model of the concept `CastingTraits_2`.
 *
 * \param[in] pgn the input polygon.
 * \param[in] i the iterator of an edge in pgn.
 * \return a pair of elements, where the first is a Boolean that indicates
 *         whether the input edge is a valid top edge, and the second
 *         is a closed range of pull-out directions represented as a pair
 *         of the extreme directions in the range. If the input edge is not
 *         a valid top edge, the range is nondeterministic.
 *         a pair of Directions is build this way [firstClockwise,secondClockwise]
 *
 * \pre `png` must be non-degenerate (has at least 3 vertices), simple, and
 * does not have three consecutive collinear vertices.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& it,
 CGAL::Orientation orientation)
{
  CastingTraits_2 traits;
  return pullout_directions(pgn, it, orientation, traits);
}

/*! Same as above with the orientation argument.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& it,
 CastingTraits_2& traits)
{
  CGAL::Orientation orientation = pgn.orientation();
  return pullout_directions(pgn, it, orientation, traits);
}

/*! Same as above with the orientation and traits arguments.
 */
template <typename CastingTraits_2>
std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
                          typename CastingTraits_2::Direction_2> >
pullout_directions
(const CGAL::Polygon_2<CastingTraits_2>& pgn,
 const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& it)
{
  CGAL::Orientation orientation = pgn.orientation();
  CastingTraits_2 traits;
  return pullout_directions(pgn, it, orientation, traits);
}

} // namespace Single_mold_translational_casting
} // namespace Set_movable_separability_2
} // namespace CGAL

#endif
