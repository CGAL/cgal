// Copyright (c) 2016 Tel-Aviv University (Israel).
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
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_PULLOUT_DIRECTIONS_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H
#define CGAL_PULLOUT_DIRECTIONS_SINGLE_MOLD_TRANSLATIONAL_CASTING_2_H

#include <CGAL/Polygon_2.h>
#include "Set_movable_separability_2/Utils.h"
#include "Set_movable_separability_2/Circle_arrangment.h"

namespace CGAL {

  namespace Set_movable_separability_2 {


    /*! Same as below with the additional traits argument.
     * \param[in] traits the traits to use.
     *
     *   algorithm:
     *   	this function implements a very simple algorithm... it just keep at any stage the current
     *   	intersection in  [firstClockwise,secondClockwise].
     *   	When a new semicircle appear the possible cases are as such:
     *   	(let f:=firstClockwise, s:=secondClockwise, a:=newSemicircleFirstClockwise , b:=newSemicircleSecondClockwise)
     *	REMEBER THAT THIS ARE SEGMENTS ON A CIRCLE! NOT ON A LINE!
     * 1. [f,s] contained in [a,b]
     *   	     f	   s	  	  *      f		s   *          f	s  *     f	s
     *   	a		b	  *      a		b   *     a		b  *     a		b
     *   	_________________     	  *      _________________  *     _________________*     _________________
     *   	     f	   s		  *      f		s   *          f	s  *     f	s
     * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     * 2.  a contained in (f,s] and b is not  /  or in other words / s in [a,b) and f is not in [a,b] (it is enough to ask if s is in [a,b] since fs+ab is less than 2*pi)
     * 		   f		s	  *	 f	s
     *   		a		b *   		a	b
     *   	_________________ 	  *   	_________________
     *   		f	s	  *   		fs
     * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     * 3. b contained in [f,s) and a is not /  or in other words / f in (a,b] and s is not in [a,b] (it is enough to ask if f is in [a,b] since fs is shorter the ab)
     * 	    	   f	    s	  * 	    	f	s
     *       a		b	  *   a		b
     *   	_________________ *   	_________________
     *   	   f	s	  *    		fs
     * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     * 4. 	no intersection between [f,s] and [a,b] / case a: or in other words /  f,s are not in [a,b]
     *   	     f	   s  	    *  	f		s
     *   	b		a   *	b		a
     *   	_________________   * 	_________________
     *   	NO INTERSECTION!    *  	NO INTERSECTION! (the only case in which this is possible is if (f,s) was not changes, and then (f,s) is an open arc)
     * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     * 5. Illegal cases
     *  	f		s *  	f		s
     *   	     a	   b	  *    	     b	  a
     *   	__________________*	__________________
     *   	THIS CASE CANT HAPPEN!! [a,b] is an semicircle, and (f,s) is a semicircle or less
     */
    template <typename CastingTraits_2>
    std::pair<bool, std::pair<typename CastingTraits_2::Direction_2,
    typename CastingTraits_2::Direction_2> >
    pullout_directions_single_mold_translational_casting_2
    (const CGAL::Polygon_2<CastingTraits_2>& pgn, const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& i, CastingTraits_2& traits)
    {
      CGAL_precondition(pgn.is_simple());
      CGAL_precondition(!internal::is_any_edge_colinear(pgn));
      CGAL_precondition(pgn.edges_end()!=i);
      CGAL::Orientation poly_orientation = pgn.orientation();

      typedef CastingTraits_2               Casting_traits_2;
      typename Casting_traits_2::Direction_2 clockFirst, clockSecond; //the returned range is [clockFirst,clockSecond]


      auto segment_outer_circle =
	  internal::get_segment_outer_circle<Casting_traits_2>(*i, poly_orientation);
      clockFirst=segment_outer_circle.first;
      clockSecond=segment_outer_circle.second;
      //well theoretically, this is a bug since the current intersection is currently (clockFirst,clockSecond)
      //and not [clockFirst,clockSecond].. but this edges will surly change since we are in a polygon

      bool isRangeSmallerThanSemicircle=false;
      auto cc_in_between = traits.counterclockwise_in_between_2_object();

      for (auto e_it = pgn.edges_begin(); e_it != pgn.edges_end(); ++e_it) {
	  if(e_it==i) continue;
	  //std::cout<<"f "<<clockFirst<<" s "<<clockSecond<<std::endl;
	  auto segment_outer_circle =
	      internal::get_segment_outer_circle<Casting_traits_2>(*e_it, poly_orientation);
	  // std::cout<<"a "<<segment_outer_circle.second<<" b "<<segment_outer_circle.first<<std::endl;

	  //notice that we are interested in the segment_inner_circle (segment_outer_circle.second,segment_outer_circle.first)
	  if(!isRangeSmallerThanSemicircle)
	    {
	      if(segment_outer_circle.first==clockSecond && segment_outer_circle.second == clockFirst)
		{
		  // std::cout<<"case 1b"<<std::endl<<std::endl;

		  // the arc is the range case 1b
		  continue;
		}
	      if(segment_outer_circle.first==clockFirst&& segment_outer_circle.second ==clockSecond  )
		{
		  //	  std::cout<<"case 4b"<<std::endl<<std::endl;

		  // the arc the opposite of the range case 4b
		  return std::make_pair(false, std::make_pair(clockFirst, clockSecond));
		}
	      isRangeSmallerThanSemicircle=true;
	    }
	  bool fBetweenAB = !cc_in_between(clockFirst,segment_outer_circle.second,segment_outer_circle.first);
	  //is true if segment_outer_circle \in [first,clockFirst,clockSecond]
	  bool sBetweenAB = !cc_in_between(clockSecond,segment_outer_circle.second,segment_outer_circle.first);
	  //is true if segment_outer_circle \in [first,clockFirst,clockSecond]
	  if (fBetweenAB && sBetweenAB)
	    {
	      // std::cout<<"case 1"<<std::endl<<std::endl;

	      //case 1 //surly not case 4b since [f,s] is less then a semicircle
	      continue;
	    }
	  if (!fBetweenAB && sBetweenAB)
	    {
	      // std::cout<<"case 2"<<std::endl<<std::endl;

	      //case 2 - return a,s
	      clockFirst = segment_outer_circle.second;
	    }
	  else if(fBetweenAB && !sBetweenAB)
	    {
	      //  std::cout<<"case 3"<<std::endl<<std::endl;

	      //case 3 - return f,b
	      clockSecond = segment_outer_circle.first;
	    }
	  else
	    {
	      //  std::cout<<"case 4a"<<std::endl<<std::endl;

	      //case 4a
	      return std::make_pair(false, std::make_pair(clockFirst, clockSecond));
	    }
      }

      return std::make_pair(true, std::make_pair(clockFirst, clockSecond));
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
    pullout_directions_single_mold_translational_casting_2
    (const CGAL::Polygon_2<CastingTraits_2>& pgn,

     const typename CGAL::Polygon_2<CastingTraits_2>::Edge_const_iterator& i)
     {

      CastingTraits_2 traits;
      return pullout_directions_single_mold_translational_casting_2(pgn, i, traits);
     }



  } // end of namespace Set_movable_separability_2
} // end of namespace CGAL

#endif
