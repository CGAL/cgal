/*
 * utils.h
 *
 *  Created on: Jan 17, 2017
 *      Author: root
 */

#ifndef SET_MOVABLE_SEPARABILITY_2_INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_2_UTILS_H_
#define SET_MOVABLE_SEPARABILITY_2_INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_2_UTILS_H_
#include <CGAL/enum.h>
#include <CGAL/Polygon_2.h>


namespace CGAL {
  namespace Set_movable_separability_2 {
    namespace internal {

      /*! \fn std::pair<typename Kernel::Direction_2,typename Kernel::Direction_2> get_segment_outer_circle(typename Kernel::Segment_2 seg, CGAL::Orientation orientation)
       * \param[in] seg the polygon segment
       * \param[in] orientation the orientation of the segment (and the polygon).
       *   if CLOCKWISE then the outer half circle is to the left.
       * \return the open outer half-circle of the edge.
       */
      template <typename Kernel>
      inline std::pair<typename Kernel::Direction_2, typename Kernel::Direction_2>
      get_segment_outer_circle(const typename Kernel::Segment_2 seg,
			       const CGAL::Orientation orientation)
			       {
	typename Kernel::Direction_2 forward( seg);
	typename Kernel::Direction_2 backward(-forward);
	return (orientation == CGAL::Orientation::CLOCKWISE) ?
	    std::make_pair(backward, forward) : std::make_pair(forward, backward);
			       }

      template <typename Kernel>
      bool is_any_edge_colinear(const CGAL::Polygon_2<Kernel>& pgn)
      {
	typedef typename CGAL::Point_2<Kernel>                Point_2;
	typedef typename CGAL::Polygon_2<Kernel>              Polygon_2;
	typedef typename Polygon_2::Vertex_const_iterator     Vertex_const_iterator;
	Vertex_const_iterator vci = pgn.vertices_begin();
	Point_2 firstVar = *(vci++);
	Point_2 secondVar = *(vci++);
	Point_2 thirdVar = *(vci++);
	for (; vci != pgn.vertices_end(); ++vci) {
	    firstVar = secondVar;
	    secondVar = thirdVar;
	    thirdVar = *vci;
	    if (CGAL::collinear(firstVar, secondVar, thirdVar)) return true;
	}
	vci = pgn.vertices_begin();
	firstVar = secondVar;
	secondVar = thirdVar;
	thirdVar = *(vci++);
	if(CGAL::collinear(firstVar, secondVar, thirdVar)) return true;

	firstVar = secondVar;
	secondVar = thirdVar;
	thirdVar = *(vci++);
	if (CGAL::collinear(firstVar, secondVar, thirdVar)) return true;

	return false;
      }
    } // end of namespace internal
  } // end of namespace Set_movable_separability_2
} // end of namespace CGAL



#endif /* SET_MOVABLE_SEPARABILITY_2_INCLUDE_CGAL_SET_MOVABLE_SEPARABILITY_2_UTILS_H_ */
