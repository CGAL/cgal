// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
#ifndef CGAL_ARR_LM_MIDDLE_EDGES_GENERATOR_H
#define CGAL_ARR_LM_MIDDLE_EDGES_GENERATOR_H

/*! \file
* Definition of the Arr_middle_edges_landmarks_generator<Arrangement> template.
*/

#include <CGAL/Arr_point_location/Arr_lm_generator.h>

CGAL_BEGIN_NAMESPACE

/*! \class
* This class is related to the Landmarks point location, and given as 
* a parameter (or template parameter) to it. 
* It inherites from Arr_lm_generator and  implements the 
* function called "void _create_point_list(Point_list &)" 
* to creates the set of landmarks, which are the middle points of the 
* arrangement edges, which must be segments !
* IMPORTANT: THIS ALGORITHM WORKS ONLY FOR SEGMENTS !!!
*/
template <class Arrangement_, 
		  class Nearest_neighbor_ 
			  = Arr_landmarks_nearest_neighbor <typename Arrangement_::Traits_2> >
class Arr_middle_edges_landmarks_generator 
	: public Arr_landmarks_generator <Arrangement_, Nearest_neighbor_>
{
public:
	typedef Arrangement_								Arrangement_2;
	typedef typename Arrangement_2::Traits_2			Traits_2;
	typedef Arr_middle_edges_landmarks_generator<Arrangement_2, Nearest_neighbor_> 	
														Self;
	typedef typename Traits_2::Point_2					Point_2;
	typedef std::vector<Point_2>						Points_set;
	typedef typename Arrangement_2::Edge_const_iterator	
													Edge_const_iterator;
	typedef typename Arrangement_2::Halfedge_const_handle	
													Halfedge_const_handle;

private:

  /*! Copy constructor - not supported. */
  Arr_middle_edges_landmarks_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

	
public: 
	  /*! Constructor. */
	  Arr_middle_edges_landmarks_generator 
		  (const Arrangement_2& arr, int lm_num = -1) : 
	      Arr_landmarks_generator<Arrangement_2, Nearest_neighbor_> (arr)
	  {
		  PRINT_DEBUG("Arr_middle_edges_landmarks_generator constructor.");

		  build_landmarks_set();
	  }

protected:
   /*!
   * create a set of middle_edges points 
   * the number of points is equal to the number of edges in the arrangement.
   */
	virtual void _create_points_set (Points_set & points)
	{
		PRINT_DEBUG("create_middle_edges_points_list");

		Edge_const_iterator    eit;
		Halfedge_const_handle  hh;
		for (eit=p_arr->edges_begin(); eit != p_arr->edges_end(); eit++)
		{
			//get 2 endpoints of edge
			hh = *eit;
			Point_2 p1 = hh.source().point();
			Point_2 p2 = hh.target().point();
			Point_2 p ((p1.x()+p2.x())/2, (p1.y()+p2.y())/2);
	
			//put in a list 
			points.push_back(p); 

			PRINT_DEBUG("mid point is= " << p);
		} 

		PRINT_DEBUG("end create_middle_edges_points_list");
	}

};

CGAL_END_NAMESPACE


#endif
