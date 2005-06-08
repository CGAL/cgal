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
#ifndef CGAL_ARR_LM_RANDOM_GENERATOR_H
#define CGAL_ARR_LM_RANDOM_GENERATOR_H

/*! \file
* Definition of the Arr_random_landmarks_generator<Arrangement> template.
*/

#include <CGAL/Arr_point_location/Arr_lm_generator.h>
#include <CGAL/Random.h>

CGAL_BEGIN_NAMESPACE

/*! \class
* This class is related to the Landmarks point location, and given as 
* a parameter (or template parameter) to it. 
* It inherites from Arr_lm_generator and  implements the 
* function called "void _create_point_list(Point_list &)" 
* to creates the set of random landmarks.
*/
template <class Arrangement_, 
		  class Nearest_neighbor_ 
			  = Arr_landmarks_nearest_neighbor <typename Arrangement_::Traits_2> >
class Arr_random_landmarks_generator 
	: public Arr_landmarks_generator <Arrangement_, Nearest_neighbor_>
{
public:
	typedef Arrangement_								Arrangement_2;
	typedef typename Arrangement_2::Traits_2			Traits_2;
	typedef Arr_random_landmarks_generator<Arrangement_2, Nearest_neighbor_> 	
														Self;
	typedef typename Traits_2::Point_2					Point_2;
	typedef std::vector<Point_2>						Points_set;
	typedef typename Arrangement_2::Vertex_const_iterator	
														Vertex_const_iterator;



protected:

	// Data members:
	int number_of_landmarks; 

private:

  /*! Copy constructor - not supported. */
  Arr_random_landmarks_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

	
public: 
	  /*! Constructor. */
	  Arr_random_landmarks_generator 
		  (const Arrangement_2& arr, int lm_num = -1) : 
	      Arr_landmarks_generator<Arrangement_2, Nearest_neighbor_> (arr), 
		  number_of_landmarks (lm_num)
	  {
		  PRINT_DEBUG("Arr_random_landmarks_generator constructor. "
			  <<"number_of_landmarks = "<< number_of_landmarks); 

		  build_landmarks_set();
	  }

protected:
  /*!
   * create a set of random points. 
   * the number of points is given as a parametr to this class' constructor.
   * the random coordinates are chosen using CGAL::Random function.
   * The coordinates are of type double, and are selected randomly in the 
   * Arrangement's bounding rectangle. This is actually the bounding rectangle 
   * of the Arrangement's vertices 
   */
	virtual void _create_points_set (Points_set & points)
	{
		PRINT_DEBUG("create_random_points_list");

		//find bounding box
		double x_min=0, x_max=0, y_min=0, y_max=0;
		double x,y;
		Vertex_const_iterator vit; 
		for (vit=p_arr->vertices_begin(); vit != p_arr->vertices_end(); vit++)
		{
			x = CGAL::to_double(vit->point().x());
			y = CGAL::to_double(vit->point().y());
			if (x < x_min) x_min = x;
			if (x > x_max) x_max = x;
			if (y < y_min) y_min = y;
			if (y > y_max) y_max = y;
		}

		CGAL::Random     random;

		// n is the number of random points.
		//if n is not given in the constructor then this number
		//is set to be the number of vertices in the arrangement.
		int n; 
		if (number_of_landmarks > 0)
			n = number_of_landmarks;
		else
			n= p_arr->number_of_vertices();

		//loop n time 
		for (int i=0; i<n; i++) 
		{
			//create the random point
			double px = random.get_double(x_min, x_max);
			double py = random.get_double(y_min, y_max);
			Point_2 p(px, py);

			//put in a list 
			points.push_back(p); 

			PRINT_DEBUG("random point "<<i<< " is= " << p);
		}
		PRINT_DEBUG("end create_random_points_list");

	}

};

CGAL_END_NAMESPACE


#endif
