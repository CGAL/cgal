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
#ifndef CGAL_ARR_LM_HALTON_GENERATOR_H
#define CGAL_ARR_LM_HALTON_GENERATOR_H

/*! \file
* Definition of the Arr_halton_landmarks_generator<Arrangement> template.
*/

#include <CGAL/Arr_point_location/Arr_lm_generator.h>

CGAL_BEGIN_NAMESPACE

/*! \class
* This class is related to the Landmarks point location, and given as 
* a parameter (or template parameter) to it. 
* It inherites from Arr_lm_generator and  implements the 
* function called "void _create_point_list(Point_list &)" 
* to creates the set of halton sequence points as landmarks.
*/
template <class Arrangement_, 
		  class Nearest_neighbor_ 
			  = Arr_landmarks_nearest_neighbor <typename Arrangement_::Traits_2> >
class Arr_halton_landmarks_generator 
	: public Arr_landmarks_generator <Arrangement_, Nearest_neighbor_>
{
public:
	typedef Arrangement_								Arrangement_2;
	typedef typename Arrangement_2::Traits_2			Traits_2;
	typedef Arr_halton_landmarks_generator<Arrangement_2, Nearest_neighbor_> 	
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
  Arr_halton_landmarks_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

	
public: 
	  /*! Constructor. */
	  Arr_halton_landmarks_generator 
		  (const Arrangement_2& arr, int lm_num = -1) : 
	      Arr_landmarks_generator<Arrangement_2, Nearest_neighbor_> (arr), 
		  number_of_landmarks (lm_num)
	  {
		  PRINT_DEBUG("Arr_halton_landmarks_generator constructor. "
			  <<"number_of_landmarks = "<< number_of_landmarks); 

		  build_landmarks_set();
	  }

protected:
   /*!
   * create a set of halton points 
   * the number of points is given as a parametr to this class' constructor.
   * We first calculate the Arrangement's bounding rectangle. 
   * This is actually the bounding rectangle of the Arrangement's vertices.
   * Then, it calculates the halton sequence for this number of points, 
   * and scales it to fit the arrangement size.
   */
	virtual void _create_points_set (Points_set & points)
	{
		PRINT_DEBUG("create_halton_points_list");

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

		// n is the number of halton points.
		//if n is not given in the constructor then this number
		//is set to be the number of vertices in the arrangement.
		int n; 
		if (number_of_landmarks > 0)
			n = number_of_landmarks;
		else
			n= p_arr->number_of_vertices();

		//calculate
		double x_trans = x_max - x_min;
		double y_trans = y_max - y_min;
		PRINT_DEBUG( "x_max= "<< x_max <<" x_min = "<< x_min );
		PRINT_DEBUG( "y_max= "<< y_max <<" y_min = "<< y_min );
		PRINT_DEBUG( "x_trans= "<< x_trans <<" y_trans = "<< y_trans );

		//create halton sequence
		double base_inv;
		int digit, i, seed2;
		int base[2], leap[2], seed[2];
		double r[2], px, py;
		int ndim = 2;
		int step = 1;
		seed[0] = seed[1] = 0; 
		leap[0] = leap[1] = 1;
		base[0] = 2;
		base[1] = 3;
		for ( step = 1; step <= n; step++ )
		{
			for ( i = 0; i < ndim; i++ )
			{
				seed2 = seed[i] + step * leap[i];
				r[i] = 0.0E+00;
				base_inv = 1.0E+00 / ( ( double ) base[i] );
				while ( seed2 != 0 )
				{
					digit = seed2 % base[i];
					r[i] = r[i] + ( ( double ) digit ) * base_inv;
					base_inv = base_inv / ( ( double ) base[i] );
					seed2 = seed2 / base[i];
				}
			}
			//i_to_halton ( ndim, step, seed, leap, base, r );
			//r[0] is the x_coord
			//r[1] is the y_coord
			px = r[0]*x_trans; 
			py = r[1]*y_trans;
			Point_2 p(px, py);

			//put in a list 
			points.push_back(p); 

			PRINT_DEBUG("halton point "<< step++ <<" is= " << p );
		}

		PRINT_DEBUG("end create_halton_points_list");
	}

};

CGAL_END_NAMESPACE


#endif
