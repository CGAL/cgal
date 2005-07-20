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
#ifndef CGAL_ARR_LANDMARKS_NEAREST_NEIGHBOR_H
#define CGAL_ARR_LANDMARKS_NEAREST_NEIGHBOR_H

/*! \file
* Definition of the Arr_landmarks_nearest_neighbor<Arrangement> template.
*/
#include <CGAL/basic.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>

CGAL_BEGIN_NAMESPACE

/*! \class
* A class that answers nearest neighbor queries.
* It recieves a set of points, and builds a kd-tree for them.
* Given a query point, it finds the closest point to the query.
*/
template <class Traits_>
class Arr_landmarks_nearest_neighbor 
{
public:
	typedef Traits_										                    Traits_2;
  typedef typename Traits_2::Approximate_number_type	  ANT;
	typedef typename Traits_2::Point_2					          Point_2;
	typedef Arr_landmarks_nearest_neighbor<Traits_2>	    Self;

	class  NN_Point_2 {
	public:
		NN_Point_2() 
		{ 
			vec[0]= vec[1]  = 0; 
		}

    NN_Point_2 (Point_2 &pnt)
		{ 
			m_point = pnt;
      Traits_2 traits; 
			vec[0]= traits.approximate_2_object()(pnt, 0);
      vec[1]= traits.approximate_2_object()(pnt, 1);
		}

		NN_Point_2 (Point_2 &pnt, Object &object) 
		{ 
			m_object = object;
			m_point = pnt;
      Traits_2 traits; 
			vec[0]= traits.approximate_2_object()(pnt, 0);
      vec[1]= traits.approximate_2_object()(pnt, 1);
		}

		Object &object() { return m_object; }

		Point_2 &point() { return m_point; }

		ANT& x()
		{
		  return vec[ 0 ];
		}
		
		ANT& y()
		{
		  return vec[ 1 ];
		}

		bool operator==(const NN_Point_2& p) const 
		{	
			return (x() == p.x()) && (y() == p.y())  ;
		}

		bool  operator!=(const NN_Point_2& p) const 
		{ 
			return ! (*this == p); 
		}

		ANT		    vec[2];
		Object		m_object;
		Point_2		m_point;
	}; //end of class


	class Construct_coord_iterator {
	public:
		const ANT* operator()(const NN_Point_2& p) const 
		{ return static_cast<const ANT*>(p.vec); }

		const ANT* operator()(const NN_Point_2& p, int)  const
		{ return static_cast<const ANT*>(p.vec+2); }
	};


	typedef CGAL::Search_traits<ANT, NN_Point_2, const ANT*,
								Construct_coord_iterator>	Traits;
	typedef CGAL::Orthogonal_k_neighbor_search<Traits>      Neighbor_search;
	typedef typename Neighbor_search::iterator              Neighbor_iterator;
	typedef typename Neighbor_search::Tree                  Tree;
	typedef std::list<NN_Point_2>                           Point_list;	
	typedef typename Point_list::const_iterator             Input_iterator;

protected:
	Tree            * tree;
	bool            b_valid_tree;

private:

  /*! Copy constructor - not supported. */
  Arr_landmarks_nearest_neighbor (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

public:
	  /*! Default constructor. */
	  Arr_landmarks_nearest_neighbor () : 
	      tree(0), 
		    b_valid_tree(false) 
	  {
	  }

	  /*! Distructor - must be because we allocated a tree */
	  ~Arr_landmarks_nearest_neighbor() 
	  {
		  if (b_valid_tree && tree){
			  delete tree;
		  }
		  b_valid_tree = false;
		  tree = 0;
	  }

	  /*! init should allocate the tree and initialize it with all points.
	  it will be better if we could create a blank tree, 
	  and afterwards add points to it
	  */
	  void init(Input_iterator begin, Input_iterator beyond) 
	  {
		  if (b_valid_tree || tree) {
			  std::cerr << "ERROR: init - tree exists" << std::endl;
			  return;
		  }

		  if (begin != beyond) {
			  tree = new Tree(begin, beyond);
			  b_valid_tree = true;
		  }
	  }

	  /*! clean - deletes the tree in order to create a new one later
	  */
	  void clean() 
	  {
		  if (b_valid_tree && tree){
			  delete tree;
		  }
		  b_valid_tree = false;
		  tree = 0;
	  }

	  /*! find the nearest point to the query (and its location object)
	  */
	  Point_2 find_nearest_neighbor(Point_2 query, Object &obj) const
	  {
		  //create NN_Point_2 from Point_2 
		  NN_Point_2 nn_query(query);

		  //use the tree to find nearest landmark
		  CGAL_assertion (b_valid_tree);
		  Neighbor_search search(*tree, nn_query, 1);

		  //get the result
		  NN_Point_2 nearest_p = search.begin()->first;

		  //get the object 
		  obj = nearest_p.object();
		  
		  return (nearest_p.point());
	  }
	 
};

CGAL_END_NAMESPACE


#endif
