// Copyright (c) 2004  Tel-Aviv University (Israel).
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
// Author(s)     : Idit Haran <haranidi@post.tau.ac.il>

#ifndef CGAL_PM_NEAREST_NEIGHBOR_H
#define CGAL_PM_NEAREST_NEIGHBOR_H

//#define CGAL_PM_WALK_DEBUG
//#define CGAL_LM_DEBUG

//----------------------------------------------------------
//Pm includes
//----------------------------------------------------------
#include <CGAL/basic.h>
//#include <CGAL/Cartesian.h>
//#include <CGAL/MP_Float.h>
//#include <CGAL/Quotient.h>
//#include <CGAL/Planar_map_2/Pm_traits_wrap_2.h>
//#include <CGAL/point_generators_2.h>
//#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Pm_default_dcel.h>
#include <iostream>
#include <stdio.h>

//#define CGAL_LM_DEBUG

CGAL_BEGIN_NAMESPACE

/*! This is a model for the concept: 
nearest neigbor for landmarks point location.
It is the nearest neighbor algorithm, that should be changed easily
*/
template <class Planar_map>   class Pm_nearest_neighbor 
{
public:
	typedef typename Planar_map::Vertex_handle                Vertex_handle;

	class  NN_Point_2 {
	public:
		NN_Point_2() { vec[0]= vec[1]  = 0; }
		NN_Point_2 (double x, double y) { vec[0]=x; vec[1]=y;  }
		NN_Point_2 (double x, double y, Vertex_handle vh) { vec[0]=x; vec[1]=y;  vert = vh; }

		double x() const { return vec[ 0 ]; }
		double y() const { return vec[ 1 ]; }
		Vertex_handle vertex() const { return vert; }

		double& x() { return vec[ 0 ]; }
		double& y() { return vec[ 1 ]; }
		Vertex_handle vertex() {return vert; }

		bool operator==(const NN_Point_2& p) const 
		{	
			return (x() == p.x()) && (y() == p.y() && vertex() == p.vertex()) ;
		}

		bool  operator!=(const NN_Point_2& p) const 
		{ 
			return ! (*this == p); 
		}

		double vec[2];
		Vertex_handle  vert;
	}; //end of class


	class Construct_coord_iterator {
	public:
		const double* operator()(const NN_Point_2& p) const 
		{ return static_cast<const double*>(p.vec); }

		const double* operator()(const NN_Point_2& p, int)  const
		{ return static_cast<const double*>(p.vec+2); }
	};

	//----------------------------------------------------------
	// Types
	//----------------------------------------------------------

	typedef CGAL::Search_traits<double, NN_Point_2, const double*, Construct_coord_iterator>  Traits;
	typedef CGAL::Orthogonal_k_neighbor_search<Traits>        Neighbor_search;
	typedef typename Neighbor_search::iterator                             Neighbor_iterator;
	typedef typename Neighbor_search::Tree                                   Tree;
	typedef std::list<NN_Point_2>                                                     Point_list;	
	typedef typename Point_list::const_iterator                             Input_iterator;
	typedef Pm_nearest_neighbor<Planar_map>						    Self;
	//----------------------------------------------------------

public:
	/*! Constructor  */
	Pm_nearest_neighbor() : 
	  tree(0), 
		  b_valid_tree(false)
	  {
#ifdef CGAL_LM_DEBUG 
		  std::cout << "constructor of nearest neighbor " << std::endl;
#endif
	  }

	  /*! Distructor - must be because we allocated a tree */
	  ~Pm_nearest_neighbor() 
	  {
		  if (b_valid_tree && tree){
#ifdef CGAL_LM_DEBUG 
			  std::cout << "~: deleting the tree. tree = " << tree ;
			   std::cout << "size = " << tree->size() << std::endl;
			   //tree->statistics(std::cout);
#endif
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

		  if (begin == beyond) {

#ifdef CGAL_LM_DEBUG 
			  std::cout << "empty list" << std::endl;
#endif
		  }
		  else {
#ifdef CGAL_LM_DEBUG 
			  std::cout << "allocating a new tree" << std::endl;
#endif
			  tree = new Tree(begin, beyond);
			  b_valid_tree = true;
		  }
	  }

	  /*! clean - deletes the tree in order to create a new one later
	  */
	  void clean() 
	  {
		  if (b_valid_tree && tree){
#ifdef CGAL_LM_DEBUG 
			  std::cout << "~: deleting the tree. tree = " << tree ;
			   std::cout << "size = " << tree->size() << std::endl;
			   //tree->statistics(std::cout);
#endif
			  delete tree;
#ifdef CGAL_LM_DEBUG 
			  std::cout << "after deleting the tree" << std::endl;
#endif
		  }
		  else {
#ifdef CGAL_LM_DEBUG 
			  std::cout << "no tree to clean" << std::endl;
#endif
		  }
		  b_valid_tree = false;
		  tree = 0;
	  }


	  //Point_2 & find_nearest_neighbor(/*Point_2 & query*/)
	  NN_Point_2 find_nearest_neighbor(const NN_Point_2 & query) const
	  {
#ifdef CGAL_LM_DEBUG 
		  std::cout << "finding nearest neighbor of "<< query.x() <<' ' <<query.y();
#endif

		  // Initialize the search structure, and search all N points
		  Neighbor_search search(*tree, query, 1); 

		  // report the N nearest neighbors and their distance
		  // This should sort all N points by increasing distance from origin
		  for (Neighbor_iterator it = search.begin(); it != search.end(); ++it)
		  {
#ifdef CGAL_LM_DEBUG 
			  std::cout << "  is "<<(it->first).x()<<' '<<(it->first).y()<< std::endl;
#endif
		  }

		  return (search.begin()->first);
	  }

	  /*! ixx: should be implemented 
	  */
	  //void insert_point(Point_2 & p) { } 

	  /*! ixx: should be implemented
	  */
	  //void insert_points(Point_list & p_list) { } 


protected:
	Tree        * tree;
	bool          b_valid_tree;

public: 
	/*! Copy Constructor. */
	Pm_nearest_neighbor(const Self & nn) :
		tree(0),
		b_valid_tree(false)
	{
		std::cout << "nn copy contructor" << std::endl;
		*this = nn;
	}

	/*! Assignment Operator */
	Self & operator=(const Self & nn)
	{
		std::cout << "nn operator =" << std::endl;
		if (this == &nn)
			return (*this);

		if (tree != 0)
			delete tree;
		tree = 0;
		b_valid_tree = false;

		if (nn.b_valid_tree)
		{
			tree = new Tree (*(nn.tree));
			b_valid_tree = true;
		}

		return (*this);
	}
};

CGAL_END_NAMESPACE


#endif //CGAL_PM_NEAREST_NEIGHBOR_H
