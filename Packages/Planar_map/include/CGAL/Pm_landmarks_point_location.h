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

#ifndef CGAL_PM_LANDMARKS_POINT_LOCATION_H
#define CGAL_PM_LANDMARKS_POINT_LOCATION_H

//#define CGAL_LM_DEBUG
#define LM_CLOCK_DEBUG

//----------------------------------------------------------
//Pm includes
//----------------------------------------------------------
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Pm_point_location_base.h>
#include <iostream>
#include <stdio.h>
#include <time.h>

////////////////////////////////////////////////////////
//    LANDMARKS STRATEGY
////////////////////////////////////////////////////////

CGAL_BEGIN_NAMESPACE

template <class Planar_map, class Nearest_neighbor> 
class Pm_landmarks_point_location : 
	public Pm_point_location_base<Planar_map> {
public:
	//----------------------------------------------------------
	// Pm Types
	//----------------------------------------------------------

	typedef typename Planar_map::Traits                       Traits;
	typedef typename Traits::Kernel                           Kernel;
	typedef typename Kernel::Segment_2                        Segment;
	typedef typename Traits::Point_2                          Point_2;
	typedef typename Traits::Curve_2                          Curve_2;
	typedef typename Traits::X_monotone_curve_2               X_monotone_curve_2;
	typedef Pm_point_location_base<Planar_map>                Base;
	typedef Pm_landmarks_point_location<Planar_map, Nearest_neighbor>    Self;
	typedef typename Planar_map::Face_iterator                Face_iterator;
	typedef typename Planar_map::Halfedge_iterator            Halfedge_iterator;
	typedef typename Planar_map::Vertex_iterator              Vertex_iterator;
	typedef typename Planar_map::Edge_iterator                Edge_iterator;
	typedef typename Planar_map::Vertex_handle                Vertex_handle;
	typedef typename Planar_map::Vertex_const_handle          Vertex_const_handle;
	typedef typename Planar_map::Halfedge_handle              Halfedge_handle;
	typedef typename Planar_map::Halfedge_const_handle        Halfedge_const_handle;
	typedef typename Planar_map::Face_handle                  Face_handle;
	typedef typename Planar_map::Halfedge_around_vertex_circulator 
		Halfedge_around_vertex_circulator;
	typedef typename Planar_map::Ccb_halfedge_circulator Ccb_halfedge_circulator;
	typedef typename Planar_map::Ccb_halfedge_const_circulator 
		Ccb_halfedge_const_circulator;
	typedef typename Base::Halfedge_handle_iterator      Halfedge_handle_iterator;
	typedef typename Planar_map::Holes_iterator Holes_iterator;
	typedef typename Planar_map::Holes_const_iterator Holes_const_iterator;

	typedef typename Base::Token                              Token;
	typedef typename Planar_map::Locate_type                  Locate_type;
	typedef typename Planar_map::Traits_wrap                  Traits_wrap;
	typedef typename Nearest_neighbor::NN_Point_2      NN_Point_2;

	typedef std::list<NN_Point_2>                                                 NN_Point_list;
	//----------------------------------------------------------

protected:
	// typedef const Self* cPLp;
	typedef const Self* const_Self_ptr;

public:
	// Constructor
	Pm_landmarks_point_location() : 
	  Pm_point_location_base<Planar_map>(),
		  pm(0),
		  traits(0),
		  updated_nn(false), 
		  verbose(false)
	  {
			#ifdef LM_CLOCK_DEBUG
				clock_ff = 0.0; 
				clock_fi= 0.0; 
				clock_ni= 0.0; 
				clock_for_nn_search = 0.0; 
				clock_for_walk = 0.0; 
				clock_fciif = 0.0;
				clock_new_alg = 0.0;
				entries_to_fi = 0;
			#endif
	  }

	  //Destructor
	~Pm_landmarks_point_location() 
	  {
		  #ifdef LM_CLOCK_DEBUG
				std::cout << "total time to walk is " << clock_for_walk <<" clocks" << std::endl;
				std::cout << "total time to nn search is " << clock_for_nn_search <<" clocks" << std::endl;
				std::cout << "total time to ff (find face) is " << clock_ff <<" clocks" << std::endl;
				std::cout << "total time to fi (find intersection) is " << clock_fi <<" clocks" << std::endl;
				std::cout << "total time to ni (nearest intersection) is " << clock_ni <<" clocks" << std::endl;
				std::cout << "total time to fciif(find closest intersection in face) is " << clock_fciif <<" clocks" << std::endl;
				std::cout << "total time to new algorithm is " << clock_new_alg <<" clocks" << std::endl;
				std::cout << "total entries to fi(find intersection) is " << entries_to_fi <<" times" << std::endl;
				getchar();
		  #endif
	  }

	  void init(Planar_map & pmp, const Traits & tr) 
	  {
		  CGAL_precondition_msg(pm == NULL,
			  "Point location instance should be uninitialized "
			  "(Do not use the same instance for more than one map).");

		  pm = &pmp;
		  traits = (Traits_wrap*)(&tr);
	  }

	  //   void init(Planar_map & pmp, Traits & tr) 
	  //   {
	  // #ifdef CGAL_LM_DEBUG
	  //     std::cout << "init PL" << std::endl;
	  // #endif

	  //     CGAL_precondition_msg(pm == NULL,
	  //     "Point location instance should be uninitialized "
	  //     "(Do not use the same instance for more than one map).");

	  //     pm = &pmp;
	  //     traits = (Traits_wrap*)(&tr);

	  //     create_landmarks_tree();
	  //   }

	  inline void insert(Halfedge_handle hh, const X_monotone_curve_2 & cv) 
	  {insert_halfedge_to_ln_tree(hh, cv); }

	  Halfedge_const_handle locate(const Point_2 & p, Locate_type & lt) const;

	  Halfedge_handle locate(const Point_2 & p, Locate_type & lt);

	  Halfedge_const_handle vertical_ray_shoot(const Point_2& p, Locate_type& lt, bool up)
		  const;
	  Halfedge_handle vertical_ray_shoot(const Point_2& p, Locate_type& lt, bool up);

	  inline void split_edge(const X_monotone_curve_2 &, Halfedge_handle, Halfedge_handle,
		  //additions by iddo for arrangement
		  const X_monotone_curve_2 &, const X_monotone_curve_2 &) 
	  {updated_nn = false; create_landmarks_tree(); }

	  inline void merge_edge(const X_monotone_curve_2 &, const X_monotone_curve_2 &, Halfedge_handle, 
		  //additions by iddo for arrangement
		  const X_monotone_curve_2 &)   
	  {updated_nn = false;  create_landmarks_tree();}

	  inline void remove_edge(Halfedge_handle) 
	  {updated_nn = false;  create_landmarks_tree();}

	  inline void remove_edge(const Halfedge_handle_iterator &,
		  const Halfedge_handle_iterator &) 
	  {updated_nn = false; create_landmarks_tree(); };

	  inline void clear() 
	  {updated_nn = false; create_landmarks_tree();}

	  inline void update(const Halfedge_handle_iterator &,
		  const Halfedge_handle_iterator &,
		  const Token& token) 
	  {updated_nn = false; create_landmarks_tree(); }

private:

	//function that updates the kd-tree for the nearest neightbor
	void create_landmarks_tree() ;

	void insert_halfedge_to_ln_tree(Halfedge_handle hh, const X_monotone_curve_2 &cv) ;

	//function that walks from the vertex to the point
	//Halfedge_const_handle 
	void
		walk(Vertex_handle vh, 
		const Point_2 & p, 
		Halfedge_const_handle& e,
		Locate_type& lt) const;

	void find_face (const Point_2 & p, 
		Vertex_handle vh,
		bool & found_vertex_or_edge, 
		bool & new_vertex, 
		bool & found_face,
		Vertex_handle & out_vertex, 
		Halfedge_handle & out_edge,
		Locate_type& lt  ) const;

	void find_intersection (const Point_2 & p, 
		//Vertex_handle vh,
		Curve_2 &seg, 
		Halfedge_handle e,
		int  & num_of_intersections, 
		bool & change_side, 
		bool & found_edge,
		Point_2 & closest_interect_point,
		bool & new_vertex, 
		Vertex_handle & out_vertex ) const;

	bool is_point_in_face (const Point_2 & p, 
		const Ccb_halfedge_circulator & face, 
		bool & found_edge,             
		bool & found_vertex,         
		Halfedge_handle  & out_edge) const;

	bool find_closest_intersection_in_face (const Point_2 & p,            
		Vertex_handle  v,     
		const Ccb_halfedge_circulator & face,                     
		Halfedge_handle  & out_edge) const ;

	bool find_real_intersection (const Point_2 & p,   
		Vertex_handle  v,        
		Halfedge_handle e,
		Point_2 & out_point) const  ;			

#ifdef CGAL_LM_DEBUG

	void debug() {}

	void debug(const Halfedge_handle& e) const
	{
		{
			if (e!=pm->halfedges_end()) 
				std::cerr << "(" << e->source()->point() << "," 
				<< e->target()->point() << ")" << std::flush;
			else std::cerr << "(oo)";
		}
	}

#endif

public:
	inline const Traits * get_traits() const {return traits;}

protected:
	Planar_map      * pm;
	Traits_wrap     * traits;
	Nearest_neighbor  nn;
	bool              updated_nn;

	#ifdef LM_CLOCK_DEBUG
		mutable double clock_ff; //find face
		mutable double clock_fi; //find intersection
		mutable double clock_ni; //nearest intersection (to left/to right)
		mutable double clock_for_nn_search ; 
		mutable double clock_for_walk ; 
		mutable double clock_fciif;
		mutable double clock_new_alg;
		mutable int entries_to_fi;
	#endif
	const bool verbose; 
	};

CGAL_END_NAMESPACE

#include <CGAL/Pm_landmarks_point_location.C>

#endif //PM_LANDMARKS_POINT_LOCATION_H
