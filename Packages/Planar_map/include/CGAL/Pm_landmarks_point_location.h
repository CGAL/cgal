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
#define LANDMARKS_CLOCK_DEBUG
//#define TRAITS_CLOCK_DEBUG

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
  //typedef typename Kernel::Segment_2                        Segment;
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

  typedef std::list<NN_Point_2>                                       NN_Point_list;
  typedef std::list<Halfedge_handle>                              Edge_list;
  typedef typename Edge_list::iterator               Std_edge_iterator;
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
    flipped_edges.clear();

#ifdef LANDMARKS_CLOCK_DEBUG
    clock_ff = 0.0; 
    clock_for_nn_search = 0.0; 
    clock_for_walk = 0.0; 
    clock_find_edge = 0.0;
    clock_new_alg = 0.0;
    clock_create_nn = 0.0;
    clock_is_point = 0.0;
    clock_check_app = 0.0;
    clock_nn_and_walk = 0.0;

    entries_find_edge = 0;
    entries_to_check_app = 0;
    entries_to_find_face = 0;
    entries_is_point_in_face = 0;
#endif

#ifdef TRAITS_CLOCK_DEBUG
    //count entries
    e_compare_distance = 0;
    e_point_equal = 0;
    e_curve_is_between_cw = 0;
    e_point_in_x_range = 0;
    e_curve_compare_y_at_x = 0;
    e_nearest_intersection_to_left = 0;
    e_nearest_intersection_to_right = 0;
    e_curve_split = 0;
    e_compare_xy = 0;
    e_curves_compare_y_at_x_left = 0;
    e_curves_compare_y_at_x_right = 0;
    e_curves_compare_cw = 0;
    //count clocks
    c_compare_distance = 0.0;
    c_point_equal = 0.0;
    c_curve_is_between_cw = 0.0;
    c_point_in_x_range = 0.0;
    c_curve_compare_y_at_x = 0.0;
    c_nearest_intersection_to_left = 0.0;
    c_nearest_intersection_to_right = 0.0;
    c_curve_split = 0.0;
    c_compare_xy = 0.0;
    c_curves_compare_y_at_x_left = 0.0;
    c_curves_compare_y_at_x_right = 0.0;
    c_curves_compare_cw = 0.0;

#endif
  }

  //Destructor
  ~Pm_landmarks_point_location() 
  {
#ifdef TRAITS_CLOCK_DEBUG
    std::cout << std::endl;
    //count entries
    //std::cout << "e_compare_distance =  " << e_compare_distance << std::endl;
    //std::cout << "e_point_equal =  " << e_point_equal << std::endl;
    std::cout << "e_curve_is_between_cw =  " << e_curve_is_between_cw << std::endl;
    //std::cout << "e_point_in_x_range =  " << e_point_in_x_range << std::endl;
    std::cout << "e_curve_compare_y_at_x =  " << e_curve_compare_y_at_x << std::endl;
    //std::cout << "e_nearest_intersection_to_left =  " << e_nearest_intersection_to_left << std::endl;
    //std::cout << "e_nearest_intersection_to_right =  " << e_nearest_intersection_to_right << std::endl;        
    std::cout << "e_curve_split =  " << e_curve_split << std::endl;
    std::cout << "e_compare_xy =  " << e_compare_xy << std::endl;
    std::cout << "e_curves_compare_y_at_x_left =  " << e_curves_compare_y_at_x_left << std::endl;
    std::cout << "e_curves_compare_y_at_x_right =  " << e_curves_compare_y_at_x_right << std::endl;
    std::cout << "e_curves_compare_cw =  " << e_curves_compare_cw << std::endl;
    std::cout << std::endl;
    //count clocks
    //std::cout << "c_compare_distance =  " << c_compare_distance << std::endl;
    //std::cout << "c_point_equal =  " << c_point_equal << std::endl;
    std::cout << "c_curve_is_between_cw =  " << c_curve_is_between_cw << std::endl;
    //std::cout << "c_point_in_x_range =  " <<  c_point_in_x_range<< std::endl;
    std::cout << "c_curve_compare_y_at_x =  " << c_curve_compare_y_at_x << std::endl;
    //std::cout << "c_nearest_intersection_to_left =  " << c_nearest_intersection_to_left << std::endl;
    //std::cout << "c_nearest_intersection_to_right =  " << c_nearest_intersection_to_right << std::endl;
    std::cout << "c_curve_split =  " << c_curve_split << std::endl;
    std::cout << "c_compare_xy =  " << c_compare_xy << std::endl;
    std::cout << "c_curves_compare_y_at_x_left =  " << c_curves_compare_y_at_x_left << std::endl;
    std::cout << "c_curves_compare_y_at_x_right =  " << c_curves_compare_y_at_x_right << std::endl;
    std::cout << "c_curves_compare_cw =  " << c_curves_compare_cw << std::endl;
    //getchar();
#endif

#ifdef LANDMARKS_CLOCK_DEBUG
    //std::cout << std::endl;
    //double seconds;
    //std::cout << "create landmarks tree  = " << clock_create_nn  << std::endl;        
    //std::cout << "nn + walk =  " << clock_nn_and_walk << std::endl;        
    //std::cout << "nn search = " << clock_for_nn_search  << std::endl;
    //seconds = (double) (clock_for_walk) / (double) CLOCKS_PER_SEC;
    //std::cout << "walk =  " << clock_for_walk <<" clocks, " << seconds << " seconds" << std::endl;
    //std::cout << "  find face =  " << clock_ff << std::endl;
    //std::cout << "  new algorithm = " << clock_new_alg  << std::endl;
    //std::cout << "    is_point_in_face =  " << clock_is_point  << std::endl;    
    //std::cout << "    find_edge_to_flip = " << clock_find_edge  << std::endl;
    //std::cout << "      check_approximate intersection = " << clock_check_app  << std::endl;
    //std::cout << std::endl;
    //std::cout << "entries_to_find_face =  " << entries_to_find_face  << std::endl;        
    //std::cout << "entries_is_point_in_face =  " << entries_is_point_in_face  << std::endl;        
    //std::cout << "entries find_edge_to_flip = " <<  entries_find_edge << std::endl; 
    //std::cout << "  entries_to_check_approximate_intersection =  " << entries_to_check_app << std::endl;
    //std::cout << "create landmarks tree  = " << clock_create_nn  << std::endl;        
    std::cout << "nn= " << clock_nn_and_walk - clock_for_walk  << std::endl;
    std::cout << "walk= " << clock_for_walk << std::endl;
    std::cout << "ff= " << clock_ff << std::endl;
    std::cout << "na= " << clock_new_alg  << std::endl;
    std::cout << "is= " << clock_is_point  << std::endl;    
    std::cout << "fe= " << clock_find_edge  << std::endl;
    std::cout << "ca= " << clock_check_app  << std::endl;
    std::cout << "eff= " << entries_to_find_face  << std::endl;        
    std::cout << "eis= " << entries_is_point_in_face  << std::endl;        
    std::cout << "efe= " <<  entries_find_edge << std::endl; 
    std::cout << "eca= " << entries_to_check_app << std::endl;
    //getchar();
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

  Halfedge_const_handle vertical_ray_shoot(const Point_2& p, Locate_type& lt,
                                           bool up) const;
  Halfedge_handle vertical_ray_shoot(const Point_2& p, Locate_type& lt,
                                     bool up);

  inline void split_edge(const X_monotone_curve_2 &, Halfedge_handle,
                         Halfedge_handle,
                         //additions by iddo for arrangement
                         const X_monotone_curve_2 &,
                         const X_monotone_curve_2 &)
  {updated_nn = false; }

  inline void merge_edge(const X_monotone_curve_2 &,
                         const X_monotone_curve_2 &, Halfedge_handle,
                         //additions by iddo for arrangement
                         const X_monotone_curve_2 &)   
  {updated_nn = false; }

  inline void remove_edge(Halfedge_handle) 
  {updated_nn = false; }

  inline void remove_edge(const Halfedge_handle_iterator &,
                          const Halfedge_handle_iterator &) 
  {updated_nn = false; };

  inline void clear() 
  {updated_nn = false; }

  inline void update(const Halfedge_handle_iterator &,
                     const Halfedge_handle_iterator &,
                     const Token& token) 
  {updated_nn = false; }

private:
  //function that updates the kd-tree for the nearest neightbor
  void create_landmarks_tree() const;

  template <class T_Kernel_Tag, class T_Point_2>
  void point_to_double_coords
  (const T_Point_2 & p, double & x, double & y, T_Kernel_Tag) const
  {
    x = CGAL::to_double(p.x());
    y = CGAL::to_double(p.y());
  }

#if defined(CEP_LEDA_RAT_KERNEL_TRAITS_H)
  template <class T_Point_2>
  void point_to_double_coords
  (const T_Point_2 & p, double & x, double & y, Leda_rational_kernel_tag) const
  {
    x = CGAL::to_double(p.xcoord());
    y = CGAL::to_double(p.ycoord());
  }
#endif

  /*! Extract the x and y coordinate of a point and convert to double. Dispatch
   * according to the parameterized kernel. This assumes that the type 'Kernel'
   * is defined by the traits, but this is NOT a requirement, which makes this
   * solution erroneous!
   */
  void point_to_double_coords(const Point_2 & p, double & x, double & y) const
  { point_to_double_coords(p, x, y, Kernel::Kernel_tag); }

  void insert_halfedge_to_ln_tree(Halfedge_handle hh,
                                  const X_monotone_curve_2 &cv) ;

  //function that walks from the vertex to the point
  //Halfedge_const_handle 
  void walk(Vertex_handle vh, 
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
                  Locate_type& lt) const;

  void new_find_face (const Point_2 & p, 
                      Vertex_handle vh,
                      bool & found_vertex_or_edge, 
                      bool & new_vertex, 
                      bool & found_face,
                      Vertex_handle & out_vertex, 
                      Halfedge_handle & out_edge,
                      Locate_type& lt) const;

  void find_intersection (const Point_2 & p, 
                          //Vertex_handle vh,
                          Curve_2 &seg, 
                          Halfedge_handle e,
                          int  & num_of_intersections, 
                          bool & change_side, 
                          bool & found_edge,
                          Point_2 & closest_interect_point,
                          bool & new_vertex, 
                          Vertex_handle & out_vertex) const;

  bool is_point_in_face (const Point_2 & p, 
                         const Ccb_halfedge_circulator & face, 
                         bool & found_edge,
                         bool & found_vertex, 
                         Halfedge_handle & out_edge) const;

  bool find_closest_intersection_in_face (const Point_2 & p,            
                                          Vertex_handle  v,     
                                          const Ccb_halfedge_circulator & face,
                                          Halfedge_handle & out_edge) const ;

  bool find_edge_to_flip (const Point_2 & p,            
                          Vertex_handle  v,
                          const Ccb_halfedge_circulator & face,
                          Halfedge_handle & out_edge) const ;

  bool find_real_intersection (const Point_2 & p,   
                               Vertex_handle  v,
                               Halfedge_handle e,
                               Point_2 & out_point) const;

  bool check_approximate_intersection (const Curve_2 & seg,  
                                       const Curve_2 & cv, 
                                       bool & intersect) const ;

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
  mutable Nearest_neighbor  nn;
  mutable bool        updated_nn;
  mutable Edge_list    flipped_edges;

#ifdef LANDMARKS_CLOCK_DEBUG
  mutable double clock_ff; //find face
  mutable double clock_for_nn_search ; 
  mutable double clock_for_walk ; 
  mutable double clock_nn_and_walk ; 
  mutable double clock_find_edge;
  mutable double clock_new_alg;
  mutable double clock_create_nn;
  mutable double clock_check_app;
  mutable double clock_is_point;

  mutable int entries_find_edge;
  mutable int entries_to_check_app;
  mutable int entries_to_find_face; 
  mutable int entries_is_point_in_face;
#endif

#ifdef TRAITS_CLOCK_DEBUG
  //count entries
  mutable int e_compare_distance;
  mutable int e_point_equal;
  mutable int e_curve_is_between_cw;
  mutable int e_point_in_x_range;
  mutable int e_curve_compare_y_at_x;
  mutable int e_nearest_intersection_to_left;
  mutable int e_nearest_intersection_to_right;
  mutable int e_curve_split;
  mutable int e_compare_xy;
  mutable int e_curves_compare_y_at_x_left;
  mutable int e_curves_compare_y_at_x_right;
  mutable int e_curves_compare_cw;
  //count clocks
  mutable double c_compare_distance;
  mutable double c_point_equal;
  mutable double c_curve_is_between_cw;
  mutable double c_point_in_x_range;
  mutable double c_curve_compare_y_at_x;
  mutable double c_nearest_intersection_to_left;
  mutable double c_nearest_intersection_to_right;
  mutable double c_curve_split;
  mutable double c_compare_xy;
  mutable double c_curves_compare_y_at_x_left;
  mutable double c_curves_compare_y_at_x_right;
  mutable double c_curves_compare_cw;
#endif

  const bool verbose; 
};

CGAL_END_NAMESPACE

#include <CGAL/Pm_landmarks_point_location.C>

#endif
